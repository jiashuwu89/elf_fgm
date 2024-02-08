import datetime
from logging import Logger
from typing import Literal, List
import numpy as np
import pandas as pd
import traceback
import sys
import requests
from pyspedas.cotrans import cotrans_lib
from . import parameter
from .function import cross_time, Bplot, igrf, preprocess, error, postprocess, output, step0, step1, detrend, wfit
from .function.coordinate import dmxl2gei, gei2obw, gei_obw_matrix
from .function.attitude import att_rot

datestr = ""

def fgm_fsp_calib_prepos(
    mission: Literal["ela", "elb"],
    starttime: datetime.datetime,
    endtime: datetime.datetime,
    fgm_cdfdata: pd.DataFrame,
    att_cdfdata: pd.DataFrame,
    pos_cdfdata: pd.DataFrame,
    ):
    """
    Note that starttime, endtime refer to the start and end of the science zone collection
    """
    # initial points exclude
    if parameter.init_secs != 0:
        starttime = starttime + datetime.timedelta(seconds = parameter.init_secs)

    df = pd.DataFrame()

    # read fgm cdf and clip
    fgm_cdfdata.set_index(f"{mission}_fgs_time", inplace=True)
    fgm_cdfdata = preprocess.clip_cdfdata(fgm_cdfdata, starttime, endtime)
    try:
        # resample att and pos to fgm time resolution
        df["att_gei"] = preprocess.resample_data(
            att_cdfdata.index, att_cdfdata[f"{mission}_att_gei"], fgm_cdfdata.index
        )
        df["pos_gei"] = preprocess.resample_data(
            pos_cdfdata.index, pos_cdfdata[f"{mission}_pos_gei"], fgm_cdfdata.index
        )
    except:
        raise error.preproc_resample_error()
    
    df["fgm_fgm"] = pd.Series(fgm_cdfdata[f"{mission}_fgs"].tolist())
    df["time"] = fgm_cdfdata.index
    df["timestamp"] = df["time"].apply(lambda ts: pd.Timestamp(ts).timestamp())

    # coordinate transformation of pos: gei -> gse -> gsm
    df["pos_gse"] = pd.Series(
        cotrans_lib.subgei2gse(
            df["timestamp"].tolist(), df["pos_gei"].tolist()
        ).tolist()
    )
    df["pos_gsm"] = pd.Series(
        cotrans_lib.subgse2gsm(
            df["timestamp"].tolist(), df["pos_gse"].tolist()
        ).tolist()
    )

    # call igrf b in gsm
    df["igrf_gsm"] = [
        igrf.get_igrf(
            df["time"][i], df["pos_gsm"][i][0], df["pos_gsm"][i][1], df["pos_gsm"][i][2]
        )
        for i in range(len(df["timestamp"]))
    ]

    # coordinate transformation of B: gsm -> gse -> gei
    df["igrf_gse"] = pd.Series(
        cotrans_lib.subgsm2gse(
            df["timestamp"].tolist(), df["igrf_gsm"].tolist()
        ).tolist()
    )
    df["igrf_gei"] = pd.Series(
        cotrans_lib.subgse2gei(
            df["timestamp"].tolist(), df["igrf_gse"].tolist()
        ).tolist()
    )

    global datestr
    datestr = df["time"][0].to_pydatetime().strftime('%Y-%m-%d/%H:%M:%S').replace("-","").replace("/","_").replace(":","")
    ############################################
    #   suyash code begin
    ############################################
    df["ctime"] = df["timestamp"] - df["timestamp"][0]
    ctime = np.array(df["ctime"])
    # get fgm data in fgm coordinate
    fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z = np.array(list(zip(*df["fgm_fgm"])))
    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z = np.array(list(zip(*df["igrf_gei"])))
    att_gei_x, att_gei_y, att_gei_z = np.array(list(zip(*df["att_gei"])))
    pos_gei_x, pos_gei_y, pos_gei_z = np.array(list(zip(*df["pos_gei"])))
    
    ctimestamp = df["timestamp"][0]

    return [ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
        fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z,
        att_gei_x, att_gei_y, att_gei_z,
        pos_gei_x, pos_gei_y, pos_gei_z]


def fgm_fsp_calib_prepos_wrapper(
        mission: Literal["ela", "elb"], 
        start_time: List[datetime.datetime], 
        end_time: List[datetime.datetime],
        f_all: List[float],
        logger: Logger) -> List[np.ndarray]:
    """Wrapper for fgm_fsp_calib_prepos_wrapper.
    Read data from all science zones and prepare for processing
    
    Args:
        mission: 'ela' or 'elb'
        start_time: List of start time of sci zones
        end_time: List of end time of sci zones
        f_all: List of rotation angle
    
    Returns:
        List of numpy array
    """
    for i in range(len(start_time)):
        sta_datestr = start_time[i].strftime("%Y%m%d")
        logger.info(f"▶️ Received {mission} collection from {start_time[i]} to {end_time[i]}")
        
        if parameter.state03 == True:
            sta_cdfpath = f"fgm_utils/test/{mission}_l1_state_defn_{sta_datestr}_v03.cdf"
        else:
            sta_cdfpath = f"fgm_utils/test/{mission}_l1_state_defn_{sta_datestr}_v02.cdf"
        fgm_cdfpath = f"fgm_utils/test/{mission}_l1_fgs_{sta_datestr}_v01.cdf" 

        if parameter.download_data == True:
            try:
                if parameter.state03 == True:
                    sta_url = f"{parameter.elfin_url}{mission}/l1/state/defn/{start_time[i].year}/{mission}_l1_state_defn_{sta_datestr}_v03.cdf"
                else:    
                    sta_url = f"{parameter.elfin_url}{mission}/l1/state/defn/{start_time[i].year}/{mission}_l1_state_defn_{sta_datestr}_v02.cdf"
                res = requests.get(sta_url)
                logger.info(f"Download file {sta_url} sucessful!")
                with open(sta_cdfpath, 'wb') as f:
                    f.write(res.content)
            except:
                logger.error(f"download error:{sta_url}")
                raise error.preproc_download_error()

            try:
                fgm_url = f"{parameter.elfin_url}{mission}/l1/fgm/survey/{start_time[i].year}/{mission}_l1_fgs_{sta_datestr}_v01.cdf"
                res = requests.get(fgm_url)
                logger.info(f"Download file {fgm_url} sucessful!")
                with open(fgm_cdfpath, 'wb') as f:
                    f.write(res.content)
            except:
                logger.error(f"download error:{fgm_url}") 
                raise error.preproc_download_error()
                

        fgm_cdfdata = pd.DataFrame(preprocess.get_cdf(fgm_cdfpath, vars=[f"{mission}_fgs_time", f"{mission}_fgs"]))
        logger.info(f"Sucessfully read cdf for {mission} from {start_time[i]} to {end_time[i]} from file {fgm_cdfpath}")
        att_cdfdata, pos_cdfdata = preprocess.get_relevant_state_data(sta_cdfpath, mission, start_time[i], end_time[i])
        logger.info(f"Sucessfully read state cdf for {mission} from {start_time[i]} to {end_time[i]} from file {sta_cdfpath}")      

        if parameter.att_rot == True:
            att_cdfdata = att_rot(att_cdfdata, parameter.att_rot_ang, parameter.att_rot_axis)    

        [
            ctime_0, ctimestamp_0,
            fgs_ful_fgm_0th_x_0, fgs_ful_fgm_0th_y_0, fgs_ful_fgm_0th_z_0, 
            fgs_igrf_gei_x_0, fgs_igrf_gei_y_0, fgs_igrf_gei_z_0,
            att_gei_x_0, att_gei_y_0, att_gei_z_0,
            pos_gei_x_0, pos_gei_y_0, pos_gei_z_0] = fgm_fsp_calib_prepos(
                mission, start_time[i], end_time[i], fgm_cdfdata, att_cdfdata, pos_cdfdata)      
        
        f_all_arry_0 = [f_all[i]] * len(fgs_ful_fgm_0th_x_0) if f_all is not None else [parameter.f] * len(fgs_ful_fgm_0th_x_0) 
        
        if i == 0: # first collection 
            [
                ctime, ctimestamp, 
                fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z,
                pos_gei_x, pos_gei_y, pos_gei_z, f_all_arry] = [
                    ctime_0, ctimestamp_0, 
                    fgs_ful_fgm_0th_x_0, fgs_ful_fgm_0th_y_0, fgs_ful_fgm_0th_z_0, 
                    fgs_igrf_gei_x_0, fgs_igrf_gei_y_0, fgs_igrf_gei_z_0, 
                    att_gei_x_0, att_gei_y_0, att_gei_z_0, 
                    pos_gei_x_0, pos_gei_y_0, pos_gei_z_0, f_all_arry_0]
            clip_start_idx = [0]
            clip_end_idx = [len(ctime)-1]

        else: # multiple collection
            ctime = np.concatenate((ctime, ctime_0 + ctime[-1] + 60))
            fgs_ful_fgm_0th_x = np.concatenate((fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_x_0))
            fgs_ful_fgm_0th_y = np.concatenate((fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_y_0))
            fgs_ful_fgm_0th_z = np.concatenate((fgs_ful_fgm_0th_z, fgs_ful_fgm_0th_z_0))
            fgs_igrf_gei_x = np.concatenate((fgs_igrf_gei_x, fgs_igrf_gei_x_0))
            fgs_igrf_gei_y = np.concatenate((fgs_igrf_gei_y, fgs_igrf_gei_y_0))
            fgs_igrf_gei_z = np.concatenate((fgs_igrf_gei_z, fgs_igrf_gei_z_0))
            att_gei_x = np.concatenate((att_gei_x, att_gei_x_0))
            att_gei_y = np.concatenate((att_gei_y, att_gei_y_0))
            att_gei_z = np.concatenate((att_gei_z, att_gei_z_0))
            pos_gei_x = np.concatenate((pos_gei_x, pos_gei_x_0))
            pos_gei_y = np.concatenate((pos_gei_y, pos_gei_y_0))
            pos_gei_z = np.concatenate((pos_gei_z, pos_gei_z_0))
            f_all_arry = np.concatenate([f_all_arry, f_all_arry_0])
            clip_start_idx.append(clip_end_idx[-1]+1) # start index of each sci zone
            clip_end_idx.append(len(ctime)-1) # end index of each sci zone
        
        return [ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z, 
                pos_gei_x, pos_gei_y, pos_gei_z,
                f_all_arry, clip_start_idx, clip_end_idx]


def fgm_fsp_calib(
    ctime: List[float],
    ctimestamp: float,
    f_all: List[float],
    fgs_ful_fgm_0th_x: List[float], 
    fgs_ful_fgm_0th_y: List[float],
    fgs_ful_fgm_0th_z: List[float],
    fgs_igrf_gei_x: List[float],
    fgs_igrf_gei_y: List[float],
    fgs_igrf_gei_z: List[float],
    att_gei_x: List[float],
    att_gei_y: List[float],
    att_gei_z: List[float],
    pos_gei_x: List[float],
    pos_gei_y: List[float],
    pos_gei_z: List[float],
    logger: Logger,
    att_loop_idx = None,
):

    if parameter.funkyfgm == True:
        try:
            preprocess.funkyfgm_check(fgs_ful_fgm_0th_x, ctime, datestr)
        except (error.funkyFGMError, error.CrossTime1Error, error.funkyFGMError_len) as e:
            if parameter.makeplot == True:
                Bplot.B_ctime_plot(ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, datestr = datestr, title = "funkyFGM")
            logger.error(e.__str__())
            return [ [] for _ in range(17) ]
        except Exception as e:
            logger.error(f"❌ funky fgm check failed. ")
            logger.error('\n'.join(traceback.format_exception(*sys.exc_info())))
            print('\n'.join(traceback.format_exception(*sys.exc_info())))


    if parameter.del_time == True:
        del_time_idx = []
        for i in range(len(parameter.del_time_idxend)):
            for j in range(parameter.del_time_idxstart[i]*10, parameter.del_time_idxend[i]*10):
                del_time_idx.append(j)
        [
            ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
            att_gei_x, att_gei_y, att_gei_z, pos_gei_x, pos_gei_y, pos_gei_z, f_all] = detrend.delete_data(
            del_time_idx, ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, att_gei_x, att_gei_y, att_gei_z,
            pos_gei_x, pos_gei_y, pos_gei_z, f_all)
        logger.info("[PREPROCESS] delete user specify data!")
    

    # check repeated ctime
    if parameter.ctime_repeat_check == True:
        ctime_idx_repeat = preprocess.ctime_check(ctime)
        if ctime_idx_repeat is not None:
            try:
                [
                    ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                    att_gei_x, att_gei_y, att_gei_z, f_all] = detrend.delete_data(
                        ctime_idx_repeat, ctime, 
                        fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                        fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                        att_gei_x, att_gei_y, att_gei_z, f_all)
            except:
                breakpoint()
            logger.info("[PREPROCESS] repeat ctime found and delete!")

    """
        # 0. step 0, ctime calibration
    """
    logger.info(f"Step 0 ctime calibration starts ... ")
    try:
        [
            fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, 
            ctime_idx, ctime_idx_time, ctime_idx_flag, ctime_idx_timediff, B_parameter0
            ] = step0.step0(
                ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z, datestr, logger, f_all,
            )
    except error.CrossTime1Error as e:
        logger.error(e.__str__())
        return [ [] for _ in range(17) ]
        logger.error(traceback.format_exception(*sys.exc_info()))
    except Exception as e:
        logger.error(f"❌ step 0 other error. Stop processing.")
        logger.error('\n'.join(traceback.format_exception(*sys.exc_info())))
        print('\n'.join(traceback.format_exception(*sys.exc_info())))
        return [ [] for _ in range(17) ]
       
    """
        # 1. step 1, B calibration
    """
    logger.info(f"Step 1 calibration starts ... ")
    try:
        if parameter.wfit_run == True:
            [
                cross_times_calib, w_syn_d_calib, T_spins_d_calib, DMXL_2_GEI_fsp,
                fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, 
                fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z,
                B_parameter1,
                ] = wfit.step1(
                    #ctime, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, 
                    ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                    att_gei_x, att_gei_y, att_gei_z,
                    datestr, logger, ctime_idx, ctime_idx_time, ctime_idx_flag, ctime_idx_timediff, f_all,
                )
        else:
            [
                cross_times_calib, w_syn_d_calib, T_spins_d_calib, DMXL_2_GEI_fsp,
                fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, 
                fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z,
                B_parameter1,
                ] = step1.step1(
                    #ctime, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, 
                    ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                    att_gei_x, att_gei_y, att_gei_z,
                    datestr, logger, ctime_idx, ctime_idx_time, ctime_idx_flag, ctime_idx_timediff, f_all,
                )

    except:
        logger.error(f"❌ step 1 other error. Stop processing.")
        logger.error('\n'.join(traceback.format_exception(*sys.exc_info())))
        print('\n'.join(traceback.format_exception(*sys.exc_info())))
        return [ [] for _ in range(17) ]
        
    if parameter.output == True:
        # full res
        fgs_res_dmxl_x = fgs_ful_dmxl_x - fgs_igrf_dmxl_x
        fgs_res_dmxl_y = fgs_ful_dmxl_y - fgs_igrf_dmxl_y
        fgs_res_dmxl_z = fgs_ful_dmxl_z - fgs_igrf_dmxl_z

        FGM_datetime = list(map(lambda ts: (datetime.datetime.fromtimestamp(ctimestamp, tz=datetime.timezone.utc) + 
                        datetime.timedelta(seconds=ts)).strftime('%Y-%m-%d/%H:%M:%S.%f'), ctime))
        output.output_txt(
            FGM_datetime, 
            [fgs_res_dmxl_x, fgs_res_dmxl_y, fgs_res_dmxl_z], 
            ['Timestamp','fsp_gei_x','fsp_gei_y','fsp_gei_z'], 
            title='ela_fgs_res_dmxl')  

        print(f"std res_x: {np.std(fgs_res_dmxl_x)}") 
        print(f"std res_y: {np.std(fgs_res_dmxl_y)}") 
        print(f"std res_z: {np.std(fgs_res_dmxl_z)}")  
        #breakpoint()  

    """
        # 2 : step 2 fgs fsp resolution
    """
    logger.info(f"Step 2 fsp data starts ... ")
    try:
        [
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z] = cross_time.fsp_igrf(
                ctime, cross_times_calib, T_spins_d_calib, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z
        )
        [
            fgs_fsp_ful_dmxl_x, fgs_fsp_ful_dmxl_y, fgs_fsp_ful_dmxl_z] = cross_time.fsp_ful(
                ctime, cross_times_calib, T_spins_d_calib, fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z
        )
    except error.postproc_fgs_igrf as e:
        logger.error(e.__str__())
        return [ [] for _ in range(17) ]
    except Exception as e:
        logger.error('\n'.join(traceback.format_exception(*sys.exc_info())))
        print('\n'.join(traceback.format_exception(*sys.exc_info())))
        
    #[
    #    fgs_fsp_ful_dmxl_x, fgs_fsp_ful_dmxl_y, fgs_fsp_ful_dmxl_z] = cross_time.fsp_igrf(
    #        ctime, cross_times_calib, T_spins_d_calib, fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z
    #)

    del_idx = np.where((fgs_fsp_ful_dmxl_x == 0) & (fgs_fsp_ful_dmxl_y == 0) & (fgs_fsp_ful_dmxl_z == 0))
    [
        cross_times_calib, DMXL_2_GEI_fsp,
        fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z, 
        fgs_fsp_ful_dmxl_x, fgs_fsp_ful_dmxl_y, fgs_fsp_ful_dmxl_z] = detrend.delete_data(
            del_idx, cross_times_calib, DMXL_2_GEI_fsp,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z, 
            fgs_fsp_ful_dmxl_x, fgs_fsp_ful_dmxl_y, fgs_fsp_ful_dmxl_z,
        )

    fgs_fsp_res_dmxl_x = fgs_fsp_ful_dmxl_x - fgs_fsp_igrf_dmxl_x
    fgs_fsp_res_dmxl_y = fgs_fsp_ful_dmxl_y - fgs_fsp_igrf_dmxl_y
    fgs_fsp_res_dmxl_z = fgs_fsp_ful_dmxl_z - fgs_fsp_igrf_dmxl_z

    """
        # 3: step 4 detrend
        this step is swap with step 4, because in postporcess, ctime_idx needs to compare with average. if big trend exists, the spike is not prominent
    """
    if parameter.fsp_detrend == False:
        fgs_fsp_res_dmxl_trend_x = [0] * len(fgs_fsp_res_dmxl_x)
        fgs_fsp_res_dmxl_trend_y = [0] * len(fgs_fsp_res_dmxl_x)
        fgs_fsp_res_dmxl_trend_z = [0] * len(fgs_fsp_res_dmxl_x)

    if parameter.fsp_detrend == True:
        [
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z] = detrend.detrend_quad(
                cross_times_calib, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z)

        if parameter.makeplot == True:
            Bplot.B_ctime_plot(
                cross_times_calib, [fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_trend_x], 
                [fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_trend_y], [fgs_fsp_res_dmxl_z, fgs_fsp_res_dmxl_trend_z], title="res_dmxl_fsp_trend", scatter = True, 
                datestr = datestr, ctime_idx_flag = ctime_idx_flag
            )
        # detrend res dmxl
        fgs_fsp_res_dmxl_x = fgs_fsp_res_dmxl_x - fgs_fsp_res_dmxl_trend_x
        fgs_fsp_res_dmxl_y = fgs_fsp_res_dmxl_y - fgs_fsp_res_dmxl_trend_y
        fgs_fsp_res_dmxl_z = fgs_fsp_res_dmxl_z - fgs_fsp_res_dmxl_trend_z

        # get ful dmxl
        fgs_fsp_ful_dmxl_x = fgs_fsp_res_dmxl_x + fgs_fsp_igrf_dmxl_x
        fgs_fsp_ful_dmxl_y = fgs_fsp_res_dmxl_y + fgs_fsp_igrf_dmxl_y
        fgs_fsp_ful_dmxl_z = fgs_fsp_res_dmxl_z + fgs_fsp_igrf_dmxl_z
    else:
        # get ful dmxl
        fgs_fsp_ful_dmxl_x = fgs_fsp_res_dmxl_x + fgs_fsp_igrf_dmxl_x
        fgs_fsp_ful_dmxl_y = fgs_fsp_res_dmxl_y + fgs_fsp_igrf_dmxl_y
        fgs_fsp_ful_dmxl_z = fgs_fsp_res_dmxl_z + fgs_fsp_igrf_dmxl_z
    
    # transform ful dmxl to ful gei
    [fgs_fsp_ful_gei_x, fgs_fsp_ful_gei_y, fgs_fsp_ful_gei_z] = dmxl2gei(
        fgs_fsp_ful_dmxl_x, fgs_fsp_ful_dmxl_y, fgs_fsp_ful_dmxl_z, DMXL_2_GEI_fsp)
    
    # get igrf gei
    [
        fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z] = cross_time.fsp_igrf(
        ctime, cross_times_calib, T_spins_d_calib, fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z)

    # get res gei
    fgs_fsp_res_gei_x = fgs_fsp_ful_gei_x - fgs_fsp_igrf_gei_x
    fgs_fsp_res_gei_y = fgs_fsp_ful_gei_y - fgs_fsp_igrf_gei_y 
    fgs_fsp_res_gei_z = fgs_fsp_ful_gei_z - fgs_fsp_igrf_gei_z 
    
    """
        # 4 : step 3 delete spike and rogue points
    """
    try:
        [
            cross_times_calib, DMXL_2_GEI_fsp, 
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
            ] = postprocess.fsp_spike_del(
            ctime, ctime_idx, ctime_idx_flag, ctime_idx_timediff, 
            cross_times_calib, w_syn_d_calib, DMXL_2_GEI_fsp,
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, 
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y,fgs_fsp_igrf_dmxl_z, 
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
            logger,
        )
    except error.fsp_spike_del_error as e:
        logger.error(e.__str__())
        return [ [] for _ in range(17) ]
    except Exception as e:
        logger.error('\n'.join(traceback.format_exception(*sys.exc_info())))
        print('\n'.join(traceback.format_exception(*sys.exc_info())))

    """
        # 5 : step 5 final data and plot
    """
    if parameter.makeplot == True:
        #Bplot.B_ctime_plot(
        #    cross_times_calib, fgs_fsp_res_dmxl_x, 
        #    fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, title="res_dmxl_fsp", scatter = True, 
        #    ctime_idx_time = ctime[ctime_idx[0]], datestr = datestr, xlimt = [ctime[ctime_idx[0]]-10, ctime[ctime_idx[0]]+10]
        #)
        if parameter.att_loop == True:
            Bplot.B_ctime_plot(
                cross_times_calib, fgs_fsp_res_dmxl_x, 
                fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, title=f"res_dmxl_fsp{att_loop_idx}", scatter = True, 
                ctime_idx_time = ctime_idx_time, datestr = datestr, ctime_idx_flag = ctime_idx_flag,
            )
            Bplot.B_ctime_plot(
                cross_times_calib, fgs_fsp_res_gei_x, 
                fgs_fsp_res_gei_y, fgs_fsp_res_gei_z, title=f"res_gei_fsp{att_loop_idx}", scatter = True, 
                ctime_idx_time = ctime_idx_time, datestr = datestr, ctime_idx_flag = ctime_idx_flag
            )
        else:
            Bplot.B_ctime_plot(
                cross_times_calib, fgs_fsp_res_dmxl_x, 
                fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, title=f"res_dmxl_fsp", scatter = True, 
                ctime_idx_time = ctime_idx_time, datestr = datestr, ctime_idx_flag = ctime_idx_flag,
            )
            Bplot.B_ctime_plot(
                cross_times_calib, fgs_fsp_res_gei_x, 
                fgs_fsp_res_gei_y, fgs_fsp_res_gei_z, title=f"res_gei_fsp", scatter = True, 
                ctime_idx_time = ctime_idx_time, datestr = datestr, ctime_idx_flag = ctime_idx_flag
            )
    
    #FGM_datetime = list(map(lambda ts: (df["time"][0].to_pydatetime() + 
    #                           datetime.timedelta(seconds=ts)).strftime('%Y-%m-%d/%H:%M:%S'), cross_times_calib))
    #breakpoint()
    FGM_timestamp = ctimestamp + cross_times_calib     
    
    if parameter.gei2obw == True:
        [pos_fsp_gei_x, pos_fsp_gei_y, pos_fsp_gei_z] = cross_time.fsp_igrf(ctime, cross_times_calib, T_spins_d_calib, pos_gei_x, pos_gei_y, pos_gei_z)
        [GEI_2_OBW, OBW_2_GEI] = gei_obw_matrix(fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z, pos_fsp_gei_x, pos_fsp_gei_y, pos_fsp_gei_z)
        [fgs_fsp_res_obw_x, fgs_fsp_res_obw_y, fgs_fsp_res_obw_z] = gei2obw(fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z, GEI_2_OBW)
        if parameter.makeplot == True:
            if att_loop_idx is not None:
                Bplot.B_ctime_plot(cross_times_calib, fgs_fsp_res_obw_x, fgs_fsp_res_obw_y, fgs_fsp_res_obw_z, title=f"res_obw_fsp{att_loop_idx}", 
                    ctime_idx_time = ctime_idx_time, datestr = datestr, ctime_idx_flag = ctime_idx_flag)
            else:
                Bplot.B_ctime_plot(cross_times_calib, fgs_fsp_res_obw_x, fgs_fsp_res_obw_y, fgs_fsp_res_obw_z, title=f"res_obw_fsp", 
                    ctime_idx_time = ctime_idx_time, datestr = datestr, ctime_idx_flag = ctime_idx_flag)

    
    return [
        FGM_timestamp,
        fgs_fsp_res_dmxl_x,
        fgs_fsp_res_dmxl_y,
        fgs_fsp_res_dmxl_z,
        fgs_fsp_igrf_dmxl_x,
        fgs_fsp_igrf_dmxl_y,
        fgs_fsp_igrf_dmxl_z,
        fgs_fsp_res_dmxl_trend_x,
        fgs_fsp_res_dmxl_trend_y,
        fgs_fsp_res_dmxl_trend_z,
        fgs_fsp_res_gei_x,
        fgs_fsp_res_gei_y,
        fgs_fsp_res_gei_z,
        fgs_fsp_igrf_gei_x,
        fgs_fsp_igrf_gei_y,
        fgs_fsp_igrf_gei_z,
        B_parameter1,
    ]