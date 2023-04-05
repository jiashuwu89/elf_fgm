import datetime
from distutils.log import Log
from typing import Literal
import numpy as np
import pandas as pd
from pyspedas.cotrans import cotrans_lib
from .. import parameter
from . import cross_time, Bplot, igrf, preprocess, error, postprocess, output, step0, step1, detrend
from .coordinate import dmxl2gei, gei2obw, gei_obw_matrix
import traceback
import sys
from ..rotation_angle.floop_plot import Gthphi_f

datestr = ""

def fgm_fsp_calib_prepos(
    mission: Literal["ela", "elb"],
    starttime: datetime.datetime,
    endtime: datetime.datetime,
    fgm_cdfdata: pd.DataFrame,
    att_cdfdata: pd.DataFrame,
    pos_cdfdata: pd.DataFrame,
    logger: Log,
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

    # resample att and pos to fgm time resolution
    df["att_gei"] = preprocess.resample_data(
        att_cdfdata.index, att_cdfdata[f"{mission}_att_gei"], fgm_cdfdata.index
    )
    df["pos_gei"] = preprocess.resample_data(
        pos_cdfdata.index, pos_cdfdata[f"{mission}_pos_gei"], fgm_cdfdata.index
    )
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



def fgm_fsp_calib_ploop(
    ctime: list,
    ctimestamp: float,
    fgs_ful_fgm_0th_x: list, 
    fgs_ful_fgm_0th_y: list,
    fgs_ful_fgm_0th_z: list,
    fgs_igrf_gei_x: list,
    fgs_igrf_gei_y: list,
    fgs_igrf_gei_z: list,
    att_gei_x: list,
    att_gei_y: list,
    att_gei_z: list,
    pos_gei_x: list,
    pos_gei_y: list,
    pos_gei_z: list,
    logger: Log,
    ploop_i: int,
):

    """phase shift fgm x
    """
    if parameter.p_loop_xvalue[ploop_i] != 0:
        if parameter.p_loop_xvalue[ploop_i] > 0:
            fgs_ful_fgm_0th_x[parameter.p_loop_xvalue[ploop_i]:] = fgs_ful_fgm_0th_x[:-parameter.p_loop_xvalue[ploop_i]] # shift forward
            [
                ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z, pos_gei_x, pos_gei_y, pos_gei_z] = detrend.delete_data(
                range(parameter.p_loop_xvalue[ploop_i]), ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, att_gei_x, att_gei_y, att_gei_z,
                pos_gei_x, pos_gei_y, pos_gei_z)
        else:
            fgs_ful_fgm_0th_x[:parameter.p_loop_xvalue[ploop_i]] = fgs_ful_fgm_0th_x[-parameter.p_loop_xvalue[ploop_i]:] # shift backward
            [
                ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z, pos_gei_x, pos_gei_y, pos_gei_z] = detrend.delete_data(
                range(parameter.p_loop_xvalue[ploop_i],0), ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, att_gei_x, att_gei_y, att_gei_z,
                pos_gei_x, pos_gei_y, pos_gei_z)
    """phase shift fgm y
    """
    if parameter.p_loop_yvalue[ploop_i] != 0:
        if parameter.p_loop_yvalue[ploop_i] > 0:
            fgs_ful_fgm_0th_y[parameter.p_loop_yvalue[ploop_i]:] = fgs_ful_fgm_0th_y[:-parameter.p_loop_yvalue[ploop_i]] # shift forward
            [
                ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z, pos_gei_x, pos_gei_y, pos_gei_z] = detrend.delete_data(
                range(parameter.p_loop_yvalue[ploop_i]), ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, att_gei_x, att_gei_y, att_gei_z,
                pos_gei_x, pos_gei_y, pos_gei_z)
        else:
            fgs_ful_fgm_0th_x[:parameter.p_loop_yvalue[ploop_i]] = fgs_ful_fgm_0th_x[-parameter.p_loop_yvalue[ploop_i]:] # shift backward
            [
                ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z, pos_gei_x, pos_gei_y, pos_gei_z] = detrend.delete_data(
                range(parameter.p_loop_yvalue[ploop_i],0), ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, att_gei_x, att_gei_y, att_gei_z,
                pos_gei_x, pos_gei_y, pos_gei_z)

    [
        ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
        fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
        att_gei_x, att_gei_y, att_gei_z, pos_gei_x, pos_gei_y, pos_gei_z] = detrend.delete_data(
        range(ploop_i), ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
        fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, att_gei_x, att_gei_y, att_gei_z,
        pos_gei_x, pos_gei_y, pos_gei_z)


    logger.info(f"Step 0 preprocess starts ... ")
    # check data sanity
    
    if parameter.funkyfgm == True:
        try:
            preprocess.funkyfgm_check(fgs_ful_fgm_0th_x, ctime, datestr)
        except (error.funkyFGMError, error.CrossTime1Error, error.funkyFGMError_len) as e:
            if parameter.makeplot == True:
                Bplot.B_ctime_plot(ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, datestr = datestr, title = "funkyFGM")
            logger.error(e.__str__())
            return [ [] for _ in range(16) ]
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
            att_gei_x, att_gei_y, att_gei_z, pos_gei_x, pos_gei_y, pos_gei_z] = detrend.delete_data(
            del_time_idx, ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, att_gei_x, att_gei_y, att_gei_z,
            pos_gei_x, pos_gei_y, pos_gei_z)


    # check repeated ctime
    if parameter.ctime_repeat_check == True:
        ctime_idx_repeat = preprocess.ctime_check(ctime)
        if ctime_idx_repeat is not None:
            [
                ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z,
                pos_gei_x, pos_gei_y, pos_gei_z] = detrend.delete_data(
                    ctime_idx_repeat, ctime, 
                    fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                    att_gei_x, att_gei_y, att_gei_z,
                    pos_gei_x, pos_gei_y, pos_gei_z)
            logger.info("[PREPROCESS] repeat ctime found and delete!")

    """
        # 0. step 0, ctime calibration
    """
    logger.info(f"Step 0 ctime calibration starts ... ")
    try:
        [
            fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, 
            ctime_idx, ctime_idx_time, ctime_idx_flag, ctime_idx_timediff,
            B_parameter0,
            ] = step0.step0(
                ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z, datestr, logger, parameter.f,
            )
    except error.CrossTime1Error as e:
        logger.error(e.__str__())
        return [ [] for _ in range(16) ]
        logger.error(traceback.format_exception(*sys.exc_info()))
    except Exception as e:
        logger.error(f"❌ step 0 other error. Stop processing.")
        logger.error('\n'.join(traceback.format_exception(*sys.exc_info())))
        print('\n'.join(traceback.format_exception(*sys.exc_info())))
        return [ [] for _ in range(16) ]

    """
        # 1. step 1, B calibration
    """
    logger.info(f"Step 1 calibration starts ... ")
    try:
        [
            cross_times_calib, w_syn_d_calib, T_spins_d_calib, DMXL_2_GEI_fsp,
            fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, 
            fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z,
            B_parameter1,
            ] = step1.step1(
                #ctime, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, 
                ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z,
                datestr, logger, ctime_idx, ctime_idx_time, ctime_idx_flag, ctime_idx_timediff, parameter.f, ploop_i = ploop_i
            )
    except:
        logger.error(f"❌ step 1 other error. Stop processing.")
        logger.error('\n'.join(traceback.format_exception(*sys.exc_info())))
        print('\n'.join(traceback.format_exception(*sys.exc_info())))
        return [ [] for _ in range(16) ]
        
    if parameter.output == True:
        # full res
        fgs_res_dmxl_x = fgs_ful_dmxl_x - fgs_igrf_dmxl_x
        fgs_res_dmxl_y = fgs_ful_dmxl_y - fgs_igrf_dmxl_y
        fgs_res_dmxl_z = fgs_ful_dmxl_z - fgs_igrf_dmxl_z

        FGM_datetime = list(map(lambda ts: (ctimestamp + 
                        datetime.timedelta(seconds=ts)).strftime('%Y-%m-%d/%H:%M:%S.%f'), ctime))
        output.output_txt(FGM_datetime, fgs_res_dmxl_x, fgs_res_dmxl_y, fgs_res_dmxl_z, title='ela_fgs_res_dmxl')  

        print(f"std res_x: {np.std(fgs_res_dmxl_x)}") 
        print(f"std res_y: {np.std(fgs_res_dmxl_y)}") 
        print(f"std res_z: {np.std(fgs_res_dmxl_z)}")  
        #breakpoint()  

    """
        # 2 : step 2 fgs fsp resolution
    """
    logger.info(f"Step 2 fsp data starts ... ")
    [
        fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z] = cross_time.fsp_igrf(
            ctime, cross_times_calib, T_spins_d_calib, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z
    )
    [
        fgs_fsp_ful_dmxl_x, fgs_fsp_ful_dmxl_y, fgs_fsp_ful_dmxl_z] = cross_time.fsp_ful(
            ctime, cross_times_calib, T_spins_d_calib, fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z
    )
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
        return [ [] for _ in range(16) ]
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
        Bplot.B_ctime_plot(
            cross_times_calib, fgs_fsp_res_dmxl_x, 
            fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, title="res_dmxl_fsp", scatter = True, 
            ctime_idx_time = ctime_idx_time, datestr = datestr, ctime_idx_flag = ctime_idx_flag
        )
        Bplot.B_ctime_plot(
            cross_times_calib, fgs_fsp_res_gei_x, 
            fgs_fsp_res_gei_y, fgs_fsp_res_gei_z, title="res_gei_fsp", scatter = True, 
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
            Bplot.B_ctime_plot(cross_times_calib, fgs_fsp_res_obw_x, fgs_fsp_res_obw_y, fgs_fsp_res_obw_z, title="res_obw_fsp", 
            ctime_idx_time = ctime_idx_time, datestr = datestr, ctime_idx_flag = ctime_idx_flag)

    # convert Bparamters to angles in degree
    G1, G2, G3, th1, th2, th3, ph1, ph2, ph3 = Gthphi_f(
        B_parameter1[0], B_parameter1[1], B_parameter1[2], 
        B_parameter1[4], B_parameter1[5], B_parameter1[6], 
        B_parameter1[8], B_parameter1[9], B_parameter1[10])
    B_parameter1_ang = [G1, G2, G3, th1, th2, th3, ph1, ph2, ph3, B_parameter1[3]/G1, B_parameter1[7]/G2, B_parameter1[11]/G3]

    return [fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, B_parameter1, B_parameter1_ang]