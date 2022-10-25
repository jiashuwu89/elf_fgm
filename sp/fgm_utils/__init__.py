import datetime
from distutils.log import Log
from typing import Literal
import numpy as np
import pandas as pd
from pyspedas.cotrans import cotrans_lib
from . import parameter
from .function import cross_time, Bplot, igrf, preprocess, error, postprocess, output, step0, step1

datestr = ""

def fgm_fsp_calib(
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

#    print(f"fgm data:{df['fgm_fgm']}")
#    print(f"att_data:{df['att_gei']}")
#    print(f"pos_data:{df['pos_gei']}")

    # df.set_index('time', inplace = True)
    # pprint(df.head())

    # iyear, idoy, ih, im, isec = cotrans_lib.get_time_parts(df['timestamp'])
    # print(f"year:{iyear}, doy:{idoy}, h:{ih}, m:{im}, sec:{isec}")

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
    # tstart = datetime.datetime(2022, 1, 12, 15, 45, 59)
    # xgsm = -2431.1245629621699
    # ygsm = 3822.9186030446831
    # zgsm = 5059.6970615621403
    # bxgsm, bygsm, bzgsm = get_igrf(tstart, xgsm, ygsm, zgsm)
    # print(bxgsm, bygsm, bzgsm)

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

    logger.info(f"Step 0 preprocess starts ... ")
    # check data sanity
    try:
        preprocess.funkyfgm_check(fgs_ful_fgm_0th_x, ctime, datestr)
    except (error.funkyFGMError, error.CrossTime1Error, error.funkyFGMError_len) as e:
        if parameter.makeplot == True:
            Bplot.B_ctime_plot(ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, datestr = datestr, title = "funkyFGM")
        logger.error(e.__str__())
        return [ [] for _ in range(16) ]
    
    # check repeated ctime
    if parameter.ctime_repeat_check == True:
        ctime_idx_repeat = preprocess.ctime_check(ctime)
        if ctime_idx_repeat is not None:
            ctime = np.delete(ctime, ctime_idx_repeat)
            fgs_ful_fgm_0th_x = np.delete(fgs_ful_fgm_0th_x, ctime_idx_repeat)
            fgs_ful_fgm_0th_y = np.delete(fgs_ful_fgm_0th_y, ctime_idx_repeat)
            fgs_ful_fgm_0th_z = np.delete(fgs_ful_fgm_0th_z, ctime_idx_repeat)
            fgs_igrf_gei_x = np.delete(fgs_igrf_gei_x, ctime_idx_repeat)
            fgs_igrf_gei_y = np.delete(fgs_igrf_gei_y, ctime_idx_repeat)
            fgs_igrf_gei_z = np.delete(fgs_igrf_gei_z, ctime_idx_repeat)
            att_gei_x = np.delete(att_gei_x, ctime_idx_repeat)
            att_gei_y = np.delete(att_gei_y, ctime_idx_repeat)
            att_gei_z = np.delete(att_gei_z, ctime_idx_repeat)
            logger.info("[PREPROCESS] repeat ctime found and delete!")

    """
        # 0. step 0, ctime calibration
    """
    logger.info(f"Step 0 ctime calibration starts ... ")
    [
        fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, 
        ctime_idx, ctime_idx_time, ctime_idx_flag, ctime_idx_timediff
        ] = step0.step0(
            ctime, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
            att_gei_x, att_gei_y, att_gei_z, datestr, logger
        )
        
    """
        # 1. step 1, B calibration
    """
    logger.info(f"Step 1 calibration starts ... ")
    [
        cross_times_calib, w_syn_d_calib, T_spins_d_calib,
        fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, 
        fgs_ful_gei_x, fgs_ful_gei_y, fgs_ful_gei_z,
        fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z,
        ] = step1.step1(
            ctime, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, 
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
            att_gei_x, att_gei_y, att_gei_z,
            datestr, logger, ctime_idx, ctime_idx_time, ctime_idx_flag, ctime_idx_timediff
        )

    if parameter.output == True:
        # full res
        fgs_res_dmxl_x = fgs_ful_dmxl_x - fgs_igrf_dmxl_x
        fgs_res_dmxl_y = fgs_ful_dmxl_y - fgs_igrf_dmxl_y
        fgs_res_dmxl_z = fgs_ful_dmxl_z - fgs_igrf_dmxl_z

        FGM_datetime = list(map(lambda ts: (df["time"][0].to_pydatetime() + 
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
        fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z] = cross_time.fsp_igrf(
            ctime, cross_times_calib, T_spins_d_calib, fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z
    )
    [
        fgs_fsp_ful_dmxl_x, fgs_fsp_ful_dmxl_y, fgs_fsp_ful_dmxl_z] = cross_time.fsp_ful(
            ctime, cross_times_calib, T_spins_d_calib, fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z
    )
    [
        fgs_fsp_ful_gei_x, fgs_fsp_ful_gei_y, fgs_fsp_ful_gei_z] = cross_time.fsp_ful(
            ctime, cross_times_calib, T_spins_d_calib, fgs_ful_gei_x, fgs_ful_gei_y, fgs_ful_gei_z
    )

    del_idx = np.where((fgs_fsp_ful_gei_x == 0) & (fgs_fsp_ful_gei_y == 0) & (fgs_fsp_ful_gei_z == 0))
    fgs_fsp_igrf_dmxl_x = np.delete(fgs_fsp_igrf_dmxl_x, del_idx)
    fgs_fsp_igrf_dmxl_y = np.delete(fgs_fsp_igrf_dmxl_y, del_idx)
    fgs_fsp_igrf_dmxl_z = np.delete(fgs_fsp_igrf_dmxl_z, del_idx)
    fgs_fsp_igrf_gei_x = np.delete(fgs_fsp_igrf_gei_x, del_idx)
    fgs_fsp_igrf_gei_y = np.delete(fgs_fsp_igrf_gei_y, del_idx)
    fgs_fsp_igrf_gei_z = np.delete(fgs_fsp_igrf_gei_z, del_idx)
    fgs_fsp_ful_dmxl_x = np.delete(fgs_fsp_ful_dmxl_x, del_idx)
    fgs_fsp_ful_dmxl_y = np.delete(fgs_fsp_ful_dmxl_y, del_idx)
    fgs_fsp_ful_dmxl_z = np.delete(fgs_fsp_ful_dmxl_z, del_idx)
    fgs_fsp_ful_gei_x = np.delete(fgs_fsp_ful_gei_x, del_idx)
    fgs_fsp_ful_gei_y = np.delete(fgs_fsp_ful_gei_y, del_idx)
    fgs_fsp_ful_gei_z = np.delete(fgs_fsp_ful_gei_z, del_idx)
    cross_times_calib = np.delete(cross_times_calib, del_idx)

    fgs_fsp_res_dmxl_x = fgs_fsp_ful_dmxl_x - fgs_fsp_igrf_dmxl_x
    fgs_fsp_res_dmxl_y = fgs_fsp_ful_dmxl_y - fgs_fsp_igrf_dmxl_y
    fgs_fsp_res_dmxl_z = fgs_fsp_ful_dmxl_z - fgs_fsp_igrf_dmxl_z

    fgs_fsp_res_gei_x = fgs_fsp_ful_gei_x - fgs_fsp_igrf_gei_x
    fgs_fsp_res_gei_y = fgs_fsp_ful_gei_y - fgs_fsp_igrf_gei_y
    fgs_fsp_res_gei_z = fgs_fsp_ful_gei_z - fgs_fsp_igrf_gei_z

    try:
        [
            cross_times_calib, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z
            ] = postprocess.fsp_spike_del(
            ctime, ctime_idx, ctime_idx_flag, ctime_idx_timediff, 
            cross_times_calib, w_syn_d_calib,
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, 
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y,fgs_fsp_igrf_dmxl_z, 
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
            logger
        )
    except error.fsp_spike_del_error as e:
        logger.error(e.__str__())
        return [ [] for _ in range(16) ]

    fgs_fsp_res_dmxl_trend_x = [0] * len(fgs_fsp_res_gei_x)
    fgs_fsp_res_dmxl_trend_y = [0] * len(fgs_fsp_res_gei_x)
    fgs_fsp_res_dmxl_trend_z = [0] * len(fgs_fsp_res_gei_x)

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
    FGM_timestamp = df["timestamp"][0] + cross_times_calib     

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
    ]
