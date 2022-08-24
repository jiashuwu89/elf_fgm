import datetime
from distutils.log import Log
from typing import Literal
import logging
import numpy as np
import pandas as pd
from pyspedas.cotrans import cotrans_lib
from . import parameter
from .function import calibration, coordinate, cross_time, Bplot, detrend, igrf, preprocess, error
from .function import output
from scipy.optimize import curve_fit

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
    B_S1_corr, B_S2_corr, B_S3_corr = np.array(list(zip(*df["fgm_fgm"])))
    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z = np.array(list(zip(*df["igrf_gei"])))
    att_gei_x, att_gei_y, att_gei_z = np.array(list(zip(*df["att_gei"])))


    # check data sanity
    try:
        preprocess.funkyfgm_check(B_S1_corr, ctime, datestr)
    except (error.funkyFGMError, error.CrossTime1Error) as e:
        logger.error(e.__str__())
        return [ [] for _ in range(16) ]
    """
        timestamp correction; TODO: automatic detection
    """
    ctime, ctime_idx = calibration.ctime_calib(ctime, B_x = B_S1_corr, B_y = B_S2_corr, B_z = B_S3_corr)
    
    if len(ctime_idx) > 0  and parameter.makeplot == True:
        #Bplot.ctimediff_plot(ctime, ctime_idx, datestr = datestr)
        Bplot.ctimediff_plot(ctime, ctime_idx, ctime_idx_zoom=ctime_idx[1])
    breakpoint()

    if parameter.del_spike_10hz == True:
        """
        # TEST: delete rogue point and interpolate 
        ctime = np.delete(ctime, [1363, 1364])
        B_S1_corr = np.delete(B_S1_corr, [1363, 1364])
        B_S2_corr = np.delete(B_S2_corr, [1363, 1364])
        B_S3_corr = np.delete(B_S3_corr, [1363, 1364])
        fgs_igrf_gei_x = np.delete(fgs_igrf_gei_x, [1363, 1364])
        fgs_igrf_gei_y = np.delete(fgs_igrf_gei_y, [1363, 1364])
        fgs_igrf_gei_z = np.delete(fgs_igrf_gei_z, [1363, 1364])
        att_gei_x = np.delete(att_gei_x, [1363, 1364])
        att_gei_y = np.delete(att_gei_y, [1363, 1364])
        att_gei_z = np.delete(att_gei_z, [1363, 1364])  
        """

        """
        # TEST: linear fit 
        fit_opt, fit_covar = curve_fit(
            calibration.linear_fit, ctime[1362:1366], B_S1_corr[1362:1366]
        )
        B_S1_corr_del = calibration.linear_fit(ctime[1362:1366], *fit_opt)
        B_S1_corr[1363:1365] = B_S1_corr_del[1:3]

        fit_opt, fit_covar = curve_fit(
            calibration.linear_fit, ctime[1362:1366], B_S2_corr[1362:1366]
        )
        B_S2_corr_del = calibration.linear_fit(ctime[1362:1366], *fit_opt)
        B_S2_corr[1363:1365] = B_S2_corr_del[1:3]

        fit_opt, fit_covar = curve_fit(
            calibration.linear_fit, ctime[1362:1366], B_S3_corr[1362:1366]
        )
        B_S3_corr_del = calibration.linear_fit(ctime[1362:1366], *fit_opt)
        B_S3_corr[1363:1365] = B_S3_corr_del[1:3]
        """

        spike_idx1s = [1363, 3238]
        spike_idx2s = [1365, 3240]
        for idx in range(len(ctime_idx)): 
            # TEST: sine fit 
            spike_idx1 = spike_idx1s[idx]
            spike_idx2 = spike_idx2s[idx]
            fit_len = 50
            fit_opt, fit_covar = curve_fit(
                calibration.sine_fit, ctime[spike_idx1-fit_len:spike_idx2+fit_len], B_S1_corr[spike_idx1-fit_len:spike_idx2+fit_len],
                p0=[1, np.max(np.abs(B_S1_corr - np.mean(B_S1_corr))), 2.2, 0, 0]
            )
            #B_S1_corr_del = calibration.sine_fit(ctime[spike_idx1-fit_len:spike_idx2+fit_len], *fit_opt)
            B_S1_corr_del = calibration.sine_fit(ctime[spike_idx1:spike_idx2], *fit_opt)
            B_S1_corr[spike_idx1:spike_idx2] = B_S1_corr_del

            fit_opt, fit_covar = curve_fit(
                calibration.sine_fit, ctime[spike_idx1-fit_len:spike_idx2+fit_len], B_S2_corr[spike_idx1-fit_len:spike_idx2+fit_len],
                p0=[1, np.max(np.abs(B_S2_corr - np.mean(B_S2_corr))), 2.2, 0, 0]
            )
            #B_S2_corr_del = calibration.sine_fit(ctime[spike_idx1-fit_len:spike_idx2+fit_len], *fit_opt)
            B_S2_corr_del = calibration.sine_fit(ctime[spike_idx1:spike_idx2], *fit_opt)
            B_S2_corr[spike_idx1:spike_idx2] = B_S2_corr_del
        
            fit_opt, fit_covar = curve_fit(
                calibration.sine_fit, ctime[spike_idx1-fit_len:spike_idx2+fit_len], B_S3_corr[spike_idx1-fit_len:spike_idx2+fit_len],
                p0=[1, np.max(np.abs(B_S3_corr - np.mean(B_S3_corr))), 2.2, 0, 0]
            )
            #B_S3_corr_del = calibration.sine_fit(ctime[spike_idx1-fit_len:spike_idx2+fit_len], *fit_opt)
            B_S3_corr_del = calibration.sine_fit(ctime[spike_idx1:spike_idx2], *fit_opt)
            B_S3_corr[spike_idx1:spike_idx2] = B_S3_corr_del
            #Bplot.B_ctime_plot(
            #    ctime[spike_idx1-fit_len:spike_idx2+fit_len], [B_S1_corr_del, B_S1_corr[spike_idx1-fit_len:spike_idx2+fit_len]], 
            #    [B_S2_corr_del, B_S2_corr[spike_idx1-fit_len:spike_idx2+fit_len]], [B_S3_corr_del, B_S3_corr[spike_idx1-fit_len:spike_idx2+fit_len]], scatter=True, title="del_rogue_10hz")

            if parameter.makeplot == True: 
                Bplot.B_ctime_plot(ctime, B_S1_corr, B_S2_corr, B_S3_corr, xlimt = [(spike_idx1-28)/10, (spike_idx2+28)/10], scatter=True, title=f"del_rogue_10hz_{idx}")
        



    #d_B_S1 = np.gradient(B_S1_corr) / np.gradient(ctime)
    #d_B_S2 = np.gradient(B_S2_corr) / np.gradient(ctime)
    #d_B_S3 = np.gradient(B_S3_corr) / np.gradient(ctime)   
    
    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, B_S1_corr, B_S2_corr, B_S3_corr)
    """
        corr - cross time determination
    """
    try:
        [
            cross_times_corr_1_select, cross_times_corr_1_mids_select, 
            T_spins_d_pad_corr_1_select, w_syn_d_corr_1_select] = cross_time.cross_time_stage_1(
            ctime, B_S3_corr,
        )
    except error.CrossTime1Error as e:
        logger.error(e.__str__())
        logger.error("cross time error return empty!")
        return [ [] for _ in range(16) ]

    [
        cross_times_corr_2_select, cross_times_corr_2_mids_select, 
        T_spins_d_pad_corr_2_select, w_syn_d_corr_2_select] = cross_time.cross_time_stage_2(
        ctime, B_S3_corr, cross_times_corr_1_select, T_spins_d_pad_corr_1_select,
    )
    
    [
        cross_times_corr_3_select, T_spins_d_corr_3_select, w_syn_d_corr_3_select] = cross_time.cross_time_stage_3(
            ctime, B_S3_corr, cross_times_corr_2_select, T_spins_d_pad_corr_2_select
    )

    if parameter.makeplot == True:
        Bplot.B_ctime_plot(
            ctime, B_S1_corr, B_S2_corr, B_S3_corr, 
            cross_times=cross_times_corr_3_select, title = "crosstime", datestr = datestr, xlimt = [ctime[ctime_idx[1]]-10, ctime[ctime_idx[1]]+10]
        )
  
    """
        corr - phase angle integration
    """
    [
        phi_corr, cross_times_corr, w_syn_d_corr, T_spins_d_corr] = cross_time.phase_integration(
        ctime, cross_times_corr_1_select, cross_times_corr_1_mids_select, w_syn_d_corr_1_select, T_spins_d_pad_corr_1_select,
        cross_times_corr_2_select, cross_times_corr_2_mids_select, w_syn_d_corr_2_select, T_spins_d_pad_corr_2_select,
        cross_times_corr_3_select, w_syn_d_corr_3_select, T_spins_d_corr_3_select,
    )    

    """
        IGRF coorindate transformation: gei -> dmxl -> smxl -> fgm
    """
    [
        DMXL_2_GEI, GEI_2_DMXL] = coordinate.dmxl_gei_matrix(
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, att_gei_x, att_gei_y, att_gei_z
    )

    # B igrf rotate from gei to dmxl
    [
        fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z] = coordinate.gei2dmxl(
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, GEI_2_DMXL
    )

    # B igrf rotate from dmxl to smxl
    [
        fgs_igrf_smxl_x, fgs_igrf_smxl_y, fgs_igrf_smxl_z] = coordinate.dmxl2smxl(
            fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, phi_corr
    )

    # B igrf rotate from smxl to fgm
    [
        fgs_igrf_fgm_x, fgs_igrf_fgm_y, fgs_igrf_fgm_z] = coordinate.smxl2fgm(
            fgs_igrf_smxl_x, fgs_igrf_smxl_y, fgs_igrf_smxl_z
    )

    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, [B_S1_corr, fgs_igrf_fgm_x], [B_S2_corr, fgs_igrf_fgm_y], 
    #        [B_S3_corr, fgs_igrf_fgm_z], plot3 = True, title="ful_igrf_fgm_before1stcali")      


    """
        use igrf to calibrate fgs data
    """
    [fgs_ful_fgm_x, fgs_ful_fgm_y, fgs_ful_fgm_z, B_parameter] = calibration.calib_leastsquare(
        B_S1_corr, B_S2_corr, B_S3_corr, fgs_igrf_fgm_x, fgs_igrf_fgm_y, fgs_igrf_fgm_z
    )  

    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, [fgs_ful_fgm_x, fgs_igrf_fgm_x], [fgs_ful_fgm_y, fgs_igrf_fgm_y], 
    #        [fgs_ful_fgm_z, fgs_igrf_fgm_z], plot3 = True, title="ful_igrf_fgm_after1stcali") 
    """
        calib - data cross time determination
    """
    [
        cross_times_calib_1_select, cross_times_calib_1_mids_select, 
        T_spins_d_pad_calib_1_select, w_syn_d_calib_1_select] = cross_time.cross_time_stage_1(
        ctime, fgs_ful_fgm_z,
    )
    #if parameter.makeplot == True:
    #    Bplot.B_ctime_plot(
    #    ctime, fgs_ful_fgm_x, fgs_ful_fgm_y, 
    #    fgs_ful_fgm_z, title="1cross_time_ful_fgm", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-10, ctime[ctime_idx[0]]+10],
    #    gap_time = [ctime[ctime_idx[0]]], cross_times = cross_times_calib_1_select,
    #    )
    [
        cross_times_calib_2_select, cross_times_calib_2_mids_select, 
        T_spins_d_pad_calib_2_select, w_syn_d_calib_2_select] = cross_time.cross_time_stage_2(
        ctime, fgs_ful_fgm_z, cross_times_calib_1_select, T_spins_d_pad_calib_1_select,
    )
    #if parameter.makeplot == True:
    #    Bplot.B_ctime_plot(
    #    ctime, fgs_ful_fgm_x, fgs_ful_fgm_y, 
    #    fgs_ful_fgm_z, title="2cross_time_ful_fgm", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-10, ctime[ctime_idx[0]]+10],
    #    gap_time = [ctime[ctime_idx[0]]], cross_times = cross_times_calib_2_select,
    #    )
    [
        cross_times_calib_3_select, T_spins_d_calib_3_select, w_syn_d_calib_3_select] = cross_time.cross_time_stage_3(
            ctime, fgs_ful_fgm_z, cross_times_calib_2_select, T_spins_d_pad_calib_2_select
    )
    #if parameter.makeplot == True:
    #    Bplot.B_ctime_plot(
    #    ctime, fgs_ful_fgm_x, fgs_ful_fgm_y, 
    #    fgs_ful_fgm_z, title="3cross_time_ful_fgm", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-10, ctime[ctime_idx[0]]+10],
    #    gap_time = [ctime[ctime_idx[0]]], cross_times = cross_times_calib_3_select,
    #    )
  
    """
        calib - phase angle integration
    """
    # Remember that the stage 1 and 2 sample angular velocities at mid points of zero-crossings
    [
        phi_calib, cross_times_calib, w_syn_d_calib, T_spins_d_calib] = cross_time.phase_integration(
        ctime, cross_times_calib_1_select, cross_times_calib_1_mids_select, w_syn_d_calib_1_select, T_spins_d_pad_calib_1_select,
        cross_times_calib_2_select, cross_times_calib_2_mids_select, w_syn_d_calib_2_select, T_spins_d_pad_calib_2_select,
        cross_times_calib_3_select, w_syn_d_calib_3_select, T_spins_d_calib_3_select,
    )
    
    if parameter.makeplot == True:
        Bplot.phase_plot(ctime, phi_calib, cross_times_calib, datestr = datestr, 
            xlimt = [ctime[ctime_idx[1]]-10, ctime[ctime_idx[1]]+10], gap_time = [ctime[ctime_idx[1]]])

    #Bplot.phase_plot(ctime, phi_calib, cross_times_calib)
    #breakpoint()

    """
        fgs data coordinate transform
    """
    # B full rotate from fgm to smxl
    [
        fgs_ful_smxl_x, fgs_ful_smxl_y, fgs_ful_smxl_z] = coordinate.fgm2smxl(
            fgs_ful_fgm_x, fgs_ful_fgm_y, fgs_ful_fgm_z
    )

    # B full rotate from smxl to dmxl
    [
        fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z] = coordinate.smxl2dmxl(
            fgs_ful_smxl_x, fgs_ful_smxl_y, fgs_ful_smxl_z, phi_calib
    )
 
    # B full rotate from dmxl to gei
    [
        fgs_ful_gei_x, fgs_ful_gei_y, fgs_ful_gei_z] = coordinate.dmxl2gei(
            fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, DMXL_2_GEI
    )

    # 2nd calibration
    if parameter.cali_2nd == True:
        #if parameter.makeplot == True: 
        #    Bplot.B_ctime_plot(ctime, [fgs_ful_fgm_x, fgs_igrf_fgm_x], [fgs_ful_fgm_y, fgs_igrf_fgm_y], 
        #        [fgs_ful_fgm_z, fgs_igrf_fgm_z], plot3 = True, title="ful_igrf_fgm_before2ndcali") 
        #    Bplot.B_ctime_plot(ctime, [fgs_ful_dmxl_x0, fgs_igrf_dmxl_x], [fgs_ful_dmxl_y0, fgs_igrf_dmxl_y], 
        #        [fgs_ful_dmxl_z0, fgs_igrf_dmxl_z], plot3 = True, title="ful_igrf_dmxl_before2ndcali")       

        # 2nd calibration of B in dmxl 
        [fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, B_parameter] = calibration.calib_leastsquare(
            fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, init = B_parameter 
        )


        #if parameter.makeplot == True: 
        #    Bplot.B_ctime_plot(ctime, [fgs_ful_dmxl_x, fgs_igrf_dmxl_x], [fgs_ful_dmxl_y, fgs_igrf_dmxl_y], 
        #        [fgs_ful_dmxl_z, fgs_igrf_dmxl_z], plot3 = True, title="ful_igrf_dmxl_after2ndcali") 

    # full res
    fgs_res_dmxl_x = fgs_ful_dmxl_x - fgs_igrf_dmxl_x
    fgs_res_dmxl_y = fgs_ful_dmxl_y - fgs_igrf_dmxl_y
    fgs_res_dmxl_z = fgs_ful_dmxl_z - fgs_igrf_dmxl_z

    # delete rogue points
    if parameter.del_rogue == True:
        del_index = detrend.del_rogue(
            ctime, fgs_res_dmxl_x, fgs_res_dmxl_y, fgs_res_dmxl_z
        )
        ctime = np.delete(ctime, del_index)
        fgs_res_dmxl_x = np.delete(fgs_res_dmxl_x, del_index)
        fgs_res_dmxl_y = np.delete(fgs_res_dmxl_y, del_index)
        fgs_res_dmxl_z = np.delete(fgs_res_dmxl_z, del_index)
        fgs_igrf_dmxl_x = np.delete(fgs_igrf_dmxl_x, del_index)
        fgs_igrf_dmxl_y = np.delete(fgs_igrf_dmxl_y, del_index)
        fgs_igrf_dmxl_z = np.delete(fgs_igrf_dmxl_z, del_index)
        fgs_igrf_gei_x = np.delete(fgs_igrf_gei_x, del_index)
        fgs_igrf_gei_y = np.delete(fgs_igrf_gei_y, del_index)
        fgs_igrf_gei_z = np.delete(fgs_igrf_gei_z, del_index)
        fgs_ful_dmxl_x = np.delete(fgs_ful_dmxl_x, del_index)
        fgs_ful_dmxl_y = np.delete(fgs_ful_dmxl_y, del_index)
        fgs_ful_dmxl_z = np.delete(fgs_ful_dmxl_z, del_index)
        fgs_ful_gei_x = np.delete(fgs_ful_gei_x, del_index)
        fgs_ful_gei_y = np.delete(fgs_ful_gei_y, del_index)
        fgs_ful_gei_z = np.delete(fgs_ful_gei_z, del_index)
        DMXL_2_GEI = np.delete(DMXL_2_GEI, del_index, axis = 0)  

    if parameter.makeplot == True:
        Bplot.B_ctime_plot(
        ctime, fgs_res_dmxl_x, fgs_res_dmxl_y, 
        fgs_res_dmxl_z, title="res_dmxl", datestr = datestr, xlimt = [ctime[ctime_idx[1]]-10, ctime[ctime_idx[1]]+10],
        gap_time = [ctime[ctime_idx[1]]], cross_times = cross_times_calib, scatter = True
        )
    breakpoint()
    # output txt
    if parameter.detrend_fsp == True:
        [
            fgs_res_dmxl_trend_x, fgs_res_dmxl_trend_y, fgs_res_dmxl_trend_z] = detrend.detrend_linear(
                    ctime, fgs_res_dmxl_x, fgs_res_dmxl_y, fgs_res_dmxl_z
        )
        if parameter.makeplot == True:
            Bplot.B_ctime_plot(
            ctime, fgs_res_dmxl_x - fgs_res_dmxl_trend_x, fgs_res_dmxl_y - fgs_res_dmxl_trend_y, 
            fgs_res_dmxl_z - fgs_res_dmxl_trend_z, title="res_dmxl", datestr = datestr
            )

    if parameter.output == True:
        FGM_datetime = list(map(lambda ts: (df["time"][0].to_pydatetime() + 
                        datetime.timedelta(seconds=ts)).strftime('%Y-%m-%d/%H:%M:%S.%f'), ctime))
        output.output_txt(FGM_datetime, fgs_res_dmxl_x, fgs_res_dmxl_y, fgs_res_dmxl_z, title='ela_fgs_res_dmxl')  

        print(f"std res_x: {np.std(fgs_res_dmxl_x)}") 
        print(f"std res_y: {np.std(fgs_res_dmxl_y)}") 
        print(f"std res_z: {np.std(fgs_res_dmxl_z)}")  
        #breakpoint()  
        
    """
        fgs fsp resolution
    """
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
    fgs_fsp_res_dmxl_x = fgs_fsp_ful_dmxl_x - fgs_fsp_igrf_dmxl_x
    fgs_fsp_res_dmxl_y = fgs_fsp_ful_dmxl_y - fgs_fsp_igrf_dmxl_y
    fgs_fsp_res_dmxl_z = fgs_fsp_ful_dmxl_z - fgs_fsp_igrf_dmxl_z

    fgs_fsp_res_gei_x = fgs_fsp_ful_gei_x - fgs_fsp_igrf_gei_x
    fgs_fsp_res_gei_y = fgs_fsp_ful_gei_y - fgs_fsp_igrf_gei_y
    fgs_fsp_res_gei_z = fgs_fsp_ful_gei_z - fgs_fsp_igrf_gei_z

    # detrend
    if parameter.detrend_fsp == True:
        [
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z]= detrend.detrend_linear_2point(
                cross_times_calib, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z
        )

        
        fgs_fsp_res_dmxl_x = fgs_fsp_res_dmxl_x - fgs_fsp_res_dmxl_trend_x
        fgs_fsp_res_dmxl_y = fgs_fsp_res_dmxl_y - fgs_fsp_res_dmxl_trend_y
        fgs_fsp_res_dmxl_z = fgs_fsp_res_dmxl_z - fgs_fsp_res_dmxl_trend_z

        fgs_fsp_ful_dmxl_x = fgs_fsp_res_dmxl_x + fgs_fsp_igrf_dmxl_x
        fgs_fsp_ful_dmxl_y = fgs_fsp_res_dmxl_y + fgs_fsp_igrf_dmxl_y
        fgs_fsp_ful_dmxl_z = fgs_fsp_res_dmxl_z + fgs_fsp_igrf_dmxl_z

        """
            coordinate transform
        """
        DMXL_2_GEI_fsp = cross_time.fsp_matrix(ctime, cross_times_calib, T_spins_d_calib, DMXL_2_GEI)
        [
            fgs_fsp_ful_gei_x, fgs_fsp_ful_gei_y, fgs_fsp_ful_gei_z] = coordinate.dmxl2gei(
                fgs_fsp_ful_dmxl_x, fgs_fsp_ful_dmxl_y, fgs_fsp_ful_dmxl_z, DMXL_2_GEI_fsp
        )
        [
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z] = coordinate.dmxl2gei(
                fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z, DMXL_2_GEI_fsp
        )

        fgs_fsp_res_gei_x = fgs_fsp_ful_gei_x - fgs_fsp_igrf_gei_x
        fgs_fsp_res_gei_y = fgs_fsp_ful_gei_y - fgs_fsp_igrf_gei_y
        fgs_fsp_res_gei_z = fgs_fsp_ful_gei_z - fgs_fsp_igrf_gei_z

    else:    
        fgs_fsp_res_dmxl_trend_x = [0] * len(fgs_fsp_res_gei_x)
        fgs_fsp_res_dmxl_trend_y = [0] * len(fgs_fsp_res_gei_x)
        fgs_fsp_res_dmxl_trend_z = [0] * len(fgs_fsp_res_gei_x)

    if parameter.del_spike_fsp == True:
        """
        [
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, 
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z] = calibration.del_spike_fsp(
            cross_times_calib, ctime[ctime_idx], 
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, 
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z
        )
        """
        for ispike in ctime[ctime_idx]:
            index = min(range(len(cross_times_calib)), key=lambda i: abs(cross_times_calib[i] - ispike))
            if index - 1 < 0 :
                del_index =  [index, index+1]
            elif index + 1 == len(cross_times_calib) :
                del_index =  [index-1, index]
            else:
                del_index =  [index-1, index, index+1]

            cross_times_calib = np.delete(cross_times_calib, del_index)
            fgs_fsp_res_dmxl_x = np.delete(fgs_fsp_res_dmxl_x, del_index)
            fgs_fsp_res_dmxl_y = np.delete(fgs_fsp_res_dmxl_y, del_index)
            fgs_fsp_res_dmxl_z = np.delete(fgs_fsp_res_dmxl_z, del_index)
            fgs_fsp_igrf_dmxl_x = np.delete(fgs_fsp_igrf_dmxl_x, del_index)
            fgs_fsp_igrf_dmxl_y = np.delete(fgs_fsp_igrf_dmxl_y, del_index)
            fgs_fsp_igrf_dmxl_z = np.delete(fgs_fsp_igrf_dmxl_z, del_index)
            fgs_fsp_res_dmxl_trend_x = np.delete(fgs_fsp_res_dmxl_trend_x, del_index)
            fgs_fsp_res_dmxl_trend_y = np.delete(fgs_fsp_res_dmxl_trend_y, del_index)
            fgs_fsp_res_dmxl_trend_z = np.delete(fgs_fsp_res_dmxl_trend_z, del_index)
            fgs_fsp_res_gei_x = np.delete(fgs_fsp_res_gei_x, del_index)
            fgs_fsp_res_gei_y = np.delete(fgs_fsp_res_gei_y, del_index)
            fgs_fsp_res_gei_z = np.delete(fgs_fsp_res_gei_z, del_index)
            fgs_fsp_igrf_gei_x = np.delete(fgs_fsp_igrf_gei_x, del_index)
            fgs_fsp_igrf_gei_y = np.delete(fgs_fsp_igrf_gei_y, del_index)
            fgs_fsp_igrf_gei_z = np.delete(fgs_fsp_igrf_gei_z, del_index)

    if parameter.makeplot == True:
        #Bplot.B_ctime_plot(
        #    cross_times_calib, fgs_fsp_res_dmxl_x, 
        #    fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, title="res_dmxl_fsp", scatter = True, 
        #    gap_time = ctime[ctime_idx], datestr = datestr, xlimt = [ctime[ctime_idx[0]]-10, ctime[ctime_idx[0]]+10]
        #)
        Bplot.B_ctime_plot(
            cross_times_calib, fgs_fsp_res_dmxl_x, 
            fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, title="res_dmxl_fsp", scatter = True, 
            gap_time = ctime[ctime_idx], datestr = datestr,
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
