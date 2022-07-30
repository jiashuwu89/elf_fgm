import datetime
import traceback

import numpy as np
import pandas as pd
from cdflib import CDF, cdfepoch
from geopack import geopack
from pyspedas.cotrans import cotrans_lib
from scipy.interpolate import interp1d
from .function import calibration
from .function import coordinate
from .function import cross_time
from .function import Bplot
from .function import detrend
from . import parameter


def get_igrf(time: datetime.datetime, xgsm: float, ygsm: float, zgsm: float):

    """
    Input
        xgsm,ygsm,zgsm: The given position in cartesian GSM coordinate in
        Re (earth radii, 1 Re = 6371.2 km).
    Return
        bxgsm,bygsm,bzgsm: Cartesian GSM components of the internal magnetic
        field in nT.
    """

    t1 = time
    t0 = datetime.datetime(1970, 1, 1)
    ut = (t1 - t0).total_seconds()
    tilt = geopack.recalc(ut)
    Re = 6371.2  # in km
    xgsm = xgsm / Re
    ygsm = ygsm / Re
    zgsm = zgsm / Re
    bxgsm, bygsm, bzgsm = geopack.igrf_gsm(xgsm, ygsm, zgsm)
    return (bxgsm, bygsm, bzgsm)


# def get_cdf(cdfpath: str, vars: Union[list[str], None]):
def get_cdf(cdfpath, vars):
    
    try:
        cdf = CDF(cdfpath)
        cdfinfo = cdf.cdf_info()
        data = {}
        if vars is None:
            vars = cdfinfo["zVariables"]
            print(f"{cdfpath} variables: {vars}")
        for var in vars:
            val = cdf.varget(var)
            if var.endswith("_time"):
                data[var] = list(
                    map(lambda t: cdfepoch.to_datetime(t)[0], val.tolist())
                )
            elif isinstance(val, np.ndarray):
                data[var] = val.tolist()
            else:
                data[var] = val

        return data

    except Exception as e:
        return {
            "message": "Failed to open state CDF",
            "error": str(e),
            "traceback": "".join(traceback.format_exception(None, e, e.__traceback__)),
        }


def clip_cdfdata(
    df: pd.Series, starttime: pd.Timestamp, endtime: pd.Timestamp
) -> pd.Series:

    startindex = df.index.get_indexer([starttime], method="nearest")[0]
    endindex = df.index.get_indexer([endtime], method="nearest")[0]

    return df[startindex:endindex]


def resample_data(
    cur_time: pd.Timestamp, cur_data: pd.Series, target_time: pd.Timestamp
) -> pd.Series:

    cur_data_np = np.array(cur_data.to_list())
    dimens = cur_data_np.shape[1]  # x, y, z
    interp_data = np.zeros((len(target_time), dimens))

    x = (cur_time - target_time[0]).total_seconds()
    x_interp = (target_time - target_time[0]).total_seconds()
    for dimen in range(dimens):
        f = interp1d(x, cur_data_np[:, dimen])
        interp_data[:, dimen] = f(x_interp)

    return pd.Series(interp_data.tolist())


def fgm_fsp_calib(starttime_str: str, endtime_str: str, sta_cdfpath: str, fgm_cdfpath: str):

    df = pd.DataFrame()

    # time range for sci zone
    starttime = pd.to_datetime(starttime_str)
    endtime = pd.to_datetime(endtime_str)

    """
        initial points exclude
    """
    if parameter.init_secs != 0:
        starttime = starttime + datetime.timedelta(seconds = parameter.init_secs)

    # read state cdf for att
    att_cdfdata = pd.DataFrame(
        get_cdf(sta_cdfpath, vars=["ela_att_time", "ela_att_gei"])
    )
    att_cdfdata.set_index("ela_att_time", inplace=True)
    att_cdfdata = clip_cdfdata(
        att_cdfdata,
        starttime - datetime.timedelta(minutes=2),
        endtime + datetime.timedelta(minutes=2),
    )

    # read state cdf for pos; not read together with att b/c different length
    pos_cdfdata = pd.DataFrame(
        get_cdf(sta_cdfpath, vars=["ela_pos_gei"])
    )  # not read state time b/c too slow
    pos_cdfdata["ela_pos_time"] = pd.date_range(
        start=starttime_str[0:10], periods=len(pos_cdfdata["ela_pos_gei"]), freq="S"
    )
    pos_cdfdata.set_index("ela_pos_time", inplace=True)
    pos_cdfdata = clip_cdfdata(
        pos_cdfdata,
        starttime - datetime.timedelta(minutes=2),
        endtime + datetime.timedelta(minutes=2),
    )

    # read fgm cdf and clip
    fgm_cdfdata = pd.DataFrame(get_cdf(fgm_cdfpath, vars=["ela_fgs_time", "ela_fgs"]))
    fgm_cdfdata.set_index("ela_fgs_time", inplace=True)
    fgm_cdfdata = clip_cdfdata(fgm_cdfdata, starttime, endtime)

    # resample att and pos to fgm time resolution
    df["att_gei"] = resample_data(
        att_cdfdata.index, att_cdfdata["ela_att_gei"], fgm_cdfdata.index
    )
    df["pos_gei"] = resample_data(
        pos_cdfdata.index, pos_cdfdata["ela_pos_gei"], fgm_cdfdata.index
    )
    df["fgm_fgm"] = pd.Series(fgm_cdfdata["ela_fgs"].tolist())
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
        get_igrf(
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

    ############################################
    #   suyash code begin
    ############################################
    df["ctime"] = df["timestamp"] - df["timestamp"][0]
    ctime = np.array(df["ctime"])
    # get fgm data in fgm coordinate
    B_S1_corr, B_S2_corr, B_S3_corr = np.array(list(zip(*df["fgm_fgm"])))
    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z = np.array(list(zip(*df["igrf_gei"])))
    att_gei_x, att_gei_y, att_gei_z = np.array(list(zip(*df["att_gei"])))
    
    """
        exclude bad data; TODO: automatic detection
    """
    #breakpoint()
    badlim1 = 30
    badlim2 = 50
    if parameter.bad_data_correct == True:
        ctime = np.delete(ctime, np.s_[parameter.bad_data-badlim1:parameter.bad_data+badlim2])
        B_S1_corr = np.delete(B_S1_corr, np.s_[parameter.bad_data-badlim1:parameter.bad_data+badlim2])
        B_S2_corr = np.delete(B_S2_corr, np.s_[parameter.bad_data-badlim1:parameter.bad_data+badlim2])
        B_S3_corr = np.delete(B_S3_corr, np.s_[parameter.bad_data-badlim1:parameter.bad_data+badlim2])
        fgs_igrf_gei_x = np.delete(fgs_igrf_gei_x, np.s_[parameter.bad_data-badlim1:parameter.bad_data+badlim2])
        fgs_igrf_gei_y = np.delete(fgs_igrf_gei_y, np.s_[parameter.bad_data-badlim1:parameter.bad_data+badlim2])
        fgs_igrf_gei_z = np.delete(fgs_igrf_gei_z, np.s_[parameter.bad_data-badlim1:parameter.bad_data+badlim2])
        att_gei_x = np.delete(att_gei_x, np.s_[parameter.bad_data-badlim1:parameter.bad_data+badlim2])
        att_gei_y = np.delete(att_gei_y, np.s_[parameter.bad_data-badlim1:parameter.bad_data+badlim2])
        att_gei_z = np.delete(att_gei_z, np.s_[parameter.bad_data-badlim1:parameter.bad_data+badlim2])


    """
        timestamp correction; TODO: automatic detection
    """
    
    ctime, ctime_idx = calibration.ctime_calib(ctime, B_x = B_S1_corr, B_y = B_S2_corr, B_z = B_S3_corr)

    d_B_S1 = np.gradient(B_S1_corr) / np.gradient(ctime)
    d_B_S2 = np.gradient(B_S2_corr) / np.gradient(ctime)
    d_B_S3 = np.gradient(B_S3_corr) / np.gradient(ctime)   
    #Bplot.B1_ctime1_plot3(ctime, d_B_S1, d_B_S2, d_B_S3, err = ctime[ctime_idx], xlimt=[58, 66], title="dB")
    #breakpoint()
    #Bplot.B1_ctime1_plot3(ctime, B_S1_corr, B_S2_corr, B_S3_corr, err = ctime[ctime_idx], xlimt=[78, 86], title="B")
    #breakpoint()
    """
        corr - cross time determination
    """
    [
        cross_times_corr_1_select, cross_times_corr_1_mids_select, 
        T_spins_d_pad_corr_1_select, w_syn_d_corr_1_select] = cross_time.cross_time_stage_1(
        ctime, B_S3_corr,
    )

    [
        cross_times_corr_2_select, cross_times_corr_2_mids_select, 
        T_spins_d_pad_corr_2_select, w_syn_d_corr_2_select] = cross_time.cross_time_stage_2(
        ctime, B_S3_corr, cross_times_corr_1_select, T_spins_d_pad_corr_1_select,
    )

    [
        cross_times_corr_3_select, T_spins_d_corr_3_select, w_syn_d_corr_3_select] = cross_time.cross_time_stage_3(
            ctime, B_S3_corr, cross_times_corr_2_select, T_spins_d_pad_corr_2_select
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

    """
        use igrf to calibrate fgs data
    """
    [B_S1_calib, B_S2_calib, B_S3_calib] = calibration.calib_leastsquare(
        B_S1_corr, B_S2_corr, B_S3_corr, fgs_igrf_fgm_x, fgs_igrf_fgm_y, fgs_igrf_fgm_z
    )

    if parameter.makeplot == True: 
        #Bplot.phase_plot(ctime, phi_corr, cross_times_corr, ctime[ctime_idx], xlimt=[70,85])
        #breakpoint()
        #Bplot.B2_ctime1_plot3(ctime, B_S1_calib, B_S2_calib, B_S3_calib, fgs_igrf_fgm_x, fgs_igrf_fgm_y, fgs_igrf_fgm_z, ctime[ctime_idx], xlimt=[70,85])
        Bplot.B1_ctime1_plot3(ctime, B_S1_calib-fgs_igrf_fgm_x, B_S2_calib-fgs_igrf_fgm_y, 
            B_S3_calib-fgs_igrf_fgm_z, err = ctime[ctime_idx], cross_times = cross_times_corr)
        #breakpoint()
    """
        calib - data cross time determination
    """
    [
        cross_times_calib_1_select, cross_times_calib_1_mids_select, 
        T_spins_d_pad_calib_1_select, w_syn_d_calib_1_select] = cross_time.cross_time_stage_1(
        ctime, B_S3_calib,
    )

    [
        cross_times_calib_2_select, cross_times_calib_2_mids_select, 
        T_spins_d_pad_calib_2_select, w_syn_d_calib_2_select] = cross_time.cross_time_stage_2(
        ctime, B_S3_calib, cross_times_calib_1_select, T_spins_d_pad_calib_1_select,
    )

    [
        cross_times_calib_3_select, T_spins_d_calib_3_select, w_syn_d_calib_3_select] = cross_time.cross_time_stage_3(
            ctime, B_S3_calib, cross_times_calib_2_select, T_spins_d_pad_calib_2_select
    )
    
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
 

    """
        fgs data coordinate transform
    """
    # B full rotate from fgm to smxl
    [
        fgs_ful_smxl_x, fgs_ful_smxl_y, fgs_ful_smxl_z] = coordinate.fgm2smxl(
            B_S1_calib, B_S2_calib, B_S3_calib
    )

    # B full rotate from smxl to dmxl
    [
        fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z] = coordinate.smxl2dmxl(
            fgs_ful_smxl_x, fgs_ful_smxl_y, fgs_ful_smxl_z, phi_calib
    )

    if parameter.makeplot == True: 
        Bplot.B2_ctime1_plot3(ctime, fgs_ful_dmxl_x, fgs_ful_dmxl_y, 
            fgs_ful_dmxl_z, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, title="dmxl")

        #Bplot.B2_ctime1_plot3(ctime, fgs_ful_smxl_x, fgs_ful_smxl_y, 
        #    fgs_ful_smxl_z, fgs_igrf_smxl_x, fgs_igrf_smxl_y, fgs_igrf_smxl_z, title="smxl")    

    if parameter.detrend_fsp == True:
        """
            generate fsp of igrf and ful
        """
        [
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z] = cross_time.fsp_igrf(
                ctime, cross_times_calib, T_spins_d_calib, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z
        )

        [
            fgs_fsp_ful_dmxl_x, fgs_fsp_ful_dmxl_y, fgs_fsp_ful_dmxl_z] = cross_time.fsp_igrf(
                ctime, cross_times_calib, T_spins_d_calib, fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z
        )

        """
            detrend fsp_res
        """
        fgs_fsp_res0_dmxl_x = fgs_fsp_ful_dmxl_x - fgs_fsp_igrf_dmxl_x
        fgs_fsp_res0_dmxl_y = fgs_fsp_ful_dmxl_y - fgs_fsp_igrf_dmxl_y
        fgs_fsp_res0_dmxl_z = fgs_fsp_ful_dmxl_z - fgs_fsp_igrf_dmxl_z

        [
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z]= detrend.detrend_linear_2point(
                cross_times_calib, fgs_fsp_res0_dmxl_x, fgs_fsp_res0_dmxl_y, fgs_fsp_res0_dmxl_z
        )

        fgs_fsp_res_dmxl_x = fgs_fsp_res0_dmxl_x - fgs_fsp_res_dmxl_trend_x
        fgs_fsp_res_dmxl_y = fgs_fsp_res0_dmxl_y - fgs_fsp_res_dmxl_trend_y
        fgs_fsp_res_dmxl_z = fgs_fsp_res0_dmxl_z - fgs_fsp_res_dmxl_trend_z

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
        """
            fgs res data full resolution
        """
        fgs_res_dmxl_x = fgs_ful_dmxl_x - fgs_igrf_dmxl_x
        fgs_res_dmxl_y = fgs_ful_dmxl_y - fgs_igrf_dmxl_y
        fgs_res_dmxl_z = fgs_ful_dmxl_z - fgs_igrf_dmxl_z    

        """
            fgs ful field data dmxl to gei
        """
        [
            fgs_ful_gei_x, fgs_ful_gei_y, fgs_ful_gei_z] = coordinate.dmxl2gei(
                fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, DMXL_2_GEI
        )
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

        fgs_fsp_res_dmxl_trend_x = 0
        fgs_fsp_res_dmxl_trend_y = 0
        fgs_fsp_res_dmxl_trend_z = 0

    if parameter.makeplot == True:
        Bplot.B1_ctime1_plot3(cross_times_calib, fgs_fsp_res_dmxl_x, 
            fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, ctime[ctime_idx], xlimt = [10, 355],
            title="res in dmxl with detrend")
    """
            fgs ful field data dmxl to gei
    """
    #FGM_datetime = list(map(lambda ts: (df["time"][0].to_pydatetime() + 
    #                                datetime.timedelta(seconds=ts)).strftime('%Y-%m-%d/%H:%M:%S'), 
    #                    cross_times_calib))

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
