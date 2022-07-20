import calendar
import datetime
import random
import string
import traceback
from pprint import pprint
from tracemalloc import start
from typing import List, Union

import numpy as np
import pandas as pd
from cdflib import CDF, cdfepoch
from fastapi import APIRouter, Query
from geopack import geopack
from pyspedas.cotrans import cotrans_lib
from scipy import signal
from scipy.integrate import simpson, trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import function.Bplot as Bplot
import function.calibration as calibration
import function.coordinate as coordinate
import function.cross_time as cross_time
import parameter


router = APIRouter(
    prefix="/fgm_calib",
    tags=["fgm_calib"],
    responses={404: {"description": "Not found"}},
)


@router.get("/get_numbers")
def get_two_numbers():
    a = random.randint(1, 10)
    b = random.randint(1, 100)
    return (a, b)


@router.get("/get_igrf")
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


@router.get("/get_cdf")
def get_cdf(cdfpath: str, vars: Union[list[str], None] = Query(default=None)):

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


@router.get("/fgm_calib")
def fgm_fsp_calib(starttime_str: str, endtime_str: str, sta_cdfpath: str, fgm_cdfpath: str):

    df = pd.DataFrame()

    # time range for sci zone
    starttime = pd.to_datetime(starttime_str)
    endtime = pd.to_datetime(endtime_str)

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
    delta_t = np.median(ctime[1:] - ctime[:-1])
    # get fgm data in fgm coordinate
    B_S1_corr, B_S2_corr, B_S3_corr = list(zip(*df["fgm_fgm"]))

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
    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z = np.array(list(zip(*df["igrf_gei"])))
    att_gei_x, att_gei_y, att_gei_z = np.array(list(zip(*df["att_gei"])))
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

    #Bplot.B2_ctime_plot(ctime, fgs_igrf_smxl_x, fgs_igrf_smxl_y, fgs_igrf_smxl_z, fgs_ful_smxl_x, fgs_ful_smxl_y, fgs_ful_smxl_z, "igrf and ful in smxl")
    #Bplot.B_ctime_plot(ctime, fgs_ful_smxl_x-fgs_igrf_smxl_x, fgs_ful_smxl_y-fgs_igrf_smxl_y, fgs_ful_smxl_z-fgs_igrf_smxl_z, "igrf and ful in smxl")
    #Bplot.B2_ctime_plot(
    #    ctime, fgs_res0_dmxl_x, fgs_res0_dmxl_y, fgs_res0_dmxl_z, fgs_res_dmxl_trend_x, 
    #    fgs_res_dmxl_trend_y, fgs_res_dmxl_trend_z, "X1 = res_dmxl    X2 = res_dmxl_detrend")     

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
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z]= calibration.detrend_linear(
                cross_times_calib, fgs_fsp_res0_dmxl_x, fgs_fsp_res0_dmxl_y, fgs_fsp_res0_dmxl_z
        )

        #Bplot.B2_ctime_plot(
        #cross_times_calib, fgs_fsp_res0_dmxl_x, fgs_fsp_res0_dmxl_y, fgs_fsp_res0_dmxl_z, 
        #fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z, "X1 = fsp_res    X2 = fsp_res_trend")
        # 

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
        
        #Bplot.B_ctime_plot(
        #cross_times_calib, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, "X1 = fsp_res_detrend")     

    else:
        """
            fgs res data detrend
        """
        fgs_res0_dmxl_x = fgs_ful_dmxl_x - fgs_igrf_dmxl_x
        fgs_res0_dmxl_y = fgs_ful_dmxl_y - fgs_igrf_dmxl_y
        fgs_res0_dmxl_z = fgs_ful_dmxl_z - fgs_igrf_dmxl_z

        [
            fgs_res_dmxl_trend_x, fgs_res_dmxl_trend_y, fgs_res_dmxl_trend_z]= calibration.detrend_linear(
                ctime, fgs_res0_dmxl_x, fgs_res0_dmxl_y, fgs_res0_dmxl_z
        )

        fgs_res_dmxl_x = fgs_res0_dmxl_x - fgs_res_dmxl_trend_x
        fgs_res_dmxl_y = fgs_res0_dmxl_y - fgs_res_dmxl_trend_y
        fgs_res_dmxl_z = fgs_res0_dmxl_z - fgs_res_dmxl_trend_z

        #Bplot.B2_ctime_plot(
        #ctime, fgs_res0_dmxl_x, fgs_res0_dmxl_y, fgs_res0_dmxl_z, 
        #fgs_res_dmxl_trend_x, fgs_res_dmxl_trend_y, fgs_res_dmxl_trend_z, "X1 = res    X2 = res_trend")     

        """
            fgs ful field data dmxl to gei
        """
        fgs_ful_dmxl_x = fgs_res_dmxl_x + fgs_igrf_dmxl_x
        fgs_ful_dmxl_y = fgs_res_dmxl_y + fgs_igrf_dmxl_y
        fgs_ful_dmxl_z = fgs_res_dmxl_z + fgs_igrf_dmxl_z
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

        #Bplot.B_ctime_plot(
        #cross_times_calib, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, "X1 = fsp_res_detrend")  

        fgs_fsp_res_gei_x = fgs_fsp_ful_gei_x - fgs_fsp_igrf_gei_x
        fgs_fsp_res_gei_y = fgs_fsp_ful_gei_y - fgs_fsp_igrf_gei_y
        fgs_fsp_res_gei_z = fgs_fsp_ful_gei_z - fgs_fsp_igrf_gei_z


    #Bplot.B_ctime_plot(cross_times_calib,fgs_fsp_res_dmxl_x,fgs_fsp_res_dmxl_y,fgs_fsp_res_dmxl_z,"fsp_res_dmxl")

    return [
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


if __name__ == "__main__":

    starttime_str = "2022-01-12 15:45:51"
    endtime_str = "2022-01-12 15:52:04"

    sta_cdfpath = "test/ela_l1_state_defn_20220112_v01.cdf"
    fgm_cdfpath = "test/ela_l1_fgs_20220112_v01.cdf"

    #starttime_str = "2022-01-14 15:45:50"
    #endtime_str = "2022-01-14 15:52:04"

    #sta_cdfpath = "test/ela_l1_state_defn_20220114_v01.cdf"
    #fgm_cdfpath = "test/ela_l1_fgs_20220114_v01.cdf"

    #starttime_str = "2022-06-27 08:52:34"
    #endtime_str = "2022-06-27 08:58:47"

    #sta_cdfpath = "test/ela_l1_state_defn_20220627_v01.cdf"
    #fgm_cdfpath = "test/ela_l1_fgs_20220627_v01.cdf"

    #starttime_str = "2022-01-25 16:28:23"
    #endtime_str = "2022-01-25 16:34:35"

    #sta_cdfpath = "test/ela_l1_state_defn_20220125_v01.cdf"
    #fgm_cdfpath = "test/ela_l1_fgs_20220125_v01.cdf"

    #starttime_str = "2022-07-06/09:45:53"
    #endtime_str = "2022-07-06/09:52:11"

    #sta_cdfpath = "test/ela_l1_state_defn_20220706_v01.cdf"
    #fgm_cdfpath = "test/ela_l1_fgs_20220706_v01.cdf"

    [
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
    ] = fgm_fsp_calib(starttime_str, endtime_str, sta_cdfpath, fgm_cdfpath)

