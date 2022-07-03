
from fastapi import APIRouter, Query
from geopack import geopack
import random
import datetime
import traceback
import calendar
from pprint import pprint
from cdflib import CDF, cdfepoch
from typing import Union
import numpy as np
from pyspedas.cotrans import cotrans_lib 

router = APIRouter(
    prefix="/jwu_test",
    tags=["jwu_test"],
    responses={404: {"description": "Not found"}},
)


@router.get("/get_numbers")
def get_two_numbers():
    a = random.randint(1, 10)
    b = random.randint(1, 100)
    return (a, b)


@router.get("/get_igrf")
def get_igrf(
        time: datetime.datetime,
        xgsm: float,
        ygsm: float,
        zgsm: float):

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
    ut = (t1-t0).total_seconds()
    tilt = geopack.recalc(ut)
    Re = 6371.2  # in km
    xgsm = xgsm/Re
    ygsm = ygsm/Re
    zgsm = zgsm/Re
    bxgsm, bygsm, bzgsm = geopack.igrf_gsm(xgsm, ygsm, zgsm)
    return (bxgsm, bygsm, bzgsm)


@router.get("/get_cdf")
def get_cdf(
        cdfpath: str,
        vars: Union[list[str], None] = Query(default=None)):

    try:
        cdf = CDF(cdfpath)
        cdfinfo = cdf.cdf_info()
        data = {}
        if vars is None:
            vars = cdfinfo["zVariables"]

        for var in vars:
            val = cdf.varget(var)
            if var.endswith("_time"):
                data[var] = list(map(lambda t: cdfepoch.to_datetime(t)[0], val.tolist()))
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


if __name__ == "__main__":

    # read cdf
    cdfpath = "test/ela_l1_state_defn_20220112_v01.cdf"
    cdfdata = get_cdf(cdfpath, vars=["ela_att_time", "ela_att_gei", "ela_pos_gei"])

    # time range for sci zone
    starttime = datetime.datetime(2022, 1, 12, 15, 45, 59)
    endtime = datetime.datetime(2022, 1, 12, 15, 52, 4)

    # get index for starttime and endtime
    startindex = cdfdata["ela_att_time"].index(min(
        cdfdata["ela_att_time"], key=lambda t: abs(t - starttime)))
    endindex = cdfdata["ela_att_time"].index(min(
        cdfdata["ela_att_time"], key=lambda t: abs(t - endtime)))

    # clip cdf data
    data = {}
    data["ela_att_time"] = cdfdata["ela_att_time"][startindex:endindex]
    data["ela_att_gei"] = cdfdata["ela_att_gei"][startindex:endindex]  # time res 1 min
    data["ela_pos_gei"] = cdfdata["ela_pos_gei"][startindex*60:endindex*60:60]  # time res 1 s
  
    # covert datetime to utc timestamp 
    # ref: https://stackoverflow.com/questions/5067218/get-utc-timestamp-in-python-with-datetime
    data["ela_timestamp"] = [calendar.timegm(ts.utctimetuple()) for ts in data["ela_att_time"]]
    #iyear, idoy, ih, im, isec = cotrans_lib.get_time_parts(data["ela_timestamp"])
    #print(f"year:{iyear}, doy:{idoy}, h:{ih}, m:{im}, sec:{isec}")

    # coordinate transformation of pos: gei -> gse -> gsm
    data["ela_pos_gse"] = cotrans_lib.subgei2gse(data["ela_timestamp"], data["ela_pos_gei"])
    data["ela_pos_gsm"] = cotrans_lib.subgse2gsm(data["ela_timestamp"], data["ela_pos_gse"])

    # get igrf b in gsm
    data["ela_igrf_gsm"] = [get_igrf(
        data["ela_att_time"][i],
        data["ela_pos_gsm"][i][0],
        data["ela_pos_gsm"][i][1],
        data["ela_pos_gsm"][i][2])
        for i in range(len(data["ela_att_time"]))]

    # coordinate transformation of B: gsm -> gse -> gei
    data["ela_igrf_gse"] = cotrans_lib.subgsm2gse(data["ela_timestamp"], data["ela_igrf_gsm"])
    data["ela_igrf_gei"] = cotrans_lib.subgse2gei(data["ela_timestamp"], data["ela_igrf_gse"])


    pprint([data["ela_att_time"][0], data["ela_igrf_gei"][0], data["ela_att_gei"][0], data["ela_pos_gei"][0], data["ela_pos_gsm"][0], data["ela_pos_gse"][0]])

    
    # tstart = datetime.datetime(2022, 1, 12, 15, 45, 59)
    # xgsm = -2431.1245629621699
    # ygsm = 3822.9186030446831
    # zgsm = 5059.6970615621403
    # bxgsm, bygsm, bzgsm = get_igrf(tstart, xgsm, ygsm, zgsm)
    # print(bxgsm, bygsm, bzgsm)
