import logging
import os
import random

import pandas as pd
import datetime as dt
from typing import List
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from ..fgm_utils import fgm_fsp_calib
from ..fgm_utils import parameter


router = APIRouter(
    prefix="/fgm_calib",
    tags=["fgm_calib"],
    responses={404: {"description": "Not found"}},
)


class FgmCalibRequest(BaseModel):
    fgs_time: List[dt.datetime]
    fgs: List[List[float]]


# TODO: Confirm all are floats?
class FgmCalibResponse(BaseModel):
    fgs_fsp_time: List[dt.datetime]
    fgs_fsp_res_dmxl: List[List[float]]
    fgs_fsp_res_dmxl_trend: List[List[float]]
    fgs_fsp_res_gei: List[List[float]]
    fgs_fsp_igrf_dmxl: List[List[float]]
    fgs_fsp_igrf_gei: List[List[float]]


@router.get("/get_numbers")
def get_two_numbers():
    a = random.randint(1, 10)
    b = random.randint(1, 100)
    return (a, b)


@router.post("/fgm_calib")
def fgm_calib(fgm_calib_request: FgmCalibRequest) -> FgmCalibResponse:
    """
    Assumes that timestamps and data is sorted chronologically
    """
    logger = logging.getLogger("fgm_calib.fgm_calib")
    if not fgm_calib_request.fgs_time or not fgm_calib_request.fgs:
        raise HTTPException(status_code=404, detail="empty science zone")

    # TODO: Remove when we read multiple state CDFs
    if fgm_calib_request.fgs_time[0].date() != fgm_calib_request.fgs_time[-1].date():
        raise HTTPException(status_code=404, detail="start and end time should be same day")

    starttime_str = fgm_calib_request.fgs_time[0].strftime("%Y-%m-%d %H:%M:%S")
    endtime_str = fgm_calib_request.fgs_time[-1].strftime("%Y-%m-%d %H:%M:%S")

    # TODO: Avoid hardcoding ELA
    sta_datestr = fgm_calib_request.fgs_time[0].strftime("%Y%m%d")
    sta_cdfpath = os.path.join(parameter.STATE_DATA_DIR, f"ela_l1_state_defn_{sta_datestr}_v01.cdf")

    # TODO: Avoid hardcoding ELA
    fgm_data = pd.DataFrame({
        "ela_fgs_time": fgm_calib_request.fgs_time,
        "ela_fgs": fgm_calib_request.fgs,
    })

    logger.info(starttime_str, endtime_str, sta_cdfpath, fgm_data)

    try:
        [
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
        ] = fgm_fsp_calib(starttime_str, endtime_str, sta_cdfpath, fgm_data)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

    # Note: Transposing
    return FgmCalibResponse(
        fgs_fsp_time=FGM_timestamp,
        fgs_fsp_res_dmxl=list(map(list, zip(fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z))),
        fgs_fsp_res_dmxl_trend=list(map(list, zip(fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z))),
        fgs_fsp_res_gei=list(map(list, zip(fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z))),
        fgs_fsp_igrf_dmxl=list(map(list, zip(fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z))),
        fgs_fsp_igrf_gei=list(map(list, zip(fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z))),
    )
