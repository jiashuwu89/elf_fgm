import logging
import random
import traceback

import pandas as pd
import datetime as dt
from typing import List, Literal
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from ..fgm_utils import fgm_fsp_calib, fgm_fsp_calib_prepos_wrapper
from ..fgm_utils.function import error

router = APIRouter(
    prefix="/fgm_calib",
    tags=["fgm_calib"],
    responses={404: {"description": "Not found"}},
)


class FgmCalibRequest(BaseModel):
    mission_id: Literal[1, 2]
    fgs_time: List[dt.datetime]
    fgs: List[List[float]]


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
    Assumes that the data provided is a chronologically sorted list of data
    points corresponding to a single science zone collection.

    TODO: Rename to something clearer - maybe just calculate_fsp_data and GET /fsp?
    """
    logger = logging.getLogger("sp")
    if not fgm_calib_request.fgs_time or not fgm_calib_request.fgs:
        raise HTTPException(status_code=404, detail="empty science zone")

    # Reformat/Restructure input data
    mission = "ela" if fgm_calib_request.mission_id == 1 else "elb"
    start_time = fgm_calib_request.fgs_time[0]
    end_time = fgm_calib_request.fgs_time[-1]
    fgm_data = pd.DataFrame({
        f"{mission}_fgs_time": fgm_calib_request.fgs_time,
        f"{mission}_fgs": fgm_calib_request.fgs,
    })
    logger.info(f"▶️ Received {mission} collection from {start_time} to {end_time}")

    [ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
            att_gei_x, att_gei_y, att_gei_z, 
            pos_gei_x, pos_gei_y, pos_gei_z,
            f_all_arry, clip_start_idx, clip_end_idx] = fgm_fsp_calib_prepos_wrapper(mission, start_time, end_time, None, logger)

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
        B_parameter]=fgm_fsp_calib(
            ctime, ctimestamp, f_all_arry,
            fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z,
            att_gei_x, att_gei_y, att_gei_z,
            pos_gei_x, pos_gei_y, pos_gei_z,
            logger, mission
    )
        
    except Exception as e:
        traceback_msg = traceback.format_exc()
        logger.error(f"fsp calibration failed ({e}): {traceback_msg}")
        raise
    logger.info(f"⏹️ End of fsp calibration for {mission} from {start_time} to {end_time}\n")

    # Note: Transposing
    return FgmCalibResponse(
        fgs_fsp_time=list(FGM_timestamp),
        fgs_fsp_res_dmxl=list(map(list, zip(fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z))),
        fgs_fsp_res_dmxl_trend=list(map(list, zip(fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z))),
        fgs_fsp_res_gei=list(map(list, zip(fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z))),
        fgs_fsp_igrf_dmxl=list(map(list, zip(fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z))),
        fgs_fsp_igrf_gei=list(map(list, zip(fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z))),
    )
