import random

from fastapi import APIRouter
from pydantic import BaseModel

from ..fgm_utils import fgm_fsp_calib


router = APIRouter(
    prefix="/fgm_calib",
    tags=["fgm_calib"],
    responses={404: {"description": "Not found"}},
)


class FgmCalibRequest(BaseModel):
    # TODO: fill this in!
    temp: int


# TODO: Confirm all are floats?
class FgmCalibResponse(BaseModel):
    # TODO: fill this in!
    temp: int
    # fgm_fgs_fsp_dmxl: List[float]
    # fgm_fgs_fsp_igrf_dmxl: List[float]
    # fgm_fgs_fsp_gei: List[float]
    # fgm_fgs_fsp_igrf_gei: List[float]


@router.get("/get_numbers")
def get_two_numbers():
    a = random.randint(1, 10)
    b = random.randint(1, 100)
    return (a, b)


@router.get("/fgm_calib")
def fgm_calib(starttime_str: str, endtime_str: str, sta_cdfpath: str, fgm_cdfpath: str) -> FgmCalibResponse:
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
    ] = fgm_fsp_calib(starttime_str, endtime_str, sta_cdfpath, fgm_cdfpath)

    return FgmCalibResponse(temp=0)
