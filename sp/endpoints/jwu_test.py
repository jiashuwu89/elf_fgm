
from fastapi import APIRouter, Depends
import random
from geopack import geopack, t89
import datetime

router = APIRouter(
    prefix="/jwu_test",
    tags=["jwu_test"],
    responses={404: {"description": "Not found"}},
)

@router.get("/get_numbers")
def get_two_numbers():
    a = random.randint(1,10)
    b = random.randint(1,100)
    return (a,b)

@router.get("/get_igrf")
def get_igrf(xgsm: float, ygsm: float, zgsm: float):
    """
        Input
        xgsm,ygsm,zgsm: The given position in cartesian GSM coordinate in Re (earth radii, 1 Re = 6371.2 km).
        Return
        bxgsm,bygsm,bzgsm: Cartesian GSM components of the internal magnetic field in nT.
    """
    t1 = datetime.datetime(2022,1,14,17,17,58)
    t0 = datetime.datetime(1970,1,1)
    ut = (t1-t0).total_seconds()
    tilt = geopack.recalc(ut)
    Re = 6371.2  # in km
    xgsm = xgsm/Re
    ygsm = ygsm/Re
    zgsm = zgsm/Re
    bxgsm,bygsm,bzgsm = geopack.igrf_gsm(xgsm,ygsm,zgsm)
    return (bxgsm, bygsm, bzgsm)
