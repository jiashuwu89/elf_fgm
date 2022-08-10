from geopack import geopack
import datetime

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