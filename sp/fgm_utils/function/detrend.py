from typing import List
import numpy as np
from scipy.optimize import curve_fit
from .. import parameter 
from . import calibration
from . import Bplot

def detrend_linear(
    ctime: List[float], B_x: List[float], B_y: List[float], B_z: List[float]
    ):
    """detrend with linear fit 

    """
    B_x_trend = calibration.linear_fit(
                ctime,
                *curve_fit(calibration.linear_fit, ctime, B_x)[0],
    )

    B_y_trend = calibration.linear_fit(
                ctime,
                *curve_fit(calibration.linear_fit, ctime, B_y)[0],
    )
    
    B_z_trend = calibration.linear_fit(
                ctime,
                *curve_fit(calibration.linear_fit, ctime, B_z)[0],
    )

    #Bplot.B2_ctime_plot(ctime, B_x, B_y, B_z, B_x_trend, B_y_trend, B_z_trend, "res_dmxl and trend_dmxl")    
    return [B_x_trend, B_y_trend, B_z_trend]


def detrend_linear_2point(
    ctime: List[float], B_x: List[float], B_y: List[float], B_z: List[float]
    ):
    """detrend with a linear trend with only first and last two points
    one problem is sometimes start and end has large spikes need to be removed
    check derivative of first and last three points
    """
    trend = np.zeros((3, len(ctime)))
    B_all = [B_x, B_y, B_z]
    for B_i in range(3):
        B = B_all[B_i]
        d_B = np.gradient(B) / np.gradient(ctime)
        index1 = [i for i in d_B[0:3] if np.abs(i) > parameter.fsp_detrend_cutoff*np.average(np.abs(d_B))]
        index2 = [i for i in d_B[-4:-1] if np.abs(i) > parameter.fsp_detrend_cutoff*np.average(np.abs(d_B))]
        if (not index1 and not index2):
            trend[B_i, :] = B[0] + (ctime-ctime[0])*(B[-1] - B[0])/(ctime[-1]-ctime[0])
        elif (not index1 and index2):
            trend[B_i, :] = B[0] + (ctime-ctime[0])*(B[-4] - B[0])/(ctime[-4]-ctime[0])
        elif (index1 and not index2):   
            trend[B_i, :] = B[3] + (ctime-ctime[3])*(B[-1] - B[3])/(ctime[-1]-ctime[3])
        else:
            trend[B_i, :] = B[3] + (ctime-ctime[3])*(B[-4] - B[3])/(ctime[-4]-ctime[3])

    #Bplot.B2_ctime_plot(ctime, B_x, B_y, B_z, B_x_trend, B_y_trend, B_z_trend, "res_dmxl and trend_dmxl")    
    return [trend[0,:], trend[1,:], trend[2,:]]


def detrend_quad_log(
    ctime: List[float], B_x: List[float], B_y: List[float], B_z: List[float]
    ):
    """detrend with quadratic fit 

    #TODO: verify fit on log scale

    """
    x_min = np.abs(np.min(B_x))+1
    y_min = np.abs(np.min(B_y))+1
    z_min = np.abs(np.min(B_z))+1
    B_x = B_x + x_min
    B_y = B_y + y_min
    B_z = B_z + z_min

    B_x_trend = calibration.quad_fit(
                ctime,
                *curve_fit(calibration.quad_fit, ctime, np.log10(B_x))[0],
    )
    B_x_trend = 10**B_x_trend - x_min

    B_y_trend = calibration.quad_fit(
                ctime,
                *curve_fit(calibration.quad_fit, ctime, np.log10(B_y))[0],
    )
    B_y_trend = 10**B_y_trend - y_min
    
    B_z_trend = calibration.quad_fit(
                ctime,
                *curve_fit(calibration.quad_fit, ctime, np.log10(B_z))[0],
    )
    B_z_trend = 10**B_z_trend - z_min

    #Bplot.B2_ctime_plot(ctime, B_x, B_y, B_z, B_x_trend, B_y_trend, B_z_trend, "res_dmxl and trend_dmxl")
    return [B_x_trend, B_y_trend, B_z_trend]


def detrend_quad(
    ctime: List[float], B_x: List[float], B_y: List[float], B_z: List[float]
    ):
    """detrend with quadratic fit 

    """
    B_x_trend = calibration.quad_fit(
                ctime,
                *curve_fit(calibration.quad_fit, ctime, B_x)[0],
    )

    B_y_trend = calibration.quad_fit(
                ctime,
                *curve_fit(calibration.quad_fit, ctime, B_y)[0],
    )
    
    B_z_trend = calibration.quad_fit(
                ctime,
                *curve_fit(calibration.quad_fit, ctime, B_z)[0],
    )

    #Bplot.B2_ctime_plot(ctime, B_x, B_y, B_z, B_x_trend, B_y_trend, B_z_trend, "res_dmxl and trend_dmxl")    
    return [B_x_trend, B_y_trend, B_z_trend]


def detrend_cube(
    ctime: List[float], B_x: List[float], B_y: List[float], B_z: List[float]
    ):
    """detrend with cubic fit 

    """

    B_x_trend = calibration.cube_fit(
                ctime,
                *curve_fit(calibration.cube_fit, ctime, B_x)[0],
    )

    B_y_trend = calibration.cube_fit(
                ctime,
                *curve_fit(calibration.cube_fit, ctime, B_y)[0],
    )
    
    B_z_trend = calibration.cube_fit(
                ctime,
                *curve_fit(calibration.cube_fit, ctime, B_z)[0],
    )
    return [B_x_trend, B_y_trend, B_z_trend]
    #Bplot.B2_ctime_plot(ctime, B_x, B_y, B_z, B_x_trend, B_y_trend, B_z_trend, "res_dmxl and trend_dmxl")    


def del_rogue(ctime: List[float], B_x: List[float], B_y: List[float], B_z: List[float]):

    B = np.sqrt(B_x**2 + B_y**2 + B_z**2)
    dB = np.gradient(B) / np.gradient(ctime)
    dB_ave = np.average(dB)
    dB_std = np.std(dB)
    dBx = np.gradient(B_x) / np.gradient(ctime)
    dBx_ave = np.average(dBx)
    dBx_std = np.std(dBx)
    dBy = np.gradient(B_y) / np.gradient(ctime)
    dBy_ave = np.average(dBy)
    dBy_std = np.std(dBy)
    dBz = np.gradient(B_z) / np.gradient(ctime)
    dBz_ave = np.average(dBz)
    dBz_std = np.std(dBz)
    #Bplot.B_ctime_plot_single(ctime, dB)
    #Bplot.B_ctime_plot(ctime, dBx, dBy, dBz)

    index = [*range(10)] + [*range(len(dB)-10, len(dB))] if np.median(np.diff(ctime)) < 0.15 else [*range(3)] + [*range(len(dB)-3, len(dB))]
    del_index_1 = [
        i for i in index 
        if (
            dB[i] > dB_ave + parameter.eps_rogue * dB_std or 
            dB[i] < dB_ave - parameter.eps_rogue * dB_std)
        ]

    del_index_3 = [
        i for i in index 
        if (
            (dBx[i] > dBx_ave + parameter.eps_rogue * dBx_std or dBx[i] < dBx_ave - parameter.eps_rogue * dBx_std) and
            (dBy[i] > dBy_ave + parameter.eps_rogue * dBy_std or dBy[i] < dBy_ave - parameter.eps_rogue * dBy_std) and 
            (dBz[i] > dBz_ave + parameter.eps_rogue * dBz_std or dBz[i] < dBz_ave - parameter.eps_rogue * dBz_std)
            )
        ]
    return np.union1d(np.array(del_index_1, dtype=int),np.array(del_index_3, dtype=int)).tolist()

def delete_data(del_idx, *argv):
    return tuple(np.delete(arg, del_idx, axis = 0) for arg in argv)
