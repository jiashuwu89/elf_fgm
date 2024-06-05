from typing import List, Optional
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
    ctime: List[float], 
    B_x: Optional[List[float]] = None, 
    B_y: Optional[List[float]] = None, 
    B_z: Optional[List[float]] = None,
    inlier_idx_x: Optional[List[int]] = None, 
    inlier_idx_y: Optional[List[int]] = None, 
    inlier_idx_z: Optional[List[int]] = None
    ):
    """detrend with quadratic fit 
        if outlinear is not None, then remove outliner of z component
    """
    def detrend_component(component, inlier_idx):
        if component is None:
            return None
        if inlier_idx is None:
            trend = calibration.quad_fit(
                ctime,
                *curve_fit(calibration.quad_fit, ctime, component)[0],
            )
        else:
            trend = calibration.quad_fit(
                ctime,
                *curve_fit(calibration.quad_fit, ctime[inlier_idx], component[inlier_idx])[0],
            )
        return trend
    
    B_x_trend = detrend_component(B_x, inlier_idx_x)
    B_y_trend = detrend_component(B_y, inlier_idx_y)
    B_z_trend = detrend_component(B_z, inlier_idx_z)
    
    #Bplot.B2_ctime_plot(ctime, B_x, B_y, B_z, B_x_trend, B_y_trend, B_z_trend, "res_dmxl and trend_dmxl")    
    return [B_x_trend, B_y_trend, B_z_trend]


def detrend_quadcube(
    ctime: List[float], 
    B_x: Optional[List[float]] = None, 
    B_y: Optional[List[float]] = None, 
    B_z: Optional[List[float]] = None,
    inlier_idx_x: Optional[List[int]] = None, 
    inlier_idx_y: Optional[List[int]] = None, 
    inlier_idx_z: Optional[List[int]] = None
    ):
    """detrend with quadratic fit for x and y
       z use cube fit
        if outlinear is not None, then remove outliner of z component
    """
    def detrend_component(component, inlier_idx, fit_func):
        if component is None:
            return None
        if inlier_idx is None:
            trend = fit_func(
                ctime,
                *curve_fit(fit_func, ctime, component)[0],
            )
        else:
            trend = fit_func(
                ctime,
                *curve_fit(fit_func, ctime[inlier_idx], component[inlier_idx])[0],
            )
        return trend
    
    B_x_trend = detrend_component(B_x, inlier_idx_x, calibration.quad_fit)
    B_y_trend = detrend_component(B_y, inlier_idx_y, calibration.quad_fit)
    B_z_trend = detrend_component(B_z, inlier_idx_z, calibration.cube_fit)
    
    return [B_x_trend, B_y_trend, B_z_trend]


def detrend_cube(
    ctime: List[float], 
    B_x: Optional[List[float]] = None, 
    B_y: Optional[List[float]] = None, 
    B_z: Optional[List[float]] = None,
    inlier_idx_x: Optional[List[int]] = None, 
    inlier_idx_y: Optional[List[int]] = None, 
    inlier_idx_z: Optional[List[int]] = None
    ):
    """detrend with cubic fit 

    """
    def detrend_component(component, inlier_idx):
        if component is None:
            return None
        if inlier_idx is None:
            trend = calibration.cube_fit(
                ctime,
                *curve_fit(calibration.cube_fit, ctime, component)[0],
            )
        else:
            trend = calibration.cube_fit(
                ctime,
                *curve_fit(calibration.cube_fit, ctime[inlier_idx], component[inlier_idx])[0],
            )
        return trend
    
    B_x_trend = detrend_component(B_x, inlier_idx_x)
    B_y_trend = detrend_component(B_y, inlier_idx_y)
    B_z_trend = detrend_component(B_z, inlier_idx_z)

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

def remove_outliers(data, sigma = 1):
    """This function removes outliers before detrend
    """
    mean = np.median(np.abs(data))
    std_dev = np.std(np.abs(data))

    lower_bound = mean - (sigma * std_dev)
    upper_bound = mean + (sigma * std_dev)

    filter_idx = (data >= lower_bound ) & (data <= upper_bound)

    return filter_idx


def iter_detrend(ctime, 
                 fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, 
                 fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z,
                 detrend_func,):
    """iteratively determine outliers and fit the baseline
    """
    [fgs_igrf_dmxl_x_detrend, fgs_igrf_dmxl_y_detrend, fgs_igrf_dmxl_z_detrend] = detrend_func(
        ctime,
        B_x = fgs_igrf_dmxl_x, 
        B_y = fgs_igrf_dmxl_y, 
        B_z = fgs_igrf_dmxl_z)
    
    # the first iteration will use difference between ful and igrf to determine outliers
    fgs_res_dmxl_x = fgs_ful_dmxl_x-fgs_igrf_dmxl_x
    fgs_res_dmxl_y = fgs_ful_dmxl_y-fgs_igrf_dmxl_y
    fgs_res_dmxl_z = fgs_ful_dmxl_z-fgs_igrf_dmxl_z
    
    inlier_idx_x = remove_outliers(fgs_res_dmxl_x, sigma=5)
    inlier_idx_y = remove_outliers(fgs_res_dmxl_y, sigma=5)
    inlier_idx_z = remove_outliers(fgs_res_dmxl_z, sigma=2)
    # x and y only exclude outliers once. if iter too many times a lot of points will be excluded, the results will have a large trend
    [fgs_ful_dmxl_x_detrend, fgs_ful_dmxl_y_detrend, fgs_ful_dmxl_z_detrend] = detrend_func(
        ctime,
        B_x = fgs_ful_dmxl_x, 
        B_y = fgs_ful_dmxl_y, 
        B_z = fgs_ful_dmxl_z,
        inlier_idx_x = inlier_idx_x,
        inlier_idx_y = inlier_idx_y,
        inlier_idx_z = inlier_idx_z)
    
    # iter for z
    for i in range(3):
        fgs_res_dmxl_z = fgs_ful_dmxl_z - fgs_ful_dmxl_z_detrend
        inlier_idx_z = remove_outliers(fgs_res_dmxl_z, sigma=2)

        _, _, fgs_ful_dmxl_z_detrend = detrend_func(
            ctime,
            B_z = fgs_ful_dmxl_z,
            inlier_idx_z = inlier_idx_z)

    return [fgs_igrf_dmxl_x_detrend, fgs_igrf_dmxl_y_detrend, fgs_igrf_dmxl_z_detrend, fgs_ful_dmxl_x_detrend, fgs_ful_dmxl_y_detrend, fgs_ful_dmxl_z_detrend]