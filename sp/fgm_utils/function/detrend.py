from typing import List, Optional
import numpy as np
from scipy.optimize import curve_fit
from .. import parameter 
from . import calibration
from . import Bplot
import datetime


def detrend_linear(
    ctime: List[float], 
    B_x: Optional[List[float]] = None, 
    B_y: Optional[List[float]] = None, 
    B_z: Optional[List[float]] = None,
    inlier_idx_x: Optional[List[int]] = None, 
    inlier_idx_y: Optional[List[int]] = None, 
    inlier_idx_z: Optional[List[int]] = None
    ):
    """detrend with linear fit 

    """
    def linear_fit_component(component, inlier_idx):
        if component is None:
            return None
        if inlier_idx is None:
            trend = calibration.linear_fit(
                ctime,
                *curve_fit(calibration.linear_fit, ctime, component)[0],
            )
        else:
            trend = calibration.linear_fit(
                ctime,
                *curve_fit(calibration.linear_fit, ctime[inlier_idx], component[inlier_idx])[0],
            )
        return trend

    B_x_trend = linear_fit_component(B_x, inlier_idx_x)
    B_y_trend = linear_fit_component(B_y, inlier_idx_y)
    B_z_trend = linear_fit_component(B_z, inlier_idx_z)
    
    return [B_x_trend, B_y_trend, B_z_trend]


def detrend_linear_2point(
    ctime: List[float], 
    B_x: Optional[List[float]] = None, 
    B_y: Optional[List[float]] = None,
    B_z: Optional[List[float]] = None,
    ):
    """detrend with a linear trend with only first and last two points
    one problem is sometimes start and end has large spikes need to be removed
    check derivative of first and last three points
    """
    def detrend_linear_2point_component(B):
        if B is None:
            return None
        d_B = np.gradient(B) / np.gradient(ctime)
        index1 = [i for i in d_B[0:3] if np.abs(i) > parameter.fsp_detrend_cutoff*np.average(np.abs(d_B))]
        index2 = [i for i in d_B[-4:-1] if np.abs(i) > parameter.fsp_detrend_cutoff*np.average(np.abs(d_B))]
        if (not index1 and not index2):
            trend = B[0] + (ctime-ctime[0])*(B[-1] - B[0])/(ctime[-1]-ctime[0])
        elif (not index1 and index2):
            trend = B[0] + (ctime-ctime[0])*(B[-4] - B[0])/(ctime[-4]-ctime[0])
        elif (index1 and not index2):   
            trend = B[3] + (ctime-ctime[3])*(B[-1] - B[3])/(ctime[-1]-ctime[3])
        else:
            trend = B[3] + (ctime-ctime[3])*(B[-4] - B[3])/(ctime[-4]-ctime[3])

        return trend
   
    B_x_trend = detrend_linear_2point_component(B_x)
    B_y_trend = detrend_linear_2point_component(B_y)
    B_z_trend = detrend_linear_2point_component(B_z)
    
    return [B_x_trend, B_y_trend, B_z_trend]


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
        inlier_idx: index used to fit 
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

    filter_idx = (np.abs(data) >= lower_bound ) & (np.abs(data) <= upper_bound)

    return filter_idx

def select_percentile(data, percent = 50):
    """This function select the lower percentile of the data
    """
    threshold = np.percentile(data, percent)
    low_idx = data <= threshold

    return low_idx

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


def iter_detrend_xyz(cross_times, fgs_res_dmxl_x, 
                 fgs_res_dmxl_y, fgs_res_dmxl_z, detrend_func, detrend_method=3, detrend_percent=85, detrend_portion=None):
    """iteratively determine outliers and fit the baseline, try both quad and cube fit, pick the one with smaller residual
    Parameter
        - detrend_func: the function for xy detrend and initial z detrend
        - detrend_method: detrend method for z detrend
    """

    ## detrend x, y
    # the first iteration will use difference between ful and igrf to determine outliers
    inlier_idx_x = remove_outliers(fgs_res_dmxl_x, sigma=3)
    inlier_idx_y = remove_outliers(fgs_res_dmxl_y, sigma=3)
    if detrend_portion is None:
        inlier_idx_z = remove_outliers(fgs_res_dmxl_z, sigma=5)
    else:
        n = len(fgs_res_dmxl_z)
        inlier_idx_z = np.zeros(n, dtype=bool)
        for start, end in detrend_portion:
            start_idx = int(np.floor(start * n))
            end_idx = int(np.ceil(end * n))
            inlier_idx_z[start_idx:end_idx] = True
 
    # x and y only exclude outliers once. if iter too many times a lot of points will be excluded, the results will have a large trend
    [fsp_trend_x_linear, fsp_trend_y_linear, fsp_trend_z_linear] = detrend_func(
        cross_times,
        B_x = fgs_res_dmxl_x, 
        B_y = fgs_res_dmxl_y,
        B_z = fgs_res_dmxl_z,
        inlier_idx_x = inlier_idx_x,
        inlier_idx_y = inlier_idx_y,
        inlier_idx_z = inlier_idx_z,
    )

    if parameter.makeplot == True:
        Bplot.B_ctime_plot(
            cross_times, [fgs_res_dmxl_z, fsp_trend_x_linear], [fgs_res_dmxl_z, fsp_trend_y_linear], 
            [fgs_res_dmxl_z, fsp_trend_z_linear], cross_times=cross_times[~inlier_idx_z], scatter=True)
        

    fgs_res_dmxl_z = fgs_res_dmxl_z - fsp_trend_z_linear

    iter = 1 if detrend_percent == 100 else 3

    ## detrend z
    if detrend_method == 1:
        # exclude points according to gradient, and fit quad and cube
        low_idxs = np.ones(len(fgs_res_dmxl_z), dtype=bool)
        for i in range(iter):
            dy = np.diff(fgs_res_dmxl_z, i)
            gradients = np.abs(dy) 
            gradients = np.pad(gradients, (0, i), mode='edge')
            low_idx = select_percentile(gradients, percent=detrend_percent)
            low_idxs = np.logical_and(low_idxs, low_idx) # select the lower 85% in each iteration. and combine

        # fit the one with smaller resdidual
        _, _, fsp_trend_z_quad = detrend_quad(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)
        _, _, fsp_trend_z_cube = detrend_cube(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)

        fsp_trend_z_quad_res = np.abs(fgs_res_dmxl_z - fsp_trend_z_quad)
        fsp_trend_z_cube_res = np.abs(fgs_res_dmxl_z - fsp_trend_z_cube)

        if np.mean(fsp_trend_z_quad_res) < np.mean(fsp_trend_z_cube_res):
            fsp_trend_z = fsp_trend_z_quad
            low_idxs_final = low_idxs
        else:
            fsp_trend_z = fsp_trend_z_cube
            low_idxs_final = low_idxs

    elif detrend_method == 2:
        # exclude points according to magnitude of res, and fit quad and cube
        low_idxs = np.ones(len(fgs_res_dmxl_z), dtype=bool)
        # fit the one with smaller resdidual
        _, _, fsp_trend_z_quad = detrend_quad(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)
        _, _, fsp_trend_z_cube = detrend_cube(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)

        fsp_trend_z_quad_res = (fgs_res_dmxl_z - fsp_trend_z_quad)**2
        fsp_trend_z_cube_res = (fgs_res_dmxl_z - fsp_trend_z_cube)**2

        for i in range(iter):
            low_idx_quad = select_percentile(fsp_trend_z_quad_res, percent=detrend_percent)
            low_idxs_quad = np.logical_and(low_idx_quad, low_idxs)
            _, _, fsp_trend_z_quad = detrend_quad(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs_quad)
            fsp_trend_z_quad_res = (fgs_res_dmxl_z - fsp_trend_z_quad)**2

        for i in range(iter):
            low_idx_cube = select_percentile(fsp_trend_z_cube_res, percent=detrend_percent)
            low_idxs_cube = np.logical_and(low_idx_cube, low_idxs)
            _, _, fsp_trend_z_cube = detrend_cube(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs_cube)
            fsp_trend_z_cube_res = (fgs_res_dmxl_z - fsp_trend_z_cube)**2

        if np.mean(fsp_trend_z_quad_res) < np.mean(fsp_trend_z_cube_res):
            fsp_trend_z = fsp_trend_z_quad 
            low_idxs_final = low_idxs_quad
        else:
            fsp_trend_z = fsp_trend_z_cube
            low_idxs_final = low_idxs_cube

    elif detrend_method == 3: 
        # method 1  
        low_idxs = np.ones(len(fgs_res_dmxl_z), dtype=bool)
        for i in range(iter):
            dy = np.diff(fgs_res_dmxl_z, i)
            gradients = np.abs(dy) 
            gradients = np.pad(gradients, (0, i), mode='edge')
            low_idx = select_percentile(gradients, percent=detrend_percent)
            low_idxs = np.logical_and(low_idxs, low_idx) # select the lower 85% in each iteration. and combine

        # fit the one with smaller resdidual
        _, _, fsp_trend_z_quad = detrend_quad(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)
        _, _, fsp_trend_z_cube = detrend_cube(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)

        fsp_trend_z_quad_res = np.abs(fgs_res_dmxl_z - fsp_trend_z_quad)
        fsp_trend_z_cube_res = np.abs(fgs_res_dmxl_z - fsp_trend_z_cube)

        thr = 50 ## this is the threshold to check how many data points in residual are below this threshold. 
        # choose the method with most data points below this threshold 

        if np.sum(fsp_trend_z_quad_res < thr) > np.sum(fsp_trend_z_cube_res < thr):
            fsp_res_z_1 = np.sum(fsp_trend_z_quad_res < thr) 
            fsp_trend_z_1 = fsp_trend_z_quad
            low_idxs_1 = low_idxs.copy()
        else:
            fsp_res_z_1 = np.sum(fsp_trend_z_cube_res < thr)
            fsp_trend_z_1 = fsp_trend_z_cube
            low_idxs_1 = low_idxs.copy()

        # method 2
        low_idxs = np.ones(len(fgs_res_dmxl_z), dtype=bool)
        _, _, fsp_trend_z_quad = detrend_quad(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)
        _, _, fsp_trend_z_cube = detrend_cube(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)

        fsp_trend_z_quad_res = np.abs(fgs_res_dmxl_z - fsp_trend_z_quad)
        fsp_trend_z_cube_res = np.abs(fgs_res_dmxl_z - fsp_trend_z_cube)

        for i in range(iter):
            low_idx_quad = select_percentile(fsp_trend_z_quad_res, percent=detrend_percent)
            low_idxs_quad = np.logical_and(low_idx_quad, low_idxs)
            _, _, fsp_trend_z_quad = detrend_quad(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs_quad)
            fsp_trend_z_quad_res = np.abs(fgs_res_dmxl_z - fsp_trend_z_quad)

        for i in range(iter):
            low_idx_cube = select_percentile(fsp_trend_z_cube_res, percent=detrend_percent)
            low_idxs_cube = np.logical_and(low_idx_cube, low_idxs)
            _, _, fsp_trend_z_cube = detrend_cube(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs_cube)
            fsp_trend_z_cube_res = np.abs(fgs_res_dmxl_z - fsp_trend_z_cube)
 
        if np.sum(fsp_trend_z_quad_res < thr) > np.sum(fsp_trend_z_cube_res < thr):
            fsp_res_z_2 = np.sum(fsp_trend_z_quad_res < thr) 
            fsp_trend_z_2 = fsp_trend_z_quad 
            low_idxs_2 = low_idxs_quad
        else:
            fsp_res_z_2 = np.sum(fsp_trend_z_cube_res < thr) 
            fsp_trend_z_2 = fsp_trend_z_cube
            low_idxs_2 = low_idxs_cube

        if fsp_res_z_1 > fsp_res_z_2:
            fsp_trend_z =  fsp_trend_z_1
            low_idxs_final = low_idxs_1
        else:
            fsp_trend_z =  fsp_trend_z_2
            low_idxs_final = low_idxs_2

    elif detrend_method == 4: 
        # exclude points according to gradient, fit quad only quad
        if detrend_portion is None:
            low_idxs = np.ones(len(fgs_res_dmxl_z), dtype=bool)
        else:
            n = len(fgs_res_dmxl_z)
            low_idxs = np.zeros(n, dtype=bool)
            for start, end in detrend_portion:
                start_idx = int(np.floor(start * n))
                end_idx = int(np.ceil(end * n))
                low_idxs[start_idx:end_idx] = True

        for i in range(iter):
            dy = np.diff(fgs_res_dmxl_z, i)
            gradients = np.abs(dy) 
            gradients = np.pad(gradients, (0, i), mode='edge')
            low_idx = select_percentile(gradients, percent=detrend_percent)
            low_idxs = np.logical_and(low_idxs, low_idx) # select the lower 85% in each iteration. and combine

        # fit the one with smaller resdidual
        _, _, fsp_trend_z_quad = detrend_quad(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)

        fsp_trend_z_quad_res = np.abs(fgs_res_dmxl_z - fsp_trend_z_quad)

        fsp_trend_z = fsp_trend_z_quad
        low_idxs_final = low_idxs

    elif detrend_method == 5: 
        # exclude points according to magnitude of res, fit quad only quad
        if detrend_portion is None:
            low_idxs = np.ones(len(fgs_res_dmxl_z), dtype=bool)
        else:
            n = len(fgs_res_dmxl_z)
            low_idxs = np.zeros(n, dtype=bool)
            for start, end in detrend_portion:
                start_idx = int(np.floor(start * n))
                end_idx = int(np.ceil(end * n))
                low_idxs[start_idx:end_idx] = True
    
        # x and y only exclude outliers once. if iter too many times a lot of points will be excluded, the results will have a large trend
        _, _, fsp_trend_z_quad = detrend_quad(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)
        fsp_trend_z_quad_res = (fgs_res_dmxl_z - fsp_trend_z_quad)**2
        
        for i in range(iter):
            low_idx_quad = select_percentile(fsp_trend_z_quad_res, percent=detrend_percent)
            low_idxs_quad = np.logical_and(low_idx_quad, low_idxs)
            _, _, fsp_trend_z_quad = detrend_quad(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs_quad)
            fsp_trend_z_quad_res = (fgs_res_dmxl_z - fsp_trend_z_quad)**2

        fsp_trend_z = fsp_trend_z_quad 
        low_idxs_final = low_idxs_quad

    elif detrend_method == 6: 
        # exclude points according to magnitude of res, fit linear only
        low_idxs = np.ones(len(fgs_res_dmxl_z), dtype=bool)
        _, _, fsp_trend_z_quad = detrend_linear(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs)

        fsp_trend_z_quad_res = (fgs_res_dmxl_z - fsp_trend_z_quad)**2

        for i in range(iter):
            low_idx_quad = select_percentile(fsp_trend_z_quad_res, percent=detrend_percent)
            low_idxs_quad = np.logical_and(low_idx_quad, low_idxs)
            _, _, fsp_trend_z_quad = detrend_linear(cross_times, B_z = fgs_res_dmxl_z, inlier_idx_z = low_idxs_quad)
            fsp_trend_z_quad_res = (fgs_res_dmxl_z - fsp_trend_z_quad)**2

        fsp_trend_z = fsp_trend_z_quad 
        low_idxs_final = low_idxs_quad

    elif detrend_method == 7: 
        # use the two  end points only to detrend
        _, _, fsp_trend_z = detrend_linear_2point(cross_times, B_z = fgs_res_dmxl_z)
        low_idxs_final = np.ones(len(fgs_res_dmxl_z), dtype=bool)

   
    if parameter.makeplot == True:
        Bplot.B_ctime_plot(
            cross_times, [fgs_res_dmxl_z, fsp_trend_z], [fgs_res_dmxl_z, fsp_trend_z], 
            [fgs_res_dmxl_z, fsp_trend_z], cross_times=cross_times[~low_idxs_final], scatter=True)

    return fsp_trend_x_linear, fsp_trend_y_linear, fsp_trend_z + fsp_trend_z_linear


detrend_list  = {
    "2022-04-02/18:20:00": { #80
        'method': 5,
        'mission': 'elb',
        'percent': 80,
        'portion': None,
    },
    "2022-04-06/12:15:00": { #89
        'method': 2,
        'mission': 'ela',
        'percent': 80,
        'portion': None,
    },
    "2022-04-03/17:32:00": { #83
        'method': 5,
        'mission': 'elb',
        'percent': 100,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-04-02/15:50:00": { #79
        'method': 3,
        'mission': 'ela',
        'percent': 85,
        'portion': None,
    },
    "2022-04-04/18:40:00": { #86
        'method': 3,
        'mission': 'ela',
        'percent': 85,
        'portion': None,
    },
    "2022-04-03/20:35:00": { #84
        'method': 5,
        'mission': 'elb',
        'percent': 100,
        'portion': None,
    },
    "2022-04-25/07:45:00": { #106
        'method': 3,
        'mission': 'ela',
        'percent': 85,
        'portion': None,
    },
    "2022-04-01/17:35:00": { #75
        'method': 4,
        'mission': 'elb',
        'percent': 95,
        'portion': None,
    },
    "2022-04-01/19:10:00": { #76
        'method': 4,
        'mission': 'elb',
        'percent': 98,
        'portion': None,
    },
    "2022-04-02/10:40:00": { #82
        'method': 5,
        'mission': 'elb',
        'percent': 85,
        'portion': None,
    },
    "2022-04-02/21:35:00": { #77
        'method': 5,
        'mission': 'elb',
        'percent': 75,
        'portion': None,
    },
    "2021-04-13/06:15:00": { #128
        'method': 2,
        'mission': 'ela',
        'percent': 80,
        'portion': None,
    },
    "2022-01-28/06:20:00": { #135
        'method': 5,
        'mission': 'ela',
        'percent': 80,
        'portion': None,
    },
    "2022-02-14/16:15:00": { #136
        'method': 5,
        'mission': 'ela',
        'percent': 75,
        'portion': None,
    },
    "2022-04-14/03:55:00": { #137
        'method': 5,
        'mission': 'ela',
        'percent': 75,
        'portion': None,
    },
    "2022-04-23/11:05:00": { #139
        'method': 5,
        'mission': 'ela',
        'percent': 70,
        'portion': None,
    },
    "2022-05-12/06:25:00": { #140
        'method': 6,
        'mission': 'ela',
        'percent': 50,
        'portion': None,
    },
    "2022-05-28/01:00:00": { #142
        'method': 6,
        'mission': 'ela',
        'percent': 50,
        'portion': None,
    },
    "2022-02-06/18:04:00": { #153
        'method': 5,
        'mission': 'elb',
        'percent': 80,
        'portion': None,
    },
    "2022-03-06/05:10:00": { #154
        'method': 6,
        'mission': 'elb',
        'percent': 80,
        'portion': None,
    },
    "2022-03-19/19:16:00": { #155
        'method': 4,
        'mission': 'elb',
        'percent': 96,
        'portion': None,
    },
    "2022-04-07/18:58:00": { #157
        'method': 5,
        'mission': 'elb',
        'percent': 70,
        'portion': None,
    },
    "2022-01-17/00:55:00": { #158
        'method': 7,
        'mission': 'elb',
        'percent': 70,
        'portion': None,
    },
    "2022-01-16/15:40:00": { #159
        'method': 5,
        'mission': 'elb',
        'percent': 99,
        'portion': None,
    },
    "2022-01-16/22:50:00": { #160
        'method': 7,
        'mission': 'elb',
        'percent': 99,
        'portion': None,
    },
    "2022-01-15/18:58:00": { #161
        'method': 7,
        'mission': 'elb',
        'percent': 99,
        'portion': None,
    },
    "2022-08-07/18:50:00": { #164
        'method': 4,
        'mission': 'elb',
        'percent': 60,
        'portion': [(0, 0.35),(0.98,1)],
    },
    "2022-03-11/18:54:00": { #168
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.6)],
    },
    "2022-03-15/18:40:00": { #169
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.5),(0.95, 1)],
    },
    "2022-03-18/02:16:00": { #170
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.5),(0.95, 1)],
    },
    "2022-03-19/10:46:00": { #171
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.5),(0.95, 1)],
    },
    "2022-03-21/12:10:00": { #172
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.5),(0.95, 1)],
    },
    "2022-03-22/05:20:00": { #173
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2022-03-24/03:26:00": { #174
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.5),(0.98, 1)],
    },
    "2022-03-29/19:12:00": { #175
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.5),(0.98, 1)],
    },
    "2022-03-30/19:52:00": { #176
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.5),(0.98, 1)],
    },
    "2022-04-05/03:58:00": { #177
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.5),(0.98, 1)],
    },
    "2022-04-11/12:34:00": { #178
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-04-17/16:26:00": { #179
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-04-18/12:30:00": { #180
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-04-19/04:06:00": { #181
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2022-04-20/12:16:00": { #182
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-04-22/16:34:00": { #183
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-04-23/11:06:00": { #184
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-04-29/14:48:00": { #185
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.4),(0.98, 1)],
    },
    "2022-04-30/12:22:00": { #186
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-05-03/12:38:00": { #187
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-05-17/02:56:00": { #188
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2022-05-11/21:24:00": { #189
        'method': 5,
        'mission': 'ela',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2020-08-27/21:42:00": { #190
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.5),(0.98, 1)],
    },
    "2021-03-07/18:50:00": { #191
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    }, 
    "2021-03-08/19:44:00": { #192
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },  
    "2021-03-08/22:54:00": { #193
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    }, 
    "2021-03-09/00:28:00": { #194
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2021-03-13/04:02:00": { #195
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2021-06-19/13:04:00": { #196
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2022-01-03/05:26:00": { #197
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.5),(0.98, 1)],
    },
    "2022-01-03/07:00:00": { #198
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.5),(0.98, 1)],
    },
    "2022-01-09/07:06:00": { #199
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-01-09/21:00:00": { #200
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-01-12/18:44:00": { #201
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-01-13/17:58:00": { #202
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-01-14/06:24:00": { #203
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-01-15/22:50:00": { #204
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2022-01-22/18:50:00": { #205
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-01-25/12:08:00": { #206
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2022-02-03/11:22:00": { #207
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2022-02-12/04:06:00": { #208
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-02-14/04:04:00": { #209
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-02-18/15:04:00": { #210
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2022-02-20/19:24:00": { #211
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-02-21/15:34:00": { #212
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-02-22/19:22:00": { #213
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-03-08/18:58:00": { #214
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-03-15/16:24:00": { #215
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-03-16/04:48:00": { #216
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-03-16/09:26:00": { #217
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.4),(0.98, 1)],
    },
    "2022-03-16/20:10:00": { #218
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-03-22/02:56:00": { #219
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-03-23/15:58:00": { #220
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-03-24/04:22:00": { #221
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-03-25/18:56:00": { #222
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-03-27/21:54:00": { #223
        'method': 5,
        'mission': 'elb',
        'percent': 90,
        'portion': [(0, 0.6),(0.98, 1)],
    },
    "2022-03-29/04:50:00": { #224
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-03-29/17:10:00": { #225
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.4), (0.98, 1)],
    },
    "2022-03-30/11:44:00": { #226
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-04-02/13:48:00": { #227
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-04-05/18:54:00": { #228
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-04-06/04:26:00": { #229
        "method": 5,
        "mission": "elb",
        "percent": 90,
        'portion': [(0, 0.05),(0.6, 1)],
    },
    "2022-04-06/11:58:00": { #230
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-04-08/00:58:00": { #231
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-04-15/18:00:00": { #232
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.3), (0.98, 1)],
    },
    "2022-04-16/20:12:00": { #233
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.4), (0.98, 1)],
    },
    "2022-04-19/03:58:00": { #234
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.05), (0.6, 1)],
    },
    "2022-04-20/04:38:00": { #235
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.05), (0.6, 1)],
    },
    "2022-04-22/02:40:00": { #236
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-04-24/10:08:00": { #237
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-04-27/15:08:00": { #238
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-04-27/18:10:00": { #239
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-04-28/15:46:00": { #240
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.4), (0.98, 1)],
    },
    "2022-04-28/17:18:00": { #241
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.4), (0.98, 1)],
    },
    "2022-04-29/14:52:00": { #242
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.4), (0.98, 1)],
    },
    "2022-05-02/07:36:00": { #243
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.4), (0.98, 1)],
    },
    "2022-05-19/06:28:00": { #244
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-02-21/19:40:00": { #245
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-01-04/18:36:00": { #246
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.65), (0.98, 1)],
    },
    "2022-01-01/18:36:00": { #247
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-01-15/18:58:00": { #161
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.05), (0.5, 1)],
    },
    "2022-01-08/18:40:00": { #248
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-02-16/04:02:00": { #249
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-01-11/05:40:00": { #250
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-01-13/18:02:00": { #251
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-01-26/21:50:00": { #252
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-02-05/06:16:00": { #253
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-02-10/16:24:00": { #254
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.05), (0.6, 1)],
    },
    "2022-02-11/06:08:00": { #255
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-03-01/05:32:00": { #256
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-03-04/09:30:00": { #257
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.02), (0.4, 1)],
    },
    "2022-03-04/23:08:00": { #258
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-03-11/05:04:00": { #259
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.6), (0.98, 1)],
    },
    "2022-03-14/08:46:00": { #260
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-03-14/10:18:00": { #261
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-03-15/15:36:00": { #262
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-03-21/04:36:00": { #263
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.02), (0.6, 1)],
    },
    "2022-04-30/15:24:00": { #264
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-05-07/02:56:00": { #265
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.02), (0.6, 1)],
    },
    "2022-05-25/23:00:00": { #266
        "method": 5,
        "mission": "ela",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2022-04-06/18:04:00": { #267
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
    "2021-10-03/03:46:00": { #268
        "method": 5,
        "mission": "elb",
        "percent": 90,
        "portion": [(0, 0.5), (0.98, 1)],
    },
}