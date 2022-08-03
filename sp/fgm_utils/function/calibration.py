from typing import List
import numpy as np
from scipy.optimize import curve_fit
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import lsqr
from .. import parameter 
from . import Bplot

def linear_fit(x, m, c):
    return m * x + c

def quad_fit(x, a, b, c):
    return a * x**2 + b * x + c

def cube_fit(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

def sine_fit(x, alpha, A, w, p, k):
    return alpha * A * np.sin(w * x + p) + k

def cosine_fit(x, alpha, A, w, p, k):
    return alpha * A * np.cos(w * x + p) + k

def mask_neg(L):
    L_mask = np.zeros(len(L))
    for i in range(len(L)):
        if L[i] > 0:
            L_mask[i] = L[i]
        else:
            L_mask[i] = np.nan
    return L_mask


def moving_average(time, signal, T=10, func=quad_fit):
    signal_avg = np.zeros(len(signal))
    for i in range(len(signal)):
        if time[i] - time[0] < T / 2:
            low_lim = time[i]
            high_lim = time[i] + T
        elif time[-1] - time[i] < T / 2:
            low_lim = time[i] - T
            high_lim = time[i]
        else:
            low_lim = time[i] - T / 2
            high_lim = time[i] + T / 2
        idx = (time >= low_lim) & (time <= high_lim)
        fit_opt, fit_covar = curve_fit(func, time[idx], signal[idx])
        signal_avg[i] = func(time[i], *fit_opt)
    return signal_avg


def running_spline(time, sample_time, sample_signal, T=10, func=cube_fit):
    signal_fit = np.zeros(len(time))
    for i in range(len(time)):
        if time[i] - time[0] < T / 2:
            low_lim = time[0]
            high_lim = time[0] + T
        elif time[-1] - time[i] < T / 2:
            low_lim = time[-1] - T
            high_lim = time[-1]
        else:
            low_lim = time[i] - T / 2
            high_lim = time[i] + T / 2

        idx_sample = (sample_time >= low_lim) & (sample_time <= high_lim)
        fit_opt, fit_covar = curve_fit(
            func, sample_time[idx_sample], sample_signal[idx_sample]
        )
        signal_fit[i] = func(time[i], *fit_opt)
    return signal_fit


def running_filter(time, signal, eps=1, T=20):
    filter_idx = []
    for i in range(len(time)):
        if time[i] - time[0] < T / 2:
            low_lim = time[0]
            high_lim = time[0] + T
        elif time[-1] - time[i] < T / 2:
            low_lim = time[-1] - T
            high_lim = time[-1]
        else:
            low_lim = time[i] - T / 2
            high_lim = time[i] + T / 2
        idx = (time >= low_lim) & (time <= high_lim)
        med = np.median(signal[idx])
        std = np.std(signal[idx])
        if signal[i] >= med - eps * std and signal[i] <= med + eps * std:
            filter_idx.append(i)
    return np.array(filter_idx)


def calib_leastsquare(
    B_S1_corr, B_S2_corr, B_S3_corr, fgs_igrf_fgm_x, fgs_igrf_fgm_y, fgs_igrf_fgm_z
):

    """use B igrf to calibrate fgs data in fgm coordinate 
       least sqaure fitting 
    """
    n = len(B_S1_corr)
    b = np.concatenate((B_S1_corr, B_S2_corr, B_S3_corr))
    A = np.zeros((3 * n, 12))
    A[0:n, 0] = fgs_igrf_fgm_x
    A[0:n, 1] = fgs_igrf_fgm_y
    A[0:n, 2] = fgs_igrf_fgm_z
    A[0:n, 3] = np.ones(n)
    A[n : 2 * n, 4] = fgs_igrf_fgm_x
    A[n : 2 * n, 5] = fgs_igrf_fgm_y
    A[n : 2 * n, 6] = fgs_igrf_fgm_z
    A[n : 2 * n, 7] = np.ones(n)
    A[2 * n : 3 * n, 8] = fgs_igrf_fgm_x
    A[2 * n : 3 * n, 9] = fgs_igrf_fgm_y
    A[2 * n : 3 * n, 10] = fgs_igrf_fgm_z
    A[2 * n : 3 * n, 11] = np.ones(n)
    A = csc_matrix(A)
    x = lsqr(A, b, atol=1e-10, btol=1e-10)[0]

    orth = np.array([[x[0], x[1], x[2]], [x[4], x[5], x[6]], [x[8], x[9], x[10]]])
    offsets = np.array([x[3], x[7], x[11]])
    calib = np.linalg.inv(orth)

    B_S1_calib = (
        calib[0, 0] * (B_S1_corr - offsets[0])
        + calib[0, 1] * (B_S2_corr - offsets[1])
        + calib[0, 2] * (B_S3_corr - offsets[2])
    )
    B_S2_calib = (
        calib[1, 0] * (B_S1_corr - offsets[0])
        + calib[1, 1] * (B_S2_corr - offsets[1])
        + calib[1, 2] * (B_S3_corr - offsets[2])
    )
    B_S3_calib = (
        calib[2, 0] * (B_S1_corr - offsets[0])
        + calib[2, 1] * (B_S2_corr - offsets[1])
        + calib[2, 2] * (B_S3_corr - offsets[2])
    ) 

    return [B_S1_calib, B_S2_calib, B_S3_calib]   

def ctime_calib(ctime, B_x = None, B_y = None, B_z = None):
    """degap when ctime has large jumps
       multipolar jumps: 
       unipolar jumps:  
    """
    delta_t = np.median(ctime[1:]-ctime[:-1])
    ctime_adj = ctime[1:]-ctime[:-1] - delta_t
    ctime_idx = np.where(np.abs(ctime_adj) > parameter.ctime_thrhld)[0]
    Bplot.ctimediff_plot(ctime, ctime_idx)

    err_idx_1 = []
    err_idx_2 = []
    err_num = 0

    i = 0
    while i < len(ctime_idx):
        if i == len(ctime_idx)-1:
            err_idx_1.append(ctime_idx[i])
            err_idx_2.append(0)
            err_num += 1

            i += 1
        else:         
            if (ctime_adj[ctime_idx[i]]*ctime_adj[ctime_idx[i+1]] < 0 # multipolar jumps
                and np.abs(np.abs(ctime_adj[ctime_idx[i]]) - np.abs(ctime_adj[ctime_idx[i+1]])) < 1e-3): # pos and neg jumps similar magnitude
                err_idx_1.append(ctime_idx[i])
                err_idx_2.append(ctime_idx[i+1])
                err_num += 1
                
                mid_idx = [i for i in range(ctime_idx[i]+1, ctime_idx[i+1])]
                ctimediff_mid = np.median(ctime[mid_idx][1:]-ctime[mid_idx][:-1])
                shift_1 = np.abs((ctime[1:]-ctime[:-1])[ctime_idx[i]] - ctimediff_mid)
                shift_2 = np.abs((ctime[1:]-ctime[:-1])[ctime_idx[i+1]] - ctimediff_mid)
                shift_idx_act = [i for i in range(ctime_idx[i]+1, ctime_idx[i+1]+1)]
                shift = np.mean([shift_1, shift_2])
                ctime[shift_idx_act] -= shift

                i += 2
            else: # unipolar jumps
                err_idx_1.append(ctime_idx[i])
                err_idx_2.append(0)
                err_num += 1

                if parameter.ctime_correct == True:
                    # TODO: not correct for unipolar spike
                    """shift all the way to the end
                    if np.abs(ctime_adj[ctime_idx[i]]) < 0.1:   # small jumps not gaps
                        mid_idx = [i for i in range(ctime_idx[i]+1, len(ctime))]
                        shift = np.median(np.abs((ctime[1:]-ctime[:-1])[ctime_idx[i]] - (ctime[1:]-ctime[:-1])[ctime_idx[i]-1]))
                        shift_idx_act = [i for i in range(ctime_idx[i]+1, len(ctime))]
                        ctime[shift_idx_act] -= shift
                    """
                    """
                    if np.abs(ctime_adj[ctime_idx[i]]) < 0.1:   # small jumps not gaps
                        fit_range = range(ctime_idx[i]-20, ctime_idx[i]+20)
                        B_x_fit = sine_fit(
                            fit_range, *curve_fit(
                            sine_fit, fit_range, B_x[fit_range])[0])
                        B_y_fit = sine_fit(
                            fit_range, *curve_fit(
                            sine_fit, fit_range, B_y[fit_range])[0])
                        B_z_fit = sine_fit(
                            fit_range, *curve_fit(
                            sine_fit, fit_range, B_z[fit_range])[0])  
                    """
                i += 1

    #Bplot.ctimediff_plot(ctime, ctime_idx)
    return ctime, ctime_idx    