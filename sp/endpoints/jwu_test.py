
import string
from tracemalloc import start
from fastapi import APIRouter, Query
from geopack import geopack
import random
import datetime
import traceback
import calendar
from pprint import pprint
from cdflib import CDF, cdfepoch
from typing import Union, List
import numpy as np
from scipy.interpolate import interp1d
from pyspedas.cotrans import cotrans_lib 
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import lsqr
from scipy.integrate import trapezoid, simpson
from scipy.optimize import curve_fit
from scipy import signal

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
            print(f"{cdfpath} variables: {vars}")
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

def clip_cdfdata(
    df:pd.Series, 
    starttime:pd.Timestamp, endtime:pd.Timestamp) -> pd.Series:

    startindex = df.index.get_indexer([starttime], method='nearest')[0]
    endindex = df.index.get_indexer([endtime], method='nearest')[0]

    return df[startindex:endindex]
    

def resample_data(
    cur_time:pd.Timestamp, cur_data:pd.Series, 
    target_time:pd.Timestamp) -> pd.Series:

    cur_data_np = np.array(cur_data.to_list())
    dimens = cur_data_np.shape[1]  # x, y, z
    interp_data = np.zeros((len(target_time), dimens))

    x = (cur_time - target_time[0]).total_seconds()
    x_interp = (target_time - target_time[0]).total_seconds()
    for dimen in range(dimens):
        f = interp1d(x, cur_data_np[:,dimen])
        interp_data[:,dimen] = f(x_interp)

    return pd.Series(interp_data.tolist())

def linear_fit(x, m, c):
    return m*x + c

def quad_fit(x, a, b, c):
    return a*x**2 + b*x + c

def cube_fit(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

def sine_fit(x, alpha, A, w, p, k):
    return alpha*A*np.sin(w*x+p) + k

def cosine_fit(x, alpha, A, w, p, k):
    return alpha*A*np.cos(w*x+p) + k

def mask_neg(L):  
    L_mask = np.zeros(len(L))
    for i in range(len(L)):   
        if(L[i]>0):
            L_mask[i] = L[i]
        else:
            L_mask[i] = np.nan        
    return L_mask    


def moving_average(time, signal, T=10, func=quad_fit):    
    signal_avg = np.zeros(len(signal)) 
    for i in range(len(signal)): 
        if(time[i]-time[0]<T/2):
            low_lim = time[i]
            high_lim = time[i]+T
        elif(time[-1]-time[i]<T/2):
            low_lim = time[i]-T
            high_lim = time[i]
        else:   
            low_lim = time[i]-T/2
            high_lim = time[i]+T/2
        idx = (time>=low_lim)&(time<=high_lim)
        fit_opt, fit_covar = curve_fit(func, time[idx], signal[idx])
        signal_avg[i] = func(time[i], *fit_opt)
    return signal_avg


def running_spline(time, sample_time, sample_signal, T=10, func=cube_fit):  
    signal_fit = np.zeros(len(time))   
    for i in range(len(time)):   
        if(time[i]-time[0]<T/2):
            low_lim = time[0]
            high_lim = time[0]+T
        elif(time[-1]-time[i]<T/2):
            low_lim = time[-1]-T
            high_lim = time[-1]
        else:   
            low_lim = time[i]-T/2
            high_lim = time[i]+T/2
        idx = (time>=low_lim)&(time<=high_lim)
        idx_sample = (sample_time>=low_lim)&(sample_time<=high_lim)
        fit_opt, fit_covar = curve_fit(func, sample_time[idx_sample], sample_signal[idx_sample])
        signal_fit[i] = func(time[i], *fit_opt)   
    return signal_fit


def running_filter(time, signal, eps=1, T=20):  
    filter_idx = []    
    for i in range(len(time)):   
        if(time[i]-time[0]<T/2):
            low_lim = time[0]
            high_lim = time[0]+T
        elif(time[-1]-time[i]<T/2):
            low_lim = time[-1]-T
            high_lim = time[-1]
        else:   
            low_lim = time[i]-T/2
            high_lim = time[i]+T/2          
        idx = (time>=low_lim)&(time<=high_lim)        
        med = np.median(signal[idx])
        std = np.std(signal[idx])       
        if (signal[i]>=med-eps*std and signal[i]<=med+eps*std):
            filter_idx.append(i)         
    return np.array(filter_idx)  


@router.get("/fgm_calib")
def fgm_calib(
    starttime_str:str, endtime_str:str,
    sta_cdfpath:str, fgm_cdfpath:str):

    df = pd.DataFrame()

    # time range for sci zone
    starttime = pd.to_datetime(starttime_str)
    endtime = pd.to_datetime(endtime_str)

    # read state cdf for att
    att_cdfdata = pd.DataFrame(get_cdf(sta_cdfpath, vars=['ela_att_time', 'ela_att_gei']))
    att_cdfdata.set_index('ela_att_time', inplace=True)
    att_cdfdata = clip_cdfdata(att_cdfdata, 
        starttime - datetime.timedelta(minutes=2), endtime + datetime.timedelta(minutes=2))

    # read state cdf for pos; not read together with att b/c different length
    pos_cdfdata = pd.DataFrame(get_cdf(sta_cdfpath, vars=['ela_pos_gei']))  # not read state time b/c too slow
    pos_cdfdata['ela_pos_time'] = pd.date_range(start='2022-01-12', periods=len(pos_cdfdata['ela_pos_gei']), freq='S')
    pos_cdfdata.set_index('ela_pos_time', inplace=True)
    pos_cdfdata = clip_cdfdata(pos_cdfdata, 
        starttime - datetime.timedelta(minutes=2), endtime + datetime.timedelta(minutes=2))
 
    # read fgm cdf and clip
    fgm_cdfdata = pd.DataFrame(get_cdf(fgm_cdfpath, vars=['ela_fgs_time', 'ela_fgs']))
    fgm_cdfdata.set_index('ela_fgs_time', inplace=True)
    fgm_cdfdata = clip_cdfdata(fgm_cdfdata, starttime, endtime)

    # resample att and pos to fgm time resolution
    df['att_gei'] = resample_data(att_cdfdata.index, att_cdfdata['ela_att_gei'], fgm_cdfdata.index)
    df['pos_gei'] = resample_data(pos_cdfdata.index, pos_cdfdata['ela_pos_gei'], fgm_cdfdata.index)
    df['fgm_fgm'] = pd.Series(fgm_cdfdata['ela_fgs'].tolist())
    df['time'] = fgm_cdfdata.index
    df['timestamp'] = df['time'].apply(lambda ts: pd.Timestamp(ts).timestamp()) 

    #df.set_index('time', inplace = True)
    #pprint(df.head())

    #iyear, idoy, ih, im, isec = cotrans_lib.get_time_parts(df['timestamp'])
    #print(f"year:{iyear}, doy:{idoy}, h:{ih}, m:{im}, sec:{isec}")

    # coordinate transformation of pos: gei -> gse -> gsm
    df['pos_gse'] = pd.Series(cotrans_lib.subgei2gse(df['timestamp'].tolist(), df['pos_gei'].tolist()).tolist())
    df['pos_gsm'] = pd.Series(cotrans_lib.subgse2gsm(df['timestamp'].tolist(), df['pos_gse'].tolist()).tolist())

    # call igrf b in gsm
    df['igrf_gsm'] = [get_igrf(
        df['time'][i], df['pos_gsm'][i][0], df['pos_gsm'][i][1], df['pos_gsm'][i][2])
        for i in range(len(df['timestamp']))]
    # tstart = datetime.datetime(2022, 1, 12, 15, 45, 59)
    # xgsm = -2431.1245629621699
    # ygsm = 3822.9186030446831
    # zgsm = 5059.6970615621403
    # bxgsm, bygsm, bzgsm = get_igrf(tstart, xgsm, ygsm, zgsm)
    # print(bxgsm, bygsm, bzgsm)

    # coordinate transformation of B: gsm -> gse -> gei
    df['igrf_gse'] = pd.Series(cotrans_lib.subgsm2gse(df['timestamp'].tolist(), df['igrf_gsm'].tolist()).tolist())
    df['igrf_gei'] = pd.Series(cotrans_lib.subgse2gei(df['timestamp'].tolist(), df['igrf_gse'].tolist()).tolist())

    ############################################
    #   suyash code begin
    ############################################
    df['ctime'] = df['timestamp']-df['timestamp'][0]
    ctime = np.array(df['ctime'])
    delta_t = np.median(ctime[1:]-ctime[:-1])
    # get fgm data in fgm coordinate
    B_S1_corr, B_S2_corr, B_S3_corr  = list(zip(*df['fgm_fgm']))

    proper_pad = False  # fails when data have gaps
    fit_running_spline = False
    relative_integrate = True
    bidirectional_integrate = False
    func = cube_fit

    ############################################
    #   1.1.1 stage one - crossing times - choose
    ############################################
    d_B_S3_corr = np.gradient(B_S3_corr)/np.gradient(ctime)
    cross_times_d_B_S3_corr_1 = []
    for i in range(1, len(ctime)-2):
        if(d_B_S3_corr[i-1]>0 and d_B_S3_corr[i]>0 and d_B_S3_corr[i+1]<0 and d_B_S3_corr[i+2]<0 
        and B_S3_corr[i-1]>0 and B_S3_corr[i]>0 and B_S3_corr[i+1]>0 and B_S3_corr[i+2]>0):  
        #jwu: when gap exits, dB can jump from positive to negative 
            y1 = d_B_S3_corr[i]
            y2 = d_B_S3_corr[i+1]
            x1 = ctime[i]
            x2 = ctime[i+1]
            cross_times_d_B_S3_corr_1.append((y2*x1 - y1*x2)/(y2 - y1))      
    # List of crossing times
    cross_times_d_B_S3_corr_1 = np.array(cross_times_d_B_S3_corr_1)
    # List of middle points of crossing times, helpful for interpolation purposes earlier
    cross_times_d_B_S3_corr_1_mids = .5*(cross_times_d_B_S3_corr_1[1:]+cross_times_d_B_S3_corr_1[:-1])
    # Spin-periods computed as difference between consecutive zero-crossings
    T_spins_d_corr_1 = cross_times_d_B_S3_corr_1[1:]-cross_times_d_B_S3_corr_1[:-1]
    # Corresponding angular velocities
    w_syn_d_corr_1 = 2*np.pi/T_spins_d_corr_1
    ############################################
    #   1.1.2 stage one - crossing times - select
    ############################################
    eps_corr_1 = 1e5
    # Get indices of valid spin-periods
    valid_idx_corr_1 = running_filter(cross_times_d_B_S3_corr_1_mids, w_syn_d_corr_1, eps_corr_1)

    # Select the corresponding crossing time mid-points, synodic angular velocities and spin-periods
    cross_times_d_B_S3_corr_1_mids_select = cross_times_d_B_S3_corr_1_mids[valid_idx_corr_1]
    w_syn_d_corr_1_select = w_syn_d_corr_1[valid_idx_corr_1]
    T_spins_d_corr_1_select = T_spins_d_corr_1[valid_idx_corr_1]

    # Reconstruct the selected crossing times themselves
    cross_times_d_B_S3_corr_1_select = cross_times_d_B_S3_corr_1_mids_select-T_spins_d_corr_1_select/2
    cross_times_d_B_S3_corr_1_select = np.concatenate((cross_times_d_B_S3_corr_1_select, np.array([cross_times_d_B_S3_corr_1_mids_select[-1]+T_spins_d_corr_1_select[-1]/2])))

    ############################################
    #   1.1.3 stage one - crossing times - pad
    ############################################
    if (proper_pad == True):    
        if (fit_running_spline == True):
            T_spins_d_pad_corr_1_select = running_spline(cross_times_d_B_S3_corr_1_select, 
                                              cross_times_d_B_S3_corr_1_mids_select,
                                              T_spins_d_corr_1_select, T=30)
        else:        
            T_spins_d_pad_corr_1_select = func(cross_times_d_B_S3_corr_1_select, 
                                           *curve_fit(func, 
                                                      cross_times_d_B_S3_corr_1_mids_select, 
                                                      T_spins_d_corr_1_select)[0])
    else:
            T_spins_d_pad_corr_1_select = np.append(T_spins_d_corr_1, T_spins_d_corr_1[-1])

    ############################################
    #   1.2.1 stage two - crossing times - choose
    ############################################ 
    # For noting phases for determining zero-crossings
    phase_corr = np.zeros(len(ctime))
    phase_corr[:] = np.nan

    for i in range(len(cross_times_d_B_S3_corr_1_select)):
        # Zero-crossing from stage 1
        t0 = cross_times_d_B_S3_corr_1_select[i]
        # Isolate a period around the zero-crossing
        idx = ((ctime-t0)>=-T_spins_d_pad_corr_1_select[i]/2) & ((ctime-t0)<=T_spins_d_pad_corr_1_select[i]/2) 
        # Use the arcsine function to get phase angle around the zero-crossing
        phase_corr[idx] = np.arcsin(d_B_S3_corr[idx]/np.max(np.abs(d_B_S3_corr[idx])))
    # Record zero crossings as locations of the phase passing over from positive to negative
    cross_times_d_B_S3_corr_2 = []

    for i in range(1, len(ctime)-2):
        if(phase_corr[i-1]>0 and phase_corr[i]>0 and phase_corr[i+1]<0 and phase_corr[i+2]<0 and
            B_S3_corr[i-1]>0 and B_S3_corr[i]>0 and B_S3_corr[i+1]>0 and B_S3_corr[i+2]>0):
        #if(phase_corr[i-1]>0 and phase_corr[i]>0 and phase_corr[i+1]<0 and phase_corr[i+2]<0):
            y1 = phase_corr[i]
            y2 = phase_corr[i+1]
            x1 = ctime[i]
            x2 = ctime[i+1]
            cross_times_d_B_S3_corr_2.append((y2*x1 - y1*x2)/(y2 - y1))
        
    # Obtaining synodic angular velocity samples - similar to stage 1
    cross_times_d_B_S3_corr_2 = np.array(cross_times_d_B_S3_corr_2)
    cross_times_d_B_S3_corr_2_mids = .5*(cross_times_d_B_S3_corr_2[1:]+cross_times_d_B_S3_corr_2[:-1])
    T_spins_d_corr_2 = cross_times_d_B_S3_corr_2[1:]-cross_times_d_B_S3_corr_2[:-1]
    w_syn_d_corr_2 = 2*np.pi/T_spins_d_corr_2       

    ############################################
    #   1.2.2 stage two - crossing times - select
    ############################################ 
    eps_corr_2 = 1e+5
    valid_idx_corr_2 = running_filter(cross_times_d_B_S3_corr_2_mids, w_syn_d_corr_2, eps_corr_2)

    cross_times_d_B_S3_corr_2_mids_select = cross_times_d_B_S3_corr_2_mids[valid_idx_corr_2]
    w_syn_d_corr_2_select = w_syn_d_corr_2[valid_idx_corr_2]

    T_spins_d_corr_2_select = T_spins_d_corr_2[valid_idx_corr_2]

    cross_times_d_B_S3_corr_2_select = cross_times_d_B_S3_corr_2_mids_select-T_spins_d_corr_2_select/2
    cross_times_d_B_S3_corr_2_select = np.concatenate((cross_times_d_B_S3_corr_2_select, np.array([cross_times_d_B_S3_corr_2_mids_select[-1]+T_spins_d_corr_2_select[-1]/2])))

    ############################################
    #   1.2.3 stage two - crossing times - pad
    ############################################
    if (proper_pad == True):  
        if (fit_running_spline == True):
            T_spins_d_pad_corr_2_select = running_spline(cross_times_d_B_S3_corr_2_select, 
                                              cross_times_d_B_S3_corr_2_mids_select,
                                              T_spins_d_corr_2_select, T=30)
        else:
            T_spins_d_pad_corr_2_select = func(cross_times_d_B_S3_corr_2_select, 
                                           *curve_fit(func, 
                                                      cross_times_d_B_S3_corr_2_mids_select, 
                                                      T_spins_d_corr_2_select)[0])
    else:
        T_spins_d_pad_corr_2_select = np.append(T_spins_d_corr_2, T_spins_d_corr_2[-1])

    ############################################
    #   1.3.1 stage three - crossing times - choose
    ############################################
    N_spins_fit = 4
    # For detecting peaks by fitting B_S3 itself instead of 
    #fitting its derivative and computing zero-crossings
    peak_detect = False
    # Stage three - refine crossing times again
    cross_times_d_B_S3_corr_3 = []
    w_syn_d_corr_3 = []

    for i in range(len(cross_times_d_B_S3_corr_2_select)):
        # Get the crossing time, synodic period, and angular velocity from stage 2
        t0 = cross_times_d_B_S3_corr_2_select[i]
        T_syn = T_spins_d_pad_corr_2_select[i]
        w_syn = 2*np.pi/T_syn
    
        # The total time over which the fit will be computed
        T_avg = N_spins_fit*T_syn
    
        # Dealing with edge cases around the beginning and the end of time
        # Similar idea to moving_average, running_spline, and running-filter
        if(t0-ctime[0]<T_avg/2):
            low_lim = ctime[0]
            high_lim = ctime[0]+T_avg
        elif(ctime[-1]-t0<T_avg/2):
            low_lim = t0-T_avg
            high_lim = t0
        else:   
            low_lim = t0-T_avg/2
            high_lim = t0+T_avg/2

        # Get time instances within N_spins_fit
        idx = (ctime>=low_lim) & (ctime<=high_lim)
    
        # Initialize time for this segment at the zero-crossing
        ctime_slice = ctime[idx]-t0

        # Slice the signal itself
    
        # In case you are trying to find the maxima of B_S3 directly
        if(peak_detect == True):
            signal_slice = B_S3_corr[idx]
            spin_func = lambda x, A, w, p: cosine_fit(x, -1, A, w, p, 0)

        # In case you are trying to go the derivative/ zero-crossing route
        else:        
            signal_slice = d_B_S3_corr[idx]
            spin_func = lambda x, A, w, p: sine_fit(x, 1, A, w, p, 0)

        # Fit the curve you want to work with!
        # Good initial guesses p0 really help
        spin_opt, spin_covar = curve_fit(spin_func, 
                                     ctime_slice, 
                                     signal_slice, 
                                     p0=[np.max(np.abs(signal_slice-np.mean(signal_slice))), 
                                         w_syn, 0])

        # Using the zero-phase and the angular velocity, computing the deviation to the crossing time
        delta_t0 = -spin_opt[2]/spin_opt[1]

        # Contextualize the computing perturbation in terms of the relative time for the science zone
        cross_times_d_B_S3_corr_3.append(t0+delta_t0)
        # Also save the fitted value for the angular velocity
        w_syn_d_corr_3.append(spin_opt[1])

    cross_times_d_B_S3_corr_3 = np.array(cross_times_d_B_S3_corr_3)
    w_syn_d_corr_3 = np.array(w_syn_d_corr_3)

    ############################################
    #   1.3.2 stage three - crossing times - select
    ############################################
    eps_corr_3 = 2
    valid_idx_corr_3 = running_filter(cross_times_d_B_S3_corr_3, w_syn_d_corr_3, eps_corr_3, T=50)
    cross_times_d_B_S3_corr_3_select = cross_times_d_B_S3_corr_3[valid_idx_corr_3]
    w_syn_d_corr_3_select = w_syn_d_corr_3[valid_idx_corr_3]
    zero_crossing_method = 3

    # Remember that the stage 1 and 2 sample angular velocities at mid points of zero-crossings
    if (zero_crossing_method == 1):
        cross_times_d_B_S3_corr = cross_times_d_B_S3_corr_1_select
        cross_times_d_B_S3_corr_mids = cross_times_d_B_S3_corr_1_mids_select
        w_syn_d_corr = w_syn_d_corr_1_select
    
    elif (zero_crossing_method == 2):
        cross_times_d_B_S3_corr = cross_times_d_B_S3_corr_2_select
        cross_times_d_B_S3_corr_mids = cross_times_d_B_S3_corr_2_mids_select
        w_syn_d_corr = w_syn_d_corr_2_select    
    
    # But stage 3 samples angular velocities at the zero-crossings themselves
    else:
        cross_times_d_B_S3_corr = cross_times_d_B_S3_corr_3_select
        w_syn_d_corr = w_syn_d_corr_3_select


    # Appropriately generate angular velocity samples at all the FGM sample times
    if (fit_running_spline == True):   
        if(zero_crossing_method == 1 or zero_crossing_method == 2):
            w_syn_corr = running_spline(ctime, cross_times_d_B_S3_corr_mids, w_syn_d_corr, T=30)
        else:
            w_syn_corr = running_spline(ctime, cross_times_d_B_S3_corr, w_syn_d_corr, T=50)
    
    else:
        if(zero_crossing_method == 1 or zero_crossing_method == 2):
            w_syn_corr = func(ctime, *curve_fit(func, cross_times_d_B_S3_corr_mids, w_syn_d_corr)[0])
        else:
            w_syn_corr = func(ctime, *curve_fit(func, cross_times_d_B_S3_corr, w_syn_d_corr)[0])    

    if (relative_integrate == True): 
        # Use multiple reference points for integration
        idx0_corrs = np.array([np.where(ctime<=t0_corr)[0][-1] for t0_corr in cross_times_d_B_S3_corr])
        if (fit_running_spline == True):
            if(zero_crossing_method == 1 or zero_crossing_method == 2):
                w_t0_corrs = running_spline(cross_times_d_B_S3_corr, cross_times_d_B_S3_corr_mids, w_syn_d_corr, T=30)
            else:
                w_t0_corrs = running_spline(cross_times_d_B_S3_corr, cross_times_d_B_S3_corr, w_syn_d_corr, T=50)
        else:
            if(zero_crossing_method == 1 or zero_crossing_method == 2):
                w_t0_corrs = func(cross_times_d_B_S3_corr, *curve_fit(func, cross_times_d_B_S3_corr_mids, w_syn_d_corr)[0])
            else:
                w_t0_corrs = func(cross_times_d_B_S3_corr, *curve_fit(func, cross_times_d_B_S3_corr, w_syn_d_corr)[0])      
    else:
        # Use just one reference point for integration  
        t0_corr = cross_times_d_B_S3_corr[0]
        idx0_corr = np.where(ctime<=t0_corr)[0][-1]
        if (fit_running_spline == True):
            if(zero_crossing_method == 1 or zero_crossing_method == 2):
                w_t0_corr = running_spline([t0_corr], cross_times_d_B_S3_corr_mids, w_syn_d_corr, T=30)[0]
            else:
                w_t0_corr = running_spline([t0_corr], cross_times_d_B_S3_corr, w_syn_d_corr, T=50)[0]
        else:
            if(zero_crossing_method == 1 or zero_crossing_method == 2):
                w_t0_corr = func(t0_corr, *curve_fit(func, cross_times_d_B_S3_corr_mids, w_syn_d_corr)[0])
            else:
                w_t0_corr = func(t0_corr, *curve_fit(func, cross_times_d_B_S3_corr, w_syn_d_corr)[0])


    phi_corr = np.zeros(len(ctime))
    for i in range(len(phi_corr)):
        if (relative_integrate == True):       
            if(bidirectional_integrate == True):            
                # In this case, find the closest zero-crossings in absolute sense
                t0_corr_idx = np.argmin(np.abs(ctime[i]-cross_times_d_B_S3_corr))

            else:           
                # In this case, find the closest zero crossing that comes before a given time
                # Except for times before the zero crossing, then use the first zero crossing           
                if(np.sum((ctime[i]-cross_times_d_B_S3_corr)>0) == 0):
                    t0_corr_idx = 0
                else:        
                    t0_corr_idx = np.nanargmin(mask_neg(ctime[i]-cross_times_d_B_S3_corr))
            t0_corr = cross_times_d_B_S3_corr[t0_corr_idx]
            idx0_corr = idx0_corrs[t0_corr_idx]
            w_t0_corr = w_t0_corrs[t0_corr_idx]

        # Do the integral based on the relative position of the time and the relevant zero-crossing  
        if (i<idx0_corr):
            phi_corr[i] = -simpson(w_syn_corr[i:idx0_corr+1], x=ctime[i:idx0_corr+1], dx=delta_t)
        elif (i>idx0_corr):
            phi_corr[i] = simpson(w_syn_corr[idx0_corr:i+1], x=ctime[idx0_corr:i+1], dx=delta_t)
        else:
            phi_corr[i] = 0
        
        # Correct for the zero-crossing not being in the time array  
        phi_corr[i] -= .5*(w_t0_corr+w_syn_corr[idx0_corr])*(t0_corr-ctime[idx0_corr])  

    B_IGRF_GEI_x, B_IGRF_GEI_y, B_IGRF_GEI_z  = np.array(list(zip(*df['igrf_gei'])))    
    att_GEI_x, att_GEI_y, att_GEI_z = np.array(list(zip(*df['att_gei'])))
    B_hat_IGRF_GEI_x = B_IGRF_GEI_x/np.sqrt(B_IGRF_GEI_x**2 + B_IGRF_GEI_y**2 + B_IGRF_GEI_z**2)
    B_hat_IGRF_GEI_y = B_IGRF_GEI_y/np.sqrt(B_IGRF_GEI_x**2 + B_IGRF_GEI_y**2 + B_IGRF_GEI_z**2)
    B_hat_IGRF_GEI_z = B_IGRF_GEI_z/np.sqrt(B_IGRF_GEI_x**2 + B_IGRF_GEI_y**2 + B_IGRF_GEI_z**2)

    DMXL_2_GEI = np.zeros((len(ctime), 3, 3))
    for i in range(len(ctime)):
        u_hat = np.array([att_GEI_x[i], att_GEI_y[i], att_GEI_z[i]])
        b_hat = np.array([B_hat_IGRF_GEI_x[i], B_hat_IGRF_GEI_y[i], B_hat_IGRF_GEI_z[i]])
    
        DMXL_2_GEI[i,:,0] = np.cross(b_hat, u_hat)
        DMXL_2_GEI[i,:,1] = np.cross(u_hat, np.cross(b_hat, u_hat))
        DMXL_2_GEI[i,:,2] = u_hat
    
        DMXL_2_GEI[i,:,0] /= np.linalg.norm(DMXL_2_GEI[i,:,0])
        DMXL_2_GEI[i,:,1] /= np.linalg.norm(DMXL_2_GEI[i,:,1])
        DMXL_2_GEI[i,:,2] /= np.linalg.norm(DMXL_2_GEI[i,:,2]) 

    B_x_DMXL_true = np.zeros(len(ctime))
    B_y_DMXL_true = np.zeros(len(ctime))
    B_z_DMXL_true = np.zeros(len(ctime))
    for i in range(len(ctime)):
        GEI_2_DMXL = np.linalg.inv(DMXL_2_GEI[i])  
        B_x_DMXL_true[i] = GEI_2_DMXL[0,0]*B_IGRF_GEI_x[i] + GEI_2_DMXL[0,1]*B_IGRF_GEI_y[i] + GEI_2_DMXL[0,2]*B_IGRF_GEI_z[i]
        B_y_DMXL_true[i] = GEI_2_DMXL[1,0]*B_IGRF_GEI_x[i] + GEI_2_DMXL[1,1]*B_IGRF_GEI_y[i] + GEI_2_DMXL[1,2]*B_IGRF_GEI_z[i]
        B_z_DMXL_true[i] = GEI_2_DMXL[2,0]*B_IGRF_GEI_x[i] + GEI_2_DMXL[2,1]*B_IGRF_GEI_y[i] + GEI_2_DMXL[2,2]*B_IGRF_GEI_z[i]
    B_x_SMXL_corr_respin = np.cos(-phi_corr)*B_x_DMXL_true - np.sin(-phi_corr)*B_y_DMXL_true
    B_y_SMXL_corr_respin = np.sin(-phi_corr)*B_x_DMXL_true + np.cos(-phi_corr)*B_y_DMXL_true
    B_z_SMXL_corr_respin = B_z_DMXL_true

    f = 44*np.pi/180

    B_S1_IGRF = np.cos(f)*B_x_SMXL_corr_respin + np.sin(f)*B_z_SMXL_corr_respin
    B_S2_IGRF = -np.sin(f)*B_x_SMXL_corr_respin + np.cos(f)*B_z_SMXL_corr_respin
    B_S3_IGRF = -B_y_SMXL_corr_respin

    n = len(B_S1_corr)
    b = np.concatenate((B_S1_corr, B_S2_corr, B_S3_corr))
    A = np.zeros((3*n, 12))
    A[0:n,0] = B_S1_IGRF
    A[0:n,1] = B_S2_IGRF
    A[0:n,2] = B_S3_IGRF
    A[0:n,3] = np.ones(n)
    A[n:2*n,4] = B_S1_IGRF
    A[n:2*n,5] = B_S2_IGRF
    A[n:2*n,6] = B_S3_IGRF
    A[n:2*n,7] = np.ones(n)
    A[2*n:3*n,8] = B_S1_IGRF
    A[2*n:3*n,9] = B_S2_IGRF
    A[2*n:3*n,10] = B_S3_IGRF
    A[2*n:3*n,11] = np.ones(n)
    A = csc_matrix(A)
    x = lsqr(A, b, atol=1e-10, btol=1e-10)[0]    

    orth = np.array([[x[0],x[1],x[2]], [x[4], x[5], x[6]], [x[8], x[9], x[10]]])
    offsets = np.array([x[3], x[7], x[11]])
    calib = np.linalg.inv(orth)

    B_S1_calib = calib[0, 0]*(B_S1_corr-offsets[0]) + calib[0, 1]*(B_S2_corr-offsets[1]) + calib[0, 2]*(B_S3_corr-offsets[2])
    B_S2_calib = calib[1, 0]*(B_S1_corr-offsets[0]) + calib[1, 1]*(B_S2_corr-offsets[1]) + calib[1, 2]*(B_S3_corr-offsets[2])
    B_S3_calib = calib[2, 0]*(B_S1_corr-offsets[0]) + calib[2, 1]*(B_S2_corr-offsets[1]) + calib[2, 2]*(B_S3_corr-offsets[2])


    ############################################
    #   2.1.1 stage one - crossing times - choose
    ############################################
    d_B_S3_calib = np.gradient(B_S3_calib)/np.gradient(ctime)

    cross_times_d_B_S3_calib_1 = []
    for i in range(len(ctime)-1):
        #if(d_B_S3_calib[i]*d_B_S3_calib[i+1]<0 and d_B_S3_calib[i]>0):
        if(d_B_S3_corr[i-1]>0 and d_B_S3_corr[i]>0 and d_B_S3_corr[i+1]<0 and d_B_S3_corr[i+2]<0 
        and B_S3_corr[i-1]>0 and B_S3_corr[i]>0 and B_S3_corr[i+1]>0 and B_S3_corr[i+2]>0): 
            y1 = d_B_S3_calib[i]
            y2 = d_B_S3_calib[i+1]
            x1 = ctime[i]
            x2 = ctime[i+1]
            cross_times_d_B_S3_calib_1.append((y2*x1 - y1*x2)/(y2 - y1))
        
    cross_times_d_B_S3_calib_1 = np.array(cross_times_d_B_S3_calib_1)
    cross_times_d_B_S3_calib_1_mids = .5*(cross_times_d_B_S3_calib_1[1:]+cross_times_d_B_S3_calib_1[:-1])
    T_spins_d_calib_1 = cross_times_d_B_S3_calib_1[1:]-cross_times_d_B_S3_calib_1[:-1]
    w_syn_d_calib_1 = 2*np.pi/T_spins_d_calib_1

    ############################################
    #   2.1.2 stage one - crossing times - select
    ############################################
    eps_calib_1 = eps_corr_1
    valid_idx_calib_1 = running_filter(cross_times_d_B_S3_calib_1_mids, w_syn_d_calib_1, eps_calib_1)
    cross_times_d_B_S3_calib_1_mids_select = cross_times_d_B_S3_calib_1_mids[valid_idx_calib_1]
    w_syn_d_calib_1_select = w_syn_d_calib_1[valid_idx_calib_1]

    T_spins_d_calib_1_select = T_spins_d_calib_1[valid_idx_calib_1]

    cross_times_d_B_S3_calib_1_select = cross_times_d_B_S3_calib_1_mids_select-T_spins_d_calib_1_select/2
    cross_times_d_B_S3_calib_1_select = np.concatenate((cross_times_d_B_S3_calib_1_select, np.array([cross_times_d_B_S3_calib_1_mids_select[-1]+T_spins_d_calib_1_select[-1]/2])))

    ############################################
    #   2.1.3 stage one - crossing times - pad
    ############################################
    if (proper_pad == True):   
        if (fit_running_spline == True):
            T_spins_d_pad_calib_1_select = running_spline(cross_times_d_B_S3_calib_1_select, 
                                              cross_times_d_B_S3_calib_1_mids_select,
                                              T_spins_d_calib_1_select, T=30)      
        else:        
            T_spins_d_pad_calib_1_select = func(cross_times_d_B_S3_calib_1_select, 
                                           *curve_fit(func, 
                                                      cross_times_d_B_S3_calib_1_mids_select, 
                                                      T_spins_d_calib_1_select)[0])
    else:
        T_spins_d_pad_calib_1_select = np.append(T_spins_d_calib_1, T_spins_d_calib_1[-1])

    ############################################
    #   2.2.1 stage two - crossing times - pad
    ############################################
    phase_calib = np.zeros(len(ctime))

    for i in range(len(cross_times_d_B_S3_calib_1_select)): 
        t0 = cross_times_d_B_S3_calib_1_select[i]  
        idx = ((ctime-t0)>=-T_spins_d_pad_calib_1_select[i]/2) & ((ctime-t0)<=T_spins_d_pad_calib_1_select[i]/2)
        phase_calib[idx] = np.arcsin(d_B_S3_calib[idx]/np.max(np.abs(d_B_S3_calib[idx])))
    
    cross_times_d_B_S3_calib_2 = []

    for i in range(len(ctime)-1):
        #if(phase_calib[i]*phase_calib[i+1]<0 and phase_calib[i]>0):
        if(phase_corr[i-1]>0 and phase_corr[i]>0 and phase_corr[i+1]<0 and phase_corr[i+2]<0 and
            B_S3_corr[i-1]>0 and B_S3_corr[i]>0 and B_S3_corr[i+1]>0 and B_S3_corr[i+2]>0):
            y1 = phase_calib[i]
            y2 = phase_calib[i+1]
            x1 = ctime[i]
            x2 = ctime[i+1]
            cross_times_d_B_S3_calib_2.append((y2*x1 - y1*x2)/(y2 - y1))
        
    cross_times_d_B_S3_calib_2 = np.array(cross_times_d_B_S3_calib_2)
    cross_times_d_B_S3_calib_2_mids = .5*(cross_times_d_B_S3_calib_2[1:]+cross_times_d_B_S3_calib_2[:-1])
    T_spins_d_calib_2 = cross_times_d_B_S3_calib_2[1:]-cross_times_d_B_S3_calib_2[:-1]
    w_syn_d_calib_2 = 2*np.pi/T_spins_d_calib_2

    ############################################
    #   2.2.2 stage two - crossing times - select
    ############################################
    eps_calib_2 = eps_corr_2

    valid_idx_calib_2 = running_filter(cross_times_d_B_S3_calib_2_mids, w_syn_d_calib_2, eps_calib_2)

    cross_times_d_B_S3_calib_2_mids_select = cross_times_d_B_S3_calib_2_mids[valid_idx_calib_2]
    w_syn_d_calib_2_select = w_syn_d_calib_2[valid_idx_calib_2]

    T_spins_d_calib_2_select = T_spins_d_calib_2[valid_idx_calib_2]

    cross_times_d_B_S3_calib_2_select = cross_times_d_B_S3_calib_2_mids_select-T_spins_d_calib_2_select/2
    cross_times_d_B_S3_calib_2_select = np.concatenate((cross_times_d_B_S3_calib_2_select, np.array([cross_times_d_B_S3_calib_2_mids_select[-1]+T_spins_d_calib_2_select[-1]/2])))

    ############################################
    #   2.2.3 stage two - crossing times - pad
    ############################################
    if (proper_pad == True):   
        if (fit_running_spline == True):
            T_spins_d_pad_calib_2_select = running_spline(cross_times_d_B_S3_calib_2_select, 
                                              cross_times_d_B_S3_calib_2_mids_select,
                                              T_spins_d_calib_2_select, T=30)       
        else:
        
            T_spins_d_pad_calib_2_select = func(cross_times_d_B_S3_calib_2_select, 
                                           *curve_fit(func, 
                                                      cross_times_d_B_S3_calib_2_mids_select, 
                                                      T_spins_d_calib_2_select)[0])
    else:
        T_spins_d_pad_calib_2_select = np.append(T_spins_d_calib_2, T_spins_d_calib_2[-1])

    ############################################
    #   2.3.1 stage three - crossing times - choose
    ############################################
    R2_filter = True
    R2_limit = 0.8

    # Stage three - refine crossing times again
    cross_times_d_B_S3_calib_3 = []
    w_syn_d_calib_3 = []
    R2score = []

    for i in range(len(cross_times_d_B_S3_calib_2_select)):
        # Get the crossing time
        t0 = cross_times_d_B_S3_calib_2_select[i]
        T_syn = T_spins_d_pad_calib_2_select[i]
        w_syn = 2*np.pi/T_syn
        T_avg = N_spins_fit*T_syn
        
        # fitting will not work if a big gap exist   
        if(t0-ctime[0]<T_avg/2):
            low_lim = ctime[0]
            high_lim = ctime[0]+T_avg
        elif(ctime[-1]-t0<T_avg/2):
            low_lim = t0-T_avg
            high_lim = t0
        else:   
            low_lim = t0-T_avg/2
            high_lim = t0+T_avg/2

        # Get time points within N_spins_fit
        idx = (ctime>=low_lim) & (ctime<=high_lim)
    
        # Initialize in this segment
        ctime_slice = ctime[idx]-t0

        # Slice the signal
        if(peak_detect == True):
            signal_slice = B_S3_calib[idx]
            spin_func = lambda x, A, w, p: cosine_fit(x, -1, A, w, p, 0)
        else:        
            signal_slice = d_B_S3_calib[idx]
            spin_func = lambda x, A, w, p: sine_fit(x, 1, A, w, p, 0)
            #spin_func = sine_fit_trial
    
        spin_opt, spin_covar = curve_fit(spin_func, 
                                     ctime_slice, 
                                     signal_slice, 
                                     p0=[np.max(np.abs(signal_slice-np.mean(signal_slice))), 
                                         w_syn, 0])
    
        # use R2 to exclude bad results (usually because gaps in the data)
        y_model = spin_func(ctime_slice, *spin_opt)
        R2score.append(1-sum((signal_slice-y_model)**2)/sum((signal_slice-np.mean(signal_slice))**2))
    
        delta_t0 = -spin_opt[2]/spin_opt[1] 

        cross_times_d_B_S3_calib_3.append(t0+delta_t0)
        w_syn_d_calib_3.append(spin_opt[1])

    if (R2_filter == True):
        # if R2 too low, replace with the results from stage two 
        for i in range(len(cross_times_d_B_S3_calib_3)):
            if (R2score[i] < R2_limit):
                cross_times_d_B_S3_calib_3[i] = cross_times_d_B_S3_calib_2_select[i]
                w_syn_d_calib_3[i] = T_spins_d_pad_calib_2_select[i]

    cross_times_d_B_S3_calib_3 = np.array(cross_times_d_B_S3_calib_3)
    w_syn_d_calib_3 = np.array(w_syn_d_calib_3)

    ############################################
    #   2.3.2 stage three - crossing times - select
    ############################################
    eps_calib_3 = eps_corr_3

    valid_idx_calib_3 = running_filter(cross_times_d_B_S3_calib_3, w_syn_d_calib_3, eps_calib_3, T=50)
    cross_times_d_B_S3_calib_3_select = cross_times_d_B_S3_calib_3[valid_idx_calib_3]
    w_syn_d_calib_3_select = w_syn_d_calib_3[valid_idx_calib_3]

    if (zero_crossing_method == 1):
        cross_times_d_B_S3_calib = cross_times_d_B_S3_calib_1_select
        cross_times_d_B_S3_calib_mids = cross_times_d_B_S3_calib_1_mids_select
        w_syn_d_calib = w_syn_d_calib_1_select   
    elif (zero_crossing_method == 2):
        cross_times_d_B_S3_calib = cross_times_d_B_S3_calib_2_select
        cross_times_d_B_S3_calib_mids = cross_times_d_B_S3_calib_2_mids_select
        w_syn_d_calib = w_syn_d_calib_2_select     
    else:
        cross_times_d_B_S3_calib = cross_times_d_B_S3_calib_3_select
        w_syn_d_calib = w_syn_d_calib_3_select

    if (fit_running_spline == True):   
        if(zero_crossing_method == 1 or zero_crossing_method == 2):
            w_syn_calib = running_spline(ctime, cross_times_d_B_S3_calib_mids, w_syn_d_calib, T=25)
        else:
            w_syn_calib = running_spline(ctime, cross_times_d_B_S3_calib, w_syn_d_calib, T=50)   
    else:
        if(zero_crossing_method == 1 or zero_crossing_method == 2):
            w_syn_calib = func(ctime, *curve_fit(func, cross_times_d_B_S3_calib_mids, w_syn_d_calib)[0])
        else:
            w_syn_calib = func(ctime, *curve_fit(func, cross_times_d_B_S3_calib, w_syn_d_calib)[0])

    ############################################
    #   integrate phi
    ############################################
    if (relative_integrate == True):
        idx0_calibs = np.array([np.where(ctime<=t0_calib)[0][-1] for t0_calib in cross_times_d_B_S3_calib])

        if (fit_running_spline == True):
            if(zero_crossing_method == 1 or zero_crossing_method == 2):
                w_t0_calibs = running_spline(cross_times_d_B_S3_calib, cross_times_d_B_S3_calib_mids, w_syn_d_calib, T=30)
            else:
                w_t0_calibs = running_spline(cross_times_d_B_S3_calib, cross_times_d_B_S3_calib, w_syn_d_calib, T=50)
        else:
            if(zero_crossing_method == 1 or zero_crossing_method == 2):
                w_t0_calibs = func(cross_times_d_B_S3_calib, *curve_fit(func, cross_times_d_B_S3_calib_mids, w_syn_d_calib)[0])
            else:
                w_t0_calibs = func(cross_times_d_B_S3_calib, *curve_fit(func, cross_times_d_B_S3_calib, w_syn_d_calib)[0])           
    else:   
        t0_calib = cross_times_d_B_S3_calib[0]
        idx0_calib = np.where(ctime<=t0_calib)[0][-1]
        if (fit_running_spline == True):
            if(zero_crossing_method == 1 or zero_crossing_method == 2):
                w_t0_calib = running_spline([t0_calib], cross_times_d_B_S3_calib_mids, w_syn_d_calib, T=30)[0]
            else:
                w_t0_calib = running_spline([t0_calib], cross_times_d_B_S3_calib, w_syn_d_calib, T=50)[0]
        else:
            if(zero_crossing_method == 1 or zero_crossing_method == 2):
                w_t0_calib = func(t0_calib, *curve_fit(func, cross_times_d_B_S3_calib_mids, w_syn_d_calib)[0])
            else:
                w_t0_calib = func(t0_calib, *curve_fit(func, cross_times_d_B_S3_calib, w_syn_d_calib)[0])

    phi_calib = np.zeros(len(ctime))
    for i in range(len(phi_calib)):  
        if (relative_integrate == True):     
            if(bidirectional_integrate == True):
                t0_calib_idx = np.argmin(np.abs(ctime[i]-cross_times_d_B_S3_calib))
            else:
                if(np.sum((ctime[i]-cross_times_d_B_S3_calib)>0) == 0):
                    t0_calib_idx = 0
                else:        
                    t0_calib_idx = np.nanargmin(mask_neg(ctime[i]-cross_times_d_B_S3_calib))
        
            t0_calib = cross_times_d_B_S3_calib[t0_calib_idx]
            idx0_calib = idx0_calibs[t0_calib_idx]
            w_t0_calib = w_t0_calibs[t0_calib_idx]

        if (i<idx0_calib):
            phi_calib[i] = -simpson(w_syn_calib[i:idx0_calib+1], x=ctime[i:idx0_calib+1], dx=delta_t)
        elif (i>idx0_calib):
            phi_calib[i] = simpson(w_syn_calib[idx0_calib:i+1], x=ctime[idx0_calib:i+1], dx=delta_t)
        else:
            phi_calib[i] = 0
        
        phi_calib[i] -= .5*(w_t0_calib+w_syn_calib[idx0_calib])*(t0_calib-ctime[idx0_calib])


    B_x_SMXL_calib = np.cos(f)*B_S1_calib - np.sin(f)*B_S2_calib
    B_y_SMXL_calib = -B_S3_calib
    B_z_SMXL_calib = np.sin(f)*B_S1_calib + np.cos(f)*B_S2_calib

    B_x_DMXL_calib = np.cos(phi_calib)*B_x_SMXL_calib - np.sin(phi_calib)*B_y_SMXL_calib
    B_y_DMXL_calib = np.sin(phi_calib)*B_x_SMXL_calib + np.cos(phi_calib)*B_y_SMXL_calib
    B_z_DMXL_calib = B_z_SMXL_calib

    if (zero_crossing_method == 1): 
        T_spins_d_calib = T_spins_d_pad_calib_1_select    
    elif (zero_crossing_method == 2):   
        T_spins_d_calib = T_spins_d_pad_calib_2_select 
    else:  
        T_spins_d_calib = 2*np.pi/w_syn_d_calib_3_select

    B_res_smooth_x = np.zeros(len(cross_times_d_B_S3_calib))
    B_res_smooth_y = np.zeros(len(cross_times_d_B_S3_calib))
    B_res_smooth_z = np.zeros(len(cross_times_d_B_S3_calib))

    for i in range(0, len(cross_times_d_B_S3_calib)):
        t0 = cross_times_d_B_S3_calib[i]
        T_syn = T_spins_d_calib[i]
        w_syn = 2*np.pi/T_syn
        idx = ((ctime-t0)>=-.5*T_syn) & ((ctime-t0)<=.5*T_syn)   
        ctime_slice = ctime[idx] 
        B_x_DMXL_resid_slice = (B_x_DMXL_calib-B_x_DMXL_true)[idx]
        B_y_DMXL_resid_slice = (B_y_DMXL_calib-B_y_DMXL_true)[idx]
        B_z_DMXL_resid_slice = (B_z_DMXL_calib-B_z_DMXL_true)[idx]
    
        FAC_func = lambda x, A, k: sine_fit(x, 1, A, w_syn, -w_syn*t0, k)
    
        B_res_smooth_x[i] = curve_fit(FAC_func, ctime_slice, B_x_DMXL_resid_slice)[0][1]
        B_res_smooth_y[i] = curve_fit(FAC_func, ctime_slice, B_y_DMXL_resid_slice)[0][1]
        B_res_smooth_z[i] = curve_fit(FAC_func, ctime_slice, B_z_DMXL_resid_slice)[0][1]

    IGRF_x_gei = np.zeros(len(ctime))
    IGRF_y_gei = np.zeros(len(ctime))
    IGRF_z_gei = np.zeros(len(ctime))

    for i in range(len(ctime)):
        IGRF_x_gei[i] = DMXL_2_GEI[i][0,0]*B_x_DMXL_true[i] + DMXL_2_GEI[i][0,1]*B_y_DMXL_true[i] + DMXL_2_GEI[i][0,2]*B_z_DMXL_true[i]
        IGRF_y_gei[i] = DMXL_2_GEI[i][1,0]*B_x_DMXL_true[i] + DMXL_2_GEI[i][1,1]*B_y_DMXL_true[i] + DMXL_2_GEI[i][1,2]*B_z_DMXL_true[i]
        IGRF_z_gei[i] = DMXL_2_GEI[i][2,0]*B_x_DMXL_true[i] + DMXL_2_GEI[i][2,1]*B_y_DMXL_true[i] + DMXL_2_GEI[i][2,2]*B_z_DMXL_true[i]    

    B_x_gei = np.zeros(len(ctime))
    B_y_gei = np.zeros(len(ctime))
    B_z_gei = np.zeros(len(ctime))

    for i in range(len(ctime)):
        B_x_gei[i] = DMXL_2_GEI[i][0,0]*B_x_DMXL_calib[i] + DMXL_2_GEI[i][0,1]*B_y_DMXL_calib[i] + DMXL_2_GEI[i][0,2]*B_z_DMXL_calib[i]
        B_y_gei[i] = DMXL_2_GEI[i][1,0]*B_x_DMXL_calib[i] + DMXL_2_GEI[i][1,1]*B_y_DMXL_calib[i] + DMXL_2_GEI[i][1,2]*B_z_DMXL_calib[i]
        B_z_gei[i] = DMXL_2_GEI[i][2,0]*B_x_DMXL_calib[i] + DMXL_2_GEI[i][2,1]*B_y_DMXL_calib[i] + DMXL_2_GEI[i][2,2]*B_z_DMXL_calib[i]    

    IGRF_smooth_x = np.zeros(len(cross_times_d_B_S3_calib))
    IGRF_smooth_y = np.zeros(len(cross_times_d_B_S3_calib))
    IGRF_smooth_z = np.zeros(len(cross_times_d_B_S3_calib))

    IGRF_smooth_x_gei = np.zeros(len(cross_times_d_B_S3_calib))
    IGRF_smooth_y_gei = np.zeros(len(cross_times_d_B_S3_calib))
    IGRF_smooth_z_gei = np.zeros(len(cross_times_d_B_S3_calib))
    for i in range(0, len(cross_times_d_B_S3_calib)):
    
        t0 = cross_times_d_B_S3_calib[i]
        T_syn = T_spins_d_calib[i]
        w_syn = 2*np.pi/T_syn
    
        idx = ((ctime-t0)>=-.5*T_syn) & ((ctime-t0)<=.5*T_syn)
    
        ctime_slice = ctime[idx]
    
        IGRF_smooth_x[i] = np.average(B_x_DMXL_true[idx])
        IGRF_smooth_y[i] = np.average(B_y_DMXL_true[idx])
        IGRF_smooth_z[i] = np.average(B_z_DMXL_true[idx])
    
        # gei
        IGRF_smooth_x_gei[i] = np.average(IGRF_x_gei[idx])
        IGRF_smooth_y_gei[i] = np.average(IGRF_y_gei[idx])
        IGRF_smooth_z_gei[i] = np.average(IGRF_z_gei[idx]) 

    B_smooth_x = np.zeros(len(cross_times_d_B_S3_calib))
    B_smooth_y = np.zeros(len(cross_times_d_B_S3_calib))
    B_smooth_z = np.zeros(len(cross_times_d_B_S3_calib))

    B_smooth_x_gei = np.zeros(len(cross_times_d_B_S3_calib))
    B_smooth_y_gei = np.zeros(len(cross_times_d_B_S3_calib))
    B_smooth_z_gei = np.zeros(len(cross_times_d_B_S3_calib))

    for i in range(0, len(cross_times_d_B_S3_calib)):
    
        t0 = cross_times_d_B_S3_calib[i]
        T_syn = T_spins_d_calib[i]
        w_syn = 2*np.pi/T_syn
    
        idx = ((ctime-t0)>=-.5*T_syn) & ((ctime-t0)<=.5*T_syn)
    
        ctime_slice = ctime[idx]
    
        B_x_DMXL_slice = B_x_DMXL_calib[idx]
        B_y_DMXL_slice = B_y_DMXL_calib[idx]
        B_z_DMXL_slice = B_z_DMXL_calib[idx]
    
        FAC_func = lambda x, A, k: sine_fit(x, 1, A, w_syn, -w_syn*t0, k)
    
        B_smooth_x[i] = curve_fit(FAC_func, ctime_slice, B_x_DMXL_slice)[0][1]
        B_smooth_y[i] = curve_fit(FAC_func, ctime_slice, B_y_DMXL_slice)[0][1]
        B_smooth_z[i] = curve_fit(FAC_func, ctime_slice, B_z_DMXL_slice)[0][1]
    
        # gei
        B_x_gei_slice = B_x_gei[idx]
        B_y_gei_slice = B_y_gei[idx]
        B_z_gei_slice = B_z_gei[idx]
    
        FAC_func = lambda x, A, k: sine_fit(x, 1, A, w_syn, -w_syn*t0, k)
    
        B_smooth_x_gei[i] = curve_fit(FAC_func, ctime_slice, B_x_gei_slice)[0][1]
        B_smooth_y_gei[i] = curve_fit(FAC_func, ctime_slice, B_y_gei_slice)[0][1]
        B_smooth_z_gei[i] = curve_fit(FAC_func, ctime_slice, B_z_gei_slice)[0][1]
    

    return [B_smooth_x, B_smooth_y, B_smooth_z, 
        IGRF_smooth_x, IGRF_smooth_y, IGRF_smooth_z, 
        B_smooth_x_gei, B_smooth_y_gei, B_smooth_z_gei, 
        IGRF_smooth_x_gei, IGRF_smooth_y_gei, IGRF_smooth_z_gei]


if __name__ == "__main__":   
    
    starttime_str = "2022-01-12 15:45:51"
    endtime_str = "2022-01-12 15:52:04"

    # read state cdf for att
    sta_cdfpath = "test/ela_l1_state_defn_20220112_v01.cdf"   
    # read fgm cdf and clip
    fgm_cdfpath = "test/ela_l1_fgs_20220112_v01.cdf"

    [B_smooth_x, B_smooth_y, B_smooth_z, 
    IGRF_smooth_x, IGRF_smooth_y, IGRF_smooth_z, 
    B_smooth_x_gei, B_smooth_y_gei, B_smooth_z_gei, 
    IGRF_smooth_x_gei, IGRF_smooth_y_gei, IGRF_smooth_z_gei]=fgm_calib(starttime_str, endtime_str,sta_cdfpath, fgm_cdfpath)

    print(B_smooth_x)
    