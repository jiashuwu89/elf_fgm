from typing import List
import numpy as np
from scipy.optimize import curve_fit
from .. import parameter
from scipy.integrate import simpson
from . import calibration
from .Bplot import B_ctime_plot, B_ctime_plot_single, omega_fit, omega_stage123
from .error import CrossTime1Error, postproc_fgs_igrf
func = calibration.quad_fit

def cross_time_stage_1(
    ctime: List[float], B_S3: List[float], ctime_idx = None, ctime_idx_flag = None,
):
    """cross time determination stage 1
    1. use dB/dt = 0 find the peak of B_z
    2. use eps to select cross time
    3. pad
    """
    #--------------------------------------------
    #   1 - dB/dt
    #--------------------------------------------
    d_B_S3 = np.gradient(B_S3) / np.gradient(ctime)
    cross_times_1 = []
    for i in range(1, len(ctime) - 2):
        if (
            d_B_S3[i - 1] > 0
            and d_B_S3[i] > 0
            and d_B_S3[i + 1] < 0
            and d_B_S3[i + 2] < 0
            and B_S3[i - 1] > 0
            and B_S3[i] > 0
            and B_S3[i + 1] > 0
            and B_S3[i + 2] > 0
        ):
            # jwu: when gap exits, dB can jump from positive to negative
            
            y1 = d_B_S3[i]
            y2 = d_B_S3[i + 1]
            x1 = ctime[i]
            x2 = ctime[i + 1]
            if x2 - x1 < 5:  # if x1 and x2 seperate more than 1 spin period. there is a gap in data. this cross time should be disgard
                cross_times_1.append((y2 * x1 - y1 * x2) / (y2 - y1))
    # List of crossing times
    cross_times_1 = np.array(cross_times_1)
    # List of middle points of crossing times, helpful for interpolation purposes earlier
    cross_times_1_mids = 0.5 * (
        cross_times_1[1:] + cross_times_1[:-1]
    )
    # Spin-periods computed as difference between consecutive zero-crossings
    T_spins_d_1 = cross_times_1[1:] - cross_times_1[:-1]
    # Corresponding angular velocities
    w_syn_d_1 = 2 * np.pi / T_spins_d_1

    #--------------------------------------------
    #   2 - select
    #--------------------------------------------
    if len(cross_times_1) < 5:
        raise CrossTime1Error(1)

    # Get indices of valid spin-periods
    valid_idx_1 = calibration.running_filter(
        cross_times_1_mids, w_syn_d_1, parameter.eps_1,
    )

    # Select the corresponding crossing time mid-points, synodic angular velocities and spin-periods
    cross_times_1_mids_select = cross_times_1_mids[
        valid_idx_1
    ]
    w_syn_d_1_select = w_syn_d_1[valid_idx_1]
    T_spins_d_1_select = T_spins_d_1[valid_idx_1]

    # Reconstruct the selected crossing times themselves
    cross_times_1_select = (
        cross_times_1_mids_select - T_spins_d_1_select / 2
    )
    cross_times_1_select = np.concatenate(
        (
            cross_times_1_select,
            np.array(
                [
                    cross_times_1_mids_select[-1]
                    + T_spins_d_1_select[-1] / 2
                ]
            ),
        )
    )

    #--------------------------------------------
    #    3 - pad
    #--------------------------------------------
    if parameter.proper_pad == True:
        if parameter.fit_running_spline == True:
            T_spins_d_pad_1_select = calibration.running_spline(
                cross_times_1_select,
                cross_times_1_mids_select,
                T_spins_d_1_select,
                T=30,
            )
        else:
            T_spins_d_pad_1_select = func(
                cross_times_1_select,
                *curve_fit(
                    func, cross_times_1_mids_select, T_spins_d_1_select
                )[0],
            )
    else:
        T_spins_d_pad_1_select = np.append(T_spins_d_1, T_spins_d_1[-1])

    return [cross_times_1_select, cross_times_1_mids_select, T_spins_d_pad_1_select, w_syn_d_1_select]


def cross_time_stage_2(
    ctime: List[float], B_S3: List[float],
    cross_times_1_select: List[float], T_spins_d_pad_1_select: List[float],
):
    """cross time determination stage 2
    1. use phase to find the peak of B_z
    2. use eps to select cross time
    3. pad
    """
    #--------------------------------------------
    #   1 - phase angle
    #--------------------------------------------
    d_B_S3 = np.gradient(B_S3) / np.gradient(ctime)
    phase = np.zeros(len(ctime))
    phase[:] = np.nan

    for i in range(len(cross_times_1_select)):
        # Zero-crossing from stage 1
        t0 = cross_times_1_select[i]
        # Isolate a period around the zero-crossing
        idx = ((ctime - t0) >= -T_spins_d_pad_1_select[i] / 2) & (
            (ctime - t0) <= T_spins_d_pad_1_select[i] / 2
        )
        
        # Use the arcsine function to get phase angle around the zero-crossing
        phase[idx] = np.arcsin(d_B_S3[idx] / np.max(np.abs(d_B_S3[idx])))  
        
        # if gap exists, idx could be []; example starttime_str = "2022-01-01 14:57:42"
             
    # Record zero crossings as locations of the phase passing over from positive to negative
    cross_times_2 = []

    for i in range(1, len(ctime) - 2):
        if (
            phase[i - 1] > 0
            and phase[i] > 0
            and phase[i + 1] < 0
            and phase[i + 2] < 0
            and B_S3[i - 1] > 0
            and B_S3[i] > 0
            and B_S3[i + 1] > 0
            and B_S3[i + 2] > 0
        ):
            # if(phase_corr[i-1]>0 and phase_corr[i]>0 and phase_corr[i+1]<0 and phase_corr[i+2]<0):
            y1 = phase[i]
            y2 = phase[i + 1]
            x1 = ctime[i]
            x2 = ctime[i + 1]
            cross_times_2.append((y2 * x1 - y1 * x2) / (y2 - y1))

    # Obtaining synodic angular velocity samples - similar to stage 1
    cross_times_2 = np.array(cross_times_2)
    cross_times_2_mids = 0.5 * (
        cross_times_2[1:] + cross_times_2[:-1]
    )
    T_spins_d_2 = cross_times_2[1:] - cross_times_2[:-1]
    w_syn_d_2 = 2 * np.pi / T_spins_d_2

    #--------------------------------------------
    #   2 - select
    #--------------------------------------------
    if len(cross_times_2) < 3 :
        raise CrossTime1Error(2)

    valid_idx_2 = calibration.running_filter(
        cross_times_2_mids, w_syn_d_2, parameter.eps_2
    )

    cross_times_2_mids_select = cross_times_2_mids[
        valid_idx_2
    ]
    w_syn_d_2_select = w_syn_d_2[valid_idx_2]

    T_spins_d_2_select = T_spins_d_2[valid_idx_2]

    cross_times_2_select = (
        cross_times_2_mids_select - T_spins_d_2_select / 2
    )
    cross_times_2_select = np.concatenate(
        (
            cross_times_2_select,
            np.array(
                [
                    cross_times_2_mids_select[-1]
                    + T_spins_d_2_select[-1] / 2
                ]
            ),
        )
    )

    #--------------------------------------------
    #   3 - pad
    #--------------------------------------------
    if parameter.proper_pad == True:
        if parameter.fit_running_spline == True:
            T_spins_d_pad_2_select = calibration.running_spline(
                cross_times_2_select,
                cross_times_2_mids_select,
                T_spins_d_2_select,
                T=30,
            )
        else:
            T_spins_d_pad_2_select = func(
                cross_times_2_select,
                *curve_fit(
                    func, cross_times_2_mids_select, T_spins_d_2_select
                )[0],
            )
    else:
        T_spins_d_pad_2_select = np.append(T_spins_d_2, T_spins_d_2[-1])

    return [cross_times_2_select, cross_times_2_mids_select, T_spins_d_pad_2_select, w_syn_d_2_select]


def cross_time_stage_3(
    ctime: List[float], B_S3: List[float],
    cross_times_2_select: List[float], T_spins_d_pad_2_select: List[float], 
    ctime_idx = None, ctime_idx_flag = None, ctime_idx_timediff = None
):
    """cross time determination stage 3
    1. curvefit to find B peak
    2. select
    3. TODO: pad? 
    """
    d_B_S3 = np.gradient(B_S3) / np.gradient(ctime)
    cross_times_3 = []
    w_syn_d_3 = []
    R2score = []
    dt0 = []
    for i in range(len(cross_times_2_select)):
        # Get the crossing time, synodic period, and angular velocity from stage 2
        t0 = cross_times_2_select[i]
        T_syn = T_spins_d_pad_2_select[i]
        w_syn = 2 * np.pi / T_syn

        # The total time over which the fit will be computed
        T_tot = parameter.N_spins_fit * T_syn
        w_avg = 2 * np.pi / np.median(T_spins_d_pad_2_select)

        # Dealing with edge cases around the beginning and the end of time
        # Similar idea to moving_average, running_spline, and running-filter
        if t0 - ctime[0] < T_tot / 2:
            low_lim = ctime[0]
            high_lim = ctime[0] + T_tot
        elif ctime[-1] - t0 < T_tot / 2:
            low_lim = t0 - T_tot
            high_lim = t0
        else:
            low_lim = t0 - T_tot / 2
            high_lim = t0 + T_tot / 2

        # Get time instances within N_spins_fit
        idx = (ctime >= low_lim) & (ctime <= high_lim)

        # Initialize time for this segment at the zero-crossing
        ctime_slice = ctime[idx] - t0

        """
        if i > 0 and i < len(cross_times_2_select) -  1:
            cross_time1 = cross_times_2_select[i-1]
            cross_time2 = cross_times_2_select[i+1]
        elif i == 0:
            cross_time1 = ctime[0]
            cross_time2 = cross_times_2_select[i+1]
        elif i == len(cross_times_2_select) - 1:
            cross_time1 = cross_times_2_select[i-1]
            cross_time2 = ctime[-1]
        """
        # delete spins with 2.5s gap, they may result in spikes in fsp 
        if parameter.cross0_spike_del == True and ctime_idx is not None and ctime_idx_flag is not None and ctime_idx_timediff is not None:
            flag = 0
            for ctime_idx_idx, ctime_idx_val in enumerate(ctime_idx):
                if ctime_idx_flag[ctime_idx_idx] >= 2:
                    spike_time1 = ctime[ctime_idx_val]
                    spike_time2 = ctime[ctime_idx_val] + ctime_idx_timediff[ctime_idx_idx]
                    if (spike_time2 < low_lim or spike_time1 > high_lim) == 0 : 
                        if low_lim < spike_time2 < cross_times_2_select[i]:
                            low_lim = spike_time2 + np.pi/w_avg
                            idx = (ctime >= low_lim) & (ctime <= high_lim)
                            if len(ctime[idx]) <= 5 : # if too little date for fit after delete spike
                                flag = 1
                            else:
                                ctime_slice = ctime[idx] - t0
                        elif cross_times_2_select[i] < spike_time1 < high_lim: 
                            high_lim = spike_time1 - np.pi/w_avg
                            idx = (ctime >= low_lim) & (ctime <= high_lim)
                            if len(ctime[idx]) <= 5 : # if too little date for fit after delete spike
                                flag = 1
                            else:    
                                ctime_slice = ctime[idx] - t0
                        elif spike_time1 < cross_times_2_select[i] < spike_time2 : # if current zero crossing is inside spike
                            flag = 1
                            break

                    """
                    # use previous result if this fit include spike
                    if spike_time2 < low_lim or spike_time1 > high_lim : 
                        pass
                    else:
                        flag = 1 # if slice include spike, then skip fitting and use cross time stage 2 result
                    """
                    """
                    if spike_time2 < cross_time1 or spike_time1 > cross_time2 : 
                        pass
                    else:
                        flag = 3
                        breakpoint()
                        break
                    """

            if flag == 1:
                cross_times_3.append(t0)
                w_syn_d_3.append(w_syn)
                R2score.append(1)
                continue
            elif flag == 3:
                continue

        # Slice the signal itself
        # In case you are trying to find the maxima of B_S3 directly
        if parameter.peak_detect == True:
            signal_slice = B_S3[idx]
            spin_func = lambda x, A, w, p, k: calibration.cosine_fit(x, -1, A, w, p, k)

        # In case you are trying to go the derivative/ zero-crossing route
        else:
            signal_slice = d_B_S3[idx]
            spin_func = lambda x, A, w, p, k: calibration.sine_fit(x, 1, A, w, p, k)

        # Fit the curve you want to work with!
        # Good initial guesses p0 really help
        # JWu:we need B_max of dB = 0 with negative slope, this find -sin(phi) where phi = 0
        try:
            spin_opt, spin_covar = curve_fit(
                spin_func,
                ctime_slice,
                signal_slice,
                p0=[-np.max(np.abs(signal_slice - np.mean(signal_slice))), w_avg, 0, 0],
                bounds=((-np.inf, -np.inf, -np.inf, -100), (0, np.inf, np.inf, 100)),
            )
            # Using the zero-phase and the angular velocity, computing the deviation to the crossing time
            delta_t0 = -spin_opt[2] / spin_opt[1]
            signal_slice_fit = spin_func(ctime_slice, *spin_opt)
            R2score.append(1-sum((signal_slice - signal_slice_fit)**2)/sum((signal_slice - np.mean(signal_slice))**2))
            # Contextualize the computing perturbation in terms of the relative time for the science zone
            if t0 + delta_t0 < 0:
                delta_t0 = -(spin_opt[2] - 2*np.pi) / spin_opt[1]
            cross_times_3.append(t0 + delta_t0)
            dt0.append(delta_t0)
            # Also save the fitted value for the angular velocity
            w_syn_d_3.append(spin_opt[1])

        except:
            #breakpoint()
            print(f"cross time 3 determination fitting fail first try: {i}\n")
            try:
                spin_opt, spin_covar = curve_fit(
                    spin_func,
                    ctime_slice,
                    signal_slice,
                    p0=[np.max(np.abs(signal_slice - np.mean(signal_slice))), w_avg, 0, 0]
                )
                # Using the phase == pi , computing the deviation to the crossing time
                delta_t0 = (np.pi - spin_opt[2]) / spin_opt[1]
                signal_slice_fit = spin_func(ctime_slice, *spin_opt)
                R2score.append(1-sum((signal_slice - signal_slice_fit)**2)/sum((signal_slice - np.mean(signal_slice))**2))
                # Contextualize the computing perturbation in terms of the relative time for the science zone
                if t0 + delta_t0 < 0:
                    delta_t0 = -(spin_opt[2] - 2*np.pi) / spin_opt[1]
                cross_times_3.append(t0 + delta_t0)
                dt0.append(delta_t0)
                # Also save the fitted value for the angular velocity
                w_syn_d_3.append(spin_opt[1])

            except:
                # if fit still cannot find
                delta_t0 = 0
                R2score.append(0)
                cross_times_3.append(t0 + delta_t0)
                dt0.append(delta_t0)
                # Also save the fitted value for the angular velocity
                w_syn_d_3.append(2 * np.pi / T_spins_d_pad_2_select[i])

    if parameter.R2_filter == True:
        # if R2 too low, replace with the results from stage two 
        for i in range(len(cross_times_2_select)):
            if (R2score[i] < parameter.R2_thrhld):
                cross_times_3[i] = cross_times_2_select[i]
                w_syn_d_3[i] = 2 * np.pi / T_spins_d_pad_2_select[i]

    cross_times_3 = np.array(cross_times_3)
    w_syn_d_3 = np.array(w_syn_d_3)

    """
    2.5s spike: delete zero crossings within  2.5s + halft spin
    """
    """
    if parameter.cross0_spike_del == True and ctime_idx is not None and ctime_idx_flag is not None and ctime_idx_timediff is not None:
        # delete 2.5 purple spike
        ctime_idx_time = ctime[ctime_idx[ctime_idx_flag == 3]]
        ctime_idx_timediffs = ctime_idx_timediff[ctime_idx_flag == 3]
        for ctime_idx_time_idx, ctime_idx_time_val in enumerate(ctime_idx_time):
            try:
                idx1 = np.where(cross_times_3 > ctime_idx_time_val - 1.5*np.pi/w_avg)[0][0]
                idx2 = np.where(cross_times_3 < ctime_idx_time_val + ctime_idx_timediffs[ctime_idx_time_idx] + 1.5*np.pi/w_avg)[0][-1]
                idxs = range(idx1, idx2) if idx2 + 1 >= len(w_syn_d_3) else range(idx1, idx2+1)
                cross_times_3 = np.delete(cross_times_3, idxs)
                w_syn_d_3 = np.delete(w_syn_d_3, idxs)
            except:
                continue
    """
    #B_ctime_plot_single(cross_times_3, dt0)
    #breakpoint()
    #--------------------------------------------
    #   2 select
    #--------------------------------------------
    if len(cross_times_3) < 3:
        raise CrossTime1Error(3)

    valid_idx_3 = calibration.running_filter(
        cross_times_3, w_syn_d_3, parameter.eps_3
    )
    cross_times_3_select = cross_times_3[valid_idx_3]
    w_syn_d_3_select = w_syn_d_3[valid_idx_3] 
    T_spins_d_3_select = 2 * np.pi / w_syn_d_3_select

    return [cross_times_3_select, T_spins_d_3_select, w_syn_d_3_select]


def phase_integration(
    ctime, cross_times_1_select, cross_times_1_mids_select, 
    w_syn_d_1_select, T_spins_d_pad_1_select,
    cross_times_2_select, cross_times_2_mids_select,
    w_syn_d_2_select, T_spins_d_pad_2_select,
    cross_times_3_select, w_syn_d_3_select, T_spins_d_3_select, CrossTime_Update=False,
):
    """cross time determination stage 3
    curvefit to find B peak
    w_syn_d is the measured spin period
    w_syn: use cross_times, w_syn_d to fit ctime
    w_t0s: use cross_times, w_syn_d to fit cross_times
    """
    if CrossTime_Update == True:
        CrossTime_Update = parameter.CrossTime_Update

    if parameter.zero_crossing_method == 1:
        cross_times = cross_times_1_select
        cross_times_mids = cross_times_1_mids_select
        w_syn_d = w_syn_d_1_select
        T_spins_d = T_spins_d_pad_1_select

    elif parameter.zero_crossing_method == 2:
        cross_times = cross_times_2_select
        cross_times_mids = cross_times_2_mids_select
        w_syn_d = w_syn_d_2_select
        T_spins_d = T_spins_d_pad_2_select

    # But stage 3 samples angular velocities at the zero-crossings themselves
    else:
        cross_times = cross_times_3_select
        w_syn_d = w_syn_d_3_select
        T_spins_d = T_spins_d_3_select

    # in case first cross time is before 0
    if cross_times[0] < 0:
        cross_times = np.delete(cross_times, 0)
        w_syn_d = np.delete(w_syn_d, 0)
        T_spins_d = np.delete(T_spins_d, 0)

    delta_t = np.median(ctime[1:]-ctime[:-1])
    if parameter.fit_running_spline == True:
        if parameter.zero_crossing_method == 1 or parameter.zero_crossing_method == 2:
            w_syn = calibration.running_spline(
                ctime, cross_times_mids, w_syn_d, T=30
            )
        else:
            w_syn = calibration.running_spline(
                ctime, cross_times, w_syn_d, T=50
            )

    else:
        if parameter.zero_crossing_method == 1 or parameter.zero_crossing_method == 2:
            w_syn = func(
                ctime, *curve_fit(func, cross_times_mids, w_syn_d)[0]
            )
        else:
            #breakpoint()
            w_syn = func(
                ctime, *curve_fit(func, cross_times, w_syn_d)[0] #****
            )
    if parameter.relative_integrate == True:
        # Use multiple reference points for integration
        idx0s = np.array(
                [np.where(ctime <= t0)[0][-1] for t0 in cross_times]
        )
        if parameter.fit_running_spline == True:
            if parameter.zero_crossing_method == 1 or parameter.zero_crossing_method == 2:
                w_t0s = calibration.running_spline(
                    cross_times,
                    cross_times_mids,
                    w_syn_d,
                    T=30,
                )
            else:
                w_t0s = calibration.running_spline(
                    cross_times, cross_times, w_syn_d, T=50
                )
        else:
            if parameter.zero_crossing_method == 1 or parameter.zero_crossing_method == 2:
                w_t0s = func(
                    cross_times,
                    *curve_fit(func, cross_times_mids, w_syn_d)[0],
                )
            else:
                w_t0s = func(
                    cross_times,
                    *curve_fit(func, cross_times, w_syn_d)[0], #****
                )
        w_fit = w_t0s
        if parameter.zero_crossing_method == 1 or parameter.zero_crossing_method == 2:
            cross_times_fit = 2*np.pi/w_t0s-2*np.pi/w_syn_d + cross_times_mids
            if CrossTime_Update == True:
                cross_times_mids = cross_times_fit
        else:
            cross_times_fit = 2*np.pi/w_t0s-2*np.pi/w_syn_d + cross_times
            if CrossTime_Update == True:
                cross_times = cross_times_fit
    else:
        # Use just one reference point for integration
        t0 = cross_times[0]
        idx0 = np.where(ctime <= t0)[0][-1]
        if parameter.fit_running_spline == True:
            if parameter.zero_crossing_method == 1 or parameter.zero_crossing_method == 2:
                w_t0 = calibration.running_spline(
                    [t0], cross_times_mids, w_syn_d, T=30
                )[0]
                w_t0s = calibration.running_spline(
                    cross_times,
                    cross_times_mids,
                    w_syn_d,
                    T=30,
                )
            else:
                w_t0 = calibration.running_spline(
                    [t0], cross_times, w_syn_d, T=50
                )[0]
                w_t0s = calibration.running_spline(
                    cross_times, cross_times, w_syn_d, T=50
                )
        else:
            if parameter.zero_crossing_method == 1 or parameter.zero_crossing_method == 2:
                w_t0 = func(
                    t0, *curve_fit(func, cross_times_mids, w_syn_d)[0],
                )
                w_t0s = func(
                    cross_times,
                    *curve_fit(func, cross_times_mids, w_syn_d)[0],
                )
            else:
                w_t0 = func(
                    t0, *curve_fit(func, cross_times, w_syn_d)[0]
                )
                w_t0s = func(
                    cross_times,
                    *curve_fit(func, cross_times, w_syn_d)[0],
                )
        
        w_fit = w_t0s
        if parameter.zero_crossing_method == 1 or parameter.zero_crossing_method == 2:
            cross_times_fit = 2*np.pi/w_t0s-2*np.pi/w_syn_d + cross_times_mids
            if CrossTime_Update == True:
                cross_times_mids = cross_times_fit
        else:
            cross_times_fit = 2*np.pi/w_t0s-2*np.pi/w_syn_d + cross_times
            if CrossTime_Update == True:
                cross_times = cross_times_fit

    phi = np.zeros(len(ctime))
    for i in range(len(phi)):
        if parameter.relative_integrate == True:
            if parameter.bidirectional_integrate == True:
                # In this case, find the closest zero-crossings in absolute sense
                t0_corr_idx = np.argmin(np.abs(ctime[i] - cross_times))

            else:
                # In this case, find the closest zero crossing that comes before a given time
                # Except for times before the zero crossing, then use the first zero crossing
                if np.sum((ctime[i] - cross_times) > 0) == 0:
                    t0_corr_idx = 0
                else:
                    t0_corr_idx = np.nanargmin(
                        calibration.mask_neg(ctime[i] - cross_times)
                    )
            t0 = cross_times[t0_corr_idx]
            idx0 = idx0s[t0_corr_idx]
            w_t0 = w_t0s[t0_corr_idx]

        # Do the integral based on the relative position of the time and the relevant zero-crossing
        if i < idx0:
            phi[i] = -simpson(
                w_syn[i : idx0 + 1], x=ctime[i : idx0 + 1], dx=delta_t
            )
        elif i > idx0:
            phi[i] = simpson(
                w_syn[idx0 : i + 1], x=ctime[idx0 : i + 1], dx=delta_t
            )
        else:
            phi[i] = 0

        # Correct for the zero-crossing not being in the time array
        phi[i] -= (
            0.5 * (w_t0 + w_syn[idx0]) * (t0 - ctime[idx0])
        )            

    return [phi, cross_times, w_syn_d, T_spins_d, cross_times_fit, w_fit]     


def fsp_igrf(ctime, cross_times, T_spins_d, fgs_x, fgs_y, fgs_z): 
    """generate igrf in fsp resolution
        use average directly
    """
    fgs_fsp_x = np.zeros(len(cross_times))
    fgs_fsp_y = np.zeros(len(cross_times))
    fgs_fsp_z = np.zeros(len(cross_times))
    for i in range(0, len(cross_times)):

        t0 = cross_times[i]
        T_syn = T_spins_d[i]

        idx = ((ctime - t0) >= -0.5 * T_syn) & ((ctime - t0) <= 0.5 * T_syn)
        try:
            if len(ctime[idx]) > 3:
                fgs_fsp_x[i] = np.average(fgs_x[idx])
                fgs_fsp_y[i] = np.average(fgs_y[idx])
                fgs_fsp_z[i] = np.average(fgs_z[idx])
            else:
                fgs_fsp_x[i] = 0
                fgs_fsp_y[i] = 0
                fgs_fsp_z[i] = 0 
        except:
            raise postproc_fgs_igrf()

    return [fgs_fsp_x, fgs_fsp_y, fgs_fsp_z]


def fsp_ful(ctime, cross_times, T_spins_d, fgs_x, fgs_y, fgs_z): 
    """generate full field in fsp
        use curve fit
    """
    fgs_fsp_x = np.zeros(len(cross_times))
    fgs_fsp_y = np.zeros(len(cross_times))
    fgs_fsp_z = np.zeros(len(cross_times))

    for i in range(0, len(cross_times)):

        t0 = cross_times[i]
        T_syn = T_spins_d[i]
        w_syn = 2 * np.pi / T_syn
        idx = ((ctime - t0) >= -0.5 * T_syn) & ((ctime - t0) <= 0.5 * T_syn)
        
        if len(ctime[idx]) > 3:
            ctime_slice = ctime[idx]

            fgs_x_slice = fgs_x[idx]
            fgs_y_slice = fgs_y[idx]
            fgs_z_slice = fgs_z[idx]

            FAC_func = lambda x, A, k: calibration.sine_fit(x, 1, A, w_syn, -w_syn * t0, k)
            
            fgs_fsp_x[i] = curve_fit(FAC_func, ctime_slice, fgs_x_slice)[0][1]
            fgs_fsp_y[i] = curve_fit(FAC_func, ctime_slice, fgs_y_slice)[0][1]
            fgs_fsp_z[i] = curve_fit(FAC_func, ctime_slice, fgs_z_slice)[0][1]
        else:
            # when there is a big gap and no enough data for fitting
            fgs_fsp_x[i] = 0
            fgs_fsp_y[i] = 0
            fgs_fsp_z[i] = 0 

    return [fgs_fsp_x, fgs_fsp_y, fgs_fsp_z]


def fsp_matrix(ctime, cross_times, T_spins_d, rotation_matrix): 
    """generate rotation matrix in fsp resolution
    """
    rotation_matrix_fsp = np.zeros((len(cross_times), 3, 3))

    for i in range(0, len(cross_times)):

        t0 = cross_times[i]
        T_syn = T_spins_d[i]

        idx = ((ctime - t0) >= -0.5 * T_syn) & ((ctime - t0) <= 0.5 * T_syn)

        m = np.average(rotation_matrix[idx, :, :], axis=0)
        rotation_matrix_fsp[i, :, :] = m

    return rotation_matrix_fsp   
