from audioop import maxpp
import numpy as np
from scipy.interpolate import interp1d
from .. import parameter 
from . import Bplot, calibration, error
from bisect import bisect_left
from scipy.optimize import curve_fit
import traceback
from functools import reduce

def ctime_calib(ctime, B_x, B_y , B_z, cross_times, logger = None):
    """degap when ctime has jumps
       multipolar jumps: 
       unipolar jumps:  
    """
    delta_t = np.median(ctime[1:]-ctime[:-1])
    ctime_adj = ctime[1:]-ctime[:-1] - delta_t
    ctime_idx = np.where(np.abs(ctime_adj) > parameter.ctime_thrhld)[0]

    #if parameter.makeplot == True:
    #    Bplot.ctimediff_plot(ctime, ctime_idx = ctime_idx)
    #breakpoint()

    """
    flag = 1 multiploar spike, corrected
    flag = 2 unipolar spike 1/80 s
    flag = 3 unipolar gap
    flag = 4 other spike
    ['red','orange','magenta','darkviolet']
    """
    ctime_idx_flag = np.zeros(len(ctime_idx), dtype = int)
    spike_ctime_idxs = []  # save 1/80 spike
    i = 0
    while i < len(ctime_idx):
        if i < len(ctime_idx)-1 and (ctime_adj[ctime_idx[i]]*ctime_adj[ctime_idx[i+1]] < 0 # multipolar jumps
            and np.abs(np.abs(ctime_adj[ctime_idx[i]]) - np.abs(ctime_adj[ctime_idx[i+1]])) < 1e-3): # pos and neg jumps similar magnitude
            ctime_idx_flag[i] = 1
            ctime_idx_flag[i+1] = 1
            mid_idx = [i for i in range(ctime_idx[i]+1, ctime_idx[i+1])]
            ctimediff_mid = np.median(ctime[mid_idx][1:]-ctime[mid_idx][:-1])
            shift_1 = np.abs((ctime[1:]-ctime[:-1])[ctime_idx[i]] - ctimediff_mid)
            shift_2 = np.abs((ctime[1:]-ctime[:-1])[ctime_idx[i+1]] - ctimediff_mid)
            shift_idx_act = [i for i in range(ctime_idx[i]+1, ctime_idx[i+1]+1)]
            shift = np.mean([shift_1, shift_2])
            ctime[shift_idx_act] -= shift

            i += 2
        else: # unipolar jumps
            delta_dt = ctime_adj[ctime_idx[i]]
            """
                this code only works for unipolar spike of 1/80 s, larger than that is gaps
            """
            if np.abs(np.abs(delta_dt) - 0.0125) < 0.01 :
                """ unipolar spike with 1/80 s
                    if no calibration, then mark the orginal spike locaiton
                    if calibrate, then find the actual spike before the original one and mark 
                """
                if parameter.ctime_correct != True:  # if no calibration, then mark the orginal spike locaiton
                    spike_ctime_idxs.append(ctime_idx[i]) 
                    ctime_idx_flag[i] = 2 
                else:  
                    # first everything on the right of the spike
                    ctime[ctime_idx[i]+1:] = ctime[ctime_idx[i]+1:] - delta_dt
                    try:
                        # get the index for 3 spins with spike, and idx for other spins without spike
                        avg_ctime_idx, spike_ctime_idx = sepoch_getspin(ctime_idx[i], ctime, cross_times, ctime_idx) # index for three spins
                    
                        #if parameter.makeplot == True:
                        #    Bplot.B_ctime_plot(ctime, B_x, B_y, B_z, 
                        #    ctime_idx_time = ctime[ctime_idx],xlimt = [ctime[avg_ctime_idx[0][0]], ctime[avg_ctime_idx[-1][-1]]],
                        #    scatter = True, title = "before_spike_correction")

                        # get average Bx, By, Bz
                        avg_Bx, avg_By, avg_Bz, spike_Bx, spike_By, spike_Bz = sepoch_avg(spike_ctime_idx, avg_ctime_idx, B_x, B_y, B_z, ctime)
                        #Bplot.B_3d(spike_Bx, spike_By, spike_Bz)
                        #Bplot.B_ctime_plot_single(ctime[spike_ctime_idx], spike_Bx**2+spike_By**2+spike_Bz**2, scatter = True)
                        
                        # subtract Bx, By, Bz, get spike index
                        spike_idx1, spike_idx2 = sepoch_sub(ctime, ctime_idx[i], spike_ctime_idx, avg_Bx, avg_By, avg_Bz, spike_Bx, spike_By, spike_Bz)
                        
                        spike_ctime_idx1 = spike_ctime_idx[spike_idx1]
                        spike_ctime_idx2 = spike_ctime_idx[spike_idx2]
                        
                    except (error.spikeError_t1t2, error.spikeError_spikespin, error.spikeError_spikcrosstime) as e:
                        
                        logger.error(e.__str__())
                        spike_ctime_idx1 = ctime_idx[i]

                    # move delta_t to t1
                    ctime[spike_ctime_idx1:] = ctime[spike_ctime_idx1:] + delta_dt
                    spike_ctime_idxs.append(spike_ctime_idx1)
                    ctime_idx_flag[i] = 2 

            elif np.abs(delta_dt) - 2 > 0 : # gaps
                ctime_idx_flag[i] = 3
            else:
                ctime_idx_flag[i] = 4
            i += 1


    return ctime, ctime_idx, ctime_idx_flag, spike_ctime_idxs    


def sepoch_getspin(idx, ctime, cross_times, ctime_idx):
    """
    get ctime index for spins with spike
    TODO: end index is the one spin after spike
    """

    # get index for spike
    spike_spin_idx2 = np.where(cross_times>ctime[idx])[0] + 1 # index of spin cross zero after spike
    if spike_spin_idx2.size != 0 :
        spike_spin_idx2 = spike_spin_idx2[0]
        spike_spin_idx1 = spike_spin_idx2 - parameter.spike_find_spin_num  
        if spike_spin_idx1 < 0 :
            spike_spin_idx1 = 0 
            spike_find_spin_num = spike_spin_idx2 - spike_spin_idx1
        elif spike_spin_idx2 >= len(cross_times) :
            spike_spin_idx2 = len(cross_times) - 1
            spike_find_spin_num = spike_spin_idx2 - spike_spin_idx1
        else:
            spike_find_spin_num = parameter.spike_find_spin_num 

        spike_ctime_idx1, spike_ctime_time1 = find_closest(ctime, cross_times[spike_spin_idx1])
        spike_ctime_idx2, spike_ctime_time2 = find_closest(ctime, cross_times[spike_spin_idx2])
        if np.abs(spike_ctime_time1 - cross_times[spike_spin_idx1]) < 0.5 and np.abs(spike_ctime_time2 - cross_times[spike_spin_idx2]) < 0.5:
            spike_ctime_idx = [*range(spike_ctime_idx1, spike_ctime_idx2)]
        else:
            raise error.spikeError_spikespin(ctime[idx])
        
        if round((spike_ctime_time2 - spike_ctime_time1)/2.8) != parameter.spike_find_spin_num:
            spike_find_spin_num = round((spike_ctime_time2 - spike_ctime_time1)/2.8)

        # get index for average spins
        avg_spin_idx1 = spike_spin_idx1 - 10 if spike_spin_idx1 > 10 else 0 # left end of the avg spin
        avg_spin_idx2 = spike_spin_idx1 + 10 if spike_spin_idx2 < len(cross_times)-10 else len(cross_times) - spike_find_spin_num # right end of the avg spin

        avg_ctime_idx = []
        for spin_i in range(avg_spin_idx1, avg_spin_idx2):
            if spin_i + spike_find_spin_num <= spike_spin_idx1 or spin_i >= spike_spin_idx2 : # exclude spin with spike
                avg_ctime_idx1, avg_ctime_time1 = find_closest(ctime, cross_times[spin_i])
                avg_ctime_idx2, avg_ctime_time2 = find_closest(ctime, cross_times[spin_i + spike_find_spin_num])
                if np.abs(avg_ctime_time1 - cross_times[spin_i]) < 0.1 and np.abs(avg_ctime_time2 - cross_times[spin_i + spike_find_spin_num]) < 0.1:
                    avg_ctime_idx.append([*range(avg_ctime_idx1, avg_ctime_idx2)])

        # if spin includes other spikes, delete it
        for ctime_idx_i in ctime_idx:
            if ctime_idx_i != idx:
                for avg_idx_i, avg_idx in enumerate(avg_ctime_idx):
                    if ctime[ctime_idx_i] >= ctime[avg_idx][0] and ctime[ctime_idx_i] <= ctime[avg_idx][-1]:
                        avg_ctime_idx.pop(avg_idx_i)                 
        
        return avg_ctime_idx, spike_ctime_idx
    else:
        raise error.spikeError_spikcrosstime(ctime[idx])


def sepoch_avg(spike_ctime_idx, avg_ctime_idxs, B_x, B_y, B_z, ctime):
    """
    get Bx, By, Bz average from avg_spin
    """
    spike_ctime_diff = ctime[spike_ctime_idx] - ctime[spike_ctime_idx[0]]
    avg_Bx_sum = np.zeros(len(spike_ctime_diff))
    avg_By_sum = np.zeros(len(spike_ctime_diff))
    avg_Bz_sum = np.zeros(len(spike_ctime_diff))
    avg_Bx_count = np.zeros(len(spike_ctime_diff))
    avg_By_count = np.zeros(len(spike_ctime_diff))
    avg_Bz_count = np.zeros(len(spike_ctime_diff))
    
    for avg_ctime_idx in avg_ctime_idxs:
        avg_ctime_diff = ctime[avg_ctime_idx] - ctime[avg_ctime_idx][0]
        #Bplot.B_ctime_plot(ctime, B_x, B_y, B_z, cross_times=ctime[[avg_ctime_idx[0],avg_ctime_idx[-1]]], xlimt = [ctime[avg_ctime_idx[0]]-10, ctime[avg_ctime_idx[-1]]+10])
        #breakpoint()

        for spike_diff_i, spike_diff_time in enumerate(spike_ctime_diff):
            close_pos, close_value = find_closest(avg_ctime_diff, spike_diff_time) 
            if np.abs(close_value - spike_diff_time) < 0.05:
                avg_Bx_sum[spike_diff_i] = avg_Bx_sum[spike_diff_i] + B_x[avg_ctime_idx[close_pos]]
                avg_Bx_count[spike_diff_i] = avg_Bx_count[spike_diff_i] + 1
                avg_By_sum[spike_diff_i] = avg_By_sum[spike_diff_i] + B_y[avg_ctime_idx[close_pos]]
                avg_By_count[spike_diff_i] = avg_By_count[spike_diff_i] + 1
                avg_Bz_sum[spike_diff_i] = avg_Bz_sum[spike_diff_i] + B_z[avg_ctime_idx[close_pos]]
                avg_Bz_count[spike_diff_i] = avg_Bz_count[spike_diff_i] + 1
            
    avg_Bx = avg_Bx_sum / avg_Bx_count
    avg_By = avg_By_sum / avg_By_count
    avg_Bz = avg_Bz_sum / avg_Bz_count
    

    return avg_Bx, avg_By, avg_Bz, B_x[spike_ctime_idx], B_y[spike_ctime_idx], B_z[spike_ctime_idx]


def sepoch_sub(ctime, ctime_idx, spike_ctime_idx, avg_Bx, avg_By, avg_Bz, spike_Bx, spike_By, spike_Bz):
    """
    subtract avg and spike Bx, By, Bz and find t1, t2
    """
    diff_Bx = avg_Bx - spike_Bx
    diff_By = avg_By - spike_By
    diff_Bz = avg_Bz - spike_Bz
    diff_B2 = diff_Bx**2+diff_By**2+diff_Bz**2

    B_std = np.std(np.sort(diff_B2)[:int(len(diff_B2)*0.6)])
    B_spike_idx = np.where((diff_B2 > np.median(diff_B2) + B_std*4))[0]
    #Bplot.B_ctime_plot_single(ctime[spike_ctime_idx], diff_B2, scatter=True)
    if len(B_spike_idx) == 0:
        Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
            title="x1 = B_avg, x2 = B_spike", scatter = True)
        raise error.spikeError_t1t2(ctime[ctime_idx])

    B_spike_idx_chunk = np.where(np.diff(B_spike_idx) != 1)[0] + 1
    B_spike_idx_chunk = np.insert(B_spike_idx_chunk, 0, 0)
    B_spike_idx_chunk = np.insert(B_spike_idx_chunk, len(B_spike_idx_chunk), len(B_spike_idx))

    spike_idx1 = []
    spike_idx2 = []
    for B_spike_idx_chunk_i in range(1, len(B_spike_idx_chunk)):
        B_spike_idx_current = B_spike_idx[B_spike_idx_chunk[B_spike_idx_chunk_i-1]:B_spike_idx_chunk[B_spike_idx_chunk_i]]
        if ctime[spike_ctime_idx[B_spike_idx_current[0]]] <= ctime[ctime_idx] <= ctime[spike_ctime_idx[B_spike_idx_current[-1]]]: 
            spike_idx1 = B_spike_idx_current[0] - 1 
            spike_idx2 = B_spike_idx_current[-1]+ 1 
            if parameter.makeplot == True:
                Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
                    title="x1 = B_avg, x2 = B_spike", scatter = True, ctime_idx_time=ctime[spike_ctime_idx][[spike_idx1,spike_idx2]], ctime_idx_flag = [2, 2]
                )
      
    if spike_idx1 == []:
        Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
            title="x1 = B_avg, x2 = B_spike", scatter = True)
        raise error.spikeError_t1t2(ctime[ctime_idx])
        
    return spike_idx1, spike_idx2


def find_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.

    ref: https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value/12141511#12141511
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0, myList[0]
    if pos == len(myList):
        return len(myList)-1, myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return pos, after
    else:
        return pos-1, before


def spike_sinefit(ctime, Bx, By, Bz, spike_ctime_idxs):

    spike_idx1s = [spike_ctime_idx - 0 for spike_ctime_idx in spike_ctime_idxs] 
    spike_idx2s = [spike_ctime_idx + 2 for spike_ctime_idx in spike_ctime_idxs]
    for idx in range(len(spike_ctime_idxs)): 

        spike_idx1 = spike_idx1s[idx]
        spike_idx2 = spike_idx2s[idx]
        fit_opt, fit_covar = curve_fit(
            calibration.sine_fit2, ctime[spike_idx1 - parameter.spike_10hz_fit:spike_idx2 + parameter.spike_10hz_fit], 
            Bx[spike_idx1-parameter.spike_10hz_fit:spike_idx2+parameter.spike_10hz_fit],
            p0=[np.max(np.abs(Bx - np.mean(Bx))), 2.2, 0, 0]
        )
        B_S1_corr_del = calibration.sine_fit2(ctime, *fit_opt)
        Bx_del = calibration.sine_fit2(ctime[spike_idx1:spike_idx2], *fit_opt)
        Bx[spike_idx1:spike_idx2] = Bx_del

        fit_opt, fit_covar = curve_fit(
            calibration.sine_fit2, ctime[spike_idx1 - parameter.spike_10hz_fit:spike_idx2 + parameter.spike_10hz_fit], 
            By[spike_idx1 - parameter.spike_10hz_fit:spike_idx2 + parameter.spike_10hz_fit],
            p0=[np.max(np.abs(By - np.mean(By))), 2.2, 0, 0]
        )
        B_S2_corr_del = calibration.sine_fit2(ctime, *fit_opt)
        By_del = calibration.sine_fit2(ctime[spike_idx1:spike_idx2], *fit_opt)
        By[spike_idx1:spike_idx2] = By_del

        fit_opt, fit_covar = curve_fit(
            calibration.sine_fit2, ctime[spike_idx1 - parameter.spike_10hz_fit:spike_idx2 + parameter.spike_10hz_fit], 
            Bz[spike_idx1 - parameter.spike_10hz_fit:spike_idx2 + parameter.spike_10hz_fit],
            p0=[np.max(np.abs(Bz - np.mean(Bz))), 2.2, 0, 0]
        )
        B_S3_corr_del = calibration.sine_fit2(ctime, *fit_opt)
        Bz_del = calibration.sine_fit2(ctime[spike_idx1:spike_idx2], *fit_opt)
        Bz[spike_idx1:spike_idx2] = Bz_del

        if parameter.makeplot == True: 
            Bplot.B_ctime_plot(
                ctime, [Bx, B_S1_corr_del], [By, B_S2_corr_del], [Bz, B_S3_corr_del], xlimt = [(spike_idx1-28)/10, (spike_idx2+28)/10], scatter=True, title=f"del_rogue_10hz_{idx}")

    return Bx, By, Bz
