from audioop import maxpp
from calendar import c
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
    flag = 3 unipolar gap 2.5 s
    flag = 4 unipolar large gap
    flag = 5 unipolar negative spike
    ['red','orange','magenta','darkviolet', 'green']
    """
    ctime_idx_flag = np.zeros(len(ctime_idx), dtype = int)
    ctime_idx_timediff = np.zeros(len(ctime_idx))
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
            if np.abs(delta_dt - 0.0125) < 0.001 :
                """ unipolar spike with 1/80 s
                    if no calibration, then mark the orginal spike locaiton
                    if calibrate, then find the actual spike before the original one and mark 
                """
                if parameter.ctime_correct != True:  # if no calibration, then mark the orginal spike locaiton
                    spike_ctime_idxs.append(ctime_idx[i]) 
                    ctime_idx_flag[i] = 2 
                else:  
                    """
                    spike = 1/80s procedure:
                        1. get average B for three spins
                        2. determine true spike start and end with average B - spike B
                        3. sine fit
                    """
                    ctime[ctime_idx[i]+1:] = ctime[ctime_idx[i]+1:] - delta_dt
                    try:
                        # get the index for 3 spins with spike, and idx for other spins without spike
                        avg_ctime_idx, spike_ctime_idx = sepoch_getspin_80(ctime_idx[i], ctime, cross_times, ctime_idx) # index for three spins
                    
                        #if parameter.makeplot == True:
                        #    Bplot.B_ctime_plot(ctime, B_x, B_y, B_z, 
                        #    ctime_idx_time = ctime[ctime_idx],xlimt = [ctime[avg_ctime_idx[0][0]], ctime[avg_ctime_idx[-1][-1]]],
                        #    scatter = True, title = "before_spike_correction")

                        # get average Bx, By, Bz
                        avg_Bx, avg_By, avg_Bz, spike_Bx, spike_By, spike_Bz = sepoch_avg_80(spike_ctime_idx, avg_ctime_idx, B_x, B_y, B_z, ctime)
                        #Bplot.B_3d(spike_Bx, spike_By, spike_Bz)
                        #Bplot.B_ctime_plot_single(ctime[spike_ctime_idx], spike_Bx**2+spike_By**2+spike_Bz**2, scatter = True)
                        
                        # subtract Bx, By, Bz, get spike index
                        spike_idx1, spike_idx2 = sepoch_sub_80(ctime, ctime_idx[i], spike_ctime_idx, avg_Bx, avg_By, avg_Bz, spike_Bx, spike_By, spike_Bz)
                        
                        spike_ctime_idx1 = spike_ctime_idx[spike_idx1]
                        spike_ctime_idx2 = spike_ctime_idx[spike_idx2]
                    except (error.spikeError80_t1t2, error.spikeError80_spikespin, error.spikeError80_spikcrosstime) as e:
                        
                        logger.error(e.__str__())
                        spike_ctime_idx1 = ctime_idx[i]

                    # move delta_t to t1
                    ctime[spike_ctime_idx1:] = ctime[spike_ctime_idx1:] + delta_dt
                    spike_ctime_idxs.append(spike_ctime_idx1)
                    ctime_idx_flag[i] = 2 
            elif np.abs((np.abs(delta_dt) - 0.01)) < 0.01:
                """
                a negative spike a little bit smaller than 1/80s, example:     
                                starttime_str = ["2022-01-05/13:26:12"] 
                                endtime_str = ["2022-01-05/13:32:24"]
                        this type spike starts at ctime[citme_idx], no need to find true start
                """
                ctime[ctime_idx[i]+1:] = ctime[ctime_idx[i]+1:] - delta_dt
                ctime_idx_flag[i] = 5
                ctime_idx_timediff[i] = delta_dt

            elif np.abs(delta_dt) > 2.4 and  np.abs(delta_dt) < 2.7 : # gaps
            #elif np.abs(delta_dt) > 1.4 and  np.abs(delta_dt) < 1.6 : # gaps
                ctime_idx_flag[i] = 3
                ctime_idx_timediff[i] = delta_dt
            else:
                ctime_idx_flag[i] = 4
                ctime_idx_timediff[i] = delta_dt
            i += 1


    return ctime, ctime_idx, ctime_idx_flag, ctime_idx_timediff, spike_ctime_idxs    


def sepoch_getspin_80(idx, ctime, cross_times, ctime_idx):
    """
    calibration step 1 for 1/80 s 
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
            raise error.spikeError80_spikespin(ctime[idx])
        
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
        raise error.spikeError80_spikcrosstime(ctime[idx])


def sepoch_avg_80(spike_ctime_idx, avg_ctime_idxs, B_x, B_y, B_z, ctime):
    """
    calibration step 2 for 1/80 s 
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


def sepoch_sub_80(ctime, ctime_idx, spike_ctime_idx, avg_Bx, avg_By, avg_Bz, spike_Bx, spike_By, spike_Bz):
    """
    calibration step 3 for 1/80 s 
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
        raise error.spikeError80_t1t2(ctime[ctime_idx])
    
    # B_spike_idx_chunk: index of chunk
    B_spike_idx_chunk = np.where(np.diff(B_spike_idx) != 1)[0] + 1
    B_spike_idx_chunk = np.insert(B_spike_idx_chunk, 0, 0)
    B_spike_idx_chunk = np.insert(B_spike_idx_chunk, len(B_spike_idx_chunk), len(B_spike_idx))

    spike_idx1 = []
    spike_idx2 = []
    for B_spike_idx_chunk_i in range(1, len(B_spike_idx_chunk)):
        B_spike_idx_current = B_spike_idx[B_spike_idx_chunk[B_spike_idx_chunk_i-1]:B_spike_idx_chunk[B_spike_idx_chunk_i]]
        # if chunk start < spike time < chunk end
        if ctime[spike_ctime_idx[B_spike_idx_current[0]]] <= ctime[ctime_idx] <= ctime[spike_ctime_idx[B_spike_idx_current[-1]]]: 
            spike_idx1 = B_spike_idx_current[0] - 1 
            spike_idx2 = B_spike_idx_current[-1]+ 1 
            if parameter.makeplot == True:
                Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
                    title="x1 = B_avg, x2 = B_spike", scatter = True, ctime_idx_time=ctime[spike_ctime_idx][[spike_idx1,spike_idx2]], ctime_idx_flag = [2, 2]
                )
      
    if spike_idx1 == []:

        if parameter.makeplot == True:
            Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
                title=f"{ctime_idx}_x1 = B_avg, x2 = B_spike", scatter = True)
        raise error.spikeError80_t1t2(ctime[ctime_idx])
      
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
    """
    calibration step 4 for 1/80 s 
    """
    sine_fit2 = lambda x, A, w, p, k: calibration.sine_fit(x, 1, A, w, p, k)
    spike_idx1s = [spike_ctime_idx - 0 for spike_ctime_idx in spike_ctime_idxs] 
    spike_idx2s = [spike_ctime_idx + 2 for spike_ctime_idx in spike_ctime_idxs]
    for idx in range(len(spike_ctime_idxs)): 

        spike_idx1 = spike_idx1s[idx]
        spike_idx2 = spike_idx2s[idx]
        fit_opt, fit_covar = curve_fit(
            sine_fit2, ctime[spike_idx1 - parameter.spike_10hz_fit:spike_idx2 + parameter.spike_10hz_fit], 
            Bx[spike_idx1-parameter.spike_10hz_fit:spike_idx2+parameter.spike_10hz_fit],
            p0=[np.max(np.abs(Bx - np.mean(Bx))), 2.2, 0, 0]
        )
        B_S1_corr_del = sine_fit2(ctime, *fit_opt)
        Bx_del = sine_fit2(ctime[spike_idx1:spike_idx2], *fit_opt)
        Bx[spike_idx1:spike_idx2] = Bx_del

        fit_opt, fit_covar = curve_fit(
            sine_fit2, ctime[spike_idx1 - parameter.spike_10hz_fit:spike_idx2 + parameter.spike_10hz_fit], 
            By[spike_idx1 - parameter.spike_10hz_fit:spike_idx2 + parameter.spike_10hz_fit],
            p0=[np.max(np.abs(By - np.mean(By))), 2.2, 0, 0]
        )
        B_S2_corr_del = sine_fit2(ctime, *fit_opt)
        By_del = sine_fit2(ctime[spike_idx1:spike_idx2], *fit_opt)
        By[spike_idx1:spike_idx2] = By_del

        fit_opt, fit_covar = curve_fit(
            sine_fit2, ctime[spike_idx1 - parameter.spike_10hz_fit:spike_idx2 + parameter.spike_10hz_fit], 
            Bz[spike_idx1 - parameter.spike_10hz_fit:spike_idx2 + parameter.spike_10hz_fit],
            p0=[np.max(np.abs(Bz - np.mean(Bz))), 2.2, 0, 0]
        )
        B_S3_corr_del = sine_fit2(ctime, *fit_opt)
        Bz_del = sine_fit2(ctime[spike_idx1:spike_idx2], *fit_opt)
        Bz[spike_idx1:spike_idx2] = Bz_del

        if parameter.makeplot == True: 
            Bplot.B_ctime_plot(
                ctime, [Bx, B_S1_corr_del], [By, B_S2_corr_del], [Bz, B_S3_corr_del], xlimt = [(spike_idx1-28)/10, (spike_idx2+28)/10], scatter=True, title=f"del_rogue_10hz_{idx}")

    return Bx, By, Bz


def spike_fill(ctime, cross_times, ctime_idx, ctime_idx_flag, B_x, B_y, B_z, B_x_igrf, B_y_igrf, B_z_igrf, att_x, att_y, att_z):
    """
    calibration step 0 for 2.5 s gap
    cannot put in ctime_calib because ctime_idx has to be shift
    """
    for i in range(len(ctime_idx)):
        if ctime_idx_flag[i] == 3:
            spike_ctime = np.arange(ctime[ctime_idx[i]]+0.1, ctime[ctime_idx[i]+1], 0.1)
            # step 1: get avg spin idx, this is the ctime index of 2.5 s segement in 10 spins 
            avg_ctime_idxs = sepoch_getspin_25(ctime_idx[i], cross_times, ctime, ctime_idx, spike_ctime)
            if len(avg_ctime_idxs) == 0:
                breakpoint()
            # step 2: average and fill in the gap and shift all data    
            ctime, ctime_idx, B_x, B_y, B_z, B_x_igrf, B_y_igrf, B_z_igrf, att_x, att_y, att_z = sepoch_avg_25(
                ctime_idx[i], avg_ctime_idxs, spike_ctime, ctime, ctime_idx, 
                B_x, B_y, B_z, B_x_igrf, B_y_igrf, B_z_igrf, att_x, att_y, att_z
            )
            if len(ctime) != len(B_x):
                print("not same length")
                breakpoint()
    return ctime, ctime_idx, B_x, B_y, B_z, B_x_igrf, B_y_igrf, B_z_igrf, att_x, att_y, att_z


def sepoch_getspin_25(spike_ctime_idx, cross_times, ctime, ctime_idx, spike_ctime):
    """
    calibration step 1 for 2.5 s gap
    get indexes for 10 segements of 2.5 s data 
    """
    # get spin index in ct for spike
    spike_spin_idx2 = np.where(cross_times>ctime[spike_ctime_idx])[0] # index of spin cross zero after spike
    if spike_spin_idx2.size != 0 :
        spike_spin_idx2 = spike_spin_idx2[0] # ct index right after spike
        if spike_ctime[-1] > cross_times[spike_spin_idx2]:  # spike start in one spin. 
            spike_spin_idx2 = spike_spin_idx2 + 1
            spike_find_spin_num = 2
        else:
            spike_find_spin_num = 2
        spike_spin_idx1 = spike_spin_idx2 - spike_find_spin_num # ct index before spike

        spike_ctime_idx1, spike_ctime_time1 = find_closest(ctime, cross_times[spike_spin_idx1]) # ctime idx for ct after spike
        spike_ctime_idx2, spike_ctime_time2 = find_closest(ctime, cross_times[spike_spin_idx2]) # ctime idx for ct before spike
        if np.abs(spike_ctime_time1 - cross_times[spike_spin_idx1]) < 0.5 and np.abs(spike_ctime_time2 - cross_times[spike_spin_idx2]) < 0.5:
            spike_ctime_idx = [*range(spike_ctime_idx1, spike_ctime_idx2)]
        else:
            raise error.spikeError25_spikespin(ctime[spike_ctime_idx])

        # get spike relative time in spike spin
        spike_rela_time1 = spike_ctime[0] - spike_ctime_time1
        spike_rela_time2 = spike_ctime[-1] - spike_ctime_time1

        # get index for average spins
        avg_spin_idx1 = spike_spin_idx1 - 5 if spike_spin_idx1 > 5 else 0 # left end of the avg spin start
        avg_spin_idx2 = spike_spin_idx1 + 5 if spike_spin_idx2 < len(cross_times) - 5 else len(cross_times) - spike_find_spin_num # right end of the avg spin start

        avg_ctime_idxs = []
        for spin_i in range(avg_spin_idx1, avg_spin_idx2):
            avg_ctime_idx1, avg_ctime_time1 = find_closest(ctime, cross_times[spin_i]) # ctime idx for start of this avg spin
            avg_ctime_idx2, avg_ctime_time2 = find_closest(ctime, cross_times[spin_i + spike_find_spin_num]) # ctime idx for end of this avg spin
            avg_ctime_idx_range = range(avg_ctime_idx1, avg_ctime_idx2)
            if np.abs(avg_ctime_time1 - cross_times[spin_i]) < 0.1 and np.abs(avg_ctime_time2 - cross_times[spin_i + spike_find_spin_num]) < 0.1:
                avg_rela_times = ctime[avg_ctime_idx_range] - avg_ctime_time1  # get relative time to ct
                avg_rela_time_idx1 = find_closest(avg_rela_times, spike_rela_time1)[0]  # idx for spike_rela_time1 in avg_rela_times
                avg_rela_time_idx2 = find_closest(avg_rela_times, spike_rela_time2)[0] + 1  # idx for spike_rela_time1 in avg_rela_times
                if avg_rela_time_idx2 - avg_rela_time_idx1  == len(spike_ctime): 
                    avg_ctime_idx1 = avg_ctime_idx_range[avg_rela_time_idx1] # ctime idx for 2.5s average segment
                    avg_ctime_idx2 = avg_ctime_idx_range[avg_rela_time_idx2]
                    # if spin includes other spikes, delete it
                    flag = 0
                    for ctime_idx_i in ctime_idx:
                        if ctime[avg_ctime_idx1] <= ctime[ctime_idx_i] <= ctime[avg_ctime_idx2]:
                            flag = 1

                    if flag == 0:
                        avg_ctime_idxs.append([*range(avg_ctime_idx1, avg_ctime_idx2)])              
        
        return avg_ctime_idxs
    else:
        raise error.spikeError25_spikcrosstime(ctime[spike_ctime_idx])


def sepoch_avg_25(
    spike_ctime_idx, avg_ctime_idxs, spike_ctime, ctime, ctime_idx, 
    Bx, By, Bz, Bx_igrf, By_igrf, Bz_igrf, att_x, att_y, att_z):
    """
    calibration step 2 for 2.5 s gap
    get average Bx, By, Bz 
    replace the ctime and Bx, By, Bz with fill in
    """
    avg_Bx = np.average(Bx[avg_ctime_idxs], axis = 0)
    avg_By = np.average(By[avg_ctime_idxs], axis = 0)
    avg_Bz = np.average(Bz[avg_ctime_idxs], axis = 0)
    avg_Bx_igrf = np.average(Bx_igrf[avg_ctime_idxs], axis = 0)
    avg_By_igrf = np.average(By_igrf[avg_ctime_idxs], axis = 0)
    avg_Bz_igrf = np.average(Bz_igrf[avg_ctime_idxs], axis = 0)
    avg_att_x = np.average(att_x[avg_ctime_idxs], axis = 0)
    avg_att_y = np.average(att_y[avg_ctime_idxs], axis = 0)
    avg_att_z = np.average(att_z[avg_ctime_idxs], axis = 0)

    ctime_fill = np.insert(ctime, spike_ctime_idx+1, spike_ctime)
    Bx_fill = np.insert(Bx, spike_ctime_idx+1, avg_Bx)
    By_fill = np.insert(By, spike_ctime_idx+1, avg_By)
    Bz_fill = np.insert(Bz, spike_ctime_idx+1, avg_Bz)
    Bx_igrf_fill = np.insert(Bx_igrf, spike_ctime_idx+1, avg_Bx_igrf)
    By_igrf_fill = np.insert(By_igrf, spike_ctime_idx+1, avg_By_igrf)
    Bz_igrf_fill = np.insert(Bz_igrf, spike_ctime_idx+1, avg_Bz_igrf)
    att_x_fill = np.insert(att_x, spike_ctime_idx+1, avg_att_x)
    att_y_fill = np.insert(att_y, spike_ctime_idx+1, avg_att_y)
    att_z_fill = np.insert(att_z, spike_ctime_idx+1, avg_att_x)

    #Bplot.B_ctime_plot(ctime_fill, Bx_fill, By_fill, Bz_fill, 
    #    ctime_idx_time=ctime[ctime_idx], xlimt = [ctime[spike_ctime_idx]-15, ctime[spike_ctime_idx]+15], cross_times=[ctime[avg_ctime_idxs[0]][0], ctime[avg_ctime_idxs[0]][-1]])
    #breakpoint()
    ctime_idx_fill = [ i+len(spike_ctime) if i > spike_ctime_idx else i for i in ctime_idx ]

    return [ctime_fill, ctime_idx_fill, Bx_fill, By_fill, Bz_fill, 
        Bx_igrf_fill, By_igrf_fill, Bz_igrf_fill, att_x_fill, att_y_fill, att_z_fill]

def getidx_spike_fsp_25(ctime, ctime_idx, ctime_idx_flag, cross_times):
    """
    get cross time index for 2.5s spike in fsp resolution
    """
    spike25_fsp_idx = []
    for ctime_idx_i in range(len(ctime_idx)):
        if ctime_idx_flag[ctime_idx_i] == 3 :
            idx = np.where(cross_times > ctime[ctime_idx[ctime_idx_i]])[0]
            if len(idx) != 0 :
                spike25_fsp_idx.append(idx[0])

    return spike25_fsp_idx
