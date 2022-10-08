import numpy as np
from .. import parameter 
from . import Bplot, error,  calibration
from bisect import bisect_left
from scipy.optimize import curve_fit


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