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
    new version so only check the chunk start time > ctime_idx
    """
    diff_Bx = avg_Bx - spike_Bx
    diff_By = avg_By - spike_By
    diff_Bz = avg_Bz - spike_Bz
    diff_B2 = diff_Bx**2+diff_By**2+diff_Bz**2

    B_std = np.std(np.sort(diff_B2)[:int(len(diff_B2)*0.6)])
    B_spike_idx = np.where((diff_B2 > np.median(diff_B2) + B_std*4))[0]

    #Bplot.B_ctime_plot_single(ctime[spike_ctime_idx], diff_B2, scatter=True)
    if len(B_spike_idx) == 0:
        if parameter.makeplot == True:
            Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
                title=f"{ctime_idx}_x1 = B_avg, x2 = B_spike", scatter = True)
        raise error.spikeError80_t1t2(ctime[ctime_idx])

    # B_spike_idx_chunk: index of chunk
    B_spike_idx_chunk = np.where(np.diff(B_spike_idx) != 1)[0] + 1
    B_spike_idx_chunk = np.insert(B_spike_idx_chunk, 0, 0)
    B_spike_idx_chunk = np.insert(B_spike_idx_chunk, len(B_spike_idx_chunk), len(B_spike_idx))

    spike_ctime_idx1 = []
    spike_ctime_idx2 = []
    """find chunk start < spike < chunk end
    """
    for B_spike_idx_chunk_i in range(1, len(B_spike_idx_chunk)):
        B_spike_idx_current = B_spike_idx[B_spike_idx_chunk[B_spike_idx_chunk_i-1]:B_spike_idx_chunk[B_spike_idx_chunk_i]]
        # if chunk start < spike time < chunk end
        if ctime[ctime_idx] - 2.8 <= ctime[spike_ctime_idx[B_spike_idx_current[0]]] <= ctime[ctime_idx] : 
            spike_ctime_idx1 = spike_ctime_idx[B_spike_idx_current[0]] - 1 
            if parameter.makeplot == True:
                Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
                    title=f"{ctime_idx}_x1 = B_avg, x2 = B_spike", scatter = True, ctime_idx_time=ctime[[spike_ctime_idx1, ctime_idx]], ctime_idx_flag = [5, 2]
                ) 
    if spike_ctime_idx1 == []:
        if parameter.makeplot == True:
            Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
                title=f"{ctime_idx}_x1 = B_avg, x2 = B_spike", ctime_idx_time=ctime[ctime_idx], ctime_idx_flag=[2], scatter = True)
        raise error.spikeError80_t1t2(ctime[ctime_idx])
      
    return spike_ctime_idx1, spike_ctime_idx2


def sepoch_sub_80_0(ctime, ctime_idx, spike_ctime_idx, avg_Bx, avg_By, avg_Bz, spike_Bx, spike_By, spike_Bz):
    """
    calibration step 3 for 1/80 s 
    subtract avg and spike Bx, By, Bz and find t1, t2
    old version, requires chunk start < ctime_idx < chunk end
    """
    diff_Bx = avg_Bx - spike_Bx
    diff_By = avg_By - spike_By
    diff_Bz = avg_Bz - spike_Bz
    diff_B2 = diff_Bx**2+diff_By**2+diff_Bz**2

    B_std = np.std(np.sort(diff_B2)[:int(len(diff_B2)*0.6)])
    B_spike_idx = np.where((diff_B2 > np.median(diff_B2) + B_std*4))[0]

    #Bplot.B_ctime_plot_single(ctime[spike_ctime_idx], diff_B2, scatter=True)
    if len(B_spike_idx) == 0:
        if parameter.makeplot == True:
            Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
                title="x1 = B_avg, x2 = B_spike", scatter = True)
        raise error.spikeError80_t1t2(ctime[ctime_idx])
    
    # B_spike_idx_chunk: index of chunk
    B_spike_idx_chunk = np.where(np.diff(B_spike_idx) != 1)[0] + 1
    B_spike_idx_chunk = np.insert(B_spike_idx_chunk, 0, 0)
    B_spike_idx_chunk = np.insert(B_spike_idx_chunk, len(B_spike_idx_chunk), len(B_spike_idx))

    spike_idx1 = []
    spike_idx2 = []
    """find chunk start < spike < chunk end
    """
    for B_spike_idx_chunk_i in range(1, len(B_spike_idx_chunk)):
        B_spike_idx_current = B_spike_idx[B_spike_idx_chunk[B_spike_idx_chunk_i-1]:B_spike_idx_chunk[B_spike_idx_chunk_i]]
        # if chunk start < spike time < chunk end
        if ctime[spike_ctime_idx[B_spike_idx_current[0]]] <= ctime[ctime_idx] <= ctime[spike_ctime_idx[B_spike_idx_current[-1]]]: 
            spike_idx1 = B_spike_idx_current[0] - 1 
            spike_idx2 = B_spike_idx_current[-1]+ 1 
            if parameter.makeplot == True:
                Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
                    title=f"{ctime_idx}_x1 = B_avg, x2 = B_spike", scatter = True, ctime_idx_time=ctime[spike_ctime_idx][[spike_idx1,spike_idx2]], ctime_idx_flag = [2, 2]
                ) 

    if spike_idx1 == []:
        """if chunk around spike not found, use the cloest chunk before spike
        """
        for B_spike_idx_chunk_i in range(1, len(B_spike_idx_chunk)):
            B_spike_idx_current = B_spike_idx[B_spike_idx_chunk[B_spike_idx_chunk_i-1]:B_spike_idx_chunk[B_spike_idx_chunk_i]]
            if ctime[spike_ctime_idx[B_spike_idx_current[0]]] <= ctime[ctime_idx]: 
                spike_idx1 = B_spike_idx_current[0] - 1 
                spike_idx2 = B_spike_idx_current[-1]+ 1 
                if parameter.makeplot == True:
                    Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
                        title=f"{ctime_idx}_x1 = B_avg, x2 = B_spike", scatter = True, ctime_idx_time=ctime[spike_ctime_idx][[spike_idx1,spike_idx2]], ctime_idx_flag = [2, 2]
                    ) 

    if spike_idx1 == []:
        if parameter.makeplot == True:
            Bplot.B_ctime_plot(ctime[spike_ctime_idx], [avg_Bx, spike_Bx], [avg_By, spike_By], [avg_Bz, spike_Bz], 
                title=f"{ctime_idx}_x1 = B_avg, x2 = B_spike", ctime_idx_time=[ctime[ctime_idx]], scatter = True)
        raise error.spikeError80_t1t2(ctime[ctime_idx])
    
    return spike_idx1, spike_idx2


def spike_sinefit_80(ctime, Bx, By, Bz, spike_ctime_idxs):
    """
    calibration step 4 for 1/80 s 
    """
    sine_fit2 = lambda x, A, w, p, k: calibration.sine_fit(x, 1, A, w, p, k)
    spike_idx1s = [spike_ctime_idx - 0 for spike_ctime_idx in spike_ctime_idxs] 
    spike_idx2s = [spike_ctime_idx + 2 for spike_ctime_idx in spike_ctime_idxs]
    for idx in range(len(spike_ctime_idxs)): 
        try:
            spike_idx1 = spike_idx1s[idx]
            spike_idx2 = spike_idx2s[idx]
            fit_opt, fit_covar = curve_fit(
                sine_fit2, ctime[spike_idx1 - parameter.spike_fit_len_80:spike_idx2 + parameter.spike_fit_len_80], 
                Bx[spike_idx1-parameter.spike_fit_len_80:spike_idx2+parameter.spike_fit_len_80],
                p0=[np.max(np.abs(Bx - np.mean(Bx))), 2.2, 0, 0]
            )
            B_S1_corr_del = sine_fit2(ctime, *fit_opt)
            Bx_del = sine_fit2(ctime[spike_idx1:spike_idx2], *fit_opt)
            Bx[spike_idx1:spike_idx2] = Bx_del

            fit_opt, fit_covar = curve_fit(
                sine_fit2, ctime[spike_idx1 - parameter.spike_fit_len_80:spike_idx2 + parameter.spike_fit_len_80], 
                By[spike_idx1 - parameter.spike_fit_len_80:spike_idx2 + parameter.spike_fit_len_80],
                p0=[np.max(np.abs(By - np.mean(By))), 2.2, 0, 0]
            )
            B_S2_corr_del = sine_fit2(ctime, *fit_opt)
            By_del = sine_fit2(ctime[spike_idx1:spike_idx2], *fit_opt)
            By[spike_idx1:spike_idx2] = By_del

            fit_opt, fit_covar = curve_fit(
                sine_fit2, ctime[spike_idx1 - parameter.spike_fit_len_80:spike_idx2 + parameter.spike_fit_len_80], 
                Bz[spike_idx1 - parameter.spike_fit_len_80:spike_idx2 + parameter.spike_fit_len_80],
                p0=[np.max(np.abs(Bz - np.mean(Bz))), 2.2, 0, 0]
            )
            B_S3_corr_del = sine_fit2(ctime, *fit_opt)
            Bz_del = sine_fit2(ctime[spike_idx1:spike_idx2], *fit_opt)
            Bz[spike_idx1:spike_idx2] = Bz_del

            """
            if parameter.makeplot == True: 
                Bplot.B_ctime_plot(
                    ctime, [Bx, B_S1_corr_del], [By, B_S2_corr_del], [Bz, B_S3_corr_del], xlimt = [(spike_idx1-28)/10, (spike_idx2+28)/10], scatter=True, title=f"del_rogue_10hz_{idx}")
            """
        except:
            pass
    return Bx, By, Bz
