import numpy as np
from .. import parameter 
from . import error, ctime_spike_80
from bisect import bisect_left
from datetime import datetime, timedelta

def ctime_calib(ctime, B_x, B_y , B_z, cross_times, logger = None, datestr = None):
    """This function find spike in ctime
    flag = 1 multiploar spike, correct here. 
        color:
            red

    flag = 2 unipolar spike 1/80 s, 0.0125 +/- 0.001 = 0.0115 - 0.0135
        paramter:
            ctime_correct_80
        color:
            orange
        procedure:
            2.1 get average B for three spins
            2.2 determine true spike start and end with average B - spike B
            2.3 sine fit (not here, in init.py)

    flag = 3 unipolar gap 2.5 s or larger
        color:
            magenta
        parameter:
            cross0_spike_del 
        procedure: 
            3.1 in cross_time stage 3, exclude spins with spikes when do sine fit  
            3.2 in cross_time stage 3, exclude zero crossing adjacent to spike

    flag = 4 all other unipolar spike
        color:
            darkviolet

    flag = 5 np.abs(delta_dt) < 0.0115 
        color:
            green
        parameter: 
            ctime_correct_100
        procedure:
            5.1 correct delta after spike

    """
    delta_t = np.median(ctime[1:]-ctime[:-1])
    ctime_adj = ctime[1:]-ctime[:-1] - delta_t
    ctime_idx = np.where(np.abs(ctime_adj) > parameter.ctime_thrhld)[0]

    #if parameter.makeplot == True:
    #    Bplot.ctimediff_plot(ctime, ctime_idx = ctime_idx)
    #breakpoint()

    ctime_idx_flag = np.zeros(len(ctime_idx), dtype = int)
    ctime_idx_timediff = ctime_adj[ctime_idx]
    spike_ctime_idxs = []  # save 1/80 spike
    i = 0
    while i < len(ctime_idx):
        if (i < len(ctime_idx)-1 and ctime_adj[ctime_idx[i]] > 0 and ctime_adj[ctime_idx[i+1]] < 0 # multipolar jumps
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
            if np.abs(delta_dt - 0.0125) < 0.001 or np.abs(delta_dt - 0.0125*2) < 0.001:
                """ unipolar spike with 1/80 s
                    if no calibration, then mark the orginal spike locaiton
                    if calibrate, then find the actual spike before the original one and mark 
                """
                if parameter.ctime_correct_80 != True:  # if no calibration, then mark the orginal spike locaiton
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
                        avg_ctime_idx, spike_ctime_idx = ctime_spike_80.sepoch_getspin_80(ctime_idx[i], ctime, cross_times, ctime_idx) # index for three spins
                    
                        #if parameter.makeplot == True:
                        #    Bplot.B_ctime_plot(ctime, B_x, B_y, B_z, 
                        #    ctime_idx_time = ctime[ctime_idx],xlimt = [ctime[avg_ctime_idx[0][0]], ctime[avg_ctime_idx[-1][-1]]],
                        #    scatter = True, title = "before_spike_correction")

                        # get average Bx, By, Bz
                        avg_Bx, avg_By, avg_Bz, spike_Bx, spike_By, spike_Bz = ctime_spike_80.sepoch_avg_80(spike_ctime_idx, avg_ctime_idx, B_x, B_y, B_z, ctime)
                        #Bplot.B_3d(spike_Bx, spike_By, spike_Bz)
                        #Bplot.B_ctime_plot_single(ctime[spike_ctime_idx], spike_Bx**2+spike_By**2+spike_Bz**2, scatter = True)
                        
                        #idx_beforespike = int(np.where(spike_ctime_idx == ctime_idx[i])[0])+10
                        idx_beforespike = len(spike_ctime_idx)
                        # subtract Bx, By, Bz, get spike index
                        spike_ctime_idx1, spike_ctime_idx2 = ctime_spike_80.sepoch_sub_80(
                            ctime, ctime_idx[i], spike_ctime_idx[:idx_beforespike], 
                            avg_Bx[:idx_beforespike], avg_By[:idx_beforespike], avg_Bz[:idx_beforespike], 
                            spike_Bx[:idx_beforespike], spike_By[:idx_beforespike], spike_Bz[:idx_beforespike]
                        )
    
                        if parameter.ctime_correct_80_skip is True:
                            """TODO: delete 1/80s spike that has worse result
                            """
                            spike_ctime_idx1 = spike_skip_80(ctime, ctime_idx[i], datestr, spike_ctime_idx1)
                        # move delta_t to t1
                        ctime[spike_ctime_idx1:] = ctime[spike_ctime_idx1:] + delta_dt
                        spike_ctime_idxs.append(spike_ctime_idx1)
                        ctime_idx_flag[i] = 2 
                    except (error.spikeError80_t1t2, error.spikeError80_spikespin, error.spikeError80_spikcrosstime) as e:
                        #ctime[ctime_idx[i]+1:] = ctime[ctime_idx[i]+1:] + delta_dt
                        logger.error(e.__str__())
                        ctime_idx_flag[i] = 2
                        spike_ctime_idx1 = ctime_idx[i]
                        spike_ctime_idxs.append(spike_ctime_idx1)
                        ctime[spike_ctime_idx1:] = ctime[spike_ctime_idx1:] + delta_dt

            elif np.abs(delta_dt) < 0.0115:
                """
                a negative spike a little bit smaller than 1/80s, example:     
                                starttime_str = ["2022-01-05/13:26:12"] 
                                endtime_str = ["2022-01-05/13:32:24"]
                        this type spike starts at ctime[citme_idx], no need to find true start
                """
                if parameter.ctime_correct_100 != True :
                    ctime[ctime_idx[i]+1:] = ctime[ctime_idx[i]+1:] - delta_dt
                ctime_idx_flag[i] = 5

            elif np.abs(delta_dt) > 2.4 : # gaps
                ctime_idx_flag[i] = 3
            else:
                ctime_idx_flag[i] = 4
            i += 1

    return ctime, ctime_idx, ctime_idx_flag, ctime_idx_timediff, spike_ctime_idxs    


def spike_skip_80(ctime, ctime_idx, datestr, spike_ctime_idx1):
    """These are bad correction for 1/80s 
    """
    ctime_datetime = list(map(lambda ts: datetime.strptime(datestr, '%Y%m%d_%H%M%S') + timedelta(seconds=ts), ctime))
    if (
        datetime(2022, 3, 28, 3, 12, 36, 200000) < ctime_datetime[ctime_idx] < datetime(2022, 3, 28, 3, 12, 36, 400000) or
        datetime(2022, 3, 28, 3, 9, 28, 600000) < ctime_datetime[ctime_idx] < datetime(2022, 3, 28, 3, 9, 28, 900000) or
        datetime(2022, 3, 26, 20, 16, 30, 600000) < ctime_datetime[ctime_idx] < datetime(2022, 3, 26, 20, 16, 30, 900000) or
        datetime(2022, 3, 26, 15, 41, 47, 600009) < ctime_datetime[ctime_idx] < datetime(2022, 3, 26, 15, 41, 48, 000000) or 
        datetime(2022, 3, 16, 16, 32, 52, 299991) < ctime_datetime[ctime_idx] < datetime(2022, 3, 16, 16, 32, 52, 699991)
        ):
        return ctime_idx
    elif(datetime(2022, 3, 10, 4, 22, 36, 200000) < ctime_datetime[ctime_idx] < datetime(2022, 3, 10, 4, 22, 36, 600000)):
        idx = bisect_left(ctime_datetime, datetime(2022, 3, 10, 4, 22, 35, 800000))
        return idx
    else:
        return spike_ctime_idx1
