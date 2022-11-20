import datetime as dt
from pickle import FALSE
import numpy as np
import os
from typing import Literal


def get_state_cdf_path(mission: Literal["ela", "elb"], date: dt.date) -> str:
    """Get the path to the CDF corresponding to the given mission and date.

    NOTE: This function relies on the fact that `sp-server` will run on the
    sciproc-vm, which has an up-to-date copy of all elfin data.

    Parameters
    ----------
    mission : Literal["ela", "elb"]
    date: dt.date
    """
    year_str = str(date.year)
    sta_datestr = date.strftime("%Y%m%d")
    return os.path.join("/nfs/elfin-mirror/elfin", mission, "l1/state/defn", year_str, f"{mission}_l1_state_defn_{sta_datestr}_v02.cdf")


proper_pad = True  # fails when data have gaps
fit_running_spline = True
relative_integrate = True
bidirectional_integrate = True

eps_1 = 5
eps_2 = 5
eps_3 = 4
N_spins_fit = 4
peak_detect = False # For detecting peaks by fitting B_S3 itself instead of 
    #fitting its derivative and computing zero-crossings
zero_crossing_method = 3   
f = 44 * np.pi / 180
init_secs = 0 # seconds

"""ctime spike paramter
"""
# in ctime_calib, criterion for finding ctime gaps difference between jumps and ctime median in seconds
ctime_thrhld = 0.003

# in ctime_calib, degap unipolar gap 1/80s if true
ctime_correct_80 = True  # degap unipolar gap 1/80s if true

# in ctime_calib, some events has 1/80s spike doesn't need calibration starttime_str = ["2022-03-28/03:08:11"] 
ctime_correct_80_skip = True 

# in ctime_spike_80, num of spins in which to find [t1, t2] for 2.5 s spike
spike_find_spin_num = 4 

# in spike_sinefit_80,  number of points to fit around spike
spike_fit_len_80 = 28 

# in ctime_calib, degap unipolar gap < 0.01s if true
ctime_correct_100 = True

"""zero crossing stage 3
"""
# cross_time, ctime_spike
cross0_spike_del = True # 2.5s gap: before sine fit delete spins with spike in zero crossing determine stage 3
# init, if true delete fsp data points aroud spike for both 1/80s
fsp_spike_del_type2 = True # compare std with and without spike
# init, if true delete fsp data points aroud spike for both 2.5s
fsp_spike_del_type3 = True # compare std with and without spike
# init, if true delete fsp data points around purple spike
fsp_spike_del_type4 = True # compare std with and without spike

R2_filter = True  # in cross time determination stage 3, use R2_thrhld to exclude bad fit
R2_thrhld = 0.8

"""preprocess parameter
"""
# funky fgm
Spinrate_thrhld = 0.2  # if std of spin rate > median of spin rate * Spinrate_thrhld, funky fgm skip collection
# completeness over 100%
ctime_repeat_check = True

"""post calibration parameter
"""
cali_2nd = False  # add a second calibration in dmxl coordinate to solve the trend
cali_3rd = False  # add a third calibration 
cali_4th = False  # add a third calibration 
#del_rogue = False   # del rogue points in the first and last three points         
del_rogue_fsp = True # del rogue points in fsp resolution, which has a better result than del_rogue
eps_rogue = 3 # delete points outside med-std*eps_rogue and med+std*eps_rogue
fsp_detrend = True # detrend in dmxl if true
fsp_detrend_cutoff = 6 # detrend in dmxl if true

#del_spike_fsp = False  # delete spike in fsp resolution

"""output paramter
"""
makeplot = True
savepng = True
output = False # if true output to txt file
download_data = True

