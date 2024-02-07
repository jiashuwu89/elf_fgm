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

elfin_url = "https://data.elfin.ucla.edu/"

proper_pad = False  # fails when data have gaps
fit_running_spline = False
relative_integrate = False
bidirectional_integrate = False

init_secs = 0 # inital seconds to exclude
funkyfgm = False
f = (90-55) * np.pi / 180

"""zero crossing parameter
"""
eps_1 = 1
eps_2 = 2
eps_3 = 2
N_spins_fit = 4
peak_detect = False # For detecting peaks by fitting B_S3 itself instead of 
    #fitting its derivative and computing zero-crossings
zero_crossing_method = 3   

"""ctime spike paramter
"""
# in ctime_calib, criterion for finding ctime gaps difference between jumps and ctime median in seconds
ctime_thrhld = 0.003

# in ctime_calib, degap unipolar gap 1/80s if true
ctime_correct_80 = False  # degap unipolar gap 1/80s if true

# in ctime_calib, some events has 1/80s spike doesn't need calibration starttime_str = ["2022-03-28/03:08:11"] 
ctime_correct_80_skip = False 

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
cali_2nd = False  # add a second calibration 
cali_3rd = False # add a third calibration 
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

"""specify time interval for some events
"""
del_time = False
del_time_idxstart = [260] # 2022-07-20/13:32:24
del_time_idxend = [306]

"""gei to obw
"""
gei2obw = True

"""change zero crossing location according to omege, in phase interation  
"""
CrossTime_Update = False

"""add boundary to least square fitting, in calibration
"""
fit_bound = True

"""loop f, if set ignore the f value above 
"""
f_loop = False
f_loop_value = list(map(lambda x: x * (np.pi / 180), range(0, 360, 100)))

"""att rotation by a angle around one axis
"""
att_rot = False
att_rot_ang = 0.5 # deg
att_rot_axis = 'z'

"""loop of attitude around the original attitude by an angle
"""
att_loop = False
att_loop_width = 5 # deg
att_loop_step = 1 # step when rotate from 0 to 360
att_loop_figure = False # plot attitude vector and rotated vectors in 3d plot

"""use att in txt file. this usually require output att to txt first
"""
att_csv = False
att_csv_filename = "fgm_utils/att_csv/20210317_1829_1851_att_gei.csv"

"""run mva during each run
"""
mva_fgm = False

"""use angle from mva as rotation angle, works only when mva_fgm is true
"""
mvaang_rotang = False

"""print out B parameter in each calibration
"""
Bpara_out = True

"""attitude determination
"""
att_determine = False

"""attitude split
"""
att_split = False
# if this is set, 
# it will divide interval into equal length snippets, 
# if not set, will check att_split_idx
att_split_num = None
# start time of each snippet, 
# if this is set, will use the setted ind to divided att, 
# if not set, will use ind for sci zone
#att_split_idx = [0, 13000, 26599, 46000]  # event 1
att_split_idx = None
att_split_detrend = False  # will detrend attitude before fitting

"""run batch sci zones from fgm_availablity.csv
"""
batch_run = False

"""spin rate fit
"""
wfit_run = True
wfit_run_maxiter = 1
wfit_gradient_figure = True  ##genreate the figure for gradient 
wfit_gradient_choice_lst = {
    1: 'cross_time0', 
    2: 'wfit_m', 
    3: 'wfit_c', 
    4: 't',
    }
wfit_gradient_choice = 4

"""calibration in dmxl
"""
cal_dmxl = False