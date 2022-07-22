import numpy as np

STATE_DATA_DIR = "/Users/jamesking/Desktop/elfin/OPS/science/gitlab/sp-server/sp/fgm_utils/test"

proper_pad = False  # fails when data have gaps
fit_running_spline = False
relative_integrate = True
bidirectional_integrate = False
eps_1 = 1e5
eps_2 = 1e5
eps_3 = 2
N_spins_fit = 4
peak_detect = False # For detecting peaks by fitting B_S3 itself instead of 
    #fitting its derivative and computing zero-crossings
zero_crossing_method = 3   
f = 44 * np.pi / 180
detrend_fsp = True
time_correct = True
err_idx = [1020, 1045]