import numpy as np

proper_pad = True  # fails when data have gaps
fit_running_spline = True
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
init_secs = 0 # seconds

detrend_fsp = True
detrend_cutoff = 6 # criterion of linear detrend with first and second points

ctime_correct = False # degap unipolar gap if true
ctime_thrhld = 5e-5 # criterion for finding ctime gaps
                    #difference between jumps and ctime median in seconds

bad_data_correct = False
bad_data = 794 # example event 0625

makeplot = False
savepng = True


