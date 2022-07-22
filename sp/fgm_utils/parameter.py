import datetime as dt
import numpy as np
import os
from typing import Literal


STATE_DATA_DIR = "/nfs/elfin-mirror/elfin/ela/l1/state/defn/2022/"


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
    return os.path.join("/nfs/elfin-mirror/elfin", mission, "l1/state/defn", year_str, f"{mission}_l1_state_defn_{sta_datestr}_v01.cdf")


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


