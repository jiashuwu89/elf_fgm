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
