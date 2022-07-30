import datetime as dt
import numpy as np
import os
from typing import Literal


STATE_DATA_DIR = "/nfs/elfin-mirror/elfin/ela/l1/state/defn/2022/"


def get_state_data_dir(mission_id: Literal[1, 2], data_datetime: dt.datetime) -> str:
    mission_str = "ela" if mission_id == 1 else "elb"
    year_str = str(data_datetime.year)
    sta_datestr = data_datetime.strftime("%Y%m%d")
    return os.path.join("/nfs/elfin-mirror/elfin", mission_str, "l1/state/defn", year_str, f"ela_l1_state_defn_{sta_datestr}_v01.cdf")


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
