from cdflib import CDF, cdfepoch
from typing import List, Literal, Tuple, Union
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import datetime
import os.path
from . import error 
from .. import parameter
from . import Bplot


def get_cdf(cdfpath: str, vars: Union[List[str], None]):
    """Read CDF
    """
    try:
        cdf = CDF(cdfpath)
        cdfinfo = cdf.cdf_info()
        data = {}
        if vars is None:
            vars = cdfinfo["zVariables"]
            print(f"{cdfpath} variables: {vars}")
        for var in vars:
            val = cdf.varget(var)
            if var.endswith("_time"):
                data[var] = list(
                    map(lambda t: cdfepoch.to_datetime(t)[0], val.tolist())
                )
            elif isinstance(val, np.ndarray):
                data[var] = val.tolist()
            else:
                data[var] = val

        return data
    except:
        if os.path.exists(cdfpath) == False:
            raise error.cdfError(cdfpath, "CDF not found")
        else:
            raise error.cdfError(cdfpath, "CDF not read")


def clip_cdfdata(
    df: pd.Series, starttime: pd.Timestamp, endtime: pd.Timestamp
) -> pd.Series:
    """clip data
    """
    startindex = df.index.get_indexer([starttime], method="nearest")[0]
    endindex = df.index.get_indexer([endtime], method="nearest")[0]
    if endindex == len(df) - 1 :
        endindex = len(df)

    return df[startindex:endindex]


def resample_data(
    cur_time: pd.Timestamp, cur_data: pd.Series, target_time: pd.Timestamp
) -> pd.Series:
    """resample data
    """
    cur_data_np = np.array(cur_data.to_list())
    if len(cur_data_np.shape) == 1:
        interp_data = np.zeros(len(target_time))
        x = (cur_time - target_time[0]).total_seconds()
        x_interp = (target_time - target_time[0]).total_seconds()
        f = interp1d(x, cur_data_np, fill_value=np.median(cur_data), bounds_error=False) 
        interp_data = f(x_interp)
    else:    
        dimens = cur_data_np.shape[1]  # x, y, z
        interp_data = np.zeros((len(target_time), dimens))
        x = (cur_time - target_time[0]).total_seconds()
        x_interp = (target_time - target_time[0]).total_seconds()
        for dimen in range(dimens):
            f = interp1d(x, cur_data_np[:, dimen])
            interp_data[:, dimen] = f(x_interp)

    return pd.Series(interp_data.tolist())


def get_relevant_state_data(sta_cdfpath: str, mission: Literal["ela", "elb"], 
    starttime: datetime.datetime, endtime: datetime.datetime) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # read state cdf for att
    att_cdfdata = pd.DataFrame(
        get_cdf(sta_cdfpath, vars=[f"{mission}_att_time", f"{mission}_att_gei"])
    )
    att_cdfdata.set_index(f"{mission}_att_time", inplace=True)
    att_cdfdata = clip_cdfdata(
        att_cdfdata,
        starttime - datetime.timedelta(minutes=2),
        endtime + datetime.timedelta(minutes=2),
    )

    # read state cdf for pos; not read together with att b/c different length
    pos_cdfdata = pd.DataFrame(
        get_cdf(sta_cdfpath, vars=[f"{mission}_pos_gei"])
    )  # not read state time b/c too slow
    pos_cdfdata[f"{mission}_pos_time"] = pd.date_range(
        start=starttime.date(), periods=len(pos_cdfdata[f"{mission}_pos_gei"]), freq="S"
    )
    pos_cdfdata.set_index(f"{mission}_pos_time", inplace=True)
    pos_cdfdata = clip_cdfdata(
        pos_cdfdata,
        starttime - datetime.timedelta(minutes=2),
        endtime + datetime.timedelta(minutes=2),
    )
    return att_cdfdata, pos_cdfdata


def funkyfgm_check(B_x, ctime, datestr):

    if len(B_x) < 3:
        raise error.funkyFGMError_len()
    cross_times = []
    dB_x = np.gradient(B_x) / np.gradient(ctime)
    for i in range(1, len(ctime) - 2):
        if (
            dB_x[i] < 0
            and dB_x[i+1] > 0
        ):
            # jwu: when gap exits, dB can jump from positive to negative
            y1 = dB_x[i]
            y2 = dB_x[i + 1]
            x1 = ctime[i]
            x2 = ctime[i + 1]
            cross_times.append((y2 * x1 - y1 * x2) / (y2 - y1))
              
    #if parameter.makeplot == True:
    #    Bplot.B_ctime_plot_single(ctime, dB_x, cross_times=cross_times, title = "funckyfgm", datestr = datestr)

    if len(cross_times) < 3 :
        raise error.CrossTime1Error(0)
    cross_times = np.array(cross_times)
    cross_times_diff = cross_times[1:-1] - cross_times[0:-2]
    idx = cross_times_diff < 4
    med = np.median(cross_times_diff[idx])
    std = np.std(cross_times_diff[idx])
    if med < 2.5 or med > 3.2 or std > parameter.Spinrate_thrhld * med:
        raise error.funkyFGMError(med, std)
    return


def ctime_check(ctime):
    """some collections have completeness over 100%, example 
        2022-03-01/02:25:52	2022-03-01/02:32:03
        ctime repeated with very small difference, delete repeated ones
    """
    delta_t = np.median(ctime[1:]-ctime[:-1])
    ctime_adj = ctime[1:]-ctime[:-1] - delta_t
    ctime_idx_repeat = []
    i = 0
    while i < len(ctime_adj):
        # if two adjacent ctime_idx sum up to one, delete the repeat one
        if i < len(ctime_adj)-1 and np.abs(np.abs(ctime_adj[i] + ctime_adj[i+1]) - 0.1) < 0.001:
            ctime_idx_repeat.append(i+1)
            i += 2
        else:
            i += 1
    return ctime_idx_repeat
