from . import error 
import pandas as pd
import datetime as dt

def get_betaCSV(csvpath: str):
    """Read start_time from orbits_fgm_cal_{mission}.csv
    """
    try:
        data = pd.read_csv(csvpath)
        data['start_time_datetime'] = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d %H:%M:%S.%f"), data['start_time']))
    except:
        data = []

    return data


def get_beta(df: pd.DataFrame, mid_time: dt.datetime):
    """get the linear interpolation of beta angle at certain time
    """
    mask = (df['start_time_datetime'] <= mid_time)
    before = df[mask].iloc[-1]
    after = df[~mask].iloc[0]

    before_timestamp = before['start_time_datetime'].timestamp()
    after_timestamp = after['start_time_datetime'].timestamp()
    target_timestamp = mid_time.timestamp()

    # calculate slope
    slope = (after['beta_angle'] - before['beta_angle']) / (after_timestamp - before_timestamp)

    # get beta angle
    beta = before['beta_angle'] + slope*(target_timestamp - before_timestamp)

    return beta
