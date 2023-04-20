import datetime as dt
import pandas as pd
import logging.config
from . import fgm_fsp_calib, parameter
from .function import error, preprocess, fgm_fsp_calib_floop
import requests
import numpy as np
from .function.Bplot import Gain_f
from .function.attitude import att_rot

def getCSV(csvpath: str, startdate: str, enddate: str):
    try:
        data = pd.read_csv(csvpath)
        data['Time Start_datetime'] = data['Time Start'].apply(lambda ts: dt.datetime.strptime(ts,'%Y-%m-%d/%H:%M:%S'))
        data['Time End_datetime'] = data['Time End'].apply(lambda ts: dt.datetime.strptime(ts,'%Y-%m-%d/%H:%M:%S'))
        startdate_datetime = dt.datetime.strptime(startdate,'%Y-%m-%d')
        enddate_datetime = dt.datetime.strptime(enddate,'%Y-%m-%d')
        data_select = data[(data['Time Start_datetime'] > startdate_datetime) & (data['Time End_datetime'] < enddate_datetime)]      
    except:
        raise error.SCreadError()
    else:
        if len(data_select) == 0:
            raise error.SCreadError()

    return [data_select['Time Start_datetime'].tolist(), data_select['Time End_datetime'].tolist()] 

if __name__ == "__main__":

    config_file = ('../logging.conf')
    logging.config.fileConfig(config_file)
    logger = logging.getLogger("sp")
    logger.setLevel(logging.INFO)

    mission = "elb"
    csvpath = f"fgm_utils/temp/{mission}_fgm_data_availability.csv"
    elfin_url = "https://data.elfin.ucla.edu/"
    
    """
    startdate = "2022-01-01"
    enddate = "2022-01-31"
    try:
        start_time, end_time = getCSV(csvpath, startdate, enddate)
    except error.SCreadError as e:
        logger.error(e.__str__())
    """
    if parameter.runthree == True:

        # # run ela three long collection together
        # mission = "ela"
        # starttime_str = [
        #     "2019-04-30/18:30:52",
        #     "2019-08-02/02:01:49",
        #     "2019-08-06/07:39:26",] # long collection
        # endtime_str = [
        #     "2019-04-30/20:00:04",
        #     "2019-08-02/02:56:39",
        #     "2019-08-06/08:33:53",]
        
        # run elb three collections, two of them are long
        # mission = "elb"
        # starttime_str = [
        #     "2020-04-30/04:53:05",
        #     "2020-04-30/06:16:29",] # long collection
        # endtime_str = [
        #     "2020-04-30/05:43:40",
        #     "2020-04-30/07:16:44",]
        # f_3 = [
        #     (180+55.0) * np.pi / 180, 
        #     (180+55.1) * np.pi / 180, 
        #     ]
        
        # # run elb three collections, two of them are long
        # mission = "elb"
        # starttime_str = [
        #     "2022-01-01/00:43:01",
        #     "2022-01-01/06:58:16",
        #     "2022-01-01/08:32:03",
        #     "2022-01-01/11:51:34",
        #     "2022-01-01/13:11:28",
        #     "2022-01-01/13:24:31",
        #     "2022-01-01/14:43:59",
        #     "2022-01-01/16:16:16",
        #     "2022-01-01/18:00:45",
        #     "2022-01-01/19:20:21",
        #     "2022-01-01/19:32:53",
        #     "2022-01-01/22:24:46"] # long collection
        # endtime_str = [
        #     "2022-01-01/00:49:11",
        #     "2022-01-01/07:04:29",
        #     "2022-01-01/08:38:17",
        #     "2022-01-01/11:57:47",
        #     "2022-01-01/13:17:42",
        #     "2022-01-01/13:30:45",
        #     "2022-01-01/14:50:13",
        #     "2022-01-01/16:22:28",
        #     "2022-01-01/18:06:57",
        #     "2022-01-01/19:26:33",
        #     "2022-01-01/19:39:05",
        #     "2022-01-01/22:30:58"]
        
        # # run elb three collections, two of them are long
        # mission = "elb"
        # starttime_str = [
        #     "2022-01-03/00:57:59",
        #     "2022-01-03/05:26:33",
        #     "2022-01-03/07:00:31",
        #     "2022-01-03/08:34:14",
        #     "2022-01-03/11:40:46",
        #     "2022-01-03/13:13:35",
        #     "2022-01-03/14:46:05",
        #     "2022-01-03/16:18:21",
        #     "2022-01-03/16:30:53",
        #     "2022-01-03/18:02:50",
        #     "2022-01-03/19:22:22",
        #     "2022-01-03/22:26:51",] # long collection
        # endtime_str = [
        #     "2022-01-03/01:04:09",
        #     "2022-01-03/05:32:48",
        #     "2022-01-03/07:06:43",
        #     "2022-01-03/08:40:26",
        #     "2022-01-03/11:46:58",
        #     "2022-01-03/13:19:49",
        #     "2022-01-03/14:52:19",
        #     "2022-01-03/16:24:35",
        #     "2022-01-03/16:37:07",
        #     "2022-01-03/18:09:04",
        #     "2022-01-03/19:28:36",
        #     "2022-01-03/22:33:03",]

        # mission = "elb"
        # starttime_str = [
        #     "2020-05-15/03:47:54",
        #     "2020-05-15/10:16:29",
        #     ]
        # endtime_str = [
        #     "2020-05-15/04:20:23",
        #     "2020-05-15/10:41:06",
        #     ]
        # f_3 = [
        #     (180+55.0) * np.pi / 180, 
        #     (180+55.1) * np.pi / 180, 
        #     ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-03-05/12:07:10",
        #     "2021-03-05/16:02:23",
        #     "2021-03-05/19:11:10",
        #     "2021-03-05/22:18:52",
        #     "2021-03-05/23:05:38",
        #     ]
        # endtime_str = [
        #     "2021-03-05/12:13:23",
        #     "2021-03-05/16:08:36",
        #     "2021-03-05/19:17:22",
        #     "2021-03-05/22:25:04",
        #     "2021-03-05/23:11:51",
        #     ]
        # f_3 = [
        #     (180+53.60) * np.pi / 180, 
        #     (180+53.55) * np.pi / 180, 
        #     (180+53.54) * np.pi / 180, 
        #     (180+53.55) * np.pi / 180, 
        #     (180+53.60) * np.pi / 180, 
        #     ]

        # mission = "elb"
        # starttime_str = [
        #     "2021-03-06/15:34:11",
        #     "2021-03-06/16:20:56",
        #     "2021-03-06/17:08:05",]
        # endtime_str = [
        #     "2021-03-06/15:40:24",
        #     "2021-03-06/16:27:07",
        #     "2021-03-06/17:14:08",
        #     ]
        # f_3 = [
        #     (180+53.5) * np.pi / 180, 
        #     (180+53.5) * np.pi / 180, 
        #     (180+53.5) * np.pi / 180, 
        #     ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-03-08/19:44:56",
        #     "2021-03-08/22:06:54",
        #     "2021-03-08/22:53:45",
        #     "2021-03-08/23:40:42"]
        # endtime_str = [
        #     "2021-03-08/19:51:10",
        #     "2021-03-08/22:13:06",
        #     "2021-03-08/22:59:59",
        #     "2021-03-08/23:46:45"
        #     ]
        # f_3 = [
        #     (180+53.49) * np.pi / 180, 
        #     (180+53.5) * np.pi / 180, 
        #     (180+53.5) * np.pi / 180, 
        #     (180+53.46) * np.pi / 180, 
        #     ]

        # mission = "elb"
        # starttime_str = [
        #     "2021-03-13/04:01:13",
        #     "2021-03-13/04:47:29",
        #     "2021-03-13/07:52:36",
        #     "2021-03-13/10:11:55",
        #     "2021-03-13/11:45:04",
        #     "2021-03-13/16:25:26",
        #     "2021-03-13/17:59:49",
        #     "2021-03-13/19:34:43",]
        # endtime_str = [
        #     "2021-03-13/04:07:25",
        #     "2021-03-13/04:53:42",
        #     "2021-03-13/07:58:49",
        #     "2021-03-13/10:18:08",
        #     "2021-03-13/11:51:16",
        #     "2021-03-13/16:31:39",
        #     "2021-03-13/18:06:02",
        #     "2021-03-13/19:40:47",
        #     ]
        # f_3 = [
        #     (180+53.6) * np.pi / 180, 
        #     (180+53.7) * np.pi / 180, 
        #     (180+53.55) * np.pi / 180, 
        #     (180+53.60) * np.pi / 180, 
        #     (180+53.61) * np.pi / 180, 
        #     (180+53.59) * np.pi / 180, 
        #     (180+53.62) * np.pi / 180, 
        #     (180+53.59) * np.pi / 180, 
        #     ]

        # mission = "elb"
        # starttime_str = [
        #     "2021-03-16/19:09:32",
        #     "2021-03-16/21:31:37",
        #     "2021-03-16/23:05:30",
        #     "2021-03-16/23:52:05",]
        # endtime_str = [
        #     "2021-03-16/19:15:44",
        #     "2021-03-16/21:37:49",
        #     "2021-03-16/23:11:44",
        #     "2021-03-16/23:58:09",
        #     ]
        # f_3 = [
        #     (180+53.68) * np.pi / 180, 
        #     (180+53.7) * np.pi / 180, 
        #     (180+53.7) * np.pi / 180, 
        #     (180+53.7) * np.pi / 180, 
        #     ]

        # mission = "elb"
        # starttime_str = [
        #     "2021-03-17/18:29:27",
        #     "2021-03-17/19:17:01",
        #     "2021-03-17/20:51:46",
        #     "2021-03-17/22:25:57",]
        # endtime_str = [
        #     "2021-03-17/18:35:40",
        #     "2021-03-17/19:23:13",
        #     "2021-03-17/20:57:58",
        #     "2021-03-17/22:32:00",
        #     ]
        # f_3 = [
        #     (180+53.68) * np.pi / 180, 
        #     (180+53.7) * np.pi / 180, 
        #     (180+53.7) * np.pi / 180, 
        #     (180+53.7) * np.pi / 180, 
        #     ]

        # mission = "elb"
        # starttime_str = [
        #     "2021-03-23/16:04:04",
        #     "2021-03-23/16:51:14",
        #     "2021-03-23/18:25:59",
        #     "2021-03-23/19:13:18",
        #     "2021-03-23/20:00:50",]
        # endtime_str = [
        #     "2021-03-23/16:10:16",
        #     "2021-03-23/16:57:25",
        #     "2021-03-23/18:32:12",
        #     "2021-03-23/19:19:24",
        #     "2021-03-23/20:06:53",
        #     ]
        # f_3 = [
        #     (180+53.68) * np.pi / 180, 
        #     (180+53.70) * np.pi / 180, 
        #     (180+53.70) * np.pi / 180, 
        #     (180+53.68) * np.pi / 180, 
        #     (180+53.66) * np.pi / 180, 
        #     ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-04-14/06:10:41",
        #     "2021-04-14/06:57:02",
        #     "2021-04-14/08:29:58",
        #     "2021-04-14/10:03:02",
        #     "2021-04-14/11:36:13",
        #     ]
        # endtime_str = [
        #     "2021-04-14/06:16:53",
        #     "2021-04-14/07:03:15",
        #     "2021-04-14/08:36:10",
        #     "2021-04-14/10:09:14",
        #     "2021-04-14/11:42:23",
        #     ]
        # f_3 = [
        #        (180+53.68) * np.pi / 180,
        #        (180+53.68) * np.pi / 180,
        #        (180+53.68) * np.pi / 180,
        #        (180+53.68) * np.pi / 180,
        #        (180+53.68) * np.pi / 180,
        #        (180+53.68) * np.pi / 180,]
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-07-27/00:18:26",
        #     "2021-07-27/00:49:56",]
        # endtime_str = [
        #     "2021-07-27/00:24:29",
        #     "2021-07-27/00:56:02",
        #     ]

        # mission = "elb"
        # starttime_str = [
        #     "2021-09-23/08:35:44",
        #     "2021-09-23/09:22:50",]
        # endtime_str = [
        #     "2021-09-23/08:41:56",
        #     "2021-09-23/09:29:03",
        #     ]

        # mission = "elb"
        # starttime_str = [
        #     "2021-10-02/07:39:03",
        #     "2021-10-02/08:26:09",]
        # endtime_str = [
        #     "2021-10-02/07:45:16",
        #     "2021-10-02/08:32:21",
        #     ]

        # mission = "elb"
        # starttime_str = [
        #     "2021-11-15/03:03:43",
        #     "2021-11-15/04:02:46",]
        # endtime_str = [
        #     "2021-11-15/03:09:56",
        #     "2021-11-15/04:08:59",
        #     ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-11-30/15:21:57",
        #     "2021-11-30/16:21:04",]
        # endtime_str = [
        #     "2021-11-30/15:28:10",
        #     "2021-11-30/16:27:07",
        #     ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-12-04/08:42:15",
        #     "2021-12-04/09:17:55",]
        # endtime_str = [
        #     "2021-12-04/08:48:26",
        #     "2021-12-04/09:24:08",
        #     ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-12-04/17:02:45",
        #     "2021-12-04/18:01:25",]
        # endtime_str = [
        #     "2021-12-04/17:08:53",
        #     "2021-12-04/18:07:28",
        #     ]
        
        # mission = "elb"  # can't run
        # starttime_str = [
        #     "2021-12-05/19:22:59",
        #     "2021-12-05/20:21:34",
        #     ]
        # endtime_str = [
        #     "2021-12-05/19:29:02",
        #     "2021-12-05/20:27:37",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-12-08/18:09:10",
        #     "2021-12-08/18:42:41",
        #     ]
        # endtime_str = [
        #     "2021-12-08/18:15:13",
        #     "2021-12-08/18:48:44",
        #     ] 

        # mission = "elb"
        # starttime_str = [
        #     "2021-12-16/08:07:59",
        #     "2021-12-16/09:06:33",
        #     ]
        # endtime_str = [
        #     "2021-12-16/08:14:13",
        #     "2021-12-16/09:12:47",
        #     ] 

        # mission = "elb"
        # starttime_str = [
        #     "2021-12-16/15:20:11",
        #     "2021-12-16/15:53:24",
        #     ]
        # endtime_str = [
        #     "2021-12-16/15:26:14",
        #     "2021-12-16/15:59:35",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-12-18/12:51:45",
        #     "2021-12-18/13:04:46",
        #     ]
        # endtime_str = [
        #     "2021-12-18/12:57:57",
        #     "2021-12-18/13:11:00",
        #     ] 
        
        
        # mission = "elb"
        # starttime_str = [
        #     "2022-03-08/15:55:41",
        #     "2022-03-08/16:08:02",
        #     ]
        # endtime_str = [
        #     "2022-03-08/16:01:51",
        #     "2022-03-08/16:14:14",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2022-03-09/16:39:30",
        #     "2022-03-09/16:51:55",
        #     ]
        # endtime_str = [
        #     "2022-03-09/16:45:43",
        #     "2022-03-09/16:58:09",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2022-03-14/21:49:47",
        #     "2022-03-14/22:02:31",
        #     ]
        # endtime_str = [
        #     "2022-03-14/21:55:59",
        #     "2022-03-14/22:08:42",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2022-03-23/14:27:39",
        #     "2022-03-23/14:40:00",
        #     ]
        # endtime_str = [
        #     "2022-03-23/14:33:50",
        #     "2022-03-23/14:46:13",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2022-07-12/17:42:02",
        #     "2022-07-12/17:54:32",
        #     ]
        # endtime_str = [
        #     "2022-07-12/17:48:47",
        #     "2022-07-12/18:01:16",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-04-04/04:16:38",
        #     "2021-04-04/05:02:49",
        #     ]
        # endtime_str = [
        #     "2021-04-04/04:22:50",
        #     "2021-04-04/05:08:59",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-05-02/00:13:48",
        #     "2021-05-02/00:59:58",
        #     ]
        # endtime_str = [
        #     "2021-05-02/00:20:00",
        #     "2021-05-02/01:06:10",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-05-05/16:05:33",
        #     "2021-05-05/16:52:45",
        #     ]
        # endtime_str = [
        #     "2021-05-05/16:11:44",
        #     "2021-05-05/16:58:48",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-06-10/03:26:43",
        #     "2021-06-10/04:12:52",
        #     ]
        # endtime_str = [
        #     "2021-06-10/03:32:56",
        #     "2021-06-10/04:19:03",
        #     ] 
        
        # mission = "elb"
        # starttime_str = [
        #     "2021-07-15/12:25:02",
        #     "2021-07-15/13:12:12",
        #     ]
        # endtime_str = [
        #     "2021-07-15/12:31:14",
        #     "2021-07-15/13:18:26",
        #     ] 

        # mission = "elb"
        # starttime_str = [
        #     "2022-01-01/13:11:28",
        #     "2022-01-01/13:24:31",
        #     "2022-01-01/14:43:59",
        #     "2022-01-01/16:16:16",
        #     "2022-01-01/18:00:45",
        #     "2022-01-01/19:20:21",
        #     "2022-01-01/19:32:53",
        #     "2022-01-01/22:24:46",
        #     ]
        # endtime_str = [
        #     "2022-01-01/13:17:42",
        #     "2022-01-01/13:30:45",
        #     "2022-01-01/14:50:13",
        #     "2022-01-01/16:22:28",
        #     "2022-01-01/18:06:57",
        #     "2022-01-01/19:26:33",
        #     "2022-01-01/19:39:05",
        #     "2022-01-01/22:30:58"
        #     ] 
        # f_3 = [
        #     (180+54.3) * np.pi / 180, 
        #     (180+54.0) * np.pi / 180, 
        #     (180+54.3) * np.pi / 180, 
        #     (180+54.3) * np.pi / 180, 
        #     (180+54.16) * np.pi / 180,
        #     (180+54.3) * np.pi / 180,
        #     (180+54.18) * np.pi / 180,
        #     (180+54.35) * np.pi / 180] # this has to match the number of starttime_str

        mission = "elb"
        starttime_str = [
            "2022-01-03/00:57:59",
            "2022-01-03/05:26:33",
            "2022-01-03/07:00:31",
            "2022-01-03/08:34:14",
            "2022-01-03/11:40:46",
            ]
        endtime_str = [
            "2022-01-03/01:04:09",
            "2022-01-03/05:32:48",
            "2022-01-03/07:06:43",
            "2022-01-03/08:40:26",
            "2022-01-03/11:46:58",
            ] 
        f_3 = [
            (180+54.3) * np.pi / 180, 
            (180+54.3) * np.pi / 180, 
            (180+54.3) * np.pi / 180, 
            (180+54.3) * np.pi / 180, 
            (180+54.3) * np.pi / 180,
            ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2022-03-09/15:07:43",
        #     "2022-03-09/16:39:30",
        #     "2022-03-09/16:51:55",
        #     "2022-03-09/19:43:29",
        #     "2022-03-09/22:48:59",
        #     ]
        # endtime_str = [
        #     "2022-03-09/15:13:56",
        #     "2022-03-09/16:45:43",
        #     "2022-03-09/16:58:09",
        #     "2022-03-09/19:49:41",
        #     "2022-03-09/22:55:11",
        #     ] 
        # f_3 = [
        #     (180+54.15) * np.pi / 180, 
        #     (180+54.14) * np.pi / 180, 
        #     (180+54.24) * np.pi / 180, 
        #     (180+54.15) * np.pi / 180, 
        #     (180+54.17) * np.pi / 180,
        #     ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2022-04-01/11:35:51",
        #     "2022-04-01/13:20:13",
        #     "2022-04-01/14:51:49",
        #     "2022-04-01/17:42:54",
        #     "2022-04-01/19:14:56",
        #     "2022-04-01/21:00:18",
        #     ]
        # endtime_str = [
        #     "2022-04-01/11:42:05",
        #     "2022-04-01/13:26:25",
        #     "2022-04-01/14:58:03",
        #     "2022-04-01/17:49:08",
        #     "2022-04-01/19:21:08",
        #     "2022-04-01/21:06:30",
        #     ] 
        # f_3 = [
        #     (180+54.08) * np.pi / 180, 
        #     (180+54.18) * np.pi / 180, 
        #     (180+54.16) * np.pi / 180, 
        #     (180+54.03) * np.pi / 180, 
        #     (180+54.03) * np.pi / 180,
        #     (180+54.18) * np.pi / 180,
        #     ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2022-06-14/01:00:03",
        #     "2022-06-14/07:09:09",
        #     "2022-06-14/14:45:14",
        #     "2022-06-14/22:26:14",
        #     ]
        # endtime_str = [
        #     "2022-06-14/01:06:15",
        #     "2022-06-14/07:15:21",
        #     "2022-06-14/14:51:26",
        #     "2022-06-14/22:32:26",
        #     ] 
        # f_3 = [
        #     (180+54.87) * np.pi / 180, 
        #     (180+54.89) * np.pi / 180, 
        #     (180+54.87) * np.pi / 180, 
        #     (180+54.87) * np.pi / 180, 
        #     ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2022-08-17/00:55:03",
        #     "2022-08-17/02:26:31",
        #     "2022-08-17/11:30:01",
        #     "2022-08-17/23:41:54",
        #     ]
        # endtime_str = [
        #     "2022-08-17/01:01:15",
        #     "2022-08-17/02:32:43",
        #     "2022-08-17/11:36:14",
        #     "2022-08-17/23:48:06",
        #     ] 
        # f_3 = [
        #     (180+54.72) * np.pi / 180, 
        #     (180+54.73) * np.pi / 180, 
        #     (180+54.73) * np.pi / 180, 
        #     (180+54.72) * np.pi / 180, 
        #     ]
        
        # mission = "elb"
        # starttime_str = [
        #     "2022-09-08/04:19:28",
        #     "2022-09-08/16:19:39",
        #     "2022-09-08/19:21:01",
        #     "2022-09-08/22:24:09",
        #     ]
        # endtime_str = [
        #     "2022-09-08/04:25:41",
        #     "2022-09-08/16:25:51",
        #     "2022-09-08/19:27:14",
        #     "2022-09-08/22:30:22",
        #     ] 
        # f_3 = [
        #     (180+54.61) * np.pi / 180, 
        #     (180+54.73) * np.pi / 180, 
        #     (180+54.73) * np.pi / 180, 
        #     (180+54.72) * np.pi / 180, 
        #     ]

        start_time = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d/%H:%M:%S"), starttime_str))
        end_time = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d/%H:%M:%S"), endtime_str))

        fgm_cdfdata_3 = pd.DataFrame()
        att_cdfdata_3 = pd.DataFrame()
        pos_cdfdata_3 = pd.DataFrame()
        for i in range(len(start_time)):

            sta_datestr = start_time[i].strftime("%Y%m%d")
            logger.info(f"▶️ Received {mission} collection from {start_time[i]} to {end_time[i]}")
            sta_cdfpath = f"fgm_utils/test/{mission}_l1_state_defn_{sta_datestr}_v02.cdf"
            fgm_cdfpath = f"fgm_utils/test/{mission}_l1_fgs_{sta_datestr}_v01.cdf" 

            if parameter.download_data == True:
                try:
                    sta_url = f"{elfin_url}{mission}/l1/state/defn/{start_time[i].year}/{mission}_l1_state_defn_{sta_datestr}_v02.cdf"
                    res = requests.get(sta_url)
                    open(sta_cdfpath, 'wb').write(res.content)
                except:
                    print(f"download error:{sta_url}")

                try:
                    fgm_url = f"{elfin_url}{mission}/l1/fgm/survey/{start_time[i].year}/{mission}_l1_fgs_{sta_datestr}_v01.cdf"
                    res = requests.get(fgm_url)
                    open(fgm_cdfpath, 'wb').write(res.content)
                except:
                    print(f"download error:{fgm_url}") 
                    breakpoint() 

            fgm_cdfdata = pd.DataFrame(preprocess.get_cdf(fgm_cdfpath, vars=[f"{mission}_fgs_time", f"{mission}_fgs"]))
            logger.info(f"Sucessfully read cdf for {mission} from {start_time[i]} to {end_time[i]}")
            att_cdfdata, pos_cdfdata = preprocess.get_relevant_state_data(sta_cdfpath, mission, start_time[i], end_time[i])
            logger.info(f"Sucessfully read state cdf for {mission} from {start_time[i]} to {end_time[i]}")      

            if parameter.att_rot == True:
                att_cdfdata = att_rot(att_cdfdata, parameter.att_rot_ang, parameter.att_rot_axis)    

            fgm_cdfdata_3 = pd.concat([fgm_cdfdata_3, fgm_cdfdata])
            att_cdfdata_3 = pd.concat([att_cdfdata_3, att_cdfdata])
            pos_cdfdata_3 = pd.concat([pos_cdfdata_3, pos_cdfdata])
                    
            [
                ctime_0, ctimestamp_0,
                fgs_ful_fgm_0th_x_0, fgs_ful_fgm_0th_y_0, fgs_ful_fgm_0th_z_0, 
                fgs_igrf_gei_x_0, fgs_igrf_gei_y_0, fgs_igrf_gei_z_0,
                att_gei_x_0, att_gei_y_0, att_gei_z_0,
                pos_gei_x_0, pos_gei_y_0, pos_gei_z_0] = fgm_fsp_calib_floop.fgm_fsp_calib_prepos(
                    mission, start_time[i], end_time[i], fgm_cdfdata, att_cdfdata, pos_cdfdata, logger)
            f_3_0 = [f_3[i]] * len(fgs_ful_fgm_0th_x_0)

            if i == 0:
                [
                    ctime, ctimestamp, 
                    fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                    att_gei_x, att_gei_y, att_gei_z,
                    pos_gei_x, pos_gei_y, pos_gei_z, f_3] = [
                        ctime_0, ctimestamp_0, 
                        fgs_ful_fgm_0th_x_0, fgs_ful_fgm_0th_y_0, fgs_ful_fgm_0th_z_0, 
                        fgs_igrf_gei_x_0, fgs_igrf_gei_y_0, fgs_igrf_gei_z_0, 
                        att_gei_x_0, att_gei_y_0, att_gei_z_0, 
                        pos_gei_x_0, pos_gei_y_0, pos_gei_z_0, f_3_0]

            else:
                ctime = np.concatenate((ctime, ctime_0 + ctime[-1] + 30))
                fgs_ful_fgm_0th_x = np.concatenate((fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_x_0))
                fgs_ful_fgm_0th_y = np.concatenate((fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_y_0))
                fgs_ful_fgm_0th_z = np.concatenate((fgs_ful_fgm_0th_z, fgs_ful_fgm_0th_z_0))
                fgs_igrf_gei_x = np.concatenate((fgs_igrf_gei_x, fgs_igrf_gei_x_0))
                fgs_igrf_gei_y = np.concatenate((fgs_igrf_gei_y, fgs_igrf_gei_y_0))
                fgs_igrf_gei_z = np.concatenate((fgs_igrf_gei_z, fgs_igrf_gei_z_0))
                att_gei_x = np.concatenate((att_gei_x, att_gei_x_0))
                att_gei_y = np.concatenate((att_gei_y, att_gei_y_0))
                att_gei_z = np.concatenate((att_gei_z, att_gei_z_0))
                pos_gei_x = np.concatenate((pos_gei_x, pos_gei_x_0))
                pos_gei_y = np.concatenate((pos_gei_y, pos_gei_y_0))
                pos_gei_z = np.concatenate((pos_gei_z, pos_gei_z_0))
                f_3 = np.concatenate([f_3, f_3_0])


        [
            FGM_timestamp, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z, B_parameter]=fgm_fsp_calib_floop.fgm_fsp_calib_all3(
            ctime, ctimestamp, f_3,
            fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z,
            att_gei_x, att_gei_y, att_gei_z,
            pos_gei_x, pos_gei_y, pos_gei_z,
            logger
        )

        # define a dataframe for restoring f loop result
        [G11, G12, G13, O1, G21, G22, G23, O2, G31, G32, G33, O3] = [*B_parameter]   
        G1 = (G11**2 + G12**2 + G13**2)**0.5 
        G2 = (G21**2 + G22**2 + G23**2)**0.5
        G3 = (G31**2 + G32**2 + G33**2)**0.5
        
        th1 = np.degrees(np.arccos(G13/G1))
        th2 = np.degrees(np.arccos(G23/G2))
        th3 = np.degrees(np.arccos(G33/G3))

        ph1 = np.degrees(np.arctan(G12/G11))
        ph2 = np.degrees(np.arctan(G22/G21))
        ph3 = np.degrees(np.arctan(G32/G31))

        print(f'G1: {G1} G2:{G2} G3:{G3}\n')
        print(f'th1: {th1} th2:{th2} th3:{th3}\n')
        print(f'ph1: {ph1} ph2:{ph2} ph3:{ph3}\n')
        print(f'O1/G1: {O1/G1} O2/G2:{O2/G2} O3/G3:{O3/G3}\n')

        columns = [ 
            'G1', 'G2', 'G3',
            'th1','th2','th3',
            'ph1','ph2','ph3',
            'O1/G1','O2/G2','O3/G3']
        floop_df = pd.DataFrame(columns=columns)
        new_row = [G1, G2, G3, th1, th2, th3, ph1, ph2, ph3, O1/G1, O2/G2, O3/G3]
        new_row = [dict(zip(columns, new_row))]
        floop_df = pd.concat([floop_df] + [pd.DataFrame(new_row)])
        filename = f"fgm_utils/rotation_angle/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}_run3_{mission}.csv"
        floop_df.to_csv(filename, index=False)

        
    else:
        #starttime_str = ["2019-04-30/18:30:52"] # long collection
        #endtime_str = ["2019-04-30/20:00:04"]
        #starttime_str = ["2019-08-02/02:01:49"] # long collection
        #endtime_str = ["2019-08-02/02:56:39"]
        #starttime_str = ["2019-08-06/07:39:26"] # long collection
        #endtime_str = ["2019-08-06/08:33:53"]
        #starttime_str = ["2022-03-30/19:53:44"] # ela long collection in 2022
        #endtime_str = ["2022-03-30/20:12:53"]
        starttime_str = ["2022-09-10/23:42:55"] # jiang event ela
        endtime_str = ["2022-09-10/23:49:06"] 
        #starttime_str = ["2022-03-13/19:33:28"]
        #endtime_str = ["2022-03-13/19:39:41"]
        #starttime_str = ["2022-03-31/16:11:16"] # a good event
        #endtime_str = ["2022-03-31/16:17:27"] 
        #starttime_str = ["2022-04-01/13:48:32"] # a good event
        #endtime_str = ["2022-04-01/13:54:45"] 
        #starttime_str = ["2022-04-01/15:20:07"] # a good event
        #endtime_str = ["2022-04-01/15:26:20"] 
        #starttime_str = ["2021-03-23/03:39:52"] # another event with o larger than w
        #endtime_str = ["2021-03-23/03:46:04"] 
        #starttime_str = ["2021-03-23/06:44:45"] # another event with o larger than w
        #endtime_str = ["2021-03-23/06:50:57"] 

        #starttime_str = ["2021-03-26/03:14:37"] # same day, different orbit, bad fgm
        #endtime_str = ["2021-03-26/03:20:51"] 
        #starttime_str = ["2021-03-26/02:28:33"] # same day, different orbit, good fgm
        #endtime_str = ["2021-03-26/02:34:44"] 
        #starttime_str = ["2022-01-14/18:49:45"]  # same day for ela and elb, opposite fgm ela
        #endtime_str = ["2022-01-14/18:55:58"] 
        #starttime_str = ["2022-01-14/18:46:04"] # same day for ela and elb, opposite fgm elb
        #endtime_str = ["2022-01-14/18:52:18"] 

        #starttime_str = ["2019-04-30/18:30:52"] # long collection 1h 30m
        #endtime_str = ["2019-04-30/20:00:04"] # have to make funkyfgm = false
        #starttime_str = ["2019-08-02/02:01:49"] # long collection 55m 
        #endtime_str = ["2019-08-02/02:56:39"] 
        #starttime_str = ["2019-08-06/07:39:26"] # long collection 55m
        #endtime_str = ["2019-08-06/08:33:53"] 
        #starttime_str = ["2020-06-11/13:16:55"] # elb long collection 1
        #endtime_str = ["2020-06-11/14:01:10"]
        #starttime_str = ["2022-01-14/18:46:04"] # elb opposite result
        #endtime_str = ["2022-01-14/18:52:18"]
        #starttime_str = ["2022-01-14/18:49:45"] # ela opposite result
        #endtime_str = ["2022-01-14/18:55:58"]
        #starttime_str = ["2022-01-14/20:18:08"] # elb similar result
        #endtime_str = ["2022-01-14/20:24:21"]
        #starttime_str = ["2022-01-14/20:21:46"] # elb similar result
        #endtime_str = ["2022-01-14/20:27:58"]
        #starttime_str = ["2022-04-29/07:13:57"] # elb opposite result
        #endtime_str = ["2022-04-29/07:20:10"]
        #starttime_str = ["2022-04-29/07:09:44"] # ela opposite result
        #endtime_str = ["2022-04-29/07:16:13"]
        #starttime_str = ["2021-03-03/07:13:37"] # elb good result in 2021
        #endtime_str = ["2021-03-03/07:19:51"]
        #starttime_str = ["2021-03-05/03:37:08"] # elb bad result in 2021
        #endtime_str = ["2021-03-05/03:43:20"]
        #starttime_str = ["2021-07-14/18:36:52"] # elb bad result in 2021 july
        #endtime_str = ["2021-07-14/18:43:04"]
        #starttime_str = ["2021-10-18/19:45:25"] # elb good result in 2021 oct
        #endtime_str = ["2021-10-18/19:51:38"]
        #starttime_str = ["2021-11-16/06:25:05"] # elb bad result in 2021 nov
        #endtime_str = ["2021-11-16/06:31:19"]
        #starttime_str = ["2021-12-18/06:37:45"] # elb result in 2021 dec
        #endtime_str = ["2021-12-18/06:43:59"]
        #starttime_str = ["2021-12-31/07:43:59"] # elb result in 2021 dec
        #endtime_str = ["2021-12-31/07:50:11"]
        #starttime_str = ["2022-01-08/10:59:14"] # elb result in 2021 dec
        #endtime_str = ["2022-01-08/11:05:27"]
        #starttime_str = ["2022-01-11/16:25:46"] # elb result in 2021 dec
        #endtime_str = ["2022-01-11/16:31:59"]
        #starttime_str = ["2022-01-14/15:42:02"] # elb result in 2021 dec
        #endtime_str = ["2022-01-14/15:48:16"]
        #starttime_str = ["2022-01-12/15:40:28"] # elb result in 2021 dec
        #endtime_str = ["2022-01-12/15:46:42"]
        #starttime_str = ["2022-01-13/17:59:20"] # elb result in 2021 dec
        #endtime_str = ["2022-01-13/18:05:33"]
        #starttime_str = ["2022-08-02/06:26:50"] # elb result in 2021 dec
        #endtime_str = ["2022-08-02/06:33:03"]
        #starttime_str = ["2022-02-15/03:06:56"] # ela result in 2021 dec
        #endtime_str = ["2022-02-15/03:13:09"]
        #starttime_str = ["2022-03-16/10:12:07"] # ela result in 2021 dec
        #endtime_str = ["2022-03-16/10:18:18"]
        #starttime_str = ["2022-04-02/08:08:39"] # ela result in 2021 dec
        #endtime_str = ["2022-04-02/08:14:53"]
        #starttime_str = ["2022-03-23/08:56:57"] # ela result in 2021 dec
        #endtime_str = ["2022-03-23/09:03:10"]
        #starttime_str = ["2022-03-29/08:29:15"] # ela result in 2021 dec
        #endtime_str = ["2022-03-29/08:35:28"]
        #starttime_str = ["2020-08-19/05:43:11"] # ela result in 2021 dec
        #endtime_str = ["2020-08-19/05:49:25"]
        #starttime_str = ["2021-08-13/21:21:04"] # ela result in 2021 dec
        #endtime_str = ["2021-08-13/21:27:17"]
        #starttime_str = ["2022-06-07/11:27:16"] # ela result in 2021 dec
        #endtime_str = ["2022-06-07/11:33:27"]

        start_time = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d/%H:%M:%S"), starttime_str))
        end_time = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d/%H:%M:%S"), endtime_str))

        for i in range(len(start_time)):

            sta_datestr = start_time[i].strftime("%Y%m%d")
            logger.info(f"▶️ Received {mission} collection from {start_time[i]} to {end_time[i]}")
            sta_cdfpath = f"fgm_utils/test/{mission}_l1_state_defn_{sta_datestr}_v02.cdf"
            fgm_cdfpath = f"fgm_utils/test/{mission}_l1_fgs_{sta_datestr}_v01.cdf" 

            if parameter.download_data == True:
                try:
                    sta_url = f"{elfin_url}{mission}/l1/state/defn/{start_time[i].year}/{mission}_l1_state_defn_{sta_datestr}_v02.cdf"
                    res = requests.get(sta_url)
                    open(sta_cdfpath, 'wb').write(res.content)
                except:
                    print(f"download error:{sta_url}")

                try:
                    fgm_url = f"{elfin_url}{mission}/l1/fgm/survey/{start_time[i].year}/{mission}_l1_fgs_{sta_datestr}_v01.cdf"
                    res = requests.get(fgm_url)
                    open(fgm_cdfpath, 'wb').write(res.content)
                except:
                    print(f"download error:{fgm_url}") 
                    breakpoint()   

            try: 
                fgm_cdfdata = pd.DataFrame(preprocess.get_cdf(fgm_cdfpath, vars=[f"{mission}_fgs_time", f"{mission}_fgs"]))
                logger.info(f"Sucessfully read cdf for {mission} from {start_time[i]} to {end_time[i]}")
                att_cdfdata, pos_cdfdata = preprocess.get_relevant_state_data(sta_cdfpath, mission, start_time[i], end_time[i])
                logger.info(f"Sucessfully read state cdf for {mission} from {start_time[i]} to {end_time[i]}")
            except error.cdfError as e:
                logger.error(e.__str__())
            else:
                if parameter.f_loop == False:
                    [
                        FGM_timestamp,
                        fgs_fsp_res_dmxl_x,
                        fgs_fsp_res_dmxl_y,
                        fgs_fsp_res_dmxl_z,
                        fgs_fsp_igrf_dmxl_x,
                        fgs_fsp_igrf_dmxl_y,
                        fgs_fsp_igrf_dmxl_z,
                        fgs_fsp_res_dmxl_trend_x,
                        fgs_fsp_res_dmxl_trend_y,
                        fgs_fsp_res_dmxl_trend_z,
                        fgs_fsp_res_gei_x,
                        fgs_fsp_res_gei_y,
                        fgs_fsp_res_gei_z,
                        fgs_fsp_igrf_gei_x,
                        fgs_fsp_igrf_gei_y,
                        fgs_fsp_igrf_gei_z,
                    ] = fgm_fsp_calib(mission, start_time[i], end_time[i], fgm_cdfdata, att_cdfdata, pos_cdfdata, logger)
                    logger.info(f"⏹️ End of fsp calibration for {mission} from {start_time[i]} to {end_time[i]}\n")
                else:
                    [
                        ctime, ctimestamp,
                        fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
                        fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z,
                        att_gei_x, att_gei_y, att_gei_z,
                        pos_gei_x, pos_gei_y, pos_gei_z] = fgm_fsp_calib_floop.fgm_fsp_calib_prepos(
                            mission, start_time[i], end_time[i], fgm_cdfdata, att_cdfdata, pos_cdfdata, logger)

                    # define a dataframe for restoring f loop result
                    columns = ['rotate_ang','res_dmxl_x', 'res_dmxl_y', 'res_dmxl_z', 
                        'G11', 'G12', 'G13','O1',
                        'G21','G22','G23','O2',
                        'G31','G32','G33','O3']
                    floop_df = pd.DataFrame(columns=columns)
                    for floop_i in parameter.f_loop_value:
                        [
                            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, B_parameter] = fgm_fsp_calib_floop.fgm_fsp_calib_floop(
                                ctime, ctimestamp,
                                fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z,
                                att_gei_x, att_gei_y, att_gei_z,
                                pos_gei_x, pos_gei_y, pos_gei_z,
                                logger, floop_i)
                        logger.info(f"⏹️ End of fsp calibration for {mission} from {start_time[i]} to {end_time[i]} with rotation angle {floop_i}\n")

                        """f loop output: get res std and B parameter to big 
                        """     
                        new_row = [floop_i*180/np.pi, np.std(fgs_fsp_res_dmxl_x), np.std(fgs_fsp_res_dmxl_y), np.std(fgs_fsp_res_dmxl_z), *B_parameter]
                        new_row = [dict(zip(columns, new_row))]
                        floop_df = pd.concat([floop_df] + [pd.DataFrame(new_row)])

                    filename = f"fgm_utils/rotation_angle/floop_{mission}.csv"
                    floop_df.to_csv(filename, index=False)
                    breakpoint()
                    #Gain_f(floop_df['rotate_ang'], floop_df['G11'], floop_df['G22'], floop_df['G33'])
