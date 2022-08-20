import datetime as dt
import pandas as pd
import logging.config
from . import fgm_fsp_calib
from .function import error, preprocess

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
    mission = "ela"
    csvpath = f"fgm_utils/temp/{mission}_fgm_data_availability.csv"
    
    """
    startdate = "2022-01-14"
    enddate = "2022-01-15"
    try:
        start_time, end_time = getCSV(csvpath, startdate, enddate)
    except error.SCreadError as e:
        logger.error(e.__str__())
    """

    #starttime_str = ["2022-01-14/15:45:50"]
    #endtime_str = ["2022-01-14/15:52:04"]
    #starttime_str = ["2022-06-23/04:00:07"]
    #endtime_str = ["2022-06-23/04:06:19"]
    #starttime_str = ["2022-01-14/17:17:52"]
    #endtime_str = ["2022-01-14/17:24:03"]

    starttime_str = ["2022-06-23/04:00:07"]
    endtime_str = ["2022-06-23/04:06:19"]
    #starttime_str = ["2022-01-14/23:26:43"]
    #endtime_str = ["2022-01-14/23:32:54"]
    #starttime_str = "2022-01-12 15:45:51"
    #endtime_str = "2022-01-12 15:52:04"
    #starttime_str = "2022-01-14 15:45:50"
    #endtime_str = "2022-01-14 15:52:04"
    #starttime_str = "2022-01-12 17:17:52"
    #endtime_str = "2022-01-12 17:24:05"
    #starttime_str = "2022-01-12 19:02:20"
    #endtime_str = "2022-01-12 19:08:31"
    #starttime_str = "2022-01-12 20:21:48"
    #endtime_str = "2022-01-12 20:27:57"
    start_time = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d/%H:%M:%S"), starttime_str))
    end_time = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d/%H:%M:%S"), endtime_str))

    for i in range(len(start_time)):

        sta_datestr = start_time[i].strftime("%Y%m%d")
        sta_cdfpath = f"fgm_utils/test/{mission}_l1_state_defn_{sta_datestr}_v02.cdf"
        fgm_cdfpath = f"fgm_utils/test/{mission}_l1_fgs_{sta_datestr}_v01.cdf"  
        logger.info(f"Received {mission} collection from {start_time[i]} to {end_time[i]}")

        

        try: 
            fgm_cdfdata = pd.DataFrame(preprocess.get_cdf(fgm_cdfpath, vars=[f"{mission}_fgs_time", f"{mission}_fgs"]))
            logger.info(f"Sucessfully read cdf for {mission} from {start_time[i]} to {end_time[i]}")
            att_cdfdata, pos_cdfdata = preprocess.get_relevant_state_data(sta_cdfpath, mission, start_time[i], end_time[i])
            logger.info(f"Sucessfully read state cdf for {mission} from {start_time[i]} to {end_time[i]}")
        except error.cdfError as e:
            logger.error(e.__str__())
        else:
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
            logger.info(f"End of fsp calibration for {mission} from {start_time[i]} to {end_time[i]}\n")


    #starttime_str = "2022-01-12 15:45:51"
    #endtime_str = "2022-01-12 15:52:04"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220112_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220112_v01.cdf"

    # SPIKE
    #starttime_str = "2022-01-14 15:45:50"
    #endtime_str = "2022-01-14 15:52:04"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220114_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220114_v01.cdf"

    #starttime_str = "2022-01-14 23:26:43"
    #endtime_str = "2022-01-14 23:32:54"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220114_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220114_v01.cdf"

    #starttime_str = "2022-06-27 08:52:34"
    #endtime_str = "2022-06-27 08:58:47"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220627_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220627_v01.cdf"

    #starttime_str = "2022-01-25 16:28:23"
    #endtime_str = "2022-01-25 16:34:35"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220125_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220125_v01.cdf"

    #starttime_str = "2022-07-06/09:45:53"
    #endtime_str = "2022-07-06/09:52:11"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220706_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220706_v01.cdf"

    #starttime_str = "2022-06-25 06:27:08"	
    #endtime_str = "2022-06-25/06:33:20"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220625_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220625_v01.cdf"

    # TODO:has spike
    #starttime_str = "2022-06-25 09:29:30"
    #endtime_str = "2022-06-25 09:35:44"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220625_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220625_v01.cdf"

    # TODO: has unipolar spike in the middle
    #starttime_str = "2022-06-25 18:36:12"	
    #endtime_str = "2022-06-25/18:42:25"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220625_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220625_v01.cdf"

    # TODO: no spikes but first point rogue; example asked by V
    #starttime_str = "2022-06-23 02:27:59"
    #endtime_str = "2022-06-23 02:34:12"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220623_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220623_v01.cdf"

    # TODO: two unipolar spikes, example asked by V
    #starttime_str = "2022-06-23 04:00:07"
    #endtime_str = "2022-06-23 04:06:19"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220623_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220623_v01.cdf"

    # TODO: one spikes
    #starttime_str = "2022-06-23 08:34:48"
    #endtime_str = "2022-06-23 08:41:01"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220623_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220623_v01.cdf"

    #
    #starttime_str = "2022-06-23 11:36:39"
    #endtime_str = "2022-06-23 11:42:52"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220623_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220623_v01.cdf"

    #starttime_str = "2022-06-17 01:08:33"
    #endtime_str = "2022-06-17 01:14:47"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

    #starttime_str = "2022-06-17 04:13:18"
    #endtime_str = "2022-06-17 04:19:30"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

    #starttime_str = "2022-06-17 08:48:13"
    #endtime_str = "2022-06-17 08:54:24"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

    # has gap
    #starttime_str = "2022-06-17 11:50:13"
    #endtime_str = "2022-06-17 11:56:26"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

    #starttime_str = "2022-06-17 14:52:03"
    #endtime_str = "2022-06-17 14:58:16"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

    #starttime_str = "2022-06-17 20:59:35"
    #endtime_str = "2022-06-17 21:05:48"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

