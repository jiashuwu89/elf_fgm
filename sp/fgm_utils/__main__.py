import datetime as dt
import pandas as pd
import logging.config
from . import fgm_fsp_calib, parameter
from .function import error, preprocess
import requests

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

    mission = "ela"
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

    #starttime_str = ["2022-02-08/22:24:45"] #  still has gap in 1/80s spike. because t1 t2 not found, fixed
    #endtime_str = ["2022-02-08/22:30:58"]
    #starttime_str = ["2022-02-18/16:12:38"] # 1/80 spike. fixed
    #endtime_str = ["2022-02-18/16:18:50"]
    #starttime_str = ["2022-07-18/11:29:48"]
    #endtime_str = ["2022-07-18/11:36:00"]
    #starttime_str = ["2022-01-05/23:38:30"] # no unipolar spike
    #endtime_str = ["2022-01-05/23:44:41"]
    #starttime_str = ["2022-01-04/06:26:23"] # still has jagged
    #endtime_str = ["2022-01-04/06:32:34"]
    #starttime_str = ["2022-01-05/20:33:27"]  #still has gap in 1/80s spike. t1 found but not calibrated well
    #endtime_str = ["2022-01-05/20:39:40"]  
    #starttime_str = ["2022-01-01/11:52:27"]
    #endtime_str = ["2022-01-01/11:58:42"]
    #starttime_str = ["2022-08-17/23:24:26"] # ela
    #endtime_str = ["2022-08-17/23:30:39"]
    #starttime_str = ["2022-08-12/12:56:24"] # elb
    #endtime_str = ["2022-08-12/13:02:36"]
    #starttime_str = ["2022-08-07/21:55:33"] # elb
    #endtime_str = ["2022-08-07/22:01:51"]
    #starttime_str = ["2022-08-07/18:50:59"] # elb
    #endtime_str = ["2022-08-07/18:57:15"]
    #starttime_str = ["2022-01-03/16:42:58"] # step even without time spike
    #endtime_str = ["2022-01-03/16:49:10"]
    #starttime_str = ["2022-06-25/18:36:12"] # large amplitude
    #endtime_str = ["2022-06-25/18:42:25"]
    #starttime_str = ["2019-04-30/18:30:52"] # long collection
    #endtime_str = ["2019-04-30/20:00:04"] # XYrotate = -33
    starttime_str = ["2019-08-02/02:01:49"] # long collection
    endtime_str = ["2019-08-02/02:56:39"] # XYrotate = -34
    #starttime_str = ["2019-08-06/07:39:26"] # long collection
    #endtime_str = ["2019-08-06/08:33:53"] # XYrotate = -35
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
