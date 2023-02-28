import datetime as dt
import pandas as pd
import logging.config
from . import fgm_fsp_calib, parameter
from .function import error, preprocess, fgm_fsp_calib_floop
import requests
import numpy as np
from .function.Bplot import Gain_f

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

    #starttime_str = ["2019-04-30/18:30:52"] # long collection
    #endtime_str = ["2019-04-30/20:00:04"]
    #starttime_str = ["2019-08-02/02:01:49"] # long collection
    #endtime_str = ["2019-08-02/02:56:39"]
    #starttime_str = ["2019-08-06/07:39:26"] # long collection
    #endtime_str = ["2019-08-06/08:33:53"]
    #starttime_str = ["2022-01-16/22:53:11"] # o largerthan w
    #endtime_str = ["2022-01-16/22:59:22"]

    #starttime_str = ["2021-03-21/08:27:09"] # Qi murong ela
    #endtime_str = ["2021-03-21/08:33:23"]    
    #starttime_str = ["2021-03-21/10:00:07"] # Qi murong ela
    #endtime_str = ["2021-03-21/10:06:20"]    
    #starttime_str = ["2021-03-21/11:33:13"] # Qi murong ela
    #endtime_str = ["2021-03-21/11:39:25"]
    #starttime_str = ["2021-03-22/13:59:54"] # Qi murong ela
    #endtime_str = ["2021-03-22/14:06:07"]     
    #starttime_str = ["2021-03-21/06:30:33"] # Qi murong elb
    #endtime_str = ["2021-03-21/06:36:46"]     
    #starttime_str = ["2021-03-22/12:03:10"] # Qi murong elb
    #endtime_str = ["2021-03-22/12:09:23"]  
    #starttime_str = ["2022-07-20/13:32:24"] # Vassilis event
    #endtime_str = ["2022-07-20/13:38:27"] 
    #starttime_str = ["2021-12-01/22:19:07"] # check error events
    #endtime_str = ["2021-12-01/22:25:10"]
    #starttime_str = ["2021-03-14/15:45:48"] # yangyang gap event
    #endtime_str = ["2021-03-14/15:52:02"]   
    #starttime_str = ["2022-09-07/19:13:44"] # yangyang gap event
    #endtime_str = ["2022-09-07/19:19:57"] 
    starttime_str = ["2022-09-10/23:42:55"] # jiang event ela o larger than w
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
    #starttime_str = ["2019-06-29/02:55:49"] # V first event in phase
    #endtime_str = ["2019-06-29/03:00:03"]

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
