import datetime as dt
import pandas as pd
import logging.config
from . import fgm_fsp_calib_prepos, fgm_fsp_calib, parameter
from .function import error, preprocess
from .function.postprocess import Bpara2Gthphi
import requests
import numpy as np
from scipy.interpolate import interp1d
from .function.Bplot import Gain_f
from .function.attitude import att_rot, att_loop, att_plot
from .eventlist import eventlist
from .function.output import output_txt
from .function.mva import mva

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

    """
    ====================
    settings 
    ====================
    """
    mission = "elb"
    csvpath = f"fgm_utils/temp/{mission}_fgm_data_availability.csv"
    elfin_url = "https://data.elfin.ucla.edu/"

    eventnum = 13
    starttime_str = eventlist[mission][eventnum]["starttime_str"]
    endtime_str = eventlist[mission][eventnum]["endtime_str"]
    f_all = eventlist[mission][eventnum].get("f_all", None)

    start_time = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d/%H:%M:%S"), starttime_str))
    end_time = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d/%H:%M:%S"), endtime_str))

    fgm_cdfdata_all = pd.DataFrame()
    att_cdfdata_all = pd.DataFrame()
    pos_cdfdata_all = pd.DataFrame()

    """
    ====================
    read all sci zones data 
    ====================
    """
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

        fgm_cdfdata_all = pd.concat([fgm_cdfdata_all, fgm_cdfdata])
        att_cdfdata_all = pd.concat([att_cdfdata_all, att_cdfdata])
        pos_cdfdata_all = pd.concat([pos_cdfdata_all, pos_cdfdata])
                
        [
            ctime_0, ctimestamp_0,
            fgs_ful_fgm_0th_x_0, fgs_ful_fgm_0th_y_0, fgs_ful_fgm_0th_z_0, 
            fgs_igrf_gei_x_0, fgs_igrf_gei_y_0, fgs_igrf_gei_z_0,
            att_gei_x_0, att_gei_y_0, att_gei_z_0,
            pos_gei_x_0, pos_gei_y_0, pos_gei_z_0] = fgm_fsp_calib_prepos(
                mission, start_time[i], end_time[i], fgm_cdfdata, att_cdfdata, pos_cdfdata)
        f_all_arry_0 = [f_all[i]] * len(fgs_ful_fgm_0th_x_0) if f_all is not None else [parameter.f] * len(fgs_ful_fgm_0th_x_0) 

        if i == 0: # first collection 
            [
                ctime, ctimestamp, 
                fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z,
                pos_gei_x, pos_gei_y, pos_gei_z, f_all_arry] = [
                    ctime_0, ctimestamp_0, 
                    fgs_ful_fgm_0th_x_0, fgs_ful_fgm_0th_y_0, fgs_ful_fgm_0th_z_0, 
                    fgs_igrf_gei_x_0, fgs_igrf_gei_y_0, fgs_igrf_gei_z_0, 
                    att_gei_x_0, att_gei_y_0, att_gei_z_0, 
                    pos_gei_x_0, pos_gei_y_0, pos_gei_z_0, f_all_arry_0]
            clip_start_idx = [0]
            clip_end_idx = [len(ctime)-1]

        else: # multiple collection
            ctime = np.concatenate((ctime, ctime_0 + ctime[-1] + 60))
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
            f_all_arry = np.concatenate([f_all_arry, f_all_arry_0])
            clip_start_idx.append(clip_end_idx[-1]+1) # start index of each sci zone
            clip_end_idx.append(len(ctime)-1) # end index of each sci zone

    if parameter.att_csv == True:
        """
        ====================
        att from txt: if true, replace the att with data in txt
        ====================
        """
        df = pd.read_csv(parameter.att_csv_filename) 
        if len(df) == len(att_gei_x):
            att_gei_x = df["att_gei_x"]
            att_gei_y = df["att_gei_y"]
            att_gei_z = df["att_gei_z"]

    """
    ====================
    start processing
    ====================
    """
    Bpara_out = []
    Gthphi_out = []
    f_out = []
    att_rot_out = [] 
    res_rot_out = []
    if parameter.att_loop == True:
        # att loop
        rot_len = len(np.arange(-parameter.att_loop_width, parameter.att_loop_width, parameter.att_loop_step))**2 + 1
        att_rot_out =  []
        att_gei_x_rot = np.random.rand(len(ctime), rot_len)
        att_gei_y_rot = np.random.rand(len(ctime), rot_len)
        att_gei_z_rot = np.random.rand(len(ctime), rot_len)

        """test
        """
        #att_gei_x =  np.ones_like(att_gei_x)
        #att_gei_y =  np.zeros_like(att_gei_y)
        #att_gei_z =  np.zeros_like(att_gei_z)

        # run att_gei no rotate as the first one
        att_gei_x_rot[:, 0] = att_gei_x
        att_gei_y_rot[:, 0] = att_gei_y
        att_gei_z_rot[:, 0] = att_gei_z

        # get attitude for all rot angle and concat sci zone
        for iclip_start_idx, iclip_end_idx in zip(clip_start_idx, clip_end_idx):         # loop over sci zone
            # rotate attitude vector for the start and end of each sci zone
            start_rotate_vector = att_loop(att_gei_x[iclip_start_idx], att_gei_y[iclip_start_idx], att_gei_z[iclip_start_idx], parameter.att_loop_width, parameter.att_loop_step)
            end_rotate_vector = att_loop(att_gei_x[iclip_end_idx], att_gei_y[iclip_end_idx], att_gei_z[iclip_end_idx], parameter.att_loop_width, parameter.att_loop_step)
            # loop over all rotate vectors 
            for idx, (istart_rotate_vector, iend_rotate_vector) in enumerate(zip(start_rotate_vector, end_rotate_vector)):
                interp_func = interp1d([ctime[iclip_start_idx], ctime[iclip_end_idx]], [istart_rotate_vector, iend_rotate_vector], axis = 0)
                interp_point = interp_func(ctime[iclip_start_idx:iclip_end_idx+1]) 
                att_gei_x_rot[iclip_start_idx:iclip_end_idx+1, idx+1] = interp_point[:, 0]
                att_gei_y_rot[iclip_start_idx:iclip_end_idx+1, idx+1] = interp_point[:, 1]
                att_gei_z_rot[iclip_start_idx:iclip_end_idx+1, idx+1] = interp_point[:, 2]
        
        if parameter.att_loop_figure == True:
            att_plot([att_gei_x[0], att_gei_y[0], att_gei_z[0]], list(zip(att_gei_x_rot[0,:], att_gei_y_rot[0,:], att_gei_z_rot[0,:])))
            breakpoint()

        if parameter.output == True:
            # output att for future run, you have to speicfy the number of the att to output.
            FGM_datetime = list(map(lambda ts: (
                dt.datetime.fromtimestamp(ctimestamp, tz=dt.timezone.utc) + dt.timedelta(seconds=ts)).strftime('%Y-%m-%d/%H:%M:%S.%f'), ctime))
            output_txt(
                FGM_datetime, 
                [att_gei_x_rot[:, 83], att_gei_y_rot[:, 83], att_gei_z_rot[:, 83]], 
                ['Timestamp','att_gei_x','att_gei_y','att_gei_z'], 
                title='att_gei')
            breakpoint()
        
        if parameter.mva_fgm == True:
            mva(ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z)
        
        # run calib with rotate att
        for idx in range(rot_len):
            [
                FGM_timestamp, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
                fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
                fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
                fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
                fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z, B_parameter]=fgm_fsp_calib(
                ctime, ctimestamp, f_all_arry,
                fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z,
                att_gei_x_rot[:,idx], att_gei_y_rot[:,idx], att_gei_z_rot[:,idx],
                pos_gei_x, pos_gei_y, pos_gei_z,
                logger, att_loop_idx = idx
            )
            Bpara_out.append(B_parameter)
            Gthphi_out.append(Bpara2Gthphi(B_parameter))
            f_out.append(f_all_arry[0]/ (np.pi / 180))
            res_rot = [(x**2 + y**2 + z**2)**0.5 for x, y, z in zip(fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z)]
            res_rot_out.append(np.median(res_rot))
            att_rot_out.append([att_gei_x_rot[0,idx], att_gei_y_rot[0,idx], att_gei_z_rot[0,idx]])

            
        Gthphi_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}_{mission}_attloop_Gthphi_{parameter.att_loop_width}_{parameter.att_loop_step}.csv"
        Bpara_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}_{mission}_attloop_Bpara_{parameter.att_loop_width}_{parameter.att_loop_step}.csv"
       
    elif parameter.f_loop == True:
        # f loop
        for if_loop in parameter.f_loop_value:
            f_all_arry = [if_loop] * len(fgs_ful_fgm_0th_x) 
            [
                FGM_timestamp, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
                fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
                fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
                fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
                fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z, B_parameter]=fgm_fsp_calib(
                ctime, ctimestamp, f_all_arry,
                fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z,
                att_gei_x, att_gei_y, att_gei_z,
                pos_gei_x, pos_gei_y, pos_gei_z,
                logger
            )
            Bpara_out.append(B_parameter)
            Gthphi_out.append(Bpara2Gthphi(B_parameter))
            f_out.append(if_loop/ (np.pi / 180))
            att_rot_out.append([att_gei_x[0], att_gei_y[0], att_gei_z[0]])
            res_rot = [(x**2 + y**2 + z**2)**0.5 for x, y, z in zip(fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z)]
            res_rot_out.append(np.median(res_rot))

        Gthphi_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}_{mission}_floop_Gthphi.csv"
        Bpara_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}_{mission}_floop_Bpara.csv"
        
    else:
        # no att loop, no floop, just a single run       
        [
            FGM_timestamp, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z, B_parameter]=fgm_fsp_calib(
            ctime, ctimestamp, f_all_arry,
            fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z, 
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z,
            att_gei_x, att_gei_y, att_gei_z,
            pos_gei_x, pos_gei_y, pos_gei_z,
            logger,
        )
        Bpara_out = [B_parameter]
        Gthphi_out = [Bpara2Gthphi(B_parameter)]
        f_out = [parameter.f/ (np.pi / 180)]
        att_rot_out = [-1]
        res_rot = [(x**2 + y**2 + z**2)**0.5 for x, y, z in zip(fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z)]
        res_rot_out.append(np.median(res_rot))
        att_rot_out.append([att_gei_x[0], att_gei_y[0], att_gei_z[0]])
        Gthphi_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}_{mission}_Gthphi.csv"
        Bpara_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}_{mission}_Bpara.csv"

    """
    ====================
    output
    ====================
    """

    Gthphi_columns = [ 
        'f', 'att_rot', 'res_med',
        'G1', 'G2', 'G3',
        'th1','th2','th3',
        'ph1','ph2','ph3',
        'O1/G1','O2/G2','O3/G3']
    Bpara_columns = [ 
        'f', 'att_rot', 'res_med',
        'G11', 'G12', 'G13',
        'G21','G22','G23',
        'G31','G32','G33',
        'O1','O2','O3']
    Gthphi_df = pd.DataFrame(columns = Gthphi_columns)
    Bpara_df = pd.DataFrame(columns = Bpara_columns)
  
    for iline in range(len(Bpara_out)):
        Gthphi_row = [f_out[iline], att_rot_out[iline], res_rot_out[iline], *Gthphi_out[iline]]
        Bpara_row = [f_out[iline], att_rot_out[iline], res_rot_out[iline], *Bpara_out[iline]]
        Gthphi_row = [dict(zip(Gthphi_columns, Gthphi_row))]
        Bpara_row = [dict(zip(Bpara_columns, Bpara_row))]

        Gthphi_df = pd.concat([Gthphi_df] + [pd.DataFrame(Gthphi_row)])
        Bpara_df = pd.concat([Bpara_df] + [pd.DataFrame(Bpara_row)])

    Gthphi_df.to_csv(Gthphi_filename, index=False)
    Bpara_df.to_csv(Bpara_filename, index=False)

    
   