import datetime as dt
import pandas as pd
import logging.config
from . import fgm_fsp_calib_prepos_wrapper, fgm_fsp_calib, parameter
from .function import error, preprocess
from .function.postprocess import Bpara2Gthphi
import numpy as np
import traceback
import sys
from .function.attitude import att_rot, att_loop, att_plot
from .eventlist import eventlist
from .function.output import output_txt
from .function.mva import mva
from .function.preprocess import get_fgmCSV
from .function.beta import get_betaCSV, get_beta
from scipy.interpolate import interp1d


def process_attloop(
        ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
        fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
        att_gei_x, att_gei_y, att_gei_z, 
        pos_gei_x, pos_gei_y, pos_gei_z,
        f_all_arry, clip_start_idx, clip_end_idx, logger):
    """processing code for loop of attitude
    """
    Bpara_out = []
    Gthphi_out = []
    f_out = []
    att_rot_out = [] 
    res_rot_out = []

    rot_len = len(np.arange(-parameter.att_loop_width, parameter.att_loop_width, parameter.att_loop_step))**2 + 1
    att_rot_out =  []
    att_gei_x_rot = np.random.rand(len(ctime), rot_len)
    att_gei_y_rot = np.random.rand(len(ctime), rot_len)
    att_gei_z_rot = np.random.rand(len(ctime), rot_len)

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
        if B_parameter != []:
            Bpara_out.append(B_parameter)
            Gthphi_out.append(Bpara2Gthphi(B_parameter))
            f_out.append(f_all_arry[0]/ (np.pi / 180))
            res_rot = [(x**2 + y**2 + z**2)**0.5 for x, y, z in zip(fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z)]
            res_rot_out.append(np.median(res_rot))
            att_rot_out.append([att_gei_x_rot[0,idx], att_gei_y_rot[0,idx], att_gei_z_rot[0,idx]])
        
    return Bpara_out, Gthphi_out, f_out, att_rot_out, res_rot_out


def process_floop(
        ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
        fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
        att_gei_x, att_gei_y, att_gei_z, 
        pos_gei_x, pos_gei_y, pos_gei_z,
        f_all_arry, logger):
    """processing code for loop of rotation angle
    """
    Bpara_out = []
    Gthphi_out = []
    f_out = []
    att_rot_out = [] 
    res_rot_out = []
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

    return Bpara_out, Gthphi_out, f_out, att_rot_out, res_rot_out


def process_single(
        ctime, ctimestamp, 
        fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
        fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
        att_gei_x, att_gei_y, att_gei_z, 
        pos_gei_x, pos_gei_y, pos_gei_z,
        f_all_arry, logger):
    """processing code for a single science zone
    """
    Bpara_out = []
    Gthphi_out = []
    f_out = []
    att_rot_out = [] 
    res_rot_out = []
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

    Bpara_out.append(B_parameter)
    Gthphi_out.append(Bpara2Gthphi(B_parameter)) if np.any(B_parameter) else []
    f_out.append(f_all_arry[0]/ (np.pi / 180))
    att_rot_out.append(-1)
    res_rot =  np.median([(x**2 + y**2 + z**2)**0.5 for x, y, z in zip(fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z)])
    res_rot_out.append(res_rot)
    #res_out = [(x**2 + y**2 + z**2)**0.5 for x, y, z in zip(fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z)]
    #print(f"median of residual: {np.median(res_out)}")
    
    return Bpara_out, Gthphi_out, f_out, att_rot_out, res_rot_out


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
    fgmcsvpath = f"fgm_utils/{mission}_fgm_data_availability.csv"
    betacsvpath = f"fgm_utils/orbits_fgm_cal_{mission}.csv"
    beta_df = get_betaCSV(betacsvpath)
    beta_out = []
    if parameter.batch_run == True:
        # not batch run, run one example from eventlist
        starttime_str = '2022-01-13/00:00:00'
        endtime_str = '2022-01-13/02:00:00'
        start_times, end_times = get_fgmCSV(fgmcsvpath, starttime_str, endtime_str)

        Bpara_out, Gthphi_out, f_out, att_rot_out, res_rot_out, beta_out, start_times_str, end_times_str = [], [], [], [], [], [], [], []

        for start_time, end_time in zip(start_times, end_times):
            start_time = [start_time]
            end_time = [end_time]
            """
            ====================
            read all sci zones data 
            ====================
            """
            try:
                [ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
                    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                    att_gei_x, att_gei_y, att_gei_z, 
                    pos_gei_x, pos_gei_y, pos_gei_z,
                    f_all_arry, clip_start_idx, clip_end_idx] = fgm_fsp_calib_prepos_wrapper(mission, start_time, end_time, None, logger)
                """
                ====================
                start processing
                ====================
                """
                # no att loop, no floop, just a single sci zone      
                Bpara_single, Gthphi_single, f_single, att_rot_single, res_rot_single = process_single(
                    ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
                    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                    att_gei_x, att_gei_y, att_gei_z, 
                    pos_gei_x, pos_gei_y, pos_gei_z,
                    f_all_arry, logger)
            except (error.preproc_resample_error, error.preproc_download_error) as e:
                # this error usually raise when 
                logger.error(e.__str__())
                Bpara_single, Gthphi_single, f_single, att_rot_single, res_rot_single = [], [], [], [], []
            except Exception as e:
                logger.error(f"❌ step 0 other error. Stop processing.")
                logger.error('\n'.join(traceback.format_exception(*sys.exc_info())))
                print('\n'.join(traceback.format_exception(*sys.exc_info())))
                Bpara_single, Gthphi_single, f_single, att_rot_single, res_rot_single = [], [], [], [], []
            
            Bpara_out.append(Bpara_single)
            Gthphi_out.append(Gthphi_single)
            f_out.append(f_single)
            att_rot_out.append(att_rot_single)
            res_rot_out.append(res_rot_single)

            """
            ====================
            get beta
            ====================
            """
            mid_time = start_time[0] + 0.5*(end_time[0] - start_time[0])
            beta_angle = get_beta(beta_df, mid_time)
            beta_out.append(beta_angle)

            """
            ====================
            get date array
            ====================
            """
            start_time_str = dt.datetime.strftime(start_time[0], "%Y-%m-%d/%H:%M:%S")
            end_time_str = dt.datetime.strftime(end_time[0], "%Y-%m-%d/%H:%M:%S")
            start_times_str.append(start_time_str)
            end_times_str.append(end_time_str)

        
        Gthphi_filename = f"fgm_utils/fitting_csv/{starttime_str[0:10]}_{starttime_str[11:13]}{starttime_str[14:16]}_{mission}_Gthphi.csv"
        Bpara_filename = f"fgm_utils/fitting_csv/{starttime_str[0:10]}_{starttime_str[11:13]}{starttime_str[14:16]}_{mission}_Bpara.csv"
    else:
        eventnum = 63
        starttime_str = eventlist[mission][eventnum]["starttime_str"]
        endtime_str = eventlist[mission][eventnum]["endtime_str"]
        f_all = eventlist[mission][eventnum].get("f_all", None)

        start_time = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d/%H:%M:%S"), starttime_str))
        end_time = list(map(lambda ts: dt.datetime.strptime(ts, "%Y-%m-%d/%H:%M:%S"), endtime_str))

        """
        ====================
        read all sci zones data 
        ====================
        """
        [ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
            att_gei_x, att_gei_y, att_gei_z, 
            pos_gei_x, pos_gei_y, pos_gei_z,
            f_all_arry, clip_start_idx, clip_end_idx] = fgm_fsp_calib_prepos_wrapper(mission, start_time, end_time, f_all, logger)

        if parameter.att_split == True and parameter.att_split_num is None and parameter.att_split_idx is None:
            parameter.att_split_idx = clip_start_idx

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
        if parameter.att_loop == True:
            # loop of attitude
            Bpara_out, Gthphi_out, f_out, att_rot_out, res_rot_out, Gthphi_filename, Bpara_filename = process_attloop(
                ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z, 
                pos_gei_x, pos_gei_y, pos_gei_z,
                f_all_arry, clip_start_idx, clip_end_idx) 
            
            Gthphi_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}{starttime_str[0][14:16]}_{mission}_attloop_Gthphi_{parameter.att_loop_width}_{parameter.att_loop_step}.csv"
            Bpara_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}{starttime_str[0][14:16]}_{mission}_attloop_Bpara_{parameter.att_loop_width}_{parameter.att_loop_step}.csv"

        elif parameter.f_loop == True:
            # loop of f
            Bpara_out, Gthphi_out, f_out, att_rot_out, res_rot_out, Gthphi_filename, Bpara_filename = process_floop(
                ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z, 
                pos_gei_x, pos_gei_y, pos_gei_z,
                f_all_arry, logger)
            
            Gthphi_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}{starttime_str[0][14:16]}_{mission}_floop_Gthphi.csv"
            Bpara_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}{starttime_str[0][14:16]}_{mission}_floop_Bpara.csv"

        else:
            # no att loop, no floop, just a single sci zone      
            Bpara_out, Gthphi_out, f_out, att_rot_out, res_rot_out = process_single(
                ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z,
                fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
                att_gei_x, att_gei_y, att_gei_z, 
                pos_gei_x, pos_gei_y, pos_gei_z,
                f_all_arry, logger)
            
            Gthphi_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}{starttime_str[0][14:16]}_{mission}_Gthphi.csv"
            Bpara_filename = f"fgm_utils/fitting_csv/{starttime_str[0][0:10]}_{starttime_str[0][11:13]}{starttime_str[0][14:16]}_{mission}_Bpara.csv"

        start_times_str, end_times_str = [], []
        start_times_str.append(starttime_str[0])
        end_times_str.append(endtime_str[-1])
    logger.info("Processing is done. Now output to file ...")
    """
    ====================
    get beta
    ====================
    """
    mid_time = start_time[0] + 0.5*(end_time[-1] - start_time[0])
    beta_angle = get_beta(beta_df, mid_time)
    beta_out.append(beta_angle)
    
    """
    ====================
    output
    ====================
    """

    Gthphi_columns = [ 
        'start_time', 'end_time', 'f', 'att_rot', 'beta_ang', 'res_med',
        'G1', 'G2', 'G3',
        'th1','th2','th3',
        'ph1','ph2','ph3',
        'O1/G1','O2/G2','O3/G3']
    Bpara_columns = [ 
        'start_time', 'end_time', 'f', 'att_rot', 'beta_ang', 'res_med',
        'G11', 'G12', 'G13', 'O1',
        'G21','G22','G23', 'O2',
        'G31','G32','G33', 'O3']
    Gthphi_df = pd.DataFrame(columns = Gthphi_columns)
    Bpara_df = pd.DataFrame(columns = Bpara_columns)

    for iline in range(len(Bpara_out)):
        Gthphi_row = [start_times_str[iline], end_times_str[iline], f_out[iline], att_rot_out[iline], beta_out[iline], res_rot_out[iline], *Gthphi_out[iline]]
        Bpara_row = [start_times_str[iline], end_times_str[iline], f_out[iline], att_rot_out[iline], beta_out[iline], res_rot_out[iline], *Bpara_out[iline]]
        Gthphi_row = [dict(zip(Gthphi_columns, Gthphi_row))]
        Bpara_row = [dict(zip(Bpara_columns, Bpara_row))]

        Gthphi_df = pd.concat([Gthphi_df] + [pd.DataFrame(Gthphi_row)])
        Bpara_df = pd.concat([Bpara_df] + [pd.DataFrame(Bpara_row)])

    Gthphi_df.to_csv(Gthphi_filename, index=False)
    Bpara_df.to_csv(Bpara_filename, index=False)

    logger.info("Done.")

    
   