from .. import parameter
from .ctime_spike_80 import find_closest
from .error import fsp_spike_del_error
from .detrend import del_rogue, delete_data
import numpy as np
from typing import List, Tuple

def fsp_spike_del(
    ctime, ctime_idx, ctime_idx_flag, ctime_idx_timediff, 
    cross_times_calib, w_syn_d_calib, DMXL_2_GEI_fsp,
    fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z, 
    fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y,fgs_fsp_igrf_dmxl_z, 
    fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
    fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
    fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
    logger):

    """type 3, 2.5s pink spike
    """
    cross_times_calib_del = []
    w_avg = np.median(w_syn_d_calib)
    if parameter.fsp_spike_del_type3 == True:
        # delete 2.5 purple spike if necenssary
        ctime_idx_time_2 = ctime[ctime_idx[ctime_idx_flag == 3]]
        ctime_idx_timediffs = ctime_idx_timediff[ctime_idx_flag == 3]
        fgs_fsp_res_dmxl_norm = fgs_fsp_res_dmxl_x**2 + fgs_fsp_res_dmxl_y**2 + fgs_fsp_res_dmxl_z**2
        for ctime_idx_time_idx, ctime_idx_time_val in enumerate(ctime_idx_time_2):
            try:
                idx1_del = np.where(cross_times_calib > ctime_idx_time_val - 1.5*np.pi/w_avg)[0][0]
                idx2_del = np.where(cross_times_calib < ctime_idx_time_val + ctime_idx_timediffs[ctime_idx_time_idx] + 1.5*np.pi/w_avg)[0][-1]
                idxs_del = range(idx1_del, idx2_del) if idx2_del + 1 >= len(cross_times_calib) else range(idx1_del, idx2_del+1)
                [cross_times_calib_del.append(idxs_del_i) for idxs_del_i in idxs_del]
                logger.debug(f"[POSTPROCESS] FSP spike delete success: 2.5s pink spike {ctime_idx_time_val} type 3!")
                """ delete according to std doesn't have a good performance
                # std for 10 spins
                idx1 = np.where(cross_times_calib > ctime_idx_time_val - 5*np.pi/w_avg)[0][0]
                idx2 = np.where(cross_times_calib < ctime_idx_time_val + ctime_idx_timediffs[ctime_idx_time_idx] + 5*np.pi/w_avg)[0][-1]
                idxs = range(idx1, idx2) if idx2 + 1 >= len(cross_times_calib) else range(idx1, idx2+1) # index for 10 spins 
                clip = fgs_fsp_res_dmxl_norm[idxs]
                clip_ctime = cross_times_calib[idxs]
                clip_std = np.std(clip)

                clip3 = np.delete(clip, idxs_del)
                clip3_std = np.std(clip3)
                if clip_std > clip3_std:
                    [cross_times_calib_del.append(idxs_del_i) for idxs_del_i in idxs_del]
                    logger.info(f"zero crossing for 2.5s purple spike {ctime_idx_time_val} is deleted!")
                else:
                    #print(f"2.5s purple spike {ctime_idx_time_val} clip1 std:{clip_std}")  # this is the std around orange spike
                    idx_del_25 = np.zeros(len(idxs_del))
                    for idxs_del_i, idxs_del_val in enumerate(idxs_del):
                        idx_del_25[idxs_del_i] = ctime_spike_80.find_closest(clip_ctime, cross_times_calib[idxs_del_val])[0]
                        clip2 = np.delete(clip, idx_del_25[idxs_del_i])
                        clip2_std = np.std(clip2)
                        #print(f"2.5s purple spike {cross_times_calib[idxs_del_i]} clip2 std:{clip2_std}") # this is the std around orange spike if exclude the spike 

                        if clip_std > clip2_std:
                            cross_times_calib_del.append(idxs_del_val)
                            logger.info(f"zero crossing for 2.5s purple spike {cross_times_calib[idxs_del_val]} is deleted!")
                """
            except:
                #breakpoint()
                logger.warning(f"❗[POSTPROCESS] FSP spike delete error: 2.5s pink spike {ctime_idx_time_val} type 3!")
                if len(fgs_fsp_res_dmxl_x) < 6:
                    raise fsp_spike_del_error
                continue
    
        [
            cross_times_calib, DMXL_2_GEI_fsp, 
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
            ] = delete_data(
                cross_times_calib_del, cross_times_calib, DMXL_2_GEI_fsp,
                fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
                fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
                fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
                fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
                fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
        )
        
        
    """type 2, 1/80s spike orange
    """
    cross_times_calib_del = []
    if parameter.fsp_spike_del_type2 == True:
    # delete 1/80 orange spike if necenssary
        ctime_idx_time_2 = ctime[ctime_idx[ctime_idx_flag == 2]]
        ctime_idx_timediffs = ctime_idx_timediff[ctime_idx_flag == 2]
        fgs_fsp_res_dmxl_norm = fgs_fsp_res_dmxl_x**2 + fgs_fsp_res_dmxl_y**2 + fgs_fsp_res_dmxl_z**2
        cross_times_calib_del = []
        for ctime_idx_time_idx, ctime_idx_time_val in enumerate(ctime_idx_time_2):
            try:
                idx1 = np.where(cross_times_calib > ctime_idx_time_val - 15*np.pi/w_avg)[0][0]
                idx2 = np.where(cross_times_calib < ctime_idx_time_val + ctime_idx_timediffs[ctime_idx_time_idx] + 15*np.pi/w_avg)[0][-1]
                idxs = range(idx1, idx2) if idx2 + 1 >= len(cross_times_calib) else range(idx1, idx2+1)
                clip = fgs_fsp_res_dmxl_norm[idxs]
                clip_ctime = cross_times_calib[idxs]
                clip_std = np.std(clip)
                #logger.info(f"1/80s orange spike {ctime_idx_time_val} clip1 std:{clip_std}")  # this is the std around orange spike
                
                idx = find_closest(clip_ctime, ctime_idx_time_val)[0]
                clip2 = np.delete(clip, idx)
                clip2_std = np.std(clip2)
                #logger.info(f"1/80s orange spike {ctime_idx_time_val} clip2 std:{clip2_std}") # this is the std around orange spike if exclude the spike 
                if clip_std > clip2_std:
                    idx_ctime = find_closest(cross_times_calib, ctime_idx_time_val)[0]
                    cross_times_calib_del.append(idx_ctime)
                    logger.debug(f"[POSTPROCESS] FSP spike delete success: 1/80s orange spike {ctime_idx_time_val} type 2.")
                if idx > 0 and idx < len(clip_ctime) and clip_ctime[idx] - clip_ctime[idx-1] < 3:
                    clip2 = np.delete(clip, idx-1)
                    clip2_std = np.std(clip2)
                    if clip_std > clip2_std:
                        idx_ctime = find_closest(cross_times_calib, ctime_idx_time_val)[0]
                        cross_times_calib_del.append(idx_ctime-1)
                if idx < len(clip_ctime)-1 and clip_ctime[idx+1] - clip_ctime[idx] < 3:
                    clip2 = np.delete(clip, idx+1)
                    clip2_std = np.std(clip2)
                    if clip_std > clip2_std:
                        idx_ctime = find_closest(cross_times_calib, ctime_idx_time_val)[0]
                        cross_times_calib_del.append(idx_ctime+1)
            except:
                #breakpoint()
                logger.warning(f"❗[POSTPROCESS] FSP spike delete error:  1/80s orange spike {ctime_idx_time_val} type 2! ")
                if len(fgs_fsp_res_dmxl_x) < 6:
                    raise fsp_spike_del_error
                continue

        [
            cross_times_calib, DMXL_2_GEI_fsp,
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
            ] = delete_data(
                cross_times_calib_del, cross_times_calib, DMXL_2_GEI_fsp,
                fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
                fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
                fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
                fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
                fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
        )


    """type 4, purple spike, other gaps
    """
    cross_times_calib_del = []
    if parameter.fsp_spike_del_type4 == True:
        ctime_idx_time_2 = ctime[ctime_idx[ctime_idx_flag == 4]]
        ctime_idx_timediffs = ctime_idx_timediff[ctime_idx_flag == 4]
        fgs_fsp_res_dmxl_norm = fgs_fsp_res_dmxl_x**2 + fgs_fsp_res_dmxl_y**2 + fgs_fsp_res_dmxl_z**2
        cross_times_calib_del = []
        for ctime_idx_time_idx, ctime_idx_time_val in enumerate(ctime_idx_time_2):
            try:
                
                # pick a interval around ctime_idx and calc std
                idx1 = np.where(cross_times_calib > ctime_idx_time_val - 30*np.pi/w_avg)[0][0]
                idx2 = np.where(cross_times_calib < ctime_idx_time_val + ctime_idx_timediffs[ctime_idx_time_idx] + 30*np.pi/w_avg)[0][-1]
                idxs = range(idx1, idx2) if idx2 + 1 >= len(cross_times_calib) else range(idx1, idx2+1)
                idxs= list(idxs)
                idxs_inters = list(set(cross_times_calib_del) & set(idxs))
                idxs = [i for i in idxs if i not in idxs_inters]
                clip = fgs_fsp_res_dmxl_norm[idxs]
                clip_ctime = cross_times_calib[idxs]
                clip_mean = np.mean(clip)
                clip_std = np.std(clip)
                #logger.info(f"1/80s orange spike {ctime_idx_time_val} clip1 std:{clip_std}")  # this is the std around orange spike
                
                idx = find_closest(clip_ctime, ctime_idx_time_val)[0] # position of ctime_idx in clip
                # if the ctime_idx is too large than three sigma
                if clip[idx] > clip_mean + clip_std*3 or clip[idx] < clip_mean - clip_std*3:
                    idx_ctime = find_closest(cross_times_calib, ctime_idx_time_val)[0] # position of cime_idx in cross_times_calib
                    cross_times_calib_del.append(idx_ctime)
                    logger.debug(f"[POSTPROCESS] FSP spike delete success: other purple spike {ctime_idx_time_val} type 4.")
                # check the point before ctime_idx and within a spin period, if std decrease after delete this point, then delete it   
                if idx > 0 and idx < len(clip_ctime) and clip_ctime[idx] - clip_ctime[idx-1] < 3:
                    if clip[idx-1] > clip_mean + clip_std*3 or clip[idx-1] < clip_mean - clip_std*3:
                        idx_ctime = find_closest(cross_times_calib, ctime_idx_time_val)[0] # position of cime_idx in cross_times_calib
                        cross_times_calib_del.append(idx_ctime-1)
                # check the point after ctime_idx and within a spin period, if std decrease after delete this point, then delete it   
                if idx < len(clip_ctime)-1 and clip_ctime[idx+1] - clip_ctime[idx] < 3: 
                    if clip[idx+1] > clip_mean + clip_std*3 or clip[idx+1] < clip_mean - clip_std*3:
                        idx_ctime = find_closest(cross_times_calib, ctime_idx_time_val)[0]
                        cross_times_calib_del.append(idx_ctime+1)
            except:
                #breakpoint()
                logger.warning(f"❗[POSTPROCESS] FSP spike delete error: other purple spike {ctime_idx_time_val} type 4 !")
                if len(fgs_fsp_res_dmxl_x) < 6:
                    raise fsp_spike_del_error
                continue

        [
            cross_times_calib, DMXL_2_GEI_fsp,
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
            ] = delete_data(
                cross_times_calib_del, cross_times_calib, DMXL_2_GEI_fsp,
                fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
                fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
                fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
                fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
                fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
        )

    """type 5, green spike, < 0.1 
    """
    cross_times_calib_del = []
    if parameter.fsp_spike_del_type4 == True:
        ctime_idx_time_2 = ctime[ctime_idx[ctime_idx_flag == 5]]
        ctime_idx_timediffs = ctime_idx_timediff[ctime_idx_flag == 5]
        fgs_fsp_res_dmxl_norm = fgs_fsp_res_dmxl_x**2 + fgs_fsp_res_dmxl_y**2 + fgs_fsp_res_dmxl_z**2
        cross_times_calib_del = []
        for ctime_idx_time_idx, ctime_idx_time_val in enumerate(ctime_idx_time_2):
            try:
                idx1 = np.where(cross_times_calib > ctime_idx_time_val - 15*np.pi/w_avg)[0][0]
                idx2 = np.where(cross_times_calib < ctime_idx_time_val + ctime_idx_timediffs[ctime_idx_time_idx] + 15*np.pi/w_avg)[0][-1]
                idxs = range(idx1, idx2) if idx2 + 1 >= len(cross_times_calib) else range(idx1, idx2+1)
                clip = fgs_fsp_res_dmxl_norm[idxs]
                clip_ctime = cross_times_calib[idxs]
                clip_std = np.std(clip)
                #logger.info(f"1/80s orange spike {ctime_idx_time_val} clip1 std:{clip_std}")  # this is the std around orange spike
                
                idx = find_closest(clip_ctime, ctime_idx_time_val)[0]
                clip2 = np.delete(clip, idx)
                clip2_std = np.std(clip2)
                #logger.info(f"1/80s orange spike {ctime_idx_time_val} clip2 std:{clip2_std}") # this is the std around orange spike if exclude the spike 
                if clip_std > clip2_std:
                    idx_ctime = find_closest(cross_times_calib, ctime_idx_time_val)[0]
                    cross_times_calib_del.append(idx_ctime)
                    logger.debug(f"[POSTPROCESS] FSP spike delete success: < 0.1s green spike {ctime_idx_time_val} type 5. ")
                if idx > 0 and idx < len(clip_ctime) and clip_ctime[idx] - clip_ctime[idx-1] < 3:
                    clip2 = np.delete(clip, idx-1)
                    clip2_std = np.std(clip2)
                    if clip_std > clip2_std:
                        idx_ctime = find_closest(cross_times_calib, ctime_idx_time_val)[0]
                        cross_times_calib_del.append(idx_ctime-1)
                if idx < len(clip_ctime)-1 and clip_ctime[idx+1] - clip_ctime[idx] < 3:
                    clip2 = np.delete(clip, idx+1)
                    clip2_std = np.std(clip2)
                    if clip_std > clip2_std:
                        idx_ctime = find_closest(cross_times_calib, ctime_idx_time_val)[0]
                        cross_times_calib_del.append(idx_ctime+1)
            except:
                #breakpoint()
                logger.warning(f"❗[POSTPROCESS] FSP spike delete error: < 0.1s green spike {ctime_idx_time_val} type 5!")
                continue
        
        if len(fgs_fsp_res_dmxl_x) < 6:
            raise fsp_spike_del_error

        [
            cross_times_calib, DMXL_2_GEI_fsp,
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
            ] = delete_data(
                cross_times_calib_del, cross_times_calib, DMXL_2_GEI_fsp,
                fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
                fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
                fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
                fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
                fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
        )

    """delete rogue points
    """
    if parameter.del_rogue_fsp == True:
        del_index = del_rogue(
            cross_times_calib, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z
        )
        [
            cross_times_calib, DMXL_2_GEI_fsp,
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
            ] = delete_data(
                del_index, cross_times_calib, DMXL_2_GEI_fsp,
                fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
                fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
                fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
                fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
                fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
        )
    
    # delete rogue points again in fsp data
    if parameter.del_rogue_fsp == True:
        del_index = del_rogue(
            cross_times_calib, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z
        )
        [
            cross_times_calib, DMXL_2_GEI_fsp,
            fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
            fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
            fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
            fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
            fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
            ] = delete_data(
                del_index, cross_times_calib, DMXL_2_GEI_fsp,
                fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
                fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z,
                fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
                fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
                fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
        )
        logger.debug(f"[POSTPROCESS] detrend successfully done. ")
    
    return [cross_times_calib, DMXL_2_GEI_fsp, fgs_fsp_res_dmxl_x, fgs_fsp_res_dmxl_y, fgs_fsp_res_dmxl_z,
        fgs_fsp_igrf_dmxl_x, fgs_fsp_igrf_dmxl_y, fgs_fsp_igrf_dmxl_z, 
        fgs_fsp_res_gei_x, fgs_fsp_res_gei_y, fgs_fsp_res_gei_z,
        fgs_fsp_igrf_gei_x, fgs_fsp_igrf_gei_y, fgs_fsp_igrf_gei_z,
        fgs_fsp_res_dmxl_trend_x, fgs_fsp_res_dmxl_trend_y, fgs_fsp_res_dmxl_trend_z,
    ]

def Bpara2Gthphi(B_parameter: List[float]) -> Tuple[float]:

    [G11, G12, G13, O1, G21, G22, G23, O2, G31, G32, G33, O3] = [*B_parameter]   
    G1 = (G11**2 + G12**2 + G13**2)**0.5 
    G2 = (G21**2 + G22**2 + G23**2)**0.5
    G3 = (G31**2 + G32**2 + G33**2)**0.5
    
    th1 = np.degrees(np.arccos(G13/G1))
    th2 = np.degrees(np.arccos(G23/G2))
    th3 = np.degrees(np.arccos(G33/G3))

    ph1 = np.degrees(np.arctan2(G12, G11))
    ph2 = np.degrees(np.arctan2(G22, G21))
    ph3 = np.degrees(np.arctan2(G32, G31))
          
    return [G1, G2, G3, th1, th2, th3, ph1, ph2, ph3, O1/G1, O2/G2, O3/G3]


def Gthphi2Bpara(Gthphi: List[float]) -> Tuple[float]:

    [G1, G2, G3, th1, th2, th3, ph1, ph2, ph3, O1G1, O2G2, O3G3] = [*Gthphi]  
    th1, th2, th3, ph1, ph2, ph3 = np.radians([th1, th2, th3, ph1, ph2, ph3])

    G11 = G1*np.sin(th1)*np.cos(ph1)
    G12 = G1*np.sin(th1)*np.sin(ph1)
    G13 = G1*np.cos(th1)

    G21 = G2*np.sin(th2)*np.cos(ph2)
    G22 = G2*np.sin(th2)*np.sin(ph2)
    G23 = G2*np.cos(th2)

    G31 = G3*np.sin(th3)*np.cos(ph3)
    G32 = G3*np.sin(th3)*np.sin(ph3)
    G33 = G3*np.cos(th3)

    O1 = O1G1 * G1
    O2 = O2G2 * G2
    O3 = O3G3 * G3
          
    return [G11, G12, G13, O1, G21, G22, G23, O2, G31, G32, G33, O3]


def print_parameter(parameter_set: List, Gthphi: bool = False):
    """print parameters to screen
    """
    if Gthphi == False:
        # parameters are G11, G12, G13 ...
        (G11, G12, G13, O1, G21, G22, G23, O2, G31, G32, G33, O3) = parameter_set[0] 
        print('===============calib para================')  
        print(f'G11: {G11} G12:{G12} G13:{G13}\n')
        print(f'G21: {G21} G22:{G22} G23:{G23}\n')
        print(f'G31: {G31} G32:{G32} G33:{G33}\n')
        print(f'O1: {O1} O2:{O2} O3:{O3}\n')
        print('=============calib para end==============')
    else:
        # parameters are G th phi
        (G1, G2, G3, th1, th2, th3, ph1, ph2, ph3, O1G1, O2G2, O3G3) = parameter_set[0]  
        print('===============calib para================')
        print(f'G1: {G1} G2:{G2} G3:{G3}\n')
        print(f'th1: {th1} th2:{th2} th3:{th3}\n')
        print(f'ph1: {ph1} ph2:{ph2} ph3:{ph3}\n')
        print(f'O1/G1: {O1G1} O2/G2:{O2G2} O3/G3:{O3G3}\n')
        print('=============calib para end==============')
    return
          


    return