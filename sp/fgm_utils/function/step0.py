from .. import parameter
from . import cross_time, error, coordinate, calibration, ctime_spike, Bplot
from .ctime_spike_80 import spike_sinefit_80

def step0(
    ctime, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, 
    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
    att_gei_x, att_gei_y, att_gei_z,
    datestr, logger):

    """
        0. precalibration: time calibration
    """
    """
        0.1 corr - cross time determination
    """
    try:
        [
            cross_times_0th_1, cross_times_0th_1_mids, 
            T_spins_0th_1, w_syn_0th_1] = cross_time.cross_time_stage_1(
            ctime, fgs_ful_fgm_1st_z,
        )
    except error.CrossTime1Error as e:
        raise error.CrossTime1Error(0)

    [
        cross_times_0th_2, cross_times_0th_2_mids, 
        T_spins_0th_2, w_syn_0th_2] = cross_time.cross_time_stage_2(
        ctime, fgs_ful_fgm_1st_z, cross_times_0th_1, T_spins_0th_1,
    )

    [
        cross_times_0th_3, T_spins_0th_3, w_syn_0th_3] = cross_time.cross_time_stage_3(
            ctime, fgs_ful_fgm_1st_z, cross_times_0th_2, T_spins_0th_2
    )
    logger.debug(f"[0.1] zero crossing is done.")

    """
        0.2 corr - phase angle integration
    """
    [
        phi_corr, cross_times_corr, w_syn_d_corr, T_spins_d_corr, cross_times_corr_fit, w_syn_d_corr_fit] = cross_time.phase_integration(
        ctime, cross_times_0th_1, cross_times_0th_1_mids, w_syn_0th_1, T_spins_0th_1,
        cross_times_0th_2, cross_times_0th_2_mids, w_syn_0th_2, T_spins_0th_2,
        cross_times_0th_3, w_syn_0th_3, T_spins_0th_3,
    )    
    logger.debug(f"[0.2] phase angle is done.")

    """
        0.3 IGRF coorindate transformation: gei -> dmxl -> smxl -> fgm
    """
    [
        DMXL_2_GEI, GEI_2_DMXL] = coordinate.dmxl_gei_matrix(
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, att_gei_x, att_gei_y, att_gei_z
    )

    # B igrf rotate from gei to dmxl
    [
        fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z] = coordinate.gei2dmxl(
            fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, GEI_2_DMXL
    )

    # B igrf rotate from dmxl to smxl
    [
        fgs_igrf_smxl_x, fgs_igrf_smxl_y, fgs_igrf_smxl_z] = coordinate.dmxl2smxl(
            fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, phi_corr
    )

    # B igrf rotate from smxl to fgm
    [
        fgs_igrf_fgm_x, fgs_igrf_fgm_y, fgs_igrf_fgm_z] = coordinate.smxl2fgm(
            fgs_igrf_smxl_x, fgs_igrf_smxl_y, fgs_igrf_smxl_z
    )

    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, [fgs_ful_fgm_1st_x, fgs_igrf_fgm_x], [fgs_ful_fgm_1st_y, fgs_igrf_fgm_y], 
    #        [fgs_ful_fgm_1st_z, fgs_igrf_fgm_z], plot3 = True, title="ful_igrf_fgm_before1stcali")      

    logger.debug(f"[0.3] ctime spike correction is done.")

    """
        0.4 use igrf to calibrate fgs data
    """
    [
        fgs_ful_fgm_x, fgs_ful_fgm_y, fgs_ful_fgm_z, B_parameter] = calibration.calib_leastsquare(
        fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, fgs_igrf_fgm_x, fgs_igrf_fgm_y, fgs_igrf_fgm_z
    )  

    logger.debug(f"[0.4] igrf and B full field calibrated in fgm coordinate.")       

    """
        0.5 ctime correction
    """
    fgs_res_fgm_x = fgs_ful_fgm_x - fgs_igrf_fgm_x
    fgs_res_fgm_y = fgs_ful_fgm_y - fgs_igrf_fgm_y
    fgs_res_fgm_z = fgs_ful_fgm_z - fgs_igrf_fgm_z

    ctime, ctime_idx, ctime_idx_flag, ctime_idx_timediff, spike_ctime_idxs = ctime_spike.ctime_calib(
            ctime, fgs_res_fgm_x, fgs_res_fgm_y, fgs_res_fgm_z, cross_times_corr, logger = logger, datestr = datestr
    )

    ctime_idx_time = ctime[ctime_idx]

    # fit 1/80 s spike with sine 
    if parameter.ctime_correct_80 == True:
        fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z = spike_sinefit_80(
            ctime, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, spike_ctime_idxs
    )

    logger.debug(f"[0.5] ctime spike correction is done.")
    if parameter.makeplot == True:
        #Bplot.B_ctime_plot(ctime, fgs_res_fgm_x, fgs_res_fgm_y, fgs_res_fgm_z, ctime_idx_time = ctime[ctime_idx], ctime_idx_flag = ctime_idx_flag)
        Bplot.ctimediff_plot(ctime, ctime_idx, ctime_idx_flag, datestr = datestr)

    return [fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, ctime_idx, ctime_idx_time, ctime_idx_flag, ctime_idx_timediff]