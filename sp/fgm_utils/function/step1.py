from .. import parameter
from . import cross_time, error, coordinate, calibration, ctime_spike, Bplot
from .ctime_spike_80 import spike_sinefit_80

def step1(
    ctime, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, 
    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
    att_gei_x, att_gei_y, att_gei_z,
    datestr, logger, ctime_idx, ctime_idx_time, ctime_idx_flag, ctime_idx_timediff):   

    """
        # 1. first run
    """
    [
        cross_times_1st_1, cross_times_1st_1_mids, 
        T_spins_1st_1, w_syn_1st_1] = cross_time.cross_time_stage_1(
        ctime, fgs_ful_fgm_1st_z,
    )
    [
        cross_times_1st_2, cross_times_1st_2_mids, 
        T_spins_1st_2, w_syn_1st_2] = cross_time.cross_time_stage_2(
        ctime, fgs_ful_fgm_1st_z, cross_times_1st_1, T_spins_1st_1,
    )
    [
        cross_times_1st_3, T_spins_1st_3, w_syn_1st_3] = cross_time.cross_time_stage_3(
            ctime, fgs_ful_fgm_1st_z, cross_times_1st_2, T_spins_1st_2
    )
    logger.debug(f"[1.1] zero crossing is done. ")

    """
        # 1.1 corr - phase angle integration
    """
    [
        phi_1st, cross_times_1st, w_syn_1st, T_spins_1st] = cross_time.phase_integration(
        ctime, cross_times_1st_1, cross_times_1st_1_mids, w_syn_1st_1, T_spins_1st_1,
        cross_times_1st_2, cross_times_1st_2_mids, w_syn_1st_2, T_spins_1st_2,
        cross_times_1st_3, w_syn_1st_3, T_spins_1st_3,
    )    
    logger.debug(f"[1.1] phase angle is done. ")

    """
        # 1.2 IGRF coorindate transformation: gei -> dmxl -> smxl -> fgm
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
        fgs_igrf_smxl_1st_x, fgs_igrf_smxl_1st_y, fgs_igrf_smxl_1st_z] = coordinate.dmxl2smxl(
            fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, phi_1st
    )

    # B igrf rotate from smxl to fgm
    [
        fgs_igrf_fgm_1st_x, fgs_igrf_fgm_1st_y, fgs_igrf_fgm_1st_z] = coordinate.smxl2fgm(
            fgs_igrf_smxl_1st_x, fgs_igrf_smxl_1st_y, fgs_igrf_smxl_1st_z
    )
    logger.debug(f"[1.2] igrf rotate gei -> dmxl -> smxl -> fgm. ")

    """
        # 1.3 use igrf to calibrate fgs data
    """
    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, [fgs_ful_fgm_1st_x, fgs_igrf_fgm_1st_x], [fgs_ful_fgm_1st_y, fgs_igrf_fgm_1st_y], 
    #        [fgs_ful_fgm_1st_z, fgs_igrf_fgm_1st_z], plot3 = True, title="ful_igrf_dmxl_before1stcali", xlimt = [0, 50])       

    # 1st calibration of B in dmxl 
    [
        fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z, B_parameter] = calibration.calib_leastsquare(
        fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, fgs_igrf_fgm_1st_x, fgs_igrf_fgm_1st_y, fgs_igrf_fgm_1st_z
    )
    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, [fgs_ful_fgm_2nd_x, fgs_igrf_fgm_1st_x], [fgs_ful_fgm_2nd_y, fgs_igrf_fgm_1st_y], 
    #        [fgs_ful_fgm_2nd_z, fgs_igrf_fgm_1st_z], plot3 = True, title="ful_igrf_dmxl_after1stcali", xlimt = [0, 50]) 
  
    #if parameter.makeplot == True :
    #    Bplot.B_ctime_plot(
    #        ctime, fgs_ful_fgm_2nd_x - fgs_igrf_fgm_1st_x, fgs_ful_fgm_2nd_y - fgs_igrf_fgm_1st_y, 
    #        fgs_ful_fgm_2nd_z - fgs_igrf_fgm_1st_z, xlimt = [0,50],
    #       title="fgs_res_fgm after 1st run")

    logger.debug(f"[1.3] first calibration done in fgm coordinate. ")

    """
        1.4 calib - data cross time determination
    """
    [
        cross_times_2nd_1, cross_times_2nd_1_mids, 
        T_spins_2nd_1, w_syn_2nd_1] = cross_time.cross_time_stage_1(
        ctime, fgs_ful_fgm_2nd_z,
    )
    #if parameter.makeplot == True:
    #    Bplot.B_ctime_plot(
    #    ctime, fgs_ful_fgm_x, fgs_ful_fgm_y, 
    #    fgs_ful_fgm_z, title="1cross_time_ful_fgm", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-10, ctime[ctime_idx[0]]+10],
    #    ctime_idx_time = ctime[ctime_idx[0]], cross_times = cross_times_calib_1_select,
    #    )

    [
        cross_times_2nd_2, cross_times_2nd_2_mids, 
        T_spins_2nd_2, w_syn_2nd_2] = cross_time.cross_time_stage_2(
        ctime, fgs_ful_fgm_2nd_z, cross_times_2nd_1, T_spins_2nd_1,
    )

    # if parameter.makeplot == True:
    #     Bplot.B_ctime_plot(
    #     ctime, fgs_ful_fgm_x, fgs_ful_fgm_y, 
    #     fgs_ful_fgm_z, title="2cross_time_ful_fgm", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-20, ctime[ctime_idx[0]]+20],
    #     ctime_idx_time = ctime[ctime_idx], cross_times = cross_times_calib_2_select,
    #     )

    [
        cross_times_2nd_3, T_spins_2nd_3, w_syn_2nd_3] = cross_time.cross_time_stage_3(
            ctime, fgs_ful_fgm_2nd_z, cross_times_2nd_2, T_spins_2nd_2, 
            ctime_idx = ctime_idx, ctime_idx_flag = ctime_idx_flag, ctime_idx_timediff = ctime_idx_timediff
    )
    logger.debug(f"[1.4] zero crossing calib 1-3 is done. ")

    # if parameter.makeplot == True:
    #     Bplot.B_ctime_plot(
    #         ctime, fgs_ful_fgm_x, fgs_ful_fgm_y, 
    #         fgs_ful_fgm_z, title="3cross_time_ful_fgm", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-10, ctime[ctime_idx[0]]+10],
    #         ctime_idx_time = ctime[ctime_idx], cross_times = cross_times_calib_3_select,
    #    )
    """
        1.5 calib - phase angle integration
    """
    # Remember that the stage 1 and 2 sample angular velocities at mid points of zero-crossings
    [
        phi_2nd, cross_times_2nd, w_syn_2nd, T_spins_2nd] = cross_time.phase_integration(
        ctime, cross_times_2nd_1, cross_times_2nd_1_mids, w_syn_2nd_1, T_spins_2nd_1,
        cross_times_2nd_2, cross_times_2nd_2_mids, w_syn_2nd_2, T_spins_2nd_2,
        cross_times_2nd_3, w_syn_2nd_3, T_spins_2nd_3,
    )
    logger.debug(f"[1.5] phi angle calib is done. ")

    #if parameter.makeplot == True and len(ctime_idx) != 0:
    #    Bplot.phase_plot(
    #        ctime, phi_calib, cross_times_calib, datestr = datestr, 
    #        xlimt = [ctime[ctime_idx[3]]-20, ctime[ctime_idx[3]]+20], ctime_idx = ctime_idx
    #        )

    """
        1.6 IGRF coorindate transformation : gei -> dmxl -> smxl -> fgm : 2nd calibration
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
        fgs_igrf_smxl_2nd_x, fgs_igrf_smxl_2nd_y, fgs_igrf_smxl_2nd_z] = coordinate.dmxl2smxl(
            fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, phi_2nd
    )

    # B igrf rotate from smxl to fgm
    [
        fgs_igrf_fgm_2nd_x, fgs_igrf_fgm_2nd_y, fgs_igrf_fgm_2nd_z] = coordinate.smxl2fgm(
            fgs_igrf_smxl_2nd_x, fgs_igrf_smxl_2nd_y, fgs_igrf_smxl_2nd_z
    )
    logger.debug(f"[1.6] igrf rotate gei -> dmxl -> smxl -> fgm. ")

    """
        # 1.7 use igrf to calibrate fgs data: 2nd calibration
    """
    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, [fgs_ful_fgm_1st_x, fgs_igrf_fgm_1st_x], [fgs_ful_fgm_1st_y, fgs_igrf_fgm_1st_y], 
    #        [fgs_ful_fgm_1st_z, fgs_igrf_fgm_1st_z], plot3 = True, title="ful_igrf_dmxl_before1stcali", xlimt = [0, 50])       

    # 1st calibration of B in dmxl 
    [
        fgs_ful_fgm_3rd_x, fgs_ful_fgm_3rd_y, fgs_ful_fgm_3rd_z, B_parameter] = calibration.calib_leastsquare(
        fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z, 
        fgs_igrf_fgm_2nd_x, fgs_igrf_fgm_2nd_y, fgs_igrf_fgm_2nd_z, init = B_parameter
    )
    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, [fgs_ful_fgm_2nd_x, fgs_igrf_fgm_1st_x], [fgs_ful_fgm_2nd_y, fgs_igrf_fgm_1st_y], 
    #        [fgs_ful_fgm_2nd_z, fgs_igrf_fgm_1st_z], plot3 = True, title="ful_igrf_dmxl_after1stcali", xlimt = [0, 50]) 
  
    #if parameter.makeplot == True :
    #    Bplot.B_ctime_plot(
    #        ctime, fgs_ful_fgm_2nd_x - fgs_igrf_fgm_1st_x, fgs_ful_fgm_2nd_y - fgs_igrf_fgm_1st_y, 
    #        fgs_ful_fgm_2nd_z - fgs_igrf_fgm_1st_z, xlimt = [0,50],
    #       title="fgs_res_fgm after 1st run")

    logger.debug(f"[1.7] second calibration done in fgm coordinate. ")

    """
        1.8 data cross time determination: 3nd calibration
    """
    [
        cross_times_3rd_1, cross_times_3rd_1_mids, 
        T_spins_3rd_1, w_syn_3rd_1] = cross_time.cross_time_stage_1(
        ctime, fgs_ful_fgm_3rd_z,
    )
    #if parameter.makeplot == True:
    #    Bplot.B_ctime_plot(
    #    ctime, fgs_ful_fgm_x, fgs_ful_fgm_y, 
    #    fgs_ful_fgm_z, title="1cross_time_ful_fgm", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-10, ctime[ctime_idx[0]]+10],
    #    ctime_idx_time = ctime[ctime_idx[0]], cross_times = cross_times_calib_1_select,
    #    )

    [
        cross_times_3rd_2, cross_times_3rd_2_mids, 
        T_spins_3rd_2, w_syn_3rd_2] = cross_time.cross_time_stage_2(
        ctime, fgs_ful_fgm_3rd_z, cross_times_3rd_1, T_spins_3rd_1,
    )

    # if parameter.makeplot == True:
    #     Bplot.B_ctime_plot(
    #     ctime, fgs_ful_fgm_x, fgs_ful_fgm_y, 
    #     fgs_ful_fgm_z, title="2cross_time_ful_fgm", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-20, ctime[ctime_idx[0]]+20],
    #     ctime_idx_time = ctime[ctime_idx], cross_times = cross_times_calib_2_select,
    #     )

    [
        cross_times_3rd_3, T_spins_3rd_3, w_syn_3rd_3] = cross_time.cross_time_stage_3(
            ctime, fgs_ful_fgm_3rd_z, cross_times_3rd_2, T_spins_3rd_2, 
            ctime_idx = ctime_idx, ctime_idx_flag = ctime_idx_flag, ctime_idx_timediff = ctime_idx_timediff
    )
    logger.debug(f"[1.8] zero crossing calib 1-3 is done. ")

    # if parameter.makeplot == True:
    #     Bplot.B_ctime_plot(
    #         ctime, fgs_ful_fgm_x, fgs_ful_fgm_y, 
    #         fgs_ful_fgm_z, title="3cross_time_ful_fgm", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-10, ctime[ctime_idx[0]]+10],
    #         ctime_idx_time = ctime[ctime_idx], cross_times = cross_times_calib_3_select,
    #    )
    """
        1.9 calib - phase angle integration
    """
    # Remember that the stage 1 and 2 sample angular velocities at mid points of zero-crossings
    [
        phi_3rd, cross_times_3rd, w_syn_3rd, T_spins_3rd] = cross_time.phase_integration(
        ctime, cross_times_3rd_1, cross_times_3rd_1_mids, w_syn_3rd_1, T_spins_3rd_1,
        cross_times_3rd_2, cross_times_3rd_2_mids, w_syn_3rd_2, T_spins_3rd_2,
        cross_times_3rd_3, w_syn_3rd_3, T_spins_3rd_3,
    )
    logger.debug(f"[1.9] phi angle calib is done. ")

    #if parameter.makeplot == True and len(ctime_idx) != 0:
    #    Bplot.phase_plot(
    #        ctime, phi_calib, cross_times_calib, datestr = datestr, 
    #        xlimt = [ctime[ctime_idx[3]]-20, ctime[ctime_idx[3]]+20], ctime_idx = ctime_idx
    #        )

    """
        1.10 IGRF coorindate transformation : gei -> dmxl -> smxl -> fgm : 3nd calibration
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
        fgs_igrf_smxl_3rd_x, fgs_igrf_smxl_3rd_y, fgs_igrf_smxl_3rd_z] = coordinate.dmxl2smxl(
            fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, phi_3rd
    )

    # B igrf rotate from smxl to fgm
    [
        fgs_igrf_fgm_3rd_x, fgs_igrf_fgm_3rd_y, fgs_igrf_fgm_3rd_z] = coordinate.smxl2fgm(
            fgs_igrf_smxl_3rd_x, fgs_igrf_smxl_3rd_y, fgs_igrf_smxl_3rd_z
    )
    logger.debug(f"[1.10] igrf rotate gei -> dmxl -> smxl -> fgm. ")

    """
        # 1.11 use igrf to calibrate fgs data: 3rd calibration
    """
    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, [fgs_ful_fgm_1st_x, fgs_igrf_fgm_1st_x], [fgs_ful_fgm_1st_y, fgs_igrf_fgm_1st_y], 
    #        [fgs_ful_fgm_1st_z, fgs_igrf_fgm_1st_z], plot3 = True, title="ful_igrf_dmxl_before1stcali", xlimt = [0, 50])       

    # 1st calibration of B in dmxl 
    [
        fgs_ful_fgm_4th_x, fgs_ful_fgm_4th_y, fgs_ful_fgm_4th_z, B_parameter] = calibration.calib_leastsquare(
        fgs_ful_fgm_3rd_x, fgs_ful_fgm_3rd_y, fgs_ful_fgm_3rd_z, 
        fgs_igrf_fgm_3rd_x, fgs_igrf_fgm_3rd_y, fgs_igrf_fgm_3rd_z, init = B_parameter
    )
    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, [fgs_ful_fgm_2nd_x, fgs_igrf_fgm_1st_x], [fgs_ful_fgm_2nd_y, fgs_igrf_fgm_1st_y], 
    #        [fgs_ful_fgm_2nd_z, fgs_igrf_fgm_1st_z], plot3 = True, title="ful_igrf_dmxl_after1stcali", xlimt = [0, 50]) 
  
    #if parameter.makeplot == True :
    #    Bplot.B_ctime_plot(
    #        ctime, fgs_ful_fgm_2nd_x - fgs_igrf_fgm_1st_x, fgs_ful_fgm_2nd_y - fgs_igrf_fgm_1st_y, 
    #        fgs_ful_fgm_2nd_z - fgs_igrf_fgm_1st_z, xlimt = [0,50],
    #       title="fgs_res_fgm after 1st run")

    logger.debug(f"[1.11] third calibration done in fgm coordinate. ")

    """
        1.12 rotate fgs data from fgm coordinate to dmxl
    """
    # B full rotate from fgm to smxl
    [
        fgs_ful_smxl_x, fgs_ful_smxl_y, fgs_ful_smxl_z] = coordinate.fgm2smxl(
            fgs_ful_fgm_4th_x, fgs_ful_fgm_4th_y, fgs_ful_fgm_4th_z
    )
    # if parameter.makeplot == True:
    #     Bplot.B_ctime_plot(
    #     ctime, fgs_ful_smxl_x, fgs_ful_smxl_y, 
    #     fgs_ful_smxl_z, title="3cross_time_ful_smxl", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-20, ctime[ctime_idx[0]]+20],
    #     ctime_idx_time = ctime[ctime_idx], cross_times = cross_times_calib_3_select,
    #     )

    # B full rotate from smxl to dmxl
    [
        fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z] = coordinate.smxl2dmxl(
            fgs_ful_smxl_x, fgs_ful_smxl_y, fgs_ful_smxl_z, phi_3rd
    )

    if parameter.makeplot == True :
        Bplot.B_ctime_plot(ctime, [fgs_ful_dmxl_x, fgs_igrf_dmxl_x], [fgs_ful_dmxl_y, fgs_igrf_dmxl_y], 
            [fgs_ful_dmxl_z, fgs_igrf_dmxl_z], plot3 = True, title="ful_igrf_dmxl_after3rdcali") 

    breakpoint

    logger.debug(f"[1.12] final B full field rotate from fgm to smxl to dmxl. ")
    cross_times = cross_times_3rd
    w_syn = w_syn_3rd
    T_spins = T_spins_3rd
    # if parameter.makeplot == True:
    #     for i in ctime_idx:
    #         Bplot.B_ctime_plot(
    #         ctime, fgs_ful_dmxl_x, fgs_ful_dmxl_y, 
    #         fgs_ful_dmxl_z, title=f"3cross_time_ful_dmxl_{i}", datestr = datestr, xlimt = [ctime[i]-20, ctime[i]+20],
    #         ctime_idx_time = ctime[ctime_idx], cross_times = cross_times_calib_3_select,
    #         )

    # B full rotate from dmxl to gei
    [
        fgs_ful_gei_x, fgs_ful_gei_y, fgs_ful_gei_z] = coordinate.dmxl2gei(
            fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, DMXL_2_GEI
    )


    #if parameter.makeplot == True: 
    #    Bplot.B_ctime_plot(ctime, fgs_res_dmxl_x, fgs_res_dmxl_y, fgs_res_dmxl_z, scatter = True) 

    #if parameter.makeplot == True and len(ctime_idx) != 0 :
    #    Bplot.B_ctime_plot(
    #    ctime, fgs_res_dmxl_x, fgs_res_dmxl_y, 
    #    fgs_res_dmxl_z, title="res_dmxl", datestr = datestr, xlimt = [ctime[ctime_idx[0]]-10, ctime[ctime_idx[0]]+10],
    #    ctime_idx_time = ctime[ctime_idx[0]], cross_times = cross_times_calib, scatter = True
    #    )

    return [
        cross_times, w_syn, T_spins,
        fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, 
        fgs_ful_gei_x, fgs_ful_gei_y, fgs_ful_gei_z,
        fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z,
        ]