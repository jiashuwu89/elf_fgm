"""This code is modified according to step1.py. delete 2nd, 3rd calibration. add wfit iteration
"""
from .. import parameter
from . import cross_time, coordinate, calibration, Bplot
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import curve_fit


func = calibration.quad_fit

def cal_res(phi, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, 
            f, fgs_ful_fgm_x, fgs_ful_fgm_y, fgs_ful_fgm_z):

    # B igrf rotate from dmxl to smxl
    [
        fgs_igrf_smxl_1st_x, fgs_igrf_smxl_1st_y, fgs_igrf_smxl_1st_z] = coordinate.dmxl2smxl(
            fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, phi
    )
   
    [
        fgs_igrf_fgm_1st_x, fgs_igrf_fgm_1st_y, fgs_igrf_fgm_1st_z] = coordinate.smxl2fgm(
            fgs_igrf_smxl_1st_x, fgs_igrf_smxl_1st_y, fgs_igrf_smxl_1st_z, f
    )

    # 1st calibration of B in dmxl 
    [
        fgs_ful_fgm_x_calib, fgs_ful_fgm_y_calib, fgs_ful_fgm_z_calib, B_parameter] = calibration.calib_leastsquare(
        fgs_ful_fgm_x, fgs_ful_fgm_y, fgs_ful_fgm_z, fgs_igrf_fgm_1st_x, fgs_igrf_fgm_1st_y, fgs_igrf_fgm_1st_z
    )

    res = np.sqrt((
        (fgs_ful_fgm_x_calib-fgs_igrf_fgm_1st_x)**2 + (fgs_ful_fgm_y_calib-fgs_igrf_fgm_1st_y)**2 + (fgs_ful_fgm_z_calib-fgs_igrf_fgm_1st_z)**2
        )/len(fgs_ful_fgm_x_calib))

    return np.median(res)


def step1(
    ctime, ctimestamp, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, 
    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, 
    att_gei_x, att_gei_y, att_gei_z,
    datestr, logger, ctime_idx, ctime_idx_time, ctime_idx_flag, ctime_idx_timediff, f):   

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
        phi_1st, cross_times_1st, w_syn_1st, T_spins_1st, cross_times_1st_fit, w_syn_1st_fit] = cross_time.phase_integration(
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

    [
        fgs_igrf_fgm_1st_x, fgs_igrf_fgm_1st_y, fgs_igrf_fgm_1st_z] = coordinate.smxl2fgm(
            fgs_igrf_smxl_1st_x, fgs_igrf_smxl_1st_y, fgs_igrf_smxl_1st_z, f
    )
    logger.debug(f"[1.2] igrf rotate gei -> dmxl -> smxl -> fgm. ")

    """
        # 1.3 use igrf to calibrate fgs data
    """          
    [
        fgs_fsp_ful_fgm_x, fgs_fsp_ful_fgm_y, fgs_fsp_ful_fgm_z] = cross_time.fsp_ful(
            ctime, cross_times_1st, T_spins_1st, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z
    )
    [
        fgs_fsp_igrf_fgm_x, fgs_fsp_igrf_fgm_y, fgs_fsp_igrf_fgm_z] = cross_time.fsp_ful(
            ctime, cross_times_1st, T_spins_1st, fgs_igrf_fgm_1st_x, fgs_igrf_fgm_1st_y, fgs_igrf_fgm_1st_z
    )
    
    # 1st calibration of B in dmxl 
    [
        fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z, B_parameter] = calibration.calib_leastsquare(
        fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, fgs_igrf_fgm_1st_x, fgs_igrf_fgm_1st_y, fgs_igrf_fgm_1st_z
    )

    # update zero cross after 1st calib
    [
        cross_times_2nd_1, cross_times_2nd_1_mids, 
        T_spins_2nd_1, w_syn_2nd_1] = cross_time.cross_time_stage_1(
        ctime, fgs_ful_fgm_2nd_z,
    )
    [
        cross_times_2nd_2, cross_times_2nd_2_mids, 
        T_spins_2nd_2, w_syn_2nd_2] = cross_time.cross_time_stage_2(
        ctime, fgs_ful_fgm_2nd_z, cross_times_2nd_1, T_spins_2nd_1,
    )
    [
        cross_times_2nd_3, T_spins_2nd_3, w_syn_2nd_3] = cross_time.cross_time_stage_3(
            ctime, fgs_ful_fgm_2nd_z, cross_times_2nd_2, T_spins_2nd_2, 
            ctime_idx = ctime_idx, ctime_idx_flag = ctime_idx_flag, ctime_idx_timediff = ctime_idx_timediff
    )
    logger.debug(f"[1.4] zero crossing calib 1-3 updated after calibration. ")

    """
        1.5 calib - phase angle integration
    """
    # Remember that the stage 1 and 2 sample angular velocities at mid points of zero-crossings
    [
        phi_2nd, cross_times_2nd, w_syn_2nd, T_spins_2nd, cross_times_2nd_fit, w_syn_2nd_fit] = cross_time.phase_integration(
        ctime, cross_times_2nd_1, cross_times_2nd_1_mids, w_syn_2nd_1, T_spins_2nd_1,
        cross_times_2nd_2, cross_times_2nd_2_mids, w_syn_2nd_2, T_spins_2nd_2,
        cross_times_2nd_3, w_syn_2nd_3, T_spins_2nd_3,
    )
    cross_time0 = cross_times_2nd[0]
    wfit_para = curve_fit(calibration.linear_fit, cross_times_2nd, w_syn_2nd_fit)
    wfit_m = wfit_para[0][0]
    wfit_c = wfit_para[0][1]
    phi = phi_2nd
    delta_step = 1e-5
    alpha = 5e-10
    tolerance = 0.1
    for i in range(50):
        phi = phase_integration_update(ctime, cross_time0, wfit_m, wfit_c)
        res = cal_res(phi, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)
        if i == 0:
            res_prev = res
        else:
            if np.abs(res - res_prev) < tolerance:
                break

        # calculate gradient of wfit_c
        wfit_c_left = wfit_c - delta_step*wfit_c
        wfit_c_phi_left = phase_integration_update(ctime, cross_time0, wfit_m, wfit_c_left)
        wfit_c_res_left = cal_res(wfit_c_phi_left, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)
        wfit_c_right = wfit_c + delta_step*wfit_c
        wfit_c_phi_right = phase_integration_update(ctime, cross_time0, wfit_m, wfit_c_right)
        wfit_c_res_right = cal_res(wfit_c_phi_right, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)
        

        print(f"=================={i} run================")
        print(f"res:{res}")
        print(f"wfit_c_res_left: {wfit_c_res_left}")
        print(f"wfit_c_res_right: {wfit_c_res_right}")

        # calculate gradient in parameter space
        wfit_c_grad = (wfit_c_res_right - wfit_c_res_left)/(wfit_c_right - wfit_c_left)

        # update parameter
        wfit_c = wfit_c - alpha*wfit_c_grad

        res_prev = res
        

    
    breakpoint()
    logger.debug(f"[1.5] phi angle calib is updated. ")
    

    """
        change to dmxl to see the results
    """
    # B full rotate from fgm to smxl
    [
        fgs_ful_smxl_2nd_x, fgs_ful_smxl_2nd_y, fgs_ful_smxl_2nd_z] = coordinate.fgm2smxl(
            fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z, f
    )
    # B full rotate from smxl to dmxl
    [
        fgs_ful_dmxl_2nd_x, fgs_ful_dmxl_2nd_y, fgs_ful_dmxl_2nd_z] = coordinate.smxl2dmxl(
            fgs_ful_smxl_2nd_x, fgs_ful_smxl_2nd_y, fgs_ful_smxl_2nd_z, phi_2nd
    )
    if parameter.makeplot == True :
        Bplot.B_ctime_plot(ctime, [fgs_ful_dmxl_2nd_x, fgs_igrf_dmxl_x], [fgs_ful_dmxl_2nd_y, fgs_igrf_dmxl_y], 
            [fgs_ful_dmxl_2nd_z, fgs_igrf_dmxl_z], title="fuligrf_dmxl_after1stcali") 
        Bplot.B_ctime_plot(ctime, fgs_ful_dmxl_2nd_x-fgs_igrf_dmxl_x, fgs_ful_dmxl_2nd_y-fgs_igrf_dmxl_y, 
            fgs_ful_dmxl_2nd_z-fgs_igrf_dmxl_z, title="res_dmxl_after1stcali") 
        
    if parameter.makeplot == True: 
        Bplot.B_ctime_plot(ctime, [fgs_ful_smxl_2nd_x, fgs_igrf_smxl_1st_x], [fgs_ful_smxl_2nd_y, fgs_igrf_smxl_1st_y], 
            [fgs_ful_smxl_2nd_z, fgs_igrf_smxl_1st_z], plot3 = True, title="fuligrf_sxml_after1stcali")
        [
            fgs_fsp_ful_smxl_2nd_x, fgs_fsp_ful_smxl_2nd_y, fgs_fsp_ful_smxl_2nd_z] = cross_time.fsp_ful(
                ctime, cross_times_1st, T_spins_1st, fgs_ful_smxl_2nd_x, fgs_ful_smxl_2nd_y, fgs_ful_smxl_2nd_z
        )
        [
            fgs_fsp_igrf_smxl_1st_x, fgs_fsp_igrf_smxl_1st_y, fgs_fsp_igrf_smxl_1st_z] = cross_time.fsp_ful(
                ctime, cross_times_1st, T_spins_1st, fgs_igrf_smxl_1st_x, fgs_igrf_smxl_1st_y, fgs_igrf_smxl_1st_z
        )
        Bplot.B_ctime_plot(cross_times_1st, [fgs_fsp_ful_smxl_2nd_x, fgs_fsp_igrf_smxl_1st_x], [fgs_fsp_ful_smxl_2nd_y, fgs_fsp_igrf_smxl_1st_y], 
            [fgs_fsp_ful_smxl_2nd_z, fgs_fsp_igrf_smxl_1st_z], plot3 = True, title="fuligrf_smxl_fsp_after1stcali")

    
 
    cross_times = cross_times_2nd
    w_syn = w_syn_2nd
    T_spins = T_spins_2nd
    fgs_ful_dmxl_x = fgs_ful_dmxl_2nd_x
    fgs_ful_dmxl_y = fgs_ful_dmxl_2nd_y
    fgs_ful_dmxl_z = fgs_ful_dmxl_2nd_z
    DMXL_2_GEI_fsp = cross_time.fsp_matrix(ctime, cross_times, T_spins, DMXL_2_GEI)

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
        cross_times, w_syn, T_spins, DMXL_2_GEI_fsp,
        fgs_ful_dmxl_x, fgs_ful_dmxl_y, fgs_ful_dmxl_z, 
        fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z,
        B_parameter,
        ]


def phase_integration_update(
    ctime, cross_time0, wfit_m, wfit_c
):
    """This code is modified according to func phase_integration in cross_time.py
    only keep integration to the first point, and method=3
    """

    delta_t = np.median(ctime[1:]-ctime[:-1])
    w_syn = calibration.linear_fit(ctime, wfit_m, wfit_c)
    
    # Use just one reference point for integration
    t0 = cross_time0
    idx0 = np.where(ctime <= t0)[0][-1]
    w_t0 = calibration.linear_fit(
        t0, wfit_m, wfit_c
    )

    phi = np.zeros(len(ctime))
    for i in range(len(phi)):

        # Do the integral based on the relative position of the time and the relevant zero-crossing
        if i < idx0:
            phi[i] = -simpson(
                w_syn[i : idx0 + 1], x=ctime[i : idx0 + 1], dx=delta_t
            )
        elif i > idx0:
            phi[i] = simpson(
                w_syn[idx0 : i + 1], x=ctime[idx0 : i + 1], dx=delta_t
            )
        else:
            phi[i] = 0

        # Correct for the zero-crossing not being in the time array
        phi[i] -= (
            0.5 * (w_t0 + w_syn[idx0]) * (t0 - ctime[idx0])
        )            

    return phi