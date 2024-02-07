"""This code is modified according to step1.py. delete 2nd, 3rd calibration. add wfit iteration
"""
from .. import parameter
from . import cross_time, coordinate, calibration, Bplot
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import curve_fit


func = calibration.quad_fit

def cal_res(phi, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, 
            f, fgs_ful_fgm_x, fgs_ful_fgm_y, fgs_ful_fgm_z, saveresult=False, **kwargs):

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

    result = {'res': np.median(res)}

    if saveresult == True:
        Bplot.B_ctime_plot(kwargs['ctime'], fgs_ful_fgm_x_calib - fgs_igrf_fgm_1st_x, fgs_ful_fgm_y_calib - fgs_igrf_fgm_1st_y, 
            fgs_ful_fgm_z_calib - fgs_igrf_fgm_1st_z, title=f"res_fgm_{kwargs['epoch']}") 
        result['fgs_ful_fgm_x_calib'] = fgs_ful_fgm_x_calib
        result['fgs_ful_fgm_y_calib'] = fgs_ful_fgm_y_calib
        result['fgs_ful_fgm_z_calib'] = fgs_ful_fgm_z_calib
        result['fgs_igrf_fgm_1st_x'] = fgs_igrf_fgm_1st_x      
        result['fgs_igrf_fgm_1st_y'] = fgs_igrf_fgm_1st_y      
        result['fgs_igrf_fgm_1st_z'] = fgs_igrf_fgm_1st_z      
    return result


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

    if parameter.wfit_gradient_figure == True:
        gradient_figure(cross_time0, wfit_m, wfit_c, f,
                    ctime, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, 
                    fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)
        breakpoint()

    phi = phi_2nd
    delta_wfit_c = 1e-5
    lr_wfit_c = 5e-10
    delta_wfit_m = 1e-2
    lr_wfit_m = 1e-17
    delta_ct = 1e-3
    lr_ct = 0.5
    tolerance = 1e-2

    ############################
    cross_time0 = 245.51 + 0.2
    wfit_m = wfit_m - 1e-7
    breakpoint()
    ###########################
    max_iter = parameter.wfit_run_maxiter
    for epoch in range(max_iter):
        phi = phase_integration_update(ctime, cross_time0, wfit_m, wfit_c)
        kwargs = {'ctime': ctime, 'epoch': epoch}
        result = cal_res(
            phi, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, 
            fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z, saveresult=True, **kwargs)
        res = result['res']
        fgs_ful_fgm_2nd_x = result['fgs_ful_fgm_x_calib']
        fgs_ful_fgm_2nd_y = result['fgs_ful_fgm_y_calib']
        fgs_ful_fgm_2nd_z = result['fgs_ful_fgm_z_calib']
        fgs_igrf_fgm_1st_x = result['fgs_igrf_fgm_1st_x']
        fgs_igrf_fgm_1st_y = result['fgs_igrf_fgm_1st_y']
        fgs_igrf_fgm_1st_z = result['fgs_igrf_fgm_1st_z']
        phi_2nd = phi

   
        if epoch == 0:
            res_prev = res
        else:
            if (epoch == max_iter-1) or (np.abs(res - res_prev) < tolerance):
                break

        # calculate gradient of wfit_c
        wfit_c_left = wfit_c - delta_wfit_c*wfit_c
        wfit_c_phi_left = phase_integration_update(ctime, cross_time0, wfit_m, wfit_c_left)
        wfit_c_res_left = cal_res(wfit_c_phi_left, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)['res']
        
        wfit_c_right = wfit_c + delta_wfit_c*wfit_c
        wfit_c_phi_right = phase_integration_update(ctime, cross_time0, wfit_m, wfit_c_right)
        wfit_c_res_right = cal_res(wfit_c_phi_right, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)['res']
        wfit_c_grad = (wfit_c_res_right - wfit_c_res_left)/(wfit_c_right - wfit_c_left)

        # calculate gradient of wfit_m
        wfit_m_left = wfit_m - delta_wfit_m*wfit_m
        wfit_m_phi_left = phase_integration_update(ctime, cross_time0, wfit_m_left, wfit_c)
        wfit_m_res_left =  cal_res(wfit_m_phi_left, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)['res']
        
        wfit_m_right = wfit_m + delta_wfit_m*wfit_m
        wfit_m_phi_right = phase_integration_update(ctime, cross_time0, wfit_m_right, wfit_c)
        wfit_m_res_right =  cal_res(wfit_m_phi_right, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)['res']
        wfit_m_grad = (wfit_m_res_right - wfit_m_res_left)/(wfit_m_right - wfit_m_left)

        # calculate gradient of cross_time0
        cross_time0_left = cross_time0 - delta_ct*cross_time0
        cross_time0_phi_left = phase_integration_update(ctime, cross_time0_left, wfit_m, wfit_c)
        cross_time0_res_left = cal_res(cross_time0_phi_left, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)['res']

        cross_time0_right = cross_time0 + delta_ct*cross_time0
        cross_time0_phi_right = phase_integration_update(ctime, cross_time0_right, wfit_m, wfit_c)
        cross_time0_res_right = cal_res(cross_time0_phi_right, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)['res']
        cross_time0_grad = (cross_time0_res_right - cross_time0_res_left)/(cross_time0_right - cross_time0_left)

        print(f"=================={epoch} run================")
        print(f"res:{res}")
        print(f"wfit_c_res_left: {wfit_c_res_left}")
        print(f"wfit_c_res_right: {wfit_c_res_right}")
        print(f"wfit_m_res_left: {wfit_m_res_left}")
        print(f"wfit_m_res_right: {wfit_m_res_right}")
        print(f"cross_time0_res_left: {cross_time0_res_left}")
        print(f"cross_time0_res_right: {cross_time0_res_right}")
        
        # update parameter
        wfit_c = wfit_c - lr_wfit_c*wfit_c_grad
        wfit_m = wfit_m - lr_wfit_m*wfit_m_grad
        cross_time0 = cross_time0 - lr_ct*cross_time0_grad

        print(f"wfit_c is {wfit_c}, change by {lr_wfit_c*wfit_c_grad}")
        print(f"wfit_m is {wfit_m}, change by {lr_wfit_m*wfit_m_grad}")
        print(f"cross_time0 is {cross_time0}, change by {lr_ct*cross_time0_grad}")

        res_prev = res
    
    logger.debug(f"[1.5] phi angle calib is updated. ")
    
    if parameter.makeplot == True: 
        Bplot.B_ctime_plot(ctime, [fgs_ful_fgm_2nd_x, fgs_igrf_fgm_1st_x], [fgs_ful_fgm_2nd_y, fgs_igrf_fgm_1st_y], 
            [fgs_ful_fgm_2nd_z, fgs_igrf_fgm_1st_z], plot3 = True, title="fuligrf_fgm_after1stcali")  
        Bplot.B_ctime_plot(ctime, fgs_ful_fgm_2nd_x-fgs_igrf_fgm_1st_x, fgs_ful_fgm_2nd_y-fgs_igrf_fgm_1st_y, 
            fgs_ful_fgm_2nd_z-fgs_igrf_fgm_1st_z, plot3 = True, title="res_fgm_after1stcali")  
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

    if parameter.cal_dmxl == True:
        [
            fgs_ful_dmxl_3rd_x, fgs_ful_dmxl_3rd_y, fgs_ful_dmxl_3rd_z, B_parameter] = calibration.calib_leastsquare(
            fgs_ful_dmxl_2nd_x, fgs_ful_dmxl_2nd_y, fgs_ful_dmxl_2nd_z, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, init = B_parameter)
        fgs_ful_dmxl_x = fgs_ful_dmxl_3rd_x
        fgs_ful_dmxl_y = fgs_ful_dmxl_3rd_y
        fgs_ful_dmxl_z = fgs_ful_dmxl_3rd_z
    else:
        fgs_ful_dmxl_x = fgs_ful_dmxl_2nd_x
        fgs_ful_dmxl_y = fgs_ful_dmxl_2nd_y
        fgs_ful_dmxl_z = fgs_ful_dmxl_2nd_z

 
    cross_times = cross_times_2nd
    w_syn = w_syn_2nd
    T_spins = T_spins_2nd
    DMXL_2_GEI_fsp = cross_time.fsp_matrix(ctime, cross_times, T_spins, DMXL_2_GEI)


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


def gradient_figure(cross_time0, wfit_m, wfit_c, f,
                    ctime, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, 
                    fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z):
    para = parameter.wfit_gradient_choice_lst[parameter.wfit_gradient_choice]
    ranges = {
        'cross_time0': np.arange(-1.4, 1.4, 0.01),
        'wfit_m': np.arange(-0.5e-5, 0.5e-5, 1e-7),
        'wfit_c': np.arange(-0.5e-5, 0.5e-5, 1e-7),
        't': np.arange(-2e-2, 2e-2, 5e-4),
    }
    plot_config = {
        'cross_time0': {
            'title': f'cross_time0 = {cross_time0:.2f}', 
            'xlabel': 'relative seconds', 
            'ylabel': 'root mean squared residuals', 
            'fname': 'loop_crosstime'
            },
        'wfit_m': {
            'title': f"dw = {wfit_m:.4e}", 
            'xlabel': 'relative change in dw', 
            'ylabel': 'root mean squared residuals', 
            'fname': 'loop_dw'
            },
        'wfit_c': {
            'title': f"w0 = {wfit_c:.4f}", 
            'xlabel': 'w0', 
            'ylabel': 'root mean squared residuals', 
            'fname': 'loop_w0', 
            'xlimt': [-0.002, 0.002], 
            'ylimt': [0, 150]
            },
        't': {
            'title': f"t = 0.1s",
            'xlabel': 'delta t',
            'ylabel': 'root mean squared residuals',
            'fname': 'loop_t',
        },

    }

    def integrate(dx):
        if para == 'cross_time0':
            phi = phase_integration_update(ctime, cross_time0 + dx, wfit_m, wfit_c)
        elif para == 'wfit_m':
            phi = phase_integration_update(ctime, cross_time0, wfit_m + dx, wfit_c)
        elif para == 'wfit_c':
            phi = phase_integration_update(ctime, cross_time0, wfit_m, wfit_c + dx)
        elif para == 't':
            #ctime_jitter =  ctime + np.random.uniform(-np.abs(dx), np.abs(dx), ctime.shape)
            ctime_jitter =  ctime + np.random.normal(0, np.abs(dx), ctime.shape)
            phi = phase_integration_update(ctime_jitter, cross_time0, wfit_m, wfit_c)
        result = cal_res(phi, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, 
             fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)
        return result['res']

    results = [integrate(dx) for dx in ranges[para]]

    Bplot.B_ctime_plot_single(ranges[para], results, **plot_config[para])



    # # cross time loop
    # cross_time0_range = np.arange(-1.4, 1.4, 0.01)
    # cross_time0_result = []
    # for i, dt in enumerate(cross_time0_range):
    #     phi = phase_integration_update(ctime, cross_time0 + dt, wfit_m, wfit_c)
    #     result = cal_res(
    #         phi, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, 
    #         fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)
    #     cross_time0_result.append(result['res'])
    # Bplot.B_ctime_plot_single(cross_time0_range, cross_time0_result, title='cross_time0 = {:.2f}'.format(cross_time0), xlabel='relative seconds', ylabel='root mean squared residuals', fname='loop_crosstime')
    # breakpoint()

    # wfit m loop
    # wfit_m_range = np.arange(-0.5e-5, 0.5e-5, 1e-7)
    # wfit_m_result = []
    # for i, dw in enumerate(wfit_m_range):
    #     phi = phase_integration_update(ctime, cross_time0, wfit_m + dw, wfit_c)
    #     result = cal_res(
    #         phi, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, 
    #         fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)
    #     wfit_m_result.append(result['res'])
    # Bplot.B_ctime_plot_single(wfit_m_range, wfit_m_result, title="dw = {:.4e}".format(wfit_m), xlabel='relative change in dw', ylabel='root mean squared residuals', fname='loop_dw')
    # breakpoint()

    # wfit c loop
    # wfit_c_range = np.ara(, 0.005)
    # wfit_c_result = []
    # for i, dw in enumerate():
    #     phi = phase_integra wfit_c + dw)
    #     result = cal_res(
    #         phi, fgs_igrf_dmxl_x, fgs_igrf_dmxl_y, fgs_igrf_dmxl_z, f, 
    #         fgs_ful_fgm_2nd_x, fgs_ful_fgm_2nd_y, fgs_ful_fgm_2nd_z)
    #     wfit_c_result.append(result['res'])
    # Bplot.B_ctime_plot_single(wfit_c_range, wfit_c_result, title="w0 = {:.4f}".format(wfit_c), xlabel='w0', ylabel='root mean squared residuals', fname='loop_w0', xlimt=[-0.002, 0.002], ylimt=[0, 150])
    breakpoint()

    