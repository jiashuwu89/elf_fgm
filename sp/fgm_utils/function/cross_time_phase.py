from . import Bplot, cross_time
from .ctime_spike_80 import find_closest
import numpy as np

def cross_times_phase(ctime, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, cross_times_1st_3, T_spins_1st_3):

    # zero crossing for x
    [
        cross_times_1st_1_x, cross_times_1st_1_mids_x, 
        T_spins_1st_1_x, w_syn_1st_1_x] = cross_time.cross_time_stage_1(
        ctime, fgs_ful_fgm_1st_x,
    )
    [
        cross_times_1st_2_x, cross_times_1st_2_mids_x, 
        T_spins_1st_2_x, w_syn_1st_2_x] = cross_time.cross_time_stage_2(
        ctime, fgs_ful_fgm_1st_x, cross_times_1st_1_x, T_spins_1st_1_x,
    )
    [
        cross_times_1st_3_x, T_spins_1st_3_x, w_syn_1st_3_x] = cross_time.cross_time_stage_3(
            ctime, fgs_ful_fgm_1st_x, cross_times_1st_2_x, T_spins_1st_2_x
    )
    # zero crossing for y
    [
        cross_times_1st_1_y, cross_times_1st_1_mids_y, 
        T_spins_1st_1_y, w_syn_1st_1_y] = cross_time.cross_time_stage_1(
        ctime, fgs_ful_fgm_1st_y,
    )
    [
        cross_times_1st_2_y, cross_times_1st_2_mids_y, 
        T_spins_1st_2_y, w_syn_1st_2_y] = cross_time.cross_time_stage_2(
        ctime, fgs_ful_fgm_1st_y, cross_times_1st_1_y, T_spins_1st_1_y,
    )
    [
        cross_times_1st_3_y, T_spins_1st_3_y, w_syn_1st_3_y] = cross_time.cross_time_stage_3(
            ctime, fgs_ful_fgm_1st_y, cross_times_1st_2_y, T_spins_1st_2_y
    )

    Bplot.B_ctime_plot(ctime, fgs_ful_fgm_1st_x, fgs_ful_fgm_1st_y, fgs_ful_fgm_1st_z, xlimt = [0, 5], cross_times = [cross_times_1st_3_x, cross_times_1st_3_y, cross_times_1st_3])
    cross_times_1st_3_xdiff, cross_times_1st_3_ydiff, cross_times_1st_3_select = [], [], []
    for cross_times_idx, cross_times_i in enumerate(cross_times_1st_3):
        cross_times_x_idx, cross_times_x = find_closest(cross_times_1st_3_x, cross_times_i) # cross time in x and z has difference smaller than 1s, pair 
        cross_times_y_idx, cross_times_y = find_closest(cross_times_1st_3_y, cross_times_i)
        if np.abs(cross_times_x - cross_times_i) < 1 and np.abs(cross_times_x - cross_times_y) < 1:
            cross_times_1st_3_xdiff.append(cross_times_x - cross_times_i - 0.25*T_spins_1st_3[cross_times_idx])
            cross_times_1st_3_ydiff.append(cross_times_y - cross_times_i - 0.25*T_spins_1st_3[cross_times_idx])    
            cross_times_1st_3_select.append(cross_times_i)
        
    Bplot.B_ctime_plot_single(cross_times_1st_3_select, [cross_times_1st_3_xdiff, cross_times_1st_3_ydiff], scatter = True, 
        legend=[f'$\Delta$Cross Times XZ', '$\Delta$Cross Times YZ'], title = 'time_diff_xyz', ylabel='time(s)')
    breakpoint()