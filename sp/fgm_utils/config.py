import numpy as np

default_params = {
    "proper_pad": False,
    "fit_running_spline": False,
    "relative_integrate": True,
    "bidirectional_integrate": False,
    "init_secs": 0,
    "funkyfgm": True,
    "f": {
    'elb': (-90-55) * np.pi / 180,
    'ela': 85 * np.pi / 180,
    },
    "eps_1": 1,
    "eps_2": 2,
    "eps_3": 2,
    "N_spins_fit": 4,
    "peak_detect": False,
    "zero_crossing_method": 3,
    "ctime_thrhld": 0.003,
    "ctime_correct_80": False,
    "ctime_correct_80_skip": False,
    "spike_find_spin_num": 4,
    "spike_fit_len_80": 28,
    "ctime_correct_100": True,
    "cross0_spike_del": True,
    "fsp_spike_del_type2": True,
    "fsp_spike_del_type3": True,
    "fsp_spike_del_type4": True,
    "R2_filter": True,
    "R2_thrhld": 0.8,
    "Spinrate_thrhld": 0.2,
    "ctime_repeat_check": True,
    "cali_2nd": False,
    "cali_3rd": False,
    "cali_4th": False,
    "del_rogue_fsp": True,
    "eps_rogue": 3,
    "prefsp_detrend": True,
    "prefsp_detrend_func": 'detrend_cube',
    "fsp_detrend": False,
    "fsp_detrend_cutoff": 6,
    "fsp_detrend_func": 'detrend_quad',
    "makeplot": False,
    "savepng": False,
    "output": False,
    "download_data": True,
    "del_time": False,
    "del_time_idxstart": [260],
    "del_time_idxend": [306],
    "gei2obw": True,
    "CrossTime_Update": False,
    "fit_bound": False,
    "f_loop": False,
    "att_rot": False,
    "att_rot_ang": 0.5,
    "att_rot_axis": 'z',
    "att_loop": False,
    "att_loop_width": 1,
    "att_loop_step": 0.5,
    "att_loop_figure": False,
    "att_csv": False,
    "mva_fgm": False,
    "mvaang_rotang": False,
    "Bpara_out": True,
    "att_determine": False,
    "att_split": False,
    "att_split_num": None,
    "att_split_idx": None,
    "att_split_detrend": False,
    "batch_run": False,
    "wfit_run": False,
    "wfit_run_maxiter": 1,
    "wfit_gradient_figure": True,
    "wfit_gradient_choice": 4,
    "cal_dmxl": False,
    "state03": False,
    "state03_plotatt": False,
    "nonlinear": True,
    "nonlinear_phi90": False,
    "skipfit": False,
    "beforecali_detrend": True,
    "shift_loop": False,
    "read_beta": False,
    }

def compare_parameters(logger):
    from . import parameter
    print(f"PARAM CHECK START:")
    for param in default_params:
        current_value = getattr(parameter, param, None)
        if current_value != default_params[param]:
            print(f"PARAM CHECK: {param} is set {current_value}!")
    print(f"PARAM CHECK END (other params are default)")