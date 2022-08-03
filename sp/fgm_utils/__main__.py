from . import fgm_fsp_calib


if __name__ == "__main__":

    #starttime_str = "2022-01-12 04:54:09"
    #endtime_str = "2022-01-12 05:00:21"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220112_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220112_v01.cdf"

    #starttime_str = "2022-01-12 15:45:51"
    #endtime_str = "2022-01-12 15:52:04"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220112_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220112_v01.cdf"

    #starttime_str = "2022-01-14 15:45:50"
    #endtime_str = "2022-01-14 15:52:04"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220114_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220114_v01.cdf"

    #starttime_str = "2022-01-14 23:26:43"
    #endtime_str = "2022-01-14 23:32:54"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220114_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220114_v01.cdf"

    #starttime_str = "2022-06-27 08:52:34"
    #endtime_str = "2022-06-27 08:58:47"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220627_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220627_v01.cdf"

    #starttime_str = "2022-01-25 16:28:23"
    #endtime_str = "2022-01-25 16:34:35"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220125_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220125_v01.cdf"

    #starttime_str = "2022-07-06/09:45:53"
    #endtime_str = "2022-07-06/09:52:11"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220706_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220706_v01.cdf"

    #starttime_str = "2022-06-25 06:27:08"	
    #endtime_str = "2022-06-25/06:33:20"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220625_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220625_v01.cdf"

    # TODO:has spike
    #starttime_str = "2022-06-25 09:29:30"
    #endtime_str = "2022-06-25 09:35:44"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220625_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220625_v01.cdf"

    # TODO: has unipolar spike in the middle
    #starttime_str = "2022-06-25 18:36:12"	
    #endtime_str = "2022-06-25/18:42:25"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220625_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220625_v01.cdf"

    # TODO: no spikes but first point rogue; example asked by V
    #starttime_str = "2022-06-23 02:27:59"
    #endtime_str = "2022-06-23 02:34:12"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220623_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220623_v01.cdf"

    # TODO: two unipolar spikes, example asked by V
    starttime_str = "2022-06-23 04:00:07"
    endtime_str = "2022-06-23 04:06:19"
    sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220623_v01.cdf"
    fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220623_v01.cdf"

    # TODO: one spikes
    #starttime_str = "2022-06-23 08:34:48"
    #endtime_str = "2022-06-23 08:41:01"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220623_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220623_v01.cdf"

    #
    #starttime_str = "2022-06-23 11:36:39"
    #endtime_str = "2022-06-23 11:42:52"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220623_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220623_v01.cdf"

    #starttime_str = "2022-06-17 01:08:33"
    #endtime_str = "2022-06-17 01:14:47"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

    #starttime_str = "2022-06-17 04:13:18"
    #endtime_str = "2022-06-17 04:19:30"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

    #starttime_str = "2022-06-17 08:48:13"
    #endtime_str = "2022-06-17 08:54:24"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

    # has gap
    #starttime_str = "2022-06-17 11:50:13"
    #endtime_str = "2022-06-17 11:56:26"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

    #starttime_str = "2022-06-17 14:52:03"
    #endtime_str = "2022-06-17 14:58:16"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

    #starttime_str = "2022-06-17 20:59:35"
    #endtime_str = "2022-06-17 21:05:48"
    #sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220617_v01.cdf"
    #fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220617_v01.cdf"

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
    ] = fgm_fsp_calib(starttime_str, endtime_str, sta_cdfpath, fgm_cdfpath)
