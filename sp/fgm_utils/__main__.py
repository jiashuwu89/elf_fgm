import datetime as dt
import pandas as pd

from . import fgm_fsp_calib, get_cdf, get_relevant_state_data


if __name__ == "__main__":

    mission = "ela"

    starttime_str = "2022-01-12 15:45:51"
    endtime_str = "2022-01-12 15:52:04"

    starttime = dt.datetime.strptime(starttime_str, "%Y-%m-%d %H:%M:%S")
    endtime = dt.datetime.strptime(endtime_str, "%Y-%m-%d %H:%M:%S")

    sta_cdfpath = "fgm_utils/test/ela_l1_state_defn_20220112_v01.cdf"

    fgm_cdfpath = "fgm_utils/test/ela_l1_fgs_20220112_v01.cdf"
    fgm_cdfdata = pd.DataFrame(get_cdf(fgm_cdfpath, vars=["ela_fgs_time", "ela_fgs"]))

    att_cdfdata, pos_cdfdata = get_relevant_state_data(sta_cdfpath, mission, starttime, endtime)

    #starttime_str = "2022-01-14 15:45:50"
    #endtime_str = "2022-01-14 15:52:04"

    #sta_cdfpath = "test/ela_l1_state_defn_20220114_v01.cdf"
    #fgm_cdfpath = "test/ela_l1_fgs_20220114_v01.cdf"

    #starttime_str = "2022-06-27 08:52:34"
    #endtime_str = "2022-06-27 08:58:47"

    #sta_cdfpath = "test/ela_l1_state_defn_20220627_v01.cdf"
    #fgm_cdfpath = "test/ela_l1_fgs_20220627_v01.cdf"

    #starttime_str = "2022-01-25 16:28:23"
    #endtime_str = "2022-01-25 16:34:35"

    #sta_cdfpath = "test/ela_l1_state_defn_20220125_v01.cdf"
    #fgm_cdfpath = "test/ela_l1_fgs_20220125_v01.cdf"

    #starttime_str = "2022-07-06/09:45:53"
    #endtime_str = "2022-07-06/09:52:11"

    #sta_cdfpath = "test/ela_l1_state_defn_20220706_v01.cdf"
    #fgm_cdfpath = "test/ela_l1_fgs_20220706_v01.cdf"

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
    ] = fgm_fsp_calib(mission, starttime, endtime, fgm_cdfdata, att_cdfdata, pos_cdfdata)
