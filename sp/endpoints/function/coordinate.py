import numpy as np
import parameter as parameter

def dmxl_gei_matrix(
    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, att_gei_x, att_gei_y, att_gei_z
    ):
    """generate transformation matrix between dmxl and gei 
    """
    fgs_igrf_gei_x_hat = fgs_igrf_gei_x / np.sqrt(
        fgs_igrf_gei_x**2 + fgs_igrf_gei_y**2 + fgs_igrf_gei_z**2
    )
    fgs_igrf_gei_y_hat = fgs_igrf_gei_y / np.sqrt(
        fgs_igrf_gei_x**2 + fgs_igrf_gei_y**2 + fgs_igrf_gei_z**2
    )
    fgs_igrf_gei_z_hat = fgs_igrf_gei_z / np.sqrt(
        fgs_igrf_gei_x**2 + fgs_igrf_gei_y**2 + fgs_igrf_gei_z**2
    )

    DMXL_2_GEI = np.zeros((len(fgs_igrf_gei_x), 3, 3))
    GEI_2_DMXL = np.zeros((len(fgs_igrf_gei_x), 3, 3))
    for i in range(len(fgs_igrf_gei_x)):
        u_hat = np.array([att_gei_x[i], att_gei_y[i], att_gei_z[i]])
        b_hat = np.array(
            [fgs_igrf_gei_x_hat[i], fgs_igrf_gei_y_hat[i], fgs_igrf_gei_z_hat[i]]
        )

        DMXL_2_GEI[i, :, 0] = np.cross(b_hat, u_hat)
        DMXL_2_GEI[i, :, 1] = np.cross(u_hat, np.cross(b_hat, u_hat))
        DMXL_2_GEI[i, :, 2] = u_hat

        DMXL_2_GEI[i, :, 0] /= np.linalg.norm(DMXL_2_GEI[i, :, 0])
        DMXL_2_GEI[i, :, 1] /= np.linalg.norm(DMXL_2_GEI[i, :, 1])
        DMXL_2_GEI[i, :, 2] /= np.linalg.norm(DMXL_2_GEI[i, :, 2])

        GEI_2_DMXL[i, :, :] = np.linalg.inv(DMXL_2_GEI[i, :, :])

    return [DMXL_2_GEI, GEI_2_DMXL]



def gei2dmxl(fgs_gei_x, fgs_gei_y, fgs_gei_z, GEI_2_DMXL):
    """
        input: 
            fgs_gei_x/y/z: time series of fgs data in gei
            GEI_2_DMXL: times seris of rotation matrix
        output:
            fgs_dmxl_x/y/z: times series of fgs data in dmxl
    """
    fgs_dmxl_x = np.zeros(len(fgs_gei_x))
    fgs_dmxl_y = np.zeros(len(fgs_gei_x))
    fgs_dmxl_z = np.zeros(len(fgs_gei_x))
    for i in range(len(fgs_gei_x)):
        fgs_dmxl_x[i] = (
            GEI_2_DMXL[i, 0, 0] * fgs_gei_x[i]
            + GEI_2_DMXL[i, 0, 1] * fgs_gei_y[i]
            + GEI_2_DMXL[i, 0, 2] * fgs_gei_z[i]
        )
        fgs_dmxl_y[i] = (
            GEI_2_DMXL[i, 1, 0] * fgs_gei_x[i]
            + GEI_2_DMXL[i, 1, 1] * fgs_gei_y[i]
            + GEI_2_DMXL[i, 1, 2] * fgs_gei_z[i]
        )
        fgs_dmxl_z[i] = (
            GEI_2_DMXL[i, 2, 0] * fgs_gei_x[i]
            + GEI_2_DMXL[i, 2, 1] * fgs_gei_y[i]
            + GEI_2_DMXL[i, 2, 2] * fgs_gei_z[i]
        )

    return [fgs_dmxl_x, fgs_dmxl_y, fgs_dmxl_z]


def dmxl2gei(fgs_dmxl_x, fgs_dmxl_y, fgs_dmxl_z, DMXL_2_GEI):
    """
        input: 
            fgs_gei_x/y/z: time series of fgs data in gei
            GEI_2_DMXL: times seris of rotation matrix
        output:
            fgs_dmxl_x/y/z: times series of fgs data in dmxl
    """
    fgs_gei_x = np.zeros(len(fgs_dmxl_x))
    fgs_gei_y = np.zeros(len(fgs_dmxl_x))
    fgs_gei_z = np.zeros(len(fgs_dmxl_x))

    for i in range(len(fgs_dmxl_x)):
        fgs_gei_x[i] = (
            DMXL_2_GEI[i][0, 0] * fgs_dmxl_x[i]
            + DMXL_2_GEI[i][0, 1] * fgs_dmxl_y[i]
            + DMXL_2_GEI[i][0, 2] * fgs_dmxl_z[i]
        )
        fgs_gei_y[i] = (
            DMXL_2_GEI[i][1, 0] * fgs_dmxl_x[i]
            + DMXL_2_GEI[i][1, 1] * fgs_dmxl_y[i]
            + DMXL_2_GEI[i][1, 2] * fgs_dmxl_z[i]
        )
        fgs_gei_z[i] = (
            DMXL_2_GEI[i][2, 0] * fgs_dmxl_x[i]
            + DMXL_2_GEI[i][2, 1] * fgs_dmxl_y[i]
            + DMXL_2_GEI[i][2, 2] * fgs_dmxl_z[i]
        )

    return [fgs_gei_x, fgs_gei_y, fgs_gei_z]    


def dmxl2smxl(fgs_dmxl_x, fgs_dmxl_y, fgs_dmxl_z, phi):
    """
        input: 
            fgs_dmxl_x/y/z: time series of fgs data in dmxl
            phi: rotation angle 
        output:
            fgs_smxl_x/y/z: times series of fgs data in smxl
    """
    fgs_smxl_x = (
        np.cos(-phi) * fgs_dmxl_x - np.sin(-phi) * fgs_dmxl_y
    )
    fgs_smxl_y = (
        np.sin(-phi) * fgs_dmxl_x + np.cos(-phi) * fgs_dmxl_y
    )
    fgs_smxl_z = fgs_dmxl_z    

    return [fgs_smxl_x, fgs_smxl_y, fgs_smxl_z]


def smxl2dmxl(fgs_smxl_x, fgs_smxl_y, fgs_smxl_z, phi):
    """
        input: 
            fgs_dmxl_x/y/z: time series of fgs data in dmxl
            phi: rotation angle 
        output:
            fgs_smxl_x/y/z: times series of fgs data in smxl
    """   
    fgs_dmxl_x = (
        np.cos(phi) * fgs_smxl_x - np.sin(phi) * fgs_smxl_y
    )
    fgs_dmxl_y = (
        np.sin(phi) * fgs_smxl_x + np.cos(phi) * fgs_smxl_y
    )
    fgs_dmxl_z = fgs_smxl_z

    return [fgs_dmxl_x, fgs_dmxl_y, fgs_dmxl_z]


def smxl2fgm(fgs_smxl_x, fgs_smxl_y, fgs_smxl_z):
    """
        input: 
            fgs_smxl_x/y/z: time series of fgs data in fgm
        output:
            fgs_fgm_x/y/z: times series of fgs data in smxl
    """
    fgs_fgm_x = np.cos(parameter.f) * fgs_smxl_x + np.sin(parameter.f) * fgs_smxl_z
    fgs_fgm_y = np.sin(parameter.f) * fgs_smxl_x - np.cos(parameter.f) * fgs_smxl_z
    fgs_fgm_z = fgs_smxl_y  

    return [fgs_fgm_x, fgs_fgm_y, fgs_fgm_z] 


def fgm2smxl(fgs_fgm_x, fgs_fgm_y, fgs_fgm_z):
    """
        input: 
            fgs_fgm_x/y/z: time series of fgs data in smxl
        output:
            fgs_smxl_x/y/z: times series of fgs data in fgm
    """
    fgs_smxl_x = np.cos(parameter.f) * fgs_fgm_x + np.sin(parameter.f) * fgs_fgm_y
    fgs_smxl_y = fgs_fgm_z
    fgs_smxl_z = np.sin(parameter.f) * fgs_fgm_x - np.cos(parameter.f) * fgs_fgm_y

    return [fgs_smxl_x, fgs_smxl_y, fgs_smxl_z]     