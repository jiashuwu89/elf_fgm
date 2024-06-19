import numpy as np
from .. import parameter


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

        # if i == 0 or i == len(fgs_igrf_gei_x) -1:
        #     breakpoint()

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


def smxl2fgm(fgs_smxl_x, fgs_smxl_y, fgs_smxl_z, f):
    """
        input: 
            fgs_smxl_x/y/z: time series of fgs data in fgm
        output:
            fgs_fgm_x/y/z: times series of fgs data in smxl
    """
    fgs_fgm_x = np.cos(f) * fgs_smxl_x + np.sin(f) * fgs_smxl_z
    fgs_fgm_y = np.sin(f) * fgs_smxl_x - np.cos(f) * fgs_smxl_z
    fgs_fgm_z = fgs_smxl_y  

    return [fgs_fgm_x, fgs_fgm_y, fgs_fgm_z] 


def fgm2smxl(fgs_fgm_x, fgs_fgm_y, fgs_fgm_z, f):
    """
        input: 
            fgs_fgm_x/y/z: time series of fgs data in smxl
        output:
            fgs_smxl_x/y/z: times series of fgs data in fgm
    """
    fgs_smxl_x = np.cos(f) * fgs_fgm_x + np.sin(f) * fgs_fgm_y
    fgs_smxl_y = fgs_fgm_z
    fgs_smxl_z = np.sin(f) * fgs_fgm_x - np.cos(f) * fgs_fgm_y

    return [fgs_smxl_x, fgs_smxl_y, fgs_smxl_z]     


def gei_obw_matrix(
    fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z, pos_gei_x, pos_gei_y, pos_gei_z
    ):
    """use igrf and pos to get the rotation matrix between gei and obw
        obw: 
            get b first by normalizing the model field (IGRF)
            then get w = (rxb)/|rxb|
            then get o = bxw
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
    pos_gei_x_hat = pos_gei_x / np.sqrt(
        pos_gei_x**2 + pos_gei_y**2 + pos_gei_z**2 
    )
    pos_gei_y_hat = pos_gei_y / np.sqrt(
        pos_gei_x**2 + pos_gei_y**2 + pos_gei_z**2 
    )
    pos_gei_z_hat = pos_gei_z / np.sqrt(
        pos_gei_x**2 + pos_gei_y**2 + pos_gei_z**2 
    )

    OBW_2_GEI = np.zeros((len(fgs_igrf_gei_x), 3, 3))
    GEI_2_OBW = np.zeros((len(fgs_igrf_gei_x), 3, 3))
    for i in range(len(fgs_igrf_gei_x)):
        r_hat = np.array([pos_gei_x_hat[i], pos_gei_y_hat[i], pos_gei_z_hat[i]])
        b_hat = np.array(
            [fgs_igrf_gei_x_hat[i], fgs_igrf_gei_y_hat[i], fgs_igrf_gei_z_hat[i]]
        )
        BperpWest = np.cross(r_hat, b_hat)
        BperpWest_norm = BperpWest / np.linalg.norm(BperpWest)
        BperpOut = np.cross(b_hat, BperpWest_norm)
        BperpOut_norm = BperpOut / np.linalg.norm(BperpOut)
        
        GEI_2_OBW[i, 0, :] = BperpOut_norm
        GEI_2_OBW[i, 1, :] = b_hat
        GEI_2_OBW[i, 2, :] = BperpWest_norm

        # if i == 0 or i == len(fgs_igrf_gei_x) -1:
        #     breakpoint()

        OBW_2_GEI[i, :, :] = np.linalg.inv(GEI_2_OBW[i, :, :])

    return [GEI_2_OBW, OBW_2_GEI]


def gei2obw(fgs_gei_x, fgs_gei_y, fgs_gei_z, GEI_2_OBW):
    """TODO:combine this with dmxl2gei
    """
    fgs_obw_x = np.zeros(len(fgs_gei_x))
    fgs_obw_y = np.zeros(len(fgs_gei_x))
    fgs_obw_z = np.zeros(len(fgs_gei_x))

    for i in range(len(fgs_obw_x)):
        fgs_obw_x[i] = (
            GEI_2_OBW[i][0, 0] * fgs_gei_x[i]
            + GEI_2_OBW[i][0, 1] * fgs_gei_y[i]
            + GEI_2_OBW[i][0, 2] * fgs_gei_z[i]
        )
        fgs_obw_y[i] = (
            GEI_2_OBW[i][1, 0] * fgs_gei_x[i]
            + GEI_2_OBW[i][1, 1] * fgs_gei_y[i]
            + GEI_2_OBW[i][1, 2] * fgs_gei_z[i]
        )
        fgs_obw_z[i] = (
            GEI_2_OBW[i][2, 0] * fgs_gei_x[i]
            + GEI_2_OBW[i][2, 1] * fgs_gei_y[i]
            + GEI_2_OBW[i][2, 2] * fgs_gei_z[i]
        )

    return [fgs_obw_x, fgs_obw_y, fgs_obw_z] 


def obw2gei(fgs_obw_x, fgs_obw_y, fgs_obw_z, OBW_2_GEI):

    fgs_gei_x = np.zeros(len(fgs_obw_x))
    fgs_gei_y = np.zeros(len(fgs_obw_x))
    fgs_gei_z = np.zeros(len(fgs_obw_x))

    for i in range(len(fgs_obw_x)):
        fgs_gei_x[i] = (
            OBW_2_GEI[i][0, 0] * fgs_obw_x[i]
            + OBW_2_GEI[i][0, 1] * fgs_obw_y[i]
            + OBW_2_GEI[i][0, 2] * fgs_obw_z[i]
        )
        fgs_gei_y[i] = (
            OBW_2_GEI[i][1, 0] * fgs_obw_x[i]
            + OBW_2_GEI[i][1, 1] * fgs_obw_y[i]
            + OBW_2_GEI[i][1, 2] * fgs_obw_z[i]
        )
        fgs_gei_z[i] = (
            OBW_2_GEI[i][2, 0] * fgs_obw_x[i]
            + OBW_2_GEI[i][2, 1] * fgs_obw_y[i]
            + OBW_2_GEI[i][2, 2] * fgs_obw_z[i]
        )

    return [fgs_gei_x, fgs_gei_y, fgs_gei_z] 


def geo_nec_matrix(
    pos_geo_x, pos_geo_y, pos_geo_z
    ):
    """use pos to get the rotation matrix between geo and nec. Code from Sasha (anton) 
            
    """
    #e_n, e_e, e_c - NEC basis
    e_c=-np.stack([pos_geo_x, pos_geo_y, pos_geo_z]).T
    e_c=e_c/np.sqrt(np.sum(e_c**2,axis=1)).reshape(-1,1)  #normalization

    e_e=np.cross(e_c,np.array([0,0,1]))
    e_e[np.all(e_e==0,axis=1)]=[0,1,0]  #make e_e zero vector into [0, 1, 0], zero vector can appar when e_c is parallel to [0, 0, 1]
    e_e=e_e/np.sqrt(np.sum(e_e**2,axis=1)).reshape(-1,1) #normaliation

    e_n=np.cross(e_e,e_c)

    NEC_2_GEO=np.zeros([e_c.shape[0],3,3]) #NEC -> ITRF transformation matrix 

    NEC_2_GEO[:,:,0]=e_n
    NEC_2_GEO[:,:,1]=e_e
    NEC_2_GEO[:,:,2]=e_c

    GEO_2_NEC = np.transpose(NEC_2_GEO, axes=(0, 2, 1))

    return [GEO_2_NEC, NEC_2_GEO]

def geo2nec(B_x, B_y, B_z, GEO_2_NEC):
    B_nec=np.squeeze(np.matmul(GEO_2_NEC,np.column_stack((B_x, B_y, B_z))[:,:,np.newaxis]))
    return B_nec