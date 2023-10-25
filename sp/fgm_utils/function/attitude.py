import numpy as np
import pandas as pd
from typing import List, Iterable
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from .. import parameter
from .Bplot import B_ctime_plot_single

def rotate_vector(vector, angle_radians, axis='x'):

    if axis == 'x':
        rotation_matrix = np.array([[1, 0, 0],
                                    [0, np.cos(angle_radians), -np.sin(angle_radians)],
                                    [0, np.sin(angle_radians), np.cos(angle_radians)]])
    elif axis == 'y':
        rotation_matrix = np.array([[np.cos(angle_radians), 0, np.sin(angle_radians)],
                                    [0, 1, 0],
                                    [-np.sin(angle_radians), 0, np.cos(angle_radians)]])
    elif axis == 'z':
        rotation_matrix = np.array([[np.cos(angle_radians), -np.sin(angle_radians), 0],
                                    [np.sin(angle_radians), np.cos(angle_radians), 0],
                                    [0, 0, 1]])
    else:
        raise ValueError("Invalid rotation axis. Must be 'x', 'y', or 'z'.")

    return np.dot(rotation_matrix, vector)


def att_rot(att_cdf: pd.DataFrame, att_rot_ang: float, rot_axis: str) -> pd.DataFrame:

    att_gei = att_cdf['elb_att_gei'].to_list()
    att_rot_gei = [rotate_vector(att, att_rot_ang, axis = rot_axis) for att in att_gei]
    att_cdf['elb_att_gei'] = att_rot_gei

    return att_cdf


def rodrigues_rotation(v, k, theta):
    return v * np.cos(theta) + np.cross(k, v) * np.sin(theta) + k * np.dot(k, v) * (1 - np.cos(theta))

def cart2sphere(x, y, z):
    """
    Convert cartesian to spherical coordinate

    Parameter
        x, y, z: in cartesian coordinate
    Return
        r, theta, phi: theta and phi are in rad
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r) # polar angle in radians
    phi = np.arctan2(y, x) # azimuthal angle, in radians

    return r, theta, phi

def sphere2cart(r, theta, phi):
    """
    Convert spherical to cartesian coordinate

    Parameter
        r, theta, phi: theta and phi are in rad
    Return
        x, y, z: same unit as r
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z


def att_loop(
        att_x: float, 
        att_y: float, 
        att_z: float,
        width: float,
        step: float,) -> Iterable[List[float]]:

    """
    Generate a grid around original attitude
    
    Args:
        att_x (float): attitude x component
        att_y (float): attitude y component
        att_z (float): attitude z component
        width (float): degree, width of attitude loop in both lat/lon direction
        step (float): degree, step of atttitude loop in both lat/lon direction 
    
    Returns:
        Tuple of the three compoents of rotated vectors
    """
    v = np.array([att_x, att_y, att_z]) # Your original vector
    v = v / np.linalg.norm(v) # Normalize the vector
    [r_0, theta_0, phi_0] = cart2sphere(*v)

    # rotate attitude to lat = 0 lon = 0 
    #rot_vector1 = rotate_vector(v, -phi_0, axis='z')
    #rot_vector2 = rotate_vector(rot_vector1, np.pi/2 - theta_0, axis='y')
    
    rotated_points = []
    for i in np.arange(-width, width, step):
        for j in np.arange(-width, width, step):
            # in lat/lon = 0 grid
            vec_latlon = sphere2cart(1, np.pi/2 + np.deg2rad(i), 0 + np.deg2rad(j)) 
            # rot vector to gei
            rot_vector1 = rotate_vector(vec_latlon,  theta_0 - np.pi/2, axis='y')
            rot_vector2 = rotate_vector(rot_vector1, phi_0, axis='z')
            rotated_points.append(list(rot_vector2))
            #rotated_points.append(list(vec_latlon))
    
    return rotated_points


def att_plot(v, rotated_points,):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Add the attitude vector to the plot
    ax.quiver(0, 0, 0, v[0], v[1], v[2], color='red', label='Attitude vector', linewidth=1.5)

    # Add the generated points on the circle to the plot
    for points in rotated_points:
        ax.plot([0, points[0]], [0, points[1]], [0, points[2]], color='blue')

    # Add a unit sphere to the plot
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z, color='gray', alpha=0.1)

    # Set axis labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Add a legend
    ax.legend()

    # Show the plot
    plt.show()


def att_determine_func(x, Oper, Gper, Opar, Gpar, th, ph, igrf_x, igrf_y, igrf_z):
    """
    FIT Bper and Bpar to get five parameters 

    Parameters
        params: Oper, Gper, Opar, Gpar, thS, phS
    Returns
        Bmodel(np.ndarray): output of the fitted result
    """
    #Gper, Opar, Gpar, th, ph = params
    #igrf_x, igrf_y, igrf_z = args
    sth = np.sin(th)
    cth = np.cos(th)
    sph = np.sin(ph)
    cph = np.cos(ph)
 
    Bpar = igrf_x*sth*cph + igrf_y*sth*sph + igrf_z*cth
    Bparvec = np.array([Bpar*sth*cph,Bpar*sth*sph,Bpar*cth])
    Bparout = Gpar * Bpar - Opar 
    Bpervec = np.array([igrf_x, igrf_y, igrf_z]) - Bparvec
    Bper = np.sqrt(np.sum(Bpervec**2, axis=0))
    Bperout = Gper * Bper - Oper   
    Bmodel = np.append(Bparout, Bperout)

    return Bmodel

def ang_twovec(v1, v2):
    """
    Get angle between two vectors

    Parameter
        v1, v2: two vectors in the same coordinate
    Return
        ang: rad
    """
    # Convert the input lists to numpy arrays
    v1_arr = np.array(v1)
    v2_arr = np.array(v2)

    # Calculate the dot product of the two vectors
    dot_product = np.dot(v1_arr, v2_arr)

    # Calculate the magnitudes of the vectors
    v1_magnitude = np.linalg.norm(v1_arr)
    v2_magnitude = np.linalg.norm(v2_arr)

    # Calculate the cosine of the angle using the dot product and magnitudes
    cos_ang = dot_product / (v1_magnitude * v2_magnitude)

    # Calculate the angle in radians using the arccosine function
    ang = np.arccos(cos_ang)

    return ang


def att_determine(
        fgs_ful_smxl_x: np.ndarray, fgs_ful_smxl_y: np.ndarray, fgs_ful_smxl_z: np.ndarray, 
        fgs_igrf_gei_x: np.ndarray, fgs_igrf_gei_y: np.ndarray, fgs_igrf_gei_z: np.ndarray, 
        att_gei_x: np.ndarray, att_gei_y: np.ndarray, att_gei_z: np.ndarray, datestr: str):
    """
    Fit Bper and Bpar to IGRF to refine attitude th, phi in GEI coordinate. 
    
    This function is similar to idl att determination 

    Parameters
        fgs_ful_fgm_x, fgs_ful_fgm_y, fgs_ful_fgm_z: fgm calibrated data in smxl coordinate, x,y are in spin plane, z is along spin axis
        fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z:  igrf in gei coordiate
        att_gei_x, att_gei_y, att_gei_z:  original attitude list in gei coordiate 
    
    Returns 
        att_gei_refine_x, att_gei_refine_y, att_gei_refine_z: refined attitude list in gei coorindate
    """

    
    # separate fgm measurement into Bpar and Bper
    Bper_measure = np.sqrt(fgs_ful_smxl_x**2 + fgs_ful_smxl_y**2)
    Bpar_measure = fgs_ful_smxl_z
    
    if parameter.att_split == False:  #fit only one attitude
        # use the attitude in the middle as initial guess
        att_gei_init_x =  np.median(att_gei_x)
        att_gei_init_y =  np.median(att_gei_y)
        att_gei_init_z =  np.median(att_gei_z)
        r, th_init, ph_init = cart2sphere(att_gei_init_x, att_gei_init_y, att_gei_init_z)

        B_measure = np.concatenate((Bpar_measure, Bper_measure))
        # define fitting func
        x = range(len(B_measure))
        att_determine_func_args = lambda x, Oper, Gper, Opar, Gpar, th, ph: att_determine_func(x, Oper, Gper, Opar, Gpar, th, ph, fgs_igrf_gei_x, fgs_igrf_gei_y, fgs_igrf_gei_z)
        params_opt, params_cov = curve_fit(
            att_determine_func_args, x, B_measure, 
            p0=[0, 1, 0, 1, th_init, ph_init], 
            #method = 'dogbox'
            #bounds = (
            #    [-np.inf, -np.inf, -np.inf, np.deg2rad(np.rad2deg(th_init)-5), np.deg2rad(np.rad2deg(ph_init)-5)],
            #    [np.inf, np.inf, np.inf, np.deg2rad(np.rad2deg(th_init)+5), np.deg2rad(np.rad2deg(ph_init)+5)]),
            )
        signal_fit = att_determine_func_args(x, *params_opt)
        signal_fit_init = att_determine_func_args(x, params_opt[0],params_opt[1],params_opt[2], params_opt[3], th_init, ph_init)
        Bpar_fit = signal_fit[0:len(fgs_ful_smxl_x)]
        Bper_fit = signal_fit[len(fgs_ful_smxl_x):]
        th_final = params_opt[4]
        ph_final = params_opt[5]
        # get new attitude in gei
        att_gei_refine_x_mid, att_gei_refine_y_mid, att_gei_refine_z_mid = sphere2cart(1, th_final, ph_final)
        att_ang_diff = ang_twovec([att_gei_init_x, att_gei_init_y, att_gei_init_z], [att_gei_refine_x_mid, att_gei_refine_y_mid, att_gei_refine_z_mid])

        # update all att according to the middle point changes
        att_gei_refine_x = []
        att_gei_refine_y = []
        att_gei_refine_z = []
        for x, y, z in zip(att_gei_x, att_gei_y, att_gei_z):
            r_i, th_i, ph_i = cart2sphere(x, y, z)
            th_i = th_i + th_final - th_init
            ph_i = ph_i + ph_final - ph_init
            x_refine, y_refine, z_refine = sphere2cart(1, th_i, ph_i)
            att_gei_refine_x.append(x_refine)
            att_gei_refine_y.append(y_refine)
            att_gei_refine_z.append(z_refine)
        
        # output 
        print('===============att determine==============')
        print("initial att in gei: [{:.5f}, {:.5f}, {:.5f}]".format(att_gei_init_x, att_gei_init_y, att_gei_init_z))
        print("initial th and ph(deg): {:.5f}, {:.5f}".format(np.rad2deg(th_init), np.rad2deg(ph_init)))
        print("final att in gei: [{:.5f}, {:.5f}, {:.5f}]".format(att_gei_refine_x_mid, att_gei_refine_y_mid, att_gei_refine_z_mid))
        print("final th and ph(deg): {:.5f}, {:.5f}".format(np.rad2deg(th_final), np.rad2deg(ph_final)))
        print("att ang difference(deg): {:.5f}".format(np.rad2deg(att_ang_diff)))
        print('==============att determine end===========')

        att_split_idx = None
    else:  #fit multiple attitude
        att_split_idx1 = []  # start idx for att split snippet
        att_split_idx2 = []  # end idx for att split snippet
        att_split_len = [] # length of att split snippet
        if parameter.att_split_num is not None: # divide interval into equal length snippets
            idx2 = 0
            att_split_num = parameter.att_split_num
            for i in range(att_split_num):
                idx1 = idx2 
                idx2 = idx2 + int(len(fgs_igrf_gei_x)/att_split_num)        
                if i == att_split_num - 1:  #if this is the last interval, include all remaining points
                    idx2 = len(fgs_igrf_gei_x)
                length = idx2 - idx1
                att_split_idx1.append(idx1)
                att_split_idx2.append(idx2)
                att_split_len.append(length)
                att_split_idx = att_split_idx1
        else:  # divide according to att_split_idx
            att_split_num = len(parameter.att_split_idx)
            for i, _ in enumerate(parameter.att_split_idx):
                idx1 = parameter.att_split_idx[i]
                idx2 = parameter.att_split_idx[i+1] if i < len(parameter.att_split_idx)-1 else len(fgs_igrf_gei_x)
                att_split_idx1.append(idx1)
                att_split_idx2.append(idx2)
                length = idx2 - idx1
                att_split_len.append(length)
                att_split_idx = parameter.att_split_idx
 
        Bpar_fit = []
        Bper_fit = []
        Bpar_fit2 = []  #mapped attitude
        Bper_fit2 = []
        att_gei_refine_x = []
        att_gei_refine_y = []
        att_gei_refine_z = []

        for i in range(att_split_num):
            idx1 = att_split_idx1[i]
            idx2 = att_split_idx2[i]
            length = att_split_len[i]

            # use the attitude in the middle as initial guess
            att_gei_init_x =  np.median(att_gei_x[idx1:idx2])
            att_gei_init_y =  np.median(att_gei_y[idx1:idx2])
            att_gei_init_z =  np.median(att_gei_z[idx1:idx2])
            r, th_init, ph_init = cart2sphere(att_gei_init_x, att_gei_init_y, att_gei_init_z)

            B_measure = np.concatenate((Bpar_measure[idx1:idx2], Bper_measure[idx1:idx2]))

            # define fitting func
            x = range(2*length)
            att_determine_func_args = lambda x, Oper, Gper, Opar, Gpar, th, ph: att_determine_func(
                x, Oper, Gper, Opar, Gpar, th, ph, fgs_igrf_gei_x[idx1:idx2], fgs_igrf_gei_y[idx1:idx2], fgs_igrf_gei_z[idx1:idx2])
            params_opt, params_cov = curve_fit(
                att_determine_func_args, x, B_measure, 
                p0=[0, 1, 0, 1, th_init, ph_init], 
                #method = 'dogbox'
                #bounds = (
                #    [-np.inf, -np.inf, -np.inf, np.deg2rad(np.rad2deg(th_init)-5), np.deg2rad(np.rad2deg(ph_init)-5)],
                #    [np.inf, np.inf, np.inf, np.deg2rad(np.rad2deg(th_init)+5), np.deg2rad(np.rad2deg(ph_init)+5)]),
                )
            signal_fit = att_determine_func_args(x, *params_opt)
            Bpar_fit = np.append(Bpar_fit, signal_fit[:length])
            Bper_fit = np.append(Bper_fit, signal_fit[length:])
            #params_opt[4] = th_init
            #params_opt[5] = ph_init
            #breakpoint()
            th_final = params_opt[4]
            ph_final = params_opt[5]
            
            # get new attitude in gei
            att_gei_refine_x_mid, att_gei_refine_y_mid, att_gei_refine_z_mid = sphere2cart(1, th_final, ph_final)
            att_ang_diff = ang_twovec(
                [att_gei_init_x, att_gei_init_y, att_gei_init_z], 
                [att_gei_refine_x_mid, att_gei_refine_y_mid, att_gei_refine_z_mid])

            # output
            print('===============att determine==============')
            print("initial att in gei: [{:.5f}, {:.5f}, {:.5f}]".format(att_gei_init_x, att_gei_init_y, att_gei_init_z))
            print("initial th and ph(deg): {:.5f}, {:.5f}".format(np.rad2deg(th_init), np.rad2deg(ph_init)))
            print("final att in gei: [{:.5f}, {:.5f}, {:.5f}]".format(att_gei_refine_x_mid, att_gei_refine_y_mid, att_gei_refine_z_mid))
            print("final th and ph(deg): {:.5f}, {:.5f}".format(np.rad2deg(th_final), np.rad2deg(ph_final)))
            print("att ang difference(deg): {:.5f}".format(np.rad2deg(att_ang_diff)))
            print('==============att determine end===========')

            # update all att according to the middle point changes
            for x, y, z, igrf_x, igrf_y, igrf_z in zip(
                att_gei_x[idx1:idx2], att_gei_y[idx1:idx2], att_gei_z[idx1:idx2], 
                fgs_igrf_gei_x[idx1:idx2], fgs_igrf_gei_y[idx1:idx2], fgs_igrf_gei_z[idx1:idx2]):
                r_i, th_i, ph_i = cart2sphere(x, y, z)
                th_i = th_i + th_final - th_init
                ph_i = ph_i + ph_final - ph_init
                x_refine, y_refine, z_refine = sphere2cart(1, th_i, ph_i)
                att_gei_refine_x.append(x_refine)
                att_gei_refine_y.append(y_refine)
                att_gei_refine_z.append(z_refine)            
                signal_fit2 = att_determine_func(0, params_opt[0], params_opt[1], params_opt[2], params_opt[3], th_i, ph_i, igrf_x, igrf_y, igrf_z)
                Bpar_fit2 = np.append(Bpar_fit2, signal_fit2[0])
                Bper_fit2 = np.append(Bper_fit2, signal_fit2[1])
  
    # plot Bper Bpar fit
    if parameter.makeplot == True:
        #B_ctime_plot_single(
        #    range(len(fgs_ful_smxl_x)), 
        #    [Bpar_measure, signal_fit_init[0:len(fgs_ful_smxl_x)], signal_fit[0:len(fgs_ful_smxl_x)]],
        #    legend=['Bpar','Bpar_init','Bpar_fit'],title="Bpar")
        # B_ctime_plot_single(
        #     np.array(range(len(fgs_ful_smxl_x))), 
        #     [Bpar_measure, Bpar_fit, Bpar_fit2],
        #     legend=['Bpar','Bpar_fit','Bpar_fit2'],title="Bpar", datestr=datestr)
        B_ctime_plot_single(
            np.array(range(len(fgs_ful_smxl_x)))/10., 
            [Bpar_measure, Bpar_fit],
            #legend=['Bpar','Bpar_fit'],title="Bpar", datestr=datestr, cross_times=att_split_idx)
            legend=['Bpar','Bpar_fit'],title="Bpar", datestr=datestr)
        #B_ctime_plot_single(
        #    range(len(fgs_ful_smxl_x)), 
        #    [Bper_measure, signal_fit_init[len(fgs_ful_smxl_x):], signal_fit[len(fgs_ful_smxl_x):]], 
        #    legend=['Bper','Bper_init','Bper_fit'],title="Bper")
        # B_ctime_plot_single(
        #    np.array(range(len(fgs_ful_smxl_x))), 
        #    [Bper_measure, Bper_fit, Bper_fit2], 
        #    legend=['Bper','Bper_fit','Bper_fit2'],title="Bper", datestr=datestr)
        B_ctime_plot_single(
           np.array(range(len(fgs_ful_smxl_x)))/10., 
           [Bper_measure, Bper_fit], 
        #   legend=['Bper','Bper_fit'],title="Bper", datestr=datestr, cross_times=att_split_idx)
            legend=['Bper','Bper_fit'],title="Bper", datestr=datestr)
        breakpoint()


    return att_gei_refine_x, att_gei_refine_y, att_gei_refine_z
    
    