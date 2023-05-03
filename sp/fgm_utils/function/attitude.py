import numpy as np
import pandas as pd
from typing import List, Iterable
import math
import matplotlib.pyplot as plt

def rotate_vector(vector, angle_degrees, axis='x'):

    angle_radians = np.radians(angle_degrees)

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
    
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.acos(z / x) # polar angle in radians
    phi = math.atan2(y, x) # azimuthal angle, in radians

    return r, theta, phi

def sphere2cart(r, theta, phi):

    x = r * math.sin(theta) * math.cos(phi)
    y = r * math.sin(theta) * math.sin(phi)
    z = r * math.cos(theta)

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

    rotated_points = [sphere2cart(r_0, theta_0 + np.deg2rad(i), phi_0 + np.deg2rad(j)) 
                      for i in np.arange(-width/2, width/2, step)
                      for j in np.arange(-width/2, width/2, step)]   
    
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
