import numpy as np
import pandas as pd
from typing import List, Iterable

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

def att_loop(
        att_x: float, 
        att_y: float, 
        att_z: float,
        angle: float,
        step: float) -> Iterable[List[float]]:

    """
    Generate a cone of vectors around the original vector. check out attitude_cone.ipynb
    
    Args:
        att_x (float): attitude x component
        att_y (float): attitude y component
        att_z (float): attitude z component
        angle (float): angle between the rotated vector and the original vector, degree
        step (float): step from rotate from 0 to 360 deg. 
    
    Returns:
        Tuple of the three compoents of rotated vectors
    """
    v = np.array([att_x, att_y, att_z]) # Your original vector
    v = v / np.linalg.norm(v) # Normalize the vector

    angle_delta = np.radians(angle) # 2 degrees in radians

    # Find an orthogonal vector to the original vector v
    ortho = np.array([-v[1], v[0], 0])
    if np.linalg.norm(ortho) < 1e-8: # exclude situation like v= [0, 0, 1]
        ortho = np.array([0, -v[2], v[1]])
    ortho = ortho / np.linalg.norm(ortho)

    rotated_vector = rodrigues_rotation(v, ortho, angle_delta)

    # Rotate each point on the circle around the given vector v
    rotated_points = [rodrigues_rotation(rotated_vector, v, np.deg2rad(ang)) for ang in range(0, 360, step)]

    return rotated_points

