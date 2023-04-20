import numpy as np
import pandas as pd

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