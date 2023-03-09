from ..function import coordinate
import numpy as np
from .. import parameter
import pandas as pd
import datetime
from ..function import preprocess

def mva_spinaxis_angdiff(mvamin_fgm_x, mvamin_fgm_y, mvamin_fgm_z):
    [x, y, z] = coordinate.fgm2smxl(mvamin_fgm_x, mvamin_fgm_y, mvamin_fgm_z, 0.5*np.pi-parameter.f)  
    vec1 = [x, y, z]
    vec2 = [0, 0, 1]
    cos_theta = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    theta = np.arccos(cos_theta)
    angle_degrees = np.degrees(theta)
    
    return angle_degrees

def mva_spinaxis_ang(mvamin_fgm_x, mvamin_fgm_y, mvamin_fgm_z):
    [x, y, z] = coordinate.fgm2smxl(mvamin_fgm_x, mvamin_fgm_y, mvamin_fgm_z, 0.5*np.pi-parameter.f)
    alpha_x = np.arccos(x) *180 / np.pi
    alpha_y = np.arccos(y) *180 / np.pi
    alpha_z = np.arccos(z) *180 / np.pi

    #return [alpha_x, alpha_y, alpha_z]
    return [x, y, z]


