#from pytplot import store_data, tplot, get_data
#from pyspedas.cotrans.minvar_matrix_make import minvar_matrix_make
import numpy as np
import matplotlib.pyplot as plt


def mva(ctime, ctimestamp, fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z):
    """
    Get minimum variance vector direction. Because it use the mva func from pyspedas, time is required to create a tplot variable

    Parameter
        ctime: List of time
        ctimestampe: timestampe of the first time
        fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z: list of fgm data

    Return
        mva_ang: angle in deg of minimum variance vector w/r to fgm_x

    """
    # TODO: apply mva without pytplot package
    print("This function is pending... run wihtout mva for now!!!")
    return 0
    """
    store_data('fgs_ful_fgm', data={'x': ctime+ctimestamp, 'y': np.transpose(np.array([fgs_ful_fgm_0th_x, fgs_ful_fgm_0th_y, fgs_ful_fgm_0th_z]))})
    #tplot('fgs_ful_fgm')
    minvar_matrix_make('fgs_ful_fgm')
    mva_time, mva_vec = get_data('fgs_ful_fgm_mva_mat')
    #fgs_maxvarvec = -mva_flip(mva_vec[:, 0 ,:])
    fgs_maxvarvec = mva_vec[:, 0 ,:]
    #fgs_midvarvec = mva_flip(mva_vec[:, 1 ,:])
    fgs_midvarvec = mva_vec[:, 1 ,:]
    #fgs_minvarvec = -mva_flip(mva_vec[:, 2 ,:])
    fgs_minvarvec = mva_vec[:, 2 ,:]

    fgs_maxvarvec_ang = np.rad2deg(np.arccos(fgs_maxvarvec))
    fgs_midvarvec_ang = np.rad2deg(np.arccos(fgs_midvarvec))
    fgs_minvarvec_ang = np.rad2deg(np.arccos(fgs_minvarvec))

    if fgs_minvarvec_ang[:,0] > 90:
        fgs_maxvarvec = -fgs_maxvarvec
        fgs_midvarvec = -fgs_midvarvec
        fgs_minvarvec = -fgs_minvarvec
        fgs_maxvarvec_ang = np.rad2deg(np.arccos(fgs_maxvarvec))
        fgs_midvarvec_ang = np.rad2deg(np.arccos(fgs_midvarvec))
        fgs_minvarvec_ang = np.rad2deg(np.arccos(fgs_minvarvec))


    print(f"Minimum variance vector dirction:")
    print(f"MinVar median angle to fgm_x [deg]: {np.median(fgs_minvarvec_ang[:,0])}")
    print(f"MinVar median angle to fgm_y [deg]: {np.median(fgs_minvarvec_ang[:,1])}")
    print(f"MinVar median angle to fgm_z [deg]: {np.median(fgs_minvarvec_ang[:,2])}")

    print(f"Middle variance vector dirction:")
    print(f"MidVar median angle to fgm_x [deg]: {np.median(fgs_midvarvec_ang[:,0])}")
    print(f"MidVar median angle to fgm_y [deg]: {np.median(fgs_midvarvec_ang[:,1])}")
    print(f"MidVar median angle to fgm_z [deg]: {np.median(fgs_midvarvec_ang[:,2])}")

    print(f"Maxium variance vector dirction:")
    print(f"MaxVar median angle to fgm_x [deg]: {np.median(fgs_maxvarvec_ang[:,0])}")
    print(f"MaxVar median angle to fgm_y [deg]: {np.median(fgs_maxvarvec_ang[:,1])}")
    print(f"MaxVar median angle to fgm_z [deg]: {np.median(fgs_maxvarvec_ang[:,2])}")

    return np.median(fgs_minvarvec_ang[:,0])
    """

def mva_flip(vec):
    """check the sign of adjacent vectors, flip if the signs are opposite
    """
    for i in range(1, len(vec)):
        if np.dot(vec[i-1], vec[i]) < 0:
            vec[i] = -vec[i]

    return vec

def mva_plot(minvec, midvec, maxvec):

    minvec = np.array(minvec)
    midvec = np.array(midvec)
    maxvec = np.array(maxvec)

    if len(minvec) == 3:
        if isinstance(minvec[0], list):
            num = 3
        else:
            num = 1
    else:
        num = len(minvec)
    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    line_styles = ['-', '--', '-.', ':']
    # Plot the vectors
    for i in range(num):
        if i == 0:
            ax.plot([0, minvec[i, 0]], [0, minvec[i, 1]], [0, minvec[i, 2]], color='r', label='MinVec', linestyle=line_styles[i % len(line_styles)])
            ax.plot([0, midvec[i, 0]], [0, midvec[i, 1]], [0, midvec[i, 2]], color='b', label='MidVec', linestyle=line_styles[i % len(line_styles)])
            ax.plot([0, maxvec[i, 0]], [0, maxvec[i, 1]], [0, maxvec[i, 2]], color='g', label='MaxVec', linestyle=line_styles[i % len(line_styles)])
        else:
            ax.plot([0, minvec[i, 0]], [0, minvec[i, 1]], [0, minvec[i, 2]], color='r', label='', linestyle=line_styles[i % len(line_styles)])
            ax.plot([0, midvec[i, 0]], [0, midvec[i, 1]], [0, midvec[i, 2]], color='b', label='', linestyle=line_styles[i % len(line_styles)])
            ax.plot([0, maxvec[i, 0]], [0, maxvec[i, 1]], [0, maxvec[i, 2]], color='g', label='', linestyle=line_styles[i % len(line_styles)])
 
    # Set axis labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Add a legend
    ax.legend()

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])

    # Display the plot
    plt.show()


if __name__ == "__main__":

    """how to run: poetry run python -m sp.fgm_utils.function.mva
    """
    # minvec_1 = [-0.59340312, 0.80455562,  -0.02372729]
    # minvec_2 = [8.06379100e-01,  -5.91398973e-01, -2.93239046e-05]
    # midvec_1 = [0.01808225, -0.01614588, -0.99970613]
    # midvec_2 = [-0.00331838, -0.00447507, -0.99998448]
    # maxvec_1 = [0.80470228, 0.59365778, 0.00496716]
    # maxvec_2 = [-0.59138966,  -0.80636668, 0.00557109]

    minvec_1 = [0.59365421, -0.80438989,  0.02305605]
    minvec_2 = [-8.05641836e-01,  5.92402917e-01,  1.27893616e-04]

    midvec_1 = [0.0222197 , -0.01225501, -0.999678]
    midvec_2 = [-0.07012253, -0.09514917, -0.99299016]

    maxvec_1 = [-0.80441343, -0.59397535, -0.01059806]
    maxvec_2 = [0.5882381 ,  0.80000339, -0.11819694]

    mva_plot([minvec_1, minvec_2], [midvec_1, midvec_2], [maxvec_1, maxvec_2])
