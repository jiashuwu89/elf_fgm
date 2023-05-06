from matplotlib import pyplot as plt
import pandas as pd
from typing import Union
import numpy as np
from scipy.interpolate import griddata
from .. import parameter

def parameter_plot(df: pd.DataFrame, title: str):

    fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(6, 8))
    axes[0].plot(df[["G1", "G2", "G3"]])
    axes[0].legend(["G1", "G2", "G3"])
    axes[0].set_ylabel("Gain")
    axes[0].set_title(title)

    axes[1].plot(df[["th1", "th2", "th3"]])
    axes[1].set_xlim([-1, 358])
    axes[1].legend(["th1", "th2", "th3"])
    axes[1].set_ylabel("theta")

    axes[2].plot(df[["ph1", "ph2", "ph3"]])
    axes[2].legend(["ph1", "ph2", "ph3"])
    axes[2].set_ylabel("phi")

    axes[3].plot(df[["O1/G1", "O2/G2", "O3/G3"]])
    axes[3].legend(["O1/G1", "O2/G2", "O3/G3"])
    axes[3].set_ylabel("offset")

    axes[4].plot(df[["res_med"]])
    axes[4].set_ylabel("res median")

    plt.tight_layout()
    plt.show()
    
def parameter_plot_3d(df: pd.DataFrame, title: str):
    
    att_rot = [
        (i,j) for i in np.arange(-5, 5, 0.5)
        for j in np.arange(-5, 5, 0.5)]
    att_rot.insert(0, (0, 0))

    x, y = zip(*att_rot)
    # Create a regular grid for interpolation
    grid_x, grid_y = np.mgrid[min(x):max(x):50j, min(y):max(y):50j]

    # Interpolate the data to the regular grid
    grid_z = griddata(att_rot, df['res_med'], (grid_x, grid_y), method='cubic')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    ax.plot_surface(grid_x, grid_y, grid_z, cmap='viridis', linewidth=0, antialiased=False)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('3D surface plot')
    plt.show()

    return

def parameter_plot_2d(df: pd.DataFrame, title: str, columns: Union[str, list]):

    if isinstance(columns, str):
        columns = [columns]

    # generate geographic grid
    att_rot = [
        (i,j) for i in np.arange(-parameter.att_loop_width, parameter.att_loop_width, parameter.att_loop_step)
        for j in np.arange(-parameter.att_loop_width, parameter.att_loop_width, parameter.att_loop_step)]
    att_rot.insert(0, (0, 0))

    if len(att_rot) != len(df):
        print("length of data and grid doesn't match!")
        return 

     # Create a regular grid for interpolation
    x, y = zip(*att_rot)
    grid_x, grid_y = np.mgrid[min(x):max(x):50j, min(y):max(y):50j]

    # create a figure
    n_columns = len(columns)
    fig, axes = plt.subplots(ncols=n_columns, nrows=1, figsize=(5 * n_columns, 4 * 1), sharex=True, sharey=True)
    
    for idx, column in enumerate(columns):
        ax = axes[idx] if n_columns >1 else axes

        # Interpolate the data to the regular grid
        grid_z = griddata(att_rot, df[column], (grid_x, grid_y), method='cubic')

        contour = ax.contourf(grid_x, grid_y, grid_z, levels=20, cmap='viridis')
        fig.colorbar(contour, ax=ax)

        ax.set_xticks(np.unique(x))
        ax.set_yticks(np.unique(y))
        ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.5)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title(f"{title}  {column}")


    plt.tight_layout()
    plt.show()

    return

def get_csv(path: str) -> pd.DataFrame:

    df = pd.read_csv(path)

    return df

if __name__ == "__main__":
    

    path = "fgm_utils/fitting_csv/2021-04-14_06_elb_attloop_Gthphi_10_2.csv"
    #path = "2021-03-06_15_elb_attloop_Gthphi.csv"
    #path = "2021-03-17_18_elb_attloop_Gthphi.csv"

    df = get_csv(path)
    #parameter_plot(df, "att_rot", ["G1", "G2", "G3"], "2021-03-06_15-elb")
    #parameter_plot(df, "2021-03-17_18_elb")
    O1_lab = -558
    O2_lab = -277
    O3_lab = -511

    df["O1/G1"] = df["O1/G1"] - O1_lab
    df["O2/G2"] = df["O2/G2"] - O2_lab
    df["O3/G3"] = df["O3/G3"] - O3_lab
    #parameter_plot_2d(df, "2021-03-17_18_elb", ["O1/G1","O2/G2","O3/G3"])
    parameter_plot_2d(df, "2021-03-17_18_elb", "res_med")
