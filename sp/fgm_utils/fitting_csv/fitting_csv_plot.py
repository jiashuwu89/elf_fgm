from matplotlib import pyplot as plt
import pandas as pd
from typing import List
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

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

def parameter_plot_2d(df: pd.DataFrame, title: str, column: str):

    att_rot = [
        (i,j) for i in np.arange(-10, 10, 2)
        for j in np.arange(-10, 10, 2)]
    att_rot.insert(0, (0, 0))

    x, y = zip(*att_rot)
    # Create a regular grid for interpolation
    grid_x, grid_y = np.mgrid[min(x):max(x):50j, min(y):max(y):50j]

    # Interpolate the data to the regular grid
    grid_z = griddata(att_rot, df[column], (grid_x, grid_y), method='cubic')

    # Create the contour plot
    plt.contourf(grid_x, grid_y, grid_z, levels=20, cmap='viridis')
    plt.colorbar()  # Display a colorbar with the scale
    plt.xticks(np.unique(x))
    plt.yticks(np.unique(y))
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.5)

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f"{title}  {column}")
    plt.show()

    return

def get_csv(path: str) -> pd.DataFrame:

    df = pd.read_csv(path)

    return df

if __name__ == "__main__":
    
    #path = "2021-03-06_15_elb_attloop_Gthphi.csv"
    #path = "2021-03-17_18_elb_attloop_Gthphi.csv"
    path = "2021-03-17_18_elb_attloop_Gthphi.csv"
    df = get_csv(path)
    #parameter_plot(df, "att_rot", ["G1", "G2", "G3"], "2021-03-06_15-elb")
    #parameter_plot(df, "2021-03-17_18_elb")
    parameter_plot_2d(df, "2021-03-17_18_elb", "O3/G3")
