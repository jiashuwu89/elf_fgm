import pandas as pd
import matplotlib.pyplot as plt

def Gain_f(f, Gain_x, Gain_y, Gain_z, mission = "", filename = None):
    """plot Gain x, y, z as a function of rotation angle
    """
    fig, ax = plt.subplots(1, figsize=(15,7))

    ax.scatter(f, Gain_x, label='Gain x')
    ax.scatter(f, Gain_y, label='Gain y')
    ax.scatter(f, Gain_z, label='Gain z')
    ax.set_xlabel('rotation angle (deg)')
    ax.set_ylabel('Gain')
    ax.legend() 
    plt.title(f"{mission}")
    filename = "Gain_f" if filename is None else filename + "_Gain_f"
    plt.savefig(f"{filename}_{mission}") if mission != "" else plt.savefig(f"{filename}")
    plt.close()


def STD_f(f, Gain_x, Gain_y, Gain_z, mission = "", filename = None):
    """plot STD x, y, z as a function of rotation angle
    """
    fig, ax = plt.subplots(1, figsize=(15,7))

    ax.scatter(f, Gain_x, label='STD x')
    ax.scatter(f, Gain_y, label='STD y')
    ax.scatter(f, Gain_z, label='STD z')
    ax.set_xlabel('rotation angle (deg)')
    ax.set_ylabel('STD')
    ax.legend() 
    plt.title(f"{mission}")
    filename = "STD_f" if filename is None else filename + "_STD_f"
    plt.savefig(f"{filename}_{mission}") if mission != "" else plt.savefig(f"{filename}")
    plt.close()

if __name__ == "__main__":
    
    mission = "elb"
    date = "20220429_071357"
    rotAng_df = pd.read_csv(f"rotAng_loop_{mission}_360_{date}.csv")

    Gain_f(rotAng_df['rotate_ang'], rotAng_df['G11'], rotAng_df['G22'], rotAng_df['G33'], mission = mission, filename = date)
    STD_f(rotAng_df['rotate_ang'], rotAng_df['res_dmxl_x'], rotAng_df['res_dmxl_y'], rotAng_df['res_dmxl_z'], mission = mission, filename = date)
    
    breakpoint()