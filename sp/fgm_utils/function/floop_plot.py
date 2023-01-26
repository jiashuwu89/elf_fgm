import pandas as pd
import matplotlib.pyplot as plt

def Gain_f(f, Gain_x, Gain_y, Gain_z, mission = None):
    """plot Gain x, y, z as a function of rotation angle
    """
    fig, ax = plt.subplots(1, figsize=(15,7))

    ax.scatter(f, Gain_x, label='Gain x')
    ax.scatter(f, Gain_y, label='Gain y')
    ax.scatter(f, Gain_z, label='Gain z')
    ax.set_xlabel('rotation angle (deg)')
    ax.set_ylabel('Gain')
    ax.legend() 
    plt.savefig(f"Gain_f_{mission}") if mission is not None else plt.savefig(f"Gain_f")
    ax.title(f"{mission}")
    plt.close()

if __name__ == "__main__":
    
    rotAng_A = pd.read_csv("rotAng_loop_ela.csv")
    rotAng_B = pd.read_csv("rotAng_loop_elb.csv")

    Gain_f(rotAng_A['rotate_ang'], rotAng_A['G11'], rotAng_A['G22'], rotAng_A['G33'], mission = 'ela')
    Gain_f(rotAng_B['rotate_ang'], rotAng_B['G11'], rotAng_B['G22'], rotAng_B['G33'], mission = 'elb')
    breakpoint()