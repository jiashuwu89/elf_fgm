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

def Allpara_f(f, G11, G12, G13, O1, G21, G22, G23, O2, G31, G32, G33, O3, mission = "", filename = None):
    """plot Gain x, y, z as a function of rotation angle
    """
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7,7))

    ax1.plot(f, G11, label='Gain x')
    ax1.plot(f, G22, label='Gain y')
    ax1.plot(f, G33, label='Gain z')
    ax1.set_xlabel('rotation angle (deg)')
    ax1.set_ylabel('Gain')
    ax1.legend() 
    #plt.title(f"{mission}")

    ax2.plot(f, O1, label='Offset x')
    ax2.plot(f, O2, label='Offset y')
    ax2.plot(f, O3, label='Offset z')
    ax2.set_xlabel('rotation angle (deg)')
    ax2.set_ylabel('Offset')
    ax2.legend() 

    ax3.plot(f, G12, label='G12')
    ax3.plot(f, G13, label='G13')
    ax3.plot(f, G21, label='G21')
    ax3.plot(f, G23, label='G23')
    ax3.plot(f, G31, label='G31')
    ax3.plot(f, G32, label='G32')
    ax3.set_xlabel('rotation angle (deg)')
    ax3.set_ylabel('off-diagonal')
    ax3.legend() 

    filename = "AllPara_f" if filename is None else filename + "_AllPara_f"
    #plt.savefig(f"{filename}_{mission}") if mission != "" else plt.savefig(f"{filename}")
    plt.show()
    #plt.close()

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
    
    mission = "ela"
    #date = "20190430_183052" # ela long collection
    #date = "20190802_020149" # elb long collection
    date = "20190806_073926" # elb long collection
    rotAng_df = pd.read_csv(f"rotAng_loop_{mission}_360_{date}.csv")

    #Gain_f(rotAng_df['rotate_ang'], rotAng_df['G11'], rotAng_df['G22'], rotAng_df['G33'], mission = mission, filename = date)
    #STD_f(rotAng_df['rotate_ang'], rotAng_df['res_dmxl_x'], rotAng_df['res_dmxl_y'], rotAng_df['res_dmxl_z'], mission = mission, filename = date)
    Allpara_f(rotAng_df['rotate_ang'], rotAng_df['G11'], rotAng_df['G12'], rotAng_df['G13'], rotAng_df['O1'], 
        rotAng_df['G21'], rotAng_df['G22'], rotAng_df['G23'], rotAng_df['O2'],
        rotAng_df['G31'], rotAng_df['G32'], rotAng_df['G33'], rotAng_df['O3'], mission = mission, filename = date)
    breakpoint()