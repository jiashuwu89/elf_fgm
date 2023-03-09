import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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


def Gthphi_f(G11, G12, G13, G21, G22, G23, G31, G32, G33):

    G1 = (G11**2 + G12**2 + G13**2)**0.5 
    G2 = (G21**2 + G22**2 + G23**2)**0.5
    G3 = (G31**2 + G32**2 + G33**2)**0.5
    
    th1 = np.degrees(np.arccos(G13/G1))
    th2 = np.degrees(np.arccos(G23/G2))
    th3 = np.degrees(np.arccos(G33/G3))

    ph1 = np.degrees(np.arctan(G12/G11))
    ph2 = np.degrees(np.arctan(G22/G21))
    ph3 = np.degrees(np.arctan(G32/G31))

    return G1, G2, G3, th1, th2, th3, ph1, ph2, ph3


def Gthphi_f_plot(f, G11, G12, G13, O1, G21, G22, G23, O2, G31, G32, G33, O3, mission = "", filename = None):
    """plot Gain x, y, z as a function of rotation angle
    """
    f_len = len(f)
    G1 = np.zeros(f_len)
    G2 = np.zeros(f_len)
    G3 = np.zeros(f_len)
    th1 = np.zeros(f_len)
    th2 = np.zeros(f_len)
    th3 = np.zeros(f_len)
    ph1 = np.zeros(f_len)
    ph2 = np.zeros(f_len)
    ph3 = np.zeros(f_len)
    
    for i, _ in enumerate(f):
        G1[i], G2[i], G3[i], th1[i], th2[i], th3[i], ph1[i], ph2[i], ph3[i] =  Gthphi_f(
            G11[i], G12[i], G13[i], G21[i], G22[i], G23[i], G31[i], G32[i], G33[i])
    
    idx_0 = [i for i in range(len(f)-1) if ph2[i] > 0 and ph2[i+1] < 0][0] + 1
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(5,8))

    ax1.plot(f, G1, label='G1')
    ax1.plot(f, G2, label='G2')
    ax1.plot(f, G3, label='G3')
    ax1.set_xlabel('rotation angle (deg)')
    ax1.set_ylabel('|G|')
    ax1.legend() 
    ax1.set_title(f'G1: {round(np.average(G1), 1)}, G2: {round(np.average(G2), 1)}, G3: {round(np.average(G3), 1)}', fontsize=10)

    ax2.plot(f, th1, label='th1')
    ax2.plot(f, th2, label='th2')
    ax2.plot(f, th3, label='th3')
    ax2.set_xlabel('rotation angle (deg)')
    ax2.set_ylabel('elevation angle (deg)')
    ax2.legend()
    ax2.set_title(f'th1: {round(np.average(th1), 1)}, th2: {round(np.average(th2), 1)}, th3: {round(np.average(th3), 1)}', fontsize=10) 

    ax3.plot(f, ph1, label='ph1')
    ax3.plot(f, ph2, label='ph2')
    ax3.plot(f, ph3, label='ph3')
    ax3.set_xlabel('rotation angle (deg)')
    ax3.set_ylabel('azimuthal angle (deg)')
    ax3.legend(loc='upper right') 
    ax3.axvline(x=f[idx_0], color='red', linestyle='--')
    ax3.set_title(f'rotAng = {f[idx_0]} deg ph1 = {round(ph1[idx_0], 1)} deg when ph2 = 90 deg')

    fig.subplots_adjust(hspace=0.4, wspace=0.4) # inch

    filename = "Gthph_f" if filename is None else filename + "_Gthph_f"
    plt.savefig(f"{filename}_{mission}") if mission != "" else plt.savefig(f"{filename}")
    #plt.show()
    #breakpoint()
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
    
    mission = "ela"
    date = "20190430_183052" # ela long collection
    #date = "20190802_020149" # elb long collection
    #date = "20190806_073926" # elb long collection
    rotAng_df = pd.read_csv(f"rotAng_loop_{mission}_360_{date}.csv")

    #Gain_f(rotAng_df['rotate_ang'], rotAng_df['G11'], rotAng_df['G22'], rotAng_df['G33'], mission = mission, filename = date)
    #STD_f(rotAng_df['rotate_ang'], rotAng_df['res_dmxl_x'], rotAng_df['res_dmxl_y'], rotAng_df['res_dmxl_z'], mission = mission, filename = date)
    Allpara_f(rotAng_df['rotate_ang'], rotAng_df['G11'], rotAng_df['G12'], rotAng_df['G13'], rotAng_df['O1'], 
        rotAng_df['G21'], rotAng_df['G22'], rotAng_df['G23'], rotAng_df['O2'],
        rotAng_df['G31'], rotAng_df['G32'], rotAng_df['G33'], rotAng_df['O3'], mission = mission, filename = date)
        
    Gthphi_f_plot(rotAng_df['rotate_ang'], rotAng_df['G11'], rotAng_df['G12'], rotAng_df['G13'], rotAng_df['O1'], 
        rotAng_df['G21'], rotAng_df['G22'], rotAng_df['G23'], rotAng_df['O2'],
        rotAng_df['G31'], rotAng_df['G32'], rotAng_df['G33'], rotAng_df['O3'], mission = mission, filename = date)
    #breakpoint()