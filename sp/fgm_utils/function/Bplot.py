import matplotlib.pyplot as plt
from typing import List
import numpy as np
from .. import parameter

def phase_plot(ctime, phi, cross_time = None, ctime_idx = None, xlimt = None, title = "phase", datestr = None):

    filename = datestr + "_" + title if datestr is not None else title
    fig, ax = plt.subplots(1, figsize=(12,7))

    ax.plot(ctime, phi)
    ax.scatter(ctime, phi)
    if cross_time is not None: ax.scatter(cross_time, np.zeros(len(cross_time)), color='r')
    if ctime_idx is not None: 
        [ax.axvline(ctime[ts], color='r') for ts in ctime_idx]
    if xlimt is not None: ax.set_xlim(xlimt)

    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{filename}.png") 
    plt.close()


def ctimediff_plot(
    ctime, ctime_idx, ctime_idx_flag, ctime_idx_zoom = None, title = "ctimediff", datestr = None,
):

    filename = datestr + "_" + title if datestr is not None else title
    fig, ax = plt.subplots(1, figsize = (12, 7))
    ctime_adj = ctime[1:]-ctime[:-1]
    #ax.plot(ctime_adj, color='blue')
    ax.scatter(ctime[:-1], ctime_adj, color='blue')
    colors = ['red','orange','magenta','darkviolet','green']
    [ax.scatter(ctime[ctime_idx[i]], ctime_adj[ctime_idx[i]], color=colors[ctime_idx_flag[i]-1]) for i in range(len(ctime_idx))]
    ax.set_title('Differences between consecutive time steps')
    ax.set_xlabel('ctime')
    ax.set_ylabel('$t_{i+1} - t_i$')  
    if ctime_idx_zoom is not None:
        ax.set_xlim([ctime[ctime_idx_zoom]-5*2.8, ctime[ctime_idx_zoom]+5*2.8])
    ax.ticklabel_format(useOffset=False)
    
    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{filename}1.png") 
    ax.set_ylim([0, 0.2])
    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{filename}2.png")  
    plt.close()


def B_ctime_plot(
    ctime: List[float], B_x: List[float],
    B_y: List[float], B_z: List[float],
    ctime_idx_time = None, plot3 = True, scatter = False,
    title = "B_ctime_plot", xlimt = None, cross_times = None, 
    datestr = None, ylimt = None, ctime_idx_flag = None
):
    """plot function
    Parameter
        plot3: if true, plot 3 subplots
    """
    filename = datestr + "_" + title if datestr is not None else title
    if np.array(B_x).ndim == np.array(B_y).ndim == np.array(B_z).ndim:
        if np.array(B_x).ndim == np.array(ctime).ndim == 1:
            # one set of data and time
            dim = 1
            ctime = [ctime, []]
            B = [[B_x, []], [B_y, []],[B_z, []]]
            labels = [['x',''],['y',''],['z','']]
            y_labels = ['X Field (nT)', 'Y Field (nT)', 'Z Field (nT)']
        elif np.array(B_x).ndim == np.array(ctime).ndim == 2:
            # two sets of data and time
            dim = 2
            B = [B_x, B_y, B_z]
            labels = [['x1','x2'],['y1','y2'],['z1','z2']]
            y_labels = ['X Field (nT)', 'Y Field (nT)', 'Z Field (nT)']
        elif np.array(B_x).ndim == np.array(ctime).ndim + 1:
            # two sets of data and one set of time
            dim = 2
            ctime = [ctime, ctime]
            B = [B_x, B_y, B_z]
            labels = [['x1','x2'],['y1','y2'],['z1','z2']]
            y_labels = ['X Field (nT)', 'Y Field (nT)', 'Z Field (nT)']
    else:
        print("B_ctime_plot: not same length!")
        return 

    ctime_idx_time = [ctime_idx_time] if isinstance(ctime_idx_time,np.float64) else ctime_idx_time
    ctime_idx_time = np.array(ctime_idx_time) if isinstance(ctime_idx_time,list) else ctime_idx_time
    if plot3 == True: # three subplots
        fig, ax = plt.subplots(3, figsize=(12,7))
        for i in range(dim):
            for j in range(3):           
                ax[j].plot(ctime[i], B[j][i], label=labels[j][i], alpha=.5)
                ax[j].scatter(ctime[i], B[j][i], label=[], alpha=.5) if scatter == True else None
                if ctime_idx_time is not None and len(ctime_idx_time) != 0:
                    if ctime_idx_flag is not None: 
                        colors = ['red','orange','magenta','darkviolet','green']
                        for k in range(ctime_idx_time.size):
                            ax[j].axvline(ctime_idx_time[k], color=colors[ctime_idx_flag[k]-1]) 
                    else:
                        for k in range(len(ctime_idx_time)):
                            ax[j].axvline(ctime_idx_time[k], color='black') 
                ax[j].set_title(filename) if j == 0 else None
                ax[j].set_xlim(xlimt) if xlimt is not None else None
                ax[j].set_ylim(ylimt[j]) if ylimt is not None else None
                if cross_times is not None:
                    [ax[j].axvline(k, linestyle='--', color='black') for k in cross_times]
                ax[j].set_xlabel('Relative Time (seconds)')
                ax[j].set_ylabel(y_labels[j])
                ax[j].legend()
    else: # all in one plot
        fig, ax = plt.subplots(1, figsize=(12,7))
        for i in range(dim):
            for j in range(3):
                ax.plot(ctime[i], B[j][i], label=labels[j][i], alpha=.5)
        ax.set_title(filename)
        ax.set_xlim(xlimt) if xlimt is not None else None
        ax.set_ylim(ylimt) if ylimt is not None else None
        ax.set_xlabel('Relative Time (seconds)')
        ax.set_ylabel('Field (nT)')
        ax.legend()

    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{filename}") 
    plt.close()


def B_ctime_plot_single(
    ctime: List[float], B: List[float], scatter = False,
    title = "B_ctime_plot_single", xlimt = None, ylimt = None, fname = None,
    cross_times = None, datestr = None, legend = None, xlabel = None, ylabel=None,
):  
    """
    Plot time vs B in a single plot

    Parameter
        ctime: List of time
        B: List of B or list of B lists
        scatter: if true, plot scatter plot instead of line plot
        title: plot title, use together with datestr to save the plot
        datestr: date string, use together with title to save the plot
        xlimt: list of two numbers for xlimt
        ylimt: list of two numbers for ylimt
        cross_times: list of cross times to plot verticle lines
        legend: list of legend, must be same length as B      
    """
    title = datestr + "_" + title if datestr is not None else title
    filename = fname if fname is not None else title

    title = title + ""
    dim = np.array(B).ndim
    fig, ax = plt.subplots(1, figsize=(12,7))

    if dim == 1:
        ax.plot(ctime, B, alpha=.5)
        ax.scatter(ctime, B, alpha=.5) if scatter == True else None
    else:
        if legend is None:
            [ax.plot(ctime, B[i], alpha=.5, label = f"X_{i}") for i in range(len(B))]
        else:
            [ax.plot(ctime, B[i], alpha=.5, label = f"{legend[i]}") for i in range(len(B))]
    if cross_times is not None:
        [ax.axvline(k, linestyle='--') for k in cross_times if cross_times is not None]
    ax.set_xlim(xlimt) if xlimt is not None else None
    ax.set_ylim(ylimt) if ylimt is not None else None
    ax.set_xlabel(xlabel) if xlabel is not None else ax.set_xlabel('Relative Time (seconds)')
    ax.set_ylabel(ylabel) if ylabel is not None else ax.set_ylabel('B (nT)')
    ax.set_title(title)
    ax.legend() if legend is not None else None
    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{filename}") 
    plt.close()


def B_3d(B_x: List[float],B_y: List[float], B_z: List[float]):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(B_x, B_y, B_z, 'gray')
    plt.show()
    plt.close()


def omega_stage123(cross_times_1, w_syn_1, 
    cross_times_2=None,  w_syn_2=None, 
    cross_times_3=None,  w_syn_3=None,
    cross_times_fit=None,  w_syn_fit=None,
    title = "omega_stage_fit", datestr = None, xlimt = None, ylimt = None):
    filename = datestr + "_" + title if datestr is not None else title
    fig, ax = plt.subplots(1, figsize=(15,7))

    ax.scatter(cross_times_1, w_syn_1, label='Stage 1')
    ax.scatter(cross_times_2, w_syn_2, label='Stage 2') if cross_times_2 is not None else None
    ax.scatter(cross_times_3, w_syn_3, label='Stage 3') if cross_times_3 is not None else None
    ax.scatter(cross_times_fit, w_syn_fit, label='Fit') if cross_times_fit is not None else None
    ax.set_xlabel('Relative Time (s)')
    ax.set_ylabel(r'$\omega_{synodic}$')

    ax.legend() 
    if ylimt is not None: ax.set_ylim(ylimt)
    if xlimt is not None: ax.set_xlim(xlimt)
    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{filename}") 
    plt.close()


def omega_fit(cross_times_org, w_syn_org, 
    cross_times_fit = None,  w_syn_fit = None, 
    title = "omega_fit", datestr = None, xlimt = None, ylimt = None):
    """plot spin period in stage 1, 2, 3 and fitted
    """
    filename = datestr + "_" + title if datestr is not None else title
    fig, ax = plt.subplots(1, figsize=(15,7))

    ax.scatter(cross_times_org, w_syn_org, label='org')
    ax.scatter(cross_times_fit, w_syn_fit, label='fit') if cross_times_fit is not None else None
    ax.set_xlabel('Relative Time (s)')
    ax.set_ylabel(r'$\omega_{synodic}$')

    #ax.set_ylim(2.217, 2.2225)
    ax.legend() 
    if ylimt is not None: ax.set_ylim(ylimt)
    if xlimt is not None: ax.set_xlim(xlimt)
    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{filename}") 
    plt.close()


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
    plt.savefig(f"fgm_utils/temp/Gain_f_{mission}") if mission is not None else plt.savefig(f"fgm_utils/temp/Gain_f")
    plt.close()


def att_compare(att_cdfdata, column, att_cdfdata2 = None, datestr = None):
    fig, ax = plt.subplots(3, 1, figsize=(10, 15))

    df = att_cdfdata.copy()
    df['x'] = df[column].apply(lambda x: x[0])
    df['y'] = df[column].apply(lambda x: x[1])
    df['z'] = df[column].apply(lambda x: x[2])
    
    if att_cdfdata2 is not None:
        df2 =  att_cdfdata2.copy()
        df2['x'] = df2[column].apply(lambda x: x[0])
        df2['y'] = df2[column].apply(lambda x: x[1])
        df2['z'] = df2[column].apply(lambda x: x[2])

    ax[0].plot(df.index, df['x'], label='v03', linestyle='-', marker='*')
    ax[0].plot(df2.index, df2['x'], label='v02', linestyle='-', marker='*') if att_cdfdata2 is not None else None
    ax[0].set_title('x')
    ax[0].set_xlabel('time')
    ax[1].plot(df.index, df['y'], label='v03', linestyle='-', marker='*')
    ax[1].plot(df2.index, df2['y'], label='v02', linestyle='-', marker='*') if att_cdfdata2 is not None else None
    ax[1].set_title('y')
    ax[1].set_xlabel('time')
    ax[2].plot(df.index, df['z'], label='v03', linestyle='-', marker='*')
    ax[2].plot(df2.index, df2['z'], label='v02', linestyle='-', marker='*') if att_cdfdata2 is not None else None
    ax[2].set_title('z')
    ax[2].set_xlabel('time')
    plt.legend()

    def get_angle(vec1, vec2):
        dotprod = np.dot(vec1, vec2)
        norm_v1 = np.linalg.norm(vec1)
        norm_v2 = np.linalg.norm(vec2)
        cos_angle = dotprod / (norm_v1 * norm_v2)
        ang_rad = np.arccos(np.clip(cos_angle, -1.0, 1.0))
        return np.degrees(ang_rad)
    
    if att_cdfdata2 is not None:
        angles = df.apply(lambda row: get_angle(np.array(row[column]), np.array(df2.loc[row.name, column])), axis=1)
    print(f"Max angle difference: {np.max(angles)}")
    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{datestr}_att") 