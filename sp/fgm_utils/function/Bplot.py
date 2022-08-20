from socketserver import DatagramRequestHandler
import matplotlib.pyplot as plt
from typing import List
import numpy as np
from .. import parameter

def phase_plot(ctime, phi, cross_time = None, gap_time = None, xlimt = None, title = "phase", datestr = None):

    filename = datestr + "_" + title if datestr is not None else title
    fig, ax = plt.subplots(1, figsize=(12,7))

    ax.plot(ctime, phi)
    ax.scatter(ctime, phi)
    if cross_time is not None: ax.scatter(cross_time, np.zeros(len(cross_time)), color='r')
    if gap_time is not None: 
        [ax.axvline(ts, color='r') for ts in gap_time]
    if xlimt is not None: ax.set_xlim(xlimt)

    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{filename}.png") 
    plt.close()

def ctimediff_plot(
    ctime, ctime_idx, ctime_idx_zoom = None, title = "ctimediff", datestr = None,
):

    filename = datestr + "_" + title if datestr is not None else title
    fig, ax = plt.subplots(1, figsize = (12, 7))
    ctime_adj = ctime[1:]-ctime[:-1]
    #ax.plot(ctime_adj, color='blue')
    ax.scatter(ctime[:-1], ctime_adj, color='blue')
    ax.scatter(ctime[ctime_idx], ctime_adj[ctime_idx], color='orange')
    ax.set_title('Differences between consecutive time steps')
    ax.set_xlabel('ctime')
    ax.set_ylabel('$t_{i+1} - t_i$')  
    if ctime_idx_zoom is not None:
        ax.set_xlim([ctime[ctime_idx_zoom]-5*2.8, ctime[ctime_idx_zoom]+5*2.8])
    ax.ticklabel_format(useOffset=False)
    
    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{filename}.png")  
    plt.close()

def B_ctime_plot(
    ctime: List[float], B_x: List[float],
    B_y: List[float], B_z: List[float],
    gap_time = None, plot3 = True, scatter = False,
    title = "B_ctime_plot", xlimt = None, cross_times = None, 
    datestr = None, ylimt = None,
):

    filename = datestr + "_" + title if datestr is not None else title
    if np.array(B_x).shape == np.array(B_y).shape == np.array(B_z).shape:
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

    if plot3 == True: # three subplots
        fig, ax = plt.subplots(3, figsize=(12,7))
        for i in range(dim):
            for j in range(3):
                ax[j].plot(ctime[i], B[j][i], label=labels[j][i], alpha=.5)
                ax[j].scatter(ctime[i], B[j][i], label=[], alpha=.5) if scatter == True else None
                if gap_time is not None:
                    [ax[j].axvline(k, color='r') for k in gap_time]
                ax[j].set_title(filename) if j == 0 else None
                ax[j].set_xlim(xlimt) if xlimt is not None else None
                ax[j].set_ylim(ylimt[j]) if ylimt is not None else None
                if cross_times is not None:
                    [ax[j].axvline(k, linestyle='--') for k in cross_times]
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
    title = "B_ctime_plot_single", xlimt = None, ylimt = None, 
    cross_times = None, datestr = None
):

    filename = datestr + "_" + title if datestr is not None else title
    dim = np.array(B).ndim
    fig, ax = plt.subplots(1, figsize=(12,7))
    if dim == 1:
        ax.plot(ctime, B, alpha=.5)
        ax.scatter(ctime, B, alpha=.5) if scatter == True else None
    else:
        [ax.plot(ctime, B[i], alpha=.5, label = f"X_{i}") for i in range(len(B))]
    if cross_times is not None:
        [ax.axvline(k, linestyle='--') for k in cross_times if cross_times is not None]
    ax.set_xlim(xlimt) if xlimt is not None else None
    ax.set_ylim(ylimt) if ylimt is not None else None
    ax.set_xlabel('Relative Time (seconds)')
    ax.set_ylabel('B (nT)')
    #ax.legend()
    plt.show() if parameter.savepng is False else plt.savefig(f"fgm_utils/temp/{filename}") 
    plt.close()
