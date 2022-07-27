import matplotlib.pyplot as plt
from typing import List
import numpy as np
from .. import parameter

def phase_plot(ctime, phi, cross_time = None, err = None, xlimt = None):

    fig, ax = plt.subplots(1, figsize=(10,5))

    ax.plot(ctime, phi)
    ax.scatter(ctime, phi)
    if cross_time is not None: ax.scatter(cross_time, np.zeros(len(cross_time)), color='r')
    if err is not None: ax.axvline(err.all())
    if xlimt is not None: ax.set_xlim(xlimt)

    plt.show() if parameter.savepng is False else plt.savefig("fgm_utils/temp/phase_plot.png") 


def B1_ctime1_plot3(
    ctime: float, B_x: List[float],
    B_y: List[float], B_z: List[float],
    err: List[float] = None,
    title = "", xlimt = None, cross_times = None
):

    fig, ax = plt.subplots(3, figsize=(12,7))

    ax[0].plot(ctime, B_x, label='x', alpha=.5)
    #ax[0].scatter(ctime, B_x, label='', alpha=.5)
    ax[0].set_title(title)
    ax[0].set_xlabel('Relative Time (seconds)')
    ax[0].set_ylabel('X Field (nT)')
    ax[1].plot(ctime, B_y, label='y', alpha=.5)
    #ax[1].scatter(ctime, B_y, label='', alpha=.5)
    ax[1].set_xlabel('Relative Time (seconds)')
    ax[1].set_ylabel('Y Field (nT)')
    ax[2].plot(ctime, B_z, label='z', alpha=.5)
    #ax[2].scatter(ctime, B_z, label='', alpha=.5)
    ax[2].set_xlabel('Relative Time (seconds)')
    ax[2].set_ylabel('Z Field (nT)')

    if err is not None:
        ax[0].axvline(err.all())
        ax[1].axvline(err.all())
        ax[2].axvline(err.all())

    if xlimt is not None:
        ax[0].set_xlim(xlimt)
        ax[1].set_xlim(xlimt)
        ax[2].set_xlim(xlimt)

    if cross_times is not None:
        ax[0].scatter(cross_times, np.zeros(len(cross_times)), color='r')
        ax[1].scatter(cross_times, np.zeros(len(cross_times)), color='r')
        ax[2].scatter(cross_times, np.zeros(len(cross_times)), color='r')    
    
    plt.show() if parameter.savepng is False else plt.savefig("fgm_utils/temp/B1_ctime1_plot3.png") 


def B1_ctime1_plot1(
    ctime: float, B_x: List[float],
    B_y: List[float], B_z: List[float],
    title: str
):

    fig, ax = plt.subplots(1, figsize=(10,5))

    ax.plot(ctime, B_x, label='x', alpha=.5)
    ax.scatter(ctime, B_x, label='', alpha=.5)
    ax.set_title(title)
    ax.plot(ctime, B_y, label='y', alpha=.5)
    ax.scatter(ctime, B_y, label='', alpha=.5)
    ax.plot(ctime, B_z, label='z', alpha=.5)
    ax.scatter(ctime, B_z, label='', alpha=.5)
    ax.set_xlabel('Relative Time (seconds)')
    ax.set_ylabel('Field (nT)')
    ax.legend()

    plt.show()
    #plt.savefig("fgm_utils/temp/B1_ctime1_plot1.png")


def B2_ctime1_plot3(
    ctime: float, B_x: List[float],
    B_y: List[float], B_z: List[float],
    B_x2: List[float], B_y2: List[float],
    B_z2: List[float], err: List[float] = None, 
    title: str = "", xlimt = None, cross_times_corr = None
):

    fig, ax = plt.subplots(3, figsize=(12,8))

    ax[0].plot(ctime, B_x, label='x1', alpha=.5)
    ax[0].scatter(ctime, B_x, label='')
    ax[0].set_title(title)
    ax[0].set_xlabel('Relative Time (seconds)')
    ax[0].set_ylabel('X Field (nT)')
    #ax.scatter(ctime, B_x, label='', alpha=.5)
    ax[0].plot(ctime, B_x2, label='x2', alpha=.5)
    ax[0].legend()
    ax[1].plot(ctime, B_y, label='y', alpha=.5)
    #ax.scatter(ctime, B_y, label='', alpha=.5)
    ax[1].plot(ctime, B_y2, label='y2', alpha=.5)
    ax[1].scatter(ctime, B_y, label='')
    ax[1].legend()
    ax[1].set_xlabel('Relative Time (seconds)')
    ax[1].set_ylabel('Y Field (nT)')
    ax[2].plot(ctime, B_z, label='z', alpha=.5)
    #ax.scatter(ctime, B_z, label='', alpha=.5)
    ax[2].plot(ctime, B_z2, label='z2', alpha=.5)
    ax[2].scatter(ctime, B_z, label='')
    ax[2].legend()
    ax[2].set_xlabel('Relative Time (seconds)')
    ax[2].set_ylabel('Z Field (nT)')

    #plt.savefig("fgm_utils/temp/B2_ctime1_plot3.png")
    if err is not None:
        for i in range(len(err)):
            ax[0].axvline(err[i])
            ax[1].axvline(err[i])
            ax[2].axvline(err[i])

    if xlimt is not None:
        ax[0].set_xlim(xlimt)
        ax[1].set_xlim(xlimt)
        ax[2].set_xlim(xlimt)

    if cross_times_corr is not None:
        ax[0].scatter(cross_times_corr, np.zeros(len(cross_times_corr)), color='r')
        ax[1].scatter(cross_times_corr, np.zeros(len(cross_times_corr)), color='r')
        ax[2].scatter(cross_times_corr, np.zeros(len(cross_times_corr)), color='r')     

    plt.show()


def B2_ctime1_plot1(
    ctime: float, B_x: List[float],
    B_y: List[float], B_z: List[float],
    B_x2: List[float], B_y2: List[float],
    B_z2: List[float], title: str
):

    fig, ax = plt.subplots(1, figsize=(12,8))

    ax.plot(ctime, B_x, label='x1', alpha=.5)
    ax.set_title(title)
    #ax.scatter(ctime, B_x, label='', alpha=.5)
    ax.plot(ctime, B_x2, label='x2', alpha=.5)

    ax.plot(ctime, B_y, label='y', alpha=.5)
    #ax.scatter(ctime, B_y, label='', alpha=.5)
    ax.plot(ctime, B_y2, label='y2', alpha=.5)

    ax.plot(ctime, B_z, label='z', alpha=.5)
    #ax.scatter(ctime, B_z, label='', alpha=.5)
    ax.plot(ctime, B_z2, label='z2', alpha=.5)
    ax.legend()
    ax.set_xlabel('Relative Time (seconds)')
    ax.set_ylabel('Field (nT)')

    plt.savefig("fgm_utils/temp/B2_ctime1_plot1.png")


def B2_ctime2_plot3(
    ctime: float, B_x: List[float],
    B_y: List[float], B_z: List[float],
    ctime2: float, B_x2: List[float], 
    B_y2: List[float], B_z2: List[float], title: str
):

    fig, ax = plt.subplots(3, figsize=(12,8))

    ax[0].plot(ctime, B_x, label='x', alpha=.5, color='r')
    ax[0].scatter(ctime2, B_x2, label='x2', alpha=.5, color='r')
    ax[0].legend()
    ax[0].set_title(title)
    ax[0].set_xlabel('Relative Time (seconds)')
    ax[0].set_ylabel('Field (nT)')
    ax[1].plot(ctime, B_y, label='y', alpha=.5, color='b')
    ax[1].scatter(ctime2, B_y2, label='y2', alpha=.5, color='b')
    ax[1].legend()
    ax[1].set_xlabel('Relative Time (seconds)')
    ax[1].set_ylabel('Field (nT)')
    ax[2].plot(ctime, B_z, label='z', alpha=.5, color='g')
    ax[2].scatter(ctime2, B_z2, label='z2', alpha=.5, color='g')
    ax[2].legend()
    ax[2].set_xlabel('Relative Time (seconds)')
    ax[2].set_ylabel('Field (nT)')

    plt.savefig("temp/B2_ctime2_plot3.png")   


def ctimediff_plot(ctime, ctime_idx):
    fig, ax = plt.subplots(1, figsize = (12, 7))
    ctime_adj = ctime[1:]-ctime[:-1]
    ax.plot(ctime_adj, color='blue')
    ax.scatter(ctime_idx, ctime_adj[ctime_idx], color='orange')
    ax.set_title('Differences between consecutive time steps')
    ax.set_xlabel('index')
    ax.set_ylabel('$t_{i+1} - t_i$')  

    plt.savefig("fgm_utils/temp/ctimediff_plot.png")       


