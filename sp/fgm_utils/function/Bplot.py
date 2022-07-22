import matplotlib.pyplot as plt
from typing import List

def B_ctime_plot3(
    ctime: float, B_x: List[float],
    B_y: List[float], B_z: List[float],
    title: str
):

    fig, ax = plt.subplots(3, figsize=(10,5))

    ax[0].plot(ctime, B_x, label='x', alpha=.5)
    ax[0].scatter(ctime, B_x, label='', alpha=.5)
    ax[0].set_title(title)
    ax[0].set_xlabel('Relative Time (seconds)')
    ax[0].set_ylabel('X Field (nT)')
    ax[1].plot(ctime, B_y, label='y', alpha=.5)
    ax[1].scatter(ctime, B_y, label='', alpha=.5)
    ax[1].set_xlabel('Relative Time (seconds)')
    ax[1].set_ylabel('Y Field (nT)')
    ax[2].plot(ctime, B_z, label='z', alpha=.5)
    ax[2].scatter(ctime, B_z, label='', alpha=.5)
    ax[2].set_xlabel('Relative Time (seconds)')
    ax[2].set_ylabel('Z Field (nT)')

    plt.show()


def B_ctime_plot1(
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


def B2_ctime_plot(
    ctime: float, B_x: List[float],
    B_y: List[float], B_z: List[float],
    B_x2: List[float], B_y2: List[float],
    B_z2: List[float], title: str
):

    fig, ax = plt.subplots(3, figsize=(12,8))

    ax[0].plot(ctime, B_x, label='x1', alpha=.5)
    ax[0].set_title(title)
    ax[0].set_xlabel('Relative Time (seconds)')
    ax[0].set_ylabel('X Field (nT)')
    #ax.scatter(ctime, B_x, label='', alpha=.5)
    ax[0].plot(ctime, B_x2, label='x2', alpha=.5)
    ax[0].legend()
    ax[1].plot(ctime, B_y, label='y', alpha=.5)
    #ax.scatter(ctime, B_y, label='', alpha=.5)
    ax[1].plot(ctime, B_y2, label='y2', alpha=.5)
    ax[1].legend()
    ax[1].set_xlabel('Relative Time (seconds)')
    ax[1].set_ylabel('Y Field (nT)')
    ax[2].plot(ctime, B_z, label='z', alpha=.5)
    #ax.scatter(ctime, B_z, label='', alpha=.5)
    ax[2].plot(ctime, B_z2, label='z2', alpha=.5)
    ax[2].legend()
    ax[2].set_xlabel('Relative Time (seconds)')
    ax[2].set_ylabel('Z Field (nT)')

    plt.show() 

def B2_ctime_plot1(
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

    plt.show() 


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

    plt.show()        


