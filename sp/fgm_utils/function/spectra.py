import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import window_hanning, window_none
from scipy.signal import welch, windows, hann

def fgm_spectra(fgm_time: np.ndarray, fgm_z: np.ndarray):
    """calculate fgm spectra
    """
    # calculate sample frequency
    dt = np.mean(np.diff(fgm_time))
    fs = 1/dt

    window = np.hamming(len(fgm_z))
    fgm_z_windowed = fgm_z * window

    # calculate the FFT
    fft_result = np.fft.fft(fgm_z_windowed)

    fft_freq = np.fft.fftfreq(len(fft_result), 1/fs)

    psd = np.sqrt((np.abs(fft_result)**2) / len(fft_result))
    
    # Return the positive part of the spectrum
    mask = fft_freq > 0
    return fft_freq[mask], psd[mask]


def plot_fgm_spectra(fft_freq: np.ndarray, psd: np.ndarray):
    """plot spectra
    """
    fig, ax = plt.subplots()
    ax.plot(fft_freq, psd)
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('nT/sqrt(Hz)')
    #ax.set_xlim([0.1, 100])
    #ax.set_xticks([0.1, 1, 10, 100])
    #ax.set_ylim([1e-3, 1e1])
    #ax.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax.yaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax.set_title('ela')
    plt.show()


def plot_fgm_psd(fgm_time: np.ndarray, fgm_z: np.ndarray, title:str, ytitle=None, xlim=None, ylim=None):
    figure0, ax1 = plt.subplots()
    dt = np.mean(np.diff(fgm_time))
    fs = 1/dt
    Pxx, freqs = plt.psd(fgm_z, Fs=fs, NFFT=len(fgm_z), visible=False, window=window_hanning)
    #Pxx, freqs = plt.psd(fgm_z, Fs=fs, NFFT=len(fgm_z), visible=False, window=window_none)
    ax1.plot(freqs, np.sqrt(Pxx))
    ax1.set_xlabel('frequency [Hz]')
    ax1.set_ylabel('nT/sqrt(Hz)') if ytitle is None else ax1.set_ylabel(ytitle)
    ax1.set_xlim([0.05, 100]) if xlim is None else ax1.set_xlim(xlim)
    ax1.set_ylim([1e-1, 1.5e3]) if ylim is None else ax1.set_ylim(ylim)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.xaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax1.yaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax1.set_title(title)
    plt.show()
    

def plot_fgm_periodogram(fgm_time: np.ndarray, fgm_z: np.ndarray, title:str, ytitle=None, xlim=None, ylim=None):
    """
    Calculate FFT with scipy.signal.periodogram

    :param window: e.g. 'hann' or 'hamming'
    :param nfft: if none, length of x will be used
    :param scaling: 'spectrum' for power spectrum with unit of V**2, or 
                    'density' for power spectrum density with unit of V**2/Hz

    """
    from scipy.signal import periodogram
    dt = np.mean(np.diff(fgm_time))
    fs = 1/dt
    freqs, Pxx = periodogram(fgm_z, fs=fs, detrend=False, scaling='density', window='hann', nfft=len(fgm_z))
    figure0, ax1 = plt.subplots()
    ax1.plot(freqs, np.sqrt(Pxx))
    ax1.set_xlabel('frequency [Hz]')
    ax1.set_ylabel('nT/sqrt(Hz)') if ytitle is None else ax1.set_ylabel(ytitle)
    ax1.set_xlim([0.05, 100]) if xlim is None else ax1.set_xlim(xlim)
    ax1.set_ylim([1e-1, 1.5e3]) if ylim is None else ax1.set_ylim(ylim)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.xaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax1.yaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax1.set_title(title)
    plt.show()


def plot_fgm_spectrogram(
        fgm_time: np.ndarray,
          fgm_z: np.ndarray, 
          title:str, 
          ytitle=None, 
          xlim=None, 
          ylim=None,
          win='hann',
          nperseg=None,
          detr=None):
    """
    Calculate FFT with scipy.signal.spectrogram

    :param window: e.g. 'hann' or 'hamming'
    :param nfft: if none, length of x will be used
    :param scaling: 'spectrum' for power spectrum with unit of V**2, or 
                    'density' for power spectrum density with unit of V**2/Hz

    """
    from scipy.signal import spectrogram
    dt = np.mean(np.diff(fgm_time))
    fs = 1/dt
    nperseg = len(fgm_z) if nperseg is None else nperseg
    detrend = 'constant' if detr is None else detr
    window_functions = {
        'boxcar': windows.boxcar,
        'hann': windows.hann
    }
    if win in window_functions:
        window = window_functions[win](nperseg)
    freqs, times, Sxx = spectrogram(fgm_z, fs=fs, window=window, nperseg=nperseg, noverlap=nperseg // 2, nfft=nperseg, scaling='density', detrend=detrend)
    
    figure0, ax1 = plt.subplots()
    for time_idx in range(len(times)):
        ax1.plot(freqs, np.sqrt(Sxx[:, time_idx]))
    ax1.set_xlabel('frequency [Hz]')
    ax1.set_ylabel('nT/sqrt(Hz)') if ytitle is None else ax1.set_ylabel(ytitle)
    ax1.set_xlim([0.05, 100]) if xlim is None else ax1.set_xlim(xlim)
    ax1.set_ylim([1e-1, 1.5e3]) if ylim is None else ax1.set_ylim(ylim)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.xaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax1.yaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax1.set_title(title)
    plt.show()


def plot_fgm_welch(
        fgm_time: np.ndarray, 
        fgm_z: np.ndarray, 
        title:str, 
        ytitle=None, 
        xlim=None, 
        ylim=None, 
        win='hann',
        nperseg=None,
        detr=None,):
    """
    Calculate FFT with scipy.signal.periodogram

    :param window: e.g. 'hann' or 'hamming'
    :param nfft: if none, length of x will be used
    :param scaling: 'spectrum' for power spectrum with unit of V**2, or 
                    'density' for power spectrum density with unit of V**2/Hz

    """
    dt = np.mean(np.diff(fgm_time))
    fs = 1/dt
    nperseg = len(fgm_z) if nperseg is None else nperseg
    detrend = 'constant' if detr is None else detr
    window_functions = {
        'boxcar': windows.boxcar,
        'hann': windows.hann
    }
    if win in window_functions:
        window = window_functions[win](nperseg)
    freqs, psd = welch(fgm_z, fs, window=window, nperseg=nperseg, noverlap=nperseg // 2, nfft=nperseg, scaling='density', detrend=detrend)
    figure0, ax1 = plt.subplots()
    ax1.plot(freqs, np.sqrt(psd))
    ax1.set_xlabel('frequency [Hz]')
    ax1.set_ylabel('nT/sqrt(Hz)') if ytitle is None else ax1.set_ylabel(ytitle)
    ax1.set_xlim([0.05, 100]) if xlim is None else ax1.set_xlim(xlim)
    ax1.set_ylim([1e-1, 1.5e3]) if ylim is None else ax1.set_ylim(ylim)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.xaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax1.yaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax1.set_title(title)
    plt.show()

