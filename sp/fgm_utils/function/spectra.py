import numpy as np
import matplotlib.pyplot as plt

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
    figure, ax1 = plt.subplots()
    dt = np.mean(np.diff(fgm_time))
    fs = 1/dt
    Pxx, freqs = plt.psd(fgm_z, Fs=fs, NFFT=2048, visible=False)
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
    
    # figure, (ax1, ax2) = plt.subplots(nrows=2)
    # dt = np.mean(np.diff(fgm_time))
    # fs = 1/dt
    # Pxx, freqs = plt.psd(fgm_z, Fs=fs, NFFT=2048, visible=False)
    # ax1.plot(freqs, np.sqrt(Pxx))
    # ax1.set_xlabel('frequency [Hz]')
    # ax1.set_ylabel('PSD')
    # ax1.set_xlim([0.1, 100])
    # ax1.set_xticks([0.1, 1, 10, 100])
    # ax1.set_ylim([1e-3, 1e1])
    # ax1.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')
    # ax1.xaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    # ax1.yaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    # ax1.set_title('ela sqrt')

    # Pxx, freqs = plt.psd(fgm_z, Fs=fs, NFFT=1024, visible=False)
    # ax2.plot(freqs, Pxx)
    # ax2.set_xlabel('frequency [Hz]')
    # ax2.set_ylabel('PSD')
    # ax2.set_xlim([0.1, 100])
    # ax2.set_xticks([0.1, 1, 10, 100])
    # ax2.set_ylim([1e-3, 1e1])
    # ax2.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
    # ax2.set_xscale('log')
    # ax2.set_yscale('log')
    # ax2.xaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    # ax2.yaxis.grid(True, which='major', linestyle='--', linewidth=0.5)
    # ax2.set_title('ela')
    # plt.show()
