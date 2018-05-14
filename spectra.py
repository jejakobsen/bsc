import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def get_utilde(filename):
    # The function reads utilde, u0 and t from \Calculated\utilde and return them in arrays
    # utilde, u and t. Which file should be read, has to specified by filename and should end with .dat
    os.chdir(r'C:\Users\jensj\OneDrive\Skrivebord\Beregninger\Calculated\utilde')
    txt = []; utilde = []; u = []; t = []
    infile = open(filename,'r')
    for line in infile:
        txt = line.split()
        utilde.append(float(txt[0]))
        u.append(float(txt[1]))
        t.append(float(txt[2]))
    infile.close()
    return(utilde,u,t)

def vizualise_u0(u0,u02,u03,u04,fs=10000):
    """
    The function takes in the different u0, and convert the u0s by using signal.welch
    such that they can be plotted in a energy cascade. The function then plots these converted
    values in a logarithmic plot.
    """
    f, Pxx_den   = signal.welch(u0, fs)
    f2, Pxx_den2 = signal.welch(u02,fs)
    f3, Pxx_den3 = signal.welch(u03,fs)
    f4, Pxx_den4 = signal.welch(u04,fs)
    plt.loglog(f, Pxx_den,'-r',f2,Pxx_den2,'-g',f3,Pxx_den3,'-k',f4,Pxx_den4,'-y')
    plt.legend(['Kommersiell probe','Probe 1','Probe 2','Probe 3'], fontsize = 18, loc = 3)
    plt.xlabel('f [Hz]', fontsize = 18)
    plt.ylabel('E(k)', fontsize = 18)
    plt.show()
    return('Done')

utilde_exp,u0_exp,t0_exp = get_utilde('utilde_expk.dat')
utilde_exp2,u0_exp2,t0_exp2 = get_utilde('utilde_exp2k.dat')
utilde_exp3,u0_exp3,t0_exp3 = get_utilde('utilde_exp3k.dat')
utilde_exp4,u0_exp4,t0_exp4 = get_utilde('utilde_exp4k.dat')


vizualise_u0(u0_exp,u0_exp2,u0_exp3,u0_exp4)

def lowpass_filter(u):
    """
    The function takes in one array of u0 and lets each value pass through the same low
    pass filter as mentioned in urms.py. The function plots how the low pass filter works,
    and also the effect of filtering represented in the energy cascade.
    """
    from scipy.signal import butter, lfilter, freqz
    def butter_lowpass(cutoff, fs, order):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        return b, a
    def butter_lowpass_filter(data, cutoff, fs, order):
        b, a = butter_lowpass(cutoff, fs, order)
        y = lfilter(b, a, data)
        return y
    order = 4
    fs = 10000.0       # sample rate, Hz
    cutoff = 1200.0  # desired cutoff frequency of the filter, Hz
    b, a = butter_lowpass(cutoff, fs, order) # Get the filter coefficients so we can check its frequency response
    w, h = freqz(b, a,worN=10000) # Plot the frequency response.
    plt.subplot(2, 1, 1)
    plt.axvline(cutoff, color='y')
    plt.plot(0.5*fs*w/np.pi, np.abs(h), 'g')
    plt.legend(['Cut-off frekvens for filter','%s. ordens filter' %str(order)], fontsize = 18)
    plt.plot(cutoff, 0.5*np.sqrt(2), 'yo')
    plt.xlim(0, 0.5*fs)
    plt.grid()
    plt.subplot(2, 1, 2)
    f2, Pxx_den2 = signal.welch(u,fs)
    y = butter_lowpass_filter(u, cutoff, fs, order)
    f3, Pxx_den3 = signal.welch(y,fs)
    plt.loglog(f2,Pxx_den2,'-r')
    plt.loglog(f3,Pxx_den3,'-g')
    plt.legend(['Probe 1: Ikke-filtrert','Probe 1: Filtrert'], fontsize = 18)
    plt.xlabel('f [Hz]',fontsize = 18)
    plt.ylabel('E(k)', fontsize = 18)
    plt.grid()
    plt.show()
    return('Done')

lowpass_filter(u0_exp2)
