# Calculate u_rms
import os
import numpy as np

def get_calibs(filename):
    """
    The function reads the calibration points from the folder \Calculated\Calibs and returns the
    two columns as two arrays (U_calib,V_calib).
    """
    os.chdir(r'C:\Users\jensj\OneDrive\Skrivebord\Beregninger\Calculated\Calibs')
    txt = []; U_calib = []; V_calib = [];
    infile = open(filename,'r')
    for line in infile:
        txt = line.split()
        U_calib.append(float(txt[0]))
        V_calib.append(float(txt[1]))
    infile.close()
    return(U_calib,V_calib)

def get_raw(directory,i):
    """
    The same function as the first function in meanwriter.py, but this function only returns
    the voltage and time.
    """
    V = []; t0 = [];
    os.chdir(directory)
    j = -1;
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        print(filename)
        if filename.endswith(".txt"):
            infile = open(filename,'r')
            infile.readline()
            infile.readline()
            txt = []; v = [];
            j = j + 1
            for line in infile:
                txt = line.split()
                v.append(float(txt[i]))
                if j == 0:
                    t0.append(float(txt[-1]))
        else:
            print(filename)
            continue
        infile.close()
        V.append(v);
    return(V,t0)

def comp_rms(list,U_calib,V_calib):
    """
    The function takes in the voltage list, and uses the U_calib and V_calib inputs
    in order to interpolate each value voltage value and convert it to a velocity.
    It then calculates and returns the mean velocity and rms velocity. The function also
    returns all velocities when the distance from the center of the pipe is zero (u0).
    It also subtracts the mean velocity from u0 and returns this value (utilde).
    """
    meanU = []; Urms = []
    V = np.array(list)
    i = 0
    for v in V:
        u = np.interp(v,V_calib,U_calib)
        meanU.append(np.mean(u))
        u1 = (u-meanU[i])
        if i == 0:
            utilde = u1
            u0 = u
        u2 = u1**2
        u3 = np.mean(u2)
        Urms.append(np.sqrt(u3))
        i = i + 1
    return(meanU,Urms,utilde,u0)

def comp_rms_fil(list,U_calib,V_calib):
    """
    The function does exactly the same as the comp_rms-function, but once
    the interpolation has happend all velocities pass through a low pass filter
    with cutoff 1200 and order 4. Credits to Petter Vollestad for helping me implement
    it. The function returns the same as comp_rms-function except utilde.
    """
    meanU_fil = []; Urms_fil = []
    V = np.array(list)
    def lowpass_filter(u):
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
        w, h = freqz(b, a) # Plot the frequency response.
        u_corr = butter_lowpass_filter(u, cutoff, fs, order)
        return(u_corr)
    i = 0
    for v in V:
        u = np.interp(v,V_calib,U_calib)
        u_fil = lowpass_filter(u)
        if i == 0:
            u0 = u_fil
        meanU_fil.append(np.mean(u_fil))
        u1_fil = (u_fil-meanU_fil[i])
        u2_fil = u1_fil**2
        u3_fil = np.mean(u2_fil)
        Urms_fil.append(np.sqrt(u3_fil))
        i = i + 1
    return(meanU_fil,Urms_fil,u0)

def write_rms(Urms,meanU,filename):
    """
    The function writes the inputted mean and rms velocity to a file in the folder
    Calculated\RMS. Filename also have to be given and should end with .dat.
    """
    os.chdir(r'C:\Users\jensj\OneDrive\Skrivebord\Beregninger\Calculated\RMS')
    outfile = open(filename,'w')
    for mu,urms in zip(meanU,Urms):
        outfile.write('%24s %25s \n' % (mu,urms))
    outfile.close()
    return('Done')

def write_utilde(utilde,u0,time,filename):
    #The function writes the inputted velocity and fluctuating velocity to a file in the folder
    #Calculated\utilde. Filename also have to be given and should end with .dat.
    os.chdir(r'C:\Users\jensj\OneDrive\Skrivebord\Beregninger\Calculated\utilde')
    outfile = open(filename,'w')
    for ut,u,t in zip(utilde,u0,time):
        outfile.write('%24s %24s %24s \n' % (ut,u,t))
    outfile.close()
    return('Done')

def read_comp_write(calibfile,rawpath,channel,outfile,outfile2):
    """
    The function runs all of the functions above by adding the calibration file, the path to the
    measurements-folder and which channel to read from. Also two outfile-names should be given.
    """
    U_calib, V_calib = get_calibs(calibfile)
    exp_dir = os.fsencode(rawpath)
    exp_data, t0 = get_raw(exp_dir,channel)
    meanU_exp,Urms_exp, utilde_exp, u0_exp = comp_rms(exp_data,U_calib,V_calib)
    meanU_exp_fil, Urms_exp_fil,u0 = comp_rms_fil(exp_data, U_calib,V_calib)
    write_rms(Urms_exp_fil,meanU_exp_fil,outfile)
    write_utilde(utilde_exp,u0_exp,t0,outfile2)
    return('Done')

read_comp_write('calib1.dat',r'C:\Users\jensj\Desktop\Raw\Maaling',1,'urms1.dat','utilde_exp.dat')
read_comp_write('calib3.dat',r'C:\Users\jensj\Desktop\Raw\Maaling2',0,'urms2.dat','utilde_exp2.dat')
read_comp_write('calib4.dat',r'C:\Users\jensj\Desktop\Raw\Maaling3',1,'urms3.dat','utilde_exp3.dat')
read_comp_write('calib5.dat',r'C:\Users\jensj\Desktop\Raw\Maaling4',0,'urms4.dat','utilde_exp4.dat')
