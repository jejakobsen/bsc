import os
import numpy as np
from math import pi, sqrt, cos
import matplotlib.pyplot as plt

def mean_reader(filename):
    """
    The function reads and returns the mean voltage, mean flow and distance from center from a file
    in the \Calculated\Mean directory. The input file should be from the experiments when the distance is varying.
    End filename with .dat as always.
    """
    os.chdir(r'C:\Users\jensj\OneDrive\Skrivebord\Beregninger\Calculated\Mean')
    txt = []; meanV = []; meanM = []; y = []
    infile = open(filename,'r')
    for line in infile:
        txt = line.split()
        meanV.append(float(txt[0]))
        meanM.append(float(txt[1]))
        y.append(float(txt[2]))
    infile.close()
    return(meanV,meanM,y)

def calib_reader(filename):
    """
    The function read and returns the calibration points from a file in the Calculated\Calibs directory.
    The input filename should as always end with .dat.
    """
    os.chdir(r'C:\Users\jensj\OneDrive\Skrivebord\Beregninger\Calculated\Calibs')
    txt = []; U_cal = []; V_cal = [];
    infile = open(filename,'r')
    for line in infile:
        txt = line.split()
        U_cal.append(float(txt[0]))
        V_cal.append(float(txt[1]))
    infile.close()
    return(U_cal,V_cal)

def rms_reader(filename):
    """
    The function read and returns the mean and rms velocities from a file in the directory
    \Calculated\RMS. The input file should be from the experiments when the distance is varying.
    End filename with .dat as always.
    """
    os.chdir(r'C:\Users\jensj\OneDrive\Skrivebord\Beregninger\Calculated\RMS')
    txt = []; meanU = []; urms = [];
    infile = open(filename,'r')
    for line in infile:
        txt = line.split()
        meanU.append(float(txt[0]))
        urms.append(float(txt[1]))
    infile.close()
    return(meanU,urms)

def get_sim_data(filename):
    """
    The function reads and returns the mean and rms velocities from the DNS-file. The file was given from
    Anis Ayati, and processed by Petter Vollestad to a .txt-file.
    The function also read and returns the distance from the center in the pipe.
    """
    ls = []; us = []; Us = []; ys = []
    a = []; b = []; c = [];
    os.chdir(r'C:\Users\jensj\OneDrive\Skrivebord\Beregninger\Simulering')
    infile = open(filename,'r')
    for line in infile:
        ls.append(line)
    us, Us, ys = ls[3:305], ls[310:612], ls[619:921]
    for i in range(301):
        a.append(float(us[i]))
        b.append(float(Us[i]))
        c.append(float(ys[i]))
    usn,Usn,ysn = a,b,c
    return(usn,Usn,ysn)

# Commercial probe
meanV,meanM,y = mean_reader('exp_mean.dat')
U_calib,V_calib = calib_reader('calib1.dat')
meanU, urms = rms_reader('urms1.dat')

# Hand made probes
meanV2,meanM2,y2 = mean_reader('exp_mean2.dat')
U_calib2,V_calib2 = calib_reader('calib3.dat')
meanU2, urms2 = rms_reader('urms2.dat')


meanV3,meanM3,y3 = mean_reader('exp_mean3.dat')
U_calib3,V_calib3 = calib_reader('calib4.dat')
meanU3, urms3 = rms_reader('urms3.dat')

meanV4,meanM4,y4 = mean_reader('exp_mean4.dat')
U_calib4,V_calib4 = calib_reader('calib5.dat')
meanU4, urms4 = rms_reader('urms4.dat')

# Simulations
usim, Usim, ysim = get_sim_data('re44k.txt')

def meanplot1(meanV1,meanM1,y1,meanV2,meanM2,y2,y_just):
    """
    The function plots two subplots of the DNS',commercial probe' and probe 1's mean velocity profile.
    A correction factor is being multiplied in order to standardize the height measurements. Also, the bulk
    velocity (Ub) is used to standardize U. Another correction factor of 2mm is added in order to
    'center' probe 1. This makes up the second subplot, where the first subplot is without this correction.
    """
    def comp_UbRe(meanMn):
        """
        Same function as in calib_writer
        """
        meanMnp = np.array(meanMn)
        rho = 1.19; r = 0.05; nu = 15.0*10**(-6)
        A = pi*(r)**2; L = 2*r
        Ub = meanMnp/(rho*A); Re = (Ub*L)/nu
        return(Ub,Re)
    corr_fac = 10**(-3)/0.05
    Ub1,Re1 = comp_UbRe(meanM1)
    Ub2,Re2 = comp_UbRe(meanM2)
    U_std1 = meanV1/Ub1
    U_std2 = meanV2/Ub2
    ynp1 = np.array(y1);
    ynp2 = np.array(y2); ynp4 = np.array(y2)+y_just
    y_std1 = ynp1*corr_fac;  y_std2 = ynp2*corr_fac; y_std4 = ynp4*corr_fac
    plt.subplot(2,2,1)
    plt.plot(U_std1,y_std1,'xr',U_std2,y_std2,'og',Usim,ysim,'m',Usim,-np.array(ysim),'m')
    plt.grid('on')
    plt.ylabel('y/r',fontsize=18)
    plt.ylim(-1.5,1.5)
    plt.legend(['Kommersiell probe','Probe 1','DNS'],fontsize=16)
    plt.subplot(2,2,3)
    plt.plot(U_std1,y_std1,'xr',U_std2,y_std4,'og',Usim,ysim,'m',Usim,-np.array(ysim),'m')
    plt.grid('on')
    plt.ylabel('y/r',fontsize=18)
    plt.xlabel('U/Ub',fontsize=18)
    plt.ylim(-1.5,1.5)
    plt.legend(['Kommersiell probe','Probe 1, y+2mm','DNS'],fontsize=16)
    return('Done')

def rmsplot1(Urms1,y1,Urms2,y2):
    """
    The function plots two subplots of the DNS',commercial probe' and probe 1's rms velocity profile.
    A correction factor is being multiplied in order to standardize the height measurements. Also, the
    rms velocity is standardized by utau. Another correction factor of 2mm is added in order to
    'center' probe 1. This makes up the second subplot, where the first subplot is without this correction.
    """
    corr_fac = 10**(-3)/0.05
    Utau = 0.36821
    from math import cos
    U_rms1 = np.array(Urms1)/Utau
    U_rms2 = np.array(Urms2)/Utau
    y_std1 = np.array(y1)*corr_fac
    y_std2 = np.array(y2)*corr_fac; y_std4 = (np.array(y2)+2.0)*corr_fac
    usimnp = np.array(usim)
    f = sqrt(0.002695)
    plt.subplot(2,2,2)
    plt.plot(U_rms1,y_std1,'xr',U_rms2,y_std2,'og',usimnp/f,np.array(ysim),'m',usimnp/f,-np.array(ysim),'m')
    plt.grid('on')
    plt.ylim(-1.5,1.5)
    plt.subplot(2,2,4)
    plt.plot(U_rms1,y_std1,'xr',U_rms2,y_std4,'og',usimnp/f,np.array(ysim),'m',usimnp/f,-np.array(ysim),'m')
    plt.grid('on')
    plt.xlabel('urms/utau',fontsize=18)
    plt.ylim(-1.5,1.5)
    return('Done')

def meanplot2(meanV1,meanM1,y1,meanV2,meanM2,y2,meanV3,meanM3,y3,y_just1,y_just2,y_just3,fs):
    """
    The function plots two subplots of the DNS and probe 1,2,3's mean velocity profile.
    A correction factor is being multiplied in order to standardize the height measurements. Also, the bulk
    velocity (Ub) is used to standardize U. Another correction factor of 2mm is added in order to
    'center' probe 1,2,3. This makes up the second subplot, where the first subplot is without this correction.
    """
    def comp_UbRe(meanMn):
        meanMnp = np.array(meanMn)
        rho = 1.19; r = 0.05; nu = 15.0*10**(-6)
        A = pi*(r)**2; L = 2*r
        Ub = meanMnp/(rho*A); Re = (Ub*L)/nu
        return(Ub,Re)
    corr_fac = 10**(-3)/0.05
    Ub1,Re1 = comp_UbRe(meanM1)
    Ub2,Re2 = comp_UbRe(meanM2)
    Ub3,Re3 = comp_UbRe(meanM3)
    U_std1 = meanV1/Ub1
    U_std2 = meanV2/Ub2
    U_std3 = meanV3/Ub3
    ynp1 = np.array(y1); ynp11 = np.array(y1)+y_just1
    ynp2 = np.array(y2); ynp21 = np.array(y2)+y_just2
    ynp3 = np.array(y3); ynp31 = np.array(y3)+y_just3
    y_std1 = ynp1*corr_fac;  y_std2 = ynp2*corr_fac; y_std3=ynp3*corr_fac; y_std4 = ynp11*corr_fac; y_std5 = ynp21*corr_fac;  y_std6 = ynp31*corr_fac
    plt.subplot(2,2,1)
    plt.plot(U_std1,y_std1,'.g',U_std2,y_std2,'+k',U_std3,y_std3,'+y',Usim,ysim,'m',Usim,-np.array(ysim),'m')
    plt.grid('on')
    plt.ylabel('y/r',fontsize=18)
    plt.ylim(-1.5,1.5)
    plt.legend(['Probe 1','Probe 2','Probe 3','DNS'],fontsize=fs)
    plt.subplot(2,2,3)
    plt.plot(U_std1,y_std4,'.g',U_std2,y_std5,'+k',U_std3,y_std6,'+y',Usim,ysim,'m',Usim,-np.array(ysim),'m')
    plt.grid('on')
    plt.ylabel('y/r',fontsize=18)
    plt.xlabel('U/Ub',fontsize=18)
    plt.ylim(-1.5,1.5)
    plt.legend(['Probe 1, y+2mm','Probe 2, y+2mm','Probe 3, y+2mm','DNS'],fontsize=fs)
    return('Done')

def rmsplot2(Urms1,y1,Urms2,y2,Urms3,y3,y_just1,y_just2,y_just3):
    """
    The function plots two subplots of the DNS' and probe 1,2,3's rms velocity profile.
    A correction factor is being multiplied in order to standardize the height measurements. Also, the
    rms velocity is standardized by utau. Another correction factor of 2mm is added in order to
    'center' probe 1. This makes up the second subplot, where the first subplot is without this correction.
    """
    corr_fac = 10**(-3)/0.05
    Utau = 0.36821
    from math import cos
    U_rms1 = np.array(Urms1)/Utau
    U_rms2 = np.array(Urms2)/Utau
    U_rms3 = np.array(Urms3)/Utau
    ynp1 = np.array(y1); ynp11 = np.array(y1)+y_just1
    ynp2 = np.array(y2); ynp21 = np.array(y2)+y_just2
    ynp3 = np.array(y3); ynp31 = np.array(y3)+y_just3
    y_std1 = ynp1*corr_fac;  y_std2 = ynp2*corr_fac; y_std3=ynp3*corr_fac; y_std4 = ynp11*corr_fac; y_std5 = ynp21*corr_fac;  y_std6 = ynp31*corr_fac
    usimnp = np.array(usim)
    f = sqrt(0.002695)
    plt.subplot(2,2,2)
    plt.plot(U_rms1,y_std1,'.g',U_rms2,y_std2,'+k',U_rms3,y_std3,'+y',usimnp/f,np.array(ysim),'m',usimnp/f,-np.array(ysim),'m')
    plt.grid('on')
    plt.ylim(-1.5,1.5)
    plt.subplot(2,2,4)
    plt.plot(U_rms1,y_std4,'.g',U_rms2,y_std5,'+k',U_rms3,y_std6,'+y',usimnp/f,np.array(ysim),'m',usimnp/f,-np.array(ysim),'m')
    plt.grid('on')
    plt.xlabel('urms/utau',fontsize=18)
    plt.ylim(-1.5,1.5)
    return('Done')

meanplot1(meanU,meanM,y,meanU2,meanM2,y2,2.0)
rmsplot1(urms,y,urms2,y2)
plt.show()

meanplot2(meanU2,meanM2,y2,meanU3,meanM3,y3,meanU4,meanM4,y4,2.0,2.0,2.0,14)
rmsplot2(urms2,y2,urms3,y3,urms4,y4,2.0,2.0,2.0)
plt.show()


def comp_UbRe(meanMn):
    """
    Same function as in calib_writer
    """
    meanMnp = np.array(meanMn)
    rho = 1.19; r = 0.05; nu = 15.0*10**(-6)
    A = pi*(r)**2; L = 2*r
    Ub = meanMnp/(rho*A); Re = (Ub*L)/nu
    return(Ub,Re)
print(meanU[0],meanU2[0],meanU3[0],meanU4[0]) # Prints mean velocities at y=0.
print(comp_UbRe(meanM[0]),comp_UbRe(meanM2[0]),comp_UbRe(meanM3[0]),comp_UbRe(meanM4[0])) # Prints Reynolds number and bulk speed at y=0.

"""
8.772593863366643 8.35141708413993 8.471884068044336 8.4814954556256
(6.867000536175182, 45780.003574501214) (6.8040107933201215, 45360.07195546749) (6.84428770142332, 45628.58467615547) (6.856233080448134, 45708.2205363209)
"""
