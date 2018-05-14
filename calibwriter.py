# Note: At the bottom of the script there are four lines in bulks of two. This has to be run in sequence in
# order for the script to work.

import os
import numpy as np
from math import pi,exp
import matplotlib.pyplot as plt

def read_mean(filename):
    """
    The function read filename from the \Calculated\Mean folder and returns mean
    voltage, mass flow and distance from center.
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

def comp_UbRe(meanM):
    """
    The function calculate and returns bulk velocity (Ub) and Reynolds number (Re) with
    the help of the inputted mean mass flow.
    """
    meanMnp = np.array(meanM)
    rho = 1.19; r = 0.05; nu = 1.5*10**(-5)
    A = pi*(r)**2; L = 2*r
    Ub = meanMnp/(rho*A); Re = (Ub*L)/nu
    return(Ub,Re)

def Nikuradse(Ub,Re):
    """
    Credits to Anis Ayati. The function interpolates the inputted Reynolds number (Re)
    and retrieves the corresponding n value from a curve. This curve represents
    the relationship between Re and n. It then calculates and returns the center velocity (U_center)
    with the help of n and inputted the bulk velocity.
    """
    from scipy.interpolate import CubicSpline
    R = np.array([4e3,2.3e4,1.1e5,1.1e6,2.0e6,3.2e6])
    exp = np.array([6.0,6.6,7,8.8,10,10])
    from scipy.interpolate import CubicSpline
    cs = CubicSpline(R,exp)
    n = cs(Re)
    U_center = Ub*(n+1)*((2.0*n)+1)/(2.0*n*n)
    return(U_center)

def calib_poly(U_c,meanV,order):
    """
    This function finds a polynomial fit between the inputted points (U_c,meanV). The order of the
    polynomial is also an input. The function returns points who describe the polynomial (U_calib,V_calib), as well
    as the function (p) and the inverse function (pp).
    """
    U_c = np.sort(U_c); meanV = np.sort(meanV)
    U_c[0] = 0.0
    p = np.polyfit(U_c,meanV,order)
    pp = np.polyfit(meanV,U_c,order)
    U_calib = np.linspace(0.0,11.0,1001)
    V_calib = np.polyval(p,U_calib)
    return(U_calib,V_calib,p,pp)

def calib_kings(U_c,meanV,f,max):
    """
    The function finds the best values for A and B from Kings law. It uses the inputted
    points (U_c,meanV) and the least square metod in order to do so. f is the exponential number "n" in Kings law and is also a input value.
    Then the function calculate and returns points (U_calib,V_calib) from the squared Kings law. It also returns
    A,B and the squared Kings law.
    """
    U_c = np.array(U_c,dtype=float); meanV = np.array(np.power(meanV,2),dtype=float)
    U_c = np.sort(U_c); meanV = np.sort(meanV);
    Unew = U_c**(f)
    Unew[0] = 0.0
    from scipy.optimize import curve_fit
    def func(U, A, B):
        return(A+B*U)
    popt, pcov = curve_fit(func,Unew,meanV)
    def func1(U):
        return(np.sqrt(popt[0]+popt[1]*U**f))
    U_calib = np.linspace(0,max,max*1000+1)
    V_calib = func1(U_calib)
    return(U_calib,V_calib,popt,func1)

def piecewise_calib(p,f,max,limit):
    """
    This function takes the polynomial and kings law function-expressions, and creates
    a piecewise function composed of points. limit sets the limit for where
    the points stop following the polynomial function, and Kings law begin. max sets the upper
    boundary for the piecewise function. The function returns the points that makes up the piecewise function.
    """
    u1 = np.linspace(0,limit,limit*1000+1); u2 = np.linspace(limit,max,limit*1000+1)
    v1 = np.polyval(p,u1); v2 = f(u2)
    U_pw = np.append(u1,u2)
    V_pw = np.append(v1,v2)
    return(U_pw,V_pw)

def visualize_calib(U_calibk,V_calibk,U_calibp,V_calibp,U_c,meanV,popt,n,pp,U_pw,V_pw,ucut,i):
    """
    This function plots the generated calibration points, and writes the expressions of the function used in the
    plots. It also shows the cut for the piecewise function. This is done for two calibration curves
    at the same time, which should describe the i-input. The return is a plot.
    """
    U_c = np.sort(U_c); meanV = np.sort(meanV)
    U_c[0] = 0.0
    plt.subplot(2,2,1+i)
    plt.plot(U_c,meanV,'og',U_calibk,V_calibk,'--k',U_calibp,V_calibp,'--r')
    plt.axvline(x=ucut)
    plt.legend(['Datapunkter','E^2=%.6s+%.6s U^{%s}' %(popt[0],popt[1],n),'E=(%.6s)+(%.6s) U+(%.6s) U^2+(%.6s) U^3+(%.6s) U^4' %(pp[0],pp[1],pp[2],pp[3],pp[4]),'Ux=%.5s'%ucut],fontsize = 7, loc = 4)
    if i == 0:
        plt.ylabel('E [V]',fontsize = 18)
    plt.grid('on')
    plt.subplot(2,2,3+i)
    plt.plot(U_c,meanV,'og',U_pw,V_pw,'-m')
    if i == 0:
        plt.ylabel('E [V]',fontsize = 18)
    plt.xlabel('U [m/s]', fontsize = 18)
    plt.legend(['Datapunkter','Epw'],fontsize = 20, loc = 4)
    plt.grid('on')
    return('Done')

def calib_writer(U_calib,V_calib,filename):
    """
    The function writes the points (U_calib,V_calib) which makes up the piecewise function to the folder
    \Calculated\Calibs. This is done in two columns. The filename have to be given and end with .dat.
    """
    os.chdir(r'C:\Users\jensj\OneDrive\Skrivebord\Beregninger\Calculated\Calibs')
    outfile = open(filename,'w')
    for uc,vc in zip(U_calib, V_calib):
        outfile.write('%24s %25s \n' % (uc,vc))
    outfile.close()
    return('Done')

def read_calc_write(infile,order,n,max,pw_lim,outfile):
    """
    The function run all of the of the functions above, except the visualize_calib-function due to error.
    Therefore, the visualize_calib-function is put outside this function.
    Give the name of the infile, the order of polynomial desired, n in Kings law, the piecewise limit
    and the name of the outfile.
    """
    meanV_cal,meanM_cal,Y_cal = read_mean(infile)
    Ub_cal,Re_cal = comp_UbRe(meanM_cal)
    Uc_cal = Nikuradse(Ub_cal,Re_cal)
    U_calibp, V_calibp,p,pp = calib_poly(Uc_cal,meanV_cal,order)
    U_calibk, V_calibk,popt,func1 = calib_kings(Uc_cal,meanV_cal,n,max)
    U_pw, V_pw = piecewise_calib(p,func1,max,pw_lim)
    calib_writer(U_pw,V_pw,outfile)
    return(U_calibk,V_calibk,U_calibp,V_calibp,Uc_cal,meanV_cal,popt,n,pp,U_pw,V_pw,pw_lim)

U_calibk,V_calibk,U_calibp,V_calibp,Uc_cal,meanV_cal,popt,n,pp,U_pw,V_pw,pw_lim = read_calc_write('calib_mean.dat',4,0.498,12.5,1.73,'calib1.dat')
visualize_calib(U_calibk,V_calibk,U_calibp,V_calibp,Uc_cal,meanV_cal,popt,n,pp,U_pw,V_pw,pw_lim,0)
U_calibk,V_calibk,U_calibp,V_calibp,Uc_cal,meanV_cal,popt,n,pp,U_pw,V_pw,pw_lim = read_calc_write('calib_mean3.dat',4,0.525,12.5,2.4,'calib3.dat')
visualize_calib(U_calibk,V_calibk,U_calibp,V_calibp,Uc_cal,meanV_cal,popt,n,pp,U_pw,V_pw,pw_lim,1)
plt.show()

U_calibk,V_calibk,U_calibp,V_calibp,Uc_cal,meanV_cal,popt,n,pp,U_pw,V_pw,pw_lim = read_calc_write('calib_mean4.dat',4,0.54,11.5,1.9,'calib4.dat')
visualize_calib(U_calibk,V_calibk,U_calibp,V_calibp,Uc_cal,meanV_cal,popt,n,pp,U_pw,V_pw,pw_lim,0)
U_calibk,V_calibk,U_calibp,V_calibp,Uc_cal,meanV_cal,popt,n,pp,U_pw,V_pw,pw_lim = read_calc_write('calib_mean5.dat',4,0.56,11.5,1.66,'calib5.dat')
visualize_calib(U_calibk,V_calibk,U_calibp,V_calibp,Uc_cal,meanV_cal,popt,n,pp,U_pw,V_pw,pw_lim,1)
plt.show()
