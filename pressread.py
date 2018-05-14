import numpy as np
import os

def conv_meas(directory):
    """
    The function needs the directory path and will return every dP-value as
    one array.
    """
    dP = []
    os.chdir(directory)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        print(filename)
        if filename.endswith('.txt'):
            infile = open(filename,'r')
            infile.readline()
            infile.readline()
            for line in infile:
                txt = line.split()
                dp = txt[-2]
                dp = dp.replace(",",".")
                dP.append(float(dp))
        else:
            print(filename)
            continue
        infile.close()
    return(dP)

def utau(dp):
    """
    The function takes in the array of dP, computes the mean of dP
    and calculates and returns tau and utau.
    """
    from math import sqrt
    dp = np.array(dp); R = 0.05; rho = 1.19
    print(np.mean(dp)) # print mean dP
    x = 12.281
    tau = (np.mean(dp)/x)*R/2;
    utau = sqrt(tau/rho)
    return(tau,utau)

press_dir = os.fsencode(r'C:\Users\jensj\Desktop\Raw\Trykkfall')
dP = conv_meas(press_dir)
tau, utau = utau(dP)
print(tau,utau) # print tau and utau
