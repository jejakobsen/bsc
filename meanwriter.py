# Note: Before running the script, one has to #-out some of the meanwriter-functions due to memory
# loss if your computer isn't strong enough. Therefore the script have to be run in segments.
import os
import numpy as np

def conv_meas(directory,i):
    """
    Add directory path and which channel to read lines from.
    The directory should consist of n experiments.
    i = 0 = channel 1, i = 1 = channel 2. The function will go through every file
    in the directory which ends with ".txt", and read column i, 3 and the last (-1).
    The function returns 1 big list composed of 3 bigger sublists of respectively voltage, mean flow and time.
    Each sublist has n sublists, which represent the nth experiment.
    """
    V = []; M = []; T = [];
    os.chdir(directory)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        print(filename)
        if filename.endswith(".txt"):
            infile = open(filename,'r')
            infile.readline()
            infile.readline()
            txt = []; v = []; m = []; t = []
            for line in infile:
                txt = line.split()
                v.append(float(txt[i]))
                m.append(float(txt[3]))
                t.append(float(txt[-1]))
        else:
            print(filename)
            continue
        infile.close()
        V.append(v);M.append(m);T.append(t)
    list = [V,M,T]
    return(list)

def comp_mean(list):
    """
    Computes and returns the mean of mass flow and voltage of each sublist in the
    inputted list.
    """
    meanV = []; meanM = [];
    V = np.array(list[0])
    M = np.array(list[1])
    i = 0
    for v in V:
        meanV.append(np.mean(v))
        i = i + 1
    for m in M:
        meanM.append(np.mean(m))
    return(meanV,meanM)

def write_means(V,M,Y,filename):
    """
    The function takes in the mean voltage and massflow of n experiments.
    Y has to be put in manually and represents the distance from center at where
    each of the n experiment was done. Y has to be in the same order as each
    experiment was in the directory folder. All values are being written
    orderly in 3 columns in a file. Therefore the filename has to be an input.
    It is recommended to use .dat at the end of the filename. The files are
    conducted in a mutual folder called \Calculated\Mean.
    """
    os.chdir(r'C:\Users\jensj\OneDrive\Skrivebord\Beregninger\Calculated\Mean')
    outfile = open(filename,'w')
    for v, m, y in zip(V, M, Y):
        outfile.write('%24s %25s %5s \n' % (v,m,y))
    outfile.close()
    return(zip(V,M,Y))

def meanwriter(Y,path,channel,filename):
    """
    Add Y, path and which channel to read from in order to run all of the functions above.
    """
    dir = os.fsencode(path)
    data = conv_meas(dir,channel)
    meanV,meanM = comp_mean(data)
    write_means(meanV,meanM,Y,filename)
    return('Done')

#Y_calib = 13*[0.0]
#meanwriter(Y_calib,r'C:\Users\jensj\Desktop\Raw\Kalibrering1',1,'calib_mean.dat')
#meanwriter(Y_calib,r'C:\Users\jensj\Desktop\Raw\Kalibrering2',0,'calib_mean2.dat')
#meanwriter(Y_calib,r'C:\Users\jensj\Desktop\Raw\Kalibrering3',0,'calib_mean3.dat')
#Y_calib2 = 8*[0.0]
#meanwriter(Y_calib2,r'C:\Users\jensj\Desktop\Raw\Kalibrering4',1,'calib_mean4.dat')
#Y_calib3 = 9*[0.0]
#meanwriter(Y_calib3,r'C:\Users\jensj\Desktop\Raw\Kalibrering5',0,'calib_mean5.dat')
#Y_exp = [0.0,5.0,10.0,15.0,20.0,25.0,30.0,33.0,36.0,39.0,42.0,-5.0,-10.0,-15.0,-20.0,-25.0,-30.0,-33.0,-36.0,-39.0,-42.0,-44.0,-46.0,-47.0,-47.5]
#meanwriter(Y_exp,r'C:\Users\jensj\Desktop\Raw\Maaling',1,'exp_mean.dat')
#Y_exp2 = [0.0,5.0,10.0,15.0,20.0,25.0,30.0,33.0,36.0,39.0,42.0,45.0,47.0,-5.0,-10.0,-15.0,-20.0,-25.0,-30.0,-33.0,-36.0,-39.0,-42.0,-45.0,-47.0]
#meanwriter(Y_exp2,r'C:\Users\jensj\Desktop\Raw\Maaling2',0,'exp_mean2.dat')
Y_exp3 = [0.0,-1.0,-2.0,-10.0,-20.0,-30.0]
#meanwriter(Y_exp3,r'C:\Users\jensj\Desktop\Raw\Maaling3',1,'exp_mean3.dat')
meanwriter(Y_exp3,r'C:\Users\jensj\Desktop\Raw\Maaling4',0,'exp_mean4.dat')
