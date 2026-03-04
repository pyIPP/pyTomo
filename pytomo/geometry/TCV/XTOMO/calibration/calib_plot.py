#!/usr/bin/env python 
# -*- coding: utf-8 -*-

from numpy import *
import matplotlib 
#from mpl_toolkits.mplot3d import Axes3D
#matplotlib.rcParams['backend'] = 'Agg'
from scipy.signal import fftconvolve

from matplotlib.pyplot import * 
from scipy.stats import nanmedian,nanstd, nanmean






#   první skok . změna všech diod, 2. taky, 3. změna filtrů (nyní by měli mít mnohem přesněji 75nm, 4. chyba H2 kamery 






SXR_detectors=['A','B','C','D','E','F','G','H','I','J']


#for i in [0,2,5,8, 11, 12, 13]:
    #print SXR_detectors[i]
#exit()

if __name__ == "__main__":
    
    shots = []
    print 
    calibs = []
    import os
    for file in os.listdir("./"):
        if file.endswith(".txt"):
            #if int(file[10:-4]) > 26000:
                #print file[10:-4]
            
                #continue
            cams,calb_ = loadtxt(file, 
                dtype={'names': ('cam', 'calib'),'formats': ('S4', 'd')},
                unpack=True)
            #print len(calb_)
            try:
                if all(calb_==1):
                    continue
                shots.append( int(file[12:-4]))

                calibs.append(calb_)
            except:
                pass
            #if len(calibs[-1]) == 14: #missing F camera
                #calibs[-1] = r_[nan, calibs[-1]]
            #print len(calibs[-1])

            #SaveAsher(  int(file[10:-4]))
    ind = argsort(shots)
    shots = array(shots)[ind]
    
    calibs = array(calibs)[ind,:]
    
    print 'shots:', len(shots)
    
    ind = ones_like(shots, dtype='bool')
    
    #for s in [27790,28750,27791 ]: 
        #ind[shots == s] = False
    
    shots = shots[ind]
    calibs = calibs[ind,:]

    
    
    

    
    calib_mean = nanmean(calibs,axis=1)
    #calib_mean[(shots>30420)&(shots<31900)] = nanmedian(calibs[(shots>30420)&(shots<31900),:][:,[1,6,9,12,13,14]],axis=1)
    calibs/= calib_mean[:,None]

    nplot = len(SXR_detectors)
    nrow = int(sqrt(nplot))
    ncol = int(ceil(float(nplot)/nrow))
    
    f, axis = subplots( nrow,ncol, sharex=True, sharey=True)
    
    f.subplots_adjust(hspace=0.05, wspace = 0.12)

    i_plot = 0
    axis_flat = axis.ravel()
    for  i in range(calibs.shape[1]):


        #print shots, calibs[:,i]
        #axis_flat[i_plot].set_title(SXR_detectors[i])
        axis_flat[i_plot].plot(shots, calibs[:,i],'.')
        #show()
        axis_flat[i_plot].set_ylim(0.8, 1.2)
        axis_flat[i_plot].axhline(y=1,c='k',ls='--')
        axis_flat[i_plot].text(0.1, 0.1, SXR_detectors[i] , transform=axis_flat[i_plot].
                transAxes,backgroundcolor='w')
        
        #axis_flat[i_plot].set_xticks([26000,28000,30000,32000])
        #if SXR_detectors[i][0] == 'H':
        #axis_flat[i_plot].axvline(x = 30420,c='k')
     
        #axis_flat[i_plot].axvline(x = 27350,c='k')
        #axis_flat[i_plot].axvline(x = 28523,c='k')
        #axis_flat[i_plot].axvline(x = 30135,c='k')
        #axis_flat[i_plot].axvline(x = 31700,c='k')

        axis_flat[i_plot].grid('on')

        i_plot += 1
        

    print shots
    axis_flat[0].set_xlim(50275, 50325)
    axis_flat[0].set_ylim(0,2)
    
    for s,f in   zip(SXR_detectors,calibs[shots > 50200].mean(0)):
        print '%.3f,'%f,
        
    print

    #axis_flat[0].text(27350,0.9,27350,backgroundcolor='w',rotation=90,fontsize=10)
    #axis_flat[0].text(28523,0.9,28523,backgroundcolor='w',rotation=90,fontsize=10)
    #axis_flat[0].text(30135,0.9,30135,backgroundcolor='w',rotation=90,fontsize=10)
    #axis_flat[0].text(30620,0.9,30620,backgroundcolor='w',rotation=90,fontsize=10)
    #axis_flat[0].text(31700,0.9,31700,backgroundcolor='w',rotation=90,fontsize=10)

    savefig('SXR_calibration.pdf')
    show()
        
        
        
        
        
    