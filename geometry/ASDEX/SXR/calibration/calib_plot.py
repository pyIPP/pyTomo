#!/usr/bin/env python 
# -*- coding: utf-8 -*-

from numpy import *
import matplotlib 
#from mpl_toolkits.mplot3d import Axes3D
#matplotlib.rcParams['backend'] = 'Agg'
from scipy.signal import fftconvolve

from matplotlib.pyplot import * 
#from scipy.stats import ,nanstd, nanmean






#   první skok . změna všech diod, 2. taky, 3. změna filtrů (nyní by měli mít mnohem přesněji 75nm, 4. chyba H2 kamery 






SXR_detectors=('F','G','H1','H2','H3','I1','I2','I3','J1','J2','J3','K1','K2','L','M')


#for i in [0,2,5,8, 11, 12, 13]:
    #print SXR_detectors[i]
#exit()

if __name__ == "__main__":
    
    shots = []
    #print 
    calibs = []
    import os
    for file in os.listdir("./"):
        if file.endswith(".txt") and file[-5]!= '_':
            print file, file[:-4]
            if int(file[:-4]) > 26000:
                #print file[10:-4]
                try:
                    calibs.append(loadtxt(file))
                    shots.append( int(file[:-4]))

                except:
                    #continue
#if 

                    cams,calb_ = loadtxt(file, 
                        dtype={'names': ('cam', 'calib'),'formats': ('S4', 'd')},
                        unpack=True)

                    #print len(calb_)
                    shots.append( int(file[:-4]))

                    calibs.append(calb_)
                if len(calibs[-1]) == 14: #missing F camera
                    calibs[-1] = r_[nan, calibs[-1]]
                #print len(calibs[-1])

                #SaveAsher(  int(file[10:-4]))
    ind = argsort(shots)
    shots = array(shots)[ind]
    
    calibs = array(calibs)[ind,:]
    
    print 'shots:', len(shots)
    
    ind = ones_like(shots, dtype='bool')
    
    for s in [27790,28750,27791 ]: 
        ind[shots == s] = False
    
    shots = shots[ind]
    calibs = calibs[ind,:]

    ind = (shots > 30420)&(shots < 30447)
    calibs[ind,2] = nan
    calibs[ind,4] = nan
    calibs[(shots>30420)&(shots<31900),3] = nan

    #calibs[:,5] = nan
    #calibs[:,7] = nan
    #calibs[shots==30898,2] = nan
    #calibs[shots==29756,6] = nan
    #calibs[shots==30904,2] = nan
    #calibs[(shots>=25995)&(shots<=27392) ,2] = nan
    #calibs[shots==27752,4] = nan
    #calibs[shots==26137,4] = nan
    #calibs[shots==28581,(-3,-4)] = nan
    #calibs[shots==28318,(-1,-2)] = nan
    ##calibs[shots==28313,(-1,-2)] = nan
    #calibs[shots==28317,(-1,-2)] = nan
    #calibs[shots==28319,(-1,-2)] = nan
    #calibs[shots==28320,(-1,-2)] = nan
    #calibs[shots==29795,(-2,)] = nan
    #calibs[shots==31563,(-2,)] = nan
    #calibs[shots==31585,(-2,)] = nan
    #calibs[shots==31426,(-1,)] = nan
    #calibs[shots==31318,(-1,)] = nan

    
    

    
    calib_mean = nanmean(calibs[:,[0,3,6,9,12,13,14]],axis=1)
    calib_mean[(shots>30420)&(shots<31900)] = nanmedian(calibs[(shots>30420)&(shots<31900),:][:,[1,6,9,12,13,14]],axis=1)
    
    #plot(shots,calib_mean ,'o' )
    #show()
    
    calibs/= calib_mean[:,None]

    nplot = len(SXR_detectors)-2
    nrow = int(sqrt(nplot))
    ncol = int(ceil(float(nplot)/nrow))
    
    f, axis = subplots( nrow,ncol, sharex=True, sharey=True)
    
    f.subplots_adjust(hspace=0.05, wspace = 0.12)
    from matplotlib.ticker import MaxNLocator
    
    
    savetxt('all_calibrations', c_[shots,calibs ], header='shot '+'\t'.join(SXR_detectors), fmt='%.d '+('\t%.4f'*calibs.shape[1]))

    i_plot = 0
    axis_flat = axis.ravel()
    for  i in range(calibs.shape[1]):

        ind1 = (shots>26000)&(shots<=27350)
        ind2 = (shots>27350)&(shots<=28523)
        ind3 = (shots>28523)&(shots<=30135)
        ind4 = (shots>30135)&(shots<=30420)
        ind5 = (shots>30420)&(shots<=31700)
        ind6 = (shots>31700)

        print i, ',%.3f,%.3f,%.3f,%.3f,%.3f,%.3f'%(nanmedian(calibs[ind1,i]), nanmedian(calibs[ind2,i]), nanmedian(calibs[ind3,i]), nanmedian(calibs[ind4,i]), nanmedian(calibs[ind5,i]), nanmedian(calibs[ind6,i]))

        if SXR_detectors[i] in ('I1','I3'):
            continue
        
        #print shots, calibs[:,i]
        #axis_flat[i_plot].set_title(SXR_detectors[i])
        axis_flat[i_plot].plot(shots, calibs[:,i],',')
        #show()
        axis_flat[i_plot].set_ylim(0.8, 1.2)
        axis_flat[i_plot].axhline(y=1,c='k',ls='--')
        axis_flat[i_plot].text(0.1, 0.1, SXR_detectors[i] , transform=axis_flat[i_plot].transAxes,backgroundcolor='w') 
        
        axis_flat[i_plot].set_xticks([26000,28000,30000,32000,34000])
        axis_flat[i_plot].set_yticks([.8,.9,1,1.1,1.2])

        #if SXR_detectors[i][0] == 'H':
        axis_flat[i_plot].axvline(x = 30420,c='k')
     
        axis_flat[i_plot].axvline(x = 27350,c='k')
        axis_flat[i_plot].axvline(x = 28523,c='k')
        axis_flat[i_plot].axvline(x = 30135,c='k')
        axis_flat[i_plot].axvline(x = 31779,c='k')
        axis_flat[i_plot].axvline(x = 33800,c='k')

        axis_flat[i_plot].grid('on')

        i_plot += 1
        #axis_flat[i_plot].yaxis.set_major_locator(MaxNLocator(6))


    print shots
    axis_flat[0].set_xlim(shots.min(), shots.max())
    axis_flat[0].text(27350,0.9,27350,backgroundcolor='w',rotation=90,fontsize=10)
    axis_flat[0].text(28523,0.9,28523,backgroundcolor='w',rotation=90,fontsize=10)
    axis_flat[0].text(30135,0.9,30135,backgroundcolor='w',rotation=90,fontsize=10)
    axis_flat[0].text(30620,0.9,30620,backgroundcolor='w',rotation=90,fontsize=10)
    axis_flat[0].text(31779,0.9,31779,backgroundcolor='w',rotation=90,fontsize=10)
    axis_flat[0].text(33800,0.9,33800,backgroundcolor='w',rotation=90,fontsize=10)

    #savefig('SXR_calibration.pdf')
    show()
        
        
        
        
        
    