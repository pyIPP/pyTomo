#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pmds
from numpy import *
import matplotlib
from  matplotlib.pyplot import *
matplotlib.rcParams['backend'] = 'Qt4Agg' 

def load_mag_surfaces(shot=0,time=0,psisur=0, number_of_pol_points=0):
    shot   =       61372
    time_vec    =       [60.5]   # Default time
    number_of_pol_points = 60   # Number of grid points
    nsurf = 30.0
    psisur  = (arange(nsurf)/nsurf)**1.5  # from Ernesto (good spacing of flux surfaces)
    psisur[0] = 0.02 #avoid singularity
    
    rsurf   =   zeros((number_of_pol_points ,len(psisur), len(time_vec)))
    zsurf   =   zeros(shape(rsurf))
    pmds.mdsconnect('mdsplus.jet.efda.org')
    
    for j in range(len(time_vec)):
        for i in range(len(psisur)):
            print(i,j, psisur[i])
            try:
                pmds.mdsvalue('flushinit(15,' +str(shot)+ ', ' +str(time_vec[j])+ ' ,0,"JETPPF","EFIT",_err)')
            except:
                pass
            pmds.mdsvalue('flupx(_err)')  # run flush
            print('flusur( ' +str(psisur[i])+ ',' +str(number_of_pol_points) +',_xp,_yp,_br,_bz,1) ')
            pmds.mdsvalue('flusur( ' +str(psisur[i])+ ',' +str(number_of_pol_points) +',_xp,_yp,_br,_bz,1) ')  # retrieve data

            rsurf[:,i,j] =  pmds.mdsvalue('_xp')
            zsurf[:,i,j] =   pmds.mdsvalue('_yp')

    pmds.mdsdisconnect()

    #plot(squeeze(rsurf), squeeze(zsurf))
    #show()


