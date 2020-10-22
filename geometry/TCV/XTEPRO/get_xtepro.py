#!/usr/bin/env python
# -*- coding: utf-8 -*-


#ssh todstrci@lac911.epfl.ch -L 8002:tcv1.epfl.ch:8000

from numpy import *
#import pmds
#from matplotlib.pylab import *
#import ././tqdm
 
import sys,os
from scipy.io import loadmat
#print os.path.abspath('../')
#exit()
parent_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

sys.path.append(os.path.abspath(parent_path))





def dget_xtepro(connect, shot,tmin=-infty,tmax=infty):
    

    path = os.path.dirname(os.path.realpath(__file__))
    
    geom = loadmat(path+'/xtepro_geometry.mat')  # load geometry parameters
    aomega = geom['aomega']


    # calculation of time base
    connect.openTree('tcv_shot',shot);   


    #**************************** get data **********************************
        
    #xtepro = [
    #(r'\atlas::dt100_southwest_001', 'channel_%.3d', range(1,33)),
    #(r'\atlas::dt100_southwest_002', 'channel_%.3d', range(1,17))]
    

    from mds_data import mds_load
    #data = []
    nodes = [r'\atlas::dt100_southwest_001', r'\atlas::dt100_southwest_002']
    chns = list(range(1,33)), list(range(1,17))
    chform = 'channel_%.3d'
    #for node, chform, chn in xtepro:
    tvec, data = mds_load(connect, nodes,chns, chform,remove_offset=False,
                            tmin=tmin,tmax=tmax)
    
        #data.append(sig)
        
    #data = hstack(data)
    # calculation of the power detected by the diode in W/m2

    #gains=1.5e6; # this is the gain actually used on this card 
    gains=1e6;
    curr_pow=0.24; # diode responsivity in in A/W
                # the theoretical value is 0.2754A/W, but the value
                # above is used in BOIVIN...   
    corr=aomega.flatten('F')/(4*pi); # correction of 4pi for the AOMEGA factor
                                                                                                                                    

    saturated = data > 10
    data-= data[tvec<0].mean(0)
    data /= corr*gains*curr_pow; # result expressed in W/m2 
    

    return tvec,data,saturated,geom


