#!/usr/bin/env python
# -*- coding: utf-8 -*-

#ssh todstrci@lac911.epfl.ch -L 8002:tcv1.epfl.ch:8000

#import pmds
from matplotlib.pylab import *
#import ././tqdm
import os.path 
from numpy import *

import sys,os
#print os.path.abspath('../../')
#sys.path.append(os.path.abspath('../../'))
#import tqdm
from scipy.io import loadmat
#BUG axuv_detectorpos
#axuv_calc_los
#axuv_get_calibration

parent_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

sys.path.append(os.path.abspath(parent_path))





def get_fir(connect, shot):
    
    #pmds.mdsopen('tcv_shot',shot)
    #fir_data = r'\results::fir_lin_int_dens_array'
    #tvec = pmds.mdsvalue('_t=dim_of('+fir_data+') ')
    #fir_data = pmds.mdsvalue('_x=data('+fir_data+') ')
    #R = pmds.mdsvalue('\diagz::fir_array:radii')
    #Z = -0.75, 0.75
    #pmds.mdsclose('tcv_shot',shot)
    
    
    connect.openTree('tcv_shot',shot)
    fir_data = r'\results::fir_lin_int_dens_array'
    tvec = asarray(connect.get('_t=dim_of('+fir_data+') '))
    fir_data = asarray(connect.get('_x=data('+fir_data+') '))
    R = asarray(connect.get('\diagz::fir_array:radii'))
    Z = -0.75, 0.75
    connect.closeTree('tcv_shot',shot)
    
    if shot > 45179:
        fringe = 12.078e18
    else:
        fringe = 10.367e18


    #FRIGES CORRECTION
    #NOTE!!  only experimental, not extensively tested
    for k in arange(len(fir_data)):
        sig = cumsum(fir_data[k])
        from scipy.signal import argrelmax

        #remove phase jumps and some low frequency information
        knots,vals = [],[]
        ind_knots = r_[list(range(0,len(tvec),20)),-1]
        peaks = infty
        ind = slice(None,None)
        dsig,tvec_ = sig,tvec
        for i in range(7):
            if not any(peaks>11):
                break
            knots.extend(tvec_[ind_knots])
            vals.extend(dsig[ind_knots])
            sort_ind = argsort(knots)
            retrofit = interp(tvec_, array(knots)[sort_ind], array(vals)[sort_ind])
            peaks = abs(dsig-retrofit)
            ind = peaks> fringe*1.5
            ind_knots = argrelmax(peaks[ind],0,100)[0]
            dsig = dsig[ind]
            tvec_ = tvec_[ind]

        sort_ind = argsort(knots)
        retrofit = interp(tvec, array(knots)[sort_ind], array(vals)[sort_ind])
        #plot(diff(retrofit));show()
        
        fringes = diff(retrofit)-unwrap(diff(retrofit)/1e19*2*pi, pi*1.5)*1e19/(2*pi)
        #plot(copy(fir_data[k]))
        #plot(diff(retrofit)-unwrap(diff(retrofit)/1e19*2*pi, pi*1.5)*1e19/(2*pi) )
        #plot(diff(retrofit) )
        
        #show()
        
        fir_data[k]-= r_[fringes,0]
        #plot(fir_data[k])
        #show()

    
    #remove vibration at 46-75 Hz!!
    background = fir_data.mean(0)
    n = len(background)/2*2
    fback = fft.rfft(background[:n])
    fvec = linspace(0,.5/(tvec[1]-tvec[0]), n/2+1)
    fback[slice(*fvec.searchsorted([46,75]))] = 0

    background_ = fft.irfft(fback)

    fir_data = fir_data[:,:n]-(background[:n]-background_)
    tvec = tvec[:n]
        
    
    
    
    #import IPython
    #IPython.embed()
    return tvec,fir_data,R,Z



