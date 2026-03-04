#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import *
#import pmds
#from matplotlib.pylab import *
import os.path
import sys,os
parent_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.abspath(parent_path))
sys.path.append(os.path.abspath(parent_path)+'/AXUV')
from scipy.io import loadmat
from get_axuv import  axuv_calc_los,geometric_data_class,axuv_detectorpos

"""
python version of axuv_get_bolo.m

rewriten by Tomas Odstrcil   tom@cbox.cz




"""



class bolostruct_class:
    def __init__(self):
        class empty:  pass
        self.geometry = empty()
        self.raw = empty()
        self.special = empty()




def get_bolo(connect,shotno,**kwarg):
    #[bolostruct]=et_bolo(shotno,varargin)
    #
    #e.g.: bolostruct=get_bolo(44697,'filtertype','gottardi')
    #
    #       varargin:               default
    #
    #       refit_tau               1
    #       offset_removal          1
    #       filtertype              'gottardi', 'bessel', 'golay', 'interpos', 'butter'
    #
    #Bessel filter parameters:
    #       filtercutoff            30
    #Gottardi parameters (see axuv_gottardi_simplified.m):
    #       knoise                  20
    #       ibave                   80
    #       ieave                   10
    #       alevel                  0.16
    #       time_noise              [-0.04 0]

    print ('Bolo data download has started.')
    
    refit_tau = kwarg['refit_tau'] if 'refit_tau' in kwarg else False
    offset_removal = kwarg['offset_removal'] if 'offset_removal' in kwarg else True
    filtercutoff = kwarg['filtercutoff'] if 'filtercutoff' in kwarg else 30
    filtertype = kwarg['filtertype'] if 'filtertype' in kwarg else 'bessel'
    knoise = kwarg['knoise'] if 'knoise' in kwarg else 20
    ibave = kwarg['ibave'] if 'ibave' in kwarg else  80
    ieave = kwarg['ieave'] if 'ieave' in kwarg else  10
    alevel = kwarg['alevel'] if 'alevel' in kwarg else  0.16
    time_noise = kwarg['time_noise'] if 'time_noise' in kwarg else [-.04, 0]

    bolostruct = bolostruct_class()

    if refit_tau:
        bolostruct.raw.tau=axuv_calc_tau_for_bolo(44703,1.16,1.16+0.1)
    

    #---properties---
    bolostruct.name='BOLO'
    bolostruct.shotno=shotno
    #---properties---
     
    #pmds.mdsopen('tcv_shot',shotno)
    connect.openTree('tcv_shot',shotno)

    #---raw data---
    #bolostruct.raw.data=pmds.mdsvalue(r'\base::bolo:source').T
    bolostruct.raw.data=asarray(connect.get(r'\base::bolo:source')).T

        
    #import IPython
    #IPython.embed()
    
    #bolostruct.raw.gains=pmds.mdsvalue(r'\base::bolo:gains')
    #bolostruct.raw.calibration=pmds.mdsvalue(r'\base::bolo:calibration')
    #etend=pmds.mdsvalue(r'\base::bolo:geom_fact')
    #if not hasattr(bolostruct.raw,'tau'):
        #bolostruct.raw.tau=pmds.mdsvalue(r'\base::bolo:tau')
            
    bolostruct.raw.gains=asarray( connect.get(r'\base::bolo:gains'))
    bolostruct.raw.calibration = asarray(connect.get(r'\base::bolo:calibration'))
    etend=asarray(connect.get(r'\base::bolo:geom_fact'))
    if not hasattr(bolostruct.raw,'tau'):
        bolostruct.raw.tau=asarray(connect.get(r'\base::bolo:tau'))
            
                 
            
    diag_path = os.path.dirname(os.path.realpath(__file__))+'/'
    bolostruct.raw.convfact=9.4*0.0015*0.004
        #this is the excitation voltage of the bridges
        #1/(excitation_voltage*detector_area)=1/(9.4*0.0015*0.004)=1.78e4
    #---raw data---

    #---geometry---
    bolostruct.geometry=axuv_calc_los(axuv_detectorpos('BOLO'))
    #---geometry---
    if os.path.isfile(diag_path+'axuv_finiteetend_bolo.mat'):
        geom = loadmat(diag_path+'axuv_finiteetend_bolo.mat')
        bolostruct.geometry.etend=squeeze(geom['etendstruct']['etend'].item())
    
    bolostruct.geometry.etend=etend
    #bolostruct.time=pmds.mdsvalue(r'dim_of(\base::bolo:signals)')
    bolostruct.time=asarray(connect.get(r'dim_of(\base::bolo:signals)'))
    start=float(connect.get('timing("401")'))
    #start=pmds.mdsvalue('timing("401")')        #this is the point when the offset removing mechanism switched off
                                            #it is worth waiting a small time after it
    if start < 0:
        offsetindexes = (bolostruct.time<=0)  &  (bolostruct.time>=0.7*start)
    else:
        raise Exception('Bad adjustment!')
    
    #plot(bolostruct.raw.data.T)
    #show()     

    bolostruct.raw.offset = mean(bolostruct.raw.data[offsetindexes],0) #these points are acquired before t=0
    if offset_removal:
            bolostruct.special.static_voltage = bolostruct.raw.data - bolostruct.raw.offset
    else:
            bolostruct.special.static_voltage = bolostruct.raw.data
    

    calib = bolostruct.raw.gains*bolostruct.raw.calibration*\
                bolostruct.geometry.etend*bolostruct.raw.convfact
            
    bolostruct.special.static_los = bolostruct.special.static_voltage/calib

    bolostruct.data=zeros_like(bolostruct.raw.data,dtype='single') 
    bolostruct.special.dynamic_voltage=zeros_like(bolostruct.data)
    #plot(bolostruct.raw.data)
    #show()


    f_nq=1/(bolostruct.time[1]-bolostruct.time[0])/2

    if  filtertype == 'bessel':
        #---filter preparation---
        #filtercutoff=30 #Hz This value is suggested by the literature.
        from scipy.signal import bessel, bilinear,lfilter
   
        banalderivfordigit,aanalderivfordigit=bessel(4,filtercutoff/f_nq,analog=False)
        banalderivfordigit=r_[banalderivfordigit[1:],0]
        bderivdigit,aderivdigit=bilinear(banalderivfordigit,aanalderivfordigit, f_nq*2)
        #bderivdigit,aderivdigit = banalderivfordigit,aanalderivfordigit
        #---filter preparation---

        for m in range(size(bolostruct.raw.data,1)):
            #BUG instability in the filter!!!!
            dynamic_part=lfilter(bderivdigit,aderivdigit,bolostruct.special.static_voltage[:,m])*bolostruct.raw.tau[m]
            bolostruct.special.dynamic_voltage[:,m]=bolostruct.special.static_voltage[:,m]+dynamic_part
            bolostruct.data[:,m]=bolostruct.special.static_los[:,m]
            bolostruct.data[:,m]+=dynamic_part/(bolostruct.raw.gains[m]*bolostruct.raw.calibration[m]\
                                                *bolostruct.geometry.etend[m]*bolostruct.raw.convfact)
            
            
        
    elif filtertype == 'gottardi':
        #BUG neot finished axuv_gottardi_simplified
        for m in range(size(bolostruct.raw.data,1)):
            ts,ys,ds=axuv_gottardi_simplified(bolostruct.time,bolostruct.special.static_voltage[:,m],
                'knoise',knoise,'ibave',ibave,'ieave',ieave,'alevel',alevel,'time_noise',time_noise)
            bolostruct.special.dynamic_voltage[:,m]=ys+bolostruct.raw.tau[m]*ds
            bolostruct.data[:,m] = (ys+bolostruct.raw.tau[m]*ds)\
                    /(bolostruct.raw.gains[m]*bolostruct.raw.calibration[m]
                    *bolostruct.geometry.etend[m]*bolostruct.raw.convfact)
        
        bolostruct.time=ts
        
    elif filtertype == 'interpos':
        #not finished interpos
        for m in range(size(bolostruct.raw.data,1)):
            yout,youtp = interpos(13,bolostruct.time,bolostruct.raw.data[:,m],bolostruct.time,1e-8)
            bolostruct.data[:,m] = (yout+bolostruct.raw.tau[m]*youtp)\
                /(bolostruct.raw.gains[m]*bolostruct.raw.calibration[m]
                   *bolostruct.geometry.etend[m]*bolostruct.raw.convfact)
                
        
    elif filtertype == 'golay': #Savitzky-Golay filter   

        from scipy.signal import savgol_coeffs
        N = 4                 # Order of polynomial fit
        F = 21                # Window length
          # Calculate S-G coefficients
        g = [savgol_coeffs(F,N, i) for i in range(3)]
        
        dx = mean(diff(bolostruct.time))
        x = bolostruct.time
        
        #bolostruct.data = zeros(size(bolostruct.raw.data), dtype=float)
        for m in range(size(bolostruct.raw.data,1)):            
            y = bolostruct.raw.data[:,m]          
            HalfWin  = ((F+1)/2) -1            
            SG0 = convolve(g[0][::-1],y ,mode='same')
            SG1 = convolve(g[1][::-1],y ,mode='same')/dx # Turn differential into derivative
            #SG2 = convolve(g[2][::-1],y ,mode='same')/dx**2  # and into 2nd derivative
            calib=bolostruct.raw.gains[m]*bolostruct.raw.calibration[m]*\
                bolostruct.geometry.etend[m]*bolostruct.raw.convfact
            bolostruct.data[:,m] = (SG0-SG1*bolostruct.raw.tau[m])/calib
        
    elif filtertype == 'butter':
        from scipy.signal import butter, filtfilt

        facq = 1./mean(diff(bolostruct.time))
        b,a = butter(5,200./facq) # Butterworth filter at 100 Hz
        
        for m in range(size(bolostruct.raw.data,1)):
            rp=gradient(bolostruct.raw.data[:,m])
            rp=filtfilt(b,a,rp)
            calib=bolostruct.raw.gains[m]*bolostruct.raw.calibration[m]*\
                bolostruct.geometry.etend[m]*bolostruct.raw.convfact
            bolostruct.data[:,m]=(bolostruct.raw.data[:,m]+rp*bolostruct.raw.tau[m])/calib

    else:
        raise Exception('Bad evaluation method')
    

    #pmds.mdsclose('tcv_shot',shotno)

    connect.closeTree('tcv_shot',shotno)

    print ('Bolo data download has finished.')

    return bolostruct









#liší se raw.tau


