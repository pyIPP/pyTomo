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





def dget_xtomo(connect, shot,tmin=-infty,tmax=infty):
    
    from mds_data import mds_load

        #
    #                                                                  -
    #[sig,t,Sat_channel,cat_default] = XTOMOGetData(shot,tmin,tmax,fans,tag)
    #
    #       INPUT:
    #       shot:   TCV shot
    #       tmin:   start time (optional, def: 0.5)
    #       tmax:   stop time (optional, def: 0.6)
    #       fans:   camera switch (optional, def. depends on shot number)
    #       tag:    'full' or 'brief' or 'liuqe' or 'thomson'
    #               or 1 (50kHz) or 2 (25kHz) (optional, def. depends on shot number)
    # OUTPUTS:
    #       sig:    calibrated xtomo signals, size: [20xsum(fans) x length(t)]
    #       tvec:   time vector
    #       Sat_channel: Saturated channels. Must be removed prior inversion
    #       cat_default: A structure with some default parameters (depends on shot number)
    # CALLS:
    #    [sig,t,Sat_channel,cat_default] = XTOMOGetData(18128,0.6,0.63,[1 1 1 1 1 1 1 1 0 1],'brief')
    #    [sig,t,Sat_channel,cat_default] = XTOMOGetData(18128,[],[],[],[])
    #    [sig,t,Sat_channel,cat_default] = XTOMOGetData(45339,0.6,0.63,[1 1 1 1 1 1 1 1 0 1],[])
    #    [sig,t,Sat_channel,cat_default] = XTOMOGetData(45339,0.6,0.63,[],'liuqe')

    # Created on 2012/02/21. CRPP/EPFL. Benoit Labit
    # Modified on 2012/02/27: If isnan(tmin) returns only cat_default
    # Modified on 2012/03/30: Ends properly even if some data are missing (shot>34800)
    # Modified on 2013/02/07: Add the possibility to download at 50kHz or 25kHz

    ##################################################################
    #
    #                        INITIALIZATION
    #
    ##################################################################



    satlevel=10  # Level in volts for saturated signals
    minlevel=3e-4

    if shot >= 13836 and shot <= 13848:
        offset_start=-0.02
        offset_stop=-0.01           #  this  one is to be used only for
    else:                                           #  shot=13836 to shot=13848
        offset_start=-0.04
        offset_stop=-0.01
    
    G,A,C = sxr_gains(connect, shot)
    #import IPython
    #IPython.embed()
    path = os.path.dirname(os.path.realpath(__file__))

    if shot<34800:
        geom=loadmat(path+'/cat_defaults2001.mat')
    else:
        geom=loadmat(path+'/cat_defaults2008_2.mat')

        #raise Exception('old geometry!! not implemented')
    


    
    print('  acquisition with 3 96-channels Dtacq cards')
    
    ##  Check that all data are there
    #pmds.mdsopen('tcv_shot',shot)
    connect.openTree('tcv_shot',shot)
    
    # Correspondance between diodes and acquisition channels due to incorrect cabling
    if shot <= 34800:
        INDEX=r_[0:180,180+array([2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,20,19])]
    else:     
        INDEX=r_[0:20, 39:19:-1,59:39:-1,60:80, 80:100,119:99:-1,139:119:-1,140:160,179:159:-1,199:179:-1]
        
    if len(INDEX)!= len(unique(INDEX)):
        raise Exception('wronmg index')


    #xtomo = [
        #(r'\atlas::dt196_xtomo_001', 'channel_%.3d', range(1,97)),
        #(r'\atlas::dt196_xtomo_002', 'channel_%.3d', range(1,97)),
        #(r'\atlas::dt196_xtomo_003', 'channel_%.3d', range(1,9))]
    
    nodes = [r'\atlas::dt196_xtomo_%.3d'%i for i in range(1,4)] 
    chns = list(range(1,97)), list(range(1,97)),list(range(1,9))
    chform = 'channel_%.3d'
    #all_data = []
    #for node, chform, chn in xtomo:
        #print 'loading from ', node
    tvec, all_data = mds_load(connect, nodes,chns, chform,remove_offset=False,
                              tmin=tmin,tmax=tmax)
        #print data.shape, tvec.shaep, tmin, tmax
        #all_data.append(data.T)
    #print 'loaded'
                            
    #import IPython
    #IPython.embed()
    all_data = all_data[:,INDEX]
    #data_list = [all_data[i/96][i%96] for i in INDEX]
    #all_data = vstack(data_list).T
    

    # Check for saturation
    saturated = all_data == satlevel# |(all_data <- satlevel)

    #print 'odecist offset!!!'
    

    ind = slice(tvec.searchsorted(offset_start), tvec.searchsorted(offset_stop))
    if ind.start == ind.stop:  #late start of measuremenst
        ind = slice(0,300)
    #print all_data.shape
    #exit()
    all_data-= all_data[ind].mean(0)

                            
    #import IPython
    #IPython.embed()
    #print 'Jak zkontrolovat že je to dobře?? nějaký hodně protažený výboj? '
    eta = TOMOCal(shot)
    calib = ones(all_data.shape[1])
    
    calib *= A   #1/etendue cm^2 * sterad'
    calib /= G   #gain
    calib *= 1/eta*4*pi*3.63/2.2e5   #absolute calibration, W/mm^2 sterad --> W/m^3  
    calib *= 1/C # correction to larger edge cells surface in detector array 

    all_data *= calib  #apply corrections
    
    #import IPython
    #IPython.embed()

    print('data from ATLAS::DT_196 loaded')
    

    return tvec,all_data,saturated,geom











def sxr_gains(connect, shot):
    # [G,A]=sxr_gains(shot,camera,diodes)
    # loads a SXR camera with a single line command
    # Works for any shot 
    #
    # INPUT:1) shot = shot number
    # OUTPUT:       G = gains specified cameras and diodes (1,10,100,1000)
    #
    #
    # Created: March 2009
    # Modified: June 2009
    #               I realized that G was incorrect because of the wrong
    #               cabling! But A is correct.
    #
    # Benoit Labit - CRPP/EPFL
    
 



    #cat_default=loadmat('/home/labit/matlab/XTOMO/cat_defaults2008_2.mat')

    
    path = os.path.dirname(os.path.realpath(__file__))
    if shot<=34800:
        cat_default = loadmat(path+'/cat_defaults2001.mat')
    else:
        cat_default = loadmat(path+'/cat_defaults2008_2.mat')
    
                

        
    angfact = cat_default['angfact']
    corfact = cat_default['corfact']



    n_cam = 10
    n_diod = 20

    # Correspondance between diodes and acquisition channels due to incorrect cabling
    if shot<=34800:
        INDEX=r_[0:180,180+array([2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,20,19])]
    else:     
        INDEX=r_[0:20, 39:19:-1,59:39:-1,60:80, 80:100,119:99:-1,139:119:-1,140:160,179:159:-1,199:179:-1]

    #pmds.mdsopen('tcv_shot',shot)
    connect.openTree('tcv_shot',shot)

    #gain
    G = empty((n_diod, n_cam))
    for ic in arange(n_cam):
        for id in arange(n_diod):
            com=r'\vsystem::tcv_publicdb_i["XTOMO_AMP:%.3d_%.3d"]'%(ic+1,id+1)
            #G[id,ic]=10**pmds.mdsvalue(com)
     
            G[id,ic]=10**int(connect.get(com))

            #A[id,ic]=angfact[d,c]
            
        

    #pmds.mdsclose('tcv_shot',shot)
    connect.closeTree('tcv_shot',shot)


    G = G.flatten(order='F')[argsort(INDEX)]
    A = angfact.flatten(order='F')#[argsort(INDEX)]
    C = corfact.flatten(order='F') 
 
    return G, A, C





def TOMOCal(shot):

    # time optimized version of MA's  xtomo_calibrate.m (located in
    # /NoTivoli/labit/IFmatlab/xtomo/)
    #
    # output:       calibration factors eta

    # G Turri Jan 2008
    # Modified by Benoit Labit

    kT=0.6       # used in calibration, to be changed if more realistic 
    texp=0.8     # calibration is necessary !! 
    nexp=0

    ###### some parameters #######################################################

    model=2        # own extension of Kingston's formula


    if shot<=34800:
        diode=2
        INDEX=r_[0:180,180+array([2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,20,19])]
    else:     
        diode=1
        INDEX=r_[0:20, 39:19:-1,59:39:-1,60:80, 80:100,119:99:-1,139:119:-1,140:160,179:159:-1,199:179:-1]
        


    path = os.path.dirname(os.path.realpath(__file__))

    if diode==1:
        #disp(' Using IRD diode geometry')
        cat_default = loadmat(path+'/cat_defaults2008_2.mat')

    else:
    #    disp(' Using CENTRONIC diode geometry')
        cat_default = loadmat(path+'/cat_defaults2001.mat')




    if diode==1:
        w=10
        D=450
        dd=6e-3       # for IRD choose SiO2 filter
        dp=0
        dp0=0.5    # 'standard' dead layer thickness
        L0=60              # 'standard' diffusion length
        filt=9      # only SiO2
        
    else:
        # data of centronic LD20-5T from Anton's get_detector.m
        w=3        # Width of depletion zone
        D=380      # absorber thickness
        dd=0.055   # Si3N4 thickness
        dp0=0.5    # 'standard' dead layer thickness
        L0=200             # 'standard' diffusion length
        filt=4      # only Si3N2
    
    fans=ones(10)
    actn=r_[0:200]

    
    th=arctan(diff(cat_default['ychord'].T)/diff(cat_default['xchord'].T))*180/pi
    nl=len(th)
    # determine angle of incidence (referred to detector surface normal)
    # for every detector

    thp=squeeze(th)
    thp[thp < 0]=180+thp[thp < 0]
    
    thdet = repeat(cat_default['vangle'][0], 20)

    angles=thdet-thp
    mist= (angles<0) & (abs(angles)>90)
    angles[mist]+= 180
    th_inc=angles*pi/180

    # --- now get diode parameters -----------------------------------------------

    if diode==1:
        result=loadtxt(path+'/new_xtomo_calibration_2008.res')
    else:
        #  disp('use 2001 calibration')
        result=fopen(path+'/cat_calibs_30kv.res') #Old caibration file
    
    #fid=fopen('CALIBRATION/eff_calibs.res')


    arraynr=result[actn,1]
    diodenr=result[actn,2]
    dp=result[actn,0] #??
    L=result[actn,4] #??
            

    dp/=cos(th_inc)**0.8  #??

    ###### get/set spectral distributions #######################################

    # According to M Anton's get_spectrum.m and get_filters.m (spec=3,filt=4):
    Epts=100         # number of points of the energy vector
    Eup = min(10*kT,50)
    KEV = linspace(0.01,Eup,Epts)
    al = TOMOLinAbs(['SI',],KEV)       

    if diode==1:
        fmu=TOMOLinAbs(['O ','SI','SI'],KEV)     # IRD
    else:
        fmu=TOMOLinAbs(['N ','SI','SI'],KEV)    #filt = 4    CENTRONIC
    
    linbe=TOMOLinAbs(['BE',],KEV)

    r=linspace(0,0.99,50)[:,None]
    psi=ones_like(r)*(1-r**2)
    Te=psi**texp
    N2=psi**nexp
    kT1=1/(kT*Te)
    

    EDIST=sum(sqrt(kT1)*exp(maximum(-kT1*KEV,-100)),0)
    EDIST=EDIST/sum(EDIST)


    tau = exp(-47*linbe) 
    EDIST = tau*EDIST/sum(tau*EDIST)

    #EDIST=ones(nl,1)*EDIST

    ####### now calculate spectrum averaged efficiencies and corrections #########

    #  according to M.Anton's eta_spec_av.m

    if diode==1:
        fd=c_[1162*dd*ones(nl), dd/3*ones(nl), dp]       # IRD
    else:
        fd=c_[1638*dd*ones(nl),0.6355*dd*ones(nl),dp] # centronic
    



    COS1=(1/cos(th_inc))[:,None]
    woft=w*COS1
    Doft=D*COS1
    FD=fd*COS1   
    weigh=EDIST*exp(dot(-FD,fmu))
    LD1=(cos(th_inc)/L)[:,None]

    aeta=1-exp(-woft*al)+al/(al+LD1)*(exp(-woft*al)-exp(-Doft*al)*exp(-(Doft-woft)*LD1))
    eta=sum((aeta*weigh),1)

    ###############  

    nil = eta < 1e-3*max(eta)
    eta[nil]=1
    # Return eta according to cabling. 2012.02.24
    eta=eta[INDEX]


    return eta





def TOMOLinAbs(filt,KEV):

    #---------[ANTON.EFFICIENCY]
    #
    #   function LINA=LINABS(filt,KEV)
    #
    #   arguments    FILT:    column vector of elements like ['BE'..]
    #                KEV:     vector of photon energies in keV
    #   returns      LINABS:  a matrix with as much rows as FILT and as much
    #                         columns as KEV 
    #   calculates linear absorption coefficients (unit 1/mu)
    #   in the energy range given by KEV
    #   cross section calculations are based on
    #   SX-ray cross-section data in file [ANTON.MATLAB]CROS2.M
    #                                M. Anton Jun 1993
    #   the whole thing is extracted from MICHAEL DUTCH's diode_resp.m
    #   modif for matlab4 23.6.94
    #
    #-----M.ANTON-----------------------------------------------------------------

    #determine number of filters used (nfilt)

    nfilt = len(filt)
    nener = argmax(KEV)+1

    #initialise LINABS

    LINA=zeros((nfilt,nener))


    #convert array of energies from keV to eV and to array of logarithms

    ep=KEV*1000         
    E=ep
    LOGE=log10(E)

    #----------LOOP OVER FILTERS--------------
    for nf,f in enumerate(filt):

        #initialise cross-section
        excr=100*ones_like(ep)

        #get cross-section data
        [RHO,NER,E1,E2,A1,A2,A3,A4] = TOMOCros2(f)

        #-- energies lower than first E1 --
        ind = KEV < E1[0]
        if any(ind):
            #**disp('below first E1')
            excr[ind] = ((A4[0]*LOGE[ind]+A3[0])*LOGE[ind]+A2[0])*LOGE[ind]+A1[0]
        
        
        for n in range(NER): # loop over each energy interval
            #-- find energies in the current range --
            ind=(KEV>=E1[n])&(KEV<=E2[n])
            if any(ind):
                excr[ind]=((A4[n]*LOGE[ind]+A3[n])*LOGE[ind]+A2[n])*LOGE[ind]+A1[n]
        
            #-- energies in-between this range + next --
            if n<NER-1:
                ind=(KEV>E2[n])&(KEV<E1[n+1])
                if any(ind):
                    #**disp(['in-between range ',num2str(n),' and ',num2str(n+1)])
                    ea=log10(E2[n]*1000)# low-E side of gap
                    eb=log10(E1[n+1]*1000)# high-E side of gap
                    excra=((A4[n]*ea+A3[n])*ea+A2[n])*ea+A1[n]
                    excrb=((A4[n+1]*eb+A3[n+1])*eb+A2[n+1])*eb+A1[n+1]
                    excr[ind]=excra+(KEV[ind]-E2[n]*ones(size(ind)))*(excrb-excra)/(E1[n+1]-E2[n])
     
        #check that excr is defined at all energies
        ind = excr > 99
        if any(ind):
            raise Exception('ERROR: unable to calculate cross-section at some energies')
            #str='       E(j),j='
            #for i=1:length(ind)
                #str=[str,num2str(ind(i)),',']end
                #disp(str)
            #return
        #end

        cr=10.0**(excr+3)
        #print LINA.shape,nf,nfilt,nener
        LINA[nf,:]=1.0e-4*RHO*cr

        ind = LINA>50.0
        if any(ind):
            LINA[ind]=LINA[ind]*0+50.0


    return LINA
#--------END LOOP OVER FILTERS-----------




def TOMOCros2(filter):

    # ----[anton.efficiency]
    #
    # function [RHO,NER,E1,E2,A1,A2,A3,A4]=cros2(filter)
    #
    # CROS.DAT, 23/6/1993
    # routine to provide cross-section data originally 
    # stored in Fortran database file CROS.DAT
    #               M. J. Dutch Feb 1992
    #
    # SYNTAX:       function [RHO,NER,E1,E2,A1,A2,A3,A4]=cros2(filter)
    #
    #  INPUT:
    #               FILTER = string containing name of filter material
    #                      e.g. 'BE' , 'AL' , 'FE' etc 'O ', 'N '
    # OUTPUTS:
    #               RHO = density in g/cm^3
    #               NER = number of energy ranges
    #               E1,E2 = limits of energy range (keV)
    #               A1,A2,A3,A4 = coeffs of polynomial fit
    #                     to log10(XS) vs log10(E)

    #casesen off
    # TOTAL CROSS-SECTION DATA FOR AL   ((HENKE , VEIGELE), INPUT FOR TRCO1
    if filter=='AL':
        RHO=2.6940
        NER=3
        E1=[0.109,0.1177,1.5610]
        E2=[0.1177,1.5590,50.00]
        A1=[3.8484,0.0003,0.0001]
        A2=[-0.8928,3.1349,4.7428]
        A3=[0.00,-1.0377,-1.9560] 
        A4=[0.00,0.00,0.1678]

    # TOTAL CROSS-SECTION DATA FOR AU    (HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='AU': 
        RHO=19.2900 
        NER=5
        E1=[0.1180,0.3418,0.5789,0.7690,2.206]
        E2=[0.3404,0.5775,0.7595,2.206,10.00]
        A1=[0.0000,0.0000,0.0000,0.000,-0.0001]
        A2=[-3.7147,-0.1215,-14.8078,-0.1987,-0.7892]
        A3=[3.2407,0.9054,11.3242,0.9428,1.0066]
        A4=[-0.6254,-0.2637,-2.1094,-0.2654,-0.2169]

    # TOTAL CROSS-SECTION DATA FOR B  (HENKE,SCHATTENBURG),INPUT FOR TRCO1
    elif filter=='B ': 
        #  RHO=2.535 NER=2 
        #  E1=[0.109,0.188]
        #  E2=[0.188,1.487]
        #  A1=[-5.7441,0.7169] A2=[11.4904,4.1999]
        #  A3=[-5.3495,-2.0451] A4=[0.6609,0.1918]

        # conv_cross_new: B 
        # converted from VEIGELE 1973
        RHO=2.5350e+00  
        NER=1 
        E1=[ 1  ]
        E2=[ 1000  ]
        A1=[ 33.68  ]
        A2=[ -19.41  ]
        A3=[ 3.343  ]
        A4=[ -0.1932  ]



    # TOTAL CROSS-SECTION DATA FOR BE   ((HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='BE': 
        RHO=1.8450 
        NER=2 
        E1=[0.109,0.1121]
        E2=[0.1121,50.00]
        A1=[3.5889,-0.0005]
        A2=[-1.4088,5.9586] 
        A3=[0.0000,-3.0330] 
        A4=[0.00,0.3381]

    # TOTAL CROSS-SECTION DATA FOR  C  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='C ': 
        RHO=2.2500 
        NER=2
        E1=[0.109,0.2838]
        E2=[0.2838,10.00]
        A1=[12.8122,0.0001]
        A2=[-12.3355,4.7466]
        A3=[4.8952,-2.1186] 
        A4=[-0.7900,0.1917]

    # TOTAL CROSS-SECTION DATA FOR  CO  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='CO': 
        RHO=8.9000 
        NER=4
        E1=[0.109,0.7936,0.9256,7.71]
        E2=[0.7786,0.9256,7.708,15.00]
        A1=[-5.0159,9.0624,0.0001,0.0000] 
        A2=[9.1499,-2.7201,4.1816,1.2259]
        A3=[-3.6783,0.000,-1.6922,-0.0717]
        A4=[0.3968,0.000,0.1361,-0.0696]

    # TOTAL CROSS-SECTION DATA FOR CR   ((HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='CR': 
        RHO=7.1800 
        NER=4
        E1=[0.109,0.5837,0.6946,5.99]
        E2=[0.5745,0.6946,5.988,50.00]
        A1=[-2.1244,8.2965,0.0001,0.0001]
        A2=[5.6865,-2.5053,4.2526,3.5562]
        A3=[-2.3605,0.000,-1.756,-1.2615]
        A4=[0.2295,0.000,0.1450,0.0800]

    # TOTAL CROSS-SECTION DATA FOR  CU  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='CU': 
        #  RHO=8.9400 NER=5
        #  E1=[0.109,0.1198,0.9510,1.0961,8.9800]
        #  E2=[0.1198,0.9311,1.0961,8.9780,20.000]
        #  A1=[4.1270,0.1362,9.6776,0.000,0.000]
        #  A2=[-1.1549,2.9734,-2.8823,4.6066,0.2483]
        #  A3=[0.000,-1.1733,0.000,-1.9188,0.3653]
        #  A4=[0.000,0.0603,0.000,0.1680,-0.1169]
        # in the following a fit for E>20keV has been added...
        #conv_cross_new: CU
        #RHO=8.9400e+00  NER=1 
        #E1=[ 20  ]
        #E2=[ 1000  ]
        #A1=[ 36.02  ]
        #A2=[ -14.59  ]
        #A3=[ 1.502  ]
        #A4=[ -3.1164e-02  ]
        RHO=8.9400 
        NER=6
        E1=[0.109,0.1198,0.9510,1.0961,8.9800,20]
        E2=[0.1198,0.9311,1.0961,8.9780,20.000,1000]
        A1=[4.1270,0.1362,9.6776,0.000,0.000,36.02 ]
        A2=[-1.1549,2.9734,-2.8823,4.6066,0.2483,-14.59 ]
        A3=[0.000,-1.1733,0.000,-1.9188,0.3653,1.502 ]
        A4=[0.000,0.0603,0.000,0.1680,-0.1169,-3.1164e-02]

    # TOTAL CROSS-SECTION DATA FOR F  (HENKE,SCHATTENBURG),INPUT FOR TRCO1
    elif filter=='F ': 
        RHO=1.579E-3
        NER=2
        E1=[0.109,0.6854]
        E2=[0.6854,1.487]
        A1=[-2.1974,8.6921] 
        A2=[7.9743,-3.9884]
        A3=[-3.9246,0.8677] 
        A4=[0.4855,-0.1398]

    # TOTAL CROSS-SECTION DATA FOR  FE  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='FE': 
        RHO=7.8600
        NER=4
        E1=[0.109,0.7211,0.8461,7.1130]
        E2=[0.7081,0.8461,7.111,50.00]
        A1=[-3.9912,5.9994,0.0001,0.000]
        A2=[7.9782,-1.6779,4.0975,3.8360]
        A3=[-3.2537,0.000,-1.6493,-1.3865]
        A4=[0.3461,0.00,0.1298,0.0950]

    # TOTAL CROSS-SECTION DATA FOR  GA  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='GA': 
        RHO=5.9180
        NER=5
        E1=[0.109,0.1581,1.1423,1.2977,10.3710]
        E2=[0.1581,1.1154,1.2977,10.3690,30.00]
        A1=[-11.0044,-9.2714,12.5976,0.0001,-0.0001]
        A2=[19.6463,13.6901,-3.7934,4.1836,6.6275]
        A3=[-9.6498,-5.1954,0.000,-1.6653,-2.6775]
        A4=[1.5177,0.5626,0.000,0.1322,0.2464]

    # TOTAL CROSS-SECTION DATA FOR  GE  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    # energy range 7 calculated from VEIGELE 1973 ,conv_cross_new.m
    elif filter=='GE': 
        RHO=5.3080
        NER=7
        E1=[0.109,0.1208,0.1800,1.2478,1.4143,11.101,30]
        E2=[0.1208,0.1800,1.2167,1.4143,11.099,30.00,1000]
        A1=[2.5223,-4.7922,-7.9763,9.1806,0.0001,0.000,62.99]
        A2=[-0.3729,6.2178,12.058,-2.6787,4.6493,4.5597,-29.52]
        A3=[0.000,-1.2333,-4.509,0.000,-1.9124,-1.7059,4.261]
        A4=[0.000,-0.1153,0.4688,0.000,0.1654,0.1327,-0.2013]



    # TOTAL CROSS-SECTION DATA FOR H  (HENKE,SCHATTENBURG),INPUT FOR TRCO1
    elif filter=='H ': 
        RHO=8.381E-5 
        NER=1
        E1=[0.109]
        E2=[1.487]
        A1=[6.9219] 
        A2=[-2.7457]
        A3=[-0.0803]
        A4=[-0.0044]

    # TOTAL CROSS-SECTION DATA FOR HE (HENKE,SCHATTENBURG),INPUT FOR TRCO1
    elif filter=='HE': 
        RHO=1.664E-4
        NER=1
        E1=[0.109]
        E2=[1.487]
        A1=[2.5115] 
        A2=[2.8776] 
        A3=[-2.1088]
        A4=[0.2445]

    # TOTAL CROSS-SECTION DATA FOR MG (HENKE,SCHATTENBURG),INPUT FOR TRCO1
    elif filter=='MG': 
        RHO=1.735 
        NER=2
        E1=[0.109,1.305]
        E2=[1.305,1.487]
        A1=[-8.1446,9.3449]
        A2=[13.7949,-2.7466]
        A3=[5.5844,0.000] 
        A4=[0.6288,0.000]

    # TOTAL CROSS-SECTION DATA FOR  MN  ((HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='MN': 
        RHO=7.3000 
        NER=4
        E1=[0.109,0.6514,0.769,6.54]
        E2=[0.6403,0.769,6.538,15.00]
        A1=[-3.6611,7.7481,0.0001,0.000] 
        A2=[7.5229,-2.3014,3.8449,-0.3190]
        A3=[-3.083,0.000,-1.5094,0.6956]
        A4=[0.3254,0.000,0.1092,-0.1663]

    elif filter=='MO':
    # TOTAL CROSS-SECTION DATA FOR  MO  (VEIGELE,AT.DATA TABLES 5 (1973))
    # fit calculated by [anton.matlab]conv_cross_new: 
        RHO=1.0200e+01  
        NER=6 
        E1=[ 1,2.52 , 2.625 , 2.866 , 17.48 , 20   ]
        E2=[ 2.52 , 2.625 , 2.866 , 17.48  ,20 , 1000   ]
        A1=[ -1.908 , 9.781 , 9.928  ,8.012 , 9.573  ,-15.71   ]
        A2=[ 7.147,  -2.773 , -2.774,  -1.075 , -2.662 , 16.11   ]
        A3=[ -3.05  ,0,  0  ,-0.4868,  0,  -4.443   ]
        A4=[ 0.3211 , 0  ,0 , 4.6994e-02,  0 , 0.3463   ]




    elif filter=='N ':
    # TOTAL CROSS-SECTION DATA FOR  N  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    #  RHO=0.0012 NER=2
    #  E1=[0.109,0.402]
    #  E2=[0.402,20.00]
    #  A1=[5.736,0.0001]
    #  A2=[-2.1165,4.8138]
    #  A3=[0.1510,-2.1299]
    #  A4=[-0.0601,0.1942]
    # total cross section data for N (Veigele 1973) FILE:CROSS_RAW.M
    # converted by CONV_CROSS_RAW.M
        RHO=1.1650e-03  
        NER=1  
        E1=[1.000  ]
        E2=[1000.000  ]
        A1=[23.8384  ]
        A2=[-11.8003  ]
        A3=[1.5448  ]
        A4=[-0.0587  ]


    # TOTAL CROSS-SECTION DATA FOR  NI  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='NI': 
        RHO=8.8780 
        NER=5
        E1=[0.109,0.1118,0.8719,1.0081,8.334]
        E2=[0.1118,0.8547,1.0081,8.332,50.00]
        A1=[3.8225,-4.1088,9.4278,0.0001,0.0001]
        A2=[-1.0203,8.000,-2.8185,4.1576,3.1790]
        A3=[0.000,-3.165,0.000,-1.665,-1.0741]
        A4=[0.00,0.3221,0.000,0.1318,0.0592]

    elif filter=='O ': 
    # TOTAL CROSS-SECTION DATA FOR O  (HENKE,SCHATTENBURG),INPUT FOR TRCO1
    #  RHO=1.331E-3 NER=2
    #  E1=[0.109,0.532]
    #  E2=[0.532,1.487]
    #  A1=[2.0547,-0.600] A2=[2.8253,5.3799]
    #  A3=[-1.9383,-2.2924] A4=[0.2326,0.2132]
    # total cross section data for O (Veigele 1973) FILE:CROSS_RAW.M
    # converted by CONV_CROSS_RAW.M
        RHO=1.3310e-03  
        NER=1  
        E1=[1.000  ]
        E2=[1000.000  ]
        A1=[19.2997  ]
        A2=[-8.4228  ]
        A3=[0.7659  ]
        A4=[-0.0016  ]


    # TOTAL CROSS-SECTION DATA FOR SC  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='SC': 
        RHO=3.0000 
        NER=4
        E1=[0.109,0.4067,0.5004,4.494]
        E2=[0.4022,0.5004,4.492,10.00]
        A1=[-10.8303,8.4591,0.0002,0.000]
        A2=[16.0399,-2.6385,4.0693,-0.1961]
        A3=[-6.5481,0.000,-1.6704,0.7102]
        A4=[0.7954,0.000,0.1314,-0.1812]

    # TOTAL CROSS-SECTION DATA FOR SI  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='SI': 
        RHO=2.3200
        NER=4
        E1=[0.109,0.1487,1.8400,1.8400]
        E2=[0.1487,1.838,1.838,50.00]
        A1=[61.3204,0.0003,0.0001,0.0001]
        A2=[-87.4143,3.0875,0.6117,4.8124]
        A3=[43.3307,-1.0054,-0.1925,-1.9686]
        A4=[-7.2176,0.000,0.000,0.1681]

    # 1.TEST MIT CR-DATEN FUER SI (VEIGELE), SELBST KONVERTIERT
    #RHO=2.3200e+00  NER=2
    #E1=[1.000  1.839  ]
    #E2=[1.839  50.000  ]
    #A1=[0.0000  -13.9975  ]
    #A2=[5.5710  15.6411  ]
    #A3=[-2.7448  -4.7379  ]
    #A4=[0.3034  0.4023  ]
    # 2.TEST MIT CR-DATEN FUER SI (VEIGELE), SELBST KONVERTIERT
    #conv_cross_new: SI
    #RHO=2.3200e+00  NER=2 
    #E1=[ 1  1.839   ]
    #E2=[ 1.839  1000   ]
    #A1=[ 9.437  9.087   ]
    #A2=[ -3.437  -0.4497   ]
    #A3=[ 0.1194  -1.084   ]
    #A4=[ 0  0.1325   ]


    # TOTAL CROSS-SECTION DATA FOR TI  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='TI': 
        RHO=4.5400 
        NER=4
        E1=[0.109,0.4615,0.5637,4.967]
        E2=[0.4555,0.5637,4.965,10.00]
        A1=[-5.2607,7.7316,0.0002,0.000]
        A2=[9.1100,-2.3478,3.9229,0.1739]
        A3=[-3.669,0.00,-1.5758,0.5023]
        A4=[0.3985,0.00,0.1177,-0.1515]

    # TOTAL CROSS-SECTION DATA FOR  V  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='V ': 
        RHO=6.1000
        NER=4
        E1=[0.109,0.5205,0.6282,5.466]
        E2=[0.5129,0.6282,5.464,20.00]
        A1=[-4.446,5.4826,0.0001,0.00]
        A2=[8.25,-0.7750,3.954,3.305]
        A3=[-3.3577,-0.2689,-1.5884,-1.1317]
        A4=[0.3619,0.000,0.1199,0.0623]


    # TOTAL CROSS-SECTION DATA FOR  W  (VEIGELE 1973)
    # eigener fit mit conv_cross_new
    #conv_cross_new: W 
    elif filter=='W ':
        RHO=1.9300e+01
        NER=11
        E1=[1,1.809,1.872,2.281,2.575,2.82,10.21,11.54,12.1,59.32,69.53]
        E2=[1.809,1.872,2.281,2.575,2.82,10.21,11.54,12.1,50,69.53,1000]
        A1=[8.866,9.308,-2.244,9.438,9.316,7.121,
                10.28,10.45,3.03,-71.42,-38.85]
        A2=[-2.812,-2.682,4.3,-2.642,-2.594,
                -0.6036,-2.715,-2.722,2.786,31.19,28.99]
        A3=[3.3677e-02,0,-1.038,0,0,-0.5887,0,0,-1.351,-3.506,-6.685]
        A4=[0,0,0,0,0,5.7917e-02,0,0,0.1107,0,0.4693]



    # TOTAL CROSS-SECTION DATA FOR  ZN  (V(HENKE , VEIGELE), INPUT FOR TRCO1
    elif filter=='ZN': 
        RHO=7.1150
        NER=5
        E1=[0.109,0.1359,1.0428,1.1936]
        E2=[0.1359,1.0197,1.1936,9.658,20.00]
        A1=[3.7437,-6.4551,9.7221,0.0001,0.000]
        A2=[-1.0203,10.5707,-2.8823,4.1237,1.803]
        A3=[0.0341,-4.0562,0.000,-1.6373,-0.3813]
        A4=[0.00,0.4237,0.00,0.1284,-0.0268]

    else:
        raise Exception('CROS2: Unknown material '+filter)


    return RHO,NER,E1,E2,A1,A2,A3,A4





#TOMOCal(40412)




