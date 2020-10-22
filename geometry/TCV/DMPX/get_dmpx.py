#!/usr/bin/env python
# -*- coding: utf-8 -*-

#ssh todstrci@lac911.epfl.ch -L 8002:tcv1.epfl.ch:8000

from numpy import *
#import pmds
#from matplotlib.pylab import *
#import ././tqdm
import os.path 

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

#parent_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
#print parent_path
#sys.path.append(parent_path)
#import tqdm
from tqdm import tqdm, trange




#mpx.top.signal
class mpx_:
    def __init__(self):
        class top_:
            #__init__(self):
                #class geom:
                    #__init__(self):

            pass
        self.top = top_()
        class bot_:pass
        self.bot = bot_()
        
        
    
class empty_class :pass
    
#class also_: pass

# MPXDATA allows you to load the DMPX calibrated signal, detectors high voltage, geometry,
# filters, an estimation of rhopsi and rhovol corresponding to the lines-of-sight as well as
# the detectors efficiency. MPXDATA works only for shots > 20030.
#
# SYNTAX
# [mpx] = mpxdata(shot,action,varargin)
#
# INPUTS
# shot          = Shot number
# action    = String containing the key character for the actions you want to execute. The actions are:
#       's':    Loads the MPX calibrated signal.
#       Channels are sorted from HFS to LFS.
#       Warning:
#       * For shots < 27127, the calibration does not include the voltage applied to the detector.
#       * For shots > 27128, it does and the signal of two shots with different voltage can be compared.
#       'v':    Loads the detector high voltage value [V] and corresponding gain.
#       'g':    Loads the detector geometry (slit and wires position from HFS to LFS) [m].
#               Creates also a dcd object to be used with the psitbx.
#       'f':    Loads the filters in front of the detectors (thicknesses in [m]).
#               + the position of the mobile absorber holder
#               + the gas in the detector
#               Use the 'pos' option to force the position of the mobile absorber holder.
#   
#       'e':    Loads the detector efficiency (probability that the detection gas absorbs an incoming photon)
#               as a function of the energy of incoming photons. The MPX output signal is proportional to
#               the detector efficiency if the escape peak is neglected.
# varargin      = Parameters for the choosen actions. You can specify:
#       'time'          Time interval (ex: [0.1 0.9]) for the signal. Default = full shot.
#               'freq'          Acquisition frequency: 
#                               0: slow 20kHz
#                               1: fast 200kHz (if available, else slow)
#                               2: ask if both are available, else slow (default)
#               'detec'         Detector selection:
#                               'top':  top detector
#                               'bot':  bottom detector
#                               'both': both detectors (default)
#               'chords'        list of chords to consider (ex: [32 48 63]). Default = all chords considered.
#                               [1:64] -> top detector
#                               [65:96] -> bottom detector
#                               NB: The 'chords' option overrules the 'detec' option.
#               'pos'           To force the position of the mobile absorber holder for the efficiency calculation.
#               'rtime'         Time interval for the rho_psi calculation (ex: [0.1 0.9]). It then calculates
#                               a value of rho on each psitbx time. Default = value used for action 's'.
#               'liuqe'         Liuqe version (1, 2, 3) for the rho_psi calculatin. Default = 1.
#               'psi'           Flux description in cylindrical coordinates ('01' psitbxpsi object) for the rho_psi
#                               calculation. Optional, useful if you already have psi loaded in your workspace.
#                               Ex: psi=psitbxtcv(shot,time,'01').
#               'opt'           Option for the rho_psi calculation
#                               -1: to get rho from -1 (HFS) to 1 (LFS). Default.
#                               +1: to get rho from 0 to 1.
#               'vol'           Non zero value to get also rho_vol. Default = 0.
#
# OUTPUT
# mpx           = structure containing the desired informations
#
# EXAMPLES
# [mpx] = mpxdata(27991,'svgf')
# [mpx] = mpxdata(27991,'sr','time',[0.9 1.1],'freq',1,'detec','bot','vol',1)
# [mpx] = mpxdata(24848,'g','detec','top')
# [mpx] = mpxdata(36989,'s','time',[0.9 0.91],'chords',[1 32 36 65 94])
#
# Written by Y. Camenen, CRPP, EPFL.
# Modified by L. Curchod, CRPP, EPFL:   15.07.2008      Correction of the error when loading rho_psi between -1 and 1.
# Modified by L. Curchod, CRPP, EPFL:   14.08.2008      Use a new mpx_calib_dec05_bis.mat file calculated with data for 2450 V as well. 
# Modified by L. Curchod, CRPP, EPFL:   15.08.2008      Corrected to sort the gain coefficients for the oct. 2005 calibration case.
# Modified by L. Curchod, CRPP, EPFL:   03.03.2009      New path to the calibration files on PC70
# Modified by Y. Camenen, Warwick, UK:  10.07.2009      Add option to load selected chords
# Modified by J. Kamleitner, CRPP, EPFL:    3.2013  moved calibration files to crpptbx to be independent of user home directories





def mpxdata(connect, shot,action,**kwarg):


    #
    # Check for non-stored shots
    #

    if shot<20030:
        raise Exception('This program works only for shot >= #20030.')
    elif shot>=34988 and shot<=35438:
        raise Warning('DMPX signal has been lost for shots #34988 to #35428.')
        ##action(strfind(action,'s'))='' # Withdraw the action of signal loading
    elif shot>=39113 and shot<=39123:
        #pass
        raise Warning('DMPX signal has been lost for shots #39113 to #39123.')
            #action(strfind(action,'s'))='' # Withdraw the action of signal loading
    
    diag_path = os.path.dirname(os.path.realpath(__file__))+'/'
    print(('\n MPXDATA for shot %d.'%shot))

    #
    # Define defaults values
    #
    definitions = empty_class()

    #mpxpath='/home/matlab/crpptbx-7.6.0/dmpx/'
    definitions.action = ['s',]
    definitions.time='**'
    definitions.rtime=[0,10]
    definitions.det='both'
    definitions.freq=2 # (0=slow, 1=fast if available, 2=ask, if fast is available)
    definitions.chords=r_[0:97] 
    definitions.liuqe=1
    definitions.opt=-1
    definitions.vol=0
    definitions.top_nchords=64
    definitions.bot_nchords=32
    abso = empty_class()
    
    
    abso.material= ['Al','Be','Be','nothing']
    abso.thickness=r_[308,550,125,NaN]*1e-6


    mpx = mpx_()
    #
    # Default values
    #
    if not 'detec' in kwarg:
        detec = definitions.det
        
    #if not 'action' in kwarg:
        #action = definitions.action
    action = list(action)
#
    for a in kwarg:                    #Loads all prepared data
        exec(a+" = kwarg['"+a+"']")
    
    #
    # Check 'chord' and 'detec' consistency
    #
    
    if 'chords' in kwarg: # If 'chords' is specified, it overrules 'detec'
        if 'detec' in kwarg:
            old_detec = detec
        
        if any(chords<0) or any(chords>definitions.top_nchords+definitions.bot_nchords):
            raise Exception('Option ''chords'' must only include numbers between 1 and %d'%(definitions.top_nchords+definitions.bot_nchords))
        
        istop=any(chords<=definitions.top_nchords)
        isbot=any(chords>definitions.top_nchords)
        
        if istop+isbot*2 == 1: detec='top'
        if istop+isbot*2 == 2: detec='bot'
        if istop+isbot*2 == 3: detec='both'

        if 'old_detec' in kwarg and old_detec!= detec:
            raise Warning('Option "detec" has been overruled from "'+old_detec+'" to "'+detec+'" by the "'+chords+'" selection.')
        
    else:
        if detec == 'both':
            chords=definitions.chords
        elif detec ==  'top':
            chords=r_[0:definitions.top_nchords]
        elif detec == 'bot':
            chords=r_[0:definitions.bot_nchords]+definitions.top_nchords
        
    
    if shot<26609:
        if  detec == 'both': 
            if shot<26555:
                print('- Load data for top detector only: Duplex MPX not installed for shot < #26555.\n')
            else:
                print('- Load data for top detector only: No bottom detector HV source for shot < #26609.\n')
            
            detec='top'
        elif detec == 'bot':
            if shot<26555:
                raise Exception('No bottom detector: Duplex MPX not installed for shot < #26555.')
            else:
                raise Exception('No bottom detector: No bottom detector HV source for shot < #26609.')
            
        
    
    #
    # To load the MPX geometry before computing rho
    #
    i_r = ''.join(action).find('r')
    if i_r != -1:
        i_g=''.join(action).find('g')
        if i_g == -1:
            action += 'g'
            i_g=len(action)-1
        
        if i_g>i_r:
            action[i_r]='g'
            action[i_g]='r'
        
    
    #
    # To load the MPX high voltage value before calibrating the signal
    #
    i_s = ''.join(action).find('s')
    if i_s != -1:
        i_v=''.join(action).find('v')
        if i_v == -1:
            action+= 'v'
            i_v=len(action)-1
        
        if i_v>i_s:
            action[i_s]='v'
            action[i_v]='s'
       
    
    #
    # To load the MPX filters before computing the energy response
    #
    i_e = ''.join(action).find('e')

    if i_e != -1:
        i_f=''.join(action).find('f')

        if i_f == -1:
            action += 'f'
            i_f=len(action)-1
        
        if i_f>i_e:
            action[i_e]='f'
            action[i_f]='e'
        
    #print 'action', action
    #############################################################
    #                     LOAD DATA                             #
    #############################################################
    for act in action:
        
        ###############################################
        #       Load the DMPX calibrated signal       # 
        ############################################### 
        if act == 's':
            #print('\ncase s\n')
            #pmds.mdsopen('tcv_shot',shot)
            connect.openTree('tcv_shot',shot)

            
            #pmds.mdsvalue('reset_public()')
            #        
            # Some initalisations and tests
            #
            print('- Load the MPX calibrated signal.\n')
            if not 'freq' in kwarg:
                freq=definitions.freq
            
            if 'time' in kwarg:
                if shot<20939: # time base was incorrect
                    time[0]=(time[0]+0.04)/0.9697-0.04
                    time[1]=(time[1]+0.04)/0.9697-0.04
                
            else:
                time=definitions.time
            
            #BUG použít k loudění tu rychlou metodu!!
            
            #echant=['[' num2str(time(1)) ':' num2str(time(2)) ':*]']
            
            echant = '[%s:%s:*]'%(time[0],time[1])
            if  shot<24087 or  shot>24725:
                DTNE1=r'\ATLAS::DT100_NORTHEAST_001:CHANNEL_0'
                DTNE2=r'\ATLAS::DT100_NORTHEAST_002:CHANNEL_0'
                DTNE3=r'\ATLAS::DT100_NORTHEAST_003:CHANNEL_0'
                DTNE1_fast=r'\ATLAS::DT100_NORTHEAST_001:SELECTED:CHANNEL_0'
                DTNE2_fast=r'\ATLAS::DT100_NORTHEAST_002:SELECTED:CHANNEL_0'
                DTNE3_fast=r'\ATLAS::DT100_NORTHEAST_003:SELECTED:CHANNEL_0'
            else:
                DTNE1=r'\ATLAS::DT100_SOUTHWEST_001:CHANNEL_0'
                DTNE2=r'\ATLAS::DT100_SOUTHWEST_002:CHANNEL_0'
                DTNE1_fast=r'\ATLAS::DT100_SOUTHWEST_001:SELECTED:CHANNEL_0'
                DTNE2_fast='r\ATLAS::DT100_SOUTHWEST_002:SELECTED:CHANNEL_0'
            
            #     
            # Test the DTACQ trigering mode
            #
            #mode1=pmds.mdsvalue(DTNE1[:-9]+'MODE')
            #mode2=pmds.mdsvalue(DTNE2[:-9]+'MODE')
            mode1=int(connect.get(DTNE1[:-9]+'MODE'))
            mode2=int(connect.get(DTNE2[:-9]+'MODE'))
            
            if detec in ['bot', 'both']:
                #mode3=pmds.mdsvalue( DTNE3[:-9]+'MODE') 
                mode3=int(connect.get( DTNE3[:-9]+'MODE'))
            else:
                mode3=4
            
            mode=mode1*mode2*mode3 #mode=64 ->ok
            # if shot>=24087&shot<25095
            if shot>=24087 and mode!=64:
                print('Random temporal gap (5 to 25 us) between the two or three MPX acquisition cards.')
                print('Random temporal gap (0.1 to 0.2ms) between DTACQ and TCV.')
                raise Warning('DTACQ not in mode 4')

            #
            # Acquisition frequency
            #
            #if shot<34988: # After TCV BIG OPENING 2007-2008, fast data are in standard nodes
                #if freq!=0:
                    #test_fast=pmds.mdsvalue(r'node_size("\\'+DTNE1_fast+'01")') #test whether the first channel of the first card has fast acquisition
                #if isscalar(test_fast):
                    #if test_fast>0:
                        #flag1=0
                        #if freq==2,
                            #freq=input(['Do you want standard acquisition f = 20kHz (0) or fast acquisition at f = 200 kHz (1) ?'])
                        #end
                    #else
                        #flag1=1
                    #end
                #else
                    #flag1=1
                #end
                #if(flag1)
                                #print('Only standard acquisition (f=20kHz) available')
                                #freq=0
                #end
                #end
            #else # After TCV BIG OPENING 2007-2008, fast data is always available in standard nodes
                #if freq==2,
                        #freq=input(['Do you want decimated data at f = 20kHz (0) or full data at f = 200 kHz (1) ?'])
                #end
            
            freq = 1 #fastest!!
            
            #        
            # Load, sort and calibrate TOP data
            #
            if detec in ['top','both']:
                if freq==1:
                    if shot<34988: # Before TCV BIG OPENING 2007-2008, fast data are in selected nodes
                        DTNE1=DTNE1_fast
                        DTNE2=DTNE2_fast
                    
                else:
                    if shot>=34988: # After TCV BIG OPENING 2007-2008, fast data are in standard nodes
                        echant = '[%f:%f:0.00005]'%time

                 
                #
                # Load timebase and check if node is full
                #
                
                datatest = []
                #BUG!!!
                #dmpx_files = [DTNE1_fast
                #for node in 
 
                #tvec = 0
                #tvec,data = mds_load(node,chn, chnfmt=None, tmin=-infty,tmax=infty,step=1,raw=True,
                    #remove_offset=False,offset_time_min = -infty, offset_time_max=0)
                #tvec = pmds.mdsvalue('_t=dim_of('+DTNE1+'01'+')')
                tvec = asarray(connect.get('_t=dim_of('+DTNE1+'01'+')'))
                ind_min = tvec.searchsorted(time[0])
                ind_max = tvec.searchsorted(time[1])-1
                echant = '[%d: %d]'%(ind_min,ind_max-1)

                if len(tvec) == 0:
                    print('No TOP detector signal.')
                    load_opt = False
                else:
                    mpx.top.tvec = tvec[ind_min:ind_max]
                    mpx.top.ch = r_[0:definitions.top_nchords]
                    load_opt = True
                
                
                

                
                if load_opt:
                    if shot>19923 and shot<20939:  # time base incorrect, need to rescale
                        mpx.top.tvec=(mpx.top.tvec+0.04)*0.9697-0.04

                    l_t=len(mpx.top.tvec)
                    #tmpdata=NaN*zeros(l_t,length(mpx.top.signal.dim{2}))
                    tmpdata = empty((l_t, len(mpx.top.ch)),dtype='float32')
                    
                    ch = ['%.2d'%k for k in r_[1:33]]
                    #ch=int2str(II')
                    #ch(ch==' ')='0'        
                    #
                    # Sort the channels from HFS to LFS
                    #
    
                    
                    II=zeros(64,dtype=int)
                    if shot<26555: #MPX, calib_coeff is already sorted for these shots
                        II[0:63:2]=r_[32:64]
                        II[1:64:2]=r_[32:0:-1]-1
                    else: #DMPX, calib_coeff is not sorted for these shots
                        II[0:63:2]=r_[32:64]
                        II[1:64:2]=r_[15:-1:-1,31:15:-1]
                        calib_coeff_t=calib_coeff_t[II]        # BUG Ordering calib_coeff for the top detector
                        #mpx.top.gainC=mpx.top.gainC[II]      # Ordering gains for the top detector (useful only for calibration of oct. 05
                        mpx.top.gainC = mpx.top.gainC[II]
                        #mpx.top.gainR = mpx.top.gainR[II]

                    #
                    # Load the raw data, not sorted, not calibrated   
                    # NB: only load the data for the channels specified in 'chords'
                    #
                 
     
                raw_calib = (1./(2**15-1)*10)
                for iii in tqdm(list(range(32))):
                          
                    
                    if any(where(II-iii==0)[0]-chords==0):
                        #print 'raw_of('+DTNE1+ch[iii]+')'+echant
                        #print tmpdata.shape, pmds.mdsvalue('raw_of('+DTNE1+ch[iii]+')'+echant).shape
                        #tmpdata[:,iii]=pmds.mdsvalue('raw_of('+DTNE1+ch[iii]+')'+echant)*raw_calib
                        tmpdata[:,iii]=connect.get('raw_of('+DTNE1+ch[iii]+')'+echant)*raw_calib

                        #pmds.mdsvalue(DTNE1+ch[iii]+echant)
                    
                    
                    if any(where(II-iii-32==0)[0]-chords==0):
                        #tmpdata[:,iii+32]=pmds.mdsvalue('raw_of('+DTNE2+ch[iii]+')'+echant)*raw_calib
                        tmpdata[:,iii+32]=connect.get('raw_of('+DTNE2+ch[iii]+')'+echant)*raw_calib

                    # Calibrate
                    #
                    mpx.top.signal = tmpdata[:,II]
                    mpx.top.saturated = mpx.top.signal == 10 
                    mpx.top.signal*= calib_coeff_t/mpx.top.gainC
                    #
                    # Some special cases
                    #
                    if shot>=25102 and shot<25933:
                        mpx.top.signal[:,35]=0
                        print('Channel 36 was missing for this shot') 
                        
                    elif shot>=27127 and shot<=28124:
                        #missing channels, problem with cable 2, repaired by DF in Dec 2004
                        mpx.top.signal[:,[2,63,61,59,57,55]]=0 
                        if shot>=27185: #one more channel missing !...
                            mpx.top.signal[:,43]=0
                            
                    
                    if shot>=28219 and shot<31446:
                        mpx.top.signal[:,[18, 20]]=0
                        print('Channel 19 and 21 were missing for this shot') 
                    
                
              
            #        
            # Load and calibrate BOTTOM data (no need to sort)
            #
            if detec in ['bot', 'both']:
                if freq==1:
                    if shot<34988: # Before TCV BIG OPENING 2007-2008, fast data are in selected nodes
                        DTNE3=DTNE3_fast
                else:
                    if shot>=34988: # After TCV BIG OPENING 2007-2008, fast data are in standard nodes
                        echant = '[%f:%f:0.00005]'%time # Need to decimate for 20 kHz

                    
                #BUG!!!!!!!!!!!!!
                #datatest = tdi([DTNE3,'01',echant])
                #tvec = pmds.mdsvalue('_t=dim_of('+DTNE1+'01'+echant+')')
                
                #tvec = pmds.mdsvalue('_t=dim_of('+DTNE3+'01'+')')
                tvec = asarray(connect.get('_t=dim_of('+DTNE3+'01'+')'))

                ind_min = tvec.searchsorted(time[0])
                ind_max = tvec.searchsorted(time[1])-1
                echant = '[%d: %d]'%(ind_min,ind_max-1)

                if len(tvec) == 0:
                    raise Warning('No BOTTOM detector signal.')
                    load_opt = False
                else:
                    mpx.bot.tvec = tvec[ind_min:ind_max]
                    mpx.bot.ch = r_[0:definitions.bot_nchords]
                    load_opt = True
                
                if load_opt:
                    l_t=len(mpx.bot.tvec)
                    tmpdata=empty((l_t,len(mpx.bot.ch)),dtype=single)
                    II=r_[0:32]
                    #ch = [str(k+1) for k in II]
                    ch = ['%.2d'%(k+1) for k in II]

                    #ch=int2str(II')
                    #ch(ch==' ')='0'        
                    #
                    # Load the raw data, sorted from HFS to LFS, non calibrated   
                    #
                    raw_calib = (1./(2**15-1)*10)
           
                    for iii in tqdm(II):
                        #tqdm
                        if any(definitions.top_nchords+where(II-iii==0)[0]-chords==0):
                            #tmpdata[:,iii]=pmds.mdsvalue('raw_of('+DTNE3+ch[iii]+')'+echant)*raw_calib                        
                            tmpdata[:,iii]=connect.get('raw_of('+DTNE3+ch[iii]+')'+echant)*raw_calib                        

                    # 
                    # Calibrate
                    #
                    
     
                    mpx.bot.signal = tmpdata
                    mpx.bot.saturated = mpx.bot.signal == 10
                    mpx.bot.signal*= calib_coeff_b/mpx.bot.gainC
                    
                    #mpx.bot.signal=tmpdata[:,II]*single(calib_coeff_b/mpx.bot.gainC)
                    #
                    # Some special cases
                    #
                    if shot>=30759 and shot<31446:
                        mpx.bot.signal[:,23]=0
                        print('Calibration problem for channel 24.') 
                    
                
            
            #        
            # merge top and bottom data
            # the bottom data is therefore rescaled by R(2kV)/2=0.2848 according to
            # coeff.m
            #
            if detec in ['top','both']:
                mpx.tvec=mpx.top.tvec
            else:
                mpx.tvec=mpx.bot.tvec
            
            mpx.dim2=chords
            if detec == 'both':
                mpx.signal=hstack((mpx.top.signal[:,chords[chords<definitions.top_nchords]], 
                            mpx.bot.signal[:,chords[chords>definitions.top_nchords]-definitions.top_nchords-1]))
                mpx.saturated=hstack((mpx.top.saturated[:,chords[chords<definitions.top_nchords]], 
                            mpx.bot.saturated[:,chords[chords>definitions.top_nchords]-definitions.top_nchords-1]))
                
                
                
            elif detec == 'top':
                mpx.signal=mpx.top.signal[:,chords[chords<definitions.top_nchords]]
                mpx.saturated=mpx.top.saturated[:,chords[chords<definitions.top_nchords]]

            else:
                mpx.signal=mpx.bot.signal[:,chords[chords>definitions.top_nchords]-definitions.top_nchords-1]
                mpx.saturated=mpx.bot.saturated[:,chords[chords>definitions.top_nchords]-definitions.top_nchords-1]


            #pmds.mdsclose('tcv_shot',shot)
            connect.closeTree('tcv_shot',shot)
            print('- MDS disconnect.\n')
            # mdsdisconnect
            #####################################################
            #       Load the detectors high voltage value       #   
            #####################################################
        
        elif act == 'v':
            
     
            
            #print('\ncase v\n')
            #pmds.mdsopen('tcv_shot',shot)
            connect.openTree('tcv_shot',shot)
            #pmds.mdsvalue('reset_public()')
            print('- Load the detectors high voltage value\n')
            
            #After #26575, the real voltage is indicated in the Vista window, before it was the reference voltage
            mm = 500 if shot<26765 else 1
    
                
            if detec in ['top','both']:
            #if strcmp(detec,'top') || strcmp(detec,'both'),
                #mpx.top.voltage=pmds.mdsvalue(r'\VSYSTEM::TCV_PUBLICDB_R["ILOT:DAC_V_01"]')*mm
                mpx.top.voltage=float(connect.get(r'\VSYSTEM::TCV_PUBLICDB_R["ILOT:DAC_V_01"]'))*mm

                if shot==32035: #special case were mpx voltage was switched of before the acquisition was finished
                    mpx.top.voltage=2000
                
                mpx.top.gainC=zeros(definitions.top_nchords)
                mpx.top.gainR=zeros(definitions.top_nchords)

                # Better than "if...else...end": the array can grow
                I=where(shot>=r_[20030,23323,26555,27127,29921,30759,31446])[0][-1] 
                if I == 0:
                    print('Detector gain dependence on the high voltage value not included in the signal calibration')
                    calib = loadmat(diag_path+'mpx_calib_first.mat')
                    calib_coeff_t=squeeze(calib['calib_coeff'])
                    mpx.top.gainC[:] = 1
                if I == 1:
                    print('Detector gain dependence on the high voltage value not included in the signal calibration')
                    calib = loadmat(diag_path+'mpx_calib_sept02.mat')
                    calib_coeff_t=squeeze(calib['calib_coeff'])
                    mpx.top.gainC[:]=1
                if I == 2:
                    print('Detector gain dependence on the high voltage value not included in the signal calibration')
                    warning('There were leaks in the top detector wire chamber for 26554<shot<27128')
                    print('Calibration is not very meaningful')

                    calib = loadmat(diag_path+'mpx_calib_may04.mat')
                    calib_coeff_t=squeeze(calib['calib_coeff'])
                    mpx.top.gainC[:]=1
                if I == 3:
                    print('Same gain dependence on the high voltage value taken for each channel')
                    
                    calib = loadmat(diag_path+'mpx_calib_july04.mat')
                    calib_coeff = squeeze(calib['calib_coeff'])
                    R = squeeze(calib['R'])

                    calib_coeff_t=calib_coeff
                    calib = loadmat(diag_path+'mpx_calib_may05.mat')
                    C = squeeze(calib['C'])
                    V = squeeze(calib['V'])

                    C=mean(C[:,:definitions.top_nchords],1) # Use the same gain for each channel
                    mpx.top.gainC[:] = exp(interp(mpx.top.voltage,V,log(C)))
                    mpx.top.gainR[:]=R
                        
                if I == 4:
                    print('Same gain dependence on the high voltage value taken for each channel')
                    calib = loadmat(diag_path+'mpx_calib_may05.mat')
                    calib_coeff = squeeze(calib['calib_coeff'])
                    C = squeeze(calib['C'])
                    V = squeeze(calib['V'])

                    calib_coeff_t=calib_coeff
                    C=mean(C[:,:definitions.top_nchords],1) # Use the same gain for each channel
                    mpx.top.gainC[:] =exp(interp(mpx.top.voltage, V,log(C)))
                    
                    R = loadmat(diag_path+'mpx_calib_july04.mat')['R'] #use the previous relative calibration
                    mpx.top.gainR[:]=R
                if I == 5:
                    #
                    # In this case, the different behaviour of the wires is contained in the matrix of gains.
                    # The calibration coefficients are in a vector: one value per wire, same value for all tensions.
                    #
                    print('Gain dependence on the high voltage value calibrated for each channel')
                    print('Leaks in the bottom detector, no relative calibration of the two detectors')
                    
                    calib = loadmat(diag_path+'mpx_calib_oct05.mat')
                    calib_coeff = squeeze(calib['calib_coeff'])
                    C = squeeze(calib['C'])
                    V = squeeze(calib['V'])
                    
                    calib_coeff_t=calib_coeff
  
                    # Interpolation to get the proper gains wrt to the high tension value
                    mpx.top.gainC[:]=[interp(mpx.top.voltage, V,log(C[:,jj])) for jj in range(definitions.top_nchords)]
                    mpx.top.gainR[:]=NaN
                if I == 6:
                    #
                    # In this case, the different behaviour of the wires is contained in the matrix of calibration coefficients.
                    # The gains are in a vector: one value per tension, same value for all wires.
                    #
                    print('Gain dependence on the high voltage value calibrated for each channel')
                    #load(sprintf('#mpx_calib_dec05_bis.mat',[mpxpath 'calibration_used/']),'calib_coeff_top','C_top_av','V_top','R')
                    calib = loadmat(diag_path+'mpx_calib_dec05_bis.mat')
                    calib_coeff_top = squeeze(calib['calib_coeff_top'])
                    C_top_av = squeeze(calib['C_top_av'])
                    V_top = squeeze(calib['V_top'])
                    R = squeeze(calib['R'])
                    calib_coeff_t = []
                    for jj in range(definitions.top_nchords): # Interpolation to get the proper calibration coefficient wrt the high tension value
                        calib_coeff_t.append(interp(mpx.top.voltage, V_top,calib_coeff_top[:,jj]))
   
                    mpx.top.gainC[:]=exp(interp(mpx.top.voltage, V_top,log(C_top_av)))
                    mpx.top.gainR =R               
                
            if detec in ['bot', 'both']:
                #mpx.bot.voltage=pmds.mdsvalue('\VSYSTEM::TCV_PUBLICDB_R["ILOT:DAC_V_02"]')*mm
                mpx.bot.voltage=float(connect.get('\VSYSTEM::TCV_PUBLICDB_R["ILOT:DAC_V_02"]'))*mm

                if shot==32035: #special case were mpx voltage was switched of before the acquisition was finished
                    mpx.top.voltage=2000
                
                mpx.bot.gainC=zeros(definitions.bot_nchords)
                mpx.bot.gainR=zeros(definitions.top_nchords)

                I=where(shot>=r_[26555,27127,29921,30759,31446])[0][-1]    
                if I == 0:
                    print('Detector gain dependence on the high voltage value not included in the signal calibration')
                    print('There were leaks in the bottom detector wire chamber for 26554<shot<27128')
                    print('Calibration is not very meaningful')
                    calib = loadmat(diag_path+'mpx_calib_may04.mat')
                    calib_coeff = squeeze(calib['calib_coeff'])
              
                    #load(sprintf('#mpx_calib_may04.mat',[mpxpath 'calibration_used/']),'calib_coeff')
                    calib_coeff_b=calib_coeff[definitions.top_nchords:definitions.top_nchords+definitions.bot_nchords]
                    mpx.bot.gainC[:]=1
                if I == 1:
                    print('Same gain dependence on the high voltage value taken for each channel')
                    
                    calib = loadmat(diag_path+'mpx_calib_july04.mat')
                    calib_coeff = squeeze(calib['calib_coeff'])
                    R = squeeze(calib['R'])

                    #load(sprintf('#mpx_calib_july04.mat',[mpxpath 'calibration_used/']),'calib_coeff','R') 
                    calib_coeff_b=calib_coeff[definitions.top_nchords:definitions.top_nchords+definitions.bot_nchords]
                    #load(sprintf('#mpx_calib_may05.mat',[mpxpath 'calibration_used/']),'C','V')
                    calib = loadmat(diag_path+'mpx_calib_may05.mat')
                    C = squeeze(calib['C'])
                    V = squeeze(calib['V'])

                    C=mean(C[:,definitions.top_nchords:definitions.top_nchords+definitions.bot_nchords],1) #use the same gain for each channel
                    #tmp=interpos(13,V,log(C),[mpx.bot.voltage mpx.bot.voltage])
                    #mpx.bot.gainC(:)=exp(tmp(1))
                    mpx.bot.gainC[:]=exp(interp(mpx.bot.voltage, V,log(C)))
                    mpx.bot.gainR[:]=R
                if I == 2:
                    print('Same gain dependence on the high voltage value taken for each channel')
                    
                    calib = loadmat(diag_path+'mpx_calib_may05.mat')
                    calib_coeff = squeeze(calib['calib_coeff'])
                    C = squeeze(calib['C'])
                    V = squeeze(calib['V'])
                    #load(sprintf('#mpx_calib_may05.mat',[mpxpath 'calibration_used/']),'calib_coeff','C','V')
                    calib_coeff_b=calib_coeff[definitions.top_nchords:definitions.top_nchords+definitions.bot_nchords]
                    C=mean(C[:,definitions.top_nchords:definitions.top_nchords+definitions.bot_nchords],1) #use the same gain for each channel
                    #tmp=interpos(13,V,log(C),[mpx.bot.voltage mpx.bot.voltage])
                    #mpx.bot.gainC(:)=exp(tmp(1))
                    mpx.bot.gainC[:]=exp(interp(mpx.bot.voltage, V,log(C)))
                    R = loadmat(diag_path+'mpx_calib_july04.mat')['R'] #use the previous relative calibration
                    mpx.bot.gainR[:]=R
                    
                if I == 3:
                    #
                    # In this case, the different behaviour of the wires is contained in the matrix of gains.
                    # The calibration coefficients are in a vector: one value per wire, same value for all tensions.
                    #
                    print('Gain dependence on the high voltage value calibrated for each channel')
                    print('Leaks in the bottom detector, no relative calibration of the two detectors')
                    calib = loadmat(diag_path+'mpx_calib_oct05.mat')
                    calib_coeff =squeeze( calib['calib_coeff'])
                    C = squeeze(calib['C'])
                    V = squeeze(calib['V'])

                    #load(sprintf('#mpx_calib_oct05.mat',[mpxpath 'calibration_used/']),'calib_coeff','C','V')
                    calib_coeff_b=calib_coeff[definitions.top_nchords:definitions.top_nchords+definitions.bot_nchords]
                    for jj in arange(definitions.bot_nchords):
                        mpx.bot.gainC[jj] = exp(interp(mpx.bot.voltage, V,log(C[:,jj+definitions.top_nchords])))
                    
                    mpx.bot.gainR[:]=NaN
                if I == 4:
                    #
                    # In this case, the different behaviour of the wires is contained in the matrix of calibration coefficients.
                    # The gains are in a vector: one value per tension, same value for all wires.
                    #
                    
               
                    print('Gain dependence on the high voltage value calibrated for each channel')
                    #load(sprintf('#mpx_calib_dec05_bis.mat',[mpxpath 'calibration_used/']),'calib_coeff_bot','C_bot_av','V_bot','R')
                    calib = loadmat(diag_path+'mpx_calib_dec05_bis.mat')
                    calib_coeff_bot = squeeze(calib['calib_coeff_bot'])
                    C_bot_av = squeeze(calib['C_bot_av'])
                    V_bot = squeeze(calib['V_bot'])
                    R = squeeze(calib['R'])

                    calib_coeff_b = []
                    for jj in range(definitions.bot_nchords):
                            #tmp=interpos(13,V_bot,calib_coeff_bot(:,jj),[mpx.bot.voltage mpx.bot.voltage])
                            #calib_coeff_b(jj)=tmp(1)
                        calib_coeff_b.append(interp(mpx.bot.voltage, V_bot,calib_coeff_bot[:,jj]))
                    
                    #tmp=interpos(13,V_bot,log(C_bot_av),[mpx.bot.voltage mpx.bot.voltage],1e6)
                    #mpx.bot.gainC(:)=exp(tmp(1))    
                    mpx.bot.gainC[:] = exp(interp(mpx.bot.voltage, V_bot,log(C_bot_av)))
                    mpx.bot.gainR[:]=NaN
                
            calib_coeff_t = asarray(calib_coeff_t)
            calib_coeff_b = asarray(calib_coeff_b)

            #pmds.mdsclose('tcv_shot',shot)
            connect.closeTree('tcv_shot',shot)

            # print('- MDS disconnect.\n')
            # mdsdisconnect
            ###############################################
            #       Load the detectors geometry [m]       #
            ###############################################
        #case 'g'
        elif act == 'g':
            #print('\ncase g\n')
            print('- Load the detectors geometry.\n')
            I=where(shot>=r_[20030,23806,26555])[0][-1] # Better than "if...else...end": the array can grow
            mpx.top.geom = {}
            mpx.bot.geom = {}

            if I == 0:
                print('MPX geometry from shot 20030 to 23806')
                mpx.top.geom['wiresR'] = r_[0.943:0.815: -0.002]
                mpx.top.geom['wiresZ'] = ones(64)*-1.558
                mpx.top.geom['length'] = 32e-3
                mpx.top.geom['slitR'] =  0.88
                mpx.top.geom['slitZ'] =  -1.
                mpx.top.geom['slit_dR'] = 2e-3
                mpx.top.geom['slit_len'] = 40e-3

                #mpx.top.geom.wires.R=r_[0.943:0.815: -0.002]
                #mpx.top.geom.wires.Z=ones(64)*-1.558
                #mpx.top.geom.wires.length=32e-3 #check it !
                #mpx.top.geom.slit.R=0.88
                #mpx.top.geom.slit.Z=-1.1 
                #mpx.top.geom.slit.dR=2e-3 #slit width
                #mpx.top.geom.slit.length=40e-3 #slit length (toroidal extension)
            if I == 1:
                print('MPX geometry from shot 23807 to 26554')
                
                mpx.top.geom['wiresR'] = r_[0.943:0.815: -0.002]
                mpx.top.geom['wiresZ'] = ones(64)*-1.14921
                mpx.top.geom['length'] = 32e-3
                mpx.top.geom['slitR'] =  0.88
                mpx.top.geom['slitZ'] =  -0.92261
                mpx.top.geom['slit_dR'] = 2e-3 #slit width
                mpx.top.geom['slit_len'] = 32e-3  #slit length (toroidal extension)
    
            if I == 2:
                print('DMPX geometry since shot 26555')
                
                mpx.top.geom['wiresR'] = r_[0.943:0.815: -0.002]
                mpx.top.geom['wiresZ'] = ones(64)*-1.15621
                mpx.top.geom['length'] = 32e-3
                mpx.top.geom['slitR'] =  0.88
                mpx.top.geom['slitZ'] =  -0.92261
                mpx.top.geom['slit_dR'] = 2e-3 #slit width
                mpx.top.geom['slit_len'] = 32e-3  #slit length (toroidal extension)
    
    
                mpx.bot.geom['wiresR'] = r_[0.942:0.816: -0.004]
                mpx.bot.geom['wiresZ'] = ones(32)*-1.17421
                mpx.bot.geom['length'] = 32e-3
                mpx.bot.geom['slitR'] =  0.88
                mpx.bot.geom['slitZ'] =  -0.92261
                mpx.bot.geom['slit_dR'] = 2e-3 #slit width
                mpx.bot.geom['slit_len'] = 32e-3  #slit length (toroidal extension)

            
            #
            # Create the dcd object
            #
            #nd=1000
            #if  detec in ['top', 'both']:
                #rs=repmat(mpx.top.geom.slit.R,1,definitions.top_nchords) #slit position R
                #zs=repmat(mpx.top.geom.slit.Z,1,definitions.top_nchords) #slit position Z
                #rd=mpx.top.geom.wires.R #wire position R
                #zd=mpx.top.geom.wires.Z #wire position Z
                #phid=zeros(definitions.top_nchords) #toroidal location angle of wires
                #tvd=zeros(definitions.top_nchords) #toroidal angles of lines of sight
                #pvd=atan((zs-zd)./(rs-rd)) #poloidal angles of lines of sight
                #pvd(pvd<0)=pi+pvd(pvd<0)
                #pvd=pi-pvd #angle between the horizontal plane and the chord, sens trigonometrique
                #mpx.top.geom.dcd=psitbxdcd(rd,zd,phid,pvd,tvd,nd)
                #if (strcmp(detec,'top'))
                    #rd=rd(chords)
                    #zd=zd(chords)
                    #phid=phid(chords)
                    #pvd=pvd(chords)
                    #tvd=tvd(chords)
                #end
            #end
            
            #if detec in ['bot', 'both']:

                #rs=repmat(mpx.bot.geom.slit.R,1,definitions.bot_nchords)
                #zs=repmat(mpx.bot.geom.slit.Z,1,definitions.bot_nchords)
                #rd=mpx.bot.geom.wires.R
                #zd=mpx.bot.geom.wires.Z
                #phid=zeros(1,definitions.bot_nchords)
                #tvd=zeros(1,definitions.bot_nchords)
                #pvd=atan((zs-zd)./(rs-rd))
                #pvd(pvd<0)=pi+pvd(pvd<0)
                #pvd=pi-pvd
                #mpx.bot.geom.dcd=psitbxdcd(rd,zd,phid,pvd,tvd,nd)
                #if detec == 'bot'
                    #rd=rd(chords-definitions.top_nchords)
                    #zd=zd(chords-definitions.top_nchords)
                    #phid=phid(chords-definitions.top_nchords)
                    #pvd=pvd(chords-definitions.top_nchords)
                    #tvd=tvd(chords-definitions.top_nchords)
                
            
            #if(strcmp(detec,'both'))
                    #rs=[repmat(mpx.top.geom.slit.R,1,definitions.top_nchords) repmat(mpx.bot.geom.slit.R,1,definitions.bot_nchords)]
                    #zs=[repmat(mpx.top.geom.slit.Z,1,definitions.top_nchords) repmat(mpx.bot.geom.slit.Z,1,definitions.bot_nchords)]
                #rd=[mpx.top.geom.wires.R mpx.bot.geom.wires.R]
                    #zd=[mpx.top.geom.wires.Z mpx.bot.geom.wires.Z]
                    #phid=zeros(1,definitions.top_nchords+definitions.bot_nchords)
                    #tvd=zeros(1,definitions.top_nchords+definitions.bot_nchords)
                    #pvd=atan((zs-zd)./(rs-rd))
                    #pvd(pvd<0)=pi+pvd(pvd<0)
                    #pvd=pi-pvd
                ## restrict to selected chords
                #rd=rd(chords)
                #zd=zd(chords)
                #phid=phid(chords)
                #pvd=pvd(chords)
                #tvd=tvd(chords)
            #end
            #mpx.geom.dcd=psitbxdcd(rd,zd,phid,pvd,tvd,nd)


            #####################################
            #       Load the DMPX filters       #
            #####################################
        #case 'f'
        elif act == 'f':
            #print('\ncase f\n')
            #pmds.mdsopen('tcv_shot',shot)
            connect.openTree('tcv_shot',shot)
            #pmds.mdsvalue('reset_public()')
            print('- Load the MPX filters and detection gas.\n')
            mpx.top.filters = {}
            mpx.bot.filters = {}
            #abso = []

            if shot<26765:
                print('No mobile absorber holder for shot<26765')
                mpx.top.filters['mobile_holder_pos']=NaN
                mpx.bot.filters['mobile_holder_pos']=NaN
                I_pos=4
            else:
                if 'pos' in kwarg and pos in [1,2,3,4]:
                    print(( 'Mobile absorber position forced to '+str(pos)))
                    I_pos=pos
                else:
                    print('Load the mobile absorber holder position')
                    
                    #print 'BUGGGG!!'
                    
                    #a = [pmds.mdsvalue('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A%d"]'%i).strip() == 'OFF' for i in  arange(4)]
                    a = [connect.get('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A%d"]'%i).strip() == 'OFF' for i in  arange(4)]

                    
                                #connect.openTree('tcv_shot',shot)

                    #a0=pmds.mdsvalue('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A0"]').strip() == 'OFF'

                    #a0=strcmp(deblank(pmds.mdsvalue('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A0"]')),'OFF')
                    #a1=strcmp(deblank(pmds.mdsvalue('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A1"]')),'OFF')
                    #a2=strcmp(deblank(pmds.mdsvalue('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A2"]')),'OFF')
                    #a3=strcmp(deblank(pmds.mdsvalue('\VSYSTEM::TCV_PUBLICDB_B["ILOT_NEB:A3"]')),'OFF')
                    I_pos=4*a[0]+3*a[1]+2*a[2]+1*a[3]
                
                mpx.top.filters['mobile_holder_pos'] = I_pos
                mpx.bot.filters['mobile_holder_pos'] = I_pos
                abso.material=abso.material[I_pos-1]
                abso.thickness=abso.thickness[I_pos-1]
            
            
            mpx.top.gas = {}
            mpx.bot.gas = {}

            I=where(shot>=r_[20030,21545,22043,23578,31942,32039])[0][-1] #gas in the wire chamber
            if I == 0:    
                
                mpx.top.gas['type']='ArCH'
                mpx.top.gas['thickness']=6e-3
            if I == 1:     
                mpx.top.gas['type']='KrCH'
                mpx.top.gas['thickness']=6e-3
            if I == 2:     
                mpx.top.gas['type']='XeCH'
                mpx.top.gas['thickness']=6e-3
            if I in  [3,5]:
                if shot<26555:
                    mpx.top.gas['type']='KrCH'
                    mpx.top.gas['thickness']=6e-3
                else:
                    if  detec in ['top','both']:
                        mpx.top.gas['type']='KrCH'
                        mpx.top.gas['thickness']=8e-3
                    
                    if  detec in ['bot','both']:
                        mpx.bot.gas['type']='KrCH'
                        mpx.bot.gas['thickness']=7.6e-3 # bottom detector is slightly thiner (construction flaw...)
                    
            if I == 4:
                if  detec in ['top','both']:
                    mpx.top.gas['type']='ArCH'
                    mpx.top.gas['thickness']=8e-3
                
                if  detec in ['bot','both']:
                    mpx.bot.gas['type']='ArCH'
                    mpx.bot.gas['thickness']=7.6e-3 # bottom detector is slightly thiner (construction flaw...)
                 
            
            I=where(shot>=r_[20030,23806,26555])[0][-1]
            if I == 0:
                mpx.top.filters['material']=['He','Be']
                mpx.top.filters['thickness']=[526.8,350e-6]
                if shot>=21946 and shot<22085: #Al absorber inserted manually in the He tube
                        mpx.top.filters['material'].append('Al')
                        mpx.top.filters['thickness'].append(308e-6)
                 
            if I == 1:
                mpx.top.filters['material']=['He','Be']
                mpx.top.filters['thickness']=[223.6e-3 ,150e-6]     
            if I == 2:
                if detec in ['top','both']: 
                    mpx.top.filters['material']=['He','Be']
                    mpx.top.filters['thickness']=[229.6e-3,100e-6]
                    #
                    # Add the filter of the mobile absorber
                    #
                    if abso.material != 'nothing':
                        I = None
                        for iii in range(len(mpx.top.filters['material'])):
                            if abso.material == mpx.top.filters['material'][iii]:
                                I=iii;break
                            
                        
                        if not I is None:
                            mpx.top.filters['thickness'][I]=mpx.top.filters['thickness'][I]+abso.thickness
                        else:
                            mpx.top.filters['material'].append(abso.material)
                            mpx.top.filters['thickness'].append(abso.thickness)
                
                
                if detec in ['bot','both']: 
                    mpx.bot.filters['material']=['He','Be','Air',mpx.top.gas['type']]
                    mpx.bot.filters['thickness']=[229.6e-3,200e-6,10e-3,mpx.top.gas['thickness']]
                    #
                    # Add the filter of the mobile absorber
                    #
                    if abso.material!= 'nothing':
                        I=None
                        for iii in range(len(mpx.bot.filters['material'])):
                            if abso.material == mpx.bot.filters['material'][iii]:
                                I=iii;break
                            
                        
                        if not I is None:
                            mpx.bot.filters['thickness'][I]=mpx.bot.filters['thickness'][I]+abso.thickness
                        else:
                            mpx.bot.filters['material'].append(abso.material)
                            mpx.bot.filters['thickness'].append(abso.thickness)
                        
                    
                    #
                    # Add the additional absorber between the two detectors (if present)
                    #
                    # Insert shots by pair, additional absorber betweem the two shots numbers
                    I=where(shot>=r_[28751,28834,28862,28891,28943,28973,30681,30690])[0][-1]
                    addi_material=['Al','Al','Al','Al']
                    addi_thick=r_[308,308,308]*1e-6
                    
                    #import IPython
                    #IPython.embed()
                    
                    if I%2 == 0: # Means there was an additional absorber between the 2 detectors
                        II=None
                        for iii in range(len(mpx.bot.filters['material'])):
                            if addi_material[I%2] == mpx.bot.filters['material'][iii]:
                                II=iii;break
                            
                        
                        if not II is None:
                            mpx.bot.filters['thickness'][II]=mpx.bot.filters['thickness'][II]+addi_thick[I%2]
                        else:
                            mpx.bot.filters['material'].append(addi_material[I%2])
                            mpx.bot.filters['thickness'].append(addi_thick[I%2])
     
     
            #pmds.mdsclose('tcv_shot',shot)
            connect.closeTree('tcv_shot',shot)
            print('- MDS disconnect.\n')
            # mdsdisconnect
            #############################################
            #       Load the detectors efficiency       #
            #############################################

        elif act == 'e':
            #print('\ncase e\n')
            
            if detec in ['bot','both']: 
                print('- Load the top detector efficiency.\n')
                st=''
                for iii in range(len(mpx.top.filters['material'])):
                    st+= mpx.top.filters['material'][iii]+str(mpx.top.filters['thickness'][iii]*1e6)+', '
                st+= ' $'
                st+= ' '+mpx.top.gas['type']+'abs' +str(mpx.top.gas['thickness']*1e6)

   
                out = XraySignals(st)
                
                out.response[isnan(out.response)]=0
                mpx.top.detector = st
                mpx.top.eff_ev=out.ev
                mpx.top.eff_response=out.response
            
            
            
            if detec in ['bot','both']: 
                print('- Load the bottom detector efficiency.\n')
                sb=''
                for iii in range(len(mpx.bot.filters['material'])):
                    sb+= mpx.bot.filters['material'][iii]+str(mpx.bot.filters['thickness'][iii]*1e6)+', '
                sb+= ' $'


                sb+= ' '+mpx.bot.gas['type']+'abs' +str(mpx.bot.gas['thickness']*1e6)

   
                out=XraySignals(sb)
                out.response[isnan(out.response)]=0
                mpx.bot.eff_ev=out.ev
                mpx.bot.eff_response=out.response       
                mpx.bot.detector = sb

           
        else:
            raise Exception('mpxdata:UnknownAction  '+act )

        
        #error('mpxdata:UnknownAction',['''' act '''' ' is an unknown action'])
         # End of 'for' loop
     # End of 'switch' cases
    #print(['Closing MDS connection to shot ' num2str(shot) '.'])
    try:
        #pmds.mdsclose('tcv_shot',shot)
        connect.closeTree('tcv_shot',shot)
    except:
        pass
    #
    
    #import IPython
    #IPython.embed()

    print('\n+ MPXDATA finished\n\n')
    
    
    return mpx
    





def XraySignals(det):

    #########################################################################
    #
    #  Function which calculates the efficiency (percentage of incoming energy detected) of a 
    #  detector, as a function of the energy of the incoming photons
    #  
    #  SYNTAX
    #
    #  out=XraySignals(det)
    #
    #  INPUT
    #
    #  det  Detector string:
    #       Everything before $ symbol is a filter, everything after
    #       is an absorber. Material symbol (case is not important) 
    #       is followed by thickness in MICRONS.
    #       
    #       LIST of materials for filters:
    #       Be = beryllium
    #       Al = aluminium
    #       Si = silicium
    #       He = helium
    #       KrCH = 90# krypton + 10# methan
    #       ArCH = 90# argon + 10# methan
    #       XeCH = 90# xenon + 10# methan
    #       AIR = air (mixture 78.1# N, 21# O, 0.015# C, 0.47# Ar)
    #       
    #       LIST of materials for absorbers:
    #       
    #       KrCHabs = 90# krypton + 10# methan
    #       ArCHabs = 90# argon + 10# methan
    #       XeCHabs = 90# xenon + 10# methan
    #
    #  EXAMPLES
    #  
    #  * two filters system with 5000 micron Ar detector  
    #  det = 'al200, Be200 $ Ar5000'
    #        
    #  * DMPX top detector with KrCH absorber
    #  det = 'He229600, Be100 $ KrCHabs8000'
    #        
    #  * DMPX bottom detector with KrCH absorber
    #  det = 'He229600, Be200, KrCH8000, AIR10000 $ KrCHabs7600'
    #
    #
    #  OUTPUT
    #   
    #  out.response=out.tr*out.abso        system efficiency 
    #  out.abso                             absorber efficency
    #  out.tr                               filter transmision 
    #  out.ev                               energy vector used
    # 
    #  From Alexei Zabolotsky routine
    #########################################################################


    filt,abso,error=XraySignals_Strings(det);
        
    if error==1:
        print('XraySignals ERROR: cannot process detectors string') 
        return
           
        
    E=linspace(200,2e5,3000);   # energy vector (eV)

#***************  get filter transmission characteristics *********************
    tr = [] 
    for F in filt:
        [tx,error]=XraySignals_SxCurves(F['material'],F['thick'],E)
        tr.append(tx)      # each filter transmission as a function of frequency
        
    tr=prod(tr,0)  # total filter transmission as a function of frequency

#***************  get detector absorbtion characteristics *********************
    diode_eff = []
    for A in abso:
        [tx,error]=XraySignals_SxCurves(A['material'],A['thick'],E);
        diode_eff.append(1-tx)  # each diode efficency as a function of frequency
    
    diode_eff=prod(diode_eff,0) # total diode efficency as a function of frequency


#****************************** Form output ***************************************    
    class out_: pass
    out = out_()
    out.abso=diode_eff;          # diode efficency
    out.tr=tr;                    # filter transmision 
    out.ev=E;                     # energy vector used
    out.response=tr*diode_eff;   # system energy respond 
    #import IPython
    #IPython.embed()
        
    return out





def XraySignals_Strings(det):

    ##############################################################
    # 
    #  Analyse string contaning the information about the XRAY 
    #  diods. Prepares data for XraySignals.m
    # 
    #  SYNTAX
    #   function [filt,abso,error]=XraySignals_Strings(det);
    #     
    #  INPUT string
    #    EXAMPLE det=['Be50,  Si0.7, O1500, N4500 $ Be200']
    #
    #from Alexei Zabolotsky routine
    #########################################################################
    error=0;
    det = det.replace(' ','') # remove blanc spaces  
         
    ind=det.find('$')
    
    if ind == -1:
        raise Exception('Strings ERROR: No detector found in det string')   

    if sum([d=='$' for d in det]) > 1:
        raise Exception('Strings ERROR: More than one detector found in det string')   
        
    
    filtersstr=det[:ind-1] # take only symbols before "$" sign
    detectorstr=det[ind+1:]
    

    #*************************** Translate detector string ************************************* 


    detectors = detectorstr.split(',')
    abso = [{} for d in detectors]


    print('      ****************  Absorbers: ********************')
    for i, det in enumerate(detectors):
        #letters = ~array([d.isdigit() for d in det])
        #tmp2=detectorstr(index(i)+1:index(i+1)-1);
        abso[i]['thick']=float(''.join([d for d in det if not d.isalpha()])) # take numbers as thickness 
        abso[i]['material']=''.join([d for d in det if d.isalpha()])   # take letters as absorber
        print('      *         '+abso[i]['material']+'  thickness '+str(abso[i]['thick'])+'  micron')
       

    #*************************** Translate filter string ************************************* 
    filters = filtersstr.split(',')
    filt = [{} for f in filters]
    print('      ****************  Filters: *********************')

    for i, fil in enumerate(filters):
        #letters = ~array([d.isdigit() for d in fil])
        #tmp2=detectorstr(index(i)+1:index(i+1)-1);
        filt[i]['thick']=float(''.join([d for d in fil if not d.isalpha()]))  # take numbers as thickness 
        filt[i]['material']=''.join([d for d in fil if d.isalpha()])   # take letters as absorber
        print('      *         '+filt[i]['material']+'  thickness '+str(filt[i]['thick'])+'  micron')
       

    

    #index=[0,findstr(filtersstr,','),length(filtersstr)+1]; 
    #nfilters=length(findstr(filtersstr,','))+1;       # number of filters
    #for i=1:nfilters
        #tmp2=filtersstr(index(i)+1:index(i+1)-1);
        #filt(i).thick=str2num(tmp2(~isletter(tmp2)));  # take numbers as thickness 
        #filt(i).material=tmp2(isletter(tmp2));             # take letters as absorber
        #print(['      *         ',filt(i).material,'  thickness ',int2str(filt(i).thick),'  micron'])
       
    print('      ************************************************')
    return filt,abso,error
 


def XraySignals_SxCurves(abso,thick,e):

    error=0
    tr=[];
    thick=thick/1e4;  # from microns to cm

    abso=abso.upper()

    
    directory = os.path.dirname(os.path.realpath(__file__))+'/DATA/'
          
    #BUG rho is density at atmospheric pressure and room temperature!!! why?? 
    try:
        if abso=='BE':
            ev,murho=loadtxt(directory+'be.dat',unpack=True);     

            # mass attenuation coefficient (cm2 g-1)
            rho=1.845;            # Nominal density:    (g cm-3)
            tr=exp(-(murho*rho*thick))
            
        elif abso=='AL':
            #out=loadmat(directory+'al.dat'); 
            ev,murho=loadtxt(directory+'al.dat',unpack=True);     
            rho=2.6941;
            tr=exp(-(murho*rho*thick));
            
        elif abso=='SI':
            #out=loadmat(directory +'si.dat'); 
            ev,murho=loadtxt(directory+'si.dat',unpack=True);     

            rho=2.3200;
            tr=exp(-(murho*rho*thick));            
            
        elif abso=='HE':
            #out=loadmat(directory +'he.dat');    
            ev,murho=loadtxt(directory+'he.dat',unpack=True);     
            rho=1.6640E-04;
            tr=exp(-(murho*rho*thick)); 
            
        elif abso=='KRCH':
            #out=loadmat(directory+ 'kr.dat');  
            ev,murhokr=loadtxt(directory+'kr.dat',unpack=True);     

            rhokr=3.4840E-03;
            murhokr=interp(e,ev*1000,murhokr);

            #out=loadmat(directory +'ch4.dat');
            #ev=(out(:,1));
            ev,murhoch4=loadtxt(directory+'ch4.dat',unpack=True);     

            rhoch4=0.424;
            #murhoch4=(out(:,2));
            murhoch4=interp(e,ev*1000,murhoch4);

            #Mixture 90# Kr and 10# CH4 
            tr=exp(-(murhoch4*rhoch4*thick*0.1+murhokr*rhokr*thick*0.9));     #only KR generates photoelectron detected by MPX wires
        
        elif abso=='KRCHABS':
            ##out=loadmat(directory +'kr.dat');  
            ev,murhokr=loadtxt(directory+'kr.dat',unpack=True);     

            #ev=(out(:,1));
            rhokr=3.4840E-03;
            #murhokr=(out(:,2));
            murhokr=interp(e,ev*1000,murhokr);

            #Mixture 90# Kr and 10# CH4 
            tr=exp(-(murhokr*rhokr*thick*0.9));     
            
        elif abso=='ARCH':
            #out=loadmat(directory +'ar.dat');  
            ev,murhoar=loadtxt(directory+'ar.dat',unpack=True);     

            #ev=(out(:,1));
            
            rhoar=3.4840E-03;
            #murhoar=(out(:,2));
            murhoar=interp(e,ev*1000,murhoar);

            #out=loadmat(directory+ 'ch4.dat');
            #ev=(out(:,1));
            ev,murhoch4=loadtxt(directory+'ch4.dat',unpack=True);     

            rhoch4=0.424;
            #murhoch4=(out(:,2));
            murhoch4=interp(e,ev*1000,murhoch4);

            #Mixture 90# Ar and 10# CH4 
            tr=exp(-(murhoch4*rhoch4*thick*0.1+murhoar*rhoar*thick*0.9));                

        elif abso=='ARCHABS':
            #out=loadmat(directory +'ar.dat');     
            ev,murhoar=loadtxt(directory+'ar.dat',unpack=True);     

            #ev=(out(:,1));
            rhoar=3.4840E-03;
            #murhoar=(out(:,2));
            murhoar=interp(e,ev*1000,murhoar);

            #Mixture 90# Ar and 10# CH4 
            tr=exp(-(murhoar*rhoar*thick*0.9));                
                
        elif abso=='XECH':
            #out=loadmat(directory +'xe.dat');     
            ev,murhoxe=loadtxt(directory+'xe.dat',unpack=True);     

            #ev=(out(:,1));
            rhoxe=5.4580E-03;
            #murhoxe=(out(:,2));
            murhoxe=interp(e,ev*1000,murhoxe);

            #out=loadmat(directory +'ch4.dat');
            #ev=(out(:,1));
            ev,murhoch4=loadtxt(directory+'ch4.dat',unpack=True);     

            rhoch4=0.424;
            #murhoch4=(out(:,2));
            murhoch4=interp(e,ev*1000,murhoch4);

            #Mixture 90# Xe and 10# CH4 
            tr=exp(-(murhoch4*rhoch4*thick*0.1+murhoxe*rhoxe*thick*0.9));     
            
        elif abso=='XECHABS':
            ##out=loadmat(directory +'xe.dat');      
            ev,murhoxe=loadtxt(directory+'xe.dat',unpack=True);     

            #ev=(out(:,1));
            rhoxe=5.4580E-03;
            #murhoxe=(out(:,2));
            murhoxe=interp(e,ev*1000,murhoxe);

            #Mixture 90# Xe and 10# CH4 
            tr=exp(-(murhoxe*rhoxe*thick*0.9));     

        elif abso=='AIR'  :     
            #out=loadmat(directory +'n.dat'); 
            ev,murhon=loadtxt(directory+'n.dat',unpack=True);     

            #ev=(out(:,1));
            rhon=1.1650E-03;
            #murhon=(out(:,2));
            murhon=interp(e,ev*1000,murhon);

            #out=loadmat(directory+ 'o.dat');  
            ev,murhoo=loadtxt(directory+'o.dat',unpack=True);     

            #ev=(out(:,1));
            rhoo=1.3310E-03;
            #murhoo=(out(:,2));
            murhoo=interp(e,ev*1000,murhoo);

            #out=loadmat(directory+ 'c.dat');  
            ev,murhoc=loadtxt(directory+'c.dat',unpack=True);     

            #ev=(out(:,1));
            rhoc=2.2600;
            #murhoc=(out(:,2));
            murhoc=interp(e,ev*1000,murhoc);

            #out=loadmat(directory+ 'ar.dat');  
            ev,murhoar=loadtxt(directory+'ar.dat',unpack=True);     

            #ev=(out(:,1));
            rhoar=3.4840E-03;
            #murhoar=(out(:,2));
            murhoar=interp(e,ev*1000,murhoar);
            
            #Mixture 78.1# N, 21# O, 0.015# C, 0.47# Ar 
            tr=exp(-(murhon*rhon*thick*0.781+murhoo*rhoo*thick*0.21+\
                murhoc*rhoc*thick*0.00015+murhoar*rhoar*thick*0.0047));     
            
        else:
            raise Exception('SxCurves ERROR: Data for ' +abso+ ' does not exist')
            
    except:
        raise Exception('SxCurves ERROR: cannot load data file for '+abso)
       
        
    'HE'
    if len(tr)!= len(e):
        tr=interp(e,ev*1000,tr) 
           
    
    return tr,error
    

##import IPython
##IPython.embed()
#server =  'localhost:8003'
#pmds.mdsconnect( server)

#out = mpxdata(40412,'sgve')
##out = mpxdata(40412,'s')
#pmds.mdsdisconnect()
#import IPython
#IPython.embed()
#exit()
#out.top.geom[]
#'wiresR', 'slit_dR', 'wiresZ', 'length', 'slitR', 'slit_len', 'slitZ'



#server =  'localhost:8003'
#pmds.mdsconnect( server)



#pmds.mdsdisconnect()



#import IPython
#IPython.embed()


#mdsconnect( 'tcv1.epfl.ch')  
#[mpx] = mpxdata(40412,'sgvre'); 

#plot(mean(mpx.signal.data,1))