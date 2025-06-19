#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
    

  

import os,sys 
from collections import OrderedDict
from numpy import *
from matplotlib.pyplot import *
#from annulus import get_bd_mat
import MDSplus as mds
from scipy.signal import fftconvolve, medfilt
import time
from multiprocessing import  cpu_count
from multiprocessing.pool import Pool
from IPython import embed
import numpy

def savitzky_golay(y, window_size, order, deriv=0, rate=1,firstvals=None, lastvals=None):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = list(range(order+1))
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    if firstvals is None:
        firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    else:
        firstvals = ones(half_window)*firstvals
    if lastvals is None:
        lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
   
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval(window+'(window_len)')

    y=convolve(w/w.sum(),s,mode='valid')
    return y

"""
def mds_load(params):
    (mds_server,   TDI, tree, shot) = params
    print(TDI)
    MDSconn = mds.Connection(mds_server )
    MDSconn.openTree(tree, shot)
    output = [MDSconn.get(tdi).data() for tdi in TDI]
    MDSconn.closeTree(tree, shot)

    return output

def mds_par_load(mds_server,   TDI, tree, shot,  numTasks):

    #load a junks of a single vector

    TDI = array_split(TDI, min(numTasks, len(TDI)))

    args = [(mds_server,  tdi, tree, shot) for tdi in TDI]

    pool =  Pool(len(args))
    

    out = pool.map(mds_load,args)
    pool.close()
    pool.join()

    output = []
    for o in out:output+= o

    
    return  output

"""

def richardson_lucy(tmp):
    """\
    Performs a Richardson-Lucy deconvolution.
    """
    (img, psf, numiter,offset,relax) = tmp
    import scipy.ndimage as nd
    from scipy.signal import fftconvolve
    

    epsilon = 1e-8 
    PSF_inv = np.flipud(psf)

    psf/= sum(psf)
    img_ = maximum(img+offset,0)
    cur_est = img_.copy()
    #Wupd = 1/(1+sigma^2) sigma - detector noise
    fpsf  = fft.rfft(psf)
    fipsf = fft.rfft(PSF_inv)

    for it in range(numiter):
        fcur_est = fft.rfft(cur_est)
        img_model = fft.irfft( fcur_est*fpsf)  #convolution
        #img_model =  fftconvolve(cur_est, psf, 'full')[:-len(psf)+1]
        img_model = maximum(img_model+offset,0)+1e-8
        Wupd = 1/(1+1/sqrt(img_model/mean(img_model)))
        Wupd/=amax(Wupd)
        fcorrection = fft.rfft(1 - img_ / img_model)
        update = fft.irfft( fcorrection*fipsf)*cur_est*Wupd
        #update *= cur_est
        #update = Wupd*fftconvolve( 1 - img_ / img_model, PSF_inv, 'full')[len(PSF_inv)-1:]  * cur_est
        cur_est = maximum(0, cur_est -  relax* update);
 
    return cur_est


def medfilt_deconv(xxx_todo_changeme2):
    (data,drtau,tdel,tau,gamma, scrfact ) = xxx_todo_changeme2
    m         = max(2*int(int(drtau/tdel)/2)+1, 3)
    k         = arange(m)-((m-1)/2.)
    
    k[k != 0] = 1.0/k[k != 0]
    k        = k/tdel/(int(m/2.)*2)
    temp      = data-mean(data[:20],0)
    temp      = temp[:len(temp)//4*4].reshape(-1,4).mean(1)
    temp      = medfilt(temp,m)
    drdt      = np.convolve( temp,k, mode='same')
    drdt      = medfilt(drdt,m)
    
    return (gamma*temp-tau*drdt)/scrfact


    
class loader_BOLO():

    
    """
    class for the loading of the foil bolometer data from AUG

    """
    
    
    dt = 0.001
    sigma = 0.00
    min_error = 0.02
    
    def __init__(self, shot, geometry_path, MDSconn):
        
        """

        :var int shot:  Number of selected shot
        :var str geometry_path: Path to saved geometry, boundary or data
        :var object dd: loading library for AUG data
        :var str experiment: name of AUG  experiment
        :var int edition: edition of AUG  experiment

        """
        
        # Handle invalid shot range
        if shot <= 91000:
            raise Exception("Bolometer invalid for shots < 91000")
        
        self.tree = 'BOLOM'
        self.MDSconn = MDSconn
        self.geometry_path = geometry_path
                 
        self.shot = shot
        self.wrong_dets_damaged = []#BUG
        #self.PFM = tok.PFM
        #tok.get_psi = self.get_psi
        #tok.get_mag_contours = self.get_mag_contours

        #Relevant references are
        #/c/idl/source/efitview/diagnoses/DIII-D/bolometerpaths.pro
        #/c/4dlib/DIAGS/BOLOM/get_bolo.pro

        # Channels are [BOL_U01_P-BOL_U24_P, BOL_L01_P-BOL_L24_P]
        
        self.nchan = 24
        self.cams = 'U','L'
        


        # Radius of detectors (m)
        RU = [2.34881]*15 + [2.31206]*9
        RL = [2.31894]*11 + [2.34932]*13

        # Eleveation of detectors (m)
        ZU = [0.72817]*15 + [0.82332]*9
        ZL = [-0.77254]*11 + [-0.66881]*13

        # Toroidal angle of detectors
        phiU = [75.]*self.nchan
        phiL = [60.]*self.nchan

        angleU=[269.4, 265.6, 261.9, 258.1, 254.4, 250.9,\
        247.9, 244.9, 241.9, 238.9, 235.9, 232.9,\
        228.0, 221.3, 214.5, 208.2, 201.1, 194.0,\
        187.7, 182.2, 176.7, 171.2, 165.7, 160.2]

        angleL=[213.7, 210.2, 206.7, 203.2, 199.7, 194.4,\
        187.4, 180.4, 173.4, 166.4, 159.4, 156.0,\
        149.2, 142.4, 135.8, 129.6, 126.6, 123.6,\
        120.6, 117.6, 114.6, 111.6, 108.6, 101.9]
        
        # Create second point to define a line
        RU2 = array(RU) + 2.0 * cos(deg2rad(angleU))
        RL2 = array(RL) + 2.0 * cos(deg2rad(angleL))

        ZU2 = array(ZU) + 2.0 * sin(deg2rad(angleU))
        ZL2 = array(ZL) + 2.0 * sin(deg2rad(angleL))

        self.Phi = {'U':phiU, 'L':phiL}           
        self.R_start= {"U":RU, 'L':RL}
        self.R_end  = {"U":RU2,'L':RL2}
        self.z_start= {"U":ZU, 'L':ZL}
        self.z_end  = {"U":ZU2,'L':ZL2}
        self.theta=   {"U":rad2deg(angleU),"L":rad2deg(angleL)}
        self.etendue = {}
        self.etendue['U'] = [3.0206e4,2.9034e4,2.8066e4,2.7273e4,2.6635e4,4.0340e4,\
        3.9855e4,3.9488e4,3.9235e4,3.9091e4,3.9055e4,3.9126e4,\
        0.7972e4,0.8170e4,0.8498e4,0.7549e4,0.7129e4,0.6854e4,\
        1.1162e4,1.1070e4,1.1081e4,1.1196e4,1.1419e4,1.1761e4]

        self.etendue['L'] = [2.9321e4,2.8825e4,2.8449e4,2.8187e4,2.8033e4,0.7058e4,\
        0.7140e4,0.7334e4,0.7657e4,0.8136e4,0.8819e4,0.7112e4,\
        0.6654e4,0.6330e4,0.6123e4,2.9621e4,2.9485e4,2.9431e4,\
        2.9458e4,2.9565e4,2.9756e4,3.0032e4,3.0397e4,0.6406e4]


      
   
        self.detectors_dict = OrderedDict()
        self.cam_ind = OrderedDict()
        self.nl = 0
        self.all_los = []
        for cam in self.cams:
            self.detectors_dict[cam] = ['BOL_%s%.2d'%(cam,ch+1) for ch in range(self.nchan)]
            self.cam_ind[cam] = arange(self.nl,self.nchan+self.nl)
            self.nl+= self.nchan
            self.all_los.append(self.detectors_dict[cam])

        self.calb_0 = ones(len(self.cams))
        self.dets = arange(self.nl)
        self.dets_index = [v for k,v in self.cam_ind.items()]
        
        
      
        try:
          #  ss
            self.MDSconn.openTree(self.tree,self.shot)
            TDIcall = "\\BOLOM::TOP.PRAD_01:TIME"
            self.tvec = self.MDSconn.get(TDIcall).data()/1e3
            self.MDSconn.closeTree(self.tree,self.shot)
            self.real_time_data = False
        except:
            import time
            print('Standart bolometer data are not availible')
            for i in range(10):
                self.tvec = self.MDSconn.get(f'dim_of(PTDATA2("DGSDPWRU01",{self.shot},1))').data()/1e3  
                if len(self.tvec) > 1:
                    break
                print('No realtime bolometer data, waiting 10s...')
                time.sleep(10)

            if len(self.tvec) == 1:
                raise Exception('No Realtime bolo data!!')
            
            self.real_time_data = True
            nt = len(self.tvec)
            dt = (self.tvec[-1]-self.tvec[0]) / (nt+1)
            self.downsample = max(int(1e-2/dt), 1)
            #downsample
            self.tvec = self.tvec[:nt // self.downsample * self.downsample].reshape(-1, self.downsample).mean(1)
                
            
            

    #def get_psi(self, tvec,xgridc, ygridc):
        #from .. import map_equ
        ##import map_equ

        #eqm = map_equ.equ_map(None)
        #eqm.pfm   = self.PFM['pfm']
        #eqm.t_eq  = self.PFM['pfm_tvec']
        #eqm.Rmesh = self.PFM['pfm_R']
        #eqm.Zmesh = self.PFM['pfm_Z']
        #eqm.eq_open=True        
 
        #R,Z = meshgrid(xgridc,ygridc)
        #Psi = eqm.rz2rho( R[None] , Z[None] , tvec, coord_out='Psi')  
        
        #return squeeze(sqrt(maximum(0,Psi)))
    
    #def get_mag_contours(self,tvec,rho):
        ##import map_equ
        #from .. import map_equ

        #eqm = map_equ.equ_map(None)
        #eqm.pfm   = self.PFM['pfm']
        #eqm.t_eq  = self.PFM['pfm_tvec']
        #eqm.Rmesh = self.PFM['pfm_R']
        #eqm.Zmesh = self.PFM['pfm_Z']
        #eqm.psi0 = eqm.t_eq*nan
        #eqm.eq_open=True        
        
        #return eqm.rho2rz(rho**2, tvec,'Psi',True)  
        
        
    def get_total_rad(self):
        
         
        try:
            assert not self.real_time_data 
            try:
                powers = list(np.load(geometry_path+'/power_%d.npz'%self.shot, allow_pickle=True)['power'])
            except:
                #try:
                powers = []
                print('get_total_rad')
                signals = 'PRAD_TOT','PRAD_DIVL','PRAD_DIVU'
                self.MDSconn.openTree(self.tree,self.shot)
                for sig in signals:
                    TDIcall = "\\BOLOM::TOP.PRAD_01.PRAD:"
                    prad = self.MDSconn.get(TDIcall+sig).data()
                    self.powers.append((sig, self.tvec, prad))
                
                self.MDSconn.closeTree(self.tree,self.shot)
                for tree in ('EFIT02','EFIT01'):
                    try:
                        self.MDSconn.openTree(tree,self.shot)
                        
                        ptot = zeros_like(self.tvec)
                        TDIcall = "\\%s::TOP.RESULTS.CONFINEMENT.POWER:PNBI"%tree
                        data = self.MDSconn.get('_x ='+TDIcall).data()
                        tvec = self.MDSconn.get('dim_of(_x)').data()/1000
                        dt = mean(diff(tvec))
                        decay = 150e-3#s
                        from scipy.signal import lfilter
                        data = lfilter((dt/decay ,), [1,-1+dt/decay], data)

                        ptot+=interp(self.tvec,tvec,data )
                        
                        TDIcall = "\\%s::TOP.RESULTS.CONFINEMENT.POWER:PECH"%tree
                        data = self.MDSconn.get('_x ='+TDIcall).data()
                        tvec = self.MDSconn.get('dim_of(_x)').data()/1000
                        data = lfilter((dt/decay ,), [1,-1+dt/decay], data)
                        ptot+=interp(self.tvec,tvec,data )
                        
                        TDIcall = "\\%s::TOP.RESULTS.CONFINEMENT.POWER:POH"%tree
                        data = self.MDSconn.get('_x ='+TDIcall).data()
                        tvec = self.MDSconn.get('dim_of(_x)').data()/1000
                        ptot+=interp(self.tvec,tvec,data )
                
                        
                        self.powers.append(("%s:PTOT"%tree,self.tvec,ptot ))
                        self.MDSconn.closeTree(tree,self.shot)
                        break
                    except:
                        pass
           
                np.savez_compressed(geometry_path+'/power_%d.npz'%self.shot, power=self.powers)
        except:
            powers = {}
            
       
        return powers
      
      
      
      
    def get_data_(self,tmin=-infty,tmax=infty):
        
        imin,imax = self.tvec.searchsorted((tmin, tmax))
        imax += 1
        
        if hasattr(self,'cache_fast'):
            data, data_err = self.cache_fast
            return self.tvec[imin:imax],data[imin:imax],data_err[imin:imax]
        
        
        #_, power, _ = self.get_data()
        #BOLOM::TOP.RAW:BOL_L23_V
        TDIcall = "\\BOLOM::TOP.RAW:"
        name = self.detectors_dict['U'][1]

        #self.MDSconn.openTree(self.tree,self.shot)
        #tvec = self.MDSconn.get('dim_of('+TDIcall+name+'_V)').data()/1000

              
        #data = zeros((len(tvec), self.nl),dtype=single)
        #data_err = zeros((len(tvec), self.nl),dtype=single)
        #tau = zeros(self.nl)
        #scrfact = zeros(self.nl)
        #gamma = zeros(self.nl)
        #kappa =  self.MDSconn.get( "\\BOLOM::TOP.PRAD_01.PRAD:KAPPA").data()[:self.nl]
        print('loading data...')
        #dt = mean(diff(tvec))
        #import IPython
        #IPython.embed()
        TDI_data = ['dim_of('+TDIcall+name+'_V)', ]
        TDI_params = []
        for cam, index in self.cam_ind.items():
            for i in index:
                name = self.detectors_dict[cam][i-index[0]]
                TDI_data.append(TDIcall+name+'_V')
                TDI_params.append(TDIcall+name+'_V:TAU')
                TDI_params.append(TDIcall+name+'_V:GAM')
                TDI_params.append(TDIcall+name+'_V:SCRFACT')

        self.MDSconn.openTree(self.tree,self.shot)
        params = self.MDSconn.get('['+','.join(TDI_params)+']').data()
        tau, gamma,scrfact =  params.reshape(-1,3).T
        #t = time.time()
        data = self.MDSconn.get('['+','.join(TDI_data)+']').data()
        #print(time.time()-t)

        #t = time.time()
        
 
        tvec = data[0]/1000 
        data   = data[1:].T
        #tau    = hstack(out[1::4])
        #gamma  = hstack(out[2::4])
        scrfact= scrfact**2#BUG in order to match the official bolometers data!!
        
        print('done')


        offset = (tvec<0)&(tvec  > -0.7)
        noise = std(data[offset],axis=0)

        data-= data[offset].mean(0)
        it = tvec.searchsorted(-.7)//2*2
        data= data[it:]
        tvec= tvec[it:]
        fs = 1/diff(tvec).mean()

        data -= outer(tvec,data[tvec>10].mean(0))/tvec[tvec>10].mean()
        import scipy.signal
        notch_filt = ((7.1,1),(6,1))
        for f0,width in notch_filt:
            b,a = scipy.signal.iirnotch(f0,f0/width,fs)
            data = scipy.signal.filtfilt(b, a, data,axis=0)


        #noise = std(data[offset],axis=0)
        drtau = 0.05*maximum(1,sqrt(noise/median(noise))) #more smoothing for noise channels
        
  
        method = 'richardson_lucy'
        print('Deconvolution using %s method'%method)
        t,x = self.deconvolve_bolo(tvec, data,tau,gamma,scrfact, method, drtau=drtau )
        print('Done')
        #import IPython
        #IPython.embed()
        from shared_modules import MovingAveradge
        downsample = 40
        
        self.tvec = t.reshape()
        #downsample = round(mean(diff(self.tvec))/mean(diff(t)))
        #x =  MovingAveradge(copy(x.T),downsample).T
        #t =  MovingAveradge(copy(t),downsample)
        pwr = zeros((len(self.tvec), self.nl), dtype='single')
        pwr_err = zeros((len(self.tvec), self.nl), dtype='single')
        pwr_err+= std(x[t<0])*sqrt(noise/median(noise))[None,:]*5+pwr.mean(1)[:,None]*0.10 #just guess 

        for i in range(self.nl):
            pwr[:,i] = interp(self.tvec,t,x[:,i])

        for cam, index in self.cam_ind.items():
            pwr[:,index] *= array(self.etendue[cam])[index-index[0]]*1e4 #W/m^2
            pwr_err[:,index] *= array(self.etendue[cam])[index-index[0]]*1e4 #W/m^2

        self.cache_fast =  pwr,pwr_err

        #import IPython
        #IPython.embed()
        ##pwr[self.tvec<3] = 0
        #plot(mean(pwr[(self.tvec>2)&(self.tvec<3)],0)/mean(power[(self.tvec>2)&(self.tvec<3)],0)/scrfact)
        #plot(mean(power[(self.tvec>2)&(self.tvec<3)],0))

        return self.tvec[imin:imax], pwr[imin:imax],pwr_err[imin:imax]
     
     
     
     
        
        t1,x1 = self.deconvolve_bolo(tvec, data,tau,gamma,scrfact, 'median', drtau=drtau )
        t2,x2 = self.deconvolve_bolo(tvec, data,tau,gamma,scrfact, 'golay', drtau=drtau )
        #t3,x3 = self.deconvolve_bolo(tvec, data,tau,gamma,scrfact, 'richardson_lucy', drtau=drtau )
        t4,x4 = self.deconvolve_bolo(tvec, data,tau,gamma,scrfact, 'wiener', drtau=drtau )
        
        
        
        pwr = zeros((len(self.tvec), self.nl), dtype='single')
        pwr_err = zeros((len(self.tvec), self.nl), dtype='single')
        
        import IPython
        IPython.embed()
        
        from shared_modules import MovingAveradge
        downsample = round(mean(diff(self.tvec))/mean(diff(t3)))
        x3 =  MovingAveradge(copy(x3.T),downsample).T
        t3 =  MovingAveradge(copy(t3),downsample)

        for i in range(self.nl):
            pwr[:,i] = interp(self.tvec,t3,x3[:,i])
            pwr_err[:,i]+= (interp(self.tvec, t2,x2[:,i])-pwr[:,i])**2
            pwr_err[:,i]+= (interp(self.tvec, t4,x4[:,i])-pwr[:,i])**2
            pwr_err[:,i] = MovingAveradge(sqrt(pwr_err[:,i]),20)

            #pwr_err[:,i] = std(vstack((interp(self.tvec, t2,x2[:,i]),interp(self.tvec, t4,x4[:,i]),pwr[:,i])),0)  #BUG UGLY
        for cam, index in self.cam_ind.items():
            pwr[:,index] *= array(self.etendue[cam])[index-index[0]]*1e4 #W/m^2
            pwr_err[:,index] *= array(self.etendue[cam])[index-index[0]]*1e4 #W/m^2

        self.cache_fast =  pwr,pwr_err

        import IPython
        IPython.embed()
        
        return self.tvec[imin:imax], pwr[imin:imax],pwr_err[imin:imax]
     
     
        
        f,ax = subplots(2,2, sharex=True, sharey=True)
        ax[0,0].plot(t1, x1)
        ax[0,0].set_title('median')
        ax[0,1].plot(t2, x2)
        ax[0,1].set_title('golay')
        #ax[1,0].plot(t3, x3)
        #ax[1,0].set_title('richardson_lucy')
        ax[1,1].plot(t4[:len(t4)//16*16].reshape(-1,16).mean(1), x4[:len(t4)//16*16].reshape(-1, 16,self.nl).mean(1))
        ax[1,1].set_title('wiener')
        show()
        from scipy.signal import fftconvolve, medfilt
        
        imshow(x1[:len(t4)//16*16].reshape(-1, 16,self.nl).mean(1), aspect='auto')
        imshow(x1 , aspect='auto')


        import IPython
        IPython.embed()
        
        
#imshow(x1, interpolation='nearest', aspect='auto',vmin=0)
        
        ch = 36
        f,ax = subplots(2,sharex=True)
        ax[0].plot(self.tvec,power[:,ch],'--' )
        ax[0].plot(t1,x1[:,ch])
        ax[0].plot(t2,x2[:,ch])
        ax[0].plot(t3,x3[:,ch],'k')
        ax[0].plot(t4,x4[:,ch],':')
        tau_ = (tau/gamma)[ch]
        h = exp(-(t4-t4[0])/tau_)
        dt= mean(diff(t4))
        retro4 = fftconvolve( x4[:,ch]*scrfact[ch]*dt/tau[ch],h,mode='full')[:-len(h)+1]
        
        h = exp(-(t3-t3[0])/tau_)
        dt= mean(diff(t3))
        retro3 = fftconvolve( x3[:,ch]*scrfact[ch]*dt/tau[ch],h,mode='full')[:-len(h)+1]
        
        ax[1].plot(tvec,data[:,ch])
        ax[1].plot(t4,retro4)
        ax[1].plot(t3,retro3,'k')

        show()
        
        
        
                
        #pwr[:,11]*=1.095
        #pwr[:,35]*=1.54
        #pwr[:,36]*=1.43
        #pwr[:,37]*=1.32
        #pwr[:,38]*=1.20


        #self.MDSconn.closeTree(self.tree,self.shot)
        
        
        
        
        
        
        
        
        
        
        
        
  
        #plot(time_,pwr[:,16] )
        #plot(self.tvec, power[:,16] )
        #plot(tvec, pwr2[:,16] )
        #plot(tvec, pwr3[:,16] )
        #plot(tvec, pwr4[:,16] )
        #plot(tvec, pwr5[:,16] )

        
        #ch = 33
        #f,ax=subplots(1,2, sharex=True)
        ##subplot(121)
        #ax[0].plot(tvec,data[:,ch])
        #tau_ = (tau/gamma)[ch]
        #h = exp(-(tvec-tvec[0])/tau_)
        #x_ = fftconvolve( pwr6[:,ch]*scrfact[ch]*dt/tau[ch],h,mode='full')[:-len(h)+1]
        #ax[0].plot(tvec, x_,'--')
        #ax[0].plot(tvec, x_- data[:,ch],'--')

        #y_ =  pwr6[:,ch]
        #ax[1].plot(tvec, y_)
        #ax[1].plot(self.tvec,power[:,ch])
        #ax[0].set_title(ch)
        #ax[1].grid(True)
        #show()
        
        #plot(tvec, pwr6[:,ch])
        #plot(self.tvec,power[:,ch]) 
        ##plot(time_, pwr[:,ch])

        
        #f,ax=subplots(1,2)
        #vmax = mean(pwr6*10)
        #ax[0].imshow(pwr6, interpolation='nearest', aspect='auto',vmin=0,vmax=vmax)
        #ax[1].imshow(power, interpolation='nearest', aspect='auto',vmin=0,vmax=vmax)
        
        #plot(self.tvec, power)
        
        
        
        #x = data[:,ch]
        #tau_ = (tau/gamma)[ch]
        #h = exp(-(tvec-tvec[0])/tau_)
        #x_ = fftconvolve( x,h,mode='full')[:-len(h)+1]*dt/tau[ch]
        #plot(x)
        #plot(x_*gamma[ch]+gradient(x_)*tau[ch]/dt)
        
        
        #plot(tvec, medfilt(x*gamma[ch]+gradient(x,dt)*tau[ch],301)/scrfact[ch])
        #plot(tvec,  pwr6[:,ch])
        
        

        #plot(pwr.mean(0)/power.mean(0))

        #y = np.fft.irfft(fy)
        #plot(y)
        
        
        #plot(tvec, data)
        
        #from scipy.optimize import curve_fit
        #f = lambda x,a,b,c: exp(-(x)/a)*b+c
        #p0 = 0.6, 2e4, 0
        #ioff()
        #tau_new = zeros(self.nl)
        #for i in range(self.nl):
            #ind = (tvec > 7.5)&(tvec < 10)
            #try:
                #popt,pcov =  curve_fit(f,tvec[ind],data[:,i][ind],p0)
            #except:
                #pass
            #tau_new[i] = popt[0]
            #title(i)
            #plot(tvec[tvec>7],data[:,i][tvec>7])
            #plot(tvec[tvec>7],f(tvec[tvec>7], *popt))
            #show()
            
           


        
        
        
        from skimage.restoration import denoise_tv_bregman
        
        
        #from sklearn.grid_search import GridSearchCV
        #for penalty in ('l1', 'tv1d'):
            #clf = FistaClassifier(penalty=penalty)
            #clf.alpha = 1
            #gs = GridSearchCV(clf, {'alpha': np.logspace(-3, 3, 10)})
            #gs.fit(X, y)
            #coefs = gs.best_estimator_.coef_
            #plt.plot(coefs.ravel(), label='%s penalty' % penalty, lw=3)
                
                
        #plot(tau_new)
        #plot(tau/gamma)
        #data = medfilt(data[:,28],3)
        
        #plot(medfilt(data[:,28],3))
        
        #kouknoiut n kana 14,16!!!
        #potlacovat gradient v log scale? 

        #plot(x)
        
        
        
        
        
        from scipy.stats.mstats import mquantiles
        fdrdt = np.fft.rfft(drdt[:,29],axis=0)
        fmax =  mquantiles(abs(fdrdt), 0.99)
        fdrdt[abs(fdrdt) > fmax] = 0
        
        drdt = np.fft.irfft(fdrdt,axis=0)
        
        t = arange(1000)/100.*2*pi
        t = tvec
        a = sin(t)

        plot(gradient(a, mean(diff(t)))+a)
        plot(a)
        
        
        
        N = len(a)
        k = 2*pi/mean(diff(t))/N*arange(N/2+1)
        
        da = np.fft.irfft( k*np.fft.rfft(a))
        plot(da,'--')
        

        pwr = (gamma*data+tau*drdt)/scrfact
        pwr = pwr.reshape(-1, 256,self.nl).mean(1)
        tvec_ = tvec.reshape(-1, 256).mean(1)
        pwr[:,[0,12,15,21]] = 0
        
        plot(pwr)
        figure()
        plot(self.tvec,power)
        show()
        
        
        imshow(pwr, interpolation='nearest', aspect='auto');show()
        
        
        plot(tvec_, pwr[:,29])
        plot(self.tvec,power[:,29])



        
        ion()       

        import IPython
        IPython.embed()
        
        
 
        
        
    def deconvolve_bolo(self,tvec, data,tau,gamma,scrfact, filtertype, drtau ):

        pwr = zeros((len(tvec), self.nl))
        dt = mean(diff(tvec))
        #import IPython
        #IPython.embed()
        
        if filtertype == 'richardson_lucy': #Richardson-Lucy deconvolution.

            data = medfilt(data,(3,1)) #remove spikes
            s = 8
            offset = std(data[tvec<0],0)*2            
            tvec = tvec[:len(tvec)//s*s].reshape(-1,s).mean(1)
            n = 2**int(ceil(log2(sum(tvec>0))))
            data = data[:len(data)//s*s].reshape(-1,s,self.nl).mean(1)[-n:]
            tvec = tvec[-n:]
            dt*=8
            #multiprocessing!!!
            
            args = []
            numiter = 1000
            relax = 1.6

            for m in range(self.nl): 
                tau_ = (tau/gamma)[m]
                h = exp(-(tvec-tvec[0])/tau_)
                nitr = int(numiter*drtau[m]/0.05)
                args.append((data[:,m],h,nitr,offset[m],relax))

        
            pool = Pool(cpu_count())
            out = map(richardson_lucy,args)
            pool.close()
            pool.join()
            
            
            pwr = vstack(list(out)).T*(tau/dt/scrfact/array([sum(a[1]) for a in args]))
  
            
        
        
        
                
            
        if filtertype == 'median': #median filter  

            tvec     = tvec[:len(tvec)//4*4].reshape(-1,4).mean(1)
            tdel      = dt*4
            #args = [(data[:,i],drtau[i],tdel,tau[i],gamma[i], scrfact[i] ) for i in range(self.nl)]
            one = ones(len(tau))
            args = zip(data.T,drtau*one,tdel*one,tau,gamma,scrfact)

            pool = Pool(cpu_count())
            out = pool.map(medfilt_deconv,args)
            pool.close()
            pool.join()
            
            pwr = vstack(out).T
  
            
        

            
      
     
        
        elif filtertype == 'golay': #Savitzky-Golay filter   

            
            from scipy.signal import savgol_coeffs
            N = 2                 # Order of polynomial fit
            F = int_(drtau/dt/2)*2+1                # Window length
            # Calculate S-G coefficients
            for m in range(self.nl):            
                g = [savgol_coeffs(F[m],N, i) for i in range(2)]
                y = data[:,m]          
                SG0 = convolve(g[0][::-1],y, mode='same')
                SG1 = convolve(g[1][::-1],y, mode='same')/dt # Turn differential into derivative
                pwr[:,m] = (SG0*gamma[m]-SG1*tau[m])/scrfact[m]
            

        
        elif filtertype == 'butter':
            from scipy.signal import butter, filtfilt
            facq = 1./ mean(diff(tvec))
            for m in range(self.nl):            
                fmax = 1./drtau[m]
                b,a = butter(5,fmax/facq) # Butterworth filter at 100 Hz
                rp=gradient(data[:,m])/dt
                rp=filtfilt(b,a,rp)
                pwr[:,m] =  (data[:,m]*gamma[m]+rp*tau[m])/scrfact[m]

        
            
        elif filtertype == 'spline':
            from scipy.interpolate import UnivariateSpline

            for m in range(self.nl):            
                s =UnivariateSpline(tvec.reshape(-1,16).mean(-1), data[:,m].reshape(-1,16).mean(-1),s=.06,k=3)
                yout = s(tvec, nu=0)
                youtp = s(tvec, nu=1)
                pwr[:,m] =  (yout*gamma[m]+youtp*tau[m])/scrfact[m]

        elif  filtertype == 'bessel':
            #---filter preparation---
      

                        
            filtercutoff=30 #Hz This value is suggested by the literature.
            f_nq=1/(tvec[1]-tvec[0])/2
                         

            from scipy.signal import bessel, bilinear,lfilter

            banalderivfordigit,aanalderivfordigit=bessel(3,filtercutoff/f_nq,analog=False)
            banalderivfordigit=r_[banalderivfordigit[1:],0]
            bderivdigit,aderivdigit=bilinear(banalderivfordigit,aanalderivfordigit, f_nq*2)
            #bderivdigit,aderivdigit = banalderivfordigit,aanalderivfordigit
            #---filter preparation---
            #print  (abs(roots(aanalderivfordigit))>=1)
            #print  (abs(roots(aderivdigit))


            for m in range(self.nl):            
                #BUG instability in the filter!!!!
                dynamic_part=lfilter(bderivdigit,aderivdigit,data[:,m] ) 
                pwr[:,m] =  (data[:,m]*gamma[m]+dynamic_part*tau[m])/scrfact[m]

            
        elif  filtertype == 'wiener':

            
            pwr = zeros((len(tvec), self.nl))
            for ch in range(self.nl):
                x = medfilt(data[:,ch],3)
                fx = np.fft.rfft(x)
                fx[2000:]=0
                tau_ = (tau/gamma)[ch]
                h = exp(-(tvec-tvec[0])/tau_)
                H = np.fft.rfft(h)
                N = abs(np.fft.rfft(x*(tvec<0)))
                S=2e-4/exp(drtau[ch]/0.05)
                smoothN = savitzky_golay(N, 101,1,firstvals=mean(N[:200]))
                G = conj(H)*S/(H*conj(H)*S+smoothN)
                pwr[:,ch]=  np.fft.irfft(G*fx)*tau[ch]/dt/scrfact[ch]
           #embed()
        return tvec, pwr
       
    def get_data_realtime(self,tmin=-infty,tmax=infty):
        print('Loading realtime bolometer data')
       
        data = []
 
       
        for cam, index in self.cam_ind.items():
            for i in index:
                name = self.detectors_dict[cam][i-index[0]]
                
                TDIcall = f'PTDATA2("DGSDPWR{name[4:]}",{self.shot},1)'

                data.append( self.MDSconn.get(TDIcall).data()) #W
                data[i] *= self.etendue[cam][i-index[0]]*1e4 #W/m^2
        
        data =  np.array(data).T
 
        
        data = data[:len(data) // self.downsample * self.downsample].reshape(-1, self.downsample, data.shape[1]).mean(1)
      
        data_err = zeros_like(data,dtype=single)
        #down
       
        shot_data = data#[(self.tvec > 0)&(self.tvec < 5)]
  
        from scipy.signal import butter, sosfiltfilt
        sos = butter(4, 0.02, output='sos')
        data_smooth = sosfiltfilt(sos, shot_data, axis=0)
   
        #Mean absolute deviation
        data_err = 1.22* np.exp(sosfiltfilt(sos,sosfiltfilt(sos, np.log(np.abs(shot_data-data_smooth+1e-6)), axis=0), axis=0))
        data_err += np.abs(data_smooth) * 0.05 #5% calibration error 
        data_err = np.single(data_err)
  
        self.cache = data,data_err
        
                
        mdata = np.mean((shot_data),0)
        likely_invalid =  (data_err[0]  <1 ) |(mdata < np.median(mdata)/10)#(mdata <  data_err[0] *  
        
        try:
            import config
            config.wrong_dets_pref = np.where(likely_invalid)[0]
        except: 
            pass
        
        
       
    def get_data(self,tmin=-infty,tmax=infty):
         
        #try load realtime data 
        if self.real_time_data and not hasattr(self,'cache'):
            self.get_data_realtime(tmin,tmax)
        
        imin,imax = self.tvec.searchsorted((tmin, tmax))
        imax += 1
        
        if hasattr(self,'cache'):
            data, data_err = self.cache

         
            return self.tvec[imin:imax],data[imin:imax],data_err[imin:imax]
        
        

    

        
        ntim = len(self.tvec)

        data = zeros((ntim, self.nl),dtype=single)
        data_err = zeros((ntim, self.nl),dtype=single)
        
        TDIcall = "\\BOLOM::TOP.PRAD_01.POWER:"
        self.MDSconn.openTree(self.tree,self.shot)
        for cam, index in self.cam_ind.items():
            for i in index:
                name = self.detectors_dict[cam][i-index[0]]
                data[:,i] = self.MDSconn.get(TDIcall+name+'_P').data() #W
                data[:,i] *= self.etendue[cam][i-index[0]]*1e4 #W/m^2

        #kappa =  self.MDSconn.get( "\\BOLOM::TOP.PRAD_01.PRAD:KAPPA").data()  #factors to calculate total radiated power 

        #data[:,0] =0
     
        self.MDSconn.closeTree(self.tree,self.shot)

        #data_err += std(data[self.tvec < 0],0)*2
        #data_err += 0.1*mean(data,1)[:,None]
        
        
        shot_data = data#[(self.tvec > 0)&(self.tvec < 5)]
        #data_err += np.median(np.abs(np.diff(shot_data[::10],axis=0)),0)/2
        
         
        
        import matplotlib.pylab as plt
        from scipy.signal import butter, sosfiltfilt
        sos = butter(4, 0.02, output='sos')
        data_smooth = sosfiltfilt(sos, shot_data, axis=0)
        
        
        
        init_err = data[self.tvec<0].std(0)/2
        #Mean absolute deviation
        data_err = 1.22* np.exp(sosfiltfilt(sos,sosfiltfilt(sos, np.log(np.abs(shot_data-data_smooth+1e-6)+init_err), axis=0), axis=0))
        data_err += np.abs(data_smooth) * 0.05 #5% calibration error 
        data_err = np.single(data_err)
        
        #median absolute deviation
        #data_err += 1.22*np.mean(abs(shot_data-data_smooth), 0)
        
       # for i in range(50):
        #    plt.plot(shot_data[:,i])
         #   plt.plot(data_smooth[:,i])
         #   plt.title(i)
         #   plt.show()
        
        
       # for i in range(len(data)):
        #    plt.title(i)
        #    plt.errorbar(self.tvec, data[:,i], data_err[:,i])
        #    plt.plot(self.tvec[(self.tvec > 0)&(self.tvec < 5)], data_smooth[:,i])
        #   plt.show()
 
 
        
        self.cache = data,data_err
        
        
        mdata = np.mean((shot_data),0)
        likely_invalid =  (data_err[0]  <1 ) |(mdata < np.median(mdata)/10)#(mdata <  data_err[0] * 3) |
       # likely_invalid = (mdata <  data_err[0] * 3) | (data_err[0]  <1 )|(mdata < np.median(mdata)/10)
        #np.mean(np.abs(np.diff(shot_data[::10],axis=0)),0)*2 
        
       # plt.plot(mdata)
        #plt.plot( data_err[0] * 3)
        #plt.plot(np.arange(len(mdata))[likely_invalid], mdata[likely_invalid],'o')
        
       # plt.plot(np.median(np.abs(np.diff(shot_data[::10],0)),0), '--'  )
        #plt.plot(  np.mean(np.abs(np.diff(shot_data[::10],axis=0)),0)   , '--'  )
        
       # plt.show()
        
        try:
            import config
            config.wrong_dets_pref = np.where(likely_invalid)[0]
        except: 
            pass
        
        #import config
        #config.wrong_dets_pref = np.where(likely_invalid)[0]
        
        #print(config.wrong_dets_pref)
       # data_err[:,likely_invalid] = np.inf
 
        
        return self.tvec[imin:imax],data[imin:imax],data_err[imin:imax]

    
   
   
   
   
       
 
    def load_geom(self,path):
 
        self.geometry_version = 0
        
        
        
        #f,axis = subplots(1,sharex=True,sharey=True)
        #f.subplots_adjust(hspace=0.10, wspace = 0.1)
        #from matplotlib.collections import PolyCollection
        ###prepare files with cordinates

        for i,det in enumerate(self.cams):

            R1 = self.R_start[det]
            Z1 = self.z_start[det]
            R2 = self.R_end[det] 
            Z2 = self.z_end[det]

            
            m_alpha = abs(abs(arctan2(Z2-Z1,R2-R1))-pi/2)
            alphas = arctan2(Z2-Z1,R2-R1)
            alphas = unwrap(alphas)
            m_alpha = mean(m_alpha)
            
            
            xfile = open(self.geometry_path+'detector_%s_x.txt'%det,'w')
            yfile = open(self.geometry_path+'detector_%s_y.txt'%det,'w')
            
            
            alpha = arctan2(Z2-Z1,R2-R1) #to same jako alphas? 
            
            theta = ones_like(alpha)*median(abs(diff(alphas)))/2 #exactly nonoverlaping LOS 
            theta/= 2  #actual LOS are narrower...
            #theta/= 10  #BUG 

            L = hypot(R2-R1, Z2-Z1)

            verts = []

            if m_alpha < pi/4:       
                Lr = L*abs(sin(alpha)/sin(pi-theta-alpha))
                Ll = L*abs(sin(pi-alpha)/sin(alpha-theta))
                R21 = R1+Ll*cos(alpha-theta)
                R22 = R1+Lr*cos(alpha+theta)
                
                savetxt(self.geometry_path+'detector_%s_x.txt'%det, c_[R21,R22,R1], fmt='%5.4f')
                savetxt(self.geometry_path+'detector_%s_y.txt'%det, c_[Z2,Z1], fmt='%5.4f')
                for r1,r21,z1,r22,z2 in zip(R1,R21,Z1,R22,Z2):
                    verts.append([[r1,z1],[r21,z2 ],[r22,z2]])

            else:
                
                
                #ind = tan(pi-abs(alpha)-theta)*tan(pi-abs(alpha)+theta) >= 0
                 ###solve special case for almost exactly vertical LOS
                #Z21 = copy(Z2)
                #Z22 = copy(Z2) +1e-2
                
                #Z21[ind] = (Z1-(R2-R1)*tan(pi-abs(alpha)-theta)*sign(alpha))[ind]
                #Z22[ind] = (Z1-(R2-R1)*tan(pi-abs(alpha)+theta)*sign(alpha))[ind]
                
                Z21 = (Z1-(R2-R1)*tan(pi-abs(alpha)-theta)*sign(alpha))
                Z22 = (Z1-(R2-R1)*tan(pi-abs(alpha)+theta)*sign(alpha))
                
                
                
                

                savetxt(self.geometry_path+'detector_%s_x.txt'%det, c_[R2,R1], fmt='%5.4f')
                savetxt(self.geometry_path+'detector_%s_y.txt'%det, c_[Z21,Z22,Z1], fmt='%5.4f')  
   
                for z1,z21,r1,z22,r2 in zip(Z1,Z21,R1,Z22,R2):
                    verts.append([[r1,z1],[r2,z21],[r2,z22]])


            #color = ('r', 'b', 'g', 'k', 'm', 'y', 'k', 'gray')

            #from matplotlib.collections import PolyCollection
            #title(det)
            #verts = array(verts)
            #ax = gca()
            #coll = PolyCollection(verts, color=color[i],edgecolors='none',alpha = 0.1)
            #ax.add_collection(coll)
            #ax.autoscale_view()
            #ax.axis('equal')
            #ax.set_xlim(1.1, 2.15)
        #show()
    


 

def main():
    import os
    import os,sys
    #, shot, geometry_path,MDSconn
   # c = mds.Connection('localhost' )
    c = mds.Connection('atlas.gat.com' )
    
    bolo = loader_BOLO(200380,'/home/tomas/tomography/geometry/DIIID/BOLO/',c )
    T = time.time()
    bolo.get_data(1.9,3)
    #bolo.get_total_rad(4,6)

    #print('loaded in ',time.time()-T)
    #bolo.load_geom()


    

if __name__ == "__main__":
    main()
