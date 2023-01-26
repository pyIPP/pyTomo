#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,sys
from numpy import *
import time
import config
from collections import OrderedDict
from scipy.interpolate import interp1d
from scipy.stats.mstats import mquantiles
import threading
from shutil import copyfile


from tokamak import Tokamak
from shared_modules import in1d,fast_svd,read_config, warning


sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib/')

global loader


class ASDEX(Tokamak):
    """
    Representation of class of tokamak. Contain basic information about geometry, borde, magnetic field and basic setting used for recontruction and postprocessing.

    ASDEX(input_path, 'SXR-slow',  shot, nx, ny,virt_chord )


    :var str name:  Name of tokamak
    :var int index: Index of the tokamak
    :var double xmin, xmax, ymin, ymax: boundaries of recontructed area of tokamak emissivity
    :var str geometry_path: Path to saved geometry, boundary or data
    :var double sigma:  Expected noise in data  += ``sigma``* sqrt(data) * sqrt(max(data))  (Poisson like noise)
    :var double min_error: Expected noise in data  += ``min_error`` * sqrt(max(data))
    :var int shot:  NUmber of selected shot
    :var int norm: Used to normalize input parameters to meters
    :var str t_name:  Physical units used to measure time for the tokamak
    :var int npix:  total number of pixels (nx*ny)
    :var double dx: Width of one cell
    :var double dy: Height of one cell
    :var bool local_network: True if it is possible to load data from network

    :var array data: Input signals from main diagnostics
    :var array tvec: Time vector that belongs to the ``data``
    :var array tsteps: Lenght of tvec /  data
    :var array error: Expected errors of data,
    :var array dets: Array of all detectors in diagnostics
    :var array wrong_dets: Array of deleted detectors in diagnostics
    #:var array magx, magy: Arrays of magnetic lines possitions
    :var spmatrix Tmat: Matrix of geometry
    :var array Xchords, Ychords: Arrays of chords coordinates, used in ``make_plots``
    :var array boundary: Coordinates of points on boundary of tokamak
    :var array emiss_0: First guess of plasma emissivity, improve convergence
    """


    def __init__(self, input_diagn, input_parameters,  only_prepare = False, load_data_only=False ):

        """ Initialize new tokamak and set default values.

            :param dict input_parameters:   List of selected paramerers
            :param str input_diagn:  Name of main diagnostics
            :param array vmax: artificial shift of magnetic field

        """
    

        name  = 'ASDEX'
        
        shot = input_parameters['shot']
        self.mag_diag = input_parameters['mag_diag']
        self.mag_exp  = input_parameters['mag_exp']
        self.mag_ed   = input_parameters['mag_ed']
        try:
            from . import dd
            self.dd = dd.shotfile()
        except:
            self.dd = None
        
        self.local_path =  input_parameters['local_path']
        self.program_path =  input_parameters['program_path']
        


        path = 'geometry/'+name
        if not os.path.exists(self.local_path+path): os.makedirs(self.local_path+path)

        path += '/'+ input_diagn
        if not os.path.exists(self.local_path+path): os.makedirs(self.local_path+path)

        
        if input_diagn in ["SXR",'SXR_fast'] :
            coord= [1.05, 2.2, -1.02, 0.94]  
            rad_profile = 'Emiss'

            if input_diagn == 'SXR_fast':
                path = 'geometry/'+name+'/'+ 'SXR'

        elif input_diagn in ['BOLO', 'AXUV']:               
            coord= [1, 2.25, -1.3, 1.25]   
            rad_profile = 'Emiss'

        #self.zoom_coord = [1.1,1.46,2,2,1 ],[-.5,-.83,-.7,-1.3,-1.3 ]
        
        self.geometry_path_program = self.program_path+path
        path = self.local_path+path
        
        if self.geometry_path_program+'/border.txt' != path+'/border.txt':
            copyfile(self.geometry_path_program+'/border.txt',path+'/border.txt' )
   

        self.input_diagn = input_diagn
        
        self.min_tvec = 0
        self.max_tvec = 9

        sigma = 0.02              
        min_error = 0.01           
        norm = 1.0           
        t_name = 's'              
        
        Tokamak.__init__(self, input_parameters, input_diagn, 
                        coord,  name, norm, sigma, min_error, path,t_name)
        
        
        
        #print self.geometry_path,input_diagn,name
        
        #exit()
        self.Tokamak = Tokamak
        self.allow_negative = False
        self.dets_index = (1,) #default
        Tokamak.rad_profile = rad_profile
        self.rad_profile = rad_profile

        #share "common wisdom", 
        calib_path = path+'/calibration/'
        if not os.path.exists(calib_path):  
            try:
                os.mkdir(calib_path)
            except:
                print('error:os.mkdir(calib_path)')
                raise
        wch_path = path+'/wrong_channels/'
        if not os.path.exists(wch_path): 
            try:
                os.mkdir(wch_path)
            except:
                print('error: os.mkdir(wch_path)')
                raise
        
            
        calib_path_orig = self.geometry_path_program+'/calibration/%d.txt'%self.shot
        if (not os.path.isfile(calib_path+'%d.txt'%self.shot) or not config.useCache) and os.path.isfile(calib_path_orig):
            if calib_path_orig!= calib_path+'%d.txt'%self.shot:
                copyfile(calib_path_orig,calib_path+'%d.txt'%self.shot)
            
        wch_path_orig = self.geometry_path_program+'/wrong_channels/%d'%self.shot
        if (not os.path.isfile(wch_path+'%d'%self.shot) or not config.useCache) and os.path.isfile(wch_path_orig):
            if wch_path_orig!= wch_path+'%d'%self.shot:
                copyfile(wch_path_orig,wch_path+'%d'%self.shot)

        if os.path.exists(self.geometry_path_program+'/virt_chord_profile.txt'):
            if self.geometry_path_program+'/virt_chord_profile.txt'!= path+'/virt_chord_profile.txt':
                copyfile(self.geometry_path_program+'/virt_chord_profile.txt',path+'/virt_chord_profile.txt')
        

        #print self.Tokamak.rad_profile
        

        if not load_data_only:
     
            self.load_mag_equilibrium()

            self.min_tvec = 0
            self.max_tvec = self.tsurf[-1]+0.1

            self.VesselStructures()

            self.load_LBO()
            self.load_ICRH()
            
            if self.input_diagn in ('SXR','SXR_fast'):
                self.load_others()
    
  
        if input_diagn in ("SXR",'SXR_fast') :
            if self.shot < 20807 and self.shot > 1000 :
                raise Exception('Old SXR system is not suppoted!!')
            
            self.SXR()
            self.index = 9 if input_diagn == 'SXR' else 10 
            
            

            #else:
                #config.wrong_dets_defined = []
            #print ' config.wrong_dets_defined', config.wrong_dets_defined
            
            
        if input_diagn == 'BOLO':
            self.BOLO()
            self.index = 11

        if input_diagn == "AXUV":
            self.AXUV()
            self.index = 12 
            
            
        self.max_tvec = min(self.max_tvec,  self.tvec[-1])
        try:
            self.load_geom()
        except:
            print('load_geom failured')
            #pass
        
        Tokamak._post_init(self)

        if not only_prepare:
            self.prepare_tokamak()
            
        #if input_diagn == 'SXR_fast':
            #global loader
            #wrong_dets_pref = loader.all_los[loader.cam_ind['F']]
            #wrong_dets_pref_ind = where(in1d(loader.all_los,wrong_dets_pref))[0]
            #config.wrong_dets_pref = unique(r_[config.wrong_dets_pref,wrong_dets_pref_ind])
            


            
        
    def __del__(self):
        if not self.dd is None: self.dd.Close()
   
     
     
     
     
    def SXR(self):
        
        
        self.default_calb = 'variable'   
        self.allow_self_calibration = False
        self.sigma = 0.006 #just guess!!
        self.min_error = 0.004#just guess!!
        
        
        self.shift = 13
        self.boundary_lcfs = True
        

        from .SXR.get_sxr import loader_SXR
        global loader

        fast_data = self.input_diagn == 'SXR_fast'
        loader = loader_SXR(self.shot,self.geometry_path, self.dd, fast_data, 
                                 experiment='AUGD',edition=0)

        
        self.detectors_dict= loader.detectors_dict
        
        
   
        self.calb_0 = loader.calb_0
        self.nl =  loader.nl
 
        self.dets_index = loader.dets_index
        self.dets = loader.dets
        self.tvec = loader.tvec
        self.Phi  = loader.Phi
        self.SampFreq  = loader.SampFreq

        if self.input_diagn == 'SXR_fast':
            self.sample_freq = 5e5
        else:
            self.sample_freq = 5e3
  
        
        
    def BOLO(self):
        
        self.allow_self_calibration = True
        self.use_median_filter = False
        #self.weight_power = 1
        self.beta = 0
        self.shift = 5
        self.boundary_lcfs = False


        from .BOLO.get_bolo import loader_BOLO
        
        global loader

        loader = loader_BOLO(self.shot,self.geometry_path, self.dd, experiment='AUGD',edition=0)
                
        self.detectors_dict= loader.detectors_dict
   
        self.calb_0 = loader.calb_0
        self.nl =  loader.nl
        self.tvec = loader.tvec

        self.dets_index = loader.dets_index
        self.dets =  loader.dets
        self.sample_freq =  1/loader.dt

        self.total_rad = loader.get_total_rad()
        
    def AXUV(self):
        
        
        
        self.allow_self_calibration = False
        self.boundary_lcfs = False
        self.shift = 16
        self.sigma = 0.01
        self.min_error = 0.01
        
        
        from .AXUV.get_axuv import loader_AXUV
        global loader

        
        loader = loader_AXUV(self.shot,self.geometry_path, self.dd, experiment='AUGD',edition=0)
                  
        self.detectors_dict= loader.detectors_dict
   
        self.calb_0 = loader.calb_0
        self.nl =  loader.nl
        
        self.dets_index = loader.dets_index
        self.dets = loader.dets
        self.tvec = loader.tvec
        self.tvec = linspace(loader.tvec[0],loader.tvec[-1], loader.tvec.size/loader.n_smooth)
        self.sample_freq =  500e3/loader.n_smooth

        
    def load_geom(self):
        global loader

        try:
        #if not self.dd is None:
            loader.load_geom(self.geometry_path_program)
        except Exception as e:
        #else:
            print('Geometry could not be loaded:: ')
            
        #print loader.geometry_version
        #exit()
        if hasattr(loader,'geometry_version'):
            self.geometry_version =  loader.geometry_version


        
    def get_data(self, failsave = True,tmin = -infty,tmax = infty):
        
        #print( 'tmin, tmax',tmin, tmax)
        global loader

        #exit()
        t = time.time()

        tvec, data, error = loader.get_data(tmin,tmax)


 
        
        wrong_dets_damaged = where(in1d(loader.all_los,loader.wrong_dets_damaged))[0]
        
        self.wrong_dets_damaged  = unique(r_[self.wrong_dets_damaged,wrong_dets_damaged])
                
       
        ind_t = slice(None,None)
        if  len(tvec) > 5000:  ind_t = random.randint(len(tvec),size=1000)
        
        mean_data = nanmean(data[ind_t],axis=0)
        #print 'asdex config.wrong_dets_defined ', config.wrong_dets_defined 
        

       #include also "nonstatistical noise"  by min_error and sigma
        error =  hypot(error, single(mean_data*self.sigma),out=error)
        error += abs(mean_data)*self.min_error

        print('loading time: %.1f s'%( time.time()-t))
        
        #print data.shape, self.nl
        return  tvec,data, error 
         
        
        
        
  
    def load_LBO(self):
        #print 'load_LBO'
        lbo_path = self.local_path +"/geometry/ASDEX/LBO/"
        
        if not os.path.exists(lbo_path):
            try:
                os.mkdir(lbo_path)
            except:
                print('error: os.mkdir(lbo_path)')
 #               except
                
        lbo_path += "LBO_"+str(self.shot)+".txt"
        
        #import IPython
        #IPython.embed()
        
        try:
            self.impur_inject_t = atleast_1d(loadtxt(lbo_path))             
            if all(isnan(self.impur_inject_t)):
                self.impur_inject_t = None
        except Exception as e:
            #print e
            #try:
            if not self.dd is None and self.dd.Open('LBO',self.shot):
                    
                #20Hz synchronization with the time 0
                LBOtvec = self.dd.GetTimebase('Time_DAC')
                LBO = self.dd.GetSignalGroup('DAC-Sig').T
                self.dd.Close()
 
                impur_inject_t = linspace(0,8,8*20+1)
                LBO_ind = zeros_like(impur_inject_t, dtype='bool')
                for i in where(any(LBO>1,axis=0))[0]:  
                    LBO_ind[(impur_inject_t>=LBOtvec[i]+1e-3)&(impur_inject_t<=LBOtvec[i+1]-1e-3)]=True

                lbo_times = impur_inject_t[LBO_ind]
                savetxt(lbo_path, lbo_times,fmt='%3.5f' )
                self.impur_inject_t = lbo_times
                #debug( 'LBO time:'+str(lbo_times))
            else:
            #except Exception as e:
                #print e

                self.impur_inject_t = None
                savetxt(lbo_path, (nan,))
         
        #print  self.impur_inject_t


    def load_others(self):
        #load centrifugal asymmetry profile and bremsstrahlung level 
        
 
        data_path = self.geometry_path+'/Bremsstrahlung_%d.npz'%self.shot
        
        if os.path.isfile(data_path) or self.dd is None:
            return 
        
   
    
        diag = 'EQH' if self.mag_diag == 'TRA' else self.mag_diag 
        stat = self.eqm.eq_open
        if not stat:
            stat = self.eqm.Open(self.shot, diag=diag, exp=self.mag_exp, ed=self.mag_ed)
        if not stat:
            print('Warning: equlibrium for shot:%d diag:%s  exp:%s  ed:%d  was not found!! FPP will be used'%(self.shot,diag,self.mag_exp,self.mag_ed))
            self.mag_diag = 'FPP'
            self.mag_exp = 'AUGD'
            self.mag_ed = 0
            stat = self.eqm.Open(self.shot, diag=self.mag_diag, exp=self.mag_exp, ed=self.mag_ed)
            if not stat:
                raise Exception('equlibrium for shot:%d diag:%s  exp:%s  ed:%d  was not found!! use another one'%(self.shot,self.mag_diag,self.mag_exp,self.mag_ed))
            
            

    
        if self.dd.Open('IDA',self.shot):
            
            
            rhop = self.dd.GetAreabase('rhop')
            tvec = self.dd.GetTimebase('time')
            ne = self.dd.GetSignalGroup('ne')
            Te = self.dd.GetSignalGroup('Te')
            
            ne_low = self.dd.GetSignalGroup('ne_lo')
            Te_low = self.dd.GetSignalGroup('Te_lo')
            ne_high = self.dd.GetSignalGroup('ne_up')
            Te_high = self.dd.GetSignalGroup('Te_up')

            self.dd.Close()
            

            
            ind_t = slice(*tvec.searchsorted((self.min_tvec, self.max_tvec)))
            rhop = rhop[ind_t].mean(0)
            ind_r = slice(0,rhop.searchsorted(1),3)
            ne = ne[ind_r,ind_t]
            Te = Te[ind_r,ind_t]
            ne_low = ne_low[ind_r,ind_t]
            Te_low = Te_low[ind_r,ind_t]
            ne_high = ne_high[ind_r,ind_t]
            Te_high = Te_high[ind_r,ind_t]
      
            rhop = rhop[ind_r]
            tvec = tvec[ind_t]



        elif self.dd.Open('VTA',self.shot):
            


            tvec = self.dd.GetTimebase('TIM_CORE')
            ind_t = slice(*tvec.searchsorted((self.min_tvec, self.max_tvec)))

            tvec = tvec[ind_t]

            R = self.dd.GetSignalGroup('R_core', nbeg=ind_t.start+1,nend=ind_t.stop)
            z = self.dd.GetSignalGroup('Z_core')
            z,R = meshgrid(z,R)
            rho = self.eqm.rz2rho(R,z, tvec)
    

            Te_     = self.dd.GetSignalGroup('Te_c',   nbeg=ind_t.start+1,nend=ind_t.stop)
            Ne_     = self.dd.GetSignalGroup('Ne_c',   nbeg=ind_t.start+1,nend=ind_t.stop)
            Ne_up_  = self.dd.GetSignalGroup('Neupp_c',nbeg=ind_t.start+1,nend=ind_t.stop)
            Ne_low_ = self.dd.GetSignalGroup('Nelow_c',nbeg=ind_t.start+1,nend=ind_t.stop)
            Te_up_  = self.dd.GetSignalGroup('Teupp_c',nbeg=ind_t.start+1,nend=ind_t.stop)
            Te_low_ = self.dd.GetSignalGroup('Telow_c',nbeg=ind_t.start+1,nend=ind_t.stop)
            
            self.dd.Close()
            
            rhop = linspace(0,1,100)
            
            valid = (Te_>10)&(Ne_>10)
            valid_ind = any(valid,axis=1)
            tvec = tvec[valid_ind]
            
            Ne_low_[Ne_low_<=0] = nan
            Te_low_[Te_low_<=0] = nan

            Te = zeros((len(rhop), sum(valid_ind)))
            ne = zeros((len(rhop), sum(valid_ind)))
            ne_high = zeros((len(rhop), sum(valid_ind)))
            ne_low = zeros((len(rhop), sum(valid_ind)))
            Te_high = zeros((len(rhop), sum(valid_ind)))
            Te_low = zeros((len(rhop), sum(valid_ind)))

            for i,ii in enumerate(where(valid_ind)[0]):
                ind = valid[ii]
                sind = arange(rho.shape[1])[ind][argsort(rho[ii,ind])]
                Te[:,i] = interp(rhop, rho[ii,sind], Te_[ii,sind],left=nan)
                ne[:,i] = interp(rhop, rho[ii,sind], Ne_[ii,sind],left=nan)
                ne_high[:,i] = interp(rhop, rho[ii,sind], Ne_up_[ii,sind],left=nan)
                ne_low[:,i] = interp(rhop, rho[ii,sind], Ne_low_[ii,sind],left=nan)
                Te_high[:,i] = interp(rhop, rho[ii,sind], Te_up_[ii,sind],left=nan)
                Te_low[:,i] = interp(rhop, rho[ii,sind], Te_low_[ii,sind],left=nan)

            
        else:
            savez(data_path)
            return 
            
        
        
        
        
        
        self.dd.Open('JOU',self.shot)
        He = ''.join([g.decode('utf-8') for g in  self.dd.GetParameter('FILLING','GAS_He')])
        self.dd.Close()
          
        Zeff = 1.1  #BUG just guess!!
        if He == 'He':  Zeff = 2.1

        from .radiation.Bremsstrahlung import CalcBackgroundRadiation,CalcRadiation,CalcZ
        
        
        BR = ones_like(Te)*nan
        BR_low = ones_like(Te)*nan
        BR_up = ones_like(Te)*nan
        W = ones_like(Te)*nan
        W_low = ones_like(Te)*nan
        W_up = ones_like(Te)*nan
        
        ind = isfinite(Te)&isfinite(ne)
        BR[ind] = CalcBackgroundRadiation(ne[ind],Te[ind], Zeff)
        W[ind]  = CalcRadiation(ne[ind],Te[ind],1e-4)
        ind = isfinite(Te_low)&isfinite(ne_low)
        BR_low[ind] = CalcBackgroundRadiation(ne_low[ind],Te_low[ind], Zeff/1.1)
        W_low[ind]  = CalcRadiation(ne_low[ind],Te_low[ind],1e-4)
        ind = isfinite(Te_high)&isfinite(ne_high)
        BR_up[ind] = CalcBackgroundRadiation(ne_high[ind],Te_high[ind], Zeff*1.3)
        W_up[ind]  = CalcRadiation(ne_high[ind],Te_high[ind],1e-4)



        from scipy.constants import e,m_p,m_u, mu_0


        m_i = Zeff*2
        m_z = 183.
        Z_i = 1.

        
        
        if self.dd.Open('CEZ',self.shot) or self.dd.Open('COZ',self.shot):
            
            if 'R_time' in self.dd.GetNames():
                R_cxrs = self.dd.GetAreabase('R_time')
                z_cxrs = self.dd.GetAreabase('z_time')
                ind = all(R_cxrs > 0,0)
            else:
                R_cxrs = self.dd.GetAreabase('R')
                z_cxrs = self.dd.GetAreabase('z')
                ind = R_cxrs > 0

                
            tvec_cxrs = self.dd.GetTimebase('time')

            R_cxrs = R_cxrs[...,ind]
            z_cxrs = z_cxrs[...,ind]
            Ti = self.dd.GetSignal('Ti_c')[:,ind]+1
            vtor = self.dd.GetSignal('vrot')[:,ind]
            self.dd.Close()

            self.eqm.read_ssq()

            R0 = interp(tvec_cxrs,self.eqm.t_eq, self.eqm.ssq['Rmag'])[:,None]
            
            rho_asym = self.eqm.rz2rho(R_cxrs,z_cxrs,tvec_cxrs)
       
            mach = sqrt(2*m_u/e*vtor**2/(2*Ti))
      
            Te_ = interp1d(tvec, Te,axis=1,bounds_error=False)(tvec_cxrs)
            Te_ = array([interp(r, rhop, te) for r, te in zip(rho_asym, Te_.T)])

            Te_[~isfinite(Te_)] = Ti[~isfinite(Te_)]
            Z_w  = CalcZ(Te_)


            M02 = (mach**2*m_z/m_i*(1-(m_i/m_z*Z_w*Zeff)*Te_/(Ti+Zeff*Te_)))/R_cxrs**2*R0**2
            #print('mach', nanmax(M02), sum(~isfinite(M02)), R0.min())

            asym_CF = exp((M02/R0**2)*(R0**2-R_cxrs**2))
            asym_CF = abs(1-asym_CF)/(asym_CF+1) #not consistent with tomography definition for large Mach! 


            asym=array([interp(rhop,r_[0,sort(r)],r_[0,a[argsort(r)]], right=0) for r,a in zip(rho_asym,asym_CF)])
                    
            asym = interp1d( tvec_cxrs, asym,axis=0,bounds_error=False)(tvec).T
 
            
        else:
            asym = zeros_like(BR )
            
            
                    
        savez_compressed(data_path,
                    tvec=single(tvec),rho=single(rhop),
                    BR=single(BR.T),
                    BR_low=single(BR_low.T),
                    BR_high=single(BR_up.T),
                    W=single(W.T),
                    W_low=single(W_low.T),
                    W_high=single(W_up.T),
                    asym=single(asym.T))
            
    
        
        
        
    
    def load_ICRH(self):
        
        icrh_path = self.local_path +"/geometry/ASDEX/ICRH_data/"
        if not os.path.exists(icrh_path):
            try:
                os.mkdir(icrh_path)
            except:
                print('error: os.mkdir(icrh_path)')
                raise
        ICRH_pos_file = icrh_path+"/ICRH_pos_"+str(self.shot)+".npz"

        try:
            assert config.useCache, 'do not use cache'
            ICRH_rezonance = load(ICRH_pos_file,allow_pickle=True,encoding='latin1')
            if 'tvec' in ICRH_rezonance:
                self.ICRH_rezonance = {k:v for k,v in ICRH_rezonance.items()}
        except:
            print('loading ICRH')
      
         
            try:
                assert not self.dd is None, 'No DD library'
                assert self.dd.Open('ICP', self.shot), 'No ICP shotfile'

                freqs = array([self.dd.GetParameter('Frequenz', 'Freq%s'%d) for d in range(1,5)])
                ICRH_tvec = self.dd.GetTimebase('T-B')
                
                PICRF = c_[[self.dd.GetSignal( "pnet%s"%d) for d in range(1,5)]]
                ICRH_on = PICRF > 1e5
                
                if self.shot > 30160:
                    self.dd.Open('MBI', self.shot)
                    mag_tvec = self.dd.GetTimebase('BTFABB')
                    Bt = self.dd.GetSignal('BTFABB')
                else:
                    self.dd.Open('MAI', self.shot)
                    mag_tvec = self.dd.GetTimebase('BTF')
                    Bt = self.dd.GetSignalCalibrated('BTF')
                self.dd.Close()

                Bt = interp(ICRH_tvec, mag_tvec[:len(Bt)-1],Bt[:-1])

                from scipy import constants
                R0 = 1.65
                BtR = abs(Bt)*R0 #T*m
                #He3 minority discharges!!!
                m_q_minor = 3/2. if self.shot in r_[28313:28321,28390,28391,31552:31556,31562:31564,31585,31587,31566,34694:34705 ] else 1
                Bc = (2*pi*freqs[None,:]*m_q_minor*constants.m_u/constants.e)
                R = BtR[:,None]/Bc

                R[~ICRH_on.T] = nan
                self.ICRH_rezonance = {'tvec':single(ICRH_tvec),'R': single(R)}


                savez_compressed(ICRH_pos_file , **self.ICRH_rezonance)
            except Exception as e:
                savez_compressed(ICRH_pos_file , None)
                print('no ICRH', e)
                
                
   

   
    def load_mag_equilibrium(self):
        
         #BUG!!!!  eqm do not support TRANSP CDF files with equlibrium
         
        eq_path = self.local_path +'/geometry/ASDEX/equilibrium/'

        if not os.path.exists(eq_path):
            try:
                os.mkdir(eq_path)
            except:
                print('error: os.mkdir(eq_path)')
                raise
        eq_path += 'MagField_fast_%d.npz'%self.shot
        
        from . import map_equ
        if not hasattr(self,'eqm'):
            self.eqm = map_equ.equ_map()   
            
     
        try:
            #zz
            data = load(eq_path,allow_pickle=True,encoding='latin1')
     
            self.tsurf = data['tsurf']
            self.surf_coeff = data['surf_coeff']
            #if 
            mag_diag = data['diag'].item()
            if isinstance(mag_diag,bytes): mag_diag = mag_diag.decode('utf-8')
            if mag_diag != self.mag_diag and not (self.mag_diag == 'EQI' and mag_diag == 'EQH'):
                warning('Warning: requested equilibrium differs from the stored: %s vs %s'%(mag_diag, self.mag_diag))

            self.mag_axis = {'tvec':data['tvec_fast'], 'Rmag':data['Rmag'], \
                'Zmag':data['Zmag'],'ahor':data['ahor'],'bver':data['bver']}
  
       
        except Exception as e:
            print('Magnetic equilibrium could not be loaded: '+str(e))
        
            
            diag = 'EQH' if self.mag_diag == 'TRA' else self.mag_diag 
            
            eq_diags = ((diag,self.mag_exp,self.mag_ed),('EQI','AUGD',0),('FPP','AUGD',0))
            
            stat = False
            
            for diag,mag_exp,mag_ed  in eq_diags:
                stat = self.eqm.Open(self.shot, diag=diag, exp=mag_exp, ed=mag_ed)
                if stat: break
                warning('Warning: equlibrium for shot:%d diag:%s  exp:%s  ed:%d  was not found!! other will be used'%(self.shot,diag,mag_exp,mag_ed))

            if not stat:
                raise Exception('equlibrium for shot:%d diag:%s  exp:%s  ed:%d  was not found!! use another one'%(self.shot,self.mag_diag,self.mag_exp,self.mag_ed))

            
            mag_diag = self.mag_diag
            from .mag_equ import  Equlibrium

            EQU = Equlibrium(self.eqm, self.shot, mag_diag,self.mag_exp,self.mag_ed)
            
            if mag_diag == 'TRA':
                output = EQU.getTranspEquilibrium()
            else:
                output = EQU.getStandartEquilibrium()
            output[mag_diag] = mag_diag
            savez_compressed(eq_path,diag=mag_diag,exp=self.mag_exp,ed=self.mag_ed, **output)

            
            self.tsurf = output['tsurf']
            self.surf_coeff = output['surf_coeff']
            self.mag_axis = output
            self.mag_axis['tvec'] = output['tvec_fast']




        self.mag_dt = amax(diff(self.mag_axis['tvec']))

        if mag_diag == 'TRA':
            self.mag_axis['Rmag'] = copy(self.surf_coeff[:,-1,0,0])
            self.mag_axis['Zmag'] = copy(self.surf_coeff[:,-1,0,1])
            self.surf_coeff[:,-1,0,:2] = 0

        
     
        

       
    def mag_equilibrium(self, tvec, preferCache = True, dryRun = False, return_mean = False,
                  surf_slice=slice(None,None),n_rho=60,n_theta=100,rho=None, radial_coordinate=None):
    
    
        """ Method to load magnetic field (contours). Try to load from local cache, if it is not possible try 
        to load from JET network and store it localy. Loaded field is only with limited precision (max 100 timeslices and it is interpolated)
         
        :param array tvec: Time vector  of requested  magnetic field
        :param bool dryRun: Do not return interpolated field
        :param bool preferCache: Use cached data if possible, it is faster

        :var array magx, magy: Arrays of magnetic lines possitions
        :var array tsurf: Time vector that belongs to ``magx``, ``magy``
        :raises Asseration: No connection to JET network"
        :rtype:  Coordinates of magnetic field interpolated to the requested time vector ``tvec``
        """
        
                
        def surf_polyval(rho,theta, p_rcos,p_zcos,p_rsin,p_zsin):
            #evaluate equilibrium polynom
            nmom = size(p_rcos,-1)
            rho = atleast_1d(rho)

            P = double(dstack((p_rcos,p_zcos,p_rsin,p_zsin)))
            moments = zeros((nmom, 4, size(rho)))
                
            #horner scheme
            for p in P:
                moments *= rho[None,None]
                moments += p[:,:,None]
            

            angle = outer(arange(nmom),theta )
            C = cos(angle)
            S = sin(angle)
            r_plot = tensordot(moments[:,0].T,C,axes=([-1,0])) #rcos
            r_plot+= tensordot(moments[:,2].T,S,axes=([-1,0])) #rsin
            z_plot = tensordot(moments[:,1].T,C,axes=([-1,0])) #zcos
            z_plot+= tensordot(moments[:,3].T,S,axes=([-1,0])) #zsin
        
            
            return r_plot.T, z_plot.T


        if rho is None:
            rho = linspace(0,1,n_rho)#[surf_slice]
        tvec = atleast_1d(tvec).astype('double')
        

      
        tvec  =  copy(tvec)


        tvec[tvec<self.tsurf[0]] = self.tsurf[0]
        tvec[tvec>self.tsurf[-1]] = self.tsurf[-1]

        
        
        ind  = slice(max(0,searchsorted( self.tsurf, tvec.min()-self.mag_dt*2)-1), 
                     min(len(self.tsurf),searchsorted( self.tsurf, tvec.max()+self.mag_dt*2)+1))

        from shared_modules import MovingAveradge

        theta = linspace(-pi,pi,n_theta)
        fast_tvec = self.mag_axis['tvec']

        #calculate separatrix
        if len(tvec) > 1:
            n_smooth = max(int(mean(diff(tvec))/mean(diff(fast_tvec))),1)
        else:
            n_smooth = sum((fast_tvec<tvec[-1]+self.mag_dt)&(fast_tvec>tvec[0]-self.mag_dt))

        surf_coeff = copy(self.surf_coeff[ind])
            
        tsurf = copy(self.tsurf[ind])

        surf_coeff = interp1d(tsurf, surf_coeff,axis=0,copy=False, bounds_error=False, fill_value=nan)(tvec)

        if return_mean:
            ind_fast = slice(fast_tvec.searchsorted(tvec[0]-self.mag_dt),
                             fast_tvec.searchsorted(tvec[-1]+self.mag_dt)+1)
            ahor = atleast_1d(median(self.mag_axis['ahor'][ind_fast]))
            bver = atleast_1d(median(self.mag_axis['bver'][ind_fast]))
            R0   = atleast_1d(median(self.mag_axis['Rmag'][ind_fast]))
            Z0   = atleast_1d(median(self.mag_axis['Zmag'][ind_fast]))
            
        else:
            fast_tvec = fast_tvec

            Dt = (tvec[-1]-tvec[0])/len(tvec)+(fast_tvec[-1]-fast_tvec[0])/fast_tvec.size
            tmin = max(tvec[0],fast_tvec[0 ])- Dt*(n_smooth+1)
            tmax = min(tvec[-1],fast_tvec[-1])+Dt*(n_smooth+1)
            imin = len(fast_tvec)-(-fast_tvec[::-1]).searchsorted(-tmin)
            imax = fast_tvec.searchsorted(tmax)
            
            
            
            ahor = MovingAveradge(double(copy(self.mag_axis['ahor'][imin:imax])),n_smooth)
            bver = MovingAveradge(double(copy(self.mag_axis['bver'][imin:imax])),n_smooth)
            R0   = MovingAveradge(double(copy(self.mag_axis['Rmag'][imin:imax])),n_smooth)
            Z0   = MovingAveradge(double(copy(self.mag_axis['Zmag'][imin:imax])),n_smooth)
            mag_tvec = fast_tvec[imin:imax]

            ahor = interp(tvec,mag_tvec, ahor)
            bver = interp(tvec,mag_tvec, bver)
            R0   = interp(tvec,mag_tvec,   R0)
            Z0   = interp(tvec,mag_tvec,   Z0)
                

             
        #renormalize almost circular  flux surfaces by fast equilibrium data
        surf_coeff[...,0::2] *= ahor[:,None,None,None]
        surf_coeff[...,1::2] *= bver[:,None,None,None]

        surf_coeff[:,-1,0,0] += R0
        surf_coeff[:,-1,0,1] += Z0
        
        
        
             
        #apply shift from the config file
        input_parameters = read_config(self.local_path+'/tomography.cfg')
        self.magfield_shift_core = input_parameters['magfield_shift_core']
        self.magfield_shift_lcfs = input_parameters['magfield_shift_lcfs']

        
        #shift only the center, not the separatrix            
        surf_coeff[:,-1,0,:2] += self.magfield_shift_core
        surf_coeff[:,-2,0,:2] -= self.magfield_shift_core
        surf_coeff[:,-2,0,:2] += self.magfield_shift_lcfs


   
        
        if  return_mean:      
            p_rcos,p_zcos,p_rsin,p_zsin = mean(surf_coeff, 0).T
            magx,magy = surf_polyval( rho ,theta, p_rcos.T,p_zcos.T,p_rsin.T,p_zsin.T)
            
     
        
        else:
            p_rcos,p_zcos,p_rsin,p_zsin = surf_coeff.T
            p_rcos,p_zcos,p_rsin,p_zsin = p_rcos.T,p_zcos.T,p_rsin.T,p_zsin.T
            
            magx = empty((size(theta), size(rho),size(tvec)),dtype='single')
            magy = empty((size(theta), size(rho),size(tvec)),dtype='single')
            

            for it, t in enumerate(tvec): #slowest
                magx[:,:,it],magy[:,:,it] = surf_polyval( rho ,theta, p_rcos[it],p_zcos[it],p_rsin[it],p_zsin[it])

   
        if any(isnan(magx)):
            raise Exception('nans in mag. surfaces!!')
        
        
        if radial_coordinate is None: radial_coordinate = self.radial_coordinate
        
        if radial_coordinate == 'r_a'  and size(rho)>1:
            self.convert_rho_2_r_a(rho, magx,magy)
        
        if radial_coordinate == 'r_V'  and size(rho)>1:
            self.convert_rho_2_r_V(rho, magx,magy)
        
        

        
        return rho,magx, magy




    def VesselStructures(self):
        try:
            from .map_equ import get_gc
            comp_r, comp_z = get_gc(self.shot)
            self.struct_dict = {c:array((comp_r[c], comp_z[c])) for c in list(comp_r.keys())}
        except:
            
            vessel_file = open('geometry/ASDEX/vessel.coords')

            struct_dict = {}
            self.struct_dict = {}
            struct_name = None
            for line in vessel_file:
                if len(line)==1:
                    continue
                elif line[0] == '#':
                    struct_name = line[1:].strip()
                    struct_dict[struct_name] = []
                else:
                    struct_dict[struct_name].append([float(n) for n in line.split()])
                            
            vessel_file.close()
            for k, items in struct_dict.items():    
                if len(items) > 0: 
                    self.struct_dict[k] = array(items).T
 
    
#https://www.aug.ipp.mpg.de/cgibin/local_or_sfread/cview.cgi?shot=29022
    def ShotOverview(self):
        return

        if not os.path.exists(self.local_path +'/geometry/ASDEX/overview_plots'):
            os.makedirs(self.local_path +'/geometry/ASDEX/overview_plots')
            
        name  = self.local_path +'/geometry/ASDEX/overview_plots/overview_%d.png'%self.shot

        if os.path.isfile(name):
            return 
        
        print('loading overview plot')
        url = 'https://www.aug.ipp.mpg.de/cgibin/local_or_sfread/cview.cgi?shot=%d'

        url_backup = 'https://www.aug.ipp.mpg.de/aug/local/aug_only/plotsigweb/plotsigweb.%d.png'
        u = 'dG9kc3RyY2k=\n'
        p = 'cGl0b21laGVzbG8=\n'

        import urllib.request, urllib.error, urllib.parse,base64
        
        
        

        try:
            password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
            password_mgr.add_password(None, "https://www.aug.ipp.mpg.de/cgibin/local_or_sfread/", u.decode('base64'), p.decode('base64'))

            handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
            opener = urllib.request.build_opener(handler)

            timeout  = 120
            try:
                imgData = opener.open(url%self.shot, timeout=timeout).read()
            except:
                imgData = opener.open(url_backup%self.shot, timeout=timeout).read()

            output = open(name,'wb')
            output.write(imgData)
            output.close()
        except Exception as e:
            print(e, 'ShotOverview plot was not loaded')
            pass
        
