 #!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
import sys,os
from matplotlib.pyplot import *
from tokamak import Tokamak
from collections import OrderedDict
import config
#AUTOR: Tomas Odstrcil  tomas.odstrcil@ipp.mpg.de
from shutil import copyfile
from scipy.interpolate import interp1d
from shared_modules import read_config, warning, debug

import time 
#from matplotlib.pylab import *


global loader
global loader2

#BUG xtomo ignuruje ten čas začátku a konce!!
 
class DIIID(Tokamak):
    """
    Representation of class of tokamak. Contain basic information about 
geometry, 
    borde, magnetic field and basic setting used for recontruction and 
postprocessing.

    >>> DIIID(input_path, 'SXR',  shot, nx, ny,virt_chord )


    :var str name:  Name of tokamak
    :var int index: Index of the tokamak
    :var double xmin, xmax, ymin, ymax: boundarys of recontructed area of 
tokamak emissivity
    :var str geometry_path: Path to saved geometry, boundary or data
    :var double sigma:  Expected noise in data  += ``sigma``* sqrt(data) * 
sqrt(max(data))  (Poisson like noise)
    :var double min_error: Expected noise in data  += ``min_error`` * 
sqrt(max(data))
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
    :var array magx, magy: Arrays of magnetic lines possitions
    :var spmatrix Tmat: Matrix of geometry
    :var array Xchords, Ychords: Arrays of chords coordinates, used in 
``make_plots``
    :var array boundary: Coordinates of points on boundary of tokamak
    """


    def __init__(self, input_diagn, input_parameters,  only_prepare = False,load_data_only=False  ):
        debug('__ninitDIIID')
        """ Initialize new tokamak and set default values.

            :param dict input_parameters:   List of selected paramerers
            :param str input_diagn:  Name of main diagnostics
            :param array magfield_shift: artificial shift of magnetic field

        """
 
 
        name  = 'DIIID'

        self.mds_server = input_parameters['mds_server']
        self.rad_profile = 'Gaussian'
        self.mag_diag = input_parameters['mag_diag']
        self.mag_exp  = input_parameters['mag_exp']
        self.mag_ed   = input_parameters['mag_ed']
        self.input_diagn = input_diagn
        fast_data = False
        if input_diagn.endswith('fast'):
            fast_data = True
            input_diagn = input_diagn[:-5]


        path = 'geometry/'+name+'/'+ input_diagn

        self.local_path =  input_parameters['local_path']
        self.program_path =  input_parameters['program_path']
        self.geometry_path_program = self.program_path+path
        path = self.local_path +path

        if not os.path.exists(path):
             os.makedirs(path)

        print(input_diagn)

        
        self.shot = input_parameters['shot']

        coord= [0.95,2.45, -1.38, 1.38]

        sigma = 0.001
        min_error = 0.001
        norm = 1.0
        t_name = 's'
        
        
        self.min_tvec = 0.
        self.max_tvec = 10

        self.wrong_dets_damaged  = []


 
        calib_path = path+'/calibration/'
        if not os.path.exists(calib_path):  os.mkdir(calib_path)
            
        wch_path = path+'/wrong_channels/'
        if not os.path.exists(wch_path):  os.mkdir(wch_path)
            
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
        
        if os.path.exists(self.geometry_path_program+'/border.txt'):
            if self.geometry_path_program+'/border.txt'!= path+'/border.txt':
                copyfile(self.geometry_path_program+'/border.txt',path+'/border.txt')
        if os.path.exists(self.geometry_path_program+'/border_sas.txt'):
            if self.geometry_path_program+'/border_sas.txt'!= path+'/border_sas.txt':
                copyfile(self.geometry_path_program+'/border_sas.txt',path+'/border_sas.txt')
        
        
        border_file = 'border_sas.txt' if self.shot > 168847 else 'border.txt' 
        self.vessel_boundary = loadtxt(path+'/'+border_file)

        Tokamak.__init__(self, input_parameters, input_diagn, coord,
                         name, norm, sigma, min_error, path,t_name)




        if input_diagn == 'SXR':
            self.SXR(fast_data,toroidal=False)
            self.index = 21+fast_data
        if input_diagn == 'BOLO':
            self.BOLO()   
            self.index = 23
        if input_diagn == 'DISRAD':
            self.DISRAD()
            self.index = 24

        self.VesselStructures()
        
        self.load_geom()

        self.allow_negative = False


        Tokamak._post_init(self)
        self.Tokamak = Tokamak
        self.input_diagn =input_diagn 

 
        if not only_prepare:
            self.prepare_tokamak()
        
        if load_data_only:
            return
        
        if input_diagn == 'SXR':
            self.load_others()
        self.load_LBO()
        self.load_mag_equilibrium()

        
        
    def SXR(self, fast_data=False,toroidal=False):
        
        
        self.default_calb = 'variable'   
        self.allow_self_calibration = False
        self.sigma = 0.01 #just guess!!
        self.min_error = 0.005#just guess!!
        self.allow_negative = True

        from .SXR.load_SXR import loader_SXR
        try:
            import MDSplus as mds
            c = mds.Connection(self.mds_server )
        except:
            c = None
            print('MDS connection with server "%s" was not established'%self.mds_server) 
        global loader
        loader = loader_SXR(self.shot,self.geometry_path,c, fast_data, toroidal)
        
        from copy import deepcopy
        self.detectors_dict = deepcopy(loader.detectors_dict)

        self.calb_0 = loader.calb_0
        self.nl =  loader.nl
        
        self.dets_index = loader.dets_index
        self.dets = loader.dets
        self.tvec = loader.tvec
        self.Phi = [loader.Phi[k] for k in list(self.detectors_dict.keys())]
    
        if self.input_diagn in ['SXR', 'SXR fast'] :
            global loader2
            try:
                loader2 = loader_SXR(self.shot,self.geometry_path,c, fast_data,~toroidal)
                for k, i in loader2.detectors_dict.items():
                    self.detectors_dict[k] = i
                
                self.calb_0 = r_[self.calb_0,loader2.calb_0]
                
                self.dets_index+=[ind+self.nl for ind in  loader2.dets_index]
                self.Phi += [loader2.Phi[k] for k in list(loader2.detectors_dict.keys())]
              
                self.dets = r_[ self.dets, self.nl+loader2.dets]
                self.tvec2 = loader2.tvec
                self.nl +=  loader2.nl
            except:
                loader2 = None
            

        self.sample_freq = len(self.tvec)/(self.tvec[-1]-self.tvec[0])
        
        if not fast_data: 
            self.sample_freq = 1e3  #data are downsampled
            
        self.Phi = deg2rad(hstack(self.Phi))

    def BOLO(self):
        
        self.default_calb = 'variable'   
        self.allow_self_calibration = False
        self.sigma = 0.005 #just guess!!
        self.min_error = 0.005#just guess!!
        global loader

        from .BOLO.BOLO import loader_BOLO
        
        import MDSplus as mds
        c = mds.Connection(self.mds_server )

        loader = loader_BOLO(self.shot,self.geometry_path,c)
        
        self.detectors_dict= loader.detectors_dict

        self.calb_0 = loader.calb_0
        self.nl =  loader.nl
        
        self.dets_index = loader.dets_index
        self.dets = loader.dets
        self.tvec = loader.tvec
        
        self.sample_freq = len(self.tvec)/(self.tvec[-1]-self.tvec[0])
    
        
        self.total_rad = loader.get_total_rad()
        self.use_pfm = True


    def init_equ(self):
        from . import map_equ
        eqm = map_equ.equ_map(None)
        eqm.pfm   = self.PFM['pfm']
        eqm.t_eq  = self.PFM['pfm_tvec']
        eqm.Rmesh = self.PFM['pfm_R']
        eqm.Zmesh = self.PFM['pfm_Z']
        eqm.eq_open=True  
        return eqm

    def get_psi(self, tvec,xgridc, ygridc):
        eqm = self.init_equ()
        R,Z = meshgrid(xgridc,ygridc)
        Psi = eqm.rz2rho(R[None],Z[None],tvec,coord_out='Psi')  
      
        return squeeze(sqrt(maximum(0,Psi)))
    
    def get_mag_contours(self,tvec,rho):
        eqm = self.init_equ()
        eqm.psi0 = eqm.t_eq*nan
        eqm.eq_open=True        
        
        return eqm.rho2rz(rho**2, tvec,'Psi',True)  
            
  
                
    def load_geom(self):
 
        global loader

        loader.load_geom(self.geometry_path_program)
        if self.input_diagn in ['SXR', 'SXR fast'] :
            global loader2
            if loader2 is not None:
                loader2.load_geom(self.geometry_path_program)

  
            
    def get_data(self, failsave = True,tmin = -infty,tmax = infty):

        t = time.time()
        global loader

    
        tvec, data, error = loader.get_data(tmin,tmax)
        wrong_dets_damaged = where(in1d(loader.all_los,loader.wrong_dets_damaged))[0]
      

        if self.input_diagn in ['SXR', 'SXR fast']:
            global loader2
            if loader2 is not None:
                tvec2, data2, error2 = loader2.get_data(tmin,tmax)
                nt = min(tvec.size, tvec2.size)
                data = hstack((data[:nt], data2[:nt]))
                error = hstack((error[:nt],error2[:nt]))
                tvec = tvec[:nt]
                #BUG funguje to i pro rychla data?

                wrong_dets_damaged2 = where(in1d(loader2.all_los,loader2.wrong_dets_damaged))[0]
                wrong_dets_damaged = r_[wrong_dets_damaged,wrong_dets_damaged2+loader.nl]
             
                O = zeros((loader.nl,loader2.nl))
              
            
        self.wrong_dets_damaged  = unique(r_[self.wrong_dets_damaged,wrong_dets_damaged])
       
        ind_t = slice(None,None)
        if  len(tvec) > 5000:  ind_t = random.randint(len(tvec),size=1000)
        
        mean_data = nanmean(data[ind_t],axis=0)

       #include also "nonstatistical noise"  by min_error and sigma
        error =  hypot(error, single(mean_data*self.sigma),out=error)
        error += abs(mean_data)*self.min_error

        print('loading time: %.1f s'%( time.time()-t))

        return  tvec, data, error 
         
        
        
        
        
    def load_mag_equilibrium(self):
        
         #BUG!!!!  eqm do not support TRANSP CDF files with equlibrium
         
        eq_path = self.local_path +'/geometry/DIIID/equilibrium/'

        if not os.path.exists(eq_path):
            os.mkdir(eq_path)

        eq_path += 'MagField_fast_%d.npz'%self.shot

        try:
            data = load(eq_path,allow_pickle=True)
            self.tsurf = data['tsurf']
            self.surf_coeff = data['surf_coeff']
            mag_diag = data['diag'].item()
            if mag_diag != self.mag_diag and not (self.mag_diag == 'EQI' and mag_diag == 'EQH'):
                warning('Warning: requested equilibrium differs from the stored: %s vs %s'%(mag_diag, self.mag_diag))

            self.mag_axis = {'tvec':data['tvec_fast'], 'Rmag':data['Rmag'], \
            'Zmag':data['Zmag'],'ahor':data['ahor'],'bver':data['bver']}
            
            if 'pfm' in data:
                self.PFM = {'pfm':data['pfm'],'pfm_tvec':data['pfm_tvec'],'pfm_R':data['pfm_R'],'pfm_Z':data['pfm_Z']}

        
        except Exception as e:
            print('eqerr', e)
            import MDSplus as mds
            MDSconn = mds.Connection(self.mds_server )

            from . import map_equ
            if not hasattr(self,'eqm'):
                self.eqm = map_equ.equ_map(MDSconn)   
            

            
            eq_diags = ((self.mag_diag,self.mag_exp,self.mag_ed),('EFIT01','DIIID',0),('EFITRT1','DIIID',0) )
            
            stat = False
            
            for diag,mag_exp,mag_ed  in eq_diags:
                stat = self.eqm.Open(self.shot, diag=diag, exp=mag_exp, ed=mag_ed)
                #print size(self.eqm.t_eq)
                if stat and size(self.eqm.t_eq) >2 : break
                warning('Warning: equlibrium for shot:%d diag:%s  exp:%s  ed:%d  was not found!! other will be used'%(self.shot,diag,mag_exp,mag_ed))

            if not stat:
                raise Exception('equlibrium for shot:%d diag:%s  exp:%s  ed:%d  was not found!! use another one'%(self.shot,self.mag_diag,self.mag_exp,self.mag_ed))


            mag_diag = self.mag_diag
            from .mag_equ import  Equlibrium

            EQU = Equlibrium(MDSconn,self.eqm, self.shot, mag_diag,self.mag_exp,self.mag_ed)
            
            if 'TRA' in mag_diag:
                output = EQU.getTranspEquilibrium()
            else:
                output = EQU.getStandartEquilibrium()
            output[mag_diag] = mag_diag

            
            self.tsurf = output['tsurf']
            self.surf_coeff = output['surf_coeff']
            self.mag_axis = output
            self.mag_axis['tvec'] = output['tvec_fast']
          
            try:
                    
                pfm = copy(self.eqm.pfm)
                pfm-= self.eqm.psi0
                pfm/= (self.eqm.psix-self.eqm.psi0)
                pfm_tvec = self.eqm.t_eq
                pfm_R = self.eqm.Rmesh
                pfm_Z = self.eqm.Zmesh
                self.PFM = {'pfm':pfm,'pfm_tvec':pfm_tvec,'pfm_R':pfm_R,'pfm_Z':pfm_Z}
                output.update(self.PFM)
            except:
                pass
                
                
            savez_compressed(eq_path,diag=mag_diag,exp=self.mag_exp,
                             ed=self.mag_ed,**output)


            MDSconn.closeAllTrees()


        self.mag_dt = amax(diff(self.mag_axis['tvec']))
        self.mag_axis['Rmag'] = copy(self.surf_coeff[:,-1,0,0])
        self.mag_axis['Zmag'] = copy(self.surf_coeff[:,-1,0,1])
        self.surf_coeff[:,-1,0,:2] = 0
        if size(self.mag_axis['ahor']) == 1 and (self.mag_axis['ahor'] is None or self.mag_axis['ahor'].item() is None):
            self.mag_axis['ahor'] = ones_like(self.mag_axis['tvec'])
            self.mag_axis['bver'] = ones_like(self.mag_axis['tvec'])

    def VesselStructures(self):
        self.struct_dict = {}
        
        from .map_equ import get_gc
        comp_r,comp_z = get_gc(self.shot)
        self.struct_dict['inner boundary'] = comp_r['vessel'],comp_z['vessel']
        
        


       
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
            
            #print self.mag_axis['ahor']
            ahor = atleast_1d(median(self.mag_axis['ahor'][ind_fast]))
            bver = atleast_1d(median(self.mag_axis['bver'][ind_fast]))
            R0   = atleast_1d(median(self.mag_axis['Rmag'][ind_fast]))
            Z0   = atleast_1d(median(self.mag_axis['Zmag'][ind_fast]))
            
        else :
            fast_tvec = fast_tvec

            Dt = (tvec[-1]-tvec[0])/len(tvec)+(fast_tvec[-1]-fast_tvec[0])/fast_tvec.size
            tmin = max(tvec[0],fast_tvec[0 ])#- Dt*(n_smooth+1)
            tmax = min(tvec[-1],fast_tvec[-1])#+Dt*(n_smooth+1)
            imin = len(fast_tvec)-(-fast_tvec[::-1]).searchsorted(-tmin)
            imax = fast_tvec.searchsorted(tmax)
            imin = max(0, imin-n_smooth-1 )
            imax = min(len(fast_tvec), imax+n_smooth+1 ) #BUG can couse problems at the begining and end

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
      
        if hasattr( config, 'magfield_shift') and  config.magfield_shift != (0,0):
            self.magfield_shift_core = config.magfield_shift
            print('config shift')
        else:
            self.magfield_shift_core = input_parameters['magfield_shift_core']
    
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
        
        

        
    def bremsstrahlung_power(self,n_e, T_e,Zeff):
        '''
        Calculate the radiation power density from bremsstrahlung radiation.
        This calcuation uses the popular formula in cgs/eV/A units giving power in W/cm**3
        or equivalently MW/m**3
        power = [1.89e-28 * gff * n_e**2 * Zeff / (sqrt(T_e)*lambda**2)] * exp[-hc/(T_e*lambda)]
        here gff = 1.35 T_e**0.15 and n_e in (cm**-3), T_e in (eV) and lambda in (A).
        for these units hc=1.2400e4
        Total bremsstrahlung (NRL Formulary) power in W/cm**3 is
        power = 1.69e-32 ne**2 sqrt(Te) Zeff

        inputs
        (none) : Returns total bremsstrahlung radiation.
        wavelength: scalar or array where scalar used for a VB filterscope system
            and array used for SXR or other system with wavelength dependent response.
            Must be ascending in wavelength.
        transmission: Filter transmisson function of (ascending) wavelength.
        outputs
        Radiated power in (time,space) for total bremsstrahlung power (W/cm**3 or MW/m**3)
        visible bremsstrahlung at wavelength (W/cm**3/A) or filtered,
        wavelength integrated emissivity for SXR (W/cm**3)
        '''
        
        h = 6.6261e-34
        c = 2.9979e8
        q = 1.6022e-19
        
        filter_dict = loadtxt(self.geometry_path+'/radiation/adas414_adf35_Be125_Si300.dat')
        # Energy in eV
        energy = filter_dict[::-1,0]
        # Filter response on energy(ascending) and wavelength(descending) grid
        transmission = filter_dict[::-1,1]
        # Wavelength in A
        wavelength = h*c/(energy*q)*1.e10
        from scipy import  integrate
        T_e = maximum(T_e,1)
        # Free-free gaunt factor
        gff = 1.35 * T_e**0.15
        # hc for the quantity [hc/(Te*lambda)] for T_e in (eV) and lambda in (A)
        hc = 1.24e4

        coeff = 1.89e-28 * (n_e*1e-6)**2 * Zeff * gff / np.sqrt(T_e)

        vb_wl = outer(coeff , transmission/wavelength**2) * exp(-hc / outer(T_e,wavelength))
        vb_wl_int = integrate.simps(vb_wl, x=wavelength, axis=1)
        
        return vb_wl_int*1e6
        
        
    def load_LBO(self):

        lbo_path = self.local_path +"/geometry/DIIID/LBO/"
        
        if not os.path.exists(lbo_path):
            os.mkdir(lbo_path)

        lbo_path += "LBO_"+str(self.shot)+".txt"
    
        
        try:
            self.impur_inject_t = atleast_1d(loadtxt(lbo_path))             
            if all(isnan(self.impur_inject_t)):
                self.impur_inject_t = None
        except Exception as e:
            
            
            try:
   
    
    
                import MDSplus as mds
                MDSconn = mds.Connection(self.mds_server )
                            
                #LBOSYNC  = MDSconn.get('PTDATA("LBOSYNC", %d) '%self.shot).data()
                LBOSHUTT = MDSconn.get('_x = PTDATA("LBOSHUTT", %d)'%self.shot).data()
                LBOQSWCH = MDSconn.get('_x = PTDATA("LBOQSWCH", %d)'%self.shot).data()
                assert len(LBOSHUTT) != 1, 'no lbo data'

                tvec = MDSconn.get('dim_of(_x)').data()
                
                LBOQSWCH  = LBOQSWCH >= 1
                LBOSHUTT  = LBOSHUTT >= 1 #BUG 
                    
                assert any(LBOQSWCH&LBOSHUTT), 'no LBO'
                
                lbo_times = tvec[where((diff(double(LBOQSWCH&LBOSHUTT)) > 0 ))]/1e3

                savetxt(lbo_path, lbo_times,fmt='%3.5f' )
                self.impur_inject_t = lbo_times
                print('LBO times: '+ str(lbo_times))
          
                    
            except:
                #raise
                self.impur_inject_t = None
                savetxt(lbo_path, (nan,))

            

    def load_others(self):
        #load centrifugal asymmetry profile and bremsstrahlung level 
        loaded = False
        
        data_path = self.geometry_path+'/Bremsstrahlung_%d.npz'%self.shot

        if os.path.isfile(data_path) and config.useCache:
            return 
            #pass
        print('Loading radiation data ')

        try:
            
            kin_data = load( self.geometry_path+'/kin_data_%d.npz'%self.shot)
            
            
            ne = kin_data['ne'].item()
            ne_data = ne['data']
            ne_rho  = ne['rho']
            ne_tvec = ne['tvec']
            ne_err  = ne['err']
            
            
            Ti = kin_data['Ti'].item()
            Ti_data = Ti['data']
            Ti_rho  = Ti['rho']
            Ti_tvec = Ti['tvec']
            Ti_err  = Ti['err']
            
            Te = kin_data['Te'].item()
            Te_data = Te['data']
            Te_rho  = Te['rho']
            Te_tvec = Te['tvec']
            Te_err  = Te['err']
            
            omega =  kin_data['omega'].item()
            omega_data = omega['data']
            omega_rho  = omega['rho']
            omega_tvec = omega['tvec']
            omega_err  = omega['err']
            
            if 'Zeff' in kin_data:
                Zeff = kin_data['Zeff'].item()
                Zeff_data = Zeff['data']
                Zeff_rho  = Zeff['rho']
                Zeff_tvec = Zeff['tvec']
                Zeff_err  = Zeff['err']
            else:
                nimp = kin_data['nimp'].item()
                nimp_data = nimp['data']
                Zeff_rho  = nimp['rho']
                Zeff_tvec = nimp['tvec']
                nimp_err  = nimp['err']
                Zimp = 6  #carbon
                Zmain = 1  #carbon
                Zeff_data  = (Zimp**2*nimp - Zimp*Zmain*nimp + Zmain*ne)/ne

            print('Loaded from kindata')



            
            loaded = True
            print('kinfit')
        except:
            pass
   
        
 
        if not loaded:
            try:
                #assert not loaded, 'already loaded' 

                
                try:
                    import MDSplus as mds
                    MDSconn = mds.Connection(self.mds_server )
                except:
                    
                    #c = None
                    print('MDS connection with server "%s" was not established'%self.mds_server) 
                    return
                

                
                
                MDSconn.openTree('IONS', self.shot)
                Ti_data = MDSconn.get('_x=\\IONS::TOP.PROFILE_FITS.ZIPFIT.ITEMPFIT').data()*1e3
                Ti_rho  = MDSconn.get('dim_of(_x,0)').data()
                Ti_tvec = MDSconn.get('dim_of(_x,1)').data()/1e3
                Ti_err  = abs(MDSconn.get('error_of(_x)').data())*1e3
                
                omega_data = MDSconn.get('_x=\\IONS::TOP.PROFILE_FITS.ZIPFIT.TROTFIT').data()*1e3
                omega_rho  = MDSconn.get('dim_of(_x,0)').data()
                omega_tvec = MDSconn.get('dim_of(_x,1)').data()/1e3
                omega_err  = abs(MDSconn.get('error_of(_x)').data())*1e3

                nimp = MDSconn.get('_x=\\IONS::TOP.PROFILE_FITS.ZIPFIT.ZDENSFIT').data()*1e19
                Zeff_rho  = MDSconn.get('dim_of(_x,0)').data()
                Zeff_tvec = MDSconn.get('dim_of(_x,1)').data()/1e3
                nimp_err  = abs(MDSconn.get('error_of(_x)').data())*1e19
                Zimp = 6  #carbon
                Zmain = 1  #carbon

                MDSconn.closeTree('IONS', self.shot)
                
                MDSconn.openTree('ELECTRONS', self.shot)

                ne_data = MDSconn.get('_x=\\ELECTRONS::TOP.PROFILE_FITS.ZIPFIT.EDENSFIT').data()*1e19
                ne_rho  = MDSconn.get('dim_of(_x,0)').data()
                ne_tvec = MDSconn.get('dim_of(_x,1)').data()/1e3
                ne_err  = abs(MDSconn.get('error_of(_x)').data())*1e19

                Te_data = MDSconn.get('_x=\\ELECTRONS::TOP.PROFILE_FITS.ZIPFIT.ETEMPFIT').data()*1e3
                Te_rho  = MDSconn.get('dim_of(_x,0)').data()
                Te_tvec = MDSconn.get('dim_of(_x,1)').data()/1e3
                Te_err  = abs(MDSconn.get('error_of(_x)').data())*1e3
                
                ne_ = interp1d(ne_tvec, ne_data, axis=0, bounds_error=False)(Zeff_tvec)
                Zeff_data  = (Zimp**2*nimp - Zimp*Zmain*nimp + Zmain*ne_)/ne_

                MDSconn.closeTree('ELECTRONS', self.shot)
            except Exception as e:
                print('loading of ZIPFIT failed: '+str(e)) 
                
                return 
                
            
        tvec = unique(r_[Te_tvec,ne_tvec,Zeff_tvec,omega_tvec,Ti_tvec])
 
        rind = ne_rho<=1
        
      

        
        rhot = ne_rho[rind]
        Ti = interp1d(Ti_tvec, Ti_data[:,rind], axis=0, bounds_error=False)(tvec)
        Ti_err = interp1d(Ti_tvec, Ti_err[:,rind], axis=0, bounds_error=False)(tvec)

        omega = interp1d(omega_tvec, omega_data[:,rind], axis=0, bounds_error=False)(tvec)
        omega_err = interp1d(omega_tvec, omega_err[:,rind], axis=0, bounds_error=False)(tvec)

        Zeff = interp1d(Zeff_tvec, Zeff_data[:,rind], axis=0, bounds_error=False,fill_value="extrapolate")(tvec)
        
        Te = interp1d(Te_tvec, Te_data[:,rind], axis=0, bounds_error=False)(tvec)
        Te_err = interp1d(Te_tvec, Te_err[:,rind], axis=0, bounds_error=False)(tvec)

        ne = interp1d(ne_tvec, ne_data[:,rind], axis=0, bounds_error=False)(tvec)
        ne_err = interp1d(ne_tvec, ne_err[:,rind], axis=0, bounds_error=False)(tvec)



        sys.path.append( self.geometry_path)
        from .SXR.radiation.Bremsstrahlung import CalcBackgroundRadiation,CalcRadiation,CalcZ, LoadAtomData,LoadFile
        
        
        BR = ones_like(Te)*nan
        BR_low = ones_like(Te)*nan
        BR_up = ones_like(Te)*nan
        W = ones_like(Te)*nan
        W_low = ones_like(Te)*nan
        W_up = ones_like(Te)*nan
            
        ind = isfinite(Te)&isfinite(ne)
        ind2 = isfinite(Te)&isfinite(ne)&isfinite(Zeff)

        filt = 14
        BR[ind2 ] = CalcBackgroundRadiation(ne[ind2],Te[ind2], Zeff[ind2],filt=filt)
        BR_low[ind2]= CalcBackgroundRadiation(maximum(1,ne-ne_err)[ind2],maximum(1,Te-Te_err)[ind2], Zeff[ind2]/1.1,filt=filt)
        BR_up[ind2] = CalcBackgroundRadiation((ne+ne_err)[ind2],(Te-Te_err)[ind2], Zeff[ind2]*1.3,filt=filt)
    
    
        W[ind]  = CalcRadiation(ne[ind],Te[ind],1e-4,filt=filt)
        W_low[ind]  = CalcRadiation(maximum(1,ne-ne_err)[ind],maximum(1,Te-Te_err)[ind],1e-4,filt=filt)
        W_up[ind]   = CalcRadiation((ne+ne_err)[ind],(Te+Te_err)[ind],1e-4,filt=filt)


        from scipy.constants import e,m_p,m_u, mu_0


        m_i = 2
        m_z = 183.
        Z_i = 1.


        R0 = 1.7
        a0 = 0.5  
        R_cxrs = R0+a0*rhot
        vtor = omega*R_cxrs#/(2*pi)   
        mach = sqrt(abs(2*m_u/e*vtor**2/(2*Ti)))

        Z_w  = CalcZ(Te)


        M02 = (mach**2*m_z/m_i*(1-(m_i/m_z*Z_w*Zeff)*Te/(Ti+Zeff*Te)))/R_cxrs**2*R0**2
        
        prof_LFS = exp((M02/R0**2)*(R_cxrs**2-R_cxrs**2)) # = 1
        prof_HFS = exp((M02/R0**2)*((R0-a0*rhot)**2-R_cxrs**2))
        asym = abs(prof_LFS-prof_HFS)/(prof_HFS+prof_LFS) #not consistent with tomography definition for large Mach! 
  
        savez_compressed(data_path,
                    tvec=single(tvec),
                    rho=single(rhot), #rho toroidal 
                    BR=single(BR),
                    BR_low=single(BR_low),
                    BR_high=single(BR_up),
                    W=single(W),
                    W_low=single(W_low),
                    W_high=single(W_up),
                    asym=single(asym))
         
