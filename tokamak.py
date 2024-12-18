#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,sys
from numpy import *
import matplotlib
from scipy import sparse
from shared_modules import in1d, debug,MovingAveradge,smooth
from orthogonal_trans import *
import time
from scipy.stats.mstats import mquantiles
from scipy.signal import medfilt
#from copy import deepcopy
import config
from annulus import get_bd_mat
import gc
from scipy.interpolate import splprep, splev
from scipy.interpolate import interp1d

#AUTOR: Tomas Odstrcil  tomas.odstrcil@ipp.mpg.de



class Tokamak(object):
    """
    Representation of class of tokamak. Contain basic information about geometry, 
    borde, magnetic field and basic setting used for recontruction and postprocessing.

    >>> JET(input_path, 'SXR-slow',  shot, nx, ny,virt_chord )


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
    :var array magx, magy: Arrays of magnetic lines possitions
    :var spmatrix Tmat: Matrix of geometry
    :var array Xchords, Ychords: Arrays of chords coordinates, used in ``make_plots``
    :var array BdMat: Matrix that is nonzero inside the tokamak
    :var array boundary: Coordinates of points on boundary of tokamak
    :var array emiss_0: First guess of plasma emissivity, improve convergence
    """


        
    def __init__(self, input_parameters,  input_diagn, coord,  name, norm, sigma, min_error, path,t_name):
        """ Initialize new tokamak and set default values.

            :param str input_diagn:  Name of main diagnostics
            :param int shot:
            :param int nx:
            :param int ny:
            :param int virt_chord:
            :param array magfield_shift: artificial shift of magnetic field

        """

        for a in ['shot',  'nx', 'ny', 'virt_chord',  'magfield_shift_core',
                  'magfield_shift_lcfs','transform_order_a', 'radial_coordinate',
                  'transform_order_r',  'boundary_width', 'transform_index',
                'prune_dets',    
                "boundary",'cos_com','sin_com','use_median_filter' ]:
            exec( "self." + a+" = input_parameters['"+a+"']")
        print("=========  TOKAMAK: "+ name+ " ====================")


        # BASIC PARAMETERS
        
        self.name = name
        self.xmin = coord[0]  #NOTE left side of most left pixel!!!
        self.xmax = coord[1]  #NOTE right side of right most pixel!!!
        self.ymin = coord[2]
        self.ymax = coord[3]
        self.plot_coord = coord # default expectation => coord contains whole vessel 
        self.geometry_path = os.path.normpath(path)+os.sep
        self.tmp_folder = input_parameters['tmp_folder']
 
        self.allow_self_calibration = False  #  use special recontruction to obtain optimal ratios of detectords

        self.wrong_dets_damaged = []
        self.wrong_dets = []

        self.rad_profile = 'Emiss' 
        self.set_resolution(self.nx, self.ny)
        self.sigma = sigma
        self.min_error = min_error
        self.norm = norm
        self.t_name = 's'   # time in seconds
        self.nl = NaN   # default number of chord is unknown
        self.default_calb  = 'const'  # const, none, variable
        self.camera = False  # in the case of tangential  camera => True
        self.allow_negative = True  # negative values are most often unphysical but not always 
        self.beta = 0
        self.use_pfm = False  #  use manetic equlibrium directly from poloidal fux matrix

        if self.radial_coordinate == 'r_a':
            self.rho_label = '$r/a$'
        elif self.radial_coordinate == 'rho_pol':
            self.rho_label = r'$\rho_\theta$'
        elif self.radial_coordinate == 'r_V':
            self.rho_label = r'$\rho_V$'
        else:
            raise Exception('Coordinate is not implemented')
        
        
        
        self.expect_zero = False
        self.allow_reflections = False # use modeling of background reflections for 2D cameras 
        self.impur_inject_t = None   #times of laser blow off or other kinds of impurity innjections
        self.shift = 10    #initial guess of the MFI method    (non necessery for linear-like methods)
        self.geometry_version=None #if more different version of geometry was avalible (due to diagnostics changes etc..)
        self.boundary_lcfs = True# the boundary is given by the last closed flux surface
        if not hasattr(self,'vessel_boundary'):
            print('loading vessel')
            self.vessel_boundary = loadtxt(self.geometry_path+'border.txt')
        
        if all(self.vessel_boundary[0,:] != self.vessel_boundary[-1,:]):   # first and last point are different            
            self.vessel_boundary = r_[self.vessel_boundary, self.vessel_boundary[(0,),:]]
            
        
        
        # transformations =========================================================

        self.wide_boundary = False
        self.boundaryWidth = 10


    def _post_init(self):
        """ 
        Commands that needs to be run after the data loading 
        """
        self.Ndets = len(self.dets_index)   #  number of standalone detectors 
        
        if os.path.isfile(os.path.join(self.geometry_path,'wrong_channels',str(self.shot))):
            with open(os.path.join(self.geometry_path,'wrong_channels',str(self.shot)),'r') as file:
                wrong_dets_stored = [int(a) for a in file]
                config.wrong_dets_pref  = unique(r_[wrong_dets_stored,config.wrong_dets_pref])
 

        
    def set_resolution(self,nx, ny):
        self.nx = nx
        self.ny = ny
        self.npix = nx*ny
        
        #width of the pixels
        self.dx =(self.xmax-self.xmin)/double(self.nx)
        self.dy =(self.ymax-self.ymin)/double(self.ny)

        #NOTE coordinates of the lower left corner of the pixel! 
        self.xgrid = linspace(self.xmin,self.xmax,self.nx,endpoint=False)
        self.ygrid = linspace(self.ymin,self.ymax,self.ny,endpoint=False)

        
        self.area_axis = r_[self.xmin,self.xmax, self.ymin, self.ymax]
        # axis for geometry generation (camera)
        self.geometry_axis = r_[self.xmin,self.xmax, self.ymin, self.ymax, -self.xmax, self.xmax]


    def generate_transform(self, tvec):


        tvec = asarray(tvec.mean())
        debug('Transformation : ' + str(self.transform_index) ) 
        params = [self,tvec, self.transform_order_r,self.transform_order_a,
                         self.cos_com, self.sin_com]
        
        if   self.transform_index == 0: # None
            Transform = sparse.eye(self.npix,self.npix)
        elif self.transform_index == 1: # Abel
            Transform= generalized_abel_transform(*params)
            self.wide_boundary = True
            
        elif self.transform_index == 2: # Cormack
            Transform = global_basis_transform('zernike',*params)
            self.wide_boundary = True
            
        elif self.transform_index == 3: # Fourier-Bessel
            Transform = global_basis_transform('bessel',*params)
            self.wide_boundary = True
        else:
            raise Exception("missing transform")
            return        

        self.npix = size(Transform,1)

        debug('Transform done')
        
        return  Transform


    def get_boundary(self,N=100, coords = None,time=None,k=1):

        if  coords is None:
            if self.boundary_lcfs and not time is None and self.boundary > 0:
                #get the separatrix
                rho,magr,magz=self.mag_equilibrium(time,surf_slice=[-1,],rho=1)
                coords = hstack((magr,magz))[:,:,0]
            else:
                coords = self.vessel_boundary
      
        if any(coords[-1,:] != coords[0,:]) and linalg.norm(coords[-1,:]-coords[0,:])<1e-10:
            #rounding error
            coords[0,:] = coords[-1,:]
        elif any(coords[-1,:] != coords[0,:]):
            coords = r_[coords, coords[(0,),:]]
        
        try:
            tckp,u = splprep(coords.T,s=1e-5,k=k, per=1)
            boundary_p = array(splev(linspace(0,1,N),tckp)).T
        except:
            return coords

        return boundary_p



        
    def prepare_tokamak(self):
        debug('prepare_tokamak')

        if True or not 'Tmat' in vars(self) or size(self.Tmat,1) != self.nx*self.ny :
            from geom_mat_setting import geom_mat_setting
            self.Tmat , self.Xchords, self.Ychords = geom_mat_setting(self,self.nx,  self.ny, self.virt_chord)
        
        self.correction_matrix = ones(self.npix)  # correction of artifacts in creation of the weigthing matrix 

        # First guess of emissivity, can improve first step and convergence ??
        print("===========create tokamak===========")
            
        # used for normalization in MFI
        self.normTmat = self.Tmat.sum()/self.Tmat.nnz
 
        debug( "object tokamak done")


    def get_data(self, failsave = True,tmin = -infty,tmax = infty,return_err=False,tind=None):
        """
        Returns loaded data from cache -> prevents memory issues
        """        
        #BUG very ugly solution!!!
        
  
        try:
            
            
            debug( 'load from cache')
            d = load(self.geometry_path+'data_cache.npz')

            if d['shot'] != self.shot:
                raise  Exception('wrong shot number')
           
            dets_pruned = d['dets_pruned']

            if return_err:
           
                saturated = d['saturated'] if 'saturated' in d else None 
                return d['err'], saturated
            
            
            data = d['data']
            tvec = d['tvec']

            assert shape(data)[0] == len(tvec), "Wrong size of saved data,\
                regenerate data !!" + str(shape(data)) + " "+str(shape(self.dets))
                
            t_ind = slice(*tvec.searchsorted([tmin,tmax]))
   

            
            errors = self.generate_errors(data[t_ind],t_ind,all_types=True)
            return tvec[t_ind],data[t_ind], errors
    
            print('data loaded from cache')
        except Exception as e:
            print('get_data '+str(e))

            if failsave:
                time.sleep(1)
                print('get_data self.prepare_tokamak()')
                return self.get_data( False, tmin = tmin,tmax = tmax,return_err=return_err)
            else:
                print('failsave '+str(failsave))
                raise
 
    def generate_errors(self, data,t_ind,all_types=True):
        """
        Try to generate expected errors according to values of min_error and sigma in
        class tokamak. The expected errors are partially constant and partialy depends on emissivity
        """
        N,ndet = data.shape
        error  = fabs(data)+1e-10
        std_offset,saturated = self.get_data(return_err=True)
        std_offset = std_offset[t_ind]
        
        if std_offset is None or size(std_offset) == 1: #not defined
            std_offset = zeros_like(data[0])

        if saturated is None or size(saturated) == 1: #not defined
            saturated = zeros_like(data,dtype='bool')
        else:
            saturated = saturated[t_ind]
  
        saturated |= data < - abs(mean(data))*0.1  #sometimes are saturated channels negative
        std_offset[isnan(std_offset)] = 0   #error from the remooving of the offset
        
   

        ind = slice(None,None)   if N < 1000 else random.randint(N,size=1000)
        tmp_error = error[ind]
        if N>10:  m = mquantiles(tmp_error, 0.9,axis=0)
        else:     m = nanmean(tmp_error,axis=0)

        
        Nsmooth = max(1,N//20)


        error = (1/error[:(N/Nsmooth)*Nsmooth].reshape(N//Nsmooth, Nsmooth,ndet)).mean(1)
        error[error!= 0] = error[error!= 0]
        error[error== 0] = infty
        error = 1/repeat(c_[error.T, error.T[:,-1]], Nsmooth, 1)[:,:N].T        

        error[(data.mean(0) <= 1e-6)|isnan(error)] = infty
        error[saturated] = infty

        error = sqrt(error, out = error)   # ugly fix for old  numpy
 
        error *= self.sigma*sqrt(abs(m))

        error += self.min_error * nanmean(data.ravel())
        error += std_offset

        if any(isnan(error)):
            print('nans in errors')

        return error
    


    def mag_equilibrium(self, tvec, preferCache = True, dryRun = False, return_mean = False,
                      surf_slice=slice(None,None),n_rho=None,n_theta=None,rho=None):
    
        """ Method to load magnetic field (contours). Try to load from local cache,
        if it is not possible try to load from JET network and store it localy. 
        Loaded field is only with limited precision (max 100 timeslices and it is interpolated)

        :param array tvec: Time vector  of requested  magnetic field
        :param bool dryRun: Do not return interpolated field
        :param bool preferCache: Use cached data if possible, it is faster

        :var array magx, magy: Arrays of magnetic lines possitions
        :var array tsurf: Time vector that belongs to ``magx``, ``magy``
        :raises Asseration: No connection to JET network"
        :rtype:  Coordinates of magnetic field interpolated to the requested time vector ``tvec``
        :int n_rho   number of returned flux surfaces (not avaible for all tokakamks!!) 
        :int n_theta   number of returned flux surfaces angles (not avaible for all tokakamks!!) 

        """
        
        mag_path = os.path.join(self.geometry_path,'equilibrium','MagField_'+str(self.shot)+'.npz')

        
        if not hasattr(self, 'magx'):
            try:
                assert  preferCache, 'Load regulary data'   #  (self.name == "JET" and self.local_network)
                
                d = load( mag_path )
                

                self.tsurf = d['tsurf']
                magx = d['magx']
                magy = d['magy']
                
    

                ntheta,nrho,nt = magx.shape
                

                
                if self.name == 'JET':   
                    #Correction of wrong equilibria from JET!!!

                    R = mean(std(magx, axis = 0), axis = 1)
                    ind  = argsort(R)
                    self.magx = zeros((ntheta,nrho+1,nt))
                    self.magy = zeros((ntheta,nrho+1,nt))

                    self.magx[:,1:,:] = magx[:,ind,:]
                    self.magy[:,1:,:] = magy[:,ind,:]

                    xc = magx[:,0,:].mean(0)[None,:]
                    yc = magy[:,0,:].mean(0)[None,:]
                    
                    self.magx[:,0,:] = (magx[:,ind[0],:]-xc)/10+xc
                    self.magy[:,0,:] = (magy[:,ind[0],:]-yc)/10+yc
                
                else:
                    self.magx = magx
                    self.magy = magy


            except Exception as e:
                print("preferCache", preferCache, str(e))
                #download the equilbrium
                self.load_mag_equilibrium()
                
        
        rho = linspace(0,1,self.magx.shape[1])[surf_slice]
        magx = self.magx[:,surf_slice,:]
        magy = self.magy[:,surf_slice,:]


        if not dryRun:
            # FIXME in case when the data are out of intepolation range

            tvec = copy(tvec)
            tvec = reshape(tvec, (-1))

            tvec[tvec < self.tsurf[0]] = self.tsurf[0]
            tvec[tvec > self.tsurf[-1]] = self.tsurf[-1]
            
            
            from scipy.interpolate import  interp1d
            if not return_mean or len(tvec) == 1:
                magx=interp1d(self.tsurf,magx,assume_sorted=True,copy=False)(tvec)
                magx = magx.astype('float32')
                magy=interp1d(self.tsurf,magy,assume_sorted=True,copy=False)(tvec)
                magy = magy.astype('float32')

            else:
                # save memory usage
                ind = slice(max(0,self.tsurf.searchsorted(amin(tvec))-1), 
                            self.tsurf.searchsorted(amax(tvec)))

                magx = mean(magx[:,:,ind],2)   #take mean field during tvec and omitted the center of field [:,0]
                magy = mean(magy[:,:,ind],2)
        
            from shared_modules import read_config

            input_parameters = read_config('tomography'+".cfg")
        
            if config.magfield_shift != (0,0):
                self.magfield_shift = config.magfield_shift
            else:
                self.magfield_shift = input_parameters['magfield_shift_core']
           
            self.magfield_shift_lcfs = input_parameters['magfield_shift_lcfs']

            if self.magfield_shift != [0,0]:
                x_shift = self.magfield_shift[0]+rho*(self.magfield_shift_lcfs[0]-self.magfield_shift[0])
                y_shift = self.magfield_shift[1]+rho*(self.magfield_shift_lcfs[1]-self.magfield_shift[1])

                if ndim(magx) == 3:
                    magx += x_shift[:,None]*self.norm
                    magy += y_shift[:,None]*self.norm
                else:
                    magx += x_shift*self.norm
                    magy += y_shift*self.norm
                    
        if self.radial_coordinate == 'r_a'  and size(rho)>1:
            rho,magx, magy = self.convert_rho_2_r_a(rho, magx,magy)
        
        if self.radial_coordinate == 'r_V'  and size(rho)>1:
            rho,magx, magy = self.convert_rho_2_r_V(rho, magx,magy)
 
                    
        return rho,magx, magy
   
                 
    def convert_rho_2_r_a(self, rho, magx,magy ):
        #convert from rho_pol to r/a 
        r_a  = hypot(magx-magx[:,(0,)], magy-magy[:,(0,)]).mean(0)
        lfs = amax(magx,axis=0)
        hfs = amin(magx,axis=0)
        r_a = (lfs-hfs)/2

        i_rho_sep = rho.searchsorted(1)
        r_a /= r_a[(i_rho_sep,)]
        if magx.ndim == 2:
            r_a = squeeze(r_a)
            magx[:] = interp1d(r_a, magx,axis=1,assume_sorted=True)(minimum(rho, r_a[-1]))
            magy[:] = interp1d(r_a, magy,axis=1,assume_sorted=True)(minimum(rho, r_a[-1]))
        else:
            for it in range(magx.shape[2]): #slowest
                magx[:,:,it] = interp1d(r_a[:,it], magx[:,:,it],axis=1,assume_sorted=True)(minimum(rho, r_a[-1,it]))
                magy[:,:,it] = interp1d(r_a[:,it], magy[:,:,it],axis=1,assume_sorted=True)(minimum(rho, r_a[-1,it]))
  
    def convert_rho_2_r_V(self,rho, magx,magy, ):
        #convert from rho_pol to r_V
        i_rho_sep = rho.searchsorted(1)

        V = sum(2*pi*((magx[1:]+roll(magx[1:],1,0))/2)**2*(magy[1:]-roll(magy[1:],1,0)),0)/2
        r_V = sqrt(maximum(0,V/V[i_rho_sep]))
       
        if magx.ndim == 2:
            magx[:] = interp1d(r_V, magx,axis=1,assume_sorted=True)(minimum(rho, r_V[-1]))
            magy[:] = interp1d(r_V, magy,axis=1,assume_sorted=True)(minimum(rho, r_V[-1]))
        else:
            for it in range(magx.shape[2]): #slowest
                magx[:,:,it] = interp1d(r_V[:,it], magx[:,:,it],axis=1,assume_sorted=True)(minimum(rho, r_V[-1,it]))
                magy[:,:,it] = interp1d(r_V[:,it], magy[:,:,it],axis=1,assume_sorted=True)(minimum(rho, r_V[-1,it]))



    def mag_theta_star(self,time,rho,magr,magz,rz_grid=True, extrapolate = 1.1, nx=None,ny=None):

        if self.radial_coordinate  == 'r_a':
            raise Exception('for r_a coordinate is not posssible to calculate mag_theta_star!!')
            
        n_theta, n_rho = magr.shape
        theta = linspace(0,2*pi,n_theta,endpoint=False)

        #is it OK?
        magr, magz = copy(magr.T), copy(magz.T)
        
        r0,z0 = magr[0].mean(), magz[0].mean()
        #calculate gradient of Phi with resoect to R and z
        magr[0] +=  (magr[1]-magr[0])/100  #coorrection of the singularity in the  center
        magz[0] +=  (magz[1]-magz[0])/100


        drdrho,drtheta = gradient(magr)
        dzdrho,dztheta = gradient(magz)
        dpsidrho,dpsitheta = gradient(tile(rho**2, (n_theta,1)).T )

        grad_rho = dstack((drdrho,dzdrho,dpsidrho ))
        grad_theta = dstack((drtheta,dztheta,dpsitheta))
        normal = cross(grad_rho,grad_theta,axis=-1)

                

        dpsi_dr = -normal[:,:,0]/(normal[:,:,2]+1e-8) #Bz
        dpsi_dz = -normal[:,:,1]/(normal[:,:,2]+1e-8) #Br

    #WARNING not defined on the magnetics axis

        dtheta_star = ((magr-r0)**2+(magz-z0)**2)/(dpsi_dz*(magz-z0)+dpsi_dr*(magr-r0))/magr
        theta = arctan2(magz-z0,-magr+r0)

        
        theta = unwrap(theta-theta[:,(0,)],axis=1)


        
        from scipy.integrate import cumtrapz

        #definition of the thetat star by integral
        theta_star = cumtrapz(dtheta_star,theta,axis=1,initial=0)
        correction = (n_theta-1.)/n_theta
        if all(magr[:,0]==magr[:,-1]) and all(magz[:,0]==magz[:,-1]):
            correction = 1
        theta_star/= theta_star[:,(-1,)]/(2*pi)/correction     #normalize to 2pi
            
            
        if not rz_grid:
            return theta_star

        
        from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
        if not (self.boundary_lcfs and self.boundary) or rho.max() < 1:
            bdx,bdy = self.get_boundary(size(magr,1),time=time).T
        else:
            bdx,bdy = magr[-1], magz[-1]
        
        magx_out = (bdx-mean(bdx))*extrapolate/rho[-1]+mean(bdx)
        magy_out = (bdy-mean(bdy))*extrapolate/rho[-1]+mean(bdy)
        
        Ninterp = NearestNDInterpolator(c_[magr[-1],magz[-1]],theta_star[-1])
        magz = c_[magz.T, magy_out].T
        magr = c_[magr.T, magx_out].T
        
        theta_out = Ninterp(c_[ magx_out,magy_out ])
        theta_star_ = c_[theta_star.T, theta_out].T


        cos_th, sin_th = cos(theta_star_), sin(theta_star_)
        

        Linterp = LinearNDInterpolator(c_[magr.ravel(),magz.ravel()], cos_th.ravel())
        
        
        if nx is None: nx = self.nx
        if ny is None: ny = self.ny

        #centers of the pixels
        xgridc = linspace(self.xmin,self.xmax,nx+1)
        xgridc = (xgridc[1:]+xgridc[:-1])/2
        ygridc = linspace(self.ymin,self.ymax,ny+1)
        ygridc = (ygridc[1:]+ygridc[:-1])/2
        
        BdMat = get_bd_mat(self, nx, ny ,boundary=c_[bdx,bdy])

        X,Y = meshgrid(xgridc,ygridc)  #centers of pixels
        X,Y = X.flatten('F')[~BdMat],Y.flatten('F')[~BdMat]
        cos_grid = Linterp(c_[X.ravel(),Y.ravel()]).reshape(X.shape)
        Linterp.values[:,0] = sin_th.ravel()#trick save a some  computing time
        sin_grid = Linterp(c_[X.ravel(),Y.ravel()]).reshape(X.shape)  
       
        theta_star_ = arctan2(sin_grid, cos_grid)
        
    
        theta_star_rz = zeros((nx,ny))
        theta_star_rz.flat[~BdMat] =  theta_star_
        
       
        return theta_star_rz, theta_star
        
        
    def prepare_emiss0(self, time):
        
        from phantom_generator import phantom_generator
        self.emiss_0 = phantom_generator(self,asarray(time), self.nx, self.ny, 
                                        profile = self.rad_profile )[3]

        self.emiss_0  = reshape(self.emiss_0, (-1,1))
        
        
        

    def prepare_data(self, tmin, tmax, data_smooth, data_undersampling, error_scale, 
                        prune_dets = 1,  detsCutOff = True, allowed_dets = None):
        """
        Main function that prepare data.

        *   Crops data to slected range ``tmin``, ``tmax``
        *   Use under sampling to decrease number of reconstructed/plotted timeslices.
        *   Use moving average to smooth the data.
        *   Remove data that belongs to ``wrong_dets``
        *   Remove NaN detectors
        *   Renormalize data form sensors according to tokamak.calb


        :param class self: ``Class`` of  tokamak with all important data.
        :param double tmin, tmax: Time range reconstruction
        :param int data_smooth: Number bins used for regularization (moving average)
        :param int data_undersampling:    Use every `n`-th snapshot
        :param double error_scale:  Increase size of expected errobars
        :param bool detsCutOff:  Remove wrong detectors from data
        :var array data:  Final array of data
        :var array tvec:  Final array of time vecor
        :var array dets:  Final array of dets
        :var int normTmat:    Norm of geometry matrix
        :var array normData:  Maximum of data in each timeslice
        """

   
        errorSmooth = max(data_smooth * data_undersampling, 10)  # errobars should be always smoothed
        
        dt = data_smooth/self.sample_freq/2*1.1

        tvec,data,error = self.get_data(tmin=tmin-dt,tmax=tmax+dt)

        data_smooth = min(max(int(data_smooth),1), len(tvec))
        
        data,error = copy(data.T), copy(error.T) 
     

        Tmat = self.Tmat if hasattr(self,'Tmat') else None

        debug('\nGenerate errors')
   

        dets = self.dets if allowed_dets is None else allowed_dets
            
    
        self.nl = len(dets)  #number of LOS
        mean_data = mean(data)
        
     
        if self.use_median_filter and data_smooth < 1000:
            debug('Median filter of data ')
            for j in dets:   
                data[j]    = medfilt(data[j], (data_smooth//2)*2+1)
                if not all(isfinite(error[j])):   continue
                error[j] = 1/medfilt(1/error[j], (data_smooth//2)*2+1)

                
        elif size(data,1) > 1 and data_smooth > 1:
            for i in dets: 
                ind = (isfinite(error[i])) & (error[i] > 0)
                if sum(ind) > data_smooth:
                    error[i,ind] =  1./MovingAveradge(1./error[i,ind],data_smooth,-1)
                error[i,~ind] = infty
                data[i] = MovingAveradge(copy(data[i]),data_smooth)
        
        imin,imax = tvec.searchsorted((tmin,tmax))
        t_ind = slice(imin,imax+1, data_undersampling)

        nerr =  max(min(5,data_smooth ),2)

        du = max(1,int(round(min(data_smooth, len(tvec))/float(nerr))))

        data_reshaped = data[:,:len(tvec)//du*du].reshape(-1,len(tvec)//du,du)
        data_reshaped = data_reshaped.mean(-1)
        nt = data_reshaped.shape[-1]
        du2 = min(nerr, nt)

        data_reshaped = data_reshaped[:,:nt//du2*du2].reshape(-1,nt//du2,du2)
        data_reshaped = diff(data_reshaped, axis=2)
        std_data = einsum('ijk,ijk->ij', data_reshaped,data_reshaped)/(du2-1)-data_reshaped.mean(2)**2

        std_tvec = tvec[:len(tvec)//du*du].reshape(len(tvec)//du,du)
        std_tvec = std_tvec.mean(-1)
        std_tvec = std_tvec[:nt//du2*du2].reshape(nt//du2,du2).mean(-1)



        
        std_data[std_data<0] = 0  #numerical errors
        std_data[~isfinite(std_data)] = 0  #overflow
        
        std_data = sqrt(std_data)
        std_data/= sqrt(du2)


        data  = data[:,t_ind]   #use every n-th snapshot
        error = error[:,t_ind]
        tvec  =  tvec[t_ind]

        error /= sqrt(sqrt(data_smooth)) #just guess!!
        std_tvec[[0,-1]] =  tvec[0],  tvec[-1]
        if len(std_tvec) > 1 or len(std_tvec) != len(tvec):
            error+= interp1d( std_tvec, std_data,axis=1,assume_sorted=True)(tvec)
        else:
            error+= std_data

        if any(isnan(data)): debug( 'data prepare nan!!')
        ind_correct = self.get_correct_dets(data)

        dets_index = self.dets_index
        self.calb = array(self.get_calb())
  
        
        if not all(self.calb==1):
            if len(self.calb) == self.Ndets:
                for i in range(self.Ndets):
                    data[dets_index[i], :] *= mean(self.calb[i])   # if case of different time vector
                    error[dets_index[i], :] *= mean(self.calb[i]) 
            elif len(self.calb) == self.nl:
                data*= self.calb[:,None]
                error*= self.calb[:,None]
                
     

        dets = dets[ind_correct]

        if detsCutOff:
            if not any(ind_correct):
                raise Exception('All data channels was removed!! ')

            if Tmat!= None: Tmat= Tmat[dets,:]
            error=error[ind_correct,:]
            data= data[ind_correct,:]
        
    
            
        if len(tvec) > 2: #remove "blocked detectors"
            error[:,all(abs(data-nanmean(data))<1e-3,0)]  = infty
        
        normData = nanmax(data, axis=0) + nanmax(data)*1e-6

        gc.collect()

        error *= error_scale
        
        
        if data.dtype != float32 or error.dtype != float32:
            print('Warning Data or error are not in a single precision', data.dtype, error.dtype)
            data = data.astype(float32, copy=False)
            error = error.astype(float32, copy=False)

        debug( "prepare data done")
       
        return data, error, tvec, dets, Tmat, normData


    def get_calb(self):
        
        calib_file = self.geometry_path+'calibration'+os.sep+'%d.txt'%(self.shot)
        calb = self.calb_0

            
        if self.default_calb != 'none' and self.Ndets > 1:# and  tokamak.allow_self_calibration:
            if not os.path.exists(self.geometry_path+'calibration'):
                os.makedirs(self.geometry_path+'calibration')
                
            if not os.path.isfile(calib_file): 
                from os import listdir
    
                files = listdir(self.geometry_path+'calibration')
                shots = []
                for f in files:
                    if not f.endswith('.txt'): continue
                    num =  f.split('.')[0]
                    if num.isdigit(): shots.append(int(num))
                    
                if len(shots):
                    shots.sort()
                    ifile = argmin(abs(array(shots)-self.shot))
                    calib_file = self.geometry_path+'calibration'+os.sep+'%d.txt'%(shots[ifile])
 
            #BUG OLD AND NEW FORMAT OF THE CALIBRATION FILES
            try:
                cams,calb_ = loadtxt(calib_file, dtype={'names': ('cam', 'calib'),'formats': ('S4', 'd')}, unpack=True)
            except:
                #old format
                try:
                    calb_ = asarray(loadtxt(calib_file))
                    cams = asarray(list(self.detectors_dict.keys()))
                    f = open(calib_file, 'w')
                    for k,i in zip(cams, calb_):
                        f.write('%s  %.3f\n'%( k,i ))
                    f.close()
                except:
                    calb_ = None
      
                if not calb_ is None:
                    calb_ = {c:c_ for c,c_ in zip(cams,calb_)}
                    try:
                        calb_ = [calb_[c] for c in list(self.detectors_dict.keys())]
                        if size(calb_) in (self.Ndets,self.nl):    calb = calb_ 
                        else:
                            print('wrong length of the calibration file', size(calb_), size(self.dets_index), self.Ndets)
                    except Exception as e:
                        print(e)
                                    
            if not calb_ is None and len(calb_) == len(self.dets_index):
                calb = calb_ 

        #save calibration
        if self.default_calb != 'none':
            if not os.path.exists(self.geometry_path+'calibration'): 
                os.makedirs(self.geometry_path+'calibration')
            cams = list(self.detectors_dict.keys())
            f = open(calib_file, 'w')
            for k,i in zip(cams, calb):  f.write('%s  %.3f\n'%( k,i ))
            f.close()    
            

        calb = asarray(calb)
        calb[(calb == 0)|~isfinite(calb)] = 1  
   
        return calb

    def get_correct_dets(self, data = None, wrong_dets = [] , include_pref  = True):

        ind_correct = ones_like(self.dets, dtype=bool)
        if include_pref:
            ind_correct &= ~in1d(self.dets, config.wrong_dets_pref)         # is in "wrong_dets_pref"       
   
        ind_correct &= ~in1d(self.dets, self.wrong_dets_damaged)         # is in damaged detectors,
        ind_correct &= ~in1d(self.dets, wrong_dets)
        
  

        if not data is None:
            if any(isnan(data)):
                debug( 'NaNs in the data detected!!')
        
                ind_correct &= ~any(isnan(data),1)             # have NaN values

            if not self.camera:
                ind_correct[ind_correct]= ~all(data[ind_correct]<=1e-5,1)                 # selected range values are zeros

        return ind_correct


