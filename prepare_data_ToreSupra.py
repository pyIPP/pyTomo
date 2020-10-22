#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from numpy import *
import matplotlib
from scipy import sparse
import sys
from scipy.io import loadmat, savemat
#from matplotlib.pyplot import *
from numpy.matlib import repmat
from tokamak import Tokamak
import config
from shared_modules import blur_image, nanmedian

class ToreSupra(Tokamak):
    def __init__(self, input_diagn, input_parameters,  only_prepare = False  ):
	#(self, input_path, input_diagn,  shot, nx, ny, virt_chord, MagField_shift, only_prepare= False):
        # BASIC PARAMETERS
        #print "inicializace"
        name = 'ToreSupra'

        #coord = [1700, 3100, -700, 700]
	#coord = [1.600, 3.300, -.800, .800]
	coord = [2.2, 2.9, -.4, .65] # camera 1  # 12.5
	#coord = [2.3, 3., -.5, .6] # camera 1  # 12.5

	#coord = [2.75, 3.1, -.05, .4] # camera 2
	#coord = [2.65, 3.1, -.05, .4] # camera 2

    
        sigma = 0.00
        min_error = 0.1
        #self.shot = shot
        #norm = 1000.0
	norm = 1.0

        t_name = 'ms'   # time in seconds
        virt_chord = 1
        
	path = 'geometry_toresupra/'+input_diagn
	Tokamak.__init__(self, input_parameters, input_diagn, coord,  name, norm, sigma, min_error, path,t_name)

	#self.vessel_coord = [1.600, 3.300, -.800, .800]
	#self.plot_coord = coord
	self.plot_coord =  [2.3, 2.8, -.35, .55]  # camera 1
	
	
	
        #self.nx = nx
        #self.ny = ny
        #self.npix = nx*ny
        #self.dx =(self.xmax-self.xmin)/double(self.nx)
        #self.dy =(self.ymax-self.ymin)/double(self.ny)
        #self.xgrid=arange(self.xmin,self.xmax, self.dx) # !! linspace is less precise !!!
        #self.ygrid=arange(self.ymin,self.ymax, self.dy)
        #self.boundary = loadtxt(self.geometry_path+'/border.txt')

        #self.MagField_shift = MagField_shift

	if input_diagn == "SXR":
	    self.SXR()
	    self.index = 5
	if input_diagn == "Camera":
	    self.Camera(input_parameters)
	    self.index = 18

        try:
            self.mag_equilibrium(self.tvec,True, True)
        except:
            raise

        #CREATE GEOMETRIC MATRIX
        from annulus import get_bd_mat
        self.BdMat = get_bd_mat(self)     #detection of boundary, g_model(sum(T)==0) = 0 is not correct

        from geom_mat_setting import geom_mat_setting
        self.Tmat , self.Xchords, self.Ychords = geom_mat_setting(self,self.nx,  self.ny, self.virt_chord)
        self.tsteps = len(self.tvec)

        #self.dets = arange(size(self.Tmat, 0))

        self.normTmat = self.Tmat.sum()/self.Tmat.nnz

        from emissivity_generator import emissivity_generator


        # First guess of emissivity, can improve first step and convergence ??

        try:
            tmp, tmp, self.emiss_0 = emissivity_generator(self,array([mean(self.tvec)]), self.nx, self.ny)
        except:
            print("emissivity failed")
            self.emiss_0 = ones((self.nx*self.ny))
        self.emiss_0  = reshape(self.emiss_0, (-1,1))

        print("dokoncena inicializace")

    #def generate_errors(self, data):
        #"""
        #Try to generate expected errors according to values of min_error and sigma in class tokamak. The expected errors are partially constant and partialy depends on emissivity
        #"""
        #N = size(data,1)
        #error  = data.copy()
        #correct = ((error > 0) & ~isnan(error))
        #error[~correct] = amin(error[correct])
        #m = amax(error, 0)
        #error = sqrt(error)
        #error =  error* self.sigma * sparse.spdiags(amax(error,0), 0, N,N)
        #error += self.min_error * m

        #return error

    def SXR(self):


        try:
            Input = loadmat(os.path.normpath(self.geometry_path+'/TSdata_example.mat'), struct_as_record=True)
        except:
            raise IOError("Can't load", None, os.path.normpath(self.geometry_path+'/TSdata_example.mat'))
        self.data = Input['Y'].T

        #self.wrong_dets = []

        #print "uuuuuuuuuuuuuuuuuuuuu"

        #print shape(Input['dets']), shape(self.data)
        self.dets = arange(amax(squeeze(Input['dets']-1)))
        self.wrong_dets_damaged = delete(self.dets,squeeze( Input['dets'])-1)

        self.tvec = Input['tsel']

	self.save_data( self.data)


    #else:
        ##generate artificial data
        #tvec = arange(50,52.5, 0.01)
        #[data, tvec, dets, wrong_dets] = createdata(sigma,min_error, xmin, xmax, ymin, ymax, shot, tvec, tokamak,shifts,alfa0 )

    #jet = loadtxt(geometry_path+'/dataJET.txt').T
    ##jet[isnan(jet)] = 0
    ##data[isnan(data)] = 0
    ##clf()
    ##plot(jet)
    ##plot(data.T)
    ##show()
    ##print shape(jet),  shape(data)
    ##exit()
    #data = jet

    #RIGID_BODY = shifts+1
    #tsteps = len(tvec)
    #ndets = len(dets)
    #print  'pred rigid body', RIGID_BODY
    #if  RIGID_BODY>1:
        #print '========rigid body data TORE============= ', shape(data),  RIGID_BODY, tsteps, ndets, RIGID_BODY
        #print tvec, shape(tvec)
        #data_tmp = zeros((RIGID_BODY*ndets , tsteps+RIGID_BODY))
        #if len(tvec) > RIGID_BODY:
            #for i in range(RIGID_BODY):
                #print i
                #print shape(data)
                #print shape( data[:,i:(i+RIGID_BODY)])
                #print 'tmp', shape(data_tmp)

                #data_tmp[i*ndets:(i+1)*ndets, i:(tsteps+i)]  = data

        #data = data_tmp[:,RIGID_BODY:-RIGID_BODY]
        #tvec = tvec[:-RIGID_BODY]
        #dets = arange(size(data,0))
        ##data[isnan(data)] = 0
        ##pcolor(data)
        ##show()
        ##print where(isnan(data))
        #print 'OUTPUT !!!!', shape(data), shape(dets), shape(tvec)



        ##plot(data[~isnan(data)])
        ##show()
        ##exit()
    #return data, tvec, dets, wrong_dets

            ##del dY

#def createdata(sigma,min_error,xmin, xmax, ymin, ymax, shot, tvec, tokamak, shifts,alfa0):
    ##sigma = 0.0001
    #wrong_dets = []
    #tsteps = len(tvec)
    #alpha = 2*pi/(shifts+1) *arange(tsteps)+pi/4
    ##alpha = 2*pi/(tsteps+1)*arange(tsteps)+alfa0


    #from geom_mat_setting import geom_mat_setting
    #To,X_tmp, Y_tmp =  geom_mat_setting(50,50,xmin, xmax, ymin, ymax,'geometry_toresupra', 1,  shot, tvec, tokamak,0,alfa0)
    #print "alfa0===========", alfa0, size(To,0)
    #dets = arange(size(To,0))
    #from generator_dat import generator_dat
    #[data, dY, G_orig]=generator_dat(To, 50, 50, tvec,sigma,min_error, shifts, alpha)
    #return data, tvec, dets, wrong_dets

    #47760_1
    
    
    def Camera2(self, input_parameters):
        print("============= camera =============")
        input_parameters['boundary'] = False
	if not self.artif_data:
	    input_parameters['error_scale'] = 10
	    
        self.camera = True
        self.allow_self_calibration = False
        self.wrong_dets_pref = []
        self.min_error = 0.05
        self.sigma = 0.00
        #self.allow_divertor = False
        #self.smooth_boundary = False
        self.virt_chord = 1
        self.SOL_pos = 0.9
        self.allow_negative = True
        self.weight_power = 1
        self.expect_zero = True
	self.geometry_axis = array([self.xmin,self.xmax, self.ymin, self.ymax, 0.1*self.xmax , self.xmax])

	#N_steps = 200

	if True:  # not self.artif_data:
	    print("loading data")
	    data = list()
	    #from scipy.signal import medfilt2d
	    data_0 =loadmat(self.geometry_path + "/47760_1_sub.mat")['img_sub']
	    print("shape", shape(data_0))
	    
	    for i in arange(size(data_0,2)):
		d = double(data_0[:,:,i])
		if nanmedian(abs(d)) < 1:
		    continue
		#d =  blur_image(d, 1)
		d = d[::2, ::2]
		d =  rot90(rot90(rot90(d)))
		#imshow(d)
		#show()
	
		#d = medfilt2d(d, (3,3))
		#d = d[::3,::3]
		
		res = shape(d)
		d = ravel(d)[:,None]
		data.append( d )
		
		
	    data = hstack(data)
	    #data = data[:,:,220:]
	    #for i in range(size(data,1)):
		#imshow(reshape(data[:,i], res), interpolation="nearest")
		#show()
		
	    #data = data[::2, ::2, :]

	    #res = shape(data)[:2]
	    N_steps = shape(data)[1]
	    #print N_steps
	    #exit()
	    #data = reshape(data, (-1, N_steps) )
	    #imshow(isnan(data))
	    #show()
	    #exit()
	    
	    self.tore_data = reshape(data, (res[0], res[1], -1))

	    self.save_data( data )
	    #del data 
	    
	if  self.artif_data:
	    res = (30, 50)
	    N_steps = 200
        # artificial !!!!
        tvec = linspace(0, 10, N_steps)

	self.camera_res = res


        nl = prod(res)
        dets = arange(nl)

	# camera 2
	coord =   [2.7, 3.1, -0.05, .45] #
	cam_coord = [3.00, 0.0, 3.0]
	
	
        distorsion = 0
	alpha_rot =  0
	center = array([0,0]) # distorsion center
	self.cam_zoom = 1 # 4
	self.pixel_zoom = 1 # 4

	a1 = arctan2(coord[1] - cam_coord[0], cam_coord[2])
	a2 = arctan2(coord[0] - cam_coord[0], cam_coord[2])
	xgrid = cam_coord[2]*tan(linspace(a1, a2, res[0] * self.cam_zoom)) + cam_coord[0]

	a1 = arctan2(coord[3] - cam_coord[1], cam_coord[2])
	a2 = arctan2(coord[2] - cam_coord[1], cam_coord[2])
	ygrid = cam_coord[2]*tan(linspace(a1, a2, res[1] * self.cam_zoom)) + cam_coord[1]

	from distorsion import  coord_distort
	[xend, yend]  = coord_distort(xgrid, ygrid, center, distorsion)
        xmean = mean(xend); ymean = mean(yend)
	xend, yend = [ (xend-xmean)*cos(alpha_rot) + (yend-ymean)*sin(alpha_rot)+ xmean, (xmean- xend)*sin(alpha_rot) + (yend-ymean)*cos(alpha_rot)+ymean ]

	self.cam_coord = cam_coord;

        zend = zeros(nl * self.cam_zoom**2)  # mriz je v ose

        xstart = ones(nl * self.cam_zoom**2) * cam_coord[0] #+ cam_x
        ystart = ones(nl * self.cam_zoom**2) * cam_coord[1] #+ cam_y
        zstart = ones(nl * self.cam_zoom**2) * cam_coord[2] #+ cam_z

	##################################################
	xend  = xend[::-1, ::-1]
	yend  = yend[::-1, ::-1]
	##################################################
	
	xend = reshape(xend, (-1), order="F")
	yend = reshape(yend, (-1), order="F")
	zend = reshape(zend, (-1), order="F")

	print("saving geometry")
	i = 0
	save(self.geometry_path+'/detector_'+str(i)+'_x.npy', array([xstart, xend ]).T)
	save(self.geometry_path+'/detector_'+str(i)+'_y.npy', array([ystart, yend ]).T)
	save(self.geometry_path+'/detector_'+str(i)+'_z.npy', array([zstart, zend ]).T)
	print("done")
	
	
        self.calb = ones(len(dets))
        self.dets_index = [arange(i*res[0],(i+1)*res[0] ) for i in range(res[1]) ]

        self.tvec = tvec
        self.dets = dets
	self.nl = len(dets)




    def Camera(self, input_parameters):
        print("============= camera =============")
        input_parameters['boundary'] = False
        if not self.artif_data:
	    input_parameters['error_scale'] = 10
	else:
	    input_parameters['error_scale'] = 5

        self.camera = True
        self.allow_self_calibration = False
        self.wrong_dets_pref = []
        self.min_error = 0.05
        self.sigma = 0.00
        #self.allow_divertor = False
        #self.smooth_boundary = False
        self.virt_chord = 1
        self.SOL_pos =  0.3
        self.allow_negative = True
        self.weight_power = 1
        self.expect_zero = True
	self.MagField_shift = [-0.25,0]
	print("posunuji magneticke pole !!!!!!!!!!!!!!!!!")

	N_steps = 200
        tvec = linspace(50, 52, N_steps)

	#if not self.artif_data:
	    #xxxx
	print("loading data")
	data = list()
	from scipy.signal import medfilt2d

	for i in 1+arange(N_steps):
	    d = double(loadmat(self.geometry_path + '/data_new/data_sub_'+str(i)+'.mat')['data'])
	    d = d[100:, :150]

	    d =  rot90(rot90(rot90(d)))
    
	    d = medfilt2d(d, (3,3))
	    d = d[::3,::3]
	    
	    res = shape(d)
	    d = ravel(d)[:,None]
	    data.append( d )
	data = hstack(data)
	
	
	self.save_data( data )
	    #del data 
	self.tore_data = reshape(data, (res[0], res[1], -1))

	    
	if self.artif_data:
	    res = (32, 59)
        # artificial !!!!

	
	

	self.camera_res = res


        nl = prod(res)
        dets = arange(nl)

	# camera 1
	coord =   [2.2, 3, -.35, .7] #
	cam_coord = [3./sqrt(2), 0.1, 4]
	

	
        distorsion = 0
	alpha_rot =  0
	center = array([0,0]) # distorsion center
	self.cam_zoom = 1 # 4
	self.pixel_zoom = 1 # 4

	a1 = arctan2(coord[1] - cam_coord[0], cam_coord[2])
	a2 = arctan2(coord[0] - cam_coord[0], cam_coord[2])
	xgrid = cam_coord[2]*tan(linspace(a1, a2, res[0] * self.cam_zoom)) + cam_coord[0]

	a1 = arctan2(coord[3] - cam_coord[1], cam_coord[2])
	a2 = arctan2(coord[2] - cam_coord[1], cam_coord[2])
	ygrid = cam_coord[2]*tan(linspace(a1, a2, res[1] * self.cam_zoom)) + cam_coord[1]

	from distorsion import  coord_distort
	[xend, yend]  = coord_distort(xgrid, ygrid, center, distorsion)
        xmean = mean(xend); ymean = mean(yend)
	xend, yend = [ (xend-xmean)*cos(alpha_rot) + (yend-ymean)*sin(alpha_rot)+ xmean, (xmean- xend)*sin(alpha_rot) + (yend-ymean)*cos(alpha_rot)+ymean ]

	self.cam_coord = cam_coord;

        zend = zeros(nl * self.cam_zoom**2)  # mriz je v ose

        xstart = ones(nl * self.cam_zoom**2) * cam_coord[0] #+ cam_x
        ystart = ones(nl * self.cam_zoom**2) * cam_coord[1] #+ cam_y
        zstart = ones(nl * self.cam_zoom**2) * cam_coord[2] #+ cam_z

	##################################################
	xend  = xend[::-1, ::-1]
	yend  = yend[::-1, ::-1]
	##################################################
	
	xend = reshape(xend, (-1), order="F")
	yend = reshape(yend, (-1), order="F")
	zend = reshape(zend, (-1), order="F")

	print("saving geometry")
	i = 0
	save(self.geometry_path+'/detector_'+str(i)+'_x.npy', array([xstart, xend ]).T)
	save(self.geometry_path+'/detector_'+str(i)+'_y.npy', array([ystart, yend ]).T)
	save(self.geometry_path+'/detector_'+str(i)+'_z.npy', array([zstart, zend ]).T)
	print("done")
	
	
        self.calb = ones(len(dets))
        self.dets_index = [arange(i*res[0],(i+1)*res[0] ) for i in range(res[1]) ]

        self.tvec = tvec
        self.dets = dets
	self.nl = len(dets)


    def mag_equilibrium(self, tvec, preferCache = True,DryRun=False, return_mean = False):
        REAL_FIELD = True
        if REAL_FIELD:
            Input = loadmat(os.path.normpath(self.geometry_path+'/TSdata_example.mat'), struct_as_record=True)
            magnet_field = array((Input['xc'][0], Input['yc'][0])).T
	    magnet_field /=  1000
            from prepare_data import convert_mag_field
            magx, magy = convert_mag_equilibrium(magnet_field)
        else:
            bd = loadtxt(geometry_path+'/border.txt')
            magx = zeros((len(bd),20))
            magy = zeros((len(bd),20))
            c = mean(bd,0)
            lin = linspace(0,1,20)
            for i in range(20):
                magx[:,i] = (bd[:,0]-c[0])*lin[i]+c[0]
                magy[:,i] = (bd[:,1]-c[1])*lin[i]+c[1]
	    
	self.tsurf = linspace(amin(self.tvec)-0.1, amax(self.tvec)+0.1, 123)

        self.magx = reshape(magx, (size(magx,0), size(magx,1), 1))
        self.magy = reshape(magy, (size(magy,0), size(magy,1), 1))
        if not return_mean:
	    self.magx = self.magx.repeat(len(tvec), axis = 2)
	    self.magy = self.magy.repeat(len(tvec), axis = 2)
	    
	self.magx += self.MagField_shift[0]
	self.magy += self.MagField_shift[1]

	#self.magy -= 0.05

        return self.magx, self.magy


    #def get_data(self):
        #"""
        #Returns loaded data from cache -> prevents memory issues
        #"""
        #data = load(self.geometry_path+'/data_cache.npy')
        #return data

    #def save_data(self, data):
        #print "saving to cache"
        #save(self.geometry_path+'/data_cache.npy', data)
        #print "done"
