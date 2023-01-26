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
from scipy.stats.mstats import mquantiles
import config
from emissivity_generator import emissivity_generator
from tokamak import Tokamak
from scipy.signal import  fftconvolve

class COMPASS(Tokamak):
    def __init__(self, input_diagn, input_parameters,  only_prepare = False  ):

        """ Initialize new tokamak and set default values.

            :param dict input_parameters:   List of selected paramerers
            :param str input_diagn:  Name of main diagnostics
            :param array MagField_shift: artificial shift of magnetic field

        """
        # BASIC PARAMETERS
        name = 'COMPASS'
        #coord = [0.3,0.9,-0.45,0.45]
        coord = [0.2,0.9,-0.45,0.45]

	sigma = 0.0
        min_error = 0.02
        norm = 1.0
        t_name = 'ms'
        self.sigma = 0
	self.index = 2
        self.min_error = 0.01

        #input_diagn = "camera+vis"

	path = 'geometry_compass/'+input_diagn
	Tokamak.__init__(self, input_parameters, input_diagn, coord,  name, norm, sigma, min_error, path,t_name)



	if input_diagn == "SXR":
	    self.SXR()
	    self.index = 2
	if input_diagn == "BOLO":
	    self.Bolo()
	    self.index = 2
	if input_diagn == "camera":
	    self.Camera()
	    self.index = 13
	if input_diagn == "CamVIS":
	    self.CamVIS()
	    self.index = 14
	if input_diagn == "VIS":
	    self.VIS()
	    self.index = 17
	    
	self.tsteps = len(self.tvec)
	self.ports = loadtxt(self.geometry_path + 'ports.txt')
	    
	if not only_prepare:
	    self.prepare_tokamak()




        #print "dokoncena inicializace"
        #exit()

        
    def VIS(self):
	from geom_mat_setting import geom_mat_setting
	self.artif_data = True
	self.geometry_path =os.path.normpath( self.input_path+'geometry_compass/SXR')

	self.Tmat , self.Xchords, self.Ychords = geom_mat_setting(self,self.nx,  self.ny, self.virt_chord)
	self.geometry_path =os.path.normpath( self.input_path+'geometry_compass/VIS')
	tvec = arange(200)/10.0           #time vector
	self.tvec  = tvec
	data , tmp, G_orig = emissivity_generator(self, tvec, self.nx, self.ny, 'Hollow')


	self.dets = arange(size(self.Tmat, 0))

	self.dets_index = list()
	self.dets_index.append(arange(shape(self.Tmat)[0]))

	self.sigma = ones(len(self.dets))*0.05
	self.error = 0.01*amax(data)
	data +=  self.error*random.randn(size(data,0),size(data,1) )
	data = data.T

	self.calb = ones(len(self.dets_index))
	self.save_data( data.T )


	
    def CamVIS(self):
	from geom_mat_setting import geom_mat_setting
	self.artif_data = True

	usecache = config.useCache
	config.useCache = False
	self.geometry_path =os.path.normpath( self.input_path+'geometry_compass/SXR')
	Tmat_VIS , Xchords_VIS, Ychords_VIS = geom_mat_setting(self,self.nx,  self.ny, self.virt_chord)
	config.useCache = usecache

	#exit()
	
	self.geometry_path =os.path.normpath( self.input_path+'geometry_compass/camera')
	self.Camera()
	Tmat_cam , Xchords_cam, Ychords_cam = geom_mat_setting(self,self.nx,  self.ny, self.virt_chord)

############################
	self.Tmat =  sparse.csc_matrix( sparse.vstack([Tmat_cam, Tmat_VIS]) )
	#print shape( Xchords_cam), shape(Xchords_VIS )

	Nc = size(Xchords_cam,0)  # number of points used for description of chord
	from scipy import interpolate
		
	Xchords_VIS = interpolate.interp1d(linspace(0,1,size(Xchords_VIS,0)), Xchords_VIS.T)(linspace(0,1,Nc)).T
	Ychords_VIS = interpolate.interp1d(linspace(0,1,size(Ychords_VIS,0)), Ychords_VIS.T)(linspace(0,1,Nc)).T

	#print shape( Xchords_cam), shape(Xchords_VIS )

	self.Xchords = vstack([ Xchords_cam.T, Xchords_VIS.T ]).T
	self.Ychords = vstack([ Ychords_cam.T, Ychords_VIS.T ]).T

	#self.Xchords = self.Xchords_VIS
	#self.Ychords = self.Ychords_VIS

	#savez('chords', Xchords = self.Xchords, Ychords = self.Ychords )
	#print "saved"

	#exit()
	
	##img_smaller_36.png

	#print self.Tmat, shape(self.Tmat), type(self.Tmat)
	#exit()
	
	#print self.Tmat
	self.camera = False # neni to cista kamera ...

	#print shape(self.Tmat ), shape(Tmat_cam), shape(Tmat_VIS)
	self.geometry_path =os.path.normpath( self.input_path+'geometry_compass/CamVIS')

	tvec = arange(200)/10.0           #time vector
	self.tvec  = tvec
	data , tmp, G_orig = emissivity_generator(self, tvec, self.nx, self.ny, 'Hollow')
	


	self.dets = arange(size(self.Tmat, 0))

	self.dets_index = list()
	self.dets_index.append(arange(shape(Tmat_cam)[0]))
	self.dets_index.append(1+amax(self.dets_index[-1])+ arange(shape(Tmat_VIS)[0] ))

	sigma = ones(len(self.dets))
	sigma[self.dets_index[0]] = 0.05
	sigma[self.dets_index[1]] = 0.02
	#print sigma
	#exit()
	self.sigma = sigma
	self.error = 0.01*amax(data)
	data +=  self.error*random.randn(size(data,0),size(data,1) )
	data = data.T
	
	self.calb = ones(len(self.dets_index))
	self.save_data( data.T)
	self.tvec = tvec

	#print self.nx, self.ny

	#exit()
	
	
    def SXR(self):


	if self.artif_data:
	    #generate artificial data
	    tvec = arange(200)/10.0           #time vector
	    self.tvec  = tvec
	    data , tmp, G_orig = emissivity_generator(self, tvec, self.nx, self.ny, 'Hollow')

	    self.error = 0.02*amax(data)
	    data +=  self.error*random.randn(size(data,0),size(data,1) )
	    data = data.T

	else:
	    data_path = self.geometry_path+'/SXR_'+str(self.shot)+'.npz'

	    try:
		d = load(data_path)
		data = d['data']
		tvec = d['tvec']
		dets = arange(size(data, 1))
		del d
		print("loading done")
	    except:
		from scipy import interpolate
		from scipy.io import loadmat, savemat
		d = list()
		for i in ['SAA','SBA']:
		    print(i)
		    d.append(loadmat(self.input_path + 'geometry_compass/Data/'+str(self.shot)+'/'+i+'raw.mat')['sig'])

		data = array(d, copy=False)
		tvec = squeeze(loadmat(self.input_path + 'geometry_compass/Data/'+str(self.shot)+'/'+i+'raw.mat')['time'])

		data = reshape(data, (-1, len(tvec))).T

		#pcolor(data[::10000,:])
		#show()
	
		print("smoothing")
		dets = arange(size(data, 1))
		decimate = 100
		for i in dets:
		    data[:,i] = fftconvolve(data[:,i], ones(decimate)/decimate, mode='same')
		data = data[::decimate,:]
		tvec = tvec[::decimate]

		min_data = zeros(len(dets))
		for i in range(len(dets)):
		    min_data[i] = mquantiles(data[(tvec < 0.93) & ~isnan(data[:,i]),i], 0.2)

		#plot(min_data)
		#show()

		#exit()
		
   
		ind_wrong = where( (min_data < 4) | (min_data > 5) )[0]
		ind_wrong = int_(array( [ ind_wrong, [48] ]))
		
		data -= min_data
		data[:,ind_wrong] = NaN


		savez(data_path , data = single(data), tvec = tvec, dets = dets)



	
	self.dets = arange(size(data, 1))

	self.dets_index = list()
	self.dets_index.append(arange(35))
	self.dets_index.append(1+amax(self.dets_index[-1])+ arange(35))

	self.calb = ones(len(self.dets_index))
	self.save_data( data.T)
	self.tvec = tvec


    def Camera(self):
        print("============= camera =============")
        self.camera = True
        self.allow_self_calibration = False
        self.wrong_dets_pref = []
        self.min_error = 0.01
        self.sigma = 0.00
        #self.allow_divertor = False
        #self.smooth_boundary = False
        self.virt_chord = 1
        self.expect_zero = True
        self.weight_power  = 4
        self.allow_reflections = True

	self.boundary_full = loadtxt( self.geometry_path  +'/border_all.txt')
	self.boundary_simple = loadtxt( self.geometry_path  +'/border_simple.txt')

        #self.geometry_path =  os.path.normpath(self.input_path+'geometry_camera')

        #for i in os.listdir(self.geometry_path):
          #if 'detector_' in i:
           #os.remove(self.geometry_path+'/'+i)

        #print "artificial data"
        #self.artif_data = True
        # artificial !!!!

        tvec = linspace(50, 52, 10)

	def rot270(x):
	    return rot90(rot90(rot90(x)))

	if not self.artif_data:
	    print("************ loading data from file **************")
	    data  = imread( self.geometry_path + '/img_smaller_36.tif')
	    #imshow(data)
	    #show()
	    
	    #shot = 2858
	    #frame = 93
	    shot = 3487
	    frame = 13
	    
	    #2752    54
	    #2753
	    #2754 45

	    #shot =  2858  # 94 95 

	    #3068 63

	    #3135 - 48  plasma i s okrajema 

	    #3496  90 - D plasma 

	    #img_small   9735
	    #shot = 3496 ; frame =  90

	    data_all  =  loadmat(self.geometry_path+'/data/'+str(self.shot)+"_small.mat")['imgs']
	    
	    [H,W,N] = shape(data_all)
	    res = [W,H]
	    
	    tvec =  arange(size(data_all,2))
	    
	    data = zeros((W*H,N))
	    for i in tvec:
		#print i
		d = ((data_all[:,:,i]))
		#imshow(d, interpolation="nearest", cmap="gray")
		#savefig(self.geometry_path+'/data/'+str(shot)+"_"+str(i)+".png")
		data[:,i] = ravel(rot270(d))
		
		
		
	    
	    
	    #imshow(data)
	    #show()
	    
	    #data =  double(rot270(data))

	    FOV = rot270(load( self.geometry_path + '/fov.npy' ))
	    ports_area = rot270(load( self.geometry_path + '/ports.npy' ))
	    limiter = rot270(load( self.geometry_path + '/limiter.npy' ))

	    #FOV[:50,:] = 0

	    #imshow(FOV, interpolation = "nearest")
	    #show()
	    
	    background = imread(self.geometry_path + '/img_background.png')
	    self.background = double(ravel(rot270(background[::2,::2])))
	    
	    #print shape(background)
	    #imshow(background, interpolation = "nearest", cmap="gray")
	    #show()
	    
	    

	    #print shape(FOV)
	    #exit()
	    

	    self.FOV = FOV != 0
	    self.ports_area = ports_area != 0
	    self.limiter = limiter != 0

	    #res = shape(data)
	    
	    #print self.ports
	    #eeee
	    #data = load('data_compass.npy')
	    #data = reshape(data[:,0], res)
	    
	    #print shape(data)
	    #print 
	    #d = reshape(data, res)
	    #print shape(d), shape(data)
	    ##d[FOV] = nan
	    #subplot(1,2,1)
	    #data[ports] = 0
	    #data[limiter] = 0
	    #imshow(data)
	    #show()
	    
	    #colorbar()
	    #subplot(1,2,2)
	    #imshow(FOV)
	    #colorbar()
	    ##subplot(1,3,2)
	    ##bg = load('diff_model.npy')
	    ##imshow(bg)
	    ##colorbar()
	    #show()
	    
	    config.wrong_dets_pref = where(ravel(FOV) == 0)[0]
	    
	    #print "============================"
	    
	    #exit()
	    
	    
	    config.camera_res = res
	    
	    #data = double(ravel(data))
	    #data[isnan(data)]  = 0
	    
	    #data = data[:,None].repeat(len(tvec), axis=1)
	    self.save_data( data )
	    del data 

        # artificial !!!!


        #coord = [-0.2,-0.9,-0.45,0.45]# real coordinates

        # fungují  paramerery !!!!!!!!!!!!!!!!!!!!
        #coord = [-0.7, -0.2 ,-0.12,0.4]# camera area, z = 0


        cam_size = [5,5]
        
	try:  # use the value in config if availible
	    res = array(config.camera_res)
	except:
	    #res = array( (124, 113) ) # (124, 113) )  # mozna to rozliseni musi byt opacene ??
	    
	    res = array( (40,40) )
	    
	    
	#res = array( (992/4, 900/4)) # (124, 113) )  # mozna to rozliseni musi byt opacene ??

	self.camera_res = res


        nl = prod(res)
	dets = arange(nl)
	    #xxx

        #xgrid = linspace(coord[1], coord[0], res[0] * self.cam_zoom)
        #ygrid = linspace(coord[3], coord[2], res[1] * self.cam_zoom)

        #coord = [-0.66, -0.16 ,-0.12,0.38]# camera area, z = 0
        #coord = [-0.87, -0.0 ,-0.25,0.40]# camera area, z = 0  # 4.2.2013
	#coord = [-0.7, -0.12 ,-0.13,0.33]  # 19.4.2013, asi nefunguje 
	
        #coord =  [-0.7, -0.11 ,-0.13,0.33]    # camera area, z = 0  # new 8.5 
        
	coord =  [-0.695, -0.105 ,-0.185,0.335]    # camera area, z = 0  # new

	#cam_coord = [-0.75/sqrt(2), 0.18, 0.75/sqrt(2)]  #  !!  funguje
	#cam_coord = [-0.8/sqrt(2), 0.17, 0.8/sqrt(2)]  #  2. verzze co funguje 
	#cam_coord = [-0.75/sqrt(2), 0.19, 0.75/sqrt(2)]  # x,y,z # 4.2.2013
	#cam_coord = [-0.75/sqrt(2), 0.15, 0.75/sqrt(2)]  # 19.4.2013, asi nefunguje 
	cam_coord = [-0.75/sqrt(2), 0.20, 0.75/sqrt(2)]  # x,y,z # 4.2.2013


        distorsion = 0.5 # 19.4.2013, asi nefunguje 
	distorsion = 0.0 # 
	
	#alpha_rot =  -0.12
	alpha_rot =  -0.02  # new 
	alpha_rot =  -0.1  # new 

	center = array([-0.4491 ,  0.25]) # distorsion center
	self.cam_zoom = 1  # 4
	self.pixel_zoom = 1 # 4



	
	#dist = cam_coord[2]   # distance of camera from FOV place
	#FOV_center = [-0.4, 0.15] #[-0.47, 0.1]

	
	#xgrid = linspace(coord[1], coord[0], res[0] * self.cam_zoom)

	#print xgrid
	a1 = arctan2(coord[1] - cam_coord[0], cam_coord[2])
	a2 = arctan2(coord[0] - cam_coord[0], cam_coord[2])
	#print a1, a2
	xgrid = cam_coord[2]*tan(linspace(a1, a2, res[0] * self.cam_zoom)) + cam_coord[0]
	#print xgrid2
	#plot(xgrid)
	#plot(xgrid2)
	#show()

        #ygrid = linspace(coord[3], coord[2], res[1] * self.cam_zoom)

	a1 = arctan2(coord[3] - cam_coord[1], cam_coord[2])
	a2 = arctan2(coord[2] - cam_coord[1], cam_coord[2])
	ygrid = cam_coord[2]*tan(linspace(a1, a2, res[1] * self.cam_zoom)) + cam_coord[1]

	#print a1, a2
	#print ygrid
	#print ygrid2

	#plot(ygrid)
	#plot(ygrid2)
	#show()

	#exit()
	
        #print xgrid
        #print ygrid
	#alpha_max = 35. / 180.*pi
        #xgrid = FOV_center[0]  + dist * sin(linspace(alpha_max, -alpha_max, res[0] * self.cam_zoom))
	#alpha_max = 35. / 180.*pi
	#ygrid = FOV_center[1]  + dist * sin(linspace(alpha_max, -alpha_max, res[1] * self.cam_zoom))
	#print xgrid
	#print ygrid
	#exit()


        # !!! !!!!!! udelat to otoceni korektne okolo stredu toho FOV !!!! 
        #[xend, yend] = meshgrid(xgrid, ygrid)
        #alpha_rot = -0.15 #  pootoceni kamery !!  funguje !!!!!!


	#center = array([mean(xgrid), mean(ygrid)])
	#print center
	#exit()

	from distorsion import  coord_distort
	[xend, yend]  = coord_distort(xgrid, ygrid, center, distorsion)
	#  pootoceni kamery !!  funguje !!!!!!
	#plot(xend, yend, 'b.')
        #xend, yend = [ xend*cos(alpha_rot) + yend*sin(alpha_rot), - xend*sin(alpha_rot) + yend*cos(alpha_rot) ]
        
        
        xmean = mean(xend); ymean = mean(yend)
	xend, yend = [ (xend-xmean)*cos(alpha_rot) + (yend-ymean)*sin(alpha_rot)+ xmean, (xmean- xend)*sin(alpha_rot) + (yend-ymean)*cos(alpha_rot)+ymean ]

	#plot(xend, yend, 'r,')
	#show()


	self.cam_coord = cam_coord;

        zend = zeros(nl * self.cam_zoom**2)  # mriz je v ose

        xstart = ones(nl * self.cam_zoom**2) * cam_coord[0] #+ cam_x
        ystart = ones(nl * self.cam_zoom**2) * cam_coord[1] #+ cam_y
        zstart = ones(nl * self.cam_zoom**2) * cam_coord[2] #+ cam_z

	##################################################
	#print shape(xend)
	xend  = xend[::-1, ::-1]
	yend  = yend[::-1, ::-1]
	#zend  = zend[::-1, ::-1]
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


	#self.FOV = coord
	
        self.tvec = tvec
        self.dets = dets
	self.nl = len(dets)

	#plot(array([xstart, xend ]), array([zstart, zend ]), 'ko-')
	#show()
	#plot(array([ystart, yend ]), array([zstart, zend ]), 'ko-')
	#show()
	#exit()


	    
    def Bolo(self):

	#self.artif_data = False
	
	if self.artif_data:
	    #generate artificial data

	    tvec = arange(200)/10.0           #time vector
	    self.tvec = tvec
	    
	    #data , tmp, G_orig = emissivity_generator(self, tvec , self.nx, self.ny, 'Hollow')
	    self.nl = 120
	    dets = arange(self.nl)
	    
	    #self.error = 0.02*amax(data)
	    #data +=  self.error*random.randn(size(data,0),size(data,1) )
	    #data = data.T

	    #exit()
	    

	else:

	    data_path = self.geometry_path+'/Bolo_'+str(self.shot)+'.npz'

	    try:
		#xxx
		d = load(data_path)
		data = d['data']
		tvec = d['tvec']
		dets = arange(size(data, 1))
		del d
		print("loading done")
	    except:

		from scipy.io import loadmat, savemat
		d = list()
		for i in ['BAC','BBN', 'BCC', 'BDN', 'BEC','BFC']:
		    print(i)
		    d.append(loadmat(self.input_path + 'geometry_compass/Data/'+str(self.shot)+'/'+i+'raw.mat')['sig'])

		data = array(d, copy=False)
		tvec = squeeze(loadmat(self.input_path + 'geometry_compass/Data/'+str(self.shot)+'/'+i+'raw.mat')['time'])

		data = reshape(data, (-1, len(tvec))).T

		print("smoothing")
		dets = arange(size(data, 1))
		decimate = 100
		for i in dets:
		    data[:,i] = fftconvolve(data[:,i], ones(decimate)/decimate, mode='same')
		data = data[::decimate,:]
		tvec = tvec[::decimate]

		min_data = zeros(len(dets))
		for i in range(len(dets)):
		    min_data[i] = mquantiles(data[(tvec < 0.93) & ~isnan(data[:,i]),i], 0.2)
		ind_wrong = where((min_data < 0.2) | (min_data > 0.5))[0]
		ind_wrong = int_(concatenate( [ ind_wrong, array([60,96,12,75,25,32,99,119 ]) ]))
		#ind_wrong = int_(concatenate( [ ind_wrong, arange(80,100)  ]))
		#ind_wrong = int_(concatenate( [ ind_wrong, arange(60,80)  ]))


		
		data -= min_data
		data[:,ind_wrong] = NaN

		
		savez(data_path , data = single(data), tvec = tvec, dets = dets)

	    self.save_data( data.T )

	#pcolor(data)
	#show()

	#plot(data.T)
	#show()
	#exit()

	#print shape(data)
	

	self.tvec = tvec 
        self.dets = dets

        #exit()
        

	self.dets_index = list()
        self.dets_index.append(arange(20))
        self.dets_index.append(1+amax(self.dets_index[-1])+ arange(20))
        self.dets_index.append(1+amax(self.dets_index[-1])+ arange(20))
        self.dets_index.append(1+amax(self.dets_index[-1])+ arange(20))
        self.dets_index.append(1+amax(self.dets_index[-1])+ arange(20))
        self.dets_index.append(1+amax(self.dets_index[-1])+ arange(20))

	#print shape(self.dets_index)
	#print self.dets_index

	self.calb = [1,1,1,0.68,0.12, 0.1]

	self.default_calb = "none"
	


    def load_mag_equilibrium(self):
        magnet_field =  loadtxt(os.path.normpath(self.geometry_path+'/magnet_field.txt'))
        from prepare_data import convert_mag_field
        magx, magy = convert_mag_equilibrium(magnet_field)
        #magx = magx[:,::-1]
	#magy = magy[:,::-1]

	self.tsurf = linspace(amin(self.tvec)-0.1, amax(self.tvec)+0.1, 100)
	
        self.magx = reshape(magx, (size(magx,0), size(magx,1), 1))
        self.magy = reshape(magy, (size(magy,0), size(magy,1), 1))
        self.magx = self.magx.repeat(len(self.tsurf), axis = 2)
        self.magy = self.magy.repeat(len(self.tsurf), axis = 2)

