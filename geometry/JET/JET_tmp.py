#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
import os, re
from numpy import *
import matplotlib
from scipy import sparse
import sys
from scipy.io import loadmat, savemat
from matplotlib.pyplot import *
from numpy.matlib import repmat
from shared_modules import in1d
from orthogonal_trans import *
import time
from tokamak import Tokamak
from scipy.stats.mstats import mquantiles
import config

class JET(Tokamak):
    """
    Representation of class of tokamak. Contain basic information about geometry, borde, magnetic field and basic setting used for recontruction and post_processing.

    >>> JET(input_path, 'SXR-slow',  shot, nx, ny,virt_chord )


    :var str name:  Name of tokamak
    :var int index: Index of the tokamak
    :var double xmin, xmax, ymin, ymax: boundarys of recontructed area of tokamak emissivity
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
    """


    def __init__(self, input_diagn, input_parameters,  only_prepare = False  ):

        """ Initialize new tokamak and set default values.

            :param dict input_parameters:   List of selected paramerers
            :param str input_diagn:  Name of main diagnostics
            :param array magfield_shift: artificial shift of magnetic field

        """
	
	
        name  = 'JET'
	path_tmp = re.sub( '(.+)[-]+(.+)', r'\1',  input_diagn)
        if input_diagn in ["SXR-slow","SXR-fast", 'camera'] :
	    coord= [185.0,400, -145.0,195.0]
	    coord= [180.0,400, -145.0,195.0]
	    #path ='geometry_JET'
	elif input_diagn in ['bolo']:
	    coord = [180.0,397.2,-190.0,210.0]
	    #path = 'geometry_JET_bolo'
	elif input_diagn in ['neutrons']:
	    coord = [202.2,392.2,-140.0,190.0]

	#path = '/work/mimrisek/geometry/JET/SXR'
	path = 'geometry/JET/'+path_tmp


	sigma = 0.0
        min_error = 0.005
        norm = 100.0
        t_name = 's'
        self._new_calb_shot = 80000  # new SXR diagn. is used from this shot number 
        self._new_diagn_shot = 83900  # new SXR diagn. is used from this shot number 

        #######################################
        self.local_network = os.system(' [ `domainname` == "jaclinux" ] ') == 0  # !!!! True if you are in JET local network !!!
        #######################################

	Tokamak.__init__(self, input_parameters, input_diagn, coord,  name, norm, sigma, min_error, path,t_name)


	#if self.shot > self._new_diagn_shot:
        if self.shot > self._new_calb_shot:
            if input_diagn == "SXR-slow":
                self.SXR('new_slow')
                self.index = 4
            if input_diagn == "SXR-fast":
                self.SXR('new_fast')
                self.index = 5
        else:
            if input_diagn == "SXR-slow":
                self.SXR('slow')
                self.index = 4

            if input_diagn == "SXR-fast":
                self.SXR('fast')
                self.index = 5

	if input_diagn == "neutrons":
	    self.Neutrons('DD_ppf') #DD,DD_ppf,DT,DT_pf
	    self.index = 7
	if input_diagn == "bolo":
	    self.Bolo()
	    self.index = 11
	if input_diagn == "camera":
	    self.Camera()
	    self.index = 12

        #print 'JET 121'
        #import code
        #code.interact(local=dict(globals(), **locals()))

	Tokamak._post_init(self)

	self.divertor_coord  = [210,330,-200,-80]


	if not only_prepare:
	    self.prepare_tokamak()

	print("object JET done")



    def SXR(self, SOURCE):
        """ Main diagnostics is Soft X Rays, try to load data from local cache and if it is not possible, try to load dta from JET network.

        :param str SOURCE: Choose if the data should be loaded from fast or slow data acquisition system

        :var list calb_0: Effeciency calibration for each sensor, it is list of integers or list of vectors with length == tsteps
        :var int time_interval: Time block used to load data from JET from fast DAQ
        :raises Asseration: No connection to JET network"
        """

        #self.xmin=190.0
        #self.xmax=397.2
        #self.ymin=-145.0
        #self.ymax=195.0

        ##self.xmin=100.0
        #self.xmax=400.2
        #self.ymin=-200.0
        #self.ymax=250.0

        #self.geometry_path =os.path.normpath(self.input_path+'geometry_JET')

        #self.calb_0 = [1.1, 1, 0.78]  # relative calibration, used in prepare_data()
        self.default_calb = 'variable'

        if not self.artif_data:
	    time_interval = 0

	    if self.shot < self._new_calb_shot:
		calbV = loadtxt(self.geometry_path+'/calbV.txt').T  #*1.1
		calbS4 =loadtxt(self.geometry_path+'/calbS4.txt').T  #*1; #0.78 # 0.78
		calbT =loadtxt(self.geometry_path+'/calbT.txt').T   #*0.78 *NaN#0.71
	    else:
		print("!!!!!!!! nepoužívám vzájemnou kalibraci !!!!!!!!!!!")
		calbV = ones(35)
		calbS4 = ones(17)
		calbT = ones(35)

	    calb_all = concatenate([calbV,calbS4,calbT ])
	    self.dets_index = list()
	    self.dets_index.append(arange(len(calbV)))
	    self.dets_index.append(1+amax(self.dets_index[-1])+ arange(len(calbS4)))
	    self.dets_index.append(1+amax(self.dets_index[-1])+ arange(len(calbT)))

	    python_data = self.geometry_path+'/SXR_Data_'+str(self.shot)+'_'+str(time_interval)+'_'+SOURCE+'.npz'
	    matlab_data =  self.geometry_path+'/'+SOURCE+str(self.shot)+'.mat'

	    print(python_data)
	    print(matlab_data)

	    try:
		print("loading", SOURCE)
		if os.path.isfile(python_data):
		    d = load(python_data)
		    data = d['data']
		    tvec = d['tvec']
		    dets = arange(size(data, 0))
		    del d
		    print("loading done")
		elif os.path.isfile(matlab_data):
		    dets, tvec, data =  self.get_mat_data(matlab_data)
		else:
		    assert False,  'Missing data'

	    except:
		#raise
		try:
		    assert  self.local_network, "cannot load data"
		    from .mds_load_sxr import mds_load_sxr
		    tvec, data = mds_load_sxr(self.shot, self.geometry_path, SOURCE, time_interval)
		except:
		    raise

		data = data.T
		dets = arange(size(data, 0))


		try:
		    savez(python_data , data = single(data), tvec = tvec, dets = dets)
		except:
		    print("data saving failed !!!!")


	    data = data.T
	    data *= sparse.spdiags(calb_all, 0, len(calb_all),len(calb_all) )
	    data = data.T


	    #print 'ODSTRANENA KALIBRACE SENZORŮ !!!'
	    #self.wrong_dets_damaged = array([8,46, 57,59, 78])-1        # problems ???
	    if self.shot < self._new_calb_shot:
		self.wrong_dets_damaged = array([44,53,55, 57, 78,81,83,85,87])-1        # problems ???
                if self.shot > 69541 & self.shot < 69943:
		    self.wrong_dets_damaged = array([12,44,53,55, 57, 78,81,83,85,87])-1        # problems ???

	    else:
		self.wrong_dets_damaged = array([33,34, 35, 47])-1



	    #########################################################################################
	    # !!! SETUP OF THE WRONG DETECTORS !!!!
	    config.wrong_dets_pref =  self.dets_index[2]   # detektor T0
	    #self.wrong_dets_pref = concatenate([self.wrong_dets_pref , self.dets_index[0]])   # detektor V0


	    #########################################################################################


	    d = nansum(data[~in1d(dets,config.wrong_dets_pref ),:], 0)
	    emiss_index = d > 0.01*amax(d)
	    data = data[:,emiss_index]

	    tvec = tvec[emiss_index]

	    self.save_data( data,tvec)
	    #self.dets_index = [arange(32), 32+arange(24)]
            self.dets_index = [arange(35), 35+arange(17), 35+17 + arange(35)]

	else:
	    # artificial !!!!
	    dets = arange(87)
	    tvec = linspace(50, 52, 100)
	    #self.dets_index = [arange(32), 32+arange(24)]
            self.dets_index = [arange(35), 35+arange(17), 35+17 + arange(35)]
	    
	
	self.calb_0 = [1,1, 1]

	    # artificial !!!!




        self.tvec = tvec
        self.dets = dets


    def get_mat_data(self, path):


        from scipy import interpolate
        from scipy.io import loadmat, savemat

        d = loadmat(path)

        print("Loading data from .mat file")



        sigS4 = interpolate.interp1d(squeeze(d['trawsig2']) , d['rawsigS'].T, bounds_error=False, fill_value=0)( squeeze(d['trawsig']))

        sigT = d['rawsigT'].T
        sigV = d['rawsigV'].T

        data = vstack([sigV, sigS4, sigT])
        tvec = squeeze(d['trawsig'])

        dets = arange(size(data, 0))

        return dets, tvec, data

    def Camera(self):
        print("============= camera =============")
        self.camera = True
        self.allow_self_calibration = False
        self.wrong_dets_pref = []
        self.min_error = 0.0001
        self.sigma = 0.00
        self.allow_divertor = False
        #self.smooth_boundary = False
        self.virt_chord = 1

        #self.geometry_path =  os.path.normpath(self.input_path+'geometry_camera')

        #for i in os.listdir(self.geometry_path):
          #if 'txt' in i:
           #os.remove(self.geometry_path+'/'+i)

        print("artificial data")
        self.artif_data = True
        # artificial !!!!

        tvec = linspace(50, 52, 100)
        # artificial !!!!

	#coord= [190.0,397.2,-145.0,195.0]
	#coord = [180.0,397.2,-190.0,210.0]  # real coordinates

        coord= [100.,430.,-120.0,140.0]
        cam_coord = [400/sqrt(2), 0, 400/sqrt(2)]  # x,y,z
        cam_size = [5,5]
	self.cam_coord = cam_coord;

        res = array( [60,60] )
	self.camera_res = res

	self.cam_zoom = 1 #8
	self.pixel_zoom = 1 # 4

        nl = prod(res)
        dets = arange(nl)

        xgrid = linspace(coord[0], coord[1], res[0] * self.cam_zoom)
        ygrid = linspace(coord[2], coord[3], res[1] * self.cam_zoom)
        
        [xend, yend] = meshgrid(xgrid, ygrid)
        alpha = 0.005   # pootoceni kamery !! 
        xend, yend = [ xend*cos(alpha) + yend*sin(alpha), - xend*sin(alpha) + yend*cos(alpha) ]

        zend = zeros(nl * self.cam_zoom**2)  # mriz je v ose

	print("POZOR UPGRADE ")
        xstart = ones(nl * self.cam_zoom**2) * cam_coord[0] #+ cam_x
        ystart = ones(nl * self.cam_zoom**2) * cam_coord[1] #+ cam_y
        zstart = ones(nl * self.cam_zoom**2) * cam_coord[2] #+ cam_z

	xend = ravel(xend)
	yend = ravel(yend)
	zend = ravel(zend)

	i = 0
	savetxt(self.geometry_path+'/detector_'+str(i)+'_x.txt', array([xstart, xend ]).T, fmt='%5.4f')
	savetxt(self.geometry_path+'/detector_'+str(i)+'_y.txt', array([ystart, yend ]).T, fmt='%5.4f')
	savetxt(self.geometry_path+'/detector_'+str(i)+'_z.txt', array([zstart, zend ]).T, fmt='%5.4f')

          #savetxt(self.geometry_path+'/dist_'+str(i)+'.txt', ravel(xgrid[i])  , fmt='%5.4f')

        self.calb_0 = ones(len(dets))
        self.dets_index = [arange(i*res[0],(i+1)*res[0] ) for i in range(res[1]) ]

        #print dets, self.dets_index

        #exit()

        self.tvec = tvec
        self.dets = dets
	self.nl = len(dets)



    def Bolo(self):
        """ Main diagnostics are Bolometers,
        """

	print("="*20, "loading bolometry ", "="*20)



        #self.virt_chord = 1

        self.allow_self_calibration = False


        #self.wrong_dets = []
        self.wrong_dets_pref = []
        self.min_error = 0.015
        self.sigma = 0.03

        self.allow_divertor = True
        self.smooth_boundary = False


	###############################
        #self.artif_data = False
        ###############################


        python_data =  self.geometry_path+ '/BOLO_'+str(self.shot)+'.npz'
        matlab_data =  self.geometry_path+'/BOLO'+str(self.shot)+'.mat'


        if not self.artif_data:
	    #try:
	    print("loading")
	    if os.path.isfile(python_data):
		d = load(python_data)
		data = d['data']
		tvec = d['tvec']
		dets = arange(size(data, 1))
		del d
		print("loading done")
	    elif os.path.isfile(matlab_data):
		dets, tvec, data =  self.get_mat_data(matlab_data)
	    else:
		#assert False,  'Missing data'
	    #except:
		from .mds_load_bolo import mds_load_bolo
		tvec, data = mds_load_bolo(self.shot, self.geometry_path, 'KB5', "")

		dets = arange(size(data, 1))
		try:
		    savez(python_data , data = single(data), tvec = tvec, dets = dets)
		except:
		    print("data saving failed !!!!")



	    #print "==================================================================="
	    #print shape(data), shape(~in1d(dets,self.wrong_dets_pref ))

	    #print shape(data[:,~in1d(dets,self.wrong_dets_pref )])

	    # remove too low emissivity data
	    d = data[:, ~in1d(dets, config.wrong_dets_pref )]
	    d = nansum(d , axis = 1)

	    from scipy.signal import medfilt

	    Nfilt = max(1,2*floor(len(d)/500.0)+1)

	    ##d = medfilt(d, Nfilt)
	    ##d = convolve(d, ones(Nfilt)/Nfilt, mode='same')
	    #emiss_index = (d > 0.01*mquantiles(d[~isnan(d)],0.99))
	    #data = data[emiss_index,:]

	    #tvec = tvec[emiss_index]

	    #plot(d)
	    #show()

	    calbV = loadtxt(self.geometry_path+'/calbV.txt').T  #*1.1
	    calbH = loadtxt(self.geometry_path+'/calbH.txt').T  #*1; #0.78 # 0.78

	    calb_all = concatenate([calbV,calbH])

	    #data *= sparse.spdiags(calb_all, 0, len(calb_all),len(calb_all) )


	    #ind  = concatenate((arange(14) , 16 +arange(56-16)))
    	    #ind  = ((arange(18, 32)))

	    #pcolor(data[::100, :  ])
	    #colorbar()
	    #show()
	    #exit()

	    self.calb_error = 1/calb_all / median(1/calb_all)

	    ind_ok = (tvec > 45) & (tvec < 60)  #   & (d >= 0)

	    data -=  median(data[ tvec < 45,:],0)   # remove systematic bias


	    #plot(std(data[ tvec > 65,:],0))
	    #plot(mean(data[ tvec > 65,:],0))

	    #print tvec

	    tvec = tvec[ind_ok]
	    data = data[ind_ok,:]

	    #plot(mean(data,0))
	    #plot(std(data,0))

	    #plot(data[:,50:])


	    #show()
	    #exit()



	    std_data = std(data, axis = 0)

	    self.wrong_dets_damaged =  where( (mean(isnan(data), axis=0) == 1) |  (std_data >  3*mquantiles(std_data,0.9) ) )[0].tolist()

	    self.wrong_dets_damaged.append(7)  # divne hodnoty !!

	    data = data.T  # +     nanmin(data)

	    self.save_data( data,tvec)
	else:
	    print("artificial data")
	    # artificial !!!!
	    dets = arange(56)
	    tvec = linspace(50, 52, 100)
        # artificial !!!!


        self.calb_0 = [1,1]
	self.dets_index = [arange(32), 32+arange(24)]

        self.tvec = tvec
        self.dets = dets
	self.nl = len(dets)


    def Neutrons(self,SOURCE):
        """ Main diagnostics are Neutrons, try to load data only from local cache, loading from network is not supported.
        """


        #self.xmin = 202.2     #must be in the same units as geometry matrix
        #self.xmax = 392.2
        #self.ymin = -140
        #self.ymax = 190
        #geometry = 'geometry_JET_neutrons'
        #self.geometry_path =os.path.normpath( self.input_path+geometry)


        #self.default_calb = 'variable'

        if not self.artif_data:
	    time_interval = 0

            if SOURCE in ['DD_ppf','DT_ppf']:
                calb_all = ones(19)*1e-9
	    else:
                detangle = loadtxt(self.geometry_path+'/detangle.txt').T
                attens = loadtxt(self.geometry_path+'/attens.txt').T
                effic = loadtxt(self.geometry_path+'/effic.txt')
                calb_all = 1/(effic*attens*detangle)*1e-9

            #import code
            #code.interact(local=dict(globals(), **locals()))
		
	    self.dets_index = list()
	    self.dets_index.append(arange(len(calb_all)))
	
	    python_data = self.geometry_path+'/NEUT_Data_'+str(self.shot)+'_'+str(time_interval)+'_'+SOURCE+'.npz'
	    matlab_data = self.geometry_path+'/'+SOURCE+str(self.shot)+'.mat'

	    print(python_data)
	    print(matlab_data)

	    try:
		print("loading", SOURCE)
		if os.path.isfile(python_data):
		    d = load(python_data)
		    data = d['data']
		    tvec = d['tvec']
		    dets = arange(size(data, 0))
		    del d
		    print("loading done")
		elif os.path.isfile(matlab_data):
		    dets, tvec, data =  self.get_mat_data(matlab_data)
		else:
		    assert False,  'Missing data'

	    except:
		#raise
		try:
		    assert  self.local_network, "cannot load data"
		    from mds_load_neut import mds_load_neut
		    tvec, data = mds_load_neut(self.shot, self.geometry_path, SOURCE, time_interval)
		except:
		    raise

                #import code
                #code.interact(local=dict(globals(), **locals()))

		data = data.T
		dets = arange(size(data, 0))


		try:
		    savez(python_data , data = single(data), tvec = tvec, dets = dets)
		except:
		    print("data saving failed !!!!")

           

	    data = data.T
	    data *= sparse.spdiags(calb_all, 0, len(calb_all),len(calb_all) )
	    data = data.T


            self.wrong_dets_damaged = array([]) #array([8,10,11])-1



	    #d = nansum(data[~in1d(dets,config.wrong_dets_pref ),:], 0)
	    #emiss_index = d > 0.01*amax(d)

	    #data = data[:,emiss_index]

	    #tvec = tvec[emiss_index]

	    self.save_data( data,tvec)
	    #self.dets_index = [arange(32), 32+arange(24)]
            self.dets_index = [arange(10),10+arange(9)]

	else:
	    # artificial !!!!
	    dets = arange(19)
	    tvec = linspace(50, 52, 100)
	    #self.dets_index = [arange(32), 32+arange(24)]
            self.dets_index = [arange(10),10+arange(9)]
	    
	
	self.calb_0 = [1,1]

	    # artificial !!!!

        self.tvec = tvec
        self.dets = dets



        #try:
            #Input = loadmat(os.path.normpath(self.geometry_path+'/'+str(self.shot)+'/JETdata.mat'), struct_as_record=True)
        #except:
            #raise IOError("Can't load", None, os.path.normpath(self.geometry_path+'/'+str(self.shot)+'/JETdata.mat'))
        #self.data = Input['Y'].T
        #self.dets = arange(amax(  squeeze(Input['dets']-1)))
        #self.wrong_dets = delete(self.dets,squeeze(Input['dets'])-1)
        #self.tvec = Input['tsel']
        #self.min_error = 0.001

        #self.calb_0 = [1]
        #self.dets_index = [self.dets]



    def load_mag_field(self):
	#  LOAD FIELD FROM NETWORK

	print("load_mag_field")

	try:
	    shot = self.shot if not self.artif_data else 65670 
	    d  = load(self.geometry_path+'/MagField_'+str(shot)+'.npz')
	    self.tsurf = d['tsurf']
	    self.magx = d['magx']
	    self.magy = d['magy']
	except Exception as msg:
            print(msg)
	    try:
                

		assert  self.local_network ,'cannot load mag from local cache'
		import pmds
		N_points = 100   # Number of grid points
		nsurf = 30.0
		psisur  = (arange(0.1, nsurf)/nsurf)**1.5  # from Ernesto (good spacing of flux surfaces)
		#psisur[0] = 0.02 #avoid singularity

		N_mag = 500;
		if len(self.tvec) > N_mag:
		    tsurf_new = linspace(self.tvec[0], self.tvec[-1], N_mag)
		else:
		    tsurf_new = self.tvec

		rsurf   =   zeros((N_points ,len(psisur), len(tsurf_new)))
		zsurf   =   zeros(shape(rsurf))
		pmds.mdsconnect('mdsplus.jet.efda.org')

		for j in range(len(tsurf_new)):
		    print("loading field", tsurf_new[j])
		    for i in range(len(psisur)):
                        #import code
                        #code.interact(local=dict(globals(), **locals()))
			try:
			    pmds.mdsvalue('flushinit(15,' +str(self.shot)+ ', ' +str(tsurf_new[j])+ ' ,0,"JETPPF","EFIT",_err)')
			except:
			    pass

			try:
			    pmds.mdsvalue('flupx(_err)')  # run flush
			except:
			    print("pmds mag fiels loading failed (maybe missing efit data ??)")
			    raise

			pmds.mdsvalue('flusur( ' +str(psisur[i])+ ',' +str(N_points) +',_xp,_yp,_br,_bz,1,_ier) ')  # retrieve data

			rsurf[:,i,j] =  pmds.mdsvalue('_xp')
			zsurf[:,i,j] =  pmds.mdsvalue('_yp')

		pmds.mdsdisconnect()
		magx_new = rsurf
		magy_new = zsurf
	    except Exception as msg:
                print(msg)
		print("load mag field failed <= loading mean field from file")
		#raise
		Input = loadmat(os.path.normpath(self.geometry_path+'/'+str(self.shot)+'/magsurf'+str(self.shot)+'.mat'), struct_as_record=True)
		tsurf_new = Input['tsurf']
		magx_new = Input['xx1'].transpose(1,2,0)
		magy_new = Input['yy1'].transpose(1,2,0)


	    if 'self.magx' in locals():           # ADD NEW DATA TO LOCAL CACHE !!!!

		print('ADD NEW DATA TO LOCAL CACHE !!!!')
		print("pred", shape(magx))
		#print shape(magx), shape(magx_new), shape(tsurf), shape(tsurf_new)
		self.magx = concatenate([self.magx, magx_new], axis=2)
		self.magy = concatenate([self.magy, magy_new], axis=2)
		self.tsurf = vstack([self.tsurf, tsurf_new])
		s = argsort(self.tsurf, axis=0)
		self.magx = self.magx[:,:,s]
		self.magx = self.magx[:,:,s]
		self.tsurf = self.tsurf[s]
		u, indices = unique(self.tsurf, return_index=True)
		self.magx = self.magx[:,:,indices]
		self.magy = self.magy[:,:,indices]
		self.tsurf = self.tsurf[indices]
		#print "pouzgzufzu", shape(self.magx)
		#exit()

	    else:
		print("tvorim tvec ...")
		self.tsurf =tsurf_new ; self.magx=magx_new; self.magy=magy_new

	    savez(self.geometry_path+'/MagField_'+str(self.shot)+'.npz', tsurf=self.tsurf, magx = self.magx, magy = self.magy)
	    #print "ulozeno"

	#self.nl = len(dets)


