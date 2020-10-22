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
from scipy.interpolate import interpolate
from scipy.signal import medfilt, convolve
import urllib.request, urllib.parse, urllib.error
from copy import deepcopy
from tokamak import Tokamak
from scipy.stats.mstats import mquantiles


class Golem(Tokamak):
    def __init__(self, input_diagn, input_parameters,  only_prepare = False, Color = 0  ):

        """ Initialize new tokamak and set default values.

            :param dict input_parameters:   List of selected paramerers
            :param str input_diagn:  Name of main diagnostics
            :param array MagField_shift: artificial shift of magnetic field

        """
        # BASIC PARAMETERS

        artif_data = False

        name = 'Golem'
        coord = [0.28,  0.52, -0.12, 0.12]
        path = 'geometry_golem'


	sigma = 0.0
        min_error = 0.08
        norm = 1.0
        t_name = 's'

	path = 'geometry_golem/'+input_diagn
	
	Tokamak.__init__(self, input_parameters, input_diagn, coord,  name, norm, sigma, min_error, path,t_name)


        if input_diagn == "SXR":
            self.Bolometry()
            self.index = 1
        if input_diagn == "camera":
            self.Camera(Color)

            self.index = Color + 7
        self.boundary_coord = loadtxt(self.geometry_path+'/border.txt')

	self.tsteps = len(self.tvec)

	if not only_prepare:
	    self.prepare_tokamak()
	#self.load_mag_equilibrium()
	self.mag_equilibrium([], preferCache = False, dryRun = True)

#mag_equilibrium(

    def Bolometry(self):
        """
        Priprava dat pro tokamak Golem !!!!!
        """
        #self.geometry_path =  os.path.normpath(self.input_path+'geometry_golem')


        try:
            Nidatap = loadtxt(self.input_path+'/data_golem/'+str(self.shot))
            Nidatap = abs(Nidatap)  #TODO !!!!!! opravit !!!!!!!!
        except:
            print("Can't load data from local file")
            try:
                address = "http://golem.fjfi.cvut.cz/operation/shots/"+str(self.shot)+"/das/Radiation/ON/1210Bolometry/data"
                data = urllib.request.urlopen(address)
            except:
                print("Can't load the pulse "+str(self.shot)+" from internet")
                raise

            data = data.read()
            if 'Not Found' in data:
                print("Can't load "+address)
                raise IOError("Can't load", None, address)
            else:
                output = open('tmp/data.txt','w')
                output.write(data)
                output.close()
                try:
                    Nidatap = loadtxt('tmp/data.txt')   #convert text to array of numbers
                except:
                    print("Can't load data from file (from internet)")
                    raise


        first = 0.006
        smooth = 2
        self.tvec = Nidatap[:,0]
        self.dets = arange(size(Nidatap,1)-1)
        last = self.tvec[-1]
        for i in range(len(self.tvec)-10) :        #detekce konce pulzu, (signál je < 0.01)
            if amax(mean((Nidatap[i+arange(10),1:]),axis=1)) < 0.01 and self.tvec[i] > 0.01:
                last = self.tvec[i]
                break

        #pulse = arange(len(self.tvec))
        pulse = where((self.tvec > first) & (self.tvec < last))[0]
        self.data = Niself.datap[pulse,1:]
        self.tvec = self.tvec[pulse]
        for j in self.dets:
            self.data[:,j] = medfilt(self.data[:,j], 5)
        self.data = self.data.T
        self.data[self.data < 0 ] = 0    #remove negatives

	self.save_data( self.data)


    def Camera(self, color):
        # high speed visible camera

        print("loading data ")

        #self.geometry_path =  os.path.normpath(self.input_path+'geometry_golem_camera')

        self.wrong_dets_pref = []
        self.wrong_dets = []
        self.corrected_dets = []

        self.make_geometry()

        ##Nahradit chybejici detektory za NaN
        #====================
        time_shift = 0.0# 0.1 #ms
        self.calb = [1,1]
        #==================

        global WIDTH, NUM_CAM, gap_size, frame_len
        WIDTH = 336
        NUM_CAM=0
        # orez
        gap_size = 16  #TODO ! zkontrolovat !!!
        frame_len = 96

	tvec_all = list()
	data_all  = list()
	#calibration = [1.0, 1.0]#   [1.0,1.5]

	ds = DataSource(self.geometry_path)


	# LOAD PHOTODIODE  -- TOTAL EMISSIVITY, TIME SYNCHRONIZATION
	#try:
	    #self.pDiode = loadtxt(self.geometry_path+'/pDiode_'+str(self.shot)+'.txt')
	    #plas
	#except:
	    #try:
	base = "http://golem.fjfi.cvut.cz/cgi-bin/data/"+str(self.shot)

	plasma_start = float(ds.open(base+"/plasma_start").read())
	plasma_end = float(ds.open(base+"/plasma_end").read())
	
	self.pDiode = loadtxt(ds.open(base+"/photodiode"))
	self.pDiode = self.pDiode[(self.pDiode[:,0] > plasma_start -  1.5e-3) & (self.pDiode[:,0] < plasma_end + 2e-3), :]
	self.pDiode[:,1] /= mquantiles(self.pDiode[:,1], 0.95)

	#print ds.open(base+"/plasma_start")

	#print pDiode

	#self.pDiode = loadtxt(self.geometry_path+'/pDiode_'+str(self.shot)+'.txt')
	#plasma_start = 0.008
	#plasma_end = 0.021

	#for i in ['plasma_start', 'plasma_end']:
	    #address = "http://golem.fjfi.cvut.cz/cgi-bin/data/"+str(self.shot)+"/"+i
	    #exec(i + '= loadtxt(urllib.urlopen(address))')
	    

	##pDiode[:,1] = abs(pDiode[:,1])**(3.45)
	#pDiode[:,1] /= amax(pDiode[:,1])
	##pDiode[:,0] -= plasma_start
	##plot(pDiode[:,0], pDiode[:,1])
	##show()
	##exit()
	
	#self.pDiode = pDiode

	#savetxt(self.geometry_path+'/pDiode_'+str(self.shot)+'.txt', pDiode)
	    #except:
		##raise 
		#print "LOADING PHOTODIODE FROM WEB FAILED"


	self.dets_index = list()

	self.emiss_colors = list()
	
	for CAMERA in [1,2]:
	    print("loading camera "+ str(CAMERA))

	    try:
		data = imread(ds.open( "http://golem.fjfi.cvut.cz/operation/shots/"+str(self.shot)+"/diagnostics/Radiation/0211FastCamera.ON/"+str(CAMERA)+"/CorrectedRGB.png"))
		NUM_CAM +=1
		if CAMERA == 1:
		    self.dets_index.append(arange(WIDTH))
		else:
		    self.dets_index.append(arange(WIDTH)+amax(self.dets_index)+1)

	    except:
		#raise 
		assert not NUM_CAM == 0 and CAMERA == 2,  "Can't load the pulse "+str(self.shot) + " " + str(NUM_CAM)
		continue

	    print("================== CAMERA "+  str(CAMERA) + " loaded")

	    data = double(data)

	    data = data*255


	    #převod intenzity do lineární škály
	    a = 43000#+-600
	    b = 140#+-3
	    c = 51.0/1460
	    data = data/(a-b*data)/c   #normované od 0 do 1

	    # it should be now automatic ...
	    #data *= self.calb[CAMERA-1] # !!!!!!!!!!!!!!!!!!§

	    from shared_modules import nanmedian
	    self.emiss_colors.append( sum(data,0) )

	    self.emiss_colors[CAMERA-1], ind_wrong = self.fix_signal(self.emiss_colors[CAMERA-1], axis = 0)
	    
	    #  !!!!!!!!find time vector !!!!!!!!!! 
	    N = len(self.emiss_colors[CAMERA-1] )
	    tvec = linspace(plasma_start, plasma_start + N*(4.8/624) * 1e-3, N)

	    cam =  self.emiss_colors[CAMERA-1] / nanmax( self.emiss_colors[CAMERA-1] )
	    from scipy.optimize import fmin
	    from scipy.interpolate import interp1d
	    cam_tmp  = cam[:,0]*7 + cam[:,1]*0 + cam[:,2]*0.8
	    ind=~isnan(cam_tmp)

	    cam_tmp = cam_tmp[ind]
	    tvec_tmp = tvec[ind]
	    
	    #plot( tvec_tmp,  cam_tmp)
	    #show()

	    def find_tshift(tshift, return_result = False):
		D = interp1d(tvec_tmp+tshift, cam_tmp, bounds_error=False, fill_value=0 )(self.pDiode[:,0]) - self.pDiode[:,1] 
		if return_result: return D
		return nansum(D**2)
	    
	    tshift = fmin(find_tshift, 0)
	    print("tshift"+str(  tshift)) 
	    if abs(tshift) > 0.0015:
		tshift = 0  # neco selhalo 
		
	    tvec += tshift 
	    #cam_new = find_shift(shift, return_result = True)
	    
	    #plot(tvec+shift, cam )
	    #plot( tvec+shift, cam[:,0]*7 + cam[:,1]*0 + cam[:,2]*0.8 , 'k')

	    #plot(self.pDiode[:,0], self.pDiode[:,1])
	    #show()
	    #exit()
	    
	    if color == 3:  # sum all colors
		data = sum(data,2)
	    else:
		data = data[:,:,color]
		
	    #show()
	    #imshow(data)
	    #show()
	    
	    data = self.fix_signal(data, axis = 1, ind_wrong = ind_wrong)
	    
	    #plot(sum(data,0))
	    #show()
	    #imshow(data)
	    #show()
	    
	    


	    gaps_filling = True
	    

	    if gaps_filling:

	
		print("======== GAPS FILLING ===============")

		data[:,[0,-1]] = nanmin(data[data>0])

		ind_ok = sum(data,0) > 0
		data = interpolate.interp1d( array(where(ind_ok)[0]), data[:,ind_ok],  kind='linear')(arange(size(data,1)))


	    for i  in arange( size(data,1)):
		data[:,i] = medfilt(data[:,i],21)  #remove outliers in on timeslice  >>   reflections ...

	    data_all.append(data)
	    tvec_all.append(tvec)



	if NUM_CAM == 2:
	    tmin = max(amin(tvec_all[0]), amin(tvec_all[1]))
	    tmax = min(amax(tvec_all[0]), amax(tvec_all[1]))
	    for i in range(2):
		ind =  (tvec_all[i] >= tmin ) & (tvec_all[i] <= tmax )
		data_all[i] = data_all[i][:,ind]
		tvec_all[i] = tvec_all[i][:,ind]
		self.emiss_colors[i] = self.emiss_colors[i][ind,:]

	    data = zeros((2*WIDTH, max( len(tvec_all[0]), len(tvec_all[1])) ))
	    data[:WIDTH,:len(tvec_all[0])] = data_all[0]
	    data[WIDTH:,:len(tvec_all[1])] = data_all[1]
	    self.tvec = tvec_all[0]
	    
	else:
	    tmin = amin(tvec_all[0])
	    tmax = amax(tvec_all[0])
	    self.calb = [1]

	    self.tvec = arange(ceil(tmax-tmin)+1e-6)*(5.0/624) * 1e-3

	self.dets = arange(size(data, 0))
	self.data = data

	save('tmp/emiss_color', self.emiss_colors)
	

        self.data = self.data / amax(self.data)
        self.data += 0.0001


	self.save_data( self.data)

    
    def fix_signal(self,signal, axis = 0, ker = 35, ind_wrong = None):
	""" remove damages signal from cameras """
	
	if axis == 1:
	    signal = signal.T
	if ind_wrong is not None:
	    #subplot(1,2,1)
	    #imshow(copy(signal))

	    signal[ind_wrong,:] = nan
	    #subplot(1,2,2)

	    #imshow(signal)
	    #show()
	    signal = signal[where(~ind_wrong)[0][0]:where(~ind_wrong)[0][-1], :]
	    if axis == 1:
		return  signal.T
	    return signal
	
	ind_wrong = sum(signal, 1)  < 0.2*median(sum(signal, 1) )  # remove weak signal 
	#ind_wrong = convolve(ind_wrong, ones(5), mode="same") > 0
	
	#plot(sum(signal, 1))
	#show()
	
	s = std(signal, 0) 
	signal /= s
	d = sum(signal, 1) 
	

	
	
	signal_tmp = copy(signal)

	for i in range(size(signal, 1)):
	    signal_tmp[~ind_wrong,i] = medfilt( signal_tmp[~ind_wrong,i] , ker )
	sdiff = amax((signal_tmp - signal)**2,1)


	ind_diff = sdiff > mquantiles(sdiff, 0.90)
	ind_diff = convolve(ind_diff, ones(8), mode="same") > 0

	ind_wrong |= ind_diff
	
	#plot(ind_wrong)
	# close the image
	ind_wrong = convolve(~ind_wrong, ones(12), mode="same") < 0.5
	ind_wrong = convolve(ind_wrong, ones(14), mode="same") > 0.5
	#plot(1.1+ind_wrong)
	#show()

	signal[ind_wrong, :] = nan

	sl = slice(where(~ind_wrong)[0][0],where(~ind_wrong)[0][-1])

	signal = signal[sl, :]#  cutoff the dark beginning andend of the plasma

	#plot(signal)
	#show()
	
	signal *= s
	

	
	if axis == 1:
	    signal = signal.T
 
	return signal, ind_wrong

    def load_mag_equilibrium(self):
	print("standart mag field failed")
	
	#xxxx
	
	try:
	    data = loadtxt(self.geometry_path+'/plasma_parameters_'+str(self.shot))
	    boundary = loadtxt(self.geometry_path+'/border.txt')

	    mag_tvec=data[:,0]
            xmass =  data[:,1]
            ymass =  data[:,2]
            power =  data[:,3]

            tsurf =  (mag_tvec-mean(mag_tvec))*1.1+mean(mag_tvec)   # extent time range for interpolation

            minpower = 0.001*median(power)
            corr = (power > minpower) & ~isnan(xmass) & ~isnan(ymass)

            ymass_new= polyval(polyfit(mag_tvec[corr], medfilt(ymass[corr], 1+2*floor(sum(corr)/200)), 6), tsurf)
            xmass_new = polyval(polyfit(mag_tvec[corr], medfilt(xmass[corr],1+2*floor(sum(corr)/200)), 6), tsurf)

            nbord = len(boundary)
            tsteps = len(tsurf)
            nsurf = 60; #number of magnetic surfaces

            bx = boundary[:,0]
            by = boundary[:,1]

            #plot(xmass)
            #plot(xmass_new, '.')
            #plot(ymass)
            #plot(ymass_new, '.')

            #show()

            #exit()


	    ymass = ymass_new
	    xmass = xmass_new

            magx = zeros((nbord, nsurf, tsteps))
            magy = zeros((nbord, nsurf, tsteps))


            bx = reshape(bx,(-1,1))
            by = reshape(by,(-1,1))
	    lin=linspace(0.01,1,nsurf)*1.1;
            lin = reshape(lin, (1,-1))

            for i in arange(tsteps):  #mag. field is circular
                magx[:,:,i] =((bx-mean(bx))*lin+xmass[i])	#artificial magnetic field
                magy[:,:,i] =((by-mean(by))*lin+ymass[i])

            self.tsurf = tsurf
            self.magx = magx
            self.magy = magy

	except:
	    #print "failed "
	    #exit()
	    #raise
	    magnet_field =  loadtxt(os.path.normpath(self.geometry_path+'/magnet_field.txt'))
	    from prepare_data import convert_mag_field
	    magx, magy = convert_mag_equilibrium(magnet_field)

	    magx = reshape(magx, (size(magx,0), size(magx,1), 1))
	    magy = reshape(magy, (size(magy,0), size(magy,1), 1))
	    self.magx = magx.repeat(len(self.tvec), axis = 2)
	    self.magy = magy.repeat(len(self.tvec), axis = 2)
	    self.tsurf = self.tvec


    def make_geometry(self):


        n_chord = 336

        # bocni fotak

        x0_start = 0.4*ones((n_chord, 1))
        y0_start = reshape(linspace(0.09, -0.08, n_chord), (-1,1))

        y0_end_0 = -0.025*ones((n_chord, 1))
        y0_end_1 = 0.03*ones((n_chord, 1))
        x0_end   = 0.81*ones((n_chord, 1))


        # horni fotak
        x1_start = reshape(linspace(0.50, 0.32, n_chord), (-1,1))
        y1_start = 0*ones((n_chord,1))

        x1_end_0 = 0.39*ones((n_chord, 1))
        x1_end_1 = 0.44*ones((n_chord, 1))
        y1_end = 0.315*ones((n_chord, 1))


        x_0 = hstack([x0_end, x0_start])
        y_0 = hstack([y0_end_0, y0_end_1, y0_start])
        x_1 = hstack([x1_end_0, x1_end_1, x1_start])
        y_1 = hstack([y1_end, y1_start])

        savetxt(self.geometry_path+'/detector_0_x.txt',x_0, fmt='%5.4f')
        savetxt(self.geometry_path+'/detector_0_y.txt', y_0, fmt='%5.4f')
        savetxt(self.geometry_path+'/detector_1_x.txt',x_1 , fmt='%5.4f')
        savetxt(self.geometry_path+'/detector_1_y.txt', y_1, fmt='%5.4f')



    def data_correction(self, data):



		
        global WIDTH, NUM_CAM, gap_size, frame_len

        n_pic =  size(data,1)/(frame_len+gap_size)

        pic= zeros((336,96,3,n_pic))

        pic = double(pic)*255

        #převod intenzity do lineární škály
        a = 43000#+-600
        b = 140#+-3
        c = 51.0/1460

        pic = pic/(a-b*pic)/c   #normované od 0 do 1

        SumOfPics = sum(sum(pic,2),2)   #sloučení obrázků a barev

        mask = zeros(shape(SumOfPics))


        for i in range(336):
            mask[i,:] =  SumOfPics[i,:]/(mean(SumOfPics[i,:]+0.001))

        mask[mask<1] = 1  #vylepšit?

        maskRGB = zeros((336,96,3))
        for i in range(3):
            maskRGB[:,:,i] = mask

        #originál pro porovnání
        OrigRGB = zeros((336,gap_size/2,3))
        CorrectedRGB = zeros((336,gap_size/2,3))
        for i in range(size(pic,3)):
            OrigRGB = hstack([OrigRGB,pic[:,:,:,i],zeros((336,gap_size,3))])
            CorrectedRGB =hstack([CorrectedRGB,pic[:,:,:,i]/maskRGB,zeros((336,gap_size,3))])

        ## přenormování intenzity na nelineární škálu
        OrigRGB         = OrigRGB*a/(1/c+b*OrigRGB)
        CorrectedRGB    = CorrectedRGB*a/(1/c+b*CorrectedRGB)

        return OrigRGB

