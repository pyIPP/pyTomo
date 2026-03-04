#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pmds
from numpy import *
import matplotlib
from scipy import sparse
import sys
from scipy.io import loadmat, savemat
import os

GETDAT=False  #PMDS vesion or GETDAT version


#vyresit multithreading

def loader(dets,pulse ):
    if GETDAT:
        from getdat import getdat

    print("Connecting")
    pmds.mdsconnect('mdsplus.jet.efda.org')
    print('loading data...')
    n_dets = 0
    for j in range(len(dets)):
        n_dets +=  dets[j]['N']

    k = 0

    for j in range(len(dets)):
        print(j)
        dets_tmp = dets[j]
        #print dets
        #print len(dets)
        #print len(dets_tmp)
        gap = dets_tmp['gap']
        name = dets_tmp['name']
        N = dets_tmp['N']
        dtb = dets_tmp['dtb']
        print(N)
        python_bolo_data='bolo_data'

        if N == 0:
            name2 = [ '%(dtb)s/%(name)s' % {'dtb':dtb, 'name':name }]
        else:
            if gap==2:
                name2 = [ '%(dtb)s/%(name)s%(N)03d' % {'dtb':dtb, 'name':name , 'N':i} for i in arange(N)+1]
            elif gap==1:
                name2 = [ '%(dtb)s/%(name)s%(N)02d' % {'dtb':dtb,'name':name , 'N':i} for i in arange(N)+1]
            elif gap == 0:
                name2 = [ '%(dtb)s/%(name)s%(N)d' % {'dtb':dtb,'name':name , 'N':i} for i in arange(N)+1]

        #print name
        #return

        #     another option how to make it

        for i in arange(len(name2)):
            print(name2[i])

	    #print "exit"
	    #exit()
	    

            try:
                if GETDAT:
                    d, tvec_tmp, nwords, title, units, ierr = getdat(name2[i], pulse)
                else:
                    name_tmp = '_sig=jet("'+name2[i]+'",'+str(pulse)+')'
                    d= array( pmds.mdsvalue(name_tmp) ) 

                    if type(d) is not str and len(d)!=0:
                        tvec_tmp = array(pmds.mdsvalue('dim_of(_sig,0)'))
			try:
			    tvec_tmp = array(pmds.mdsvalue('dim_of(_sig,1)'))
			except:
			    pass
			
                    else:
                        print(d)

                if j == 0 and i == 0:    # first step !!!!
                    tvec = tvec_tmp
                    ts = len(tvec)
                    tvec = tvec[int(ts/40):-int(ts/40)]    #fixes the end and beginning in interpolation
                    datalen = len(tvec)
                    data = ones((datalen, 0))
                    #exit()

                if type(d) is not str and len(d)!=0:
                    from scipy.interpolate import interpolate
		    
                    d=interpolate.interp1d(tvec_tmp , d.T, bounds_error=False, fill_value=0)(tvec)
                    #data[:,k]  = d
		    		    
		    data = hstack([ data, d.T ])
		    
		    
                #else:
                    #print "loading error detector:", i, d
            except:
                raise
                print('error channel', i)

            k +=1

            sys.stdout.write("\r %2.1f %%" %(i/double(N)*100))
            sys.stdout.flush()
        sys.stdout.write("\r")
	
    print("diconnecting")
    pmds.mdsdisconnect()
    print("saving")
    savez(python_bolo_data , data = single(data), tvec = tvec, dets = dets)
    return tvec, data




def mds_load_bolo(pulse, data_path, SOURCE,time_interval,  UPLOAD=False):

    """
    Beta version allowing to read SXR data for JET tokamak class
    """

    #TODO zařidat načítání neutronu + mag pole
    #zkusit savez na ukladani
    #Ndets = 87
    #UPLOAD=False

    #str(time_interval) = '4'

    #try:
        #print 'loading'
        #source = load(data_path+'/Data_'+str(pulse)+'_'+str(time_interval)+'_'+SOURCE+'.npz')
        #data = source['arr_0']
        #tvec = source['arr_1']
        #print 'loaded'
        ##if UPLOAD:
            ##try:
                ##datasend(data_path, 'Data_'+str(pulse)+'_'+str(time_interval)+'_'+SOURCE+'.npz')
            ##except:
                ##raise
    #except:
        #raise


    #print name2
    dets = list()
    if SOURCE == 'KB':
        dets.append( {'name':'DB/B5HR-UBOL<RAW:', 'dtb':'jpf', 'N':24, 'gap':2})     #gap is number of 0 before integer
        dets.append( {'name':'DB/B5VR-UBOL<RAW:','dtb':'jpf', 'N':24, 'gap':2})
        dets.append( {'name':'DB/B3-UBOL<RAW:', 'dtb':'jpf','N':24, 'gap':2})
        dets.append( {'name':'DB/B4-UBOL<RAW:', 'dtb':'jpf','N':24, 'gap':2})
    elif SOURCE == "KB5":
	#BOLO/KB5V 
	#BOLO/KB5H
	dets.append( {'name':'BOLO/KB5V', 'dtb':'ppf','N':0, 'gap':0})
	dets.append( {'name':'BOLO/KB5H', 'dtb':'ppf','N':0, 'gap':0})

	pass
    else:
	print(SOURCE, SOURCE =="KB")
        print('Source undefined')
        exit()

    tvec, data = loader(dets ,pulse )

    if SOURCE == 'new':
        #data_new = zeros((size(data,0), 87))+NaN
        #data_new[:,arange(87) != 52] = data
        #data_new[:,arange(87) != 35] = data
        data = insert(data, 35, ones(len(data))*NaN, 1)
        #data = data_new


    #print tvec
    #print data
    #exit()

    #print 'saving'
    #datafile = data_path+'/Data_'+str(pulse)+'_'+str(time_interval)+'_'+SOURCE
    #savez(data_path+'/Data_'+str(pulse)+'_'+('FAST' if FAST else 'SLOW')+'.npz', data, tvec)   #caching
    #savez(datafile, data, tvec)
    #savemat(datafile, {'data':data, 'tvec':tvec})
    #print 'saved'
    #if UPLOAD:
        #import multiprocessing
        ##p = multiprocessing.Process(target=datasend, args=(data_path, 'Data_'+str(pulse)+'_'+SOURCE+'.npz'))
        ##p.start()
        ##p2 = multiprocessing.Process(target=datasend, args=(data_path, 'Data_'+str(pulse)+'_'+SOURCE+'.mat'))
        ##p2.start()

        #datasend(data_path, 'Data_'+str(pulse)+'_'+str(time_interval)+'_'+SOURCE+'.npz')
        ##datasend(data_path, 'Data_'+str(pulse)+'_'+SOURCE+'.mat')

        #try:
            #os.remove(datafile+'.npz')
            #os.remove(datafile+'.mat')
        #except:
            #raise


    return [tvec,data]


if __name__ == "__main__":

    for PULSE in [68546]:
        for SOURCE in [ 'GAMMA', 'DENSITY', 'HA', 'SLOW_SXR', 'FAST_SXR' ]:
            mds_load_sxr(PULSE, './data', SOURCE, True)
    print('done')


