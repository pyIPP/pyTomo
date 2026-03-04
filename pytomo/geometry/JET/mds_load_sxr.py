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

    k = -1

    for j in range(len(dets)):
        dets_tmp = dets[j]
        #print dets
        #print len(dets)
        #print len(dets_tmp)
        gap = dets_tmp['gap']
        name = dets_tmp['name']
        N = dets_tmp['N']
        time = dets_tmp['time']
        dtb = dets_tmp['dtb']

        if N == 0:
            name2 = [ '%(name)s%(t)s' % {'name':name , 't':time} ]
        else:
	    if gap==3: # H detector starts with channel 02
		name2 = [ '%(dtb)s/%(name)s%(N)02d%(t)s' % {'dtb':dtb,'name':name , 'N':i, 't':time} for i in arange(N)+2]
            elif gap==2:
                name2 = [ '%(dtb)s/%(name)s%(N)03d%(t)s' % {'dtb':dtb, 'name':name , 'N':i, 't':time} for i in arange(N)+1]
            elif gap==1:
                name2 = [ '%(dtb)s/%(name)s%(N)02d%(t)s' % {'dtb':dtb,'name':name , 'N':i, 't':time} for i in arange(N)+1]
            elif gap == 0:
                name2 = [ '%(dtb)s/%(name)s%(N)d%(t)s' % {'dtb':dtb,'name':name , 'N':i, 't':time} for i in arange(N)+1]

        #print name
        #return

        #     another option how to make it



        from tqdm import trange
        #for seq in trange(sequence,desc='Main cycle: '):
        for i in trange(N,desc='Loading SXR dat'+str(pulse)+': '):
            #print name2[i]
            k +=1

            try:
                if GETDAT:
                    d, tvec_tmp, nwords, title, units, ierr = getdat(name2[i], pulse)
                else:
                    name_tmp = '_sig=jet("'+name2[i]+'",'+str(pulse)+')'
                    #print name_tmp
                    d=pmds.mdsvalue(name_tmp)

                    #print d

                    if type(d) is not str and len(d)!=0:
                        tvec_tmp = array(pmds.mdsvalue('dim_of(_sig,0)'))
                    else:
			print(name_tmp)
                        print(d)
                        continue


                if j == 0 and i == 0:    # first step !!!!
                    tvec = tvec_tmp
                    ts = len(tvec)
                    tvec = tvec[int(ts/40):-int(ts/40)]    #fixes the end and beginning in interpolation
                    datalen = len(tvec)
                    data = ones((datalen, n_dets))*NaN
                    #exit()

                if type(d) is not str and len(d)!=0:
                    from scipy.interpolate import interpolate
                    d=interpolate.interp1d(tvec_tmp,d, bounds_error=False, fill_value=0)(tvec)
                    data[:,k]  = d
                #else:
                    #print "loading error detector:", i, d
            except:
                raise
                print('error channel', i)


            #sys.stdout.write("\r %2.1f %%" %(i/double(N)*100))
            #sys.stdout.flush()
        #sys.stdout.write("\r")
    pmds.mdsdisconnect()
    return tvec, data




def mds_load_sxr(pulse, data_path, SOURCE,time_interval,  UPLOAD=False):

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

    print("==========  mds_load_sxr ============")

    #print name2
    dets = list()
    if SOURCE == 'slow':
        dets.append( {'name':'DB/J3-SXR<V', 'dtb':'jpf', 'N':35, 'gap':0, 'time':'/1'})     #gap is number of 0 before integer
        dets.append( {'name':'DB/J5-SXR<S4:','dtb':'jpf', 'N':17, 'gap':2, 'time':'/1'})
        dets.append( {'name':'DB/J3-SXR<T', 'dtb':'jpf','N':35, 'gap':0, 'time':'/1'})
    elif SOURCE == 'fast':
        dets.append( {'name':'DB/J4-CATS<V', 'dtb':'jpf','N':35, 'gap':0, 'time':'/'+str(time_interval)})     #gap is number of 0 before integer
        dets.append( {'name':'DD/J5-RTVS<S4:', 'dtb':'jpf','N':17, 'gap':2, 'time':'/'+str(time_interval)})    #time cant be 'all' => Memory Error
        dets.append( {'name':'DB/J4-CATS<T', 'dtb':'jpf','N':35, 'gap':0, 'time':'/'+str(time_interval)})
    elif SOURCE in ['new_slow', 'new_fast']:
        if pulse > 83900 :
            dets.append( {'name':'SXR/V', 'dtb':'ppfgud','N':35, 'gap':1, 'time':'?JETPPF'})     #gap is number of 0 before integer
            dets.append( {'name':'SXR/H', 'dtb':'ppfgud','N':17, 'gap':3, 'time':'?JETPPF'})     #here gap is used to distinguish H
            dets.append( {'name':'SXR/T', 'dtb':'ppfgud','N':35, 'gap':1, 'time':'?JETPPF'})     #time cant be 'all' => Memory Error
        else:
            dets.append( {'name':'SXR/V', 'dtb':'ppfgud','N':35, 'gap':1, 'time':'?chain1'})     #gap is number of 0 before integer
            dets.append( {'name':'SXR/H', 'dtb':'ppfgud','N':17, 'gap':3, 'time':'?chain1'})     #here gap is used to distinguish H
            dets.append( {'name':'SXR/T', 'dtb':'ppfgud','N':35, 'gap':1, 'time':'?chain1'})     #time cant be 'all' => Memory Error
                
                
    else:
        print('Source undefined')
        exit()

    tvec, data = loader(dets ,pulse )

    #from signaltools import decimate


    if SOURCE == 'new_slow':
        #downsample the data

        from shared_modules import MovingAveradgeFast
        data = MovingAveradgeFast(data, 10, axis=0)
        #introduce to tvec the same shift
        tvec = MovingAveradgeFast(double(tvec), 10) 
        
       
        data = data[::100, :]
        tvec = tvec[::100]



    if SOURCE in ['new_slow', 'new_fast']:
        #data_new = zeros((size(data,0), 87))+NaN
        #data_new[:,arange(87) != 52] = data
        #data_new[:,arange(87) != 35] = data
        #data = insert(data, 35, ones(len(data))*NaN, 1)
        #data = delete(data, 36, 1)
        data[:,35] = NaN
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



