#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from numpy import *
#from  matplotlib.pyplot import *
from scipy import sparse 
from numpy import linalg
import time
from scipy.interpolate import RectBivariateSpline, interp1d
from scipy.io import loadmat
import pickle
import os,os.path
import socket
from scipy.sparse import spdiags, eye
from copy import deepcopy, copy
import gc
import multiprocessing 
from prepare_data import * 
from shared_modules import debug
import  config
from shared_modules import make_postprocessing
from asymmetries import CalcAsymNew

    

try:
    from multiprocessing import Process, Pool
    import threading
    threading._DummyThread._Thread__stop = lambda x:40
except:
    pass


##AUTOR: Tomas Odstrcil  tomas.odstrcil@ipp.mpg.de

def tomography(inputs, tokamak, progress = None):
    """
    The heart of the tomography. Load setting from dictionary, call function to prepare data and perform reconstruction. Can use parallel version. This function is called by ``GUI.py``


    :param dict inputs:  dictionary with all setting generated by `main` and `loadSetting`
    :param QObject progress:    Object used to move with progressbar in GUI
    :var int regularizatalgorotihmsion: TODO remove
    :var int numSteps: Number of steps of recontruction, depends on tsteps and number of rapid blocks
    :var array G: Final reconstructed emissivity, shape (npix, tsteps)
    :var array sequence: Sequence of all tasks used for solving
    :var int numcores:  Number of threads for multiprocessing
    :var object pool:   Object for multiprocessing
    :var int numTasks:  Number of tasks for one calling of main(..)
    :var array CHI2:    Final chi2

    """

    print("="*10+" RECONSTRUCTION STARTED "+"="*10)

    try:
        from sksparse.cholmod import cholesky
    except:
        print("Can't find Cholesky decomposition for sparse matrices !!! \nFalling back to slow version (~6x slower)\n Install scikits.sparse package")

    tokamak.rgmin = inputs['rgmin']
    inputs['main_run'] = False #flag that this is not the main run yet
    # ================== PREPROCESSING,CALIBRATION ===========

    #special option for the reconstruction of the MHD modes
    if inputs['rem_back'] and inputs['presolver'] == 0 and  not inputs['substitute_back']:
        tokamak.allow_negative = True
    if tokamak.allow_negative:
        inputs['positive_constrain'] = False
        
  
    debug( "tokamak.default_calb", tokamak.default_calb)

    tokamak.prepare_emiss0( (inputs['tmax']+inputs['tmin'])/2.)
    
    data_all, error_all, tvec, dets, Tmat_all, normData = tokamak.prepare_data(inputs['tmin'], inputs['tmax'], 
                        inputs['data_smooth'], inputs['data_undersampling'], inputs['error_scale'], inputs['prune_dets'], False)

    tsteps = len(tvec)

    data  =  data_all[dets]
    error = error_all[dets]
    Tmat  =  Tmat_all[dets]

    if tokamak.default_calb != 'none':# and  tokamak.allow_self_calibration:
        if inputs['ratiosolver'] != 0:
            calb = pre_calibrate(tokamak, inputs, inputs['ratiosolver'], inputs['ifishmax'], 
                                    data, error, data_all,error_all, tvec, dets,Tmat, normData)
            cams = list(tokamak.detectors_dict.keys())
            calib_file = tokamak.geometry_path+'calibration'+os.sep+'%d.txt'%tokamak.shot

            f = open(calib_file, 'w')
            for k,i in zip(cams, calb): f.write('%s  %.3f\n'%( k,i ))
            f.close()
            debug(  'calibration factors:'+str(calb))
    
    tsteps = len(tvec)


    if inputs['solver'] == 0 or tsteps == 1 or 2*inputs['rapid_blocks'] >= tsteps:
        numSteps = tsteps
    else:
        numSteps = inputs['rapid_blocks']
    nx,ny = inputs['nx'], inputs['ny']
    
    num_iter = sum(r_[(inputs['ratiosolver']!=0), inputs['plot_all']*(1+inputs['plot_svd']),inputs['postprocessing'],inputs['rem_back'],\
            inputs['presolver']!= 0 and inputs['solver'] == 0,inputs['asymmetry'],
        inputs['sawtooths'],inputs['impurities'],inputs['plot_poloidal_spectrum'],2])
    if not progress is None:   #try to determine progress of reconstruction
        progress.setNumSteps(num_iter)
        progress.setNumSubSteps(1)
        progress.iterateStep()  #progress for data loading 


    
    ##==================== Phantom Generator==========================
    use_phantom_data = inputs['use_phantom_data']
    if isinstance(inputs['use_phantom_data'],  ndarray):
        use_phantom_data = array_str(use_phantom_data)

    if str(use_phantom_data) != "None" or inputs['solver'] == 7:

        if str(use_phantom_data) == "None": use_phantom_data = 'Asymmetric'

        from phantom_generator import phantom_generator
        nx_new = max(min(100, 4*nx),nx)
        ny_new = max(min(150, 4*ny),ny)
        nx_new = 200
        ny_new = 300
        
        scale = nanmean(data.ravel())*5#/nanmean(data_.ravel())
 
            
        data_,error_, tvec_new, G0 = phantom_generator(tokamak, tvec, nx_new=nx_new, 
                        ny_new=ny_new, profile = use_phantom_data, hollowness = 0.98,
                        edge = 0.7,scale=scale)

        print('generted')
  
        data = data_[dets,:]

        G_phantom = G0
        data_ = copy(data)
        
        if use_phantom_data == 'moving_circle':
            error[:] = 1
            
        elif inputs['solver'] != 7 and inputs['error_scale'] >  0.01:

            data[isfinite(error)] += error[isfinite(error)]*random.randn(sum(isfinite(error)))
            

        tokamak.emiss_0 = ones_like(G0)
    else:
        try:
            os.remove(tokamak.tmp_folder+'/Emiss0.npz')
        except:
            pass                

 

    
    #  ================== PRESOLVE CYCLE ========================



    TotTime = time.time()
    
    subst_ind = slice(None,None)
    if inputs['rem_back']:
        #substract averadge of the first 10% of points! 
        subst_ind = slice(0,int(ceil(tsteps*inputs['bcg_subst_fract'])))
        
    ifishmax = inputs['ifishmax']
    lam0 = 3*ones(numSteps)
    if inputs['solver']!= 7 and (inputs['solver'] == 0 or inputs['rem_back'] or inputs['rotation_tomography']):
        print("="*10+" PRESOLVE" +"="*10)
 
        G0,retro_fit, lam0 = presolve(tokamak, data[:,subst_ind], error[:,subst_ind], tvec[subst_ind],
                    Tmat, dets,  normData, inputs['regularization'], ifishmax, inputs, inputs['presolver'],inputs['rapid_blocks'])
        if not progress is None: progress.iterateStep()

        if inputs['presolver'] != 0:  ifishmax = 1 #will increase the stability
        debug( 'ifishmax in presolve is set to 1!!'  )
        debug('--------Prepare TIME - time  %g' % (time.time()-TotTime))
    elif inputs['solver'] != 7 or str(inputs['use_phantom_data']) == "None": #None solver 
        try:
            G0 = single(tokamak.emiss_0*mean(data)/mean(Tmat*tokamak.emiss_0)) #make a size of profile more similar 
        except:
            print('error single(tokamak.emiss_0*mean(data)/mean(Tmat*tokamak.emiss_0)) ')
            print(tokamak.emiss_0.shape, Tmat.shape)

    
    #===================== Remove background (if requested) ==================================  

    if inputs['rem_back'] and (inputs['presolver']!= 0 or (tokamak.allow_negative and tsteps > 1)):
        debug( 'rem back')
        if hasattr( tokamak,'mag_axis') and tokamak.name != 'DIIID' :

            
            #try to substract the position peturbation of the plasma core
            #if the plasma is moving, the movement of the bacground profile+LBO profile together
            #will create large false assymetries close to the core
    


            
            from  scipy.interpolate import RectBivariateSpline, UnivariateSpline
          
            G0_ = mean(G0[:,subst_ind],axis=1).reshape((ny,nx), order='F')

            interp2 = RectBivariateSpline(tokamak.ygrid, tokamak.xgrid,G0_)
            ind = (tokamak.mag_axis['tvec']>=tvec.min()-0.05)&(tokamak.mag_axis['tvec']<=tvec.max()+0.05)

            Zmag = tokamak.mag_axis['Zmag'][ind]
            Rmag = tokamak.mag_axis['Rmag'][ind]
            mag_tvec = copy(tokamak.mag_axis['tvec'][ind])

            R,z = interp(tvec,mag_tvec, Rmag),interp(tvec,mag_tvec, Zmag)
            R0,z0 = mean(R[subst_ind]), mean(z[subst_ind])
            dR,dz = R-R0, z-z0
            
            data0 = median(data[:,subst_ind],axis=1)
            data -= data0[:,None]

            G0_ = G0_.ravel(order='F')
            data_diff = zeros_like(data)
            peak_perturb = G0.max(0)/G0_.max()
            if G0.ndim==2:
               peak_perturb =  nanmean(data,0)/nanmean(data0)
            
            #remove movements of the plasma during LBO 
            for it in range(tsteps): 
                G0_i = interp2(tokamak.ygrid-dz[it],tokamak.xgrid-dR[it])
                data_diff[:,it] = Tmat*(G0_i.ravel(order='F')-G0_)*peak_perturb[it]
             

            data -= data_diff
            G0    = mean(G0[:,subst_ind],axis=1)[:,None]
            if not tokamak.allow_negative:
                data += Tmat*G0_[:,None]
            
            if str(use_phantom_data) != "None" or inputs['solver'] == 7:
                G_phantom-= mean(G_phantom[:,subst_ind],axis=1)[:,None]
            corrupted = isinf(error)
            error = tile(data[:,subst_ind].std(1)*inputs['error_scale'], (len(tvec),1)).T
            error[corrupted] = infty
            
        else:         

            data -= mean(data[:,subst_ind],axis=1)[:,None]
            error = tile(data[:,subst_ind].std(1)*inputs['error_scale'], (len(tvec),1)).T
            error[:] = error.mean()  #assume the same error for all LOS 
            corrupted = isinf(error)
            
            error[corrupted] = infty

            G0    = mean(G0[:,subst_ind],axis=1)[:,None]
            if not tokamak.allow_negative:
                data += Tmat*G0

            if str(use_phantom_data) != "None" or inputs['solver'] == 7:
                G_phantom-= mean(G_phantom[:,subst_ind],axis=1)[:,None]
                
                G0 = G_phantom
    else:
        #some channels can be due to a noise below zero
        pass

    ##==============MAIN CYCLE =======================

    if inputs['solver'] in r_[0:8]:
        print("="*10+" MAIN SOLVE "+"="*10)
        time_0 = time.time()

        import multiprocessing
        inputs['main_run'] = True #flag that this is the main run
        if inputs['solid_parallel'] or tokamak.npix < 1000:    #faster, but unbreakable, without progressbar, for nx*ny<1000 is parallel processing too slow
            numTasks = 1
            debug( 'Uses solid parallel solve => Cancel button will not work')
        else:
            numTasks = int(ceil(numSteps/(inputs['n_cpu']*4.)))
        numTasks = min(100,numTasks) #small speed improvement for very large number of blocks
        
        if config.DEBUG or config.no_multiprocessing:
            numTasks = numSteps

        ind = [slice(i*tsteps//numSteps,(i+1)*tsteps//numSteps) for i in range(numSteps)]
        


        if size(G0,1)== 1:
            G0 = [G0 for i in range(numSteps)]
        elif size(G0,1)== tsteps:
            G0 = [G0[:,ii] for ii in ind]

        
 
        if size(lam0) != numSteps:
            lam0 = ones(numSteps)*mean(lam0)
        
        #MDSplus object is not pickable
        if hasattr(tokamak, 'MDSconn' ):
            del tokamak.MDSconn
        if hasattr(tokamak, 'mds_server' ):
            del tokamak.mds_server
        mds_server = inputs.pop('mds_server', None)
             
            
        postprocessing = True
        sequence = array([(i,data[:,ii],error[:,ii],tvec[ii],Tmat,dets, normData[ii], 
                    G0[i],lam0[i],  numSteps, postprocessing,inputs, tokamak,
                    inputs['solver'], inputs['ifishmax'], time_0,postprocessing)  for i,ii in enumerate(ind)],copy=False,dtype=object)    
            
        sequence = array_split(sequence, numTasks)
   
        if config.DEBUG or numSteps==1 or config.no_multiprocessing:
            method =  map 
            pool = None
        else:
            pool = multiprocessing.Pool(min(inputs['n_cpu'],numSteps))
            method = pool.map  
            
                  
        results = []
        try:
            from PyQt5.QtCore import currentThread
        except:
            currentThread = None
     
        try:
            from tqdm import tqdm
            if not progress is None: progress.setNumSubSteps(len(sequence))
            for seq in tqdm(sequence,desc='Main cycle: '):
                if currentThread is not None:
                    if currentThread().isInterruptionRequested():
                        raise Exception('Externally interrupted!')
           
                out = method(main_cycle, seq)  #   pool.imap(main_cycle, seq)  #!!!open main cycle !!!
                gc.collect()
                results.extend(out)
                if not progress is None:  progress.iterateSubStep()
        except Exception as e:
            print('Error '+str( e))     #try to find the error
            out = list(map(main_cycle, seq)) 
            results.extend(out)
            
        finally:
            if pool is not None:
                pool.close()
                pool.join()
        
        sys.stdout.write("\r MAIN SOLVE DONE\n")
        if not progress is None: progress.iterateStep()
        inputs['mds_server'] = mds_server


        if numSteps == 1:
            output = results[0]

        else:
            keys = list(results[0].keys())
            output = {}
            keys.remove('bnd')

            for key in keys:   
                output[key] = concatenate([d[key] for d in results],0)
                for d in results: del d[key]
            output['bnd'] = []
            for d in results: output['bnd'].extend(d['bnd'])

        debug('--------Prepare TIME 3 - time  %g' % (time.time()-TotTime))

    else:
        raise Exception('solver %d not implemented yet'%inputs['solver'])
    
    
    TotTime = time.time()-TotTime
    print('Reconstruction time %.1fs' % TotTime )
    #=======================================
    print('DONE ..... ')
    gc.collect()        # free unused RAM memory !!!


    mean_chi2 = exp(nansum(log(output['chi2']))/tsteps)
    mean_lam = median(output['lam0'])
    
    print('Statistics: \n Mean time: %.3gs \n Mean Chi2: %.3g \n Mean Lam: %.3g \n Time slices: %d\n Blocks: %d'\
            %(TotTime/tsteps,mean_chi2,mean_lam,tsteps ,numSteps ))
    if 'mag_rho' in output:
        output['mag_rho'] =  output['mag_rho'][0]
    output['data'] = data.T
    output['tvec'] = tvec
    output['error']= error.T
    output['dets'] = dets
    output['data_all']  = data_all.T
    output['error_all'] = error_all.T
    output['Tmat_all']  = Tmat_all


    if str(use_phantom_data) != "None":

        
        Emiss0 = load(tokamak.tmp_folder+'/Emiss0.npz')
        
        xgridc = Emiss0['xgridc']
        ygridc = Emiss0['ygridc']
        ygridc_out = Emiss0['Zc']
        xgridc_out = Emiss0['Rc']
        G0 = Emiss0['G0']#.flatten(order='F')
        from scipy.interpolate import RectBivariateSpline

        emissivity_out = zeros((tsteps,len(xgridc), len(ygridc) ))
        for it in range(tsteps):
            Rint = RectBivariateSpline(xgridc_out,ygridc_out,output['g'][it].reshape(nx, ny,order='C'),kx=1,ky=1)                        
            emissivity_out[it] = Rint(xgridc,ygridc)
                
  
        dif  = linalg.norm(output['g'].T-G_phantom)/linalg.norm(G_phantom)
        dif2 = linalg.norm(emissivity_out-G0.T)/linalg.norm(G0)
        
     
        print('Comparism with phantom %.3f%%  %.3f%%'%(dif*100,dif2*100 ))
        dif3 = linalg.norm(output['retro']-data_.T)/linalg.norm(data_)
        dif4 = linalg.norm(output['data'] -data_.T)/linalg.norm(data_)

        f=open('out','a')
        f.write('%d %d %.5f %.5f %.5f %.5f\n'%(tokamak.nx,tokamak.ny, dif ,dif2, dif3,dif4  ))
        

    return inputs, tokamak, progress,output



def pre_calibrate(tokamak, inputs, solver,ifishmax, data, error, 
                        data_all,error_all, tvec, dets, Tmat,normData):
    
    """
    estimate relative calibration of the detectors by the method 

    :param class tok: ``Class`` of  tokamak with all important data.
    :param dict inputs:  dictionary with all setting generated by `main` and `loadSetting`
    :param int solver:  Type of presolving method
                       0 - Use default values
                       1 - Rapid
                       2 - SVD
                       ....
    :param int solver:     
    
    """

    tmin = inputs['tmin']
    tmax = inputs['tmax']
    data_smooth= inputs['data_smooth'] * 4  # slower changes !!!   
    data_undersampling= inputs['data_undersampling']
    error_scale= inputs['error_scale']
    prune_dets= inputs['prune_dets']
  
    Ndets = tokamak.Ndets
   
    
    regularization = 1 if inputs['regularization']!=0 else 0


    if Ndets == 1:
        return [1,]


    G0,retro_fit,lam = presolve(tokamak, data, error, tvec, Tmat, dets,  normData, 
                      regularization, ifishmax,inputs,solver,inputs['rapid_blocks'])
 
    new_calib = ones(Ndets)  #multiplicative modification of the orig calibration 

    det_ind = where([any(in1d(ind, dets)) for ind in tokamak.dets_index])[0]

        

    R = retro_fit.T
    D = data
    E = (1/error).mean(1)
    Eind = E == 0
    E[Eind ] = infty
    E[~Eind]= 1/E[~Eind]

    calib = ones(len(det_ind))
    for i,ind in enumerate(det_ind):
        ind = tokamak.dets_index[ind]
        ind = in1d(dets,ind)
        if any(ind):
            calib[i] = mean(D[ind]*R[ind]/E[ind,None])/mean(D[ind]**2/E[ind,None])+1e-5

    new_calib[det_ind] = calib
    new_calib[isnan(new_calib)] = 1
    
    if len(new_calib) == len(tokamak.get_calb()):
        new_calib *= tokamak.get_calb()
 

    mean_corr = exp(mean([log(abs(new_calib[i])) for i in det_ind]))
    new_calib /= mean_corr
    

    for i,ii in enumerate(det_ind):
        ind = tokamak.dets_index[ii]
        dind = in1d(dets,ind)
        data[dind]*= calib[i]
        error[dind]*= calib[i]
        data_all[ind]*= calib[i]
        error_all[ind]*= calib[i]

    return new_calib



def presolve(tokamak, data, error, tvec, Tmat, dets,  normData, 
             regularization,ifishmax, inputs, solver,rapid_blocks,G0=None):
    """
    Function containing all usable presolve methods : Rapid, small  resolution, SVD

    :param bool solver:  Type of presolving method
                       0 - no presolve
                       1,2,3,4,5,6 -  MFI rapid, SVD,SVD2, QR, GEV,GSVD
                       7 - SMALL RESOLUTION PRESOLVE

    :param array data:      Measured data
    :param array error:     Expected errors of data
    :param array tvec:      Time vector
    :param spmatrix Tmat:      Geometry matrix
    :param array dets:      Array of allowed detectors
    :param array normData:  Maximum of data in each timeslice
    :param dict inputs:  dictionary with all setting generated by `main` and `loadSetting`
    :param int numSteps:  Number of steps of recontruction, depends on tsteps and number of rapid blocks

    :var double lam0: Initial value of parameter lambda in min. fisher
    """
    from multiprocessing import Pool
    inputs = inputs.copy()
    inputs['regularization'] = regularization
    inputs['postprocessing'] = False
    solid_parallel = inputs['solid_parallel']
    nx = tokamak.nx
    ny = tokamak.ny
    tsteps = len(tvec)

    if solver == 0 or tsteps == 1 or  2*rapid_blocks > tsteps:
        numSteps = tsteps
    else:
        numSteps = rapid_blocks

    lam0 = 3 *ones(numSteps)
    retro = None
    postprocessing=False
    
    #MDSplus object is not pickable
    if hasattr(tokamak, 'MDSconn' ):
        del tokamak.MDSconn
        
    if solver in r_[1:7]:  # MFI rapid, SVD,SVD2, QR, GEV,GSVD   #(numSteps > 10 and solver == 1) or solver == 2:
        if G0 == None:
            G0 = single(tokamak.emiss_0)
        time_0 = time.time()
        nt = size(tvec)
        
        ind = [slice(i*nt//numSteps,(i+1)*nt//numSteps) for i in range(numSteps)]
        sequence = array([(i,data[:,ii],error[:,ii],tvec[ii],Tmat,dets, normData[ii],G0.mean(-1),lam0[i], 
                         numSteps, False,inputs, tokamak,
                        solver, ifishmax, time_0, postprocessing)  for i,ii in enumerate(ind)],copy=False,dtype=object) 



        out = []
        if  config.DEBUG:
            out = list(map(main_cycle, sequence))
        else:
            numTasks = 1 if solid_parallel else int(ceil(numSteps/inputs['n_cpu']))

            sequence = array_split(sequence, numTasks)

            from tqdm import tqdm
            pool = Pool(min(inputs['n_cpu'],numSteps) )

            for seq in tqdm(sequence,desc='Presolver: '):

                result = pool.map(main_cycle, seq)
                gc.collect()
                out.extend(result)
                
            pool.close()
            pool.join()
    
    
        if numSteps == 1:
            G0 =  out[0]['g'].T
        else:
            G0 = concatenate([d['g'] for d in out],0).T
        lam0 = concatenate([d['lam0'] for d in out],0)
        retro = concatenate([d['retro'] for d in out],0)

        gc.collect()
        sys.stdout.write("\r PRESOLVE DONE\n")
        


    if solver == 0:
        print("NO PRESOLVE")
        G0 = tokamak.emiss_0
        G0*= mean(data)/mean(Tmat*G0) #make a size of profile more similar 
        retro_fit = Tmat*G0
        
    if solver == 7:
        print("SMALL RESOLUTION PRESOLVE")

        
        tokamak_tmp = deepcopy(tokamak)
        tokamak_tmp.nx = nx//2
        tokamak_tmp.ny = ny//2
        
  

        Tmat = array(Tmat.todense())
        Tmat = sparse.csr_matrix(Tmat.reshape(-1,2,ny//2,2,nx//2,order='F').sum(3).sum(1).reshape(size(Tmat,0),-1,order='F'))

        from scipy.ndimage import zoom
        G0 = zoom(tokamak.emiss_0.reshape(ny,nx,order='F'),0.5).reshape(-1,1,order='F')


        
        tokamak_tmp.normTmat = tokamak.normTmat
        tokamak_tmp.npix = tokamak.npix//4
        if tokamak.transform_index != 0:
            Warning('Presolving with nontrivial transformation is not yet supported')
            
        tokamak_tmp.Transform = ones(1)  
        tokamak_tmp.transform_index = 0  
     
     
        lam0 = 3*ones(numSteps)
        time_0 = time.time()
        
        nt = size(tvec)
        ind = [slice(i*nt//numSteps,(i+1)*nt//numSteps) for i in range(numSteps)]
        
        sequence = [(i,data[:,ii],error[:,ii],tvec[ii],Tmat,dets, normData[ii], G0,lam0[i],0, numSteps,
                     False,inputs, tokamak_tmp,0, ifishmax, time_0,postprocessing)  for i,ii in enumerate(ind)]
    
   
        try:
            pool = Pool(min(inputs['n_cpu'],numSteps))
            out =  pool.map(main_cycle, sequence)
        except:
            print('multiprocessing in preprocessing has failured, turning it off')
            out =  list(map(main_cycle, sequence))
        finally:
            pool.close()
            pool.join() 

        sys.stdout.write("\r DONE\n")

            
        lam0 = concatenate([d['lam0'] for d in out],0)
        G0 = concatenate([d['g'] for d in out],0).T
        retro = concatenate([d['retro'] for d in out],0)

        G0 = G_interpolate(G0, nx/2, ny/2,  nx, ny)
        sys.stdout.write("\r\n")

    gc.collect()
    


    return G0,retro, lam0





def main_cycle(input):
    """
        Cycle iterating over the numSteps, (rapid_blocks), and reconstruct blocks/timeslices.

        :param int step: Number of current step from  range numSteps
        :param array data:      Measured data
        :param array error:     Expected errors of data
        :param array tvec:      Time vector
        :param array dets:      Array of allowed detectors
        :param array Tmat:      Geometry matrix
        :param int normTmat:    Norm of geometry matrix
        :param array normData:  Maximum of data in each timeslice
        :param array G0:        Initial guess of emissivity
        :param int numSteps:  Number of steps of recontruction, depends on tsteps and number of rapid blocks
        :param dict inputs:  dictionary with all setting generated by `main` and `loadSetting`
        :param double time_0:   Initial time of recontruction, it is used to guess remaining time
    """
    debug('Starter main cycle')
    if os.name != 'nt':
        os.nice(3)

    

    step,data,error,tvec, Tmat,dets,  normData, G0,lam0,  numSteps, \
            _postprocessing, inputs, tokamak, solver, _ifishmax,time_0,postprocessing = input


    boundary = inputs['boundary']
    regularization = inputs['regularization']
    
    reg_params = {}
    reg_params['danis'] =  inputs['danis']
    reg_params['pfx_weight'] = inputs.get('pfx_weight', 4)
    reg_params['sol_weight'] = inputs.get('sol_weight', 2)
    data = data.T
    error = error.T

    
    if G0.ndim >= 2 and solver!= 7:
        G0 = mean(G0,-1)
        
     
 
    lam0 = median(lam0)
    output = {}

    if solver in [0,1]:  #direct inversion
        from minfisher import minfisher
        tvec, g,retro, chi2, lam0,bnd,SigmaGsample = minfisher(data, error,tvec,Tmat,dets,
                        normData, G0,lam0, tokamak, reg_params, boundary,
                        regularization,_ifishmax, _postprocessing)
    elif solver in [2,3,4,5,6]: #SVD,SVD2, QR, GEV,GSVD solver
        from linear_methods import linear_methods
        tvec, g , chi2, lam0,retro,bnd,SigmaGsample  =  linear_methods(tokamak,inputs,
                data, error, tvec, Tmat,dets, normData, G0,  reg_params, boundary,regularization,
                solver,  True,_ifishmax, _postprocessing)
    elif solver in [7,]: # no solver - return phantom. Useful for testing of the postprocessing methods
        if (G0.shape[-1] == 1 and  len(tvec)!= 1) or G0.ndim == 1:
            g = tile(G0, (1,len(tvec)))
        else: g = G0
        retro = (Tmat*G0).T
        chi2 = linalg.norm((retro-data)/error,axis=1)**2/sum(isfinite(error),1)
        bnd = tokamak.get_boundary(N=100,time=mean(tvec))
        bnd = [bnd]*len(tvec)
        SigmaGsample = G0

    else:
        raise Exception( "\n\n Missing solver %d"%solver)

    if isscalar(chi2): chi2 = ones_like(tvec)*chi2
    if isscalar(lam0): lam0 = ones_like(tvec)*lam0

    if tokamak.transform_index!=0 and solver!= 7:
        g = tokamak.Transform*g
        if SigmaGsample is not None:
            SigmaGsample = tokamak.Transform*SigmaGsample

    output['g'] = g.T
    if SigmaGsample is not None:
        output['g_samples'] = SigmaGsample.T[None]
        
    output['retro'] = retro
    output['lam0'] = lam0
    output['bnd']  = bnd
    output['chi2'] = chi2
    output['chi2_real'] = linalg.norm((retro-data)/error,axis=1)**2/sum(isfinite(error),1)


    #place for the varius post-processing evaluated in multithreading
    if inputs['postprocessing'] and postprocessing:
        output = make_postprocessing(output,tokamak,g,data, error,retro, tvec, 
                                      dets,SigmaGsample,inputs)


    if inputs['asymmetry'] and postprocessing:
        output = CalcAsymNew(output,tokamak,g,data, error,retro, tvec, dets,
                             SigmaGsample,inputs)

    return output




def G_interpolate(g0, nx, ny,  nx_new, ny_new):
    """
    Interpolate solution on higher resolution
    Used for solving profile of function and for presolve in lower resolution
    """
    

    from scipy.ndimage import zoom
    tsteps = g0.shape[-1]
    g0 = reshape(g0, (ny,nx, tsteps), order='F')


    g_new = empty((ny_new, nx_new, tsteps))
    for i  in range(tsteps):
        g_new[:,:,i] = zoom(squeeze( g0[:,:,i]),2)

        
    g_new = g_new.reshape((nx_new*ny_new, tsteps), order='F')

    return g_new







