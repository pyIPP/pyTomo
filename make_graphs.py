#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
import config
import time
import sys,os
import scipy as sp
from scipy.stats.mstats import mquantiles
from tqdm import trange
from shutil import copyfile
from scipy.interpolate import interp1d
from matplotlib.ticker import MaxNLocator,FormatStrFormatter,ScalarFormatter,NullFormatter
from matplotlib.transforms import Bbox
from matplotlib.figure import Figure
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvas as FigureCanvas
from matplotlib import colors,cbook,cm
from matplotlib.colors import Normalize

from shared_modules import fsvd,fast_svd,debug, get_rho_tan
from annulus import get_bd_mat, get_rho_field_mat
try:
    from multiprocessing import Process, Pool, cpu_count
    import threading
    threading._DummyThread._Thread__stop = lambda x:40
except:
    pass

import matplotlib


try:
    from scipy.stats import nanmedian,nanstd
except:
    pass

global my_cmap
from matplotlib import colors,cbook,cm

cdict = {
'red': ((0, 0.1, 0.1),(0.04, 0.15, 0.15),(0.07, 0.25, 0.25), (0.12, 0.4, 0.4),\
    (0.2, 0.6, 0.6),(0.35, 0.9, 0.9),(0.4, 0.9, 0.9),(0.6, 1, 1),(0.8, 1, 1),(0.9, 1, 1),(1, 1, 1)),
'green': ((0, 0.1, 0.1),(0.04, 0.15, 0.15),(0.07, 0.2, 0.2), (0.12, 0.15,0.15 ),\
    (0.2, 0.1, 0.1),(0.35, 0.1, 0.1),(0.4, 0.1, 0.1),(0.6, 0.6, 0.6),(0.8, 1, 1),(0.9, 1, 1),(1, 1, 1)),
'blue' : ((0, 0.3, 0.3),(0.04, 0.5, 0.5),(0.07, 0.6, 0.6), (0.12, 0.6,0.6),\
    (0.2, 0.6, 0.6),(0.35, 0.6, 0.6),(0.4, 0.1, 0.1), (0.6, 0, 0),(0.8, 0, 0),(0.9, 0.8, 0.8),(1, 1, 1))
}
my_cmap_ = colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
my_cmap = my_cmap_

     

class LogFormatterTeXExponent(FormatStrFormatter, object):
    """Extends pylab.LogFormatter to use 
    tex notation for tick labels."""
    
    def __init__(self, *args, **kwargs):
        super(LogFormatterTeXExponent, 
              self).__init__(*args, **kwargs)
        
    def __call__(self, *args, **kwargs):
        """Wrap call to parent class with 
        change to tex notation."""
        label = super(LogFormatterTeXExponent, 
                      self).__call__(*args, **kwargs)
        
        import re
        label = label.replace('--','-')
        x = float(label)
        label = '%.3g'%x

     
        if 'e' in label:
            if label[:2] == '1e':
                label_ = '10^{%d}'%int(label[2:])
            else:
                label_ = re.sub(r'e(\S)0?(\d+)',r'\\cdot 10^{\1\2}',str(label))
            label_.replace('.','\!.\!')

            label_ = "$" + label_ + "$"
            label_ = label_.replace('+','')
            #print label_
        else:
            label_ = '$%g$'%x

        return label_
        

class symlogNorm(Normalize):    
    def __init__(self,linthresh,vmin=None,vmax=None,clip=False):
        Normalize.__init__(self,vmin,vmax,clip)
        self.linthresh=float(linthresh)
        self.vmax=float(vmax)
        self.vmin=float(vmin)


    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax = self.vmin, self.vmax          
        if vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)
            # ma division is very slow; we can take a shortcut
            resdat = result.data;   resdat/=self.linthresh
            resdat = arcsinh(resdat,out=resdat)
            vmin = arcsinh(vmin/self.linthresh)
            vmax = arcsinh(vmax/self.linthresh)
            resdat-= vmin
            resdat/= vmax-vmin
  
            result = ma.array(resdat, mask=result.mask, copy=False)

        if is_scalar:
            result = result[0]

        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax = self.vmin, self.vmax
        vmin = arcsinh(vmin/self.linthresh)
        vmax = arcsinh(vmax/self.linthresh)
        val = asarray(value)
        try:
            val*=(vmax-vmin);  val+= vmin
            val = sinh(val,out= val)
            val*= self.linthresh
            return ma.asarray(val)
        except:
            print('problem!')
            return val



class MyFormatter(ScalarFormatter):   # improved format of axis
    def __call__(self, x, pos=None):
        self.set_scientific(True)
        self.set_useOffset(True)
        self.set_powerlimits((-3, 3))
        return ScalarFormatter.__call__(self, x, pos)


#help share large global variables between processes on Windows
def glob_initializer( plot_details_):
    global  plot_details
    plot_details = plot_details_


def make_graphs(input_data, plot_svd = False):
    """
    Prepare graphs for GUI and export. This is main plotting module. It can plot via Gnuplot (fast) and matplotlib (slow, better quality).

    Load and prepare data for plotting, save important values on disk, and plot ...

    :param bool _plot_svd: If `True` substract first mode of SVD
    :param dict inputs:  dictionary with all setting generated by `main` and `loadSetting`
    :param class tokamak: ``Class`` of  tokamak with all important data.
    :var array G:    Array of all results
    :var array chi2:  Array of final chi2 used during reconstruction


    """
    
     
    global inputs, tokamak, plot_details, my_cmap

    inputs, tokamak, progress,results = input_data
    output_path = inputs['output_path']

    tmp_folder = inputs['tmp_folder']
     
    shot = inputs['shot']

    my_cmap = my_cmap_ if inputs['cmap'] == 'ingesson' else cm.get_cmap(inputs['cmap'])
    my_cmap._init()

    #NOTE for transparent cmaps 
    #from matplotlib.colors import ListedColormap
    ## Get the colormap colors
    #my_cmap = my_cmap_(arange(my_cmap_.N))

    ## Set alpha
    #my_cmap[0,-1] = 0

    ## Create new colormap
    #my_cmap = ListedColormap(my_cmap)
        
    if inputs['blacken_negative']: #regions of the nagative values will be black
        my_cmap.set_under('k')


    tmin = inputs['tmin']
    tmax = inputs['tmax']
    print('PLOTTING ..'+ str(shot))

    dets = tokamak.dets
    data = results['data_all']
    error = results['error_all']
    Tmat = results['Tmat_all']
    tvec = results['tvec']
    G = results['g'].T
    G*= tokamak.norm    #convert from W/m^2/cm na W/m^3
    G_samples = results.get('g_samples',None)
    if G_samples is not None:
        G_samples = G_samples.T
        G_samples*= tokamak.norm    #convert from W/m^2/cm na W/m^3
        
    

    data[:,results['dets']]  = results['data']  #include data correction (bacground substraction) done during preprocessing  
    error[:,results['dets']] = results['error']  #include data correction (bacground substraction) done during preprocessing  
    retro = zeros_like(data)
    retro[:,results['dets']] = results['retro']
    
    chi2 = results['chi2']
    lam = results['lam0']

    rem_dets = where(~in1d(arange(data.shape[1]), results['dets']))[0]
 
    try:
        retro[:,rem_dets] = (Tmat[rem_dets]*G).T
    except:
        print('error retro[:,rem_dets] = (Tmat[rem_dets]*G).T', Tmat.shape, rem_dets,shape, G.shape) 
        import IPython
        IPython.embed()
    
    

    #=========================================
    dets = dets[tokamak.get_correct_dets(data.T)]   
    #=========================================


    tsteps = len(tvec)

    inputs['tsteps'] = tsteps   # used in GUI.py
    base = '_SVD_' if plot_svd else '_rec_'


    if plot_svd and tsteps < 4:
        print("Not enough time slices for SVD")
        return



    t_plot = time.time()

    rem_fsa = inputs['rem_fsa']
    rem_fsa &= tokamak.transform_index == 1
    symemtric_cmap = False
    if (inputs['rem_back'] or plot_svd or rem_fsa) and not inputs['cmap_neg'] == 'none':
        my_cmap = cm.get_cmap(inputs['cmap_neg'])
        my_cmap._init()
        symemtric_cmap = True
 
        
    if tokamak.transform_index == 1: 
        if rem_fsa:
            raise Exception('implementation not finished!!')

   
    if inputs['rem_back']:  # remove background
        print("REMOVING BACKGROUND !!!")

        subst_ind = slice(None,int(ceil(tsteps*inputs['bcg_subst_fract'])))
        data_0 = average(data[subst_ind,:],weights=1/error[subst_ind,:]**2+1e-6,axis=0)
        dets_ = results['dets']
        try:
            #data from the phantom
            Emiss0 = load(tmp_folder+'/Emiss0.npz')
            G0 = Emiss0['G'].reshape(tokamak.ny*tokamak.nx,tsteps, order='F')
            #normalize it to one!
            vmax = G0.max()
            G0/= vmax
            data/= vmax
            retro/= vmax
            G/= vmax
            data0 = (tokamak.Tmat*G0).T

            print('remove phantom')
        except:
            G0 = mean( G[:,subst_ind],axis=1)[:,None]
            data0 = median(data[subst_ind],axis=0)[None,:]
            
            
        retro0 =  median(retro[subst_ind,:],axis=0)[None,:]


 
        data -= data0
        retro -= retro0
        retro[:,dets_] += (retro0-data0)[:,dets_]
        G  -= G0
        savez(tmp_folder+'/'+'data0_'+str(shot),data0=data_0)
        savez(tmp_folder+'/'+'G0_'+str(shot),  G0=G0)
    #from IPython import embed
    #embed()
    gres=G.reshape(tokamak.ny,tokamak.nx,tsteps, order='F')
    if G_samples is not None:
        gres_samples=G_samples.reshape(tokamak.ny,tokamak.nx,G_samples.shape[1],-1, order='F')
 
        
    if plot_svd:
        n_svd_rem = 1
        ur,sr,vr = fast_svd(G,n_svd_rem)

        G0 = dot(ur*sr,vr)
        retro0 = (tokamak.Tmat*G0).T
        G_svd = G-G0
        data -= retro0
        retro -= retro0
        gres=reshape(G_svd,(tokamak.ny,tokamak.nx,tsteps), order='F')
        plot_autoscale=False
        

    dets_dict = channels = None
    if hasattr(tokamak, 'detectors_dict'):
        detectors_dict = tokamak.detectors_dict
        channels = [j for cam, ch in detectors_dict.items() for j in ch]
  
    savez_compressed(tmp_folder+'/'+'data_'+str(shot),tvec=tvec, retro=retro, 
                data=data, err=error,dets=dets,dets_dict=dets_dict,channels=channels)
    savez_compressed(tmp_folder+'/'+'convergence_'+str(shot), tvec=tvec,chi2=chi2,lam=lam)
    
    if inputs['enable_output']:
        copyfile(tmp_folder+'/'+'data_'+str(shot)+'.npz',output_path+'/data_'+str(shot)+'.npz')
        copyfile(tmp_folder+'/'+'convergence_'+str(shot)+'.npz',output_path+'/convergence_'+str(shot)+'.npz')
        
        

    #what was this step for???  it is usually doing nothing
    xslice  = slice( max(0, int(round((tokamak.plot_coord[0] - tokamak.xmin)/tokamak.dx)) ), min(tokamak.nx, int(round((tokamak.plot_coord[1]-tokamak.xmin)/tokamak.dx)) ) )
    yslice  = slice( max(0, int(round((tokamak.plot_coord[2] - tokamak.ymin)/tokamak.dy)) ), min(tokamak.ny, int(round((tokamak.plot_coord[3]-tokamak.ymin)/tokamak.dy)) ) )
    gres = gres[yslice, xslice]
    if G_samples is not None:
        gres_samples = gres_samples[yslice, xslice]

    
    ny, nx, nt = shape(gres)
    xmin, xmax, ymin, ymax = tokamak.plot_coord
    
    xpix=linspace(xmin,xmax, nx)/tokamak.norm
    ypix=linspace(ymin,ymax, ny)/tokamak.norm
    rhop,magx_all, magy_all = tokamak.mag_equilibrium(tvec,return_mean=True)


    ##============store equilibrium =================
    #print 'save HR equilibrium ' fc

    #rho = linspace(0,1.005, 200)
    #magx, magy = tokamak.mag_equilibrium(tvec,return_mean=True,n_rho=200,n_theta=200,rho=rho)

    #theta_star = tokamak.mag_theta_star(tvec.mean(), rho, magx, magy, rz_grid=False )
    #theta_star_rz = tokamak.mag_theta_star(tvec.mean(), rho[:-1], magx[:,:-1], magy[:,:-1], rz_grid=True )

    #M = get_rho_field_mat(tokamak,  mean(tvec) )

    #savez_compressed(local_path+'tmp/'+'mag_'+str(shot),magx=magx[:,:-1],magy = magy[:,:-1], 
          #theta_star=theta_star[:-1], rho_pol=M,theta_star_rz=theta_star_rz,BdMat=BdMat)
    #print 'done save HR equilibrium '
    
    #==============================================

    BdMat = get_bd_mat(tokamak,time=tvec.mean())

    #save the data 
    try:
        gresnorm = 1
        if inputs['save_profiles'] or inputs['enable_output']:
            #print('saving emissivity ... ')
            print('saving tomography ...   ', end = '')

            BdMat = BdMat.reshape((tokamak.ny, tokamak.nx), order="F")
            
            xgrid = (tokamak.xgrid+tokamak.dx/2)/tokamak.norm  #centers of the pixels
            ygrid = (tokamak.ygrid+tokamak.dy/2)/tokamak.norm  #centers of the pixels
            gresnorm = maximum((1+gres.max(0).max(0)),(1-gres.min(0).min(0)))/65504
            gres/= gresnorm[None,None,:]

            if G_samples is not None:
                gres_samples_norm = maximum((1+gres_samples.max(0).max(0)),(1-gres_samples.min(0).min(0)))/65504
                gres_samples /= gres_samples_norm[None, None]
                gres_samples = gres_samples.astype(float16)
            else:
                gres_samples_norm = None
                gres_samples  = None
                
            tokamak_tmp = inputs.pop('tokamak_tmp') #do not save it!
            inputs['impur_inject_t'] = tokamak.impur_inject_t
            name = 'Emissivity_%.3f-%.3f'%(tvec[0], tvec[-1])+base+str(shot)
            savez_compressed(tmp_folder+'/'+name,gres=gres.astype(float16),tvec=tvec,\
                        rvec=xgrid,zvec=ygrid,inputs=inputs,gres_norm=gresnorm,BdMat=BdMat,
                        gres_samples=gres_samples, gres_samples_norm=gres_samples_norm)#,Rho=Rho)
            if inputs['enable_output']:
                copyfile(tmp_folder+'/'+name+'.npz',output_path+'/'+name+'.npz')
                if inputs['rem_back']:
                    copyfile(tmp_folder+'/'+'G0_'+str(shot)+'.npz',output_path+'/G0_'+str(shot)+'.npz')
                    copyfile(tmp_folder+'/'+'data0_'+str(shot)+'.npz',output_path+'/data0_'+str(shot)+'.npz')
            else:
                inputs['tokamak_tmp'] = tokamak_tmp
  
            print('done')
    except Exception as e:
        print('saving of radiation profile has failured', e)
    finally:
        if inputs['save_profiles'] or inputs['enable_output']:
            gres*= gresnorm[None,None,:]
            
  

    #add padding to to graphs
    nx_new = (nx*(ymax-ymin)/(xmax-xmin))
    padding_size = int(max(round((nx_new-nx+1)/2), 0))
    padding = zeros((ny, padding_size),dtype=gres.dtype)
    padding*= nan
    
    nx += padding_size*2
    #correction for the padding
    xmin -= padding_size*tokamak.dx
    xmax += padding_size*tokamak.dx
    xpix = linspace(xmin,xmax, nx)/tokamak.norm
    ypix = linspace(ymin,ymax, ny)/tokamak.norm
    


    #dictionary with all setings for plotting
    plot_details = {}
    plot_details['geometry'] = array((xmin,xmax,ymin,ymax))/tokamak.norm
    plot_details['mat_size'] = array((nx,ny))
    plot_details['output_path'] = inputs['output_path']
    plot_details['local_path'] = inputs['local_path']
    plot_details['enable_output'] = inputs['enable_output']
    plot_details['geometry_path'] = tokamak.geometry_path
    plot_details['shot'] = tokamak.shot
    plot_details['detectors_dict'] = getattr(tokamak,'detectors_dict',None)
    plot_details['dets_index'] = tokamak.dets_index
    plot_details['input_diagn'] = tokamak.input_diagn
    plot_details['name'] = tokamak.name
    plot_details['transform_index'] = tokamak.transform_index
    plot_details['Xchords'] = tokamak.Xchords
    plot_details['Ychords'] = tokamak.Ychords
    plot_details['dx'] = tokamak.dx
    plot_details['dy'] = tokamak.dy
    plot_details['norm'] = tokamak.norm
    plot_details['struct_dict'] = getattr(tokamak,'struct_dict',None)
    plot_details['transform_index'] = getattr(tokamak,'transform_index',None)
    plot_details['xmin'] = tokamak.xmin
    plot_details['xmax'] = tokamak.xmax
    plot_details['ymin'] = tokamak.ymin
    plot_details['ymax'] = tokamak.ymax
    plot_details['nx'] = tokamak.nx
    plot_details['ny'] = tokamak.ny
    plot_details['ygrid'] = tokamak.ygrid
    plot_details['xgrid'] = tokamak.xgrid

    plot_details['plot_autoscale'] = inputs['plot_autoscale']
    plot_details['plot_svd'] = plot_svd
    plot_details['plot_contours'] = inputs['plot_contours']
    plot_details['n_contours'] = inputs['n_contours']
    plot_details['boundary'] =  results['bnd']
    plot_details['blacken_negative'] =  inputs['blacken_negative']
    plot_details['show_fit_residuum'] = inputs['show_fit_residuum']
    plot_details['dpi'] = inputs['dpi']
    plot_details['img_size'] = inputs['img_size']
    
    plot_details['rem_back'] = inputs['rem_back']
    plot_details['rem_fsa'] = inputs['rem_fsa']
    plot_details['output_type'] = inputs['output_type']
    plot_details['tmp_folder']  = inputs['tmp_folder']
    plot_details['plot_chords'] = inputs['plot_chords']

    if hasattr(tokamak, 'ICRH_rezonance'):
        plot_details['ICRH_position'] = tokamak.ICRH_rezonance

    savez_compressed(tmp_folder+'/'+'mag_'+str(shot),magx=single(magx_all),magy = single(magy_all), tvec=tvec)
    
    

   
    

    if inputs.get('plot_surfaces',False) and inputs['plot_theta_star']:
        THETA = []
        MAGX = []
        MAGY = []

            #============store equilibrium =================
        

        rho = linspace(0,1.005, 200)**3
        n_rho = sum(rho<=1.)
        theta_star_contours_R = zeros((tsteps,n_rho, 20))
        theta_star_contours_Z = zeros((tsteps,n_rho, 20))
        from tqdm import tqdm,trange
        for it in trange(tsteps,desc='Calculate theta star: '):
            rhop_,magx_, magy_ = tokamak.mag_equilibrium(tvec[it],return_mean=True,n_theta=200,rho=rho)

            theta_star = tokamak.mag_theta_star(tvec[it], rhop_, squeeze(magx_),squeeze(magy_), rz_grid=False )
            
            THETA.append(theta_star)
            MAGX.append(magx_)
            MAGY.append(magy_)

            theta = linspace(0,2*pi, 20, endpoint=False)
            for i in range(n_rho):
                theta_star_contours_R[it,i] = interp(theta,theta_star[i],squeeze(magx_[:,i]))
                theta_star_contours_Z[it,i] = interp(theta,theta_star[i],squeeze(magy_[:,i]))

        plot_details['theta_star'] =  theta_star_contours_R,theta_star_contours_Z
        savez_compressed(tmp_folder+'/theta_star', theta = single(THETA),
                         magx=single(MAGX), magy = single(MAGY), rhop = rhop_, tvec=tvec )
    
    if (inputs['plot_all'] or tsteps == 1) and inputs['plot_surfaces']:

        
        if 'magx' in results:
            skip_mag = max(size(results['magx'], 1)//10, 1)    #use at most 10 lines in the graph
            ind = slice(skip_mag-1,None,skip_mag)
            magx_all = results['magx'][:,ind].T
            magy_all = results['magy'][:,ind].T
            rhop = results['mag_rho'][ind]
        else:
            ind = slice(6-1,None,6)
            rhop,magx_all,magy_all = tokamak.mag_equilibrium(tvec,surf_slice=ind,n_rho=11)

        plot_details['mag_field'] = True
        plot_details['mag_field_rhop'] = rhop

        if 'plot_surfaces_outer' in inputs and inputs['plot_surfaces_outer'] and tokamak.use_pfm:

            rho_out = linspace(1,1.1,5)
            magr_out, magz_out = tokamak.get_mag_contours(tvec,rho_out)
            plot_details['mag_out'] = magr_out, magz_out
  
    elif 'mag_field' in plot_details:
        plot_details.pop('mag_field')
        

  
    
    

    prof_norm = None
    data_norm = None
    if inputs['plot_loglin']: 
        #estimate a normalization
        ind = slice(None,None) if tsteps < 1000 else random.randint(tsteps,size=1000)       
        data_norm = median(abs(data[:,results['dets']]).mean(0))/2/inputs['loglin_threshold']
        prof_norm = std(G[:,ind])/inputs['loglin_threshold']

    vmin,vmax = mquantiles(gres.max(0).max(0), (0.02, 0.98))
    vmin = min(0, vmin)
    
    if symemtric_cmap:
        vmax = max(-vmin,vmax)
        vmin = min(vmin,-vmax)
    
    if inputs['plot_all']  or tsteps == 1:
  
        try:
            n_cpu = cpu_count()
            #use initializer to share 
            p = Pool(n_cpu)
        except Exception as e:
            print('Pool initialization error: ', e)
            n_cpu = 1

        if tsteps == 1:
            name = ['previewA'+str(shot)+base, 'previewB'+str(shot)+base]
            namesA,namesB = [name[0],],[name[1],]
        else:
            namesA = ['brightness_%d_%.4d'%(shot,ts)+base for ts in range(tsteps)]
            namesB = ['emissivity_%d_%.4d'%(shot,ts)+base for ts in range(tsteps)]

        try:
            if not rem_back:
                vmax = load(tmp_folder+'/Emiss0.npz')['G0'].max()
        except: 
            pass 
            
        data_min,data_max = mquantiles(data[:,dets][isfinite(error[:,dets])],[.01,.99])
        glim = vmin, vmax
        limit = min(0,data_min), max(data_max,amax(retro))
  


        if 'mag_field' in plot_details: #BUG is ot working if False
            mag_field = array((magx_all,magy_all),copy=False)
            mag_field/= tokamak.norm
        else:
            mag_field = nan*ones(tsteps)
            
        #high quality matplotlib output
        n_split = 1 if tsteps < 2*n_cpu else 2*n_cpu
        ind = [slice(i*tsteps//n_split,(i+1)*tsteps//n_split) for i in range(n_split)]

        args = [(r_[ii],gres[...,ii],padding,data[ii,:], error[ii,:],retro[ii,:],dets
                ,tvec,my_cmap,mag_field[...,ii],chi2[ii],lam[ii],tsteps,tokamak.rho_label,
                base,(data_norm,prof_norm),(limit,glim),plot_details) for ii in ind]
            

        if tokamak.input_diagn in ['BOLO', ]:
            inputs['plot_data_rho_tg'] = False
        plot_method = matplotlib_data_tg if inputs['plot_data_rho_tg'] else matplotlib_data


        if rcParams['backend'].lower() == 'agg' and not config.DEBUG:
            p.imap(plot_method,args,1 )
            p.imap(matplotlib_image,args,1 )
    

        else:
            from tqdm import tqdm
            for a in tqdm(args,desc='single thread plotting: '):
                plot_method(a)
                matplotlib_image(a)
          
        try:
            p.close()
            p.join()
        except:
            pass
      

    #convert *1_rec_.png -background white -alpha remove %04d.png 
    #mencoder "mf://*.png" -o ./movie_30579.avi -ovc lavc -lavcopts vcodec=mpeg4:mbd=1:vbitrate=2800
    #ffmpeg -r 60 -f image2 -s 1920x1080 -i %2d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test.mp4
    #mencoder "mf://*.png" -o ./movie_30579.avi -ovc lavc -lavcopts vcodec==mpeg2:tsaf:vbitrate=8000
        #print tsteps > 50  , make_movie , plot_all
        if tsteps > 50 and inputs['make_movie'] and inputs['plot_all'] :
            debug( '\nmake a movie')

            #BUG size of image is hardcoded :( 
            movie_name = output_path+'movie_%d%s.mp4'%(shot,base)
            if os.path.exists(movie_name):
                os.remove(movie_name)
            
            #os.remove(output_path+'/_movie_%d_%s.mp4'%(shot,base))
            os.system('ffmpeg   -f lavfi -i color=c=white:s=544x451 -i '+tmp_folder+'/emissivity_'+str(shot)+'_%04d'+base+'.png -filter_complex "[0:v][1:v]overlay=shortest=1,format=yuv420p[out]" -map "[out]" '+movie_name)

            #from subprocess import Popen, PIPE, STDOUT

            #import IPython
            #IPython.embed()
            #from subprocess import Popen, PIPE, STDOUT
            #mencoder "mf://'+tmp_folder+'/*1_rec_.png" -o movie_29624.avi -ovc lavc -lavcopts vcodec=mjpeg
            #os.system('mencoder "mf://'+tmp_folder+'/*1'+base+'.png" -o '+output_path\
                #+'/movie_%d_%.4f-%.4f%s.avi -ovc lavc -lavcopts vcodec=mpeg4:mbd=1:vbitrate=2000'%(shot,tvec[0], tvec[-1],base))
            #try:
                 #Popen(['ffmpeg',' -f lavfi -i color=c=white:s=544x451 -i '+tmp_folder+'/emissivity_'+str(shot)+'_%04d'+base+
                        #'.png -filter_complex "[0:v][1:v]overlay=shortest=1,format=yuv420p[out]" -map "[out]" '+
                        #output_path+'/movie_%d_%s.mp4'%(shot,base)],stdin=PIPE, stderr=PIPE,stdout=PIPE)
                 
                            #Popen(['ffmpeg',' -f lavfi -i color=c=white:s=544x451 -i '+tmp_folder+'emissivity_'+str(shot)+'_%04d'+base+'.png -filter_complex "[0:v][1:v]overlay=shortest=1,format=yuv420p[out]" -map "[out]" '+ movie_name],stdin=PIPE,)
                 
                #Popen(['mencoder',' "mf://'+tmp_folder+'emissivity_'+str(shot)+'_%04d'+base+'.png" -o '   +output_path+'/movie_%d.avi -ovc lavc -lavcopts vcodec=mjpeg'])
                 
                 #os.system('mencoder "mf:/'+tmp_folder+'emissivity_*'+base+'.png" -o '   +output_path+'/movie_%d.avi  -ovc lavc -lavcopts vcodec=mpeg4:mbd=1:vbitrate=2800')
                 
            #except subprocess.CalledProcessError as e:
                #print 'return code :',e.returncode
                #print 'error message:',e.output
                
                
                #os.system( 'mencoder "mf://'+tmp_folder+'/*1_rec_.png" -o movie_29624.avi -ovc lavc -lavcopts vcodec=mjpeg')
   #mencoder "mf://*.png" -o ./movie_30579.avi  -ovc x264 -x264encopts pass=1:preset=veryslow:fast_pskip=0:tune=film:frameref=10:bitrate=3000:threads=auto 

    if  tsteps > 1:
        #try: #cut throught magnetics axis
        Rc,Zc = magx_all[:,0].mean(), magy_all[:,0].mean()
        xfract = int_((Rc-(xmin+padding_size*tokamak.dx))/(xmax-padding_size*tokamak.dx-(xmin+padding_size*tokamak.dx))*(nx-2*padding_size))
        yfract = int_((Zc-ymin)/(ymax-ymin)*ny)

        g_profilx = gres[yfract,:,arange(tsteps)]     
        g_profily = gres[:,xfract,arange(tsteps)].T
    
        converged = '     Converged: %d/%d'%(sum(abs(chi2) < 0.1), tsteps)
        mean_chi =  '     Mean Chi2: %.2f' % exp(nansum(log(chi2))/len(chi2))
        name = ('previewA'+str(shot)+base, 'previewB'+str(shot)+base)
        norm = Normalize()
        if inputs['plot_loglin']:
            norm = symlogNorm(prof_norm,vmin=vmin,vmax=vmax)
            
        if symemtric_cmap:
            if inputs['plot_loglin']:
                norm = symlogNorm(prof_norm,vmin=vmin,vmax=vmax)
            else:
                norm = Normalize(vmin=vmin,vmax=vmax)
  
        matplotlib_preview(g_profilx, g_profily,padding,data, dets, 'Profile R'+converged, 
                    'Profile Z'+mean_chi, tokamak.t_name, tvec , my_cmap,norm,  name,inputs['plot_loglin'],
                    symemtric_cmap)



    print('\rPlotting time %.1f s' % (time.time()-t_plot))


def update_errorbar(err_plot, x,y,yerr):
    
    plotline, caplines, barlinecols = err_plot

    # Replot the data first
    plotline.set_data(x,y)

    # Find the ending points of the errorbars
    error_positions = (x,y-yerr), (x,y+yerr)

    # Update the caplines
    if len(caplines) > 0:
        for j,pos in enumerate(error_positions):
            caplines[j].set_data(pos)

    # Update the error bars
    barlinecols[0].set_segments(list(zip(list(zip(x,y-yerr)), list(zip(x,y+yerr))))) 
        



 

def matplotlib_preview(g_profilx, g_profily, padding, data,dets, titleX, titleY,
                       t_name, tvec , my_cmap,norm, filename, plot_loglin,sym_colorbar  ):
    """
    Make preview of recontruction using the Matplotlib library, it is slower than GNUPLOT but it has higher quality
    """

    global  plot_details

    enable_output = plot_details['enable_output']
    geometry = plot_details['geometry']


    shot = plot_details['shot']
    xgrid = plot_details['xgrid']/plot_details['norm']
    ygrid = plot_details['ygrid']/plot_details['norm']
    
    under = max(1,len(tvec)//500)
    tvec = tvec[:,None]

    g_profilx = g_profilx[::under,:]
    g_profily = g_profily[::under,:]
    padding = tile(padding[0,:], (g_profilx.shape[0],1))
    g_profilx = concatenate([padding, g_profilx, padding], axis=1)
    
    if all(g_profilx==0) or all(g_profily==0):
        raise Exception('Zero cut of the emissivity profile')
    
    gmax = max(nanmax(g_profilx), nanmax(g_profily))
    if gmax > 3e5:
        fact,pre = 1e6, 'M'
    elif gmax> 3e2:
        fact,pre = 1e3, 'k'
    else:
        fact,pre = 1e0, ''
        
        
    if plot_loglin:
        fact,pre = 1e0, ''

    plot_2D_adv(xgrid,tvec,g_profilx/fact,filename[0],title="Emissivity cut in R-plane",
            xlabel='Position [m]',ylabel='Time ['+t_name+']' ,sym_colorbar=sym_colorbar,
            clabel='Emissivity [%sW/m$^3$]'%pre, cmap=my_cmap,norm=norm)
    
    plot_2D_adv(xgrid,tvec,g_profily/fact,filename[1],title="Emissivity cut in Z-plane",
            xlabel='Position [m]',ylabel='Time ['+t_name+']',sym_colorbar=sym_colorbar,
            clabel='Emissivity [%sW/m$^3$]'%pre, cmap=my_cmap,norm=norm)
    
    if enable_output:
        plot_2D_adv(xgrid,tvec,g_profilx/fact,'g_profilX_'+str(shot),title="Emissivity cut in R-plane",
                xlabel='Position [m]',ylabel='Time ['+t_name+']',sym_colorbar=sym_colorbar, 
                clabel='Emissivity [%sW/m$^3$]'%pre, cmap=my_cmap,norm=norm)
        
        plot_2D_adv(xgrid,tvec,g_profily/fact,'g_profilY_'+str(shot),title="Emissivity cut in Z-plane",
                xlabel='Position [m]',ylabel='Time ['+t_name+']',sym_colorbar=sym_colorbar,
                clabel='Emissivity [%sW/m$^3$]'%pre, cmap=my_cmap,norm=norm)

  
def matplotlib_data_tg(params):
    #advacent plotting method for the showing  the measured brighntess and back calculated values at function of tangential radius

    try:
        #clean catch in foked processes
        matplotlib.font_manager._get_font.cache_clear()
    except:
        pass
    
    (ntvec,G,padding, Data,Error,Fc,dets,tvec,my_cmap,
                    mag_field,chi2,lam,tsteps,rho_label,base,linthresh,limits,plot_details) = params
    T = time.time()
    from collections import OrderedDict
    from scipy.ndimage.interpolation import map_coordinates
    from shared_modules import MovingAveradge,extrap1d

    if os.name != 'nt':
        os.nice(3)
    
    shot = plot_details['shot']
    enable_output=   plot_details['enable_output']
    geometry=   plot_details['geometry']
    mat_size=   plot_details['mat_size']
    output_path=   plot_details['output_path']
    local_path=   plot_details['local_path']
    boundary=   plot_details['boundary']
    plot_contours=   plot_details['plot_contours']
    n_contours=   plot_details['n_contours']
    dpi=   plot_details['dpi']
    img_size=   plot_details['img_size']


    magx_all = mag_field[0]
    magy_all = mag_field[1]
    rhop = plot_details['mag_field_rhop']

    plot_loglin = linthresh[0]!= None
    dmin,dmax = limits[0]
    output_type = plot_details['output_type']
    tmp_folder = plot_details['tmp_folder']
    
    if plot_details['rem_back'] or plot_details['plot_svd']:
        plot_details['show_fit_residuum'] = False

            

    dets_dict = plot_details['detectors_dict']
    cam_ind = OrderedDict()

    for det, key in zip(plot_details['dets_index'],list(dets_dict.keys())):
        if len(key)==2 and key[-1].isdigit(): key =  key[ 0]  #convect SXR different cameras from AUG together 
        if not key in  cam_ind: cam_ind[key] = []
        cam_ind[key].append(det)
    #jont some comeraras from the same poloidal port
    if plot_details['name'] == 'ASDEX' and plot_details['input_diagn'][:3] == 'SXR' and 'F' in cam_ind:
        cam_ind_ = OrderedDict()
        cam_ind_['F+G'] = cam_ind.pop('F')+cam_ind.pop('G')
        cam_ind_.update(cam_ind)
        cam_ind = cam_ind_
        
    if plot_details['name'] == 'ASDEX' and plot_details['input_diagn'] =='AXUV':
        cam_ind_ = OrderedDict()
        cam_ind_['D13+DVC'] = cam_ind.pop('D13')+cam_ind.pop('DVC')
        cam_ind_.update(cam_ind)
        cam_ind = cam_ind_        
 
    if plot_details['name'] == 'DIIID' and '45R1' in cam_ind and '195R1' in cam_ind:
        cam_ind_ = OrderedDict()
        cam_ind_['45R1+195R1'] = cam_ind.pop('45R1')+cam_ind.pop('195R1')
        cam_ind_.update(cam_ind)
        cam_ind = cam_ind_
        
    if plot_details['name'] == 'DIIID' and '90RP1a' in cam_ind and '90RP1b' in cam_ind :
        cam_ind_ = OrderedDict()
        cam_ind_['90RP1a+90RP1b'] = cam_ind.pop('90RP1a')+cam_ind.pop('90RP1b')
        cam_ind_['90RM1a+90RM1b'] = cam_ind.pop('90RM1a')+cam_ind.pop('90RM1b')
        cam_ind_.update(cam_ind)
        cam_ind = cam_ind_ 
        
    n_col = len(cam_ind)
    n_row = 1
    resid0 = infty
    for i_col in range(int(ceil(sqrt(len(cam_ind))))+1):
        for i_row in range(i_col+1):
            resid = (i_col+1)*(i_row+1) - len(cam_ind)
            if resid < 0: continue
            if resid<resid0:
                resid0,n_col,n_row = resid,i_col+1,i_row+1
            

    
    #prepare plots and axis 
    
    f = Figure((img_size*1.5,img_size))
    FigureCanvas(f)
    f.subplots_adjust(hspace=0.05, wspace = 0.10)
    axis = OrderedDict()
    ax = None
    c = 'k','b','g'

    for i,(det, inds) in enumerate(cam_ind.items()):
        ax = f.add_subplot(n_row, n_col,i+1,sharex=ax,sharey=ax )
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.get_yaxis().set_tick_params(direction='in')
        ax.get_xaxis().set_tick_params(direction='in')
        axis[ax] = {}
        ax.xaxis.grid(True)
        ax.set_xlim(-1.1,1.1)
        if plot_details['input_diagn'] == 'AXUV' and plot_details['name']=='ASDEX':
            ax.set_xlim(-1.2,1.2)

        fcmax = abs(Fc).max()  
        if  fcmax > 2e8:
            fac,pre = 1,''
        elif fcmax > 2e5:
            fac,pre = 1e6,'M'
        elif fcmax > 2e2:
            fac,pre = 1e3,'k'
        else:
            fac,pre = 1.,''

        ax.axhline(y=0,c='k',zorder=0)

        if plot_loglin:
            ax.set_yscale('symlog', linthreshy=linthresh[0]*2)
            fac,pre = 1.,''
            
        
        if not plot_details['plot_autoscale']:
            ax.set_ylim((-amax(Fc[:,dets])*.05+dmin)/fac ,dmax/fac)
            
        if plot_details['rem_back'] or (plot_details['rem_fsa'] and plot_details['transform_index']==1) or plot_details['plot_svd']:
            ax.set_ylim(dmin/fac,dmax/fac)

        if i < len(cam_ind)-n_col:
            ax.xaxis.offsetText.set_visible(False)
            for label in ax.get_xticklabels():
                label.set_visible(False)
        else:
            ax.set_xlabel('tangent %s'%rho_label)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        if i%n_col != 0:
            ax.yaxis.offsetText.set_visible(False)
            for label in ax.get_yticklabels():
                label.set_visible(False)
        else:
            ax.set_ylabel('Brightness [%sW/m$^2$]'%pre)

        axis[ax]['retro_spline'], = ax.plot([],[],'k')

        axis[ax]['retro_line'] = []
        axis[ax]['retro_points'] = []
        axis[ax]['data_points'] = []
        axis[ax]['data_points_removed'] = []
        axis[ax]['resid'] = []
        axis[ax]['retro_spline'] = ax.plot([],[],'k')[0]

        for j,ind in enumerate(inds):
            axis[ax]['retro_points'].append(ax.plot([],[],c[j%3]+'x')[0])
            axis[ax]['data_points'].append(ax.errorbar(0,0,0,fmt=c[j%3]+'_'))
            axis[ax]['data_points_removed'].append(ax.errorbar(0,0,0,fmt='.',c='r'))
            axis[ax]['resid'].append(ax.errorbar(0,0,0,fmt=c[j%3]+'_',capsize=0))
            
        if det.find('+') != -1: 
            cameras = det.split('+')
        elif plot_details['input_diagn'][:3] == 'SXR' and plot_details['name']=='ASDEX':
            cameras = [k for k in list(dets_dict.keys()) if k[0]==det]
            if det.find('+')!= -1: cameras = det.split('+')
        else:
            cameras = det,
                    
        for k,(cam, ic) in enumerate(zip(cameras, c)):
            ax.text(0.3, 0.2-k/15., cam, transform=ax.transAxes,color=ic)#,backgroundcolor='w')

    G = G.reshape(plot_details['ny'], plot_details['nx'],-1, order='F')

    #prepare data for plotting and the plotting  
    for it, tt in enumerate(ntvec):
        
        t = time.time()
 
        rho_tangent,chordsx,chordsy = get_rho_tan(rhop,magx_all[:,:,it],
                            magy_all[:,:,it],plot_details['Xchords'], plot_details['Ychords']) #40ms
  
        fc = Fc[it]
        data = Data[it]
        error = Error[it]
        

        scaling = array([plot_details['dx']  ,plot_details['dy']  ])
        offset  = array([plot_details['xmin'],plot_details['ymin']])+scaling/2
        
        ax = list(axis.keys())[0]
        vmax = 0

        for i,(cam, inds) in enumerate(cam_ind.items()):
            unsample=5
            ax = list(axis.keys())[i]
            rho_tg_cam = []
            retro_interp_cam = []
            ind_ = hstack(inds)
            #orientation = sign(mean(tokamak.Ychords[-1,ind_]-tokamak.Ychords[0,ind_]))
            #orientation = 1d
            for j,ind in enumerate(inds):
                orientation = sign(mean(plot_details['Ychords'][-1,ind]-plot_details['Ychords'][0,ind]))

                correct = in1d(ind, dets)&isfinite(error[ind])
                rhotg = orientation*rho_tangent[ind]
                axis[ax]['retro_points'][j].set_data(rhotg, fc[ind]/fac)
            
                #calculate interpolated retrofits in between the measurements
                x = (arange(unsample*len(ind))-unsample/2)*1./unsample
                X = extrap1d(arange(len(ind)),chordsx[:,ind].T)(x) #0.8ms
                Y = extrap1d(arange(len(ind)),chordsy[:,ind].T)(x) #0.8ms

                coords = c_[X.ravel(),Y.ravel()]#0.1ms
                
                idx = (coords-offset)/scaling#0.6ms
                
                interpG = map_coordinates(G[...,it].T,idx.T,order=2).reshape(X.shape)  #2.4ms

                ifc = (hypot(gradient(X)[1], gradient(Y)[1])*interpG).sum(1) #0.8ms #make line integration
                

                retro_interp = MovingAveradge(copy(ifc),unsample) #broadenning due a to a finite width of LOS
                
                #retrofit will go now exactly through real retrofit values!
                retro_interp-= interp(ind[0]+x,ind,interp(ind,ind[0]+x,retro_interp))
                retro_interp+= interp(ind[0]+x,ind,Fc[it,ind])
     
                

                #estimate rho tangental for retrofit
                rhotg_ = extrap1d(ind,rhotg)(ind[0]+x)
                rho_tg_cam.append(rhotg_)
                retro_interp_cam.append(retro_interp)
        
                x = orientation*rho_tangent[ind[correct]]
                y = data[ind[correct]]/fac
                yerr = error[ind[correct]]/fac
                update_errorbar(axis[ax]['data_points'][j], x,y,yerr) #0.4ms
                if len(y):
                    vmax = max(vmax, y.max()*fac )

                x = orientation*rho_tangent[ind[~correct]]
                y = data[ind[~correct]]/fac
                yerr = error[ind[~correct]]/fac
                update_errorbar(axis[ax]['data_points_removed'][j], x,y,yerr)#0.4ms


                if plot_details['show_fit_residuum']:
                    x = orientation*rho_tangent[ind[correct]]
                    y = (data[ind[correct]]-fc[ind[correct]])/fac
                    yerr = error[ind[correct]]/fac
                    update_errorbar(axis[ax]['resid'][j], x,y,yerr)#0.4ms
                    
            #averadge profile over almost identical cameras
            Rho = linspace(amin(hstack(rho_tg_cam)),amax(hstack(rho_tg_cam)),200)
            W = zeros_like(Rho)
            Retro = zeros_like(Rho)
            for rho,retro in zip(rho_tg_cam,retro_interp_cam ):
                sort_ind = argsort(rho)
                w = interp(Rho,rho[sort_ind],ones_like(rho), left=0,right=0)
                W+= w
                Retro+= interp(Rho,rho[sort_ind],retro[sort_ind])*w
            #W==0 => gap between detector arrays
            Retro[W!=0]/= W[W!=0]            
            axis[ax]['retro_spline'].set_data(Rho[W!=0], Retro[W!=0]/fac)
            vmax = max(vmax, Retro.max() )
       
        if plot_details['plot_autoscale']:
            ax.set_ylim(-vmax/fac*.05 ,vmax/fac*1.1)
            if plot_details['rem_back'] or plot_details['plot_svd']:
                vmax = amax(abs(fc[dets]))/fac*1.1
                ax.set_ylim(-vmax ,vmax)
    
            if plot_loglin:
                ax.set_ylim(0,vmax/fac*1.1)
                    

        name = 'brightness_%d_%.4d'%(shot,ntvec[it])+base
        if tsteps==1:
            name = 'previewA'+str(shot)+base 
            
        filename = tmp_folder + '/' + name + '.png'
        s = img_size/6
        f.savefig(filename,transparent=True, dpi=dpi,bbox_inches=Bbox([[.3*s,.1*s],[8.3*s,5.5*s]]))
        
      
        if enable_output:
            name = 'brightness_%d_%.4d'%(shot,ntvec[it])+base
            f.savefig(output_path+ name + '.'+output_type,transparent=True, dpi=dpi,bbox_inches='tight')

    
    try:
        sys.stdout.write("\r %2.1f %% fps:%2.1f" %((ntvec[-1])*100./tsteps/2, size(ntvec)/(time.time()-T)))
        sys.stdout.flush()  # maybe problem with paralelization ??
    except:
        pass
        
 
def matplotlib_data(params):
    
    try:
        #clean catch in foked processes
        matplotlib.font_manager._get_font.cache_clear()
    except:
        pass
    

    (ntvec,gprof,padding, data,error,fc,dets,tvec,my_cmap,
                      mag_field,chi2,lam,tsteps,rho_label,base,linthresh,limits,plot_details) = params
    
    if os.name != 'nt':
        os.nice(3)
    


    enable_output=   plot_details['enable_output']
    geometry=   plot_details['geometry']
    mat_size=   plot_details['mat_size']
    output_path=   plot_details['output_path']
    dpi=   plot_details['dpi']
    img_size=   plot_details['img_size']


    nx,ny = mat_size
    xmin, xmax, ymin, ymax = geometry

    dxm=plot_details['dx']/plot_details['norm']
    dym=plot_details['dy']/plot_details['norm']
    shot = plot_details['shot']
    output_type = plot_details['output_type']
    tmp_folder  = plot_details['tmp_folder']

    dmin,dmax = limits[0]

    T = time.time()


    Ndets = size(data,1)
    fig = Figure((1.33*img_size,img_size))
    FigureCanvas(fig)
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_locator(MaxNLocator(10))

    errorbar1 = ax.errorbar(0,0,0,fmt='+b',label='input data')
    errorbar2 = ax.errorbar(0,0,0,fmt='.r',label='removed data')

    ax.axhline(0, c='k')

    for di in plot_details['dets_index'][:-1]:
        ax.axvline(x=0.5+amax(di), ls='--')


    #upper axis with detectors names 
    if plot_details['detectors_dict'] is not None:
        xlabels = [median(ind) for ind in plot_details['dets_index']]
        labels = list(plot_details['detectors_dict'].keys())
        ax2 = ax.twiny()
        ax2.set_xticks(xlabels)
        ax2.set_xticklabels(labels)
        ax2.set_xlim(0,Ndets)
        ax2.tick_params(axis='x', which='major', pad=0) 
        
        
    retrofit,= ax.plot([],[],'k-+', label='retrofit',zorder=99)
    leg = ax.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    from matplotlib.ticker import AutoMinorLocator
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    if linthresh[0]!= None:
        ax.set_yscale('symlog', linthreshy=linthresh[0]*2)

    ax.set_xlabel('detector')
    ylabel = ax.set_ylabel('')
    ax.set_xlim(0-0.5, Ndets-0.5)
    ax.set_ylim(dmin, dmax)

    for i,ts in enumerate(ntvec):
        if plot_details['plot_autoscale']:
            dmin = min(0,dmin)
            dmax = fc[i].max()*1.1
        
        if dmax > 2e5:
            fact,pre = 1e6, 'M'
        elif dmax > 2e2:
            fact,pre = 1e3, 'k'
        else:
            fact,pre = 1e0, ''
        if linthresh[0]!= None:
            fact,pre = 1e0, ''

        ylabel.set_text('Brightness [%sW/m$^2$] '%pre)
        ax.set_ylim(dmin/fact,dmax/fact)

        ind = isfinite(error[i,dets])
        x = dets[ind]
        y = data[i,x]/fact
        yerr = error[i,x]/fact
        
        update_errorbar(errorbar1, x,y,yerr)

        #do the same for wrong detectors
        wrong_dets = ~in1d(arange(Ndets), dets)
        wrong_dets|= ~isfinite(error[i,:])
        x = where(wrong_dets)[0]
        y = data[i,wrong_dets]/fact
        yerr = error[i,wrong_dets]/fact
        
        update_errorbar(errorbar2, x,y,yerr)

        retrofit.set_data(arange(Ndets),fc[i,:]/fact)
        
        if tsteps > 1 and False:
            time_line2.set_xdata([tvec[ts],tvec[ts]])
            time_line1.set_xdata([tvec[ts],tvec[ts]])
        name = 'brightness_%d_%.4d'%(shot,ts)+base

        if tsteps==1:
            name = 'previewA'+str(shot)+base 

        filename = tmp_folder + '/' + name + '.png'
        s = img_size/6
        fig.savefig(filename,transparent=True, dpi=dpi,bbox_inches=Bbox([[.2*s,.1*s],[7.5*s,5.7*s]]))
 

        if enable_output:
            fig.savefig(output_path+'/brightness_%.4d_%d.%s'%(ts,shot,output_type))
            

    try:
        sys.stdout.write("\r %2.1f %% fps:%2.1f" %((ntvec[-1])*100./tsteps, size(ntvec)/(time.time()-T)))
        sys.stdout.flush()  # maybe problem with paralelization ??
    except:
        pass

   
        

 
 

    
    
def matplotlib_image(params):    
    """
    Make preview of recontruction using the Matplotlib library, it is slower than GNUPLOT but it has higher quality
    """   
    
    try:
        #clean catch in foked processes
        matplotlib.font_manager._get_font.cache_clear()
    except:
        pass

    (ntvec,gprof,padding, data,error,fc,dets,tvec,my_cmap,
                        mag_field,chi2,lam,tsteps,rho_label,base,linthresh,limits,plot_details) = params
    
    if os.name != 'nt':
        os.nice(3)
    
    try:
        import matplotlib._cntr as cntr
    except:   #slower option from new matplolib      
        import matplotlib._contour as _contour

    
    enable_output=   plot_details['enable_output']
    geometry=   plot_details['geometry']
    mat_size=   plot_details['mat_size']
    output_path=   plot_details['output_path']
    plot_contours=   plot_details['plot_contours']
    n_contours=   plot_details['n_contours']
    dpi=   plot_details['dpi']
    img_size=   plot_details['img_size']
        



    nx,ny = mat_size
    xmin, xmax, ymin, ymax = geometry
    dxm=plot_details['dx']/plot_details['norm']
    dym=plot_details['dy']/plot_details['norm']
    shot = plot_details['shot']
    output_type = plot_details['output_type']
    tmp_folder  = plot_details['tmp_folder']

    gmin,gmax = limits[1]

    T = time.time()


    Ndets = size(data,1)

    fig = Figure((8*img_size/6.,img_size))
    FigureCanvas(fig)

    if plot_details['plot_chords']:  # !! allow plotting of chords !!! 
        from geom_mat_setting import loadgeometry
        Xchord, Ychord, distance, nl,virt_chord  = loadgeometry(plot_details['geometry_path'],
                                                list(plot_details['detectors_dict'].keys()), 1)   
        Xchord = Xchord[:,dets]
        Ychord = Ychord[:,dets]



    area = (xmin, xmax, ymin,ymax)

    area_tight = (plot_details['xmin']/plot_details['norm'],plot_details['xmax']/plot_details['norm']\
        ,(plot_details['ymin'])/plot_details['norm'],(plot_details['ymax'])/plot_details['norm'])
 
    sym_colorbar = False
    white_bg = True

    if plot_details['rem_back'] or (plot_details['rem_fsa'] and plot_details['transform_index']==1) or plot_details['plot_svd']:
        sym_colorbar = True
        white_bg = True
        

    fact,unit = 1,''
    
    if linthresh[1] != None:
        norm = symlogNorm(linthresh[1]/8.5,vmin=gmin,vmax=gmax)
        unit = 'W'
        
        ticks = [10**ceil(log10(gmax))]
        for i in range(3):
            ticks = [float('%.1g' % (ticks[0]/10))]+ticks
        if sym_colorbar:
            ticks = ticks[::-1]+[0,]+ticks
        else:
            ticks = [0,]+ticks

    else:
        if gmax > 1e6:
            fact, unit = 1e-6,'MW'
        elif gmax > 1e3:
            fact, unit = 1e-3,'kW'
        else: fact, unit = 1,'W'
            
        gmin,gmax = gmin*fact,gmax*fact 
        norm = Normalize(vmin=gmin,vmax=gmax)
        ticks = None


    if plot_details['blacken_negative'] and not sym_colorbar:
        extend = 'min'
    else: extend = 'neither'


    if enable_output:
        plt_list = { 'color': [my_cmap, 'gray', ''] }#,  'bw': ['Greys', 'black', 'bw'],
    else:
        plt_list = { 'color': [my_cmap, 'gray', '']  }
        

    for typ,cmap_typ in plt_list.items():
  
        fig.set_size_inches(8,6)
        ax = fig.add_subplot(111)
        ax.xaxis.set_major_locator( MaxNLocator(5) )
        ax.yaxis.set_major_locator( MaxNLocator(5) )


        ax.axis(geometry)
        if plot_contours:
            gmin = 0
            levels = norm.inverse(linspace(norm(gmin),norm(gmax),n_contours))
            prof_img = ax.contourf(gprof[:,:,0]*fact,norm=norm,#vmin=gmin,vmax=gmax,
                                   levels=levels,extent=area_tight,cmap=cmap_typ[0],
                                   extend=extend)

            CS2 = ax.contour(prof_img, levels=prof_img.levels,colors = 'k',linewidths=.2,  origin='lower')#, hold='on')
        else:
            
            prof_tmp = concatenate([padding, gprof[:,:,0]*fact, padding], axis=1)
            prof_img= ax.imshow(prof_tmp,cmap= cmap_typ[0],extent=area,norm=norm,#vmin=gmin,vmax=gmax,
                interpolation='nearest',aspect = str((ymax-ymin)/(xmax-xmin)), origin='lower')
            prof_img.set_clim((gmin,gmax)) 


        if linthresh[1] != None:
            c = fig.colorbar(prof_img,format=LogFormatterTeXExponent('%.2e'),ticks=ticks,extend=extend)
        else:
            c = fig.colorbar(prof_img,format=LogFormatterTeXExponent('%.2e'),extend=extend)
            tick_locator = MaxNLocator(nbins=7)
            c.locator = tick_locator
            c.update_ticks()
        
        #BUG hardcoded :(
        c.ax.set_position((0.76, 0.10, 0.03, 0.8))
        c.set_label('Emissivity [%sm$^{-3}$]'%unit, labelpad=4)

        if 'mag_field' in plot_details: #BUG is not working if False
            mag_surfs=ax.plot(mag_field[0,:,:,0],mag_field[1,:,:,0],'--',c='0.75',lw=.5)
            
        if 'mag_out' in plot_details:
            magr_out,magz_out = plot_details['mag_out']
            mag_surfs_out = [ax.plot(r,z,':',c='0.75',lw=.5)[0] for r,z in zip(magr_out[0],magz_out[0])]

        if 'theta_star' in plot_details:
            theta_star=ax.plot(plot_details['theta_star'][0][0],plot_details['theta_star'][0][0],'--',c='0.75',lw=.5)


        if 'boundary' in plot_details:
            bound, = ax.plot([],[],cmap_typ[1], lw=2)

        if 'ICRH_position' in  plot_details:
            if 'R2D' in plot_details["ICRH_position"]:
                v_line, = ax.plot([],[],ls='--',c='0.5',lw=1.5)
            else:
                Nlines = size(plot_details["ICRH_position"]['R'],1)
                v_line = [ax.axvline(x=0,ls='--',c='0.5',lw=1.5) for i in range(Nlines)]


        ax.set_xlabel('R [m]')
        ax.set_ylabel('z [m]', labelpad=-5)
        tlt = ax.set_title('')
        if plot_details['plot_chords']:
            ax.plot(Xchord/plot_details['norm'], Ychord/plot_details['norm'], 'gray', lw=.25)

        if plot_details['struct_dict'] is not None:
            for _,(xstruct,ystruct) in plot_details['struct_dict'].items():
                if len(xstruct) > 1:
                    color = '.5' if typ =='bw' or white_bg else 'w'
                    ax.plot(xstruct,ystruct, color,lw=0.5)
                    

        if tsteps> 1:
            dt = abs(mean(diff(tvec)))
            n_digit = max(0, -int(floor(log10(dt))))
        else:
            n_digit = 3
        
        for i,ts in enumerate(ntvec):
            
            
            if plot_details['plot_autoscale']:
                vmax = amax(gprof[:,:,i])
                vmin = min(0,amin(gprof[:,:,i] ))
                gmax = vmax*fact
                gmin = vmin*fact
  
            if sym_colorbar:
                gmin = min(gmin,-gmax)
                gmax = max(-gmin,gmax)
                
            if plot_details['blacken_negative']:
                gmin = -gmax/100
            if plot_contours:

                
                prof = gprof[:,:,i]*fact
                prof[prof>gmax*(1-1e-5)] = gmax*(1-1e-5)
                if not plot_details['blacken_negative'] and not sym_colorbar:
                    prof[prof<gmin] = gmin
  
                gmin -= gmax*1e-10
                
                if linthresh[1] != None:
                    norm = symlogNorm(linthresh[1]/8.5,vmin=gmin,vmax=gmax)
                else:
                    norm = Normalize(vmin=gmin,vmax=gmax)

                for coll in prof_img.collections+CS2.collections:
                    try:    ax.collections.remove(coll)
                    except ValueError: pass#Everything is not removed for some reason!    
                
                centerng = array((0,0, 0, -dym))

                dlev = (norm(gmax)-norm(gmin))/n_contours
                levels = norm.inverse(linspace(norm(gmin),norm(gmax),n_contours))
 

                prof_img = ax.contourf( prof,extent=array(area_tight),
                            norm=norm,vmin=gmin,vmax=gmax,levels=levels,origin='lower',extend=extend,cmap=cmap_typ[0])
                            #aspect = str((ymax-ymin)/(xmax-xmin)))
                
                CS2 = ax.contour(prof_img, levels=prof_img.levels[1:],colors = 'k',linewidths=.2,
                                        origin='lower',norm=norm,vmin=gmin,vmax=gmax)

                c._boundaries = asarray(levels)
                prof_img.set_clim(gmin,gmax)
                c.vmin = gmin
                c.vmax = gmax
                c.update_ticks()

            else:
                prof_tmp = concatenate([padding, gprof[:,:,i]*fact, padding], axis=1)
                prof_img.set_data(prof_tmp)
                if plot_details['plot_autoscale']:
                    prof_img.set_clim((gmin,gmax))

            title=('#%d, t=%3.'+str(n_digit)+'f, $\chi^2$= %4.2f $\lambda$=%.2f')%(shot,tvec[ts],chi2[i],lam[i])

            tlt.set_text(title)
            

            if 'mag_field' in plot_details:
                for j,surf_plt in enumerate(mag_surfs):    
                    surf_plt.set_data(mag_field[0,:,j,i],mag_field[1,:,j,i])
                    
            if 'mag_out' in plot_details:
                for j,surf_plt in enumerate(mag_surfs_out):
                    surf_plt.set_data(magr_out[ts][j], magz_out[ts][j])

            if 'boundary' in plot_details:
                bound.set_data(plot_details['boundary'][ts][:,0],plot_details['boundary'][ts][:,1])
            
            
            if 'theta_star' in plot_details:
                for j,plt_theta in enumerate(theta_star): 
                    plt_theta.set_data(plot_details['theta_star'][0][ts,:,j],
                                    plot_details['theta_star'][1][ts,:,j])
        
            if 'ICRH_position' in  plot_details:
                ind = argmin(abs(tvec[ts]-plot_details["ICRH_position"]['tvec']))
                ricrh = plot_details["ICRH_position"]['R'][ind]
                if 'R2D' in plot_details["ICRH_position"] and any(isfinite(ricrh)):
                    #take in account plasma diamagnetics
                    ind = argmin(abs(tvec[ts]-plot_details["ICRH_position"]['timeR2D']))
                    rhoR2D = plot_details["ICRH_position"]['rhoR2D'][ind][::-1]
                    R2D = plot_details["ICRH_position"]['R2D'][:,ind][::-1]
                    rho = linspace(0,1,mag_field.shape[2])
                    B_Bc = interp(rho, rhoR2D,R2D)/mag_field[0,:,:,i]
                            
                    try: 
                        countour = cntr.Cntr(mag_field[0,:,:,i],mag_field[1,:,:,i],B_Bc )
                        nlist = countour.trace(level0=1,level1=1,nchunk=0)

                    except: #slower option from new matplolib 
                        gen = _contour.QuadContourGenerator(mag_field[0,:,:,i],mag_field[1,:,:,i],B_Bc, False, 0)
                        nlist = gen.create_contour(1)

                    lines = nlist[:len(nlist)//2]
                    lines = vstack(lines)
                    lines = lines[argsort(lines[:,1])]
                    v_line.set_data( lines[:,0], lines[:,1])
                elif  any(isfinite(ricrh)):
                    for r,v in zip(ricrh,v_line): v.set_xdata([r,r])


            name = 'emissivity_%d_%.4d'%(shot,ts)+base
                        
            if tsteps==1: name = 'previewB'+str(shot)+base 
            filename = tmp_folder+'/' + name + '.png'
            ax.set_aspect('equal', 'datalim')

            
            if typ != 'bw':   
                s = img_size/6.
                fig.savefig(filename,transparent=True, dpi=dpi,bbox_inches=Bbox([[.3*s,.1*s],[7.1*s,5.74*s]]))

            ax.axis(geometry)

            if enable_output:
                ax_bbox = ax.get_position()
                ax.axis(area_tight)
                ax.set_aspect('equal', 'box')
                #BUG bbox_inches cause jumping of the plot axis in time, it is bad for movies 
                fig.savefig(output_path+'/emissivity_%.4d_%d_%s%s.'%(ts,shot,cmap_typ[2],base)
                            +output_type,bbox_inches='tight' )
                 
                ax.axis(geometry)
                ax.set_position(ax_bbox)

    try:
        sys.stdout.write("\r %2.1f %% fps:%2.1f" %((ntvec[-1])*100./tsteps/2+50, size(ntvec)/(time.time()-T)))
        sys.stdout.flush()  # maybe problem with paralelization ??
    except:
        pass



    
def make_svd(inputs, tokamak,progress,results ):
    """
    Create subplots from first N modes of SVD and its time evolution
    """
  
        
    print('PLOTTING ..SVD')
    tmin = inputs['tmin']
    tmax = inputs['tmax']

    tvec = results['tvec']
    n_svd_show = inputs['n_svd_show']
    

    tsteps = len(tvec)
    if tsteps < n_svd_show:
        return 

    geometry = tokamak.area_axis
    boundary = tokamak.get_boundary(100,time=mean(tvec))

    import matplotlib.gridspec as gridspec

    print('SVD ..')
    ind = (tsteps//10000+1)
    G = results['g']
    G_tmp = G[::ind]
    tvec_tmp = tvec[::ind]


    ur,sr,vr = fsvd(G_tmp, n_svd_show,i=2)
    v_sign = sign(mean(vr[:,0]))
    ur[:,0]*= v_sign
    vr[:,0]*= v_sign
    vr = vr.T

    savez_compressed(inputs['tmp_folder']+'/svd',u=ur,s=sr,v=vr)


    print("SVD done")
    # Singular Value Decomposition of the 2D emissivity evolution

    tsteps = len(tvec_tmp)

    rhop,magx, magy  = tokamak.mag_equilibrium(tvec,n_rho=11, return_mean=True)
    if len(rhop)!= 11:
        skip_mag = max(size(magx, 1)/10, 1)    #use maximaly 10 lines in the graph
        ind = slice(skip_mag-1,None,skip_mag)
        rhop,magx, magy = rhop[ind],magx[:,ind], magy[:,ind]
    
    # plotting    
    fig = Figure( (1.3*inputs['img_size'],inputs['img_size']))
    FigureCanvas(fig)
    gs = gridspec.GridSpec(2, n_svd_show, width_ratios=[1,]*n_svd_show,height_ratios=[2,1])
    gs.update(left=0.1, right=0.90, hspace=0.01)
    fig.subplots_adjust(hspace=0.0, wspace=0.05)

    xmin,xmax,ymin,ymax = geometry
    
    area_tight = (tokamak.xmin/tokamak.norm,tokamak.xmax/tokamak.norm\
        ,(tokamak.ymin)/tokamak.norm,(tokamak.ymax)/tokamak.norm)
    

    for i in range(n_svd_show):
        ax = fig.add_subplot(gs[i])
        ax.plot(magx, magy, '0.4',lw=.2)
        vmax = max(-vr[i,:].min(),vr[i,:].max())
        vmin = min(-vr[i,:].max(),vr[i,:].min())

        ax.plot(boundary[:,0]/tokamak.norm, boundary[:,1]/tokamak.norm, '0.6', linewidth=1)
        ax.imshow(reshape(vr[i,:],(tokamak.ny,tokamak.nx),order='F'),cmap=inputs['cmap_neg'],extent=geometry
                ,origin='lower',vmin=vmin,vmax=vmax)
        
        if i: ax.yaxis.set_major_formatter( NullFormatter() )
        ax.xaxis.set_major_locator( MaxNLocator(3) )
        
        if i == 0:  ax.set_ylabel('topos',labelpad=-5)
        for xlabel_i in ax.axes.get_xticklabels():
            xlabel_i.set_fontsize(10)
        for ylabel_i in ax.axes.get_yticklabels():
            ylabel_i.set_fontsize(10)
   
        ax.axis(area_tight)

    M_norm = sp.linalg.norm((G_tmp.flatten()))

    ur = dot(vr, G.T).T
    from scipy.signal import welch

    ax = None
    ymax = 0
    ymin = infty
    for i in range(n_svd_show):
        ax =  fig.add_subplot(gs[n_svd_show+i],sharey=ax)

        ax.set_title('Energy=%.3f%%'%(sr[i]**2/M_norm**2*100),fontsize=7)
  
        if inputs['svd_freq_domain']:
            f_s = len(tvec)/(tvec[-1]-tvec[0])
            fvec, Pxx = welch(ur[:,i],f_s, nperseg = min(len(tvec), 32), scaling='spectrum')
       
            ax.loglog(fvec[1:], Pxx[1:]  )
            ymin = min(ymin,Pxx.min())
            ymax = max(ymax,Pxx.max())


            if i == 0:
                ax.set_ylabel('FT(chronos)')
            else:
                ax.yaxis.set_major_formatter( NullFormatter())
                
            ax.set_ylim(ymin,ymax)
            ax.set_xlim(fvec[1],fvec[-1])
            ax.grid(True)

            ax.set_xlabel('f [Hz]',fontsize=9, labelpad=2)

        else:
            if i == 0:
                ax.set_ylabel('chronos')
            
            ax.plot(tvec,ur[:,i]*sr[i])       # first order chronos in relative signal energy units 
            ax.axhline(y=0,color='k')
            ax.set_xlim(tvec.min(), tvec.max())
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(5))
            
        for xlabel_i in ax.axes.get_xticklabels():
            xlabel_i.set_fontsize(10)
        for ylabel_i in ax.axes.get_yticklabels():
            ylabel_i.set_fontsize(10)

    name = inputs['tmp_folder']+'/SVD'+'_'+str(inputs['shot'])+'.png'
    fig.savefig(name, transparent=True, dpi=inputs['dpi'], bbox_inches='tight')

    if inputs['enable_output']:
        fig.savefig(inputs['output_path']+'/SVD'+'_'+str(inputs['shot'])+'.'+inputs['output_type'],
                    bbox_inches='tight',transparent=False, dpi=150)


    print("SVD end")

def plot_adv(tvec, data, name, save_data = True,**kwarg):
    
    global inputs, plot_details

    img_size = plot_details['img_size']
    dpi = plot_details['dpi']
    output_path = inputs['output_path']
    #local_path = inputs['local_path'] 
    enable_output =  inputs['enable_output']
    output_type=  inputs['output_type']
    
    sub = 0;  subplots = 0 


    fig = Figure((1.3*img_size,img_size))
    FigureCanvas(fig)
    fig.subplots_adjust(hspace=0.1, wspace = 0.07)


    for d in data: 
        if 'subplt' in d:  subplots += 1
        
    axes = []
    ax = None
    for d in data:
        if 'subplt' in d:
            sub += 1
            ax = fig.add_subplot(subplots, 1, sub)
            d.pop('subplt')
            axes.append(ax)
        elif ax is None:
            ax = fig.add_subplot(1, 1, 1)
            #ax = gca()
        fmt = d.pop('fmt') if 'fmt' in d else '' 
        xvec = d.pop('tvec') if 'tvec' in d else tvec
        if not 'yerr' in d:
            ax.plot(xvec, d.pop('data'),fmt, label=d.pop('label'))
        else:
            ax.errorbar(xvec,d.pop('data'),label=d.pop('label'),yerr=d.pop('yerr'))
            
        ax.set_xlim(amin(xvec), amax(xvec))
        ax.yaxis.set_major_locator(MaxNLocator(6))
        if xvec is tvec and not tokamak.impur_inject_t is None:
            for lbo in tokamak.impur_inject_t:
                if lbo>tvec.min() and lbo<tvec.max():
                    ax.axvline(x=lbo,color='k',lw=.5,ls='--')
            
        ax.set(**d)
    ax.set(**kwarg)
    
    if axes == []: axes.append(ax)
    for ax in axes:
        try:
            leg = ax.legend(loc='best', fancybox=True)
            for l in leg.legendHandles:  l.set_linewidth(5)
            leg.get_frame().set_alpha(0.8)
        except:
            pass

    if enable_output:
        fig.savefig(output_path+'/'+name+'.'+output_type, bbox_inches='tight',transparent=False, dpi=120)

    name = inputs['tmp_folder']+'/'+name+'.png'

    fig.savefig(name, transparent=True, dpi=inputs['dpi'], bbox_inches='tight')
    
    
    

def plot_2D_adv(yvec,xvec,data,name,plot_type=0,cmap=my_cmap_,norm=Normalize(),sym_colorbar=False,**kwarg):
    
    global inputs,tokamak

        
    img_size = plot_details['img_size']
    dpi = plot_details['dpi']
    output_path = inputs['output_path']
    tmp_folder = inputs['tmp_folder'] 
    enable_output =  inputs['enable_output']
    output_type =  inputs['output_type']
    inputs['blacken_negative']
    
    
    vmin,vmax = mquantiles(data[isfinite(data)],(0.01, 0.99))
    vmin =  min(0,vmin) 
        
    if sym_colorbar:
        vmax = max(-vmin,vmax)
        vmin = min(vmin,-vmax)
    elif inputs['blacken_negative']:
        vmin = (vmin-vmax)/1e3
        
    
        
 
    fig = Figure((1.3*img_size,img_size))
    FigureCanvas(fig)

    ax = fig.add_subplot(111)
    extent = (amin(yvec),amax(yvec), amin(xvec), amax(xvec))
    c = None
    if plot_type == 0:
        im = ax.imshow(data, interpolation='bilinear', origin='lower', cmap =cmap, extent=extent,norm=norm)

        if isinstance(norm, symlogNorm):
            ticks = norm.inverse(linspace(0,1,7,endpoint=True))
            ticks = [float('%.1g' % (t)) for t in ticks] 
        else:
            ticks = None

        c = fig.colorbar(im,format=LogFormatterTeXExponent('%.2e'),ticks=ticks)

  
        if ticks is None:
            c.locator = MaxNLocator(nbins=7)
            c.update_ticks()
 

        
        if inputs['rem_back'] or (inputs['rem_fsa'] and tokamak.transform_index==1):  #set symetrical color scale with 0 in the middle
            im.set_clim([min(vmin,-vmax), max(-vmin, vmax)])
        else:
            im.set_clim([vmin,vmax])
        levels = linspace(vmin,nanmax(data),30)

        ax.contour( data, levels=levels, linewidths=.1,  colors='k',extent=extent)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.axis(extent)
        
    elif plot_type == 1:
        im = ax.imshow( data, interpolation='bilinear', origin='lower', cmap ='binary', extent=extent,norm=norm)
  
        c = fig.colorbar(im,format=LogFormatterTeXExponent('%.2e'))

        ax.contour( data, 20, linewidths=1,  colors='k',extent=extent)
    elif plot_type == 2:
        ax.contour( data, 20, linewidths=0.5,  colors='k',extent=extent)

    if 'clabel' in kwarg:
        clabel = kwarg.pop('clabel')
        if c != None: c.set_label(clabel)

    
    ax.axis('tight')
    ax.set(**kwarg)
    
    if not tokamak.impur_inject_t is None:
        c = 'k' if inputs['rem_back'] else 'w'
        for lbo in tokamak.impur_inject_t:
            if ax.get_xlabel().find('Time')!= -1:
                ax.axvline(lbo,color=c,lw=.5)
            if ax.get_ylabel().find('Time')!= -1:
                ax.axhline(lbo,color=c,lw=.5)

                 


    if enable_output:
        fig.savefig(output_path+name+'.'+output_type, bbox_inches='tight',transparent=False, dpi=120)

    fig.savefig(tmp_folder+'/'+name+'.png', bbox_inches='tight',transparent=True, dpi=dpi)




def postprocessing_plot(input_data):
    """
    Determine and plot evolution of Chi2, output power, center position, 
    emissivity profile. 
    """
              


    global inputs, tokamak, plot_details, my_cmap
    

    inputs, tokamak, progress,results = input_data

    print('PLOTTING postprocessing ..')
    tmp_folder = inputs['tmp_folder']
    shot = str(inputs['shot'])

    tvec = results['tvec'] 
    t = time.time()
    tsteps = len(tvec)
    t_label = 'Time ['+tokamak.t_name+']'
    

    G = results['g'].T
    power = results['power']
    profile = results['fsa_emiss']


    m_negative = sum(G<0)
    if  m_negative > 0:   print('negative values in the G  %.3f%%'%(m_negative*100./size(G)))
      
    resid = results['retro']-results['data'] 
    resid/= results['error']
    
    #
    xcpix=linspace(tokamak.xmin+tokamak.dx/2,tokamak.xmax+tokamak.dx/2, tokamak.nx)/tokamak.norm
    ycpix=linspace(tokamak.ymin+tokamak.dy/2,tokamak.ymax+tokamak.dy/2, tokamak.ny)/tokamak.norm   #raise


    gres = G.reshape(tokamak.ny, tokamak.nx, tsteps, order='F')


    # LAMBDA
    ylim = (0,1) if all((results['lam0']>=0)&(results['lam0']<=1)) else  (None,None)
     
    plot_adv(tvec, [{'data':results['lam0'],'label':'$\lambda$','ylim':ylim}],
                'lam'+'_'+shot,
                title="Regularization level",
                xlabel=t_label) 
                #,ylabel='log $\lambda $')
 

    #CHI2  # chi2
    
    try:
        plot_adv(tvec, [{'tvec':tvec, 'data':results['chi2_real'], 
                            'ylabel':"$\chi^2/doF$ residuum",
                            'xlabel':t_label,
                            'subplt':True,
                            'yscale':'log' , 
                            'label':'Real'},
        {'tvec':tvec, 'data':results['chi2'], 
                            'ylabel':"$\chi ^2/doF$",
                            'xlabel':t_label,
                            'yscale':'log',
                            'label':'Guess'},
        {'tvec':results['dets'], 'data':mean(resid,0), 'yerr': std(resid,0) ,
                            'label':'mean residuum for each detector' ,
                            'ylabel':"Normalised residuum",\
                            'xlabel':'Detector' ,
                            'ylim': (-1.5*amax(abs(mean(resid,0))),1.5*amax(abs(mean(resid,0)))), 
                            'subplt':True} ], 
            'chi2'+'_'+shot)

    except:
        print(results['chi2_real'])
        


        

    if inputs['post_proces_equi']:
        
        try:
            import matplotlib._cntr as cntr
        except: #slower option        
            import matplotlib._contour as _contour
            
        magx, magy = results['magx'].T,results['magy'].T

        R,Z = meshgrid(tokamak.xgrid+tokamak.dx/2,tokamak.ygrid+tokamak.dy/2)

        shaf_shift_mag = zeros(tsteps)
        shaf_shift_sxr = zeros(tsteps)
        elongation_sxr = zeros(tsteps)
        elongation_mag = zeros(tsteps)
        elong_all,Rc_all,Zc_all = [],[],[]
        lim = [.15**2,.8]
        lim = [0.1,1]
        fun = lambda x:(x)

        n_mag = magx.shape[1]
        
     
        
        for it in trange(tsteps,desc='Computing emiss. contours: '): 


            r_a = hypot(magx[:,:,it]-magx[:,(0,),it],magy[:,:,it]-magy[:,(0,),it]).mean(0)
            r_a/= r_a[-1]
            r_a_eq = linspace(1./n_mag ,1-1./n_mag,n_mag)

            levels = interp(fun(r_a_eq),r_a, profile[it])

            #contour searching routine from matplotlib     
            try: 
                c = cntr.Cntr(R,Z, gres[:,:,it])
            except: #slower option  
                gen = _contour.QuadContourGenerator(R,Z, gres[:,:,it],bool_(Z*0), False, 0)
 

            rho_contours = []
            for i,lev in enumerate(levels):
                try:
                    nlist = c.trace(level0=lev, level1=lev, nchunk=0)
                    nlist = nlist[:len(nlist)/2]
                except: #slower option  
                    nlist = gen.create_contour(lev)

                lines = nlist[:len(nlist)//2]
                if len(lines) == 0:
                    lines = [empty((0,2)),]
                line = []
                #choose the longest line
                for l in lines:
                    if len(l)>=len(line):
                        line = l
                rho_contours.append(line)    
                
            from shared_modules import   fitEllipse
   
            elong = ones(len(rho_contours))*nan
            Rc = ones(len(rho_contours))*nan
            Zc = ones(len(rho_contours))*nan
            for ic, C in enumerate(rho_contours):
                try:
                    x,y = C.T
                    elips = fitEllipse(x,y)
                    Rc[ic],Zc[ic] = elips.ellipse_center()
                    Rc[ic]*= 0.999  #correction estimated from the comparism with the phantoms
                    a,b = elips.ellipse_axis_length()
                    elong[ic] = b/a*0.99  #correction estimated from the comparism with the phantoms
                except:
                    continue

            a,b = polyfit(fun(r_a_eq), Rc,1)
            shaf_shift_sxr[it] = -a

            ind = (fun(r_a)<lim[1])&(fun(r_a)> lim[0])
            a,b = polyfit(r_a[ind]**2,magx[:,ind,it].mean(0),1)
            shaf_shift_mag[it] = -a

            
            elongation_sxr[it] = nanmedian(elong)
            elongation_mag[it] = nanmedian(magy[:,ind,it].std(0)/magx[:,ind,it].std(0))
            
            
            elong_all.append(elong)
            Rc_all.append(Rc)
            Zc_all.append(Zc)

                


        elong_all = vstack(elong_all)
        Rc_all = vstack(Rc_all)
        Zc_all = vstack(Zc_all)

        Rm = ones(n_mag)*nan
        Zm = ones(n_mag)*nan
        elongm = ones(n_mag)*nan
        for ir in range(1,n_mag):
            try:
                a = fitEllipse(double(magx[:,ir].mean(-1)),double(magy[:,ir].mean(-1)))
                Rm[ir],Zm[ir] = ellipse_center(a)
                Rm[ir]*= 0.999  #correction estimated from the comparism with the phantoms
                a,b = ellipse_axis_length(a)
                elongm[ir] = b/a*0.99  #correction estimated from the comparism with the phantoms
            except:
               pass

        fig = Figure((1.6*plot_details['img_size'],)*2)
        FigureCanvas(fig)
        
        fig.subplots_adjust(hspace=0.1, wspace = 0.07)
        fig.suptitle('SXR/equilibrium comparism, Discharge: %s, t = %.3s'%(shot, mean(tvec)))

        ax = fig.add_subplot(221)
        ax.plot( fun(r_a_eq), nanmean(elong_all,axis=0))
        ax.fill_between(fun(r_a_eq), nanmean(elong_all,axis=0) - nanstd(elong_all,axis=0),
                        nanmean(elong_all,axis=0) + nanstd(elong_all,axis=0)  ,alpha=.2, facecolor='b', edgecolor='None')
        ax.set_ylabel('Elongation')
        ax.text(0.1, 0.9, 'SXR', transform=ax.transAxes,color='b')
        ax.text(0.1, 0.83, 'equilibrium', transform=ax.transAxes,color='g')
        ax.plot( r_a ,elongm )

        ax.set_ylim(1,2)

        ax = fig.add_subplot(222)

        ax.plot( fun(r_a_eq), nanmean(Rc_all,axis=0))
        ax.fill_between(fun(r_a_eq), nanmean(Rc_all,axis=0) - nanstd(Rc_all,axis=0),
                        nanmean(Rc_all,axis=0) + nanstd(Rc_all,axis=0),alpha=.2,facecolor='b', edgecolor='None')
        ax.text(0.7, 0.9, 'SXR', transform=ax.transAxes,color='b')
        ax.text(0.7, 0.83, 'equilibrium', transform=ax.transAxes,color='g')

        ax.plot( r_a , Rm)

        ax.set_ylabel('R shift [m]')
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")

        ax = fig.add_subplot(223)

        ax.plot( fun(r_a_eq), nanmean(Zc_all,axis=0))
        ax.fill_between(fun(r_a_eq), nanmean(Zc_all,axis=0) - nanstd(Zc_all,axis=0),
                        nanmean(Zc_all,axis=0) + nanstd(Zc_all,axis=0),alpha=.2, 
                        facecolor='b', edgecolor='None')
        ax.text(0.1, 0.2, 'SXR', transform=ax.transAxes,color='b')
        ax.text(0.1, 0.13, 'equilibrium', transform=ax.transAxes,color='g')
            
        print('je ot blbe!!!')

     
        ax.plot( r_a , Zm)

        ax.set_ylabel('Z shift [m]')
        ax.set_xlabel('r/a')

        ax = fig.add_subplot(224)
        ax.set_aspect('equal', 'datalim')
        ax.set_ylabel('Z [m]')
        ax.set_xlabel('R [m]')
        ax.text(0.7, 0.07, 'SXR', transform=ax.transAxes,color='b')
        ax.text(0.7, 0.13, 'equilibrium', transform=ax.transAxes,color='g')
        from scipy.interpolate import interp1d

        magx_red = interp1d(r_a, magx[:,:,it],axis=1,assume_sorted=True)(r_a_eq[4::5])
        magy_red = interp1d(r_a, magy[:,:,it],axis=1,assume_sorted=True)(r_a_eq[4::5])

        ax.plot(magx_red, magy_red,'g--',lw=.5)

        ax.plot(magx[:,-1,it], magy[:,-1,it],'k-.',lw=.5)
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        #BUG dělat contury z průměru?? 
        for A in rho_contours[4::5]: ax.plot(*A.T,color='b',lw=.5)

        fig.savefig(tmp_folder+'/emiss_contours.pdf')



        #Shafranov shift
        plot_adv(tvec, [{'data': shaf_shift_sxr*100,'label':'SXR'} ,{'data': shaf_shift_mag*100,
                'label':'magnetics'}],'shafr_shift'+'_'+shot,title="Shafranov shift (only w/o asymmetries)",
                xlabel=t_label,ylabel='$\Delta(0)$ [cm]',ylim=(0,median(shaf_shift_mag*100)*2))
        
        plot_adv(tvec, [{'data': elongation_sxr,'label':'SXR'} ,{'data': elongation_mag,
                'label':'magnetics'}],'elongation'+'_'+shot,title="Elongation (only w/o asymmetries)",
                xlabel=t_label,ylabel='$\kappa$ [-]',ylim=(1,1.7))
        
 
    cmap = my_cmap
    if inputs['rem_back'] and inputs['cmap_neg']!= 'none':  # remove background  
        cmap = cm.get_cmap(inputs['cmap_neg'])
        
        
        
    

    bcg_subst_fract = inputs['bcg_subst_fract']
    subst_ind = slice(None,int(ceil(tsteps*bcg_subst_fract))) 
    n_mag = profile.shape[1]
    profile0 = zeros(n_mag)
    if  inputs['rem_back']:
        power-= power[subst_ind].mean()
        profile0 = profile[subst_ind].mean(0)
        profile-= profile0

    if nanmax(power) > 1e5:
        fact,pre = 1e6, 'M'
    elif nanmax(power) > 1e2:
        fact,pre = 1e3, 'k'
    else:
        fact,pre = 1e0, ''
    
    # POWER
    #gres is in W/m^3
    #dV = 2*pi*xcpix*dx*dy  #volume of the single pixel
    power_list = []

    if hasattr(tokamak,'total_rad'):
        for label, tvec_rad, rad in tokamak.total_rad:
            fmt = 'o--' if len(tvec_rad) < 100 else '--'
            power_list.append({'tvec':tvec_rad, 'data':rad/fact, 'label':label,'fmt':fmt})

    
    
    if 'power_div' in results:
        power_div = results['power_div']
        power_list.append({'data':power_div/fact, 'label':'Divertor Power'}, )
        power_list.append({'data':(power-power_div)/fact, 'label':'Core Power'}, )
        
    power_list.append({'data':power/fact, 'label':'Tomography','fmt':'k'})
    
    if hasattr(tokamak, 'pDiode'):
        diode = tokamak.pDiode
        plot_adv(tvec, [{'tvec':diode[:,0], 'data':diode[:,1]*mquantiles(power[isfinite(power)], 0.95),\
            'label':'Photodiode'}, { 'data':power, 'label':'Camera'}],'power'+'_'+shot,\
                title="Radiated power",xlabel=t_label,ylabel='Power [%sW]'%pre, ylim=(0,None))
    else:
        plot_adv(tvec, power_list, 'power'+'_'+shot,  title="Radiated power", 
                 xlabel=t_label, ylabel='Power [%sW]'%pre,
                 ylim=(0,amax(power/fact)*1.2), xlim=(tvec[0], tvec[-1]))


    if nanmax(profile) > 1e5:
        fact,pre = 1e6, 'M'
    elif nanmax(profile) > 1e2:
        fact,pre = 1e3, 'k'
    else:
        fact,pre = 1e0, ''
    

    units = '[%sW/m$^3$]'%pre
    norm = Normalize()
    if inputs['plot_loglin']:
        vmax = amax(abs(profile))
        vmin = 0 if not inputs['rem_back'] and not  inputs['plot_svd'] and not inputs['rem_fsa'] else -vmax
        norm = symlogNorm(std(profile)/10,vmin=vmin,vmax=vmax)
        units = '[W/m$^3$]'
        fact = 1

        

    profile_name = 'profile_%s_%.3f-%.3f'%(shot,tvec[0],tvec[-1])
   
    plot_2D_adv(tvec,[0,1],profile.T/fact,profile_name, 
                title='Emissivity profile',
                cmap=cmap,norm=norm,xlabel=t_label,ylabel=tokamak.rho_label,
                clabel='Emissivity %s'%units,plot_type=0)
    
    
    copyfile(tmp_folder+'/'+profile_name+'.png',tmp_folder+'/profile_%s.png'%shot)

    
    emiss_x_dicts = [{'data':results['xmass'][:,2],'label':'Emissivity 2/3'},\
                     {'data':results['xmass'][:,1],'label':'Emissivity 1/3'},\
                     {'data':results['xmass'][:,0],'label':'Emissivity 0/3'}]
                
    emiss_y_dicts = [{'data':results['ymass'][:,2],'label':'Emissivity 2/3'},\
                     {'data':results['ymass'][:,1],'label':'Emissivity 1/3'},\
                     {'data':results['ymass'][:,0],'label':'Emissivity 0/3'}]
                         
              
    mag_x_dicts = [{'data':results['xmass_eq'][:,2],'label':'Magnetic 2/3'},\
                   {'data':results['xmass_eq'][:,1],'label':'Magnetic 1/3'},\
                   {'data':results['xmass_eq'][:,0],'label':'Magnetic 0/3'}]
                
    mag_y_dicts = [{'data':results['ymass_eq'][:,2],'label':'Magnetic 2/3'},\
                   {'data':results['ymass_eq'][:,1],'label':'Magnetic 1/3'},\
                   {'data':results['ymass_eq'][:,0],'label':'Magnetic 0/3'}]
                         
        
    if hasattr(tokamak, 'mag_axis'):
        if any(tokamak.mag_axis['tvec']<tvec.min()):
            ind_min = where(tokamak.mag_axis['tvec']<tvec.min())[0][-1]
        else:  ind_min = None
            
        
        if any(tokamak.mag_axis['tvec']>tvec.max()):
            ind_max = where(tokamak.mag_axis['tvec']>tvec.max())[0][0]+1 
        else: ind_max = None

        mag_ind=slice(ind_min,ind_max)
        mag_tvec=tokamak.mag_axis['tvec'][mag_ind]
        Zmag = tokamak.mag_axis['Zmag'][mag_ind]
        Rmag = tokamak.mag_axis['Rmag'][mag_ind]
    else:
        Rmag = results['magx'][:,0,:].mean(1)
        Zmag = results['magy'][:,0,:].mean(1)
        mag_tvec = tvec

    rmin1,rmax1 = mquantiles(results['xmass'], (.05, .95))
    rmin4,rmax4 = mquantiles(Rmag,   (.05, .95))
    rmin = min(rmin1,rmin4)-.01
    rmax = max(rmax1,rmax4)+.01

    zmin1,zmax1 = mquantiles(results['ymass'], (.05, .95))
    zmin4,zmax4 = mquantiles(Zmag,   (.05, .95))
    zmin = min(zmin1,zmin4)-.01
    zmax = max(zmax1,zmax4)+.01


    
    plot_adv(tvec,[{'tvec':mag_tvec,'data':Rmag+mean(inputs['magfield_shift_core'][0]),\
        'label':'Magnetics'}]+emiss_x_dicts+mag_x_dicts,'xmass'+'_'+shot,\
            ylim=(rmin,rmax),title="Center of mass in R coordinate",\
        xlabel=t_label,ylabel='R [m]',xlim=(tvec.min(),tvec.max()))
    plot_adv(tvec,[{'tvec':mag_tvec,'data':Zmag+mean(inputs['magfield_shift_core'][1]),\
        'label':'Magnetics'}]+emiss_y_dicts+mag_y_dicts,'ymass'+'_'+shot,\
            ylim=(zmin,zmax),title="Center of mass in Z coordinate",\
        xlabel=t_label,ylabel='Z [m]',xlim=(tvec.min(),tvec.max()))


    name = tmp_folder+'/plasma_parameters_'+shot+'.txt'
    tvec_digits = int(max(0,-ceil(log10(mean(diff(tvec))))))+1

    output = {}

    if inputs['post_proces_equi']:
        eq_out = c_[tvec,results['xmass'][:,1],results['ymass'][:,1],results['xmass_eq'][:,1],results['ymass_eq'][:,1],
                    elongation_sxr,elongation_mag,shaf_shift_sxr,shaf_shift_mag]
        savetxt(tmp_folder+'/equilibrium_'+shot+'.txt', eq_out,
                    fmt=[ '%1.'+str(tvec_digits)+'e',]+['%.5e']*8,
        header='time [s]\tR 1/3 [m]\tz 1/3 [m]\t Rmag 1/3 [m]\tzmag 1/3 [m]\tSXR_elong\tmag_elong\tshaf_shioft_sxr\t shaf_shift_mag')
        
        output['elongation_sxr'] = elongation_sxr
        output['elongation_mag'] = elongation_mag
        output['shaf_shift_sxr'] = shaf_shift_sxr
        output['shaf_shift_mag'] = shaf_shift_mag

    savetxt(name, c_[tvec,results['xmass'],results['ymass'],results['xmass_eq'],
                    results['ymass_eq'] ,results['position'],results['power'],results['power']-results['power_div'] ],
                    fmt=[ '%1.'+str(tvec_digits)+'e',]+['%.5e']*16,
            header='time [s]\tR 0/3 [m]\tR 1/3 [m]\tR 2/3 [m]\tz 0/3 [m]\tz 1/3 [m]\tz 2/3[m]'\
                    +'\t R_eq 0/3 [m]\tR_eq 1/3 [m]\tR_eq 2/3 [m]\tz_eq 0/3 [m]\t'\
                    +'z_eq 1/3 [m]\tz_eq 2/3[m]\t rmag [m]\tzmag\[m]\tpower_tot[W]\tpower_core[W]')


    savetxt(tmp_folder+'/emiss_profile0_%s.txt'%shot,profile0)
    
    savetxt(tmp_folder+'/emiss_profile_%s.txt'%shot, c_[tvec, profile],
            fmt=['%1.'+str(tvec_digits)+'e']+['%1.4e',]*n_mag)
    
    output['xmass'] = results['xmass']
    output['ymass'] = results['ymass']
    output['xmass_eq'] = results['xmass_eq']
    output['ymass_eq'] = results['ymass_eq']
    output['ymass_eq'] = results['ymass_eq']
    output['power'] = single(results['power'])
        
    if 'power_div' in results:
        output['power_div'] = results['power_div']
    output['position'] = single(results['position'])
    output['rho'] = results['mag_rho']
    output['tvec'] = tvec
    output['fsa_emiss0'] = profile0
    output['fsa_emiss'] = profile

    
    
    if 'fsa_emiss_err' in results:
        output['fsa_emiss_err'] = results['fsa_emiss_err']
        output['fsa_emiss_cov'] = results['fsa_emiss_covar'].mean(0)

    post_process_path = tmp_folder+'/postprocessing_%s_%.3f-%.3f.npz'%(shot, tvec[0], tvec[-1])
    savez_compressed(post_process_path,rho_lbl=tokamak.radial_coordinate,**output )


    if inputs['enable_output']:
        output_path = inputs['output_path']
        copyfile(post_process_path,output_path+'/postprocessing_'+shot+'.npz')
        copyfile(tmp_folder+'/'+'emiss_profile_'+shot+'.txt',output_path+'/emiss_profile_'+shot+'.txt')
        copyfile(tmp_folder+'/'+'emiss_profile0_'+shot+'.txt',output_path+'/emiss_profile0_'+shot+'.txt')
        copyfile(tmp_folder+'/'+'plasma_parameters_'+shot+'.txt',output_path+'/plasma_parameters_'+shot+'.txt')
        if inputs['post_proces_equi']:
            copyfile(tmp_folder+'/'+'equilibrium_'+shot+'.txt',output_path+'/equilibrium_'+shot+'.txt')

 
    print('Postprocessing time %.1fs' % (time.time()-t)) 

  




def CalcPoloidalModeSpectrum(input_data):
    
    #-m number spectrogram


    inputs, tokamak,progress,results = input_data
    
    global plot_details 

    
    img_size = plot_details['img_size']
    dpi = plot_details['dpi']

    tvec = results['tvec']
    nx = inputs['nx']
    ny = inputs['ny']
    tsteps = len(tvec)
    tmp_path = inputs['tmp_folder']
    

    cmap_neg = inputs['cmap_neg']

    gres = results['g'].T.reshape(ny, nx, tsteps, order='F')

    n_theta = 64
    rhop,magr, magz = tokamak.mag_equilibrium(tvec, return_mean=True,n_theta=n_theta, radial_coordinate='rho_pol')
    rhop,magr, magz = rhop[:-1], magr[:,:-1], magz[:,:-1]

    theta   = tokamak.mag_theta_star(tvec.mean(),rhop,magr,magz,rz_grid=False)

    n_mag = size(magr,1)
    n_modes = 7

    scaling = array([tokamak.dx,tokamak.dy])
    offset = array([tokamak.xmin,tokamak.ymin])+scaling/2
    
    from scipy.ndimage.interpolation import map_coordinates
    
    gres_mean = gres.mean(-1)
    theta0 = linspace(0,2*pi, theta.shape[1], endpoint=False)
    coords = c_[magr.ravel(),magz.ravel()].T
    idx = (coords-offset[:,None])/ scaling[:,None]

    complex_profile = zeros((n_mag, n_modes-1,tsteps),dtype=complex64)

    
    map_prof0 = map_coordinates(gres_mean.T,idx,order=2)
    map_prof0 = map_prof0.reshape(magz.shape)

    for j in range(n_mag): map_prof0[:,j] = interp(theta0,theta[j],map_prof0[:,j])
    
    
  
    A = exp(1j*outer(arange(1, n_modes),theta0) )
    for it in trange(tsteps,desc='Poloidal modes: '): 
        map_prof = map_coordinates((gres[...,it]/(gres_mean+1e-6)).T-1,idx,order=2)
        map_prof = map_prof.reshape(magz.shape)
        for j in range(n_mag):
            map_prof[:,j] = interp(theta0,theta[j,:],  map_prof[:,j])
        complex_profile[:,:,it] = fft.rfft(map_prof,axis=0)[1:n_modes].T

    
    profile = mean(abs(complex_profile/n_theta)**2, -1)
    phase = sign(median(diff(angle(complex_profile),axis=-1),axis=-1))

            
    fig = Figure((img_size,img_size))
    FigureCanvas(fig)
    fig.clf()
    vmax = profile[:-n_mag//10].max()
    ax = fig.add_subplot(111)
    im = ax.imshow(phase.T*profile.T, aspect='auto',interpolation='nearest', 
                extent=(0-.5/(n_mag),1+.5/(n_mag),.5,n_modes-.5),origin='lower',
                cmap=cmap_neg,vmin=-vmax, vmax=vmax)

    ax.set_xlabel(r'$\rho_{pol}$')
    ax.set_ylabel('poloidal mode number $m$')
    cb = fig.colorbar(im)
    cb.set_label('Relative aplitude')
    fig.savefig(tmp_path+'/PoloidalModeNumber.png', transparent=True, dpi=dpi)




def plot3Dprofile(G, tokamak):
    #3D plot 

    G = G.mean(-1)/100000.

    from mpl_toolkits.mplot3d import Axes3D
    fig = Figure()
    FigureCanvas(fig)
    ax = fig.add_subplot(111, projection='3d')
    G[G==0] = nan
    R,Z = meshgrid(tokamak.xgrid, tokamak.ygrid)
    fig.set_size_inches(8,6)
    ax = fig.gca(projection='3d')
    ax.pbaspect = [2.0, 0.6, 0.9]
 
    ax.plot_surface(R, Z, G, rstride=1, cstride=1,facecolor='w',lw=.1,color='w')

    ax.set_xlabel('R [m]')
    ax.set_ylabel('z [m]')
    ax.set_zlim(0,None)

    ax.auto_scale_xyz([1,2.2], [-1,1], [0, G.max()])
    ax.set_aspect('equal')

    show()
    
    
    import IPython
    IPython.embed()


def add_plot_description(shot,time,tok,diag,ax=None,txt=None):

    description = '%s %s #%d t= %.4fs'%(tok,diag,shot, time)
    
    if not txt is None: 
        txt.set_text(description)
        return 
    
    
    if ax is None: ax = gca()


    txt = text(1.02,.06,description,rotation='vertical', 
                transform=ax.transAxes,verticalalignment='bottom',
                size='xx-small',backgroundcolor='none')
    
    return txt




def colorLOS_plot(time,data,error, magx_all,magy_all,extent,G,show_background=False):
    
    import matplotlib.pylab as plt
    global inputs, tokamak, plot_details

    xmin, xmax, ymin, ymax = plot_details['geometry']
    
    line_c = 'w' if show_background else 'k'
    ax = plt.subplot(111)
    if show_background:ax.imshow(G, origin='lower', interpolation='nearest', extent=extent)

    Rc,Zc = magx_all[:,0].mean(), magy_all[:,0].mean()
    
    if hasattr(tokamak,'struct_dict'):
        for _,struct in tokamak.struct_dict.items():
            if size(struct)>1:
                ax.plot(struct[0], struct[1], line_c,lw=.5)
       
    from geom_mat_setting import loadgeometry
    xchords, ychords, distance, nl,virt_chord  = loadgeometry(tokamak.geometry_path, list(tokamak.detectors_dict.keys()), 1)   
  
        
    lengths = hypot(xchords[-1]-xchords[0],ychords[-1]-ychords[0])

    data[isinf(error)| (error < 1)]= nan
    data = data/lengths
    data = data/amax(data[isfinite(data)])

    dxm=tokamak.dx/tokamak.norm
    dym=tokamak.dy/tokamak.norm
    for d,x,y in zip(data, xchords.T,ychords.T):
        #YlOrRd_r
        if isnan(d): continue
        plt.plot(x,y,c=cm.jet(d),alpha=.9,lw=3)

    plt.axes().set_aspect('equal')
    ax.plot(magx_all[:,-1],magy_all[:,-1],lw=1,c=line_c)
    ax.plot(Rc,Zc,'+',c=line_c,ms=10)
    ax.axis(extent)
    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    add_plot_description(tokamak.shot,time,tokamak.name,tokamak.input_diagn,ax=ax)
    plt.show()
  




if __name__ == "__main__":

    global inputs, plot_details,tokamak
    input = load('input.npz')['input']
    inputs= load('input.npz')['inputs'].item()
    plot_details= load('input.npz')['plot_details'].item()
    tokamak= load('input.npz')['tokamak'].item()

    matplotlib_data_tg2(input)

