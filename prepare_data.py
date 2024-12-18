#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
#from matplotlib.pyplot import *
from scipy.io import loadmat 
import os
import warnings
import time
from scipy import sparse
from scipy.sparse import spdiags, eye
from geom_mat_setting import geom_mat_setting
#from  prepare_data_golem import Golem
#from prepare_data_JET import JET
#from prepare_data_compass import COMPASS
#from prepare_data_ToreSupra import ToreSupra
#from prepare_data_ASDEX import ASDEX
from geometry import *
from matplotlib import rcParams
import gc
from shared_modules import in1d, debug
from matplotlib import colors
import config
#from matplotlib.ticker import  MaxNLocator 
 

def loaddata(inputs, useCache = True, prepare_tokamak = False):
    """
    Prepare and create object Tokamak containing all useful information for reconstruction and postprocessing.
    If object Tokamak already exists in ``inputs`` and input parameters are the same, the function return 
    ``inputs``. Otherwise new object is created.

    :param dict inputs: dictionary with all setting generated by `main` and `loadSetting`
    :rtype: dict
    """

    try:   # prevent data reloading !!! And save them when possible
        #assert useCache

        tok = inputs['tokamak_tmp']
        tok_changed = inputs['tok_index'] != tok.index
        
        shot_changed = inputs['shot'] != tok.shot
        resolution_changed  = inputs['nx'] != tok.nx or inputs['ny'] != tok.ny

        tok.set_resolution(inputs['nx'], inputs['ny'])
        config.wrong_dets_pref = unique(r_[config.wrong_dets_defined,config.wrong_dets_pref])
        if not tok_changed and not shot_changed and not resolution_changed:
            return tok
        
        if tok_changed:
            print("================== Tokamak/Diagnotics changed ===============")
            inputs['wrong_dets'] = []
            config.wrong_dets_pref = []
            config.wrong_dets_defined = []
     
            
        if shot_changed:
            print("================== Shot changed ===============")
            tok.shot = inputs['shot']
            inputs['wrong_dets'] = []
            config.wrong_dets_pref = []
            config.wrong_dets_defined = []
            
        #if resolution_changed:
            #print('!!!!!!!!!!Regenerate matrix of geometry?? when is it called??')
            #tok.Tmat , Xchords, Ychords = geom_mat_setting(tok, tok.nx,tok.ny,tok.virt_chord )
        
    except Exception as e:
        pass

    only_prepare =True
 
    tokamak_list =  [
        [Golem, 'bolometry'], [Golem, 'camera'],
        [COMPASS, 'BOLO'],  [COMPASS, 'SXR'],   [COMPASS, 'camera'],
        [JET, 'SXR-slow'], [JET, 'SXR-fast'],  [JET, 'BOLO'],  [JET, 'neutrons'],
        [ASDEX, 'SXR'], [ASDEX,'SXR_fast'], [ASDEX, 'BOLO'],[ASDEX, 'AXUV'],
        [ToreSupra, 'SXR'], [ToreSupra, 'Camera'],
        [TCV, 'XTOMO'], [TCV,'AXUV'],[TCV,'BOLO'],[TCV,'DMPX'],[TCV,'FIR'],[TCV,'XTEPRO'],
        [DIIID, 'SXR'], [DIIID,'SXR fast'],[DIIID,'BOLO'],[DIIID,'DISRAD'],
        ]
    
    tok_ = tokamak_list[inputs['tok_index']][0](tokamak_list[inputs['tok_index']][1],  inputs,  only_prepare )
    
    inputs['tokamak_tmp'] = tok_

    # TEST OF USER SELECTED WRONG DETECTORS

    if  ((not size(tok_.wrong_dets)==0) and all([i not in arange(size(tok_.Tmat,0)) for i in tok_.wrong_dets])):
        tok_.wrong_dets = []
        print("Incorrect wrong detectors!")

    gc.collect()
    
    if prepare_tokamak:
        tok_.prepare_tokamak()
 
    return tok_


   


def preview(fig, inputs, tokamak, plot_2D = True, show_prew= False):
    """
    Make preview for GUI, If the GUI is not enabled (``show_prew`` is ``True`` it create its own GUI with plot. Output is saved to ./tmp file and loaded by GUI.

    :param inputs dict:  dictionary with all setting generated by `main` and `loadSetting`
    :param bool plot_2D:  The output plot is 2D pcolor graph if ``True``, if ``False`` use 1D plot and if tsteps == 1 plot expected errorbars.
    :param bool show_prew:  If ``True`` create new windows with plot. Cannot be used together with GUI.
    """

    from make_graphs import my_cmap_

    
   
    data, error, tvec, dets, Tmat, normData = tokamak.prepare_data(inputs['tmin'],
                inputs['tmax'], inputs['data_smooth'], inputs['data_undersampling'],
                 inputs['error_scale'],tokamak.prune_dets,False)
    
  
    tsteps = sum((tvec >= inputs['tmin']) & (tvec <= inputs['tmax']))
    data_undersampling_tmp = max(inputs['data_undersampling'], int(tsteps/500.0))
    
    
    inputs['tsteps'] = tsteps
 
    wrong_dets = where( all(~isfinite(error), 0) )[0]
 
 
    ind_correct = tokamak.get_correct_dets(data, wrong_dets = wrong_dets)
    data = copy(data)
    nch = data.shape[0]
    
    from scipy.stats.mstats import mquantiles

    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    fig.clf()

    fig.subplots_adjust(left=.15-plot_2D*.05,top=.95,right=.90+plot_2D*.05,bottom=.1)
    from matplotlib.ticker import MaxNLocator

    font_size = 10
    ax = fig.add_subplot(111)
    
    if len(tvec) > 1 :
        data[~ind_correct, :] = NaN

        ind = slice(None,None, max(1, int(len(tvec)/1000)))

        data = data[:,ind]
        tvec = tvec[ind]
        vmax = mquantiles(data[isfinite(data)], 0.99)[0]

    if nanmax(data) > 2e5:
        fact,pre = 1e6, 'M'
    elif nanmax(data) > 2e2:
        fact,pre = 1e3, 'k'
    else:
        fact,pre = 1e0, ''

    ax.point_label = ax.text(0,0,'', fontsize=8 ,zorder=99)
    if len(tvec) > 1 and plot_2D == True:
        
        extent = (0.5, nch+.5,tvec[0],tvec[-1])
        #TODO remove in future
        try:
            im = ax.imshow(data.T/fact,cmap = my_cmap_,extent=extent,vmin=0,vmax=vmax/fact
                        ,interpolation='nearest', origin='lower', pickradius=1)
        except:
            im = ax.imshow(data.T/fact,cmap = my_cmap_,extent=extent,vmin=0,vmax=vmax/fact
                ,interpolation='nearest', origin='lower', picker=True)  
            ax.xaxis.set_pickradius(1)
            ax.yaxis.set_pickradius(1)  
        ax.axis('tight')

        ax.set_title('Input data, number of time steps: %i' % tsteps,fontsize=font_size)
        #xlabel('Channel')
        ax.set_ylabel('Time ['+tokamak.t_name+']',fontsize=font_size)
        from make_graphs import LogFormatterTeXExponent
        cbar = fig.colorbar(im,format=LogFormatterTeXExponent('%.2e'))
        ax.axis((0.5, nch+.5,tvec[0],tvec[-1]))
        cbar.ax.tick_params(labelsize=font_size ) 
        #from matplotlib import ticker
        cbar.locator = MaxNLocator(nbins=5)
        cbar.update_ticks()
        cbar.set_label('Intensity [%sW/m$^2$]'%pre,fontsize=10)

        if not tokamak.impur_inject_t is None:
            for lbot in tokamak.impur_inject_t:
                ax.axhline(y=lbot, lw=1, color='w') 
        
        ax.set_xlabel('Channel',fontsize=font_size)
        rcParams['xtick.direction'] = 'out'
        rcParams['ytick.direction'] = 'out'
        
    else:  
        if len(tvec) == 1:
            """
            try:#TODO remove in future
                ax.errorbar(arange(len(data[:,0]))[ind_correct]+1, data[ind_correct,0]/fact
                        ,yerr=error[ind_correct,0]/fact,fmt='-xb',lw=.5, pickradius=2)
                ax.plot(arange(len(data[:,0]))[~ind_correct]+1, data[~ind_correct,0]/fact,'or', picker=2)
            except: """
            ax.errorbar(arange(len(data[:,0]))[ind_correct]+1, data[ind_correct,0]/fact
                    ,yerr=error[ind_correct,0]/fact,fmt='-xb',lw=.5,   picker=True)
            ax.plot(arange(len(data[:,0]))[~ind_correct]+1, data[~ind_correct,0]/fact,'or' , picker=True)    
            ax.xaxis.set_pickradius(5)
            ax.yaxis.set_pickradius(5.)   
                
            ax.set_title('Input data, time: %g %s' % (tvec, tokamak.t_name),fontsize=font_size)

            if not any(ind_correct): ind_correct[:] = True
            ax.axis([0,nch+1, 0 , nanmax(data[ind_correct,0])/fact])
            ax.set_xlabel('Channel',fontsize=font_size)

        else:
            ax.plot(tvec, data[dets,:].T/fact ,lw=.2)
            ax.set_xlim(tvec.min(),tvec.max())
            ax.set_xlabel('t [s]',fontsize=font_size)
            ax.set_title('Input data, number of time steps: %i' % tsteps,fontsize=font_size)
            ax.set_ylim( -vmax*0.1/fact,vmax*1.2/fact)
            ax.axhline(0,c='k')
            ax.xaxis.set_major_locator(MaxNLocator(5))

            if not tokamak.impur_inject_t is  None:
                for lbot in tokamak.impur_inject_t:
                    ax.axvline(x=lbot, lw=1, color='k') 
        
        ax.set_ylabel('Intensity  [%sW/m$^{2}$]'%pre,fontsize=font_size)
        
        
    if plot_2D or len(tvec) == 1:
        for ind in tokamak.dets_index[:-1]:
            ax.axvline(x=1.5+amax(ind), linestyle='-' ,color='k')
            ax.axvline(x=1.5+amax(ind), linestyle='--',color='w')


    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(font_size) # Size here overrides font_prop
 

    ax.xaxis.set_major_locator(MaxNLocator(12))





