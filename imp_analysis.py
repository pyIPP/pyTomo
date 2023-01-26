#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
from scipy.interpolate import RectBivariateSpline, interp1d
from matplotlib.colors import LogNorm,NoNorm
import os, sys
from matplotlib.ticker import MaxNLocator,NullFormatter
from matplotlib.transforms import Bbox
from scipy.stats.mstats import mquantiles
from matplotlib.backends.backend_agg import FigureCanvas as FigureCanvas
from matplotlib.figure import Figure
from multiprocessing import Process
from time import sleep
from matplotlib.pylab import *    
from IPython import embed 



def update_fill_between(fill,x,y_low,y_up,min,max ):
    paths, = fill.get_paths()
    nx = len(x)
    y_low = maximum(y_low, min)
    y_low[y_low==max] = min
    y_up = minimum(y_up,max)
    paths.vertices[1:nx+1,1] = y_up
    paths.vertices[nx+1,1] = y_up[-1]
    paths.vertices[nx+2:-1,1] = y_low[::-1]
    paths.vertices[0,1] = y_up[0]
    paths.vertices[-1,1] = y_up[0]
    

    

#only for SXR AUG
def imp_analysis( input_data):

    print('Impurity analysis...')
    
    input_parameters, tokamak, progress,results = input_data
  
    tvec = results['tvec']
    fsa_emiss = results['fsa_emiss']
    shot = tokamak.shot
    rho = linspace(0,1,fsa_emiss.shape[1])
    G = results['g'].reshape(-1, tokamak.ny, tokamak.nx, order='F')
    magx = results['magx']
    magy = results['magy']
    nt, nr, n_theta = magx.shape
    indR = argmax(magx.reshape(-1,n_theta),1)
    LFS_R = magx.reshape(-1,n_theta)[arange(nr*nt),indR].reshape(nt,nr)
    LFS_z = magy.reshape(-1,n_theta)[arange(nr*nt),indR].reshape(nt,nr)
    indR = argmin(magx.reshape(-1,n_theta),1)
    HFS_R = magx.reshape(-1,n_theta)[arange(nr*nt),indR].reshape(nt,nr)
    HFS_z = magy.reshape(-1,n_theta)[arange(nr*nt),indR].reshape(nt,nr)
    
    plot_autoscale = input_parameters['plot_autoscale']
    geometry_path = tokamak.geometry_path
    tmp_folder = tokamak.tmp_folder
    output_path = input_parameters['output_path']
    dpi = input_parameters['dpi']

    HFS_emiss = zeros((nt,nr),dtype=single)
    LFS_emiss = zeros((nt,nr),dtype=single)

    scaling = array([tokamak.dx,tokamak.dy])
    offset = array([tokamak.xmin,tokamak.ymin])#+scaling/2
    
    from scipy.ndimage.interpolation import map_coordinates    

    for i in range(nt): 
        coords = c_[LFS_R[i].ravel(),LFS_z[i].ravel()].T
        idx = (coords-offset[:,None])/ scaling[:,None]
        LFS_emiss[i] = map_coordinates(G[i].T,idx,order=2)
        coords = c_[HFS_R[i].ravel(),HFS_z[i].ravel()].T
        idx = (coords-offset[:,None])/ scaling[:,None]
        HFS_emiss[i] = map_coordinates(G[i].T,idx,order=2)
        
            
    try:
        
        BRdata = load(geometry_path+'Bremsstrahlung_%d.npz'%shot)
        tvec_data =   double(BRdata['tvec'])
        rho_data_ = tile(BRdata['rho'], (len(tvec_data), 1))
        BR_ = double(BRdata['BR'])
        BR_low_ =double( BRdata['BR_low'])
        BR_high_ = double(BRdata['BR_high'])
        BR_high_[isinf(BR_high_)] = amax(BR_high_[isfinite(BR_high_)])
        W_ =  double(BRdata['W'])
        try:
            W_low_ =   double(BRdata['W_low']) #optimistics guess
            W_high_ =    double(BRdata['W_high'])
        except:
            W_low_ =  W_*.8
            W_high_ = W_*1.2
            
            
            
    except:

        try:
            data = load('./asymmetry/asymmetry_emiss_%d_.npz'%shot)
            tvec_data = data['tvec']
            rho_data_ = data['rho_pol']

            BR_ = data['Bremsstrahlung']
            W_ = data['W_Radiation']
        
            data_low  = load('./asymmetry/asymmetry_emiss_%d_low.npz'%shot) 
            data_high = load('./asymmetry/asymmetry_emiss_%d_up.npz'%shot )
            
            BR_low_ = data_low['Bremsstrahlung']
            W_low_ = data_low['W_Radiation']
            BR_high_  = data_high['Bremsstrahlung']
            W_high_  = data_high['W_Radiation']
        
            
        except:
            print('no radiation data!!!')
            return 
        
            
    
    

    BR_[isinf(BR_)] = nan
    W_[isinf(W_)] = nan
    
        
    BR = zeros((len(tvec), BR_.shape[1]))
    BR_low = zeros((len(tvec), BR_low_.shape[1]))
    BR_high = zeros((len(tvec), BR_high_.shape[1]))
    W = zeros((len(tvec), W_.shape[1]))
    W_low = zeros((len(tvec), W_low_.shape[1]))
    W_high = zeros((len(tvec), W_high_.shape[1]))
    rho_data = zeros((len(tvec), rho_data_.shape[1]))
    dt = mean(diff(tvec)) if size(tvec) > 1  else 0

    for it,t in enumerate((tvec)):
        ind = slice(tvec_data.searchsorted(t-dt),tvec_data.searchsorted(t+dt)+1)
        #embed()
        BR[it] = nanmean(BR_[ind],0)
        BR_low[it] = nanmean(BR_low_[ind],0)
        BR_high[it] = nanmean(BR_high_[ind],0)
        W[it] = nanmean(W_[ind],0)
        W_low[it] = nanmean(W_low_[ind],0)
        W_high[it] = nanmean(W_high_[ind],0)
        rho_data[it] = nanmean(rho_data_[ind],0)

    
    if tokamak.name == 'ASDEX':
        zeff = input_parameters['zeff']
        print('Assume Zeff = ', zeff)
        BR*= zeff/1.1
        BR_low*= zeff/1.1  #default zeff used to calc Bremsstrahlung is 1.1!!
        BR_high*= zeff/1.1
    else:
        zeff = 1

    


    BR = array([interp(rho,r, p ) for r,p in zip( rho_data, BR)])
    BR_low = array([interp(rho,r, p ) for r,p in zip( rho_data, BR_low)])
    BR_high = array([interp(rho,r, p ) for r,p in zip( rho_data, BR_high)])
    W = array([interp(rho,r, p ) for r,p in zip( rho_data, W)])
    W_low = array([interp(rho,r, p ) for r,p in zip( rho_data, W_low)])
    W_high = array([interp(rho,r, p ) for r,p in zip( rho_data, W_high)])
    W_high[isnan(W_high)] = infty
     
    W_high[isnan(W_high)] = infty
    BR_high[isnan(BR_high)] = infty




 
    BR_low[BR_low == 0] = BR_low[BR_low!=0].min()
    
    gmax = fsa_emiss.max(1)
    
    if gmax.max() > 1e6:
        fact, unit = 1e-6,'MW'
    elif gmax.max() > 1e3:
        fact, unit = 1e-3,'kW'
    else: fact, unit = 1,'W'
    if input_parameters['plot_loglin']:
        fact, unit = 1,'W'
        

    fig = Figure((17,6))
    FigureCanvas(fig)
    ax = []
    ax.append(fig.add_subplot(1,3,1))
    ax.append(fig.add_subplot(1,3,2,sharex=ax[0]))
    ax.append(fig.add_subplot(1,3,3,sharex=ax[0]))
    ax[0].set_xlim(0,1)

    from scipy.signal import order_filter

    gmax = order_filter(gmax,ones(11),10)
    
    cw  = ones_like(fsa_emiss )*nan
    Zeff = ones_like(fsa_emiss )*nan
    Zeff_up = ones_like(fsa_emiss )*nan
    Zeff_low = ones_like(fsa_emiss )*nan

    for it,t in enumerate(tvec):
        if all(isnan(W[it])|isnan(BR[it])): continue
        br = copy(BR[it])
        br[isnan(br)] = br[isfinite(br)][0]
        w = copy(W[it])
        w[isnan(w)] = w[isfinite(w)][0]
        cw[it]  = (fsa_emiss[it]-br)/w
        Zeff[it] = (HFS_emiss[it])/(br/zeff)
        Zeff_up[it] = (HFS_emiss[it]-W_low[it]/2)/BR_low[it]
        Zeff_up[it,w/br>0.3] = 6
        Zeff_low[it] = (HFS_emiss[it]-W_high[it]*2)/(BR_high[it]/1.3/zeff*1.1)

    cw_max = order_filter(nanmax(cw[:,rho<0.4],1),ones(11),10)*2
    
    if not plot_autoscale:
        cw_max[:] = nanmax(cw_max)
    cw_max = maximum(cw_max, 1e-4)

    BR_err_plot = ax[0].fill_between(rho, -rho-1,1+ rho,alpha=.3,facecolor='b', edgecolor='None')
    BR_plot, = ax[0].plot(rho, BR[0],'b',label='Bremsstrahlung')

    W_err_plot = ax[0].fill_between(rho, -rho,rho,alpha=.3,facecolor='r', edgecolor='None')
    W_plot, = ax[0].plot(rho, W[0],'r',label='W radiation $c_w=10^{-4}$')
    
    Emiss_plt, = ax[0].plot(rho, fsa_emiss[0],'k',label='FSA emissivity')
    Emiss_lfs_plt, = ax[0].plot(rho, LFS_emiss[0],'k--')
    Emiss_hfs_plt, = ax[0].plot(rho, HFS_emiss[0],'k--')

    ax[0].set_xlabel(tokamak.rho_label,fontsize=15)
    ax[0].set_ylabel('Emissivity [%sm$^{-3}$]'%unit,fontsize=13)
    leg=ax[0].legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.7)
    
    for l in leg.legendHandles:  l.set_linewidth(10)
    ax[0].yaxis.set_major_locator(MaxNLocator(5))
        
    ax[1].set_ylabel('$10^4\ c_W$',fontsize=13)
    ax[1].set_xlabel(tokamak.rho_label ,fontsize=15)
    ax[1].yaxis.set_major_locator(MaxNLocator(5))
    
    
    cw_high = (fsa_emiss[0]-BR_low[0])/(W_low[0]+1e-1)
    cw_low  = maximum(0, (fsa_emiss[0]-BR_high[0])/W_high[0])
    cw_low[isnan(cw_low)]
    cw_err_plot = ax[1].fill_between(rho,-rho,rho,alpha=.3,facecolor='k', edgecolor='None')
    cw_plot, = ax[1].plot(rho, cw[it],'k')
    ax[1].set_ylim(0,max(1e-4,cw_max[it]))
        
    ax[2].set_ylabel('$Z_{eff}$',fontsize=13)
    ax[2].set_xlabel(tokamak.rho_label ,fontsize=15)
    ax[2].yaxis.set_major_locator(MaxNLocator(5))
    zeff_err_plot = ax[2].fill_between(rho,-rho,rho,alpha=.3,facecolor='k', edgecolor='None')
    zeff_plot, = ax[2].plot(rho, cw[it],'k')
    ax[2].set_ylim(0,3)
    ax[2].set_yticks([0,1,2,3])
    ax[2].axhline(1,c='k')

                    
    txt = ax[1].text(.05, .9,  '' ,transform=ax[1].transAxes)

                 
    plot_loglin = input_parameters['plot_loglin']
                 
    ax[0].set_ylim(0,max(nanmax(gmax), BR.max())*fact*1.2)
           

    if input_parameters['plot_all'] or len(tvec)== 1:
        from tqdm import tqdm,trange
        for it in trange(len(tvec),desc='Plotting impurity analysis: '):
            txt.set_text("time %.3fs"%tvec[it])
            
            if plot_autoscale:
                ax[0].set_ylim(0,max(gmax[it], BR[it].max())*fact*1.2)

            update_fill_between(BR_err_plot,rho,BR_low[it]*fact, BR_high[it]*fact, 0, ax[0].get_ylim()[1])
            BR_plot.set_ydata(BR[it]*fact)
    
            update_fill_between(W_err_plot,rho,W_low[it]*fact, W_high[it]*fact, 0, ax[0].get_ylim()[1])
            W_plot.set_ydata( W[it]*fact)

            Emiss_plt.set_ydata(fsa_emiss[it]*fact)
            Emiss_lfs_plt.set_ydata(LFS_emiss[it]*fact)
            Emiss_hfs_plt.set_ydata( HFS_emiss[it]*fact)

            if not all(isnan(W[it])): 
                ind = (fsa_emiss[it]-BR_low[it])/fsa_emiss[it] <.5
                cw_high = (fsa_emiss[it]-BR_low[it])/(W_low[it]+1e-6)
                cw_low  = (fsa_emiss[it]-BR_high[it])/(W_high[it]+1e-6)
                if isfinite(cw_max[it]):
                    ax[1].set_ylim(0,cw_max[it])
           
                update_fill_between(cw_err_plot,rho,cw_low,cw_high,0,cw_max[it])
                cw[it, ind] = nan
                cw_plot.set_ydata(cw[it])
            
            
            ax[1].set_ylabel('$10^4\ c_W$',fontsize=13)
            ax[1].set_xlabel(tokamak.rho_label,fontsize=15)
            ax[1].yaxis.set_major_locator(MaxNLocator(5))
            if plot_loglin:
                ax[0].set_yscale('symlog', linthreshy= BR[it].max()*fact)

            update_fill_between(zeff_err_plot,rho,Zeff_low[it],Zeff_up[it],0,6)
            zeff_plot.set_ydata( Zeff[it])

            try:
                fig.savefig(tmp_folder+'imp_plot_%.4d.png'%it, bbox_inches=Bbox([[1,0.05],[16,5.7]]),transparent='True',dpi=dpi)
            except Exception as e:
                print(e)
                

            if input_parameters['enable_output']:
                fig.savefig(output_path+'impurities_%.4d_%d.%s'%(it,shot,input_parameters['output_type']),
                            dpi=dpi, bbox_inches='tight', transparent=True)


    savez_compressed(tmp_folder+'/imp_analysis_%d'%shot, tvec=tvec,rho=rho,
                     BR=single(BR), 
                     W=single(W),
                     emiss=single(fsa_emiss), 
                     LFS_emiss=single(LFS_emiss), 
                     HFS_emiss=single(HFS_emiss),
                     cw=single(cw*1e-4) )
                
    if len(tvec) == 1:
        return 

    fig.clf()
    fig.set_size_inches((16,5.5))

    
    ax = []
    ax.append(fig.add_subplot(1,3,1))
    ax.append(fig.add_subplot(1,3,2,sharex=ax[0],sharey=ax[0]))
    ax.append(fig.add_subplot(1,3,3,sharex=ax[0],sharey=ax[0]))
    for a in ax[1:]:
        for label in a.get_yticklabels():
            label.set_visible(False)
        a.yaxis.offsetText.set_visible(False)

       

    fig.subplots_adjust(hspace=0.13, wspace = 0)

    
    vmax = nanmax(cw[:,rho<.4],1)
    vmax = mquantiles(vmax[isfinite(vmax)],0.9)[0]
    vmax = max(vmax,1e-5)
    cw[isnan(cw)] = 0
    try:
        levels=linspace(0,vmax,15)
    except:
        vmax = 1
        levels=linspace(0,vmax,15)

        
    CM = ax[0].contourf(tvec, rho, cw.T,levels, vmin=0,vmax=vmax, extend='both',cmap='jet')
    CM.cmap.set_over('brown')
    CM.cmap.set_under('black')
    from matplotlib import colorbar

    cax,kw = colorbar.make_axes(ax[0])
    cb = fig.colorbar(CM, cax=cax, **kw)
    tick_locator = MaxNLocator(nbins=4)   
    ax[0].set_title(r'W concentration [$10^{-4}$]')
    cb.locator = tick_locator
    cb.update_ticks()
    ax[0].set_xlabel('time [s]')
    ax[0].set_ylabel(tokamak.rho_label,fontsize=15)
    ax[0].set_xlim(tvec[0],tvec[-1])
    ax[0].set_ylim(0,1)

    
    vmax = nanmax(fsa_emiss,1)
    vmax = mquantiles(vmax,0.95)[0]
    levels=linspace(0,vmax,15)
    CM = ax[1].contourf(tvec, rho, fsa_emiss.T*fact,levels*fact, vmin=0,vmax=vmax*fact, extend='max',cmap='jet')
    CM.cmap.set_over('brown')
    cax,kw = colorbar.make_axes(ax[1])
    cb = fig.colorbar(CM, cax=cax, **kw)
    tick_locator = MaxNLocator(nbins=4)   
    ax[1].set_title(r'Emissivity [%s/m$^3$]'%unit)
    cb.locator = tick_locator
    cb.update_ticks()
    ax[1].set_xlabel('time [s]')

    
    vmax = nanmax(W,1)
    vmax = mquantiles(vmax[isfinite(vmax)],0.95)[0]
    levels=linspace(0,vmax,15)
    CM = ax[2].contourf(tvec, rho, W.T*fact,levels*fact, vmin=0,vmax=vmax*fact, extend='max',cmap='jet')
    CM.cmap.set_over('brown')
    cax,kw = colorbar.make_axes(ax[2])
    cb = fig.colorbar(CM, cax=cax, **kw)
    tick_locator = MaxNLocator(nbins=4)   
    ax[2].set_title(r'10$^{-4}$ W [%s/m$^3$]'%unit)
    cb.locator = tick_locator
    cb.update_ticks()
    ax[2].set_xlabel('time [s]')

    fig.savefig(tmp_folder+'impurities_2D_%d.png'%shot, dpi=dpi, bbox_inches='tight', transparent=True)


    if input_parameters['enable_output']:
        fig.savefig(output_path+'impurities_2D_%d.%s'%(shot,input_parameters['output_type']),
                    dpi=dpi, bbox_inches='tight', transparent=True)




