#!/usr/bin/env python 
# -*- coding: utf-8 -*-

from numpy import *
from matplotlib.pylab import *
from scipy.signal import argrelmin, find_peaks_cwt, medfilt, argrelmax

from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d
from scipy.stats.mstats import mquantiles ,skew
from tqdm import trange
from shared_modules import  fsvd, poly_min
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvas as FigureCanvas
from matplotlib.ticker import MaxNLocator,NullFormatter
from multiprocessing import Process

#29460 - doctr blbě, ani nic neidentifikovalo 
#identifikovat sawotooths z profilu s pak počítat FSA jen z těch v okkolí? 

def MAD(x):     return median(abs(x-median(x)))*1.3

def win_diff(x,w=1):
    dx = zeros_like(x)
    dx[w:] = x[:-w]
    dx[:-w]-= x[w:]
    dx/= w
    dx[:w] = 0
    dx[-w:] = 0
    return -dx




def inverse_radius_tcv(tokamak,rho, events,sawtooths):
    #vybrat ten nejvnitřnější co obepíná střed, dát tam nějaký rozsah, minimální a maximální voliskoti? 
    #kreslit ten invezní poloměr do tomografie? 
    #less sensitive to the error in equilibrium
    import matplotlib._cntr as cntr
    R,Z = meshgrid(tokamak.xgrid,tokamak.ygrid )
    _,rc,zc = tokamak.mag_equilibrium(sawtooths, surf_slice=[0,])
    rc = rc.mean(0)[0]
    zc = zc.mean(0)[0]
    
    
    rhop,magx, magy = tokamak.mag_equilibrium( sawtooths,surf_slice=slice(6,None,6),n_rho=11)
    rho_inv = []
    
    
    for i, t  in enumerate(sawtooths):
        u,s,v = fsvd(events[i],2,i=0)

        
        c = cntr.Cntr(R,Z, u[:,1].reshape(tokamak.ny,tokamak.nx,order='F') )
        nlist = c.trace(level0=0,level1=0,nchunk=0)
        lines = nlist[:len(nlist)//2]
        lines = [l for l in lines if len(l)>10 ]
        dist = []
        for l in lines:  dist.append(mean(hypot(l[:,0]-rc[i], l[:,1]-zc[i])))
        line = lines[argmin(dist)]

        #works also if the plasma position is not right! 
        r_inv = hypot(line[:,0]-mean(line[:,0]), line[:,1]-mean(line[:,1])).mean()
        r_rho = hypot(magx[:,:,i]-rc[i],magy[:,:,i]-zc[i]).mean(0)
        rho_inv.append(interp(r_inv, r_rho,rhop ))
        
    
    return sawtooths, rho_inv
    
    

def inverse_radius(rho, events,sawtooths):
    #alternative!
    #https://crpplocal.epfl.ch/wiki/SoftX_analysis#Hints_about_.CF.81inv_computation

    
    #compute an inversion radius for the given sawtooth careshes (events) 
    std_emiss = r_[[e.std(0)/e.mean(0) for e in events]]

    t,r_ind = argrelmin(std_emiss,axis=1, order=sum(rho<.05))

    rho = linspace(0,1,std_emiss.shape[1])


    rho_inv = []
    rho_inv_tvec = []
    dx = mean(diff(rho))

    #subpixel precision of the minimum 
    for i in range(len(r_ind)):
        it = t[i]
        x = rho[r_ind[i]]
        if x > .6 or x < 0.05:
            continue
    
        y1,y2,y3 = std_emiss[it,r_ind[i]-1],std_emiss[it,r_ind[i]],std_emiss[it,r_ind[i]+1]

        rho_inv.append(poly_min(x, dx, (y1,y2,y3)))
        rho_inv_tvec.append(sawtooths[it])

    return rho_inv_tvec, rho_inv
    


def sawtooths_detection(input_data):
    #BUG a lot of hardcoded constants!!
    print('============== Sawtooths detection ===========')
    #load data and all necessary variables 
    global   plot_details, my_cmap

    input_parameters, tokamak, progress,results = input_data

    from shared_modules import read_config
    
    local_path = input_parameters['local_path']
    tmp_folder = input_parameters['tmp_folder']
    output_path = input_parameters['output_path']
    output_type = input_parameters['output_type']
    input_parameters = read_config(local_path+'/tomography.cfg')

    tvec = results['tvec']
    G = results['g'].T

    threshold = input_parameters['threshold']
    min_dt = input_parameters.get('min_dt',1e-3)
    max_dt = input_parameters['max_dt']
    crash_time = input_parameters['crash_time']

    tsteps = len(tvec)
    
    
    if len(tvec) < 20:
        raise Exception('Too short signal')

    
    #calculate a FSA emissivity (duplicated by postprocessing!!)
    tind = int_(linspace(0,len(tvec)-1,20))
    rhop,magx, magy = tokamak.mag_equilibrium(tvec[tind],return_mean=False)

    n_mag = size(magx,1)
    emiss = zeros(( tsteps,n_mag))
    rho = linspace(0,1,n_mag)

    scaling = array([tokamak.dx,tokamak.dy])
    offset = array([tokamak.xmin,tokamak.ymin])+scaling/2
    
    from scipy.ndimage.interpolation import map_coordinates

    
    gres = G.reshape(tokamak.ny,tokamak.nx, tsteps, order='F')

    for i in trange(tsteps,desc='Mapping Emiss to rho: '): 
        it = argmin(abs(tvec[tind]-tvec[i]))
        coords = c_[magx[:,:,it].ravel(),magy[:,:,it].ravel()].T
        idx = (coords-offset[:,None])/ scaling[:,None]
        map_prof = map_coordinates(gres[...,i].T,idx,order=1)
        map_prof = map_prof.reshape(magy.shape[:2])
        emiss[i]   = mean(map_prof,0)

    #calculate a detection signal 
    u,s,v = fsvd(emiss,2)
 
    sig = u[:,1]/(u[:,0]+1e-7)/mean(u[:,1]/(u[:,0]+1e-7))
   

    dt = (tvec[-1]-tvec[0])/len(tvec)
 
    #simple high pass filter
    sig = -win_diff(sig,max(1,int(crash_time/dt)))
    
    #remove background
    n_smooth = minimum(len(sig)//2,  2*(int(max_dt/dt)//2)+1)
    

    for i in range(len(sig)//n_smooth+1):
        ind = slice(i*n_smooth, max((i+1)*n_smooth,len(sig)))
        sig[ind]-= median(sig[ind])

    #normalize 
    n = 6*(int(max_dt/dt)//2)+1

    for i in range(len(sig)//n+1):
        ind = slice(i*n, min((i+1)*n,len(sig)))
        sig[ind]/= 2*median(abs(sig[ind]))
 
    if not all(isfinite(sig)):
        print('Sawtooth indicator is not finite')
        return

    noise_perc = sum(abs(sig-median(sig))/MAD(sig)<3)/float(len(sig))*100


    orientation = sign(skew(sig))
    sig*= orientation


    #identify positions of the sawtooths
    try:
        i_sawtooth = array(argrelmax(sig, order= int(2*crash_time/dt)))
        i_sawtooth = i_sawtooth[sig[i_sawtooth]> threshold   ]
    except:
        i_sawtooth = []
    

    
    
    fig = Figure((15,6))
    FigureCanvas(fig)
    
    #plot results
    axis = empty((2,2), dtype=object)
    for i in range(axis.size):
        axis.flat[i] = fig.add_subplot(2,2,i+1, sharex=axis[0,0])

    fig.subplots_adjust(hspace=0.07, wspace = 0.07)


    sawtooths = tvec[i_sawtooth]
    sawtooths = [poly_min(tvec[i],(tvec[i+1]-tvec[i-1])/2,sig[i-1:i+2]) for i in i_sawtooth]
    sawtooths = array(sawtooths)
    
   
    [axis[0,0].axvline(x=t,c='k',lw=.3) for t in sawtooths]
    axis[0,0].plot(tvec, sig,lw=.5)
    axis[0,0].axhline(y=MAD(sig)  )
    axis[0,0].axhline(y=-MAD(sig) )
    axis[0,0].axhline(y=mquantiles(sig, 0.98)[0]*threshold,ls='--' )
    axis[0,0].plot(sawtooths, sig[ i_sawtooth],'o')

    axis[0,0].set_ylabel('Sawtooth signal')
    axis[0,0].set_xlim(tvec[0],tvec[-1])
    axis[0,0].set_ylim(mquantiles(sig,[0.02,1]))
    

    dt = diff(sawtooths)*1e3
    axis[1,0].plot((sawtooths[1:]+sawtooths[:-1])/2,dt,'x')
    axis[1,0].set_xlabel('time [s]')
    axis[1,0].set_ylabel('$\Delta$ t [ms]')
    axis[1,0].grid(True)
    axis[1,0].set_ylim(0.,mquantiles(diff(tvec[i_sawtooth])*1e3,0.95)[0]*1.3)
    axis[1,0].set_xlim(tvec[0],tvec[-1])
    if size(dt) > 2:
        axis[1,0].axhline(y=median(dt))
        axis[1,0].text(.2,.8,str('%.1fms'%median(dt)),
                    transform=axis[1,0].transAxes)

    #analyze a sigle sawtooth crashes to identify an inverion radius
    sawtooths = sawtooths[(sawtooths< tvec[-1]-4e-3)&(sawtooths> tvec[0]+2e-3)]
    events_fsa = [emiss[tvec.searchsorted(s-2e-3):tvec.searchsorted(s+4e-3)-1] for s in sawtooths]
    events = [G[:,tvec.searchsorted(s-2e-3):tvec.searchsorted(s+4e-3)-1] for s in sawtooths]
 
    for label in axis[0,0].get_xticklabels():   label.set_visible(False)
    axis[0,0].xaxis.offsetText.set_visible(False)
    
    for label in axis[0,1].get_xticklabels():   label.set_visible(False)
    axis[0,1].xaxis.offsetText.set_visible(False)
    
    


    axis[0,1].grid(True)
    axis[0,1].set_ylabel('inversion radius '+tokamak.rho_label)

    axis[0,1].set_ylim(0,.5)
    axis[0,1].set_xlim(tvec[0],tvec[-1])
    axis[0,1].yaxis.tick_right()
    axis[0,1].yaxis.set_label_position("right")
 
    
    if len(events_fsa) > 0:  #evaluate only if some crashes was detected
        rho_inv_tvec, rho_inv = inverse_radius(rho,  events_fsa,sawtooths)

        axis[1,1].axhline(y=median(rho_inv),c='w')
        axis[0,1].plot(rho_inv_tvec,rho_inv,'b.')

        axis[0,1].axhline(y=median(rho_inv))
        axis[0,1].text(.1,.8,str('$\\rho_{inv}$=%.2f'%median(rho_inv))
                   ,transform=axis[0,1].transAxes)
        
    axis[1,1].set_xlabel('time [s]')
    #axis[1,1].set_ylim(0,.5)
    [axis[1,1].axvline(x=t+mean(diff(tvec)),c='w',lw=.5) for t in sawtooths]
    [axis[1,1].axvline(x=t+mean(diff(tvec)),c='k',lw=.5,ls='--') for t in sawtooths]

    axis[1,1].yaxis.tick_right()
    axis[1,1].yaxis.set_label_position("right")
    axis[1,1].set_ylabel('emissivity/max(emissivity)')
    axis[1,1].set_xlim(tvec[0],tvec[-1])

   
    extent=(tvec[0],tvec[-1],0,1)
    axis[1,1].imshow(   (emiss/amax(emiss,1)[:,None]).T, vmin=0,
               vmax=mquantiles(emiss/amax(emiss,1)[:,None],.999), 
               extent=extent,origin='lower',aspect='auto')
    
    axis[1,1].set_ylim(0,1)
    fig.savefig(tmp_folder+'/sawtooths'+str(tokamak.shot)+'.png',
                bbox_inches='tight',transparent=True, dpi=82)

    if input_parameters['enable_output']:
        fig.savefig(output_path+'/sawtooths'+str(tokamak.shot)+'.'+output_type,
                bbox_inches='tight',transparent=True, dpi=82)

    
    savetxt(tmp_folder+'/sawtooths_%d.txt'%tokamak.shot,sawtooths,
                fmt='%.5g', header='sawtooth crashes [s]')
    if len(events_fsa) > 0: 
        savetxt(tmp_folder+'/sawtooth_inversion_radius_%d.txt'%tokamak.shot,
                c_[rho_inv_tvec,rho_inv] ,fmt='%.5g',
                header='sawtooth crashes [s], inversion radius')
        
        
        
        
        
