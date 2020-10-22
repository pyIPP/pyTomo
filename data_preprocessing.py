#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os, random
import matplotlib

try:
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

except:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *

        
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    try:
        from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
    except:
        from matplotlib.backends.backend_qt4agg import NavigationToolbar as NavigationToolbar
        
from matplotlib.figure import Figure
import numpy as np
import config


#import matplotlib.pyplot as plt
try:
    from signal_analysis.stft  import stft
    from signal_analysis.sfft  import sfft
    from signal_analysis.sstft import sstft
except:
    print('pyfftw is missing')


from matplotlib.widgets import Slider, RadioButtons,RectangleSelector
from matplotlib.ticker import NullFormatter, ScalarFormatter
from scipy.signal import argrelmax
import scipy.sparse as sp
from colorsys import hls_to_rgb
from matplotlib.ticker import MaxNLocator
import fconf

import os
from numpy import *
#from matplotlib import Widget
from scipy import sparse
import sys
#from matplotlib.pyplot import *
import time
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp2d,RectBivariateSpline
import os,sys
#from BlockConv import BlockConv
from scipy.signal import firwin2, freqz, lfilter, lfiltic,get_window,fftconvolve,firwin
from numpy.linalg import qr
import scipy.signal as signal
from scipy.fftpack import ifft,fft,fftfreq
from matplotlib.widgets import MultiCursor,Widget



#AUTOR: Tomas Odstrcil  tomas.odstrcil@ipp.mpg.de




#TODO okno s residuem - dát tam nahoru i jeden grafík s retrofitem vs orig data?? a měnit tanm obsh podle toho na který kanál kliknu, a spojit zoomování 
#TODO normální ono na výběr kanálů, zapínání/vypínánání kanálů prostředním tlačítkem, 


# TODO udělat aby to dovedlo trackovat frekvenci módu a pak to vzít v úvahu při friltraci
# TODO aby to nenačítalo data pořád zvova, 
# TODO proč se to 2x přepočte při změně kanálu? 
#TODO jak vyřešit dvojítý přepočet při změně velikosti? 

#proč je výsledek filtru, - komplexní fáze náhodné číslo? Zafixovat počáteční fázi? ale to nejde, výdyť to musí být fixní!! zafixovat počáteční fázi co od toho odčítám? 

#vzít v úvahu vzorkovací frekvenci kanálu
#udělat korekci aplitudy u těch starých kanálů!!! ipp report o SXRstr 13


class MyFormatter(ScalarFormatter):   # improved format of axis
    def __call__(self, x, pos=None):
        self.set_scientific(True)
        self.set_useOffset(True)
        self.set_powerlimits((-3, 3))
        return ScalarFormatter.__call__(self, x, pos)





def colorize(z):
    r = np.abs(z)
    arg = np.angle(z) 

    h = (arg + pi)  / (2 * pi) + 0.5
    l = 1.0 - 1.0/(1.0 + r**0.3)
    s = 0.8

    c = np.vectorize(hls_to_rgb) (h,l,s) # --> tuple
    c = np.array(c)  # -->  array of (3,n,m) shape, but need (n,m,3)
    c = c.T.swapaxes(0,1)
    return c


def comb_filter(win_t, dt, f0,n=5,window='hamming',pass_DC=True):

    harm =  arange(not pass_DC,n+1)

    tvec_win = arange(0,win_t, dt )

    if len(tvec_win)%2 == 0:
        tvec_win = tvec_win[:-1]

    win = exp(1j*2*pi*f0*tvec_win[None,:])**harm[:,None]

    win = r_[win.real, win.imag]
    win_fun = get_window(window, len(tvec_win))

    q,r = qr(win.T)

    filter = dot(q[len(tvec_win)//2,:], q.T)

    filter*= win_fun
    if pass_DC:
        filter/= sum(filter)
    
    return  filter


def comb_filter2(win_t, dt, f0,n=5,window='hamming',pass_DC=True):
    
    length = int(win_t/dt)
    
    harm =  arange(not pass_DC,n+1)
    
    f_pass = (harm[:,None] +array([0,1e-5]) ).ravel()*f0
    if pass_DC:
        f_pass = f_pass[1:]

    filter = firwin(length, f_pass, pass_zero=pass_DC,nyq=1/dt/2, window = window)
    return filter



def comb_filter3(win_t, dt, f0,n=5,window='hamming',pass_DC=True):
    
    n = 2*(int(win_t/dt)//2)+1
    filter = zeros(n)
    dn = 1/(f0*dt*pi)
    i = 0
    while(int(i*dn)<n/2):
        filter[n//2+int(i*dn)]=1
        filter[n//2-int(i*dn)]=1
        i+=1
    filter/= sum(filter)
    return filter



def filter_bank(  win_t, dt, f0,n=5,window='hamming' ):
    
    length = int(win_t/dt)

    filter = zeros((n+1,length))
    for i in arange(n+1):
        filter[i,:] = firwin(length, [(i+1e-3)*f0,(i+2e-3)*f0], pass_zero=False,
                        nyq=1/dt/2, window = window)

    return filter



class MultiCursor(Widget):
    """
    Provide a vertical line cursor shared between multiple axes

    Example usage::

        from matplotlib.widgets import MultiCursor
        from pylab import figure, show, np

        t = np.arange(0.0, 2.0, 0.01)
        s1 = np.sin(2*np.pi*t)
        s2 = np.sin(4*np.pi*t)
        fig = figure()
        ax1 = fig.add_subplot(211)
        ax1.plot(t, s1)


        ax2 = fig.add_subplot(212, sharex=ax1)
        ax2.plot(t, s2)

        multi = MultiCursor(fig.canvas, (ax1, ax2), color='r', lw=1)
        show()

    """
    def __init__(self, canvas, axes, useblit=True,vertical=True,horizontal=False, **lineprops):

        self.canvas = canvas
        self.axes = axes
        xmin, xmax = axes[-1].get_xlim()
        xmid = 0.5*(xmin+xmax)
        
        ymin, ymax = axes[-1].get_ylim()
        ymid = 0.5*(ymin+ymax)
        if useblit:
            lineprops['animated'] = True
            
        self.lines = []
        for ax in axes:
            self.lines.append({})
            if vertical:
                self.lines[-1]['x'] =  ax.axvline(xmid, visible=False, **lineprops)
            if horizontal:
                self.lines[-1]['y'] =  ax.axhline(ymid, visible=False, **lineprops)

        self.visible = True
        self.useblit = useblit
        self.background = None
        self.needclear = False

        self.canvas.mpl_connect('motion_notify_event', self.onmove)
        self.canvas.mpl_connect('draw_event', self.clear)


    def clear(self, event):
        'clear the cursor'
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(
                                            self.canvas.figure.bbox)
        for lines_dict in self.lines:
            for line in list(lines_dict.values()):
                line.set_visible(False)


    def onmove(self, event):
        if event.inaxes is None:
            return
        if not self.canvas.widgetlock.available(self):
            return
        self.needclear = True
        if not self.visible:
            return

        for lines_dict in self.lines:
            if 'x' in lines_dict:
                lines_dict['x'].set_xdata((event.xdata, event.xdata))
                lines_dict['x'].set_visible(self.visible)
            if 'y' in lines_dict:
                lines_dict['y'].set_ydata((event.ydata, event.ydata))
                lines_dict['y'].set_visible(self.visible)
            

        self._update()

    def _update(self):

        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for ax, lines_dict in zip(self.axes, self.lines):
                for line in list(lines_dict.values()):
                    ax.draw_artist(line)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()
            
            
            

class SVDFilter():
    #algorithm for quasiperiodic multichannel signal filtering

    def __init__(self,shot,fig,ch0,SpecWin=None,statusbar=None,slider_svd=None,slider_harm=None,ind=slice(None,None)):

        self.SXR_detectors=('G','H1','H2','H3','I1','I2','I3','J1','J2','J3','K1','K2','L','M')

        self.shot = shot
        self.ind = ind  #BUG optializovat to ať se to zbytečně nenačítá celé!!
        self.fig = fig
        self.SpecWin=SpecWin
        self.det_ind = {}

        self.window = 'blackman'
        self.ch0 = ch0#72 #BUG
        self.n_svd = 3#BUG
        self.n_harm = 4#BUG
        print('loading data')
        self.f0 = nan
        self.TSS = nan
        self.RSS = nan
        self.fig2 = None
        
        if os.path.isfile('data/tvec_%d.npz'%self.shot):
            self.path = 'data'
        else:
            self.path = '../geometry/ASDEX/SXR_fast/'
        
        self.slider_harm = slider_harm
        self.slider_svd = slider_svd
        self.statusBar = statusbar

        if not slider_harm is None:
            slider_harm.setRange(2, 10)
            slider_harm.setValue(self.n_harm)
            slider_harm.setTracking(True)
            slider_harm.setTickPosition(QSlider.TicksBelow)
            slider_harm.setTickInterval(1)
            
        if not slider_svd is None:
            slider_svd.setRange(2, 10)
            slider_svd.setValue(self.n_svd)
            slider_svd.setTracking(True)
            slider_svd.setTickPosition(QSlider.TicksBelow)
            slider_svd.setTickInterval(1)

        self.actual = False
            
    def apply_slider(self):
        
        if not self.slider_harm is None and self.n_harm!= self.slider_harm.value():
            self.n_harm =  self.slider_harm.value()
            self.statusBar().showMessage('N HARM = %d' % self.n_harm , 2000)

        if not self.slider_svd is None and self.n_svd!= self.slider_svd.value():
            self.n_svd  =   self.slider_svd.value()
            self.statusBar().showMessage('N SVD = %d' % self.n_svd , 2000)


        self.actual = False
        

        
        
            
    def load(self):
        
         #TODO spojis s tomogrefií !!

        #path = '../geometry/ASDEX/SXR_fast/'

        tvec = load(self.path+'/tvec_%d.npz'%(self.shot))['tvec'][self.ind]
        t_range = self.SpecWin.ax.get_xlim()

        iimin = tvec.searchsorted(t_range[0])
        iimax = tvec.searchsorted(t_range[1])
        ind = slice(iimin,iimax )
        self.tvec = tvec[ind]
        
        
        self.data_list,self.wrong_det = [],[]
        self.dets = [0,]
        self.ndets = 0


        for d in self.SXR_detectors:   
            #print self.ind
            self.data_list.append(load(self.path+'/%s_%d.npy'%(d,self.shot), mmap_mode='r')[self.ind])
            self.det_ind[d] = slice(self.ndets,self.ndets+self.data_list[-1].shape[1])

            try:
                detector_stat = load(self.path+'/'+d+"_"+str(self.shot)+"_stat.npy")
                self.wrong_det.append(where(~detector_stat)[0]+self.ndets)
            except Exception as e:
                print(e)
                pass
            self.ndets+= self.data_list[-1].shape[1]
            self.dets.append(self.ndets)



        self.data = empty((iimax-iimin,self.ndets ), dtype=single)
        from tqdm import tqdm
        for det_data,i,det in tqdm(list(zip(self.data_list,self.dets,self.SXR_detectors))):
            self.data[:,i:i+det_data.shape[1]] = det_data[ind]

        print('data loaded')
        
        #BUG načítat z tomogrefie!!
        self.wrong_det = hstack(self.wrong_det+ [where(mean(self.data,0)<0)[0],where(any(isnan(self.data),0))[0]])
        self.wrong_det = r_[self.wrong_det,80]
        self.wrong_det = unique(self.wrong_det)

        #print self.ndets, data.shape[1]
        #exit()
        #ndets = data.shape[1]
        
        #BUG!!
        self.err = ones(self.ndets)
        self.err[155:172] = 1e4
        self.err[172:] = 11
        self.err[:16] = 4
        self.err[28:30] = 4
        self.err[74] = 3
        self.err[87] = 3

        self.data[:,self.wrong_det] = 0
        #self.ch0 -= sum(self.wrong_det < self.ch0) 

        #t_win, frange =  0.0007,(32737, 12832)
        #t_range = (0,10)
        
        #t_win = abs(diff(self.SpecWin.t_range))
        #t_range = self.SpecWin.ax.get_xlim()
                
                
                
        #iimin = tvec.searchsorted(t_range[0])
        #iimax = tvec.searchsorted(t_range[1])
        #ind = slice(iimin,iimax )

        #self.tvec = copy(tvec[ind])
        #self.data = copy(data[ind,:])
        
        
    def run_filter(self):

        nt = len(self.tvec)
        dt = (self.tvec[-1]-self.tvec[0])/(len(self.tvec)-1)
        f =  fftfreq(nt, dt)[:nt//2+1]
        frange = self.SpecWin.f_range
        t_win = ptp(self.SpecWin.t_range)
        t_range = self.SpecWin.ax.get_xlim()
        if t_win > ptp(t_range)/2:
            raise Exception('Please select smaller window')
        
        
               

        self.f0 = self.OptimizeF0(self.tvec, self.data, mean(frange),
                            df0 = max(frange)-mean(frange),ch = self.ch0,n_steps=1e2)
        nmax = int(self.f0*dt*nt)

        print((t_win, dt, self.f0))
        f_bank = filter_bank(  t_win, dt, self.f0,n=self.n_harm,window=self.window)
        b = f_bank[0,:]  #basis frequency filter
        self.b_comb =  comb_filter(t_win, dt, self.f0,n=self.n_harm,window=self.window,pass_DC=True)



        cmpl_exp = zeros(2**int(ceil(log2(nt))),dtype='single')
        cmpl_exp[(nmax*len(cmpl_exp))//nt] = 1
        cmpl_exp = ifft(cmpl_exp)[:nt]*len(cmpl_exp)


        offset = mean(self.data,axis=0)
        #from matplotlib.mlab import detrend_linear
        #offset = data-detrend_linear(data,axis=0 )
        self.data -= offset[None,:] 


        n_fft = 2**int(ceil(log2(nt+len(b)-1)))
        fsig = fft(self.data,axis=0,n=n_fft)#použít RFFT? 
        fb   = fft(single(b),axis=0,n=n_fft)




        #import IPython
        #IPython.embed()
        #print self.err, offset,mean(offset)

        self.retrofit = zeros((nt-f_bank.shape[1]//2, self.ndets ),dtype=single)
        weight = 1/single(self.err/(offset+mean(offset)/100))  #weight
        #plot(offset)
        #show()
        #plot(std(diff(fsig,axis=0),axis=0)/std(fsig))

        #plot(offset/mean(offset))
        #show()N=
        #plot(mean(weight)/weight/100)
        #show()
        downsample = int(ceil(1/dt/(self.f0*2*self.n_harm))) *2
        downsample = max(downsample, nt/100)
        #retrofit_unmodulate = zeros(((nt-f_bank.shape[1]/2)/downsample+1, ndets ),dtype=single)
        #phi = ones((nt-f_bank.shape[1]/2)/downsample+1)

        #fig = self.fig_spec
        self.fig.clf()
        self.fig.subplots_adjust(hspace=0.051, wspace = 0.05,left=0.05,top=0.95,right=0.99,bottom=0.04)
        #self.fig.subplots_adjust(left=0.1,top=0.85,right=0.95,bottom=0.1)

        axes = []
        ax2 = None
        for i in range(self.n_harm):
            ax1 = self.fig.add_subplot(self.n_harm,2,2*i+1,sharex=ax2,sharey=ax2 )
            ax2 = self.fig.add_subplot(self.n_harm,2,2*i+2,sharex=ax1,sharey=ax1)
            
            axes.append((ax1,ax2))
            ax1.set_ylabel('%d. harm.'%i,fontsize=10)
            ax2.yaxis.tick_right()
            for label in ax1.get_yticklabels():
                label.set_fontsize(10) # Size here overrides font_prop
            #for label in (ax2.get_xticklabels()):
                #label.set_fontsize(10) # Size here overrides font_prop
            ax1.yaxis.set_major_formatter( NullFormatter() )
            ax2.yaxis.set_major_formatter( NullFormatter() )
            
            
            #for label in ax1.get_xticklabels():
                #label.set_visible(False)
            
            for label in ax2.get_xticklabels()+ax2.get_yticklabels()+ax1.get_xticklabels():
                label.set_visible(False)
            ax2.xaxis.offsetText.set_visible(False)
            ax1.xaxis.offsetText.set_visible(False)
            ax2.yaxis.offsetText.set_visible(False)

        for label in ax1.get_xticklabels()+ax2.get_xticklabels():
            label.set_visible(True)
            label.set_fontsize(10) # Size here overrides font_prop

        ax1.xaxis.offsetText.set_visible(True)
        ax2.xaxis.offsetText.set_visible(True)
            
        #ax1.set_xlabel('detector',fontsize=10)
        #ax2.set_xlabel('detector',fontsize=10)
        #ax1.yaxis.set_major_formatter( MyFormatter() )
        len_b = len(b)//2
        self.len_b = len_b
        axes = array(axes)
        extent=(0,self.ndets,self.tvec[0],self.tvec[-len_b])
        len_red = len_b//downsample

        cmpl_exp_i = ones_like(cmpl_exp)
        phi = ones((nt-(f_bank.shape[1]+1)//2)//downsample+1)
        t = time.time()
        from tqdm  import trange
        for i_harm in trange(self.n_harm):
        
            #filter out the choosen harmonics
            C1 = ifft(fsig*fb[:,None],axis=0,overwrite_x=1 )[:nt+len(b)-1]

            filtered_harmonic = C1[(len(C1)-nt+f_bank.shape[1])//2:nt+len(b)-1-(len(C1)-nt+1)//2]
            
            
            #import IPython
            #IPython.embed()

            #shift array to the next harmonics
            fsig = self.shift_array(fsig,(nmax*n_fft)//nt)


            #complex SVD filtration of the harmonics
            filtered_harmonic_low = copy(filtered_harmonic[:-len_red:downsample,:])
            filtered_harmonic_low*=weight[None,:]  #weight channels according their error

            U,s,V = linalg.svd(filtered_harmonic_low/phi[:,None]**i_harm,full_matrices=0)
            
            fact = 1 if i_harm == 0 else 2
            svd_filtered_harm = dot(dot(conj(V[:max(1,self.n_svd-i_harm)]*weight[None,:]),
                                        filtered_harmonic.T).T, V[:max(1,self.n_svd-i_harm)]/(weight[None,:]/fact))

            #svd_filtered_harm*=weight[None,:]
            fun = colorize

            if i_harm == 1:  #global phase of the mode
                phi =  U[:,0]/abs(U[:,0])
                

            #plotting
            vmax = mquantiles(abs(svd_filtered_harm[:-len_red:downsample,:]),0.95)[0]
            axes[i_harm,0].imshow(fun(filtered_harmonic[:-len_red:downsample,:]/phi[:,None]**i_harm/vmax),
                                    origin='lower',aspect='auto',interpolation='nearest',
                                    vmin=0,vmax=1,extent=extent)
            axes[i_harm,1].imshow(fun( svd_filtered_harm[:-len_red:downsample,:]/phi[:,None]**i_harm/vmax),
                                    origin='lower',aspect='auto',interpolation='nearest',
                                    vmin=0,vmax=1,extent=extent)
                
            axes[i_harm,0].axis('tight')
            axes[i_harm,1].axis('tight')

            for det in self.dets:
                axes[i_harm,0].axvline(x=det-.5,c='w')
                axes[i_harm,1].axvline(x=det-.5,c='w')
                
            
            #print phi.shape , svd_filtered_harm[::downsample,:].shape, retrofit_unmodulate.shape
            #retrofit_unmodulate+= real(svd_filtered_harm[::downsample,:]/(phi[:,None]**i_harm))
            #retrofit_unmodulate+= einsum('ij,i->ij',svd_filtered_harm[::downsample,:].imag,(angle**i_harm).imag)
            
            #calculate a retrofit of the original signal 
            self.retrofit+= einsum('ij,i->ij',svd_filtered_harm.real,cmpl_exp_i.real[f_bank.shape[1]//2:])
            self.retrofit+= einsum('ij,i->ij',svd_filtered_harm.imag,cmpl_exp_i.imag[f_bank.shape[1]//2:])
            #if not remove_rotation:
            cmpl_exp_i*= cmpl_exp
            #angle**i_harm
            
        print('computing time: %f'%(time.time()-t))
        #exit()
        multi = MultiCursor(self.fig.canvas, axes.flatten(),horizontal=True, color='k', lw=1)

        axes[0,0].set_title('Raw signal',fontsize=10)
        axes[0,1].set_title('SVD filtered',fontsize=10)
        self.AxZoom = fconf.AxZoom()
        self.phase_cid = self.fig.canvas.mpl_connect('button_press_event', self.AxZoom.on_click)
        #self.fig.canvas.mpl_connect('scroll_event',self.WheelInteraction)
        #ylim
        #show()
        #close()
        #close()
        #imshow(fun( retrofit_unmodulate /vmax),
                                    #origin='lower',aspect='auto',interpolation='nearest',
                                    #vmin=0,vmax=1,extent=(0,ndets, 0,nt/downsample))

        self.actual = True
        self.TSS = linalg.norm(self.data[len_b:-len_b])**2
        self.RSS = linalg.norm(self.data[len_b:-len_b]-self.retrofit[:-len_b])**2
        #data-retrofit
        
        save('data', self.data[len_b:-len_b])
        save('resid', self.data[len_b:-len_b]-self.retrofit[:-len_b])


        self.data += offset[None,:]
        self.retrofit+= offset[None,:]
        #retrofit_unmodulate+= offset[None,:]

        self.retrofit[self.retrofit<=0] = 1e-3  #make it nonzero
        
        #len_b = len(b)/2

        #correctly estimated errorbars for gaussian noise!!
        errors = std(self.data[len_b:-len_b]-self.retrofit[:-len_b],0)*linalg.norm(self.b_comb) 
        #errors = std(data[len_b:-len_b,:]-retrofit[:-len_b,:],0)*linalg.norm(b_comb)   #correctly estimated errorbars for gaussian noise!!
        errors = hypot(errors, self.data.mean()*0.005)
        
        
        #return 
    
    
    
    
    
    
    
    
            
        #data += offset[None,:]
        #retrofit+= offset[None,:]
        ##retrofit_unmodulate+= offset[None,:]

        #retrofit[retrofit<=0] = 1e-3  #make it nonzero
        

        #save filtered signals
        #if remove_rotation:
            #errors = 
        #else:
        #errors = std(data[len_b:-len_b,:]-retrofit[:-len_b,:],0)*linalg.norm(b_comb)   #correctly estimated errorbars for gaussian noise!!
        #errors = hypot(errors, data.mean()*0.005)
        #for i,d in enumerate(SXR_detectors):
            #save('results/%s_%d_err.npy'%(d,shot),single(errors[dets[i]:dets[i+1]] ))
            #save('results/'+d+"_"+str(shot)+"_stat.npy",all(abs(retrofit[:,dets[i]:dets[i+1]])>1,axis=0))
            
            ##if remove_rotation:
                ##savez('results/tvec_%d.npz'%(shot),tvec=tvec[:-len_b:downsample])
                ##save('results/%s_%d.npy'%(d,shot),single(retrofit_unmodulate[:,dets[i]:dets[i+1]] ))
                ##print retrofit_unmodulate[:,dets[i]:dets[i+1]].shape, tvec[:-len_b:downsample].shape

            ##else:
            #savez('results/tvec_%d.npz'%(shot),tvec=tvec[:-len_b])
            #save('results/%s_%d.npy'%(d,shot),single(retrofit[:,dets[i]:dets[i+1]] ))



            
        ##ang = unwrap(angle((svd_filtered_harm[:,:,1])),axis=0)
        #ang -= ang.mean(0)

        #plot(tvec[:-f_bank.shape[1]/2], ang);show()
        #plot(tvec[:-f_bank.shape[1]/2], ang);show()
        #embed()
        #plot(ang.mean(0))
            
        #for det in dets:
            #axvline(x=det-.5,c='k')
            
        #embed()
            
        #ang = unwrap(angle((svd_filtered_harm[:,:,2])),axis=0)
        #ang -= ang.mean(0)
        #plot(tvec[:-f_bank.shape[1]/2], ang);show()

        #ang = unwrap(angle((svd_filtered_harm[:,:,3])),axis=0)
        #ang -= ang.mean(0)
        #plot(tvec[:-f_bank.shape[1]/2], ang);show()

        #print 'DONE plotting'





            




        #BUG neobtížnější kanál skrz střed ostrova 


        ##f3 = figure(str(shot)+' retrofit')
        ##ax_f3 = f3.add_subplot(111)
        ##N = 71
        ###plot(retrofit2[:,N]-mean(data[:-len_b,N]),label='time filtered')
        ##ind = slice(None,None )
        ##iimin = tvec.searchsorted(tvec[0])
        ##iimax = tvec.searchsorted(tvec[-1])+1
        ###print iimin, iimax, len(tvec)
        ##ind = slice(iimin,iimax )
        ##mean_data = mean(data[ind,N])
        ###print data[ind,N][len_b:].shape, retrofit[ind,N].shape,lfilter(b_comb,1,data[ind,N]-mean_data)[len_b:].shape, tvec[ind][len_b:].shape
        ##plt_data, = ax_f3.plot(data[:,N][len_b:],label='data')
        ##plt_retro, = ax_f3.plot( retrofit[:,N] ,'--',label='svd filtered')
        ##plt_lfilt, = ax_f3.plot( lfilter(b_comb,1,data[:,N]-mean_data)[2*len_b:] +mean_data,label='comb filter')
        ##ax_f3.legend()
        ###show()


        ##def update_single_plot(tmin,tmax,N):
            ###return 
            ###print 'update_single_plot'
            ##iimin = tvec.searchsorted(tmin)
            ##iimax = tvec.searchsorted(tmax)+1
            ###ind = slice(iimin,min(iimax, len(tvec)-len_b) )
            ##ind2 = slice(max(iimin,len_b),iimax)
            ###ind3 = slice(max(iimin,len_b),iimax-len_b)

            ##mean_data = mean(data[ind,N])
            ###print tvec[ind2].shape, data[ind2,N].shape, retrofit[ind,N].shape, lfilter(b_comb,1,data[max(0,iimin-len_b):iimax,N])[len_b:].shape
            ##plt_data.set_data( tvec[iimin:iimax],data[iimin:iimax,N])
            ###plt_retro.set_data(tvec[iimin+len_b:iimax ],retrofit[ iimin:iimax-len_b,N])
            ##plt_retro.set_data(tvec[max(0,iimin-len_b)+len_b:iimax ],retrofit[ max(0,iimin-len_b):iimax-len_b,N])


            ###plt_lfilt.set_data(tvec[ind2][:len_b],lfilter(b_comb,1,data[max(0,iimin-2*len_b):iimax,N]-mean_data)[ iimin-max(0,iimin-2*len_b):]+mean_data)
            ###print tvec[iimin:iimax][len_b:-len_b].shape, lfilter(b_comb,1,data[iimin:iimax,N]-mean_data)[2*len_b:].shape
            ###plt_lfilt.set_data(tvec[iimin:iimax][len_b:-len_b],lfilter(b_comb,1,data[iimin:iimax,N]-mean_data)[2*len_b:]+mean_data)
            ##plt_lfilt.set_data(tvec[:][len_b:-len_b],lfilter(b_comb,1,data[:,N]-mean_data)[2*len_b:]+mean_data)

            ###print 
            ###plt_lfilt.set_data(tvec[max(len_b,iimin):min(len(tvec)-len_b,iimax)],lfilter(b_comb,1,data[max(0,iimin-len_b):min(iimax, len(tvec)-len_b),N]-mean_data)[2*len_b:]+mean_data)
            ##if ax_f3.get_xlim() != (tvec[ind2][0],tvec[ind2][-1]):
                ##ax_f3.set_xlim(tvec[ind2][0],tvec[ind2][-1])
            ##ax_f3.set_ylim(data[ind2,N].max(),data[ind2,N].min())
            ##ax_f3._N = N
            ###exit()

            ##f3.canvas.draw()
        ##update_single_plot(0,infty, N)
        #show()

        ##plot(lfilter(b_comb,1, data[:,71]-mean(data[:-len_b,71]))[f_bank.shape[1]/2:])
        #plot(lfilter(b_comb,1, data[:,71]-mean(data[:-len_b,71]))[f_bank.shape[1]-1:],'--')

        #show()

        #plot(data[:-len_b,72]-mean(data[:-len_b,72]))
        #plot(retrofit2[:,72]-mean(data[:-len_b,72]),'--')

        #plot(lfilter(b_comb,1, data[:,71]-mean(data[:-len_b,71]))[f_bank.shape[1]/2:])
        #show()



        #print 'differences %.2f%% %.2f%% '%(100*linalg.norm(data[len_b:,:]-retrofit)/linalg.norm(data[len_b:,:]), 100*linalg.norm(data[len_b:,:]-retrofit2)/linalg.norm(data[len_b:,:]))


        self.fig2.clf()
        
        
        f = self.fig2
        #f = figure()
        import matplotlib.gridspec as gridspec
        
        gs0 = gridspec.GridSpec(2, 1,   height_ratios=[1,5] )
        gs00 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[0])
        gs01 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs0[1])
  
        

        
        
        #ax1 = Subplot( f, gs01[0])
        #f.add_subplot(ax1)
        
        #ax2 = Subplot( f, gs01[1])
        #f.add_subplot(ax2)
        
        #ax3 = Subplot( f, gs01[2])
        #f.add_subplot(ax3)
        
        
        #ax0.xaxis.set_ticks_position('top')
        

        nplots = 3
        axarr = np.empty(nplots, dtype=object)
        ax1 = Subplot(f, gs01[0])
        f.add_subplot(ax1)

        

        axarr[0] = ax1

        for i in range(1, nplots):
            axarr[i] = Subplot(f, gs01[i],sharex=ax1,sharey=ax1)
            f.add_subplot(axarr[i])

        for ax in axarr[..., 1:].flat:
            for label in ax.get_yticklabels():
                label.set_visible(False)
            ax.yaxis.offsetText.set_visible(False)

        ax0 = Subplot( f, gs00[0])
        f.add_subplot(ax0)
        ax0.xaxis.set_ticks_position('top')
        ax0.yaxis.offsetText.set_visible(False)
        ax0.xaxis.offsetText.set_visible(False)
        for label in ax0.get_xticklabels()+ax0.get_yticklabels():
            label.set_visible(False)
            
        (ax1, ax2, ax3) = axarr

        self.retrofit[:,self.wrong_det] = 0
        vmax = 50

        
        
        
        for ax in (ax0,ax1,ax2,ax3):
            for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontsize(10) # Size here overrides font_prop
        #ax2.yaxis.set_major_formatter( NullFormatter() )
        #ax3.yaxis.set_major_formatter( NullFormatter() )
        #ax0.set_text()
        
        #TODO aby to psalo čilso signálu co by se měl vyhodit ale i kód LOS
        #signal_name = ax0.text(0.05,0.7, 'detektor',transform=ax0.transAxes,fontsize=10, backgroundcolor='w')


        f.subplots_adjust(hspace=0.1, wspace = 0.05, left=0.1,top=0.93,right=0.98,bottom=0.05)



        ax2.set_title('filtered',fontsize=10)
        extent=(0,self.ndets,self.tvec[0],self.tvec[-1])
        im1 = ax2.imshow(-(self.retrofit-self.data.mean(0))/sqrt(self.retrofit.std(0)+1),aspect='auto'
                        ,vmin=-vmax,vmax=vmax,cmap='seismic',extent=extent,
                        interpolation='nearest',origin='lower');
        ax1.set_title('original',fontsize=10)
        im2 = ax1.imshow(-(self.data[len_b:,:]-self.data.mean(0))/sqrt(self.retrofit.std(0)+1),
                        aspect='auto',vmin=-vmax,vmax=vmax,cmap='seismic',
                        extent=extent,interpolation='nearest',origin='lower') ;
        ax3.set_title('original-filtered',fontsize=10)
        im3 = ax3.imshow(-(self.data[len_b:,:]-self.retrofit)/sqrt(self.retrofit.std(0)+1),
                        aspect='auto',vmin=-vmax,vmax=vmax,cmap='seismic',
                        extent=extent,interpolation='nearest',origin='lower');
        #cbar = self.fig2.colorbar(im3) 
        #cbar.ax.tick_params(labelsize=10) 

        for ax in (ax1,ax2,ax3):
            for det in self.dets:
                ax.axvline(x=det,c='k')
                ax.axvline(x=det,ls='--',c='w')
                
        self.multi2 = MultiCursor(self.fig2.canvas, (ax1,ax2,ax3),horizontal=True, color='k', lw=1)
        self.AxZoom2 = fconf.AxZoom()
        
        def update_single_plot(tmin,tmax,N):
            iimin = self.tvec.searchsorted(tmin)
            iimax = self.tvec.searchsorted(tmax)+1
            ind2 = slice(max(iimin,len_b),iimax)
            #mean_data = mean(self.data[ind,N])
            
            self.plt_data.set_data( self.tvec[iimin:iimax],self.data[iimin:iimax,N])
            self.plt_retro.set_data(self.tvec[max(0,iimin-len_b)+len_b:iimax ],self.retrofit[ max(0,iimin-len_b):iimax-len_b,N])

            if ax0.get_xlim() != (self.tvec[ind2][0],self.tvec[ind2][-1]):
                ax0.set_xlim(self.tvec[ind2][0],self.tvec[ind2][-1])
            ax0.set_ylim(self.data[ind2,N].max(),self.data[ind2,N].min())
            ax0._N = N
            self.fig2.canvas.draw()

        def mouse_interaction(event):
            #print 'mouse_interaction'
            if not hasattr(event,'button'):
                tmin,tmax = event.get_ylim()
                if ax0.get_xlim() == event.get_ylim():
                    return 
                update_single_plot(tmin,tmax,ax0._N)
                
            elif hasattr(event,'button') and event.button is 3:
                self.AxZoom2.on_click(event)
            #elif hasattr(event,'button') and event.button is 1 and event.inaxes in [ax1,ax2,ax3]:
                #N = int(event.xdata)
                
            elif hasattr(event,'button') and event.button is 2 and event.inaxes in [ax1,ax2,ax3]:
                tmin,tmax = event.inaxes.get_ylim()
                N = int(event.xdata)

                update_single_plot(tmin,tmax,N)
            
            #if not hasattr(event,'button'):

                #tmin,tmax = event.get_ylim()
                #xmin,xmax = event.get_xlim()
                #axis = [xmin,xmax,tmin,tmax]
                #if axis != ax1.axis(): ax1.axis([xmin,xmax,tmin,tmax])
                #if axis != ax2.axis(): ax2.axis([xmin,xmax,tmin,tmax])
                #if axis != ax3.axis(): ax3.axis([xmin,xmax,tmin,tmax])
                #self.fig2.canvas.draw()

                #update_single_plot(tmin,tmax,N)
        self.plt_data,  = ax0.plot([],[],label='data')
        self.plt_retro, = ax0.plot([],[] ,'--',label='svd filtered') 
        ax0.set_xlim(self.tvec[len_b],self.tvec[-1])
        
        def mouse_interaction2(event):
            tmin,tmax = event.get_xlim()
            ax1.set_ylim(tmin,tmax)
            ax2.set_ylim(tmin,tmax)
            ax3.set_ylim(tmin,tmax)
            if tmin < self.plt_data.get_xdata()[0] or tmax > self.plt_data.get_xdata()[-1]:
                update_single_plot(tmin,tmax,ax0._N)
            #if tmin < self.plt_data.get_xdata()[0] or tmax > self.plt_data.get_xdata()[-1]:
                #pass
            self.fig2.canvas.draw()
            
        update_single_plot(0,infty, self.ch0)


      
        self.cid0 = self.fig2.canvas.mpl_connect('button_press_event', mouse_interaction)
        #ax_f3.callbacks.connect('xlim_changed',  mouse_interaction2)
        self.cid1 =  ax0.callbacks.connect('xlim_changed',  mouse_interaction2)

        ax1.callbacks.connect('ylim_changed',  mouse_interaction)
        ax2.callbacks.connect('ylim_changed',  mouse_interaction)
        ax3.callbacks.connect('ylim_changed',  mouse_interaction)
        #ax1.callbacks.connect('xlim_changed',  mouse_interaction)
        #ax2.callbacks.connect('xlim_changed',  mouse_interaction)
        #ax3.callbacks.connect('xlim_changed',  mouse_interactio
        
        
        
        #f3 = figure(str(shot)+' retrofit')
        #ax_f3 = f3.add_subplot(111)
        #plot(retrofit2[:,N]-mean(data[:-len_b,N]),label='time filtered')
        #ind = slice(None,None )
        #iimin = tvec.searchsorted(tvec[0])
        #iimax = tvec.searchsorted(tvec[-1])+1
        #print iimin, iimax, len(tvec)
        #ind = slice(iimin,iimax )
        #mean_data = mean(data[ind,N])
        #print data[ind,N][len_b:].shape, retrofit[ind,N].shape,lfilter(b_comb,1,data[ind,N]-mean_data)[len_b:].shape, tvec[ind][len_b:].shape
                #N = 71
    #self.ch0

        #plt_lfilt, = ax_f3.plot( lfilter(b_comb,1,data[:,N]-mean_data)[2*len_b:] +mean_data,label='comb filter')
        #ax_f3.legend()
        
        
        
        
        
        
        



    def Save(self):
        
        len_b = self.len_b
        #correctly estimated errorbars for gaussian noise!!
        errors = std(self.data[len_b:-len_b]-self.retrofit[:-len_b],0)*linalg.norm(self.b_comb) 
        errors = hypot(errors, self.data.mean()*0.005)
        

        #save filtered signals

        for i,d in enumerate(self.SXR_detectors):
            save('results/%s_%d_err.npy'%(d,self.shot),single(errors[self.dets[i]:self.dets[i+1]] ))
            save('results/'+d+"_"+str(self.shot)+"_stat.npy",all(abs(self.retrofit[:,self.dets[i]:self.dets[i+1]])>1,axis=0))


            savez('results/tvec_%d.npz'%(self.shot),tvec=self.tvec[:-self.len_b])
            save('results/%s_%d.npy'%(d,self.shot),single(self.retrofit[:,self.dets[i]:self.dets[i+1]] ))



        print('data saved')
    



    def OptimizeF0(self, tvec, sig, f0,df0 = 200,ch=None,n_steps=100):
        
        #simple and naive algorithm to fine tune basic frequency
        ind = slice(None,None) if ch is None else slice(ch,ch+1)

        sig = sig[:,ind]-sig.mean(0)[None,ind]
        difference = zeros(n_steps)
        
        test_fun = exp(1j*2*pi*(f0-df0)*tvec)
        test_fun/= linalg.norm(test_fun)
        dtest_fun = exp(1j*2*pi*(2*df0/float(n_steps))*tvec)
        
        for i in arange(n_steps):

            retro = outer(test_fun, dot(conj(test_fun), sig))
            difference[i] = linalg.norm(retro.real-sig)
            
            test_fun*= dtest_fun
            
        #fig = figure()
        ##fig.clf()
        #ax = fig.add_subplot(111)
        #ax.plot(difference)
        #ax.plot(linspace(-df0,df0,n_steps)+f0,difference/linalg.norm(sig))
        #ax.axvline(x=f0)
        #fig.savefig('f0.pdf')
        #fig.clf()
        savetxt('data.txt', c_[linspace(-df0,df0,n_steps)+f0,difference/linalg.norm(sig)])
            
        print('F0  %f %f %d'%( f0+(argmin(difference)*2.-n_steps+1)/(n_steps)*df0, f0,ch))
            
        return f0+(argmin(difference)*2.-n_steps+1)/(n_steps)*df0
    

    def shift_array(self,x,shift):
        #roll along first axis - much faster than numpy.roll
        if shift < 0:
            raise Exception('not implemented yet')
    
        tmp = copy(x[-shift:])
        x[shift:] = x[:-shift]
        x[:shift] = tmp
        return x







class SpectraViewer(object):
    def __init__(self,sgamma,snfft,method,show_raw=True,allow_selector=True,fig_name=None,fig=None):
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        
        #matplotlib.rcParams['xtick.direction'] = 'out'
        #matplotlib.rcParams['ytick.direction'] = 'out'

        matplotlib.rcParams['xtick.major.size'] = 4
        matplotlib.rcParams['xtick.minor.size'] = 2

        matplotlib.rcParams['ytick.major.size'] = 4
        matplotlib.rcParams['ytick.minor.size'] = 2


    

        axcolor = 'lightgoldenrodyellow'
        self.show_raw = show_raw
        if fig is None:
            self.fig, self.ax = plt.subplots(num=fig_name)
        else:
            self.fig, self.ax = fig, fig.gca()
        #if show_raw:
        self.fig.subplots_adjust(left=0.11,top=0.85,right=0.98,bottom=0.1)
        self.uax = self.fig.add_axes([0.11, 0.85, 0.87, 0.13])
        self.uax.xaxis.set_major_formatter(NullFormatter())
        self.uax.yaxis.set_major_formatter(NullFormatter())
        #else:
            #self.fig.subplots_adjust(left=0.2,top=0.95,right=0.95,bottom=0.1)
        self.font_size = 10
        #self.rax = self.fig.add_axes([0.05, 0.7, 0.1, 0.15], axisbg=axcolor)
        self.ax_gamma = None# self.fig.add_axes([0.05, 0.2, 0.1, 0.05], axisbg=axcolor)
        self.ax_nfft = None#self.fig.add_axes([0.05, 0.1, 0.1, 0.05], axisbg=axcolor)
        print('allow_selector %s'%allow_selector)

        self.sdpi = 1000.

        self.f0_plot, = self.ax.plot([],[],'k',zorder=99)
        self.ax.set_xlabel('Time [s]',fontsize=self.font_size)
        self.ax.set_ylabel('Frequency [kHz]',fontsize=self.font_size)
        
        #import matplotlib.font_manager as fm
        #font = fm.FontProperties(fname=fontPath, size=10)
        #self.ax.xaxis.get_label().set_fontproperties(font)
        #self.ax.yaxis.get_label().set_fontproperties(font)
        
        for label in (self.ax.get_xticklabels() + self.ax.get_yticklabels()):
            label.set_fontsize(self.font_size) # Size here overrides font_prop


        self.methods = ('time', 'freq.','sparse')
        #i_radio = [i for i in range(3) if methods[i]==stft_disp.method][0]
        #self.radio = RadioButtons(self.rax, methods,active=0)
        
        self.sgamma  = sgamma
        self.snfft  = snfft
        self.method = method
        
        sgamma.setRange(0.01*self.sdpi, 1*self.sdpi)
        sgamma.setValue(.5*self.sdpi)
        sgamma.setTracking(True)
        sgamma.setTickPosition(QSlider.NoTicks)
        
        #self.sgamma = Slider(self.ax_gamma, 'Gamma',0.01*self.sdpi,
                             #1*self.sdpi, valinit=.5,valfmt=u'%0.1f')
        #self.snfft = Slider(self.ax_nfft, 'NFFT',  6*self.sdpi, 
                            #16*self.sdpi, valinit=13*self.sdpi,valfmt=u'%1.0f')

        
        snfft.setRange(6*self.sdpi, 16*self.sdpi)
        snfft.setValue(13*self.sdpi)
        snfft.setTracking(True)
        snfft.setTickPosition(QSlider.NoTicks)
        
        
        self.method.addItem("STFT")
        self.method.addItem("SFFT")
        self.method.addItem("Sparse")
        self.method.setCurrentIndex(0)
        
        
        sgamma.setSingleStep(.1)
        snfft.setSingleStep(1)

        #method.currentIndex()
                #self.method  = QComboBox(self.main_frame)

        #self.sgamma.label.set_position((.9,1.4))
        #self.snfft.label.set_position((0.7,1.4))
        #self.sgamma.valtext.set_position((0.35,.5))
        #self.snfft.valtext.set_position((0.25,.5))
        #self.sgamma.poly.set_facecolor('lightblue')
        #self.snfft.poly.set_facecolor('lightblue')
        

        matplotlib.rcParams['xtick.direction'] = 'in'
        matplotlib.rcParams['ytick.direction'] = 'in'

        if allow_selector:
            rectprops = dict(facecolor='gray', edgecolor = 'black',alpha=0.5, fill=True,zorder=99)
            self.RS1 = RectangleSelector(self.ax, self.line_select_callback,
                                        drawtype='box', useblit=True,
                                        button=[1,], # don't use middle button
                                       minspanx=5, minspany=5,rectprops=rectprops,
                                       spancoords='pixels')
        

            
    def init_plot(self, tvec,signal,window = 'gauss',tmin=None,tmax=None,
                  fmin0=1e3,fmax0=2e5,cmap= 'gnuplot2' ):
        
        
        #object of the spectrogram
        
        self.stft_img =  STFTImage(self.ax,gamma=.5,cmap=cmap)
        if tmin is None: tmin =  tvec[0]
        if tmax is None: tmax =  tvec[-1]

        self.t_range = tmin,tmax
        self.f_range = fmin0,fmax0
        #print self.t_range, self.f_range,'elf.t_range, self.f_range'
        #clear the image
        for im in self.ax.get_images():
            im.remove()
            del im

        
        if self.show_raw:
            
            #clear it before use
            for l in self.uax.get_lines():
                l.remove()
                del l
            
            #object of the plot with raw data
            self.data_plot = DataPlot(self.uax,lw=.5, c='k')
            self.data_plot.prepare(tvec,signal )
            self.data_plot.update_plot( tmin,tmax)

            #object of the both + STFT algorithms 
            self.stft_disp = STFTDisplay(signal,tvec,self.stft_img,self.ax,
                                        self.data_plot,win=window)
            self.cid3 = self.uax.callbacks.connect('xlim_changed', 
                                            self.stft_disp.plot_ax_update)

        else:
            self.stft_disp = STFTDisplay(signal,tvec,self.stft_img,self.ax, 
                                        win=window)
    
        self.window = window
        #print fmin0, fmax0, 'initplot', tmin, tmax
        x,y,Z = self.stft_disp(tmin, tmax, fmin0, fmax0)
        self.stft_img.prepare(x,y,Z)
        self.ax.set_ylim(y[0],y[-1])#;print 'set_ylim(y[0],y[-1])'
        self.ax.set_xlim(x[0],x[-1])#;print 'set_ylim(x[0],x[-1])'

        #print 'self.ax.set_ylim',y[0],y[-1]

        self.cid1 = self.ax.callbacks.connect('xlim_changed', self.stft_disp.ax_update)
        self.cid2 = self.ax.callbacks.connect('ylim_changed', self.stft_disp.ax_update)
        #does not work with Qt4!!
        self.cid4 = self.fig.canvas.mpl_connect('resize_event', self.stft_disp.ax_update)
        #self.radio.on_clicked(self.select_DFT_backend)
        #self.sgamma.on_changed(self.apply_slider)
        
        #majorLocator   = MaxNLocator(6)
        #minorLocator   = MaxNLocator(30)

        #self.ax.xaxis.set_major_locator(majorLocator)
        self.ax.xaxis.set_minor_locator(MaxNLocator(30))
        self.ax.yaxis.set_major_locator(MaxNLocator(6))
        self.ax.yaxis.set_minor_locator(MaxNLocator(30))

        #self.snfft.on_changed(self.apply_slider)
        self.fig.canvas.mpl_connect('scroll_event',self.WheelInteraction)
        self.apply_slider(0)

        
    def line_select_callback(self,eclick, erelease):
        'eclick and erelease are the press and release events'
        print('line_select_callback')
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        self.t_range = x1,x2
        self.f_range = y1,y2
        
        self.RS1.to_draw.set_visible(True)
        self.RS1.canvas.draw()

    def select_DFT_backend(self,label):
        
        self.stft_disp.method = self.methods[label]
        
        xstart, xend, ystart, yend = self.stft_img.ax.axis()
        x,y,Z = self.stft_disp(xstart, xend, ystart, yend)
        #print ystart, yend, y.min(), y.max(),'select_DFT_backend'
        #print xstart, xend, x.min(), x.max()
        #print

        self.stft_img.z = Z
        self.stft_img.x = x
        self.stft_img.y = y
        self.stft_img.im.set_extent((x[0],x[-1],y[0],y[-1]))


        self.apply_slider(0)
        
        
    def updateGammaSlider(self,val):
        self.sgamma.setValue(val*self.sdpi)
        #self.sgamma.value = val
        #poly = self.sgamma.poly.get_xy()
        #poly[2:4,0] = val
        #self.sgamma.poly.set_xy(poly)  
        #self.sgamma.valtext.set_text(self.sgamma.valfmt%val)

    def updateNFFTSlider(self,val):
        self.snfft.setValue(val*self.sdpi)

        #self.snfft.value = val
        #poly = self.snfft.poly.get_xy()
        #poly[2:4,0] = val
        #self.snfft.poly.set_xy(poly) 
        #self.snfft.valtext.set_text(self.snfft.valfmt%(2**val))


    def apply_slider(self,val):
        gamma = self.sgamma.value()/self.sdpi
        nfft = self.snfft.value()/self.sdpi
        nfft = 2**nfft if self.window == 'gauss' else 2**int(nfft)
            
        xstart, xend, ystart, yend = self.stft_img.im.get_extent()
        if nfft!= self.stft_disp.nfft:
            self.stft_disp.nfft = nfft
            #print ystart, yend, 'apply_slider'

            x,y,Z = self.stft_disp(xstart, xend, ystart, yend)
            self.stft_img.z = Z
        if gamma!= self.stft_img.gamma:
            self.stft_img.gamma = gamma
                
        #print ystart, yend
        self.stft_img.update_image(xstart, xend, ystart, yend,self.stft_img.z)     
        #self.snfft.valtext.set_text(self.snfft.valfmt%( 2**int(self.snfft.value())))

    def set_f0(self,f0):
        print('set_f0',f0,self.t_range,'PROC to nefunguje??')
        self.f0_plot.set_data(self.t_range,[f0/1000,f0/1000])
        self.ax.plot(self.t_range,[f0/1000,f0/1000],'k')

        self.ax.figure.canvas.draw_idle()

    def WheelInteraction(self,event):
        ##print('WheelInteraction',event.inaxes)
        #interaction with mouse wheel 
        if event.inaxes == self.ax_gamma:
            #set gamma of the spectrogram
            new_val = event.step/10.+self.sgamma.value()/self.sdpi
            new_val = max(self.sgamma.minimum/self.sdpi,min(self.sgamma.valmax/self.sdpi,new_val))

            self.updateGammaSlider(new_val)
            self.sgamma.setValue()
        
        if event.inaxes == self.ax_nfft:
            #set ratio between time/frequency resolution
            new_val = event.step+int(self.snfft.value()/self.sdpi)
            new_val = max(self.snfft.minimum/self.sdpi,min(self.snfft.maximum/self.sdpi,new_val))

            self.updateNFFTSlider(new_val)
        
        if event.inaxes == self.ax:
            #  zoom by mouse wheel
            factor = 0.9
            curr_xlim = self.ax.get_xlim()
            curr_ylim = self.ax.get_ylim()
            new_width = (curr_xlim[1]-curr_xlim[0])*factor**event.step
            new_height= (curr_ylim[1]-curr_ylim[0])*factor**event.step

            relx = (curr_xlim[1]-event.xdata)/(curr_xlim[1]-curr_xlim[0])
            rely = (curr_ylim[1]-event.ydata)/(curr_ylim[1]-curr_ylim[0])

            self.ax.set_ylim([event.ypdate_imadata-new_height*(1-rely),event.ydata+new_height*(rely)])
            self.ax.set_xlim([event.xdata- new_width*(1-relx),event.xdata+ new_width*(relx)])

        if self.show_raw and event.inaxes == self.uax:
            #  zoom by mouse wheel in upper plot 
            factor = 0.9
            curr_xlim = self.uax.get_xlim()
            new_width = (curr_xlim[1]-curr_xlim[0])*factor**event.step

            relx = (curr_xlim[1]-event.xdata)/(curr_xlim[1]-curr_xlim[0])
            xstart = event.xdata- new_width*(1-relx)
            xend = event.xdata+ new_width*(relx)
            self.ax.set_xlim([xstart,xend])

        self.apply_slider(0)
        
    def show(self):
        self.fig.show()

class STFTImage():
    def __init__(self,ax,gamma=.5,cmap='gnuplot2'):
        self.ax = ax
        self.gamma = gamma
        self.cmap = cmap

    def prepare(self,x,y,z):
        #show spectrogram
        
        #print x,y
        self.x = x
        self.y = y
        self.z = np.abs(z)
        #TODO use PCOLORMESH? 
        z = self.GammaTransform(self.z)
        self.im = self.ax.imshow(z, origin='lower', extent=(self.x.min(), 
            self.x.max(), self.y.min(), self.y.max()),aspect='auto',cmap=self.cmap)
        n_harm = 10
        #tvec,fmax = find_base_freq(x[0],x[-1],y[0],y[-1],self.z)
        #self.harm_plots = self.ax.plot(tvec,(fmax*np.arange(1,n_harm)[:,None]).T,'0.5')

        yticks = self.ax.get_yticks()

        #plot yticks in kHz
        ndig = max(0,int(np.ceil(-np.log10(np.mean(np.diff(yticks)))))+3)
        self.ax.set_yticklabels([('%.'+str(ndig)+'f')%f for f in self.ax.get_yticks()/1000]) 
       
    def GammaTransform(self,z,log=True):
        if log:
            #transformation by log
            return np.log(1+z/np.tan(np.pi/2.00001*self.gamma)/z.mean())
        else:#standart definition of the gamma
            return z**self.gamma

        
        
    def update_image(self,xstart, xend, ystart, yend,z):
        #update spectrogram image
        self.z = z
        z = self.GammaTransform(self.z)
        
        
        np.savez('spectrogram',spec=z,tmin=xstart,tend=xend, fmin=ystart,fmax= yend)
        self.im.set_data(z)
        self.im.set_extent((xstart, xend, ystart, yend))
        #print('update_image', self.im.get_extent())
        #tvec,fmax = find_base_freq(xstart, xend, ystart, yend,self.z)
        #for i,plot in enumerate(self.harm_plots):
            #plot.set_data(tvec, fmax*(i+1))
        #self.harm_plots
        
        self.ax.autoscale_view(tight=True)
        self.im.autoscale()
        self.ax.figure.canvas.draw_idle()
        
class DataPlot():
    def __init__(self,ax=None,**kwarg):
        self.kwarg = kwarg
        self.ax = ax

    def prepare(self,x,y):
        self.x = x
        self.y = y
        self.data_plot, = self.ax.plot([],[],antialiased=True, **self.kwarg)
        self.update_plot(x[0],x[-1])
        self.xstart = x[0]
        self.xend = x[-1]

        
    def update_plot(self,xstart, xend):
        print('update_plot')
        dims = self.ax.axesPatch.get_window_extent().bounds
        width = int(dims[2] + 0.5)
        
        xlim =  self.ax.get_xlim()
        self.xstart = max(self.x[0],xstart)
        self.xend = min(self.x[-1],xend)
        istart = self.x.searchsorted(self.xstart)
        iend   = self.x.searchsorted(self.xend)
        
        #random sampling
        #ind = np.sort(np.random.randint( istart, iend,width*20))
        ind = np.linspace(istart, iend,width*20)+np.random.rand(width*20)-.5
        #print ind
        ind = np.int_(np.unique(np.round(ind)))
        
        
        x = self.x[ind]
        y = self.y[ind]

        self.data_plot.set_data(x,y)
        self.ax.set_xlim(self.xstart, self.xend)
        self.ax.set_ylim(y.min(), y.max())




class STFTDisplay():
    def __init__(self,signal,tvec,image,im_ax,data_plot=None,
                 method='time', h=500, w=500,  nfft=2**13,win='gauss'):
        self.height = h
        self.width = w
        self.nfft = nfft
        self.signal = signal
        self.tvec = tvec
        self.x = None
        self.y = None
        self.image = image
        self.data_plot = data_plot
        self.im_ax = im_ax
        self.method = method
        self.dt = (tvec[-1]-tvec[0])/(len(tvec)-1)
        self.window = win
    

    def __call__(self, xstart, xend, ystart, yend):
        #print xstart, xend, ystart, yend
        #compute the signal transformation
        #BUG during zooming/unzooming is spectrogram computed two times, 
        #because of change of xlim and than chynge of ylim. 
        #import IPython
        #IPython.embed()
        import time
        t = time.time()
        #print self.method
        #print self.method 
        if self.method == 'time':
            A,fv,tv = stft(self.tvec, self.signal, self.nfft,resolution=self.width,
                    dt=self.dt,window=self.window,tmin=xstart,tmax=xend,fmin=ystart, 
                    fmax=yend,pass_DC=False)
        elif self.method == 'freq.':
            A,fv,tv = sfft(self.tvec, self.signal, self.nfft,resolution=self.height,
                        window=self.window,dt=self.dt,fmin=ystart,fmax=yend,tmin=xstart, 
                        tmax=xend,pass_DC=False)
        elif self.method == 'sparse':
            A,fv,tv = sstft(self.tvec, self.signal, self.nfft,xstart,xend,ystart,
                           yend,self.dt,self.width,self.height,zoom=4)

  
        print('stft time:%f'%( time.time()-t))
        self.x = tv
        self.y = fv

        return tv,fv,A.T

    def ax_update(self, ax):
        print('ax_update')
        #called when xlim or ylim was changed
        
        ax = self.im_ax

        ax.set_autoscale_on(False) # Otherwise, infinite loop

        #Get the number of points from the number of pixels in the window
        dims = ax.axesPatch.get_window_extent().bounds
        self.width = int(dims[2] + 0.5)
        self.height = int(dims[2] + 0.5)

        #Get the range for the new area
        xstart,ystart,xdelta,ydelta = ax.viewLim.bounds
        xend = xstart + xdelta
        yend = ystart + ydelta
     
        #set margins given by time and frequency
        if xstart < self.tvec[0] or xend > self.tvec[-1]:
            xstart = max(xstart, self.tvec[0] )
            xend   = min(xend  , self.tvec[-1])
            ax.set_xlim(xstart,xend)
            
        if ystart < 0 or yend > .5/self.dt:
            ystart = max(ystart, 0 )
            yend   = min(yend  , .5/self.dt)
            ax.set_ylim(ystart,yend)

        yticks = ax.get_yticks()
        
        #plot yticks in kHz
        ndig = max(0,int(np.ceil(-np.log10(np.mean(np.diff(yticks)))))+3)
        ax.set_yticklabels([('%.'+str(ndig)+'f')%f for f in ax.get_yticks()/1000]) 

        # Update the image object with our new data and extent
        extent = self.image.im.get_extent()
        if not self.data_plot is None and (self.data_plot.xstart != xstart or self.data_plot.xend != xend):
            #print self.data_plot.ax
            #print self.data_plot.ax.get_xlim() != (xstart, xend),self.data_plot.ax.get_xlim() , (xstart, xend)

            self.data_plot.update_plot(xstart, xend)
        #skip when nothing importnat was changed
        if self.method == 'time' and ystart >= extent[2] and \
             yend <= extent[3] and xstart==extent[0] and xend==extent[1]:
            return 
  
        if self.method == 'freq.' and ystart == extent[2] and \
             yend == extent[3] and xstart>=extent[0] and xend<=extent[1]:
            return 
        
        if self.method == 'sparse' and ystart == extent[2] and \
             yend == extent[3] and xstart==extent[0] and xend==extent[1]:
            return 

        x,y,z = self.__call__(xstart, xend, ystart, yend)
        self.image.update_image(xstart, xend, ystart, yend,z)
        print('ax_update done')

        
    def plot_ax_update(self, ax):
        xstart,xend = ax.get_xlim()
        self.im_ax.set_xlim(xstart,xend)
        
        





class DataSettingWindow(QMainWindow):
    def __init__(self, parent):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Data preprocessing tool')
        #self.ch0 = 72
        #self.shot = 31113
        #self.shot = 29222
        #self.shot = 31533

        #jak načítat data? povolit načtení celého signálu pro jeden výboj? 
        self.parent = parent
        #self.SXR_detectors=('G','H1','H2','H3','I1','I2','I3','J1','J2','J3','K1','K2','L','M')
        self.det_num = []

        #for det in self.SXR_detectors:
            #path = './geometry/ASDEX/SXR_fast/'
            #N = load(path+'%s_%d.npy'%(det,self.shot), mmap_mode='r').shape[1]
            #self.det_num.append(N)
            
        self.dpi = 100

        self.load_data()

        self.create_menu()
        self.create_status_bar()

        self.create_main_frame()
        self.create_data_table()
        #self.create_spec_table()
        #self.create_SVD_table()
        #self.create_residuum_table()
        self.main_tab.setCurrentIndex(0)
        #self.main_tab.setTabEnabled(0,False)
        
        
        
                #BUG 



        
        #self.textbox.setText('1 2 3 4')
        #self.on_draw()
    def load_data(self):
        return 
        #BUG vyřešit pořádně!!
        #self.tvec = np.load('data/tvec_%d.npz'%(self.shot))['tvec']
        #self.signal = np.load('data/%s_%d.npy'%('I2',self.shot), mmap_mode='r')[:,8]#BUG
        
        ch = self.ch0
        
        for idet, det in enumerate(self.SXR_detectors):
            if ch - self.det_num[idet] < 0:
                break
            ch-= self.det_num[idet]
            
        
                
                
        self.tvec   = np.load('signal_analysis/data/tvec_%d.npz'%(self.shot))['tvec']
        self.signal = np.load('signal_analysis/data/%s_%d.npy'%( det,self.shot), mmap_mode='r')[:,ch]
        self.sig_name = det+' %.3d'%(ch+1)
        #print det,  ch+1, self.ch0
        #exit()

        
    def save_plot(self):
        file_choices =  "PDF (*.pdf)|*.pdf"
        
        path = str(QFileDialog.getSaveFileName(self, 
                        'Save file', '', file_choices))
        
        canvas = None
        #print self.main_tab.currentIndex()
        if self.main_tab.currentIndex() == 0:
            canvas = self.canvas_data
   
        if self.main_tab.currentIndex() == 1:
            canvas = self.canvas_spec
        
        if self.main_tab.currentIndex() == 2:
            canvas = self.canvas_svd

        
        if path:
            canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
    
    def save_data(self):
        file_choices =  "NPZ (*.npz)|*.npz"
        
        path = str(QFileDialog.getSaveFileName(self, 
                        'Save file', '', file_choices))
        
        #canvas = None
        #print self.main_tab.currentIndex()
        #if self.main_tab.currentIndex() == 0:
            #canvas = self.canvas_data
   
        #if self.main_tab.currentIndex() == 1:
            #canvas = self.canvas
        
        #if self.main_tab.currentIndex() == 2:
            #canvas = self.canvas_svd

        
        #if path:
            #canvas.print_figure(path, dpi=self.dpi)
        try:
            savez(path, tvec = self.SVDF.tvec, data = self.SVDF.data, 
                        filtered=self.SVDF.retrofit)
            #TODO musí se to uložit aby to bylo správně posunuté 
            self.statusBar().showMessage('Saved to %s' % path, 2000)
        except:
            pass
    
    
    
    def on_about(self):
        msg = """ This is a tool for the preparation of the input signals for the tomography. 
       
        Basic data processing 
        *** Selecting proper time interval
        *** Averadging by box car method and downsampling
        *** Removing of the corupted channels
        
        Advaced data processing 
        *** Tool for design of the FIR filter
        *** Compressed coding by complex SVD
        *** Tool for verification of the results
        """
        QMessageBox.about(self, "About the demo", msg.strip())
    
    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        # 
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        
        QMessageBox.information(self, "Click!", msg)
    
    #def on_draw(self):
        #pass
        #""" Redraws the figure
        #"""
        ##str = unicode(self.textbox.text())
        ##self.data = map(int, str.split())
        
        ##x = range(len(self.data))

        ## clear the axes and redraw the plot anew
        ##
        ##self.image.set_data(np.random.randn(100,100))
        ##self.axes.clear()        
        ##self.axes.grid(self.grid_cb.isChecked())
        
        ##self.axes.bar(
            ##left=x, 
            ##height=self.data, 
            ##width=self.slider.value() / 100.0, 
            ##align='center', 
            ##alpha=0.44,
            ##picker=5)
        
        #self.canvas.draw()
    
    def create_main_frame(self):
        
        
        
        


        self.cWidget = QWidget(self)

        self.gridLayout = QGridLayout(self.cWidget)

        self.setWindowTitle('Data selection / Preprocessing')

        self.Expand = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.Fixed  = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        self.main_tab = QTabWidget(self)
        self.main_tab.setTabPosition(QTabWidget.North)

        self.gridLayout.addWidget(self.main_tab, 0, 0, 1,2)
        
        self.setCentralWidget(self.cWidget)
        
        
        self.resize(600, 650)
        self.setCenter()
        
        #self.connect(self.main_tab, SIGNAL('currentChanged(int)'),self.switch_tab)
        self.main_tab.currentChanged.connect(self.switch_tab)
    def onpick_data(self,event):
        #mouse inetraction, remove corrupted channels
        thisline = event.artist
        ax = event.mouseevent.inaxes

        if hasattr( thisline,'get_xdata') and ax is not None:
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            ch = xdata[ind][0]
            ax.point_label.set_position((ch+.5,ydata[ind][0]))
            name = self.LOS_names[ch-1]
            ax.point_label.set_text(name)
        else:
            ch = int(round(event.mouseevent.xdata))

        self.canvas_data.draw()
        
        if event.mouseevent.dblclick:
            wrong = list(self.ReadWrongDets())
            #print 'wrong', wrong
            
            if event.mouseevent.button == 1: #delete
                wrong.append(ch-1)
            
            if event.mouseevent.button == 3 and ch-1 in wrong: #delete
                 wrong.remove(ch-1)

            self.SetWrongDets(wrong)

            self.RefreshEvent()

    def create_data_table(self):

        #print 'create_data_table',self.parent.tokamak.index
        self.tokamak = self.parent.tokamak
        self.setting = self.parent.setting
        #self.ch_label
        
        self.tab_widget_data = QWidget(self)
        
        #self.parent = tok()

        self.main_tab.insertTab(0,self.tab_widget_data, 'Data')

        self.tab_widget_data.setSizePolicy(self.Expand)
        self.tab_widget_data.setMinimumSize(700,700)

        self.verticalLayout_data = QVBoxLayout(self.tab_widget_data)
        self.horizontalLayout_data = QHBoxLayout()


        #self.centralwidget = QWidget(self)
        #self.vLayout = QVBoxLayout(self.centralwidget)


        self.dpi = 100
        self.fig_data = Figure((7.0, 5.0), dpi=self.dpi)
        self.tab_widget_data.setToolTip('Use left/right moise button to remove/return point')

        self.canvas_data = FigureCanvas(self.fig_data)
        self.canvas_data.setParent(self.cWidget)
        
        c =  self.cWidget.palette().color(QPalette.Base)
        self.fig_data.patch.set_facecolor((c.red()/255., c.green()/255., c.blue()/255.))

        self.canvas_data.setSizePolicy(self.Expand)

        #self.horizontalLayout_data.addWidget(self.canvas )
        self.verticalLayout_data.addWidget(self.canvas_data )

        
        #self.verticalLayout_data.addLayout(self.horizontalLayout_data)
        
        #fig_data

        self.LOS_names = hstack([i for k,i in  list(self.tokamak.detectors_dict.items())])

        self.link = self.fig_data.canvas.mpl_connect('pick_event', self.onpick_data)
        self.wrong_dets_mouse = []

        #prepare first data
        self.groupBox = QGroupBox("")
        self.gridLayout = QGridLayout(self.groupBox)
        self.labelDetectors = QLabel( "Wrong detectors:")
        self.Edit_wrongDetectors = QLineEdit()
        self.refreshButton=QPushButton("Refresh")
        self.gridLayout.addWidget(self.labelDetectors, 0, 0)
        self.gridLayout.addWidget(self.Edit_wrongDetectors, 0, 1)
        self.gridLayout.addWidget(self.refreshButton, 0, 2)
        self.Edit_wrongDetectors.editingFinished.connect(self.RefreshEvent)

        self.verticalLayout_data.addWidget(self.groupBox)

        #norm = parent.tokamak.norm
        t_name = self.parent.t_name
        #parent.tokamak.t_name


        self.labelInterval = QLabel( "Time range ["+t_name+"]: ")
        self.labelFrom= QLabel( "from")
        self.tmin_spin = QDoubleSpinBox()
        self.tmin_spin.setKeyboardTracking(False)

        self.labelTo = QLabel( "to")

        self.tmax_spin = QDoubleSpinBox()
        self.tmax_spin.setKeyboardTracking(False)

        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, 
                                       QSizePolicy.Minimum)

        self.OkButton=QPushButton("OK")
        #self.tvec = self.tokamak.tvec

        self.hLayout = QHBoxLayout()
        self.hLayout.addWidget(self.labelFrom)
        self.hLayout.addWidget(self.tmin_spin)
        self.hLayout.addWidget(self.labelTo)
        self.hLayout.addWidget(self.tmax_spin)
        self.hLayout.addItem(spacerItem)
        
        self.gridLayout.addWidget(self.labelInterval, 1, 0)
        self.gridLayout.addLayout(self.hLayout, 1, 1)
        self.gridLayout.addWidget(self.OkButton, 1, 2)

        self.labelSmooth = QLabel( "Smoothing over")
        self.labelSmooth2 = QLabel( "bins")
        self.smoothSpin = QSpinBox()
        self.smoothSpin.setMaximum(999)
        self.smoothSpin.setMinimum(1)
        self.smoothSpin.setValue(self.parent.data_smooth)

        self.labelUndersampling = QLabel( "   Reconstruct every")
        self.labelUndersampling2 = QLabel( "-th snapshot")
        self.UndersamplingSpin = QSpinBox()
        self.UndersamplingSpin.setMaximum(999)
        self.UndersamplingSpin.setMinimum(1)
        self.UndersamplingSpin.setValue(self.parent.data_undersampling)

        self.gridLayout.addWidget(self.labelSmooth, 2, 0)

        #self.labelSmooth = QLabel( )
        self.check_plot_2D = QCheckBox("Use 2D plot")
        self.gridLayout.addWidget(self.check_plot_2D, 3, 0)
        self.check_plot_2D.setChecked(True)
        
        #self.data_error_averadge = QCheckBox("Mean error")
        #self.gridLayout.addWidget(self.data_error_averadge, 4, 0)
        #self.data_error_averadge.setChecked(self.parent.error_averadge)
        

        
        

        self.hLayout2 = QHBoxLayout()
        self.hLayout2.addWidget(self.smoothSpin)
        self.hLayout2.addWidget(self.labelSmooth2)
        self.hLayout2.addItem(spacerItem)
        self.hLayout2.addWidget(self.labelUndersampling)
        self.hLayout2.addWidget(self.UndersamplingSpin)
        self.hLayout2.addWidget(self.labelUndersampling2)
        self.hLayout2.addItem(spacerItem)

        self.gridLayout.addLayout(self.hLayout2, 2, 1)

        #tooltips
        self.Edit_wrongDetectors.setToolTip('Insert wrong detectors, use ":" for the intervals - 1,,3,5:40,...')
        self.smoothSpin.setToolTip('Moving average will be over [i-smooth ; i+smooth]')
        self.UndersamplingSpin.setToolTip('Undersampling to every n-th snapshot')

        #self.setCentralWidget(self.centralwidget)
        #self.resize(500, 500)
        self.resize(400, 400)

        frect = QDesktopWidget.frameGeometry(self)
        frect.moveCenter(QDesktopWidget().availableGeometry(self.parent).center());
        self.move(frect.topLeft())

        self.tmin_spin.setRange(self.tokamak.min_tvec, self.tokamak.max_tvec)
        self.tmax_spin.setRange(self.tokamak.min_tvec, self.tokamak.max_tvec)
        self.parent.tmin_spin.setRange(self.tokamak.min_tvec, self.tokamak.max_tvec)
        self.parent.tmax_spin.setRange(self.tokamak.min_tvec, self.tokamak.max_tvec)
        self.tmin_spin.setDecimals(4)
        self.tmax_spin.setDecimals(4)
        self.tmin_spin.setSingleStep(0.1)
        self.tmax_spin.setSingleStep(0.1) 

        #self.connect(self.tmin_spin, SIGNAL('valueChanged(double)'), self.parent.tminChanged)
        self.tmin_spin.valueChanged.connect(self.parent.tminChanged)
        #self.connect(self.tmax_spin, SIGNAL('valueChanged(double)'), self.parent.tmaxChanged)
        self.tmax_spin.valueChanged.connect(self.parent.tmaxChanged)
        self.tmin_spin.setValue(double(self.parent.tmin))
        self.tmax_spin.setValue(double(self.parent.tmax))

        self.RefreshEvent(True)

        #self.connect(self.tmin_spin, SIGNAL('valueChanged(double)'), self.tminChanged)
        self.tmin_spin.valueChanged.connect(self.tminChanged)
        #self.connect(self.tmax_spin, SIGNAL('valueChanged(double)'), self.tmaxChanged)
        self.tmax_spin.valueChanged.connect(self.tmaxChanged)

        #self.connect(self.OkButton, SIGNAL('clicked()'),self.CloseEvent)
        self.OkButton.clicked.connect(self.CloseEvent)
        
        #self.connect(self, SIGNAL('triggered()'), self.CloseEvent)
        self.closeEvent = self.CloseEvent
        
        #self.btnExit.clicked.connect(self.CloseEvent)
        #self.actionExit.triggered.connect(self.CloseEvent)
        #self.connect(self.refreshButton, SIGNAL('clicked()'),self.RefreshEvent)
        self.refreshButton.clicked.connect(self.RefreshEvent)
        #return 

    def tminChanged(self,value):
        if self.tmax_spin.value() < value:
            self.tmax_spin.setValue(value)
            self.parent.tmax_spin.setValue(value)
        self.parent.tmin_spin.setValue(value)


    def tmaxChanged(self,value):
        if self.tmin_spin.value() > value:
            self.tmin_spin.setValue(value)
            self.parent.tmin_spin.setValue(value)
        self.parent.tmax_spin.setValue(value)

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            self.close()



    def ReadWrongDets(self,init=False):
        dets = str(self.Edit_wrongDetectors.text())
        if len(dets.strip()) > 0:
            try:
                dets_list = eval('r_['+dets+']')
                wrong_dets_pref = array(int_(dets_list)-1, ndmin=1)
            except:
                wrong_dets_pref = []
                QMessageBox.warning(self,"Input problem",
                    "Wrong detectors are in bad format, use ...,10,11,13:15,...",QMessageBox.Ok)
                raise
        elif init:     # do not remove on init
            wrong_dets_pref = config.wrong_dets_pref
        else:
            wrong_dets_pref = []
        return  wrong_dets_pref
    
    def SetWrongDets(self, wrong):
        
        wrong = unique(wrong)

        wrong_dets_str = ''
        i = 0
        while i < len(wrong):
            wrong_dets_str = wrong_dets_str+str(wrong[i]+1)
            if i+1 < len(wrong) and wrong[i+1] == wrong[i]+1:
                while i+1 < len(wrong) and wrong[i+1] == wrong[i]+1:
                    i += 1
                wrong_dets_str = wrong_dets_str+':'+str(wrong[i]+2)
            i += 1
            wrong_dets_str = wrong_dets_str+','
        wrong_dets_str = wrong_dets_str[:-1]
            
        self.Edit_wrongDetectors.setText(wrong_dets_str)


            

    def SaveValues(self, init=False):

        """
        Save values from Data Settings dialog
        """
        self.setting['shot'] = int(self.parent.lineEdit_Shot.text())
        self.parent.tmin = self.tmin_spin.value()
        self.setting['tmin'] = self.parent.tmin
        self.parent.tmax = self.tmax_spin.value()
        self.setting['tmax'] = self.parent.tmax
        self.parent.data_smooth = self.smoothSpin.value()
        self.setting['data_smooth'] = self.parent.data_smooth
        self.parent.data_undersampling = self.UndersamplingSpin.value()
        self.setting['data_undersampling'] = self.parent.data_undersampling

            
        wrong_dets_pref =  self.ReadWrongDets(init)
        wrong = self.tokamak.dets[~self.tokamak.get_correct_dets(include_pref  = False)]
        config.wrong_dets_pref  = unique(setdiff1d( wrong_dets_pref,wrong))
        
        
                    
        if not os.path.exists(self.tokamak.geometry_path+'/wrong_channels/'):
            os.mkdir(self.tokamak.geometry_path+'/wrong_channels/')

        savetxt(self.tokamak.geometry_path+'/wrong_channels/%d'%self.setting['shot'], config.wrong_dets_pref,fmt='%d')

        self.parent.setting = self.setting
        self.parent.tokamak = self.tokamak

    def RefreshEvent(self, init=False):
        """
        Reload preview when Refresh button pressed
        """
        #self.message = QMessageBox()
        #self.message.setText("Loading ...")
        #self.message.setWindowTitle("Loading")
        #self.message.show()
        self.plot_2D = self.check_plot_2D.isChecked()
        self.SaveValues(init)
    
        try:
            ##===========make preview ============
            self.setting['tokamak_tmp'] = self.tokamak
            from prepare_data import loaddata, preview
            self.tokamak = loaddata(self.setting)
            preview(self.fig_data, self.setting, self.tokamak, self.plot_2D)
            self.canvas_data.draw()

        except:
            QMessageBox.warning(self,"Loading problem", "Data couldn't be loaded\nIf you changed tokamak now\ntry to remove tomography.npy file",QMessageBox.Ok)
            raise
        #finally:
            
            #self.message.close()


        #compress the wrong detectors to the easilly redable form
        wrong = self.tokamak.dets[~self.tokamak.get_correct_dets(include_pref=True)]
        
        self.SetWrongDets(wrong)
        

    def CloseEvent(self, *event):
        self.SaveValues()
        self.close()

        

        
        
        
    def create_spec_table(self):

        
        self.tab_widget = QWidget(self)
        self.main_tab.insertTab(1,self.tab_widget, 'Spectrogram')
        #self.tables_dict['Main'] = i_tab
        #i_tab+=1
        self.tab_widget.setSizePolicy(self.Expand)
        self.verticalLayout = QVBoxLayout(self.tab_widget)
        self.horizontalLayout = QHBoxLayout()


        self.verticalLayout.addLayout(self.horizontalLayout)


        #self.tab_widget.setLayout(vbox)
        #self.setCentralWidget(self.cWidget)
        
        self.resize(600, 650)
        self.setCenter()


        #self.dpi = 100
        self.fig_spec = Figure(dpi=self.dpi)
        self.canvas_spec = FigureCanvas(self.fig_spec)
        self.canvas_spec.setParent(self.cWidget)
        
        
        self.canvas_spec.setSizePolicy(self.Expand)
        self.canvas_spec.updateGeometry()
        
        #self.fig.patch.set_alpha(0.5)
        #self.fig.patch.set_facecolor('white')
        c =  self.cWidget.palette().color(QPalette.Base)
        
        self.fig_spec.patch.set_facecolor((c.red()/255., c.green()/255., c.blue()/255.))
        
        #self.Expand = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        #self.Fixed  = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)



        #self.cWidget = QWidget(self)
        #self.gridLayout = QGridLayout(self.cWidget)
        ##self.main_frame = QWidget()
        #self.main_frame = QTabWidget(self)
        #self.main_frame.setTabPosition(QTabWidget.North)
        ##QVBoxLayout
        #self.gridLayout.addWidget(self.main_frame, 0, 0, 1,2)

        #self.tab_widget = QWidget(self)
        #i_tab = 0
        #self.main_frame.insertTab(i_tab,self.tab_widget, 'Main')
  
        #self.tab_widget.setSizePolicy(self.Expand)
        #self.verticalLayout = QVBoxLayout(self.tab_widget)
        #self.horizontalLayout = QHBoxLayout()
        #self.verticalLayout.addLayout(self.horizontalLayout)
        #self.verticalLayout = QVBoxLayout(self.tab_widget)

        
        # Create the mpl Figure and FigCanvas objects. 
        #

        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        
        #from view_spectra import SpectraViewer
        self.snfft = QSlider(Qt.Horizontal)
        self.sgamma = QSlider(Qt.Horizontal)
        self.method  = QComboBox(self.cWidget)

        
        self.SpecWin =  SpectraViewer( self.sgamma, self.snfft ,self.method,
                    show_raw=True,allow_selector=True,fig= self.fig_spec)

        
        self.SpecWin.init_plot( self.tvec,self.signal,window = 'gauss')
        

        #self.axes = self.SpecWin.ax

        #self.axes = self.fig.add_subplot(111)
        
        #self.image = self.axes.imshow(np.random.randn(100,100))
        
        # Bind the 'pick' event for clicking on one of the bars
        #
        self.canvas_spec.mpl_connect('pick_event', self.on_pick)
        
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas_spec, self.cWidget)
        
        
        
        self.detector  = QComboBox(self.cWidget)
        self.detector.setToolTip('Select signal for the spectrogram')
        for N,det in zip(self.det_num, self.SXR_detectors):
            for i in arange(N):
                self.detector.addItem(det+'_%.3d'%(i+1)) 
                
                
        #[self.detector.addItem('det_%d'%det) for det in range(200)]
        self.detector.setCurrentIndex(self.ch0)
        #print self.ch0

        det_label2 = QLabel('Detector:')
        
        #self.connect(self.detector, SIGNAL('currentIndexChanged(int)'), 
            #self.change_det)
        self.detector.currentIndexChanged.connect(elf.change_det)
        #det_label2.setFixedWidth(3)

        #self.detector.addItem("STFT")
        #self.detector.addItem("SFFT")
        #self.detector.addItem("Sparse")
        #self.detector.setCurrentIndex(0)
        self.detector.setFixedWidth(80)

        
                #self.tab_widget_solve = QWidget(self)

        # Other GUI controls
        # 
        #self.textbox = QLineEdit()
        #self.textbox.setMinimumWidth(200)
        #self.connect(self.textbox, SIGNAL('editingFinished ()'), self.on_draw)
        self.method.setFixedWidth(70)

        #self.connect(self.method, SIGNAL('currentIndexChanged(int)'), 
                    #self.SpecWin.select_DFT_backend)
        self.method.currentIndexChanged.connect(self.SpecWin.select_DFT_backend)

        self.run_button = QPushButton("&Apply Filter")
        #self.connect(self.run_button, SIGNAL('clicked()'), self.run_svd_filt)
        self.run_button.clicked.connect(self.run_svd_filt)
        self.run_button.setFixedWidth(100)
        self.run_button.setToolTip('Evaluate SVD filter')
        self.method.setToolTip('Choose algorithm')
        self.snfft.setToolTip('Set width of window for FFT')
        self.sgamma.setToolTip('Set contrast')

        #self.grid_cb = QCheckBox("&Something..")
        
        #self.grid_cb.setChecked(False)
        #self.connect(self.grid_cb, SIGNAL('stateChanged(int)'), self.on_draw)
        
        slider_label2 = QLabel('NFFT:')
        slider_label1 = QLabel('GAMMA:')
        slider_label3 = QLabel('Method:')
        
        #spacer = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        #self.connect(sgamma, SIGNAL('valueChanged(int)'), self.on_draw)
        #self.connect(self.sgamma, SIGNAL('valueChanged(int)'),self.SpecWin.apply_slider)
        self.sgamma.valueChanged.connect(self.SpecWin.apply_slider)
        #self.connect(self.snfft , SIGNAL('valueChanged(int)'),self.SpecWin.apply_slider)
        self.snfft.valueChanged.connect(self.SpecWin.apply_slider)
        #
        # Layout with box sizers
        # 
        hbox1 = QHBoxLayout()
        hbox1.addWidget(self.mpl_toolbar)
        hbox1.setAlignment( self.mpl_toolbar, Qt.AlignVCenter)
        hbox2 = QHBoxLayout()

        for w in [ self.run_button,det_label2,self.detector,# self.grid_cb,
                  slider_label3, self.method]:
            hbox2.addWidget(w)
            hbox2.setAlignment(w, Qt.AlignVCenter)
        
        #hbox2.setAlignment(self.run_button, Qt.AlignVCenter)

        hbox2.setAlignment( slider_label3, Qt.AlignRight)
        hbox2.setAlignment( det_label2, Qt.AlignRight)

        for w in [ slider_label1, self.sgamma,slider_label2,self.snfft]:
            hbox1.addWidget(w)
            hbox1.setAlignment(w, Qt.AlignRight)
        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas_spec)

        vbox.addLayout(hbox1)
        vbox.addLayout(hbox2)

        #self.tab_widget.addLayout(vbox)
        #self.tab_widget.setLayout(vbox)
        self.horizontalLayout.addLayout(vbox)

        #self.horizontalLayout_rec.setLayout(vbox)
        #self.setCentralWidget(self.tab_widget)
        
        #self.resize(1050, 650)
        #self.setCenter()
        
        
        
    def create_SVD_table(self):
        
        
        
        
        
        self.tab_widget_svd = QWidget(self)
        self.main_tab.insertTab(2,self.tab_widget_svd, 'SVD filter')
        self.tab_widget_svd.setSizePolicy(self.Expand)
        self.verticalLayout_svd = QVBoxLayout(self.tab_widget_svd)
        self.horizontalLayout_svd = QHBoxLayout()
        #preview right
        #self.Picture_recons_1 = QLabel()
        #self.Picture_recons_1.setText('Preview profile X')
        #self.Picture_recons_1.setAlignment(QtCore.Qt.AlignCenter)
        #self.Picture_recons_1.setScaledContents(True)
        #self.Picture_recons_1.setSizePolicy(self.Fixed)
        #self.horizontalLayout_rec.addWidget(self.Picture_recons_1)



        self.verticalLayout_svd.addLayout(self.horizontalLayout_svd)


        #self.tab_widget.setLayout(vbox)
        #self.setCentralWidget(self.cWidget)
        
        #self.resize(600, 650)
        #self.setCenter()


        self.fig_svd = Figure((6.0, 5.0), dpi=self.dpi)
        self.canvas_svd = FigureCanvas(self.fig_svd)
        self.canvas_svd.setParent(self.cWidget)
        
        #self.fig_svd.patch.set_alpha(0.5)
        #self.fig_svd.patch.set_facecolor('white')
        c =  self.cWidget.palette().color(QPalette.Base)
        
        self.fig_svd.patch.set_facecolor((c.red()/255., c.green()/255., c.blue()/255.))
        #ax = self.fig_svd.add_subplot(111)
        #ax.plot([],[])
        
        self.mpl_toolbar_svd = NavigationToolbar(self.canvas_svd, self.cWidget)

        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas_svd)
        
        hbox = QHBoxLayout()
        
                
        
        slider_label1 = QLabel('N SVD:')
        slider_label2 = QLabel('N HARM:')
        
        self.nsvd = QSlider(Qt.Horizontal)
        self.nharm = QSlider(Qt.Horizontal)
        self.nsvd.setToolTip('Set number of used SVD components')
        self.nharm.setToolTip('Set number of filtered harmononics')
        
        self.run_button2 = QPushButton("&Run")
        #self.connect(self.run_button2, SIGNAL('clicked()'), self.run_svd_filt)
        self.run_button2.clicked.connect(self.run_svd_filt)
        self.run_button2.setFixedWidth(50)
        self.run_button2.setToolTip('Evaluate SVD filter')
        
        self.SVDF = SVDFilter( self.shot, self.fig_svd,self.ch0,self.SpecWin,self.statusBar ,self.nsvd,self.nharm)

        #hbox.addWidget(selfself.mpl_toolbar_svdelf.mpl_toolbar_svd, Qt.AlignVCenter)
        #self.SVDF.ch0 = self.ch0 
        
        self.run_button3 = QPushButton("&Save")
        #self.connect(self.run_button3, SIGNAL('clicked()'), self.SVDF.Save)
        self.run_button3.clicked.connect(self.SVDF.Save)
        self.run_button3.setFixedWidth(50)
        self.run_button3.setToolTip('Save filtered data')
        
        
        for w in [ self.mpl_toolbar_svd,slider_label1, self.nsvd,slider_label2,
                  self.nharm,self.run_button2,self.run_button3]:
            hbox.addWidget(w)
            hbox.setAlignment(w, Qt.AlignRight)
            

        
        #hbox.addWidget(self.mpl_toolbar_svd)
        #hbox.setAlignment( self.mpl_toolbar_svd, Qt.AlignVCenter)

        vbox.addLayout(hbox)

        #vbox.addLayout(hbox1)
        #vbox.addLayout(hbox2)

        #self.tab_widget.addLayout(vbox)
        #self.tab_widget.setLayout(vbox)
        self.horizontalLayout_svd.addLayout(vbox)
    
        #self.connect(self.nsvd, SIGNAL('valueChanged(int)'),self.SVDF.apply_slider)
        self.nsvd.valueChanged.connect(self.SVDF.apply_slider)
        #self.connect(self.nharm , SIGNAL('valueChanged(int)'),self.SVDF.apply_slider)
        self.nharm.valueChanged.connect(self.SVDF.apply_slider)
        

    def create_residuum_table(self):
        
        
        
        
        
        self.tab_widget_res = QWidget(self)
        self.main_tab.insertTab(3,self.tab_widget_res, 'Residuum')
        self.tab_widget_res.setSizePolicy(self.Expand)
        self.verticalLayout_res = QVBoxLayout(self.tab_widget_res)
        self.horizontalLayout_res = QHBoxLayout()
 

        self.verticalLayout_res.addLayout(self.horizontalLayout_res)

        self.fig_res = Figure((6.0, 5.0), dpi=self.dpi)
        self.canvas_res = FigureCanvas(self.fig_res)
        self.canvas_res.setParent(self.cWidget)
        
        c =  self.cWidget.palette().color(QPalette.Base)
        
        self.fig_res.patch.set_facecolor((c.red()/255., c.green()/255., c.blue()/255.))
 
        
        #self.mpl_toolbar_res = NavigationToolbar(self.canvas_res, self.cWidget)

        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas_res)
        
        #hbox = QHBoxLayout()
        
        self.SVDF.fig2 = self.fig_res
        
        self.mpl_toolbar_res = NavigationToolbar(self.canvas_res, self.cWidget)
        vbox.addWidget(self.mpl_toolbar_res)

        
        #slider_label1 = QLabel('N SVD:')
        #slider_label2 = QLabel('N HARM:')
        
        #self.nsvd = QSlider(Qt.Horizontal)
        #self.nharm = QSlider(Qt.Horizontal)
        #self.nsvd.setToolTip('Set number of used SVD components')
        #self.nharm.setToolTip('Set number of filtered harmononics')
        
        #self.run_button2 = QPushButton("&Run")
        #self.connect(self.run_button2, SIGNAL('clicked()'), self.run_svd_filt)
        #self.run_button2.setFixedWidth(50)
        #self.run_button2.setToolTip('Evaluate SVD filter')
        
        #for w in [ self.mpl_toolbar_svd,slider_label1,
                        #self.nsvd,slider_label2,self.nharm,self.run_button2]:
            #hbox.addWidget(w)
            #hbox.setAlignment(w, Qt.AlignRight)
            
        #self.SVDF = SVDFilter(self.fig_svd,self.SpecWin,self.statusBar, self.nsvd,self.nharm)

        #hbox.addWidget(selfself.mpl_toolbar_svdelf.mpl_toolbar_svd, Qt.AlignVCenter)
        #self.SVDF.ch0 = self.ch0 
        
        #hbox.addWidget(self.mpl_toolbar_svd)
        #hbox.setAlignment( self.mpl_toolbar_svd, Qt.AlignVCenter)
        
        #vbox.addLayout(hbox)

        #vbox.addLayout(hbox1)
        #vbox.addLayout(hbox2)

        #self.tab_widget.addLayout(vbox)
        #self.tab_widget.setLayout(vbox)
        self.horizontalLayout_res.addLayout(vbox)
    
        #self.connect(self.nsvd, SIGNAL('valueChanged(int)'),self.SVDF.apply_slider)
        #self.connect(self.nharm , SIGNAL('valueChanged(int)'),self.SVDF.apply_slider)
        
        

        

    
    def create_status_bar(self):
        self.status_text = QLabel("")
        self.statusBar().addWidget(self.status_text, 1)
        
    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        load_file_action = self.create_action("&Save plot",
            shortcut="Ctrl+S", slot=self.save_plot, 
            tip="Save the plot")
        
                
        load_data_action = self.create_action("Save &data",
            shortcut="Ctrl+D", slot=self.save_data, 
            tip="Save the data")
        
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (load_file_action,load_data_action, None, quit_action))
        
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About this tool')
        
        self.add_actions(self.help_menu, (about_action,))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)
                
                
    def setCenter(self):
        frect = QDesktopWidget.frameGeometry(self)
        frect.moveCenter(QDesktopWidget().availableGeometry(self).center());
        self.move(frect.topLeft())

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            action.triggered.connect(slot)
            #self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action
    
    def run_svd_filt(self):
        self.SVDF.load()
        self.SVDF.run_filter()
        #f0 = self.SVDF.f0
        self.SpecWin.set_f0(self.SVDF.f0 )
        self.canvas_svd.draw()
        self.canvas_res.draw()

        self.switch_tab(2)
        
    #def select_channel(self, ch0):
        #self.ch0 = ch0
        #self.SVDF.ch0  = ch0
        
        
    def switch_tab(self, tab):
        if tab == 0:
            pass
        if tab == 1:
            self.status_text.setText("Select time range by zooming, specify 1.harmonics and filter window width by selection")
            
        if tab == 2:
            self.status_text.setText("f0: %.4gkHz,\t RSS/TSS: %.3g%%"%(self.SVDF.f0/1000,
                                             self.SVDF.RSS/self.SVDF.TSS*100))
        if tab == 3:
            self.status_text.setText("Mouse: center button - select channel, right - full view")  
            
        
    def change_det(self,num):
        #print 'Not implemented yet!!'
        self.ch0=num
        #self.SVDF.ch0  = num-sum(self.SVDF.wrong_det < num) 
        #self.ch0 -= sum(self.wrong_det < self.ch0) 
        self.load_data()
        #t_range = self.SpecWin.t_range
        #f_range = self.SpecWin.f_range
        #print t_range, f_range, 't_range, f_range'
        xlims = self.SpecWin.ax.get_xlim()
        ylims = self.SpecWin.ax.get_ylim()

        
        self.SpecWin.init_plot( self.tvec,self.signal,window = 'gauss', 
                        tmin=xlims[0],tmax=xlims[1],fmin0=ylims[0],fmax0=ylims[1])
        
    
        #BUG  misí se načíst data
        



#def main():
    
    ##shot = 29222
    ##fig = figure()
    
    ##svdf = SVDFilter(shot,fig,ch0, SpecWin,statusbar,slider_svd=None,slider_harm=None,ind=slice(None,None))
    
    
    #app = QApplication(sys.argv)
    #form = SettingWindow()
    #form.show()
    #app.exec_()


#if __name__ == "__main__":
    #main()
