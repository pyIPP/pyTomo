#!/usr/bin/env python
# -*- coding: utf-8 -*-
from matplotlib.pylab import *
from numpy import *
from shutil import copyfile
from scipy.interpolate import LSQUnivariateSpline, UnivariateSpline

from scipy.interpolate import RectBivariateSpline, interp1d
from scipy.stats.mstats import mquantiles
from matplotlib.colors import LogNorm,NoNorm
from matplotlib import colors
import os, sys
from shared_modules import MovingAveradge
from time import time
import config
import fconf
from matplotlib.ticker import MaxNLocator,LogFormatter,ScalarFormatter,NullFormatter
from matplotlib.transforms import Bbox
from shared_modules import debug
import matplotlib.cm as cm
from matplotlib.figure import Figure
from matplotlib import rcParams,colorbar
from matplotlib.backends.backend_agg import FigureCanvas as FigureCanvas



try:
    from multiprocessing import Process, Pool, cpu_count
    import threading
    threading._DummyThread._Thread__stop = lambda x:40
except:
    pass


from make_graphs  import my_cmap_
my_cmap1 = my_cmap_
my_cmap2 = cm.get_cmap('RdBu')
my_cmap2._init()
my_cmap2._lut[:-3,:] = my_cmap2._lut[-4::-1,:]
my_cmap2._lut[-3:-1,:] = my_cmap2._lut[-3:-1,:][::-1]
my_cmap2._lut/= my_cmap2._lut.max(0)


class LogFormatterTeXExponent(LogFormatter, object):
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

        label = label.replace('--','-')
        x = float(label)
        label = '%.1e'%x
        
        if abs(log10(abs(x)+1e-10)) >2.001 and abs(x)>1e-10:  
            label_ = re.sub(r'e(\S)0?(\d+)',r'\\cdot 10^{\1\2}',str(label))

            label_ = "$" + label_ + "$"
            label_ = label_.replace('+','')
        elif abs(x)>1e-10:
            n = max(0,-int(log10(abs(x))-1))
            label_ = ('$%.'+str(n)+'f$')%x
        else:
            label_ = '$0$'
 
        return label_
        


class MyFormatter(ScalarFormatter): 
    def __call__(self, x, pos=None): 
        if pos==0: 
            return '' 
        else: return ScalarFormatter.__call__(self, x, pos)  



def Asymmetry2dplot(tokamak,tvec,rho,G0,GS,GC,plot_up_down_asym,rho_icrh=None,rem_back=False):


    if len(tvec) < 2: return 
    nrows = 3 if plot_up_down_asym else 2
    
    
    fig = Figure((13,6.5))
    FigureCanvas(fig)
    axes = []
    ax = None
    for i in range(nrows):
        ax = fig.add_subplot(nrows,1,i+1, sharex=ax)
        if i < nrows-1:
            for label in ax.get_xticklabels():
                label.set_visible(False)
            ax.xaxis.offsetText.set_visible(False)

        axes.append(ax)
    
    fig.subplots_adjust(hspace=0.00, wspace = 0.05)
    nonlinG0 = arcsinh(G0/std(G0)*10)
    nonlinG0 /= amax(abs(nonlinG0))
 
    GS,GC = copy(GS), copy(GC)
    G0[~isfinite(G0)| (G0==0)] = nan
    GS/=G0
    GC/=G0
 

    extent=(tvec.min(),tvec.max(),rho.max(),rho.min())
    limG = amax(abs(mquantiles(G0[isfinite(G0)], (0.001,0.999))))

    if not rem_back:
        limC = amax(abs(mquantiles(GC[isfinite(GC)&(G0 > .1*limG)], (.01,.995))))
        limS = amax(abs(mquantiles(GS[isfinite(GS)&(G0 > .1*limG)], (.01,.995))))
    else:
        limC, limS = 0.3, 0.3
    
    im1 = axes[0].imshow(GC,vmax=limC,vmin=-limC, aspect='auto', cmap=my_cmap2, extent=extent)
    if plot_up_down_asym: im2 = axes[1].imshow(GS,vmax=limS,vmin=-limS, aspect='auto', cmap=my_cmap2, extent=extent)
    im3 = axes[-1].imshow(G0, aspect='auto', cmap=my_cmap2 if rem_back else my_cmap1
                         , extent=extent, norm=LogNorm(vmin=.002*limG,vmax=limG))

    try:
        CS = axes[-1].contour(tvec,rho, G0, 20, colors='k',linewidths=.05,alpha=.2 )
    except:
        pass

     
    for ax in axes:
        ax.set_ylabel(tokamak.rho_label)

        if not all([r is None for r in array(rho_icrh,ndmin=1)]) and any(isfinite(rho_icrh)):
            ax.plot(tvec, rho_icrh,'k-', lw=1.5 ) #LFS
            ax.plot(tvec,-rho_icrh,'b--',lw=1.5 ) #HFS
    
        ax.set_ylim(1,0)
        ax.yaxis.set_major_formatter( MyFormatter() )

                
        if  hasattr(tokamak, 'impur_inject_t') and not tokamak.impur_inject_t is None:
            for lbo in tokamak.impur_inject_t:
                if lbo>=amin(tvec) and lbo<=amax(tvec):
                    ax.axvline(x=lbo,color='w',lw=1,ls='-')
                    ax.axvline(x=lbo,color='k',lw=1,ls='--')
            
    
    axes[-1].set_xlabel('time [s]')
    axes[0].yaxis.set_major_formatter( ScalarFormatter() )
 
    tick_locator = MaxNLocator(nbins=5)

   
    cax,kw = colorbar.make_axes(axes[0])
    cb = fig.colorbar(im1, cax=cax, **kw)
    cb.set_label(r'in/out asymmetry')
    cb.locator = tick_locator
    cb.update_ticks()

    if plot_up_down_asym:
        cax,kw = colorbar.make_axes(axes[1])
        cb = fig.colorbar(im2, cax=cax, **kw)
        cb.set_label(r'up/down asymmetry')
        cb.locator = tick_locator
        cb.update_ticks()
    cax,kw = colorbar.make_axes(axes[-1])
    cb = fig.colorbar(im3, cax=cax, **kw)
    cb.set_label(r'zero order')


    savez_compressed(tokamak.tmp_folder+'/asym_data_%d_%.3f'%(tokamak.shot,mean(tvec) ),
                     tvec=tvec,GC = single(GC),GS = single(GS), G0=single(G0),rho=rho,rho_icrh=rho_icrh)
    
    fig.savefig(tokamak.tmp_folder+'asymmetry_2D_%d.png'%tokamak.shot ,bbox_inches='tight',transparent=True, dpi=80 )



    
def CalcAsymmetry(tmp):
    os.nice(20)



    (tvec,jtvec, shot, magx,magy,R,Z,G,rc,zc,rho_icrh,rho_pol,ymin,
                   ymax,lin_treshold,n_nods,log_lin,auto_scale, make_plot,
                   enable_output,tsteps,plot_up_down,BdMat,rho_label,tmp_folder,geom_path,output_path) = tmp
    rho_mag = outer(ones(magx.shape[0]), linspace(0,1,magx.shape[1]))
    from scipy.interpolate import  NearestNDInterpolator

    
    rho = zeros(shape(R)+shape(tvec),dtype=single)
    BdMat = BdMat.reshape(shape(R), order='F')
    #map main grid on the rho_pol
    for jt,t in enumerate(tvec):
        NNDI = NearestNDInterpolator(vstack((magx[...,jt].ravel(),magy[...,jt].ravel())).T, rho_mag.ravel())            
        rho[...,jt] = NNDI(vstack((R.ravel(),Z.ravel())).T).reshape(R.shape)
        rho[...,jt]+= (random.randn(*rho[...,jt].shape)-.5)/magx.shape[1]/2.
        
    HFS = (R[...,None]-rc < -abs(Z[...,None]-zc))&~BdMat[:,:,None]   #isfinite(rho)
    LFS = (R[...,None]-rc >  abs(Z[...,None]-zc))&~BdMat[:,:,None]   #isfinite(rho)

    Up =  (Z[...,None]-zc >  abs(R[...,None]-rc))&~BdMat[:,:,None]   #isfinite(rho)
    Down= (Z[...,None]-zc < -abs(R[...,None]-rc))&~BdMat[:,:,None]   #isfinite(rho)

    Left  = R[...,None]<rc
    Right = R[...,None]>rc
    
    rho[abs(rho)>1] = 1


    
    up_down_rho = copy(rho)
    left_right_rho = copy(rho)
    left_right_rho[HFS]*=-1
    up_down_rho[Down]*=-1
    

    G0,Gc,Gs = [], [],[]
    
    #NOTE correction of the method in order to make it consistent
    #with sin/cos decomposition
    corr = 1.1
  
    
    if make_plot:

        fig = Figure()
        FigureCanvas(fig)

        nplots = 3 if plot_up_down else 2
        fig = Figure((5*nplots+1,6))
        FigureCanvas(fig)

        fig.subplots_adjust(hspace=0, wspace = 0)
        ax1 = fig.add_subplot(1,nplots,1)
        if plot_up_down: 
            ax2 = fig.add_subplot(1,nplots,2)
            ax3 = fig.add_subplot(1,nplots,3)
        else:
            ax3 = fig.add_subplot(1,nplots,2)

        ax3.yaxis.tick_right()
        ax3.yaxis.set_label_position("right")

        if log_lin:
            if plot_up_down: ax2.set_yscale('symlog', linthreshy=lin_treshold)
            ax1.set_yscale('symlog', linthreshy=lin_treshold)
            ax1.axhline(y=0,color='k')
            if plot_up_down: ax2.axhline(y=0,color='k')

        points1_lfs, = ax1.plot([],[],'b,')
        points1_hfs, = ax1.plot([],[],'r,')

        line_lfs, = ax1.plot(rho_pol,rho_pol,'b',label='LFS')
        line_hfs, = ax1.plot(-rho_pol,rho_pol,'b',label='HFS')
        line_hfs2, = ax1.plot(rho_pol,rho_pol,'r',label='mirror HFS')
        ax1.set_ylim(ymin,ymax)
        
        if plot_up_down:  ax2.set_ylim(ymin,ymax)
        ax1.set_xlim(-1,1)
        
        
        if plot_up_down:
            ax2.set_xlim(-1,1)
            points2_left, = ax2.plot([],[], 'r,')
            line_up_left, = ax2.plot(rho_pol,rho_pol,'r',label='HFS UP')
            line_down_left, = ax2.plot(-rho_pol,rho_pol,'r',label='HFS DOWN')
            line_down2_left, = ax2.plot(rho_pol,rho_pol,'m',lw=.5,label='mirror HFS DOWN')

            points2_right, = ax2.plot([],[], 'b,')
            line_Up_right, = ax2.plot(rho_pol,rho_pol,'b--',label='LFS UP')
            line_down_right,  = ax2.plot(-rho_pol,rho_pol,'b--',label='LFS DOWN')
            line_down2_right, = ax2.plot(rho_pol,rho_pol,'y--',lw=.5,label='LFS DOWN')
            ax2.yaxis.set_major_locator(MaxNLocator(5))

        r_plot  = ax1.axvline(x=-1,c='0.5')
        r_plot_ = ax1.axvline(x=-1,ls='--',c='0.5')
        tlt = ax1.set_title('Midplane profile cut  shot:%d'% shot)

        ax1.text(0.25,0.8, 'HFS' ,backgroundcolor='w',transform=ax1.transAxes)
        ax1.text(0.75,0.8, 'LFS' ,backgroundcolor='w',transform=ax1.transAxes)
        ax1.set_xlabel(r'%s (negative for HFS)'%rho_label)

        ax1.set_ylabel(r' $\epsilon$ [kW/m$^3$] ')
        ax1.axvline(x=0,ls="-.")
        ax1.yaxis.set_major_locator(MaxNLocator(5))
        ax3.yaxis.set_major_locator(MaxNLocator(5))

        
        try:
            Emiss0 = load(tmp_folder+'Emiss0.npz')#BUG
            G_ = Emiss0['G']
            R,Z = Emiss0['Rc'], Emiss0['Zc']
            ir = argmin(abs(R-mean(rc)))
            iz = argmin(abs(Z-mean(zc)))
            r = linspace(-1,1,100)
            srho = rho.mean(2)[iz,:]*sign(R-mean(rc))
            ax1.plot(sort(srho)   ,  G_.mean(2)[iz,argsort(srho)]/1e3,'k' )

            crho = rho.mean(2)[:,ir]*sign(Z-mean(zc))
            ax2.plot(sort(crho),G_.mean(2)[argsort(crho),ir]/1e3,'k'  )


        except Exception as e:
            pass

        
        if plot_up_down:
            tlt = ax2.set_title('')
            ax2.text(0.25,0.8, 'Down',backgroundcolor='w',transform=ax2.transAxes)
            ax2.text(0.75,0.8, 'Up'  ,backgroundcolor='w',transform=ax2.transAxes)
            ax2.yaxis.set_major_formatter( NullFormatter() )
            ax2.xaxis.set_major_formatter( MyFormatter() )
            ax2.set_xlabel(r'%s (negative for Down)'%rho_label)
            
            ax2.axvline(x=0,ls="-.")
            asym_sin, = ax3.plot(rho_pol,rho_pol,'r--',label='up/down')


        
        
        ax3.axhline(y=0,ls="--",c='k',lw=.5)
        asym_cos, = ax3.plot(rho_pol,rho_pol,'b',label='in/out')

        model = False
        plt_br = None

        try:
            s,c = loadtxt(tmp_folder+'sin'), loadtxt(tmp_folder+'cos')
            ax3.plot(linspace(0,1,len(c)),c,'b:',label='in/out  model')
            if plot_up_down: ax3.plot(linspace(0,1,len(s)),s,'r:',label='up/down model')
            text3 = ax3.text(0.75,0.8,'',backgroundcolor='w',transform=ax3.transAxes)
            text4 = ax3.text(0.75,0.7,'',backgroundcolor='w',transform=ax3.transAxes)
            model = True
        except:
            pass
        
        try:
            assert not model ,' Used phantom data'
            low_asym  = load('./asymmetry/asymmetry_emiss_%d_low.npz'%shot)#BUG #uncertainty in Mach +/- 5%
            high_asym = load('./asymmetry/asymmetry_emiss_%d_up.npz'%shot)#BUG
            asym = load('./asymmetry/asymmetry_emiss_%d_.npz'%shot)#BUG

            if 'rho_pol' in low_asym:
                rho_pol_ = low_asym['rho_pol']
                rho_pol_ = interp1d(high_asym['tvec'],rho_pol_,axis=0,bounds_error=False,assume_sorted=True)(tvec)
            else:
                print('BUG rhotor')
                rho_pol_ = tile(low_asym['rhotor'], (len(tvec), 1))

            Ac_high = interp1d(high_asym['tvec'],high_asym['Ac'],axis=0,bounds_error=False,assume_sorted=True)(tvec)*.9  #correction of the step method!
            Ac_low  = interp1d( low_asym['tvec'], low_asym['Ac'],axis=0,bounds_error=False,assume_sorted=True)(tvec)*1
            Ac  = interp1d( asym['tvec'], asym['Ac'],axis=0,bounds_error=False,assume_sorted=True)(tvec)


            
            model_up,   = ax3.plot(rho_pol_[0,:],rho_pol_[0,:],'g:')
            model_down, = ax3.plot(rho_pol_[0,:],rho_pol_[0,:],'g:')
            model, = ax3.plot(rho_pol_[0,:],rho_pol_[0,:],'g-.',label='in/out model')

            try:
                data = load('./asymmetry/asymmetry_emiss_%d_.npz'%shot)#BUG
                
                
                Bremsstrahlung = interp1d(data['tvec'],data['Bremsstrahlung'],axis=0,bounds_error=False,assume_sorted=True)(tvec)
                W_Radiation = interp1d(data['tvec'],data['W_Radiation'],axis=0,bounds_error=False,assume_sorted=True)(tvec)
                plt_br,   = ax1.plot(rho_pol_[0,:],rho_pol_[0,:],'k')
                plt_w, = ax1.plot(rho_pol_[0,:],rho_pol_[0,:],'k--')
                plt_br2, = ax2.plot(rho_pol_[0,:],rho_pol_[0,:],'k')
                plt_w2,  = ax2.plot(rho_pol_[0,:],rho_pol_[0,:],'k--')
            except:
                pass

            
        except Exception as e:
            pass
                
        try:
            assert  plt_br is None
            BRdata = load(geom_path+'/Bremsstrahlung_%d.npz'%shot)
            plt_br,   = ax1.plot([],[],'k')
            plt_br2,  = ax2.plot([],[],'k')
            plt_w,   = ax1.plot([],[],'y')
            plt_w2,  = ax2.plot([],[],'y')
            
            if 'asym' in BRdata and not 'model_up' in vars():
                Ac = interp1d(BRdata['tvec'],BRdata['asym'],axis=0,bounds_error=False,assume_sorted=True)(tvec)
                Ac_low =   Ac/1.1
                Ac_high = Ac*1.1

                rho_pol_ = tile(BRdata['rho'], (len(tvec), 1))
                        
                model_up,   = ax3.plot(rho_pol_[0,:],rho_pol_[0,:],'g:')
                model_down, = ax3.plot(rho_pol_[0,:],rho_pol_[0,:],'g:')
                model, = ax3.plot(rho_pol_[0,:],rho_pol_[0,:],'g-.',label='in/out model')


        except:
            pass
 
        ax3.xaxis.set_major_formatter( MyFormatter() )

        r_plot2 = ax3.axvline(x=-1,c='0.5')
        r_plot2_ = ax3.axvline(x=-1,ls='--',c='0.5')
        ax3.set_ylabel(r' $(\varepsilon-\langle \varepsilon \rangle )/ \langle \varepsilon \rangle  $ ')
        ax3.set_xlabel(rho_label)
        ax3.set_xlim(0,1)
        ax3.set_ylim(-.6, .9)
        ax3.grid(True)

        
        
        if enable_output:
            ax3.legend(loc='upper left')
        
        

       
    T = time()

    for jt,t in enumerate(tvec):
        try:
            sys.stdout.write("\r %2.1f %%" %(jtvec[jt]*100./tsteps))
            sys.stdout.flush()
        except:
            pass
 
        r = left_right_rho[:,:,jt][HFS[:,:,jt]|LFS[:,:,jt]]
        sorted_rho,ind  = unique(r,return_index=True)
        sorted_rho = r_[-1.005,sorted_rho,1.005]
        
        g = G[:,:,jt][HFS[:,:,jt]|LFS[:,:,jt]]
        g = r_[0,g[ind],0]
        g_norm = median(g)/2
        g_norm = abs(g_norm)+ abs(mean(g))/10.
        sorted_emiss = arcsinh(g/g_norm)


        try:
            s = LSQUnivariateSpline(sorted_rho, sorted_emiss, linspace(-.95,.95,n_nods),bbox=[-1.01,1.01],k=4)
        except:
            print('LSQUnivariateSpline')
            import IPython
            IPython.embed()
            
            
        G_HFS = sinh(s(-rho_pol))*g_norm
        G_LFS = sinh(s( rho_pol))*g_norm
        
        if make_plot:
            ind_r = sorted_rho > 0
            points1_lfs.set_data(sorted_rho[ind_r],g[ind_r]/1e3)
            points1_hfs.set_data(sorted_rho[~ind_r],g[~ind_r]/1e3)

            line_lfs.set_ydata( G_LFS/1e3) 
            line_hfs.set_ydata( G_HFS/1e3)
            line_hfs2.set_ydata( G_HFS/1e3)
            if not rho_icrh[jt,:] is None and not all(isnan((rho_icrh[jt,:]))):
                r_plot.set_xdata([nanmean(rho_icrh[jt,:]),]*2)
                r_plot_.set_xdata([-nanmean(rho_icrh[jt,:]),]*2)
                
                r_plot2.set_xdata([nanmean(rho_icrh[jt,:]),]*2)
                r_plot2_.set_xdata([-nanmean(rho_icrh[jt,:]),]*2)
                
            
      
        r = up_down_rho[:,:,jt][(Up[:,:,jt]|Down[:,:,jt]) & Left[:,:,jt]]
        sorted_rho2,ind  = unique(r,return_index=True)

        sorted_rho2 = r_[-1.005,sorted_rho2,1.005]
        
        g = G[:,:,jt][(Up[:,:,jt]|Down[:,:,jt]) & Left[:,:,jt]]
        g = r_[0,g[ind],0]
        sorted_emiss = arcsinh(g/g_norm)
        
        s = LSQUnivariateSpline(sorted_rho2, sorted_emiss, linspace(-0.95,0.95, n_nods),bbox=[-1.01,1.01],k=4)
   
        
        G_Up   = sinh(s( rho_pol))*g_norm
        G_Down = sinh(s(-rho_pol))*g_norm
        

        
        if make_plot and plot_up_down:
            points2_left.set_data(sorted_rho2,g/1e3)
            line_up_left.set_ydata(G_Up/1e3)
            line_down_left.set_ydata(G_Down/1e3)
            line_down2_left.set_ydata(G_Down/1e3)
        
        
        r = up_down_rho[:,:,jt][(Up[:,:,jt]|Down[:,:,jt]) & Right[:,:,jt]]
        sorted_rho3,ind  = unique(r,return_index=True)
        sorted_rho3 = r_[-1.005,sorted_rho3,1.005]
        
        g = G[:,:,jt][(Up[:,:,jt]|Down[:,:,jt]) & Right[:,:,jt]]
        g = r_[0,g[ind],0]
        sorted_emiss = arcsinh(g/g_norm)

        s = LSQUnivariateSpline(sorted_rho3, sorted_emiss, linspace(-0.95,0.95,n_nods),bbox=[-1.01,1.01],k=4)

        
        G_Up2   = sinh(s( rho_pol))*g_norm
        G_Down2 = sinh(s(-rho_pol))*g_norm
        
        
        if make_plot and plot_up_down:
            points2_right.set_data(sorted_rho3,g/1e3)
            line_Up_right.set_ydata(G_Up2/1e3)
            line_down_right.set_ydata(G_Down2/1e3)
            line_down2_right.set_ydata(G_Down2/1e3)
            asym_sin.set_ydata((G_Up+G_Up2-G_Down-G_Down2)/(G_Up+G_Down+G_Down2+G_Up2+1e-6)*corr)

        
        if make_plot:
            
            asym_cos.set_ydata((G_LFS-G_HFS)/(G_LFS+G_HFS+1e-6)*corr)

            try:
                model_up.set_data(rho_pol_[jt,:],Ac_high[jt,:])
                model_down.set_data(rho_pol_[jt,:],Ac_low[jt,:])
            except:
                pass
            try:
                model.set_data(rho_pol_[jt,:],Ac[jt,:])
            except:
                pass
            try:
     
                plt_br.set_data( r_[-rho_pol_[jt,::-1],rho_pol_[jt]],r_[Bremsstrahlung[jt,::-1],Bremsstrahlung[jt,:]]/1e3)  
                plt_w.set_data(  r_[-rho_pol_[jt,::-1],rho_pol_[jt]],r_[   W_Radiation[jt,::-1],   W_Radiation[jt,:]]/1e3)  
                plt_br2.set_data(r_[-rho_pol_[jt,::-1],rho_pol_[jt]],r_[Bremsstrahlung[jt,::-1],Bremsstrahlung[jt,:]]/1e3)  
                plt_w2.set_data( r_[-rho_pol_[jt,::-1],rho_pol_[jt]],r_[   W_Radiation[jt,::-1],   W_Radiation[jt,:]]/1e3)  
            
            except:
                pass
            
                       
            if auto_scale:
                ax1.set_ylim(G[:,:,jt].min()/1e3,G[:,:,jt].max()/1e3*1.1)
                if plot_up_down: ax2.set_ylim(ax1.get_ylim())

        
   
            try:
                it = argmin(abs(BRdata['tvec']-t))
                it1 = BRdata['tvec'].searchsorted((tvec[max(0,jt-1)]+t)/2)
                it2 = BRdata['tvec'].searchsorted((tvec[min(len(tvec)-1,jt+1)]+t)/2)+1

                x = r_[-BRdata['rho'][::-1], BRdata['rho']]
                y = r_[nanmean(BRdata['BR'][it1:it2,::-1],0),  nanmean(BRdata['BR'][it1:it2],0)]
                y_up = r_[nanmean(BRdata['BR_low'][it1:it2,::-1],0),  nanmean(BRdata['BR_low'][it1:it2],0)]
                y_down = r_[nanmean(BRdata['BR_high'][it1:it2,::-1],0),  nanmean(BRdata['BR_high'][it1:it2],0)]
                
                x = r_[-BRdata['rho'][::-1], BRdata['rho']]
                y = r_[BRdata['BR'][it,::-1],  BRdata['BR'][it]]
                y_up = r_[BRdata['BR_low'][it,::-1],  BRdata['BR_low'][it]]
                y_down = r_[BRdata['BR_high'][it,::-1],  BRdata['BR_high'][it]]
                ymin,ymax = ax1.get_ylim()
                
                
                plt_br.set_data(x,y/1e3)
                plt_br2.set_data(x,y/1e3)  
                
                y = r_[BRdata['W'][it,::-1],  BRdata['W'][it]]
                plt_w.set_data(x,y/1e3)
                plt_w2.set_data(x,y/1e3)  
            except:
                pass
                
                
 
            tlt.set_text("Time : %.4fs"%t)
            fig.savefig(tmp_folder+'asymmetry_plot_%.4d.png'%jtvec[jt],
                        bbox_inches=Bbox([[1.3,0.05],[15.3,5.7]]), dpi=80)
     
            if enable_output:
                fig.savefig(output_path+'asymmetry_plot_%.4d.pdf'%jtvec[jt],
                                bbox_inches='tight', dpi=80 )

            try:
                text3.set_text('$\Delta_c$ = %.2f'%std(asym_cos.get_ydata\
                    -interp(rho_pol,linspace(0,1,len(c)),c)))
                text4.set_text('$\Delta_s$ = %.2f'%std(asym_sin.get_ydata\
                    -interp(rho_pol,linspace(0,1,len(s)),s)))
            except:
                pass
            
      
            
        Gs.append((G_Up+G_Up2-G_Down-G_Down2)/4)
        Gc.append((G_LFS-G_HFS)/2*corr)
        G0.append((G_HFS+G_LFS)/2)
    
    return  G0, Gc, Gs 
    
    
    
    
#aby to bralo plot_symlog, auto_scale, 
#kreslit paraelně 

def EvalAsymmetry(input_data):
    print('EvalAsymmetry')

    input_parameters, tokamak, progress,results = input_data


    G = results['g'].T
    tvec=results['tvec']
    

    #estimate an averadge asymmetry by Dider'ss approach 
    #r,z coordinates of the middles in the pixel grid
    #dR,dZ - shift EQH equilibrium by dR,dZ
    #emiss - array(nr,nz,nt) or array(nr,nz) 
    tsteps = len(tvec)
    
    rem_back = input_parameters['rem_back']

    if rem_back:
            #substract averadge of the first 10% of points! 
        subst_ind = slice(0,int(ceil(len(tvec)*input_parameters['bcg_subst_fract'])))
        G = copy(G)-G[:,subst_ind].mean(1)[:,None]
    
        #position during this time interval was already corrected
        rhop,magx, magy = tokamak.mag_equilibrium(tvec, n_rho=100, return_mean = True)
        magx = rollaxis(tile(magx, (tsteps,1,1)).T,1)
        magy = rollaxis(tile(magy, (tsteps,1,1)).T,1)

    else:
        rhop,magx, magy = tokamak.mag_equilibrium(tvec, n_rho=100,n_theta=50, return_mean = False)


    R,Z = meshgrid(tokamak.xgrid+tokamak.dx/2,tokamak.ygrid+tokamak.dy/2) #position of pixel centers


    rc = magx[:,0,:].mean(0)
    zc = magy[:,0,:].mean(0)
        
 

    mag_tvec = tvec
    try:
        G = G.reshape((tokamak.ny, tokamak.nx,-1), order='F')
    except:
        print((tokamak.ny, tokamak.nx, G.shape))
        raise 
    rho_icrh, r_icrh = Icrh_rho(tokamak,tvec,magx,magy)

    rho_pol = linspace(0,.95,200)
    if rho_icrh is None:
        rho_icrh = array((None,)*len(tvec),ndmin=2,dtype=double).T
    

    try:
        n_cpu = cpu_count()
        n_split = 1 if tsteps < 2*n_cpu else 2*n_cpu
        p = Pool(min(n_cpu,n_split ))
    except:
        n_split = 1

    ind = [slice(i*tsteps//n_split,(i+1)*tsteps//n_split) for i in range(n_split)]

    n_nods = min(int(sqrt(tokamak.nx*tokamak.ny)/7)+10, tokamak.nx*tokamak.ny/8)

        
    if rcParams['backend'].lower() != 'agg':  #multithreading is  working only with AGG
        ind = [slice(0,tsteps),]
    ymax = amax(G)/1e3
    ymin = amin(G)/1e3 if rem_back else 0

    plot_loglin = input_parameters['plot_loglin']
    plot_all = input_parameters['plot_all'] or tsteps == 1
    auto_scale = input_parameters['plot_autoscale']
    enable_output = input_parameters['enable_output']
    plot_up_down_asym = input_parameters['plot_up_down_asym']
    radial_coordinate = input_parameters['radial_coordinate']
    output_path = input_parameters['output_path']


    lin_treshold = mean(abs(G))/10
    from annulus import get_bd_mat
    BdMat = get_bd_mat(tokamak,time=tvec.mean())
    

    args = [(tvec[ii],r_[ii],tokamak.shot,magx[:,:,ii],magy[:,:,ii],R,Z,
             G[:,:,ii],rc[ii],zc[ii],rho_icrh[ii],rho_pol,ymin,ymax,lin_treshold,n_nods,
             plot_loglin,auto_scale,plot_all,enable_output,tsteps,plot_up_down_asym,
             BdMat,tokamak.rho_label,tokamak.tmp_folder,tokamak.geometry_path,output_path) for ii in ind]

    try:
        #pp
        if plot_all and rcParams['backend'].lower() != 'agg' or config.DEBUG:
            raise Exception()
        out = list(p.imap(CalcAsymmetry,args,1 ))
    except Exception as e:
        print(e)
        out = list(map(CalcAsymmetry,args ))
    finally:
        try:
            p.close()
            p.join() 
        except:pass


    G0 = vstack([A[0] for A in out])
    Gc = vstack([A[1] for A in out])
    Gs = vstack([A[2] for A in out])
    plot_up_down_asym = input_parameters['plot_up_down_asym']


        
    nr = len(rho_pol)
    
 
    if len(tvec) > 1:
        try:
            Asymmetry2dplot(tokamak,tvec,rho_pol,G0.T,Gs.T,Gc.T,plot_up_down_asym,rho_icrh,rem_back)
        except:
            pass

    if True:
        Gc_ = (Gc/(abs(G0)+1e-3))[:,int(nr*input_parameters['rho_asym'])]
        Gs_ = (Gs/(abs(G0)+1e-3))[:,int(nr*input_parameters['rho_asym'])]

        asym_ICRH = vstack(( Gs_,Gc_ )).T
        rho_icrh_mean = input_parameters['rho_asym']
    else:
        rho_icrh2 = nanmean(rho_icrh,1)
        rho_icrh_mean = nanmean(rho_icrh2)
        if isnan(rho_icrh_mean): rho_icrh_mean = input_parameters['rho_asym']


        rho_icrh2[isnan(rho_icrh2)] = rho_icrh_mean
        rho_icrh2[(rho_icrh2>.99)|(rho_icrh2<-.99)] = 0
        Gc_ = (Gc/(abs(G0)+1e-3))[arange(len(tvec)),int_(rho_icrh2*nr)]
        Gs_ = (Gs/(abs(G0)+1e-3))[arange(len(tvec)),int_(rho_icrh2*nr)]

        asym_ICRH = vstack(( Gs_,Gc_ )).T



    #Asymmetry1dplot(tvec,asym_ICRH[:,None,:],rho_icrh_mean,rho_icrh,   tokamak , enable_output,plot_up_down_asym )
    
    
    

    rho_p = linspace(0,1,results['asym_error_sin'].shape[1])
    

    profiles = {'tvec':tvec,'rho_p':rho_p,'FSA':single(results['asym_0']) ,
                'Fc': single(results['asym_cos']),'Fs':single(results['asym_sin']) }
    
    profiles['FSA_err'] = single(results['asym_error_0'])
    profiles['Fc_err'] = single(results['asym_error_cos'])
    profiles['Fs_err'] = single(results['asym_error_sin'])

    savez_compressed(tokamak.tmp_folder+'asymmetry_results%d_%.1f.npz'%(tokamak.shot, tvec[0]),**profiles)

    return

 


def CalcAsymNew(output,Tok, G,data, error,retro, tvec, dets,SigmaGsample,input_parameters):
    
    

    tsteps = len(tvec)
    
    BdMat = output['bnd'][0]
    tcalc = time()

    from scipy.ndimage.interpolation import map_coordinates

    if 'magx' in output:
        magx = output['magx'].T
        magy = output['magy'].T
        rhop = output['mag_rho']

    else:
        rhop,magx, magy = Tok.mag_equilibrium(tvec, return_mean=False)
        output['magx'] = magx.T
        output['magy'] = magy.T
        output['mag_rho'] = rhop
    
    #WARNING it is not exact flax surface averadge!!
    n_mag = size(magx,1)
    n_theta = size(magx,0)

    profile_cos = zeros((tsteps,n_mag),dtype=single)
    profile_sin = zeros((tsteps,n_mag),dtype=single)
    profile_0   = zeros((tsteps,n_mag),dtype=single)
    rc = magx[:,0,:].mean(0)
    zc = magy[:,0,:].mean(0)
        
    scaling = array([Tok.dx,Tok.dy])
    offset  = array([Tok.xmin,Tok.ymin])+scaling/2
    theta = arctan2(magy-zc,magx-rc)
    theta = single(unwrap(theta, axis=0))
    
    wrong = abs(theta[0,1]- theta[-1,1]-2*pi) >1e-3
    theta[:,1,wrong] =  theta[:,2,wrong]  #correction in case of poor magnetic data
    G = G.reshape(Tok.ny,Tok.nx,tsteps,order='F')
    theta0 = linspace(-pi,pi,n_theta, endpoint=False)

    from scipy.linalg import lstsq,qr_multiply,solve_triangular
    
    
    A = ones((n_theta,3))
    for i in range(tsteps): 
        coords = c_[magx[:,:,i].ravel(),magy[:,:,i].ravel()].T
        idx = (coords-offset[:,None])/ scaling[:,None] 
        map_prof = map_coordinates(G[...,i].T,idx,order=1)
        map_prof = map_prof.reshape((n_theta,n_mag))

        for j in range(n_mag):
            A[:,1] = cos(theta[:,j,i])
            A[:,2] = sin(theta[:,j,i])
            x,r,r,s=lstsq(A,map_prof[:,j],overwrite_a=True,overwrite_b=False, check_finite=False)
            profile_cos[i,j] = x[1]/x[0]
            profile_sin[i,j] = x[2]/x[0]
            profile_0[i,j] = x[0]
    profile_cos[:,0] = 0
    profile_sin[:,0] = 0
    profile_0[:,0] = profile_0[:,1]

    
    output['asym_0'] = profile_0
    output['asym_cos'] = profile_cos
    output['asym_sin'] = profile_sin
    


    asym_error_0   = zeros((tsteps, n_mag),dtype=single)
    asym_error_cos = zeros((tsteps, n_mag),dtype=single)
    asym_error_sin = zeros((tsteps, n_mag),dtype=single)
    
    if SigmaGsample is None:
        SigmaGsample = zeros((Tok.ny*Tok.nx, 1))
    SigmaGsample = copy(SigmaGsample)

    n_samples =  SigmaGsample.shape[-1]
    SigmaGsample = SigmaGsample.reshape((Tok.ny, Tok.nx,n_samples), order='F')
    map_prof = zeros((n_samples,n_theta,n_mag))
    for it,t in enumerate(tvec): 

        err_cos = zeros((n_mag, n_samples),dtype=single)
        err_sin = zeros((n_mag, n_samples),dtype=single)
        err_0   = zeros((n_mag, n_samples),dtype=single)
        eq_err_r,eq_err_z = input_parameters['expected_mag_eq_err']
        for i in range(n_samples): 
            dr = eq_err_r*random.randn()
            dz = eq_err_z*random.randn()
            coords = c_[magx[:,:,it].ravel()+dr,magy[:,:,it].ravel()+dz].T
            idx = (coords-offset[:,None])/ scaling[:,None] 
            map_prof[i] = map_coordinates((SigmaGsample[...,i]+G[:,:,it]).T,idx,order=1).reshape((n_theta,n_mag))

        for j in range(n_mag-1):
            A[:,1] = cos(theta[:,j,it])
            A[:,2] = sin(theta[:,j,it])
            x,r,r,s=lstsq(A,map_prof[:,:,j].T,overwrite_a=True,overwrite_b=False, check_finite=False)
            #BUG!!! use the integral, proper scalar product!!

            err_cos[j] = x[1]/x[0]
            err_sin[j] = x[2]/x[0]
            err_0[j] = x[0]


        edge_err = 1e-3*profile_0[it].max()/(abs(profile_0[it])+1e-5) #edge accuracy error - it should not be more accurate than 0.1%of max

        asym_error_0[it]   = std(err_0,1)*input_parameters['error_scale']
        asym_error_cos[it] = std(err_cos,1)*input_parameters['error_scale']+edge_err
        asym_error_sin[it] = std(err_sin,1)*input_parameters['error_scale']+edge_err
        asym_error_cos[it][0] = 1e-3
        asym_error_sin[it][0] = 1e-3

  
    output['asym_error_0']   = asym_error_0
    output['asym_error_cos'] = asym_error_cos
    output['asym_error_sin'] = asym_error_sin
    

    debug('asymmetry calc time %.3f'%(time()-tcalc ))
    
    return output
  
    


def Icrh_rho(tokamak,tvec,magx=None,magy=None ):
    
    rho_icrh,r_icrh = None,None

    if hasattr(tokamak, 'ICRH_rezonance'):
        r_icrh = interp1d(tokamak.ICRH_rezonance['tvec'],tokamak.ICRH_rezonance['R'],
                          axis=0,bounds_error=False,assume_sorted=True)(tvec)
        if magx is None or magy is None:
            rhop,magx, magy = tokamak.mag_equilibrium(tvec)
            
        n_mag = size(magy,1)
        
        r_in = magx.min(0)
        r_out = magx.max(0)
        rho_mag = linspace(-1,1,n_mag*2)
        rho_icrh = empty_like(r_icrh)  #can be negative!!
        for i in range(len(tvec)): 
            r = r_[r_in[::-1,i],r_out[:,i]]
            rho_icrh[i,:] = interp(r_icrh[i,:],r[r!=0], rho_mag[r!=0] )
        
        rho_icrh[isnan(r_icrh)] = nan
    return rho_icrh, r_icrh

            
            
