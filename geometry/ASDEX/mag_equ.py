#!/usr/bin/env python 
# -*- coding: utf-8 -*-



import numpy as np
from numpy import *
from scipy.interpolate import interp1d
from time import time
from scipy.signal import  medfilt



import os,sys
from . import dd
dd = dd.shotfile()



from .descur import DESCUR
D = DESCUR()


from multiprocessing import  Pool, cpu_count




class MOM2RZ:
    #object for the fast evaluation of the equilibrium for arbitrary rho and theta
    def __init__(self,rho,rcos,rsin,zcos,zsin,order=3,R0=0,Z0=0,regularization=1e5):
        nmom = np.size(rcos,-1)

        self.p_rcos = zeros((order+1+nmom,nmom))
        self.p_zcos = zeros((order+1+nmom,nmom))
        self.p_rsin = zeros((order+1+nmom,nmom))
        self.p_zsin = zeros((order+1+nmom,nmom))

        order_range = r_[:order:2]

        #do not use zero order term for zero order mode
        #self.p_rcos[-order-1:,0] = polyfit_reg(rho,rcos[:,0],order_range+1)
        #self.p_zcos[-order-1:,0] = polyfit_reg(rho,zcos[:,0],order_range+1)
        #self.p_rsin[-order-1:,0] = polyfit_reg(rho,rsin[:,0],order_range+1)
        #self.p_zsin[-order-1:,0] = polyfit_reg(rho,zsin[:,0],order_range+1)

        #use n to n+order  polynom order for higher modes ( it will have n zero derivations rho zero)  and zero first derivative
        

        for i,n in enumerate(r_[1,arange(1,nmom)]):
            
            #regularization is important, else rounding errors can be huge
            self.p_rcos[-n-order_range-1,i],r = self.polyfit_reg(rho,rcos[:,i],order_range+n,regularization)
            self.p_zcos[-n-order_range-1,i],r = self.polyfit_reg(rho,zcos[:,i],order_range+n,regularization)
            self.p_rsin[-n-order_range-1,i],r = self.polyfit_reg(rho,rsin[:,i],order_range+n,regularization)
            self.p_zsin[-n-order_range-1,i],r = self.polyfit_reg(rho,zsin[:,i],order_range+n,regularization)



        #shift it by centers
        self.p_rcos[-1,0]+= R0
        self.p_zcos[-1,0]+= Z0


     
    def __call__(self,rho,theta):
        nmom = np.size(self.p_rcos,-1)

        rcos,rsin,zcos,zsin = self.polyval(rho)

        angle = np.outer(np.arange(nmom),theta)
        cos = np.cos(angle)
        sin = np.sin(angle)
        r_plot = np.tensordot(rcos,cos,axes=([-1,0]))
        r_plot+= np.tensordot(rsin,sin,axes=([-1,0]))
        z_plot = np.tensordot(zcos,cos,axes=([-1,0]))
        z_plot+= np.tensordot(zsin,sin,axes=([-1,0]))

        return r_plot, z_plot

    def polyval(self,rho):
        nmom = np.size(self.p_rcos,-1)

        rcos = np.empty((np.size(rho),nmom))
        rsin = np.empty((np.size(rho),nmom))
        zcos = np.empty((np.size(rho),nmom))
        zsin = np.empty((np.size(rho),nmom))

        for i in range(nmom):
            rcos[:,i] = np.polyval(self.p_rcos[:,i],rho)
            zcos[:,i] = np.polyval(self.p_zcos[:,i],rho)
            rsin[:,i] = np.polyval(self.p_rsin[:,i],rho)
            zsin[:,i] = np.polyval(self.p_zsin[:,i],rho)
            
        return rcos,rsin,zcos,zsin


        
    def polyfit_reg(self, x,y,order,regularization=1,weight=1):
        ind = ~any(isnan(y.reshape(len(x),-1)),1)
        x = x[ind]
        y = array(y.T,ndmin=2).T
        
        y = y[ind,...]
        
        y/= x[:,None]**(order.min()-1)
        order-= order.min()-1
        
        N_order = amax(order)+1

        if isscalar(order):
            order  = arange(order)
        order = N_order - order-1
            
        lhs = vander(x, N_order)[:,order]
    
        scale = sqrt((lhs*lhs).sum(axis=0))
        lhs /= scale
        rhs = y
        rcond = len(x)*finfo(x.dtype).eps*regularization

        c, resids, rank, s = linalg.lstsq(lhs*x[:,None] , rhs*x[:,None], rcond)
            
        c = (c.T/scale).T # broadcast scale coefficients

        #C = zeros((N_order,)+y.shape[1:])
        #C[order,...] =  c

        return squeeze(c),resids




def Descur_fit_core(params):

    (t_fract, i, rho, R_contour, Z_contour, R0, Z0,n_fourier) = params
    sys.stdout.write("\r calculate Fourier coefficients: %3.0f %%    N: %d        " %(t_fract*100,i))
    sys.stdout.flush()  
    n_rho = len(rho)
    
    


    moments_all= empty((n_rho, n_fourier,4))
    for nr, (cr,cz) in enumerate(zip(R_contour.T, Z_contour.T)):
        #don't use to high orrder in the core - regularization
        n_fourier_ = 2+int(ceil((nr+1)/float(n_rho)*(n_fourier-2)))
        moments_all[nr,:n_fourier_,:] =  D.descur_fit(cr-R0,cz-Z0,n_fourier_)
        moments_all[nr,n_fourier_:,:] = 0

        

    rcos = moments_all[:,:,0]
    rsin = moments_all[:,:,1]
    zcos = moments_all[:,:,2]
    zsin = moments_all[:,:,3]
    poly_order = 30
    
    try:
        mom_poly = MOM2RZ(rho,rcos,rsin,zcos,zsin,order=poly_order,R0=R0,Z0=Z0)
        coeff = dstack((mom_poly.p_rcos, mom_poly.p_zcos, mom_poly.p_rsin, mom_poly.p_zsin ))
    except Exception as e:
        print('Descur_fit_core',e)
        coeff = ones((poly_order+1+n_fourier,n_fourier, 4))*nan
        

    
    
    return single(coeff)

    


#different methods for correctig of the signals perturbated by elms 
def ElmCorrection(tvec_fast, sig_fast, tvec_slow, sig_slow, t_elms, dt_elms,nsmooth=11,mode=None):
    #sig_fast fast is the signal with the correct high frequancy part and wrong low frequency
    #sig_slow is opposite 
 
    
    #make a high time resolution signal
    sig = sig_fast-interp(tvec_fast,tvec_slow,interp(tvec_slow,tvec_fast,sig_fast))+interp(tvec_fast,tvec_slow,sig_slow)
    tvec = tvec_fast
    
    if mode== None:
        return tvec,sig

    #prepare intervals of the elms which will be removed from the analysis
    n_elm  = size(t_elms)
    tvec_elm = vstack((t_elms-3e-4, t_elms,t_elms+ dt_elms, t_elms+ dt_elms+3e-4)).T.ravel()
    elm = vstack((zeros(n_elm),ones(n_elm),ones(n_elm), zeros(n_elm))).T.ravel()
  
    elms_interp = interp1d(tvec_elm, elm,bounds_error = False,fill_value= 0 )(tvec)
    elms_ind = elms_interp != 0
    #just replace ELM region by interpolation
    if mode == 'interpolate':
        sig[elms_ind] = interp1d(tvec[~elms_ind],sig[~elms_ind],bounds_error=False,fill_value=0)(tvec[elms_ind])
        return tvec, sig
    
    #replace low frequency part of fast data by low frequency part of the slow data
    if mode== 'replace': #replace elms by FPG data
        dsig = sig-sig_fast
        dsig = interp1d(tvec[~elms_ind], dsig[~elms_ind],bounds_error = False,fill_value= 0 )(tvec_slow)
        
        dsig = medfilt(dsig,nsmooth)
        dsig = interp1d(tvec_slow, dsig,bounds_error = False,fill_value= 0 )(tvec)

        sig = dsig+sig_fast
        
        return  tvec, sig
    
    if mode== 'remove': #remove elms
        
        for t,dt in zip(t_elms,dt_elms) :
            sig_fast[~((tvec_fast<t-dt*0.1)|(tvec_fast>t+dt*1.1))] = nan
            sig_slow[~((tvec_slow<t-dt*0.1)|(tvec_slow>t+dt*1.1))] = nan
            
        
                    
        
        sig_fast[isnan(sig_fast)] = interp(tvec_fast[isnan(sig_fast)],tvec_fast[~isnan(sig_fast)], sig_fast[~isnan(sig_fast)] )
        sig_slow[isnan(sig_slow)] = interp(tvec_slow[isnan(sig_slow)],tvec_slow[~isnan(sig_slow)], sig_slow[~isnan(sig_slow)] )

        #sig = interp1d(tvec[~elms_ind], sig[~elms_ind],bounds_error = False,fill_value= 0 )(tvec_slow)
    
        
        dt_elms = median(diff(t_elms))
        if isnan(dt_elms): dt_elms = 0
        n_smooth = int(dt_elms/mean(diff(tvec_slow)))*2+1
        #print 'n_smooth',n_smooth
        
        dsig = medfilt(sig_slow,n_smooth)-medfilt(interp(tvec_slow,tvec_fast, sig_fast ),n_smooth)
        
        
        sig_fast+= interp(tvec_fast, tvec_slow, dsig)
        
        sig = medfilt(sig_fast,5)

        

        return tvec,sig

    
    #averadge over the elms 
    if mode == 'averadge':
        #TODO vylepšit!!!!! moc to průměruje 26061
        
        sig = interp1d(tvec[~elms_ind], sig[~elms_ind],bounds_error = False,fill_value= 0 )(tvec_slow)
        dt_elms = median(diff(t_elms))
        if isnan(dt_elms): dt_elms = 0
        n_smooth = int(dt_elms/mean(diff(tvec_slow)))*2+1
        #print 'n_smooth',n_smooth
        
        sig = medfilt(sig,n_smooth)

        
        sig = interp1d(tvec_slow, sig,bounds_error = False,fill_value= 0 )(tvec)

        return tvec, sig 






def help_fun(tmp):
    (eqm, rho, theta, time) = tmp
    return eqm.rhoTheta2rz(rho, theta,time, n_line=100)


class Equlibrium:
    def __init__(self,eqm,shot,  diag='EQI',exp='AUGD', ed=0):
        self.eqm = eqm
        self.shot = shot
        self.diag = diag
        

        
        
        
        
        
    def getTranspEquilibrium(self):
        
        
        from netCDF4 import Dataset

        cdf = None
        paths = 'm/CDF','TRANSP','CDF' 
        
        flag = False
        for i in range(100,0,-1):
            if flag: break
            for path in paths:
                try:
                    cdf_file = os.path.expanduser('~/'+path+'/%dA%.2d.CDF'%(self.shot,i))
                    cdf = Dataset(cdf_file, 'r', format='NETCDF4')
                except:
                    continue
                flag = True
                break
            
            

        
        if cdf is None:
            raise Exception( 'loading of cdf_file file has failured')


        #load TRANSP shotfile
        ntim,n_rho=cdf.variables['RMC00'].shape

        tvec=cdf.variables['TIME3'][:]
        rho_p = cdf.variables['PLFLX'][:]
        rho_p = sqrt(rho_p/rho_p[:,(-1,)])
        rho_t = cdf.variables['XB'][:]

    
        R0 = cdf.variables['RAXIS'][:]/100
        Z0 = cdf.variables['YAXIS'][:]/100

        nmom=0
        rc='RMC%.2d'%nmom
        while rc in iter(cdf.variables.keys()):
            nmom += 1
            rc='RMC%.2d'%nmom

        rcos = zeros((ntim,n_rho, nmom))
        rsin = zeros((ntim,n_rho, nmom))
        zsin = zeros((ntim,n_rho, nmom))
        zcos = zeros((ntim,n_rho, nmom))
        

        for jmom in range(nmom):
            rc='RMC%.2d'%jmom
            zc='YMC%.2d'%jmom
            rcos[:,:,jmom] = cdf.variables[rc][:]/100
            zcos[:,:,jmom] = cdf.variables[zc][:]/100


        for jmom in range(1,nmom):
            rs='RMS%.2d'%jmom
            zs='YMS%.2d'%jmom
            rsin[:,:,jmom] = cdf.variables[rs][:]/100
            zsin[:,:,jmom] = cdf.variables[zs][:]/100


        zcos[:,:,0]-= Z0[:,None]
        rcos[:,:,0]-= R0[:,None]

        coeffs = []
        for it in range(ntim):
            #BUG sometimes, very rarelly it can failure if rho_p is deviating from straight line in the core!!! 
            #during core ECRH current drive!
            mom_poly = MOM2RZ(rho_p[it],rcos[it],rsin[it],zcos[it],zsin[it],order=30,R0=R0[it],Z0=Z0[it],regularization=0.1)
            coeff = dstack((mom_poly.p_rcos, mom_poly.p_zcos, mom_poly.p_rsin, mom_poly.p_zsin ))
            coeffs.append(coeff)
 
            
            
        coeffs = array(list(coeffs))

        tvec = tvec_surf
        Rmag_ = zeros_like(tvec_surf)
        Zmag_ = zeros_like(tvec_surf)
        ahor_ = ones_like(tvec_surf)
        bver_ = ones_like(tvec_surf)
        
        
        output = {'tvec_surf':tvec_surf,
        'surf_coeff':coeff ,  
        'tvec_fast':tvec , 
        'Rmag':Rmag_,
        'Zmag':Zmag_,
        'ahor':ahor_,
        'bver':bver_  
        }
        
        
        return output

        
        
        
    def getStandartEquilibrium(self):
    

        self.eqm.read_ssq()
        self.eqm._read_scalars()
        self.eqm._read_profiles()
        self.eqm._read_pfm()


        os.nice(3)
        n_rho = 40
        n_theta = 50
        rho = linspace(0.01,0.998,n_rho)
        theta = linspace(0,2*pi,n_theta, endpoint=False)
        
        
        
        corrupted_eq = abs(self.eqm.PFxx[3]) > 1000
        
        ncpu = cpu_count()
        pool = Pool(ncpu)
        
        t_eq = self.eqm.t_eq[~corrupted_eq]
        nti_kk = len(t_eq)

        t_sequence = array_split(t_eq, min(ncpu,nti_kk))

        print('calculate contours ')
        t1 = time()
        args = [(self.eqm,rho, theta, t) for t in t_sequence]
        
        
        import config
        #from shared_modules import debug

        if config.DEBUG:
            out = list(map(help_fun,args ))
        else:
            out = pool.map(help_fun,args )
        R_cont,z_cont = hstack(out)
        print('Calc. contours: %.1f s'%( time()-t1))


        R0 = self.eqm.ssq['Rmag'][~corrupted_eq]
        Z0 = self.eqm.ssq['Zmag'][~corrupted_eq]

                
       
        
        mom_order = 6
        #mom_order = 10 #for DIII-D

        args = [(jt/float(nti_kk),jt,rho,R_cont[jt],z_cont[jt],R0[jt], Z0[jt], mom_order) for jt in arange(nti_kk)]
        
        

        t1 = time()
        
        print('\n fitting contours  fit:')
        if config.DEBUG:
            coeffs = list(map(Descur_fit_core,args  ))

        else:
            #coeffs = [ Descur_fit_core(a) for a in args]
            coeffs = pool.map(Descur_fit_core,args  )

            
        print('\nContouts fit: %.1f s'%(time()-t1))

        pool.close()
        pool.join()
        
        coeffs = array(list(coeffs))
        
        ind = all(isfinite(coeffs.reshape(nti_kk, -1)),1)
        surf_coeff = coeffs[ind]
        tvec_surf = t_eq[ind]


        tvec_slow  = t_eq
        Rc_slow = self.eqm.ssq['Rmag'][~corrupted_eq]
        Zc_slow = self.eqm.ssq['Zmag'][~corrupted_eq]
        bver_slow = self.eqm.ssq['bver'][~corrupted_eq]
        ahor_slow = self.eqm.ssq['ahor'][~corrupted_eq]
        
      

        
        if not dd.Open('FPG',self.shot):
            raise Exception('no FPG!!')

        tvec_fast  = dd.GetTimebase('TIMEF')
        Rc_fast = dd.GetSignal('Rmag')
        Zc_fast = dd.GetSignal('Zmag')
        tvec_fast = tvec_fast[:len(Zc_fast)]
        Zc_fast = dd.GetSignal('Zmag')
        bver_fast = dd.GetSignal('bver')
        ahor_fast = dd.GetSignal('ahor')
        dd.Close()

        if dd.Open('ELM',self.shot):
            tELM = dd.GetTimebase('dt_ELM')
            dt_ELM= dd.GetSignal('dt_ELM')
            dd.Close()
        else:
            tELM,dt_ELM = array((0,)),array((1,))




        mode_r = 'remove' if  size(tELM) >1 else 'replace'

        tvec,Zmag_ = ElmCorrection(tvec_fast, Zc_fast,   tvec_slow, Zc_slow,   
                                   tELM, dt_ELM,nsmooth=11,mode='replace')
        tvec,Rmag_ = ElmCorrection(tvec_fast, Rc_fast,   tvec_slow, Rc_slow,   
                                   tELM, dt_ELM,nsmooth=11,mode=mode_r)
        tvec,bver_ = ElmCorrection(tvec_fast, bver_fast, tvec_slow, bver_slow, 
                                   tELM, dt_ELM,nsmooth=11,mode='replace')
        tvec,ahor_ = ElmCorrection(tvec_fast, ahor_fast, tvec_slow, ahor_slow, 
                                   tELM, dt_ELM,nsmooth=11,mode='replace')

        bver_slow = interp(tvec_surf,tvec_slow, bver_slow)
        Zc_slow   = interp(tvec_surf,tvec_slow,   Zc_slow)
        Rc_slow   = interp(tvec_surf,tvec_slow,   Rc_slow)
        ahor_slow = interp(tvec_surf,tvec_slow, ahor_slow)



        #renormalize shape of the surfaces
        surf_coeff[:,-1,0,0] -= Rc_slow
        surf_coeff[:,-1,0,1] -= Zc_slow
        surf_coeff[...,0::2] /= ahor_slow[:,None,None,None]
        surf_coeff[...,1::2] /= bver_slow[:,None,None,None]

        if self.diag != 'EQI':
            #remove ELM perturbated timepoints 
            elm_ind = zeros_like(tvec_surf, dtype=bool)
            for t,dt in zip(tELM,dt_ELM):
                it1 = searchsorted(tvec_surf,t )
                it2 = searchsorted(tvec_surf,t+dt )
                elm_ind[it1-1: it2]  = True
                        
        
            tvec_surf = tvec_surf[~elm_ind]
            surf_coeff = surf_coeff[~elm_ind,...]

        output = {'tsurf':single(tvec_surf),
                'surf_coeff':single(surf_coeff),  
                'tvec_fast':single(tvec) , 
                'Rmag':single(Rmag_),
                'Zmag':single(Zmag_),
                'ahor':single(ahor_),
                'bver':single(bver_ ) 
              }

        return output
    
 


  
        
#EQU = Equlibrium(33974,diag='EQH')
#EQU.getStandartEquilibrium()


