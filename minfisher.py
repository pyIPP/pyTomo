#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
from scipy import sparse
#from  matplotlib.pyplot import *
import gc
from numpy import linalg
import time
from scipy.sparse import spdiags, eye
from scipy.sparse.linalg import factorized
import os
from scipy.interpolate import *
#import config 
from shared_modules import prepare_inputs, create_derivation_matrix, debug


#AUTOR: Tomas Odstrcil  tomas.odstrcil@ipp.mpg.de


try:
    from sksparse.cholmod import  cholesky,CholmodWarning,CholmodError, analyze, cholesky_AAt
except:
    print("Can't find Cholesky decomposition for sparse matrices !!! \nFalling back to slow version (~6x slower)\n Install scikits.sparse package")

    pass
    """
    import cvxopt ;import cvxopt.cholmod
        # Solve system 
        K = cvxopt.spmatrix(K.data,K.row.astype(np.int),K.col.astype(np.int))
        B = cvxopt.matrix(f[free,0])
        cvxopt.cholmod.linsolve(K,B)
        u[free,0]=np.array(B)[:,0] 
    """



def  minfisher(data, error,tvec,Tmat,dets, normData, G0, lam0, Tok, reg_params, boundary, regularization, ifishmax,
                postprocessing = False, progress=0):
  

    """
    Main algorithm for pixel tomography based on Minimum Fisher Regularisation with rapid version and anisotropic smoothing.
    This algorithm identify boundaries and create virtual sensors on boundary, prepare matrices of derivation, rescale and renormalize all inputs (Tmat, data, error)
    This file also contains outer cyclus of MFI (nonlinear iteration)


    :param array data:      Measured data
    :param array error:     Expected errors of data
    :param array tvec:      Time vector
    :param spmatrix Tmat:      Geometry matrix
    :param int normTmat:    Norm of geometry matrix
    :param array normData:  Maximum of data in each timeslice
    :param array G0:        Initial guess of emissivity
    :param double lam0:     Initial value of searched lambda
    :param class tok: ``Class`` of  tokamak with all important data.
    :param double danis:    Ratio of anisotropic matrices, B = n*Bperp+1/n*Bpar  => The smaller number the more will reconstruction follow the field For MFI is recommended 0.2, for MDIFF < 0.1
    :param double boundary:    If non-zero allow smoothing with boundary => The reconstruction will be zero on the boundary,  Set the pressure forcing reconstruction to be zero on boundary. Precision of the virtual senzors (boundary) is min(sigma)/X.
    :param int regularization: Apriory smoothing used for minimalization
    :param int ifishmax: Number of Minimum Fisher loops (outer loops), recommended is 3
    :param bool postprocessing:   Perform postprocessing --  entropy, smoothness
    :param object progress: link to QObject that move with the progressbar, do not work with multiprocessing
    :var array g: Searched final reconstruction
    :var array lam:  Array of all smoothing parameters lambda used during reoconstruction
    :var array chi2:  Array of all chi2 used during reconstruction
    :var int nl:  Number of detectors
    :var int iimax: Maximal iteration in inner cycle
    :var int BdMat:  Array nonzero outside of tokamak

    :var spmatrix diam,diar,dial,diao,diau: Sparse matrices of derivation (elementar directions)

    :var spmatrix Bx, By, Bperp, Bpar: Sparse matrices of derivation


    :var spmatrix Err:  Square root of variation matrix
    :var spmatrix fs:   Normalised data
    :var spmatrix Ts:   Normalised geometry matrix
    :var spmatrix TT:   Ts.T*Ts

    :var array boundary_vec:  Vector nonzero on boundary
    :var array lambdab:    Initial values of lambda for new step in inner cycle
    :var array G_all:    Array of all results
    :var array Flags:   Array of all flags from reconstruction

        flag = 0  -- OK

        flag = 1  -- Solution doesn't exist

        flag = 2  -- Didn't find solution

    :var int min_change:  Minimal change between steps in outer step. If the change is lower then min_change, reconstruction ends.


    :var spmatrix H:    Final version of derivation matrix
    :var spmatrix W:    Weighting matrix
    :var array entropy: Entropy
    :var array smoothnest_izo:  Smoothness izotropic
    :var array smoothnest_anizo: Smoothness anizotropic



    .. sectionauthor:: Michal Odstrčil <michalodstrcil@gmail.com>
    .. sectionauthor:: Tomáš Odstrčil <todstrci@ipp.mpg.de>
    .. sectionauthor:: Jan Mlynar


    """

    if any(isnan(error)):
        raise Exception('NaNs in the  errorbars!')
    if any(isnan(data)) :
        raise Exception('NaNs in the data !')
    if all(~isfinite(error)) :
        raise Exception('all errors are infinity')

    reconstruct = True
    tsteps = len(tvec)

    time_all = time.time()

    debug('\n=========Minimum Fisher Information Regularisation ============')
 

    Bmat,BdMat,bnd,Ts, fs,err,tvec, G0 = prepare_inputs(Tok,data, error, dets,
                    tvec, Tmat, normData, G0, reg_params, boundary, regularization,
                    reconstruct, postprocessing)


    
    gc.collect()        # free unused RAM memory !!!
    npix = Tok.npix
    
    TT = Ts.T*Ts
    
    Ts = Ts.tocsr()
    danis = reg_params['danis']


    debug('Matrix of geometry - density %2.1f %%  tsteps: %i'  %(TT.nnz*100./npix**2, tsteps))
    debug('Outer cycle')

    
    ndets = Tmat.shape[0]
    n_virt = fs.shape[0]-ndets


    g, chi2last, lam_last,chi2all = outer_cycle(Tok, Bmat,danis, TT,Ts, fs, G0 , 
            BdMat,tsteps, regularization, ifishmax,npix,lam0,progress, time_all,boundary)
    gc.collect()        # free unused RAM memory !!!

    debug('Outer cycle end')
    debug('----------------------------------------')



    if Tok.transform_index in [0,4]:
        g[BdMat] = 0
        
        
    g*=normData[None,:]/Tok.normTmat
    debug('------------ALL - time  %g' % (time.time()-time_all))


    Tmat = (Tok.Transform.T*Tmat.T).T  #Tok.Transform can have problems with __rmul__ 
    retro = (Tmat*g).T

    return tvec, g, retro, chi2last, lam_last,[single(bnd)]*tsteps,None



def outer_cycle(Tok, Bmat, danis, TT,Ts, fs, g, BdMat,tsteps, 
                regularization, ifishmax,npix,lam0,  progress, time_all,boundary ):
    """
    Outer loop of MFI recontruction, calls Netwon method to find chi2 == 1 root

    :param list Bmat: list of derivation matrices for smoothing
    """

    debug('===============Fast MFR================')

    iimax=15       # number of chi2 iterations loops (inner loops)
    lam=ones((ifishmax+5,iimax+4))*NaN
    chi2=ones((ifishmax+5,iimax+4))*NaN
    
    lam[:,0] = 30  #lambda larger than lam_max
    chi2[:,0] =  log(mean(sum(fs**2, axis=0)/sum(fs!=0,axis=0)))
    lam[:,1] = -10
    mfs = mean(fs, 1)
    pinvG,istop,itn,r1norm,r2norm,anorm,acond,arnorm,xnorm,var=sparse.linalg.lsqr(Ts, mfs,atol=1e-9, btol=1e-9)
    chi2[:,1] = log(sum((mfs-Ts*(abs(pinvG)+pinvG)/2)**2,axis=0)/sum(mfs!=0))
    
    lambdab = zeros((ifishmax+5+1, 1))

    Flags = ones((ifishmax+5, 1))*10
    lambdab[0]  = 3 if lam0 is None else  lam0  #first lambda
    
    
    dG = 100
    flag = 1
    i=0
    if tsteps == 1:
        min_change =  2 #% # minimal change in recontruction to break cyclus
        precision_num = 40  # variable accuracy of chi2 search, the higher the more precise
    else:
        precision_num = 10
        min_change =  1 #%  # minimal change in recontruction to break cyclus

       
    debug('------------init - time  %g' % (time.time()-time_all))

    debug('Outer cycle prepared')

    factor = False
    g1 = mean(g,axis=1)[:,None]
    while  ((dG > min_change and not regularization in (2,3)) or ((dG > min_change/2  or flag != 0) 
                            and regularization in (2,3,))) and i < ifishmax  :
        debug('%d    loop'%(i+1))
        H = create_derivation_matrix(g1, Bmat,BdMat, regularization, danis, Tok)
        precision  = 0.5**i/precision_num         #variable accuracy for chi2tar
        #debug(('------------init 2 - time  %g' % (time.time()-time_all)))
        g0 = g1
        del g

                
        chi2,lam,g,chi2all,flag,lambdab[i+1]=root_finder(lam,chi2,i,iimax,lambdab[i], 
                            TT,Ts,  H, fs, g0, BdMat, precision, factor,  Tok, lam0 is not None )
 
        #flag = 0  - OK
        #flag = 1  - Solution doesn't0 exist
        #flag = 2  - Didn't find solution

        debug('------------root_finder - time  %g' % (time.time()-time_all))



        g = reshape(g, (npix , tsteps), order='F')
        Flags[i] = flag
        
        ##=========Statistics =============
        g1 = mean(g,axis=1)[:,None]
        dG = mean(fabs(g0 - g1))/mean(fabs(g0) + fabs(g1))*100
        
        debug('difference in reconstruction: %5.2f%% ' % dG)
        if progress != 0:
            progress.iterateSubStep()
        i += 1
  
  
    #============= End of the fast MinFisher loop

    last = i-1    
    
    lam_last = ( lam[last,where(~isnan(lam[last,:]))[0][-1]])
    chi2last = exp(chi2[last,where(~isnan(lam[last,:]))[0][-1]])

    debug('solution existed, last chi2= %g' %exp(mean(chi2last))+' steps:%d'%i)

    
    return g, chi2last, lam_last,chi2all





def root_finder(lam,chi2, i, iimax,lambdab0,  TT,Ts,  H, fs,  g, BdMat,precision, factor, Tok, allow_search):
    """
    Search for optimal value of lambda. 
    1. option find chi2 = 1
    2. if is not possible find min(chi2)
    
    based on spline fitting of avalible points and root/fmin searching at the spline. 
    It is faster and it need less steps than any other method. 

    :param array lam:  Array of all smoothing parameters lambda used during reoconstruction
    :param array chi2:  Array of all chi2 used during reoconstruction
    :param int iimax: Maximal iteration in inner cycle
    :param int lambdab0:    Initial value of lambda
    :param spmatrix fs:   Normalised data
    :param spmatrix Ts:   Normalised geometry matrix
    :param spmatrix TT:   Ts.T*Ts
    :param spmatrix H:    Final version of derivation matrix
    :param array g: Searched final reconstruction
    :param int BdMat:  Array nonzero outside of tokamak
    :param double precision: Precision of root searching of chi2.  Precision is increased with outer iteration

    :var double lam_min, lam_max:  Range of possible values of lambda
    :var int flag: State of algorithm
    :var object factor:  Object of cholesky decomposition
    :var array weight: Array of empirical weights used to determine the best solution, if no solution exists.



    """
    
    
    lam_min = -5   #empirical constants, solution should be within this range
    lam_max = 20  
    debug('----------------------------------------')
    flag = 2
    ii = 0
    shift = Tok.shift +-log( Tok.nx*Tok.ny)*2+18


    lam[i,ii]=min(max(lambdab0,lam_min),lam_max)        #it is recommended to choose lam[i,0] > lam[i,1]
 
    try:
        try:
            t = time.time()
            factor = cholesky(1/exp(lam[i,0]+shift)*TT+H)
            #BUG inplace is slow, why???
            print('cholesky', time.time()-t)
                 

        
        #matrix is probably singular => add additional regularization by L2 norm
        except CholmodError :
            factor = cholesky(1/exp(lam[i,0]+Tok.shift -log(Ts.shape[1])*2+18)*TT+H,beta=exp(lam[i,0]-5))
        except CholmodWarning :
            factor = cholesky(1/exp(lam[i,0]+Tok.shift -log(Ts.shape[1])*2+18)*TT+H,beta=exp(lam[i,0]-5))

    except Exception as e:
        print('Exception: ',e, 'basic cholesky will by done')
    

        try:
            #problems with rounding errors??
            factor = cholesky(1/exp(lam[i,0]+Tok.shift -log(Ts.shape[1])*2+18)*TT+H)      

        except: 
            print('cholmod failured')
            factor= False
   
    from scipy.interpolate import UnivariateSpline


    del g
    

    g, chi2_all, lam_prew,neg=solve(TT,Ts,  H, fs,  lam[i, ii], BdMat, lam[i,0], factor, Tok)
    chi2[i,ii] = log(mean(chi2_all))
    
    if  abs(chi2[i,ii])<precision and not neg:
        flag =0
        lam[i, ii+1] = nan
        return chi2,lam, g,chi2_all, flag, lam[i,ii]
    
    if (neg or chi2[i,ii]<0) and lam[i,0] != lam_max:
       lam[i,ii+1] = min(lam[i,ii]+5, lam_max)
    else:
       lam[i,ii+1] = max(lam[i,ii]/2,lam_min)
       
       
    ii+=1
    del g
    g, chi2_all, lam_prew,neg=solve(TT,Ts,  H, fs,  lam[i, ii], BdMat, lam_prew, factor, Tok)
    chi2[i,ii] = log(mean(chi2_all))
    

    
    
    flag = 2
    fval = 0

    while ii  < iimax-2:
        #stop when diference between chi2[i,ii] nad prediceted fval is small enought and if even for largest 
        #lambda is still chi2<1
        try:
            if abs(chi2[i,ii]-fval)< precision and (not neg or ii>2) and fval<4\
                or ((lam[i, ii]==lam_max) and chi2[i,ii]<0) :
                flag = 0 if fval==0 else 1
                break
        except Exception as e:
            print(abs(chi2[i,ii]-fval),  precision, neg, fval, lam[i, ii], lam_max,  chi2[i,ii])
            raise e
        
        ii += 1 # ii is changed also inside the while cycle !!

        ind = argsort(lam[i, :ii])
        
        #do not print useless warnings
        import warnings
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            try:
            # Call some code that triggers a custom warning.
                S = UnivariateSpline(lam[i, ind],chi2[i,ind],k=min(2,ii-1),s=1e-5,bbox=[None,None])
            except:
                print('UnivariateSpline failured, why? ')
                import IPython
                IPython.embed()
            # ignore any non-custom warnings that may be in the list
            w = [i for i in w if issubclass(i.category, UserWarning)]

        L = linspace(lam_min,lam_max,1000)
        SL = S(L)
        

        
        
        
        xmin = L[argmin(SL)]
        fvalmin = S(xmin)
        xmax = L[argmax(SL)]
        fvalmax = S(xmax)
        xroots = []
   
        if fvalmin > 0 and fvalmax > 0 and not xmin in lam[i]:  #no root
            lam[i,ii] = xmin
            fval = fvalmin
        elif fvalmin < 0 and fvalmax < 0 and not xmax in lam[i]:  #no root
            lam[i,ii] = xmax
            fval = fvalmax
        else: #there must be some root
            S = UnivariateSpline(L,SL,k=3,s=1e-5,bbox=[lam_min,lam_max])
            xroots = S.roots()
            deriv = S(xroots,1)
            fval = 0
            if any(deriv>0):   #the expected root
                #choose the rightest root in the examined interval
                #if there any, use the righest extrapolated  root 
                roots = xroots[deriv>0]
                inner_ind = (roots>lam[i, :ii].min())&(roots<lam[i, :ii].max())
                if any(inner_ind):
                    lam[i,ii] = amax(roots[inner_ind]) 
                elif any(roots>lam[i, :ii].max()):
                    lam[i,ii] = amin(roots[roots>lam[i, :ii].max()])
                else:
                    lam[i,ii] = amax(roots)
                    
            elif ii < 3 and len(xroots):    #at the begining can be used even the wrong root
                lam[i,ii] = amax(xroots)

        #highly overestimated errorbars? 
        if xmin == lam_max and fvalmin>0 and ii  > 7:  #minimum is most probably out of the range => next step must be done
            lam[i,ii] = lam_max

            if lam_max == lam[i,ii-1]:
                flag = 1
                break
            elif lam_max in lam[i,:ii-2]:
                fval = chi2[i,lam[i,:ii-1]==lam_max]
            else:
                fval = S(lam_max)
        
        if isnan(lam[i,ii]):   #no expected root
                #place it in the  less explored region.
            # use ratio according expected distance from chi2=1
            fval = 0
            L = unique(sort(r_[lam[i,:ii], lam_max,xmin, lam_min]))
            L = L[L<=lam_max]
         
            j = argmax(diff(L))
            if L[j] in lam[i,:ii] and L[j+1] in lam[i,:ii]:  
                CH1, CH2 = chi2[i,:][lam[i,:]==L[j]][0], chi2[i,:][lam[i,:]==L[j+1]][0]
            else:
                CH1, CH2  = 1,1 
            lam[i,ii] = (abs(CH2)*L[j]+abs(CH1)*L[j+1])/(abs(CH1)+abs(CH2))   
            debug( 'no expected root, try: \tlambda=%.3f\tchi2=%.3f'%(lam[i,ii], chi2[i,ii-1]))
       

        if isnan(lam[i, ii]): #should not be possible
            print("what is wrong???")
            lam[i, ii] = random.rand(1)*(lam_max-lam_min)+lam_min
            
            
        if ii == iimax-2:  #last step - choose the best up to now known  solution

            imin = argmin(abs(chi2[i,:ii]))
            if imin == iimax-3:
                flag = 1
                break
            else:
                lam[i, ii] = lam[i, :ii][imin]
            
            
        del g  

        g, chi2_all, lam_prew,neg=solve(TT,Ts,  H, fs,  lam[i, ii], BdMat, lam_prew, factor, Tok,fval)
        chi2[i,ii] = log(mean(chi2_all))
        
        
    if isnan(lam[i,ii]):
        print('lambda is nan!!')
        import IPython
        IPython.embed() 
        
        

    return chi2,lam, g,chi2_all, flag, lam[i,ii]
        
    
    


def solve(TT, Ts,  H, fs, lam0,  BdMat, lam_prew, factor, Tok,chi2pred=nan):
    """
    Solve  of Tichononov regularization, (TT+lam0*H)*g = fs
    If it is possible,try to use sparse cholesky decomposition with updating/downdating. Otherwise use sparse solve if tsteps == 1 and ordinary solve if tsteps > 1. Finally, make emissivity zero out of boundaries and compute chi2.

    :param int lam0:    Solved value of lambda
    :param int lam_prew:    Previos value of lambda, used for cholmod updating
    :param spmatrix fs:   Normalised data
    :param spmatrix Ts:   Normalised geometry matrix
    :param spmatrix TT:   Ts.T*Ts
    :param spmatrix H:    Final version of derivation matrix
    :param array g: Searched final reconstruction
    :param int BdMat:  Array nonzero outside of tokamak

    :var double FreeMemory: Value of free RAM, used as workarounf of a bug in cholmod (MEMORY LEAK)
    :var spmatrix A:  Total matrix used to solve Ag=f

    """
    shift = Tok.shift +-log( Tok.nx*Tok.ny)*2+18
    lam = exp(double(lam0)+shift)
    lam_prew = exp(lam_prew+shift)
    nt = size(fs, 1)
    ind = all(fs!=0, axis=1)  #boundary and wrong detectors
    
    try:
        #!TODO scikits.sparse.cholmod module has important memory leaks !!!
    
        try:
            FreeMemory=os.popen("free -m").readlines()[1].split()[3]
        except:
            FreeMemory = 0
            print("failed FreeMemory")
   
   
        if not factor:
            t = time.time()

            factor = cholesky(1/lam*TT+H)
            debug('cholesky %f'%( time.time()-t))
            lam_prew = lam

        if int(FreeMemory) > 300:# and shape(Ts)[0]*4 < shape(Ts)[1] :  # low memory or not a low rank matrix !!! 
            if lam != lam_prew:
                param= lam_prew<lam
                t = time.time()
                #WARNING do not make an update of the boundary and zero channels
                factor.update_inplace(sqrt(abs(1/lam_prew-1/lam))*Ts.T, subtract=param)
                #factor.cholesky_inplace(1/lam*TT+H  )

                debug('update_inplace %f'%(  time.time()-t))

        else:
            #WORKAROUND OF BUG IN scikits.sparse.cholmod (MEMORY LEAK)
            raise Exception('not enought free memory:', int(FreeMemory) )
        
        
        
        if nt >  sum(ind):
            #it is faster to invert the geometry matrix, optimized for huge g
              #Tpsinv is sparse dense!
            Tpsinv=factor(Ts[where(ind)[0],:].T).tocsr()
            Tpsinv.data = single(Tpsinv.data/lam) 
            g = Tpsinv*fs[ind,:]
             
        else:
            g = single(factor(Ts[where(ind)[0],:].T*fs[ind,:]))/lam
        

    except Exception as e:
        debug('\nWarning: slow spsolve used -  '+str(e))
        
        #BACKUP CONFIG
        factor = factorized(TT+lam*H)   #slower way, ~3-30 times slower compared to cholesky for sparse

        if nt >  sum(ind):
            Tpsinv = empty(Ts.shape,dtype=fs.dtype)
            for j in where(ind)[0]:
                Tpsinv[j,:] = factor(array(Ts[j,:].todense())[0,:])
            g = dot(Tpsinv.T, fs)
        else:
            g = empty((Ts.shape[1],nt),dtype=fs.dtype)
            for j in range(nt):
                g[:,j] = factor(Ts[where(ind)[0],:].T*fs[ind,j])
                

    

    g = asarray(reshape(g, (len(g), -1)))
    
    neg = True
    t1 = time.time()
    chi2 = zeros(nt)
    neg_val = 0
    
    

    if not Tok.transform_index in [0,4]:
        ##fast way how to partially remove negative part of the profiles
        g_ = zeros(g.shape[0])
        if not Tok.allow_negative:
            g_ = Tok.Transform*g.mean(1)
            
            neg_val = sum(g_< 0)*nt
            g_[g_>0]=0
            g_ = Tok.Transform.T*g_
            
        chi2=sum((Ts*(g-g_[:,None])-fs)**2,0)/sum(fs!=0,axis=0)
   

    else:  # nothing 
        if Tok.transform_index == 4: g = Tok.Transform*g
        
        step = 100
        for ind in [slice(i*step,(i+1)*step) for i in range(nt//step+1)]:
            g_tmp = copy(g[:,ind])
            g_tmp[BdMat,:]=0

            if not Tok.allow_negative:
                neg_val+= sum(g_tmp<0)
                g_tmp[g_tmp<0] = 0
            chi2[ind]= sum((Ts*g_tmp-fs[:,ind])**2, axis=0)/sum(fs[:,ind]!=0,axis=0)
        
    if neg_val<0.05*g.size: neg  = False

    debug('lambda: %3.3f\t \tchi2: %2.3g\t\tchi2pred: %2.3g\t negative: %2.1f%%'\
            % (lam0, mean(chi2),exp(chi2pred),neg_val*100./nt/Tok.Transform.shape[0]))
  
      

  
    return g, chi2, lam0,neg






def plotG_chi2(TT, Ts,  H, fs,  BdMat, lam_prew, factor, Tok,lam_,chi2,g):
    #plot chi2 as function of gamma
    
    
    M = 100
    L = r_[linspace(-20,10,M ),20]
    
    CH2 = zeros(M+1)*nan
    CH2_neg = zeros(M+1)*nan
    CH2_ = zeros(M+1)*nan
    CH2_2 = zeros(M+1)*nan

    import  scipy
    import sys
    t = time.time()
    T = time.time()
    from tqdm import trange


    D =  []
    
    for k in trange(M, desc='Calculate chi2 vs.lambda: '):

        
        Tok.allow_negative = False
        g,chi2_all, lam_prew,neg = solve(TT, Ts,  H, fs, 
                        L[k],  BdMat, lam_prew, factor, Tok)
        CH2[k]= mean(chi2_all)
        Tok.allow_negative = True
        D.append(factor.D())

        g,chi2_all, lam_prew,neg = solve(TT, Ts,  H, fs, 
                        L[k],  BdMat, lam_prew, factor, Tok)

        CH2_neg[k]= mean(chi2_all)
    
        shift = Tok.shift +-log( Tok.nx*Tok.ny)*2+18
        lam0 = double(L[k])
        lam = exp(lam0+shift)

        factor_ = cholesky(TT+lam*H)
        g = factor_(Ts.T*fs)
        CH2_[k]=sum((Ts*g-fs)**2)/sum(fs!=0)
        
        g = Tok.Transform*g
        g[g<0] = 0
        g = Tok.Transform.T*g

        CH2_2[k]=sum((Ts*g-fs)**2)/sum(fs!=0)


        
    CH2[-1] =     mean(sum(fs**2, axis=0)/sum(fs!=0,axis=0))

    CH2_neg[-1] = mean(sum(fs**2, axis=0)/sum(fs!=0,axis=0))
    
    mfs = mean(fs, 1)
    
    t = time.time()
    pinvG,istop,itn,r1norm,r2norm,anorm,acond,arnorm,xnorm,var=sparse.linalg.lsqr(Ts, mfs,atol=1e-9, btol=1e-9)
    axhline(y=sum((mfs-Ts*(abs(pinvG)+pinvG)/2)**2,axis=0)/sum(mfs!=0))
    axhline(y=sum((mfs-Ts*((pinvG)+pinvG)/2)**2,axis=0)/sum(mfs!=0),ls='--')
        
                       
    #Tok.allow_negative = False
    ind = argsort(L)
    print('total time',time.time()-T)
    semilogy( L[ind],(CH2[ind]),'-',label='chi2 pos')
    semilogy( L[ind],(CH2_neg[ind]),'x',label='chi2')
    semilogy( L[ind],(CH2_[ind]),'--',label='chi2_')
    semilogy( L[ind],(CH2_2[ind]),'.-',label='chi2_ pos')

    
    semilogy(lam_,exp(chi2),'o-')
    semilogy(lam_[isfinite(lam_)][-1],exp(chi2[isfinite(lam_)][-1]),'ok')

    axhline(y=1,ls='--',c='k')
    #axhline(y=CH2[-1])

    xlabel('lambda')
    ylabel('chi2')
    legend()
        
                

    show()
    
    import IPython
    IPython.embed()
    
