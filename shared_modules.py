    #!/usr/bin/env python
# -*- coding: utf-8 -*-
#from geom_mat_setting import loadgeometry

from numpy import *
import numpy as np
from scipy import sparse
#from  matplotlib.pyplot import *
import gc, re,sys
from numpy import linalg
import time
from scipy.sparse import spdiags, eye
from scipy.sparse.linalg import spsolve
import os.path
from scipy.io import loadmat, savemat
import config
from copy import deepcopy
from scipy.signal import convolve, fftconvolve, medfilt2d, convolve2d,fftconvolve
from scipy.optimize import minimize
from annulus import get_bd_mat
from scipy.linalg import solve_triangular, sqrtm
from matplotlib.pylab import *

from scipy.interpolate import RectBivariateSpline

def prepare_inputs(Tok, data, error_all,dets, tvec, Tmat, normData, G0, danis,boundary,
                   regularization, reconstruct, postprocessing):    
    
    """
     This algorithm identify boundaries and create virtual sensors on boundary,
     prepare matrices of derivation, rescale and renormalize all its (Tmat, data, error)


    :param class tok: ``Class`` of  tokamak with all important data.
    :var double boundary:    If non-zero allow smoothing with boundary => The reconstruction 
            will be zero on the boundary,  Set the pressure forcing reconstruction to be zero on boundary. Precision of the virtual senzors (boundary) is min(sigma)/X.
    :param int regularization: Apriory smoothing used for minimalization
    :param bool reconstruct: Perform tomographic reconstruction
    :param bool postprocessing: Perform postprocessing --  entropy, smoothness
    :var int rgmin:  Small constant
    :var int noemiss:  Array nonzero outside of tokamak
    :var double danis:    Ratio of anisotropic matrices, B = n*Bperp+1/n*Bpar  => The smaller 
            number the more will reconstruction follow the field For MFI is recommended 0.2, for MDIFF < 0.1
    :param object progress: link to QObject that move with the progressbar, do not work with multiprocessing
    :var spmatrix diam,diar,dial,diao,diau: Sparse matrices of derivation (elementar directions)


    :var spmatrix Err:  Square root of variation matrix
    :var spmatrix fs:   Normalised data
    :var spmatrix Ts:   Normalised geometry matrix

    :var spmatrix H:    Final version of derivation matrix
    :var list Bmat:     List of smoothing matrices



    """


    tsteps = len(tvec)

    nx = Tok.nx
    ny = Tok.ny
    #SCALING !!! :
    normTmat = Tok.normTmat

    Tmat = deepcopy(Tmat)
    Tmat.data/=normTmat
    
 
    #remove elements of the projection matrix ot of the boundary
    from geom_mat_setting import prune_geo_matrix
    
    bd_len = int(2*pi*sqrt(nx*nx))
    bnd = Tok.get_boundary( bd_len,time=tvec.mean())
    BdMat = get_bd_mat(Tok,nx, ny,boundary=bnd)
    

 
    Tmat = prune_geo_matrix(Tmat,BdMat)

    normData = single(normData)
    data = (data/normData[:,None]).T
  
    error = ones(size(error_all,1))*infty
     #BUG all or  any???? !!!
    ind = all(isfinite(error_all)&(error_all!= 0),0)

    error[ind]= mean(normData)/mean(normData[:,None]/error_all[:,ind],0)  #in rapid version use mean error

    mag_axis = None
    #compute averadge position of the magnetic axis
    if hasattr(Tok, 'mag_axis') and Tok.mag_axis != None:
        T = linspace(tvec.min(), tvec.max(),100)
        rc = ones_like(tvec)*mean(interp(T,Tok.mag_axis['tvec'],Tok.mag_axis['Rmag']))
        zc = ones_like(tvec)*mean(interp(T,Tok.mag_axis['tvec'],Tok.mag_axis['Zmag']))
        mag_axis = c_[rc,zc].T
        

    error /= mean(normData)
 
    Tok.Transform  = Tok.generate_transform(tvec)
 
    try:
        BdMat = abs((Tok.Transform.T*BdMat.T)) > .5   #Tok.Transform can have problems with __rmul__
    except:
        print('BdMat = (Tok.Transform.T*BdMat.T).T > .1', shape(Tok.Transform),Tok.npix,Tok.nx*Tok.ny)
        
        import IPython
        IPython.embed()
        raise

    Tmat = (Tok.Transform.T*Tmat.T).T  #Tok.Transform can have problems with __rmul__ 

    #inverse transformation of F
    G0 = Tok.Transform.T*G0


    G0 = reshape((G0*normTmat*mean(1/normData)), (Tok.npix, 1), order='F')                     #initial guess of G

    

    from mat_deriv_B import mat_deriv_B

    Bmat, diag_mat = mat_deriv_B(Tok, tvec,regularization, danis)
    error = single(error)


    if boundary>0 and Tok.transform_index==0:  #if boundary is nonzero, uses the value as precision, recommended value is 50-100

        BdMat_ = BdMat.reshape(Tok.ny, Tok.nx, order="F")
        from scipy.ndimage.morphology import binary_dilation
        
        border = binary_dilation(~BdMat_)&BdMat_

        ind_border = where(border.flatten(order='F'))[0]
        n_bord = len(ind_border)
        Tadd = sparse.coo_matrix((ones(n_bord), ( ind_border ,arange(n_bord))), (nx*ny, n_bord))

        Tadd = Tok.Transform.T*Tadd

        Tmat = sparse.vstack((Tmat,Tadd.T))                # add artificial (virtual) sensors at the tokamak boundary
        data = vstack((data,zeros((n_bord,tsteps),dtype=data.dtype)))
        error = hstack((error,amin(error)*ones(n_bord,dtype=data.dtype)/sinh(boundary))) 
        error = single(error)

    nl,npix = shape(Tmat)
     

    iErr = spdiags(1/error, 0, nl, nl)
    Err = spdiags(error, 0, nl, nl)

    data /= error[:,None]
    fs = data
    Ts = iErr*Tmat
    
    return Bmat,BdMat,bnd, Ts, fs,Err, tvec,G0






def create_derivation_matrix(g, Bmat,BdMat, regularization, danis, Tok, compute_2nd_deriv = True ):
    """
    Prepare matrix of derivation according to choosed regularization method
    """
 

    # prefere stability , robusness        
    g_tmp = mean(asarray(g),1)
    g_tmp = g_tmp if Tok.transform_index==0 else Tok.Transform * g_tmp
    
    npix = len(g_tmp)
    max_g = nanmax(abs(g_tmp))

    g_tmp /= max_g
    
    #print('Tok.boundary and Tok.transform_index==0',Tok.boundary , Tok.transform_index==0)
    if Tok.boundary and Tok.transform_index==0:
        g_tmp[BdMat] = Tok.rgmin

    g_tmp[g_tmp < Tok.rgmin] = Tok.rgmin
    

    g_tmp_inv = 1/g_tmp #**Tok.weight_power

    if Tok.camera and Tok.allow_reflections:
        g_tmp_inv *= Tok.correction_matrix  # only in case of camera try to remove artifacts near to the camera automatically

    W=spdiags(g_tmp_inv,0, npix,npix,format='csr')

    if regularization  in [0,5]:  # izotropic MFI
        H = Bmat[0].T*(W*Bmat[0]) +  Bmat[1].T*(W*Bmat[1])
        
    elif regularization in [1,6]:# anizotropic MFI
        H = sigmoid(-danis)*( Bmat[0].T*(W*Bmat[0]) + Bmat[2].T*(W*Bmat[2]))\
            +sigmoid(danis) *(Bmat[1].T*(W*Bmat[1]) +Bmat[3].T*(W*Bmat[3]))
        

    elif regularization  in [2]:  # izotropic MDIFF
        if compute_2nd_deriv:
            Bmat = [b.T*b for b in Bmat]
        H = Bmat[0].T*(W*Bmat[0]) +  Bmat[1].T*(W*Bmat[1])

    elif regularization in [3]:# anizotropic MDIFF
        if compute_2nd_deriv:
            Bmat = [b.T*b for b in Bmat]
        H = sigmoid(-danis)*( Bmat[0].T*(W*Bmat[0]) + Bmat[2].T*(W*Bmat[2]))\
            +sigmoid(danis) *(Bmat[1].T*(W*Bmat[1]) +Bmat[3].T*(W*Bmat[3]))

    elif regularization in [4]:  # ENERGY
        H = Bmat[0]
        
    elif regularization in [7]: #INGESSON anisotropic
        H = Bmat[0].T*(W*Bmat[0])+Bmat[1]
        
    else:
        raise Exception('regularization not found ')
    
    if regularization in [5,6]:  #squared
        H = H.T*H

  
    if Tok.beta!= 0 :  # press reconstruction to zero 
        H = H + Tok.beta* spdiags(ones(npix),0, npix,npix)* W.data.min()
        

    if Tok.transform_index==0:
        return H
    else:
        #slow!!
        return (Tok.Transform.T*(Tok.Transform.T * H).T).T



def make_postprocessing(output,Tok, G,data, error,retro, tvec, dets,SigmaGsample,input_parameters):
    """
    Perform postprocessing of reconstructed emissivity 

    """
    debug('postprocessing started')
    T = time.time()
    


    resid = retro-data
    resid/= error  #BUG do not correcr when covariance matrix is used
    
    
    xcpix=linspace(Tok.xmin+Tok.dx/2,Tok.xmax+Tok.dx/2, Tok.nx)/Tok.norm
    ycpix=linspace(Tok.ymin+Tok.dy/2,Tok.ymax+Tok.dy/2, Tok.ny)/Tok.norm 
    
    
    tsteps = len(tvec)

    G = G.reshape(Tok.ny, Tok.nx, tsteps, order='F')
    chi2_real = mean(resid**2,1)
    
    output['chi2_real'] = chi2_real
    
    dV = 2*pi*xcpix*Tok.dx*Tok.dy  #volume of the single pixel
    output['power'] = einsum('ijk,j ->k',G, dV)
    
    rho,magr,magz=Tok.mag_equilibrium(tvec.mean(0),surf_slice=[-1,],rho=1)
    coords = hstack((magr,magz))[:,:,0]

    BdMat = get_bd_mat(Tok,Tok.nx, Tok.ny,boundary=coords)
    BdMat = BdMat.reshape(Tok.ny, Tok.nx, order='F')
    output['power_div'] = einsum('ijk,j ->k',G*BdMat[:,:,None], dV)

    position = zeros((tsteps,2))


    
    if 'magx' in output:
        magx = output['magx'].T
        magy = output['magy'].T
    else:
        rho,magx, magy = Tok.mag_equilibrium(tvec, return_mean=False) #SLOW!
        output['magx'] = magx.T
        output['magy'] = magy.T
        output['mag_rho'] =  rho[None]
    
    if hasattr(Tok,'mag_axis'):
        position[:,0] = interp(tvec,Tok.mag_axis['tvec'],Tok.mag_axis['Rmag'])+Tok.magfield_shift_core[0]
        position[:,1] = interp(tvec,Tok.mag_axis['tvec'],Tok.mag_axis['Zmag'])+Tok.magfield_shift_core[1]
        
    else:        
        position[:,0] = mean(magx[:,0,:],0)/Tok.norm
        position[:,1] = mean(magy[:,0,:],0)/Tok.norm
        Tok.mag_axis = {'tvec':tvec,'Rmag':position[:,0],'Zmag':position[:,1]}
        
    output['position'] = position 

    #compute emisivity "moments" from magnetic
    dS = zeros_like(magx)  #volume of every elementary cell in 
    for it in range(tsteps): 
        #usage of for loop will save a lot of memory
        #it should be correct with subpixel precision
        #BUG is it OK? was toroidal expansion included?? 
        diag_r_x   = magx[:,2:,it]-magx[:,:-2,it]
        diag_r_y   = magy[:,2:,it]-magy[:,:-2,it]

        diag_phi_x = r_[-magx[(0,),1:-1,it]+magx[(2,),1:-1,it],magx[2:,1:-1,it]-magx[:-2,1:-1,it],magx[(1,),1:-1,it]-magx[(-2,),1:-1,it]]
        diag_phi_y = r_[-magy[(0,),1:-1,it]+magy[(2,),1:-1,it],magy[2:,1:-1,it]-magy[:-2,1:-1,it],magy[(1,),1:-1,it]-magy[(-2,),1:-1,it]]

        norm_r_phi = hypot(diag_r_y,diag_r_x)*hypot(diag_phi_x,diag_phi_y)
        cos_phi = (diag_r_x*diag_phi_x+diag_r_y*diag_phi_y)/norm_r_phi
    
        dS[:,1:-1,it]  = .5*sqrt(fabs(1-cos_phi**2))*norm_r_phi/2
        
        
    dS[:,0,:]  = dS[:,1,:]/10 #aproximate the core 
    dS[:,-1,:]  = 2*dS[:,-2,:]-dS[:,-3,:] #aproximate the edge 
    dV = 2*pi*dS*magx*dS

    # EMISSIVITY PROFILE 
    n_mag = size(magx,1)
    fsa_emiss = zeros((n_mag, tsteps),dtype=single)
    
    scaling = array([Tok.dx,Tok.dy])
    offset = array([Tok.xmin,Tok.ymin])#+scaling/2
    
    from scipy.ndimage.interpolation import map_coordinates    
    
    for i in range(tsteps): 
        coords = c_[magx[...,i].ravel(),magy[...,i].ravel()].T
        idx = (coords-offset[:,None])/ scaling[:,None]
        map_prof = map_coordinates(G[...,i].T,idx,order=1)
        map_prof = map_prof.reshape(magy[...,i].shape)
        fsa_emiss[:,i]   = average(map_prof,0,dV[...,it])

    if SigmaGsample is not None:
        n_sample =  SigmaGsample.shape[-1]
        SigmaGsample = SigmaGsample.reshape(Tok.ny, Tok.nx,n_sample, order='F')
        fsa_emiss_err = zeros((n_mag, n_sample),dtype=single)
        coords = c_[magx.mean(-1).ravel(),magy.mean(-1).ravel()].T
        idx = (coords-offset[:,None])/ scaling[:,None]
        dS_m = dV.mean(-1)
        for i in range(n_sample): 
            map_prof = map_coordinates(SigmaGsample[...,i].T,idx,order=1)
            map_prof = map_prof.reshape(magy[...,0].shape)
            fsa_emiss_err[:,i] = average(map_prof,0,dS_m)

        output['fsa_emiss_covar'] = cov(fsa_emiss_err)[None]

        fsa_emiss_err = std(fsa_emiss_err,1)
        output['fsa_emiss_err'] = tile(fsa_emiss_err,(tsteps,1))

    
    profmax = fsa_emiss.max(0)[None,:]

    profile_0_3 = copy(fsa_emiss)
    profile_0_3/= profmax

    #devide the profile in the the 3 parts accrding teh radiated power with smooth transition
    f0 = lambda x: interp(x, (0,.2,.4),(1,1,0)) 
    f1 = lambda x: interp(x, (.2,.4,.6,.8),(0,1,1,0)) 
    f2 = lambda x: interp(x, (.6,.8,1),(0,1,1))
    
    profile_2_3 = f2(profile_0_3)
    profile_1_3 = f1(profile_0_3)
    profile_0_3 = f0(profile_0_3)


    eps = 1e-10
    norm_0_3 = einsum('ijk,jk->k',dV,profile_0_3)+eps
    Rmag_0_3 = einsum('ijk,ijk,jk,k->k',dV,magx, profile_0_3,1/norm_0_3)
    Zmag_0_3 = einsum('ijk,ijk,jk,k->k',dV,magy, profile_0_3,1/norm_0_3)

    norm_1_3 = einsum('ijk,jk->k',dV,profile_1_3)+eps
    Rmag_1_3 = einsum('ijk,ijk,jk,k->k',dV,magx, profile_1_3,1/norm_1_3)
    Zmag_1_3 = einsum('ijk,ijk,jk,k->k',dV,magy, profile_1_3,1/norm_1_3)

    norm_2_3 = einsum('ijk,jk->k',dV,profile_2_3)+eps
    Rmag_2_3 = einsum('ijk,ijk,jk,k->k',dV,magx, profile_2_3,1/norm_2_3)
    Zmag_2_3 = einsum('ijk,ijk,jk,k->k',dV,magy, profile_2_3,1/norm_2_3)

    output['fsa_emiss'] = fsa_emiss.T
    output['xmass_eq'] = vstack((Rmag_0_3,Rmag_1_3,Rmag_2_3)).T
    output['ymass_eq'] = vstack((Zmag_0_3,Zmag_1_3,Zmag_2_3)).T
    

    Gmax = amax(amax(G,0),0)[None,None,:]

    xmass = empty((3,tsteps),dtype=single) 
    ymass = empty((3,tsteps),dtype=single) 

    n_chunk = 100
    for i,fun in enumerate([f0,f1,f2]):
        for it in range(tsteps//n_chunk+1):
            ind = slice(it*n_chunk, (it+1)*n_chunk)
            G_ = fun(G[:,:,ind]/(Gmax[:,:,ind]+eps))
            g_sum = einsum('ijk->k',G_)+eps
            #BUG!!!correction factors dx/2 was found after detail comparism with phantom, 
            #but cause of this subpixel error is unknown
            xmass[i,ind] = einsum('ijk,j,k->k',G_,xcpix-Tok.dx/2,1/g_sum) 
            ymass[i,ind] = einsum('ijk,i,k->k',G_,ycpix-Tok.dy/2,1/g_sum)

    output['xmass'] = xmass.T
    output['ymass'] = ymass.T



    debug('postprocessing done at %.2fs'%(time.time()-T))

    return output
    
    
def in1d(ar1, ar2, assume_unique=False):
    """
    Test whether each element of a 1D array is also present in a second array.
    """
    if not assume_unique:
        ar1, rev_idx = unique(ar1, return_inverse=True)
        ar2 = unique(ar2)

    ar = concatenate( (ar1, ar2) )
    # We need this to be a stable sort, so always use 'mergesort'
    # here. The values from the first array should always come before
    # the values from the second array.
    order = ar.argsort(kind='mergesort')
    sar = ar[order]
    equal_adj = (sar[1:] == sar[:-1])
    flag = concatenate( (equal_adj, [False] ) )
    indx = order.argsort(kind='mergesort')[:len( ar1 )]

    if assume_unique:
        return flag[indx]
    else:
        return flag[indx][rev_idx]


def unique(ar, return_index=False, return_inverse=False):
    """
    Find the unique elements of an array.

    Returns the sorted unique elements of an array. There are two optional
    outputs in addition to the unique elements: the indices of the input array
    that give the unique values, and the indices of the unique array that
    reconstruct the input array.
    """
    try:
        ar = ar.flatten()
    except AttributeError:
        if not return_inverse and not return_index:
            items = sorted(set(ar))
            return np.asarray(items)
        else:
            ar = np.asanyarray(ar).flatten()

    if ar.size == 0:
        if return_inverse and return_index:
            return ar, np.empty(0, np.bool), np.empty(0, np.bool)
        elif return_inverse or return_index:
            return ar, np.empty(0, np.bool)
        else:
            return ar

    if return_inverse or return_index:
        perm = ar.argsort()
        aux = ar[perm]
        flag = np.concatenate(([True], aux[1:] != aux[:-1]))
        if return_inverse:
            iflag = np.cumsum(flag) - 1
            iperm = perm.argsort()
            if return_index:
                return aux[flag], perm[flag], iflag[iperm]
            else:
                return aux[flag], iflag[iperm]
        else:
            return aux[flag], perm[flag]

    else:
        ar.sort()
        flag = np.concatenate(([True], ar[1:] != ar[:-1]))
        return ar[flag]

def debug(*string):
    string = [str(s) for s in string]
    if config.DEBUG:
        
        print(''.join(string))
        
        
def warning(*warning ):
    sys.stderr.write('\x1b[1;33m' + ''.join(warning)  + '\x1b[0m\n')
        
 

def sigmoid(x):    return 1 / (1 + exp(-x))

from scipy import signal

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    ker = max(int(size),1)
    if not sizey:
        sizey = size
        kery = ker
    else:
        kery =  int(sizey)
        sizey = float(sizey)
    x, y = np.mgrid[-ker:ker+1, -kery:kery+1]
    g = np.exp(-(x**2/float(size) + y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    #improc = fftconvolve(im, g, mode='same') 

    improc = convolve2d(im, g, mode='same', boundary="symm") 
    return improc



def delete_sparse(A, ind, axis = 0):

    if dtype(ind[0]) == "bool":
        ind = where(ind)[0]

    [indx, indy] = A.nonzero()

    for i in where(in1d(A.nonzero()[axis], ind))[0]:
        A[indx[i],indy[i]] = 0

    return A

def MAD(x):
    return mean(abs(x-mean(x)))
    
 
def nandot(X,Y):       

    a,b = shape(X)
    c,d = shape(Y)
    if b != c:
        print('wrong size',X.shape,Y.shape)
        exit()
    
    Z = empty((a,d))
    
    for i in range(a):
        Z[i,:] = nanmean(Y*X[i,:][:,newaxis],axis=0)*b
    return Z

    
def fsvd(A, k=10, i=3, usePowerMethod=False):
    #BUG!! sometimes are not working!!!
## FSVD Fast Singular Value Decomposition 
# 
#   [U,S,V] = FSVD(A,k,i,usePowerMethod) computes the truncated singular
#   value decomposition of the input matrix A upto rank k using i levels of
#   Krylov method as given in [1], p. 3.
# 
#   If usePowerMethod is given as true, then only exponent i is used (i.e.
#   as power method). See [2] p.9, Randomized PCA algorithm for details.
#   can be slighly faster and lower memory consumption, but quality will be worse
#   
#  WARNING higher i can lead to problems with singledouble precison
# 
#   [1] Halko, N., Martinsson, P. G., Shkolnisky, Y., & Tygert, M. (2010).
#   An algorithm for the principal component analysis of large data sets.
#   Arxiv preprint arXiv:1007.5510, 0526. Retrieved April 1, 2011, from
#   http://arxiv.org/abs/1007.5510. 
#   
#   [2] Halko, N., Martinsson, P. G., & Tropp, J. A. (2009). Finding
#   structure with randomness: Probabilis algorithms for constructing
#   approximate matrix decompositions. Arxiv preprint arXiv:0909.4061.
#   Retrieved April 1, 2011, from http://arxiv.org/abs/0909.4061.
# 
#   See also SVD.
# 
#   Copyright 2011 Ismail Ari, http://ismailari.com.

    from scipy.sparse import issparse
    import numpy as np
    from  scipy.linalg import qr,svd

    isSparse =  issparse(A) 
    if not isSparse:  
        A = np.asmatrix(A)


    # Take (conjugate) transpose if necessary. It makes H smaller thus
    # leading the computations to be faster
    if A.shape[0] > A.shape[1]:
        A = A.T 
        isTransposed = True 
    else:
        isTransposed = False 


    m,n = A.shape
    l = k + 5 
    maxA = A.max()
    A = A/maxA

    # Form a real n×l matrix G whose entries are iid Gaussian r.v.s of zero
    # mean and unit variance
    G = np.random.randn(n,l).astype(A.dtype,copy=False)


    if usePowerMethod:
        # Use only the given exponent
        H = A*G 
        for j in range(i):
            H = A * (A.T*H) 

    else:
        # Compute the m×l matrices H^{(0)}, ..., H^{(i)}
        # Note that this is done implicitly in each iteration below.
        
        #if isSparse:
            #H = csc_matrix((m,l * (i+1)))
        #else:
        H = np.zeros((m, l * (i+1)), dtype=A.dtype) 
        
        H[:,:l] = A*G

        for j in range(2,i+2):
            H[:,(j-1)*l: j*l] = A *(A.T*H[: , (j-2)*l:(j-1)*l])

        #if  isSparse and H.size  > 0.3*np.prod(H.shape):  # sparsity
            #H = H.todense()


        # Form the m×((i+1)l) matrix H
       

    # Using the pivoted QR-decomposiion, form a real m×((i+1)l) matrix Q
    # whose columns are orthonormal, s.t. there exists a real
    # ((i+1)l)×((i+1)l) matrix R for which H = QR.  
    # XXX: Buradaki column pivoting ile yapılmayan hali.

    Q,_ = qr(H,mode='economic', overwrite_a=True, check_finite=False)

    del H
    

    # Compute the n×((i+1)l) product matrix T = A^T Q
    T = A.T*Q
    
    # Form an SVD of T
    Vt, St, W = svd(T,full_matrices=False,overwrite_a=True,check_finite=False)

    # Compute the m×((i+1)l) product matrix
    Ut = np.dot(Q,W) 
    

    # Retrieve the leftmost m×k block U of Ut, the leftmost n×k block V of
    # Vt, and the leftmost uppermost k×k block S of St. The product U S V^T
    # then approxiamtes A. 

    if isTransposed:
        V = Ut[:,:k]
        U = Vt[:,:k]    
    else:
        U = Ut[:,:k] 
        V = Vt[:,:k] 
        

    return U,St[:k]*maxA,V

    
    
    
def fast_svd(M, k):  # random projection SVD 
    import numpy.linalg as linalg
    from  scipy.linalg import qr,svd


    transpose = False
    if 3*M.shape[1]<M.shape[0]:  #Choose a optimal side in order to make it faster
        transpose = True
        M = M.T
 
    p = k+5
    R = np.random.normal(size=(p,M.shape[1])).astype(M.dtype)

    
    from scipy.sparse  import issparse
    
    if issparse(M): 
        Y = M.dot(R.T)
    else:
        M = ascontiguousarray(M)
        Y = dot(R ,M.T).T
    
    

    Q,r = qr(Y,mode='economic', overwrite_a=True, check_finite=False)
    
    if issparse(M): 
        B = M.T.dot(Q).T
    else:
        M = asfortranarray(M)
        B  = np.dot(Q.T,M)
        
    Uhat, s, v = svd(B,full_matrices=False,overwrite_a=True,check_finite=False)
    U = np.dot(Q, Uhat)
    
    if  transpose:
        U,s,v = v[:k,...].T, s[:k], U[...,:k].T
    else:
        U,s,v = U[...,:k], s[:k], v[:k,...] 
        
    return U,s,v



def svd_qr(A, tol = 1e-8, N = None): #rank revealing SVD
    #TODO otestovat vliv C nebo F poradi c matici A na rychlost vypoctu
    from  scipy.linalg import svd,qr

    swap = lambda x:x
    if size(A,0) < size(A,1):
        swap = transpose
        
    A = swap(A)
    #if size(A,0) < size(A,1):
        #transpose = True
        #A = A.T
    #else:
        #transpose = False

    q,r = qr(A,overwrite_a=False, mode='economic', pivoting=False)  #BUG jde lépe, bez výpočetu Q!!?
    u, s, vt = svd(r,full_matrices = False, overwrite_a = True )

    #r, = qr(A, overwrite_a=False,mode='r',  check_finite=False)
    #r = r[:,:k]
    #u,s,v = svd(r, full_matrices=False, overwrite_a=True, check_finite=False)
    #invR = dot(v.T, u.T/s[:,None] )
    #vt = dot(invR,A )
    
    

    index = slice(0, N)
    if N is None:
        index = s/s[0] > tol



    #if transpose:
        #V = (dot(q, u)[:,index]).T
        #U = vt[index,:].T
    #else:
    U = swap(dot(q, u)[:,index])
    V = swap(vt[index,:])
        
    S = s[index]


    return U,S,V


def MovingAveradge(sig, n,axis=-1):
    #Fast algorithm for calculation of the moving averadge
    #can failure due to rounding errors!!
    #print axis
    if n <= 1: return sig
    sig = sig.swapaxes(axis, -1)
    n = int(n)


    sig.cumsum(axis=-1,out=sig)
    right_pad  = ones(n//2, dtype=sig.dtype)*sig[...,-1][...,None]
    left_pad = zeros(shape(sig[...,0])+((n+1)//2,), dtype=sig.dtype)
    

    cs = concatenate((sig[...,n//2:],right_pad), axis=-1)
    cs-= concatenate((left_pad,sig[...,:-n//2]), axis=-1)
    cs *=1./n
    edge = 1-arange(n//2+1.)/n
    cs[...,:(n+1)//2] /= edge[-2+n%2::-1]
    cs[...,-(n+1)//2:]/= edge
    cs.swapaxes( -1,axis)
    

    return cs

def MovingAveradgeFast(sig, n,axis=-1):
    #Fast inplace algorithm for calculation of the moving averadge
    #return only data from n/2 to end-n/2!!!
    #can failure due to rounding errors!!

    sig = sig.swapaxes(axis, -1)

    sig.cumsum(axis=-1,out=sig)

    sig[...,:-n]-= sig[...,n:]
    sig = sig[...,:-n]

    sig *=-1./n

    sig =  sig.swapaxes( -1,axis)
            
    
    return sig


def denoise(img, sigma, distr= 'gauss2', boundary = False, boundary_mask= None):

    N = 3*sigma

    if boundary_mask is None:
        boundary_mask = ones(shape(img))

    if 'lorenz' ==  distr:
        N= ceil(6*sigma)
    

    x = arange(-N,N)/sigma
    [X,Y] = meshgrid(x,x)

    if 'gauss' == distr:
        a = exp(x)
    elif 'gauss2' == distr:
        a = exp( -x**2 )
    elif 'lorenz'== distr:
        mask = 1/pi / (1 + X**2 + Y**2 )
    elif 'laplace'== distr:
        a = exp( -x**2 )
        L =  -ones((3,3))
        L[1,1] = 8
    elif 'sobel'== distr:
        a = exp( -x**2 )
        L1 = array([-1,0,1-2,0,2-1,0,1]) /8.
        L2 = L1.T
        mask = L1
    else:
        print('unknown smoothing')



    if not 'mask' in locals():
        [X,Y] = meshgrid(a,a)
        mask = X*Y

    mask = mask/sum(mask[:])

    I_size = array(shape(img))
    M_size = array(shape(mask))
    outsize = M_size + I_size - 1

    # figure out when to use spatial domain vs. freq domain
    conv_time = time_conv2(I_size,M_size) # 1 conv2
    fft_time = 3*time_fft2(outsize) # 2 fft2 + 1 ifft2


    if (conv_time < fft_time):
         conv = convolve
    else:
         conv = fftconvolve
    

    if not 'sobel' ==  distr:
        img_new = conv(img, mask, 'same')
    

    if boundary:
        b = conv(boundary_mask, mask, 'same')
        img_new = img_new / b
        img_new[boundary_mask == 0] = 0


    if 'laplace' ==  distr:
        img_new = conv(img_new,L, 'same')
        b = conv(ones(shape(img)), L, 'same')
        img_new[ abs(b) > 1e3*median(b[:])]  = 0


    if 'sobel' == distr:
        a = conv(img,L1, 'same')
        b = conv(img,L2, 'same')
        img_new = sqrt(a**2 + b**2)

    return img_new


#####x other ==================================
#def nanmean(x, axis=None):
    #return nansum(x, axis = axis) / sum(~isnan(x), axis = axis)

#def nanstd(x, axis=None):
    #return std(x[~isnan(x)])


#def nanmedian(x, axis=None):
    #return median(x[~isnan(x)], axis = axis)
 
 
def read_config(path, file = ""):
    """Read data configuration from the fileobject and return it as a dictionary"""
    try:
        import configparser
    except:
        import ConfigParser as configparser

    #from collections import OrderedDict
    config = configparser.RawConfigParser()
    config_dict = dict()
    config.readfp(open(path+file, 'r'))
    for data_type in config.sections():
        config_dict.update(dict(config.items(data_type)))
        
    for key,item in config_dict.items():
        if 'True' in item: item = True
        elif 'False' in item: item = False
        elif '[]' in item: item = []
        elif 'None' in item: item = None
        elif isscalar(item): 
            item = re.sub("[\[\]]", "", item) 
            item = item.split(",")
            typ = [int, float, str]
            for t in typ:
                try:
                    item = [t(i) for i in item ]
                except:
                    continue
                break
     

            if len(item) == 1: item = item[0]
        config_dict[key] = item
    return config_dict




#TODO dodělat bolometrii centra plazmatu
def smooth(x,window_len=11,window='hanning',axis=-1):
    
    from scipy import signal 

    if not x.ndim in (1,2):
        raise ValueError("smooth only accepts 1 or 2 dimension arrays.")
    if x.shape[axis] < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
        
    x = atleast_2d(x)
    x = swapaxes(x, -1, axis)
    y = empty_like(x)
    
    
    if isinstance(window, str) or type(window) is tuple:
        w = signal.get_window(window, window_len)
    else:
        w = asarray(window)
        if len(w.shape) != 1:
            raise ValueError('window must be 1-D')
        if  w.shape[0] > x.shape[-1]:
            raise ValueError('window is longer than x.')
        window_len = w.shape[0]
    
    if window_len>20:
        conv = signal.fftconvolve
    else:
        conv = signal.convolve
    for i in range(x.shape[0]):
        s=r_[2*x[i,0]-x[i,window_len-1::-1],x[i,:],2*x[i,-1]-x[i,-1:-window_len:-1]]
        
        y[i,:]=conv(s,w/w.sum(),mode='same')[window_len:-window_len+1]     
    
    y = swapaxes(y, -1, axis)

    return squeeze(y)








def find_crossection(lineA,lineB,extrapolate=0):
    #find crosssection of two lines
    #extrapolate extend the second line 
    if extrapolate > 0:
        lineB = concatenate((((lineB[0]-lineB[1])*extrapolate+lineB[0])[None,:], 
                lineB, ((lineB[-1]-lineB[-2])*extrapolate+lineB[-1])[None,:]))

    lineA_ = lineA[None,:,:]
    lineB_ = lineB[:,None,:]
    DA = diff(lineA_,axis=1)
    DB = diff(lineB_,axis=0)
    B = (lineB_-lineA_)[:-1,:-1]
    q  = DA[:,:,0]*DB[:,:,1]-DA[:,:,1]*DB[:,:,0]
    x1 = (DA[:,:,1]*B[:,:,0]-DA[:,:,0]*B[:,:,1])/q
    x2 = (DB[:,:,1]*B[:,:,0]-DB[:,:,0]*B[:,:,1])/q
    ind = array(where(((x1>=0)&(x1<=1)&(x2>=0)&(x2<=1)))).T[:,::-1]
    crosssection = [x2[iB,iA]*DA[0,iA]+lineA[iA] for iA,iB in ind]
    positionA = [x2[iB,iA]+iA for iA,iB in ind]
    positionB = [x1[iB,iA]+iB for iA,iB in ind]

    return crosssection,positionA,positionB
    


def inside_convex_hull(boundary,x):
    """x.shape = (...,2), boundary.shape = (:,2), boundary must have the clock wise direction
    it is !!not!! working for nonconvex boundaries"""
    def cross(o, a, b): 
        return (a[:,0] - o[:,0]) * (b[...,1,None] - o[:,1]) - (a[:,1] - o[:,1]) * (b[...,0,None] - o[:,0])
    inlayers = np.all(cross(boundary[:-1,:], boundary[1:,:],x) > 0,axis=-1)
    return inlayers



 
#http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html

class fitEllipse:
    def __init__(self,x,y):
        x = x[:,np.newaxis]
        y = y[:,np.newaxis]
        D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
        S = np.dot(D.T,D)
        C = np.zeros([6,6])
        C[0,2] = C[2,0] = 2; C[1,1] = -1

        E, V =  linalg.eig(linalg.solve(S,C))
        self.a = real(V[:,3])


    def ellipse_center(self):
        a = self.a
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        num = b*b-a*c
        x0=(c*d-b*f)/num
        y0=(a*f-b*d)/num
        return np.array([x0,y0])

    def ellipse_angle_of_rotation(self ):
        a = self.a
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        return 0.5*np.arctan(2*b/(a-c))

    def ellipse_axis_length( self ):
        a = self.a
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
        down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        res1=np.sqrt(abs(up/down1))
        res2=np.sqrt(abs(up/down2))
        return np.array([res1, res2])
    
    
    
    

 
 
from scipy.interpolate import interp1d
def poly_min(x,dx,y): #find extrem of polynom 
    return x+(y[0]-y[2])*dx/(2*(y[0]+y[2]-2*y[1]))

def min_fine(x,y): #find extrem with subpixel precision
    i = argmin(y)
    if i==0 or i == len(y)-1: return y[i],x[i]    #not working at the edge! 
    A = x[i-1:i+2]**2, x[i-1:i+2], ones(3)
    a,b,c = linalg.solve(array(A).T,y[i-1:i+2])
    
    ymin =  c-b**2/(4*a)
    xmin = -b/(2*a) 
    return  ymin,xmin


    
def extrap1d(xs,ys):

    interpolator = interp1d(xs,ys,axis=0,assume_sorted=True, copy=False)
    def extrapol(xo):
        #just for 1D x and 0. axis of Y
        
        yo = empty(xo.shape+ys.shape[1:])

        # Values lower than the minimum "x" are extrapolated at the same time
        low = xo < xs[0]
        yo[low] =  squeeze(ys[0] + (xo[low,None]-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0]))

        # Values higher than the maximum "x" are extrapolated at same time
        high = xo > xs[-1]
        yo[high] = squeeze(ys[-1] + (xo[high,None]-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2]))

        # Values inside the interpolation range are interpolated directly
        inside = (xo >= xs[0])&(xo <= xs[-1])
        yo[inside] = interpolator(xo[inside])
        return yo
  

    return  extrapol



from scipy.interpolate import NearestNDInterpolator, LinearNDInterpolator

def get_rho_tan(rhop, magx, magy,XC, YC, npoints=100):

    Rc,Zc = magx[:,0].mean(), magy[:,0].mean()

    magx_ext = (magx[:,-1]-Rc)*10+Rc
    magy_ext = (magy[:,-1]-Zc)*10+Zc

    MagXY = c_[r_[Rc,c_[magx,magx_ext].ravel()] ,r_[Zc,c_[magy,magy_ext].ravel()]].T

        
    rho_mag = outer(ones(magx.shape[0]),r_[rhop,10])

    LNDI = LinearNDInterpolator(MagXY.T, r_[0,rho_mag.ravel()],fill_value=10) #14ms            

    imin = minimum(argmin(hypot(XC-Rc, YC-Zc),0), len(XC)-2)
    
    XC1 = XC[imin  , arange(len(imin))]
    XC2 = XC[imin+1, arange(len(imin))]
    YC1 = YC[imin  , arange(len(imin))]
    YC2 = YC[imin+1, arange(len(imin))]
    
    
    xmax = magx[:,-1].max()
    xmin = magx[:,-1].min()
    ymax = magy[:,-1].max()
    ymin = magy[:,-1].min()
    x = array((xmin,xmax))[:,None]
    y = (YC2-YC1)/(XC2-XC1)*(x-XC1)+YC1
    y = minimum(maximum(ymin,y ), ymax)
    x = (XC2-XC1)/(YC2-YC1)*(y-YC1)+XC1
    x = (x[0]+x[1])/2+.5*linspace(-1,1,npoints)[:,None]*(x[1]-x[0])*sign(XC2-XC1)
    y = (YC2-YC1)/(XC2-XC1)*(x-XC1)+YC1
    
    
    rho = LNDI(array((x,y)).T)

    cos_theta = zeros(len(rho))
    rho_tangent = zeros(len(rho))
    Rtan = zeros(len(rho))
    Ztan = zeros(len(rho))

    irho = arange(len(rho.T))
    for ilos,r  in enumerate(rho):
        rhomin,i_min = min_fine(irho,r)
        rho_tangent[ilos] = rhomin
        Rtan[ilos] = interp(i_min, arange(len(x)),x[:,ilos])
        Ztan[ilos] = interp(i_min, arange(len(x)),y[:,ilos]) 
        cos_theta[ilos] = (Rtan[ilos]-Rc)*(YC2[ilos]-YC1[ilos])-( XC2[ilos]-XC1[ilos])*(Ztan[ilos]-Zc)


    rho_tangent*=sign(cos_theta)


    return rho_tangent, x,y
        
    
