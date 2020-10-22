#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import *
#from  matplotlib.pyplot import *
from  scipy import sparse
from annulus import get_rho_field_mat
from mat_deriv_B import calculate_theta_angle, get_coord
from  scipy.linalg import *
from scipy import sparse
from shared_modules import delete_sparse, blur_image, debug
from scipy.interpolate import interp2d
import time
from scipy.linalg import qr

class QR_sparseT():  # it is much faster than to evaluate it directly 
        def __init__(self,A,B):
            self.A = A
            self.B = B
        def __mul__(self,x):
            return self.A*(self.B*x)

class QR_sparse:  # it is much faster than to evaluate it directly 
    """ primitive class with basic operations simulationg a sparse matrix"""
    def __init__(self,r, A,permut=None):

        
        
        self.r = sparse.csc_matrix(r)
        
        corr_ind = where(abs(self.r.diagonal())>1e-2)[0]

        self.r = self.r[corr_ind][:,corr_ind]
                                                                                                                                                                                                            
        if not permut is None:
            permut = permut[corr_ind ]

        ir = solve_triangular(self.r.todense(),eye(*self.r.shape),check_finite=False,overwrite_b=True)

        ind = argsort(permut)
        self.A = sparse.csc_matrix(A)[:,sort(permut)]

        self.ir = ir[ind,:]
        
        self.permut = permut
        self.int_permut = ind

        self.ir[abs(self.ir)<1e-10] = 0 #remove 'noise'
        self.ir = sparse.csc_matrix(self.ir)
        
        self.A.data = self.A.data.astype(single,copy=False)
        self.r.data = self.r.data.astype(single,copy=False)
        self.ir.data = self.ir.data.astype(single,copy=False)

        
        self.T = QR_sparseT(self.ir.T,self.A.T)
        self.shape = self.A.shape
        
        
 
        
        
    def __mul__(self,x):
        return self.A*(self.ir*x)

    def todense(self):
        return self.A*self.ir


def generalized_abel_transform(tokamak,tvec,N_r,N_a,name,cos_com=True,sin_com=True):

    
    extrapolate = 1.01
    M = get_rho_field_mat(tokamak, mean(tvec), extrapolate=extrapolate)

    R,Z = meshgrid(tokamak.xgrid, tokamak.ygrid)
    rhop,rc, zc = tokamak.mag_equilibrium(tvec,return_mean=True,surf_slice=slice(0,1))
    rc,zc = mean(rc),mean(zc)
    theta = arctan2(Z-zc,R-rc)
    theta = theta.ravel(order='F')[:,None]
    #ensure that in the middle point will be only N pixels
    M_min = sort(M.ravel())[N_a]  #else the problems with singularity of the transformation will rise. 
    M-= M_min
    M[M<0] = 0
    

    M[M==M.max()] = M[M!=M.max()].max()

    # interpolate on abel pixels
 
    M /= amax(M)
    c = 1.3
    M = tan(M*c)/tan(c)
    M**= 0.7


    M *= N_r

    M = M.ravel(order='F')


    M_0 = floor(M)
    M_1 = 1- M + M_0
    M_2 = M - M_0

    n_mag = int(amax(M_0))
    npix = tokamak.nx*tokamak.ny

    ##build a directly sparse transformation matrix 
    i_ind = []
    j1_ind,j2_ind = [],[]
    data1,data2 = [],[]
    
    
    for  i in range(n_mag):
        #function constant in angle but with triangular profile in radial direction
        ind = where(M_0 == i)[0]
        i_ind.append(ind)
        j1_ind.append(ones_like(ind)*max(0,i-1))
        j2_ind.append(ones_like(ind)*i)

        data1.append(M_1[ind])
        data2.append(M_2[ind])
        


    Transform = sparse.csc_matrix((hstack(data1), c_[hstack(i_ind),hstack(j1_ind)].T),shape=(npix, n_mag))
    Transform = Transform+sparse.csc_matrix((hstack(data2), c_[hstack(i_ind),hstack(j2_ind)].T),shape=(npix, n_mag))
    Transform_N = [Transform, ]
    
    for  n in range(N_a):
        if sin_com:Transform_N.append(sparse.spdiags(sin((n+1)*theta.T),(0,),npix, npix)*Transform)
        if cos_com:Transform_N.append(sparse.spdiags(cos((n+1)*theta.T),(0,),npix, npix)*Transform)
 

    Transform = sparse.hstack(Transform_N)
  
    try:
        from spqr.spqr import qr_sparse_r 
        r,e = qr_sparse_r(Transform)
    except Exception as e :
        debug( 'Warning sparseQR is not working! dense QR will be used')
        debug(str(e))
        if sparse.issparse(Transform):
            Transform = Transform.todense()
        r, = qr(Transform,overwrite_a=True,mode='r',check_finite=False)
        r = single(r[:r.shape[1],:])
        T = sparse.csc_matrix(single(Transform))
        e = arange(r.shape[1])
        
    T = QR_sparse(r,Transform,e)
    return T





def global_basis_transform(tokamak, tvec, N_r,N_a, basis,cos_com=True,sin_com=True):

    from zernike import zernike, bessel, polynom
    from annulus import get_bd_mat
    transform = eval(basis)


    extrapolate = 1.01

    nx = tokamak.nx
    ny = tokamak.ny

    rho = get_rho_field_mat(tokamak, tvec.mean(), extrapolate=extrapolate)
    rho /= amax(rho)
    
    R,Z = meshgrid(tokamak.xgrid, tokamak.ygrid)
    rhop,rc, zc = tokamak.mag_equilibrium(tvec,return_mean=True,surf_slice=slice(0,1))
    rc,zc = mean(rc),mean(zc)
    theta = arctan2(Z-zc,R-rc)
    
    BdMat = get_bd_mat(tokamak, time = tvec.mean())
    mask = isnan(rho) | array(rho == amax(rho)) | BdMat.reshape(ny,nx,order='F')


    Transform = []
    for ia in range(N_a+1):
        for isgn in [-1,1]:
            if isgn == -1 and ia == 0:
                continue
            if not cos_com and isgn == 1 and ia!= 0:
                continue
            if not sin_com and isgn == -1:
                continue
            
            for ir in range(N_r):

                if basis == 'zernike' and ((ir-ia*isgn)%2!=0 or abs(ia)>ir):
                    continue
                    
                z = asarray(transform(ir,ia*isgn,rho,theta,N_r))
                z[mask] = 0

                z = z.ravel(order='F')
                Transform.append(z)              


    Transform = vstack(Transform).T
   
     # orthogonalize base, slow.. 
    Transform_,R = qr(Transform,overwrite_a=True,mode='economic',check_finite=False,pivoting=False)  


    Transform_ = Transform_.astype('single')
    Transform_[:,abs(diag(R))>1e-6*R[0,0]]  #remove very colinear vectors
    Transform = sparse.csc_matrix(Transform_,copy=False)

    return Transform
 
