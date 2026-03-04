#!/usr/bin/env python
# -*- coding: utf-8 -*-

#cython:language_level=2
#cython:boundscheck=False
#cython:wraparound=False
#cython:infer_types=True
#cython:cdivision=True
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

from numpy import *
import numpy as np
#import matplotlib.cm as cm
#import matplotlib.mlab as mlab
from scipy import sparse

from libc.math cimport  sqrt, floor,hypot
cimport cython

cimport numpy as np
import time
#from numpy.matlib import repmat
#import time

from multiprocessing import  Pool, cpu_count

@cython.boundscheck(False)
@cython.wraparound(False)






def geom_mat_gen_cython(np.ndarray[double, ndim=2] xchord, np.ndarray[double, ndim=2] ychord,
                        np.ndarray[double] distance,Tok, int nx,int ny, np.ndarray[double] chord_profile):

    """
    Primitive but fast cython version of generator of the geometric matrix
    :param array distance: If the chords are curved, this is perpendicular distance of chord from center of tokamak
    :param class tokamak: ``Class`` of  tokamak with all important data.
    :param int nl: Number of all virtual detectors 
    :param int virt_chord: Number of virtual chords used in geometry matrix
    :param array xchords, xchords: Arrays of loaded geometry data

    """

    cdef double dx, dy
    cdef double d,t,c,L,r1,r2,z1,z2,tgalpha
    dx =(Tok.xmax-Tok.xmin)/double(nx)
    dy= (Tok.ymax-Tok.ymin)/double(ny)
    cdef int i,j,k, j_tmp

    cdef int nl = size(xchord,1)
    cdef np.ndarray[double] xgrid, ygrid,a,b

    xgrid= linspace(Tok.xmin, Tok.xmax, nx,endpoint=False)
    ygrid= linspace(Tok.ymin, Tok.ymax, ny,endpoint=False)

    a = (ychord[1,:]-ychord[0,:])/(xchord[1,:]-xchord[0,:]+1e-6) #slope
    b = ychord[1,:]-a*xchord[1,:] #intercept 

    cdef np.ndarray[double, ndim=2]  X, Y, A, B
    X = zeros((nx, nl))
    Y = zeros((nx, nl))
    A = zeros((nx, nl))
    B = zeros((nx, nl))
    
    


    for k in range(nl):
        L =xchord[0,k]*hypot(1,(ychord[1,k]-ychord[0,k])/(xchord[1,k]-xchord[0,k]))
        r1 = xchord[1,k]  #slit
        r2 = xchord[0,k]
        z1 = ychord[1,k]  #slit
        z2 = ychord[0,k]
        
        #keep the same direction
        if r1 < r2:
            for j in range(nx):
                X[j,k] = hypot(xgrid[j], distance[k]*(xgrid[j]-r1)/r1)
                Y[j,k] = a[k] *xgrid[j]+b[k]
        else:
            for j in range(nx):
                X[nx-j-1,k] = hypot(xgrid[j], distance[k]*(xgrid[j]-r1)/r1)
                Y[nx-j-1,k] = a[k] *xgrid[j]+b[k]
    
        #print k/50, hypot(X[0,k]-r1,Y[0,k]-z1),  hypot(X[0,k]-r2,Y[0,k]-z2),hypot(X[0,k]-r1,Y[0,k]-z1) > hypot(X[0,k]-r2,Y[0,k]-z2)
        #if hypot(X[0,k]-r1,Y[0,k]-z1) > hypot(X[0,k]-r2,Y[0,k]-z2):
            #X[:,k], Y[:,k] = X[::-1,k], Y[::-1,k]
    
    
                  

    direction = sign(xchord[1,:]-xchord[0,:])

    
    ncpu = cpu_count()
    pool = Pool(min(ncpu, 6))
    nvirt = len(chord_profile)
    nLOS = nl/nvirt

    II = [slice(i*nvirt,(i+1)*nvirt) for i in range(nLOS)]
    
    args = [(nx,ny,nvirt,a[ii],b[ii], xgrid, ygrid, dx,dy, distance[ii],
                   xchord[1,ii],direction[ii], Tok.xmin, Tok.xmax,chord_profile) for ii in II]
    

    TT = pool.map(help_fun, args)
    pool.close()
    pool.join()   


    TT = sparse.vstack(TT,format='csr')




    return TT,X,Y


def help_fun(input):
    nx,ny,nvirt,a,b, xgrid, ygrid, dx,dy, distance,xchord,direct,xmin,xmax,chord_profile = input
                   
    return generate(nx,ny,nvirt,a,b, xgrid, ygrid, dx,dy, distance,
                   xchord,direct,xmin, xmax,chord_profile)



def generate(int nx,int ny,int nvirt, np.ndarray[double] a,np.ndarray[double] b, \
    np.ndarray[double] xgrid, np.ndarray[double] ygrid,double dx,double dy,\
    np.ndarray[double] distance,np.ndarray[double] r1,np.ndarray[double] direct, double xmin,\
        double xmax,np.ndarray[double] chord_profile):
    
    cdef np.ndarray[float] TT_line = np.zeros(nx*ny,dtype=np.single)
    cdef np.ndarray[double] mod_xgrid ,mod_ygrid
    mod_xgrid = np.zeros(nx*ny)
    mod_ygrid = np.zeros(nx*ny)
    cdef double x1, x2, y1, y2, a_, b_, d_
    cdef float dist
    cdef int i, j

    with nogil:
        for i in range(nx*ny):
            mod_xgrid[i] =  xgrid[i//ny]
            mod_ygrid[i] =  ygrid[i%ny]

    for  i in range(nvirt):
        b_ = b[i]
        for  j in range(nx*ny):
            d_ = distance[i]*(mod_xgrid[j]-r1[i])/r1[i]
            
            #do not propagate LOS backwards behind the slit 
            if (direct[i] > 0.) and (mod_xgrid[j] > r1[i]+dx):
                continue
                
            if (direct[i] < 0.) and (mod_xgrid[j] < r1[i]-dx):
                continue
            
            

            y1 = mod_ygrid[j]
            y2 = y1+dy
            if mod_xgrid[j] >= d_:
                x1 = sqrt(mod_xgrid[j]**2-d_**2)      #nonlinearity of the chords
                x2 = sqrt((mod_xgrid[j]+dx)**2-d_**2)     # !!! pozor  to má vliv i na ten horní port na        geometrii !!!
                a_ = a[i]
            else:
                
                x1 =  sqrt(-mod_xgrid[j]**2+d_**2)      #nonlinearity of the chords
                x2 =  sqrt(-(mod_xgrid[j]+dx)**2+d_**2)     # !!! pozor  to má vliv i na ten horní port na  geometrii !!!
                a_ = -a[i]
                j =  abs( j//ny - 2*int( (d_-xmin)/dx ) ) * ny  + j%ny  # shift the x-axis
                if j < 0 or j > nx*ny-1:
                    continue

            l_1 = a_*x1+b_
            l_2 = a_*x2+b_
            l_3 = (y1-b_)/a_
            l_4 = (y2-b_)/a_
            dist = 0

            if (( l_1 < y2 and  l_2 > y1 ) or (  l_2 < y2 and  l_1 > y1 )):
                    if ( l_1 > y2 ):        # goes through the top side
                        if(l_2 > y1):       # goes through the right side
                            dist = hypot(y2-l_2, l_4-x2)
                        else  :                     # goes through the bottom side
                            dist = hypot(y2 -y1,l_4-l_3)
                    elif ( l_1 < y1) :      # goes through the bottom side
                        if(l_2 < y2) :      # goes through the right side
                            dist = hypot(y1-l_2, l_3-x2)
                        else :                       #goes through the top side
                            dist = hypot(y1-y2, l_3-l_4)
                    else :                  # goes through the center
                        if(l_2 > y2):       # goes through the top side
                            dist = hypot(l_1-y2, x1-l_4)
                        elif (l_2 < y1) :   # goes through the bottom side
                            dist = hypot(l_1-y1, x1-l_3)
                        else:
                            dist =  hypot(l_1-l_2, x1-x2 )

            if dist  > 0:
                TT_line[j] += dist*chord_profile[i]
            


    return sparse.csc_matrix(TT_line)


