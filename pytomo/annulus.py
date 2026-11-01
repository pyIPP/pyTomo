#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import time
from  scipy.interpolate import griddata as sc_interpolate
from scipy.ndimage import binary_dilation, binary_closing,binary_opening
import traceback

def get_bd_mat(tokamak, nx=None, ny=None, boundary=None,time=None, bbox=None):
    """
    Find pixels of the main grid inside and outside of boundary. 
    :param array boundary:  array of coordinates point on boundary.
    :param class tokamak: ``Class`` of  tokamak with all important data.
    :param double time: time when the boundary should be calculated
    :param int nx: horizontal pixel resolutin of the grid 
    :param int ny: vertical pixel resolutin of the grid 

     This submodule is used to identify area with no emissivity
     Boundary do not have to be convex 

    """

    if nx is None or ny is None:
        nx = tokamak.nx
        ny = tokamak.ny


    
    if boundary is None:
        bd_len = int(2*np.pi*np.sqrt(nx*nx))
        bd = tokamak.get_boundary( bd_len,time=time)
    else:
        bd = boundary
        bd_len = boundary.shape[0]
        
    if not bbox is None:
        xmin,xmax,ymin,ymax = bbox
    else:
        xmin = tokamak.xmin
        xmax = tokamak.xmax
        ymin = tokamak.ymin
        ymax = tokamak.ymax



 
    if not np.all(bd[0,:]==bd[-1,:]):   
        bd = np.r_[bd[(-1,),:],bd]

    np.random.seed(0)
    bd+= np.random.randn(*bd.shape)/1e6  #BUG problem with exactly vertical lines

    BdMat = np.zeros((ny,nx),dtype=int)
    dx = (xmax-xmin)/float(nx)
    dy = (ymax-ymin)/float(ny)
    ygridc = np.linspace(ymin,ymax,ny)#+dy/2 


    # algorithm based on Jordan curve theorem 
    bd_left  = np.minimum(bd[:-1,0],bd[1:,0])
    bd_right = np.maximum(bd[:-1,0],bd[1:,0])
    bd_down  = np.minimum(bd[:-1,1],bd[1:,1])
    bd_up    = np.maximum(bd[:-1,1],bd[1:,1])
    dby = np.diff(bd[:,1])
    dbx = np.diff(bd[:,0])
    
    k =  dby/dbx
    k[abs(k)<1e-6] = 1e-6
    m = bd[:-1,1]-k*bd[:-1,0]
    for iy in range(ny):
        y0 = ygridc[iy]
        x_cross = (y0-m)/k
        
        #find  if this xcross is relevant
        ind  = x_cross >= bd_left
        ind &= x_cross <= bd_right
        ind &=   y0    >= bd_down
        ind &=   y0    <= bd_up

        n_cross = np.maximum(0,np.around((x_cross[ind]-xmin)/dx))
        n_cross = n_cross[n_cross<nx]
        #main idea of Jordan curve theorem 
        n_cross,count = np.unique(np.int_(n_cross), return_counts=True) #if a single pixel is crossed more times
        BdMat[iy,n_cross]+= count

    #main idea of Jordan curve theorem 
    np.cumsum(BdMat, out=BdMat,axis=1)
    
    BdMat = BdMat%2==0

    #import matplotlib.pylab as plt
    #plt.title('anulus')
    #plt.imshow(BdMat, interpolation='nearest',aspect='equal',origin='lower',
           #extent=(xmin,xmax, ymin,ymax))
    #plt.plot(bd[:,0], bd[:,1],'w')
    #plt.xlim(xmin, xmax)
    #plt.ylim(ymin, ymax)
    #plt.show()
 

    return BdMat.ravel(order='F') 



	
def get_rho_field_mat(tokamak, tvec,nx=None, ny=None,   extrapolate = 0,  boundary = None):
    """
    Load all annulus for given magnetic field and create matrix of size nx*ny

    :param class tokamak: ``Class`` of  tokamak with all important data.
    :param int nx, ny:  orizuntal and verticla grid resolution
    :param array magx, magy: Coordinates of magnetic field
    :param float extrapolate: extrapolate magnetic flux surface outside of the fluxsurfaces
    :param array boundary:  array of coordinates point on boundary.

    """

    if nx is None: nx = tokamak.nx
    if ny is None: ny = tokamak.ny

    #centers of the pixels
    xgridc = np.linspace(tokamak.xmin,tokamak.xmax,nx+1)
    xgridc = (xgridc[1:]+xgridc[:-1])/2
    ygridc = np.linspace(tokamak.ymin,tokamak.ymax,ny+1)
    ygridc = (ygridc[1:]+ygridc[:-1])/2

         
    try:
        if tokamak.use_pfm:
            return tokamak.get_psi(tvec,xgridc, ygridc)
    except:
        pass

    mag_rhop,fieldx, fieldy = tokamak.mag_equilibrium(tvec, return_mean = False)

    n_theta, n_fields,n_tvec = np.shape(fieldx)
    
        

    if extrapolate:
        #C = [0.99, 1, 1.01,1.02] if getboundary else extrapolate
        bdx,bdy = tokamak.get_boundary(n_theta, boundary,time=tvec.mean()).T  
        if extrapolate == np.inf:  #make sure that it will circumscribe whole grid
            min_dist = np.amin(np.hypot(bdx-np.mean(bdx),bdy-np.mean(bdy)))
            max_dist_x = max(abs(tokamak.xmax-np.mean(bdx)),abs(tokamak.xmin-np.mean(bdx)))
            max_dist_y = max(abs(tokamak.ymax-np.mean(bdy)),abs(tokamak.ymin-np.mean(bdy)))
            extrapolate = max(1,np.hypot(max_dist_x, max_dist_y)/min_dist)

         
        #used to extrapolated boundary for smoothing with boundary
        xch0 = np.single(np.outer(bdx-np.mean(bdx),extrapolate)+np.mean(bdx)) 
        ych0 = np.single(np.outer(bdy-np.mean(bdy),extrapolate)+np.mean(bdy))

        n_bor = len(xch0)

        fieldx = np.append(fieldx, np.tile(xch0[:,None], n_tvec),axis=1)
        fieldy = np.append(fieldy, np.tile(ych0[:,None], n_tvec),axis=1)
        mag_rhop = np.append(mag_rhop,extrapolate)
        n_fields+= np.size(extrapolate)

    num_mag = np.ones_like(fieldx,dtype=int)
    num_mag *= np.arange(n_fields)[:,None]
    mag_rhop = np.tile(mag_rhop,(n_theta,1))

    M = np.empty((n_tvec, ny, nx),dtype='single')
    Xgrid, Ygrid = np.meshgrid(xgridc, ygridc)
    NUM = np.hstack((mag_rhop[1:,1:].ravel(),mag_rhop[0,0]))
 
    for i in range(n_tvec):

        #remove central mag. surface because it is degenerated and 0 = 2*pi value
        X = np.hstack((fieldx[1:,1:,i].ravel(),fieldx[0,0,i]))
        Y = np.hstack((fieldy[1:,1:,i].ravel(),fieldy[0,0,i]))
   
        try:
            #it is working only for very carefuly prepared equilibrium, but it is much faster
            assert safe_interpolation, 'safe_interpolation'
            #this grid data is fast but often failures
            import delaunay  #originaly library from matplolib, but they have depricated it

            tri = delaunay.Triangulation(X,Y) #6.5ms 
            interp = tri.linear_interpolator(NUM)
            M[i] =  interp[ygridc[0]:ygridc[-1]:complex(0,ny),xgridc[0]:xgridc[-1]:complex(0,nx)]#5ms 
            
        except:
            M[i]= sc_interpolate(np.c_[X,Y],NUM,np.c_[Xgrid.ravel(),Ygrid.ravel()],fill_value=1.05).reshape(ny,nx)

    M[np.isnan(M)] = 1.0+extrapolate
 
    return np.squeeze( M)




