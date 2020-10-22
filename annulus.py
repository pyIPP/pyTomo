#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
from matplotlib.pyplot import *
import time
from  scipy.interpolate import griddata as sc_interpolate
from scipy.ndimage.morphology import binary_dilation, binary_closing,binary_opening
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
        bd_len = int(2*pi*sqrt(nx*nx))
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



 
    if not all(bd[0,:]==bd[-1,:]):   
        bd = r_[bd[(-1,),:],bd]

    random.seed(0)
    bd+= random.randn(*bd.shape)/1e6  #BUG problem with exactly vertical lines

    BdMat = zeros((ny,nx),dtype=int)
    dx = (xmax-xmin)/float(nx)
    dy = (ymax-ymin)/float(ny)
    ygridc = linspace(ymin,ymax,ny)#+dy/2 


    # algorithm based on Jordan curve theorem 
    bd_left  = minimum(bd[:-1,0],bd[1:,0])
    bd_right = maximum(bd[:-1,0],bd[1:,0])
    bd_down  = minimum(bd[:-1,1],bd[1:,1])
    bd_up    = maximum(bd[:-1,1],bd[1:,1])
    dby = diff(bd[:,1])
    dbx = diff(bd[:,0])
    
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

        n_cross = maximum(0,around((x_cross[ind]-xmin)/dx))
        n_cross = n_cross[n_cross<nx]
        #main idea of Jordan curve theorem 
        n_cross,count = unique(int_(n_cross), return_counts=True) #if a single pixel is crossed more times
        BdMat[iy,n_cross]+= count

    #main idea of Jordan curve theorem 
    cumsum(BdMat, out=BdMat,axis=1)
    
    BdMat = BdMat%2==0

    #from matplotlib.pylab import *
    #title('anulus')
    #imshow(BdMat, interpolation='nearest',aspect='equal',origin='lower',
           #extent=(xmin,xmax, ymin,ymax))
    #plot(bd[:,0], bd[:,1],'w')
    #xlim(xmin, xmax)
    #ylim(ymin, ymax)
    #show()
 

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
    xgridc = linspace(tokamak.xmin,tokamak.xmax,nx+1)
    xgridc = (xgridc[1:]+xgridc[:-1])/2
    ygridc = linspace(tokamak.ymin,tokamak.ymax,ny+1)
    ygridc = (ygridc[1:]+ygridc[:-1])/2

         

    if tokamak.use_pfm:
        return tokamak.get_psi(tvec,xgridc, ygridc)
        

    mag_rhop,fieldx, fieldy = tokamak.mag_equilibrium(tvec, return_mean = False)

    n_theta, n_fields,n_tvec = shape(fieldx)
    
        

    if extrapolate:
        #C = [0.99, 1, 1.01,1.02] if getboundary else extrapolate
        bdx,bdy = tokamak.get_boundary(n_theta, boundary,time=tvec.mean()).T  
        if extrapolate == infty:  #make sure that it will circumscribe whole grid
            min_dist = amin(hypot(bdx-mean(bdx),bdy-mean(bdy)))
            max_dist_x = max(abs(tokamak.xmax-mean(bdx)),abs(tokamak.xmin-mean(bdx)))
            max_dist_y = max(abs(tokamak.ymax-mean(bdy)),abs(tokamak.ymin-mean(bdy)))
            extrapolate = max(1,hypot(max_dist_x, max_dist_y)/min_dist)

         
        #used to extrapolated boundary for smoothing with boundary
        xch0 = single(outer(bdx-mean(bdx),extrapolate)+mean(bdx)) 
        ych0 = single(outer(bdy-mean(bdy),extrapolate)+mean(bdy))

        n_bor = len(xch0)

        fieldx = append(fieldx, tile(xch0[:,None], n_tvec),axis=1)
        fieldy = append(fieldy, tile(ych0[:,None], n_tvec),axis=1)
        mag_rhop = append(mag_rhop,extrapolate)
        n_fields+= size(extrapolate)

    num_mag = ones_like(fieldx,dtype=int)
    num_mag *= arange(n_fields)[:,None]
    mag_rhop = tile(mag_rhop,(n_theta,1))

    M = empty((n_tvec, ny, nx),dtype='single')
    Xgrid, Ygrid = meshgrid(xgridc, ygridc)
    NUM = hstack((mag_rhop[1:,1:].ravel(),mag_rhop[0,0]))
 
    for i in range(n_tvec):

        #remove central mag. surface because it is degenerated and 0 = 2*pi value
        X = hstack((fieldx[1:,1:,i].ravel(),fieldx[0,0,i]))
        Y = hstack((fieldy[1:,1:,i].ravel(),fieldy[0,0,i]))
   
        try:
            #it is working only for very carefuly prepared equilibrium, but it is much faster
            assert safe_interpolation, 'safe_interpolation'
            #this grid data is fast but often failures
            import delaunay  #originaly library from matplolib, but they have depricated it

            tri = delaunay.Triangulation(X,Y) #6.5ms 
            interp = tri.linear_interpolator(NUM)
            M[i] =  interp[ygridc[0]:ygridc[-1]:complex(0,ny),xgridc[0]:xgridc[-1]:complex(0,nx)]#5ms 
            
        except:
            M[i]= sc_interpolate(c_[X,Y],NUM,c_[Xgrid.ravel(),Ygrid.ravel()],fill_value=1.05).reshape(ny,nx)

    M[isnan(M)] = 1.0+extrapolate
 
    return squeeze( M)




