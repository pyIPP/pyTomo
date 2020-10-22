#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
#from matplotlib.pyplot import *
from shared_modules import  debug,sigmoid, blur_image
from scipy.sparse import spdiags, eye
import time
from annulus import get_rho_field_mat, get_rho_field_mat, get_bd_mat
import config

class coord():
    pass



def mat_deriv_B(Tok, tvec, regularization,danis, nx=None, ny=None):
    """
    Prepare  matrix of derivation for unisotropic smoothing. Load magnetic field, prepare angles and determine directions of derivations. The result can be 1. or 2. derivation.

    :param class Tok: ``Class`` of  tokamak with all important data.
    :param array tvec: time vector of selected time range, matrix uses mean field during this time
    :param int regularization: Determine what type of derivation will be the result
    :var array magx,magy: Coordinater of magnetic field contours
    :var array atan2: Arctan of magnetic lines, "main result", slowest step
    :var spmatrix Bper, Bpar: Matrices of derivation
    :var danis - anisotropic ratio  - used here only by Imgesson smoothing matrix
    :param int nx,ny: horizontal, verticaůl resolution

    """


    if nx is None or ny is  None:
        nx = Tok.nx
        ny = Tok.ny

    coord = get_coord(Tok, nx, ny)
    
    
    
    t = time.time()
    


    npix = coord.npix
    
    diam = eye(npix,npix,    format='csc')     # diagonal
    diar = eye(npix,npix, ny,format='csc')# to reference right nearest neighbors
    dial = eye(npix,npix,-ny,format='csc')# reference left nearest neighbors
    diao = eye(npix,npix,  1,format='csc')# reference upper nearest neighbors
    diau = eye(npix,npix, -1,format='csc') # reference lower nearest neighbors

   
    diag_mat = diam, diar, dial, diao, diau
    
    Bmat = None
    

    if regularization in (-1,1,3,6):
        try:
    
            t = time.time()
            debug('This is MFR Anisotropic. Set smoothing along the flux surfaces...')

            rhop,magx, magy = Tok.mag_equilibrium(tvec, return_mean = True)
      

            magx = squeeze(magx)   #take mean field during tvec and omitted the center of field [:,0]
            magy = squeeze(magy)
  
            atan2 = calculate_theta_angle(Tok, mean(tvec), magx, magy, coord,extrapolate=1.2)
            rho = get_rho_field_mat(Tok, mean(tvec),extrapolate = 1.2)

            
            Bper1, Bpar1 = generate_aniso_matrix(coord, atan2,rho,  1,  1)
            Bper2, Bpar2 = generate_aniso_matrix(coord, atan2,rho, 1,  2)
            Bmat = [Bper1,Bpar1,Bper2,Bpar2]
          
            # correction of boundary conditions => zero derivation on borders !! + ensure positive definitness
            for i in range(4):
                Bmat[i] = Bmat[i] - spdiags(Bmat[i].sum(axis=1).T, 0, npix, npix)
                ind = array(abs(Bmat[i]).sum(axis=1).T == 0)
                corr = spdiags(ind,0,npix, npix)
                Bmat[i] = Bmat[i] + corr
                
      
        except IOError:
            raise
        except Exception as e:
            print('Anisotropic smoothing unavailable - Error: ' + str(e))
            regularization -= 1
            raise

    if regularization in (0,2,5):
        # find B for the rectangular mesh (normal smoothing)
        W = 1
        debug('This is MFR Classics (rectangular smoothing).')
        Bmat = rect_smoothing(Tok, coord, diam,diar, diao ,W)

    if regularization  == 4:
        debug('This is tikhonov Classics  (no smoothing)')
        Bmat = [eye(npix, npix)]
        
    if regularization  == 7:
        in_out_frac = 1 #BUG NOT FINISHED

        rhop,magx, magy = Tok.mag_equilibrium(tvec, return_mean = True,surf_slice=slice(-1,None),rho=1)

        rho = get_rho_field_mat(Tok, mean(tvec),extrapolate = infty)
        BdMat = get_bd_mat(Tok, boundary=squeeze(c_[magx, magy]))

        BdMat = blur_image(BdMat.reshape(Tok.ny, Tok.nx,order='F'),Tok.nx/50.)
        BdMat = BdMat.flatten('F')

        Bmat = generate_aniso_imgesson(coord, rho, BdMat, in_out_frac, danis)
        

    
    return Bmat, diag_mat



def rect_smoothing(Tok, coord, diam,diar, diao,W ):
    Bx =-diam+diar
    By =-diam+diao
    By /= coord.dx/coord.dy                      #to get the same size for different resolutions

    By = By*W
    Bx = Bx*W
    Bx = Bx - spdiags(Bx.sum(axis=1).T, 0, coord.npix, coord.npix)

    #correction of boundary conditions => zero derivation on borders !!
    return [Bx, By]






def get_coord(Tok, nx = None, ny = None):
    if nx is None: nx = Tok.nx
    if ny is None: ny = Tok.ny

    coord.nx = nx
    coord.ny = ny
    coord.npix =  nx*ny
    coord.dx = (Tok.xmax-Tok.xmin)/double(nx)
    coord.dy = (Tok.ymax-Tok.ymin)/double(ny)

    xgridc = linspace(Tok.xmin,Tok.xmax,nx+1)
    coord.xgrid = (xgridc[1:]+xgridc[:-1])/2
    ygridc = linspace(Tok.ymin,Tok.ymax,ny+1)
    coord.ygrid = (ygridc[1:]+ygridc[:-1])/2
    
    
    return coord



def calculate_theta_angle(Tok,time, magx, magy, coord, extrapolate = 10):
    """
    Intepolate/extrapolate data on x,y grid. The result is arctan in all points of the grid.
    Last magnetic line is used for extrapolation => smoothing with boundary.
    """

    
    bdx,bdy = Tok.get_boundary(size(magy,0),time=time).T

     #used to extrapolated boundary for smoothing with boundary
    magx_out = (bdx[::-1]-mean(bdx))*extrapolate+mean(bdx)
    magy_out = (bdy[::-1]-mean(bdy))*extrapolate+mean(bdy)
    
    magy = c_[magy, magy_out]
    magx = c_[magx, magx_out]

    xt   = zeros_like(magx)
    yt   = zeros_like(magx)
    cos1 = zeros_like(magx)
    sin1 = zeros_like(magx)

    dx =  diff(magx,axis=0)   #last point in boundary must be the same as the first
    dy =  diff(magy,axis=0)
    
     #positions at the place of the derivative: 
    xt[1:]= magx[:-1]+dx/2+coord.dx*.2 #BUG artificial correction of the observed subpixel error
    yt[1:]= magy[:-1]+dy/2-coord.dy*.2

    cos1[1:] = dx/(hypot(dx,dy)+1e-10)
    sin1[1:] = dy/(hypot(dx,dy)+1e-10)

    from scipy.interpolate import LinearNDInterpolator

    X = r_[xt[1:,1:].ravel(),xt[1,0]]
    Y = r_[yt[1:,1:].ravel(),yt[1,0]]
    Z = r_[cos1[1:,1:].ravel(),cos1[1,0]]

    try:
        assert safe_interpolation, 'safe_interpolation'


        import delaunay  #oroginaly library from matplolib, but they have depricated it
        tri = delaunay.Triangulation(X,Y)
        interp = tri.linear_interpolator(Z)
        cos_grid = interp[coord.ygrid[0]:coord.ygrid[-1]:complex(0,coord.ny),
                    coord.xgrid[0]:coord.xgrid[-1]:complex(0,coord.nx)]
    
        Z = r_[sin1[1:,1:].ravel(),sin1[1,0]]
        interp = tri.linear_interpolator(Z)
        sin_grid = interp[coord.ygrid[0]:coord.ygrid[-1]:complex(0,coord.ny),
                    coord.xgrid[0]:coord.xgrid[-1]:complex(0,coord.nx)]

    except:
        #slow but stable
        xi,yi = meshgrid(coord.xgrid,coord.ygrid)  #centers of pixels
        Linterp = LinearNDInterpolator(c_[X.ravel(),Y.ravel()], Z.ravel())
        cos_grid = Linterp(c_[xi.ravel(),yi.ravel()]).reshape(coord.ny,coord.nx)
        Z = r_[sin1[1:,1:].ravel(),sin1[1,0]]
        Linterp.values[:,0] = Z.ravel()#trick save a some  computing time        
        sin_grid = Linterp(c_[xi.ravel(),yi.ravel()]).reshape(coord.ny,coord.nx)  

    atan2 = array(arctan2(coord.dx*sin_grid,coord.dy*cos_grid))
    return atan2




def generate_aniso_matrix(coord, atan2,rho,weight, derivation):
    """
    Write prepared directions atan2 into the derivation matrix. Main magic 
    of this algorithm. Rewrite arctan into the directions 
    and decompose directions to parallel and oblique direction.

    :var array atan2: Arctan of magnetic lines, "main result", slowest step
    :param int derivation: Type of derivation that will be returned
        1. 2 are forward and backward derivatives
    
    """

    nx = coord.nx
    ny = coord.ny
    npix = coord.npix

    atan2 = atan2.ravel(order='F')
    rho = rho.ravel(order='F')/amax(rho)

    n_diag = 9
    Bper_tmp = zeros((n_diag,npix))
    Bpar_tmp = zeros((n_diag,npix))

    
    ind = isfinite(atan2)
    atan2 = atan2[ind]
    N_ind = len(atan2)



    
    for k in [0,1] :
        #decomposition of direction parallel and oblique
        direction =  int_(mod(floor(atan2 / (pi/4)+ k) , 8 ))
        dir = array((-1,2,3,4,1,-2,-3,-4),dtype=int)  #=> (-1,2,3,4,1,-2,-3,-4) is equal (-1, ny-1, ny, ny+1, 1, -ny+1, -ny, -ny-1) is compressed form
        next = squeeze(dir[direction])

        Arelativ = abs(pi/4 - mod(atan2+pi/4,pi/2))       #saw function, MAIN MAGIC

        if derivation in (1,2): #for second derivation use /2 instead of /sqrt(2)
            norm = sqrt(2)
        elif derivation == 3:
            norm  = 2
        else:
            print("bad regularization number, choose 2 or 4")
            exit()
            
            
        k1 = sin(2*Arelativ)/norm# obligue direction 
        k2 = cos(2*Arelativ)        # direction paralel with axes

        a = zeros(N_ind)

        ind_mod2 = bool_(mod(direction,2))
        a[ind_mod2] = k1[ind_mod2]
        a[~ind_mod2] = k2[~ind_mod2]

        if derivation == 1:
            Bper_tmp[n_diag//2+next,ind] = a
        if derivation == 2:
            Bper_tmp[n_diag//2-next,ind] = a             #allow second derivation
        if derivation == 3:
            Bper_tmp[n_diag//2+next,ind] = a
            Bper_tmp[n_diag//2-next,ind] = a             #allow second derivation



    ind = arange(8)
    pt1 = dir[ind%8]; pt2 = dir[(ind+2)%8]
    Bpar_tmp[n_diag//2+ pt2[ind]] = Bper_tmp[n_diag//2+ pt1[ind]]     #rotate directions by 90°

    #reduce paraell smoothing proportionaly to the circumference of the fluxsurface
    # by applying artificial weighting (for example for the divertor)

    Bpar_tmp/= sum((Bpar_tmp), 0)+1e-6
    Bper_tmp/= sum((Bper_tmp), 0)+1e-6

    min_smooth = .2
    Bper_tmp*= ((rho+min_smooth)*weight)
    Bpar_tmp*= ((rho+min_smooth)*weight)
    
    
    Bper_tmp[n_diag//2] = -Bper_tmp.sum(0) 
    Bpar_tmp[n_diag//2] = -Bpar_tmp.sum(0) 

    Bpar = spdiags(Bpar_tmp, (ny+1, ny, ny-1, 1,0,-1, -ny+1, -ny, -ny-1), npix, npix,format='csc').T
    Bper = spdiags(Bper_tmp, (ny+1, ny, ny-1, 1,0,-1, -ny+1, -ny, -ny-1), npix, npix,format='csc').T



    return Bper, Bpar





def generate_aniso_imgesson(coord, rho, BdMat, in_out_frac, danis):
    

    Psi = rho**2
    in_out_frac = 1
    Dparin  = sqrt(danis)*sigmoid( in_out_frac)*10
    Dparout = sqrt(danis)*sigmoid(-in_out_frac)*2
    Dperpin = 1/sqrt(danis)*sigmoid( in_out_frac)*10
    Dperpout= 1/sqrt(danis)*sigmoid(-in_out_frac)

    nx = coord.nx       # number of horisontal pixels
    ny = coord.ny       # number of vertical pixels
    npix  = coord.npix
    deltax= coord.dx
    deltay= coord.dy
    xmesh = coord.xgrid
    ymesh = coord.ygrid
    n_diag = 9  #result will by nona-diagonal matrix 

    #---finding the inidces of edges, corners, borders---

    # indices of edges
    i_top=slice(ny-1,npix,ny)                #or equivalently: i_top=[i_al,i_ba,i_ar]
    i_right=slice((nx-1)*ny,npix)        #or equivalently: i_right=[i_ur,i_br,i_ar]
    i_bottom=slice(0,(nx-1)*ny+1,ny)
    i_left=slice(0,ny)

    i_left_in=slice(ny, 2*ny)
    i_right_in=slice((nx-2)*ny,(nx-1)*ny)
    i_top_in=slice(ny-2,nx*ny-1,ny)                #or equivalently: i_top=[i_al,i_ba,i_ar]
    i_bottom_in=slice( 1,(nx-1)*ny+2,ny)

    #---finding the inidces of edges, corners, borders---
    

    #---diagonal and side diagonal matrices to construct differential operators---

    

    #---creating differential operators---
    D = eye(npix,format='csr')
    
    
    Dx = zeros((3,npix))
    Dx[0] =  -1/2.
    Dx[2] = 1/2.  #central derivative => /2 factor
    #calculating the gradients at the edges
    Dx[1,i_right] = -1
    Dx[0,i_right_in] = 1
    Dx[1,i_left]  = -1
    Dx[2,i_left_in]  = 1
    Dx/= deltax
    Dx = spdiags( Dx, (-ny,  0, ny), npix, npix,format='csr')



    #calculating the gradients inside the area (in this step the gradients at the edges are bad)
    Dy = zeros((3,npix))
    Dy[0] =  -1/2.
    Dy[2] = 1/2.  #central derivative => /2 factor
    #calculating the gradients at the edges
    Dy[1,i_top] = 1
    Dy[0,i_top_in] = -1
    Dy[1,i_bottom]  = -1
    Dy[2,i_bottom_in]  = 1
    Dy[0,i_top]  = 0
    Dy[2,i_bottom]  = 0
    Dy/= deltay
    Dy = spdiags( Dy, (-1,  0, 1), npix, npix,format='csr')
   
    Dxx = zeros((3,npix))
    Dxx[0] =  1
    Dxx[1] = -2
    Dxx[2] =  1  #central derivative => /2 factor
    #calculating the gradients at the edges
    Dxx[1,i_left] = -1
    Dxx[1,i_right] = -1
    Dxx/=  deltax**2
    Dxx = spdiags( Dxx, (-ny,0, ny), npix, npix,format='csr')


    ##inside
    Dyy = zeros((3,npix))
    Dyy[0] =  1
    Dyy[1] = -2
    Dyy[2] =  1  #central derivative => /2 factor
    #calculating the gradients at the edges
    Dyy[1,i_top] = -1
    Dyy[1,i_bottom] = -1
    Dyy[0,i_top] = 0
    Dyy[2,i_bottom] = 0
    Dyy/=  deltay**2
    Dyy = spdiags( Dyy, (-1,0, 1), npix, npix,format='csr')
   
    
    Dxy = Dy*Dx

    


    #---producing the diffusion matrix---
    xmeshgrid,_ = meshgrid(xmesh,ymesh)
    
    # we can define diffusion coefficient different for different parts of the plasma 
    # the original method caused to strong smoothing in the core 
    Dpar  = (Psi.flatten('F')**2+1e-1)*Dparin*(1-BdMat)  +Dparout *(BdMat)
    Dperp = (Psi.flatten('F')**2+1e-1)*Dperpin*(1-BdMat) +Dperpout*(BdMat)

    #---producing the diffusion matrix---



    #---creating c coeffitients---

  
    DyM, DxM = gradient(Psi, deltay,deltax)
    

    #NOTE it is faster then applying these derivation operators 
    DyyM = diff(Psi,2,0)/deltay**2
    DyyM = vstack((DyyM[0]/2, DyyM, DyyM[-1]))
    DxxM = diff(Psi,2,1)/deltax**2
    DxxM = hstack((DxxM[:,(0,)], DxxM, DxxM[:,(-1,)]))
    DxyM = gradient(DyM,deltay,deltax)[1]
    
    #WARNING due to linear interpolation can be the second order differences quite noisy!

    DxxM = DxxM.flatten(order='F')
    DyyM = DyyM.flatten(order='F')
    DxyM = DxyM.flatten(order='F')
    DyM  = DyM.flatten(order='F')
    DxM  = DxM.flatten(order='F')
    

    
    gradRho=DxM**2+DyM**2

    cxx=(Dperp*DxM**2+Dpar*DyM**2)/gradRho
    cyy=(Dperp*DyM**2+Dpar*DxM**2)/gradRho
    cxy=(Dperp-Dpar)*(DxM*DyM)/gradRho
    
    #NOTE it is faster then applying these derivation operators 
    dyDperp,dxDperp = gradient(Dperp.reshape(ny,nx,order='F'),deltay,deltax)
    dyDpar,  dxDpar = gradient( Dpar.reshape(ny,nx,order='F'),deltay,deltax)
    dyDperp = dyDperp.flatten(order='F')
    dxDperp = dxDperp.flatten(order='F')
    dyDpar  =  dyDpar.flatten(order='F')
    dxDpar  =  dxDpar.flatten(order='F')


    dDtermcx=DxM**2*dxDperp+DyM**2*dxDpar+DxM*DyM*(dyDperp-dyDpar)
    dDtermcy=DyM**2*dyDperp+DxM**2*dyDpar+DxM*DyM*(dxDperp-dxDpar)
    
    dlongtermcx=-2*((Dperp*DxM**2+Dpar*DyM**2)*(DxM*DxxM+DyM*DxyM)+
            (Dperp-Dpar)*(DxM*DyM)*(DxM*DxyM+DyM*DyyM))/gradRho

    dlongtermcy=-2*((Dperp*DyM**2+Dpar*DxM**2)*(DxM*DxyM+DyM*DyyM)+
            (Dperp-Dpar)*(DxM*DyM)*(DxM*DxxM+DyM*DxyM))/gradRho


    toroidaltermcx=cxx*gradRho/xmeshgrid.flatten(order='F')
    toroidaltermcy=cxy*gradRho/xmeshgrid.flatten(order='F')

    cx= ( 2*Dperp*DxxM*DxM+2*Dpar*DxyM*DyM+(Dperp-Dpar)*(DxyM*DyM+DyyM*DxM)+
            dDtermcx+dlongtermcx+toroidaltermcx)/gradRho
    
    cy= ( 2*Dperp*DyyM*DyM+2*Dpar*DxyM*DxM+(Dperp-Dpar)*(DxyM*DxM+DxxM*DyM)+
            dDtermcy+dlongtermcy+toroidaltermcy)/gradRho
    
    cx = spdiags(cx*deltax*deltay ,0,npix, npix)
    cy = spdiags(cy*deltax*deltay ,0,npix, npix)
    cxx= spdiags(cxx*deltax*deltay,0,npix, npix)
    cyy= spdiags(cyy*deltax*deltay,0,npix, npix)
    cxy= spdiags(cxy*deltax*deltay,0,npix, npix)
    
    B = cx*Dx+cy*Dy+cxx*Dxx+cyy*Dyy+2*cxy*Dxy

    return B, 



