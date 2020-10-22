#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
from matplotlib.pyplot import *
from mat_deriv_B import mat_deriv_B
from scipy.stats.mstats import mquantiles



def graph_derivation(time, tok, smoothing_matrix=None):

    xgrid = tok.xgrid   
    ygrid = tok.ygrid
    ny = tok.ny
    nx = tok.nx
    npix = nx*ny

    dx=tok.dx
    dy=tok.dy
    vessel_boundary = tok.get_boundary(N= 200,time=time)


    if smoothing_matrix is None:
        Bmat, diag_mat = mat_deriv_B(tok, time,3 ,4)
        B = Bmat[1] #paralel matrix 
    else:
        B = smoothing_matrix

    fields = []
    M = asarray(B.todense().T )      
    ind = where(M!= 0)
    
    
    for i, j in zip(ind[1],ind[0]): 
        x = xgrid[i//ny]+ dx/2
        y = ygrid[mod(i,ny)]+ dy/2
        dX = (xgrid[j//ny]+ dx/2 - x)*M[j,i]
        dY = (ygrid[mod(j,ny)]+ dy/2 - y)*M[j,i]
        fields.append(array([x,y,dX,dY]))
        
    field = vstack(fields).T


    field = field.T
    field[abs(field[:,3])+abs(field[:,2]) > 30, 2:4] = 0
    field2 = array([field[:,0], field[:,1],field[:,0]+field[:,2],field[:,3]+field[:,1]]).T
    field2 /= mean(field2)


    field[abs(field[:,2]) > 2*mquantiles(abs(field[:,2]),0.99),2:4] = 0
    field[abs(field[:,3]) > 2*mquantiles(abs(field[:,3]),0.99),2:4] = 0


    quiver(field[:,0], field[:,1],field[:,2],field[:,3], scale = 3*tok.norm/10.)  # 500
    rhop,magx, magy = tok.mag_equilibrium(time, return_mean = True )
    rho = linspace(0,2,100)
    magx, magy = tok.get_mag_contours(time,rho)
    for r,z in zip(magx[0], magy[0]):
        plot(r,z, 'r',lw=.2)
    plot(vessel_boundary[:,0], vessel_boundary[:,1], 'b')
    axis([tok.xmin,tok.xmax,tok.ymin,tok.ymax])
    gca().set_aspect('equal')
    show()


