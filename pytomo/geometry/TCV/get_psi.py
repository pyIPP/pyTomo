#!/usr/bin/env python
# -*- coding: utf-8 -*-

#ssh todstrci@lac911.epfl.ch -L 8002:tcv1.epfl.ch:8000

from numpy import *
#import pmds
from matplotlib.pylab import *
import os.path 
import sys,os
sys.path.append(os.path.abspath('../../'))
import tqdm
from scipy.io import loadmat
from shared_modules import inside_convex_hull



class psi_dow_class:
    pass


def get_psi(connect, shotnumber,liuqenumber=None,return_contours=False):
    #function [psi_dow]=get_psi(shotnumber,liuqenumber,gridsize)
    #gridsize: interpolation to this gridsize (not neccessary)


    #psi_dow = psi_dow_class()
    
    #pmds.mdsopen('tcv_shot',shotnumber)
    connect.openTree('tcv_shot',shotnumber)

    if not liuqenumber is None:
        if liuqenumber == 1:
            numberstr=''
        elif liuqenumber == 2:
            numberstr='_2'
        elif liuqenumber == 3:
            numberstr='_3'
        else:
            raise Exception('wrong liuqenumber')
    else:
        numberstr=''

    #time = pmds.mdsvalue(r'dim_of(\results::psi'+numberstr+',2)')
    #lcfs = [pmds.mdsvalue(r'\results::r_contour'+numberstr),
           #pmds.mdsvalue(r'\results::z_contour'+numberstr)]
    #mag_axis = [pmds.mdsvalue(r'\results::r_axis'+numberstr),
               #pmds.mdsvalue(r'\results::z_axis'+numberstr)]
    #psi_axis = pmds.mdsvalue(r'\results::psi_axis'+numberstr)


    #xmesh = pmds.mdsvalue(r'dim_of(\results::psi'+numberstr+',0)')
    #ymesh = pmds.mdsvalue(r'dim_of(\results::psi'+numberstr+',1)')
    #Psi   = pmds.mdsvalue(r'\results::psi'+numberstr)
    ##psi_dow.q=pmds.mdsvalue(r'\results::q_psi'+numberstr)
    #pmds.mdsclose('tcv_shot',shotnumber)
    
    
    time = asarray(connect.get(r'dim_of(\results::psi'+numberstr+',2)'))
    lcfs = [asarray(connect.get(r'\results::r_contour'+numberstr)),
           asarray(connect.get(r'\results::z_contour'+numberstr))]
    mag_axis = [asarray(connect.get(r'\results::r_axis'+numberstr)),
               asarray(connect.get(r'\results::z_axis'+numberstr))]
    psi_axis = asarray(connect.get(r'\results::psi_axis'+numberstr))


    xmesh = asarray(connect.get(r'dim_of(\results::psi'+numberstr+',0)'))
    ymesh = asarray(connect.get(r'dim_of(\results::psi'+numberstr+',1)'))
    Psi   = asarray(connect.get(r'\results::psi'+numberstr))
    #psi_dow.q=connect.get(r'\results::q_psi'+numberstr)
    connect.closeTree('tcv_shot',shotnumber)
    
    
    
    
    
    

    if return_contours:
        rho,cnt = GetRhoContours(xmesh,ymesh, Psi,psi_axis,lcfs, 60)
        return time, mag_axis,lcfs,rho,cnt
    
    

    return time, mag_axis,lcfs,xmesh,ymesh, Psi


def GetRhoContours(Rmesh,Zmesh, Psi,psi_axis,lcfs, N):
    import matplotlib._cntr as cntr
    from scipy.ndimage.interpolation import map_coordinates
    from scipy.ndimage.interpolation import zoom

    contours = []
    #tin = np.atleast_1d(tin)
    
    nr = len(Rmesh)
    nz = len(Zmesh)
    
    zoom_fac = 4
    R_= linspace(Rmesh[0], Rmesh[-1], nr*zoom_fac)
    Z_= linspace(Zmesh[0], Zmesh[-1], nz*zoom_fac)
    R_,Z_ = np.meshgrid(R_,Z_)
    rho = linspace(0,1,N+1)
    psi_edge = Psi.max(1).max(1)

    

    #plot(psi_axis);plot(Psi.min(1).min(1));show()
                
    #import IPython
    #IPython.embed()
    for it in range(len(Psi)):  
        #it = 50

        psi_c = -(rho[1:]**2-1)*psi_axis[it]  #position of the contours
        #psi_c[psi_c < Psi[it].min()] = Psi[it].min() 
        #psi_c = -(rho[1:]**2-1)*Psi[it].min()    #position of the contours

        #contour searching routine from matplotlib
        #reduce discretization of the contours 
        Psi_zoom = zoom(Psi[it],zoom_fac,order=3)
        c = cntr.Cntr(R_,Z_, Psi_zoom)

        it_contours = []

        #plot(lcfs[0][it], lcfs[1][it]);show()
        ind = isfinite(lcfs[1][it])
        boundary = c_[lcfs[0][it][ind],lcfs[1][it][ind]]
        #inside_convex_hull(boundary[::-1] , lines[1])

        for f in psi_c:
            nlist = c.trace(level0=f,level1=f,nchunk=0)
            lines = nlist[:len(nlist)//2]
            if len(lines) == 0:
                lines = [np.empty((0,2)),]
            line = []
            #choose the longest line
            for l in lines:
                if len(l)>=len(line) and any(inside_convex_hull(boundary[::-1],l)):
                    line = l
            it_contours.append(line)  
            #for l in lines: plot(l[:,0], l[:,1])
            #xlim(Rmesh[0], Rmesh[-1])
            #ylim(Zmesh[0], Zmesh[-1])
            #show()
            #close()


        
 
        contours.append(it_contours)
        #for c in it_contours:  plot(c[:,0],c[:,1])
        #show()
    

    return rho, contours


def main():
    server =  'localhost:8002'
    pmds.mdsconnect( server)
    shotno = 40412
    psi_dow,cnt  = get_psi(shotno)
    
    import IPython
    IPython.embed()

    it = psi_dow.time.searchsorted(1.)
    from matplotlib.patches import Rectangle
    rec = Rectangle([psi_dow.vessel_x[0],psi_dow.vessel_y[0]],
              diff(psi_dow.vessel_x), diff(psi_dow.vessel_y),facecolor='none')
    contour(psi_dow.xmesh, psi_dow.ymesh, -psi_dow.data[it],20,lw=.5,colors='r')
    for c in cnt[it]: plot(c[:,0],c[:,1],lw=.5,c='k')
    plot(psi_dow.lcfsx[it],psi_dow.lcfsy[it],'b')
    plot( psi_dow.axisx[it], psi_dow.axisy[it],'x' )
    
    gca().add_patch(rec)
    axis('equal')
    show()
    
    
if __name__ == "__main__":
    main()

    