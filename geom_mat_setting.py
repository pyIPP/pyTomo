#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
from scipy.interpolate import interpolate
from scipy import sparse
import sys, os
from numpy import *
from scipy.sparse import spdiags, eye
import gc
from annulus import get_bd_mat
from scipy.linalg import solve
from shared_modules import debug
import config
from scipy.signal import medfilt


    
def geom_mat_setting(tokamak,nx, ny, virt_chord, path=None):
    """
    Data settings for algorithm generating the geometric matrix. Try to load matrix from local cache, if it is not possible, try to generate new.


    :param class tok: ``Class`` of  tokamak with all important data.
    :param int virt_chord: Number of virtual chords used in geometry matrix
    :var array Xchords, Ychords: Arrays of chords coordinates, used in ``make_plots``
    :var spmatrix Tmat: Matrix of geometry

    """

    xmin = tokamak.xmin
    xmax = tokamak.xmax
    ymin = tokamak.ymin
    ymax = tokamak.ymax
    
    if path is None:
        path = tokamak.geometry_path
        
        
    ver = '' if tokamak.geometry_version is None else str(tokamak.geometry_version)
    name = path+'/Geom_matrix_%dx%d-%d-%1.2f-%1.2f-%1.2f-%.2f_%s_%s.npz'\
        %(nx,ny,virt_chord,xmin,xmax,ymin,ymin,str(tokamak.nl),ver)

    try:
  
        assert tokamak.name!='DIIID' ,  "Don't use cache at DIII-D"
        assert config.useCache ,  "Don't use cache"

        d = load(name, allow_pickle=True,encoding='latin1')
        Tmat = d['T'].item();Xchords = d['X'];Ychords = d['Y'];

    except Exception as e:
 
        if hasattr(tokamak, 'load_geom'):
            tokamak.load_geom()
   
        print('Generating geom. matrix')
        if tokamak.camera:
            xchords, ychords, zchords, nl = load_geom_3D(tokamak.geometry_path,list(tokamak.detectors_dict.keys()) )
            Tmat,Xchords,Ychords = generate_3D_matrix(xchords, ychords, zchords,tokamak)
           
        else:
            xchords, ychords, distance, nl,virt_chord  = loadgeometry(tokamak.geometry_path, list(tokamak.detectors_dict.keys()),   virt_chord)   

            #plot_projection_space(tokamak,xchords, ychords, virt_chord,nl)
       
            assert isnan(tokamak.nl) or nl == tokamak.nl, 'wrong number of LOS  '+str(nl)+' '+str(tokamak.nl)
         
            Tmat,Xchords,Ychords = generate_matrix(xchords, ychords, distance,virt_chord,nl,nx, ny, tokamak)
            tokamak.Xchords = Xchords
            tokamak.Ychords = Ychords
        
        #just for "ilustration" it calculated for the midlle of the tvec
        try:
            BdMat = get_bd_mat(tokamak, nx, ny, time=(tokamak.min_tvec+tokamak.max_tvec)/2)    # nx can be differnet of tokamk.nx !!!
        except:
            BdMat= None

        savez_compressed(name, T=Tmat,X=Xchords, Y=Ychords, xgrid = tokamak.xgrid,
                  ygrid = tokamak.ygrid, BdMat = BdMat)


    Tmat = sparse.csr_matrix(Tmat)
    #import IPython
    #IPython.embed()

    return Tmat,Xchords,Ychords


def prune_geo_matrix(Tmat,BdMat):
    #delete point ot of the boundary
    BdMat = BdMat.flatten('F')
    for i in range(Tmat.shape[0]): 
        data_slice = slice(Tmat.indptr[i],Tmat.indptr[i+1])
        Tmat.data[data_slice] *= ~BdMat[Tmat.indices[data_slice]]
    Tmat.eliminate_zeros()
    return Tmat




def loadgeometry(geometry_path,cameras, virt_chord):
    """
    Load geometry from ``geometry_path`` and prepare for processing.

    *Expected names of detectors are detector_[num]_[x,y].txt. Each file has 2-3
    columns. If there are 2, then the first is beginning of chord and the second
    is end, it is only X or Y coordinate. If there are 3 columns, the the first 
    two are range for virtual chords and the last is the other side (position of pinhole)*
    """


     #========================  Locations of detectors  ========================
    ychords = []
    xchords = []
    distance = []
    nl = 0
    import os.path

    for i,c in enumerate(cameras):
        xchord= loadtxt(geometry_path+'/detector_'+c+'_x.txt')
        ychord= loadtxt(geometry_path+'/detector_'+c+'_y.txt')
    
        ychord = atleast_2d(ychord).T
        xchord = atleast_2d(xchord).T
        
        ndet_los = size(xchord,1)
        if len(xchord) == 2 and len(ychord) == 2 and virt_chord > 1 :
            print("cannot use virtual chords !!")
            virt_chord = 1
        try:
            dist_tmp = loadtxt(geometry_path+'/dist_'+c+'.txt')
            if size(dist_tmp) == 1:
                dist_tmp = ones(ndet_los)*dist_tmp
                
        except:
            dist_tmp = zeros(ndet_los) ;  #distance of chords from the main axes, usually = 0
            
        nl += ndet_los
     
        xchord, ychord, dist  = virt_chords(xchord, ychord, virt_chord, dist_tmp)
        ychords.append(ychord), xchords.append(xchord), distance.append(dist)

    if nl == 0:
        raise NameError("Program couldn't import any geometric chords, check \
            input path and format of geometric data\n"+geometry_path+'/detector_%%s_x.txt')

    #=====================  Generate the geometric matrix ======================
    return hstack(xchords), hstack(ychords), hstack(distance), nl, virt_chord



def virt_chords(achord, bchord, virt_chord, dist):
    """
    Prepare virtual chords from loaded data if there are three columns in file.
    """
    n_chord = size(achord,1)
    n_senz = virt_chord*n_chord

    distance = repeat(dist,virt_chord)


    swap = False
    if len(achord) == 3:     #swap x,y if expansion is not in bchord
        swap = True
        achord ,bchord  = bchord , achord 
    elif len(bchord) == 2:
        bchord = vstack((bchord[0,:], bchord[1,:]))
        return achord, bchord, distance


    achords = tile(achord, (virt_chord,1)).reshape(2,-1,order='F')
 
    bchords =  zeros((2, n_senz))

    bchords[1,:] = tile(bchord[2,:], (virt_chord,1)).flatten(order='F')


    tmp_chord = zeros(n_senz)

    for i in range(n_chord):
        da = achord[1,i] - achord[0,i]
        db2 = bchord[1,i]- bchord[2,i]
        db1 = bchord[0,i]- bchord[2,i]

        alpha_1 = arctan2(bchord[0,i]-bchord[2,i], da)
        alpha_2 = arctan2(bchord[1,i]-bchord[2,i], da)
        
    
        dalpha = (alpha_1-alpha_2+pi)%(2*pi)-pi
  
        for j in range(virt_chord):
            bchords[0,i*virt_chord+j]=tan(alpha_2+dalpha/(virt_chord-1+1e-6)*(j+0.5e-6))*da +bchord[2,i]

    if swap:
        achords,bchords = bchords, achords  

    return achords, bchords, distance


def plot_projection_space(tokamak,xchords, ychords, virt_chord,nl):
    
    #plot detectors in the projection space, tested only for Asdex!!!
    import matplotlib.pylab as plt

    xchords = xchords.reshape(2,nl,virt_chord)
    ychords = ychords.reshape(2,nl,virt_chord)
    
    dets_dict = tokamak.detectors_dict if hasattr(tokamak, 'detectors_dict') else None

    from collections import OrderedDict
    cam_ind = OrderedDict()
    colors = ('b','k','r','g','y')

    for det_index, key in zip(tokamak.dets_index,list(dets_dict.keys())):
        if len(key)==2 and key[1].isdigit() : key =  key[0]  #convect different cameras from AUG together :(
        if tokamak.name=='DIIID' and len(key)==6 and key[1].isdigit() : key =  key[:-1]  #convect different cameras from AUG together :(
        
        if tokamak.name=='DIIID' and key in ['45R1','195R1'] : 
            key =  'SXR-TA'  #convect different cameras from AUG together :(

        if not key in  cam_ind: cam_ind[key] = []
        
        color = colors[len(cam_ind[key])]
        if key == 'SXR-TA':
            color = colors[len(cam_ind[key])+2]

        cam_ind[key].append((det_index,color))
        
    
    def Transformation(x,y,x0,y0):
 
        A,C = x
        B,D = y

        E,F = x0,y0  #center of the plasma 
        xi = arctan2( (D-B),(C-A))

        p = ((B-F)/(D-B)+(E-A)/(C-A))/hypot(1/(C-A), 1/(D-B))
        return -p*sign(xi+pi/2)*sign(xi-pi/2),pi-(xi+pi)%pi
    
    def TransformBoundary(bx,by,x0,y0,theta):
        
        proj_boundary = zeros((len(theta), 2))
        bx = bx - x0
        by = by - y0
        for i,t in enumerate(theta):
            bx_ = bx*cos(t)-by*sin(t)
            by_ = bx*sin(t)+by*cos(t)
            proj_boundary[i] = by_.min(), by_.max()
            
        return -proj_boundary
        
    time = tokamak.tvec.mean()
    
    rhop,X0,Y0 = tokamak.mag_equilibrium(tokamak.tvec.mean(), return_mean=True,surf_slice=[0,])
    X0,Y0 = X0[:,0].mean(),Y0[:,0].mean()
    
    theta = linspace(0,pi,100)
    
    
    xbnd, ybnd  = loadtxt(tokamak.geometry_path+'/border.txt',unpack=True)


    chamber_bnd = TransformBoundary(xbnd, ybnd,X0,Y0,theta)


    fig,ax = plt.subplots()
    rho_max = amax(abs(chamber_bnd))
    ax.fill_between(rad2deg(theta),chamber_bnd[:,0] , rho_max, interpolate=True, color='.9',edgecolor='k')
    ax.fill_between(rad2deg(theta),chamber_bnd[:,1], -rho_max, interpolate=True, color='.9',edgecolor='k')
    ax.set_ylim(-rho_max,rho_max)
    xbnd, ybnd = tokamak.get_boundary(1000, time=time).T
    from matplotlib.ticker import MaxNLocator

    chamber_bnd = TransformBoundary(xbnd, ybnd,X0,Y0,theta)
    ax.plot(rad2deg(theta),chamber_bnd ,'k--')
    proj_space_coord = []    
    for i,(det, ind_colors) in enumerate(cam_ind.items()):
        
        ncam = len(ind_colors)
        X = xchords.mean(-1)[:,ind_colors[ncam//2][0]]
        Y = ychords.mean(-1)[:,ind_colors[ncam//2][0]]
        

        ndets = X.shape[1]
        p,xi = Transformation(X[:,12*ndets//20],Y[:,12*ndets//20],X0,Y0)

        ax.text(rad2deg(xi.T),p.T,det,backgroundcolor='w',fontsize=13)

        
        for ind,c in ind_colors:
            c= 'k' if ncam == 1 else c

            p,xi = Transformation( xchords[:,ind],ychords[:,ind],X0,Y0)
            
            ax.plot(rad2deg(xi.T),p.T,','+c)
            ax.plot(rad2deg(xi[:,virt_chord//2]),p[:,virt_chord//2],'o'+c)
           
            proj_space_coord.append([rad2deg(xi[:,virt_chord//2]),p[:,virt_chord//2]])
    
            
    ax.set_xlabel('$\\xi $ [deg]')
    ax.set_ylabel('p [m]')
    ax.set_ylim(-rho_max,rho_max)
    ax.set_xlim(0,180)
    ax.axhline(0,ls='--',c='k')
    ax.text(1.01/pi*180,0.02,'Plasma core')
    ax.text(.39/pi*180,-.73,'LCFS')
    ax.text(.69/pi*180,-0.90,'Borders')
    ax.xaxis.set_major_locator(MaxNLocator(6))

    plt.show()


    from phantom_generator import phantom_generator
    print('ready')

    _,_,_,G0 = phantom_generator(tokamak, array(4), profile = 'Gaussian', hollowness = 0.98,edge = 0.7)
    G0 =  G0.mean(1).reshape(tokamak.ny, tokamak.nx,order='F')
    rvec, zvec = tokamak.xgrid+tokamak.dx/2, tokamak.ygrid+tokamak.dy/2
    try:
        #load a specific profile
        name = 'Emissivity_4.513-4.514_rec_30812.npz'
        data = load('./tmp/'+name)
        G0 = (single(data['gres'])*data['gres_norm']).mean(-1)
        rvec, zvec = data['rvec'], data['zvec']
    except:
        pass

    theta = linspace(0,pi,200)
    p = linspace(ax.get_ylim()[0], ax.get_ylim()[1],200)
    sinogram = zeros((len(theta), len(p)),dtype=single)

    
    
    scaling = array([rvec[1]-rvec[0],zvec[1]-zvec[0]])
    offset = array([rvec[0],zvec[0]])
    
    from scipy.ndimage.interpolation import map_coordinates
    #WARNING it is not exact flax surface averadge!!
    X,Y = meshgrid(p,p)


    
    from scipy.integrate import trapz

    for i,t in enumerate(theta):
        X_ = X*cos(t)-Y*sin(t)+X0
        Y_ = X*sin(t)+Y*cos(t)+Y0

        coords = c_[X_.ravel(),Y_.ravel()].T
        idx = (coords-offset[:,None])/ scaling[:,None]
        map_prof = map_coordinates(G0.T,idx,order=2)
        map_prof = map_prof.reshape(X.shape)
        
        sinogram[i] =  trapz(map_prof, p, axis=1)
    

    sinogram/= amax(sinogram)
    
    ax.contour(theta/pi*180, p, sinogram[::-1].T,linspace(0.01,1,10),vmin=0,vmax=1,colors='0.5', linewidths = 0.5)

    plt.figure()

    proj_space = hstack(proj_space_coord)
    
    
    data = load('./tmp/data_%d.npz'%tokamak.shot)
    dets = data['dets']
    err = data['err'][:,dets].mean(0)
    data = data['data'][:,dets].mean(0)
    ind = isfinite(err)
    proj_space = proj_space[:,dets][:,ind]
    data = data[ind]

    
    from scipy.interpolate import griddata
   
    grid_x, grid_y = mgrid[-1.2:1.2:1000j, 0:180:1000j]
    
    figure()

    pointsx = hstack((proj_space[0],proj_space[0]-180,180+proj_space[0]))
    pointsy = hstack((proj_space[1],-proj_space[1],-proj_space[1]))
    values = hstack((data,)*3)
    D = griddata( (pointsx,pointsy), values,  (grid_y, grid_x), method='nearest',rescale=True)
    plt.pcolormesh(grid_y.T, grid_x.T,sqrt(D).T,vmin=0,vmax=sqrt(amax(data)),cmap='nipy_spectral')
    plt.plot(proj_space[0],proj_space[1],'wo',markeredgecolor='k')
    plt.ylim(-.7, .7)
    plt.xlabel('$\\xi$ [deg]')
    plt.ylabel('p [m]')
    plt.xlim(0,180)
    c = colorbar()
    c.set_label('Intensity [W/m$^2$]')
    
    
    I,J = 6,-2 #IL
    I,J = 3,-2 #HL
    I,J = 9,-1 #MJ
    I,J = 1,-3 #GK
    I,J = -4,-1 #K1M
    I,J = -3,-1 #K2M
    I,J = 6,-1 #IM

    

    
    pairs = {'IL':(6,-2),\
    'HL':(3,-2),\
    'MJ':(9,-1),\
    'GK':(1,-3),\
    'K1M':(-4,-1),\
    'K2M':(-3,-1),\
    'IM':(6,-1)}

    for name,(I,J) in list(pairs.items()):
        stop=False
        for i in range(len(proj_space_coord[I][0])-1):
            if stop:break
            for j in range(len(proj_space_coord[J][0])-1):
                Q0 = c_[proj_space_coord[I]][:,i]
                P0 = c_[proj_space_coord[J]][:,j]
                Q1 = diff(c_[proj_space_coord[I]][:,i:i+2],axis=1)
                P1 = diff(c_[proj_space_coord[J]][:,j:j+2],axis=1)
                t = linalg.solve(c_[P1,-Q1],Q0-P0)
                #print i,j,t
                if all((t>0)&(t<1)) and abs(Q1[0])< 10 and abs(P1[0])< 10 :
                    i = i+t[1]+sum([len(a[0]) for a in proj_space_coord[:I]])
                    j = j+t[0]+sum([len(a[0]) for a in proj_space_coord[:J]])
                    stop=True
                    pairs[name] = (i,j)
                    break

    
    data = load('./tmp/data_%d.npz'%tokamak.shot)
    tvec = data['tvec']
    dets_dict = data['dets_dict']
    data = data['data']
    
    

    from scipy.interpolate import interp1d

    figure()
    proj_space = hstack(proj_space_coord)
    ticks = []
    c= 'r','g','b','k','y','m','0.5'
    for i,(name,(I,J)) in enumerate(pairs.items()):
        x = arange(proj_space.shape[1])
        if floor(I) in dets and ceil(I) in dets and floor(J) in dets and ceil(J) in dets:
            y1 = exp(interp1d(x, log(data),axis=1,kind='linear',assume_sorted=True)(I))
            y2 = exp(interp1d(x, log(data),axis=1,kind='linear',assume_sorted=True)(J))
            errorbar(len(ticks),nanmean(y2/y1), std(y2/y1),c='r',fmt='x')
            ticks.append(name)
    xticks( arange(len(ticks)),ticks, rotation='horizontal')
    xlim(-.5,len(ticks)-.5)
    ylim(.8,1.2)
    axhline(1)

    
    import IPython
    IPython.embed()
  
    
    


def generate_matrix(xchords, ychords, distance,virt_chord,nl,nx, ny,  tokamak):

    """
    From loaded geometry and tokamak parameters prepare data for matrix geenrators. Try to load cython, if it is not possible, try to use slower python version.

    :param array distance: If the chords are nonlinear, this is perpendicular distance if detector (chords) from center of tokamak
    :param class tokamak: ``Class`` of  tokamak with all important data.
    :param int nl: Number of detectors
    :param int virt_chord: Number of virtual chords used in geometry matrix
    :param array xchords, xchords: Arrays of loaded geometry data
    :var spmatrix Tmat: Matrix of geometry

    :var array Xchords, Ychords: Arrays of chords coordinates, used in ``make_plots``
    :var spmatrix Tmat: Matrix of geometry
    :var array fk: Profile of virtual chords
    """


    geometry_path = tokamak.geometry_path
    
    #===================== use custom profile of chords =======================

    fk = ones(virt_chord)/virt_chord
        
    if virt_chord > 1:
        if os.path.isfile(geometry_path+'/virt_chord_profile.txt'):
            fk =  loadtxt(geometry_path+'/virt_chord_profile.txt').T
            fk  = interp(linspace(0,1,virt_chord), linspace(0,1,len(fk)),fk)
            fk = fk/sum(fk)
        else:
            debug( 'virt_chord_profile not availible')
        
    try:
        Tmat,Xchords,Ychords = gen_cython(xchords, ychords, distance,nx, ny, tokamak,fk)
    except:
        raise
        if any(distance != 0) or alpha != 0: #distance of chords from the main axes, usually = 0
            print("Check your Cython/Pyrex/Numpy instalation!!")
            raise NameError("Nonlinear chord are not supported in non-Cython version")
        
        else:
        
            from geom_mat_gen import geom_mat_gen
            Tmat = geom_mat_gen(xchords,ychords,tokamak, nx, ny).T
            Tmat = sparse.csc_matrix(Tmat)

    return Tmat ,Xchords,Ychords


def gen_cython(xchords, ychords, distance,nx, ny, Tok,chord_profile):
    """
    Load and compile cython version of generator for matrix of derivation. Only this simple cython version supports nonlinear chords from JET

    *Cython version, stupid algorithm but 5x faster than better matrix algorithm*

    """

    import pyximport

    os.environ['CPATH'] = get_include()
    if sys.platform == 'win32':
        mingw_setup_args = { 'options': { 'build_ext': { 'compiler': 'mingw32' }} }
        pyximport.install(setup_args=mingw_setup_args)
    else:
        pyximport.install()

    # troubles with JET old python !!!!!
    from geom_mat_gen_cython import geom_mat_gen_cython
    #import IPython
    #IPython.embed()

    
    Tmat,Xchords,Ychords = geom_mat_gen_cython(xchords,ychords,distance, Tok,  nx, ny,chord_profile)

    
    virt_chord = len(chord_profile)
    Xchords = Xchords[:,virt_chord//2::virt_chord]
    Ychords = Ychords[:,virt_chord//2::virt_chord]
   
    return Tmat,Xchords,Ychords




































def generate_3D_matrix(xchord, ychord, zchord, tokamak):
    """
    From loaded geometry and tokamak parameters prepare data for matrix geenrators. Try to load cython, if it is not possible, try to use slower python version.

    :param class tokamak: ``Class`` of  tokamak with all important data.
    :param int nl: Number of detectors
    :param array xchords, xchords: Arrays of loaded geometry data
    :var spmatrix Tmat: Matrix of geometry

    :var array Xchords, Ychords: Arrays of chords coordinates, used in ``make_plots``
    :var spmatrix Tmat: Matrix of geometry
    :var array fk: Profile of virtual chords
    """

    print("mage 3D geom")

    #nx = 100  # rezliseni te interpolace 

    from geom_mat_gen_3D import geom_mat_gen_3D, geom_mat_gen


    cam_zoom = tokamak.cam_zoom

    xmin = tokamak.xmin
    xmax = tokamak.xmax
    ymin = tokamak.ymin
    ymax = tokamak.ymax
    zmin = tokamak.ymin
    zmax = tokamak.ymax
    coord = [-xmax, xmax, ymin, ymax, -xmax, xmax]

    # a) pole X a Y musí obsahovat vsechny hodnoty od xmin do xmax a ymin do ymax

    nl = size(xchord,1)
   
    print("generuji chordy")
    ## ============== chordy ==============================
    a = zeros(nl)
    b = zeros(nl)
    c = zeros(nl)
    d = zeros(nl)
    for k in range(nl):
        a[k] = (xchord[1,k]-xchord[0,k]+1e-8)/(zchord[1,k]-zchord[0,k]+1e-8)
        b[k] = (ychord[1,k]-ychord[0,k]+1e-8)/(zchord[1,k]-zchord[0,k]+1e-8)
        c[k] = xchord[0,k] - a[k]*zchord[0,k]
        d[k] = ychord[0,k] - b[k]*zchord[0,k]


    print("matice chord")
    Nsteps = 500
    
    X = zeros((Nsteps, nl), dtype=single)
    Y = zeros((Nsteps, nl), dtype=single)
    Z = zeros((Nsteps, nl), dtype=single)

    # perform interpolation to a better coordinatet
    for k in range(nl):
        #X[:,k] = linspace(xmax, -xmax, nx)
        #Y[:,k] = a[k]*X[:,k] + c[k]
        #Z[:,k] = b[k]*X[:,k] + d[k]
        if tokamak.cam_coord[2] < 0:
            Z[:,k] = linspace(-xmax, xmax, Nsteps)
        else:
            Z[:,k] = linspace(xmax, -xmax, Nsteps)
        X[:,k] = a[k]*Z[:,k] + c[k]
        Y[:,k] = b[k]*Z[:,k] + d[k]

    del a,b,c,d
    
    R = sqrt(X**2 + Z**2)
    

    ## ==============pixely ==========================================
    zoom = tokamak.pixel_zoom
    tokamak.nl /= zoom**2
    npix = tokamak.nx*tokamak.ny*zoom**2
    tok_ny = double(tokamak.ny*zoom)
    tok_nx = double(tokamak.nx*zoom)
    xgrid = linspace(xmin, xmax, tok_nx)
    ygrid = linspace(ymin, ymax, tok_ny)
    dx = xgrid[1] - xgrid[0]
    dy = ygrid[1] - ygrid[0]
    
    bd = tokamak.get_boundary()

    print("BD MAT ")
    
    nx_tmp =  max(200, tok_nx)
    ny_tmp =  max(200, tok_ny)

    from annulus import get_bd_mat
    BdMat_tmp = get_bd_mat(tokamak, nx_tmp , ny_tmp,time=(tokamak.tmin+tokamak.tmax)/2) 
    BdMat_tmp = reshape(BdMat_tmp, (ny_tmp, nx_tmp), order='F')
    

    def get_inside(X,Y,Z,k):
        #r = sqrt(X[:,k]**2 + Z[:,k]**2)
        r = hypot(X[:,k],Z[:,k] )
        y = Y[:,k]
        i = int_((r - xmin)/(xmax-xmin)*nx_tmp)
        j = int_((y - ymin)/(ymax-ymin)*ny_tmp)
        #print i,j
        i[ (i < 0) | (i >= nx_tmp) ] = 0
        j[ (j < 0) | (j >= ny_tmp) ] = 0

 
        inside = ~BdMat_tmp[j,i]
          
        inside = medfilt(inside, 5)  # remove some problems with multiple bound. crossing
           
        inside = where(inside)[0]
        if len(inside) == 1:  # avoid problems with 1px crossbourders
            inside = concatenate([inside, [inside[-1] + 1]])
                
        wrong = (i == 0) | (j == 0)


        inside = in1d( arange(len(X)), inside) & ~wrong

        return where(inside)[0]
                
    def clean_chords(X,Y,Z):
        """ detect chords in chamber and remove parts that shouldn't be visible
        """
        nl = size(X,1)
        Xch = ones((2,nl))*nan
        Ych = ones((2,nl))*nan
        Zch = ones((2,nl))*nan
        
        for k in range(nl):
            inside = get_inside(X,Y,Z,k)   
            
            #plot(bd[:,0], bd[:,1])
            #plot(sqrt(X[:,k]**2+Z[:,k]**2), Y[:,k], ':')
            #plot(sqrt(X[inside,k]**2+Z[inside,k]**2), Y[inside,k], 'k')
                
            end = where(abs(diff( inside )) > 1)[0]
            end  = setdiff1d( end, [0])

            if len(end) > 0:
                inside = inside[:end[0]]
            
            if len(inside) > 1:
                Xch[:,k] =  [X[inside[0],k], X[inside[-1], k]]
                Ych[:,k] =  [Y[inside[0],k], Y[inside[-1], k]]
                Zch[:,k] =  [Z[inside[0],k], Z[inside[-1], k]]
        
        return Xch, Ych, Zch
    
    print("clean chords")
    
    Xch, Ych, Zch = clean_chords(X,Y,Z)

    from scipy.interpolate import interp1d

    N = 100

    #Xch =  interp1d([0,1], Xch.T)(linspace(0,1,N)).T
    #Ych =  interp1d([0,1], Ych.T)(linspace(0,1,N)).T
    #Zch =  interp1d([0,1], Zch.T)(linspace(0,1,N)).T
    
    ####!!  cut off geometry 
    #ax = tokamak.geometry_axis
    #ind = (Xch > ax[0]) & (Xch < ax[1]) & (Ych > ax[2]) & (Ych < ax[3]) & (Zch > ax[4]) & (Zch < ax[5])
    #Xch[~ind] = nan
    #Ych[~ind] = nan
    #Zch[~ind] = nan


    #diam=eye(npix,npix).tocsr()                                     # diagonal
    #diar=spdiags((ones((1,npix))), tok_ny,npix,npix,format='csc')             # to reference right nearest neighbors
    #dial=spdiags((ones((1,npix))),-tok_ny,npix,npix,format='csc')             # reference left nearest neighbors
    #diao=spdiags((ones((1,npix))), 1,     npix,npix,format='csc')              # reference upper nearest neighbors
    #diau=spdiags((ones((1,npix))),-1,npix,npix,format='csc')              # reference lower nearest neighbors
    
   
    diam=eye(npix,npix,format='csc')          # diagonal
    diar=eye(npix,npix, tok_ny,format='csc')# to reference right nearest neighbors
    dial=eye(npix,npix, -tok_ny,format='csc')# reference left nearest neighbors
    diao=eye(npix,npix, 1,format='csc')# reference upper nearest neighbors
    diau=eye(npix,npix, -1,format='csc') # reference lower nearest neighbors

    diag_mat = array([diam, diar, dial, diao, diau])/5

    bd_vec = (get_bd_mat(tokamak,time=(tokamak.tmin+tokamak.tmax)/2) != 0) * sum(diag_mat,0)**2
    bd_vec = reshape(bd_vec, (tok_ny, tok_nx), order='F')
    
    bd = tokamak.get_boundary()
    #imshow(bd_vec,  interpolation="nearest", extent=tokamak.area_axis); colorbar()
    #plot(bd[:,0], -bd[:,1])
    #show()
    
    #exit()
    
    boundary = (bd_vec <1) & (bd_vec >0)
    #boundary = (bd_vec <1) & (bd_vec >0.5)
    boundary[[0,-1], :] = 0
    boundary[:,[0,-1]] = 0

    #pcolor(boundary)
    ##show()
    #imshow(boundary,  interpolation="nearest", extent=tokamak.area_axis); colorbar()
    #plot(bd[:,0], -bd[:,1])
    #show()
    
    #exit()

    #print "done"

    def get_pix(BdMat):
        #  toroidal coordinates 
        N = 50


        bd_vec = (1-BdMat)  * sum(diag_mat,0)
        bd_vec = reshape(bd_vec, (tok_ny, tok_nx), order='F')

        Xpix = zeros((N, npix))
        Ypix = zeros((N, npix))
        Zpix = zeros((N, npix))

        print("generuji pixely")

        for k in range(npix):
            # vykresli jen okraj !!!!!!!!!!!!!!!!!!!!!
            #if not bd_vec[(k%tok_ny)/zoom, (k/tok_ny)/zoom] == 1:
                #continue

            r0 = xgrid[k//tok_ny]
            y0 = ygrid[k%tok_ny]
            center = (xmax+xmin)/2
            phi = linspace(0,2*pi, N)  #linspace(pi/2,pi+pi/2, N)+pi
            Q = 5
            alpha =   linspace(0,2*pi/Q, N)

            #Xpix[:,k] = r0*cos(phi)
            #Ypix[:,k] = y0
            #Zpix[:,k] = r0*sin(phi)

            R = sqrt( (r0 - center)**2 + y0**2 )
            alpha_0 = arctan2((r0 - center), y0)
            r = R * sin(alpha+alpha_0)+center 
            y = R * cos(alpha+alpha_0) 
                    
            Xpix[:,k] = r * cos(phi)
            Ypix[:,k] =  y 
            Zpix[:,k] = r * sin(phi)

        #show()
        
        return Xpix, Ypix, Zpix


    BdMat =  get_bd_mat(tokamak,time=(tokamak.tmin+tokamak.tmax)/2)
    BdMat = reshape(BdMat, (-1), order='F')

    ############################################################x
    Xpix, Ypix, Zpix = get_pix(BdMat)
    #############################################################x
    #ind  = all(Xpix != 0, 0)
    #Xpix = Xpix[:,ind]
    #Ypix = Ypix[:,ind]
    #Zpix = Zpix[:,ind]
    #plot(Ypix, Zpix)
    #show()
    #print shape(BdMat)
    #pcolor(xgrid, ygrid, reshape(BdMat, (tok_nx, tok_ny), order="F"))
    #print Xpix
    #print Zpix
    #print Ypix
    #plot(sqrt(Xpix**2 + Zpix**2), Ypix, '.')
    #show()
    
    

    ##print "plotting"
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = figure()
    #ax = fig.add_subplot(111, projections='3d')
    #for i in range(size(Xpix,1)):
        #print shape(Xpix)
        #print Xpix[:,i]
        #ax.plot( Xpix[:,i],Ypix[:,i],Zpix[:,i], 'k-', linewidth=0.5)

    #ax.set_xlim3d(-xmax, xmax)
    #ax.set_ylim3d(ymin, ymax)
    #ax.set_zlim3d(-ymax, ymax)
    #ax.set_xlabel('x')
    #ax.set_ylabel('y')
    #ax.set_zlabel('z')


    #show()

    
    #print Xpix
    #exit()
    
    

    Rch = sqrt(Xch**2 + Zch**2)

    if not tokamak.allow_reflections :
        Tmat = geom_mat_gen(Xch, Ych, Zch, tokamak,tokamak.nx*zoom,tokamak.ny*zoom, nl, tokamak.geometry_axis )
        Tmat = sparse.csc_matrix(Tmat)
        return Tmat, Rch, Ych

    
    

    cam_zoom = tokamak.cam_zoom
     
    # ============ compass ============================
     
    D = sqrt((Xch[-1,:]-Xch[0,:])**2 +(Ych[-1,:]-Ych[0,:])**2 + (Zch[-1,:]-Zch[0,:])**2)
    


    def norm(a):
        return sqrt(sum(a**2))

    
    
    #BUG obsolete!! 
    boundary_simple = tokamak.boundary_simple
       
    
    bd = tokamak.get_boundary(500, boundary_simple )
    BdMat = get_bd_mat(tokamak, boundary = boundary_simple) * sum(diag_mat,0) != 0
    

    
    print("diff model")

    #chords_all = list()
    #chords_tmp = zeros((2, nl, 3))
    #chords_tmp[:,:,0] = Xch
    #chords_tmp[:,:,1] = Ych
    #chords_tmp[:,:,2] = Zch

    def get_diff_model(r,y):
        d = (r - bd[:,0])**2 + (y - bd[:,1])**2 
        s = sort(argsort(d)[:2])

        [r,y] = bd[s[0],:]  # smoother position  
        
        Nvec  = array( [bd[s[1],0] -bd[s[0],0] , bd[s[1],1]-bd[s[0],1] ])
        Nvec = copy(Nvec)
        Nvec[0] *= -1
        Nvec[[0,1]] = Nvec[[1,0]]
        Nvec /= norm(Nvec)

        x = ((r - xmin)/(xmax-xmin)*tok_nx)
        y = ((y - ymin)/(ymax-ymin)*tok_ny)
        pos_x = (X - x)/tok_nx
        pos_y = (Y - y)/tok_ny 
        D = 1./ sqrt( pos_x ** 2 + pos_y **2 )
        D[isinf(D)] = amax(D[~isinf(D)])
        #D = blur_image(D, 1)
        
        pos = hstack([reshape(pos_x*dx, (-1,1)), reshape(pos_y*dy, (-1,1))])
        #pos = pos / repmat(sqrt(sum(pos**2,1)), 2,1).T
        pos = pos/sqrt(sum(pos**2,1))[:,None]
        proj_angle = abs(dot(pos, Nvec))
        proj_angle = reshape(proj_angle, (tok_nx, tok_ny))
        D = D*proj_angle
        
        #D[D>0.5]  = 0.5
        
        D[reshape(BdMat, (tok_ny, tok_ny), order="F")] = 0
        
        #imshow(D, extent=tokamak.area_axis, interpolation = "nearest" )
        #colorbar()
        #plot(bd[:,0], bd[:,1])
        #show()
        
        #D[D < 0 ] = 0
        return D

    X,Y = meshgrid(arange(tok_nx), arange(tok_ny) )

    try:
        diff_model_full = load( tokamak.geometry_path + "diff_model_full.npy")
    except:
        print("diff model")
        diff_model_full = zeros((npix, nl), dtype="single")

        for i in range(nl):
            #print i
            r = sqrt(Xch[1,i]**2 + Zch[1,i]**2)
            y = Ych[1,i]
            model = get_diff_model(r,y)
            diff_model_full[:,i] = reshape( model , (-1), order="F")
        diff_model_full[diff_model_full < 0] = 0
        limit = amax(diff_model_full)/5
        diff_model_full[ diff_model_full >  limit ] = limit  # remove too strong points 
        save(tokamak.geometry_path +  'diff_model_full.npy', diff_model_full)
    
    model = sum(diff_model_full, 0)

    
    

    print("generate")
    diff_model_G = zeros((npix, npix), dtype="single")


    for i in arange(tok_nx):
        for j in arange(tok_ny):

            model = get_diff_model(tokamak.xgrid[i],tokamak.ygrid[j])
        
            model[isnan(model)] = 0

            diff_model_G[:,i*tok_nx+j] = reshape( model , (-1), order="F")

    print("plotting")
    max_g = nanmax(diff_model_G)
    diff_model_G[ diff_model_G > max_g/20 ] = max_g/20
    

    print("refl model")
    
    #Nrfl = 1
    
    #print tokamak.camera_res
    #exit()
    
    #Xch_rlf = zeros((Nrfl+1, nl))
    #Ych_rlf = zeros((Nrfl+1, nl))
    #Zch_rlf = zeros((Nrfl+1, nl))
    

    #from mpl_toolkits.mplot3d import Axes3D
    #fig = figure()
    #ax = fig.gca(projection='3d')

    #for j in range(Nrfl):
        #if j == 0:
        
    try:
        xxx
        Tmat_reflections = load(tokamak.geometry_path + 'Tmat_reflections.npy').item()
        projections  = load(tokamak.geometry_path + 'projections.npy')
        Tmat_diffusion  = load(tokamak.geometry_path + 'Tmat_diffusion.npy')

    except:
        #raise
    
        t = time.time()
        #Tmat_reflections = zeros( (npix, nl) , dtype="single" )
        Tmat_reflections = sparse.lil_matrix((npix, nl), dtype="single")
        #print shape(Tmat_reflections), Tmat_reflections,  (npix, nl)
        #exit()
        
        Tmat_diffusion = zeros((npix, nl), dtype="single")

        
        
        projections = zeros(nl)
        norm_angle = zeros((3, nl))

        ##########################
        vchords = 20
        b = 0.5  # cone angle 
        alpha  = 4   # phong model !! 
        ##############################
        
        #b *= linspace(-1,1, vchords )
        #b *= random.rand(vchords, vchords)*2-1
        #theta0 = arccos(1-b*random.rand(vchords, vchords))
        #phi0 = 2*pi*random.rand(vchords, vchords)
        
        #x = linspace(1./vchords,1, vchords)
        #[X,Y] = meshgrid(x,x)
        X = random.rand(vchords, vchords)
        Y = random.rand(vchords, vchords)
        theta0 = arccos(1-b*X)
        phi0 = 2*pi*Y
        
        #from mpl_toolkits.mplot3d import Axes3D
        #fig = figure()
        #ax = fig.gca(projection='3d')

        vec = zeros((vchords**2, 3))
        vec[:,2] = sin(ravel(phi0))*sin(ravel(theta0))
        vec[:,1] = cos(ravel(phi0))*sin(ravel(theta0))
        vec[:,0] = cos(ravel(theta0))

        
        #ax.plot(vec[:,0],vec[:,1], vec[:,2] , 'bo')
        #show()
        
        
        chords_tmp = zeros((2, nl, 3), dtype="single")
        chords_tmp[:,:,0] = Xch
        chords_tmp[:,:,1] = Ych
        chords_tmp[:,:,2] = Zch

        for i in range(nl):
            #if mod(i, 121): continue
            chords = zeros(( Nsteps,vchords**2, 3 ), dtype="single") + nan
            weights_chords = zeros((vchords**2), dtype="single") + nan

            start = chords_tmp[0,i,:]
            end = chords_tmp[1,i,:]
            Vvec = end - start
            Vvec /= norm(Vvec)
            
            R = sqrt(end[0]**2 + end[2]**2)
            Y = end[1]
            d = (R - bd[:,0])**2 + (Y - bd[:,1])**2 
            phi = arctan2(end[0] , end[2])
            #print phi
            s = sort(argsort(d)[:2])
            
            #plot(R,  chords_all[j][1,i,1], 'o')
            #plot( bd[s,0],  bd[s,1], '.')
            #plot(bd[:,0], bd[:,1])
            #show()
            
            #r = sqrt(X[:,k]**2 + Z[:,k]**2)
            #y = Y[:,k]
            #tok_nx

        
            Nvec = array( [bd[s[1],0] -bd[s[0],0] , bd[s[1],1]-bd[s[0],1] ])
            Nvec[0] *= -1
            Nvec[[0,1]] = Nvec[[1,0]]
            #print Nvec[0]
            #print bd[s[0],0]+Nvec[0]
            #print ([bd[s[0],0], (bd[s[0],0]+Nvec[0]) ])
            #X = array([bd[s[0],0], (bd[s[0],0]+Nvec[0]) ])
            #Y =  array([bd[s[1],1],(bd[s[1],1]+Nvec[1])])
            #plot( X , Y , '-')



            Nvec = array([sin(phi)*Nvec[0], Nvec[1], cos(phi)*Nvec[0]])
            Nvec /= norm(Nvec)
            
            #Nvec = Nvec*abs(1+0.05*random.randn(3))
            
            end = array([sin(phi)*bd[s,0], bd[s,1],  cos(phi)*bd[s,0] ])
            end = end[:,0]
            #end = array()

            
            
            ##print end
            #Nvec = end[:,1] - end[:,0]
            
            #plot([bd[s[0],0], bd[s[0],0] ] , [bd[s[1],0],bd[s[0],0]], '-')
            #Nvec[0] *= -1  # make it perpendicular
            #Nvec[1] *= -1  # make it perpendicular

            #print Nvec
            #exit()
            #ax.plot([Nvec[0]], [Nvec[1]],[ Nvec[2]], 'b.')
            #ax.plot([Vvec[0]], [Vvec[1]],[ Vvec[2]], 'r.')

            #tmp = Nvec / 50
            #ax.plot([end[0], end[0] + tmp[0]], [end[1], end[1]+ tmp[1]], [end[2], end[2]+ tmp[2]], 'b-')
            #ax.plot([end[0]], [end[1]], [end[2]], 'b.')
            #ax.plot([Xch[-1,i]], [Ych[-1,i]], [Zch[-1,i]], 'r.')

            #plot(sin(phi), cos(phi), '.')
            #ax.plot(end[0], end[1], end[2], 'b.')
            #if mod(i, 5) == 0:
            
            #from mpl_toolkits.mplot3d import Axes3D
            #fig = figure()
            #ax = fig.gca(projection='3d')

            #ax.plot([0,Xch[-1,i]] , [0,Ych[-1,i]], [0,Zch[-1,i]],  'r:o')
            
            
            #print abs(sum(Vvec*Nvec)), norm(Vvec), norm(Nvec)
            #norm_angle[:,i]= Nvec
            #print angle[i]
            #print end
            
            projection  =   dot(Vvec, Nvec) /  ( norm(Vvec)*norm(Nvec) ) 
            projections[i] = abs(projection)**alpha
            
            #print "sum(Vvec*Nvec)", sum(Vvec*Nvec)
            Rvec = ( Vvec - 2*(sum(Vvec*Nvec)) *Nvec ) *2
            Rvec /= norm(Rvec)
            #r = sqrt( sum(Rvec **2))
            r  = 1
            phi = arctan2( Rvec[0], Rvec[1] )
            theta = arctan2( Rvec[2], sqrt(Rvec[0]**2 + Rvec[1]**2) )


            #ax.plot([end[0],end[0]+Rvec[0]],[end[1],end[1]+Rvec[1]],[end[2],end[2]+Rvec[2]], 'g-o')
            #ax.plot([0,r*sin(phi)*cos(theta)],[0,r*cos(phi)*cos(theta)],[0,r*sin(theta)], 'g:o')
            
            #show()
            #phi0[0,0] = 0
            #theta0[0,0] = 0
            vec[0,:] = zeros(3)
            vec[0,0] = 1

            
            #rot = zeros((3,3))
            #rot[0,:] = [ sin(theta)*sin(phi), r*cos(theta)*sin(phi), -r*sin(theta)*cos(phi)]
            #rot[1,:] = [ sin(theta)*cos(phi), r*cos(theta)*cos(phi), r*sin(theta)*sin(phi)]
            #rot[2,:] = [ cos(theta), -r*sin(theta), 0]

            rot = zeros((3,3))
            rot[0,:] = [ cos(theta)*sin(phi), r*sin(theta)*sin(phi), -r*cos(theta)*cos(phi)]
            rot[1,:] = [ cos(theta)*cos(phi), r*sin(theta)*cos(phi), r*cos(theta)*sin(phi)]
            rot[2,:] = [ sin(theta), -r*cos(theta), 0]

            Lvec =  dot(rot, vec.T * r).T
            
            #Lvec =  repmat(copy(Rvec), vchords**2, 1)

            #print "phi", phi, "theta", theta, "refl", Rvec, "Lvec", Lvec


            #ax.plot(end[0]+Lvec[:,0],end[1]+Lvec[:,1],end[2]+Lvec[:,2] , 'bo')
            ##ax.plot(end[0]+vec[:,0],end[1]+vec[:,1], end[2]+vec[:,2] , 'r.')

            #ax.plot([end[0],end[0]+Rvec[0]],[end[1],end[1]+Rvec[1]],[end[2],end[2]+Rvec[2]], 'g-')
            
            #show()
                
            for k in range(vchords**2):
                if dot(Lvec[k,:],Nvec) < 0:
                    continue
                
                tmp_chord = array([end,  end+Lvec[k,:]*(xmax - xmin)*2 ]) # (Vvec - Vvec*parvec)/20 ])
                
                chords[:, k, 0] = interp1d([0,1], tmp_chord[:,0],assume_sorted=True )(linspace(0,1,Nsteps)).T 
                chords[:, k, 1] = interp1d([0,1], tmp_chord[:,1],assume_sorted=True )(linspace(0,1,Nsteps)).T 
                chords[:, k, 2] = interp1d([0,1], tmp_chord[:,2],assume_sorted=True )(linspace(0,1,Nsteps)).T 
                weights_chords[k] = abs(dot(  Lvec[k,:] , Rvec ))**alpha #  / max( dot(Lvec[k,:],Nvec), abs(dot(Vvec,Nvec)), 0.5)



            out = clean_chords(chords[:,:,0],chords[:,:,1],chords[:,:,2])
            chords = zeros((2, vchords**2, 3))
            for k in range(3):
                chords[:,:,k] = out[k]
                
            chords_0 = zeros((Nsteps, vchords**2,3))
            


            # !!  SIMULATE REFLECTED CHORDS AND REPLACE ILUMINATION BY DIFF. LIGHT 
            for k in range(vchords**2):
                chords_0[:, k, 0] = interp1d([0,1], chords[:,k,0] )(linspace(0,1,Nsteps)).T 
                chords_0[:, k, 1] = interp1d([0,1], chords[:,k,1] )(linspace(0,1,Nsteps)).T 
                chords_0[:, k, 2] = interp1d([0,1], chords[:,k,2] )(linspace(0,1,Nsteps)).T 

                r = sqrt( chords[1,k,0]**2 + chords[1,k,2]**2 )
                y = chords[1,k,1]
                if isnan(r):
                    #print "nan r", r
                    continue
                r = int((r-xmin)/(xmax-xmin)*tok_nx)
                y = int((y-ymin)/(ymax-ymin)*tok_ny)
                
                #print r, y 
                
                #print r,y, r*tok_nx+y, shape(diff_model_G)
                #print Tmat_diffusion[:,i]
                #print diff_model_G[:,int(r*tok_nx+y)]
                #print "============"
                #imshow(reshape(diff_model_G[:,int(r*tok_nx+y)], (tok_ny, tok_nx), order="F"),  extent=tokamak.area_axis)
                #plot(sqrt(chords_0[:,k,0]**2+chords_0[:,k,2]**2), -chords_0[:,k,1])
                ##plot(r,y, 'o')
                #show()
                
                #print shape(
                              
                Tmat_diffusion[:,i] += diff_model_G[:,int(r*tok_nx+y)]
                
            #print shape(Tmat_diffusion[:,i]), array(Tmat_diffusion[:,i]), (tok_ny, tok_nx)
            #imshow( reshape(Tmat_diffusion[:,i], (tok_ny, tok_nx), order="F"), extent=tokamak.area_axis)
            #plot(sqrt(chords_0[:,:,0]**2+chords_0[:,:,2]**2), -chords_0[:,:,1])
            #show()
                
            #from mpl_toolkits.mplot3d import Axes3D
            #fig = figure()
            #ax = fig.gca(projection='3d')

            #ax.plot(chords[:,:,0], chords[:,:,1],chords[:,:,2],'b-')
            #show()
            
            ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Tmat_tmp = geom_mat_gen(chords[:,:,0],chords[:,:,1],chords[:,:,2],
                                tokamak,tokamak.nx*zoom,tokamak.ny*zoom, vchords**2, tokamak.geometry_axis )
            Tmat = sparse.csc_matrix(Tmat)
            Tmat_tmp = (Tmat_tmp * sparse.spdiags(weights_chords, 0, vchords**2, vchords**2)).sum(axis=1) # change weight according to reflection angle deviance 
            
            #xxxxx  # opravit projections[i] ??  
            
            Tmat_reflections[:,i] = Tmat_tmp * projections[i]  #   zkusim menit vahu podle uhlu projekce !!!!!!!
            

            #Tmat_tmp = Tmat_tmp.toarray()
            #print shape(Tmat_tmp), tokamak.camera_res[0]*tokamak.camera_res[1]
            #pcolor( reshape( sum( Tmat_tmp, 0) , (tokamak.camera_res[0]*cam_zoom,tokamak.camera_res[1]*cam_zoom)))
            ####
            #colorbar()
            #show()
            ##n = size(Tmat,1)
            #figure()
            #print (tokamak.ny, tokamak.nx)
            
            #print "pcolor"
            #imshow( reshape(  Tmat_tmp , (tokamak.ny*zoom, tokamak.nx*zoom), order="F" ), interpolation="nearest")
            #colorbar()
            #show()
                    
            print("i%d/%d"%(i, nl))

        Tmat_reflections = sparse.csc_matrix(Tmat_reflections, dtype="single")
        save(tokamak.geometry_path + '/Tmat_reflections.npy', Tmat_reflections)
        

        save(tokamak.geometry_path + '/Tmat_diffusion.npy', Tmat_diffusion)


        
        save(tokamak.geometry_path + '/projections.npy', projections)
        
        
        
        #save('weights_chords.npy', weights_chords)

        #pcolor(weights_chords)
        #colorbar()
        #show()
        print("time %f"%(t-time.time()))
        
        #chords = reshape( chords, ( Nsteps , nl* vchords**2, 3 ))
        
        #print shape(chords)
        #out = clean_chords(chords[:,:,0],chords[:,:,1],chords[:,:,2])
        #print shape(out[1])
        
        ##ind = range(j,j+2)
        #chords = zeros((2, nl* vchords**2, 3))
        #for k in range(3):
            #chords[:,:,k] = out[k]
            
        #Nsteps  = 100
            
        #for i in range(size(chords,1)):
            #if mod(i, 31): continue
            #print(shape(chords[:,i,:].T))
            #ch = interp1d([0,1], chords[:,i,:].T)(linspace(0,1,Nsteps)).T
            #print ch
            #print  chords[:,i,:]
            ##exit()
            #plot(sqrt(ch[:,0]**2+ch[:,2]**2),ch[:,1], '.')
            #plot(sqrt(chords[0,i,1]**2+chords[0,i,2]**2), chords[0,i,1], '.')
            #plot(sqrt(chords[1,i,1]**2+chords[1,i,2]**2), chords[1,i,1], '.')
            #plot(bd[:,0], bd[:,1])
            #show()
            
        #Xch_rlf = Xch_rlf[:2]
        #Ych_rlf = Ych_rlf[:2]
        #Zch_rlf = Zch_rlf[:2]

        #tokamak.diff_model = projections
        #config.diff_model = projections
        
        #save('diff_model.npy', projections )
        
        ##for i in range(3):
            ##x = reshape(norm_angle[i,:], (tokamak.camera_res[0]*cam_zoom,tokamak.camera_res[1]*cam_zoom))
            ##x = rot90(x)
            ##imshow(x, interpolation="nearest")
            ##colorbar()
            ##show()
        
        
        #projections = reshape(projections, (tokamak.camera_res[0]*cam_zoom,tokamak.camera_res[1]*cam_zoom))
        #projections = rot90(projections)
        
        ##print projections
        
        #imshow(projections, interpolation="nearest")
        #colorbar()
        #title("amgle")
        #show()
        #exit()
        
        #diff_model_full *=  ravel(angle)
        #save('diff_model_full.npy', diff_model_full )

        #imshow(diff_model_full)
        #show()
        
        #save('chords.npy', chords)
        
        
    #Xch_rlf[ind,:] = out[0]
    #Ych_rlf[ind,:] = out[1]
    #Zch_rlf[ind,:] = out[2]
    
    #chords[:,:,0] = Xch_rlf[ind,:]
    #chords[:,:,1] = Ych_rlf[ind,:]
    #chords[:,:,2] = Zch_rlf[ind,:]
    
    #chords_all.append(chords_tmp)
        
    
        
    #print shape(Tmat_reflections), shape(Tmat_reflections.sum(0)), tokamak.camera_res
    T = array(reshape(Tmat_reflections.sum(0), tokamak.camera_res))
    
    #imshow(T, interpolation="nearest")
    #show()
    
    #print projections
    #print shape(projections), tokamak.camera_res
    T = array(reshape(projections, tokamak.camera_res))
    
    #imshow(T, interpolation="nearest")
    #show()
    
    
    #return Tmat, sqrt(Xch**2+Ych**2), Ych

    
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = figure()
    #ax = fig.add_subplot(111, projections='3d')

    #for i in range(nl):
        #if mod(i, 51): continue
        #ax.plot(Xch_tmp[:,i],Ych_tmp[:,i],Zch_tmp[:,i], 'b-')
        #ax.plot(Xch[:,i],Ych[:,i],Zch[:,i], 'r-')

    #show()
    
    #show()
    #exit()
   
    #print shape(Xch_tmp)
    #Xch_tmp =  interp1d([0,1], Xch_tmp.T)(linspace(0,1,Nsteps)).T
    #Ych_tmp =  interp1d([0,1], Ych_tmp.T)(linspace(0,1,Nsteps)).T
    #Zch_tmp =  interp1d([0,1], Zch_tmp.T)(linspace(0,1,Nsteps)).T
    #Xch = vstack([Xch, Xch_tmp])
    #Ych = vstack([Ych, Ych_tmp])
    #Zch = vstack([Zch, Zch_tmp])
    
    #print shape(Xch), shape(Zch), shape(Zch_tmp), Nsteps
    #exit()
    
    
    print("geom_mat_gen")
    #from geom_mat_gen import geom_mat_gen


    #weights_chords  = reshape(weights_chords, nl* vchords**2)

    #try:
        #xxx
        #Tmat_reflections = load('Tmat_reflections.npy')
    #except:
        #Tmat_reflections = geom_mat_gen(chords[:,:,0],chords[:,:,1],chords[:,:,2], tokamak,tokamak.nx*zoom,tokamak.ny*zoom, nl*vchords**2, tokamak.geometry_axis )
        #save('Tmat_reflections.npy', Tmat_reflections)

        #Tmat_reflections = Tmat_reflections * sparse.spdiags(weights_chords, 0,  nl* vchords**2,  nl* vchords**2)
        #Tmat_reflections =   sum(reshape( Tmat_reflections.toarray() , (npix, nl, vchords**2)), 2) / vchords**2
        #Tmat_reflections = sparse.csc_matrix(Tmat_reflections, dtype="single")
        #save('Tmat_reflections.npy', Tmat_reflections)

        
        
        
    ## !!!  tahle 2D verze nefunguje protože to nepočítá správné 3D vzdálenosti !!!!
    ##Tmat = geom_mat_gen(X,Y,Z, tokamak,tokamak.nx*zoom,tokamak.ny*zoom,nl)
    
    
    Tmat = geom_mat_gen(Xch, Ych, Zch, tokamak,tokamak.nx*zoom,tokamak.ny*zoom, nl, tokamak.geometry_axis )

    

            
    
    #Tmat = Tmat_all[0] # + 0.3*Tmat_all[1] #  + 0.1*Tmat_all[2]
    
    #print shape(Tmat_all[1])
    #imshow(reshape(Tmat.sum( axis=0 ), tokamak.camera_res), interpolation="nearest")
    #show()
    #exit()
    
           
    
    #tokamak.Tmat_reflections = Tmat_reflections
    
    #print "????????????????????????????????????"
    
      #Tmat = Tmat_all[1]
    

    #Tmat = Tmat.toarray()
    #print shape(Tmat), tokamak.camera_res[0]*tokamak.camera_res[1]
    #pcolor( reshape( sum( Tmat, 0) , (tokamak.camera_res[0]*cam_zoom,tokamak.camera_res[1]*cam_zoom)))
    ####
    #colorbar()
    #show()
    ###n = size(Tmat,1)
    ##figure()
    ##print (tokamak.ny, tokamak.nx)
    #pcolor( reshape( sum( Tmat, 1 ), (tokamak.ny*zoom, tokamak.nx*zoom), order="F") )
    #colorbar()
    #show()
    ##,Ych2,Zch2, tokamak,tokamak.nx*zoom,tokamak.ny*zoom,nl)





    #fig = figure()
    #ax = fig.add_subplot(111, projections='3d')
    #ax.plot(Xch,Ych,Zch, 'b-')
    #show()
    ##exit()
    

    
    #exit()
    

    


    #plot(Rch[-1,:], Ych[-1,:], '.')
    #plot(bd[:,0], bd[:,1])
    #show()
    #exit()

    ###print tokamak.camera_res
    ###print cam_zoom, shape(D)
    #D = reshape(D, (tokamak.camera_res[0]*cam_zoom,tokamak.camera_res[1]*cam_zoom))
    #D = rot90(D)

    #imshow(D, interpolation="nearest")
    #colorbar()
    #show()

    ###print " saving "
    ###save('D', D)
    ###savemat('Xch', {'Xch': Xch[-1,:]})
    ###savemat('Ych', {'Ych': Ych[-1,:]})
    ###savemat('Zch', {'Zch': Zch[-1,:]})
    

    ###img = double(imread('geometry_compass/camera/img_2851_12.jpg'))
    ###DOWN=8/tokamak.cam_zoom
    ####img = double(imread('geometry_compass/camera/image.jpg'))
    ###img = img[::DOWN,::DOWN]
    ####img = loadmat('geometry_compass/camera/img_background.mat' )['img_background']
    ####img = img[::2,::2]
    ####img = double(imread('img_smaller_36.png'))
    ####img = img[::DOWN,::DOWN]
    
    ####img = loadtxt('geometry_compass/camera/image.txt')
    ####img = img[::DOWN,::DOWN]
    ###FOV = load('geometry_compass/camera/FOV.npy')
    ###FOV = FOV[::DOWN,::DOWN]
    ###Schema = loadmat('FOV_Schema.mat')['FOV_Schema']
    ###Schema = Schema[::DOWN,::DOWN]

    ####subplot(2,1,1)
    ####imshow(D)
    ####subplot(2,1,2)
    ####imshow(Schema)
    ####show()
    ####print shape(FOV), shape(Schema), shape(D)
    #####exit()

    ###try:
        ###D[FOV==0] = nan
    ###except:
        ###pass
    
    ####print shape(D), shape(FOV)

    ###D[isinf(D)]  = nan
    ###D[isnan(D)]  = nan

    ####imshow(D*100, interpolation="nearest")
    ####colorbar()
    ####show()
    ###save('D_2', D)

    ###from scipy.signal import medfilt2d
    ###D = medfilt2d(D, 3)
    
    #img /= nanmean(img)
    #D /= nanmean(D)
    #print D
    #img[FOV==0] = 0
    #imshow( img - D, interpolation="nearest")
    #print "plotting"
    #D[Schema != 0 ] = nanmax(D)/2
    #imshow(  D, interpolation="nearest")
    #c = colorbar()
    #c.set_label('Distance from camera [m]')
    #axis('off')
    #savefig( 'distance.png' )
    #show()

    #exit()




    ########################################xx
    ##try:
    #D = D - nanmin(D)
    ##imshow(D); colorbar()
    ##show()
    #D = log(1+1e5*abs(D))
    #D[D>0] -= mquantiles(D[D>0], 0.1)
    #img /= nanmean(abs(img))
    #D /= nanmean(abs(D))
    #subplot(2,2,1)
    #imshow(abs(D), interpolation="nearest")#; colorbar()
    #subplot(2,2,2)
    #imshow(abs(img), interpolation="nearest")#; colorbar()
    #subplot(2,2,3)
    #imshow( img - D, interpolation="nearest")#; colorbar()
    #subplot(2,2,4)
    #plot(sqrt(tokamak.cam_coord[0]**2 + tokamak.cam_coord[2]**2), tokamak.cam_coord[1], 'ow')
    #plot(tokamak.ports[:,0], tokamak.ports[:,1], '.')
    #colorbar()
    #axis('equal')
            
    #show()
    
    ##imshow( img - D)

    ##except:
        ##pass
    #show()

    #exit()

    ########################################xx


    
    #del X,Z #,Y



    ##exit()
    
    N = 100

    Xch =  interp1d([0,1], Xch.T)(linspace(0,1,N)).T
    Ych =  interp1d([0,1], Ych.T)(linspace(0,1,N)).T
    Zch =  interp1d([0,1], Zch.T)(linspace(0,1,N)).T
    Rch = hypot(Xch, Zch)
    
    
    #Tmat = Tmat.T
    #zoom = tokamak.pixel_zoom
    #cam_zoom = tokamak.cam_zoom
    if zoom > 1:
        from orthogonal_trans import zoom_matrix
        Q = zoom_matrix(tokamak.nx, tokamak.ny, zoom)
        #print shape(Q.T), shape(Tmat)
        Tmat = Q.T*Tmat
        #print shape(Tmat)
        del Q
        
    if cam_zoom > 1:
        from orthogonal_trans import zoom_matrix
        Q = zoom_matrix(tokamak.camera_res[0], tokamak.camera_res[1], cam_zoom)
        #print shape(Tmat), shape(Q)
        Tmat = Tmat*Q
        del Q
        Xch = Xch[:,::cam_zoom**2]
        Ych = Ych[:,::cam_zoom**2]
        Rch = Rch[:,::cam_zoom**2]


    #Tmat = Tmat.T

    #save('Tmat', Tmat)
    
    ### ========== vekreslit ty pcolor grafy Tmat ==============

    ##save('Tmat', Tmat)
    #print "plotting (to array)"
    ##n = size(Tmat,0)
    #figure()
    #imshow( reshape( array(Tmat.sum(1)) , (tokamak.camera_res[0],tokamak.camera_res[1])), interpolation="nearest")
    #colorbar()
    #show()
    ##n = size(Tmat,1)
    ##figure()9500
    ##print (tokamak.ny, tokamak.nx)
    #imshow( reshape( array(Tmat.sum(0)) , (tokamak.ny, tokamak.nx), order="F") , interpolation="nearest")
    #colorbar()
    #show()

    ##figure()
    #imshow( reshape( sum( Tmat.toarray()[::401,:], 0 ), (tokamak.ny, tokamak.nx) , order="F"), interpolation="nearest" )
    #show()
    ###==========  ==============
    #Tmat = Tmat.T

    #exit()
    


    #print "shape(Tmat)", shape(Tmat)

    #plot(Rch, Ych)
    #show()
    ##exit()

    #try:
        #Tmat = Tmat.toarray()
    #except:
        #pass
    
    ###########################################
    return Tmat ,Rch, Ych
    ##########################################



        

    

    #savez('plot2', X=Xpix, Y=Ypix, Z=Zpix)

    #X2= X.copy()
    #Y2= Y.copy()
    #Z2= Z.copy()


#udelat zakroucené chordy  !!!

    #save('Xch', Xch)
    #save('Ych', Ych)
    #save('Zch', Zch)
    #save('Xpix', Xpix)
    #save('Ypix', Ypix)
    #save('Zpix', Zpix)


    ###      plotting of chords
   ##=====================================================
    #print "plotting"
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = figure()
    #ax = fig.add_subplot(111, projections='3d')
    #nl = size(Xch, 1)
    #ind = random.permutation(nl)
    #ind = ind[1:min(400,nl)]
    ##range(size(Xpix,1)):
    ##ind = arange(nl)
    #for i in ind:
        #print i, nl
        ##if mod(i,231) != 0:
            ##continue
            
        ##print i
        ## plot whole line 
        ## ax.plot( X[:,i],Y[:,i],Z[:,i], 'r') # , linestyle='None', marker='o', markersize = 5, antialiased=True)
        ## plot only the internal part of the line 
        #ax.plot( Xch[:,i],Zch[:,i],Ych[:,i], 'b') 
        ## , linestyle='None', marker='o', markersize = 5, antialiased=True)

    #for i in range(size(Xpix,1)):
        ##print shape(Xpix)
        #ax.plot( Xpix[:,i],Zpix[:,i],Ypix[:,i], 'k-', linewidth=0.5)

    ##ax.set_xlim3d(-xmax, xmax)
    ##ax.set_ylim3d(ymin, ymax)
    ##ax.set_zlim3d(-ymax, ymax)
    #ax.set_xlabel('R [m]')
    #ax.set_ylabel('R [m]')
    #ax.set_zlabel('Z [m]')


    #show()

    #plot(Zpix,Xpix)
    ##print 
    #show()
    

    #print shape(Rch), shape(Ych)

    #plot(Rch, Ych)
    #show()
    #exit()
    

    #print "pliotting"
    ##from mpl_toolkits.mplot3d import Axes3D
    ##fig = figure()
    ##ax = fig.add_subplot(111, projections='3d')


    #exit()
    ##=====================================================

    #gc.collect()

    #print X,Y,Zphi2
    #print where(isnan(X)), where(isnan(Y)), where(isnan(Z))

    print("generate chords")
    
    nx = 50
    ny = 50
    nz = 50
    #npix 
    name = 'chords_' + str(nx)+str(npix) + str(zoom)+".npy"
    try:
        xxx
        Tmat_1 = load( name ).item()
    except:
        Tmat_1 = geom_mat_gen_3D(Xch,Ych,Zch, tokamak,nx,ny,nz,nl, coord)
        save(name, Tmat_1)

    if cam_zoom > 1:
        from orthogonal_trans import zoom_matrix
        Q = zoom_matrix(tokamak.camera_res[0], tokamak.camera_res[1], cam_zoom)
        Tmat_1 = Tmat_1*Q
        X = X[:,::cam_zoom]
        Y = Y[:,::cam_zoom]
        Z = Z[:,::cam_zoom]
        del Q
    gc.collect()

    
    T = array( Tmat_1.sum(axis = 0) )
    T = reshape(T, (tokamak.camera_res[1], tokamak.camera_res[0]) , order="F")

   
    
    ######===========================================================

    print("generate pixels")
    name = 'pixels_' + str(nx)+str(npix) + str(zoom)+".npy"
    try:
        xxx
        Tmat_2 = load( name).item()
    except:
        Tmat_2 = geom_mat_gen_3D(Xpix,Ypix,Zpix, tokamak,nx,ny,nz, npix, coord)
        save(name, Tmat_2)
        
    if zoom > 1:
        from orthogonal_trans import zoom_matrix
        Q = zoom_matrix(tokamak.nx, tokamak.ny, zoom)
        Tmat_2 = Tmat_2*Q
        del Q
    gc.collect()
    
    ######save('tmat', Tmat_2)
    ######exit()
    ######print "plotting"
    T = array(Tmat_2.sum(axis = 0))
    #print("=========================")
    T = reshape(T, (tok_ny, tok_nx), order="F")
    #print((shape(T), sum(T)))
    #print(T)
    #print((where(T> 0)))
    
    
    #pcolor( T  )
    #colorbar()
    #show()
    
    
    
    
    #for i in range(nz):
        #print i
        #pcolor( reshape( array( Tmat_2.sum(axis = 1) ), (nx, ny, nz))[:,i,:]  )
        #plot(X,Y, 'r')
        
        #show()
        #exit()
        ######savefig("rez_"+str(i)+'.png')



    ######exit()

    ######print shape(Tmat_1)

    ######print shape(Tmat_2)

    #####print "saving"
    #####save('Tmat_1', Tmat_1)
    #####save('Tmat_2', Tmat_2)
    #####print "done"
    
    #exit()

    #Tmat_1 = load('Tmat_1.npy').item()
    #Tmat_2 = load('Tmat_2.npy').item()

    #print "smooting"
    
    #Tmat_1 = smooth_3D_mat(Tmat_1, nx, ny, nz)
    #Tmat_2 = smooth_3D_mat(Tmat_2, nx, ny, nz)
    
    #print "smooting2"

    #Tmat_1 = smooth_3D_mat(Tmat_1, nx, ny, nz)
    #Tmat_2 = smooth_3D_mat(Tmat_2, nx, ny, nz)

       
    #print("multiplication")
    #print (Tmat_1.shape()),(Tmat_2.shape())
    #print((shape(Tmat_1), shape(Tmat_2)))
    
    Tmat = Tmat_1.T * Tmat_2
    Tmat = Tmat.astype(float64)
   
    #print((Tmat_1.sum().sum()))
    #print((Tmat_2.sum().sum()))

    #print((Tmat.sum().sum()))

    #print shape(Tmat)

    #print nx, ny, nz, nl

    #exit()
    #print shape(Tmat)

    #print("plotting")
    n = size(Tmat,0)
    #pcolor( reshape( sum( Tmat.toarray(), 1) , (tokamak.camera_res[0],tokamak.camera_res[1]) ))
    ##(tokamak.camera_res[0],tokamak.camera_res[1])))
    #colorbar()
    #show()

    #n = size(Tmat,0)
    #pcolor( reshape( sum( Tmat.toarray()[:,::401], 1) , (sqrt(n), sqrt(n))))
    ###(tokamak.camera_res[0],tokamak.camera_res[1])))
    #colorbar()
    #show()

    #print shape(Tmat), (tokamak.ny, tokamak.nx), shape( sum( Tmat.toarray(), 0) )
    n = size(Tmat,1)
    
    #pcolor( reshape( sum( Tmat.toarray(), 0) , (tokamak.ny, tokamak.nx), order="F"))
    #colorbar()
    #show()

    #n = size(Tmat,1)
    #pcolor( reshape( sum( Tmat.toarray()[::401,:], 0) , (tokamak.ny, tokamak.nx), order="F"))
    #colorbar()
    #show()
    #print('konec')
    exit()
    
    #print "plotting"
    #for i in range(nz):
        #print i
        #pcolor( tokamak.xgrid, tokamak.ygrid, reshape( array( Tmat_1.sum(axis = 1) ), (nx, ny, nz))[:,:,i]  )
        #plot(X,Y, 'r')
        #axis([tokamak.xmin,tokamak.xmax,tokamak.ymin,tokamak.ymax])
        ##show()
        ##exit()
        #savefig("rez_"+str(i)+'.png')
    #fig = figure()
    #import mpl_toolkits.mplot3d.axes3d as p3
    #ax = p3.Axes3D(fig)
    #ax.plot_wireframe(x,y,z)
    #show()
    gc.collect()


    #print("koncim")
    #exit()
    return Tmat, R, Y




def smooth_3D_mat(mat, nx, ny, nz):   #(mat, ker, nx, ny, nz):

    #mat = load('Tmat_2.npy').item()
    #nx = 100
    #ny = 100
    #nz = 100
    ker = 0.5
    N = size(mat,0)
    from time import sleep

    mat0 = mat
    print("tolil")
    #mat = sparse.lil_matrix(shape(mat))
    print("done")
    
    #print shape(mat)
    iker = ceil(2*ker)
    ind = arange(-iker,iker+1, dtype=int)

    [X,Y,Z] = mgrid[-iker:(iker+1), -iker:(iker+1), -iker:(iker+1)]

    #[X,Y] = meshgrid(ind,ind)
    kernel = exp( - (X**2 + Y**2 + Z**2) / ker**2 )
    kernel /= sum(kernel)
    kernel = ravel(kernel)

    #print kernel
    
    ind = ravel(X + Y * nz + Z * ny * nz)

    print(("spdiag", len(ind)))

    #exit()

    m = dot(ones( (N, len(ind)) ) ,   diag(kernel)  )
    m = array(m).T

    sleep(2)
    
    #print shape(m)
    #print shape(ind)
    
    smooth_mat = spdiags(  m , ind, N,N )

    sleep(2)
    
    #print shape(smooth_mat), shape(mat)

    mat = smooth_mat * mat 

    return mat




def load_geom_3D(geometry_path,cameras):
    """
    Load geometry from ``geometry_path`` and prepare for processing.

    *Expected names of detectors are detector_[num]_[x,y].txt. Each file has 2-3 columns. If there are 2, then the first is beginning of chord and the second is end, it is only X or Y coordinate. If there are 3 columns, the the first two are range for virtual chords and the last is the other side (position of pinhole)*
    """

    print("loading geom")
     #========================  Locations of detectors  ========================
    ychords = list()
    xchords = list()
    zchords = list()
    
    for i,c in enumerate(cameras):
        try:
            xchords.append(load(geometry_path+'/detector_'+c+'_x.npy'))
            ychords.append(load(geometry_path+'/detector_'+c+'_y.npy'))
            zchords.append(load(geometry_path+'/detector_'+c+'_z.npy'))
        except:
            xchords.append(loadtxt(geometry_path+'/detector_'+c+'_x.txt'))
            ychords.append(loadtxt(geometry_path+'/detector_'+c+'_y.txt'))
            zchords.append(loadtxt(geometry_path+'/detector_'+c+'_z.txt'))
        
        
    #i = 0
    #print geometry_path+'/detector_'+str(i)+'_z.txt'
    #while os.path.isfile(geometry_path+'/detector_'+str(i)+'_z.txt') or os.path.isfile(geometry_path+'/detector_'+str(i)+'_z.npy'):
        #print i

        #i += 1
    xchords = squeeze(array(xchords)).T
    ychords = squeeze(array(ychords)).T
    zchords = squeeze(array(zchords)).T
    


    print("loading geom done")

    #if i == 0:
        #raise NameError("Program couldn't import any geometric chords, check input path and format of geometric data\n"+geometry_path+'/detector_'+str(i)+'_z.txt')
    ##=====================  Generate the geometric matrix ======================

    nl = size(xchords,1)


    return xchords, ychords, zchords, nl



