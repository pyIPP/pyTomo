

from numpy import *
from matplotlib.pyplot import *
from scipy.io import loadmat, savemat
from scipy import sparse
import config
from scipy.special import erfc
import os 

def Hamming(x, hollowness = 0.9 ,edge=0.4):
    hollowness = (1-hollowness)/2
    #print x.shape
    x = copy(x)
    x[x>=edge+hollowness]= 0
    x = x - edge+hollowness
    x = (x +abs(x)) /2
    x /= amax(x)
    P =  cos(2*pi*x)
  
    return P


def Hollow(x, hollowness = 0, edge = 1):
    """
    Artificial emissivity profile
    """

    P = Gaussian(x, 0, edge )
    P[P < 0] = 0
    P*= P* 1/(1+1/exp((-.2+x)*hollowness*150))

    return P

def Gaussian(x, hollowness = 0, edge = 1):
    """
    Artificial emissivity profile
    """
    sigma =.5
    P = exp(-((x)/sigma)**2)*tanh((1-x)*20)
    P[x>1] = 0

    return P
    
def Emiss(x, hollowness = 0, edge = 1):
    """
    Artificial emissivity profile
    """

    P = exp(-((x)/0.7)**2)

    return P

def Peaked(x, hollowness = 0, edge = 1):
    """
    Artificial emissivity profile  -Wolfram
    """
    peaking = 4
    P = (x**2+peaking**-2)**(-2)*tanh((1-x)*10)
    P[x>1] = 0

    return P



def Flat(x, hollowness = 0, edge = 1):
    """
    Flat profile
    """

    return abs(tanh((1-x)*10))+(tanh((1-x)*10))


def shepp_logan(R,Z):
    
    


    #
    #   This head phantom is the same as the Shepp-Logan except 
    #   the intensities are changed to yield higher contrast in
    #   the image.  Taken from Toft, 199-200.
    #      
    #         A    a     b    x0    y0    phi
    #        ---------------------------------
    ellipse = array([[  1 ,  .69 ,  .92 ,   0   ,  0  ,   0 ], 
                    [-.8  ,.6624 ,.8740 ,  0 , -.0184  , 0 ],
                    [-.2 , .1100 ,.3100  ,.22  ,  0  ,  -18],
                    [-.2  ,.1600 ,.4100 ,-.22  ,  0 ,    18],
                    [.1 , .2100 ,.2500 ,  0   , .35 ,   0 ],
                    [.1 , .0460 ,.0460  , 0  ,  .1   ,  0 ],
                    [.1 , .0460 ,.0460  , 0  , -.1   ,  0 ],
                    [.1 , .0460 ,.0230 ,-.08 , -.605 ,  0 ],
                    [.1  ,.0230 ,.0230 ,  0  , -.606 ,  0 ],
                    [.1  ,.0230 ,.0460  ,.06 ,-.605  ,  0 ]])
    
    ellipse = array([[  1  , .69 ,  .92 ,   0  ,   0  ,   0   ],
                    [-.98 ,.6624, .8740 ,  0 , -.0184 ,  0],
                    [-.02 ,.1100, .3100 , .22 ,   0 ,   -18],
                    [-.02 ,.1600 ,.4100 ,-.22,    0 ,    18],
                    [.01, .2100 ,.2500 ,  0 ,   .35 ,   0],
                    [.01 ,.0460, .0460 ,  0  ,  .1 ,    0],
                    [.01, .0460 ,.0460 ,  0  , -.1 ,    0],
                    [.01, .0460, .0230 ,-.08,  -.605,   0 ],
                    [.01, .0230, .0230 ,  0 ,  -.606 ,  0],
                    [.01 ,.0230 ,.0460  ,.06 , -.605 ,  0   ]])
                

            
    p = zeros((len(R),len(Z)))

    xax = linspace(-1,1,len(R))
    yax = linspace(-1,1,len(Z))
    X,Y = meshgrid(xax,yax)

    for k in range(len(ellipse)):
        asq = ellipse[k,1]**2       # a^2
        bsq = ellipse[k,2]**2       # b^2
        phi = ellipse[k,5]*pi/180  # rotation angle in radians
        x0 = ellipse[k,3]          # x offset
        y0 = ellipse[k,4]          # y offset
        A = ellipse[k,0]           # Amplitude change for this ellipse
        x=X.T-x0                    # Center the ellipse
        y=Y.T-y0  
        cosp = cos(phi) 
        sinp = sin(phi)
        idx= (x*cosp + y*sinp)**2/asq + (y*cosp - x*sinp)**2/bsq <= 1
        p[idx] += A

    return p



def phantom_generator(tokamak, tvec_new, nx_new=100, ny_new=100, profile = 'Gaussian', 
                      hollowness = 0.9, edge = 0.8,scale=1):
    """
    Create artificial profile according to given `profile` function with given time vector

    :param array tvec_new: New time vector
    :param function profile: Function of expected profile
    :param float edge: norm psi of SOL
    :param int nx_new: arbitrary number of pixels (can be higher than actual resulution to study dicretization error!)
    :param int ny_new: arbitrary number of pixels in vertical direction

    
    """
    
    print('phantom_generator '+ profile)

    from annulus import get_rho_field_mat, get_bd_mat
    
    M = get_rho_field_mat(tokamak, mean(tvec_new), nx_new, ny_new )
    M = tile(M, (size(tvec_new),1,1))
    rhop,magx, magy = tokamak.mag_equilibrium(mean(tvec_new),return_mean=True)
    
    from scipy.interpolate import RectBivariateSpline
 
    n_mag = 100

    M = M/(amax(M))
    
    x = sort(r_[linspace(0,1,n_mag),linspace(1,10,n_mag)])


    emissivity = zeros((nx_new*ny_new, size(tvec_new)),dtype=single)
    BdMat = get_bd_mat(tokamak, time=mean(tvec_new))

    M = array(M,ndmin=3)
    M = reshape(M, (-1,nx_new*ny_new), order='F')

    ind = argmin(M,1)
    ncx,ncy  = unravel_index(ind, (nx_new,ny_new))


    #interpolate profile to the actual resolution
    #center of pixels
    xgrid = linspace(tokamak.xmin,tokamak.xmax,nx_new+1)
    xgridc = (xgrid[1:]+xgrid[:-1])/2
    ygrid = linspace(tokamak.ymin,tokamak.ymax,ny_new+1)
    ygridc = (ygrid[1:]+ygrid[:-1])/2
    
    xgrid_out = linspace(tokamak.xmin,tokamak.xmax,tokamak.nx+1)
    xgridc_out = (xgrid_out[1:]+xgrid_out[:-1])/2
    ygrid_out = linspace(tokamak.ymin,tokamak.ymax,tokamak.ny+1)
    ygridc_out = (ygrid_out[1:]+ygrid_out[:-1])/2
    
    rhop,magx, magy = tokamak.mag_equilibrium(tvec_new)


    x0 = xgrid[ncx]
    y0 = ygrid[ncy]
    tsteps = size(tvec_new)

    if profile == "dot":
        emissivity = zeros((ny_new,nx_new, tsteps),dtype=single)
        raws = int(sqrt(tsteps))
        cols = tsteps/raws
        for i in range(raws):
            for j in range(cols):
                    emissivity[(ny_new*i)/raws,(nx_new*j)/cols , i*cols+j] = 100
        for i in range(cols*raws, tsteps):
            emissivity[ny_new/2, nx_new/2, i] = 10
        emissivity+= sqrt(config.rgmin)
        emissivity = emissivity.reshape((nx_new * ny_new, tsteps), order="F")


    elif profile == 'Asymmetric_full':
    
        data = load('./asymmetry/%d_full_prof.npz'%tokamak.shot)
        
        from scipy.interpolate import LinearNDInterpolator
        
        points = c_[data['R'].ravel(),data['Z'].ravel() ]        
        rhop,magx_all, magy_all = tokamak.mag_equilibrium(tvec_new,return_mean=True,n_rho=40)
 
        X,Y = meshgrid(xgridc,ygridc)
        emissivity[:] = LinearNDInterpolator(points,data['n_z'].ravel() , fill_value=0)(c_[X.flatten(order='F'), Y.flatten(order='F')])[:,None]

        
        
    elif profile == "loc_asym":

        theta = arctan2(ygridc[:,None]-y0.mean(), xgridc[None,:]-x0.mean()).flatten(order='F')
        P = exp(-x**2*5)*1e5*erfc((x-1)/10)

        P[x>=1] = 0
        emissivity[:] = interp(M,x,P).T

        fwhm = linspace(0.02,0.5,len(tvec_new))
        emissivity *= (1+exp(-((M-0.4)/fwhm[:,None])**2)*cos(theta)*0.5).T
    
    
        
        

    elif profile == "Asymmetric":
        
        try:
            emiss = loadtxt('./asymmetry/emiss_profile_%d.txt'%tokamak.shot)
            tvec_emiss,emiss = emiss[:,0], emiss[:,1:]
            it = argmin(abs(tvec_emiss-mean(tvec_new)))
            P = interp(x, linspace(0,1,size(emiss,1)),emiss[it,:] )
            
        except:
            print('emissivity profile for this shot was not found')
  
            P = exp(-x**2*5)*1e5*erfc((x-1)/10)
      
        emissivity[:] = interp(M,x,P).T

        
  
        try:
            pp
            asym = load('./asymmetry/asymmetry_emiss_%d_.npz'%tokamak.shot)
            t = asym['tvec']
            it = argmin(abs(t-mean(tvec_new)))
            x = asym['rho_pol']
            c = asym['Ac'][it]
            s = asym['As'][it]

            cx,sx = r_[0,x[it],1], r_[0,x[it],1]
            c,s = r_[0, c,0], r_[0, s,0]

        except:
            cx,c = [0,.3,.4,0.5,1],[0,0,.2,0,0]    #asymmetry parameter
            cx,c = [0,.4,.47,0.53,.6,1],[0,0,.2,0.2,0,0]    #asymmetry parameter
            cx = [0.00,0.02,0.05,0.08,0.11,0.15,0.18,0.21,0.24,0.27,0.30,0.34,0.37,0.40,0.43,0.46,0.49,0.52,0.55,0.58,0.60,0.63,0.66,0.68,0.71,0.73,0.75,0.78,0.80,0.82,0.84,0.86,0.87,0.89,0.91,0.92,0.94,0.95,0.97,0.98,0.99]
            c = [0.00,-0.03,-0.11,-0.14,-0.14,-0.15,-0.18,-0.24,-0.29,-0.33,-0.31,-0.26,-0.16,-0.09,-0.06,-0.02,0.02,0.08,0.15,0.20,0.23,0.25,0.27,0.29,0.31,0.35,0.38,0.40,0.42,0.46,0.48,0.48,0.45,0.41,0.38,0.36,0.35,0.33,0.32,0.30,0.28]

            sx,s = [0,1,2],[0,0,0]
            cx,c = [0,.2,.4,0.5,.6, 1],[0,-.2,.25,.3, .2,.0]    #asymmetry parameter
  
        from scipy.interpolate import interp1d,UnivariateSpline

        
        
        theta = arctan2(ygridc[:,None]-y0.mean(), xgridc[None,:]-x0.mean()).flatten(order='F')
        s_cos = UnivariateSpline(cx, c,k=2,s=1e-3)
        s_sin = UnivariateSpline(sx, s,k=2,s=1e-3)
        asym = exp(s_cos(M.flatten()).reshape(M.shape)*cos(theta))
        asym *= 1+s_sin(M.flatten()).reshape(M.shape)*sin(theta)
        asym[asym<0] = 0

        savetxt(tokamak.tmp_folder+'/cos',s_cos(linspace(0,1,100)) )
        savetxt(tokamak.tmp_folder+'/sin',s_sin(linspace(0,1,100)) )

        emissivity*= asym.T
      
        
    elif profile == 'snake':
        
        theta0 = linspace(0,2*pi,size(tvec_new))
        theta = arctan2(ygridc[:,None]-y0.mean(), xgridc[None,:]-x0.mean()).flatten(order='F')
        
        theta_, rho_ = linspace(-pi,pi, 100), linspace(0,1,50)
        
        for it, t in enumerate(theta0):
            
            theta_prof = exp(-((theta_-t)%(2*pi)-pi)**2/0.1)
            rho_prof = exp(-(rho_-.3 )**2/0.02)
            prof = outer(theta_prof,rho_prof)
      
            Rint = RectBivariateSpline(rho_ , theta_,prof.T )
            
            emissivity[...,it] = Rint.ev(M[it].flatten(order='F'), theta.flatten(order='F'))
            
  
        P = Gaussian(x, 0, 0.99)
        P /= amax(P)*4
        
        emissivity+= interp(M,x,P).T
        


    elif profile == 'mode':
        
        rhop,magx_all, magy_all = tokamak.mag_equilibrium(tvec_new,return_mean=True,n_rho=100)
        rho = linspace(1./magx_all.shape[0],1,magx_all.shape[0])
        theta = tokamak.mag_theta_star(tvec_new.mean(),rho,magx_all,magy_all,rz_grid=True, extrapolate = 1.1, nx=nx_new,ny=ny_new)
 
        HWWNM = 0.03
        Mode_pos = .6
        Delta = 0.2
        m = 4
        n = 1
        
        theta0 = linspace(0,2*pi,size(tvec_new))
        for it, t in enumerate(theta0):
            theta_prof = cos(m*theta+t*n)
            rho_prof = 1/(1+((M[it]-Mode_pos)/HWWNM)**2)
            emissivity[...,it] = Delta*theta_prof.flatten(order='C')*rho_prof
  
        emissivity[ isnan(emissivity)]  = 0
        
        P = Gaussian(x, 0, 0.99)
        P /= amax(P)
        
        emissivity+= interp(M,x,P).T
        emissivity[emissivity<0] = 0
        
        

    elif profile == 'shepp_logan':
        emissivity = shepp_logan(xgridc, ygridc)
        
    elif profile == "ring":
        G = Gaussian(x)
        emissivity = zeros((ny_new*nx_new, tsteps),dtype=single)
        for i in range(tsteps):
            P = Hamming(x, hollowness = 0.9 ,edge=(i+10)/float(tsteps+15))
            P+= G
            emissivity[:,i] = interp(M[i],x,P)

    elif profile == 'radial_modes':
        
        omega = linspace(0,20, len(tvec_new)/2)
        M[:] = M.mean(0)[None,:]
        h = 0.0
        edge = 0.99
        for it,w in enumerate(omega):
            
            P =  Gaussian(x, h, edge)
            P /= amax(P)
            emissivity[:, it] = interp(M[ it,:],x,P*(1+sin(2*pi*w*x+it)/10.)).T
            emissivity[:,-it] = interp(M[-it,:],x,P*(1+cos(2*pi*w*x+it)/10.)).T
        savetxt(tokamak.tmp_folder +'/emiss0', c_[x,P].T)
        
    elif profile == 'poloidal_modes':
        theta = arctan2(ygrid[:,None]-y0.mean(),xgrid[None,:]-x0.mean()).flatten(order='F')

        omega = linspace(0,20, len(tvec_new)/2)
        M[:] = M.mean(0)[None,:]
        h = 0.0
        edge = 0.99
        P =  Gaussian(x, h, edge)
        P /= amax(P)
    
        savetxt(tokamak.tmp_folder +'/emiss0', c_[x,P].T)

        P = interp(M[0,:],x,P)
        for it,w in enumerate(omega):
            emissivity[:,  it] = P*(1+sin(w*theta+it)*M[it,:]*(1-M[it,:])**2)
            emissivity[:, -it] = P*(1+cos(w*theta+it)*M[it,:]*(1-M[it,:])**2)

  
    elif profile == 'moving_circle':
        
        start_r = 2.0,
        end_r = 1.1,
        start_z = 1.0,
        end_z = -0.8,
        radius = 0.1,
        max_amp = 10.0,
        
  
        rcenter = linspace(start_r,end_r,len(tvec_new))
        zcenter = linspace(start_z,end_z,len(tvec_new))

        index=hypot(xgridc[None,:,None]-rcenter[None,None],
                    ygridc[:,None,None]-zcenter[None,None])<radius
 
        emissivity[index.reshape(nx_new*ny_new,-1,order='F')] = max_amp
        
        
    elif profile == "const":
        emissivity[:]= 1
    else:
        try:
            profile_fun = eval(profile)
        except:
            print('no such phantom exist:', (profile,))
            raise
        
        
        h = 0.3
        edge = 0.99
        P = profile_fun(x, h, edge)
        P /= amax(P)
        import os
        print('tmp',  os.path.isdir(tokamak.tmp_folder ))
        import sys
        import os
        print(os.getcwd())
        sys.stdout.flush()
       # import IPython
        #IPython.embed()
        savetxt(tokamak.tmp_folder +'/emiss0', c_[x,P*scale].T)

        emissivity[:] = interp(M,x,P).T


    from geom_mat_setting import geom_mat_setting

    if nx_new == tokamak.nx and ny_new == tokamak.ny:
        Tmat = tokamak.Tmat
    else:
        Tmat, X, Y = geom_mat_setting(tokamak, nx_new, ny_new, tokamak.virt_chord)  


    # avoid dividing by zero in the first estimation 
    emissivity[emissivity < emissivity.max()*1e-6]  = emissivity.max()*1e-6
      
    BdMat = get_bd_mat(tokamak, time=mean(tvec_new),nx=nx_new, ny=ny_new)
    emissivity[BdMat] = 0
    data =  Tmat*emissivity
    
    rescale = scale/(nanmean(data.ravel())*5)
    emissivity*= rescale
    data*= rescale

    N = size(data,1)
    if isscalar(tokamak.sigma):
        error = tokamak.sigma*data+amax(abs(data))*tokamak.min_error
    else:
        error = data*tokamak.sigma[:,None]+amax(data)*tokamak.min_error

    power1,power2,power3 = [],[],[]
    if tokamak.nx == nx_new and tokamak.ny == ny_new and False:
        emissivity_out  = emissivity
    else:
        emissivity_out = zeros((tokamak.ny*tokamak.nx,emissivity.shape[1]),dtype=single)

        for it in range(emissivity.shape[1]):
            Rint = RectBivariateSpline(xgridc,ygridc,emissivity[:,it].reshape(ny_new, nx_new,order='F').T,kx=2,ky=2)
            emissivity_out[:,it] = Rint(xgridc_out,ygridc_out).T.flatten(order='F')

    BdMat = get_bd_mat(tokamak, time=mean(tvec_new),nx=tokamak.nx, ny=tokamak.ny)

    emissivity_out[BdMat] = 0
    
    
    Tmat_ , _, _ = geom_mat_setting(tokamak, tokamak.nx, tokamak.ny, tokamak.virt_chord)

    savez_compressed(tokamak.tmp_folder+'/Emiss0.npz',G = reshape(emissivity_out,(tokamak.ny,tokamak.nx,tsteps), order='F')
            ,Rc= xgridc_out,Zc= ygridc_out,tvec=tvec_new, 
            xgridc=xgridc,ygridc=ygridc,G0=emissivity.reshape(ny_new, nx_new,tsteps,order='F'))

    return single(data),single(error), tvec_new, single(emissivity_out)





