"""
Module defining the Zernike polynomials
"""

import types
from scipy import *
from numpy import  *
from math import factorial
from numpy.ma import masked_array
#from pyoptools.misc.Poly2D import *

def polar_array(Rmax=1.,DS=0.1, pr=1.):
    """
    Function that greates 2 square matrices one with rho, and the other with
    theta, to be able to calculate functions using polar coordinates.
    It is similar to the mgrid function.

    Arguments:

    Rmax
        Limit the pupil area -Rmax<=X<=Rmax -Rmax<=Y<=Rmax
    DS
        Step between pixels
    pr
        Pupil radius. Used to normalize the pupil.

    TODO: This function should be moved to a auxiliary functions module
    """

    X,Y= mgrid[-Rmax:Rmax+DS:DS,-Rmax:Rmax+DS:DS]/pr
    r = sqrt(X**2+Y**2)
    th= arccos(X/r).T
    th= where(th<2*pi,th,0)
    th= where(X<0,2*pi-th,th)
    return r,th



#plot(vrnm(50,0,linspace(0,1,1000))
     
     
def rnm(n,m, rho):
    """
    Make radial Zernike polynomial on coordinate grid **rho**.
    @param [in] m Radial Zernike index
    @param [in] n Azimuthal Zernike index
    @param [in] rho Radial coordinate grid
    @return Radial polynomial with identical shape as **rho**
    
    
    http://www.univie.ac.at/nuhag-php/janssen/data/p156.pdf
    """
    m = abs(m)
    if (mod(n-m, 2) == 1):
        return rho*0.0
    
    def Un(x,n):
        t = arccos(x)
        t[abs(t)<1e-6] = 1e-6
        return sin((n+1)*t)/sin(t)
    
    N = 2*n+1
    

    
    Rmn = zeros_like(rho)
    
    for k in range(N):
        Rmn +=  Un( rho*cos(2*pi/N*k),n)* cos(2*pi/N*m*k)
    Rmn/= N
    
    return Rmn



def rnm_(n,m, rho):
    """
    Make radial Zernike polynomial on coordinate grid **rho**.
    @param [in] m Radial Zernike index
    @param [in] n Azimuthal Zernike index
    @param [in] rho Radial coordinate grid
    @return Radial polynomial with identical shape as **rho**
    """
    
    m = abs(m)
    if (mod(n-m, 2) == 1):
        return rho*0.0
    poly = zeros(n+1)
    for k in range((n-m)/2+1):
        poly[2*k] = (-1)**k * factorial(n-k) / ( factorial(k) * v( (n+m)/2 - k ) * factorial( (n-m)/2 - k ) )

    wf = polyval(poly,rho)
 
    return wf








def zernike(n,m,rho,theta,N):
    """
    Returns the an array with the Zernike polynomial evaluated in the rho and
    theta

    Arguments:


    *n*
        n order of the Zernike polynomial

    *m*
        m order of the Zernike polynomial

    *rho*
        Matrix containing the radial coordinates.

    *theta*
        Matrix containing the angular coordinates.

    Note: For rho>1 the returned value is 0

    Note: Values for rho<0 are silently returned as rho=0
    """

    
    rho_ = linspace(0,1,100)
    theta_ = linspace(-pi,pi,40)
    

    NC=sqrt(2*(n+1))
    Rnm=rnm(n,m,rho_)*NC

        
    if m>0:
        Zmn_t = cos(m*theta_)
    elif m<0:
        Zmn_t = sin(m*theta_)
    else:
        Zmn_t = sqrt(0.5)*ones_like(theta_)

    #trick to speed up the computation
    Zmn = interp(rho,rho_,Rnm)*interp(theta,theta_,Zmn_t)

        
    return Zmn


def i2nm(i):
    """
    Return the n and m orders of the i'th zernike polynomial

    index      0  1  2  3  4  5  6  7  8  9 N 11 12 13 14
    n-order    0  1  1  2  2  2  3  3  3  3  4  4  4  4  4
    m-order    0 -1  1 -2  0  2 -3 -1  1  3 -4 -2  0  2  4
    """
    # Calculate the polynomial order
    #Order      0   1   2   3   4
    #initindex  0   1   3   6   N
    ia=array(i)
    n=(1+(sqrt(8*(ia)+1)-3)/2).astype(int)
    ni=n*(n+1)/2
    m=-n+2*(i-ni)
    m = int(m)
    n = int(n)
    return n, m



from scipy.special import jn, jn_zeros, jnp_zeros, jv

def bessel(n,m,rho,theta,N):
    #http://www.cl.cam.ac.uk/~jrh13/papers/bessel.pdf!!!!
    #large numerical errors!

    #faster version
    rho_ = linspace(0,1,100)
    theta_ = linspace(-pi,pi,40)

    if m>0:
        Bmn_r = jn(m,rho_*jn_zeros(m,n+1)[-1])
        Bmn_t = cos(m*theta_)
    elif m<0:
        Bmn_r = jn(-m,rho_*jn_zeros(-m,n+1)[-1])
        Bmn_t = sin(m*theta_)
    else:
        Bmn_r = jn(m,rho_*jn_zeros(m,n+1)[-1])
        Bmn_t = ones_like(theta_)


    Bmn = interp(rho,rho_,Bmn_r)*interp(theta,theta_,Bmn_t)

    

    return Bmn


def polynom(n,m,rho,theta,N):

        
    if n > 0:
        coeff = exp(-rho**2*n)
    else :
        coeff = ones_like(rho)
        
    if m > 0:
        coeff*= rho**abs(m)*cos(m*theta)
    elif m < 0:
        coeff*= rho**abs(m)*sin(m*theta)

    
    coeff[rho>0.96*min(1,rho.max())] *= 0#exp(-rho[rho>1]*20)

    return coeff





        
#from matplotlib.pyplot import *

if __name__ == "__main__":
    r,th = polar_array()
    b = bessel(3,-3,r,th)
    pcolor(b)
    show()
