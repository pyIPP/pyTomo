import ctypes as ct
import numpy as np
import sys,os
#from numpy.ctypeslib import ndpointer





#libdsc.f.argtypes = [ndpointer(numpy.uint8, flags="C_CONTIGUOUS"),
                       #ct.c_size_t]



#libdsc.f.argtypes = [ct.c_long,ct.c_long,ct.c_long,ct.c_long,ct.c_long,ct.c_long,ct.c_double,ct.c_double,ct.c_double,ndpointer(numpy.double, flags="C_CONTIGUOUS"), 

#c_nmom=ct.c_long(mom_order)
#c_ok=ct.c_long(0)
#c_nphi=ct.c_long(nphi)
#c_niter=ct.c_long(4000)
#c_nstep=ct.c_long(100)
#c_nfp=ct.c_long(1)
#c_ftol=ct.c_double(1.E-9)
#c_pexp=ct.c_double(4.E0)
#c_qexp=ct.c_double(4.E0)
#nphi1=nphi+1
#nphi3=nphi+3
#c_rmnaxis=(ct.c_float * nphi1)()
#c_zmnaxis=(ct.c_float * nphi1)()
#c_rbc=(ct.c_double * mom_order * nphi3)()
#c_zbc=(ct.c_double * mom_order * nphi3)()
#c_rbs=(ct.c_double * mom_order * nphi3)()
#c_zbs=(ct.c_double * mom_order * nphi3)()

#nrz=len(r_surf)
#c_nrz=ct.c_long(nrz)

#c_rin=(ct.c_double * nrz)()
#c_zin=(ct.c_double * nrz)()
#c_rin[:] = r_surf
#c_zin[:] = z_surf

#c_rin=(ct.c_double * nrz)()
#c_zin=(ct.c_double * nrz)()
#_nmom=ct.byref(c_nmom)
#_nrz =ct.byref(c_nrz)
#_nphi=ct.byref(c_nphi)
#_niter=ct.byref(c_niter)
#_nstep=ct.byref(c_nstep)
#_nfp  =ct.byref(c_nfp)
#_ftol =ct.byref(c_ftol)
#_pexp =ct.byref(c_pexp)
#_qexp =ct.byref(c_qexp)
#_rin  =ct.byref(c_rin)
#_zin  =ct.byref(c_zin)
#_ok   =ct.byref(c_ok)
#_rbc  =ct.byref(c_rbc)
#_rbs  =ct.byref(c_rbs)
#_zbc  =ct.byref(c_zbc)
#_zbs  =ct.byref(c_zbs)
#_rmnaxis=ct.byref(c_rmnaxis)
#_zmnaxis=ct.byref(c_zmnaxis)

        #libdsc.curve3d_(_nmom,_nrz,_nphi,_niter,_nstep,_nfp,_ftol,_pexp,_qexp, \
                    #_rin,_zin,_ok,_rbc,_zbs,_rbs,_zbc,_rmnaxis,_zmnaxis)





import os
import sys

class suppress_stdout_stderr(object):
    def __enter__(self):
        self.outnull_file = open(os.devnull, 'w')
        self.errnull_file = open(os.devnull, 'w')

        self.old_stdout_fileno_undup    = sys.stdout.fileno()
        self.old_stderr_fileno_undup    = sys.stderr.fileno()

        self.old_stdout_fileno = os.dup ( sys.stdout.fileno() )
        self.old_stderr_fileno = os.dup ( sys.stderr.fileno() )

        self.old_stdout = sys.stdout
        self.old_stderr = sys.stderr

        os.dup2 ( self.outnull_file.fileno(), self.old_stdout_fileno_undup )
        os.dup2 ( self.errnull_file.fileno(), self.old_stderr_fileno_undup )

        sys.stdout = self.outnull_file        
        sys.stderr = self.errnull_file
        return self

    def __exit__(self, *_):        
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr

        os.dup2 ( self.old_stdout_fileno, self.old_stdout_fileno_undup )
        os.dup2 ( self.old_stderr_fileno, self.old_stderr_fileno_undup )

        os.close ( self.old_stdout_fileno )
        os.close ( self.old_stderr_fileno )

        self.outnull_file.close()
        self.errnull_file.close()


class Silence:
    """Context manager which uses low-level file descriptors to suppress
    output to stdout/stderr, optionally redirecting to the named file(s).
    
    >>> import sys, numpy.f2py
    >>> # build a test fortran extension module with F2PY
    ...
    >>> with open('hellofortran.f', 'w') as f:
    ...     f.write('''\
    ...       integer function foo (n)
    ...           integer n
    ...           print *, "Hello from Fortran!"
    ...           print *, "n = ", n
    ...           foo = n
    ...       end
    ...       ''')
    ...
    >>> sys.argv = ['f2py', '-c', '-m', 'hellofortran', 'hellofortran.f']
    >>> with Silence():
    ...     # assuming this succeeds, since output is suppressed
    ...     numpy.f2py.main()
    ...
    >>> import hellofortran
    >>> foo = hellofortran.foo(1)
     Hello from Fortran!
     n =  1
    >>> print "Before silence"
    Before silence
    >>> with Silence(stdout='output.txt', mode='w'):
    ...     print "Hello from Python!"
    ...     bar = hellofortran.foo(2)
    ...     with Silence():
    ...         print "This will fall on deaf ears"
    ...         baz = hellofortran.foo(3)
    ...     print "Goodbye from Python!"
    ...
    ...
    >>> print "After silence"
    After silence
    >>> # ... do some other stuff ...
    ...
    >>> with Silence(stderr='output.txt', mode='a'):
    ...     # appending to existing file
    ...     print >> sys.stderr, "Hello from stderr"
    ...     print "Stdout redirected to os.devnull"
    ...
    ...
    >>> # check the redirected output
    ...
    >>> with open('output.txt', 'r') as f:
    ...     print "=== contents of 'output.txt' ==="
    ...     print f.read()
    ...     print "================================"
    ...
    === contents of 'output.txt' ===
    Hello from Python!
     Hello from Fortran!
     n =  2
    Goodbye from Python!
    Hello from stderr
    
    ================================
    >>> foo, bar, baz
    (1, 2, 3)
    >>>

    """
    def __init__(self, stdout=os.devnull, stderr=os.devnull, mode='w'):
        self.outfiles = stdout, stderr
        self.combine = (stdout == stderr)
        self.mode = mode
        
    def __enter__(self):
        import sys
        self.sys = sys
        # save previous stdout/stderr
        self.saved_streams = saved_streams = sys.__stdout__, sys.__stderr__
        self.fds = fds = [s.fileno() for s in saved_streams]
        self.saved_fds = list(map(os.dup, fds))
        # flush any pending output
        for s in saved_streams: s.flush()

        # open surrogate files
        if self.combine: 
            null_streams = [open(self.outfiles[0], self.mode, 0)] * 2
            if self.outfiles[0] != os.devnull:
                # disable buffering so output is merged immediately
                sys.stdout, sys.stderr = list(map(os.fdopen, fds, ['w']*2, [0]*2))
        else: null_streams = [open(f, self.mode, 0) for f in self.outfiles]
        self.null_fds = null_fds = [s.fileno() for s in null_streams]
        self.null_streams = null_streams
        
        # overwrite file objects and low-level file descriptors
        list(map(os.dup2, null_fds, fds))

    def __exit__(self, *args):
        sys = self.sys
        # flush any pending output
        for s in self.saved_streams: s.flush()
        # restore original streams and file descriptors
        list(map(os.dup2, self.saved_fds, self.fds))
        sys.stdout, sys.stderr = self.saved_streams
        # clean up
        for s in self.null_streams: s.close()
        for fd in self.saved_fds: os.close(fd)
        return False

class DESCUR:

  def descur_fit(self,r_surf,z_surf,mom_order):
    #libd = '/afs/ipp/home/t/transp/AUG_equilibrium/descur/@sys/lib/descur_idl.so'
    parent_path = os.path.dirname(os.path.realpath(__file__))
    #print parent_path
    libdsc = ct.cdll.LoadLibrary(parent_path+'/descur_idl.so')

    momtype=4
    nphi=1
    c_nmom=ct.c_long(mom_order)
    c_ok=ct.c_long(0)
    c_nphi=ct.c_long(nphi)
    c_niter=ct.c_long(4000)
    c_nstep=ct.c_long(100)
    c_nfp=ct.c_long(1)
    c_ftol=ct.c_double(1.E-9)
    c_pexp=ct.c_double(4.E0)
    c_qexp=ct.c_double(4.E0)
    nphi1=nphi+1
    nphi3=nphi+3
    c_rmnaxis=(ct.c_float * nphi1)()
    c_zmnaxis=(ct.c_float * nphi1)()
    c_rbc=(ct.c_double * mom_order * nphi3)()
    c_zbc=(ct.c_double * mom_order * nphi3)()
    c_rbs=(ct.c_double * mom_order * nphi3)()
    c_zbs=(ct.c_double * mom_order * nphi3)()

    moments=np.zeros((mom_order,momtype))
    nrz=len(r_surf)
    c_nrz=ct.c_long(nrz)
    if (nrz <= 0):
      return moments
  
    c_rin=(ct.c_double * nrz)()
    c_zin=(ct.c_double * nrz)()
    c_rin[:] = r_surf
    c_zin[:] = z_surf

    _nmom=ct.byref(c_nmom)
    _nrz =ct.byref(c_nrz)
    _nphi=ct.byref(c_nphi)
    _niter=ct.byref(c_niter)
    _nstep=ct.byref(c_nstep)
    _nfp  =ct.byref(c_nfp)
    _ftol =ct.byref(c_ftol)
    _pexp =ct.byref(c_pexp)
    _qexp =ct.byref(c_qexp)
    _rin  =ct.byref(c_rin)
    _zin  =ct.byref(c_zin)
    _ok   =ct.byref(c_ok)
    _rbc  =ct.byref(c_rbc)
    _rbs  =ct.byref(c_rbs)
    _zbc  =ct.byref(c_zbc)
    _zbs  =ct.byref(c_zbs)
    _rmnaxis=ct.byref(c_rmnaxis)
    _zmnaxis=ct.byref(c_zmnaxis)
    #if np.amin(z_surf) < -1:
      #import IPython
      #IPython.embed()
      #raise Exception('Unacceptable field lines')
    #else:

    with suppress_stdout_stderr():
        libdsc.curve3d_(_nmom,_nrz,_nphi,_niter,_nstep,_nfp,_ftol,_pexp,_qexp, \
                    _rin,_zin,_ok,_rbc,_zbs,_rbs,_zbc,_rmnaxis,_zmnaxis)

    moments[:,0]=c_rbc[1]
    moments[:,1]=c_rbs[1]
    moments[:,2]=c_zbc[1]
    moments[:,3]=c_zbs[1]
    return moments



  def descur_fit_fast(self,r_surf,z_surf,mom_order):
    #LSQR decoposition of the R,Z flux surface coordinates to Fourier moments
    #r_surf - (nt,nrho,nangle) or (nrho,nangle) or (nangle)
    #higher mom_order is neccessary! 
    
    nangle = np.size(r_surf,-1)
    
    rb = np.fft.rfft(r_surf)[...,:mom_order]/nangle
    zb = np.fft.rfft(z_surf)[...,:mom_order]/nangle
    
    moments=np.zeros(r_surf.shape[:-1]+(mom_order,4))

    moments[..., :mom_order,0]=np.real(rb)
    moments[..., :mom_order,1]=np.imag(rb)
    moments[..., :mom_order,2]=np.real(zb)
    moments[..., :mom_order,3]=np.imag(zb)
    moments[...,1:mom_order,:]*=2

    return moments


  def mom2rz(self,rcos,rsin,zcos,zsin,nthe=101,theta=None):
    #composition of the flux surfaces from the Fourier moments
    #rcos - jt,jrho,jmom  (or just jrho,jmom, or jmom)
    if theta==None:
        theta = -np.linspace(0,2*np.pi,nthe,endpoint=False) 
    #theta = np.linspace(-np.pi,np.pi,nthe)

    nmom = np.size(rcos,-1)
    angle = np.outer(np.arange(nmom),theta )
    cos = np.cos(angle)
    sin = np.sin(angle)
    r_plot = np.tensordot(rcos,cos,axes=([-1,0]))
    r_plot+= np.tensordot(rsin,sin,axes=([-1,0]))
    z_plot = np.tensordot(zcos,cos,axes=([-1,0]))
    z_plot+= np.tensordot(zsin,sin,axes=([-1,0]))
    #einsum 
    
    return r_plot,z_plot

def fourier2grid(rcos,rsin,zcos,zsin,n_theta ):
    
    n_fourier = rcos.shape[1]
    theta = np.linspace(-np.pi,np.pi,n_theta)
    R = np.tensordot(rcos[:,1:],np.cos(np.arange(1,n_fourier)*theta[:,None]),[1,1])
    R+= np.tensordot(rsin      ,np.sin(np.arange(1,n_fourier)*theta[:,None]),[1,1])
    R+= rcos[:,0][:,None]

    Z = np.tensordot(zcos[:,1:],np.cos(np.arange(1,n_fourier)*theta[:,None]),[1,1])
    Z+= np.tensordot(zsin      ,np.sin(np.arange(1,n_fourier)*theta[:,None]),[1,1])
    Z+= zcos[:,0][:,None]
    
    return R,Z


def main():
 
  import sys
  sys.path.append('/afs/ipp/home/t/todstrci/TRANSP')
  import matplotlib.pylab as plt
  from scipy.interpolate import splprep, splev
  import kk


  import eqi_map,dd
  n_fourier = 4

  n_rho = 40
  t = 2
  kk = eqi_map.eqi_map()
  dd = dd.shotfile()
  contours = kk.kkeqpsp(29022,t,np.linspace(0.01,0.99,n_rho),diag='EQH',exp='AUGD',
                        ed=0,rho_lbl='rho_tor')

  s=0.0 # smoothness parameter
  k=2 # spline order
  nest=-1 # estimate of number of knots needed (-1 = maximal)
  n_the_spl=201
  D = DESCUR()
  rcos,rsin,zcos,zsin = np.zeros((n_rho,n_fourier)),np.zeros((n_rho,n_fourier)),np.zeros((n_rho,n_fourier)),np.zeros((n_rho,n_fourier))
  for i,cr in enumerate(contours[0]):
    if np.size(cr)>1:
      plt.plot(cr[:,0],cr[:,1])
      
      tckp,u = splprep([cr[:,0],cr[:,1]],s=s,k=k,nest=-1)
      rspl,zspl = splev(np.linspace(0,1,n_the_spl),tckp)
      
      moments = D.descur_fit(rspl,zspl,n_fourier)
      rcos[i,:] = moments[:,0]
      rsin[i,:] = moments[:,1]
      zcos[i,:] = moments[:,2]
      zsin[i,:] = moments[:,3]
    else:
      rcos[i,:] = np.nan
      rsin[i,:] = np.nan
      zcos[i,:] = np.nan
      zsin[i,:] = np.nan
  
  #import IPython
  #IPython.embed()
  R,Z = D.mom2rz(rcos,rsin,zcos,zsin,1000)
  plt.plot(R.T,Z.T)
  plt.show()
  


if __name__ == "__main__":
  main()
