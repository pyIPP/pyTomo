#!/usr/bin/env python""" SuiteSparseQR Python wrapper """

import os.path
import ctypes
from ctypes import c_double, c_size_t, byref, pointer, POINTER,c_long,cast
import numpy as np
from numpy.ctypeslib import ndpointer
from scipy.sparse import csc_matrix, isspmatrix_csc, SparseEfficiencyWarning
#from matplotlib.pylab import *
from scipy.sparse import rand

# Assume spqr_wrapper.so (or a link to it) is in the same directory as this file
spqrlib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__) + os.path.sep + "spqr_wrapper.so")

# Function prototypes for spqr_wrapper.so
# void qr_solve(double const *A_data, double const *A_row, double const *A_col, size_t A_nnz, size_t A_m, size_t A_n, double const *b_data, double *x_data)
spqrlib.qr_solve.argtypes = [
        ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # A_data
        ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'), # A_row
        ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'), # A_col
        c_size_t,  # A_nnz
        c_size_t,  # A_m
        c_size_t,  # A_n
        ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # b_data
        ndpointer(dtype=np.float64, ndim=1, flags=('C_CONTIGUOUS', 'WRITEABLE')), # x_data
        ]
spqrlib.qr_solve.restype = None

spqrlib.qr_solve_csr.argtypes = [
        ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # A_data
        ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'), # A_row
        ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'), # A_col
        c_size_t,  # A_nnz
        c_size_t,  # A_m
        c_size_t,  # A_n
        ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # b_data
        ndpointer(dtype=np.float64, ndim=1, flags=('C_CONTIGUOUS', 'WRITEABLE')), # x_data
        ]
spqrlib.qr_solve_csr.restype = None


spqrlib.qr_sparse_csc.argtypes = [
        ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # A_data
        ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'), # A_row
        ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'), # A_col
        c_size_t,  # A_nnz
        c_size_t,  # A_m
        c_size_t,  # A_n
        POINTER(POINTER(c_double)),#Qx
        POINTER(POINTER(c_long)),#Qi
        POINTER(POINTER(c_long)),#Qp
        POINTER(c_size_t),#Q_nnz
        POINTER(c_size_t),#Q_m
        POINTER(c_size_t),#Q_n
        POINTER(POINTER(c_double)),#Rx
        POINTER(POINTER(c_long)),#Ri
        POINTER(POINTER(c_long)),#Rp
        POINTER(c_size_t),#R_nnz
        POINTER(c_size_t),#R_m
        POINTER(c_size_t),#R_n
        POINTER(POINTER(c_long)),#E
        ]
spqrlib.qr_sparse_csc.restype = c_long


spqrlib.qr_sparse_csc_r.argtypes = [
        ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # A_data
        ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'), # A_row
        ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'), # A_col
        c_size_t,  # A_nnz
        c_size_t,  # A_m
        c_size_t,  # A_n
        POINTER(POINTER(c_double)),#Rx
        POINTER(POINTER(c_long)),#Ri
        POINTER(POINTER(c_long)),#Rp
        POINTER(c_size_t),#R_nnz
        POINTER(c_size_t),#R_m
        POINTER(c_size_t),#R_n
        POINTER(POINTER(c_long)),#E
        ]
spqrlib.qr_sparse_csc_r.restype = c_long





def qr_solve(A_data, A_row, A_col, A_nnz, A_m, A_n, b_data):
    """ Python wrapper to qr_solve """
    if len(A_data) != len(A_row) != len(A_col) != A_nnz:
        raise TypeError("A_data, A_row, A_col, A_nnz must agree")
    if len(b_data) != A_m:
        raise TypeError("b_data must be A_m long")

    x_data = np.empty(A_n, dtype=np.float64)
    spqrlib.qr_solve(
            np.require(A_data, np.float64, 'C'),
            np.require(A_row, np.int64, 'C'),
            np.require(A_col, np.int64, 'C'),
            A_nnz, A_m, A_n,
            np.require(b_data, np.float64, 'C'),
            np.require(x_data, np.float64, 'C')
            )
    return x_data

def qr_solve_csr(A_data, A_i, A_p, A_nnz, A_m, A_n, b_data):
    """ Python wrapper to qr_solve """
    if len(A_data) != len(A_i) != len(A_p) != A_nnz:
        raise TypeError("A_data, A_i, A_p, A_nnz must agree")
    if len(b_data) != A_m:
        raise TypeError("b_data must be A_m long")

    x_data = np.empty(A_n, dtype=np.float64)
    spqrlib.qr_solve_csr(
            np.require(A_data, np.float64, 'C'),
            np.require(A_i, np.int64, 'C'),
            np.require(A_p, np.int64, 'C'),
            A_nnz, A_m, A_n,
            np.require(b_data, np.float64, 'C'),
            np.require(x_data, np.float64, 'C')
            )
    return x_data



def qr_sparse(A):
    """ Python wrapper to qr """
    
    #BUG it will cause a memory leak!!

    if not isspmatrix_csc(A):
        Warning( SparseEfficiencyWarning())
    
    A = csc_matrix(A)

    A_data = A.data
    A_i  = A.indices
    A_p = A.indptr
    A_nnz = A.nnz
    A_m = A.shape[0]
    A_n = A.shape[1]

    Qx = POINTER(c_double)()#Qx
    Qi =  POINTER(c_long)()#Qi
    Qp =  POINTER(c_long)()#Qp
    Q_nnz = c_size_t()
    Q_m = c_size_t()
    Q_n = c_size_t()
    Rx = POINTER(c_double)()#Rx
    Ri = POINTER(c_long)()#Ri
    Rp = POINTER(c_long)()#Rp
    R_nnz = c_size_t()

    R_m = (c_size_t)()#R_m
    R_n = (c_size_t)()#R_n
    E =  POINTER(c_long)()#E


    err = spqrlib.qr_sparse_csc(
            np.require(A_data, np.float64, 'C'),
            np.require(A_i, np.int64, 'C'),
            np.require(A_p, np.int64, 'C'),
            A_nnz, A_m, A_n,byref(Qx),byref(Qi),
            byref(Qp),byref(Q_nnz),byref(Q_m),
            byref(Q_n),byref(Rx),byref(Ri),
            byref(Rp),byref(R_nnz),byref(R_m),byref(R_n),byref(E) )

    if err:   raise Exception('Sparse QR has failured')

    Qx = np.asarray(Qx[:Q_nnz.value])
    Qi = np.asarray(Qi[:Q_nnz.value])
    Qp = np.asarray(Qp[:Q_n.value+1])

    Q = csc_matrix((Qx, Qi, Qp), shape=(Q_m.value,  Q_n.value))
    
    Rx = np.asarray(Rx[:R_nnz.value])
    Ri = np.asarray(Ri[:R_nnz.value])
    Rp = np.asarray(Rp[:R_n.value+1])
    R = csc_matrix((Rx, Ri, Rp), shape=(R_m.value,  R_n.value))
    #E = np.arange(A_n)
    #print E
    E = np.asarray(E[:A_n])
    #print E
    ind = np.argsort(E)
    #print E[ind]
    #import IPython
    #IPython.embed()
    
    #Q = Q[:,E]
    #print np.linalg.norm((A[:,E]-(Q*R)).data)

    #print np.linalg.norm((A-(Q*R)[:,ind]).data)
    #print np.linalg.norm((A-(Q[:,:]*R[:,ind])).data)
    #imshow(R[:,:].todense())
    #show()
    #exit()
    
    return Q,R,E#R[:,ind]



def qr_sparse_r(A):
    """ Python wrapper to qr """
    
    #BUG it will cause a memory leak!!

    if not isspmatrix_csc(A):
        Warning( SparseEfficiencyWarning())
    
    A = csc_matrix(A)

    A_data = A.data
    A_i  = A.indices
    A_p = A.indptr
    A_nnz = A.nnz
    A_m = A.shape[0]
    A_n = A.shape[1]

   
    Rx = POINTER(c_double)()#Rx
    Ri = POINTER(c_long)()#Ri
    Rp = POINTER(c_long)()#Rp
    R_nnz = c_size_t()

    R_m = (c_size_t)()#R_m
    R_n = (c_size_t)()#R_n
    E =  POINTER(c_long)()#E
    
    #print bool(E)
    err = spqrlib.qr_sparse_csc_r(
            np.require(A_data, np.float64, 'C'),
            np.require(A_i, np.int64, 'C'),
            np.require(A_p, np.int64, 'C'),
            A_nnz, A_m, A_n,byref(Rx),byref(Ri),
            byref(Rp),byref(R_nnz),byref(R_m),byref(R_n),byref(E) )

    if err:   raise Exception('Sparse QR has failured')
    #print bool(E)

    #Qx = np.asarray(Qx[:Q_nnz.value])
    #Qi = np.asarray(Qi[:Q_nnz.value])
    #Qp = np.asarray(Qp[:Q_n.value+1])

    #Q = csc_matrix((Qx, Qi, Qp), shape=(Q_m.value,  Q_n.value))
    
    Rx = np.asarray(Rx[:R_nnz.value])
    Ri = np.asarray(Ri[:R_nnz.value])
    Rp = np.asarray(Rp[:R_n.value+1])
    R = csc_matrix((Rx, Ri, Rp), shape=(R_m.value,  R_n.value))
    #import IPython
    #IPython.embed()
    E =  np.asarray(E[:A_n]) if bool(E) else arange(A_n)
    
    assert  all(np.sort(E)== np.arange(A_n)), 'SPQR arror'
    #if any(sort(E)!= arange(A_n)):  #strange bug!!
        #print 'strange bug'
        #E = arange(A_n)
        
    #if bool(E): 
        #E = np.asarray(E[:A_n])
    #else:
    #E = arange(A_n)

    #E = np.asarray(E[:A_n])
    #print E
    #ind = np.argsort(E)
    #print E[ind]
   
    
    #Q = 
    #print np.linalg.norm((A[:,E]-(Q*R)).data)

    #print np.linalg.norm((A-(Q*R)[:,ind]).data)
    #print np.linalg.norm((A-(Q[:,:]*R[:,ind])).data)
    #imshow(R[:,:].todense())
    #show()
    #exit()
    
    return R,E



#qr_solve_csr(double const *A_data, long const *A_i_data, long const *A_p_data, 
                  #size_t A_nnz, size_t A_m, size_t A_n, double const *b_data, double *x_data) 

def main():
    print("Testing qr_solve")
    A_data = np.array([4, 9, 25], dtype=np.float64)
    A_row = np.array([0, 1, 2])
    A_col = np.array([0, 2, 1])
    b_data = np.array([1, 1, 1], dtype=np.float64)
    from scipy.sparse import csc_matrix, coo_matrix
    from scipy.sparse.linalg import spsolve
        
    A = np.random.randn(100,200)
    A = np.dot(A,A.T)

    #A = coo_matrix(A)

    #x_data = qr_solve(A.data, A.row, A.col,A.nnz, 3, 3, b_data)
    #print(x_data)
    
    ##import IPython
    ##IPython.embed()

    ##print len(A.data),len(A.indices),len(A.indptr  )
    ##print A.indices, A.indptr, 
    #print 'qr_solve_csr'
    #x_data = qr_solve_csr(A.data, A.indices,A.indptr,A.nnz, 3, 3, b_data)

    ##from scipy.linalq import qr
    A = rand(100, 100, density=0.01, format='csc')
    
    #print(x_data)
    #print spsolve(A,b_data)
    #print A.todense()
    #E = np.ones(A.shape[1],dtype=int)
    from time import time
    t = time()
    
    print( 'qr_sparse')
    Q,R,E  = qr_sparse(A)
    print( E)

    print( 'qr_sparse_r')
    R_,E  = qr_sparse_r(A)
    print( E)
    print( time()-t)

    #exit()
    #print time()-t
    from scipy.linalg import qr
    import matplotlib.pylab as plt
    t = time()
    q,r = qr(A.todense(), mode='economic')
    #print np.linalg.norm(A-np.dot(q,r))
    
    plt.title('r-r')
    plt.subplot(221)
    plt.imshow(r)
    plt.colorbar()

    plt.subplot(222)

    plt.imshow(R.todense())
    plt.colorbar()
    
    
    plt.subplot(223)

    plt.imshow(R_.todense())
    plt.colorbar()
    
    
    plt.subplot(224)

    plt.title('difference')
    plt.imshow((R_-R).todense())
    plt.colorbar()
    plt.show()
    
    #figure()
    #imshow(q)

    
    #show()
    #print Q.shape, R.shape, q.shape, r.shape
    


    
    
    
    #print Qx[0],Qi[0],Qp[0]
         #ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # A_data
        #ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'), # A_row
        #ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS'), # A_col
        #c_size_t,  # A_nnz
        #c_size_t,  # A_m
        #c_size_t,  # A_n
        #POINTER(c_double),#Qx
        #POINTER(c_long),#Qi
        #POINTER(c_long),#Qp
        #POINTER(c_size_t),#Q_nnz
        #POINTER(c_size_t),#Q_m
        #POINTER(c_size_t),#Q_n
        #POINTER(c_double),#Rx
        #POINTER(c_long),#Ri
        #POINTER(c_long),#Rp
        #POINTER(c_size_t),#R_nnz
        #POINTER(c_size_t),#R_m
        #POINTER(c_size_t),#R_n
        #]

if __name__ == "__main__":
    main() 
    
    
#http://www.sagemath.org/doc/numerical_sage/ctypes_examples.html
#http://osdir.com/ml/python.ctypes/2007-02/msg00006.html
#http://stackoverflow.com/questions/16023614/ctypes-c-function-returns-array-of-unknown-size
