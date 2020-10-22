# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import sys
import scipy
import scipy.linalg
import scipy.sparse
import warnings
#from matplotlib.pylab import *
from time import time

# <codecell>

"""
GSVD   Generalized Singular Value Decomposition.
   [U,V,X,C,S] = GSVD(A,B) returns unitary matrices U and V,
   a (usually) square matrix X, and nonnegative diagonal matrices
   C and S so that

       A = U*C*X'
       B = V*S*X'
       C'*C + S'*S = I 

   A and B must have the same number of columns, but may have
   different numbers of rows.  If A is m-by-p and B is n-by-p, then
   U is m-by-m, V is n-by-n and X is p-by-q where q = min(m+n,p).

   SIGMA = GSVD(A,B) returns the vector of generalized singular
   values, sqrt(diag(C'*C)./diag(S'*S)).

   The nonzero elements of S are always on its main diagonal.  If
   m >= p the nonzero elements of C are also on its main diagonal.
   But if m < p, the nonzero diagonal of C is diag(C,p-m).  This
   allows the diagonal elements to be ordered so that the generalized
   singular values are nondecreasing.
   
   GSVD(A,B,0), with three input arguments and either m or n >= p,
   produces the "economy-sized" decomposition where the resulting
   U and V have at most p columns, and C and S have at most p rows.
   The generalized singular values are diag(C)./diag(S).
   
   When I = eye(size(A)), the generalized singular values, gsvd(A,I),
   are equal to the ordinary singular values, svd(A), but they are
   sorted in the opposite order.  Their reciprocals are gsvd(I,A).

   In this formulation of the GSVD, no assumptions are made about the
   individual ranks of A or B.  The matrix X has full rank if and only
   if the matrix [A B] has full rank.  In fact, svd(X) and cond(X) are
   are equal to svd([A B]) and cond([A B]).  Other formulations, eg.
   G. Golub and C. Van Loan, "Matrix Computations", require that null(A)
   and null(B) do not overlap and replace X by inv(X) or inv(X').
   Note, however, that when null(A) and null(B) do overlap, the nonzero
   elements of C and S are not uniquely determined.

   Class support for inputs A,B:
      float: double, single

   See also SVD.

   Copyright 1984-2007 The MathWorks, Inc.
   $Revision: 1.9.4.4 $  $Date: 2007/08/03 21:26:22 $
"""


#function [U,V,Z,C,S] = csd(Q1,Q2)


def csd(Q1,Q2):
    #print 'csd'
    '''CSD Cosine-Sine Decomposition
    [U,V,Z,C,S] = csd(Q1,Q2)
    
    Given Q1 and Q2 such that Q1'*Q1 + Q2'*Q2 = I, the
    C-S Decomposition is a joint factorization of the form
        Q1 = U*C*Z' and Q2=V*S*Z'
    where U, V, and Z are orthogonal matrices and C and S
    are diagonal matrices (not necessarily square) satisfying
        C'*C + S'*S = I
    '''
    m,p = Q1.shape
    n,pb = Q2.shape
    if pb != p:
        print("gsvd Matrix Column Mismatch : Matrices must have the same number of columns")
        sys.exit()
        
            

    
    if m < n:
        V,U,Z,S,C = csd(Q2,Q1)
        j = slice(p-1,None,-1)
        
        C = C[:,j]
        S = S[:,j]
        Z = Z[:,j]
        m = min(m,p)
        i = slice(m-1,None,-1)

        C[:m,:] = C[i,:] 
        U[:,:m] = U[:,i]
        n = min([n,p])
        i = slice(n-1,None,-1)
        S[:n,:] = S[i,:]
        V[:,:n] = V[:,i]
        return U,V,Z,C,S
        # Henceforth, n <= m.
 
    U,C,Z = scipy.linalg.svd(Q1,check_finite=False,full_matrices=True)
 
    C = add_zeros(C,Q1)
    q = min(m,p)
    i = slice(0,q)
    j = slice(q-1,None,-1)
  
    C[i,i]=C[j,j]
    U[:,i] = U[:,j]
    Z[:,i] = Z[:,j]
    Z = np.fliplr(np.flipud(Z)).T
    
    S = np.dot(Q2,Z)
    if q == 0:
        k = 1
    elif m < p:
        k = n
    else:
        k = np.r_[0,np.where(np.diag(C) <= 1/np.sqrt(2))[0]+1].max()
   
    #print np.diag(C)
    V, R = scipy.linalg.qr(S[:,:k],check_finite=False,mode='full')

    S = np.dot(V.T,S)
    r = min(k,m)
    S[:,:r]=diagf(S[:,:r])
    if m == 1 and p > 1:  S[1,1] = 0

    if k < min(n,p):
        r = min(n,p)
        i = slice(k,n)
        j = slice(k,r)
 
        UT, ST, VT = scipy.linalg.svd(S[i,j],check_finite=False,full_matrices=True)


        ST = add_zeros(ST,np.zeros([n-k,r-k]))

        if k > 0:    S[:k,j] = 0
        S[i,j] = ST


        C[:,j] = np.dot(C[:,j],VT)
        V[:,i] = np.dot(V[:,i],UT)
        Z[:,j] = np.dot(Z[:,j],VT)
        i = slice(k,q)

        Q,R = scipy.linalg.qr(C[i,j],check_finite=False,mode='full')

        C[i,j] = diagf(R)
        U[:,i] = np.dot(U[:,i],Q)
        
    if m < p:
        # Diagonalize final block of S and permute blocks.
        q = min(scipy.sparse.lil_matrix(abs(diagk(C,0))>10*m*np.spacing(1)).getnnz(), 
            scipy.sparse.lil_matrix(abs(diagk(S,0))>10*n*np.spacing(1)).getnnz())
  
        # At this point, S(i,j) should have orthogonal columns and the
        # elements of S(:,q+1:p) outside of S(i,j) should be negligible.
        Q,R = scipy.linalg.qr(S[q:n,m:p],check_finite=False,mode='economic')
        S[:,q:p-1] = 0
        S[q:n,m:p] = diagf(R)
        V[:,q:n] = np.dot(V[:,q:n],Q)
        if n > 1:
            i = np.r_[q:q+p-m,:q,q+p-m:n]
        else:
            i = 1
        j = np.r_[m:p,:m]
        t = slice(0,len(j))
        C[:,t] = C[:,j]
        S = S[q:n,m:p]
        Z = Z[:,j]
        V = V[:,i]

    if n < p:
        # Final block of S is negligible.
        S[:,n:p] = 0
    

    # Make sure C and S are real and positive.
    U,C = diagp(U,C,max(0,p-m))
    C = C.real
    V,S = diagp(V,S,0)
    S = S.real

    return U, V, Z, C, S

# <codecell>

def add_zeros(C,Q):
    '''ADD_ZEROS returns the vector C padded with zeros to be the same size as matrix Q. 
    The values of C will be along the diagonal.
    
    USAGE: add_zeros(C,Q)
    '''
    
    assert C.shape > 1
    assert Q.shape > 1
    m,p=Q.shape
    n = C.size
    toto = np.zeros_like(Q,dtype=C.dtype)

    toto[:n,:n]=np.diag(C)

    return toto

# <codecell>

def diagk(X,k):
    '''DIAGK K-th matrix diagonal.
    DIAGK(X,k) is the k-th diagonal of X, even if X is a vector.'''
    if min(X.shape)> 1:
        D = np.diag(X,k)
    elif 0 <= k and 1+k <= X.shape[1]:
        D = X(1+k)
    elif k < 0 and 1-k <= X.shape[0]:
        D = X(1-k)
    else:
        D = []
    return D

# <codecell>

def diagf(X):
    ''' DIAGF Diagonal force.
    X = DIAGF(X) zeros all the elements off the main diagonal of X.
    '''
    D = np.triu(np.tril(X))
    return D

# <codecell>

def diagp(Y,X,k):
    ''' DIAGP Diagonal positive.
    Y,X = diagp(Y,X,k) scales the columns of Y and the rows of X by
    unimodular factors to make the k-th diagonal of X real and positive.
    '''
    
    
    D = diagk(X,k)
    #j = [item for item, a in enumerate(D) if a.real < 0 or a.imag != 0]
    ind = (D.real < 0) | (D.imag != 0)
    

    Y[:,ind] *= (D[ind].conjugate()/abs(D[ind]))[None,:]
    X[ind,:] *= (D[ind].conjugate()/abs(D[ind]))[:,None]

    return Y, X

# <codecell>

def trim_matrix_col(Q,p):
    '''TRIM_MATRIX_COL trims the output of the 
    Q matrix output of the qr function column-wise to the number of columns of 
    the input matrices to scipy.linalg.qr to 
    match the format of the Matlab qr() function and returns Q.
    
    USAGE trim_matrix(Q,p)
    where 
    Q is the Q output matrix of scipy.linalg.gr
    p is the number of columns of the input matrices to scipy.linalg.gr
    '''
    Q=Q[:,:p]
    return Q

# <codecell>

def trim_matrix_row(R,p):
    '''TRIM_MATRIX_ROW trims the output of the 
    R matrix output of the qr function row-wise to the number of rows of 
    the input matrices to scipy.linalg.qr to 
    match the format of the Matlab qr() function and returns R.
    
    USAGE trim_matrix(R,p)
    where 
    R is the R output matrix of scipy.linalg.gr
    p is the number of columns of the input matrices to scipy.linalg.gr
    '''
    R=R[:p,:]
    return R

# <codecell>

def gsvd(A,B,economic=True):
    
    m,p = A.shape
    n,pb = B.shape


    if pb != p:
        print("gsvd Matrix Column Mismatch : Matrices must have the same number of columns")
        sys.exit()
    QA = None
    QB = None
    

    t = time()

    if economic: 
        # Economy-sized.
        if m > p:
            QA, A = scipy.linalg.qr(A,check_finite=False,mode='economic')
            QA, A = diagp(QA,A,0)
            #print np.linalg.norm(A-dot(QA,A)), np.linalg.norm(dot(QA.T,QA) - np.eye(QA.shape[1])), 
            m = p
        if n > p:
            QB, B = scipy.linalg.qr(B,check_finite=False,mode='economic')
            QB, B = diagp(QB,B,0)
            n = p
            

    AB = np.r_[A,B]

    Q, R = scipy.linalg.qr(AB,overwrite_a=True,mode='economic',check_finite=False)


    U, V, Z, C, S = csd(Q[:m,:],Q[m:m+n,:])

    X = np.dot(R.T,Z)
    if not QA is None:  U = np.dot(QA,U)
    if not QB is None:  V = np.dot(QB,V)

    return U,V,X,C,S
    




def fast_svd(M, k):  # random projection SVD 
    import numpy.linalg as linalg
    from  scipy.linalg import qr,qr_multiply,svd


    transpose = False
    if 3*M.shape[1]<M.shape[0]:
        transpose = True
        M = M.T
 
    p = k+5

    M = ascontiguousarray(M)
    R = np.random.normal(size=(p,M.shape[1])).astype(M.dtype)
    Y = dot(R ,M.T).T

    Q,r = qr(Y,mode='economic', overwrite_a=True, check_finite=False)
    
    M = asfortranarray(M)
    B  = np.dot(Q.T,asfortranarray(M))
    
    Uhat, s, v = svd(B,full_matrices=False,overwrite_a=True,check_finite=False)
    U = np.dot(Q, Uhat)
    
    if  transpose:
        U,s,v = v[:k,...].T, s[:k], U[...,:k].T
    else:
        U,s,v = U[...,:k], s[:k], v[:k,...] 
        
    return U,s,v



def svd_qr(A, tol = 1e-8, N = None): #rank revealing SVD
    #TODO otestovat vliv C nebo F poradi c matici A na rychlost vypoctu
    from  scipy.linalg import svd,qr

    if size(A,0) < size(A,1):
        transpose = True
        A = A.T
    else:
        transpose = False

    q,r = qr(A,overwrite_a=False, mode='economic', pivoting=False)  #BUG jde lépe, bez výpočetu Q!!
    u, s, vt = svd(r,full_matrices = False, overwrite_a = True )

  
    if N is not None:
        index = zeros(len(s), dtype=bool)
        index[:N] = True
    else:
        index = s/s[0] > tol

    
    if transpose:
        V = (dot(q, u)[:,index]).T
        U = vt[index,:].T
    else:
        U = dot(q, u)[:,index]
        V = vt[index,:]
        
    S = s[index]


    return U,S,V



def main():
    import numpy as np
    import scipy.sparse as sp
    npix = 1000
    k = 100
    A = np.random.rand(k,npix)
    #B = np.random.rand(npix,npix)
    B = sp.spdiags(np.ones(npix),0,npix, npix)
    B = B-sp.spdiags(np.ones(npix),1,npix, npix)
    B = B+sp.spdiags(np.ones(npix),-1,npix, npix)
    B = B.todense()
    
    U,V,X,C,S = gsvd(A,B,economic=True)
    
    
    X.shape
    U.shape
    V.shape
    print(np.linalg.norm(np.dot(np.dot(U,C),X.T)-A),np.linalg.norm(A))
    print(np.linalg.norm(np.dot(np.dot(V,S),X.T)-B),np.linalg.norm(B))

    print(np.allclose(np.dot(np.dot(U,C),X.T),A))
    print(np.allclose(np.dot(np.dot(V,S),X.T),B))
    print(np.allclose(np.dot(C.T,C)+np.dot(S.T,S),np.eye(B.shape[1])))


    C = np.diag(C[::-1,::-1])
    S = np.diag(S)[:-k-1:-1]
    
    #V se nepotřebuje!
    X = X[:,-k:]
    
    np.linalg.norm(np.dot(np.dot(U,np.diag(C)),X.T)-A),np.linalg.norm(A)
    np.linalg.norm(np.dot(np.dot(V[:,-k:],np.diag(S)),X.T)-B)

  
if __name__ == "__main__":
    main()
