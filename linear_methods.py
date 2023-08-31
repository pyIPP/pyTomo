#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from matplotlib.pylab import *

from numpy import *
from matplotlib.pyplot import *
import time
from scipy import sparse
from shared_modules import prepare_inputs, create_derivation_matrix, make_postprocessing, debug
import scipy.ndimage as ndi
from scipy.optimize import nnls, _nnls
from annulus import get_bd_mat
from sksparse.cholmod import  cholesky,CholmodWarning,CholmodError
#from scikits.sparse.cholmod import cholesky,CholmodWarning,CholmodError
from  scipy.linalg import qr, qr_multiply, solve_triangular, cholesky_banded,eigh,svd, inv,lstsq,pinv
from scipy.sparse.linalg import svds, eigsh
from gsvd import gsvd
import warnings


#datatype used to represent the solution
dtype = 'float32'

def linear_methods(Tok,input_parameters, data, error, tvec, Tmat,dets, normData,G0, danis, boundary, 
                    regularization, solver, reconstruct, ifishermax, postprocessing = False, progress=0):
    """
    Pixel tomography based on SVD with anisotropic smoothing.
    This algorithm identify boundaries and create virtual sensors on boundary,
    prepare matrices of derivation, rescale and renormalize all its (Tmat, data, error)
    
    All methods are optimized for npix >> ndets! if is it not true, the algorithms
    will not be efficients, and it is necessary to modify them. 



    :param class tok: ``Class`` of  tokamak with all important data.
    :var double boundary:    If non-zero allow smoothing with boundary => 
    The reconstruction will be zero on the boundary,  Set the pressure forcing
    reconstruction to be zero on boundary. Precision of the virtual senzors (boundary) is min(sigma)/X.
    :param int regularization: Apriory smoothing used for minimalization
    :param bool reconstruct: Perform tomographic reconstruction
    :param bool postprocessing: Perform postprocessing --  entropy, smoothness
    :var double danis:    Ratio of anisotropic matrices, B = n*Bperp+1/n*Bpar  => 
    The smaller number the more will reconstruction follow the field For MFI is recommended 0.2, for MDIFF < 0.1
    :param object progress: link to QObject that move with the progressbar, do not work with multiprocessing
    :var spmatrix diam,diar,dial,diao,diau: Sparse matrices of derivation (elementar directions)


    :var spmatrix Err:  Square root of variation matrix
    :var spmatrix fs:   Normalised data
    :var spmatrix Ts:   Normalised geometry matrix

    :var spmatrix H:    Final version of derivation matrix
    :var list Bmat:     List of smoothing matrices


    .. sectionauthor:: Tomáš Odstrčil <tomasodstrcil@gmail.com>
    .. sectionauthor:: Michal Odstrčil <michalodstrcil@gmail.com>

    
    """
    #print G0
    
    #exit()
    if any(isnan(error)):
        raise Exception('NaNs in the  errorbars!')
    if any(~isfinite(data)) :
        raise Exception('NaNs in the data !')
    if all(~isfinite(error)) :
        raise Exception('all errors are infinity')
    
    global nx, ny, tokamak
    nx = Tok.nx
    ny = Tok.ny
    tokamak  = Tok
    solver  = int(solver)


    tsteps = len(tvec)
    #rgmin = 1e-4  # very important constant, affects stability of SVD/QR

    t_all = time.time()

    Steps = ifishermax
    

    if input_parameters['main_run'] and input_parameters['rotation_tomography']:
        Steps = 2
    

    bound = False #boundary will be introduced to H matrix in order to make H regular and reduces rank of T

    Bmat,BdMat,bnd, Ts, fs,err,tvec, G  = prepare_inputs(Tok, data, error,dets, tvec, Tmat, normData, G0, danis,bound, regularization, reconstruct, postprocessing)
                        
    ndets = size(fs, 0)
   

    for step in range(Steps):  # minfisher iterations !!!
        if progress != 0:
            progress.iterateSubStep()
        t = time.time()
        
 
        debug('\n======== Linear method -  step === %i' % (step+1))

        H = create_derivation_matrix(G, Bmat,BdMat, regularization, danis, Tok )
        

        last_step = (step == Steps-1)

        H = H.tocsc()
    

        if Tok.transform_index == 0 and boundary> 0:
            H =  H + sparse.spdiags(BdMat*sinh(maximum(0,boundary)), 0, nx*ny, nx*ny,format='csc')
       
       
       
       
        solvers = {2:PresolveSVD,3:PresolveSVD3,4:PresolveQR,5:PresolveGEV,6:PresolveGEV_singular}
        
   
        debug('--------Prepare - time  %g' % (time.time()-t))

        
       

        decomposition = solvers[solver](Ts,H, Tok.nx,Tok.ny,ndets,BdMat)
        debug('--------Decomposition - time  %g' % (time.time()-t))


        if reconstruct:
            debug('Linear method solving')
            g_fract_0 = None
            
            #input_parameters
            use_gcv = False
            lam_method = input_parameters['lambda_solver']
            positive_constrain = input_parameters['positive_constrain']

            try: g_fract_0 = float(input_parameters['lambda_solver'])
            except: pass

            error_scale = input_parameters['error_scale']
            rapid_solver = input_parameters['rapid_solver']
            
         
            estimate_sigma = input_parameters['estimate_sigma']
            G,retro, chi2,g,SigmaGsample = SolveLinearMetods( decomposition,tvec,
                        fs,ndets,dets, normData/Tok.normTmat,Ts,H,BdMat, error_scale,
                        last_step,lam_method,positive_constrain,
                        g_fract_0,rapid_solver, input_parameters,estimate_sigma,
                        input_parameters['lambda_up_lim'],input_parameters['lambda_low_lim'])
           
   
            wrong_dets = decomposition['wrong_dets']
            #BUG  not compatible wit QR!!
            U = decomposition['U']
            R = 1 if not 'R' in decomposition else decomposition['R']
            
            debug('Compute emissivity ')
            #Emis = FastLinEmisivity(decomposition,  fs, normData/Tok.normTmat, tsteps)

            debug('--------Solve - time  %g' % (time.time()-t))
            
            if debug:
                resid = linalg.norm(fs[~wrong_dets]-U*(transpose(R)*(R*(U.T.dot(fs[~wrong_dets]).T).T)))/linalg.norm(fs)
                resid_ = linalg.norm(fs-(Ts*G))/ linalg.norm(fs)
                debug( '%d. pass - unexplainable data: %.1f%%, unexplained data %.1f%%'%(step+1,resid*100, resid_*100))
                
            

        else:
            post_its = Bmat, diam, diar, diao

    debug('--------Linear method solution - time  %g' % (time.time()-t_all))
    
    
    G *= normData[None,:]/Tok.normTmat   #rescaling

    if Tok.transform_index in [0,4]:
        G[BdMat] = 0  #no emissivity out of the boundary (this emissivity is not contributing to retrofits!)
        if not SigmaGsample is None:   SigmaGsample[BdMat] = 0 
    

    retro *= normData[:,None]

    #set corrupted points and points with super large errorbars (also corrupted) to zero
    
    ind = where(~array(isfinite(err.diagonal()))[0] | all(abs(retro)<1e-5,0))[0]

        
    retro = asarray(err*retro.T).T
    
    if size(ind): 
        retro[:,ind] = (Tmat[ind]*(tokamak.Transform*G)).T


    
    


    return tvec, G, chi2, g, retro, [bnd]*tsteps,SigmaGsample





#def PresolveQR_(L,C, nx,ny,ndets):
    #"""
    #Source:N Terasaki, Y Hosodab, M Teranishia and N Iwamaa, Linear algebraic algorithms for high speed and stable reconstruction of plasma image
    #Y. Hosoda and T. Torii, A direct method for ill-posed linear operator equations--truncated least-square least-norm solutions

    #:param spmatrix L: Geometry matrix
    #:param spmatrix C: Smoothing matrix



    #"""
    #from scikits.sparse.cholmod import cholesky
    #from  scipy.linalg import qr, qr_multiply, solve_triangular

    
    
    #wrong_dets = squeeze(array(L.sum(1)==0))
    #L = L[where(~wrong_dets)[0],:]
    #M = min(sum(~wrong_dets),ndets )

    ##try:
        ##from sparsesvd import sparsesvd
        ##ut, s, vt = sparsesvd(L.tocsc(),size(L,0)) 
    ##except:
        ##print 'sparsesvd has failured'

        ##from scipy.sparse.linalg import svds
        ##ut, s, vt = svds(L.tocsc(),size(L,0)-1)  #BUG reconstruct only M-1 eigenvectors!
        ##ind = where(isfinite(s) & (s>1))[0][::-1]
        ##ut,s,vt = ut[:,ind].T,s[ind],vt[ind,:]
 

    #factor = cholesky(C) #C must by  regular
    
    #D = factor.D()
    
    #invsqrtD = sparse.spdiags(1/sqrt(D),(0,),len(D),len(D))
    ##t = time.time()
    ##vt = factor.apply_Pt(factor.solve_Lt(invsqrtD*vt.T)).T
    ##vt = (invsqrtD*factor.solve_L(factor.apply_P(vt.T))).T
    ##zero_ind = all(vt==0,axis=0)
    ##vt = vt[:,~zero_ind]

    ##vt = factor(vt.T).T
    #
    #

    ##wrong_dets = sum(fabs(ut),0) < 1e-8

    ##ut = ut[:,~wrong_dets]                           #odstraní nulové sloupcem vypočte
    ##LC =  dot(vt.T,(ut*s[:,None])).T    #L*C^-1  slow 
    
    #L = L.tocsr()
    #LC = invsqrtD*factor.solve_L(factor.apply_P(L.T))
    #LC = LC.T.toarray()

    #zero_ind = all(LC==0,axis=0)
    #LC = LC[:,~zero_ind]
    #
    #

    #Q1,S1,P = qr(LC.T, mode='economic',pivoting=True)       # 1. QR
    #S1_,P_ = qr(LC.T, mode='r',pivoting=True, check_finite=False) 
    ##allclose(LC.T[:, P], dot(Q1, S1))
    ##allclose(LC.T, dot(Q1, S1[:, P]))
    #
    #
    

    #D = diag(S1)                # čísla na diagonále nsjou setřízená!


    #S1 = S1/D[:,None]
    #S1 = S1[:,argsort(P)]


    #Q2,R2 = qr(S1.T, mode='economic')   # 2. QR

    #M =  D[:,None]*R2.T/D[None,:]


    ##Q3,R = qr(M, mode='economic')       # 3. QR

    ##V = dot(Q1,Q3) 
    #V_,R = qr_multiply(M, Q1, mode='right') # 3. QR


    #U = Q2

    #DR = diag(R)

    #R  = D[:,None]*R.T/(D*DR)[None,:]   # lower triangular now!

    #D = D*DR   

    ##V = factor(V).T     
    #V = zeros((len(zero_ind),V_.shape[1]))
    
    #
    #
    #V[~zero_ind,:] = V_
    
    #V = factor.apply_Pt(factor.solve_Lt(invsqrtD*V)).T

    ##U = dot(U,linalg.inv(R.T)).T 
    #U = solve_triangular(R,U.T,lower=True,unit_diagonal=True,overwrite_b=True) 
    ##U is not ortonormal anymore! bot it is not problem in the next solution

    #debug( 'PresolveQR done ')

    ##PlotBaseVectors({'U':U , 'V':V , 'D':D, 'R':R, 'wrong_dets':wrong_dets}, 44,66)
    ##D can be also negative!!!!
    #return {'U':U , 'V':V , 'D':D, 'R':R, 'wrong_dets':wrong_dets}






def PresolveQR(L,C, nx,ny,ndets,BdMat):
    """
    Source:N Terasaki, Y Hosodab, M Teranishia and N Iwamaa, Linear algebraic algorithms for high speed and stable reconstruction of plasma image
    Y. Hosoda and T. Torii, A direct method for ill-posed linear operator equations--truncated least-square least-norm solutions

    :param spmatrix L: Geometry matrix
    :param spmatrix C: Smoothing matrix



    """
    #from scikits.sparse.cholmod import cholesky

    #debug( 'PresolveQR')
    
    
    
    wrong_dets = squeeze(array(L.sum(1)==0))
    L = L[where(~wrong_dets)[0],:]
    m = min(sum(~wrong_dets),ndets )

    #try:
        #from sparsesvd import sparsesvd
        #ut, s, vt = sparsesvd(L.tocsc(),size(L,0)) 
    #except:
        #print 'sparsesvd has failured'

        #from scipy.sparse.linalg import svds
        #ut, s, vt = svds(L.tocsc(),size(L,0)-1)  #BUG reconstruct only M-1 eigenvectors!
        #ind = where(isfinite(s) & (s>1))[0][::-1]
        #ut,s,vt = ut[:,ind].T,s[ind],vt[ind,:]
 
    

    
    #factor = cholesky(C) #C must by  regular

    try:
        factor = cholesky(C, mode = 'simplicial')
    except:
        #solve case when B is singular
        C_trace = C.diagonal().sum()
        factor = cholesky(C,C_trace*1e-18, mode = 'simplicial')
    
    D = factor.D()
    
    invsqrtD = sparse.spdiags(1/sqrt(D),(0,),len(D),len(D))

    #vt = factor.apply_Pt(factor.solve_Lt(invsqrtD*vt.T)).T
    #vt = (invsqrtD*factor.solve_L(factor.apply_P(vt.T))).T


    #LC =  dot(vt.T,(ut*s[:,None])).T    #L*C^-1  slow 
    
    L = L.tocsr()
    LC = invsqrtD*factor.solve_L(factor.apply_P(L.T))
    LC = LC.toarray()

    zero_ind = all(LC==0,axis=1)
    LC = LC[~zero_ind,:]


    S1,P = qr(LC, mode='r',pivoting=True, check_finite=False)      #20ms        # 1. QR
    
    S1 = S1[:m,:]
    #P = P[:]

    D = diag(S1)            


    S1 = S1/D[:,None]
    S1 = S1[:,argsort(P)]


    Q2,R2 = qr(S1.T, mode='economic',check_finite=False,overwrite_a=True)       # 2. QR

    M =  D[:,None]*R2.T/D[None,:]


    #Q3,R = qr(M, mode='economic')       # 3. QR

    #V = dot(Q1,Q3) 
    #V_,R = qr_multiply(M, Q1, mode='right') # 3. QR
    R, = qr(M,mode='r',check_finite=False,overwrite_a=True)                    # 3. QR


    U = Q2

    DR = diag(R)

    R  = D[:,None]*R.T/(D*DR)[None,:]   # lower triangular now!

    D = abs(D*DR)    #make it positive (there is a freedom in the definition of sign )
    
    
    U = solve_triangular(R,U.T,lower=True,unit_diagonal=True,overwrite_b=True) 
    #U is not ortonormal anymore! bot it is not problem in the next steps of the solution


    V = zeros((L.shape[1],len(D) ),dtype='single')

    V[~zero_ind,:] = dot(U/D[:,None],LC.T).T  #10ms
    
     #theoreticaly this integration step can be applied later on the final solution, it will be  faster( but just about just 20% on total)
    V = factor.apply_Pt(factor.solve_Lt(invsqrtD*V)).T 
    U = U.T
    
    #plot(D)
    #show()
    

    return {'U':matrix(U,copy=False) , 'V':matrix(V,copy=False) , 'D':D, 
                'R':matrix(R,copy=False), 'wrong_dets':wrong_dets}







def PresolveSVD(L,C, nx,ny,nl,BdMat):
    """
    Linear algebraic algorithms for high speed and stable reconstruction of plasma image

    Source: N Terasaki, Y Hosodab, M Teranishia and N Iwamaa

    :param spmatrix L: Geometry matrix
    :param spamtrix C: Smoothing matrix
    """
   
    wrong_dets = squeeze(array(L.sum(1)==0))
    L = L[where(~wrong_dets)[0],:]
    M = min(sum(~wrong_dets),nl )
   
   
    #from scikits.sparse.cholmod import cholesky

    #ut, s, vt = sparsesvd(L.tocsc(),M) # jde této znalosti ještě nějak zneužít?
    #factor = cholesky(C)
    #vt = factor(vt.T).T
    #LC = dot(ut.T, dot(diag(s), vt)) #L*C^-1
    #LC = L*inv(C.toarray())
    

    #debug('SVD start')

    #try:
    #from sparsesvd import sparsesvd
    #ut, s, vt = sparsesvd(L.tocsc(),M) 
    #except:
        #print 'sparsesvd has failured'

    from scipy.sparse.linalg import svds
    #70ms
    try:
        u, s, vt = svds(L.tocsc(),M-1)  #BUG reconstruct only M-1 eigenvalues!!  
        ind = where(isfinite(s) & (s>1))[0][::-1] #remove noise vectors
        s = s[ind]
        vt = vt[ind,:]
        ut = u[:,ind].T
    except:
        u, s, vt = svd(L.todense(),full_matrices=False, overwrite_a=True )
        ut = u.T


    #ut, s, vt = sparsesvd(L.tocsc(),M)  
    #debug('SVD end')
    #debug('cholesky start')
        
    
    
    
    #16ms
    factor = cholesky(C)   #spíš by tu měla být "odmocnina", ne?
    #debug('cholesky end')

    #http://www.mathworks.com/matlabcentral/newsreader/view_thread/32825
    D = factor.D()
    
    invsqrtD = sparse.spdiags(1/sqrt(D),(0,),len(D),len(D))
    #t = time.time()
    #41ms
    vt = (invsqrtD*factor.solve_L(factor.apply_P(vt.T))).T
    #vt = factor(vt.T).T
    
    v = vt.T*s[None,:]
    #from time import time
    #Compute only for nonzero parts of the profile 
    zero_ind = all(v==0,axis=1) 
    #vt2 = copy(vt)

    #t = time()
    #LC = dot(ut.T, vt[:,~zero_ind]*s[:,None])            #L*C^-1 
    #u2,s2,vt_ = svd(LC,full_matrices=False, overwrite_a=True)
    #vt2 = zeros(L.shape)
    #vt2[:,~zero_ind] = vt_
    #print time()-t

    
    #t = time()
    #nový upgrade rozkladu LC^-1, aby to bylo rychlejší 
    #200ms!!!
    P,T = qr(v[~zero_ind],overwrite_a=True, mode='economic')
    
    #je C symetrická matice????
    #U,s,V = svd(dot(R,T.T),full_matrices=False, overwrite_a=False)
    #V,s,U = svd(T,full_matrices=False, overwrite_a=True)
    U,s,V = svd(T.T,full_matrices=False, overwrite_a=True)

    
    ut = dot(ut.T,U)
    vt_ = copy(vt)
    #31ms
    vt_[:,~zero_ind] = dot(V,P.T)

    #print time()-t
    

    #vt = factor(vt.T).T
    #40ms
    vt = factor.apply_Pt(factor.solve_Lt(invsqrtD*vt_.T)).T
    




    #PlotBaseVectors({'U':ut , 'V':vt , 'D':s,'wrong_dets':wrong_dets}, 44,66)

    return {'U':matrix(ut) , 'V':matrix(vt) , 'D':s, 'wrong_dets':wrong_dets}






def PresolveSVD_(L,C, nx,ny,nl):
    """
    Linear algebraic algorithms for high speed and stable reconstruction of plasma image

    Source: N Terasaki, Y Hosodab, M Teranishia and N Iwamaa

    :param spmatrix L: Geometry matrix
    :param spamtrix C: Smoothing matrix
    """
   
    wrong_dets = squeeze(array(L.sum(1)==0))
    L = L[where(~wrong_dets)[0],:]
    M = min(sum(~wrong_dets),nl )
   
   
    #from scikits.sparse.cholmod import cholesky



    #from scipy.sparse.linalg import svds

    factor = cholesky(C)   #spíš by tu měla být "odmocnina", ne?
    #debug('cholesky end')

    #http://www.mathworks.com/matlabcentral/newsreader/view_thread/32825
    D = factor.D()
    
    invsqrtD = sparse.spdiags(1/sqrt(D),(0,),len(D),len(D))
    #t = time.time()
    #41ms
    
        
    L = L.tocsr()
    LC = invsqrtD*factor.solve_L(factor.apply_P(L.T))
    LC = LC.T.toarray()

    zero_ind = all(LC==0,axis=0)
    LC = LC[:,~zero_ind]
    
    #LC = (invsqrtD*factor.solve_L(factor.apply_P(L.T))).T
    
        
    
   
    #vt = factor(vt.T).T
    
    #v = vt.T*s[None,:]
    #from time import time
    #Compute only for nonzero parts of the profile 
    #zero_ind = all(v==0,axis=1) 
    #vt2 = copy(vt)

    #t = time()
    #LC = dot(ut.T, vt[:,~zero_ind]*s[:,None])            #L*C^-1 
    #u2,s2,vt_ = svd(LC,full_matrices=False, overwrite_a=True)
    #vt2 = zeros(L.shape)
    #vt2[:,~zero_ind] = vt_
    #print time()-t

    
    #t = time()
    #nový upgrade rozkladu LC^-1, aby to bylo rychlejší 
    #200ms!!!
    #P,T = qr(LC,overwrite_a=True, mode='economic')
    
    #je C symetrická matice????
    #U,s,V = svd(dot(R,T.T),full_matrices=False, overwrite_a=False)
    #V,s,U = svd(T,full_matrices=False, overwrite_a=True)
    U,s,Vt_ = svd(LC,full_matrices=False, overwrite_a=True)

    
    
    #ut = dot(ut.T,U).T
    #vt_ = copy(vt)
    #31ms
    Vt = zeros(L.shape)
    Vt[:,~zero_ind] = Vt_

    #print time()-t
    

    #vt = factor(vt.T).T
    #40ms
    vt = factor.apply_Pt(factor.solve_Lt(invsqrtD*Vt.T)).T
    




    #PlotBaseVectors({'U':ut , 'V':vt , 'D':s,'wrong_dets':wrong_dets}, 44,66)

    return {'U':matrix(U) , 'V':matrix(vt) , 'D':s, 'wrong_dets':wrong_dets}









def PresolveSVD3(K,B,nx,ny, ndets,BdMat):
    """
    :param spmatrix K: Geometry matrix
    :param spamtrix B: Smoothing matrix D.T*D
    Solve SVD by cholesky and svd decomposition. 
    
    BUG Slow due to some bug in cholmod :( ???????

    """
    #T = time.time()
    #print K.shape,ndets

    wrong_dets = squeeze(array(K.sum(1)==0))
    
    if all(wrong_dets):
        raise Exception('Something is wrong with normalized geometry matrix')
        
        
        
    K = K[where(~wrong_dets)[0],:]  #1.2ms
    k = min(sum(~wrong_dets),ndets,K.shape[1] )
    


    #print K.shape
    #exit()
    #B = B*1e4
    #K = K*1e2
    #solve generalized eigen value problem K.T*K*x = lambda*C*x


    if not debug: warnings.filterwarnings("error")

    
    
    #null = B*ones(2400)==0
    B_trace = B.diagonal().sum()
    #singular = abs(sum( B *(B.sum(1).T==0).T)/B_trace)  <1e-20
    
    
    
    #save('B',B)
    #from scikits.sparse import cholmod
        
       
    #B = load('B.npy').item()
    #from scikits.sparse.cholmod import cholesky
    #T = time.time()
    #%timeit cholesky(B, mode = 'simplicial')
    #print  time.time()-T
    
    #print K.shape
    #import IPython
    #IPython.embed()
    #print  'init', time.time()-T
    T = time.time()
    try:
        #assert ~singular, 'singular matrix'
        F = cholesky(B, mode = 'simplicial') #5ms!!
    except:# CholmodError, CholmodWarning:
        #solve case when B is singular
        F = cholesky(B,B_trace*1e-18, mode = 'simplicial')
        #pass
        
        
    #print  'cholesky',time.time()-T
    T = time.time()

        
    if not debug: warnings.resetwarnings()
    
    if sparse.isspmatrix_csc(K):   K = K.tocsr()
    LPK = F.solve_L(F.apply_P(K.T))  #8ms
    #print 'solved L', time.time()-T

    #D = 
    
    
    #invsqrtD = sparse.spdiags(1/sqrt(D),(0,),len(D),len(D))
    zero_ins = array(LPK.sum(1)==0).ravel()|(F.D()<=0)
    #pcolormesh(F.apply_Pt(double(zero_ins)).reshape(ny,nx,order='F')+BdMat.reshape(ny,nx,order='F'))
    #figure()
    #pcolormesh(BdMat.reshape(ny,nx,order='F'))
    #show()

    #import IPython
    #IPython.embed()
        
    
    if sparse.issparse(LPK): 
        LPK = LPK.toarray() #2.5ms
        
    LPK = LPK[~zero_ins].T  #.5ms
    #LPK_ = (invsqrtD*LPK_).T
    #import IPython
    #IPython.embed()
    #print  'F.solve_L',time.time()-T
    T = time.time()

    #LPK_ = LPK.toarray()  #memory consuming 
    sqrtD = sqrt(F.D()[~zero_ins])
    
    
    #print 'to array'
    
    
    LPK/= sqrtD

    
    #LPK_ = (invsqrtD*LPK_).T.toarray()  #LC^-1

    #calculate only for nonzero coloumns of LPK_ matrix

    #LPK = LPK.astype('float32')
    #del LPK_
    #LPK = LPK.toarray()

    #LPK = LPK.T

    
    #from scipy.sparse.linalg import svds

    #u, s, vt = svds(LPK,k-1)  #BUG reconstruct only M-1 eigenvalues!!  
    #ind = where(isfinite(s) & (s>1))[0][::-1] #remove noise vectors
    #s = s[ind]
    #vt = vt[ind,:]
    #ut = u[:,ind].T
    
    
    
    
    #NOTE ugly way how to calulate fast SVD
    
    #t = time.time()
    #u,s,v = svd(LPK, full_matrices=False, overwrite_a=True, check_finite=False)

    #r, = qr(LPK, overwrite_a=False,mode='r', pivoting=False, check_finite=False)
    #r = r[:,:k]
    #u,s,v = svd(r, full_matrices=False, overwrite_a=True, check_finite=False)
    #invR = dot(v.T, u.T/s[:,None] )
    ##vt_ = dot(invR,LPK )
    #print time.time()-t
    #t = time.time()
    #tol = max(dim(m))*max(D)*.Machine$double.eps
#NOTE http://hackersome.com/p/arbenson/mrtsqr

    #T = ascontiguousarray(LPK.T)
    #dot(T.T, T)
    
    #from scipy.sparse.linalg import eigsh, svds
        #from scipy.sparse.linalg import svds

    #eigsh( LPK* LPK.T, k=k-1)
    #u, s, vt_ = svds(LPK,k-1)  #BUG reconstruct only M-1 eigenvalues!!  
    #print 'BUG!!!!!\!!'
    #u,s,vt_ = svd(LPK_, full_matrices=False, overwrite_a=True, check_finite=False)
    
    def fast_svd(A):

        transposed = False
        if A.shape[0] >  A.shape[1]:
            transposed = True
            A = A.T
            
            
        #T2 = time.time()

        LL = dot(A, A.T)  #Gram matrix
        #print  'Gram matrix',time.time()-T2
        #T2 = time.time()
  
        s,u = eigh(LL,overwrite_a=True, check_finite=False,lower=True)  
        #print  'eigh',time.time()-T2
        #T2 = time.time()

        del LL
        s,u = sqrt(maximum(s,0))[::-1], u[:,::-1] 

        try:
            rank = where(cumprod((diff(log(s[s!=0]))>-5)|(s[s!=0][1:]>median(s))))[0][-1]+2
        except:
            print(diff(log(s[s!=0])))
            for r in s: print(s)
         
    
        #rank = len(s)
        S = s[:rank]
        U = u[:,:rank]
        #rank = 

        VT = dot(U.T/S[:,None], A)  #slow...
        #print  'backprojection',time.time()-T2

        
        if transposed:
            return VT.T, S,U.T, rank
        else:
            return U,S,VT,rank 
        



    #T2 = time.time()
    #import IPython
    #IPython.embed()
    
    #LPK = random.randn(200, 6000)
    
    
    #%timeit dot(LPK, LPK.T) 
    
    #print LPK.shape, type(LPK), LPK.dtype
    U,S,vt_, rank = fast_svd(LPK) #22ms
    
    #print  'fast_svd',time.time()-T2
    #T = time.time()
    
        
    
    
    
    del LPK

    #LL = dot(LPK_, LPK_.T)
    #s,u = eigh(LL,overwrite_a=True, check_finite=False,lower=True)
    ##print 'LL'

    #del LL
    #s,u = sqrt(maximum(s,0))[::-1], u[:,::-1]
    ##dot(u[:,:rank].T/D[:,None], LPK)
    ##ut = u[:,ind].T
    #
    #
    
    ##r = qr(ascontiguousarray(LPK.T),  mode='r',check_finite=False)

    ##from scipy.linalg import cholesky
    
    ##L = cholesky( dot(LPK, LPK.T),lower=False, overwrite_a=False, check_finite=False)
    ##U,S,V = svd(L, full_matrices=False, compute_uv=True, overwrite_a=True, check_finite=False)


    
    
    
    

    ##estimate the rank
    ##print s
    ##plot(s)
    ##show()  
    ##print where(cumprod((diff(log(s[s!=0]))>-5)))
    ##try:
    #rank = where(cumprod((diff(log(s[s!=0]))>-5)|(s[s!=0][1:]>median(s))))[0][-1]+2
    ##except:
        #
        #

    ##print 'BUG!!!!!!'
    ##rank-= 1
    ##valid_ind = slice(None, rank+1)
    
    #S = s[:rank]
    #U = u[:,:rank]
    
    ##if any(S == 0) or any(isnan(S)):
        #
        #
    #vt_ = dot(U.T/S[:,None], LPK_)  #slow...
    #del LPK_

    #print 'VT'

    #wrong_dets[~wrong_dets] |=~valid_ind
#print time.time()-t
    #t = time.time()

    #ordinary SVD (slow)
    #u, s_, vt_ = svd(LPK,  full_matrices=False,check_finite = False, overwrite_a=True)
    #print time.time()-t


    #valid_ind = isfinite(s)
    #valid_ind[valid_ind] &= s[valid_ind] > 0
    
    
    vt = zeros(( K.shape[1], len(S)) )
    
    
    
    #print vt_.shape, vt[:,~zero_ins].shape, vt.shape
    
    vt_ /= sqrtD
    vt[~zero_ins] = vt_.T
    #print 'VT2'
    #vt = ascontiguousarray(vt.T)

    #vt = vt.astype(single,copy=False)
    del vt_
    V = F.solve_Lt(vt) #7ms
    del vt
    V = F.apply_Pt(V).T #1.5ms
    #print 'Solve L'
    
    #print  'solve_Lt',time.time()-T
    T = time.time()


    #valid_ind = slice(0, where(s/s.max() > 1e-7)[0][-1]+1)
    #valid_ind = slice(None,None)
  
    #sort from the most important to the least
    #valid_ind = S > 1e-6
    #S =S
    #U = (K*V)[:,valid_ind].T/S[:,newaxis]
    #U = u.T
  ##V = V[:,valid_ind]
    #V = V.T

    
    
    
    
    
    #V = V[:,valid_ind]
    #define decomposition in the same form as for the other methods
    #U =  (K*V).T
    #V = V.T
    
    

    
    
    
    
    #U*S*V*B= K!!
    
    #iterative upgrade of the existing decompostion for a new B. But U is not orthogonal!? 
    #F = cholesky(B)
    #HK = asarray(F(K.T).T.todense())
    #V_ = dot(U.T/S[:,None],HK)
    #U_ = K*V_.T
    #D_ = S*linalg.norm(U,axis=0)
    
    
    #sqrt(B.diagonal().sum())/sum(S)
    #1/sqrt(linalg.inv((K*K.T).todense()).trace())/(1/sum(1/S))

    #S[0],S[-1]
    
    #dot(linalg.inv((K*K.T).todense()),K.T)
    #lmin = 1/sqrt(linalg.inv((K*K.T).todense()).trace())
    #lmax = sqrt(B.diagonal().sum())
    #fill =  B.nnz/float(nx*ny)
    #lmax = sqrt((square(K.data).sum())/B.diagonal().sum()*(nx*ny)**2)*20
    #lmin = sqrt((square(K.data).sum())/B.diagonal().sum()*nx*ny)/5
    #print 
    #print '======== =%.3e  %.3e, %.3e, %.3e'%( lmin,S[-1], lmax,S[0])#, sqrt(lmin*lmax)/median(S)
    
    #K_ = dot(dot(U,diag(1/S)),V)
    #print 'Decomposition accuracy',linalg.norm((K_* K.T)-eye(K.shape[0]))#-len(S)+rank
    #exit()


    
    
    #print   time.time()-T

    #exit()

    return {'K':K,'U':matrix(U,copy=False), 'V':matrix(V,copy=False) , 'D':S, 'wrong_dets':wrong_dets}


def PresolveSVD3_sparse(K,B,nx,ny, ndets):
    """
    :param spmatrix K: Geometry matrix
    :param spamtrix B: Smoothing matrix D.T*D
    Solve SVD by cholesky and svd decomposition. 
    
    memore effective but slow version

    """

    wrong_dets = squeeze(array(K.sum(1)==0))
    K = K[where(~wrong_dets)[0],:]
    k = min(sum(~wrong_dets),ndets,K.shape[1] )
    npix = K.shape[1] 

    #solve generalized eigen value problem K.T*K*x = lambda*C*x

    #from scikits.sparse.cholmod import cholesky
    #from  scipy.linalg import eigh,svd,qr
    

    try:
        F = cholesky(B)
    except:
        #solve case whan B is singular
        B_trace = B.diagonal().sum()
        F = cholesky(B,B_trace*1e-20/nx/ny)
    
    
    if sparse.isspmatrix_csc(K):
        K = K.tocsr()
        
    LPK = F.solve_L(F.apply_P(K.T))
    #del K
    #print 'solved L'

    D = F.D()
    
    
    invsqrtD = sparse.spdiags(1/sqrt(D),(0,),len(D),len(D))
    zero_ins = array(LPK.sum(1)==0).ravel()
    LPK = (invsqrtD*LPK)
    LPK = LPK[where(~zero_ins)[0],:].T

    #LPK_ = LPK.T.toarray()  #memory consuming 

    #print 'to array'
    #LPK_/= sqrt(D[~zero_ins])
    #LPK_ = (invsqrtD*LPK_).T.toarray()  #LC^-1

    #calculate only for nonzero coloumns of LPK_ matrix

    #LPK = LPK.astype('float32')
    #del LPK_
    #LPK = LPK.toarray()

    #LPK = LPK.T

    
    #from scipy.sparse.linalg import svds

    #u, s, vt = svds(LPK,k-1)  #BUG reconstruct only M-1 eigenvalues!!  
    #ind = where(isfinite(s) & (s>1))[0][::-1] #remove noise vectors
    #s = s[ind]
    #vt = vt[ind,:]
    #ut = u[:,ind].T
    
    
    
    
    #NOTE ugly way how to calulate fast SVD
    
    #t = time.time()
    #u,s,v = svd(LPK, full_matrices=False, overwrite_a=True, check_finite=False)

    #r, = qr(LPK, overwrite_a=False,mode='r', pivoting=False, check_finite=False)
    #r = r[:,:k]
    #u,s,v = svd(r, full_matrices=False, overwrite_a=True, check_finite=False)
    #invR = dot(v.T, u.T/s[:,None] )
    #vt_ = dot(invR,LPK )
    #print time.time()-t
    #t = time.time()
    #tol = max(dim(m))*max(D)*.Machine$double.eps
#NOTE http://hackersome.com/p/arbenson/mrtsqr

    #T = ascontiguousarray(LPK.T)
    #dot(T.T, T)
    
    #from scipy.sparse.linalg import eigsh, svds
        #from scipy.sparse.linalg import svds

    #eigsh( LPK* LPK.T, k=k-1)
    #u, s, vt = svds(LPK,k-1)  #BUG reconstruct only M-1 eigenvalues!!  
    print('calc LL')
    LL =  LPK* LPK.T
    #LL.data = LL.data.astype('single')
    LL = LL.toarray()
    
    #LPK = LPK.toarray()
    print('eigh')

    #LL = dot( LPK, LPK.T)
    s,u = eigh(LL,overwrite_a=True, check_finite=False,lower=True)
    #print 'LL'

    del LL
    s[s<0] = 0
    s,u = sqrt(s)[::-1], u[:,::-1]
    #dot(u[:,:rank].T/D[:,None], LPK)
    #ut = u[:,ind].T
    
    
    
    #r = qr(ascontiguousarray(LPK.T),  mode='r',check_finite=False)

    #from scipy.linalg import cholesky
    
    #L = cholesky( dot(LPK, LPK.T),lower=False, overwrite_a=False, check_finite=False)
    #U,S,V = svd(L, full_matrices=False, compute_uv=True, overwrite_a=True, check_finite=False)



    
    
    #estimate the rank
    #print s
    #plot(s)
    #show()  
    #print where(cumprod((diff(log(s[s!=0]))>-5)))
    #try:
    rank = where(cumprod((diff(log(s[s!=0]))>-5)|(s[s!=0][1:]>median(s))))[0][-1]+2
    #except:

    #print 'BUG!!!!!!'
    #rank-= 1
    #valid_ind = slice(None, rank+1)
    
    S = abs(s[:rank])
    U = u[:,:rank]
    print('calc V')
    indptr = list(LPK.indptr)
    for i in where(zero_ins)[0]:  indptr.insert(i,indptr[i])    
    LPK = sparse.csc_matrix((LPK.data, LPK.indices, indptr), shape=( ndets, len(zero_ins)) )

    
    vt = zeros(( rank,npix),dtype=single )
    chunk = 1000
    US = U.T/S[:,None]
    for i in range(npix/chunk+1):
        ind = slice(i*chunk,min((i+1)*chunk,npix)) ;
        vt[:,ind] = US* LPK[:,ind]  #slow...
        
        
        
    del LPK,US
    
    
    #print 'VT'

    #wrong_dets[~wrong_dets] |=~valid_ind
#print time.time()-t
    #t = time.time()

    #ordinary SVD (slow)
    #u, s_, vt_ = svd(LPK,  full_matrices=False,check_finite = False, overwrite_a=True)
    #print time.time()-t

    
    #valid_ind = isfinite(s)
    #valid_ind[valid_ind] &= s[valid_ind] > 0
    

    #vt.indptr = array(indptr)
    #vt.shape = (vt_.shape[0],len(zero_ins))
    #vt.set_shape((vt_.shape[0],len(zero_ins)))
    #vt[:,~zero_ins] = vt_.todense()   
    
    #sparse.csc_matrix(vt).indptr
    
    #V = F.apply_Pt(F.solve_Lt(invsqrtD*vt.T)).T
    
    
    #print vt_.shape, vt[:,~zero_ins].shape, vt.shape
    #vt[:,~zero_ins] = vt_
    vt/= sqrt(D)
    #print 'VT2'
    print('calc CV')

    #vt = vt.astype(single,copy=False)
    #del vt_
    chunk = 100
    for i in range(rank/chunk+1):
        ind = slice(i*chunk,min((i+1)*chunk,rank))
        vt[ind,:] =  F.apply_Pt(F.solve_Lt(vt[ind,:].T.astype('double'))).T
         #F.apply_Pt(double(V))
    #V = vt
    #del vt
    #V = F.apply_Pt(double(V)).T
    #print 'Solve L'


    #valid_ind = slice(0, where(s/s.max() > 1e-7)[0][-1]+1)
    #valid_ind = slice(None,None)
  
    #sort from the most important to the least
    #valid_ind = S > 1e-6
    #S =S
    #U = (K*V)[:,valid_ind].T/S[:,newaxis]
    #U = u.T
  ##V = V[:,valid_ind]
    #V = V.T

    
    
    
    
    
    #V = V[:,valid_ind]
    #define decomposition in the same form as for the other methods
    #U =  (K*V).T
    #V = V.T
    
    
    #K_ = dot(dot(U,diag(1/S)),V)
    #print 'Decomposition accuracy',linalg.norm((K_* K.T)-eye(K.shape[0]))-len(s)+rank
    #exit()

    
    
    
    
    #U*S*V*B= K!!
    
    #iterative upgrade of the existing decompostion for a new B. But U is not orthogonal!? 
    #F = cholesky(B)
    #HK = asarray(F(K.T).T.todense())
    #V_ = dot(U.T/S[:,None],HK)
    #U_ = K*V_.T
    #D_ = S*linalg.norm(U,axis=0)


    return {'U':matrix(U,copy=False), 'V':matrix(vt,copy=False) , 'D':S, 'wrong_dets':wrong_dets}

    


def PresolveSVD3_sparse2(K,B,nx,ny, ndets,nsqval=200):
    """
    :param spmatrix K: Geometry matrix
    :param spamtrix B: Smoothing matrix D.T*D
    Solve SVD by cholesky and svd decomposition. 
    
    memore effective but slow version

    """
    
    
    
    wrong_dets = squeeze(array(K.sum(1)==0))
    K = K[where(~wrong_dets)[0],:]
    k = min(sum(~wrong_dets),ndets,K.shape[1] )
    npix = K.shape[1] 
    K.indices = K.indices.astype('uint16',copy=True) #BUG max 2**16 detectors!
    

    #solve generalized eigen value problem K.T*K*x = lambda*C*x

    #from scikits.sparse.cholmod import cholesky
    #from  scipy.linalg import eigh,svd,qr
    

    try:
        F = cholesky(B)
    except:
        #solve case whan B is singular
        B_trace = B.diagonal().sum()
        F = cholesky(B,B_trace*1e-16/nx/ny)
        
    PK = F.apply_P(K.T.tocsc())
    LPK = F.solve_L(PK) 
    del PK
    LPK.data = LPK.data.astype('float32',copy=False)
    del K

    D = F.D().astype('float32',copy=False)
    zero_ins = array(LPK.sum(1)==0).ravel()

    chunk = 100
    for i in range(size(LPK)/chunk+1):
        ind = slice(i*chunk,min((i+1)*chunk,size(LPK)))
        LPK.data[ind]/= sqrt(D[LPK.indices[ind]])

    
    
    
    #def sparse2dense(A,ind):
        
        
        
    #BUG memory consumption!!
    LPK = LPK[where(~zero_ins)[0],:].T #udělat tohle seříznutí dřív?? 
    #invsqrtD = sparse.spdiags(1/sqrt(D[~zero_ins]),(0,),sum(~zero_ins),sum(~zero_ins)  )



    #LPK_ = LPK.toarray()  #memory consuming 

    #print 'to array'
    #LPK_/= sqrt(D[~zero_ins])
    #LPK_ = (invsqrtD*LPK_).T.toarray()  #LC^-1

    #calculate only for nonzero coloumns of LPK_ matrix

    #LPK = LPK.astype('float32')
    #del LPK_
    #LPK = LPK.toarray()

    #LPK = LPK.T


    #from scipy.sparse.linalg import svds
    #from scipy.sparse.linalg import aslinearoperator
    
    #from scipy.linalg.interpolative import svd
    #u, s, vt_ = svd(aslinearoperator(LPK),nsqval)
    print('calc svds')
    #from fsvd import fsvd
    #u, s, vt_ = fsvd(LPK)
    u, s, vt_ = svds(LPK,min(nsqval,k-1))  #BUG reconstruct only M-1 eigenvalues!!  
    #ind = where(isfinite(s) & (s>1))[0][::-1] #remove noise vectors
    s[isnan(s)] = 0
    #s = s[ind]
    #vt = vt[ind,:]
    #ut = u[:,ind].T
    
    
    
    
    #NOTE ugly way how to calulate fast SVD
    
    #t = time.time()
    #u,s,v = svd(LPK, full_matrices=False, overwrite_a=True, check_finite=False)

    #r, = qr(LPK, overwrite_a=False,mode='r', pivoting=False, check_finite=False)
    #r = r[:,:k]
    #u,s,v = svd(r, full_matrices=False, overwrite_a=True, check_finite=False)
    #invR = dot(v.T, u.T/s[:,None] )
    #vt_ = dot(invR,LPK )
    #print time.time()-t
    #t = time.time()
    #tol = max(dim(m))*max(D)*.Machine$double.eps
#NOTE http://hackersome.com/p/arbenson/mrtsqr

    #T = ascontiguousarray(LPK.T)
    #dot(T.T, T)
    
    #from scipy.sparse.linalg import eigsh, svds
        #from scipy.sparse.linalg import svds

    #eigsh( LPK* LPK.T, k=k-1)
    #u, s, vt = svds(LPK,k-1)  #BUG reconstruct only M-1 eigenvalues!!  
    #print 'calc LL'
    #LL =  LPK* LPK.T
    #LL.data = LL.data.astype('single')
    #LL = LL.toarray()
    
    #LPK = LPK.toarray()
    #print 'eigh'

    #LL = dot( LPK, LPK.T)
    #s,u = eigh(LL,overwrite_a=True, check_finite=False,lower=True)
    #print 'LL'

    #del LL
    #s[s<0] = 0
    s,u,vt_ = s[::-1], u[:,::-1], vt_[::-1]
    #dot(u[:,:rank].T/D[:,None], LPK)
    #ut = u[:,ind].T
    
    
    
    #r = qr(ascontiguousarray(LPK.T),  mode='r',check_finite=False)

    #from scipy.linalg import cholesky
    
    #L = cholesky( dot(LPK, LPK.T),lower=False, overwrite_a=False, check_finite=False)
    #U,S,V = svd(L, full_matrices=False, compute_uv=True, overwrite_a=True, check_finite=False)



    
    
    #estimate the rank
    #print s
    #plot(s)
    #show()  
    #print where(cumprod((diff(log(s[s!=0]))>-5)))
    #try:
    rank = where(cumprod((diff(log(s[s!=0]))>-5)|(s[s!=0][1:]>median(s))))[0][-1]+2
    #except:

    #print 'BUG!!!!!!'
    #rank-= 1
    #valid_ind = slice(None, rank+1)
    
    S = abs(s[:rank])
    U = u[:,:rank]
    vt_ = vt_[:rank,:]

    #print 'calc V'
    #indptr = list(LPK.indptr)
    #for i in where(zero_ins)[0]:  indptr.insert(i,indptr[i])    
    #LPK = sparse.csc_matrix((LPK.data, LPK.indices, indptr), shape=( ndets, len(zero_ins)) )

    
    vt = zeros(( rank,npix),dtype=single )
    vt[:,~zero_ins] = vt_
    #chunk = 1000
    ##US = U.T/S[:,None]
    #for i in range(npix/chunk+1):
        #ind = slice(i*chunk,min((i+1)*chunk,npix)) ;
        #vt[:,ind] = vt_[:,ind]  #slow...
        
        
        
    #del LPK,US
    
    
    #print 'VT'

    #wrong_dets[~wrong_dets] |=~valid_ind
#print time.time()-t
    #t = time.time()

    #ordinary SVD (slow)
    #u, s_, vt_ = svd(LPK,  full_matrices=False,check_finite = False, overwrite_a=True)
    #print time.time()-t

    
    #valid_ind = isfinite(s)
    #valid_ind[valid_ind] &= s[valid_ind] > 0
    

    #vt.indptr = array(indptr)
    #vt.shape = (vt_.shape[0],len(zero_ins))
    #vt.set_shape((vt_.shape[0],len(zero_ins)))
    #vt[:,~zero_ins] = vt_.todense()   
    
    #sparse.csc_matrix(vt).indptr
    
    #V = F.apply_Pt(F.solve_Lt(invsqrtD*vt.T)).T
    
    
    #print vt_.shape, vt[:,~zero_ins].shape, vt.shape
    #vt[:,~zero_ins] = vt_
    vt/= sqrt(D)
    #print 'VT2'
    print('calc CV')

    #vt = vt.astype(single,copy=False)
    #del vt_
    chunk = 100
    for i in range(rank/chunk+1):
        ind = slice(i*chunk,min((i+1)*chunk,rank))
        vt[ind,:] =  F.apply_Pt(F.solve_Lt(vt[ind,:].T.astype('double'))).T
         #F.apply_Pt(double(V))
    #V = vt
    #del vt
    #V = F.apply_Pt(double(V)).T
    #print 'Solve L'


    #valid_ind = slice(0, where(s/s.max() > 1e-7)[0][-1]+1)
    #valid_ind = slice(None,None)
  
    #sort from the most important to the least
    #valid_ind = S > 1e-6
    #S =S
    #U = (K*V)[:,valid_ind].T/S[:,newaxis]
    #U = u.T
  ##V = V[:,valid_ind]
    #V = V.T

    
    
    
    
    
    #V = V[:,valid_ind]
    #define decomposition in the same form as for the other methods
    #U =  (K*V).T
    #V = V.T
    
    
    #K_ = dot(dot(U,diag(1/S)),V)
    #print 'Decomposition accuracy',linalg.norm((K_* K.T)-eye(K.shape[0]))-len(s)+rank
    #exit()

    
    
    
    
    #U*S*V*B= K!!
    
    #iterative upgrade of the existing decompostion for a new B. But U is not orthogonal!? 
    #F = cholesky(B)
    #HK = asarray(F(K.T).T.todense())
    #V_ = dot(U.T/S[:,None],HK)
    #U_ = K*V_.T
    #D_ = S*linalg.norm(U,axis=0)


    return {'U':matrix(U,copy=False), 'V':matrix(vt,copy=False) , 'D':S, 'wrong_dets':wrong_dets}


def PresolveSVD2(K,B,nx,ny, ndets):
    """
    :param spmatrix K: Geometry matrix
    :param spamtrix B: Smoothing matrix D.T*D
    Solve SVD by cholesky and svd decomposition. 

    """

    wrong_dets = squeeze(array(K.sum(1)==0))
    K = K[where(~wrong_dets)[0],:]
    k = min(sum(~wrong_dets),ndets,K.shape[1] )
    

    #solve generalized eigen value problem K.T*K*x = lambda*C*x

    #from scikits.sparse.cholmod import cholesky
    
  
    try:
        F = cholesky(B)
    except:
        #solve case whan B is singular
        B_trace = B.diagonal().sum()
        F = cholesky(B,B_trace*1e-16/nx/ny)
        
        
    LPK_ = F.solve_L(F.apply_P(K.T.tocsc()))  

    D = F.D()

    invsqrtD = sparse.spdiags(1/sqrt(D),(0,),len(D),len(D))
    LPK_ = (invsqrtD*LPK_).T
    #LPK_ = (invsqrtD*LPK_).T.toarray()  #LC^-1
    
    #calculate only for nonzero coloumns of LPK_ matrix
    zero_ins = array(LPK_.sum(0)==0).ravel()
    LPK = LPK_[:,where(~zero_ins)[0]].todense()

    
    
    
    
    #LPK = LPK.T

    
    #from scipy.sparse.linalg import svds

    #u, s, vt = svds(LPK,k-1)  #BUG reconstruct only M-1 eigenvalues!!  
    #ind = where(isfinite(s) & (s>1))[0][::-1] #remove noise vectors
    #s = s[ind]
    #vt = vt[ind,:]
    #ut = u[:,ind].T
    
    
    
    
    #NOTE ugly way how to calulate fast SVD
    
    #t = time.time()
    #u,s,v = svd(LPK, full_matrices=False, overwrite_a=True, check_finite=False)

    #r, = qr(LPK, overwrite_a=False,mode='r', pivoting=False, check_finite=False)
    #r = r[:,:k]
    #u,s,v = svd(r, full_matrices=False, overwrite_a=True, check_finite=False)
    #invR = dot(v.T, u.T/s[:,None] )
    #vt_ = dot(invR,LPK )
    #print time.time()-t
    #t = time.time()
    #tol = max(dim(m))*max(D)*.Machine$double.eps
#NOTE http://hackersome.com/p/arbenson/mrtsqr

    #T = ascontiguousarray(LPK.T)
    #dot(T.T, T)
    
    #from scipy.sparse.linalg import eigsh, svds
        #from scipy.sparse.linalg import svds

    #eigsh( LPK* LPK.T, k=k-1)
    #u, s, vt = svds(LPK,k-1)  #BUG reconstruct only M-1 eigenvalues!!  
    
    
    s,u = eigh(dot(LPK, LPK.T),overwrite_a=True, check_finite=False,lower=True)
    s,u = sqrt(abs(s))[::-1], u[:,::-1]
    #dot(u[:,:rank].T/D[:,None], LPK)
    #ut = u[:,ind].T
    
    
    
    #r = qr(ascontiguousarray(LPK.T),  mode='r',check_finite=False)

    #from scipy.linalg import cholesky
    
    #L = cholesky( dot(LPK, LPK.T),lower=False, overwrite_a=False, check_finite=False)
    #U,S,V = svd(L, full_matrices=False, compute_uv=True, overwrite_a=True, check_finite=False)



    #estimate the rank
    #print s
    #plot(s)
    #show()  
    #print where(cumprod((diff(log(s[s!=0]))>-5)))
    #try:
    rank = where(cumprod((diff(log(s[s!=0]))>-5)|(s[s!=0][1:]>median(s))))[0][-1]+2
    #except:

    #print 'BUG!!!!!!'
    #rank-= 1
    #valid_ind = slice(None, rank+1)
    
    D = abs(s[:rank])

    U = u[:,:rank]

    vt_ = dot(U.T/D[:,None], LPK)
    
    
    
    #wrong_dets[~wrong_dets] |=~valid_ind
#print time.time()-t
    #t = time.time()

    #ordinary SVD (slow)
    #u, s_, vt_ = svd(LPK,  full_matrices=False,check_finite = False, overwrite_a=True)
    #print time.time()-t

    
    #valid_ind = isfinite(s)
    #valid_ind[valid_ind] &= s[valid_ind] > 0
    
    vt = zeros((vt_.shape[0],len(zero_ins)) )
    
    
    
    vt[:,~zero_ins] = vt_
    
    V = F.apply_Pt(F.solve_Lt(invsqrtD*vt.T)).T


    #valid_ind = slice(0, where(s/s.max() > 1e-7)[0][-1]+1)
    #valid_ind = slice(None,None)
  
    #sort from the most important to the least
    #valid_ind = D > 1e-6
    #D = D
    #U = (K*V)[:,valid_ind].T/D[:,newaxis]
    #U = u.T
  ##V = V[:,valid_ind]
    #V = V.T

    
    
    
    
    
    #V = V[:,valid_ind]
    #define decomposition in the same form as for the other methods
    #U =  (K*V).T
    #V = V.T
    
    
    #K_ = dot(dot(U,diag(1/D)),V)
    #print 'Decomposition accuracy',linalg.norm((K_* K.T)-eye(K.shape[0]))-len(s)+rank
    #exit()

    
    
    
    
    #U*D*V*B= K!!
    
    #iterative upgrade of the existing decompostion for a new B. But U is not orthogonal!? 
    #F = cholesky(B)
    #HK = asarray(F(K.T).T.todense())
    #V_ = dot(U.T/D[:,None],HK)
    #U_ = K*V_.T
    #D_ = D*linalg.norm(U,axis=0)


    return {'U':matrix(U), 'V':matrix(V) , 'D':D, 'wrong_dets':wrong_dets}

    


def PresolveSVD3D(K,B,nx,ny, ndets):
    """
    :param spmatrix K: Geometry matrix
    :param spamtrix B: Smoothing matrix D.T*D
    Solve SVD by cholesky and svd decomposition. 

    """

    wrong_dets = squeeze(array(K.sum(1)==0))
    
    
    
    K = K[where(~wrong_dets)[0],:]
    k = min(sum(~wrong_dets),ndets,K.shape[1] )
    

    #solve generalized eigen value problem K.T*K*x = lambda*C*x

    #from scikits.sparse.cholmod import cholesky, cholesky_AAt
    #from  scipy.linalg import eigh,svd,qr
    

    try:
        F = cholesky(B)
    except:
        #solve case whan B is singular
        B_trace = B.diagonal().sum()
        F = cholesky(B,B_trace*1e-16/nx/ny)
        
    
    LPK_ = F.solve_L(F.apply_P(K.T.tocsc()))  

    D = F.D()

    invsqrtD = sparse.spdiags(1/sqrt(D),(0,),len(D),len(D))
    LPK_ = (invsqrtD*LPK_).T
    #LPK_ = (invsqrtD*LPK_).T.toarray()  #LC^-1
    
    #calculate only for nonzero coloumns of LPK_ matrix
    zero_ins = array(LPK_.sum(0)==0).ravel()
    LPK = LPK_[:,where(~zero_ins)[0]]


    LPK.sort_indices()
    
    #LPKT = LPK.tocsr().T
    
    #LPKT[:k]*LPK[:k,0]
    
    
    
    from scipy.sparse import csc_matrix
    
    def slice_csc(A,k,typ='='):
        #slicing in scipy is extremaly slow :( 
        from scipy.sparse import csc_matrix
        ind = A.indptr[k]
        if typ == '=':
            ind2 = A.indptr[k+1]
            return csc_matrix((A.data[ind:ind2],A.indices[ind:ind2],[0, ind2-ind]),shape=(A.shape[0],1))
        if typ == '<':
            return csc_matrix((A.data[:ind],A.indices[:ind],A.indptr[:k+1]),shape=(A.shape[0],k))
        #if typ == '>':
            #return csc_matrix((A.data[ind:],A.indices[ind:],A.indptr[k:]-A.indptr[k]),shape=(A.shape[0],A.shape[1]-k))

    
    LPK = LPK.T.tocsc()
    LLs = sparse.vstack([ slice_csc(LPK,i+1,'<').T*slice_csc(LPK,i) for i in range(k)])
    #LLs = sparse.vstack([ LPK[:i+1,:]*LPK[i,:].T for i in arange(k)])
    #LLs = sparse.hstack([ LPK*LPK[i,:].T for i in arange(k)])

    
    LL = zeros((k,k))
    LL[tril_indices(k)] = array(LLs.todense()).flatten()
    #ind = ones((k,k), dtype='bool')
    LL.flat = array(LLs.todense()).flatten()
    LL2 = LPK.T*LPK
    
    
    #X = [r_[x.todense] for i, x in enumerate(X)]
    
    
    subplot(1,2,1);imshow(LL2.todense());subplot(1,2,2);imshow(LL);show()

    
    
    #LPK[:,:500]
    
    
    #LPK = LPK.toarray()
    
    #[for i in]

    
    #from scipy.sparse.linalg import svds

    #u, s, vt = svds(LPK,k-1)  #BUG reconstruct only M-1 eigenvalues!!  
    #ind = where(isfinite(s) & (s>1))[0][::-1] #remove noise vectors
    #s = s[ind]
    #vt = vt[ind,:]
    #ut = u[:,ind].T
    
    
    
    
    #NOTE ugly way how to calulate fast SVD
    
    #t = time.time()
    #u,s,v = svd(LPK, full_matrices=False, overwrite_a=True, check_finite=False)

    #r, = qr(LPK, overwrite_a=False,mode='r', pivoting=False, check_finite=False)
    #r = r[:,:k]
    #u,s,v = svd(r, full_matrices=False, overwrite_a=True, check_finite=False)
    #invR = dot(v.T, u.T/s[:,None] )
    #vt_ = dot(invR,LPK )
    #print time.time()-t
    #t = time.time()
    #tol = max(dim(m))*max(D)*.Machine$double.eps
#NOTE http://hackersome.com/p/arbenson/mrtsqr

    #T = ascontiguousarray(LPK.T)
    #dot(T.T, T)
    
    #from scipy.sparse.linalg import eigsh, svds
        #from scipy.sparse.linalg import svds
        


    #s,u = eigsh( LPK* LPK.T, k=k-1)
    #u, s, vt = svds(LPK,k-1)  #BUG reconstruct only M-1 eigenvalues!!  
    
    
    #from sparsesvd import sparsesvd
    #ut, s, vt = sparsesvd( LPK,k)

    #
    #imshow(LL.todense().astype('single'),interpolation='nearest')
    #show()
    

    
    from sparsesvd import sparsesvd
    #from scipy.sparse.linalg import eigsh, svds
    
    #FL = cholesky_AAt(LPK)
    #FL.apply_P(LPK)
    
    #LL = (LPK* LPK.T)
    #FLL = cholesky(LL)

    #I = eye(k)
    #I[:] = FL.apply_P(I)
    
    

    
    


    #461
    #t = time.time();eigh((LPK* LPK.T).todense(),overwrite_a=True, check_finite=False,lower=True);print time.time()-t
    #t = time.time();eigsh( LPK* LPK.T, k=k-1);print time.time()-t
    #t = time.time();svds(LPK,k-1);print time.time()-t
    #t = time.time();sparsesvd( LPK,k);print time.time()-t
    #t = time.time();svdp(LPK , k, which='L', irl_mode=False,blocksize = 16,compute_v=False);print time.time()-t


    
    
    #1120 = 160*7
    #LL[-160:,-160:]-LL[:160,:160]
    
    
    #from pypropack import svdp
    #u,s,v = svdp(LPK , k, which='L', irl_mode=False,blocksize = 16,compute_v=False )
    #imshow(LL.astype('float32'),aspect='auto', interpolation='nearest');show()
    
    
    #TODO Very ugly step!! i kan not slve it as sparse!!
    LL = (LPK* LPK.T).todense()
    s,u = eigh(LL,overwrite_a=True, check_finite=False,lower=True)
    del LL
    s,u = sqrt(abs(s))[::-1], u[:,::-1]
    #dot(u[:,:rank].T/D[:,None], LPK)
    #ut = u[:,ind].T
    #imshow(u.astype('float32'),aspect='auto', interpolation='nearest');show()
    #imshow(sqrt(abs(u.astype('float16'))),aspect='auto', interpolation='nearest');show()

    I = 1+dot(abs(u.T),arange(k))/amax(dot(abs(u.T),arange(k)))*1e-5
    ind = argsort(I*s)[::-1]
    u = u[:,ind]
    s = s[ind]
    u = sparse.csr_matrix(u)

    

    #imshow(sqrt(abs(u[:,ind].astype('float16'))),aspect='auto', interpolation='nearest');show()
    
    
    #u*arange(k)

    #r = qr(ascontiguousarray(LPK.T),  mode='r',check_finite=False)

    #from scipy.linalg import cholesky
    
    #L = cholesky( dot(LPK, LPK.T),lower=False, overwrite_a=False, check_finite=False)
    #U,S,V = svd(L, full_matrices=False, compute_uv=True, overwrite_a=True, check_finite=False)



    #estimate the rank
    #print s
    #plot(s)
    #show()  
    #print where(cumprod((diff(log(s[s!=0]))>-5)))
    #try:
    rank = where(cumprod((diff(log(s[s!=0]))>-5)|(s[s!=0][1:]>median(s))))[0][-1]+2
    #except:

    #print 'BUG!!!!!!'
    #rank-= 1
    #valid_ind = slice(None, rank+1)
    
    D = abs(s[:rank])

    U = u[:,:rank]
    
    #sparse.spdiags(1/D, 0, k,k)

    vt_ = (U.T*sparse.spdiags(1/D, 0, k,k))* LPK
    
    
    
    #wrong_dets[~wrong_dets] |=~valid_ind
#print time.time()-t
    #t = time.time()

    #ordinary SVD (slow)
    #u, s_, vt_ = svd(LPK,  full_matrices=False,check_finite = False, overwrite_a=True)
    #print time.time()-t

    
    #valid_ind = isfinite(s)
    #valid_ind[valid_ind] &= s[valid_ind] > 0
    
    #vt = zeros((vt_.shape[0],len(zero_ins)) )
    
    
    
    
    indptr = list(vt_.indptr)
    for i in where(zero_ins)[0]:  indptr.insert(i,indptr[i])    
    vt = sparse.csc_matrix((vt_.data, vt_.indices, indptr), shape=(rank, len(zero_ins)) )
    
    #vt.indptr = array(indptr)
    #vt.shape = (vt_.shape[0],len(zero_ins))
    #vt.set_shape((vt_.shape[0],len(zero_ins)))
    #vt[:,~zero_ins] = vt_.todense()   
    
    #sparse.csc_matrix(vt).indptr
    
    V = F.apply_Pt(F.solve_Lt(invsqrtD*vt.T)).T


    #valid_ind = slice(0, where(s/s.max() > 1e-7)[0][-1]+1)
    #valid_ind = slice(None,None)
  
    #sort from the most important to the least
    #valid_ind = D > 1e-6
    #D = D
    #U = (K*V)[:,valid_ind].T/D[:,newaxis]
    #U = u.T
  ##V = V[:,valid_ind]
    #V = V.T

        
 
    
    
    
    
    #V = V[:,valid_ind]
    #define decomposition in the same form as for the other methods
    #U =  (K*V).T
    #V = V.T
    
    
    #K_ = dot(dot(U,diag(1/D)),V)
    #print 'Decomposition accuracy',linalg.norm((K_* K.T)-eye(K.shape[0]))-len(s)+rank
    #exit()

    
    
    
    
    #U*D*V*B= K!!
    
    #iterative upgrade of the existing decompostion for a new B. But U is not orthogonal!? 
    #F = cholesky(B)
    #HK = asarray(F(K.T).T.todense())
    #V_ = dot(U.T/D[:,None],HK)
    #U_ = K*V_.T
    #D_ = D*linalg.norm(U,axis=0)


    return {'U':U, 'V':V , 'D':D, 'wrong_dets':wrong_dets}

    


def PresolveGEV(K,B,nx,ny, ndets,BdMat):
    """
    The Mathematics of some Tomography Algorithms used at JET  -- L.C.Ingesson


    :param spmatrix K: Geometry matrix
    :param spamtrix B: Smoothing matrix

    """

    wrong_dets = squeeze(array(K.sum(1)==0))
    #print K.shape, wrong_dets.shape
    K = K[where(~wrong_dets)[0],:]
    k = min(sum(~wrong_dets),ndets,K.shape[1] )
    
    C = K.T*K
    k = min(shape(C)[0]-1, ndets-1)

    #NOTE faster computation by cholmod (just sometimes faster!!, example - JET)
    try:
        from scipy.sparse.linalg.interface import  LinearOperator
        from scikits.sparse.cholmod import cholesky
        Minv = LinearOperator(shape(B), cholesky(B))
    except:
        debug( 'cholesky failured')
        Minv = None
    
    #slowest step
          
    D, V = eigsh(C,M = B, k=k,Minv= Minv,tol=0)  # in case that npix < ndets

    valid_ind = D > 1e-6
        
    #GEV solved, it can be used for computation 


    D = sqrt(D[valid_ind][::-1])
    V = V[:,valid_ind][:,::-1]

    #define decomposition in the same form as for the other methods
    U =  (K*V)/D[newaxis,:]
    V = V.T
    
    
    
    #NOTE alternative algorithm by An algorithm for quadratic optimization with one
    #quadratic constraint and bounds on the variables  #G C Fehmers†
    #BUt nemarically less stable and not rank revealing!!
    #fact = 1e4  #compensate a very low rank of C matrix => huge numerical instability 
    #D, V = eigsh(C,M = B*fact+C, k=k)  # in case that npix < ndets
    #D = D/(1-D)*fact
    #valid_ind = D > 0 #can be negative -  numerical error
    #D = sqrt(D[valid_ind])
    #U = (K*V)
    #V = V[:,valid_ind].T*D[:,None]
    #U = U[:,valid_ind] #reduce numerical errors? 
    
    

    
    
    
          
    #K_ = dot(dot(U.T,diag(1/D)),V)
    #print 'Decomposition accuracy',linalg.norm((K_* K.T)-eye(K.shape[0]))
    #exit()


    return {'U':matrix(U), 'V':matrix(V) , 'D':D, 'wrong_dets':wrong_dets}

    
def PresolveGEV_singular(K,B,nx,ny, ndets,BdMat):
    
    #BUG!!
    #dict =  PresolveGEV(K,B,nx,ny, ndets)
    #from scipy.sparse.linalg import eigsh
    
    
        
    wrong_dets = squeeze(array(K.sum(1)==0))
    K = K[where(~wrong_dets)[0],:]
    k = min(sum(~wrong_dets),ndets,K.shape[1] )
    
    C = K.T*K
    k = min(shape(C)[0]-1, ndets-1)

    fact = 1e4  #compensate a very low rank of C matrix => huge numerical instability 

    D, V = eigsh(C,M = B*fact+C, k=k)  # in case that npix < ndets
    D = D/(1-D)*fact
    
    
    valid_ind = D > 1e-6

    ##debug('done')

    D = sqrt(abs(D[valid_ind]))#can be negative -  numerical error

    U = (K*V).T#/sqrt(D)[:,newaxis]

    V = V[:,valid_ind].T*D[:,None]#(V*sqrt(D)[newaxis,:]).T
    U = U[valid_ind]#/linalg.norm(U, axis=1)[valid_ind,None] #reduce numerical errors? 
    D = D[::-1]
    U = U[::-1,:].T    #BUG sort the vectors from the most important to the less
    V = V[::-1,:]




    #K_ = dot(dot(U.T,diag(1/D)),V)
    #print 'Decomposition accuracy',linalg.norm((K_* K.T)-eye(K.shape[0]))


    
    return {'U':matrix(U), 'V':matrix(V) , 'D':D, 'wrong_dets':wrong_dets}




def PresolveGSVD(K,B,nx,ny, ndets):
    """
    http://www.netlib.org/lapack/explore-html/dd/db4/dggsvd_8f.html

    :param spmatrix K: Geometry matrix
    :param spamtrix B: Smoothing matrix
    
    """
    
    #dec = PresolveGEV(K,B,nx,ny, ndets)

    wrong_dets = squeeze(array(K.sum(1)==0))

    K = K[where(~wrong_dets)[0],:]
    k = min(sum(~wrong_dets),ndets,K.shape[1] )
    #from scikits.sparse.cholmod import cholesky

    try:
        F = cholesky(B)
    except:
        #solve case whan B is singular
        B_trace = B.diagonal().sum()
        F = cholesky(B,B_trace*1e-16/nx/ny)
        

    L,D = F.L_D()
    L = sqrt(D)*F.apply_Pt(L.tocsc()).T

  
    fact = 1e4  #compensate a very low rank of C matrix => huge numerical instability 

    U,V,X,C,S = gsvd(K.todense(),L.todense()*fact)  #normalizovat K a B 
    
    k = K.shape[0]
    C = diag(C[::-1,::-1])
    S = diag(S)[:-k-1:-1]
    X = inv(X, overwrite_a=True, check_finite=False)
    X = X.T[:,-k:][:,::-1]
    X*= fact/S
    V = X.T
    U = U[:,::-1]
    D = C/S*fact
    
    
    return {'U':matrix(U), 'V':matrix(V) , 'D':D, 'wrong_dets':wrong_dets}

  

def FastLinEmisivity(presolved_decomposition, S, norm,tsteps):
    """
    Algoritms that used preprocessed vectors, and estimate the emisivity
    $$Emis = \sum_{kl} E_{kl} = \sum_{kl} \sum_i \omega_i(\gamma) \frac{(u_i, S)}{D_i} V^{(i)}_{kl}$$
    $$ \omega = \frac1{1+n_{dets} 10^{g_{chi2}}/D^2}$$
    $$Emis = (b(\gamma), S) $$
    $$ b_i = \sum_j U_{ij} \left( \frac{\omega_j}{D_j}\sum_{kl}V^{(j)}_{kl}\right)$$

    :param int tsteps: Number of time steps
    :param dict presolved_decomposition: Dictionary with presolved vectors
    :param array S: measured data
    :var double g:  smoothing factor -- it can be estimated as median(D**2/ndets) or calculated.
    :param double norm:  data S normalization

    """
    wrong_dets = presolved_decomposition['wrong_dets']
    S = reshape(S, (-1,tsteps))
    
    
   
    S = S[~wrong_dets, :]
    U = presolved_decomposition['U']
    V = presolved_decomposition['V']
    D = presolved_decomposition['D']


    sum_v = sum(V,1)
    g = median(D**2)
    w = 1/(1+g/D**2)
    b = dot(sum_v*w/D,U.T)   #emision coefficients

    emis =  dot(b, S[:size(b),:])


    return multiply(tokamak,emis,norm)



def FindNegLocMin(image,lim,neighborhood_size=5):

    if all(image >= -lim):
        return []
    
    data_min = ndi.filters.minimum_filter(image, neighborhood_size)
    minima = (image == data_min)&(image < -lim)
    ind_minima = where(minima)
    active_set = ravel_multi_index(ind_minima,image.shape,order='F')

    return active_set

 
 
 
def ForcePositivity(E,W,g,V,D,prod,resid,ndets,BdMat,chi2,tokamak):
    
    from scipy.optimize import nnls
    
    #find positive solution in the subspace of solutions descrived by basis V

    #minimise |T*g-f|^2+|Dx|^2+  where x > 0
    #using nonlinear minimizer
    # not other method was working or it was too slow. 

    
    for ts, iE in enumerate(E):
    
        if tokamak.transform_index not in [0,4]:
            iE = tokamak.Transform*iE
            BdMat = False
            
     
        lim = iE[iE>0&~BdMat].mean()/3
        
        from scipy.optimize import minimize
        
        A = (W[ts]/D)*V.T/lim
        A[BdMat] = 0
        wdg = W[ts]/D*np.exp(g[ts]/2)
        args = W[ts], wdg, prod[ts], A
            
   
        def cost_fun( x, w, wdg, p, A):
            
            cost = np.sum((w*x-p)**2)+np.sum((x*wdg)**2)
            img = A.dot(x)
            neg = img < 0
            cost += np.sum(img[neg]**2)
   
            jac = 2*(w**2*x-p*w+wdg**2*x)
            jac += 2*dot(A[neg].T,img[neg])
    
            return cost, jac
        
        out = minimize(cost_fun, prod[ts],jac=True, method='L-BFGS-B', args=args )   
 
        iE = dot(W[ts]*out.x/D, V)
     
        #inverse transformation
        E[ts] = tokamak.Transform.T*iE 
        chi2[ts] =  (np.sum((prod[ts]-W[ts]*out.x)**2)+resid[ts])/ndets
        
        
    return E, chi2
        
        
        
        
def ForcePositivity2(E,W,g,V,D,prod,resid,ndets,BdMat,chi2,tokamak):
    
    from scipy.optimize import nnls
    
    #find positive solution in the subspace of solutions descrived by basis V

    #minimise |T*g-f|^2+|Dx|^2+  where x > 0

    #formulate the original problem in V basis
    #it is a trivial least squares problem which can be calculated analytically because 
    #A is double diagonal
     
    for ts, iE in enumerate(E):
    
        active_set = set()

        if tokamak.transform_index not in [0,4]:
            iE = tokamak.Transform*iE
            BdMat = False
            
        iE[BdMat] = 0
        lim = iE[iE>0].mean()/20


        D1 = np.diag(W[ts])
        D2 = np.diag(W[ts]/D*np.exp(g[ts]/2))
        A = np.vstack((D1,D2))
        y = np.hstack((prod[ts], 0*prod[ts]))
        p = prod[ts]
   
        #repeat in 10 interations
        for k in range(5):
            image = iE.reshape(ny,nx,order='F')
 
            update_active_set =  FindNegLocMin(image,lim)#[0]
            min_vals = iE[update_active_set]
            update_active_set = set(update_active_set)

            #if there is too many constrains, use only the most negative, dont use more than 50% of free dimensions
            nconst = len(update_active_set|active_set) 
            if nconst > len(prod[ts])//2:
                nconst = len(prod[ts])//2 - len(active_set)
                update_active_set = set(array(list(update_active_set))[argsort(min_vals)][:nconst])
     
            if len(update_active_set - active_set) == 0: 
                break  
            
            active_set|=update_active_set
            
        
            
            B = (W[ts] / D) * V.T[list(active_set)]

            U,S,V_ = np.linalg.svd(B)  #must be full SVD!
            
            n = len(active_set)
            pinvB = dot(U/S[None],V_[:n]).T
            perpB = V_[n:].T
            
            invB = np.hstack((perpB, pinvB))
            
            #rotate the problem in the space of negative points
            A_ = np.dot(A,invB)
            
            #https://en.wikipedia.org/wiki/Constrained_least_squares
            X1,X2 = A_[:,:-n], A_[:,-n:]

            #TODO better way to calculate np.linalg.pinv(X1)? 
            invX1 = np.linalg.pinv(X1)
            
            #ortonormal projector on space perpediculer to A1 - inequality constrained space
            P = np.eye(len(X1)) - np.dot(X1, invX1)
            
            
            #constrained projection in the rotated space
            beta2 = nnls(np.dot(P, X2),  np.dot(P,  y))[0]
            #beta2 = np.linalg.lstsq(np.dot(P, X2),  np.dot(P,  y))[0]
            
            #unconstrained projection in the rotated space
            beta1 = np.dot(invX1,  y - np.dot(X2, beta2))
            
            beta = np.hstack((beta1, beta2))
            
            
            #map back to the original projection space
            p = np.dot(invB, beta)
            
         
           # print(np.sum((prod[ts]-W[ts]*p)**2)/ np.sum(((1-W[ts])*prod[ts])**2), len(active_set))
            
            
            iE = dot(W[ts]*p/D, V)
            if tokamak.transform_index not in [0,4]:
                iE = tokamak.Transform*iE
            
            iE[BdMat] = 0
            continue
            
            #vmin  = iE.min()
            #vmax  = iE.max()
            #print( i, k, -sum(image[image<lim]) )
            #suptitle(str((i, k, -sum(image[image<lim]) )))
            #subplot(121)
            #itr = dot(single(w*delta_prod/D), V)
            #itr[BdMat] = 0
            #imshow(itr.reshape(ny,nx,order='F'),origin='lower', interpolation='nearest',vmin=vmin,vmax=vmax)
            #colorbar()
            ##figure()
            #subplot(122)
            #imshow(iE.reshape(ny,nx,order='F'),origin='lower', interpolation='nearest',vmin=vmin,vmax=vmax)
            #colorbar()
            #contour(image,0,colors='w',origin='lower')

            #show()
            
            image2 = copy(iE.reshape(ny,nx,order='F'))
            #corr = dot(single(w*delta_prod/D), V) #.5ms
            
            image2.ravel(order='F')[list(active_set)] = nan
            figure()
            #suptitle(str(i)+'  '+str(k))
            subplot(121)
            imshow( image2.reshape(ny,nx,order='F'),origin='lower', interpolation='nearest', vmin=iE.min(), vmax=None)
            colorbar()
            contour(image2,0,colors='w',origin='lower');
            subplot(122)
            imshow( image.reshape(ny,nx,order='F'),origin='lower', interpolation='nearest', vmin=iE.min(), vmax=None)
            colorbar()
            contour(image,0,colors='w',origin='lower');


            show()



            #imshow( corr.reshape(ny,nx,order='F'),origin='lower', interpolation='nearest',  vmin=iE.min(), vmax=iE.max())
            #colorbar()
            #contour(image,0,colors='w',origin='lower');
    

            

     
        #iE[(iE>-lim*2)&(iE<0)] = 0  #delete this almost zero poinsts

        #inverse transformation
        E[ts] = tokamak.Transform.T*iE 


        chi2[ts] =  (np.sum((prod[ts]-W[ts]*p)**2)+resid[ts])/ndets
        
    return E, chi2
        
            

         

def SolveLinearMetods(presolved_decomposition,tvec, S,ndets,dets, norm,L, H,BdMat,
                      error_scale,LastCycle, lam_method,positive_constrain,
                      g_fract_0,rapid_solver,plotting=False,estimate_sigma=True,lam_up=1,lam_low=0):
    """
    Algoritms that used preprocessed vectors, data and geometry matrix.

    Source: N Terasaki, Y Hosodab, M Teranishia and N Iwamaa

    :param dict presolved_decomposition: Dictionary with presolved vectors
    :param array S: measured data
    :param int tsteps: Number of time steps
    :param bool LastCycle: compute a GCV or chi2=1. In the other steps is used only guess
    :param bool full_reconstruction: calculate the emisivity profile

    """
    #t = time.time()
    
    
    #import IPython
    #IPython.embed()

    #TODO to kreslení udělat jako speciální třídu 
    global nx, ny, tokamak 
    #print 'LastCycle',LastCycle
    from os.path import expanduser
    #if LastCycle and debug: savez_compressed(expanduser('~')+'/tomography/H.npz',H=H)
    #print ' BUG!!!!!!!!   s'
    if g_fract_0 == None:
        g_fract_0 = .8   #default value 
        #use_chi2 = not use_gcv
    #else:  
        #use_chi2 = False

    t_total = time.time()

    wrong_dets = presolved_decomposition['wrong_dets']
    


    plotting=True
    #LastCycle = True
    plotting=False

    U = presolved_decomposition['U'].astype(dtype)
    V = presolved_decomposition['V'].astype(dtype)
    D = presolved_decomposition['D'].astype(dtype)
 
    #L_ = L[where(~wrong_dets)[0], :]
    #K_ = ((U*diag(1/D)*V))
    #print '\n\nDecomposition accuracy',linalg.norm((K_* L_.T)-eye(L_.shape[0]))#-len(S)+rank
    #exit()
    
    #PlotBaseVectors(V, nx,ny,Nmax=infty)

    if not all(fabs(D) > 0):
        print('undetected colinear dimensions!')
        ind = D == 0
        if not all(ind):D,V,U = D[ind],V[ind],U[:,ind]
        
        
        
        
    assert all(fabs(D) > 0), 'undetected colinear dimensions!'
    

    if 'R' in list(presolved_decomposition.keys()):   # QR decomposition
        R =  sparse.csc_matrix(presolved_decomposition['R'])
    else:
        R = 1               # it is universal for SVD and GEV
    
    tsteps = S.shape[1]
    S = S.astype(dtype)
    
 
    chi2 = zeros(tsteps)
    g = zeros(tsteps)
    doF = zeros(tsteps)
    retro = zeros((tsteps,ndets ),dtype=dtype)
    from scipy.stats.mstats import mquantiles


    g_min = 2*log(mquantiles(D,.1))-1
    g_max = 2*log(mquantiles(D,1))+1
    
    
    lam_low = min(g_fract_0, lam_low)
    lam_up  = max(g_fract_0, lam_up)

    g_lim_low,g0,g_lim_up = mquantiles(2*log(D),[lam_low,g_fract_0,lam_up])
    
    
    try:
        g0 = mquantiles(2*log(D), float(lam_method))[0]
        if float(lam_method) == 0: g0 = -50 
    except:
        pass  


    if plotting and LastCycle:

    
        from matplotlib.ticker import MaxNLocator,LogLocator
        g_tmp = logspace(g_min-3,g_max+3,100,base=e)
        g_quant = interp(log(g_tmp),2*log(D)[::-1],linspace(0,1,len(D)))
        #semilogx(g_tmp, g_quant)
        #[axvline(x,lw=.1,c='k') for x in D**2]
        #show()

        
        #ion()
        import matplotlib.pylab as plt
        #fig = plt.figure(figsize=(15,6))
        fig = plt.figure(figsize=(10,4))

        fig.clf()
        #fig.set_size_inches(10,6)

        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        try:
            G0 = load('./tmp/Emiss0.npz')['G']
            p_model_diff, = ax1.plot( g_quant , g_tmp,'k-',label='$||G-G_0||_2/||G_0||_2$')
            p_model_min = ax1.axvline(x=-100,lw=2,ls='--',c='k')
            #p_res_min   = ax1.axhline(y=1,lw=2,ls='--',c='b')
            #min_err_plot, = ax1.plot(g_quant, g_tmp*nan,'y--')
        except:
            #print 'no emiss 0'
            
            pass
        model_min_txt = ax1.text(0.05, 0.58 ,'', backgroundcolor='w', 
                                    transform=ax1.transAxes,color='k')
        
        methods = 'chi2', 'GCV','AICc','BIC','AIC', 'PRESS','Lcurve','Quasioptimum'
        methods = 'chi2', 'GCV','PRESS','AICc','AIC'

        if 'chi2' in methods: p_chi, = ax1.plot(g_quant,g_tmp,label='chi2')
        if 'GCV' in methods: p_gcv, = ax1.plot(g_quant,g_tmp,label='GCV')
        if 'chi2' in methods: p_chi_root, = ax1.plot( 1,1 ,'*b')
        if 'GCV' in methods: p_gcv_min, = ax1.plot(1,1,'go')
        if 'AICc' in methods: p_aicc_min, = ax1.plot(1,1,'kx')
        if 'AIC' in methods: p_aic_min, = ax1.plot(1,1,'vk')
        if 'BIC' in methods: p_bic_min, = ax1.plot(1,1,'dk')
        if 'PRESS' in methods: p_press, = ax1.plot(1,1,'rs')
        if 'Quasioptimum' in methods: Qopt_plt, = ax1.plot(g_quant, g_tmp*nan,'-.',label='Quasioptimum')
        if 'Lcurve' in methods: QLcurva_plt,= ax1.plot(g_quant, g_tmp*nan,'y--',label='Lcurve')
        if 'PRESS' in methods: PRESS_plt, = ax1.plot(g_quant, g_tmp*nan,':',label='PRESS')
        if 'AIC' in methods: AIC_plt, = ax1.plot(g_quant, g_tmp*nan,'k-.',label='AIC')
        if 'AICc' in methods: AICc_plt, = ax1.plot(g_quant, g_tmp*nan,'k--',label='AICc')
        if 'BIC' in methods: BIC_plt, = ax1.plot(g_quant, g_tmp*nan,'k:',label='BIC')
        if 'GCV' in methods: gcv_diff_txt  = ax1.text(0.05, 0.66, '',backgroundcolor='w',  transform=ax1.transAxes,color='k')
        if 'chi2' in methods: chi2_diff_txt = ax1.text(0.05, 0.76 , '',backgroundcolor='w',  transform=ax1.transAxes,color='k')
        if 'PRESS' in methods: press_diff_txt = ax1.text(0.05, 0.84 , '',backgroundcolor='w',  transform=ax1.transAxes,color='k')
        if 'AICc' in methods: aicc_diff_txt = ax1.text(0.05, 0.92 , '',backgroundcolor='w',  transform=ax1.transAxes,color='k')

        #p_eig_nums = [ax1.axvline(x=d,lw =.1,color='k') for d in D**2]
        
        if 'chi2' in methods: ax1.axhline(y=1,ls='--')
        variance_plot, = ax1.plot(g_quant, g_tmp*nan,'r--')
        bias_plot,     = ax1.plot(g_quant, g_tmp*nan,'g--')


        ax1.set_xlabel('$q_\lambda$',fontsize=13)
        ax1.set_ylabel('$\chi^2/N$',fontsize=13)
        leg = ax1.legend(loc='lower right',fancybox=True,ncol=2)
        leg.get_frame().set_alpha(0.9)
        #ax.xaxis.set_major_locator( MaxNLocator(nbins = 5) )
        #ax.locator_params(tight=True, nbins=4)
        ax1.xaxis.set_major_locator(LogLocator(numticks=6))
        #ax1.set_xlim(g_tmp[0],g_tmp[-1])
        ax1.set_xlim(0,1)
        ax1.set_yscale("log", nonposx='clip')
        ax1.set_xscale("linear", nonposx='clip')

        
        
        #ax1 = subplot(111)
        prod_ = array(dot(U.T, S[~wrong_dets, :]))            
        #plot(abs(prod_),'b',lw=.01)
        ax2.plot((mean(prod_**2,axis=1)),'k',label=r'$\langle |u_i b| \rangle$')
        vmax = mean(abs(prod_[:20,:]))
        #ax1.plot(D[0]*mean(abs(prod_[0,:]))/D,label='$1/D_i$')
        ax2.plot(median(D)/D*10,'--',label='$1/D_i$')
        orig_spectrum, = ax2.plot([],[],'k-',lw=.3)

  
        ax2.plot([],[],'r',label='$w_i$')
        ax2.axhline(y=1,c='k',ls='--')
        #yscale('symlog', linthreshy=10)
        ax2.set_xlabel('$i$-th singular value',fontsize=13)
        
        ax2.set_ylabel(r'${\langle (u_i b)^2 \rangle }$',fontsize=13)
        #ax2.set_ylabel(r'$U^T f$')

        ax2.set_yscale('log')
        ax2.set_ylim(.3, mean(amax(abs(prod_),0))**2*1.3)
        ax2.set_ylim(1.2e-1, mean(amax(abs(prod_),0))**2*1.3)
        ax2.set_ylim(1.2e-1, mean(linalg.norm(S,axis=0)**2))

        #ax2.set_ylim(.1, 1e6)

        ax2.text(20,0.7,'noise level')
        leg = ax2.legend(loc='upper right',fancybox=True)
        leg.get_frame().set_alpha(0.9)
        ax3 = ax2.twinx()
        ax3.set_yscale('linear')
        p_filter, = ax3.plot([],'-.r')
        ax3.set_ylim(0,1)
        ax2.set_xlim(0,len(D))
        ax3.set_xlim(0,len(D))

        ax3.set_ylabel(r'$w_i$',fontsize=13)
        for l in leg.legendHandles:  l.set_linewidth(3)
        #ioff()
        ##fig.show()
        #show() 
        #pause(1)
        
    def w_i(g):
        w = 1./(1.+exp(g)/D**2)
        w[~isfinite(w)] = 0
        return w

        
    def Press(g, prod,resid=0,G0=None):
        #TODO how will be PRESS effected by nonzero resid? 
        
        
        w = w_i(g)
        u = array(U,copy=False)

        ind = einsum('ij,ij->i',u,u) < single(u.shape[1]*.1/u.shape[0])  #  diag(U*U.T) - smaller are wrong, linearly dependent? 

        if ind.any(): u = u[~ind]
        return sum((dot(u, (1-w)*prod)/einsum('ij,ij,j->i', u,u, 1-w))**2)/ndets
    
        #einsum('ij,ij,j->i', u,u, 1-w) it is diag(W*U*U.T) - 1-leverage- sensitivity of retrofit on these measurements
    
        #sum(U**2,1)
        
        
        

    def CurvLCurv(g, prod,resid=0):
        #curvature of the Lcurve
        w = w_i(g)
        g = exp(g)        

        a = linalg.norm(R*((w-1)*prod))**2+resid**2
        b = linalg.norm((w/D*prod))**2
        
        
        
        E = sum((D*prod)**2/(D**2+g)**3)

        kappa = (a*b/E-g*(a+g*b))/(a+g**2*b)**1.5
        
        return kappa

        
    def QuasiOpt(g, prod,resid=0):
        #Quaisi-Optimum algorithm  PC Hansen Matlab Regularization Tools
        w = w_i(g)
        return 1/(linalg.norm(prod/D*w*(1-w))*ndets) #BUG inversion is there to have a minimum!


    def GCV(g, prod,resid=0):
        #generalized crossvalidation
        w = w_i(g)
        #return 1/(1-mean(w))

        return (sum((R*((w-1)*prod))**2)+resid)/ndets/(1-mean(w))**2
    
    
    def AIC(g, prod,resid=0):
        #AIC criteriom

        w = w_i(g)
        n = size(prod)
        #k = len(w)-sum((1-w)**2)  #number of "parameters" 
        k = sum(w)
        chi2 = CHI2(g,prod,resid,bias=True,log_out=False)
        LnL = -.5*chi2*ndets-n*log(2*pi)/2
        
        return 2*k-2*LnL
        
        
        #ln_chi2 = CHI2(g,prod,resid,bias=True,log_out=True)
        #ln_chi2+= log(ndets)
        
        
        

        #return ndets*ln_chi2+2*k
    
    def AICc(g, prod,resid=0):
        #AICc is AIC with a correction for finite sample sizes.
        w = w_i(g)
        #k = len(w)-sum((1-w)**2)  #number of "parameters" 
        k = sum(w)
        n = size(prod)
        chi2 = CHI2(g,prod,resid,bias=True,log_out=False)
        LnL = -.5*chi2*ndets-n*log(2*pi)/2

        return -2*LnL+2*(k+1)*n/(n-k-2)
        #return -2*LnL+2*k*n/(n-k-1)

        
        #n = size(prod)
        #w = w_i(g)
        #k = len(w)-sum((1-w)**2)  #number of "parameters" 
        #ln_chi2 = CHI2(g,prod,resid,bias=True,log_out=True)
        
        ##BUG output can be negative!!
        #return n*ln_chi2+2*k*(k+1)/(n-k-1)
    
    
    def BIC(g, prod,resid=0):
        #BIC criteriom
        
        
        
        w = w_i(g)
        n = size(prod)
        k = len(w)-sum((1-w)**2)  #number of "parameters" 
        chi2 = CHI2(g,prod,resid,bias=True,log_out=False)
        LnL = -.5*chi2*ndets#-n*log(2*pi)/2
        
        return k*log(n)-2*LnL
    
    
        #n = ndets
        #w = w_i(g)
        #k = len(w)-sum((1-w)**2)  #nomber of "parameters" 
        #ln_chi2 = CHI2(g,prod,resid,bias=True,log_out=True)
        #ln_chi2+= log(ndets)

        #return n*ln_chi2+k*log(n)  #BUG zkontroolovat
    
    
    def CHI2(g,prod,resid,bias=True,log_out=False):
        #  discrepancy principle  - chi2
        #ndets = len(D) #BUG!!
        #w = 1/(1+ndets*exp(g)/D**2)
        w = w_i(g)
        #nDoF = ndets if bias else ndets-sum(w)
        nDoF = ndets  if bias else sum((1-w)**2)
        #nDoF = len(D)
        #resid = 0 #BUG  ale pak ne vždy existuje chi2=0
        
        if  g < g_min-5:
            if log_out:
                return log(resid+ndets*sum((R*(prod/D**2))**2))+2*g
            else:
                return ndets*sum((R*(prod/D**2))**2)*exp(2*g)+resid
   
        if log_out:
            return log((linalg.norm(R*((w-1)*prod))**2+resid)/nDoF)
        else:
            return (linalg.norm(R*((w-1)*prod))**2+resid)/nDoF

    def FindMin(F, x0,dx0,prod,resid,tol=0.01):
        #stupid but robust minimum searching algorithm.

        fg = F(x0, prod,resid)

        while abs(dx0) > tol:
            fg2 = F(x0+dx0, prod,resid)
                                
            if fg2 < fg:
                fg = fg2
                x0 += dx0
                continue
            else:
                dx0/=-2.
                
        return x0, log(fg2)
    
    
        #g_gcv  = copy(g0)
        
        #dx0 = 1.0
        #asymptote = log(ndets*sum((R*(prod[ts]/D**2))**2))#+2*g_gcv

        #fg = log(GCV(g_gcv, prod[ts],resid[ts]))

        #while abs(dg) > 0.01:
            ##print dg, exp(g_gcv)
            #fg2 = log(GCV(g_gcv+dg, prod[ts],resid[ts]))
                                
            #if fg2 < fg:# and fg2 < asymptote+2*(g_gcv+dg):
                #fg = fg2
                #g_gcv += dg
                #continue
            #else:
                #dg/=-2

    #clf()
    #plot(w_i(g0))
    #show()
   
   
    #if LastCycle:
            
        
        
    gcv = 0
    g_gcv = 0
    
    S_ = S[~wrong_dets, :]  # remove not important  LOS 
    L_ = L[where(~wrong_dets)[0], :]

    #try:
        #prod = array(dot(S_.T,U))
    prod = array(U.T.dot(S_).T)

    #except:
    

    if tsteps > U.shape[1]:
        P = eye(U.shape[0])-(R*U.T).T*(R*U.T)
        resid = linalg.norm(dot(P.astype(single),S_) ,axis=0)**2 
    else:
        #print 
        
        
        resid = linalg.norm(S_-U*(transpose(R)*(R*prod.T)),axis=0)**2 
        #resid = linalg.norm(dot(P.astype(single),S_) ,axis=0)**2 

        
    resid += linalg.norm(S[wrong_dets,:],axis=0)**2 
    #resid[:] = 0


 
    for ts in range(tsteps):
        if (not LastCycle and not lam_method=='chi2') or not lam_method\
            in ['chi2','gcv','press','aic','aicc','bic','qoptimum'] :
            g[ts] = g0   

        
        if (lam_method=='gcv' and LastCycle) or plotting :
            #GCV  minimum finding algorithm 
            g_gcv, fgcv = FindMin(GCV, g0, 1,prod[ts],resid[ts])
            g[ts]  = g_gcv


        if (lam_method=='press' and LastCycle) or plotting :
            #PRESS  predicted sum of squares, a little slower than GCV
            g_press, fpress = FindMin(Press, g0,1,prod[ts],resid[ts])
            g[ts]  = g_press     
              
              
        if (lam_method=='aic' and LastCycle) or plotting :
            #AIC
            g_aic, faic = FindMin(AIC, g0,1,prod[ts],resid[ts])
            g[ts]  = g_aic

            
        if (lam_method=='aicc' and LastCycle) or plotting :
            #AICc
            g_aicc, faicc = FindMin(AICc, g0,1,prod[ts],resid[ts])
            g[ts]  = g_aicc
            

        if (lam_method=='bic' and LastCycle) or plotting:
            #BIC 
            g_bic, fbic = FindMin(BIC, g0,1,prod[ts],resid[ts])
            g[ts]  = g_bic
            
        if (lam_method=='qoptimum' and LastCycle) or plotting:
            #quasioptimum
            g_qopt, fqopt = FindMin(QuasiOpt, g0,1,prod[ts],resid[ts])
            g[ts]  = g_qopt 


            
        #g_tmp = logspace(g_min-10,g_max+10,100,base=e)
        #chi2_tmp = zeros(100)
        #chi2_tmp2 = zeros(100)
        #GCV_tmp = zeros(100)
        
        #for i in range(100):
                #w = 1/(1+g_tmp[i]/D**2)
                #chi2_tmp2[i] = CHI2(log(g_tmp[i]),prod,resid,False)
                #chi2_tmp[i] = CHI2(log(g_tmp[i]),prod,resid)
                #GCV_tmp[i]  = GCV( log(g_tmp[i]),prod,resid)
                
        #loglog(g_tmp,chi2_tmp)
        #loglog(g_tmp,chi2_tmp2)

        #loglog(g_tmp,GCV_tmp)
        #axhline(y=1)
        #axvline(exp(g_gcv))
        #show()
        
        
        if (lam_method=='chi2' and LastCycle) or plotting:

            if sum((R*prod[ts])**2)/ndets < 1+(.1) or resid[ts]/ndets > 1 :  #no root!
                g_chi2, g[ts], f_chi2 = g0,g0, CHI2(g0,prod[ts],resid[ts])
            else:   
                    g1 = g0 - .2
                    g2 = g0 + .2
                    # derivaci lze vyjádřit analyticky!
                    fg1 = CHI2(g1,prod[ts],resid[ts],log_out=True)
                    fg2 = CHI2(g2,prod[ts],resid[ts],log_out=True)
                    
                    # very simple searching for chi2 = 1
                    steps = 0
                    while abs(fg2)>0.01  and steps < 100:
                        steps+= 1
                        dg = -fg2*(g2-g1)/(fg2-fg1)
                        g1 = g2
                        g2 += dg
                        fg1 = fg2
                        #if g2 < g_min:
                            #fg2 = log(ndets*sum((R*(prod/D**2))**2))+2*g2
                            
                        if g2 > g_max:
                            g2 = g_max
                            fg2 = CHI2(g_max,prod[ts],resid[ts],log_out=True)
                        else:
                            fg2 = CHI2(g2,prod[ts],resid[ts],log_out=True)

                    if steps==100:
                        print('no root find')
                        #g2, fg2 = g0, CHI2(g0,prod[ts],resid[ts],log_out=True)
                        
                    g[ts] = g2
                    g_chi2 = g2
                    f_chi2 = fg2
                    #chi2[ts] = exp(fg2)
 
     
            
   
        if plotting and LastCycle:
            ##plotting
            #ax2 = subplot(111)

            chi2_tmp = zeros(100)
            chi2_tmp2 = zeros(100)
            Q_tmp = zeros(100)
            L_tmp = zeros(100)
            AIC_tmp = zeros(100)
            AICc_tmp = zeros(100)
            BIC_tmp = zeros(100)

            PRESS = zeros(100)
            stat_err = zeros(100)
            bias_err = zeros(100)

            GCV_tmp = zeros(100)
            model_diff = zeros(100)
            model_diff_pos = zeros(100)
            #quasi_optim = zeros(100)
            #from the phantom analysis it looks like that gcv is much more stable and accurate guess of chi2=0
            Lcurve = zeros(100)

            G0 = None
            f0 = None
            g_optim = None
            min_res = 0
            try:
                tvec_0 = load('./tmp/Emiss0.npz')['tvec']
                G0 = load('./tmp/Emiss0.npz')['G'].reshape( L.shape[1], -1,order='F')
                ts0 = argmin(abs(tvec_0-tvec[ts]))
                G0 = G0.mean(1)
                for i in range(100):
                    w = w_i(log(g_tmp[i]))
                    model_diff[i]      = linalg.norm(G0-dot(single(w*prod[ts]/D)*norm[ts], V))/linalg.norm(G0)
                                        #(linalg.norm(G0-dot(single(w*prod[ts]/D)*norm[ts], V))/linalg.norm(G0[:,ts])*100)
                    pos = lambda x: (abs(x)+x)/2
                    #positive part
                    model_diff_pos[i]  = linalg.norm( G0-pos(dot(single(w*prod[ts]/D)*norm[ts], V)))/linalg.norm(G0)
                    prod0 = dot((L*G0)[~wrong_dets] ,U).T
                    PP = abs(prod0.T)/norm[ts]/D

                
                imin = argmin(model_diff_pos)
                def poly_min(x,dx,y): #find extrem of polynom 
                    return x+(y[0]-y[2])*dx/(2*(y[0]+y[2]-2*y[1]))
                if imin > 0:
                    g_optim = exp(poly_min(log(g_tmp[imin]),log(g_tmp[imin+1]/g_tmp[imin]),model_diff_pos[imin-1:imin+2]))
                else:
                    g_optim = 1e-100
                #g_optim = g_tmp[argmin(model_diff)]
                #title('min phantom difference %.2f'%(model_diff.min()*100))
                
                
                
                g_optim_ = interp(log(g_optim),2*log(D)[::-1],linspace(0,1,len(D)))
                print('g_optim_rms', g_optim_)
                p_model_min.set_data([g_optim_,g_optim_],[1e-5,1e5])
                p_model_diff.set_ydata(model_diff_pos)

                #w = 1/(1+g_optim/D**2)
                w = w_i(log(g_optim))
                model_min_txt.set_text('$\min\ \Delta G:$ %.1f%%'%(linalg.norm(G0-\
                    dot(single(w*prod[ts]/D)*norm[ts],V))/linalg.norm(G0)*100))
                
                #ax.text(0.05, 0.7 , 'MIN dev %.2f%%'%(linalg.norm(G0[:,ts]-\
                    #dot(single(w*prod[ts]/D)*norm[ts],V))/linalg.norm(G0[:,ts])*100),
                    #backgroundcolor='w',  transform=ax.transAxes,color='k')
                
                #if size(V) > 1e6:
                    #iV = V.T  #BUG it is very slow to calculate a pinv!!
                #else:
                iV = linalg.pinv(V)
           
                min_res = linalg.norm(G0-dot(V.T,dot(iV.T,G0).T).T)/linalg.norm(G0)
                print('min residuum', min_res)
                V_bias = G0-dot(V.T,dot(iV.T,G0).T).T
                orig_spectrum.set_data(arange(len(w)),(array(dot(U.T, L_*G0.T)).T)**2/norm[ts]**2)

                #min_err_plot.set_ydata(min_res)
    
            except:
                #print e

                
                
                pass
                #pass
               #savez_compressed('./tmp/Emiss0',G = reshape(G0,(tokamak.ny,tokamak.nx,tsteps), order='F')
                        #,Rc= tokamak.xgrid+tokamak.dx/2,Zc= tokamak.ygrid+tokamak.dy/2)
            #close()             
            
            #show()
            
            
            #subplot(131)
            #imshow(((G0[:,ts]-dot(single(w*prod/D)*norm[ts], V))/G0[:,ts]).reshape(tokamak.ny, tokamak.nx,order='F'),origin='lower',vmin=-.1, vmax = .1)
            #colorbar()
            #subplot(132)
            #imshow(G0[:,ts].reshape(tokamak.ny, tokamak.nx,order='F'),origin='lower',cmap='Paired')
            #colorbar()
            #subplot(133)
            #imshow(dot(single(w*prod/D)*norm[ts], V).reshape(tokamak.ny, tokamak.nx,order='F'),origin='lower',cmap='Paired')
            #colorbar()
            #show()
            #print tsteps
            
            
            
            #PRESS[i]  =    Press( g_gcv,prod[ts],resid[ts],G0)
            G0_ = G0/norm[ts] if not G0 is None else dot(single(w_i(g[ts])*prod[ts]/D),V)
            t = time.time()
            for i in range(100):
                #w = 1/(1+g_tmp[i]/D**2)
                #chi2_tmp[i] = (resid+sum((R*((w-1)*prod))**2))/ndets
                chi2_tmp[i] = CHI2(log(g_tmp[i]),prod[ts],resid[ts])
                GCV_tmp[i]  = GCV( log(g_tmp[i]),prod[ts],resid[ts])
                Q_tmp[i]  = QuasiOpt( log(g_tmp[i]),prod[ts],resid[ts])
                PRESS[i]  = Press( log(g_tmp[i]),prod[ts],resid[ts])
                L_tmp[i]  = CurvLCurv( log(g_tmp[i]),prod[ts],resid[ts])
                AIC_tmp[i]  = AIC( log(g_tmp[i]),prod[ts],resid[ts])
                AICc_tmp[i] = AICc( log(g_tmp[i]),prod[ts],resid[ts])
                BIC_tmp[i]  = BIC( log(g_tmp[i]),prod[ts],resid[ts])
                chi2_tmp2[i]= CHI2(log(g_tmp[i]),prod[ts],resid[ts]*0,bias=False)
                #w = 1/(1+g_tmp[i]/D**2)
                w = w_i(log(g_tmp[i]))
         
                #Lcurve[i] = linalg.norm(H*dot(single(w*prod[ts]/D), V))**2
                Lcurve[i] = sum(w**2/D**2*prod**2)
                #len(w)/float(ndets)*
                #close()
                #plot(G0_.ravel());show()
                
                linalg.norm(dot(single(w_i(g[ts])*prod[ts]/D),V))
                stat_err[i] = sqrt(sum(dot((w/D)**2,asarray(V)**2)))/linalg.norm(G0_)
                
                
                 
                
                
                #bias_err[i] = linalg.norm(V_bias+dot((1-w)/D*array(L_*G0.T).T,V))/linalg.norm(G0)
                #stat_err[i] = sqrt(sum((V/D[:,None])**2))/linalg.norm(G0_)
                try:
                    bias_err[i] = linalg.norm(V_bias+dot(  (1-w)/D*array(dot(U.T, L_*G0.T))[0]  ,V))/linalg.norm(G0)
                    #bias_err[i] = hypot(bias_err[i],stat_err[i])
                except:pass
                #plt.plot(sqrt(mean(prod_**2,axis=1)),'k',label=r'$\langle |u_i b| \rangle$')
                #plt.plot(abs(array(dot(U.T, L_*G0_.T)).T))
                #plt.show()
                
                
                
                ##X = dot(w*(random.randn(len(D),1000)+prod[ts][:,None] ).T/D, V)
                    
                #stat_err2[i]  = linalg.norm(X.std(0))/linalg.norm(G0_)
                
                #sqrt(sum(dot(( w_i(g_gcv )/D)**2,V**2)))/linalg.norm(G0_)

                #quasi_optim[i] = linalg.norm(dot(single(w**2*len(D)*g_tmp[i]/D**3*prod), V))
#imshow(X.mean(0).reshape(60,40,order='F'),interpolation='nearest' ,origin='lower');colorbar()
#imshow(X.std(0).reshape(60,40,order='F'),interpolation='nearest' ,origin='lower');colorbar()
            print(time.time()-t)



            #title(str(ts)+' '+str( LastCycle))
            #close()
            da = gradient(log(chi2_tmp))/gradient(log(g_tmp))
            dda = gradient(da)/gradient(g_tmp)
            db = gradient(log(Lcurve))/gradient(log(g_tmp))
            ddb = gradient(db)/gradient(g_tmp)
            kappa = (da*ddb-dda*db)/(da**2+db**2)**1.5
            #close()
            
             
            
            
            
            #close()
            #loglog(chi2_tmp ,Lcurve)
            #figure()
            #semilogx(chi2_tmp,kappa)
                
            #show()
            
            #from scipy.stats import chi2
            #print chi2.ppf(0.10, len(D)),  chi2.ppf(0.90, len(D)),chi2.ppf(0.50, len(D)),len(D)
            #ax2.loglog(chi2_tmp, Lcurve,'b-',label='L-curve')
            #ax2.loglog(chi2_tmp2, Lcurve,'r--',label='unbiased $\chi^2/DoF$')

            #ax2.axvline(x=1,ls=':',c='k')
            ##ax2.axvline(x=1)
            

            ##ax2.plot(interp(D**2, g_tmp,chi2_tmp),interp(D**2, g_tmp, Lcurve),'y+')
            #ax2.plot(interp(exp(g_gcv), g_tmp,chi2_tmp),interp(exp(g_gcv), g_tmp, Lcurve),'o',label='GCV')
            #ax2.loglog(1,interp(1, chi2_tmp,Lcurve),'sk',label=r'$\chi^2/k = 1$')

            #ax2.set_xlabel('$||Ax-b ||_2^2/DoF$')
            #ax2.set_ylabel('$||Dx||_2^2$')
            #ax2.xaxis.set_major_locator(LogLocator(numticks=6))
            #ax2.yaxis.set_major_locator(LogLocator(numticks=6))
            #ax2.legend(loc='best')
            #if not g_optim is None:
                #ax2.plot(interp(g_optim, g_tmp,chi2_tmp),interp(g_optim, g_tmp, Lcurve),'*',label='optimal')

            #ind = argsort()
            #ax2.axvline(x=interp(1,chi2_tmp2[::-1],chi2_tmp[::-1])  ,ls='--')

            #ax2.axvline(x=interp(exp(g_gcv), g_tmp,chi2_tmp) )

            #show()
            

            #loglog(g_tmp,chi2_tmp)
            #loglog(g_tmp,chi2_tmp2,'-.')
            #loglog(g_tmp,GCV_tmp)
            
            #loglog(g_tmp,quasi_optim,':')

            #plot(exp(g2),exp(fg2) ,'*b',label = '$\chi^2/N_{det} = 1$')
            #plot(exp(g_gcv),exp(gcv),'go', label = 'min GCV')
            #plot(exp(g0) ,CHI2(g0,prod[ts],resid[ts]),'rs',label = 'initial guess')
            #plot(exp(median(log(D**2))) ,CHI2(median(log(D**2)),prod,resid),
                 #'ks',label = 'median $\gamma_0$')

            if 'chi2' in methods: p_chi.set_ydata(chi2_tmp)
            #p_chi2.set_ydata(chi2_tmp2)
            if 'GCV' in methods: p_gcv.set_ydata(GCV_tmp)
            #ax1.plot(g_tmp, Q_tmp,'-.')
            #ax1.plot(g_tmp, kappa,'--')
            #ax1.plot(g_tmp, PRESS,':')
            if 'Quasioptimum' in methods:Qopt_plt.set_ydata(Q_tmp)
            if 'Lcurve' in methods:QLcurva_plt.set_ydata(kappa)
            if 'AIC' in methods: AIC_plt.set_ydata(AIC_tmp/100)
            if 'AICc' in methods: AICc_plt.set_ydata(AICc_tmp/100)
            if 'BIC' in methods: BIC_plt.set_ydata(BIC_tmp/100)
            if 'PRESS' in methods: PRESS_plt.set_ydata(PRESS)

            #ax1.plot(g_tmp, AIC_tmp,'k-.',label='AIC')
            #ax1.plot(g_tmp, AICc_tmp,'k--',label='AICc')
            #ax1.plot(g_tmp, BIC_tmp,'k:',label='BIC')

            
            
            
            

                    #Qopt_plt, = ax1.plot(g_tmp, g_tmp,'-.')
        #QLcurva_plt,= ax1.plot(g_tmp, g_tmp,'--')
        #PRESS_plt, = ax1.plot(g_tmp, g_tmp,':')
            
            #p_res_min   = ax1.axhline(y=1,lw=2,ls='--',c='b')
            variance_plot.set_ydata(stat_err)

            try:
                bias_plot.set_ydata(bias_err)
            except:
                pass

            
            #ax1.axhline(y= linalg.norm(dot(( w_i(g_gcv)/D)**2,V**2)) ,ls='--',c='r')
            
            
            #w_i(g_gcv)/D*X
            
            #dot( w_i(g_gcv)/D*X.T,V)
            
            
            #err = std(dot( w_i(g_gcv)/D*X.T,V),0)**2
            
            
            #plot(sqrt(dot(( w_i(g_gcv)/D)**2,V**2)))
            
            #plot(dot(single(w_i(g_gcv)*prod[ts]/D), V))
            


            #ax1.plot(g_tmp, L_tmp,'--')
        #g_quant = 
        
        #interp(g_chi2,2*log(D)[::-1],linspace(0,1,len(D)))
            if 'chi2' in methods: p_chi_root.set_data(interp(g_chi2,2*log(D)[::-1],linspace(0,1,len(D))),exp(f_chi2))
            if 'GCV' in methods: p_gcv_min.set_data(interp(g_gcv,2*log(D)[::-1],linspace(0,1,len(D))),exp(fgcv))
            if 'PRESS' in methods: p_press.set_data(interp(g_press,2*log(D)[::-1],linspace(0,1,len(D))),exp(fpress))
            if 'AICc' in methods: p_aicc_min.set_data(interp(g_aicc,2*log(D)[::-1],linspace(0,1,len(D))),exp(faicc)/100)
            if 'AIC' in methods: p_aic_min.set_data(interp(g_aic,2*log(D)[::-1],linspace(0,1,len(D))),exp(faic)/100)
            if 'BIC' in methods: p_bic_min.set_data(interp(g_bic,2*log(D)[::-1],linspace(0,1,len(D))),exp(fbic)/100)

            #p_chi, = loglog([],[])
            #p_chi2, = loglog([],[],'-.')
            #p_gcv, = loglog([],[])
            
            #p_chi_root, = plot([],[] ,'*b',label = '$\chi^2/N_{det} = 1$')
            #p_gcv_min, = plot([],[],'go', label = 'min GCV')
            #p_init, = plot([],[],'rs',label = 'initial guess') 
      
            if not G0 is None:
                #w = 1/(1*exp(g_gcv)/D**2)
                w = w_i(g_gcv)

                if 'GCV' in methods: gcv_diff_txt.set_text('$\mathrm{GCV}: \Delta G $ %.1f%%'%(linalg.norm(G0\
                        -dot(single(w*prod[ts]/D)*norm[ts],V))/linalg.norm(G0)*100))
                #w = 1/(1+ndets*exp(g2)/D**2)
                w = w_i(g2)

                if 'chi2' in methods: chi2_diff_txt.set_text('$\chi^2/N: \Delta G $ %.1f%%'%(linalg.norm(G0
                    -dot(single(w*prod[ts]/D)*norm[ts],V))/linalg.norm(G0)*100))
                
                w = w_i(g_press)

                if 'PRESS' in methods: press_diff_txt.set_text('$\mathrm{PRESS}: \Delta G $ %.1f%%'%(linalg.norm(G0
                    -dot(single(w*prod[ts]/D)*norm[ts],V))/linalg.norm(G0)*100))
            
                
                w = w_i(g_aicc)

                if 'AICc' in methods: aicc_diff_txt.set_text('$\mathrm{AICc}: \Delta G $ %.1f%%'%(linalg.norm(G0
                    -dot(single(w*prod[ts]/D)*norm[ts],V))/linalg.norm(G0)*100))
            
            
                
                
                
                #ax.text(0.05, 0.75 , 'GCV dev %.2f%%'%(linalg.norm(G0[:,ts]-dot(single(w*prod[ts]/D)*norm[ts],                                                                        
                                    #V))/linalg.norm(G0[:,ts])*100),backgroundcolor='w',  transform=ax.transAxes,color='k')
                #ax.text(0.05, 0.8, '$\chi^2/doF$ dev %.2f%%'%(linalg.norm(G0[:,ts]-dot(single(w*prod[ts]/D)*norm[ts],                                                                           
                                    #V))/linalg.norm(G0[:,ts])*100), backgroundcolor='w',  transform=ax.transAxes,color='k')
            
            ax1.set_ylim(1e-3,sum(prod[ts]**2)/ndets*2)

            #p_eig_nums
            #for d in D**2/ndets:   axvline(x=d, lw = .1,color='k')
            #asymptotes
            #loglog(g_tmp,ndets*sum((R*(prod/D**2))**2)*g_tmp**2,'-.k')
            #axhline(y=sum(prod**2)/ndets, linestyle='-.', color='k')
    
            
    

            #savefig('%d.png'%ts)
            
            #show()
            #figure()
            
            

            #ax1 = subplot(111)
            #prod_ = dot(U.T, S)            
            ##plot(abs(prod_),'b',lw=.01)
            #ax1.plot(sqrt(mean(prod_**2,axis=1)),'k',label=r'$\langle |u_i b| \rangle$')
            #vmax = mean(abs(prod_[:20,:]))
            ##ax1.plot(D[0]*mean(abs(prod_[0,:]))/D,label='$1/D_i$')
            #ax1.plot(median(D)/D,label='$1/D_i$')

            #ax1.plot([],[],'r',label='$w_i$')
            #ax1.axhline(y=1,c='k',ls='--')
            ##yscale('symlog', linthreshy=10)
            #ax1.set_xlabel('$i$-th singular value')
            #ax1.set_ylabel(r'$\langle |u_i b| \rangle $')
            #ax1.set_yscale('log')
            #ax1.set_ylim(.1, mean(amax(abs(prod_),0))*1.3)
            #ax1.text(2,1.1,'noise level')
            #leg = ax1.legend(loc='upper right',fancybox=True)
            #leg.get_frame().set_alpha(0.9)
            #ax2 = ax1.twinx()
            #ax2.set_yscale('linear')
            #ax2.plot(1/(1+ndets*exp(g_gcv)/D**2),'r')
            #ax2.set_ylim(0,1)
            #ax2.set_ylabel(r'$w_i$')
            #for l in leg.legendHandles:  l.set_linewidth(3)

            #show() 
            g_ = g0
            if lam_method=='press':         g_ = g_press
            if lam_method=='gcv':           g_ = g_gcv
            if lam_method=='chi2':          g_ = g_chi2
            if lam_method=='aic':           g_ = g_aic
            if lam_method=='aicc':          g_ = g_aicc
            if lam_method=='bic':           g_ = g_bic
            if lam_method=='qoptimum':      g_ = g_qopt
            
            p_filter.set_data(arange(len(D)),w_i(g_))
            #ion()
            #savefig('./tmp/lambda_optim_%.5f.png'%tvec[ts])
            #print './tmp/lambda_optim_%.5f.png'%tvec[ts]
            
            print('chi2: %.2f'%CHI2(g_,prod[ts],resid[ts]))
            print('DoF: %.1f'%sum((1-w_i(g_))**2))
            print('gamma: %.3f'%g_)
            
            ax2.text(0.05, 0.91 , 'N-DoF: %.1f'%sum(1-(1-w_i(g_))**2),backgroundcolor='w',  transform=ax2.transAxes,color='k')


            #show()
            #draw()
            #print ts
            #fig.canvas.draw_idle()
            #pause(1)
            #ion()
            plt.show()
            #plt.savefig('/home/tom/Desktop/CLANEK/mfi_iter4/out.svg')
            
            #exit()
            #break
            
            
            
    
            print(ts)

            
            #clf()


    #exit()
    #emulate a rapid solver - use only one median value fo the whole block - it is more stable, lower noise in the timeevolution
    if rapid_solver:
        g[:] = median(g)
        
    g = minimum(maximum(g,g_lim_low), g_lim_up)
    
    
    
    
    

        #ind = Leverage > 4*std(Leverage)
        #if any(ind):
            #print '\nMost probably broken detectors: ', dets[~wrong_dets][ind]
            #title('Leverage')
            #axhline( 4*std(Leverage))
            #plot(dets[~wrong_dets], Leverage,'.');xlabel('Detector');ylabel('Leverage');ylim(0,1),show()
        
    #print '...'
    
    
    
    




        #imshow(sqrt(E[ts].reshape(150,100, order='F')))
        #show()
        
        #imshow(R.todense())
        #show()
        
        
        
        #p = array(U.T.dot(S_[ts]).T)
        #print p.shape
        #print R.shape
        #exit()
        #if LastCycle:
            ##from matplotlib.pyplot import *
            #plot(L_*E[ts,:])
            #plot(dot(asarray(U),R.T*R*p) )
            #show()
            
            #imshow(R.T*R)
                
            
        #U*R.T*R*U.T
            
        #V.T.dot(single(w*p/D)[:,None],out=ascontiguousarray(E[ts,:]))
        
        #NOTE errorbars!! err = sqrt(dot(( w_i(g_gcv)/D)**2,V**2))
    #TODO 
    #CookDistance = zeros(ndets)
    SigmaGsample = None

    if LastCycle:
        
        E = empty(( tsteps,  size(V,1) ),dtype=dtype)  #memory consuming!!
        
        
        #remove this forcycle? and the previos one too? 
        if not sparse.issparse(V):
            V = array(V,copy=False)
            
        W = np.zeros_like(prod)
        for ts,(p,r) in enumerate(zip(prod,resid)):
            W[ts] = w_i(g[ts])
            chi2[ts] = CHI2(g[ts],p,r)
            #slowest step 
            if sparse.issparse(V):
                E[ts] = V.T.dot((W[ts]*p/D).astype(dtype))
            else:
                dot((W[ts]*p/D).astype(dtype), V,out=E[ts])
            p = transpose(R)*(R*p)  #for QR decomposition
            retro[ts,~wrong_dets] = dot(asarray(U),W[ts]*p)
            

        #Leverage = einsum('ij,ij,j->i', U,U,w_i( median(g))  )        

        #p = sum(w_i( median(g)))
        #n = ndets
        #er =  mean((retro-S.T)**2,0)
        #s2 = sum(e**2)/n
        #NOTE it is just inspired by Cook's Distance!! 
        #rel_err = er/mean(er)
        #CookDistance[~wrong_dets] = rel_err[~wrong_dets]*(Leverage/(1-Leverage))
        #CookDistance[wrong_dets] = infty
        
        

        
        
        
        
        
        w = w_i(median(g))

        #import IPython
        #IPython.embed()
        if positive_constrain :
            E,chi2 = ForcePositivity(E,W,g,V,D,prod,resid,ndets, BdMat,chi2,tokamak)
            #E,chi2 = ForcePositivity_old(E,V,D,BdMat,w,chi2,nx,ny, tokamak)
            #print retro.shape, L.shape, E.shape
            retro = asarray(L*E.T).T  #retrofit was affected
            #print retro.shape, L.shape, E.shape

                
        if estimate_sigma:
            proj_err = sqrt(median(chi2))*random.randn(U.shape[0],100)
            SigmaGsample = dot(V.T,(w/D*norm.mean())[:,None]*array(dot(U.T,proj_err)))
                
            
    else:
        w = w_i(median(g))
        E = dot((w*prod.mean(0)/D), V)
        
        
        
        #f,ax = subplots(2, sharex=True)
        #ax[0].plot(CookDistance)
        #ax[0].plot(where(~wrong_dets)[0],(Leverage/(1-Leverage)**2))
        #ax[0].plot(rel_err,'--')

        ##figure()
        #ax[1].plot( median(S.T,0))
        #ax[1].plot( median(retro,0))
        ##ax[1].plot(e)

        #show()
        
    #linalg.norm(L*E.T-S,axis=0)**2/len(D)
    
    #plot(chi2)
    #plot(linalg.norm(L*E.T-S,axis=0)**2/ndets)
    #show()
    


        
    #print '...'

        
#show()

    #E2 = copy(E.mean(0))
    #imshow( E.mean(0).reshape(ny,nx,order='F'));colorbar();show()
  
    
    #figure();imshow(E[0].reshape( (ny,nx), order="F")<0,vmin=0,origin='lower')
    
    #positive_constrain = True
    #print chi2
    #solve positivchi2ity constrain !!
    #exit()
    
    #positive_constrain = True
    
    
    #TODO také by to šlo vyřešit přidáváním virtuálních detektorů do těch minim a řešit to iterativně.

        

    #print chih
    #try:
    #savez('./tmp/decompostion_%.5f_%-5f.npz'%(tvec[0],tvec[-1]), U=U,W=w_i(median(g)),
          #D=D,V=tokamak.Transform*V.T,L=L ,wrong_dets=wrong_dets,chi2=median(chi2),norm=norm )
    #except:
        #print 'decompostion_%.5f_%-5f.npz was not saved'%(tvec[0],tvec[-1])
  
    #image = E[0].reshape(ny,nx,order='F')
    #imshow( image,vmin=0,origin='lower');contour(image,(-lim,),colors='w',origin='lower');show()

  
    #figure()
    #imshow(E[0].reshape( (ny,nx), order="F")<0,vmin=0,origin='lower')
    
    
#L
    #plot(S)
    #show()
    
    
    
    ##print chi2
    
    #imshow(E[0].reshape( (ny,nx), order="F"),vmin=0,origin='lower')
    
    #subplot(121)
    #imshow( E.mean(0).reshape(ny,nx,order='F'));colorbar()#;show()
    #contour(E.mean(0).reshape(ny,nx,order='F'),(0,) )
    #subplot(122)
    #imshow( E2.reshape(ny,nx,order='F'));colorbar()#;show()
    #contour(E2.reshape(ny,nx,order='F'),(0,) )

    #show()
    
    
    #if any(isnan(E)):
    
    
    ##figure()
    #title(LastCycle)

    #imshow( image, vmin =-img_max/100, vmax = img_max/100);show()



        #wmax = amax(w)
        #ind = w/wmax>0.01
        #sort_ind = argsort(D, order=None)  # D from QR is unsorted
        #for i in sort_ind[::-1]:
            #if w[i]>wmax*0.1:
            #print sort_ind
            #print w[i]
            #print prod[i]
            #print D[i]
            ##prod = sum(U*reshape(S[:,ts], (-1,1)),0)
            #print prod
            #print reshape(S[:,ts], (-1,1))
            #print U
            
            #imshow(reshape(w[i]*prod[i]/D[i]*V[i,:], (nx,ny), order="F"))
            #colorbar()
            #show()
            #E[:,ts]  += w[i]*prod[i]/D[i]*V[i,:]

    #debug( '\nlinear solver finished, mean time: %.1fms'%((time.time()-t_total)/(ts+1)*1000))
    
    
    E = E.T
    
    
    
    
    #L*E-S
    #BUG 
    #chi2__ = (linalg.norm((L*E)[~wrong_dets]-S,axis=0)**2/ndets)
    #E[E<0] = 0
    #chi2_ = (linalg.norm((L*E)[~wrong_dets]-S,axis=0)**2/sum(~wrong_dets))
    #print E.shape
    #imshow(E[:,0].reshape( (ny,nx), order="F"),vmin=0,origin='lower')
    #show()
    
    ##plot(chi2__);plot(chi2,'--');show()
    
    
    #K_ = dot(dot(U,diag(1/D)),V)
    #K_ = dot(dot(U,diag(1/D)),V)
    #print 'Decomposition accuracy',
    #print linalg.norm((K_* L.T)-eye(L.shape[0])[~wrong_dets])
    #imshow((K_* L.T)-eye(L.shape[0])[~wrong_dets])
    #colorbar()
    #show()
    
    
    #plot(log(chi2))
    #show()
    #g0 = mquantiles(log(D**2), g_fract_0)[0]
    
    
    
    
    
    g_ = g
    #normalize gamma as quantil of the singular values
    g = interp(g, log(D**2)[::-1], linspace(0,1,len(D)))
    
    #from scipy.stats import percentileofscore
    
    #percentileofscore(log(D**2)[::-1], g_)/100
    #percentileofscore((D**2)[::-1], exp(g_))/100

    
     
      
    
    

    #if LastCycle:
        #plot_tomo_error(tvec, E, V,norm,w,D,chi2 )





    #from scikits.sparse.cholmod import cholesky
    
    #factor = cholesky( L_.T*L_ + exp(g_[0])*H )
    #E0 = factor( L_.T*S_)
            

    #print '-----------', len(E), linalg.norm( E0-E)/linalg.norm(E0)
    
    #plot
    
    
    
        
    #import matplotlib.pylab as plt
    
    
    
    
    
    
    #for i, d in enumerate(D):
        #ndof = sum((1-w_i(log(d)*2))**2)
        #print ndof/len(D),'\t\t', interp(log(d)*2, log(D**2)[::-1], linspace(0,1,len(D)))
        
        
        

    ##L, H
    


        
        
    #print linalg.norm(L_.T*(L_*E)+g_[0]*H*E-L_.T*S_)/linalg.norm(L_.T*S_)
    
    #exit()
    #plot(L.T*(L*E)+g_[0]*H*E-L.T*S)
    #plot(L.T*S)
    #show()
    #print time.time()-t_total

    #PlotBaseVectors2(V, nx,ny)
        

    return E,retro, chi2,g,SigmaGsample#,CookDistance



def signal_handler(signum, frame):
    raise Exception("SVD Timed out!")


def plot_tomo_error(tvec, E, V,norm,w,D,chi2 ):
    from fconf import usetex
    usetex()
    #from matplotlib.pyplot import *
    
    
    
    V_ = V.T*(w/D)*sqrt(median(chi2))*norm.mean()
    #print V_.dtype
    Cov = dot(V_,V_.T)
    from matplotlib.colors import LogNorm
    
    #CM = contourf((sqrt(diag(Cov))/E.mean(1)/norm.mean()).reshape(ny,nx,order='F'),levels=logspace(-3, 0, 10) ,origin='lower',vmin=1e-1,vmax=1,norm = LogNorm())
    #CM.cmap.set_over('red')
    #CM.cmap.set_under('yellow')
    #colorbar(CM);show()
    
    #imshow((sqrt(diag(Cov))/E.mean(1)/norm.mean()).reshape(ny,nx,order='F') ,origin='lower',vmin=0,vmax=0.1);show()





    global tokamak
    
    extent = tokamak.xmin,tokamak.xmax, tokamak.ymin, tokamak.ymax
    from matplotlib.ticker import MaxNLocator

    from make_graphs import my_cmap_
    f,axis = subplots(1,3, sharex=True,sharey=True)
    E[E<0] = 0
    f.set_size_inches(10,4)
    
    img1 = axis[0].imshow((E.mean(1)*norm.mean()).reshape(ny,nx,order='F')/1e3,
                          origin='lower',vmin=0,extent=extent,aspect='equal',cmap=my_cmap_)

    #levels = linspace(0,E.mean(1).max(), 10)*norm.mean()/1e3
    #img1 = axis[0].contourf(tokamak.xgrid+tokamak.dx/2,tokamak.ygrid+tokamak.dy/2,
                            #(E.mean(1)*norm.mean()).reshape(ny,nx,order='F')/1e3,levels,
                          #origin='lower',vmin=0,extent=extent,aspect='equal',cmap=my_cmap_)


    
    Sig = sqrt(diag(Cov))
    from annulus import get_bd_mat, get_rho_field_mat
    bnd = tokamak.get_boundary( 100,time=tvec.mean())
    rho = get_rho_field_mat(tokamak,  tvec.mean())
    levels = .2, .4, .6,.8, 1
    CS = axis[0].contour(tokamak.xgrid+tokamak.dx,tokamak.ygrid+tokamak.dy, rho ,levels, colors='w')
    axis[0].clabel(CS, inline=1, fontsize=10, color='w',fmt='%1.1f',inline_spacing=-2)

    BdMat = get_bd_mat(tokamak,tokamak.nx, tokamak.ny,boundary=bnd)
    
    pointx = int(nx*0.4)
    pointy = ny/2
    Sig[BdMat] = 0
    img2=axis[1].imshow((Sig/(E.mean(1)*norm.mean()+1)).reshape(ny,nx,order='F') 
                        ,origin='lower',vmin=1e-2,extent=extent,aspect='equal',
                        vmax=1,norm=LogNorm(vmin=1e-2, vmax=1),cmap=my_cmap_)
    Corr = Cov[pointy+pointx*ny].reshape(ny,nx,order='F')/ Cov[pointy+pointx*ny,pointy+pointx*ny]
    #Corr[abs(Corr)<.05] = nan
    img3=axis[2].imshow( Corr,origin='lower',vmin=-1,vmax=1,cmap='PuOr_r',extent=extent,aspect='equal')
    
    axis[0].axhline(tokamak.ygrid[pointy]+tokamak.dy/2,c='w')
    axis[0].axvline(tokamak.xgrid[pointx]+tokamak.dx/2,c='w')
    axis[1].axhline(tokamak.ygrid[pointy]+tokamak.dy/2,c='k')
    axis[1].axvline(tokamak.xgrid[pointx]+tokamak.dx/2,c='k')
    axis[2].axhline(tokamak.ygrid[pointy]+tokamak.dy/2,c='k')
    axis[2].axvline(tokamak.xgrid[pointx]+tokamak.dx/2,c='k')
    axis[0].text(0.1, 0.9 ,'a)',transform=axis[0].transAxes,color='w')
    axis[1].text(0.1, 0.9 ,'b)',transform=axis[1].transAxes,color='w')
    axis[2].text(0.1, 0.9 ,'c)',transform=axis[2].transAxes,color='k')

    
    for _,struct in list(tokamak.struct_dict.items()):
        if size(struct)>1: 
            axis[0].plot(struct[0], struct[1], 'w',lw=.5)
            axis[1].plot(struct[0], struct[1], 'k',lw=.5)
            axis[2].plot(struct[0], struct[1], 'k',lw=.5)
    #axes().set_aspect('equal')


    axis[0].axis(extent)
    axis[0].set_xlabel('R [m]')
    axis[0].set_ylabel('Z [m]', labelpad=-5)
    axis[1].set_xlabel('R [m]')
    axis[2].set_xlabel('R [m]')
    axis[0].xaxis.set_major_locator(MaxNLocator(4))
    axis[1].xaxis.set_major_locator(MaxNLocator(4))
    axis[2].xaxis.set_major_locator(MaxNLocator(4))
    axis[0].yaxis.set_major_locator(MaxNLocator(5))
    
    axis[0].set_title('Emissivity [kW/m$^3$]')
    axis[1].set_title('Relative variance')
    axis[2].set_title('Correlation')
    
    c1=colorbar(img1,ax=axis[0])
    c2=colorbar(img2,ax=axis[1])
    c3=colorbar(img3,ax=axis[2])
    tick_locator = MaxNLocator(nbins=5)
    c1.locator = tick_locator
    c1.update_ticks()
    #c2.locator = tick_locator
    #c2.update_ticks()
    c3.locator = tick_locator
    c3.update_ticks()
    show()
    
    
    
    
    
    
        
def PlotBaseVectors2(V, nx,ny):    
    #plot base vectors of the reconstruction space of the linear methods
    #from make_graphs import my_cmap2
    #V = presolved_decomposi{'U':U, 'V':V , 'D':D, 'wrong_dets':wrong_dets}tion['V']
    #ion()
    from matplotlib.ticker import NullFormatter

    
    N = size(V,0)
    #from matplotlib.pyplot import *
    

    f,ax = subplots(5,15, sharex=True, sharey=True)
    
    ax = ax.flatten()

    for i in range(len(ax)):
        solution = reshape(V[i],(ny,nx), order='F')
        vmax = amax(abs(solution))
        ax[i].imshow(solution, origin='lower', vmin = -vmax, vmax = vmax,cmap='seismic')
    
    f.subplots_adjust(hspace=0.0, wspace = 0.)
    ax[0].yaxis.set_major_formatter( NullFormatter() )
    ax[0].xaxis.set_major_formatter( NullFormatter() )
    #tight_layout()

    show()

        
    
    
def PlotBaseVectors(V, nx,ny,Nmax=infty):    
    #plot base vectors of the reconstruction space of the linear methods
    #from make_graphs import my_cmap2
    #V = presolved_decomposi{'U':U, 'V':V , 'D':D, 'wrong_dets':wrong_dets}tion['V']
    #ion()
    N = size(V,0)

    
    for i in range(min(N,Nmax)):
        solution = reshape(V[i,:],(ny,nx), order='F')
        #vmin,vmax = solution.min(),solution.max()
        vmax = amax(abs(solution))
        #vmin = min(vmin,-vmax)
        #print vmin, vmax
        imshow(solution, origin='lower', vmin = -vmax, vmax = vmax)
        colorbar()
        #im.set_clim((vmin,vmax))

        #draw()
        show()

    fig = figure()   
    fig.show()

    
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')

    im =  ax.imshow(zeros((ny,nx)), origin='lower',cmap='seismic')
    txt = ax.set_title('')
    #pause(1)
    #pause(.1)

    for i in range(min(N,Nmax)):
        print(i) 
        txt.set_text('%d/%d'%(i+1, min(N,Nmax)))

        solution = reshape(V[i,:],(ny,nx), order='F')
        vmin,vmax = solution.min(),solution.max()
        vmin = min(vmin,-vmax)
        im.set_data(solution )
        im.set_clim((vmin,vmax))

        draw()
        #show()
        pause(.1)
        #close()
        

