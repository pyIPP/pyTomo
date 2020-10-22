// C wrapper to SparseSuiteQR library et al. for Python

// We pass in the sparse matrix data in a COO sparse matrix format. Cholmod
// refers to this as a "cholmod_triplet" format. This is then converted to its
// "cholmod_sparse" format, which is a CSC matrix.

#include <stdio.h>
#include <stdlib.h>
#include "SuiteSparseQR_C.h"
// #include "SuiteSparseQR.hpp"

#define Long SuiteSparse_long
#define Int SuiteSparse_long
 
#include "cholmod_core.h"

// using namespace std;


void qr_solve(double const *A_data, long const *A_row, long const *A_col,
              size_t A_nnz, size_t A_m, size_t A_n, double const *b_data, double *x_data) {
    // Solves the matrix equation Ax=b where A is a sparse matrix and x and b
    // are dense column vectors. A and b are inputs, x is solved for in the
    // least squares sense using a rank-revealing QR factorization.
    //
    // Inputs
    //
    // A_data, A_row, A_col: the COO data
    // A_nnz: number of non-zero entries, ie the length of the arrays A_data, etc
    // A_m: number of rows in A
    // A_n: number of cols in A
    // b_data: the data in b. It is A_m entries long.
    //
    // Outputs
    //
    // x_data: the data in x. It is A_n entries long
    //
    // MAKE SURE x_data is allocated to the right size before calling this function
    //
    cholmod_common Common, *cc;
    cholmod_sparse *A_csc;
    cholmod_triplet *A_coo;
    cholmod_dense *b, *x;
    size_t k;
    // Helper pointers
    long *Ai, *Aj;
    double *Ax, *bx, *xx;

    /* start CHOLMOD */
    cc = &Common ;
    cholmod_l_start (cc) ;

    // Create A, first as a COO matrix, then convert to CSC
    A_coo = cholmod_l_allocate_triplet(A_m, A_n, A_nnz, 0, CHOLMOD_REAL, cc);
    
    
    if (A_coo == NULL) {
        fprintf(stderr, "ERROR: cannot allocate triplet");
        return;
    }
    // Copy in data
    Ai = A_coo->i;
    Aj = A_coo->j;
    Ax = A_coo->x;
    for (k=0; k<A_nnz; k++) {
        Ai[k] = A_row[k];
        Aj[k] = A_col[k];
        Ax[k] = A_data[k];
    }
    
    
    
    A_coo->nnz = A_nnz;
    // Make sure the matrix is valid
    if (cholmod_l_check_triplet(A_coo, cc) != 1) {
        fprintf(stderr, "ERROR: triplet matrix is not valid");
        return;
    }
    // Convert to CSC
    A_csc = cholmod_l_triplet_to_sparse(A_coo, A_nnz, cc);

//     for (k=0; k<A_nnz; k++) {
// //         cout<<"  %d   \n",(A_coo->i)[k];
// //         printf("  %L   \n",(Aj)[k]);      
//     }
//      for (k=0; k<A_nnz; k++) {
// //         std::cout <<A_csc->p[k]<<std::endl;
//         printf("  %3.2f   \n",(A_coo->x)[k]);      
// 
//     }
    

    
    
    
    // Create b as a dense matrix
    b = cholmod_l_allocate_dense(A_m, 1, A_m, CHOLMOD_REAL, cc);
    bx = b->x;
    for (k=0; k<A_m; k++) {
        bx[k] = b_data[k];
    }
    // Make sure the matrix is valid
    if (cholmod_l_check_dense(b, cc) != 1) {
        fprintf(stderr, "ERROR: b vector is not valid");
        return;
    }

    // Solve for x
    x = SuiteSparseQR_C_backslash_default(A_csc, b, cc);

    // Return values of x
    xx = x->x;
    for (k=0; k<A_n; k++) {
        x_data[k] = xx[k];
        
        
    }
    

    /* free everything and finish CHOLMOD */
    cholmod_l_free_triplet(&A_coo, cc);
    cholmod_l_free_sparse(&A_csc, cc);
    cholmod_l_free_dense(&x, cc);
    cholmod_l_free_dense(&b, cc);
    cholmod_l_finish(cc);
    return;
}

// void qr_solve_csr(double const *A_data, long const *A_p, long const *A_i, size_t A_nnz, size_t A_m, size_t A_n, double const *b_data, double *x_data) {
//TODO dodělat impementaci té fukce solve ale pro csr

void qr_solve_csr(double const *A_data, long const *A_i_data, long const *A_p_data, 
                  size_t A_nnz, size_t A_m, size_t A_n, double const *b_data, double *x_data) {
    // Solves the matrix equation Ax=b where A is a sparse matrix and x and b
    // are dense column vectors. A and b are inputs, x is solved for in the
    // least squares sense using a rank-revealing QR factorization.
    //
    // Inputs
    //
    // A_data, A_row, A_col: the COO data
    // A_nnz: number of non-zero entries, ie the length of the arrays A_data, etc
    // A_m: number of rows in A
    // A_n: number of cols in A
    // b_data: the data in b. It is A_m entries long.
    //
    // Outputs
    //
    // x_data: the data in x. It is A_n entries long
    //
    // MAKE SURE x_data is allocated to the right size before calling this function
    //
    cholmod_common Common, *cc;
    cholmod_sparse *A_csc;
    cholmod_dense *b, *x;
    size_t k;
    
    // Helper pointers
//     Long *A_p,*A_i;
    double *bx, *xx;
    int stype ;         
    double *Ax ;


    /* start CHOLMOD */
    cc = &Common ;
    cholmod_l_start (cc) ;
    int sorted,packed;
    sorted=1;
    packed=1;
    stype = 0;
    A_csc = cholmod_l_allocate_sparse(A_m, A_n, A_nnz, sorted,packed,stype, CHOLMOD_REAL, cc);
        
    if (A_csc == NULL) {
        fprintf(stderr, "ERROR: cannot allocate sparse");
        return;
    }
//     cholmod_dense *B ;
//     cholmod_sparse *Bsparse ;
    
//     B = cholmod_l_ones (A_m, 1, A_csc->xtype, cc) ;
//     Bsparse = cholmod_l_dense_to_sparse (B, 1, cc) ;
    
    
    
    
    
    
//     double *Ax, */*Az*/ ;
//     cholmod_sparse *A ;
    Int *Ap, *Ai ;
//     Int  n ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

//     RETURN_IF_NULL_COMMON (NULL) ;
//     Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate the matrix */
    /* ---------------------------------------------------------------------- */

//     n =  A_m<A_n?A_m:A_n ;
//     A = CHOLMOD(allocate_sparse) (nrow, ncol, n, TRUE, TRUE, 0, xtype,
//             Common) ;
// 
//     if (Common->status < CHOLMOD_OK)
//     {
//         return (NULL) ;     /* out of memory or inputs invalid */
//     }
//             


    
    /* ---------------------------------------------------------------------- */
    /* create the identity matrix */
    /* ---------------------------------------------------------------------- */

    Ap = A_csc->p ;
    Ai = A_csc->i ;
    Ax = A_csc->x ;


    for (k=0; k<A_nnz; k++) {
        Ai[k] = A_i_data[k];
        Ax[k] = A_data[k];
    }
    for (k=0; k<=A_n; k++) {
        Ap[k] = A_p_data[k];
    }
     
    

//     TODO Make sure the matrix is valid
    if (cholmod_l_check_sparse(A_csc, cc) != 1) {
        fprintf(stderr, "ERROR: sparse matrix is not valid\n");
        return;
    }
//     return;


    // Create b as a dense matrix
    b = cholmod_l_allocate_dense(A_m, 1, A_m, CHOLMOD_REAL, cc);
    bx = b->x;
    for (k=0; k<A_m; k++) {
        bx[k] = b_data[k];
    }
    // Make sure the matrix is valid
    if (cholmod_l_check_dense(b, cc) != 1) {
        fprintf(stderr, "ERROR: b vector is not valid");
        return;
    }

    // Solve for x
    x = SuiteSparseQR_C_backslash_default(A_csc, b, cc);
    // Return values of x
    xx = x->x;
    for (k=0; k<A_n; k++) {
        x_data[k] = xx[k];
    }
//     return;

    /* free everything and finish CHOLMOD */
//     cholmod_l_free_triplet(&A_coo, cc);
    cholmod_l_free_sparse(&A_csc, cc);
    cholmod_l_free_dense(&x, cc);
    cholmod_l_free_dense(&b, cc);
    cholmod_l_finish(cc);
    return;
}


int qr_sparse_csc(double const *A_data, long const *A_i_data, long const *A_p_data, size_t A_nnz, size_t A_m, size_t A_n,
    double **Qx, long **Qi, long **Qp, size_t *Q_nnz, size_t *Q_m, size_t *Q_n,
    double **Rx, long **Ri, long **Rp, size_t *R_nnz, size_t *R_m, size_t *R_n, long** E) {

    cholmod_common Common, *cc;
    cholmod_sparse *A_csc,*Q, *R, *I;
    double *Qx_, *Rx_;
    long *Ri_, *Rp_,*Qi_, *Qp_;
    size_t k;
    
    // Helper pointers
    int stype ; 

    /* start CHOLMOD */
    cc = &Common ;
    cholmod_l_start (cc) ;
    int sorted,packed;
    sorted=1;
    packed=1;
    stype = 0;
    A_csc = cholmod_l_allocate_sparse(A_m, A_n, A_nnz, sorted,packed,stype, CHOLMOD_REAL, cc);
    
    
    
    I = cholmod_l_speye (A_m, A_m, CHOLMOD_REAL, cc) ;
    
    

        
    if (A_csc == NULL) {
        fprintf(stderr, "ERROR: cannot allocate sparse");
        return 1;
    }

    Int *Ap, *Ai ;
    double *Ax ;


    Ap = A_csc->p ;
    Ai = A_csc->i ;
    Ax = A_csc->x ;


    for (k=0; k<A_nnz; k++) {
        Ai[k] = A_i_data[k];
        Ax[k] = A_data[k];
    }
    for (k=0; k<=A_n; k++) {
        Ap[k] = A_p_data[k];
    }
     
    

    if (cholmod_l_check_sparse(A_csc, cc) != 1) {
        fprintf(stderr, "ERROR: sparse matrix is not valid\n");
        return 1;
    }
    Long* Etmp = malloc(A_n*sizeof(Long));
//     Long* Etmp = NULL;

    
    Long econ = A_n;
//     Long econ = 0;

    /* number of rows of R and columns of Q 
    to return. The default is m. Usingn gives the standard 
    economy form (as in the MATLAB qr(A,0)). A value less 
    thanthe estimated rank r is set to r, so opts.econ=0 
    gives the “rank-sized” factorization,where size(R,1)==nnz(diag(R))==r.*/



    int ordering = SPQR_ORDERING_DEFAULT;
    double tol = SPQR_DEFAULT_TOL;
    int getCTX = 1; //??

//     SuiteSparseQR_C_QR(ordering,tol,econ,A_csc,&Q,&R,&Etmp,cc);
    
    SuiteSparseQR_C             // returns rank(A) estimate, (-1) if failure
(
    // inputs:
    ordering,           // all, except 3:given treated as 0:fixed
    tol,             // columns with 2-norm <= tol are treated as 0
    econ,              // e = max(min(m,econ),rank(A))
    getCTX,             // 0: Z=C (e-by-k), 1: Z=C', 2: Z=X (e-by-k)
    A_csc,      // m-by-n sparse matrix to factorize
    I,// sparse m-by-k B
    NULL, // dense  m-by-k B
    // outputs:
    &Q,   // sparse Z
    NULL,    // dense Z
    &R,     // R factor, e-by-n1
    &Etmp,               // size n column permutation, NULL if identity
    NULL,     // m-by-nh Householder vectors
    NULL,           // size m row permutation
    NULL,   // 1-by-nh Householder coefficients
    cc      // workspace and parameters
);


    double  *tmpR_x;
    
    Int *tmpR_i,*tmpR_p,*tmpQ_p,*tmpQ_i;
    
    R_nnz[0] = R->nzmax;
    R_m[0] = R->nrow;
    R_n[0] = R->ncol;
    
    
    tmpR_p = R->p ;
    tmpR_i = R->i ;
    tmpR_x = R->x ;
    
    Q_nnz[0] = Q->nzmax;
    Q_m[0] = Q->nrow;
    Q_n[0] = Q->ncol;

    tmpQ_p = Q->p ;
    tmpQ_i = Q->i ;

    double *tmpQ_x = (double *)(Q->x);

    Rp_ = malloc((*R_n+1)*sizeof(long));
    Ri_ = malloc(*R_nnz*sizeof(long));
    Rx_ = malloc(*R_nnz*sizeof(double));
    
    Qp_ = malloc((*Q_n+1)*sizeof(long));
    Qi_ = malloc(*Q_nnz*sizeof(long));
    Qx_ = malloc(*Q_nnz*sizeof(double));
    
    
    //memory inicialization has failured
    if((Rp_==NULL)|(Ri_==NULL)|(Rx_==NULL)|(Qp_==NULL)|(Qi_==NULL)|(Qx_==NULL))
    {
        fprintf(stderr, "ERROR: cannot allocate memory for Q and R");
        return 1;
    }

//     for(k = 0; k <= A_n;k++)
//     {
//         E_[k] = Etmp[k];
//     }
    
    for(k = 0; k <= *R_n;k++)
    {
        Rp_[k] = tmpR_p[k];
    }
    for(k = 0; k < *R_nnz;k++)
    {
        Rx_[k] = tmpR_x[k];
        Ri_[k] = tmpR_i[k];

    }
    for(k = 0; k <= *Q_n;k++)
    {
        Qp_[k] = tmpQ_p[k];
    }
    for(k = 0; k < *Q_nnz;k++)
    {
        Qx_[k] = tmpQ_x[k];
        Qi_[k] = tmpQ_i[k];

    }
    
    
    *Qx = Qx_;    *Qi = Qi_;    *Qp = Qp_;
    *Rx = Rx_;    *Rp = Rp_;    *Ri = Ri_;
    *E = Etmp;
    

    /* free everything and finish CHOLMOD */
    cholmod_l_free_sparse(&A_csc, cc);
    cholmod_l_free_sparse(&Q, cc);
    cholmod_l_free_sparse(&R, cc);
    cholmod_l_finish(cc);
    
    
    

    return 0;
};




// =============================================================================
// === SuiteSparseQR_C_QR ======================================================
// =============================================================================

// [Q,R,E] = qr(A), returning Q as a sparse matrix




// extern "C" {
    
    
// #include "SuiteSparseQR_C.h"




   
int qr_sparse_csc_r(double const *A_data, long const *A_i_data, long const *A_p_data, size_t A_nnz, size_t A_m, size_t A_n,
    double **Rx, long **Ri, long **Rp, size_t *R_nnz, size_t *R_m, size_t *R_n, long ** E) {

    cholmod_common Common, *cc;
    cholmod_sparse *A_csc, *R;
//     cholmod_dense * /*Bdense*/;
    double  *Rx_;
    long *Ri_, *Rp_;
    size_t k;
    Long *E_;
    
    // Helper pointers
    int stype,getCTX ; 
    
    getCTX = 0; //??

    // start CHOLMOD 
    cc = &Common ;
    cholmod_l_start (cc) ;
    int sorted,packed;
    sorted=1;
    packed=1;
    stype = 0;


    
    A_csc = cholmod_l_allocate_sparse(A_m, A_n, A_nnz, sorted,packed,stype, CHOLMOD_REAL, cc);
        
    if (A_csc == NULL) {
        fprintf(stderr, "ERROR: cannot allocate sparse");
        return 1;
    }

    Int *Ap, *Ai ;
    double *Ax ;

    Long* Etmp = malloc(A_n*sizeof(Long));

    Ap = A_csc->p ;
    Ai = A_csc->i ;
    Ax = A_csc->x ;


    for (k=0; k<A_nnz; k++) {
        Ai[k] = A_i_data[k];
        Ax[k] = A_data[k];
    }
    for (k=0; k<=A_n; k++) {
        Ap[k] = A_p_data[k];
    }
     
    

    if (cholmod_l_check_sparse(A_csc, cc) != 1) {
        fprintf(stderr, "ERROR: sparse matrix is not valid\n");
        return 1;
    }
    
//  To return full-sized results, set econ = m.  Then C and R will have m rows,
//  and C' will have m columns.
//
//  To return economy-sized results, set econ = n.  Then C and R will have k
//  rows and C' will have k columns, where k = min(m,n).
//
//  To return rank-sized results, set econ = 0.  Then C and R will have k rows
//  and C' will have k columns, where k = r = the estimated rank of A.
    
    
    Long econ = A_n;
//     econ = 0;  

    int ordering = SPQR_ORDERING_DEFAULT;  
    double tol = SPQR_DEFAULT_TOL;

//     SuiteSparseQR_C_R(ordering,tol,econ,A_csc,&R,&Etmp,cc);
 
  SuiteSparseQR_C /* returns rank(A) estimate, (-1) if failure */
(
    /* inputs: */
     ordering,               /* all, except 3:given treated as 0:fixed */
     tol,                 /* columns with 2-norm <= tol treated as 0 */
     econ,      /* e = max(min(m,econ),rank(A)) */
    getCTX,                 /* 0: Z=C (e-by-k), 1: Z=C', 2: Z=X (e-by-k) */
    A_csc,          /* m-by-n sparse matrix to factorize */
    NULL,    /* sparse m-by-k B */
    NULL,     /* dense  m-by-k B */
    /* outputs: */
    NULL,   /* sparse Z */
    NULL,    /* dense Z */
    &R,         /* e-by-n sparse matrix */
    &Etmp,       /* size n column perm, NULL if identity */
    NULL,         /* m-by-nh Householder vectors */
    NULL,   /* size m row permutation */
    NULL,       /* 1-by-nh Householder coefficients */
    cc          /* workspace and parameters */
) ;
    

    double  *tmpR_x;
    
    Int *tmpR_i,*tmpR_p;
    
    R_nnz[0] = R->nzmax;
    R_m[0] = R->nrow;
    R_n[0] = R->ncol;
    
    
    tmpR_p = R->p ;
    tmpR_i = R->i ;
    tmpR_x = R->x ;
    
 

    Rp_ = malloc((*R_n+1)*sizeof(long));
    Ri_ = malloc(*R_nnz*sizeof(long));
    Rx_ = malloc(*R_nnz*sizeof(double));
    E_ = malloc(A_n*sizeof(Long));


    //memory inicialization has failured
    if((Rp_==NULL)|(Ri_==NULL)|(Rx_==NULL))
    {
        fprintf(stderr, "sparse qr: ERROR: cannot allocate memory for  R");
        return 1;
    }


    for(k = 0; k < A_n;k++)
    {
        E_[k] = Etmp[k];
    }
    
    for(k = 0; k <= *R_n;k++)
    {
        Rp_[k] = tmpR_p[k];
    }
    for(k = 0; k < *R_nnz;k++)
    {
        Rx_[k] = tmpR_x[k];
        Ri_[k] = tmpR_i[k];

    }

    
    
    *Rx = Rx_;    *Rp = Rp_;    *Ri = Ri_;
    *E = E_;


    // free everything and finish CHOLMOD 
    cholmod_l_free_sparse(&A_csc, cc);
    cholmod_l_free_sparse(&R, cc);
    cholmod_l_finish(cc);
    
    
    

    return 0;
}
