"""Sparse LDL' factorization (CHOLMOD through scikit-sparse), P A P' = L D L'.

Works with scikit-sparse 0.4 and >= 0.5, which changed the API. The returned
factor has the 0.4 interface used by the solvers: D(), solve_L, solve_Lt,
apply_P, apply_Pt. The factorization is always computed in float64.
"""
import numpy as np
from scipy import sparse
from sksparse import cholmod

CholmodError = cholmod.CholmodError
CholmodWarning = cholmod.CholmodWarning
NEW_API = hasattr(cholmod, 'ldl_factor')  # scikit-sparse >= 0.5


class _Factor:
    """Adapter of the scikit-sparse >= 0.5 factor to the 0.4 method names."""

    def __init__(self, f):
        self.f = f
        self.perm = f.perm
        self.iperm = np.argsort(f.perm)

    def D(self):
        return self.f.D.diagonal()

    def solve_L(self, b):
        return self.f.solve(b, system='L')

    def solve_Lt(self, b):
        return self.f.solve(b, system='Lt')

    def apply_P(self, b):
        return b[self.perm]

    def apply_Pt(self, b):
        return b[self.iperm]


def ldl(A, beta=0.0, order='default'):
    """Simplicial LDL' factorization of the symmetric matrix A + beta*I.

    Raises CholmodError if A is not positive definite."""
    A = sparse.csc_matrix(A, dtype=np.float64)
    if NEW_API:
        return _Factor(cholmod.ldl_factor(A, beta, order=order))
    return cholmod.cholesky(A, beta, mode='simplicial', ordering_method=order)


def to_dense(x):
    """Solves return sparse or dense arrays depending on the version and input."""
    return x.toarray() if sparse.issparse(x) else np.asarray(x)
