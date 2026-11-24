"""LDL' adapter must give the same results with scikit-sparse 0.4 and >= 0.5."""
import numpy as np
import pytest
from scipy import sparse

from pytomo.sparse_cholesky import ldl, to_dense, CholmodError


def spd(n=60, seed=0):
    A = sparse.random(n, n, 0.1, random_state=seed)
    return (A @ A.T + sparse.eye(n)).tocsc()


def test_factorization_and_solves():
    A = spd()
    f = ldl(A)
    L_inv_PA = to_dense(f.solve_L(f.apply_P(A.toarray())))  # L^-1 P A
    D = f.D()
    # P A P' = L D L'  =>  L^-1 P A P' L^-T = D
    M = to_dense(f.solve_L(f.apply_P(L_inv_PA.T)))
    assert np.allclose(M, np.diag(D), atol=1e-10)
    b = np.random.default_rng(0).random((60, 3))
    x = f.apply_Pt(f.solve_Lt(to_dense(f.solve_L(f.apply_P(b))) / D[:, None]))
    assert np.allclose(A @ x, b)


def test_sparse_rhs_and_beta():
    A = spd()
    f = ldl(A, 0.5)
    b = sparse.random(60, 5, 0.2, format='csr', random_state=1)
    y = to_dense(f.solve_L(f.apply_P(b)))
    assert y.shape == (60, 5)
    x = f.apply_Pt(f.solve_Lt(y / f.D()[:, None]))
    assert np.allclose((A + 0.5 * sparse.eye(60)) @ x, b.toarray())


def test_float32_input_is_factorized_in_float64():
    f = ldl(spd().astype(np.float32))
    assert f.D().dtype == np.float64


def test_singular_matrix_raises():
    A = sparse.diags(np.r_[np.ones(9), 0.0]).tocsc()
    with pytest.raises(CholmodError):
        ldl(A)
