"""Pytomo interfaces used by pyspecview/roto_tomo.

The roto_tomo tests import the real downstream module; point PYSPECVIEW_PATH to a
pyspecview checkout (the directory containing roto_tomo/), otherwise they are skipped.
"""
import importlib
import os
import sys

import numpy as np
import pytest
from scipy import sparse

from conftest import D3D_SHOT
import fake_mds

T0 = 2.0
PYSPECVIEW = os.environ.get('PYSPECVIEW_PATH', '')


def flux_surfaces(eq, n_rho=60, n_theta=120):
    """(n_theta, n_rho) R, Z of the analytic flux surfaces, as eqm.rhoTheta2rz returns them."""
    rho = np.linspace(0, 1, n_rho)
    th = np.linspace(0, 2 * np.pi, n_theta)[:, None]
    R = eq.R0 + eq.shift * (1 - rho**2) + eq.a * rho * np.cos(th)
    Z = eq.Z0 + eq.kappa * eq.a * rho * np.sin(th)
    return rho, R, Z


class AuxTokamak:
    """Copy of roto_tomo.roto_tomo.tokamak (the minimal tokamak interface pytomo needs)."""
    norm = 1
    allow_negative = False
    allow_divertor = False
    transform_index = 0
    camera = False
    beta = 0
    boundary = 0
    weight_power = 1
    use_pfm = False
    smooth_boundary = True

    def __init__(self, rhop, magr, magz, xgrid, ygrid):
        self.dx = np.mean(np.diff(xgrid))
        self.dy = np.mean(np.diff(ygrid))
        self.xgrid, self.ygrid = xgrid, ygrid
        self.nx, self.ny = len(xgrid), len(ygrid)
        self.npix = self.nx * self.ny
        self.xmin, self.ymin = xgrid[0], ygrid[0]
        self.xmax, self.ymax = self.dx + xgrid[-1], self.dy + ygrid[-1]
        self.rhop, self.magr, self.magz = rhop, magr, magz

    def get_boundary(self, *var, **kwarg):
        R0, Z0 = self.magr[:, 0].mean(), self.magz[:, 0].mean()
        ext = 1. if self.boundary else 1.2
        return np.array(((self.magr[::-1, -1] - R0) * ext + R0, (self.magz[::-1, -1] - Z0) * ext + Z0)).T

    def mag_equilibrium(self, tvec, return_mean):
        if return_mean:
            return self.rhop, self.magr, self.magz
        return self.rhop, self.magr[:, :, None], self.magz[:, :, None]


@pytest.fixture(scope='module')
def aux():
    eq = fake_mds.Equilibrium()
    rho, R, Z = flux_surfaces(eq)
    xgrid = np.linspace(0.95, 2.45, 80, endpoint=False)
    ygrid = np.linspace(-1.38, 1.38, 120, endpoint=False)
    tok = AuxTokamak(rho, R, Z, xgrid, ygrid)
    X, Y = np.meshgrid(xgrid + tok.dx / 2, ygrid + tok.dy / 2)
    return tok, eq.psin(X, Y)


def test_imported_names():
    from pytomo.mat_deriv_B import mat_deriv_B  # noqa: F401
    from pytomo.shared_modules import read_config  # noqa: F401
    from pytomo.graph_derivation import graph_derivation  # noqa: F401
    from pytomo.geom_mat_setting import geom_mat_setting, prune_geo_matrix  # noqa: F401
    from pytomo.annulus import get_bd_mat, get_rho_field_mat  # noqa: F401
    from pytomo.make_graphs import LogFormatterTeXExponent  # noqa: F401
    import pytomo.config  # noqa: F401
    # roto_tomo does `import pytomo.geometry.DIIID as Tok; Tok(diag, inputs, ...)`
    import pytomo.geometry.DIIID as D3D
    import pytomo.geometry.ASDEX as AUG
    assert callable(D3D) and callable(AUG)


def test_read_config_types():
    from pytomo.shared_modules import read_config
    cfg = read_config(os.path.join(os.path.dirname(fake_mds.__file__), '..', 'pytomo', 'tomography_D3D.cfg'))
    assert cfg['shot'] == 163303 and cfg['nx'] == 80
    assert cfg['danis'] == 4 and isinstance(cfg['rgmin'], float)
    assert cfg['lambda_solver'] == 'gcv'
    assert cfg['wrong_dets'] == [] and cfg['use_phantom_data'] is None
    assert cfg['positive_constrain'] is False and cfg['reconstruct'] is True
    assert cfg['magfield_shift_lcfs'] == [0.0, 0.0]
    assert cfg['max_cook_dist'] == np.inf


def test_rho_field_and_boundary(aux):
    from pytomo.annulus import get_bd_mat, get_rho_field_mat
    tok, psin = aux
    rho = get_rho_field_mat(tok, T0).reshape(psin.shape)
    inside = psin < 0.95
    assert np.allclose(rho[inside], np.sqrt(psin[inside]), atol=0.01)
    tok.boundary = 1  # boundary exactly at the separatrix
    outside = get_bd_mat(tok, time=T0).reshape(tok.ny, tok.nx, order='F')
    tok.boundary = 0
    clear = np.abs(psin - 1) > 0.1
    assert np.mean(outside[clear] == (psin[clear] > 1)) > 0.99


def test_anisotropic_derivative_matrices(aux):
    from pytomo.mat_deriv_B import mat_deriv_B
    tok, psin = aux
    Bmat, diag_mat = mat_deriv_B(tok, 0, 1, {})
    assert len(Bmat) == 4 and len(diag_mat) == 5
    assert all(sparse.issparse(B) and B.shape == (tok.npix, tok.npix) for B in Bmat)
    f = np.exp(-3 * psin).ravel(order='F')  # flux function
    Bper, Bpar = Bmat[0], Bmat[1]
    core = (psin < 0.8).ravel(order='F')
    assert np.linalg.norm((Bpar @ f)[core]) < 0.1 * np.linalg.norm((Bper @ f)[core])


@pytest.mark.parametrize('reg,n', [(0, 2), (2, 2), (4, 1)])
def test_isotropic_derivative_matrices(aux, reg, n):
    from pytomo.mat_deriv_B import mat_deriv_B
    tok, _ = aux
    Bmat, _ = mat_deriv_B(tok, 0, reg, {})
    assert len(Bmat) == n
    const = np.ones(tok.npix)
    inner = np.zeros((tok.ny, tok.nx), dtype=bool)
    inner[1:-1, 1:-1] = True  # vertical derivative wraps around between columns at the grid edge
    inner = inner.ravel(order='F')
    if reg != 4:  # derivative of a constant vanishes
        for B in Bmat:
            assert np.abs(B @ const)[inner].max() < 1e-10


def test_prune_geo_matrix():
    from pytomo.geom_mat_setting import prune_geo_matrix
    T = sparse.random(10, 12, density=0.5, format='csr', random_state=0)
    BdMat = np.zeros((3, 4), dtype=bool)
    BdMat[0] = True
    Tp = prune_geo_matrix(T, BdMat)
    assert Tp.shape == T.shape
    dense = Tp.toarray()
    assert np.all(dense[:, BdMat.ravel(order='F')] == 0)


def test_theta_star(d3d_server, d3d_sxr, aux):
    _, eq = d3d_server
    _, tok = d3d_sxr
    rho, R, Z = flux_surfaces(eq, n_theta=100)
    th_rz, th = tok.mag_theta_star(T0, rho, R, Z, rz_grid=True)
    assert th_rz.shape == (tok.nx, tok.ny)
    assert th.shape == (len(rho), 100)
    assert np.allclose(th[1:, -1], 2 * np.pi * 99 / 100 * 100 / 99, rtol=1e-6) or np.allclose(th[1:, 0], 0)
    assert np.all(np.abs(th_rz) <= np.pi + 1e-6)


# --------------------------------------------------------- real roto_tomo

@pytest.fixture(scope='module')
def roto_tomo():
    if not os.path.isdir(os.path.join(PYSPECVIEW, 'roto_tomo')):
        pytest.skip('set PYSPECVIEW_PATH to a pyspecview checkout')
    pytest.importorskip('PyQt5')
    sys.path.insert(0, PYSPECVIEW)
    return importlib.import_module('roto_tomo.roto_tomo')


def test_roto_tomo_regularization_pipeline(roto_tomo, d3d_server, d3d_sxr):
    """Reproduces the pytomo calls of Roto_tomo.prepare_tomo with an analytic equilibrium."""
    from pytomo.annulus import get_bd_mat, get_rho_field_mat
    from pytomo.geom_mat_setting import geom_mat_setting, prune_geo_matrix
    from pytomo.mat_deriv_B import mat_deriv_B
    _, eq = d3d_server
    _, tok = d3d_sxr
    rhop, magr, magz = flux_surfaces(eq, n_rho=50, n_theta=100)
    T_full = geom_mat_setting(tok, tok.nx, tok.ny, tok.virt_chord, path=None)[0]
    xgridc = (tok.xgrid + tok.dx / 2) / tok.norm
    ygridc = (tok.ygrid + tok.dy / 2) / tok.norm
    aux = roto_tomo.tokamak(rhop, magr, magz, tok.xgrid, tok.ygrid)
    aux.vessel_boundary = tok.vessel_boundary

    theta_star_rz, theta_star = tok.mag_theta_star(T0, rhop, magr, magz, rz_grid=True)
    theta_star_rz, theta_star = theta_star_rz.T, theta_star.T
    rho_mat = get_rho_field_mat(aux, T0)
    BdMat = get_bd_mat(aux, time=T0).reshape(tok.ny, tok.nx, order='F')
    T = prune_geo_matrix(T_full, BdMat)
    Bmat, _ = mat_deriv_B(aux, 0, 1, {})
    Ht = roto_tomo.build_reg_mat_time(rho_mat, theta_star_rz, theta_star, rhop, magr, magz,
                                      xgridc, ygridc, BdMat, 2 * np.pi / 100)
    assert T.shape == T_full.shape == (tok.nl, tok.npix)
    assert len(Bmat) == 4
    assert sparse.issparse(Ht) or isinstance(Ht, np.ndarray)

    # HarmSolver (sparse LDL' + SVD) on noise-free synthetic data, run in-process
    import queue
    from conftest import phantom_on_grid
    G0 = phantom_on_grid(tok, eq, [T0])[:, 0] + 1e-3
    Hper, Hpar = roto_tomo.create_derivation_matrix(G0.reshape(-1, 1), Bmat, 4)
    b = T @ G0
    good = b > 0.02 * b.max()
    berr = np.where(good, 0.01 * b.max(), np.inf)
    npix = tok.npix
    E = sparse.spdiags(np.exp(1j * theta_star_rz.flatten('F')), 0, npix, npix)
    W = sparse.spdiags(1 / G0, 0, npix, npix, format='csc')
    for n in (0, 1):  # n = 1 exercises the complex part of the solve
        qin, qout = queue.Queue(), queue.Queue()
        h = roto_tomo.HarmSolver(qin, qout, n, E, W, Ht, Hper, Hpar, T, BdMat, b, berr, tok.nl)
        qin.put(0.5)
        qin.put(None)
        h.run()
        G, retro, lam, chi2 = qout.get()
        assert G.shape == (npix,) and np.all(np.isfinite(G))
        if n == 0:
            inside = ~BdMat.flatten('F')
            assert np.corrcoef(G.real[inside], G0[inside])[0, 1] > 0.9
