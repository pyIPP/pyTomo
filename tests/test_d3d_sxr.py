"""End-to-end characterization of the DIII-D slow SXR workflow on a synthetic shot.

Data and equilibrium are served by the fake MDSplus server (tests/fake_mds.py).
Reference results in tests/baseline were produced by the original code
(tests/make_baseline.py); refactoring must keep them within the tolerances below.
"""
import glob
import json
import os

import numpy as np
import pytest

import fake_mds
from cases import CASES
from conftest import BASELINE, reconstruct, phantom_on_grid

T0 = 2.0  # inside the flat top of fake_mds.time_trace


def rel(a, b):
    return np.linalg.norm(np.asarray(a) - b) / np.linalg.norm(b)


# ------------------------------------------------------------------ loading

def test_tokamak_setup(d3d_sxr):
    _, tok = d3d_sxr
    assert tok.name == 'DIIID' and tok.index == 21
    assert list(tok.detectors_dict) == ['90RP1a', '90RP1b', '90RM1a', '90RM1b']
    assert tok.nl == 64 and len(tok.dets) == 64
    assert tok.Tmat.shape == (64, tok.nx * tok.ny) == (64, 80 * 120)
    assert np.allclose(tok.vessel_boundary[0], tok.vessel_boundary[-1])


def test_prepare_data_float32_and_calibrated(d3d_server, d3d_sxr):
    _, eq = d3d_server
    inputs, tok = d3d_sxr
    data, error, tvec, dets, Tmat, normData = tok.prepare_data(
        T0, T0 + 0.05, 1, 1, 1.0, None, False)
    assert data.dtype == np.float32 and error.dtype == np.float32
    assert data.shape == error.shape == (64, len(tvec))
    assert np.all(error[np.isfinite(error)] > 0)
    # slow DAS digitizes only channels 1-28 of each array for new shots
    good = np.r_[0:28, 32:60]
    bright = fake_mds.sxr_brightness(eq)
    ratio = np.median(data[good], 1) / bright[good]
    # a few channels are rescaled by load_SXR.hardcoded_corrections (real detector issues)
    assert abs(np.median(ratio) - 1) < 0.015
    assert np.mean(abs(ratio - 1) < 0.02) > 0.9, ratio


def test_unused_channels_are_flagged(d3d_sxr):
    inputs, tok = d3d_sxr
    data = tok.prepare_data(T0, T0 + 0.05, 1, 1, 1.0, None, False)[0]
    ok = tok.get_correct_dets(data)
    assert not ok[np.r_[28:32, 60:64]].any()
    assert ok[np.r_[0:28, 32:60]].sum() >= 54


def test_get_data_is_repeatable(d3d_sxr):
    # regression: DIIID.get_data used to inflate the cached errors in place
    _, tok = d3d_sxr
    e1 = tok.get_data(tmin=T0, tmax=T0 + 0.05)[2].copy()
    e2 = tok.get_data(tmin=T0, tmax=T0 + 0.05)[2]
    assert np.array_equal(e1, e2)


# ------------------------------------------------------------- equilibrium

def test_flux_surfaces_match_equilibrium(d3d_server, d3d_sxr):
    _, eq = d3d_server
    _, tok = d3d_sxr
    rho, magx, magy = tok.mag_equilibrium(np.array([T0]), return_mean=True)
    assert tok.radial_coordinate == 'r_V'
    psin = eq.psin(magx, magy)  # (n_theta, n_rho)
    expected = np.interp(rho, fake_mds.r_V(eq, np.linspace(0, 1, 200)), np.linspace(0, 1, 200))
    inner = rho > 0.1
    assert np.allclose(psin.mean(0)[inner], expected[inner], atol=0.02)
    # the Fourier-moment parameterization of the surfaces deviates near the edge
    spread = psin.std(0)
    assert np.all(spread[rho < 0.7] < 0.01)
    assert np.all(spread[inner] < 0.05)


def _rho_field_error(eq, tok):
    from pytomo.annulus import get_rho_field_mat
    rho = get_rho_field_mat(tok, T0)
    X, Y = np.meshgrid(tok.xgrid + tok.dx / 2, tok.ygrid + tok.dy / 2)
    psin = eq.psin(X, Y)
    return psin, np.abs(rho.reshape(psin.shape) - fake_mds.r_V(eq, psin))


def test_rho_field_matrix(d3d_server, d3d_sxr):
    psin, err = _rho_field_error(d3d_server[1], d3d_sxr[1])
    assert err[psin < 0.5].max() < 0.015
    assert err[psin < 1.2].max() < 0.08  # current edge accuracy (parameterized surfaces)


@pytest.mark.xfail(strict=True, reason='target of direct psin use instead of surface parameterization (spec 11.2/14)')
def test_rho_field_matrix_target(d3d_server, d3d_sxr):
    psin, err = _rho_field_error(d3d_server[1], d3d_sxr[1])
    assert err[psin < 1.0].max() < 0.01


def test_boundary_matrix(d3d_server, d3d_sxr):
    from pytomo.annulus import get_bd_mat
    _, eq = d3d_server
    _, tok = d3d_sxr
    outside = get_bd_mat(tok, time=T0).reshape(tok.ny, tok.nx, order='F')
    X, Y = np.meshgrid(tok.xgrid + tok.dx / 2, tok.ygrid + tok.dy / 2)
    psin = eq.psin(X, Y)
    clear = np.abs(psin - 1) > 0.1  # ignore pixels cut by the separatrix
    assert np.mean(outside[clear] == (psin[clear] > 1)) > 0.99


# ---------------------------------------------------------------- geometry

def test_geometry_matrix_line_integrals(d3d_server, d3d_sxr):
    _, eq = d3d_server
    _, tok = d3d_sxr
    G = phantom_on_grid(tok, eq, [T0])[:, 0]
    bright = fake_mds.sxr_brightness(eq)
    retro = tok.Tmat @ G
    good = bright > 0.05 * bright.max()
    assert rel(retro[good], bright[good]) < 0.05


# ---------------------------------------------------------- reconstruction

def _baseline(name):
    return np.load(os.path.join(BASELINE, 'd3d_sxr_%s.npz' % name))


@pytest.mark.parametrize('name', list(CASES))
def test_reconstruction_matches_baseline(d3d_server, d3d_sxr, name):
    _, eq = d3d_server
    _, tok = d3d_sxr
    ref = _baseline(name)
    meta = json.load(open(os.path.join(BASELINE, 'meta.json')))[name]
    _, out = reconstruct(tok, **CASES[name])

    np.testing.assert_array_equal(out['tvec'], ref['tvec'])
    np.testing.assert_array_equal(out['dets'], ref['dets'])
    assert out['g'].shape == ref['g'].shape
    assert rel(out['data'], ref['data']) < 1e-5
    assert rel(out['g'], ref['g']) < 2e-3
    assert rel(out['retro'], ref['retro']) < 2e-3
    assert np.allclose(out['chi2'], ref['chi2'], rtol=0.02, atol=1e-3)

    err = rel(out['g'].T, phantom_on_grid(tok, eq, out['tvec']))
    assert err < meta['phantom_err'] + 0.005


def test_svd2_reconstruction_accuracy(d3d_server, d3d_sxr):
    _, eq = d3d_server
    _, tok = d3d_sxr
    _, out = reconstruct(tok, solver=3, regularization=1)
    P = phantom_on_grid(tok, eq, out['tvec'])
    assert rel(out['g'].T, P) < 0.10
    # total emitted power (sum over pixels) within 3%
    assert abs(out['g'].sum() / P.sum() - 1) < 0.03
    assert np.all(np.isfinite(out['chi2'])) and np.median(out['chi2']) < 2


def test_output_dtypes(d3d_sxr):
    """Documents the current dtypes; the spec asks for float32 emissivity (see README/report)."""
    _, tok = d3d_sxr
    _, out = reconstruct(tok, solver=3, regularization=1)
    assert out['data'].dtype == np.float32
    assert out['g'].dtype in (np.float32, np.float64)


def test_saved_profiles_file(d3d_sxr):
    from pytomo.make_graphs import make_graphs
    inputs, tok = d3d_sxr
    inp, out = reconstruct(tok, solver=3, regularization=1, plot_all=False)
    tmp = inp['tmp_folder']
    for f in glob.glob(os.path.join(tmp, 'Emissivity_*.npz')):
        os.remove(f)
    make_graphs((inp, tok, None, out))
    files = glob.glob(os.path.join(tmp, 'Emissivity_*%d.npz' % tok.shot))
    assert len(files) == 1
    d = np.load(files[0], allow_pickle=True)
    for key in ('gres', 'tvec', 'rvec', 'zvec', 'inputs', 'gres_norm', 'BdMat'):
        assert key in d
    assert d['gres'].dtype == np.float16
    assert d['gres'].shape == (tok.ny, tok.nx, len(out['tvec']))
    G = d['gres'].astype(np.float32) * d['gres_norm'][None, None]
    G = G.reshape(-1, G.shape[-1], order='F').T
    assert rel(G, out['g']) < 2e-3


# ---------------------------------------------- known limitations (targets)

def test_pedestal_profile_current_accuracy(d3d_server, d3d_sxr_pedestal):
    _, eq = d3d_server
    _, tok = d3d_sxr_pedestal
    _, out = reconstruct(tok, solver=3, regularization=1)
    err = rel(out['g'].T, phantom_on_grid(tok, eq, out['tvec'], 'pedestal'))
    assert err < 0.30  # baseline 0.28: the edge gradient is smeared by MFI smoothing


@pytest.mark.xfail(strict=True, reason='target of the pedestal-aware regularization (spec section 10)')
def test_pedestal_profile_target_accuracy(d3d_server, d3d_sxr_pedestal):
    _, eq = d3d_server
    _, tok = d3d_sxr_pedestal
    _, out = reconstruct(tok, solver=3, regularization=1)
    assert rel(out['g'].T, phantom_on_grid(tok, eq, out['tvec'], 'pedestal')) < 0.15


@pytest.mark.xfail(strict=True, reason='positivity is a soft penalty, not a hard constraint (spec section 9.2)')
def test_positivity_is_enforced(d3d_sxr_pedestal):
    _, tok = d3d_sxr_pedestal
    _, out = reconstruct(tok, solver=3, regularization=7, positive_constrain=True)
    g = out['g']
    assert g.min() >= -1e-3 * g.max()
