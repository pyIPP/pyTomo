"""Shared fixtures.

All tests run against an isolated HOME (pyTomo writes ~/tomography) and the
offline MDSplus replacement from fake_mds.py, so no database access is needed.
"""
import os
import shutil
import sys
import tempfile

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS = os.path.join(ROOT, 'tests')
for p in (ROOT, TESTS):
    if p not in sys.path:
        sys.path.insert(0, p)

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
os.environ.setdefault('MPLBACKEND', 'Agg')

# isolated home, must be set before pytomo resolves '~'
_HOME = tempfile.mkdtemp(prefix='pytomo_home_')
os.environ['HOME'] = _HOME

import fake_mds  # noqa: E402

D3D_SHOT = 163303           # peaked SXR profile
D3D_SHOT_PEDESTAL = 163304  # pedestal-like SXR profile
SXR_GEOM = os.path.join(ROOT, 'pytomo', 'geometry', 'DIIID', 'SXR')
BASELINE = os.path.join(TESTS, 'baseline')


def _prepare_local_files():
    # SXR calibration settings are read from the *local* geometry path (see report)
    dst = os.path.join(_HOME, 'tomography', 'geometry', 'DIIID', 'SXR')
    os.makedirs(dst, exist_ok=True)
    for f in os.listdir(SXR_GEOM):
        if f.endswith('.dat'):
            shutil.copy(os.path.join(SXR_GEOM, f), dst)


@pytest.fixture(scope='session')
def d3d_server():
    _prepare_local_files()
    server, eq = fake_mds.synthetic_d3d({D3D_SHOT: 'peaked', D3D_SHOT_PEDESTAL: 'pedestal'},
                                        SXR_GEOM + os.sep)
    return server, eq


def default_inputs(**kw):
    from pytomo import pytomo, config
    inputs = pytomo.loadSetting()
    inputs.update(tok_index=21, shot=D3D_SHOT, tmin=2.0, tmax=2.05, n_cpu=2,
                  enable_output=False, save_profiles=True)
    inputs.update(kw)
    for d in ('tmp_folder', 'output_path'):
        os.makedirs(inputs[d], exist_ok=True)
    config.useCache = False
    config.wrong_dets_pref = np.array([], dtype=int)
    config.wrong_dets_defined = []
    return inputs


def load_d3d_sxr(shot):
    from pytomo.prepare_data import loaddata
    np.random.seed(0)
    inputs = default_inputs(shot=shot)
    tok = loaddata(inputs)
    tok.prepare_tokamak()
    return inputs, tok


_loaded = {}


def _get_d3d_sxr(shot):
    # DIIID keeps its data loader in a module global, so only the most recently
    # created DIIID object returns correct data; reload when switching shots.
    if _loaded.get('shot') != shot:
        _loaded.update(shot=shot, value=load_d3d_sxr(shot))
    return _loaded['value']


@pytest.fixture
def d3d_sxr(d3d_server):
    """DIII-D slow SXR tokamak object with geometry matrix (shared, do not modify)."""
    return _get_d3d_sxr(D3D_SHOT)


@pytest.fixture
def d3d_sxr_pedestal(d3d_server):
    return _get_d3d_sxr(D3D_SHOT_PEDESTAL)


def reconstruct(tok, **kw):
    """Run the main reconstruction workflow (as in the CLI) and return its output dict."""
    from pytomo.main import tomography
    from pytomo import config
    inputs = default_inputs(shot=tok.shot, **kw)
    inputs['tokamak_tmp'] = tok
    config.wrong_dets_pref = np.array([], dtype=int)
    np.random.seed(0)
    _, _, _, output = tomography(inputs, tok)
    return inputs, output


def phantom_on_grid(tok, eq, tvec, profile='peaked'):
    """True emissivity (npix, nt) on the pixel centres, Fortran pixel order as in pyTomo."""
    xc = tok.xgrid + tok.dx / 2
    yc = tok.ygrid + tok.dy / 2
    X, Y = np.meshgrid(xc, yc)
    E = fake_mds.PROFILES[profile](eq.psin(X, Y)).ravel(order='F')
    return E[:, None] * fake_mds.time_trace(np.asarray(tvec))[None, :]
