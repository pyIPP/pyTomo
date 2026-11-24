"""Record reference reconstructions of the synthetic DIII-D SXR shot.

Run once with the unmodified code in the oldest supported environment:
    python -m pytest tests/make_baseline.py -s
"""
import json
import os
import time

import numpy as np

from cases import CASES
from conftest import BASELINE, reconstruct, phantom_on_grid


def test_make_baseline(d3d_server, d3d_sxr):
    _, eq = d3d_server
    _, tok = d3d_sxr
    os.makedirs(BASELINE, exist_ok=True)
    meta = {}
    for name, kw in CASES.items():
        t = time.time()
        _, out = reconstruct(tok, **kw)
        dt = time.time() - t
        P = phantom_on_grid(tok, eq, out['tvec'])
        err = np.linalg.norm(out['g'].T - P) / np.linalg.norm(P)
        np.savez_compressed(os.path.join(BASELINE, 'd3d_sxr_%s.npz' % name),
                            g=out['g'].astype(np.float32), retro=out['retro'].astype(np.float32),
                            chi2=out['chi2'], lam0=out['lam0'], tvec=out['tvec'],
                            data=out['data'], error=out['error'], dets=out['dets'])
        meta[name] = dict(phantom_err=float(err), runtime_s=round(dt, 2))
        print(name, meta[name])
    with open(os.path.join(BASELINE, 'meta.json'), 'w') as f:
        json.dump(meta, f, indent=1)
