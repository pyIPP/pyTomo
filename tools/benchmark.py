"""Time the main stages of the DIII-D SXR workflow on the synthetic shot.

    python tools/benchmark.py [out.json]

Reports wall time per stage and the peak resident memory of the process.
"""
import json
import os
import resource
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, 'tests'), ROOT]
import conftest  # noqa: E402  (isolated HOME, fake MDSplus)
import fake_mds  # noqa: E402
import numpy as np  # noqa: E402


def main(out=None):
    res = {}
    fake_mds.synthetic_d3d({conftest.D3D_SHOT: 'peaked'}, conftest.SXR_GEOM + os.sep)
    conftest._prepare_local_files()
    from pytomo.prepare_data import loaddata
    from pytomo.mat_deriv_B import mat_deriv_B
    from pytomo.make_graphs import make_graphs

    def timed(name, f, *a, **k):
        t = time.perf_counter()
        r = f(*a, **k)
        res[name] = round(time.perf_counter() - t, 3)
        return r

    inputs = conftest.default_inputs()
    tok = timed('load data + equilibrium', loaddata, inputs)
    timed('geometry matrix', tok.prepare_tokamak)
    timed('prepare data (51 slices)', tok.prepare_data, 2.0, 2.05, 1, 1, 1.0, None, False)
    timed('regularization matrix (aniso MFI)', mat_deriv_B, tok, np.array([2.0, 2.05]), 1, {'danis': 4})
    for name, kw in [('SVD2', dict(solver=3)), ('SVD2 + positivity', dict(solver=3, positive_constrain=True)),
                     ('QR', dict(solver=4))]:
        inp, o = timed('reconstruction %s (51 slices)' % name, conftest.reconstruct, tok, regularization=1, **kw)
    timed('plots + profile export', make_graphs, (inp, tok, None, o))
    res['peak RSS MB'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    res['versions'] = dict(python=sys.version.split()[0], numpy=np.__version__)
    print(json.dumps(res, indent=1))
    if out:
        json.dump(res, open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:])
