"""Record the MDSplus data used by a pyTomo run, to replay it offline in the tests.

At DIII-D (with database access):

    python tests/record_mds.py OUTDIR -- 180000 -t 21 --tmin 2 --tmax 2.1 -R 1

runs the normal pyTomo CLI with every `Connection.get(expr)` result also written to
OUTDIR (one .npy file per expression, safe with the multiprocessing data loaders).
`load_recording(OUTDIR)` returns a fake_mds.FakeServer serving these data.
"""
import hashlib
import json
import os
import sys

import numpy as np


def _key(expr):
    return hashlib.sha1(expr.encode()).hexdigest()[:16]


def install_recorder(outdir):
    import MDSplus
    os.makedirs(outdir, exist_ok=True)
    Base = MDSplus.Connection

    class RecordingConnection(Base):
        def get(self, expr, *args):
            res = Base.get(self, expr, *args)
            try:
                value = res.data()
                k = _key(expr)
                np.save(os.path.join(outdir, k + '.npy'), np.asarray(value), allow_pickle=True)
                with open(os.path.join(outdir, k + '.json'), 'w') as f:
                    json.dump({'expr': expr}, f)
            except Exception as e:  # never break the real run
                print('record_mds: could not store %s: %s' % (expr, e))
            return res

    MDSplus.Connection = RecordingConnection


def load_recording(outdir):
    import fake_mds
    signals = {}
    for f in os.listdir(outdir):
        if f.endswith('.json'):
            expr = json.load(open(os.path.join(outdir, f)))['expr']
            v = np.load(os.path.join(outdir, f[:-5] + '.npy'), allow_pickle=True)
            signals[expr] = v.item() if v.dtype == object and v.ndim == 0 else v
    return fake_mds.FakeServer(signals)


if __name__ == '__main__':
    outdir, sep, args = sys.argv[1], sys.argv[2], sys.argv[3:]
    assert sep == '--', __doc__
    install_recorder(outdir)
    from pytomo import pytomo
    sys.argv = ['pytomo'] + args
    pytomo.main()
