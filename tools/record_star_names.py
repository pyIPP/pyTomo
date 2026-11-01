"""Record where each name of `from X import *` comes from (run in the reference env:
numpy 1.26, matplotlib 3.8). Output: tools/star_names.json  {module: {name: qualified}}."""
import builtins
import json
import os

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab

ORIGINS = [('np', np), ('np.linalg', np.linalg), ('np.random', np.random), ('np.fft', np.fft),
           ('plt', plt), ('mlab', mlab)]
out = {}
for mod in ['numpy', 'matplotlib.pyplot', 'matplotlib.pylab', 'pylab', 'scipy']:
    ns = {}
    exec('from %s import *' % mod, ns)
    table = {}
    for name, obj in ns.items():
        if name.startswith('__'):
            continue
        if getattr(builtins, name, None) is obj:
            continue  # builtin (e.g. pylab restores max/min/abs/round)
        for prefix, m in ORIGINS:
            if getattr(m, name, object()) is obj:
                table[name] = prefix + '.' + name
                break
        else:
            table[name] = None  # unknown origin, must be reviewed
    out[mod] = table
json.dump(out, open(os.path.join(os.path.dirname(__file__), 'star_names.json'), 'w'), indent=0, sort_keys=True)
print({m: (len(t), sorted(k for k, v in t.items() if v is None)[:20]) for m, t in out.items()})
