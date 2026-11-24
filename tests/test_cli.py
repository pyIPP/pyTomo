"""Command line workflow: `python -m pytomo.pytomo SHOT [options]` (run in-process)."""
import glob
import os
import sys

import numpy as np

from conftest import D3D_SHOT, default_inputs


def test_cli_reconstruction_writes_profiles(d3d_server, monkeypatch):
    from pytomo import pytomo
    import conftest
    conftest._loaded.clear()  # main() creates its own DIIID object (resets the global loader)
    inputs = default_inputs()
    for f in glob.glob(os.path.join(inputs['tmp_folder'], 'Emissivity_*.npz')):
        os.remove(f)
    argv = ['pytomo', str(D3D_SHOT), '-t', '21', '--tmin', '2.0', '--tmax', '2.02',
            '-R', '1', '-S', '1', '--solver', '3', '-x', '80', '-y', '120']
    monkeypatch.setattr(sys, 'argv', argv)
    pytomo.main()
    files = glob.glob(os.path.join(inputs['tmp_folder'], 'Emissivity_*%d.npz' % D3D_SHOT))
    assert len(files) == 1
    d = np.load(files[0], allow_pickle=True)
    assert d['gres'].dtype == np.float16
    assert d['gres'].shape[:2] == (120, 80)
    assert np.isclose(d['tvec'][0], 2.0, atol=1e-3)
