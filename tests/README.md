# pyTomo test suite

The tests run fully offline. `tests/fake_mds.py` provides a fake `MDSplus`
module that serves a synthetic DIII-D shot:

- an analytic equilibrium (shifted elliptic flux surfaces) as EFIT01 gEQDSK/aEQDSK nodes,
- raw slow poloidal SXR signals (`PTDATA2("SX90R{P,M}1S{ch}")`) computed as line
  integrals of a known emissivity along the real DIII-D chords and divided by the
  calibration, so the full loader chain (filtering, offsets, calibration) is exercised.

Shot 163303 has a peaked SXR profile, shot 163304 a pedestal-like profile with a
steep edge gradient. The known phantom gives an absolute accuracy reference.

## Running

```bash
python setup.py build_ext --inplace          # Cython geometry matrix generator
export PYSPECVIEW_PATH=/path/to/pyspecview   # optional, enables the roto_tomo tests
python -m pytest                              # ~40 s
```

GUI tests use `QT_QPA_PLATFORM=offscreen` (set automatically).
Environments: `requirements-min.txt` (Python 3.12, oldest supported libraries,
strict baseline tolerances) and `requirements-current.txt`.

## Files

| file | content |
|---|---|
| `test_d3d_sxr.py` | DIII-D SXR loading, calibration, equilibrium, geometry matrix, reconstruction vs. baseline and phantom, saved profile format, known limitations (xfail) |
| `test_downstream_api.py` | functions used by `pyspecview/roto_tomo`, including the real `roto_tomo` module |
| `test_gui.py` | main window, resizing, reconstruction in the worker thread, GUI stall time |
| `test_cli.py` | `pytomo SHOT -t 21 -R 1 ...` writes the `Emissivity_*.npz` file |
| `baseline/` | reference reconstructions (`make_baseline.py`) and timings (`tools/benchmark.py`) |
| `record_mds.py` | records the MDSplus data of a real run for offline replay |

Baselines were recorded with Python 3.12 / NumPy 1.26 / SciPy 1.13 / Matplotlib 3.8.
Regenerate them only for an intended change of results:
`python -m pytest tests/make_baseline.py -s`.

## Real data

At DIII-D, record a shot once and replay it in tests:

```bash
python tests/record_mds.py rec_180000 -- 180000 -t 21 --tmin 2 --tmax 2.1 -R 1
```

`record_mds.load_recording('rec_180000')` returns a server for `fake_mds.install()`.
