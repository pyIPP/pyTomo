# Modernization log

Concise record of the refactoring (spec: "Pytomo Modernization and Refactoring").
Status legend: **done** = implemented and verified by tests, **open** = not started.

## Environments tested

| | Python | NumPy | SciPy | Matplotlib | scikit-sparse | PyQt5 |
|---|---|---|---|---|---|---|
| oldest | 3.12.3 | 1.26.4 | 1.13.1 | 3.8.4 | 0.4.16 | 5.15.11 |
| current | 3.13.18 | 2.5.4 | 1.18.1 | 3.11.2 | 0.4.16 | 5.15.11 |

Not tested: Python 3.14, PyQt6/PySide6, real DIII-D/AUG/TCV database access, DIII-D bolometers.

## Stage 0: inspection and test harness (done)

The repository had no tests. Added an offline harness (`tests/`, see `tests/README.md`):
fake MDSplus server with a synthetic DIII-D shot, 35 tests covering DIII-D SXR loading,
equilibrium, geometry, reconstruction (SVD2, QR, positivity, MFI/imgesson),
saved profile format, CLI, GUI (offscreen) and the roto_tomo interface.
Reference reconstructions in `tests/baseline/`.

Baseline (oldest environment, 80x120 pixels, 51 time slices, 2 CPUs):

| stage | time [s] |
|---|---|
| load data + equilibrium | 2.25 |
| geometry matrix | 0.24 |
| prepare data | 0.97 |
| anisotropic regularization matrix | 0.11 |
| reconstruction SVD2 / SVD2+positivity / QR | 1.80 / 1.85 / 2.12 |
| plots + profile export | 0.42 |
| peak RSS | 432 MB |

Phantom error (relative L2) on the peaked profile: SVD2 aniso MFI 8.1%, iso MFI 8.9%,
QR 8.2%, imgesson 24.5%. Pedestal profile: 28% (aniso MFI).

## Stage 1: bug fix found by the tests (done)

- `DIIID.get_data` / `ASDEX.get_data` modified the loader cache in place
  (`np.hypot(..., out=error)`), so the errors grew by ~8% on every call; repeated
  reconstructions (e.g. pressing START twice in the GUI) gave different results.
  First-call results are unchanged (verified bit-identical).

## Stage 2: library compatibility, explicit imports (done, partially)

- `matplotlib.cm.get_cmap` (removed in 3.9) -> `matplotlib.colormaps[...]`.
- `scipy.ndimage.interpolation/morphology` -> `scipy.ndimage`, `integrate.simps` -> `simpson`.
- NumPy 2 exports `max, min, abs, round, pow, bool` in `from numpy import *`, which
  silently changes the meaning of builtins in all star-importing modules
  (`annulus.get_rho_field_mat` failed with TypeError). `tools/dewildcard.py` rewrites
  star imports of numpy/pyplot/pylab to `np.`/`plt.` names using the origin table of the
  reference environment, so behaviour is preserved. Applied to 11 core/DIII-D modules;
  **46 files still use star imports** (next: GUI, make_graphs, main, tokamak, machines).

Results in the oldest environment are within 2e-3 of the baseline. In the current
environment they differ by 1-8%, traced to the DIII-D flux-surface parameterization
(`mag_equ.Descur_fit_core`, order-20 polynomial fit in float32): the fitted surfaces
move by up to 2 cm between NumPy versions. This is resolved by stage "equilibrium" below.

## Findings for the next stages

Compatibility
- scikit-sparse 0.5 changed the cholmod API (`solve_L`, `apply_P`, `mode=` removed);
  the code requires `scikit-sparse<0.5`. Either pin it or isolate the few cholmod calls
  of the kept SVD2/QR path behind one small function.
- 128 uses of `interp1d`, 5 of `scipy.fftpack` (legacy APIs, not removed yet).
- `exec` in `tokamak.py`, `eval` of user input for `--wrong_dets` in `pytomo.py`.

Correctness
- MFI solvers 0/1 (`minfisher.py`) are broken (`prepare_inputs()` signature mismatch):
  confirms they are obsolete. `PresolveSVD2` exists but is unused; solver 3 ("SVD2" in
  the docs) calls `PresolveSVD3`. GEV_singular (6) takes ~120 s.
- The positivity option is a quadratic penalty (L-BFGS-B), not a constraint: with
  imgesson, negative values down to -19% of the peak remain (`test_positivity_is_enforced`, xfail).
- Emissivity output `g` is float64 (spec: float32).
- `DIIID` keeps its data loader in a module global: only the most recently created
  DIIID object returns correct data (affects roto_tomo and any second tokamak object).
- `annulus.get_bd_mat` re-seeds the global NumPy RNG (`np.random.seed(0)`).
- `DIIID.SXR`: `loader_SXR(..., ~toroidal)` uses bitwise not on a bool (works by accident).
- DIII-D SXR calibration settings (`SXRsettings*.dat`) are read from the local
  `~/tomography/geometry` copy, which a fresh installation does not have.
- `pytomo.geometry.DIIID` is the class, not the subpackage (shadowed in
  `geometry/__init__.py`); roto_tomo relies on this.
- Vertical derivative matrices wrap around between pixel columns at the grid edge.
- Edge accuracy of flux surfaces / rho matrix: up to 0.07 in rho near the separatrix
  (parameterization error, `test_rho_field_matrix_target`, xfail).

GUI
- Creating a second `MainWindow` in one QApplication segfaults; the window must be
  deleted through Qt (`deleteLater`) before Python's GC, otherwise exit segfaults.
- START deletes the whole tmp folder (`shutil.rmtree(tmp_folder)`).
- During a reconstruction the event loop stalls for at most ~0.3 s (worker QThread).
- Plots are written to PNG and loaded into fixed-size `QLabel`s.

## Proposed order of the remaining stages

1. Explicit imports in the remaining modules (tool + pyflakes + tests), remove `exec`/`eval`.
2. Solver reduction to SVD2 + QR; remove minfisher, GEV, GSVD and dead code; float32 policy
   with float64 solver option; string identifiers for solvers/regularizations.
3. Equilibrium: use psi(R, Z) directly (removes descur parameterization and caches);
   rewrite `annulus` on psin; update `mat_deriv_B`, plotting, roto_tomo.
4. Machine/diagnostic structure (`machines/`, string identifiers like `D3D_SXR_fast`,
   machine defaults in config), remove the global loader.
5. Pedestal-aware regularization and a true positivity constraint (benchmarked with the
   pedestal phantom; xfail tests become the acceptance tests).
6. GUI: embedded Matplotlib canvases, layouts, config window, Outputs tab; then CLI/README.

Blocked without user input: real DIII-D data (record with `tests/record_mds.py`), AUG/TCV
data access, DIII-D bolometer test data, decision on PyQt5 vs Qt6 for the deployment.
