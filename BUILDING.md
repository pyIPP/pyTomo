# Building pyTomo

## Compiling Cython Extensions

pyTomo includes a Cython module (`geom_mat_gen_cython.pyx`) that must be compiled before use. This is **required for multiprocessing to work correctly**.

### Quick Start

```bash
cd /path/to/pyTomo
./build_cython.sh
```

Or manually:

```bash
cd /path/to/pyTomo
python setup.py build_ext --inplace
```

### What This Does

- Compiles `geom_mat_gen_cython.pyx` → `geom_mat_gen_cython.c` (C code)
- Compiles C code → `geom_mat_gen_cython.*.so` (shared library)
- Places the `.so` file in the same directory for easy import

### Why Pre-compilation is Required

The previous approach used `pyximport` for just-in-time compilation, but this causes issues with `multiprocessing`:

1. **Problem**: Worker processes spawned by `multiprocessing` don't inherit the dynamically compiled module
2. **Result**: `ModuleNotFoundError: No module named 'geom_mat_gen_cython'` in worker processes
3. **Solution**: Pre-compile the module with `setup.py` so it's available to all processes

### Requirements

- Python >= 3.8
- Cython >= 0.29
- NumPy
- C compiler (gcc, clang, or MSVC)

### Installation for Development

For development, install pyTomo in editable mode after building:

```bash
cd /path/to/pyTomo
python setup.py build_ext --inplace
pip install -e .
```

### Cleaning Build Artifacts

```bash
rm -f geom_mat_gen_cython.c
rm -f geom_mat_gen_cython.*.so
rm -rf build/
```

### Troubleshooting

**Error: "No module named 'geom_mat_gen_cython'"**
- Solution: Run `./build_cython.sh` or `python setup.py build_ext --inplace`

**Error: "gcc: command not found" or similar**
- Solution: Install a C compiler
  - Linux: `sudo apt-get install build-essential` (Debian/Ubuntu) or `sudo yum install gcc` (RHEL/CentOS)
  - macOS: `xcode-select --install`
  - Windows: Install Microsoft Visual C++ Build Tools

**Error: "numpy/arrayobject.h: No such file or directory"**
- Solution: Ensure NumPy is installed: `pip install numpy`

### For Conda-forge Feedstock

When creating a conda-forge feedstock, the build script should include:

```yaml
build:
  script: {{ PYTHON }} setup.py build_ext --inplace && {{ PYTHON }} -m pip install . -vv
```

This ensures the Cython module is compiled during the conda build process.
