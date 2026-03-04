#!/bin/bash
# Build script for pyTomo Cython extensions
# This script compiles the Cython module in-place so it can be imported

set -e  # Exit on error

echo "Building pyTomo Cython extensions..."

# Navigate to the pyTomo directory
cd "$(dirname "$0")"

# Clean old build artifacts first
echo "Cleaning old build artifacts..."
rm -f geom_mat_gen_cython.c
rm -f geom_mat_gen_cython*.so
rm -f geom_mat_gen_cython*.pyd
rm -rf build/

# Build the Cython extension in-place
echo "Compiling Cython extension..."
python setup.py build_ext --inplace

echo ""
echo "Build complete!"
echo "The compiled module 'geom_mat_gen_cython' is now available for import."
echo ""
echo "Compiled files created:"
ls -lh geom_mat_gen_cython*.so geom_mat_gen_cython.c 2>/dev/null || echo "  (check for .so or .pyd files)"
