#!/bin/bash
# Clean all Cython build artifacts

echo "Cleaning Cython build artifacts..."

cd "$(dirname "$0")"

# Remove compiled Cython files
rm -f geom_mat_gen_cython.c
rm -f geom_mat_gen_cython*.so
rm -f geom_mat_gen_cython*.pyd

# Remove build directories
rm -rf build/
rm -rf dist/
rm -rf *.egg-info/

# Remove Python bytecode
find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
find . -name "*.pyc" -delete 2>/dev/null || true

echo "Clean complete!"
echo ""
echo "Now run: python setup.py build_ext --inplace"
