"""
pyTomo - Advanced tomography reconstruction code for fusion plasma diagnostics

pyTomo is an advanced code for tomographical inversion from various line
integrated diagnostics like SXR or bolometers on fusion tokamaks (DIII-D,
ASDEX Upgrade, JET, TCV, etc.). The main advantages are a simple GUI,
very fast inversion, and high accuracy (depends on the quality of input data).

Key features:
- Multiple inversion methods (MFI, SVD, GEV, GSVD)
- Support for multiple tokamaks and diagnostics
- Fast Cython-accelerated geometry matrix generation
- Multiprocessing support for parallel reconstruction
"""

__version__ = "1.0.2"
__author__ = "Tomas Odstrcil"

# Import main classes and functions for easy access


def pytomo_class(*args, **kwargs):
    from .pytomo import pytomo_class as _cls
    return _cls(*args, **kwargs)
    

__all__ = [ '__version__']
