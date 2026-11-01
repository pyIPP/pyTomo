"""Reconstruction cases shared by make_baseline.py and the regression tests."""
CASES = {
    # name: overrides of the default DIII-D SXR inputs (solver 3 = SVD2, 4 = QR)
    'svd2_mfi_aniso':  dict(solver=3, regularization=1),
    'svd2_mfi_iso':    dict(solver=3, regularization=0),
    'svd2_imgesson':   dict(solver=3, regularization=7),
    'svd2_positive':   dict(solver=3, regularization=1, positive_constrain=True),
    'qr_mfi_aniso':    dict(solver=4, regularization=1),
}
