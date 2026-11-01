"""Offline stand-in for an MDSplus server, used by the test suite.

`FakeServer` answers `Connection.get(expr)` from a dict of TDI expression ->
array. The dict is either generated synthetically (`synthetic_d3d_shot`) or
recorded from a real server with `tests/record_mds.py`, so the same tests can
run against replayed real DIII-D data.

`install(server)` puts a fake `MDSplus` module into `sys.modules`, so code that
does `import MDSplus; MDSplus.Connection(host)` (including forked worker
processes) talks to the fake server.
"""
import re
import sys
import types

import numpy as np


class _Result:
    def __init__(self, value):
        self._value = value

    def data(self):
        return self._value


class FakeServer:
    def __init__(self, signals, rules=()):
        self.signals = dict(signals)
        self.rules = list(rules)  # (regex, callable(match) -> array)
        self.missing = []

    def lookup(self, expr):
        key = expr.strip()
        if key in self.signals:
            return self.signals[key]
        for pattern, func in self.rules:
            m = re.fullmatch(pattern, key)
            if m:
                return func(m)
        self.missing.append(key)
        raise KeyError('fake MDSplus: no data for %s' % key)


class Connection:
    server = None  # FakeServer shared by all connections

    def __init__(self, hostspec='fake', *args, **kwargs):
        self.hostspec = hostspec

    def openTree(self, tree, shot):
        pass

    def closeTree(self, tree, shot):
        pass

    def closeAllTrees(self):
        pass

    def get(self, expr, *args):
        return _Result(Connection.server.lookup(expr))


def install(server):
    mod = types.ModuleType('MDSplus')
    mod.Connection = Connection
    Connection.server = server
    sys.modules['MDSplus'] = mod
    return mod


# ---------------------------------------------------------------- synthetic DIII-D

class Equilibrium:
    """Shifted elliptic flux surfaces: R = R0 + d*(1-rho^2) + a*rho*cos(t), Z = Z0 + k*a*rho*sin(t)."""
    R0, Z0, a, kappa, shift = 1.70, 0.02, 0.60, 1.6, 0.05
    psi_axis, psi_sep = -0.30, 0.10  # Wb/rad

    def psin(self, R, Z):
        R, Z = np.broadcast_arrays(np.asarray(R, float), np.asarray(Z, float))
        rho2 = np.zeros_like(R)
        for _ in range(60):  # fixed point iteration, converges for small shift
            x = (R - self.R0 - self.shift * (1 - np.minimum(rho2, 1))) / self.a
            y = (Z - self.Z0) / (self.kappa * self.a)
            rho2 = x**2 + y**2
        return rho2

    def boundary(self, n=101):
        t = np.linspace(0, 2 * np.pi, n)
        return np.c_[self.R0 + self.a * np.cos(t), self.Z0 + self.kappa * self.a * np.sin(t)]


def emissivity(psin, amp=1e5):
    """Peaked emissivity profile in W/m^3, zero outside separatrix."""
    return amp * np.maximum(0, 1 - psin) ** 1.5


def emissivity_pedestal(psin, amp=1e5):
    """Broad H-mode like profile with a steep edge gradient at psin ~ 0.95."""
    return amp * (1 - 0.5 * np.minimum(psin, 1)) * 0.5 * (1 - np.tanh((psin - 0.95) / 0.015))


PROFILES = {'peaked': emissivity, 'pedestal': emissivity_pedestal}


def line_integrals(eq, R1, Z1, R2, Z2, func, n=2000):
    """Brightness [W/m^2] = integral of func(psin) along straight chords."""
    s = np.linspace(0, 1, n)[:, None]
    R = R1 + s * (R2 - R1)
    Z = Z1 + s * (Z2 - Z1)
    L = np.hypot(R2 - R1, Z2 - Z1)
    return func(eq.psin(R, Z)).mean(0) * L


def r_V(eq, psin, n=600):
    """Volume radial coordinate sqrt(V(psin)/V(1)) computed numerically on a fine grid."""
    R, Z = np.meshgrid(np.linspace(0.9, 2.5, n), np.linspace(-1.2, 1.2, n))
    p = eq.psin(R, Z).ravel()
    w = R.ravel()[np.argsort(p)]
    V = np.cumsum(w)
    ps = np.sort(p)
    return np.sqrt(np.interp(psin, ps, V) / np.interp(1, ps, V))


def efit_signals(eq, shot, tree='EFIT01', t_eq=None, nR=65, nZ=129):
    """gEQDSK/aEQDSK nodes used by geometry/DIIID/map_equ.py."""
    t_eq = np.arange(0.1, 6.0, 0.05) if t_eq is None else t_eq
    nt = len(t_eq)
    R = np.linspace(0.84, 2.54, nR)
    Z = np.linspace(-1.6, 1.6, nZ)
    psin = eq.psin(*np.meshgrid(R, Z))  # (nZ, nR)
    psi = eq.psi_axis + psin * (eq.psi_sep - eq.psi_axis)
    nrho = 65
    PSIN = np.linspace(0, 1, nrho)
    g = '\\%s::TOP.RESULTS.gEQDSK.' % tree
    a = '\\%s::TOP.RESULTS.aEQDSK' % tree
    A = '\\%s::TOP.RESULTS.AEQDSK' % tree
    ones = np.ones(nt)
    bdry = np.tile(eq.boundary(), (nt, 1, 1))
    lim = np.loadtxt(__file__.rsplit('/', 2)[0] + '/pytomo/geometry/DIIID/SXR/border.txt')
    sig = {
        g + 'R': R, g + 'Z': Z,
        g + 'GTIME': t_eq * 1e3,
        g + 'PSIRZ': np.tile(psi, (nt, 1, 1)) / (2 * np.pi),
        g + 'RMAXIS': ones * (eq.R0 + eq.shift), g + 'ZMAXIS': ones * eq.Z0,
        g + 'BDRY': bdry,
        g + 'SSIMAG': ones * eq.psi_axis / (2 * np.pi),
        g + 'SSIBRY': ones * eq.psi_sep / (2 * np.pi),
        g + 'CPASMA': ones * 1e6, g + 'BCENTR': ones * -2.0, g + 'RZERO': ones * 1.6955,
        g + 'PSIN': PSIN,
        g + 'QPSI': np.tile(1 + 3 * PSIN**2, (nt, 1)),
        g + 'RHOVN': np.tile(np.sqrt(PSIN), (nt, 1)),
        g + 'FPOL': np.tile(-2.0 * 1.6955 * np.ones(nrho), (nt, 1)),
        g + 'LIM': lim,
        '\\%s::TOP.COMMENTS' % tree: 'synthetic',
        a + '.ATIME': t_eq * 1e3,
        a + '.CHISQ': ones * 10.0,
        A + ':CHISQ': ones * 10.0, A + ':CONDNO': ones, A + ':ERROR': ones * 1e-4,
        A + ':TERROR': ones * 1e-4,
    }
    return sig


def time_trace(t):
    """Plasma 'on' between 0.2 and 5 s, smooth ramps, used for all channels."""
    return np.clip((t - 0.2) / 0.3, 0, 1) * np.clip((5.0 - t) / 0.3, 0, 1)


def sxr_slow_signals(eq, shot, geometry_path, profile='peaked', dt=1e-4, noise=0.002, seed=0):
    """Raw slow poloidal SXR signals [V] PTDATA2("SX90R{P,M}1S{ch}"), inverse of the calibration."""
    from pytomo.geometry.DIIID.SXR.load_SXR import xraypaths, get_calib_fact

    rng = np.random.default_rng(seed)
    tvec = np.arange(-0.2, 6.5, dt)
    trace = time_trace(tvec)
    r1, z1, r2, z2, _ = xraypaths(1e6, toroidal=False)
    bright = line_integrals(eq, r1, z1, r2, z2, PROFILES[profile])
    calib = get_calib_fact(shot, geometry_path, toroidal=False)
    # channel order of xraypaths: RP1 1-32, RM1 1-32
    cal = np.r_[calib['90RP1a'], calib['90RP1b'], calib['90RM1a'], calib['90RM1b']]
    volts = bright / cal
    sig = {}
    for icam, cam in enumerate(('90RP1', '90RM1')):
        for ch in range(32):
            i = 32 * icam + ch
            v = volts[i] * trace + noise * rng.standard_normal(tvec.size)
            name = 'PTDATA2("SX%sS%d",%d)' % (cam, ch + 1, shot)
            sig[name] = v.astype(np.float32)
            sig['dim_of(%s)' % name] = tvec * 1e3
    return sig, dict(tvec=tvec, bright=bright, calib=cal)


def sxr_brightness(eq, profile='peaked'):
    """Noise free poloidal SXR brightness [W/m^2] of the 64 channels at full plasma."""
    from pytomo.geometry.DIIID.SXR.load_SXR import xraypaths
    r1, z1, r2, z2, _ = xraypaths(1e6, toroidal=False)
    return line_integrals(eq, r1, z1, r2, z2, PROFILES[profile])


def synthetic_d3d(shots, geometry_path):
    """Fake DIII-D database; `shots` maps shot number -> emissivity profile name.

    All shots share the same equilibrium (EFIT node names do not contain the shot)."""
    eq = Equilibrium()
    server = FakeServer({})
    install(server)  # load_SXR imports MDSplus at module level
    server.signals.update(efit_signals(eq, min(shots)))
    for shot, profile in shots.items():
        sxr, _ = sxr_slow_signals(eq, shot, geometry_path, profile)
        server.signals.update(sxr)
    return server, eq
