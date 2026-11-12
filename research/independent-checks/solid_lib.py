"""Bounded-reflectance body with a finite number of wavelength bins (helper class used by other scripts).

Inputs   Stockman-Sharpe table, PCHIP-interpolated to the requested bin width on [lo, hi]; each channel divided by its own
         sum so that white is (1, 1, 1).
Model    Reflectances between 0 and 1 on each bin, linear cone catches.  Coordinates: t = 1 + z = (L + M + S)/1.5 and
         (x, y) = (L - h/3, M - h/3)/1.5 with h = L + M + S, so that (x, y)/(1 + z) is the chart (l - 1/3, m - 1/3).
Status   Exact for the finite-bin model, to linear-programming tolerance (absolute errors of order 1e-7 to 1e-6 in R).
Scope    A finite-bin body is a zonotope.  It is not the continuum body: see continuum_rim.py.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

class Solid:
    def __init__(self, step=1.0, lo=390.0, hi=700.0):
        wl_all, LMS_all = load_table()
        self.wl = np.arange(lo, hi + 1e-9, step)
        LMS = np.clip(PchipInterpolator(wl_all, LMS_all)(self.wl), 0, None)
        w = np.full(len(self.wl), step); w[0] = w[-1] = step/2
        gen = LMS*w[:, None]; self.gen = gen/gen.sum(0)                 # white = (1, 1, 1)
        h = self.gen.sum(1); self.N = len(self.wl)
        self.gx, self.gy, self.gz = (self.gen[:, 0] - h/3)/1.5, (self.gen[:, 1] - h/3)/1.5, h/1.5
        self.chrom = self.gen[:, :2]/h[:, None] - 1/3
        self.hull = ConvexHull(self.chrom); self.Nn, self.Cc = self.hull.equations[:, :2], self.hull.equations[:, 2]
        self.A = np.zeros((3, self.N + 1)); self.A[0, :self.N], self.A[1, :self.N], self.A[2, :self.N] = self.gx, self.gy, self.gz
        self.c = np.zeros(self.N + 1); self.c[self.N] = -1
        self.bounds = [(0, 1)]*self.N + [(0, None)]
    def alpha(self, theta):
        u = np.array([np.cos(theta), np.sin(theta)]); den = self.Nn @ u
        return np.min(-self.Cc[den > 1e-15]/den[den > 1e-15])
    def R(self, theta, z, full=False):
        s = float(self.N)                                   # rescale so the coefficients are of order one
        A = self.A*s; A[0, self.N], A[1, self.N] = -np.cos(theta)*s, -np.sin(theta)*s
        b = [0.0, 0.0, (z + 1)*s]
        for method in ("highs", "highs-ds", "highs-ipm"):
            res = linprog(self.c, A_eq=A, b_eq=b, bounds=self.bounds, method=method)
            if res.x is not None and res.status == 0: break
        else:
            raise RuntimeError(f"LP failed at theta={theta}, z={z}: {res.message}")
        return (res.x[self.N], res.x[:self.N]) if full else res.x[self.N]
    def kind(self, rho):
        b = (rho > 0.5).astype(int); edges = np.flatnonzero(np.diff(b)) ; tw = [float(self.wl[i]) for i in edges]
        name = {0: "white" if b[0] else "black", 1: "long-pass" if b[-1] else "short-pass", 2: "band-pass" if not b[0] else "band-stop"}.get(len(edges), f"{len(edges)} transitions")
        return name, tw
