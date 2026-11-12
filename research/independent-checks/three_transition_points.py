"""Three-transition rim points on the full table, and weighting by normals versus hue angle.

Inputs   Stockman-Sharpe table; PCHIP on 390-830 nm, no bins.
Checks   Rim normals with three interior sign changes: two of the transition vectors are independent and their cross product
         is the horizontal normal, so the criterion (which needs at least two independent transitions) applies.  Shares of
         each class by count of normals and by hue angle of the boundary point.
Status   Numerical, on 1,440 sampled normals.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
wl_all, LMS_all = load_table(); lo, hi = 390.0, 830.0
P = PchipInterpolator(wl_all, LMS_all); AD = P.antiderivative(); tot = AD(hi) - AD(lo)
T = np.array([[2/3, -1/3, -1/3], [-1/3, 2/3, -1/3], [1, 1, 1]])/1.5
c = lambda lam: (P(lam)/tot) @ T.T
grid = np.arange(lo, hi + 1e-9, 0.005); Cg = c(grid)
phis = np.linspace(0, 2*np.pi, 1440, endpoint=False); th, nz, three = [], [], []
for phi in phis:
    n = np.array([np.cos(phi), np.sin(phi), 0.0]); f = Cg @ n; idx = np.flatnonzero(np.sign(f[:-1])*np.sign(f[1:]) < 0)
    zs = [brentq(lambda l: float(c(l) @ n), grid[i], grid[i+1], xtol=1e-12) for i in idx]
    x = np.zeros(3); edges = [lo] + zs + [hi]
    for a, b in zip(edges[:-1], edges[1:]):
        if float(c((a + b)/2) @ n) > 0: x += ((AD(b) - AD(a))/tot) @ T.T
    th.append(np.arctan2(x[1], x[0])); nz.append(len(zs))
    if len(zs) == 3:
        v = [c(l) for l in zs]; best = 0.0
        for i in range(3):
            for j in range(i + 1, 3):
                cr = np.cross(v[i], v[j]); s = np.linalg.norm(cr)/(np.linalg.norm(v[i])*np.linalg.norm(v[j]))
                if s > best: best, par = s, abs(abs(cr @ n)/np.linalg.norm(cr) - 1)
        three.append((zs, best, par, [float(np.linalg.norm(u)) for u in v]))
th = np.unwrap(th); nz = np.array(nz); dth = np.abs(np.diff(np.append(th, th[0] + np.sign(th[-1] - th[0])*2*np.pi)))
print("three-transition rim normals:", len(three))
print(f"  transition wavelengths: first {min(t[0][0] for t in three):.1f}-{max(t[0][0] for t in three):.1f}, second {min(t[0][1] for t in three):.1f}-{max(t[0][1] for t in three):.1f}, third {min(t[0][2] for t in three):.1f}-{max(t[0][2] for t in three):.1f} nm")
print(f"  best pair of transition vectors: sine of angle {min(t[1] for t in three):.3f} .. {max(t[1] for t in three):.3f}; their cross product is parallel to the horizontal normal to within {max(t[2] for t in three):.1e}")
print("weighting: by count of normals vs by hue angle of the boundary point")
for k in sorted(set(nz)): print(f"  {k} transition(s): {np.mean(nz == k):.1%} of normals, {dth[nz == k].sum()/(2*np.pi):.1%} of hue angle")
