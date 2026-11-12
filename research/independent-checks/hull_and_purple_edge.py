"""Attainable chromaticity hull and its long mixture edge, in the September manuscript's conventions.

Inputs   Stockman-Sharpe 2-degree table (research/object-colour/inputs), 390-830 nm; stimuli restricted to 390-700 nm.
Checks   The convex hull of the 311 chromaticities in both conventions of Sec. 5 (white-balanced channels; raw channels
         centred on equal-energy white): vertex count and the edges joining non-adjacent samples.  For the 410-700 nm edge:
         distance from white, normal direction, angle subtended, and that the manuscript's edge formula (eq:edge-derivative)
         equals tan(theta - theta0).
Status   Numerical reproduction of numbers stated in the September manuscript (188 vertices; edges 410-700, 453-458,
         599-616, 616-700), plus an exact identity checked numerically.
Note     Earlier figures quoted for a 390-700 nm chord with raw channels centred at 1/3 used the historical pipeline's
         convention, which Sec. 5 rejects.  That script is not preserved.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, numpy as np
from scipy.spatial import ConvexHull
wl_all, LMS_all = load_table()
print("input grid:", wl_all.min(), "-", wl_all.max(), "nm,", len(wl_all), "samples")
E = np.trapezoid(LMS_all, wl_all, axis=0); print("E integrals:", np.round(E, 6))
m = (wl_all >= 390) & (wl_all <= 700); wl = wl_all[m]

def setup(balanced):
    X = LMS_all[m] / E if balanced else LMS_all[m].copy()
    white3 = np.ones(3) if balanced else E
    c = X / X.sum(1, keepdims=True); w = white3 / white3.sum()
    pts = c[:, :2] - w[:2]
    return X, white3, pts, ConvexHull(pts)

for name, balanced in (("E-balanced", True), ("raw, E-centred", False)):
    X, white3, pts, hull = setup(balanced)
    v = hull.vertices                                   # counter-clockwise
    jumps = [(wl[a], wl[b]) for a, b in zip(v, np.roll(v, -1)) if abs(a - b) > 1]
    print(f"\n[{name}] hull vertices: {len(v)}; edges joining non-adjacent samples (nm):", sorted((int(min(a, b)), int(max(a, b))) for a, b in jumps))
    a, b = pts[wl == 700][0], pts[wl == 410][0]         # the long mixture edge of the manuscript
    d = b - a; n = np.array([d[1], -d[0]]); n /= np.linalg.norm(n); dist = n @ a
    if dist < 0: n, dist = -n, -dist
    th0 = np.arctan2(n[1], n[0])
    ang = np.arccos(a @ b / np.linalg.norm(a) / np.linalg.norm(b))
    t = np.linspace(1e-5, 1 - 1e-5, 101); p = a[None] + t[:, None]*d[None]
    slope_paper = (p @ d) / (p[:, 0]*d[1] - p[:, 1]*d[0])            # manuscript Eq. (edge-derivative)
    slope_tan = np.tan(np.arctan2(p[:, 1], p[:, 0]) - th0)
    print(f"  410-700 nm edge: distance from white {dist:.4f}, normal at {np.degrees(th0):.1f} deg, subtends {np.degrees(ang):.1f} deg at white")
    print(f"  manuscript edge formula equals tan(theta - theta0): {np.allclose(slope_paper, slope_tan)}; slope range {slope_paper.min():.3f} .. {slope_paper.max():.3f}")
    r = np.hypot(*pts[v].T); print(f"  hull radius range {r.min():.4f} .. {r.max():.4f}")
