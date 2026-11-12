"""Two ways of identifying the manuscript's metric with a discrimination tensor.

Inputs   Stockman-Sharpe table; hull of the 390-700 nm chromaticities in both normalizations (September manuscript, Sec. 5).
Checks   1. The manuscript's test, in which only the angular coefficient B is identified with the simplex Fisher metric:
            the scale c* above which it holds.
         2. A full-block identification, in which the whole chromatic block is c G: nullity then needs G^-1(dF, dF) to be
            one constant along the boundary, so its max/min ratio would have to be 1 for any c.
Status   Numerical.  Part 1 reproduces the c* values of the September table (0.3761, 0.4028).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import numpy as np
from scipy.spatial import ConvexHull
WL, LMS = load_table()
cross = lambda a, b: a[..., 0]*b[..., 1] - a[..., 1]*b[..., 0]
for balanced in (True, False):
    sens = LMS/np.trapezoid(LMS, WL, axis=0) if balanced else LMS
    catches = np.trapezoid(sens, WL, axis=0); white = catches/catches.sum()
    p = sens[:311]/sens[:311].sum(axis=1)[:, None]
    ids = ConvexHull(p[:, :2]).vertices; p = p[ids]
    t = np.tile(np.linspace(1e-5, 1-1e-5, 101), len(p)); seg = np.repeat(np.arange(len(p)), 101)
    q = p[seg]*(1-t[:, None]) + p[(seg+1) % len(p)]*t[:, None]
    d = p[(seg+1) % len(p), :2] - p[seg, :2]
    xy = q[:, :2] - white[:2]; r = np.hypot(xy[:, 0], xy[:, 1]); th = np.arctan2(xy[:, 1], xy[:, 0])
    dlog = (xy*d).sum(axis=1)/cross(xy, d)                       # signed d log R / d theta
    rhat = xy/r[:, None]; that = np.column_stack([-rhat[:, 1], rhat[:, 0]])
    n = rhat - dlog[:, None]*that                                # conormal of F = r - R(theta)
    l, m = q[:, 0], q[:, 1]
    with np.errstate(divide='ignore', invalid='ignore'):
        beta = np.sqrt(np.sin(th)**2/q[:, 0] + np.cos(th)**2/q[:, 1] + (np.sin(th)-np.cos(th))**2/q[:, 2])
    keep = np.all(q > 0, axis=1) & np.isfinite(beta) & (beta < 10)
    # inverse simplex Fisher in (l,m) = single-sample multinomial covariance
    Sig = np.stack([np.stack([l*(1-l), -l*m], -1), np.stack([-l*m, m*(1-m)], -1)], -2)
    i0 = np.flatnonzero(keep)[0]
    G0 = np.array([[1/l[i0]+1/q[i0,2], 1/q[i0,2]], [1/q[i0,2], 1/m[i0]+1/q[i0,2]]])
    assert np.allclose(np.linalg.inv(G0), Sig[i0])
    Q = np.einsum('ni,nij,nj->n', n, Sig, n)[keep]
    c_star = np.max((np.abs(dlog[keep])/beta[keep])**2)
    name = 'E-balanced' if balanced else 'raw'
    print(f"{name:10s} B-only test: holds for c > {c_star:.4f}  (scale-dependent)")
    print(f"{'':10s} full-block test: G^-1(dF,dF) ranges {Q.min():.4f} .. {Q.max():.4f}, max/min = {Q.max()/Q.min():.1f}"
          f"  (must be 1 for any c)")
    for lo, hi in ((0, 60), (60, 120), (120, 180), (-180, -120), (-120, -60), (-60, 0)):
        sel = (np.degrees(th[keep]) >= lo) & (np.degrees(th[keep]) < hi)
        if sel.any(): print(f"{'':12s}hue {lo:5d}..{hi:4d} deg: {Q[sel].min():.4f} .. {Q[sel].max():.4f}")
