"""Width of the supporting-normal fan at edge-colour rim points, and the 550 nm contact at several scales.

Inputs   Stockman-Sharpe table; PCHIP (and, for part 2, a cubic spline) on 390-700 nm, no bins.
Checks   1. At edge-colour rim points, the angular range of normals n + e m (m perpendicular to c(lambda1) and n) that
            support the same reflectance over the whole spectrum.
         2. At 510 and 550 nm: the tangent functional is supporting nearby; one-sided second derivatives; 2 f/d^2 from
            0.25 to 8 nm.
Status   Numerical.
Correction (agreed in review)  550 nm is a supporting contact with a positive, scale-stable quadratic coefficient.  An
         earlier explanation of the 550 nm exponents ("nearly straight, so quadratic contact fails") was wrong.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline
from scipy.optimize import brentq
wl_all, LMS_all = load_table()
lo, hi = 390.0, 700.0
T = np.array([[2/3, -1/3, -1/3], [-1/3, 2/3, -1/3], [1, 1, 1]])/1.5

def model(kind):
    F = PchipInterpolator(wl_all, LMS_all) if kind == "pchip" else CubicSpline(wl_all, LMS_all)
    AD = F.antiderivative(); tot = AD(hi) - AD(lo)
    c = lambda lam: (F(lam)/tot) @ T.T
    dc = lambda lam, k=1: (F.derivative(k)(lam)/tot) @ T.T
    return c, dc

print("== 1. angular width of the set of supporting normals at edge-colour rim points (PCHIP continuum) ==")
c, dc = model("pchip"); grid = np.arange(lo, hi + 1e-9, 0.005); Cg = c(grid)
def pattern_ok(n, lam1, longpass):
    f = Cg @ n; tol = 1e-12*np.abs(f).max(); up = grid > lam1
    return bool(np.all(f[up] >= -tol) and np.all(f[~up] <= tol)) if longpass else bool(np.all(f[up] <= tol) and np.all(f[~up] >= -tol))
widths = []
for phi in np.linspace(0, 2*np.pi, 360, endpoint=False):
    n = np.array([np.cos(phi), np.sin(phi), 0.0]); f = Cg @ n; idx = np.flatnonzero(np.sign(f[:-1])*np.sign(f[1:]) < 0)
    if len(idx) != 1: continue
    lam1 = brentq(lambda l: float(c(l) @ n), grid[idx[0]], grid[idx[0]+1], xtol=1e-12); longpass = f[-1] > 0
    m = np.cross(c(lam1), n); m /= np.linalg.norm(m)
    def limit(sign):
        a, b = 0.0, 50.0
        if pattern_ok(n + sign*b*m, lam1, longpass): return b
        for _ in range(50):
            mid = (a + b)/2
            if pattern_ok(n + sign*mid*m, lam1, longpass): a = mid
            else: b = mid
        return a
    ep, em = limit(+1), limit(-1)
    na, nb = n + ep*m, n - em*m; ang = np.degrees(np.arccos(np.clip(na @ nb/np.linalg.norm(na)/np.linalg.norm(nb), -1, 1)))
    widths.append((lam1, ang))
widths = np.array(widths)
print(f"edge-colour rim points sampled: {len(widths)}; width of the normal fan: {widths[:, 1].min():.1f} .. {widths[:, 1].max():.1f} deg (median {np.median(widths[:, 1]):.1f}); every one above 1 deg: {bool(np.all(widths[:, 1] > 1))}")

print("\n== 2. the 550 nm contact ==")
for kind in ("pchip", "spline"):
    c, dc = model(kind)
    for lam0 in (510.0, 550.0):
        n = np.cross(c(lam0), dc(lam0)); W = np.array([0, 0, 2.0]); n = n/(n @ W)       # tangent plane of the cone at c(lam0), n(W) = 1
        f = lambda l: float(c(l) @ n)
        b = c(lam0)[2]/2.0                                                               # signal density, white = 1
        eps = 1e-4; q_minus = (f(lam0 - 2*eps) - 2*f(lam0 - eps) + f(lam0))/eps**2; q_plus = (f(lam0 + 2*eps) - 2*f(lam0 + eps) + f(lam0))/eps**2
        qs = ", ".join(f"{2*f(lam0 + d)/d**2:.3e}" for d in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0))
        qm = ", ".join(f"{2*f(lam0 - d)/d**2:.3e}" for d in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0))
        near = np.arange(lam0 - 8, lam0 + 8.001, 0.01); fmin = min(f(l) for l in near)
        print(f"[{kind}] {lam0:.0f} nm: supporting nearby (min of n.c within 8 nm = {fmin:.1e}); one-sided second derivatives {q_minus:.3e} / {q_plus:.3e}")
        print(f"         2 f/d^2 at d = 0.25, 0.5, 1, 2, 4, 8 nm above: {qs}")
        print(f"                                             below: {qm}")
