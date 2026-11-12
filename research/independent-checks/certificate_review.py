"""Independent review of research/object-colour/certify_continuum_rim.py.

Inputs   The certificate script and its CSV input, as committed on main.
Checks   The exact rational PCHIP agrees with SciPy's on sampled cells; the Bernstein sign certifier refuses a polynomial
         that changes sign inside a cell and one that touches zero, and certifies a positive one needing subdivision; a
         floating-point reproduction of the two roots and the support point falls inside the certified enclosures; white is
         a positive combination of the two transition vectors.
Status   Independent floating-point reproduction.  The rigorous statement is the certificate itself.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, csv, importlib.util, numpy as np
from fractions import Fraction as F
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
repo = str(REPO)
spec = importlib.util.spec_from_file_location("cert", f"{repo}/research/object-colour/certify_continuum_rim.py"); cert = importlib.util.module_from_spec(spec)
sys.modules["cert"] = cert; spec.loader.exec_module(cert)

print("== inputs ==")
rows = [r for r in csv.reader(open(f"{repo}/research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv")) if r]
try: float(rows[0][0])
except ValueError: rows = rows[1:]
tab = np.array([[float(v) for v in r] for r in rows]); wl, LMS = load_table()
print("the shared loader returns the CSV values exactly (441 rows, 390-830 nm):", tab.shape == (441, 4) and np.array_equal(tab[:, 0], wl) and np.array_equal(tab[:, 1:], LMS))

print("\n== the two core routines ==")
P = PchipInterpolator(wl, LMS); worst = 0.0
for j in range(3):
    exact = cert.pchip([F.from_float(float(v)) for v in LMS[:, j]])
    for i in (0, 1, 96, 97, 209, 210, 225, 226, 300, 438, 439):
        for t in (0.0, 0.3, 0.77, 1.0):
            v = float(sum(a*F(t)**k for k, a in enumerate(exact[i]))); worst = max(worst, abs(v - float(P(390 + i + t)[j])))
print(f"exact rational PCHIP agrees with SciPy's PCHIP on sampled cells (including both ends): worst difference {worst:.1e}")
IV = cert.IV; mk = lambda cs: [IV.point(F(c)) for c in cs]
def try_sign(coeffs, wanted):
    try: cert.certify_sign(cert.bernstein(mk(coeffs)), wanted); return "certified"
    except AssertionError: return "refused"
print("sign certifier on (t-0.4)(t-0.6), positive at both ends but negative inside:", try_sign([F(24, 100), F(-1), F(1)], 1))
print("sign certifier on (t-0.5)^2, touching zero inside:", try_sign([F(1, 4), F(-1), F(1)], 1))
print("sign certifier on (t-0.4)(t-0.6) + 0.011, positive everywhere but needing subdivision:", try_sign([F(24, 100) + F(11, 1000), F(-1), F(1)], 1))

print("\n== independent floating-point reproduction on the full 390-830 nm table ==")
AD = P.antiderivative(); tot = AD(830.0) - AD(390.0)
c = lambda l: P(l)/tot
q = lambda l: float(-c(l)[0] + 2*c(l)[1] - c(l)[2])
grid = np.arange(390.0, 830.0001, 0.002); qg = np.array([-1, 2, -1.0]) @ (P(grid)/tot).T
idx = np.flatnonzero(np.sign(qg[:-1])*np.sign(qg[1:]) < 0); roots = [brentq(q, grid[i], grid[i+1], xtol=1e-13) for i in idx]
print("sign changes of q on a 0.002 nm grid:", len(idx), "| roots:", ", ".join(f"{r:.11f}" for r in roots), " [certificate: 486.863080975890..948 and 599.369447760400..458]")
a = (AD(roots[1]) - AD(roots[0]))/tot; tot3 = a.sum()
xyz = np.array([(2/3)*(a[0] - tot3/3), (2/3)*(a[1] - tot3/3), (2/3)*tot3 - 1])
enc = [(0.11649647725622, 0.11649647725683), (0.20194541250715, 0.20194541250768), (0.09328313559491, 0.09328313559570)]
print("support point:", ", ".join(f"{v:.14f}" for v in xyz), "| inside the certified enclosures:", all(lo_ <= v <= hi_ for v, (lo_, hi_) in zip(xyz, enc)))
print(f"hue of this point: {np.degrees(np.arctan2(xyz[1], xyz[0])):.3f} deg; q is positive exactly on the band: {bool(np.all(qg[(grid > roots[0]) & (grid < roots[1])] > 0) and np.all(qg[(grid < roots[0]) | (grid > roots[1])] < 0))}")

print("\n== what 'horizontal' means without coordinates ==")
A = np.column_stack([c(roots[0]), c(roots[1])]); coef, res, *_ = np.linalg.lstsq(A, np.ones(3), rcond=None)
print(f"white = {coef[0]:.3f} c(lambda1) + {coef[1]:.3f} c(lambda2), residual {np.linalg.norm(A @ coef - 1):.1e}; both weights positive: {bool(np.all(coef > 0))}")
