"""Crease slopes at edge colours, and the near-black behaviour of finite bins.

Inputs   Stockman-Sharpe table via solid_lib.Solid, 390-700 nm.
Checks   1. At a long-pass or short-pass colour with transition lambda1, the one-sided slopes of R in z predicted by the two
            planes spanned by c(lambda1) and each end of the spectrum, against slopes measured on the 0.25 nm body.
         2. 1 - R/(alpha t) near black for 1 nm and 0.05 nm bins.
         3. Local exponents of the deficit alpha t - R at five spectral directions (0.05 nm bins).
Status   Numerical.
Corrections (agreed in review)
         * Part 3 is limited by linear-programming error once the deficit is below about 1e-6.  The exponents printed for
           550 nm are solver noise, not a failure of quadratic contact.  cubic_law_and_september_numbers.py repeats the test
           with exact integrals.
         * Part 2: exact contact with the cone at small signal is a bin artifact.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, time, numpy as np
from solid_lib import Solid

print("== 1. crease along the one-transition (edge) colours: predicted from tangent planes vs measured ==")
fine = Solid(0.25)
T = np.array([[2/3, -1/3, -1/3], [-1/3, 2/3, -1/3], [1, 1, 1]])/1.5            # (L, M, S) catches -> (x, y, t = 1 + z)
c = lambda lam: T @ fine.gen[int(np.argmin(np.abs(fine.wl - lam)))]
def slope(n, theta):                                                           # dR/dz of the plane with normal n along the hue section
    u = np.array([np.cos(theta), np.sin(theta), 0.0]); return -n[2]/(n @ u)
from scipy.optimize import minimize_scalar
for hue in (10, 20, 30, 40, -170, -150):
    t = np.radians(hue); zc = np.linspace(-0.9, 0.9, 73); i = int(np.argmax([fine.R(t, z) for z in zc]))
    z0 = minimize_scalar(lambda z: -fine.R(t, z), bounds=(zc[i-1], zc[i+1]), method="bounded", options={"xatol": 2e-5}).x
    R0, rho = fine.R(t, z0, full=True); name, tw = fine.kind(rho); lam1 = tw[0]
    d = 0.005; left, right = (R0 - fine.R(t, z0 - d))/d, (fine.R(t, z0 + d) - R0)/d
    s_red, s_vio = slope(np.cross(c(lam1), c(700.0)), t), slope(np.cross(c(lam1), c(390.0)), t)
    print(f"hue {hue:+4d}: {name} at {lam1:.1f} nm; measured slopes below/above {left:+.4f} / {right:+.4f}; "
          f"planes through c(lam1) and the red end / violet end give {s_red:+.4f} / {s_vio:+.4f}")

print("\n== 2. near black: finite bins sit on the cone, fine bins fall inside it ==")
for step in (1.0, 0.05):
    s = Solid(step); lam0 = 520.0; k = int(np.argmin(np.abs(s.wl - lam0))); th = np.arctan2(s.chrom[k, 1], s.chrom[k, 0]); a = s.alpha(th)
    print(f"bin width {step} nm, hue of {lam0:.0f} nm:", ", ".join(f"t={t}: 1 - R/(alpha t) = {1 - s.R(th, t - 1)/(a*t):.2e}" for t in (0.002, 0.008, 0.032)))

print("\n== 3. the cubic law at spectral directions (0.05 nm bins, 390-700 nm) ==")
s = Solid(0.05); ts = np.array([0.01, 0.02, 0.04, 0.08]); t0 = time.time()
for lam0 in (470.0, 490.0, 510.0, 530.0, 550.0):
    k = int(np.argmin(np.abs(s.wl - lam0))); th = np.arctan2(s.chrom[k, 1], s.chrom[k, 0]); a = s.alpha(th)
    on_hull = np.isclose(a, np.hypot(*s.chrom[k]), rtol=1e-9)
    d = []; width = []
    for t in ts:
        R, rho = s.R(th, t - 1, full=True); d.append(a*t - R); width.append((rho > 0.5).sum()*0.05)
    d = np.array(d); ex = np.log2(d[1:]/d[:-1])
    print(f"{lam0:.0f} nm (hull vertex: {on_hull}): deficit {', '.join(f'{x:.2e}' for x in d)}; exponents {', '.join(f'{x:.2f}' for x in ex)}; band widths {', '.join(f'{w:.1f}' for w in width)} nm")
print(f"({time.time() - t0:.0f}s)")
