"""Reflectance type at the widest point, and how the slope jump there scales with step size (0.25 nm bins).

Inputs   Stockman-Sharpe table via solid_lib.Solid, 390-700 nm, 1 nm and 0.25 nm bins.
Checks   In the 1 nm model R(theta, .) is concave and piecewise linear.  At 72 hues in the 0.25 nm model: the optimal
         reflectance at the widest point (number of transitions), and the jump in dR/dz across it at steps 0.04 to 0.005.
Status   Numerical evidence only.  One-transition points keep their jump as the step shrinks; two-transition points lose it
         in proportion.  The continuum statement is made in continuum_rim.py, and one point is certified in
         research/object-colour/Continuum_Rim_Certificate.md.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, time, numpy as np
from scipy.optimize import minimize_scalar
from solid_lib import Solid
coarse, fine = Solid(1.0), Solid(0.25)
print(f"bins: coarse {coarse.N}, fine {fine.N}")
# Point raised in review, checked directly on the 1 nm model: R(theta, .) is concave and piecewise linear
t = np.radians(100.0); zz = np.linspace(-0.2, 0.2, 401); Rz = np.array([coarse.R(t, z) for z in zz])
d2 = np.diff(Rz, 2); print(f"1 nm model, hue 100 deg: second differences of R in z are <= 0 everywhere: {np.all(d2 < 1e-10)}; "
                           f"exactly zero (straight pieces) at {np.mean(np.abs(d2) < 1e-11):.0%} of grid points")
ths = np.linspace(-np.pi, np.pi, 72, endpoint=False); rows = []; t0 = time.time()
for t in ths:
    zc = np.linspace(-0.9, 0.9, 37); i = int(np.argmax([coarse.R(t, z) for z in zc]))
    z0 = minimize_scalar(lambda z: -fine.R(t, z), bounds=(zc[max(i-1, 0)], zc[min(i+1, 36)]), method="bounded", options={"xatol": 2e-5}).x
    R0, rho = fine.R(t, z0, full=True); name, tw = fine.kind(rho)
    J = [(2*R0 - fine.R(t, z0 - d) - fine.R(t, z0 + d))/d for d in (0.04, 0.02, 0.01, 0.005)]      # jump in dR/dz across the widest point
    rows.append((np.degrees(t), z0, name, tw, J, fine.alpha(t), fine.alpha(t + np.pi)))
print(f"{len(rows)} hues in {time.time() - t0:.0f}s\n")
print("hue    z_widest  reflectance at the widest point            slope jump in z at steps 0.04 / 0.02 / 0.01 / 0.005    alpha(theta)+alpha(theta+pi)")
for r in rows[::3]:
    print(f"{r[0]:+5.0f}  {r[1]:+.3f}    {r[2]:10s} {str([round(x, 1) for x in r[3]]):22s}  {r[4][0]:.3f} / {r[4][1]:.3f} / {r[4][2]:.3f} / {r[4][3]:.3f}                 {r[5] + r[6]:.3f}")
J = np.array([r[4] for r in rows]); ratio = J[:, 3]/J[:, 0]; kinds = np.array([r[2] for r in rows])
persist = ratio > 0.6           # jump does not shrink with the step: a corner at this resolution
shrink = ratio < 0.25           # jump shrinks roughly in proportion to the step (0.125 expected for a smooth maximum)
print(f"\nslope jump persists as the step shrinks 8-fold (ratio > 0.6): {persist.sum()} hues; shrinks (ratio < 0.25): {shrink.sum()} hues; in between: {(~persist & ~shrink).sum()} hues")
for k in sorted(set(kinds)):
    sel = kinds == k; print(f"  {k:11s}: {sel.sum():2d} hues, of which jump persists at {np.sum(persist & sel)} and shrinks at {np.sum(shrink & sel)}; median jump ratio {np.median(ratio[sel]):.2f}")
near_end = np.array([any(abs(x - 390) < 3 or abs(x - 700) < 3 for x in r[3]) for r in rows])
print(f"hues whose widest-point reflectance has a transition within 3 nm of an end of the spectrum: {near_end.sum()}; jump persists at {np.sum(persist & near_end)} of them")
print(f"jump persists but no transition near a spectrum end: {np.sum(persist & ~near_end)} hues, at angles {[int(round(r[0])) for r, p, n in zip(rows, persist, near_end) if p and not n]}")
