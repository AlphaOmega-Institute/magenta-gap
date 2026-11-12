"""A boundary face through black: the reviewer's counterexample, verified.

Inputs   Stockman-Sharpe table; PCHIP on 390-830 nm.
Checks   S is identically zero on 616-830 nm in the table and the interpolant; half reflectance there gives the catch vector
         (0.06566, 0.01201, 0); it is a boundary point with tangent plane S = 0, which contains its own direction from black.
         At the certified point both intensity scaling and motion toward white cross the tangent plane.
Status   Numerical verification of an exact counterexample found in review.
Correction  It refutes the claim "no boundary point other than black has a tangent plane containing its direction from
         black".  That claim holds only when no supporting functional vanishes on the response vectors over a stretch of
         wavelengths.  The face exists here because S is tabulated as zero; this is a property of the table model.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, numpy as np
from scipy.interpolate import PchipInterpolator
wl, LMS = load_table()
P = PchipInterpolator(wl, LMS); AD = P.antiderivative(); tot = AD(830.0) - AD(390.0)
I = lambda a, b: (AD(b) - AD(a))/tot                                  # white-balanced catches, white = (1, 1, 1)
fine = np.arange(616.0, 830.0001, 0.01)
print("S on 616-830 nm: table max", LMS[wl >= 616, 2].max(), "| interpolant max", float(np.abs(P(fine)[:, 2]).max()))
print("smallest tabulated L, M over the whole table:", LMS[:, 0].min(), LMS[:, 1].min(), "| S is zero only from", wl[LMS[:, 2] == 0].min(), "nm")
a = 0.5*I(616.0, 830.0); print("half reflectance on 616-830 nm, none elsewhere: a =", np.round(a, 5), " [value reported in review: (0.06566, 0.01201, 0)]")
n = np.array([0.0, 0.0, -1.0])                                         # outward normal of the plane S = 0
grid = np.arange(390.0, 830.0001, 0.01); f = (P(grid)/tot) @ n
print("n.c <= 0 on the whole spectrum:", bool(np.all(f <= 0)), "| n.c = 0 on a stretch of", f"{np.mean(f == 0)*440:.0f} nm", "| n.a =", float(n @ a), "= support value 0, so a is a boundary point")
d1, d2 = P(630.0)/tot, P(760.0)/tot; rank = np.linalg.matrix_rank(np.array([d1, d2]), tol=1e-12)
print("response vectors at 630 and 760 nm are independent (the face is two-dimensional):", rank == 2, "| so half reflectance is in its relative interior and the tangent plane there is S = 0")
print("the intensity-scaling direction a lies in that plane:", float(n @ a) == 0.0, "| complement W - a:", np.round(1 - a, 5), "with normal +S: (W - a) direction toward white has S-component", float((1 - (1 - a))[2]))
# at the certified point of the rim certificate, both directions cross the tangent plane
x = I(486.86308097592, 599.36944776045); ny = np.array([-1.0, 2.0, -1.0])     # the certified normal, in catch coordinates
print(f"certified point: n.x = {ny @ x:.6f} (from black), n.(W - x) = {ny @ (1 - x):.6f} (toward white); both nonzero, so neither direction is tangent there")
