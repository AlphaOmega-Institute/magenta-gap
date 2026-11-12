"""The cubic near-black law at single contacts, and three numbers from the September note.

Inputs   Stockman-Sharpe table; part 1 on 390-700 nm (cubic spline and PCHIP), part 2 on the full 390-830 nm table (PCHIP).
Checks   1. With exact band integrals and no linear programming: deficit / (A z^3) for A = q/(24 b^3) at 510 and 550 nm,
            z from 5e-4 to 4e-2 (white = 1).  One-sided second derivatives are combined when they differ at a knot.
         2. The long contact where L/M is stationary; the short contact of the supporting functional of the purple face;
            the equal-contribution mixture tau* for four candidate signal functionals.
Status   Numerical reproduction.  Part 1 supports the coefficient formula at single contacts; the two-contact mixture
         coefficient is not tested.  Part 2 reproduces 703.417797 nm and tau* = 0.0674706704 from the September note and
         identifies its signal functional as the equal-weight sum of white-balanced channels.
Note     fsolve may warn at the smallest z; the printed ratios are unaffected.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline
from scipy.optimize import brentq, fsolve
wl_all, LMS_all = load_table()

print("== 1. single spectral contact, no bins and no LP: exact band integrals vs A = q / (24 b^3)  [September eq. (1), (12)] ==")
lo, hi = 390.0, 700.0
T = np.array([[2/3, -1/3, -1/3], [-1/3, 2/3, -1/3], [1, 1, 1]])/1.5; W = np.array([0, 0, 2.0])
for kind in ("spline", "pchip"):
    F = CubicSpline(wl_all, LMS_all) if kind == "spline" else PchipInterpolator(wl_all, LMS_all)
    AD = F.antiderivative(); tot = AD(hi) - AD(lo)
    c = lambda lam: (F(lam)/tot) @ T.T
    dc = lambda lam: (F.derivative()(lam)/tot) @ T.T
    I = lambda a, b: ((AD(b) - AD(a))/tot) @ T.T
    ell = lambda v: v[2]/2.0
    for lam0 in (510.0, 550.0):
        n = np.cross(c(lam0), dc(lam0)); n = n/(n @ W); p0 = c(lam0)/ell(c(lam0)); hfun = lambda v: np.linalg.det(np.array([W, p0, v]))
        f = lambda l: float(c(l) @ n); b = ell(c(lam0)); e = 1e-3
        qm, qp = (f(lam0 - 2*e) - 2*f(lam0 - e) + f(lam0))/e**2, (f(lam0 + 2*e) - 2*f(lam0 + e) + f(lam0))/e**2
        A = qm*qp/(6*b**3*(np.sqrt(qm) + np.sqrt(qp))**2)            # reduces to q/(24 b^3) when the two sides agree
        out = []
        for z in (5e-4, 1e-3, 2e-3, 4e-3, 8e-3, 1.6e-2, 4e-2):
            w = z/b
            sol = fsolve(lambda u: [ell(I(u[0], u[1])) - z, hfun(I(u[0], u[1]))/z], [lam0 - w/2, lam0 + w/2], xtol=1e-14)
            x = I(sol[0], sol[1]); out.append((z, float(n @ x), sol[1] - sol[0]))
        txt = ", ".join(f"z={z:g}: {d/(A*z**3):.3f} (band {w:.2f} nm)" for z, d, w in out)
        print(f"[{kind}] {lam0:.0f} nm: A = {A:.5f}; deficit/(A z^3) -> {txt}")

print("\n== 2. the two-contact (purple) face on the full table, 390-830 nm, PCHIP  [September endpoint audit] ==")
F = PchipInterpolator(wl_all, LMS_all); dF, d2F = F.derivative(), F.derivative(2); AD = F.antiderivative()
ratio_slope = lambda l: float(dF(l)[0]*F(l)[1] - F(l)[0]*dF(l)[1])                  # d/dlambda of L/M, up to a positive factor
lam1 = brentq(ratio_slope, 695.0, 715.0, xtol=1e-13)
print(f"long contact: L/M is stationary at {lam1:.6f} nm  [September note: approximately 703.417797 nm]; S there = {float(F(lam1)[2]):.1e}")
tang = lambda l: float(np.linalg.det(np.array([F(l), dF(l), F(lam1)])))              # plane through c(lam0), c'(lam0), c(lam1)
gridv = np.arange(392.0, 470.0, 0.01); tv = np.array([tang(l) for l in gridv]); idx = np.flatnonzero(np.sign(tv[:-1])*np.sign(tv[1:]) < 0)
cands = [brentq(tang, gridv[i], gridv[i+1], xtol=1e-13) for i in idx]
full = np.arange(390.0, 830.0001, 0.01); Call = F(full)
for lam0 in cands:
    n = np.cross(F(lam0), F(lam1)); Wraw = AD(830.0) - AD(390.0); n = n/(n @ Wraw); fall = Call @ n
    supporting = fall.min() > -1e-12*np.abs(fall).max()
    print(f"candidate short contact {lam0:.6f} nm: n.c >= 0 on the whole table: {supporting} (min {fall.min():.2e})")
    if not supporting: continue
    q0, q1 = float(d2F(lam0) @ n), float(d2F(lam1) @ n)
    print(f"   quadratic coefficients q0 = {q0:.6e}, q1 = {q1:.6e} (both positive: {q0 > 0 and q1 > 0})")
    E = np.trapezoid(LMS_all, wl_all, axis=0)
    for name, ell in (("equal-weight sum of white-balanced channels", (1/Wraw)/3), ("raw sum L+M+S", np.ones(3)/Wraw.sum()),
                      ("luminance 0.68990272 L + 0.34832189 M", np.array([0.68990272, 0.34832189, 0.0])/(np.array([0.68990272, 0.34832189, 0.0]) @ Wraw)),
                      ("L + M", np.array([1, 1, 0.0])/(Wraw[0] + Wraw[1]))):
        b0, b1 = float(F(lam0) @ ell), float(F(lam1) @ ell); A0, A1 = q0/(24*b0**3), q1/(24*b1**3); rho = (A0/A1)**(1/3)
        print(f"   signal functional = {name}: A0 = {A0:.6g}, A1 = {A1:.6g}, tau* = {rho/(1 + rho):.10f}")
print("[September note, Stockman-Sharpe / PCHIP: tau* = 0.0674706704]")
