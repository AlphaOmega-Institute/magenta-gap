"""The bounded-reflectance body with 1 nm bins, compared with the manuscript's bicone.

Inputs   Stockman-Sharpe table, 311 bins on 390-700 nm, white-balanced over that range.
Checks   Central symmetry; failure of the mirror form; containment in the order-interval body and in the bicone; volumes;
         height of the widest point at each hue; approach to the cone at the black pole.
Status   Numerical, for this finite-bin model and this height coordinate (cone sum).
Corrections (agreed in review)
         * Every fixed-hue section of a finite-bin body is a polygon, so its maximum is a corner or a plateau.  The slope
           columns printed here must not be read as evidence of a smooth maximum.  See widest_point_types.py and
           continuum_rim.py.
         * The ratio 1.00 near black is a property of the bins: a small signal can sit exactly on one spectral direction.
           In the continuum the body falls inside the cone by a cubic deficit (research/object-colour, September note).
         * This is the standard object-colour construction; nothing here is new, and it describes attainable signals, not a
           perceptual boundary.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, time, numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull
wl_all, LMS_all = load_table()
m = (wl_all >= 390) & (wl_all <= 700); wl, LMS = wl_all[m], LMS_all[m]
wts = np.gradient(wl); wts[0] = wts[-1] = 0.5                      # trapezoid weights
gen = LMS * wts[:, None]; gen = gen / gen.sum(0)                   # catches of each bin; white = (1, 1, 1)
h = gen.sum(1); gx, gy, gz = (gen[:, 0] - h/3)/1.5, (gen[:, 1] - h/3)/1.5, h/1.5
N = len(wl)
def R(theta, z):
    A = np.zeros((3, N + 1)); A[0, :N], A[1, :N], A[2, :N] = gx, gy, gz
    A[0, N], A[1, N] = -np.cos(theta), -np.sin(theta)
    c = np.zeros(N + 1); c[N] = -1
    res = linprog(c, A_eq=A, b_eq=[0, 0, z + 1], bounds=[(0, 1)]*N + [(0, None)], method="highs")
    return res.x[N]
# the manuscript's alpha(theta) for the same stimuli: hull of the bins' chromaticities
chrom = gen / h[:, None]; hull = ConvexHull(chrom[:, :2] - 1/3); Nn, Cc = hull.equations[:, :2], hull.equations[:, 2]
def alpha(theta):
    u = np.array([np.cos(theta), np.sin(theta)]); den = Nn @ u; return np.min(-Cc[den > 1e-15]/den[den > 1e-15])

t0 = time.time(); nth, nz = 72, 41
ths = np.linspace(-np.pi, np.pi, nth, endpoint=False); zs = np.linspace(-0.95, 0.95, nz)
Rg = np.array([[R(t, z) for z in zs] for t in ths]); print(f"{nth*nz} radial values in {time.time() - t0:.0f}s")
a = np.array([alpha(t) for t in ths]); ao = np.roll(a, nth//2)
central = np.max(np.abs(Rg - np.roll(Rg, nth//2, axis=0)[:, ::-1])/Rg)
mirror = np.abs(Rg - Rg[:, ::-1])/np.maximum(Rg, Rg[:, ::-1])
print(f"central symmetry R(theta,z) = R(theta+pi,-z): worst relative deviation {central:.1e}")
print(f"mirror symmetry R(theta,z) = R(theta,-z): median deviation {np.median(mirror):.0%}, worst {mirror.max():.0%}")
bic = a[:, None]*(1 - np.abs(zs))[None]; oi = np.minimum(a[:, None]*(1 + zs)[None], ao[:, None]*(1 - zs)[None])
print(f"inside the order-interval body everywhere: {np.all(Rg <= oi*(1 + 1e-7))}; inside the manuscript's bicone everywhere: {np.all(Rg <= bic*(1 + 1e-7))}")
vol = lambda Rr: (0.5*Rr**2).sum()
print(f"volume relative to the manuscript's bicone (same alpha): {vol(Rg)/vol(bic):.2f}; relative to the order-interval body: {vol(Rg)/vol(oi):.2f}; share of the body outside the bicone: {(0.5*np.maximum(Rg**2 - bic**2, 0)).sum()/vol(Rg):.0%}")
# where is the widest point at each hue, and is it a corner or a smooth maximum?
zf = np.linspace(-0.9, 0.9, 181); out = []
for t in ths[::6]:
    Rz = np.array([R(t, z) for z in zf]); i = int(np.argmax(Rz)); z0 = zf[i]
    sl = [(R(t, z0) - R(t, z0 - d))/d for d in (0.05, 0.02)]; sr = [(R(t, z0 + d) - R(t, z0))/d for d in (0.05, 0.02)]
    pred = (alpha(t + np.pi) - alpha(t))/(alpha(t + np.pi) + alpha(t))
    out.append((np.degrees(t), z0, Rz[i], pred, sl[0], sr[0], sl[1], sr[1], alpha(t), alpha(t + np.pi)))
out = np.array(out)
print("\nhue(deg)  z of widest point  [order-interval value]  slope dR/dz just below / just above (step 0.05 | step 0.02)   cone slopes alpha(theta), alpha(theta+pi)")
for r in out:
    print(f"{r[0]:7.0f}   {r[1]:+.2f}              [{r[3]:+.2f}]               {r[4]:+.3f} / {r[5]:+.3f}  |  {r[6]:+.3f} / {r[7]:+.3f}             {r[8]:.3f}, {r[9]:.3f}")
print(f"\nwidest point ranges over z = {out[:, 1].min():+.2f} .. {out[:, 1].max():+.2f}; correlation with the order-interval value: {np.corrcoef(out[:, 1], out[:, 3])[0, 1]:.2f}")
print(f"largest one-sided slope magnitude at the widest point, step 0.02: {np.abs(out[:, 6:8]).max():.3f}  (a corner would show slopes of order alpha, i.e. {out[:, 8].min():.2f} .. {out[:, 8].max():.2f})")
# approach to the cone at the black pole
for z in (-0.9, -0.98, -0.998):
    ratio = np.array([R(t, z)/((1 + z)*alpha(t)) for t in ths[::6]]); print(f"z = {z:+.3f}: R / [alpha(theta)(1+z)] over hue: {ratio.min():.2f} .. {ratio.max():.2f}")
