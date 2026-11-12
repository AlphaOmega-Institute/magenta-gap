"""The widest points of the continuum bounded-reflectance body, traced by horizontal supporting normals.

Inputs   Stockman-Sharpe table; PCHIP interpolation on [390, hi] nm with exact antiderivatives; hi = 700 by default, or the
         first command-line argument (830 for the full table).  No wavelength bins.
Checks   For 1,440 normals n = (cos phi, sin phi, 0): interior sign changes of n.c, the support point, its hue and height.
         Two or more transitions: independence of two transition vectors and that their cross product is the horizontal
         normal (the supporting-normal criterion).  One transition: the horizontal normal and a second, non-parallel normal
         both support the same edge colour over the whole spectrum (a crease).  Shares are weighted by hue angle of the
         boundary point, not by count of normals.
Status   Numerical, on sampled normals; not interval-certified.  The shares (about 78% / 22% on 390-700 nm) are estimates.
         The rigorous statement is the single certified point in research/object-colour.
Scope    Signal model and height coordinate as in solid_lib.py.  Attainable signals only.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
wl_all, LMS_all = load_table()
lo, hi = 390.0, (float(sys.argv[1]) if len(sys.argv) > 1 else 700.0)
P = PchipInterpolator(wl_all, LMS_all); AD = P.antiderivative()
tot = AD(hi) - AD(lo)
T = np.array([[2/3, -1/3, -1/3], [-1/3, 2/3, -1/3], [1, 1, 1]])/1.5
c = lambda lam: (P(lam)/tot) @ T.T                              # density in (x, y, t) coordinates
I = lambda a, b: ((AD(b) - AD(a))/tot) @ T.T                    # exact integral of the density over [a, b]
grid = np.arange(lo, hi + 1e-9, 0.005); Cg = c(grid)
print("white =", np.round(I(lo, hi), 12), "(should be (0, 0, 2))")

def zeros(n):
    f = Cg @ n; s = np.sign(f); idx = np.flatnonzero(s[:-1]*s[1:] < 0)
    zs = [brentq(lambda l: float(c(l) @ n), grid[i], grid[i+1], xtol=1e-12) for i in idx]
    touching = int(np.sum(np.abs(f[1:-1]) < 1e-14))
    return zs, f, touching
def support_point(n, zs, f):
    edges = [lo] + zs + [hi]; x = np.zeros(3); bands = []
    for a, b in zip(edges[:-1], edges[1:]):
        if float(c((a + b)/2) @ n) > 0: x += I(a, b); bands.append((a, b))
    return x, bands
def same_pattern(n, bands):
    """n supports the same reflectance: n.c >= 0 on its bands and <= 0 off them, over the whole spectrum."""
    f = Cg @ n; inside = np.zeros(len(grid), bool)
    for a, b in bands: inside |= (grid >= a) & (grid <= b)
    tol = 1e-12*np.abs(f).max()
    return bool(np.all(f[inside] >= -tol) and np.all(f[~inside] <= tol))

phis = np.linspace(0, 2*np.pi, 1440, endpoint=False); rec = []
for phi in phis:
    n = np.array([np.cos(phi), np.sin(phi), 0.0]); zs, f, touching = zeros(n); x, bands = support_point(n, zs, f)
    theta, z = np.arctan2(x[1], x[0]), x[2] - 1
    simple = all(abs(float(c(l, ) @ n) - 0) < 1e-9 and abs(float(P.derivative()(l)/tot @ T.T @ n)) > 1e-12 for l in zs)
    info = {"phi": phi, "theta": theta, "z": z, "nz": len(zs), "zs": zs, "simple": simple, "touch": touching}
    if len(zs) == 2:
        c1, c2 = c(zs[0]), c(zs[1]); cr = np.cross(c1, c2)
        info["sin"] = np.linalg.norm(cr)/(np.linalg.norm(c1)*np.linalg.norm(c2))
        nn = cr/np.linalg.norm(cr); info["parallel"] = abs(abs(nn @ n) - 1); info["nz_comp"] = abs(nn[2])
    if len(zs) == 1:
        c1 = c(zs[0]); ok_self = same_pattern(n, bands)
        m = np.cross(c1, n); m /= np.linalg.norm(m)              # perpendicular to c(lambda1) and to n
        tilt = [e for e in (0.05, -0.05) if same_pattern(n + e*m, bands)]
        ext = []
        for end in (lo, hi):                                      # extreme candidates: planes through c(lambda1) and a spectrum end
            ne = np.cross(c1, c(end)); ne /= np.linalg.norm(ne)
            for sgn in (1, -1):
                if same_pattern(sgn*ne, bands): ext.append(sgn*ne)
        info.update(ok_self=ok_self, tilt=len(tilt), ext=ext)
        if len(ext) == 2: info["wedge"] = np.degrees(np.arccos(np.clip(ext[0] @ ext[1], -1, 1)))
    rec.append(info)
nz = np.array([r["nz"] for r in rec]); th = np.unwrap([r["theta"] for r in rec])
print(f"interior sign changes of n.c along the rim: " + ", ".join(f"{k}: {np.sum(nz == k)} of {len(rec)} normals" for k in sorted(set(nz))))
print("all zeros simple (n.c' != 0):", all(r["simple"] for r in rec), "| grid points where n.c touches zero without crossing:", sum(r["touch"] for r in rec))
print(f"hue advances monotonically along the rim: {bool(np.all(np.diff(th) > 0) or np.all(np.diff(th) < 0))}; total hue swept {abs(th[-1] - th[0]) + abs(th[1] - th[0]):.4f} rad")
dth = np.abs(np.diff(np.append(th, th[0] + np.sign(th[-1] - th[0])*2*np.pi)))
for k in sorted(set(nz)): print(f"  share of the hue circle whose widest point has {k} interior transition(s): {dth[nz == k].sum()/(2*np.pi):.1%}")
two = [r for r in rec if r["nz"] == 2]
print(f"\ntwo-transition rim points: {len(two)}")
print(f"  independence of c(lambda1), c(lambda2): sine of their angle ranges {min(r['sin'] for r in two):.3f} .. {max(r['sin'] for r in two):.3f}")
print(f"  the unique normal c1 x c2 is parallel to the horizontal n to within {max(r['parallel'] for r in two):.1e}, height component at most {max(r['nz_comp'] for r in two):.1e}")
print(f"  transitions stay inside the spectrum: nearest approach to an end is {min(min(r['zs'][0] - lo, hi - r['zs'][1]) for r in two):.2f} nm")
print(f"  height of these widest points: z from {min(r['z'] for r in two):+.3f} to {max(r['z'] for r in two):+.3f}")
one = [r for r in rec if r["nz"] == 1]
print(f"\none-transition rim points: {len(one)}")
print(f"  the horizontal normal supports the edge colour over the whole spectrum: {sum(r['ok_self'] for r in one)} of {len(one)}")
print(f"  a second, non-parallel normal (tilted 0.05 about c(lambda1)) also supports it: {sum(r['tilt'] >= 1 for r in one)} of {len(one)} (both tilt signs: {sum(r['tilt'] == 2 for r in one)})")
print(f"  both extreme normals (planes through c(lambda1) and each spectrum end) support it: {sum(len(r['ext']) == 2 for r in one)} of {len(one)}" + (f"; wedge between them {min(r['wedge'] for r in one if 'wedge' in r):.1f} .. {max(r['wedge'] for r in one if 'wedge' in r):.1f} deg" if any('wedge' in r for r in one) else ""))
print(f"  transition wavelengths {min(r['zs'][0] for r in one):.1f} .. {max(r['zs'][0] for r in one):.1f} nm; heights z from {min(r['z'] for r in one):+.3f} to {max(r['z'] for r in one):+.3f}")
other = [r for r in rec if r["nz"] not in (1, 2)]
for r in other[:6]: print("  other:", round(np.degrees(r["theta"]), 1), r["nz"], [round(v, 2) for v in r["zs"]])
