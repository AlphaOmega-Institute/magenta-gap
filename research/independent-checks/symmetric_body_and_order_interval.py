"""The symmetric body, its two symmetries, and the order-interval body of the measured gamut.

Inputs   Stockman-Sharpe table; white-balanced hull of the 390-700 nm chromaticities (September manuscript, Sec. 5).
Checks   A. eq:baseline is -4 det of the displacement of E = ((1+z) + x sz + y sx)/2; the generator of eq:generator has rank one.
         B. 1 - E is the central symmetry (hue + pi, z -> -z).  For the cone K over the measured hull with white on the axis,
            K intersected with (white - K) is r <= min(alpha(theta)(T+z), alpha(theta+pi)(T-z)); ridge heights; volumes
            compared with the manuscript's body r <= alpha(theta)(T-|z|).
         C. The Hilbert-metric angular scale grows without bound towards the hull boundary at almost every hue.
Status   A, B (identities): exact, checked numerically.  Ridge heights and volumes: numerical, for this hull and this height
         coordinate (the cone sum, not luminance).  C: numerical.
Scope    The order-interval body is an outer bound for images of stimuli bounded by a white under a linear response.  It is
         not the object-colour solid, and nothing here says which body describes perception.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, numpy as np
from scipy.spatial import ConvexHull
rng = np.random.default_rng(5)
I2 = np.eye(2); sx = np.array([[0., 1.], [1., 0.]]); sz = np.array([[1., 0.], [0., -1.]])
mat = lambda x, y, z: ((1 + z)*I2 + x*sz + y*sx)/2          # kappa = T = 1

print("== A. eq:baseline is minus four times the determinant of the displacement ==")
ok = True
for _ in range(2000):
    dx, dy, dz = rng.normal(size=3)
    ok &= np.isclose(dx*dx + dy*dy - dz*dz, -4*np.linalg.det(mat(dx, dy, dz) - mat(0, 0, 0)))
print("g0(dE,dE) = -4 det(dE):", ok)
th = 1.1; gen = mat(np.cos(th), np.sin(th), 1.0) - mat(0, 0, 0)   # tangent of eq:generator: (cos, sin, 1)
print("generator tangent of eq:generator has rank", np.linalg.matrix_rank(gen), "and g0 =", round(float(-4*np.linalg.det(gen)), 12))

print("\n== B. the two symmetries ==")
x, y, z = 0.21, -0.13, 0.4
print("E -> 1 - E sends (x, y, z) to (-x, -y, -z):", np.allclose(I2 - mat(x, y, z), mat(-x, -y, -z)), " (hue theta -> theta + pi, z -> -z)")

wl_all, LMS_all = load_table()
Eint = np.trapezoid(LMS_all, wl_all, axis=0); m = (wl_all >= 390) & (wl_all <= 700); wl = wl_all[m]
X = LMS_all[m] / Eint; c = X / X.sum(1, keepdims=True); pts = c[:, :2] - 1/3
hull = ConvexHull(pts); N, C = hull.equations[:, :2], hull.equations[:, 2]
def alpha(theta):
    u = np.array([np.cos(theta), np.sin(theta)]); den = N @ u
    return np.min(-C[den > 1e-15] / den[den > 1e-15])
def in_hull(p): return np.all(N @ p + C <= 1e-12)
T = 1.0
# order interval of the cone over the hull, with white w = (0, 0, height 2T):  X in K and w - X in K
okc = True
for _ in range(20000):
    thq, zq = rng.uniform(-np.pi, np.pi), rng.uniform(-T, T); rq = rng.uniform(0, 0.7)
    h = zq + T; p = rq*np.array([np.cos(thq), np.sin(thq)])
    member = in_hull(p / h) and in_hull(-p / (2*T - h))
    formula = rq <= min(alpha(thq)*(T + zq), alpha(thq + np.pi)*(T - zq))
    okc &= (member == formula) or abs(rq - min(alpha(thq)*(T + zq), alpha(thq + np.pi)*(T - zq))) < 1e-9
print("order-interval body of the measured gamut is r <= min( alpha(theta)(T+z), alpha(theta+pi)(T-z) ):", okc)

ths = np.linspace(-np.pi, np.pi, 3600, endpoint=False); a = np.array([alpha(t) for t in ths]); ao = np.array([alpha(t + np.pi) for t in ths])
zr = T*(ao - a)/(ao + a)
def dominant(theta):
    u = np.array([np.cos(theta), np.sin(theta)]); p = alpha(theta)*u
    i = np.argmin(np.abs(N @ p + C)); verts = wl[hull.simplices[i]]
    return f"{verts.min():.0f}-{verts.max():.0f} nm edge"
i_hi, i_lo = np.argmax(zr), np.argmin(zr)
print(f"ridge height z/T over hue: from {zr.min():+.3f} (hue {np.degrees(ths[i_lo]):.0f} deg, boundary on {dominant(ths[i_lo])}) to {zr.max():+.3f} (hue {np.degrees(ths[i_hi]):.0f} deg, boundary on {dominant(ths[i_hi])})")
print(f"share of hue angles with |z_ridge|/T > 0.2: {np.mean(np.abs(zr) > 0.2):.0%}; alpha(theta)/alpha(theta+pi) ranges {np.min(a/ao):.2f} .. {np.max(a/ao):.2f}")
# volumes (per unit T^3): manuscript body r <= alpha(theta)(T-|z|); order-interval body as above
zs = np.linspace(-T, T, 2001); dth = ths[1] - ths[0]; dz = zs[1] - zs[0]
Rp = a[:, None]*(T - np.abs(zs))[None]; Ro = np.minimum(a[:, None]*(T + zs)[None], ao[:, None]*(T - zs)[None])
Vp, Vo, Vi = (0.5*Rp**2).sum()*dth*dz, (0.5*Ro**2).sum()*dth*dz, (0.5*np.minimum(Rp, Ro)**2).sum()*dth*dz
print(f"volumes: manuscript body {Vp:.4f}, order-interval body {Vo:.4f}, common part {Vi:.4f} -> {1 - Vi/Vp:.0%} of the manuscript body lies outside the order-interval body")

print("\n== C. Hilbert metric as a candidate for the angular scale B ==")
def reach(p, v):
    den = N @ v; num = -(N @ p + C); return np.min(num[den > 1e-15]/den[den > 1e-15])
for f in (0.9, 0.99, 0.999):
    b = []
    for t in ths[::10]:
        u = np.array([np.cos(t), np.sin(t)]); e = np.array([-np.sin(t), np.cos(t)]); p = f*alpha(t)*u
        b.append(0.5*(1/reach(p, e) + 1/reach(p, -e)))
    b = np.array(b); print(f"at {f:.3f} of the ceiling: Hilbert angular scale beta_H over hue: min {b.min():.1f}, median {np.median(b):.1f}")
