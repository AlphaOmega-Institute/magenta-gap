"""Linear span of the spectral signals, the collinear long-wavelength stretch, and the Hilbert metric.

Inputs   Stockman-Sharpe table, 390-700 nm.  Parts 2 and 4 use raw channels normalised by their sum.
Checks   1. Three hue projectors from one arc, and the spectral signal vectors, have rank 3.
         2. Wavelengths with tabulated S exactly zero are collinear in chromaticity.
         3. On the round model the halved Hilbert metric equals the Klein hyperbolic distance, is invariant under
            E -> S E S^T, and saturations add like velocities.
         4. On the measured hull, Hilbert distances are unchanged under random linear changes of cone coordinates and chart;
            the distance from an interior point along hue theta is (1/2) log[(1 + r/alpha(theta+pi))/(1 - r/alpha(theta))].
Status   Numerical.
Prior art and scope (agreed in review)
         * The Hilbert metric on a general convex gamut, and its comparison with MacAdam ellipses, is in Berthier (2020),
           Sec. 7.3, Eq. (106).  Reading it as perceptual discrimination is an additional hypothesis.
         * Part 4 takes the point (1/3, 1/3) of the raw chart as its interior reference point.  The invariance does not
           depend on that choice; the formula holds from any interior point.
         * The collinearity argument limits linear or affine encodings that keep those lights distinct.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from common import load_table, REPO
import sys, numpy as np
from scipy.spatial import ConvexHull
I2 = np.eye(2); sx = np.array([[0., 1.], [1., 0.]]); sz = np.array([[1., 0.], [0., -1.]])
P = lambda th: (I2 + np.cos(th)*sz + np.sin(th)*sx)/2
wl, LMS = load_table()
m = (wl >= 390) & (wl <= 700); LMS, wl = LMS[m], wl[m]

print("== 1. the purple gap is not a linear defect (point raised in review) ==")
arc = np.array([P(t).ravel() for t in (0.2, 0.9, 1.7)])
print("rank of three hue projectors from one short arc (space of real symmetric 2x2 is 3-dim):", np.linalg.matrix_rank(arc))
print("rank of the spectral lights' cone-signal vectors:", np.linalg.matrix_rank(LMS), "-> linear span is everything, linear defect = 0")

print("\n== 2. straight stretch of the locus (coordinate-free obstruction to roundness under linear maps) ==")
zeroS = wl[LMS[:, 2] == 0]
print(f"wavelengths with tabulated S exactly 0: {len(zeroS)} of {len(wl)}, from {zeroS.min():.0f} nm; these are exactly collinear (the line S = 0)")
tiny = wl[LMS[:, 2] / LMS.sum(1) < 1e-3]
print(f"wavelengths with S share below 1e-3: {len(tiny)}, from {tiny.min():.0f} nm")

print("\n== 3. Hilbert metric on the round model ==")
def dH(A, B):                       # Hilbert projective metric on the positive cone, halved
    w = np.linalg.eigvals(np.linalg.solve(B, A)).real
    return 0.5*np.log(w.max()/w.min())
rng = np.random.default_rng(3); ok = True
for _ in range(500):
    r1, r2, t1, t2 = rng.uniform(0, .95), rng.uniform(0, .95), rng.uniform(0, 6.28), rng.uniform(0, 6.28)
    A = (I2 + r1*(np.cos(t1)*sz + np.sin(t1)*sx))/2; B = (I2 + r2*(np.cos(t2)*sz + np.sin(t2)*sx))/2
    Smat = rng.normal(size=(2, 2))
    ok &= np.isclose(dH(A, B), dH(Smat@A@Smat.T, Smat@B@Smat.T)) and np.isclose(dH(A, I2/2), np.arctanh(r1))
    # Klein-model hyperbolic distance between the two chromaticity points
    p, q = r1*np.array([np.cos(t1), np.sin(t1)]), r2*np.array([np.cos(t2), np.sin(t2)])
    ok &= np.isclose(dH(A, B), np.arccosh((1 - p@q)/np.sqrt((1 - p@p)*(1 - q@q))))
print("invariant under every E -> S E S^T; distance from white = artanh(saturation); equals the Klein hyperbolic distance:", ok)
eta = 0.7; Sb = np.diag([np.exp(eta/2), np.exp(-eta/2)]); r1 = 0.4
B = Sb @ ((I2 + r1*sz)/2) @ Sb.T; B /= np.trace(B)
print("boost along the hue axis adds saturations like velocities:", np.isclose(B[0, 0] - B[1, 1], (r1 + np.tanh(eta))/(1 + r1*np.tanh(eta))))

print("\n== 4. Hilbert metric on the measured gamut (convex hull of the locus) ==")
class Gamut:
    def __init__(self, pts):
        h = ConvexHull(pts); self.n, self.c = h.equations[:, :2], h.equations[:, 2]; self.vertices = h.vertices
    def reach(self, p, v):          # largest t with p + t v inside
        den = self.n @ v; num = -(self.n @ p + self.c)
        return np.min(num[den > 1e-15] / den[den > 1e-15])
    def dist(self, p, q):
        v = q - p; tp, tm = self.reach(p, v), self.reach(p, -v)
        return 0.5*np.log((tm + 1)/tm * tp/(tp - 1))
    def norm(self, p, v):
        return 0.5*(1/self.reach(p, v) + 1/self.reach(p, -v))
def chart(X, ell, u, w): return np.column_stack([X@u, X@w]) / (X@ell)[:, None]
ell0, u0, w0 = np.ones(3), np.array([1., 0, 0]), np.array([0, 1., 0])
G0 = Gamut(chart(LMS, ell0, u0, w0))
mix = rng.dirichlet(np.ones(len(wl))*0.02, size=60) @ LMS + 0.02*LMS.sum(0)/len(wl)   # interior cone points
worst = 0.0
for trial in range(20):
    M = rng.normal(size=(3, 3))
    ell = np.linalg.solve(M.T, ell0 + rng.uniform(-0.3, 0.3, 3))          # still positive on the cone
    assert np.all((LMS@M.T)@ell > 0)
    u, w = rng.normal(size=3), rng.normal(size=3)
    G = Gamut(chart(LMS@M.T, ell, u, w)); a0 = chart(mix, ell0, u0, w0); a1 = chart(mix@M.T, ell, u, w)
    for i in range(0, 60, 2):
        d0, d1 = G0.dist(a0[i], a0[i+1]), G.dist(a1[i], a1[i+1]); worst = max(worst, abs(d0 - d1)/d0)
print(f"distances unchanged under 20 random linear changes of cone coordinates and chromaticity charts: worst relative change {worst:.1e}")
white = np.zeros(2) + 1/3; ok = True
for th in rng.uniform(-np.pi, np.pi, 200):
    u = np.array([np.cos(th), np.sin(th)]); ap, am = G0.reach(white, u), G0.reach(white, -u); r = rng.uniform(0.05, 0.95)*ap
    ok &= np.isclose(G0.dist(white, white + r*u), 0.5*np.log((1 + r/am)/(1 - r/ap)))
print("distance from white along hue theta = (1/2) log[(1 + r/alpha(theta+pi)) / (1 - r/alpha(theta))]:", ok)
