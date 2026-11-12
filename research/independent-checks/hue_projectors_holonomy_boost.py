"""Hue projectors, their sign holonomy, and the action of boosts on pure hues.

Inputs   None.
Checks   E = w+ P(theta) + w- P(theta + pi); the grading J(theta) = P(theta) - P(theta + pi); transport of the hue eigenvector
         once round the grey axis returns its negative; E -> S E S^T with det S = 1 maps pure hues to pure hues with
         tan(theta'/2) = exp(-eta) tan(theta/2).
Status   Exact identities, checked numerically.  The boost map is Berthier (2020), Eq. (89).
Scope    The sign holonomy is not observable in E, because the projector itself returns unchanged.  The projected
         (Grassmann) connection on this real line is flat; the hyperbolic metric of the colour literature is a different
         curvature question.
"""
import numpy as np
sx, sz, I2 = np.array([[0., 1.], [1., 0.]]), np.array([[1., 0.], [0., -1.]]), np.eye(2)
P = lambda th: (I2 + np.cos(th)*sz + np.sin(th)*sx)/2          # hue projector
J = lambda th: np.cos(th)*sz + np.sin(th)*sx                   # hue-opponent grading
rng = np.random.default_rng(0); ok = True
for _ in range(2000):
    th, z = rng.uniform(0, 2*np.pi), rng.uniform(-1, 1); r = rng.uniform(0, 1 - abs(z))
    E = ((1 + z)*I2 + r*J(th))/2
    lp, lm = (1 + z + r)/2, (1 + z - r)/2
    ok &= np.allclose(J(th)@J(th), I2) and np.allclose(J(th), P(th) - P(th + np.pi))
    ok &= np.allclose(E, lp*P(th) + lm*P(th + np.pi)) and np.isclose(4*np.linalg.det(E), (1 + z)**2 - r**2) and np.isclose(4*lp*lm, (1 + z)**2 - r**2)
print("J(theta)^2 = 1, J = P(hue) - P(opposite hue), E = weight * hue + weight * opposite hue, det = product of weights:", ok)
# holonomy of the hue eigenline around the achromatic axis (Kato parallel transport, dv = [P', P] v)
n = 20000; v = np.array([1.0, 0.0]); ths = np.linspace(0, 2*np.pi, n + 1)
for a, b in zip(ths[:-1], ths[1:]):
    v = P(b) @ v; v /= np.linalg.norm(v)
print("eigenvector after one trip around the hue circle:", np.round(v, 6), " (started at [1, 0]) -> holonomy -1")
# Lorentz transformations act on hue by Moebius maps:  E -> S E S^T with det S = 1
eta = 0.8; Sb = np.diag([np.exp(eta/2), np.exp(-eta/2)])
ok2 = True
for th in np.linspace(0.1, 6.1, 13):
    E2 = Sb @ P(th) @ Sb.T
    ok2 &= np.isclose(np.linalg.det(E2), 0)                      # null stays null
    th2 = np.arctan2(2*E2[0, 1], E2[0, 0] - E2[1, 1])
    ok2 &= np.isclose(np.tan(th2/2), np.exp(-eta)*np.tan(th/2))
print("boost maps pure hues to pure hues with tan(theta'/2) = exp(-eta) tan(theta/2):", ok2)
print("example: hues at 90 and 270 degrees move to", [round(float(np.degrees(2*np.arctan(np.exp(-eta)*np.tan(np.radians(t)/2))) % 360), 1) for t in (90, 270)], "(crowd toward the boost axis)")
