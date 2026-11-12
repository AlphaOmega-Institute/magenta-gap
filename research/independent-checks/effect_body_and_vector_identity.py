"""The symmetric bicone as an effect body, and the quadratic identity for vectors.

Inputs   None.
Checks   r <= 1 - |z| exactly when 0 <= E <= 1 for E = ((1+z) + x sx + y sz)/2; the boundary sheets are det E = 0 and
         det(1 - E) = 0; X^2 + Y^2 + N^2 = Z^2 for the four bilinears of a two-component vector.
Status   Exact identities, checked numerically.  Prior art: M. Berthier (2020), J. Math. Neurosci. 10:14, Sec. 5.1, Eq. (80).
Scope    States and effects are different objects.  The vector identity gives a cone; the bounded bicone comes from effects.
"""
import numpy as np
rng = np.random.default_rng(0)
sx, sy, sz = np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])
# 1. The symmetric bicone r <= 1 - |z| is the set of real symmetric 2x2 matrices with 0 <= E <= 1,
#    E = ((1+z) 1 + x sx + y sz)/2, and the Lorentzian form is 4 det E = (1+z)^2 - r^2.
ok = True
for _ in range(20000):
    x, y, z = rng.uniform(-1.5, 1.5, 3); r = np.hypot(x, y)
    E = ((1 + z)*np.eye(2) + x*sx + y*sz)/2
    w = np.linalg.eigvalsh(E)
    ok &= (w.min() >= 0 and w.max() <= 1) == (r <= 1 - abs(z))
    ok &= np.isclose(4*np.linalg.det(E), (1 + z)**2 - r**2) and np.isclose(4*np.linalg.det(np.eye(2) - E), (1 - z)**2 - r**2)
print("bicone == effect space {0 <= E <= 1}; boundary sheets are det E = 0 and det(1-E) = 0:", ok)
# ridge z = 0 on the boundary: rank-one projectors
th = 0.7; Pth = (np.eye(2) + np.cos(th)*sx + np.sin(th)*sz)/2
print("boundary point at z=0, r=1 is a projector:", np.allclose(Pth@Pth, Pth), "trace", np.trace(Pth).real)
# 2. LG ground plane: xi = (a, b) in C^{1,1}, J = sz, R = sx, Q = sy.  Bilinears X=<R>, Y=<Q>, N=<J>, Z=<1>.
ok2 = True
for _ in range(2000):
    xi = rng.normal(size=2) + 1j*rng.normal(size=2)
    X, Y, N, Z = (np.vdot(xi, M@xi).real for M in (sx, sy, sz, np.eye(2)))
    ok2 &= np.isclose(X**2 + Y**2 + N**2, Z**2)
print("X^2 + Y^2 + N^2 = Z^2 for every vector:", ok2)
for phi in (0.0, 1.0, 2.5):
    xi = np.array([1, np.exp(1j*phi)])/np.sqrt(2)
    X, Y, N, Z = (np.vdot(xi, M@xi).real for M in (sx, sy, sz, np.eye(2)))
    print(f"  Krein-null state with relative phase {phi}: N = {N:+.1e}, (X, Y)/Z = ({X/Z:+.4f}, {Y/Z:+.4f}) = (cos, sin) of phase: {np.allclose([X/Z, Y/Z], [np.cos(phi), np.sin(phi)])}; X^2+Y^2-Z^2 = {X*X+Y*Y-Z*Z:+.1e}")
