"""On LG's ground plane the squared Krein form of a vector is a (2,1) interval.

Inputs   None.  Basis |0>, |Omega_0>; J = diag(1, -1) (LG, eq:J-ground-plane); X, Y the expectations of R_0, Q_0.
Checks   [v, v]^2 = Z^2 - X^2 - Y^2 for 20,000 vectors; maps preserving the Krein form preserve that interval; the neutral
         circle Xi_0(theta) (eq:neutral-ground-circle) lies on the null cone at angle theta.
Status   Exact identity, checked numerically.
Scope    Pure vectors only.  For a mixed state the interval gains a determinant term, so equal sector weights do not imply
         saturation.  This is not the restriction of the Krein form to the fixed neutral graph, which vanishes.
"""
import numpy as np
from scipy.linalg import expm
rng = np.random.default_rng(8)
J = np.diag([1.0, -1.0]); sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]])
def stokes(v): return np.array([np.vdot(v, sx @ v).real, np.vdot(v, sy @ v).real, np.vdot(v, v).real]), np.vdot(v, J @ v).real
ok1 = ok2 = ok3 = True; eta = np.diag([1.0, 1.0, -1.0])
for _ in range(20000):
    v = rng.normal(size=2) + 1j*rng.normal(size=2); (X, Y, Z), N = stokes(v)
    ok1 &= np.isclose(Z*Z - X*X - Y*Y, N*N)                                   # Krein form squared = Lorentz interval
    # a random element of the group preserving the Krein form: exp of a J-anti-selfadjoint generator
    a, b, c = rng.normal(size=3)*0.7; G = 1j*a*np.diag([1.0, -1.0]) + b*sx + c*sy; g = expm(G)
    ok2 &= np.allclose(g.conj().T @ J @ g, J)
    (X2, Y2, Z2), N2 = stokes(g @ v)
    ok3 &= np.isclose(N2, N) and np.isclose(Z2*Z2 - X2*X2 - Y2*Y2, Z*Z - X*X - Y*Y)
print("for every vector of the ground plane, [v,v]^2 = Z^2 - X^2 - Y^2:", ok1)
print("maps preserving the Krein form preserve that (2,1) interval of (X, Y, Z):", ok2 and ok3)
for th in (0.0, 1.0, np.pi):
    xi = np.array([1.0, np.exp(1j*th)])/np.sqrt(2); (X, Y, Z), N = stokes(xi)
    print(f"neutral circle Xi_0({th:.2f}): [v,v] = {N:+.1e}, (X, Y, Z) = ({X:+.3f}, {Y:+.3f}, {Z:.3f}), angle = {np.arctan2(Y, X):.2f}")
print("vacuum |0> alone: ", [round(float(t), 3) for t in stokes(np.array([1.0, 0]))[0]], "| ground |Omega_0> alone:", [round(float(t), 3) for t in stokes(np.array([0, 1.0]))[0]], "(both on the axis X = Y = 0)")
