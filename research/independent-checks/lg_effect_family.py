"""The colour effect family built from LG's ground-plane operators.

Inputs   None.  R_0, Q_0 built from eq:ground-polar-operators.
Checks   E = ((1+z) + x R_0 + y Q_0)/2 is an effect exactly when sqrt(x^2 + y^2) <= 1 - |z| (200,000 points); conjugation by J
         reverses hue; the complement is the central symmetry; their product is the mirror symmetry; on the ground plane the
         family is the set of effects with E + J E J a multiple of the identity.
Status   Exact identities, checked numerically.  The effect body is Berthier (2020), Eq. (80).
Scope    (agreed in review) The characterisation by E + J E J holds only on the ground plane; on a larger carrier it admits
         arbitrary off-diagonal blocks (see two_encoding_routes.py).  The zero and unit effects are not LG's ground vectors.
"""
import numpy as np
rng = np.random.default_rng(13)
I2 = np.eye(2, dtype=complex); J = np.diag([1.0, -1.0]).astype(complex)
R0 = np.array([[0, 1], [1, 0]], complex); Q0 = -1j*(J @ R0)
print("R0, Q0 Hermitian, square to the identity, anticommute with each other and with J:",
      all(np.allclose(M, M.conj().T) and np.allclose(M @ M, I2) and np.allclose(M @ J, -J @ M) for M in (R0, Q0)) and np.allclose(R0 @ Q0, -Q0 @ R0))
E = lambda x, y, z, w=0.0: 0.5*((1 + z)*I2 + x*R0 + y*Q0 + w*J)
agree = True; margin = []
for _ in range(200000):
    x, y, z = rng.uniform(-1.3, 1.3, 3); ev = np.linalg.eigvalsh(E(x, y, z))
    is_effect = ev.min() >= 0 and ev.max() <= 1; in_bicone = np.hypot(x, y) <= 1 - abs(z)
    if is_effect != in_bicone: agree = False; margin.append(abs(np.hypot(x, y) - (1 - abs(z))))
print("0 <= E <= 1 exactly when sqrt(x^2 + y^2) <= 1 - |z| (200,000 random points):", agree)
x, y, z = 0.3, -0.2, 0.25
print("conjugating by J reverses hue at fixed z:  J E J = E(-x, -y, z):", np.allclose(J @ E(x, y, z) @ J, E(-x, -y, z)))
print("complement is the central symmetry:  1 - E = E(-x, -y, -z):", np.allclose(I2 - E(x, y, z), E(-x, -y, -z)))
print("their product is the bicone's mirror symmetry:  1 - J E J = E(x, y, -z):", np.allclose(I2 - J @ E(x, y, z) @ J, E(x, y, -z)))
# which effects on the ground plane belong to the colour family?
ok = True
for _ in range(20000):
    x, y, z, w = rng.uniform(-1, 1, 4); Ef = E(x, y, z, w); S = Ef + J @ Ef @ J
    ok &= np.allclose(S, (1 + z)*I2 + w*J) and (np.allclose(S, S[0, 0]*I2) == (abs(w) < 1e-12))
print("a general ground-plane effect has a fourth coordinate w (along J); E + J E J is a multiple of the identity exactly when w = 0:", ok)
v0, vO = np.array([1, 0], complex), np.array([0, 1], complex)
print("for the colour family, the vacuum and the luminous ground give the same probability (1+z)/2:",
      np.isclose(np.vdot(v0, E(x, y, z) @ v0).real, (1 + z)/2) and np.isclose(np.vdot(vO, E(x, y, z) @ vO).real, (1 + z)/2))
rho_mix = 0.5*(np.outer(v0, v0.conj()) + np.outer(vO, vO.conj())); xi = np.array([1, 1], complex)/np.sqrt(2); rho_coh = np.outer(xi, xi.conj())
print("equal incoherent mixture of the two grounds: <R0>, <Q0>, <J> =", [round(float(np.trace(rho_mix @ M).real), 3) for M in (R0, Q0, J)], "| coherent Xi_0(0):", [round(float(np.trace(rho_coh @ M).real), 3) for M in (R0, Q0, J)])
