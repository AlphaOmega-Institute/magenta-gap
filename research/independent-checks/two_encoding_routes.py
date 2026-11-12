"""Two ways colour could enter LG's carrier, on a finite model with a defect.

Inputs   None.  H_Psi = C^3, H_Phi = C^5, a random isometry u; structural operators as in eq:ground-polar-operators.
Checks   Route 1 (hue as intersector phase): <W_theta psi, E W_theta psi> = (1 + z + x cos theta + y sin theta)/2 for every
         psi; on the bridge-selected graph every input gives <R>, <Q>, <J> = 1, 0, 0; an effect with E + J E J = 1 outside the
         three-parameter family.  Route 2 (colour as carried content): F_Phi = u F_Psi u* gives equal probabilities on
         corresponding states, varying with the input; adding a multiple of the defect projection changes none of them.
Status   Exact statements, checked numerically on one random model.
Scope    Neither route is asserted as a colour realization.  Static paired statistics leave the defect readout
         undetermined; that is not permanent invisibility (see defect_extensions_and_altered_points.py).
"""
import numpy as np
rng = np.random.default_rng(31); n, m = 3, 5
M = rng.normal(size=(m, n)) + 1j*rng.normal(size=(m, n)); u, _ = np.linalg.qr(M)            # u* u = 1 on H_Psi
Z = lambda a, b: np.zeros((a, b), complex)
J = np.block([[np.eye(n), Z(n, m)], [Z(m, n), -np.eye(m)]]).astype(complex)
R = np.block([[Z(n, n), u.conj().T], [u, Z(m, m)]]); Q = -1j*(J @ R)
Pi = R @ R                                                                                   # projection onto the active carrier
W = lambda psi, th: np.concatenate([psi, np.exp(1j*th)*(u @ psi)])/np.sqrt(2)
E = lambda x, y, z: 0.5*((1 + z)*Pi + x*R + y*Q)
unit = lambda k: (lambda v: v/np.linalg.norm(v))(rng.normal(size=k) + 1j*rng.normal(size=k))

print("== route 1: hue as intersector phase ==")
ok = True
for _ in range(5000):
    psi, th = unit(n), rng.uniform(0, 2*np.pi); x, y, z = rng.uniform(-1, 1, 3); v = W(psi, th)
    ok &= np.isclose(np.vdot(v, E(x, y, z) @ v).real, (1 + z + x*np.cos(th) + y*np.sin(th))/2)
print("<W_theta psi, E W_theta psi> = (1 + z + x cos(theta) + y sin(theta))/2 for every psi:", ok)
agree = True
for _ in range(20000):
    x, y, z = rng.uniform(-1.3, 1.3, 3); ev = np.linalg.eigvalsh(E(x, y, z)); ev_act = ev[np.abs(ev) > 1e-12] if abs(1 + z) > 1e-9 else ev
    agree &= ((ev.min() >= -1e-12 and ev.max() <= 1 + 1e-12) == (np.hypot(x, y) <= 1 - abs(z)))
print("on the active carrier E is an effect exactly inside the bicone:", agree)
psi = unit(n); v = W(psi, 0.0); stat = [np.vdot(v, Mx @ v).real for Mx in (R, Q, J)]
print("on the bridge-selected graph (theta = 0) every input gives <R>, <Q>, <J> =", [round(float(s), 12) for s in stat])
B = rng.normal(size=(m, n)) + 1j*rng.normal(size=(m, n)); B *= 0.3/np.linalg.norm(B, 2)
Ebig = 0.5*np.eye(n + m) + 0.5*np.block([[Z(n, n), B.conj().T], [B, Z(m, m)]]); ev = np.linalg.eigvalsh(Ebig)
print("an effect with E + J E J = identity that is not in the three-parameter family (arbitrary off-diagonal block):",
      np.allclose(Ebig + J @ Ebig @ J, np.eye(n + m)) and ev.min() >= 0 and ev.max() <= 1)

print("\n== route 2: colour as content carried across the fixed bridge ==")
A0 = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)); H = A0 @ A0.conj().T; F_psi = H/np.linalg.eigvalsh(H).max()   # some physical effect
F_phi = u @ F_psi @ u.conj().T; Pdef = np.eye(m) - u @ u.conj().T
ok = True; spread = []
for _ in range(5000):
    psi = unit(n); p = np.vdot(psi, F_psi @ psi).real; spread.append(p)
    ok &= np.isclose(p, np.vdot(u @ psi, F_phi @ (u @ psi)).real)
print(f"paired effect F_Phi = u F_Psi u* gives the same probability on corresponding states: {ok}; the probability varies with the input ({min(spread):.2f} .. {max(spread):.2f})")
ok2 = all(np.isclose(np.vdot(u @ p_, (F_phi + c*Pdef) @ (u @ p_)).real, np.vdot(p_, F_psi @ p_).real) for p_ in [unit(n) for _ in range(2000)] for c in (0.0, 0.4, 1.0))
print("adding any multiple of the defect projection to F_Phi changes no paired probability:", ok2, "| F_Phi itself assigns probability", round(float(np.trace(F_phi @ Pdef).real), 12), "to every defect state")
