"""Dynamics that reaches the defect separates readout extensions; illustrative altered-colour points.

Inputs   None.
Checks   Two phenomenal effects that agree on every paired state stay identical under dynamics preserving the range and
         differ under dynamics mixing range and defect.  Three illustrative points (chromatic at the black level, at the
         white level, and beyond the rim) are not effects in a bicone whose tips are those levels; the enlargement about
         mid-grey that admits each.
Status   Exact statements, checked numerically.
Scope    (agreed in review) The points are illustrative coordinates, not measured appearances.  The incompatibility is
         conditional on three assumptions held together: ordinary black and white are the zero and unit effects, altered
         appearances have those same height values, and all appearances lie in this bicone.  It does not show that the
         Berthier-Provenzi approach excludes altered colours.
"""
import numpy as np
from scipy.linalg import expm
rng = np.random.default_rng(31); n, m = 3, 5
M = rng.normal(size=(m, n)) + 1j*rng.normal(size=(m, n)); u, _ = np.linalg.qr(M); Pran = u @ u.conj().T; Pdef = np.eye(m) - Pran
A0 = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)); H0 = A0 @ A0.conj().T; F_psi = H0/np.linalg.eigvalsh(H0).max()
F1, F2 = u @ F_psi @ u.conj().T, u @ F_psi @ u.conj().T + 0.7*Pdef          # two extensions that agree on every paired state
unit = lambda k: (lambda v: v/np.linalg.norm(v))(rng.normal(size=k) + 1j*rng.normal(size=k))
print("== qualification 2: dynamics that reaches the defect separates the two extensions ==")
X = rng.normal(size=(m, m)) + 1j*rng.normal(size=(m, m)); Hmix = (X + X.conj().T)/2                 # does not preserve the range
Hkeep = Pran @ Hmix @ Pran + Pdef @ Hmix @ Pdef                                                     # preserves the range
for name, H in (("phenomenal dynamics preserving the range", Hkeep), ("phenomenal dynamics mixing range and defect", Hmix)):
    gap = max(abs(np.vdot(v, (F2 - F1) @ v).real) for v in (expm(-1j*t*H) @ (u @ unit(n)) for t in np.linspace(0.1, 3, 30) for _ in range(20)))
    print(f"  {name}: largest difference between the two readouts on evolved paired states = {gap:.3f}")
print("\n== can a bicone readout hold the altered colours named in the July manuscript? ==")
I2 = np.eye(2); sx = np.array([[0., 1.], [1., 0.]]); sz = np.array([[1., 0.], [0., -1.]])
E = lambda x, y, z: 0.5*((1 + z)*I2 + x*sz + y*sx)
for name, (x, y, z) in (("chromatic at the black level (stygian)", (0.5, 0.0, -1.0)), ("chromatic at the white level (self-luminous)", (0.5, 0.0, 1.0)),
                        ("more saturated than the rim at mid level", (1.3, 0.0, 0.0)), ("ordinary rim colour", (1.0, 0.0, 0.0))):
    ev = np.linalg.eigvalsh(E(x, y, z)); print(f"  {name}: eigenvalues {np.round(ev, 3)} -> effect: {bool(ev.min() >= -1e-12 and ev.max() <= 1 + 1e-12)}")
# enlarge the bicone about mid-grey by a factor s: E' = 1/2 + (E - 1/2)/s must be an effect
print("  smallest enlargement of the bicone about mid-grey that admits each point:")
for name, (x, y, z) in (("stygian", (0.5, 0.0, -1.0)), ("self-luminous", (0.5, 0.0, 1.0)), ("oversaturated", (1.3, 0.0, 0.0))):
    s = np.hypot(x, y) + abs(z); print(f"    {name}: factor {s:.2f}; the enlarged body's dark tip then sits at z = {-s:.2f}, below ordinary black at z = -1: {s > 1}")
