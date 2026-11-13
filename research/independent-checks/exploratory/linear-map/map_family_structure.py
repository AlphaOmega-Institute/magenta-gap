"""Exploratory: what the condition-dependent family of linear maps

    A_c = ( I_lambda (x) diag(w_{c,i}^(-D/2)) (x) I_2 ) A,
    A |lambda> = |lambda> (x) sum_i sqrt(a_i(lambda)) |i> (x) |v_i>,

supplies as linear algebra, and where a ground correspondence runs into a
constraint.  w_c are the cone catches of the adapting white of condition c.

Checks
  1. each A_c is injective; A_c* A_c is wavelength-diagonal; the polar
     isometry depends on the condition
  2. dimension count at every wavelength: output space, range of one map,
     and the joint span of the ranges over all conditions
  3. the matching rule obtained by equating normalised colour outputs is the
     geometric cone-gain rule, which composes consistently; the linear-gain
     interpolation does not
  4. the selected observables commute after compression to the joint
     span, although the opponent observables do not preserve that span
  5. a fixed ray for every condition: a grey example inside the existing
     carrier, and a direct-sum extension with a no-absorption outcome

Usage:  python map_family_structure.py TABLE_CSV
Status: exploratory linear algebra on the binned table.  No generator,
ground or dynamics is supplied by the construction itself.
"""
import csv
import pathlib
import sys

import numpy as np

rows = [r for r in csv.reader(pathlib.Path(sys.argv[1]).read_text().splitlines()) if r]
tab = np.array([[float(v) for v in r] for r in rows[1:]])
wl, lms = tab[:, 0], tab[:, 1:]
wt = np.ones(len(wl)); wt[0] = wt[-1] = 0.5
a = lms / (lms * wt[:, None]).sum(axis=0)
N = len(wl)
ang = np.deg2rad([0.0, 120.0, 240.0])
nvec = np.column_stack([np.cos(ang), np.sin(ang)])
v = np.column_stack([np.cos(ang / 2), np.sin(ang / 2)])
vperp = np.column_stack([-np.sin(ang / 2), np.cos(ang / 2)])
D = 0.85
rng = np.random.default_rng(11)
whites = np.vstack([np.ones(3), rng.dirichlet([6, 6, 6], size=5) * 3])      # condition 0 is balanced white
np.set_printoptions(precision=4, suppress=True, linewidth=150)


def amp(c):
    """Amplitudes of A_c |lambda> on the three vectors |i>|v_i>: shape (N, 3)."""
    return np.sqrt(a * whites[c] ** (-D))


print(f"Conditions: {len(whites)} adapting whites, strength D = {D}.")
print("1. Injectivity and polar part")
for c in range(len(whites)):
    w2 = (amp(c) ** 2).sum(axis=1)                  # diagonal of A_c* A_c
    assert w2.min() > 0
u = [amp(c) / np.linalg.norm(amp(c), axis=1, keepdims=True) for c in range(len(whites))]
ov = min(float((u[0] * u[c]).sum(axis=1).min()) for c in range(1, len(whites)))
print(f"   A_c* A_c = diag(sum_i w_i^-D a_i(lambda)) is positive for every condition, so each A_c is injective.")
print(f"   The polar isometry sends |lambda> to its normalised amplitude vector; across conditions the smallest overlap")
print(f"   <u_0 lambda | u_c lambda> is {ov:.4f}, so the final support moves with the condition.")

print("2. Dimensions at every wavelength")
basis_reached = np.array([np.kron(np.eye(3)[i], v[i]) for i in range(3)])          # |i>|v_i>
basis_never = np.array([np.kron(np.eye(3)[i], vperp[i]) for i in range(3)])        # |i>|v_i perp>
rank = np.array([np.linalg.matrix_rank(np.array([amp(c)[k] for c in range(len(whites))])) for k in range(N)])
print("   output space C^3 (x) C^2: 6;  range of one A_c: 1.  Joint span of the ranges over these conditions:")
for r in sorted(set(rank), reverse=True):
    ws = wl[rank == r]
    print(f"      dimension {r} at {int(ws.min())}-{int(ws.max())} nm ({len(ws)} wavelengths): the 5-dimensional Ker A_c* splits into "
          f"{r - 1} inside the joint span and {6 - r} outside it")
print(f"   The smaller span is where the table has S = 0 ({int(wl[lms[:, 2] == 0].min())}-{int(wl[lms[:, 2] == 0].max())} nm): no gain activates a zero cone component.")
print(f"   The vectors |i>|v_i-perp> are orthogonal to every A_c (largest overlap {np.abs(basis_never @ basis_reached.T).max():.1e}).")
print("   'Inside the joint span' means a linear combination of preparations; it does not mean that one allowed")
print("   preparation reaches that direction (amplitudes of every preparation are nonnegative).")

print("3. Matching rule and its consistency")
s1 = rng.uniform(0, 1, N)                                           # an arbitrary spectrum
cat1 = a.T @ (s1 * wt)
for c1, c2 in ((0, 1), (2, 4)):
    # normalised colour weights under each condition; equate and solve for the matching catches
    p1 = cat1 * whites[c1] ** (-D); p1 /= p1.sum()
    cat2 = p1 * whites[c2] ** D
    geo = cat1 * (whites[c2] / whites[c1]) ** D
    print(f"   conditions {c1}->{c2}: catches of the match from equal colour states versus geometric gain rule, "
          f"max relative difference {np.abs(cat2 / cat2.sum() - geo / geo.sum()).max():.1e}")
g = lambda r, form: (r ** D if form == "geo" else 1 - D + D * r)
r12 = whites[1] / whites[0]; r23 = whites[2] / whites[1]; r13 = whites[2] / whites[0]
for form in ("geo", "lin"):
    back = g(r12, form) * g(1 / r12, form)
    comp = g(r12, form) * g(r23, form) / g(r13, form)
    print(f"   {'geometric' if form == 'geo' else 'linear   '} gains: there-and-back {back.round(4)}; composition 1->2->3 against 1->3 {comp.round(4)}")
print(f"   example from review, D = 1/2, ratios (2,1,1), linear gains there-and-back: "
      f"{((1 - 0.5 + 0.5 * np.array([2., 1, 1])) * (1 - 0.5 + 0.5 / np.array([2., 1, 1]))).round(4)}")

print("4. Selected observables on the joint span")
sz = np.array([[1.0, 0.0], [0.0, -1.0]]); sx = np.array([[0.0, 1.0], [1.0, 0.0]])
ops = [np.kron(np.diag(np.eye(3)[i]), np.eye(2)) for i in range(3)] + [np.kron(np.eye(3), sz), np.kron(np.eye(3), sx)]
comp = [basis_reached @ O @ basis_reached.T for O in ops]
off = max(np.abs(C - np.diag(np.diag(C))).max() for C in comp)
print(f"   compressions of the three cone projections and two opponent observables to the joint span are all diagonal")
print(f"   (largest off-diagonal entry {off:.1e}), so their statistics on every A_c output are classical.")
print(f"   On the full output space the opponent observables do not commute with each other: |[sz, sx]| = "
      f"{np.abs(ops[3] @ ops[4] - ops[4] @ ops[3]).max():.0f}.")
Pi = basis_reached.T @ basis_reached
leak = max(np.linalg.norm((np.eye(6) - Pi) @ O @ Pi, 2) for O in ops[3:])
print(f"   They also do not preserve the joint span (largest norm of the part mapped outside it: {leak:.3f}). So the classical")
print("   description covers these static statistics only; measurement-update operations, if available, could leave the span.")

print("5. A fixed ray for every condition")
def gain_op(c): return np.kron(np.diag(whites[c] ** (-D / 2)), np.eye(2))
par = lambda x, y: np.isclose(abs(x @ y), np.linalg.norm(x) * np.linalg.norm(y))
y_grey = basis_reached.sum(axis=0) / np.sqrt(3)                        # equal amplitudes on the three |i>|v_i>
ray = [gain_op(c) @ y_grey for c in range(len(whites))]
ovl = min(abs(ray[0] @ r) / (np.linalg.norm(ray[0]) * np.linalg.norm(r)) for r in ray[1:])
print(f"   A vector with equal amplitudes on the three |i>|v_i> is moved by the gains (ray overlap falls to {ovl:.4f}).")
# (a) a fixed grey ray inside the existing carrier (first reviewer's example): two wavelengths with S = 0,
#     both on the S register level, with orthogonal two-level vectors
k1, k2 = int(np.where(wl == 616)[0][0]), int(np.where(wl == 617)[0][0])
assert lms[k1, 2] == 0 and lms[k2, 2] == 0
e_S = np.eye(3)[2]
y = np.zeros((N, 6))
y[k1] = np.kron(e_S, [1.0, 0.0]) / np.sqrt(2); y[k2] = np.kron(e_S, [0.0, 1.0]) / np.sqrt(2)
fixed = all(par((y @ gain_op(c).T).ravel(), y.ravel()) for c in range(len(whites)))
M = y.reshape(N, 3, 2); rho = np.einsum("nik,nil->kl", M, M)
prep_overlap = max(abs(amp(c)[k] @ (basis_reached @ y[k])) for c in range(len(whites)) for k in (k1, k2))
pS = float((M[:, 2, :] ** 2).sum())
print(f"   (a) inside the existing carrier: (|616>|S>|e1> + |617>|S>|e2>)/sqrt 2 keeps its ray under every gain: {fixed};")
print(f"       reduced two-level state {rho.round(12).tolist()} (maximally mixed); overlap with every preparation {prep_overlap:.1e};")
print(f"       but its S-channel probability is {pS:.1f}, so it is not a no-absorption state.")
print(f"       The cone projections sum to the identity (deviation {np.abs(sum(ops[:3]) - np.eye(6)).max():.1e}), so every normalised")
print("       state of this carrier registers in some channel.")
# (b) direct-sum extension: one dark input ray and one output ray, untouched by adaptation
a0 = 0.7
worst = 0.0
for c in range(len(whites)):
    Ac = np.zeros((6 * N, N))
    for k in range(N):
        Ac[6 * k:6 * k + 6, k] = basis_reached.T @ amp(c)[k]
    Ah = np.zeros((6 * N + 1, N + 1)); Ah[0, 0] = a0; Ah[1:, 1:] = Ac
    d = np.zeros(N + 1); d[0] = 1; Om = np.zeros(6 * N + 1); Om[0] = 1
    worst = max(worst, np.abs(Ah @ d - a0 * Om).max(), np.abs(Ah.T @ Om - a0 * d).max())
    assert np.linalg.svd(Ah, compute_uv=False).min() > 0
    # cone readouts extended by zero vanish on the added output ray; the defect is unchanged
    for i in range(3):
        Fi = np.zeros(6 * N + 1); Fi[1:] = np.tile(np.kron(np.eye(3)[i], np.ones(2)), N)
        assert Om @ (Fi * Om) == 0
    assert (6 * N + 1) - np.linalg.matrix_rank(Ah) == 6 * N - np.linalg.matrix_rank(Ac)
print(f"   (b) direct sum with a dark input ray |d> and an output ray |Omega_d>, A-hat_c = a0 (+) A_c with a0 = {a0}:")
print(f"       A-hat_c |d> = a0 |Omega_d> and A-hat_c* |Omega_d> = a0 |d> for every condition (largest deviation {worst:.1e});")
print("       each A-hat_c is injective; all three cone readouts vanish on |Omega_d>; the dimension of Ker A-hat_c* equals")
print("       that of Ker A_c*, so the added pair creates no further defect.")
