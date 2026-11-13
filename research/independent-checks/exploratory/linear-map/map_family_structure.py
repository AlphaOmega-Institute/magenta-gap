"""Exploratory: what the condition-dependent family of linear maps

    A_c = ( I_lambda (x) diag(w_{c,i}^(-D/2)) (x) I_2 ) A,
    A |lambda> = |lambda> (x) sum_i sqrt(a_i(lambda)) |i> (x) |v_i>,

supplies as linear algebra, and where a ground correspondence runs into a
constraint.  w_c are the cone catches of the adapting white of condition c.

Checks
  1. each A_c is injective; A_c* A_c is wavelength-diagonal; the polar
     isometry depends on the condition
  2. dimension count per wavelength: output space, range of one map, joint
     span of the ranges over all conditions, and what no condition reaches
  3. the matching rule obtained by equating normalised colour outputs is the
     geometric cone-gain rule, which composes consistently; the linear-gain
     interpolation does not
  4. the selected observables commute on the joint span
  5. a dark vector: when its image keeps one ray for every condition

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

print("2. Dimensions per wavelength")
k = int(np.argmin(np.abs(wl - 550)))
stack = np.array([amp(c)[k] for c in range(len(whites))])
print(f"   output space C^3 (x) C^2: 6;  range of one A_c: 1;  joint span of the ranges over conditions at 550 nm: "
      f"{np.linalg.matrix_rank(stack)}")
basis_reached = np.array([np.kron(np.eye(3)[i], v[i]) for i in range(3)])          # |i>|v_i>
basis_never = np.array([np.kron(np.eye(3)[i], vperp[i]) for i in range(3)])        # |i>|v_i perp>
print(f"   the joint span is spanned by the three vectors |i>|v_i>; its complement |i>|v_i-perp> (dimension "
      f"{np.linalg.matrix_rank(basis_never)}) is orthogonal to every A_c (largest overlap "
      f"{np.abs(basis_never @ basis_reached.T).max():.1e}).")
print("   So Ker A_c* has dimension 5 per wavelength: 2 inside the joint span, which other conditions reach,")
print("   and 3 outside it, which no condition in the family reaches.")

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

print("5. A dark vector and a fixed ground ray")
y_grey = basis_reached.sum(axis=0) / np.sqrt(3)                        # equal amplitudes: central colour state
def gain_op(c): return np.kron(np.diag(whites[c] ** (-D / 2)), np.eye(2))
ray = [gain_op(c) @ y_grey for c in range(len(whites))]
ovl = min(abs(ray[0] @ r) / (np.linalg.norm(ray[0]) * np.linalg.norm(r)) for r in ray[1:])
bl = []
for r in ray:
    M = r.reshape(3, 2); rho = M.T @ M; rho /= np.trace(rho)
    bl.append(np.hypot(rho[0, 0] - rho[1, 1], 2 * rho[0, 1]))
print(f"   Image of dark with equal amplitudes on the three |i>|v_i>: its colour state is the centre for balanced white,")
print(f"   but the gains move it: ray overlap with the balanced case falls to {ovl:.4f}, colour radius up to {max(bl):.3f}.")
print("   A fixed ground ray for every condition needs the dark image to be a common eigenvector of all gain operators.")
single = np.kron(np.eye(3)[0], v[0])
par = lambda x, y: np.isclose(abs(x @ y), np.linalg.norm(x) * np.linalg.norm(y))
ok = all(par(gain_op(c) @ single, single) for c in range(len(whites))) and not all(par(gain_op(c) @ y_grey, y_grey) for c in range(len(whites)))
print(f"   With distinct gains those are the vectors on a single register level (checked: {ok}), whose colour state is a rim")
print("   state, or vectors on an extra register level that the gains leave alone.")
