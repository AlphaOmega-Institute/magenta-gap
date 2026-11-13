"""Exploratory: one injective linear map that carries (i) three positive cone
channels and (ii) a two-level opponent factor, built from the cone table.

Construction (finite wavelength bins, white-balanced catches a_i(lambda) >= 0):

    A |lambda> = |lambda> (x) sum_i sqrt(a_i(lambda)) |i> (x) |v_i>,   i in {L, M, S}

where |v_L>, |v_M>, |v_S> are three real unit vectors of a two-level system
at 60 degree spacing (Bloch angles 0, 120, 240 degrees).

Checks
  A. the map is injective and both families of observables pull back to
     wavelength-diagonal operators
  B. the cone channels are positive operators on the three-level register and
     reproduce the cone catches; no positive S-channel on a two-level factor
     alone can do so (two S = 0 lights with different L:M)
  C. the reduced two-level state of each light has Bloch vector equal to the
     standard opponent combination of balanced catches, so the receptor
     triangle is inscribed in the state disc with white at its centre
  D. where the monochromatic lights and their mixtures sit in that disc

Usage:  python cone_register_and_trine.py TABLE_CSV
Status: exploratory, numerical.  The construction is one choice; it is fixed
by the white balance and is not fitted to appearance data.
"""
import csv
import pathlib
import sys

import numpy as np
from scipy.spatial import ConvexHull

rows = [r for r in csv.reader(pathlib.Path(sys.argv[1]).read_text().splitlines()) if r]
tab = np.array([[float(v) for v in r] for r in rows[1:]])
wl, lms = tab[:, 0], tab[:, 1:]
wt = np.ones(len(wl)); wt[0] = wt[-1] = 0.5
a = lms / (lms * wt[:, None]).sum(axis=0)          # white-balanced: equal-energy white has equal catches
N = len(wl)
np.set_printoptions(precision=4, suppress=True, linewidth=150)

ang = np.deg2rad([0.0, 120.0, 240.0])              # Bloch angles of the three two-level states
v = np.column_stack([np.cos(ang / 2), np.sin(ang / 2)])     # real unit vectors, rows v_L, v_M, v_S
sz = np.array([[1.0, 0.0], [0.0, -1.0]]); sx = np.array([[0.0, 1.0], [1.0, 0.0]])

# A as an isometry-like matrix: output index (lambda, i, k) with i the register and k the two-level factor
A = np.zeros((N, 3, 2, N))
for i in range(3):
    for k in range(2):
        A[np.arange(N), i, k, np.arange(N)] = np.sqrt(a[:, i]) * v[i, k]
Am = A.reshape(N * 6, N)

print("A. Injectivity and diagonality")
sv = np.linalg.svd(Am, compute_uv=False)
print(f"   smallest singular value of A: {sv.min():.3e} (largest {sv.max():.3e}); "
      f"A*A is diagonal with entries a_L+a_M+a_S, smallest {a.sum(axis=1).min():.3e} at {int(wl[np.argmin(a.sum(axis=1))])} nm")


def pull(op_reg, op_two):
    F = np.kron(np.eye(N), np.kron(op_reg, op_two))
    return Am.T @ F @ Am


offdiag = 0.0
for i in range(3):
    P = np.zeros((3, 3)); P[i, i] = 1.0
    G = pull(P, np.eye(2))
    offdiag = max(offdiag, np.abs(G - np.diag(np.diag(G))).max())
    assert np.allclose(np.diag(G), a[:, i])
Gx, Gz = pull(np.eye(3), sx), pull(np.eye(3), sz)
offdiag = max(offdiag, np.abs(Gx - np.diag(np.diag(Gx))).max(), np.abs(Gz - np.diag(np.diag(Gz))).max())
print(f"   pulled-back cone projections equal diag(a_i); largest off-diagonal entry of any pulled-back observable: {offdiag:.1e}")

print("B. Positivity")
print("   cone channels are orthogonal projections on the three-level register, hence positive.")
s0 = lms[:, 2] == 0
r616, r700 = (lms[wl == w][0] for w in (616, 700))
print(f"   lights with tabulated S = 0: {int(wl[s0][0])}-{int(wl[s0][-1])} nm. L/(L+M) at 616 nm = {r616[0]/(r616[0]+r616[1]):.5f}, at 700 nm = {r700[0]/(r700[0]+r700[1]):.5f}.")
print("   A nonzero positive operator on two levels has a kernel of dimension at most one, so two S = 0 lights would")
print("   have to share one colour state and one L:M ratio. They do not, so the S channel cannot be a positive")
print("   operator on a two-level factor alone.")

print("C. Reduced two-level state of each light")
t = a.sum(axis=1)
bx = (a * np.cos(ang)).sum(axis=1) / t            # <sigma_z>, since <v|sz|v> = cos(angle)
by = (a * np.sin(ang)).sum(axis=1) / t
bx_num, by_num = np.diag(Gz) / t, np.diag(Gx) / t
print(f"   Bloch vector from the construction versus opponent formula: max difference {max(np.abs(bx - bx_num).max(), np.abs(by - by_num).max()):.1e}")
print("   opponent formula: x = a_L - (a_M + a_S)/2,  y = (sqrt(3)/2)(a_M - a_S),  both divided by a_L + a_M + a_S")
white = (a * wt[:, None]).sum(axis=0)
print(f"   equal-energy white: Bloch vector ({(white*np.cos(ang)).sum()/white.sum():+.1e}, {(white*np.sin(ang)).sum()/white.sum():+.1e}), i.e. the centre (maximally mixed)")

print("D. Where the lights sit")
r = np.hypot(bx, by)
vis = (wl >= 400) & (wl <= 700)
sel = [400, 420, 440, 460, 480, 490, 500, 520, 540, 560, 570, 580, 590, 600, 620, 650, 700]
print("   purity radius by wavelength: " + "  ".join(f"{w}:{r[wl == w][0]:.2f}" for w in sel))
k = np.argmin(np.where(vis, r, np.inf)); k2 = np.argmax(np.where(vis, r, -np.inf))
print(f"   on 400-700 nm: least pure light {int(wl[k])} nm (radius {r[k]:.3f}); most pure {int(wl[k2])} nm (radius {r[k2]:.3f}); "
      f"largest over the whole table {r.max():.4f} at {int(wl[np.argmax(r)])} nm")
loc = np.argmin(np.where((wl >= 540) & (wl <= 620), r, np.inf))
print(f"   local minimum on the long-wave side (540-620 nm): {int(wl[loc])} nm (radius {r[loc]:.3f})")
pts = np.column_stack([bx, by])
hull = ConvexHull(pts[wl <= 700])
tri = 3 * np.sqrt(3) / 4
print(f"   areas: disc pi = {np.pi:.4f}; receptor triangle {tri:.4f} = {100*tri/np.pi:.1f}% of the disc; "
      f"mixtures of 390-700 nm lights {hull.volume:.4f} = {100*hull.volume/tri:.1f}% of the triangle = {100*hull.volume/np.pi:.1f}% of the disc")
print("   no monochromatic light is a pure state: every light excites at least two cone classes.")
