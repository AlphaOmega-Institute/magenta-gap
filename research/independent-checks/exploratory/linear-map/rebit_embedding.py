"""Exploratory: can one fixed linear map carry spectra into a two-level
(rebit) colour factor, and what does the cone table then force?

Setting.  A stimulus is a positive combination of monochromatic lights.  A
linear map followed by expectation values of fixed colour observables gives
three functions of wavelength that are a linear image of the cone
fundamentals:  (t, x, y)(lambda) = M (L, M, S)(lambda).  The reduced colour
state of a monochromatic light is  t/2 * [I + (x/t) s1 + (y/t) s2], which is
a positive matrix exactly when (x/t, y/t) lies in the unit disc.  So a map
exists exactly when the spectral chromaticities fit inside an ellipse of the
chosen chromaticity plane; the ellipse becomes the rim of pure states.

This script
  1. confirms that the spectral chromaticities cannot all be pure states
     (three collinear distinct points cannot lie on one ellipse),
  2. computes minimum-area enclosing ellipses in two chromaticity planes,
     with free centre and with centre at equal-energy white,
  3. prints the resulting purity radius of each monochromatic light.

Usage:  python rebit_embedding.py TABLE_CSV
Status: exploratory, numerical.  The choice of plane and of ellipse is not
unique; nothing here is fitted to appearance data.
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
integ = (lms * wt[:, None]).sum(axis=0)
np.set_printoptions(precision=4, suppress=True, linewidth=150)


def mvee(P, centre=None, tol=1e-10, itmax=200000):
    """Minimum-area ellipse {(p-c)^T Q (p-c) <= 1} containing the rows of P.
    With centre given, the ellipse is constrained to that centre."""
    n, d = P.shape
    if centre is None:
        Qm = np.column_stack([P, np.ones(n)])
        D = d + 1
    else:
        Qm = P - centre
        D = d
    u = np.full(n, 1.0 / n)
    for _ in range(itmax):
        Mx = Qm.T @ (Qm * u[:, None])
        g = np.einsum("ij,jk,ik->i", Qm, np.linalg.inv(Mx), Qm)
        j = int(np.argmax(g)); k = g[j]
        if k <= D * (1 + tol):
            break
        beta = (k - D) / (D * (k - 1))
        u *= (1 - beta); u[j] += beta
    if centre is None:
        c = P.T @ u
        Q = np.linalg.inv(P.T @ (P * u[:, None]) - np.outer(c, c)) / d
    else:
        c = np.asarray(centre, float)
        Q = np.linalg.inv(Qm.T @ (Qm * u[:, None])) / d
    return c, Q


def radius(P, c, Q):
    V = P - c
    return np.sqrt(np.einsum("ij,jk,ik->i", V, Q, V))


def report(name, chrom, white, wl):
    hull = ConvexHull(chrom)
    verts = chrom[hull.vertices]
    print("=" * 96)
    print(name)
    for label, centre in (("free centre", None), ("centre at equal-energy white", white)):
        c, Q = mvee(verts, centre)
        r = radius(chrom, c, Q)
        rw = float(radius(white[None, :], c, Q)[0])
        # area of the gamut hull as a fraction of the disc, in ellipse-normalised coordinates
        Lc = np.linalg.cholesky(Q)
        Z = (verts - c) @ Lc
        frac = ConvexHull(Z).volume / np.pi
        touch = wl[r > 1 - 1e-6]
        print(f"  -- {label}: largest radius {r.max():.6f}; white sits at radius {rw:.3f}; "
              f"mixture gamut covers {100*frac:.1f}% of the disc")
        print(f"     lights on the rim (pure states), nm: {', '.join(str(int(w)) for w in touch)}")
        vis = (wl >= 420) & (wl <= 680)
        k = np.argmin(np.where(vis, r, np.inf))
        print(f"     least pure monochromatic light between 420 and 680 nm: {int(wl[k])} nm, radius {r[k]:.3f}")
        sel = [400, 440, 460, 480, 500, 520, 540, 560, 570, 580, 590, 600, 620, 650, 700]
        print("     radius by wavelength: " + "  ".join(f"{w}:{r[wl == w][0]:.2f}" for w in sel))
        # the purple line: mixtures of the two ends of the 390-700 range
        a, b = chrom[wl == 400][0], chrom[wl == 700][0]
        mid = radius(((a + b) / 2)[None, :], c, Q)[0]
        print(f"     midpoint of the 400-700 nm mixture line: radius {mid:.2f}")
        if centre is not None:
            # saturation relative to white along the spectrum: where is it smallest?
            k2 = np.argmin(np.where((wl >= 400) & (wl <= 700), r, np.inf))
            print(f"     smallest radius on 400-700 nm: {int(wl[k2])} nm ({r[k2]:.3f})")


# 1. collinearity: three distinct collinear chromaticities cannot be distinct points of an ellipse
bal = lms / integ
tsum = bal.sum(axis=1)
chromA = bal[:, :2] / tsum[:, None]                     # balanced simplex coordinates (l, m)
s0 = lms[:, 2] == 0
P0 = chromA[s0]
resid = np.abs(np.cross(P0[-1] - P0[0], P0 - P0[0])).max()
print(f"1. Wavelengths with tabulated S = 0: {int(wl[s0][0])}-{int(wl[s0][-1])} nm ({s0.sum()} samples); "
      f"they lie on one line (max deviation {resid:.1e}) and are {len(np.unique(np.round(P0, 12), axis=0))} distinct points.")
print("   An ellipse meets a line in at most two points, so no linear image puts all of them on the rim.")
print("   The same holds without the exact zeros wherever the spectral locus is locally straight to within its")
print("   own curvature; the minimum-area ellipses below quantify how far inside the rim the lights sit.")

# 2. two chromaticity planes
whiteA = np.array([1 / 3, 1 / 3])
report("2a. Plane normalised by the sum of white-balanced cone catches (equal-energy white at the centroid)", chromA, whiteA, wl)

KL, KM = 0.68990272, 0.34832189
lum = KL * lms[:, 0] + KM * lms[:, 1]
KS = (KL * integ[0] + KM * integ[1]) / integ[2]
chromB = np.column_stack([KL * lms[:, 0] / lum, KS * lms[:, 2] / lum])
whiteB = np.array([KL * integ[0] / (KL * integ[0] + KM * integ[1]), 1.0])
keep = wl <= 700     # beyond 700 nm S = 0 and the points add nothing new in this plane
report("2b. Plane normalised by luminance (MacLeod-Boynton l, s), 390-700 nm", chromB[keep], whiteB, wl[keep])
