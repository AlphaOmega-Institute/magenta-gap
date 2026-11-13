"""Exploratory: two ways an adapting condition could act on the linear map's
colour state, compared on public corresponding-colour data.

Colour state (see cone_register_and_trine.py): for white-balanced cone
catches a_i, sigma = sum_i a_i P_i / sum_i a_i, with P_i three rim states of
a two-level system at 120 degree Bloch spacing.

Rule R  "reweight the cone register": multiply the three catches by gains
        and keep the rim states fixed.  With gains 1/(catch of the adapting
        white) this is the classical von Kries rule.
Rule F  "filter the two-level factor": sigma -> T sigma T / trace.
        F-white: T sends the adapting white to the disc centre.
        F-gain:  T reproduces the same three cone gains as rule R.

A corresponding-colour pair is two stimuli that look alike under two
adapting conditions.  Each rule predicts the second stimulus's chromaticity
from the first by equating colour states, so no calibration of hue or
saturation scales is needed.  Error: distance in CIE 1976 u'v'.

Usage:  python adaptation_rules.py DATA_DIR
        DATA_DIR holds the .dat files of doi:10.5281/zenodo.3270903
        (Luo and Rhodes 1999; licence not specified, so not redistributed).

Status: exploratory first pass.  Complete adaptation is assumed for every
rule; cone signals are Hunt-Pointer-Estevez, equal-energy normalised; the
disc centre is equal-energy white.  CAT02 is included only as a familiar
benchmark.
"""
import pathlib
import sys

import numpy as np

DATA = pathlib.Path(sys.argv[1])
np.set_printoptions(precision=4, suppress=True, linewidth=150)

HPE = np.array([[0.38971, 0.68898, -0.07868],
                [-0.22981, 1.18340, 0.04641],
                [0.0, 0.0, 1.0]])
CAT02 = np.array([[0.7328, 0.4296, -0.1624],
                  [-0.7036, 1.6975, 0.0061],
                  [0.0030, 0.0136, 0.9834]])
ang = np.deg2rad([0.0, 120.0, 240.0])
nvec = np.column_stack([np.cos(ang), np.sin(ang)])
vvec = np.column_stack([np.cos(ang / 2), np.sin(ang / 2)])
P = np.array([np.outer(x, x) for x in vvec])


def uv(xyz):
    xyz = np.atleast_2d(xyz)
    d = xyz[:, 0] + 15 * xyz[:, 1] + 3 * xyz[:, 2]
    return np.column_stack([4 * xyz[:, 0] / d, 9 * xyz[:, 1] / d])


def state(lms):                      # lms: (n,3) -> (n,2,2), trace one
    s = np.einsum("ni,ijk->njk", lms, P)
    return s / np.trace(s, axis1=1, axis2=2)[:, None, None]


def weights(sig):                    # barycentric weights of a state in the fixed triangle
    b = np.column_stack([sig[:, 0, 0] - sig[:, 1, 1], 2 * sig[:, 0, 1]])
    return 1 / 3 + (2 / 3) * b @ nvec.T


def msqrt(S, power=0.5):
    w, U = np.linalg.eigh(S)
    return U @ np.diag(w ** power) @ U.T


def gain_filter(g):
    G = (4 / 3) * np.einsum("i,ijk->jk", g, P) - g.sum() / 3 * np.eye(2)
    return G, np.linalg.eigvalsh(G).min()


def load(path):
    rows = [ln.split() for ln in path.read_text().splitlines() if ln.strip()]
    w = np.array(rows[0], float)
    n = int(rows[1][0])
    d = np.array(rows[2:2 + n], float)
    return w[:3], w[3:], d[:, :3], d[:, 3:]


# ---- algebra of the two rules, independent of any data
print("Algebra of the two rules")
g0 = np.array([2.0, 1.0, 1.0]); T0 = np.diag([np.sqrt(2.0), np.sqrt(2.0 / 3.0)])
bl = lambda s: np.array([s[0, 0] - s[1, 1], 2 * s[0, 1]]) / np.trace(s)
print(f"   example: the filter diag(sqrt 2, sqrt(2/3)) gives cone gains {[round(float(x @ T0 @ T0 @ x), 12) for x in vvec]};"
      f" balanced white goes to {bl(np.einsum('i,ijk->jk', g0, P)).round(12)} under reweighting and {bl(T0 @ P.sum(axis=0) @ T0).round(12)} under the filter")
rng = np.random.default_rng(7)
worst, n_ok, n_no = 0.0, 0, 0
for _ in range(20000):
    g = rng.uniform(0.05, 3, 3)
    G, mn = gain_filter(g)
    bR = (g[:, None] * nvec).sum(axis=0) / g.sum()
    assert np.allclose([x @ G @ x for x in vvec], g)
    assert (mn >= -1e-12) == (np.hypot(*bR) <= 0.5 + 1e-12)
    if mn >= -1e-12:
        n_ok += 1
        Tm = msqrt(G) if mn > 0 else None
        worst = max(worst, np.abs(bl(G) - 2 * bR).max())
        if Tm is not None:
            pr = rng.dirichlet([1, 1, 1]); pp = g * pr / (g * pr).sum()
            moved = np.array([bl(Tm @ Pi @ Tm) for Pi in P])
            worst = max(worst, np.abs(bl(Tm @ np.einsum("i,ijk->jk", pr, P) @ Tm) - pp @ moved).max(), np.abs(np.hypot(*moved.T) - 1).max())
    else:
        n_no += 1
print(f"   over 20000 random gain triples: a positive filter with those gains exists for {n_ok} and not for {n_no};")
print("   it exists exactly when the reweighted white has radius at most 1/2, it is then unique, it sends white twice as far")
print(f"   as reweighting does, and it equals the same new weights placed on moved rim states (largest deviation {worst:.1e}).")
print()

names = sorted(p.name for p in DATA.glob("*.dat"))
print(f"{len(names)} files.  Columns: mean u'v' error of the predicted second chromaticity; 'out' = predictions that leave the receptor triangle.")
print(f"{'data set':16s} {'n':>4s} {'white shift':>11s} | {'none':>7s} {'R=vonK':>7s} {'F-white':>8s} {'out':>4s} {'F-gain':>8s} {'out':>4s} | {'CAT02':>7s}")
tot = {k: [] for k in ("none", "R", "Fw", "Fg", "cat")}
used = 0
for name in names:
    w1, w2, x1, x2 = load(DATA / name)
    shift = float(np.hypot(*(uv(w1) - uv(w2))[0]))
    if shift < 0.005:
        print(f"{name:16s} {len(x1):4d} {shift:11.4f} | same white chromaticity; not an illuminant change, skipped")
        continue
    used += 1
    obs = uv(x2)
    l1, lw1, lw2 = x1 @ HPE.T, HPE @ w1, HPE @ w2
    e = {}
    e["none"] = np.hypot(*(uv(x1) - obs).T)
    # Rule R: von Kries in this cone space
    pred = (l1 * (lw2 / lw1)) @ np.linalg.inv(HPE).T
    e["R"] = np.hypot(*(uv(pred) - obs).T)
    # Rule F-white: T_k = sigma(white_k)^(-1/2)
    S1, S2 = state(lw1[None])[0], state(lw2[None])[0]
    Mx = msqrt(S2, 0.5) @ msqrt(S1, -0.5)
    sig2 = np.einsum("ij,njk,lk->nil", Mx, state(l1), Mx)
    sig2 = sig2 / np.trace(sig2, axis1=1, axis2=2)[:, None, None]
    p2 = weights(sig2)
    out_w = int((p2.min(axis=1) < 0).sum())
    e["Fw"] = np.hypot(*(uv(p2 @ np.linalg.inv(HPE).T) - obs).T)
    # Rule F-gain: filters reproducing the von Kries gains 1/white
    G1, m1 = gain_filter(1 / lw1 * lw1.mean())
    G2, m2 = gain_filter(1 / lw2 * lw2.mean())
    if min(m1, m2) > 0:
        My = msqrt(G2, -0.5) @ msqrt(G1, 0.5)
        sg = np.einsum("ij,njk,lk->nil", My, state(l1), My)
        sg = sg / np.trace(sg, axis1=1, axis2=2)[:, None, None]
        pg = weights(sg)
        out_g = int((pg.min(axis=1) < 0).sum())
        e["Fg"] = np.hypot(*(uv(pg @ np.linalg.inv(HPE).T) - obs).T)
        fg = f"{e['Fg'].mean():8.4f} {out_g:4d}"
    else:
        e["Fg"] = None
        fg = f"{'no filter':>8s} {'':>4s}"
    # benchmark
    c1 = x1 @ CAT02.T
    e["cat"] = np.hypot(*(uv((c1 * ((CAT02 @ w2) / (CAT02 @ w1))) @ np.linalg.inv(CAT02).T) - obs).T)
    for k in tot:
        if e[k] is not None:
            tot[k].append(e[k])
    print(f"{name:16s} {len(x1):4d} {shift:11.4f} | {e['none'].mean():7.4f} {e['R'].mean():7.4f} {e['Fw'].mean():8.4f} {out_w:4d} {fg} | {e['cat'].mean():7.4f}")

print(f"\n{used} data sets with an illuminant change.")
for k, lab in (("none", "no adaptation"), ("R", "rule R (von Kries)"), ("Fw", "rule F-white"), ("Fg", "rule F-gain (where a filter exists)"), ("cat", "CAT02 benchmark")):
    allv = np.concatenate(tot[k])
    print(f"   {lab:38s} pairs {len(allv):4d}   mean {allv.mean():.4f}   median {np.median(allv):.4f}   90th pct {np.percentile(allv, 90):.4f}")
# paired comparison on the same pairs
r, f = np.concatenate(tot["R"]), np.concatenate(tot["Fw"])
print(f"   pairs where rule R is closer than F-white: {np.mean(r < f)*100:.1f}%;  data sets where its mean is lower: "
      f"{sum(a.mean() < b.mean() for a, b in zip(tot['R'], tot['Fw']))} of {len(tot['R'])}")
