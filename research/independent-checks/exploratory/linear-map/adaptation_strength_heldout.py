"""Exploratory: adaptation strength fitted on some experiments and scored on
others, with the encoding held fixed.

Same encoding, rules and matching criterion as adaptation_rules.py.  Each
rule gets one strength parameter D between no adaptation (0) and the full
transformation (1):

  R-lin   cone gains  D * (w2/w1) + (1 - D)          (the usual incomplete von Kries form)
  R-geo   cone gains  (w2/w1) ** D
  F       two-level filter  S2^(D/2) S1^(-D/2), S_k the unadapted state of white k

Held-out designs, using the viewing conditions recorded with the data:
  A. the eight phases measured at 1000 lux on a 20% background, from five
     experiments: fit D on four experiments, score the fifth
  B. within each of the three remaining experiments (different illuminance
     ranges and media): fit D on the other phases, score the left-out phase

Usage:  python adaptation_strength_heldout.py DATA_DIR
Status: exploratory.  Phases are not independent studies; several share
observers, samples or apparatus.  Chromaticity only.
"""
import pathlib
import sys

import numpy as np

DATA = pathlib.Path(sys.argv[1])
HPE = np.array([[0.38971, 0.68898, -0.07868], [-0.22981, 1.18340, 0.04641], [0.0, 0.0, 1.0]])
iH = np.linalg.inv(HPE)
ang = np.deg2rad([0.0, 120.0, 240.0])
nvec = np.column_stack([np.cos(ang), np.sin(ang)])
vvec = np.column_stack([np.cos(ang / 2), np.sin(ang / 2)])
P = np.array([np.outer(x, x) for x in vvec])
FAMILY = {"CSAJ.da": ("CSAJ-C", 1000), "helson.ca": ("Helson", 1000), "lam.da": ("Lam & Rigg", 1000),
          "lutchi": ("Lutchi", 1000), "Kuo": ("Kuo & Luo", 1000),
          "Brene": ("Breneman", None), "RIT": ("Braun & Fairchild", 129), "mcan": ("McCann", None)}


def uv(x):
    x = np.atleast_2d(x); d = x[:, 0] + 15 * x[:, 1] + 3 * x[:, 2]
    return np.column_stack([4 * x[:, 0] / d, 9 * x[:, 1] / d])


def state(l):
    s = np.einsum("ni,ijk->njk", l, P)
    return s / np.trace(s, axis1=1, axis2=2)[:, None, None]


def weights(s):
    b = np.column_stack([s[:, 0, 0] - s[:, 1, 1], 2 * s[:, 0, 1]])
    return 1 / 3 + (2 / 3) * b @ nvec.T


def mpow(S, p):
    w, U = np.linalg.eigh(S)
    return U @ np.diag(w ** p) @ U.T


def load(path):
    r = [ln.split() for ln in path.read_text().splitlines() if ln.strip()]
    w = np.array(r[0], float); n = int(r[1][0]); d = np.array(r[2:2 + n], float)
    return w[:3], w[3:], d[:, :3], d[:, 3:]


phases = []
for p in sorted(DATA.glob("*.dat")):
    w1, w2, x1, x2 = load(p)
    if np.hypot(*(uv(w1) - uv(w2))[0]) < 0.005:
        continue
    fam = next(v for k, v in FAMILY.items() if p.name.startswith(k))
    phases.append(dict(name=p.name, family=fam[0], lux=fam[1], w1=w1, w2=w2, x1=x1, x2=x2))


def errors(ph, rule, D):
    obs = uv(ph["x2"]); l1, lw1, lw2 = ph["x1"] @ HPE.T, HPE @ ph["w1"], HPE @ ph["w2"]
    if rule == "R-lin":
        pred = (l1 * (D * (lw2 / lw1) + (1 - D))) @ iH.T
    elif rule == "R-geo":
        pred = (l1 * (lw2 / lw1) ** D) @ iH.T
    else:
        S1, S2 = state(lw1[None])[0], state(lw2[None])[0]
        M = mpow(S2, 0.5 * D) @ mpow(S1, -0.5 * D)
        s2 = np.einsum("ij,njk,lk->nil", M, state(l1), M)
        s2 /= np.trace(s2, axis1=1, axis2=2)[:, None, None]
        pred = weights(s2) @ iH.T
    return np.hypot(*(uv(pred) - obs).T), int((pred[:, 2] < 0).sum())


GRID = np.round(np.linspace(0.40, 1.20, 81), 2)
RULES = ("R-lin", "R-geo", "F")
cache = {(ph["name"], r, D): errors(ph, r, D) for ph in phases for r in RULES for D in GRID}


def fit(train, rule):
    tot = [np.concatenate([cache[(ph["name"], rule, D)][0] for ph in train]).mean() for D in GRID]
    return GRID[int(np.argmin(tot))], float(np.min(tot))


print(f"{len(phases)} phases with an illuminant change, {sum(len(p['x1']) for p in phases)} pairs.\n")
print("In-sample: one strength per rule over all phases")
for r in RULES:
    D, e = fit(phases, r)
    full = np.concatenate([cache[(ph["name"], r, 1.0)][0] for ph in phases]).mean()
    print(f"   {r:6s} best D = {D:.2f}  mean error {e:.5f}   (at D = 0.85: "
          f"{np.concatenate([cache[(ph['name'], r, 0.85)][0] for ph in phases]).mean():.5f};  at full strength: {full:.5f})")


def held_out(units, label):
    """units: list of (test phases, training phases)."""
    print(f"\n{label}")
    out = {r: [] for r in RULES}; negz = 0; wins = 0; rows = []
    for test, train in units:
        line = []
        for r in RULES:
            D, _ = fit(train, r)
            e = np.concatenate([cache[(ph["name"], r, D)][0] for ph in test])
            out[r].append(e); line.append((D, e.mean()))
            if r == "F":
                negz += sum(cache[(ph["name"], r, D)][1] for ph in test)
        wins += line[0][1] < line[2][1]
        rows.append((test, line))
    for test, line in rows:
        tag = test[0]["family"] if len({t["family"] for t in test}) == 1 and len(test) > 1 else test[0]["name"]
        print(f"   held out {tag:18s} n={sum(len(t['x1']) for t in test):3d}  " +
              "  ".join(f"{r} D={D:.2f} err {m:.4f}" for r, (D, m) in zip(RULES, line)))
    for r in RULES:
        allv = np.concatenate(out[r])
        print(f"   {r:6s} held-out mean {allv.mean():.5f}  median {np.median(allv):.5f}  ({len(allv)} pairs)")
    print(f"   R-lin beats F in {wins} of {len(units)} held-out units; F predictions with negative Z at the fitted strengths: {negz}")
    return out


lux1000 = [p for p in phases if p["lux"] == 1000]
fams = sorted({p["family"] for p in lux1000})
A = held_out([([p for p in lux1000 if p["family"] == f], [p for p in lux1000 if p["family"] != f]) for f in fams],
             "A. 1000 lux, 20% background: leave one experiment out (five experiments, eight phases)")
B = {}
for fam in ("Breneman", "Braun & Fairchild", "McCann"):
    ps = [p for p in phases if p["family"] == fam]
    B[fam] = held_out([([p], [q for q in ps if q is not p]) for p in ps], f"B. {fam}: leave one phase out ({len(ps)} phases)")

print("\nAll held-out pairs together")
for r in RULES:
    allv = np.concatenate(A[r] + [x for fam in B for x in B[fam][r]])
    print(f"   {r:6s} mean {allv.mean():.5f}  median {np.median(allv):.5f}  ({len(allv)} pairs)")
