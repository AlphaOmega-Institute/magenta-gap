"""Exploratory: Experiment 1 (6500 K) of Morimoto, Numata, Fukuda & Uchikawa
(2021, J. Vision 21(13):3) against bounded-reflectance limits.

Usage
    python exp1_authors_comparison.py TABLE_CSV ZENODO_DATA_DIR AUTHORS_REPO_DIR

    TABLE_CSV         research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv
    ZENODO_DATA_DIR   unzipped Data/ of doi:10.5281/zenodo.5590120 (CC BY 4.0)
    AUTHORS_REPO_DIR  clone of github.com/takuma929/LuminosityThreshold (MIT),
                      commit f857dc7480e0a5d3c202d946f536537f8a8ec010

What it checks
    A. conventions: the authors' spectrum-locus file against our cone table,
       and their black-body locus against ours
    B. the authors' own prediction rule (nearest stored optimal colour) and the
       per-observer correlations it gives for Experiment 1
    C. our constrained optimisation against their stored optimal colours,
       point by point
    D. one multiplicative factor per observer versus an additive offset and an
       affine fit, with residuals by chromaticity

Status: exploratory and numerical.  Four observers.  Nothing here is a
certificate, a calibration of a perceptual boundary, or a test of the
manuscript's null-boundary hypothesis.
"""
import csv
import pathlib
import sys

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.optimize import linprog

TABLE, ZEN, REPO = (pathlib.Path(a) for a in sys.argv[1:4])
WHITE_CD = 35.0                       # cd/m^2 of a 100% reflector (paper; Lum_k in the authors' code)
np.set_printoptions(precision=4, suppress=True, linewidth=160)


def mat(name):
    m = loadmat(REPO / name, squeeze_me=True, struct_as_record=False)
    return m[[k for k in m if not k.startswith("__")][0]]


def lp_limit(A, target):
    """max t with A rho = t * target, 0 <= rho <= 1.  A is 3 x N (catches per
    unit reflectance in each bin); target is the catch vector of unit luminance.
    Rows are normalised by the white catch and the LP variable is the limit
    relative to the white luminance, so every coefficient is of order one."""
    white = A.sum(axis=1)
    V = white[0] + white[1]
    n = A.shape[1]
    c = np.zeros(n + 1)
    c[-1] = -1.0
    Aeq = np.hstack([A / white[:, None], -(target * V / white)[:, None]])
    for method in ("highs", "highs-ds", "highs-ipm"):
        res = linprog(c, A_eq=Aeq, b_eq=np.zeros(3),
                      bounds=[(0, 1)] * n + [(0, None)], method=method)
        if res.status == 0 and res.x is not None:
            break
    else:
        raise RuntimeError("LP failed")
    rho = res.x[:-1]
    trans = int(np.sum(np.abs(np.diff((rho > 0.5).astype(int)))))
    return res.x[-1] * V, trans


# ---------------------------------------------------------------- A. conventions
print("=" * 100)
print("A. Conventions")
rows = [r for r in csv.reader(TABLE.read_text().splitlines()) if r]
tab = np.array([[float(v) for v in r] for r in rows[1:]])
wl, lms = tab[:, 0], tab[:, 1:]
KL, KM = 0.68990272, 0.34832189
wtrap = np.ones(len(wl)); wtrap[0] = wtrap[-1] = 0.5
ee = (lms * wtrap[:, None]).sum(axis=0)
KS = (KL * ee[0] + KM * ee[1]) / ee[2]              # equal-energy white -> s = 1
ours = np.column_stack([KL * lms[:, 0], KM * lms[:, 1], KS * lms[:, 2]])   # L', M', S' with V = L'+M'

loc = mat("SpectrumLocus_MB.mat")                    # wavelength, l, s, L+M
assert np.array_equal(loc[:, 0], wl)
theirs = np.column_stack([loc[:, 1] * loc[:, 3], (1 - loc[:, 1]) * loc[:, 3], loc[:, 2] * loc[:, 3]])
vis = (wl >= 400) & (wl <= 700)
g = (theirs[vis, 0] + theirs[vis, 1]).sum() / (ours[vis, 0] + ours[vis, 1]).sum()
for name, j in (("L'", 0), ("M'", 1), ("S'", 2)):
    ok = vis & (ours[:, j] > 1e-6 * ours[:, j].max())
    ratio = theirs[ok, j] / (g * ours[ok, j])
    print(f"   authors' {name} / ours (global scale removed), 400-700 nm: "
          f"min {ratio.min():.5f}  median {np.median(ratio):.5f}  max {ratio.max():.5f}")


def planck(T, lam_nm):
    lam = lam_nm * 1e-9
    return lam ** -5 / np.expm1(1.438776877e-2 / (lam * T))


bbl = mat("MB_bbl_500to25000with500step.mat")
for T in (3000, 6500, 20000):
    c = (ours * (planck(T, wl) * wtrap)[:, None]).sum(axis=0)
    mine = (c[0] / (c[0] + c[1]), c[2] / (c[0] + c[1]))
    th = bbl[T // 500 - 1]
    print(f"   black body {T:5d} K: ours l={mine[0]:.5f} s={mine[1]:.5f} | authors' file l={th[0]:.5f} s={th[1]:.5f}"
          f" | s ratio {mine[1]/th[1]:.5f}")

# their optimal-colour set and illuminant
OP = mat("OP_MB/OP_6500_MB.mat")                     # l, s, L+M relative to white; last row is white
ill = mat("ill6500K.mat")                            # 400..700 step 10, relative energy
w3 = wl[vis]
print(f"   stored optimal colours: {len(OP)} rows = 2*C(301,2)+1 = {2*(301*300//2)+1};"
      f" stored white l={OP[-1,0]:.6f} s={OP[-1,1]:.6f}")
cands = {"linear": np.interp(w3, ill[:, 0], ill[:, 1]),
         "pchip": PchipInterpolator(ill[:, 0], ill[:, 1])(w3),
         "spline": CubicSpline(ill[:, 0], ill[:, 1])(w3),
         "planck6500": planck(6500, w3)}
best = None
for nm, e in cands.items():
    c = (theirs[vis] * e[:, None]).sum(axis=0)
    l, s = c[0] / (c[0] + c[1]), c[2] / (c[0] + c[1])
    err = np.hypot(20 * (l - OP[-1, 0]), s - OP[-1, 1])
    print(f"   white from authors' locus x illuminant ({nm:10s}), 400-700 sum: l={l:.6f} s={s:.6f}  weighted distance to stored white {err:.2e}")
    if best is None or err < best[0]:
        best = (err, nm, e)
print(f"   -> using '{best[1]}' interpolation of the authors' illuminant for the 'their conventions' optimisation")
E3 = best[2]
A_theirs = (theirs[vis] * E3[:, None]).T             # 3 x 301, plain sums as in a 1-nm sample grid
white_theirs = A_theirs.sum(axis=1)
V_theirs = white_theirs[0] + white_theirs[1]
A_ours = (ours * (planck(6500, wl) * wtrap)[:, None]).T
V_ours = A_ours.sum(axis=1)[:2].sum()


def target(l, s):
    return np.array([l, 1 - l, s])                   # L', M', S' with L'+M' = 1


# ---------------------------------------------------------------- B. authors' rule and data
print("=" * 100)
print("B. Authors' prediction rule and Experiment 1 data")
D = loadmat(REPO / "LT_Data.mat", squeeze_me=True, struct_as_record=False)["LT_Data"].Exp1
obs_mat = [str(o) for o in D.Key.observer]
dist_mat = [str(d) for d in D.Key.distribution]
test = np.asarray(D.TestChromaticity.Natural.ill6500K, float)            # 9 x 2
for dname in dist_mat:
    assert np.allclose(getattr(D.TestChromaticity, dname).ill6500K, test)
Y = np.stack([np.asarray(getattr(D.Results, d).ill6500K, float) for d in dist_mat], axis=2)   # chrom x obs x dist
SE = np.stack([np.asarray(getattr(D.SE, d).ill6500K, float) for d in dist_mat], axis=2)

# Zenodo CSV cross-check
lines = [r for r in csv.reader((ZEN / "Exp1.csv").read_text(encoding="utf-8-sig").splitlines())]
names_csv = lines[2][3:8]
zen, ob = {}, None
for r in lines[3:]:
    if len(r) < 8 or r[1] == "":
        continue
    if r[0].strip():
        ob = r[0].strip(); zen[ob] = []
    zen[ob].append([float(x) for x in r[1:13]])
zen = {k: np.array(v) for k, v in zen.items()}
csv_obs = list(zen)
diff = np.array([[np.abs(zen[co][:, 2:7] - Y[:, i, :]).max() for i in range(len(obs_mat))] for co in csv_obs])
match = diff.argmin(axis=1)                          # repository column carrying each CSV observer's numbers
assert sorted(match) == list(range(len(obs_mat))) and diff.min(axis=1).max() < 1e-9
labels = [None] * len(obs_mat)
for co, i in zip(csv_obs, match):
    labels[i] = co if co == obs_mat[i] else f"{co}*"
    assert np.abs(zen[co][:, 7:12] - SE[:, i, :]).max() < 1e-9
print(f"   distributions: repository {dist_mat} = Zenodo {names_csv} (same order, identical numbers)")
print(f"   observer columns: repository key {obs_mat}; the Zenodo CSV attaches these numbers to {[csv_obs[list(match).index(i)] for i in range(len(obs_mat))]}")
if any(l.endswith("*") for l in labels):
    print("   NOTE the two public sources label two observers differently; the numbers are identical.")
    print("        Below, observers carry the Zenodo CSV label, starred where the repository key differs.")
obs_mat = labels
print(f"   test chromaticities: CSV versus repository differ by at most {max(np.abs(zen[c][:, :2] - test).max() for c in csv_obs):.1e}")

sel = np.array([np.argmin(np.hypot(20 * (OP[:, 0] - l), OP[:, 1] - s)) for l, s in test])
pred_auth = OP[sel, 2]                                # relative to white
r_auth = np.array([[np.corrcoef(Y[:, i, j], pred_auth)[0, 1] for j in range(len(dist_mat))] for i in range(len(obs_mat))])
print("   correlation, observer settings versus authors' nearest-sample prediction (rows observers, columns distributions):")
print(r_auth)
print(f"   mean over 4 observers x 5 distributions: {r_auth.mean():.4f}   (paper reports 0.901)")

# ---------------------------------------------------------------- C. point by point
print("=" * 100)
print("C. Point-by-point: stored optimal colours versus constrained optimisation (limits relative to white)")
print("   columns: test l, s | nearest stored sample l, s, Y | LP(authors' conv.) at test | LP(authors' conv.) at the sample |"
      " LP(ours: 390-830, Planck) at test")
lp_t, lp_s, lp_o = [], [], []
for (l, s), k in zip(test, sel):
    a, _ = lp_limit(A_theirs, target(l, s)); a /= V_theirs
    b, _ = lp_limit(A_theirs, target(OP[k, 0], OP[k, 1])); b /= V_theirs
    c, tr = lp_limit(A_ours, target(l, s)); c /= V_ours
    lp_t.append(a); lp_s.append(b); lp_o.append(c)
    print(f"   {l:.4f} {s:.4f} | {OP[k,0]:.4f} {OP[k,1]:.4f} {OP[k,2]:.4f} | {a:.4f} | {b:.4f} | {c:.4f}")
lp_t, lp_s, lp_o = map(np.array, (lp_t, lp_s, lp_o))
print(f"   LP at the stored sample's own chromaticity minus the stored luminance: min {np.min(lp_s-pred_auth):+.5f} max {np.max(lp_s-pred_auth):+.5f}")
print(f"   LP(authors' conv.) at the test chromaticity / authors' nearest-sample prediction: min {np.min(lp_t/pred_auth):.4f} max {np.max(lp_t/pred_auth):.4f}")
print(f"   LP(ours) / LP(authors' conv.) at the test chromaticity:                           min {np.min(lp_o/lp_t):.4f} max {np.max(lp_o/lp_t):.4f}")
print(f"   LP(ours) / authors' nearest-sample prediction:                                    min {np.min(lp_o/pred_auth):.4f} max {np.max(lp_o/pred_auth):.4f}")

rng = np.random.default_rng(20261113)
idx = rng.choice(len(OP) - 1, size=400, replace=False)
idx = idx[OP[idx, 2] > 0.02]
d, kept, failed = [], [], 0
for k in idx:
    try:
        d.append(lp_limit(A_theirs, target(OP[k, 0], OP[k, 1]))[0] / V_theirs - OP[k, 2]); kept.append(k)
    except RuntimeError:
        failed += 1; fail_s = OP[k, 1]
idx = np.array(kept)
rel = np.array(d) / OP[idx, 2]
print(f"   (solver failures in the random-sample check: {failed}" + (f", at a stored colour with s = {fail_s:.3g})" if failed else ")"))
print(f"   random stored optimal colours (n={len(idx)}, luminance > 0.02): LP limit at their chromaticity minus stored luminance,"
      f" relative: min {rel.min():+.2e}  median {np.median(rel):+.2e}  95th pct {np.percentile(rel,95):+.2e}  max {rel.max():+.2e}")
print(f"   fraction within 1e-3 relative: {np.mean(np.abs(rel) < 1e-3):.3f};  fraction where the LP exceeds the stored value by more than 1%: {np.mean(rel > 1e-2):.3f}")
beaten = idx[rel > 1e-2]
print(f"   the {len(beaten)} stored colours beaten by more than 1% have s between {OP[beaten, 1].min():.3g} and {OP[beaten, 1].max():.3g}"
      f" and l between {OP[beaten, 0].min():.3f} and {OP[beaten, 0].max():.3f};"
      f" the test chromaticities have l {test[:, 0].min():.3f}-{test[:, 0].max():.3f}, s {test[:, 1].min():.3f}-{test[:, 1].max():.3f}")
print("   (s = 0 is the region where the S fundamental is tabulated as zero; chromaticity is then one-dimensional and a")
print("    two-transition reflectance is not in general the brightest at its chromaticity.  This does not touch Experiment 1.)")

for nm, p in (("authors' nearest sample", pred_auth), ("LP, authors' conventions", lp_t), ("LP, our conventions", lp_o)):
    r = np.array([[np.corrcoef(Y[:, i, j], p)[0, 1] for j in range(len(dist_mat))] for i in range(len(obs_mat))])
    print(f"   mean correlation over 20 observer-distributions using {nm:26s}: {r.mean():.4f}")

# ---------------------------------------------------------------- D. one factor per observer
print("=" * 100)
print("D. One multiplicative factor per observer, against an additive offset and an affine fit")
print("   y = setting (cd/m2), m = 35 x predicted limit.  Weighted by the reported SE (20 repetitions).")
print("   chi2/dof uses those SEs only, so it also absorbs between-surround differences; read it comparatively.")


def wls(X, y, w):
    Xw, yw = X * w[:, None], y * w
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    cov = np.linalg.inv(Xw.T @ Xw)
    res = y - X @ beta
    return beta, np.sqrt(np.diag(cov)), res, float(((res * w) ** 2).sum())


for nm, p in (("authors' nearest-sample predictions", pred_auth), ("LP, our conventions", lp_o)):
    print(f"\n   --- model limits: {nm}")
    m9 = WHITE_CD * p
    for i, o in enumerate(obs_mat):
        y = Y[:, i, :].ravel(order="F"); se = SE[:, i, :].ravel(order="F")
        m = np.tile(m9, len(dist_mat)); w = 1 / se; n = len(y)
        bM, sM, rM, cM = wls(m[:, None], y, w)
        bA, sA, rA, cA = wls(np.ones((n, 1)), y - m, w)
        bF, sF, rF, cF = wls(np.column_stack([m, np.ones(n)]), y, w)
        rmsM, rmsA, rmsF = (np.sqrt(np.mean(r ** 2)) for r in (rM, y - m - bA[0], rF))
        print(f"   {o}: multiplicative k={bM[0]:.3f}+-{sM[0]:.3f} chi2/dof={cM/(n-1):6.2f} rms={rmsM:5.2f} |"
              f" additive c={bA[0]:+6.2f}+-{sA[0]:.2f} chi2/dof={cA/(n-1):6.2f} rms={rmsA:5.2f} |"
              f" affine a={bF[0]:.3f}+-{sF[0]:.3f} b={bF[1]:+6.2f}+-{sF[1]:.2f} chi2/dof={cF/(n-2):6.2f} rms={rmsF:5.2f}")
        ratio = Y[:, i, :] / m9[:, None]                      # chrom x dist
        mean_c = ratio.mean(axis=1); sd_c = ratio.std(axis=1, ddof=1) / np.sqrt(ratio.shape[1])
        dev = (mean_c - bM[0]) / bM[0]
        slope = np.polyfit(test[:, 0], mean_c / mean_c.mean(), 1)[0]
        r_l = np.corrcoef(test[:, 0], mean_c)[0, 1]
        r_m = np.corrcoef(m9, mean_c)[0, 1]
        print("       setting/limit by chromaticity (mean over 5 surrounds): " + " ".join(f"{v:.2f}" for v in mean_c))
        print("       its SE across surrounds:                               " + " ".join(f"{v:.2f}" for v in sd_c))
        print(f"       deviation from k: {' '.join(f'{100*v:+.0f}%' for v in dev)};  correlation of the ratio with l: {r_l:+.2f},"
              f" with the limit itself: {r_m:+.2f};  spread of the ratio across chromaticities (max/min): {mean_c.max()/mean_c.min():.2f}")
    # per observer x surround factors
    k_all = np.array([[wls((WHITE_CD * p)[:, None], Y[:, i, j], 1 / SE[:, i, j])[0][0] for j in range(len(dist_mat))]
                      for i in range(len(obs_mat))])
    print("   factor fitted separately for each observer (rows) and surround (columns):")
    print(k_all)
