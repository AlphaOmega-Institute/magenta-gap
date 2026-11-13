"""Exploratory first pass: luminosity thresholds (Morimoto, Numata, Fukuda &
Uchikawa 2021, J. Vision 21(13):3) against the bounded-reflectance limit at
the same chromaticity, for Experiments 1 and 2.

Usage
    python first_pass.py TABLE_CSV ZENODO_DATA_DIR [S_AXIS_FACTOR]

    TABLE_CSV        research/object-colour/inputs/stockman_sharpe_2deg_390_830.csv
    ZENODO_DATA_DIR  unzipped Data/ of doi:10.5281/zenodo.5590120 (CC BY 4.0)
    S_AXIS_FACTOR    optional; multiplies the S-axis scale of the MODEL only.
                     It deliberately miscalibrates model against data to show
                     how much the correlations depend on that calibration.
                     It is not a property of the physical limit.

Conventions (inferred when this was written; confirmed afterwards against the
authors' own files by exp1_authors_comparison.py, section A)
  * cone fundamentals: Stockman-Sharpe 2 deg, energy units, 390-830 nm, 1 nm
  * luminance weights 0.68990272 L + 0.34832189 M
  * MacLeod-Boynton l = 0.6899 L / V, s = k S / V, equal-energy white at s = 1
  * illuminant: black-body radiator, Planck's law per unit wavelength
  * a 100% reflector has 35 cd/m^2 (paper, Experiments 1 and 2)

What it prints
  per experiment, observer and surround: the correlation between settings and
  the computed limit across test chromaticities, and setting / limit ratios.

Status: exploratory, numerical, finite-bin linear programme.  Not a
certificate.  Two readings first drawn from this output were later corrected;
see README.md in this folder.
"""
import csv
import pathlib
import sys

import numpy as np
from scipy.optimize import linprog

TABLE = pathlib.Path(sys.argv[1])
DATA = pathlib.Path(sys.argv[2])

rows = [r for r in csv.reader(TABLE.read_text().splitlines()) if r]
tab = np.array([[float(v) for v in r] for r in rows[1:]])
wl, lms = tab[:, 0], tab[:, 1:]
w = np.full(len(wl), 1.0)
w[0] = w[-1] = 0.5                      # trapezoid weights, 1 nm

KL, KM = 0.68990272, 0.34832189
ee = (lms * w[:, None]).sum(axis=0)
KS = (KL * ee[0] + KM * ee[1]) / ee[2] * (float(sys.argv[3]) if len(sys.argv) > 3 else 1.0)  # equal-energy white -> s = 1 (optional sensitivity factor)


def planck(T):
    lam = wl * 1e-9
    return lam ** -5 / np.expm1(1.438776877e-2 / (lam * T))


def mb(catch):
    v = KL * catch[0] + KM * catch[1]
    return KL * catch[0] / v, KS * catch[2] / v, v


def limit(l, s, illum):
    """Largest luminance, relative to the 100% reflector, of a reflectance
    0 <= rho <= 1 with MacLeod-Boynton chromaticity (l, s) under illum."""
    A = (lms * (illum * w)[:, None]).T            # 3 x N catches per unit rho
    white = A.sum(axis=1)
    target = np.array([l / KL, (1 - l) / KM, s / KS])   # catches with V = 1
    scale = 1.0 / white.max()
    n = A.shape[1]
    c = np.zeros(n + 1)
    c[-1] = -1.0
    # A rho = t * target, with tau = t * scale as the O(1) LP variable
    Aeq = np.hstack([A * scale, -target[:, None]])
    for method in ("highs", "highs-ds", "highs-ipm"):
        res = linprog(c, A_eq=Aeq, b_eq=np.zeros(3),
                      bounds=[(0, 1)] * n + [(0, None)], method=method)
        if res.status == 0 and res.x is not None:
            break
    else:
        raise RuntimeError("LP failed")
    rho, t = res.x[:-1], res.x[-1] / scale
    vwhite = KL * white[0] + KM * white[1]
    trans = int(np.sum(np.abs(np.diff((rho > 0.5).astype(int)))))
    frac = float(np.mean((rho > 1e-6) & (rho < 1 - 1e-6)))
    return t / vwhite, trans, frac


def read(path):
    """Return {observer: array rows (l, s, settings...)} and the column names."""
    lines = [r for r in csv.reader(path.read_text(encoding="utf-8-sig").splitlines())]
    head = lines[2]
    ncond = (len(head) - 3) // 2
    names = head[3:3 + ncond]
    out, obs = {}, None
    for r in lines[3:]:
        if not r or len(r) < 3 + ncond or r[1] == "":
            continue
        if r[0].strip():
            obs = r[0].strip()
            out[obs] = []
        out[obs].append([float(x) for x in r[1:3 + ncond]])
    return {k: np.array(v) for k, v in out.items()}, names


print("== convention test: black-body locus in this MacLeod-Boynton scaling")
for T in (3000, 5000, 6500, 10000, 20000):
    l, s, _ = mb((lms * (planck(T) * w)[:, None]).sum(axis=0))
    print(f"   T={T:6d} K   l={l:.4f}  s={s:.4f}")

WHITE_CD = 35.0
for fname, T in (("Exp1.csv", 6500), ("Exp2_3000K.csv", 3000),
                 ("Exp2_6500K.csv", 6500), ("Exp2_20000K.csv", 20000)):
    data, names = read(DATA / fname)
    illum = planck(T)
    lw, sw, _ = mb((lms * (illum * w)[:, None]).sum(axis=0))
    print(f"\n== {fname}: illuminant {T} K, white at l={lw:.4f} s={sw:.4f}; conditions {names}")
    first = next(iter(data.values()))
    lim = []
    for l, s in first[:, :2]:
        rel, trans, frac = limit(l, s, illum)
        lim.append(rel)
        print(f"   test l={l:.4f} s={s:.4f}  limit={rel:.4f} of white = {WHITE_CD*rel:6.2f} cd/m2"
              f"  transitions={trans} partial-bin fraction={frac:.4f}")
    lim = np.array(lim)
    for ob, arr in data.items():
        assert np.allclose(arr[:, :2], first[:, :2])
        for j, nm in enumerate(names):
            y = arr[:, 2 + j]
            ratio = y / (WHITE_CD * lim)
            r = np.corrcoef(y, lim)[0, 1]
            print(f"   {ob:3s} {nm:11s} r={r:+.3f}  setting/limit: min={ratio.min():.2f} "
                  f"median={np.median(ratio):.2f} max={ratio.max():.2f}  "
                  f"setting/white: {y.min()/WHITE_CD:.2f}-{y.max()/WHITE_CD:.2f}")
    allr = np.array([arr[:, 2:] for arr in data.values()])          # obs x chrom x cond
    mean = allr.mean(axis=(0, 2))
    print("   mean over observers and conditions, setting/limit by chromaticity:",
          " ".join(f"{m/(WHITE_CD*x):.2f}" for m, x in zip(mean, lim)))
    print(f"   correlation of that mean with the limit: {np.corrcoef(mean, lim)[0,1]:+.3f}")
