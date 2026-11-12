"""Symmetry of the image of independent bounded cone activations under linear opponent coordinates.

Inputs   None.  Uses the two formulas printed in the July manuscript (subsec:ordinary-optical-image):
         x = (L - M)/2, y = (L + M)/4 - S/2.  The manuscript does not print the white/black formula.
Checks   For 40 random linear white/black channels (with and without a surround term) acting on independent activations in
         [0, 100]: the body has R(theta, z) = R(theta + pi, -z); the mirror form R(theta, z) = R(theta, -z) fails; where the
         widest band sits.  Also: complementing the activations sends (x, y) to (-x, -y).
Status   Numerical illustration on chosen models.
Scope    (agreed in review) Central symmetry of the attainable body needs two hypotheses the manuscript does not supply:
         every attainable activation triple has an attainable complement, and the height coordinate reverses under it.
         The percentages printed describe these 40 models only and are not a physiological result.
"""
import itertools, numpy as np
from scipy.spatial import ConvexHull
rng = np.random.default_rng(21)
Mxy = np.array([[0.5, -0.5, 0.0], [0.25, 0.25, -0.5]])
def body(wf, ws):
    """Image of focal cube [0,100]^3 (and surround cube if ws is not None) under (x, y, z), z centred."""
    gens = [np.append(Mxy[:, i]*100, wf[i]*100) for i in range(3)]
    if ws is not None: gens += [np.array([0, 0, -ws[i]*100]) for i in range(3)]
    pts = np.array([sum(c*g for c, g in zip(cs, gens)) for cs in itertools.product([0, 1], repeat=len(gens))])
    pts = pts - (pts.max(0) + pts.min(0))/2 * np.array([0, 0, 1])       # centre z so the poles sit at -T and +T
    return ConvexHull(pts), pts[:, 2].max()
def radial(h, theta, z):
    N, C = h.equations[:, :3], h.equations[:, 3]
    u = np.array([np.cos(theta), np.sin(theta), 0.0]); den = N @ u; num = -(N @ np.array([0, 0, z]) + C)
    return np.min(num[den > 1e-12]/den[den > 1e-12])
for label, surround in (("fixed surround", False), ("surround varies too", True)):
    central, mirror_fail, tilt = True, [], []
    for trial in range(40):
        wf = rng.uniform(0.05, 1, 3); ws = rng.uniform(0.05, 1, 3) if surround else None
        h, T = body(wf, ws)
        for _ in range(60):
            t, z = rng.uniform(-np.pi, np.pi), rng.uniform(-0.9*T, 0.9*T)
            a, b, c = radial(h, t, z), radial(h, t + np.pi, -z), radial(h, t, -z)
            central &= np.isclose(a, b, rtol=1e-7); mirror_fail.append(abs(a - c)/max(a, c))
        zs = np.linspace(-0.999*T, 0.999*T, 801)
        for t in np.linspace(-np.pi, np.pi, 24, endpoint=False):
            Rz = np.array([radial(h, t, z) for z in zs]); top = zs[Rz >= Rz.max()*(1 - 1e-9)]
            tilt.append(((top.min() + top.max())/2)/T)
    tilt = np.array(tilt)
    print(f"[{label}] 40 random linear white/black channels: R(theta,z) = R(theta+pi,-z) always: {central}; "
          f"median relative failure of R(theta,z) = R(theta,-z): {np.median(mirror_fail):.0%}; "
          f"centre of the widest band: {tilt.min():+.2f}T .. {tilt.max():+.2f}T, within 0.02T of the equator at {np.mean(np.abs(tilt) < 0.02):.0%} of hues")
# the one printed fact that needs no assumption about z
ok = all(np.allclose(Mxy @ (100 - c), -(Mxy @ c)) for c in rng.uniform(0, 100, (1000, 3)))
print("printed formulas alone: complementing the cone activations sends (x, y) to (-x, -y):", ok)
