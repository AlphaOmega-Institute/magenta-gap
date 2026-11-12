"""A conormal that annihilates a timelike direction is spacelike.

Inputs   None.
Checks   200,000 random signature-(2,1) metrics for which a fixed direction W is timelike, almost all with cross-terms:
         g^-1(n, n) > 0 for conormals with n(W) = 0.  A null W admits a null conormal annihilating it.
Status   Numerical sanity check of a three-line proof (research/object-colour/Continuum_Rim_Certificate.md).
"""
import numpy as np
rng = np.random.default_rng(42); checked = 0; worst = np.inf; cross = 0
while checked < 200000:
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3))); g = Q @ np.diag([rng.uniform(0.1, 5), rng.uniform(0.1, 5), -rng.uniform(0.1, 5)]) @ Q.T   # any signature (2,1) metric
    W = np.array([0.0, 0.0, 1.0])
    if W @ g @ W >= 0: continue                                  # keep only metrics for which the achromatic direction is timelike
    cross += abs(g[0, 2]) > 0.05 or abs(g[1, 2]) > 0.05          # chromatic-height cross-terms present
    phi = rng.uniform(0, 2*np.pi); n = np.array([np.cos(phi), np.sin(phi), 0.0])      # conormal with n(W) = 0
    val = n @ np.linalg.inv(g) @ n; worst = min(worst, val/np.abs(np.linalg.eigvalsh(np.linalg.inv(g))).max()); checked += 1
print(f"{checked} random signature-(2,1) metrics with the achromatic direction timelike ({cross} with sizeable cross-terms):")
print(f"  smallest g^-1(n,n) for a conormal annihilating that direction, relative to the metric scale: {worst:.2e} (positive: {worst > 0})")
# the converse the note states: with the direction null or spacelike, a null conormal annihilating it can exist
g = np.diag([1.0, 1.0, -1.0]); Wn = np.array([1.0, 0.0, 1.0]); n = np.array([1.0, 0.0, -1.0])
print("null achromatic direction (1,0,1) in the flat metric: conormal (1,0,-1) annihilates it:", n @ Wn == 0, "and is null:", n @ np.linalg.inv(g) @ n == 0)
