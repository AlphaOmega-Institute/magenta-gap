"""Boosts do not preserve white, and the size of the region beyond a chord of the unit disc.

Inputs   None.
Checks   A boost moves the unit effect and pushes a legal effect outside the body; rotations keep white.  The fraction of
         the disc beyond a chord subtending phi is (phi - sin phi)/(2 pi), by Monte Carlo.
Status   Exact statements, checked numerically.
Corrections (agreed in review)
         * The region beyond a mixture chord is a convex restriction on preparations.  It is not a linear (kernel) defect:
           the spectral lights span the whole space.
         * That boosts move white shows white is extra data; it says nothing about why physiology selects it.
"""
import numpy as np
I2 = np.eye(2); sx = np.array([[0., 1.], [1., 0.]]); sz = np.array([[1., 0.], [0., -1.]])
P = lambda th: (I2 + np.cos(th)*sz + np.sin(th)*sx)/2
# 1. which determinant-preserving maps E -> S E S^T keep the white point (the unit effect)?
eta = 0.8; boost = np.diag([np.exp(eta/2), np.exp(-eta/2)])
c, s = np.cos(0.37), np.sin(0.37); rot = np.array([[c, -s], [s, c]])
print("boost keeps white:", np.allclose(boost @ boost.T, I2), "| rotation keeps white:", np.allclose(rot @ rot.T, I2))
E = 0.9*P(0.0) + 0.1*P(np.pi)        # a legal colour (eigenvalues 0.9, 0.1)
print("boosted colour eigenvalues:", np.round(np.linalg.eigvalsh(boost @ E @ boost.T), 3), "(> 1 means it left the colour body)")
# 2. spectral lights cover only an arc of the hue circle; mixtures of lights fill the convex hull of that arc.
#    Fraction of the unit-trace disc NOT reachable by mixing lights, as a function of the missing arc phi.
rng = np.random.default_rng(1)
pts = rng.uniform(-1, 1, size=(400000, 2)); pts = pts[(pts**2).sum(1) <= 1]
for phi_deg in (60, 90, 120):
    phi = np.radians(phi_deg)
    # missing arc centred on angle pi: chord is the line x = -cos(phi/2)
    unreachable = np.mean(pts[:, 0] < -np.cos(phi/2))
    print(f"missing arc {phi_deg:3d} deg: unreachable fraction MC = {unreachable:.4f}, formula (phi - sin phi)/(2 pi) = {(phi - np.sin(phi))/(2*np.pi):.4f}, max saturation on the chord midpoint = {np.cos(phi/2):.3f}")
