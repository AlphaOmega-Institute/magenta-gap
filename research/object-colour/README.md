# Object-colour geometry

This folder contains an exact-arithmetic boundary certificate and two separate research notes on small-signal cubic contact.

## Continuum boundary certificate

- [`Continuum_Rim_Certificate.md`](Continuum_Rim_Certificate.md): model, proof and attribution.
- [`certify_continuum_rim.py`](certify_continuum_rim.py): executable rational-arithmetic verification.
- [`continuum_rim_certificate.json`](continuum_rim_certificate.json): saved enclosing intervals and conclusions.
- [`inputs/`](inputs/): cone-response table, its hash, interpolation convention and license.

Run from the repository root; only the Python standard library is needed:

```sh
python research/object-colour/certify_continuum_rim.py
```

The certificate constructs a continuous PCHIP response density from the 390–830 nm table and proves that one boundary point has a unique tangent plane containing the achromatic direction. That boundary cannot be null for any Lorentzian metric keeping this direction timelike, including metrics with spatial–height cross-terms. At least two independent genuine reflectance transitions suffice for the uniqueness argument; additional transitions are allowed.

The result concerns the declared bounded-reflectance signal model. It does not identify a perceptual boundary, certify a percentage of the hue circle, or assert a smooth neighborhood around the certified point. The current [manuscript](../../manuscript/Magenta_Gap.tex) includes the argument and its scope.

## Cubic-contact research

- [Two-contact proof and prior-art assessment](Magenta_Gap_Proof_and_Prior_Art_2026-09-17.md): the conditional cubic asymptotic, endpoint qualifications, invariance and discrete-versus-continuum distinction.
- [Dual lower bound and endpoint counterexample](Magenta_Gap_Proof_Audit_and_Test_Design_2026-09-17.md): an independent lower-bound argument, a polynomial counterexample and candidate measurement calculations.

These notes contain distinct mathematics beyond the current paper. The cubic result assumes two distinct interior quadratic supporting contacts in a bounded continuous-spectrum reflectance model. Endpoints require additional conditions; fixed wavelength bins change the small-signal limit. The signal coordinate is a white-normalized linear response, not measured perceptual brightness.

The notes' auxiliary numerical scripts are not bundled here. [`recompute/`](../../recompute/) implements the paper's angular Fisher-bound illustration, not those cubic calculations. The notes establish no new empirical law or priority for the cubic scaling mechanism.
