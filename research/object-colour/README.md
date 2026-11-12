# Object-colour research

These research notes are imported from the September 17 archive for shared review. They are separate from the July JMP manuscript and the September 28 condensation; neither manuscript is changed by this import.

- [Proof and prior-art assessment](Magenta_Gap_Proof_and_Prior_Art_2026-09-17.md): the conditional continuum two-contact cubic relation, endpoint qualifications, discrete-versus-continuum warning, invariance scope, and prior-art limits.
- [Proof audit and candidate measurement design](Magenta_Gap_Proof_Audit_and_Test_Design_2026-09-17.md): a dual lower bound, an endpoint counterexample, and numerical test-design estimates in specified spectral models.

The result assumes two distinct interior quadratic supporting contacts and a bounded continuous-spectrum reflectance model. Endpoints require the stated extra conditions. A fixed wavelength grid changes the small-signal limit: exact cone contact can persist until bin capacities bind. Signal level is a white-normalized linear response, not assumed perceptual brightness.

The notes do not establish a Lorentzian perceptual metric, a new empirical law, or priority for the cubic scaling mechanism. Their original numerical-script and package references are retained; the associated code archive has not been imported here. The supplementary audit is explicitly a further assistant check, not an external referee report.

The Markdown contents are preserved as read from the archived notes. `PROVENANCE.json` records their filenames, source versions, byte lengths, and SHA-256 hashes.

## Continuum boundary certificate

The [continuum rim certificate](Continuum_Rim_Certificate.md) gives an explicit non-null boundary point in the full 390–830 nm interpolation model. It uses two independent genuine reflectance transitions to prove uniqueness of a horizontal supporting normal. The criterion needs at least two independent transitions; a third does not invalidate it.

Run from the repository root:

```sh
python research/object-colour/certify_continuum_rim.py
```

The script uses only the Python standard library. It builds the continuous PCHIP model from exact rational knot values and verifies the required signs and root enclosures with rational bounds. The [saved certificate](continuum_rim_certificate.json) records the verification. The [input notes](inputs/README.md) document the supplied table and its license.

This new certificate is separate from the still-unimported September numerical-code archive. It does not certify sampled percentages of the hue circle or identify a perceptual boundary. The manuscript files remain unchanged.
