# Cone-response table for the continuum certificate

`stockman_sharpe_2deg_390_830.csv` contains 441 samples at 1 nm spacing, with columns `wavelength_nm,L,M,S`. They are the Stockman–Sharpe 2-degree fundamentals retained in the September 28, 2026 source-and-calculations archive, originally obtained from Colour Science 0.4.7.

The CSV export preserves every value of the archived NumPy arrays exactly when parsed as binary64. No resampling, cutoff change, or rounding to fewer significant digits was applied.

CSV SHA-256: `144fc58d2f3b01cd8ce5ae765199fd1a2af658b7bea7a540aa3cea6bb9d02c40`.

For the certificate, these parsed binary64 knot values are declared exact rational numbers. The continuous response density is then constructed with exact rational PCHIP slopes and Hermite polynomials. Integrals of those continuous polynomials supply the white normalization. This is a declared interpolation model of a tabulated standard observer, not a claim of exact biological data.

The Colour Developers license supplied with the archived arrays is reproduced in `COLOUR_SCIENCE_LICENSE.txt`.
