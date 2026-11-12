# A certified non-null point of a continuum object-colour boundary

This calculation addresses the continuum smoothness question raised during the November 2026 manuscript audit. It is a mathematical counterexample to an everywhere-null identification in a specified metric class, using a declared interpolation model of standard cone-response data. It is not an empirical identification of a perceptual boundary.

## Model

The input CSV contains the 441 Stockman–Sharpe 2-degree LMS samples from 390 to 830 nm, exported from the archived Colour Science 0.4.7 arrays without changing any binary floating-point value. The accompanying Colour Developers license is retained. The certificate treats those stored binary values as exact rational knot data and constructs the continuous PCHIP interpolant with rational arithmetic. It does not rationalize independently rounded polynomial coefficients.

Write the resulting response density as `s(lambda)` and its three integrals over the full interval as `I_L, I_M, I_S`. Set

\[
c_i(\lambda)=s_i(\lambda)/I_i,\qquad
\mathcal O=\left\{\int_{390}^{830}c(\lambda)\rho(\lambda)\,d\lambda:0\leq\rho\leq1\right\}.
\]

Perfect white is therefore exactly `W=(1,1,1)`. For a catch vector `a`, use

\[
z=\tfrac23(a_L+a_M+a_S)-1,\qquad
x=\tfrac29(2a_L-a_M-a_S),\qquad
y=\tfrac29(-a_L+2a_M-a_S).
\]

These are the object-colour audit's normalized coordinates: `(x,y)/(1+z)` equals `(l-1/3,m-1/3)` away from black.

## The supporting-normal criterion

Let a full-dimensional continuum object-colour solid have continuous response density. Suppose a boundary point is represented by a reflectance with **at least two** genuine interior `0 <-> 1` transitions, and the response vectors at two of those transitions are linearly independent.

For any supporting normal `n`, the same reflectance must maximize the functional `n.a`. Its pointwise maximization rule is `rho=1` where `n.c>0`, and `rho=0` where `n.c<0`. Continuity forces

\[
n\cdot c(\lambda_1)=n\cdot c(\lambda_2)=0.
\]

Two independent vectors have a one-dimensional orthogonal complement. Full dimensionality excludes opposite normal orientations supporting the same point. Thus the outward normal is unique and the boundary is differentiable at that point. This does not assert a continuously differentiable neighborhood. Additional transitions do not invalidate the argument. The argument concerns actual reflectance transitions, not merely two tangencies of a supporting functional to the spectral curve.

## The selected point

Choose the horizontal normal `n=(0,1,0)` in the `(x,y,z)` coordinates and put

\[
q(\lambda)=-c_L(\lambda)+2c_M(\lambda)-c_S(\lambda).
\]

The certificate checks the sign on every interpolation interval, isolates exactly two simple sign-changing roots, and verifies that the positive set is the band between them. Consequently

\[
\rho_*(\lambda)=\mathbf1_{\{q(\lambda)>0\}}
\]

is an exact support optimizer. Approximate root locations are 486.863081 and 599.369448 nm; the executable certificate reports rigorous enclosing intervals.

The corresponding response point is enclosed by

\[
\begin{aligned}
0.11649647725622&\leq x_*\leq0.11649647725683,\\
0.20194541250715&\leq y_*\leq0.20194541250768,\\
0.09328313559491&\leq z_*\leq0.09328313559570.
\end{aligned}
\]

Thus the certified boundary point is strictly above the equator. Its positive `y` coordinate also gives a positive pairing between the supporting normal and the radial direction.

At a root, `2c_M=c_L+c_S`, and the response density in the transformed coordinates is

\[
\left((c_L-c_S)/3,\;0,\;c_L+c_S\right).
\]

The certificate bounds `c_L-c_S` strictly below zero at the first root and strictly above zero at the second; `c_L+c_S` is positive at both. The two transformed vectors are therefore independent. Together with any response vector inside the positive band, they span three dimensions; continuity therefore makes the solid full dimensional. The unique outward normal direction at the corresponding boundary point is horizontal.

## Consequence for the metric model

Let \(W\) denote the fixed achromatic **direction** (the white catch vector in the affine signal space), and let \(n\) be the nonzero supporting conormal at the certified point. The exact geometric fact is

\[
n(W)=0.
\]

In the displayed coordinates the linear part of the coordinate map sends this direction to \((0,0,2)\); the affine white point itself has coordinates \((0,0,1)\). The unique tangent plane therefore contains the achromatic direction.

For **any** Lorentzian metric of signature \((2,1)\) for which that direction is timelike at this point,

\[
g(W,W)<0\quad\Longrightarrow\quad g^{-1}(n,n)>0.
\]

Indeed, if \(n^\sharp=g^{-1}n\), then \(g(n^\sharp,W)=n(W)=0\). The orthogonal complement of a timelike vector in a Lorentzian space is positive definite. Since \(n^\sharp\) is nonzero, its squared norm is positive.

Thus the boundary cannot be null there. The previously stated class

\[
g=G_{ab}\,dx^a dx^b-C\,dz^2,\qquad G>0,\quad C>0,
\]

is a special case. **Adding spatial–height cross-terms cannot remove the obstruction while the achromatic direction remains timelike.** Making the boundary null would require that direction to be null or spacelike there, a necessary rather than sufficient change.

The result is independent of the choice of a valid linear height functional, as long as height variation at fixed linear chromatic coordinates remains parallel to \(W\). It is also unchanged by any mere coordinate transformation when the metric, direction, and conormal are transported together. Holding projective chromaticity fixed instead varies intensity along the current catch vector, generally not along \(W\). Declaring that different direction timelike changes the geometric identification; it is not a coordinate-only escape.

No finite-difference smoothness test, linear-programming tolerance, or sampled percentage of the hue circle enters this conclusion. A nonlinear response model or a perceptual identification requires separate analysis. The calculation neither changes nor refutes the conditional bicone compatibility theorem.

## Complementarity and classical full colours

The certified transition-vector signs prove more than approximate coplanarity. Their transformed vectors have opposite nonzero chromatic \(x\) components and positive height components. A positive combination therefore lies along the positive achromatic direction. Consequently the two spectral response directions are complementary relative to this model's white. This follows exactly from the existing interval certificate; a floating-point residual is unnecessary.

This construction belongs to the classical study of Ostwald's full colours or semichromes and Schrödinger optimal colours. Koenderink, *The prior statistics of object colors*, JOSA A **27**, 206–217 (2010), p. 213, relates the tangent plane's parallelism to the achromatic axis to complementary transition wavelengths: https://doi.org/10.1364/JOSAA.27.000206. Koenderink, van Doorn and Gegenfurtner, *RGB Colors and Ecological Optics* (2021), discuss the complementary-band definition and the semichrome locus in Sections 4–5: https://doi.org/10.3389/fcomp.2021.630370.

These references identify the classical construction; they are not used to impose a universal two-transition restriction. Modern object-colour analysis allows additional transitions. The normal-uniqueness criterion needs at least two independent genuine transitions and permits more. The particular verified interpolation-model point and its use to test the manuscript's metric identification are the scope of this certificate, with no claim that full colours or their supporting-plane geometry are new.

## Verification and scope

Run `python research/object-colour/certify_continuum_rim.py` from the repository root. The script uses only the Python standard library. Its conclusions are checked with exact rational arithmetic and outward rational interval bounds; decimal numbers in the report summarize those bounds.

The bundled input is finite, but the response model is a continuous piecewise-cubic density and permits every measurable reflectance between zero and one. It is not a finite-bin reflectance model. No assertion is made that the interpolation exactly equals biological sensitivity functions.

The published supporting-plane framework is described by Logvinenko, Funt and Bastani, *The object-color solid*, Journal of Vision 25(2):2 (2025), especially the pointwise optimization rule in Equations (5)–(7): https://doi.org/10.1167/jov.25.2.2. The normal-uniqueness argument above is an elementary consequence of that framework; no novelty claim is made.

The September two-contact cubic theorem concerns supporting tangencies and a small-signal asymptotic. This certificate concerns the normal to the continuum solid at a finite-signal reflectance with genuine transitions. The two statements answer different questions.
