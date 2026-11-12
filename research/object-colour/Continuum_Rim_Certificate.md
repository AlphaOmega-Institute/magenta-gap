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

For every metric of the form

\[
g=G_{ab}\,dx^a dx^b-C\,dz^2,\qquad G>0,\quad C>0,
\]

with no spatial–height cross-terms, the squared inverse-metric norm of this normal is

\[
g^{-1}(n,n)=(G^{-1})_{yy}>0.
\]

Hence the boundary at this point is not null. No finite-difference smoothness test, linear-programming tolerance, or sampled percentage of the hue circle enters this conclusion.

This is an obstruction for this response model in these specified coordinates and metric class. Allowing spatial–height cross-terms changes the hypothesis. A nonlinear response model or a perceptual identification requires separate analysis. The calculation neither changes nor refutes the conditional bicone compatibility theorem.

## Verification and scope

Run `python research/object-colour/certify_continuum_rim.py` from the repository root. The script uses only the Python standard library. Its conclusions are checked with exact rational arithmetic and outward rational interval bounds; decimal numbers in the report summarize those bounds.

The bundled input is finite, but the response model is a continuous piecewise-cubic density and permits every measurable reflectance between zero and one. It is not a finite-bin reflectance model. No assertion is made that the interpolation exactly equals biological sensitivity functions.

The published supporting-plane framework is described by Logvinenko, Funt and Bastani, *The object-color solid*, Journal of Vision 25(2):2 (2025), especially the pointwise optimization rule in Equations (5)–(7): https://doi.org/10.1167/jov.25.2.2. The normal-uniqueness argument above is an elementary consequence of that framework; no novelty claim is made.

The September two-contact cubic theorem concerns supporting tangencies and a small-signal asymptotic. This certificate concerns the normal to the continuum solid at a finite-signal reflectance with genuine transitions. The two statements answer different questions.
