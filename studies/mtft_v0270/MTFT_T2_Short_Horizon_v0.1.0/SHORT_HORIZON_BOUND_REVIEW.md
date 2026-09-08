# Independent review: finite-matrix short-horizon bounds

This review concerns a frozen numerical finite-dimensional model. It does not certify the original period geometry, identify physical time, or establish a particle Hamiltonian.

## Coordinate convention

The prior input `T2_blocks.npz` contains an ambient-coordinate `H`, together with numerical frames `U` and `F`. Its first eight coordinates are not the registered active sector. For an exact dyadic calculation define a new, explicit block-coordinate model

\[
 H_f=\begin{pmatrix}A&B\\B^*&D\end{pmatrix},\qquad
 H_r=\begin{pmatrix}A&BQ\\(BQ)^*&\operatorname{herm}(Q^*DQ)\end{pmatrix}.
\]

Here the matrices assembled and saved by the implementation are interpreted entry by entry as exact dyadic complex numbers. Check exact Hermiticity, including real diagonals. The original numerical frame and selected `Q` have floating-point orthogonality errors. These do not obstruct a theorem about the explicitly saved Hermitian matrices, but they prevent calling this a certified exact compression of the original ambient arithmetic operator. Agreement of the models with the parent study is a separate numerical control.

## Absolute Frobenius bound

Let `E_f` and `E_r` select the first eight standard coordinate vectors in dimensions 13 and 12 respectively. Put

\[
 M_k=E_f^*H_f^k E_f-E_r^*H_r^k E_r.
\]

For an order `m`, horizon `T >= 0`, and all real `|t| <= T`,

\[
\|E_f^*e^{-itH_f}E_f-E_r^*e^{-itH_r}E_r\|_F
\le N_m(T),
\]

where

\[
N_m(T)=\sum_{k=0}^{m}\frac{T^k}{k!}\|M_k\|_F
+\frac{T^{m+1}}{(m+1)!}
\left(\|H_f^{m+1}E_f\|_F+\|H_r^{m+1}E_r\|_F\right).
\]

Proof: the matrix exponential has the integral remainder

\[
R_m(H,t)=\frac{(-i)^{m+1}}{m!}\int_0^t(t-s)^m e^{-isH}H^{m+1}\,ds.
\]

For Hermitian `H`, `exp(-isH)` is unitary. Multiplication on the left by `E*` has operator norm one, so the Frobenius norm of the compressed integrand is at most `||H^(m+1) E||F`. Integrating yields the remainder used above. Negative time follows by reversing the integration variable, or by conjugate transposition. The finite Taylor sum follows from the triangle inequality. No exponential growth factor is necessary for real time and Hermitian matrices.

This active-column remainder is usually much sharper than bounding `||H||` by a row sum or Frobenius norm and raising that norm to a high power. Either bound is valid, and their minimum is valid. An order-12 certificate requires powers through order 13.

## Relative bound and exact decision

Write the full unitary in active/complement blocks. Unitarity gives

\[
U_{aa}^*U_{aa}=I_8-U_{ca}^*U_{ca}.
\]

The second term has rank at most five and eigenvalues in `[0,1]`. Consequently at least three singular values of `U_aa` are exactly one and `||U_aa||F >= sqrt(3)` at all real times. Therefore the relative Frobenius error is at most `N_m(T)/sqrt(3)` uniformly on `[-T,T]`.

If exact rational square-root upper enclosures make `N` rational, the one-percent decision can be made with the rational comparison **`N^2 <= 3/10000`**. Decimal approximations of `sqrt(3)` are unnecessary in the pass/fail calculation. A rounded decimal is only a display value.

An optional, separately declared improvement follows from Duhamel's formula:

\[
\|U_{ca}(t)\|_F\le |t|\|B\|_F,\qquad
\|U_{aa}(t)\|_F\ge\sqrt{\max(3,8-T^2\|B\|_F^2)}.
\]

The registered primary comparison may retain the simpler `sqrt(3)` denominator.

## Exact arithmetic checks and pitfalls

1. Convert each real and imaginary float by its exact integer ratio; conversion through a short decimal changes the certified matrix.
2. Include the full complex Frobenius norm: sum real-entry squares plus imaginary-entry squares.
3. If `H = Z / 2^s`, use integer complex multiplication and divide its kth power by `2^(s*k)`. Comparing full and reduced moments requires a common denominator.
4. A rational upper enclosure for `sqrt(n/d)` can use `ceil(sqrt(n*d))/d`, or a stated fixed rational grid. Prove the squared enclosure is at least `n/d`; do not assume floating-point rounding points outward.
5. The exact agreement of moments zero and one follows from identical active `A`. The familiar `M_2 = B_s B_s*` additionally requires an exact orthogonal retained/discarded decomposition. With floating-point `Q`, compute `M_2` directly.
6. A bound above one percent is inconclusive; it does not establish that the actual error exceeds one percent.
7. A certificate valid for all times in a short interval does not give a resolvent approximation or a longer-horizon guarantee.
8. This is a bound for active amplitude propagators in a direct sum. It is not a statement about density matrices or a tensor-product partial trace.

## Independent controls

The accompanying script evaluates the bound with 80-digit arithmetic from exact float ratios, compares it with a finite sampled evolution, verifies the dimension-based denominator, and checks Taylor identities in exact symbolic toy models. The sampled checks are numerical diagnostics; the general proof and a correct rational implementation supply the continuous-interval guarantee.

## Completed implementation review

The primary `dyadic_bounds.py` implementation matches the preceding proof for the supplied frozen matrices. Its use of a rational lower enclosure of `sqrt(3)` is conservative. The row-norm remainder and active-column remainder are both valid, so selecting their minimum preserves the bound.

The independent `verify_exact_certificate.py` implementation uses SymPy's `DomainMatrix` over the Gaussian rational field `QQ_I`, rather than the primary implementation's integer matrix pairs and shared powers of two. It verified all 14 active-moment norm enclosures and all 25 active-column tail and relative-bound rows using rational squared inequalities. A separately constructed rational bound proves the order-12, `T = 1/2` one-percent inequality directly. Eight invalid-input controls also passed.

An 80-digit calculation independently gives the order-12 relative bound at `T = 1/2` as approximately `0.0060190630922425304033`. Twelve exact symbolic toy controls and 161 sampled evolution controls passed. These sampled checks are supplementary and do not replace the exact interval inequality.

The first toy-control execution compared two unsimplified SymPy matrix expressions using structural equality, producing one false mismatch in the third-moment identity. Replacing that comparison by simplification of their difference resolved the mismatch; the experiment rules and numerical matrices were unchanged.

For future public JSON intake, `decode` should reject nonintegral numeric entries instead of coercing them with `int`. The actual frozen models contain integer strings and are unaffected by this input-hardening observation.
