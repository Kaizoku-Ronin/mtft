# T2 short-horizon approximation: frozen protocol

Date: 2026-09-08. Parent: MTFT_T2_Memory_Approximation_v0.1.0.

Question: can the previously selected four-direction B-SVD model meet a 1%
relative Frobenius active-propagator error budget on a continuous short interval?
This follow-up is motivated by the parent's observed 0.01984% error at t=0.25.
It is not a blinded model-selection experiment. No new subspace is fitted.

Freeze A,B,D from inputs/T2_blocks.npz and Q_svd_r4 from selected_models.npz.
Define H=[[A,B],[B*,D]] and Hr=[[A,BQ],[(BQ)*,herm(Q*DQ)]]. The leading
eight coordinates are the active space. Freeze these explicit Hermitian matrices
as exact dyadic rational complex entries. The certificate concerns this finite
model, not the exact MTFT period geometry or physical time. Q is a numerical
input; Hr is defined directly, without claiming Q is exactly isometric.

Before evaluating new bounds, fix:

- horizons T = 1/8, 1/4, 1/2, 3/4, 1;
- Taylor orders m = 4, 6, 8, 10, 12;
- primary error budget 1/100;
- choose the largest listed T with a proved uniform bound at most 1/100;
- display all orders and all horizons, including failures to certify;
- report measured errors separately; a failed bound is not proof of error;
- retain the parent failure on the long interval and spectral target.

For F(t)=[exp(-itH)]aa and Fr(t)=[exp(-itHr)]aa, bound the finite sum by
sum(k=0..m) T^k ||(H^k)aa-(Hr^k)aa||F/k!.
Hermitian integral remainders bound the two omitted active blocks in Frobenius
norm by sqrt(8) T^(m+1) (rho^(m+1)+rho_r^(m+1))/(m+1)!, where rho and rho_r
are certified operator-norm upper bounds. Obtain each from a Hermitian matrix's
maximum absolute row sum, enclosing each complex magnitude by rational square
root upper bounds. All matrix powers and comparisons use exact integer/rational
arithmetic. Enclose square roots outwards to denominator 10^30.

Divide the absolute bound by a rational lower bound for sqrt(3). Indeed the
full unitary's complement-to-active block has at most rank5, so F has at least
three singular values1 and ||F||F >= sqrt(3) for every real t.

Validation: exact Hermitian checks; zero-th and first moments agree exactly;
independent proof review; exact toy controls (including delayed coupling);
independent numerical expm vs spectral evaluation; saved certificate replay;
1001-point diagnostic grid on [0,1] and each specified endpoint. Grid agreement
does not establish continuity; the exact rational inequality does. Record all
input/source/script hashes and freeze old packages without modifying them.

## Pre-execution mathematical amendment

Independent review supplied a sharper unitary integral remainder before any
new bound was evaluated: replace sqrt(8) rho^(m+1) by ||H^(m+1)E||F, and likewise
for Hr, where E selects the eight active columns. Compute and disclose both
bounds; the minimum of these two valid bounds may certify a listed horizon.
The model, horizons, orders, 1% budget, and sqrt(3) denominator are unchanged.
This follows by applying the unitary integral remainder directly to E and
using that left projection cannot increase the Frobenius norm. It introduces
no fitted parameters and no measured-error-based choice.
