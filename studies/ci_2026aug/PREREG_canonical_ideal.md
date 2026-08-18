# Pre-registration — canonical ideal of X0(143)

Filed before any computation. Session: 2026-08-16.

## Object

The canonical embedding X0(143) -> P^12 via the 13-dimensional space
S_2(Gamma_0(143)) = H^0(X, K). Target: the degree-2 part I_2 of the
homogeneous ideal.

## Predictions (theory-side, fixed in advance)

| # | quantity | predicted | source |
|---|---|---|---|
| P1 | dim S_2(Gamma_0(143)) | 13 | genus, Riemann-Hurwitz |
| P2 | dim Sym^2 H^0(K) | 91 | g(g+1)/2 |
| P3 | h^0(2K) | 36 | Riemann-Roch, 3g-3 |
| P4 | dim I_2 = kernel | **55** | 91 - 36 = (g-2)(g-3)/2, Petri |
| P5 | dim Sym^3 H^0(K) | 455 | C(15,3) |
| P6 | h^0(3K) | 60 | 5g-5 |
| P7 | dim I_3 | 395 | 455 - 60, Max Noether |
| P8 | rank of V . I_2 in degree 3 | **395** | Petri generation |
| P9 | Sturm bound, weight 4, level 143 | 56 | k*mu/12 = 4*168/12 |

## Decision rules (fixed in advance)

- **Kernel = 55** -> non-hyperelliptic canonical curve, as predicted.
- **Kernel = 66** -> X0(143) IS hyperelliptic; canonical image is the
  rational normal curve of degree 12, ideal = C(12,2) = 66 quadrics.
  This falsifies the Petri route and is a reportable honest negative.
- **Kernel anything else** -> computation is wrong (insufficient
  q-precision or basis error), not a discovery. Re-run at higher NC
  before interpreting.
- **rank(V . I_2) = 395** -> quadrics generate the ideal (Petri).
- **rank(V . I_2) < 395** -> X0(143) is trigonal or a plane quintic
  (Enriques-Babbage); the 55 quadrics then cut out a scroll or Veronese,
  NOT the curve, and cubic generators are required. Reportable, not a
  failure.

## E2 routes

- **Route A** (kernel): exact rational nullspace over Q of the
  91-column product matrix, `matker` in PARI/GP. Coefficients q^0..q^NC
  with NC >= 120 > 2 x Sturm.
- **Route B** (independent cohomological prediction): rank of the
  13 x 55 = 715 cubics {x_k . Q_m} inside the 455-dim monomial space.
  Tests h^0(3K) = 60, a different Riemann-Roch value than the one
  fixing P4. Computed mod two distinct large primes; rank_Q >= rank_p
  always, and theory caps I_3 at 395, so rank_p = 395 forces equality.
- **Route C** (analytic, not coefficient linear algebra): evaluate the
  13 basis forms numerically at several points tau in H, form the
  canonical image vector, and check all 55 quadrics vanish. A wrong
  kernel would not vanish at generic points.

Route A and Route C share the q-expansion data but not the operation
(coefficient nullspace vs. point evaluation). Route B tests a distinct
Riemann-Roch number. No single route is treated as a certificate alone.

## Not claimed

- No physics reading of the quadrics. This session produces geometry
  only. Any mass / coupling interpretation is out of band and is not
  pre-registered here.
- The Schottky problem is not touched. tau is known to be a Jacobian
  point by construction, not by certification.
