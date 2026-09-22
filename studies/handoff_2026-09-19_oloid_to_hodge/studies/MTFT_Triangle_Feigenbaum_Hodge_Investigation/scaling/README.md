# Reproducible quadratic-map scaling calculation

Run `python compute_scaling.py` in this directory. Only Python's standard
library is required. The program rewrites the adjacent CSV and JSON files.

## What was computed

For `f_c(z)=z²+c`, the main real period-doubling cascade has superstable centers
`c_n` satisfying `f_c^(2^n)(0)=0`. Here `c_0=0`, `c_1=-1`, and `n=10` is period
1024. These centers are different from the parameter values at which an
attracting periodic orbit loses stability and bifurcates.

The two reported scaling quantities are

```
delta_n = (c_(n-1)-c_(n-2))/(c_n-c_(n-1)), n >= 2
d_n = f_(c_n)^(2^(n-1))(0), n >= 1
alpha_n = abs(d_(n-1)/d_n), n >= 2.
```

The signed spatial ratio is also recorded; it is negative. The absolute value
convention matches the positive constant approximately 2.502907875.

The affine conjugacy to the logistic family is

```
z = r(1/2-x)
c = r(2-r)/4
r = 1+sqrt(1-4c), taking the relevant real branch.
```

For the triangular polynomial family `F_r(n)=r n(n+1)`, the complete coordinate
change is `x=-n`, then `z=r(n+1/2)`. For `r != 0` this is invertible. The
scaling calculation applies after choosing to iterate this parameterized
family; it does not derive the parameter or the iteration from an invariant.

## Numerical evidence and its limits

Every result was computed independently at Decimal precisions of 60 and 90
digits. The roots use sign-change scanning followed by bisection. The search
moves left from the preceding center over 0.05 to 0.4 times the preceding
parameter gap, in steps of 0.0025 times that gap. It selects the first detected
sign change, thereby avoiding the old center, which remains a root of the
higher iterate. This branch selection is numerically successful for all
computed orders. The scan is not a proof that every root has been enumerated.

The halfway return is checked to be nonzero. Since every proper divisor of
`2^n` divides `2^(n-1)`, a certified nonzero halfway return at an exact root
would establish primitive period; our numerical checks support that condition.

Decimal arithmetic here is not outward-rounded interval arithmetic. Reported
brackets, residuals, and precision agreement are numerical diagnostics, not
rigorous enclosures or proofs of the infinite-order limits. The small
precision-rerun error must not be confused with the much larger remaining
finite-period error of each scaling ratio.

At period 1024:

| Quantity | Numerical value |
|---|---|
| `c_10` | -1.4011547825466178412186125377 |
| corresponding logistic `r_10` | 3.5699453554864685808928023823 |
| `delta_10` | 4.6691951560300171740211088012 |
| `alpha_10` | 2.5029075715122790772751516450 |
| critical return residual magnitude | about 3.85e-77 |
| 60/90-digit discrepancy in `c_10` | about 1.46e-49 |
| 60/90-digit discrepancy in `delta_10` | about 2.21e-42 |
| 60/90-digit discrepancy in `alpha_10` | about 5.10e-43 |

No new universality theorem, value of a physical coupling, or Hodge-theoretic
identification is asserted by these numerical results.
