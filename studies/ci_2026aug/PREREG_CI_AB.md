# Pre-registration — CI-A and CI-B

Filed 2026-08-16 before computation. Successor to
PREREG_canonical_ideal.md.

## Study CI-A — the anomalous quadric on the f1+f2 projection

S = {y1} u {y9..y12}, |S| = 5, dim Sym^2(S) = 15,
dim(I_2 n Sym^2(S)) = 1. Call the unique quadric Q*.

| # | prediction | reasoning |
|---|---|---|
| A1 | Q* lies in the **(+,+)** isotypic class | the class-(-,+) part of Sym^2(S) is y1·(linear in f2); such a product lies in I_2 only if a factor vanishes, and M_*(Gamma_0(143)) is an integral domain |
| A2 | the y1^2 coefficient is **nonzero** | otherwise Q* would lie in I_2 n Sym^2(f2), computed = 0 |
| A3 | consequence of A1+A2: Q* is an identity **a·e1^2 = -q(e9..e12)** in weight 4, a != 0 | |
| A4 | the "excess 1" shrinks under the grading: class (+,+) contributes 11 products into a 12-dim target, so the finding is rank 10 where 11 was generic | |

Open, no prediction (both outcomes reportable):

- **A5** — is the bilinear form G of q self-adjoint for the Hecke action,
  i.e. T^t G = G T with T = T_2 restricted to f2? If yes, q is a trace
  form Tr_{K2/Q}(c x^2) and the identity is Hecke-canonical. If no, it is
  a projection accident. The Hecke algebra does not act on the ideal, so
  there is no a priori reason to expect yes.
- **A6** — rank, determinant and definiteness of q.

Decision rule: if A1 or A2 fails, the domain argument is wrong and
everything downstream is void — re-derive before interpreting.

## Study CI-B — the newspace deficiency of 3

N = the 11 newform coordinates {y1, y2..y7, y9..y12}. Products of
newforms span 33 of the 36 dimensions of H^0(2K). Which 3 are missing?

**I claim this is forced, and predict it before computing.** The newspace
contains **no (-,-) coordinate** (the package constant AL_DECOMPOSITION
already records (-1,-1): 0 for newforms; the single (-,-) line of
S_2(Gamma_0(143)) is oldspace). Class-(-,+) monomials require either
(+,+) x (-,+) or (+,-) x (-,-). With no (-,-) coordinate available, the
newspace can only produce the first kind: y1 x f2, exactly 4 monomials,
against a 7-dimensional target.

| # | prediction |
|---|---|
| B1 | the entire deficiency lies in the **(-,+)** class |
| B2 | per-class deficiency vector = **(0, 0, 3, 0)** against targets (12, 6, 7, 11) |
| B3 | rank of newform products in class (-,+) is exactly **4** (domain argument: e1·e_k independent) |
| B4 | the missing 3 = a complement of e1·f2 inside H^0(2K)_{(-,+)}, supplied by the oldspace monomials y1·y8 and y13·f3 |

Open, no prediction:

- **B5** — which single oldspace direction suffices. Ranks of
  e1·f2 + {y1 y8} (old+ only) and e1·f2 + y13·f3 (old- only) decide it.

Decision rule: if B2 is not (0,0,3,0), the counting argument is wrong;
report the failure and do not repair the argument post hoc.

## Method

All ranks over Q by exact elimination on q-expansion coefficient vectors
truncated at q^140 (weight-4 Sturm bound 56). Q* obtained as an exact
rational nullspace and re-verified two ways: (i) direct substitution into
the q-expansions, residual must be exactly 0; (ii) analytic evaluation at
points tau in H.

## Not claimed

No physics reading. Nothing about mass ratios, generations or couplings
follows from these two studies and none is asserted.
