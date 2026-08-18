# The canonical ideal of X0(143)

**Session 2026-08-16 · Route A/B/C complete · pre-registered before computation**

Engine: PARI/GP 2.15.4 (`mfinit([143,2],1)`, `mfcoefs`, `matkerint`,
`mfheckemat`, `mfatkininit`) plus exact `Fraction` linear algebra in
Python. No floating point anywhere in Routes A and B.

---

## Headline

**dim I_2(X0(143)) = 55, exactly as predicted, and the 55 quadrics
generate the full canonical ideal.**

X0(143) sits in P^12 as a canonical curve of degree 2g-2 = 24, cut out
by 55 quadrics. It is neither hyperelliptic nor trigonal nor a plane
quintic. All three exclusions are consequences of the computation, not
inputs to it.

| # | quantity | predicted | computed | class |
|---|---|---|---|---|
| P1 | dim S_2(Gamma_0(143)) | 13 | 13 | EXACT |
| P2 | dim Sym^2 H^0(K) | 91 | 91 | EXACT |
| P3 | h^0(2K) = 3g-3 | 36 | 36 | EXACT |
| P4 | **dim I_2** | **55** | **55** | **EXACT** |
| P5 | dim Sym^3 H^0(K) | 455 | 455 | EXACT |
| P6 | h^0(3K) = 5g-5 | 60 | 60 | EXACT |
| P7 | dim I_3 | 395 | 395 | EXACT |
| P8 | rank(V . I_2) in degree 3 | 395 | 395 | CERTIFIED (two primes) |

## The three routes

**Route A — exact rational nullspace.** q-expansions to q^140 (Sturm
bound at weight 4, level 143 is 56, so this is 2.5x margin). The 91
products f_i f_j span a 36-dimensional space; the kernel is 55.
`matkerint` returns an LLL-reduced integral basis with max |entry| = 676.
Max residual over the kernel basis: **exactly 0**, in exact integer
arithmetic — not a tolerance.

**Route B — cubic generation, an independent Riemann-Roch number.**
The 13 x 55 = 715 cubics x_k . Q_m span a subspace of the 455-dimensional
degree-3 monomial space. Rank = **395** mod 2147483647 and mod 1000003,
recovering h^0(3K) = 60 = 5g-5. Since rank_Q >= rank_p always and theory
caps dim I_3 at 395, equality is forced. This route tests a different
cohomological value than the one fixing P4, and it is the route that
excludes trigonality (Enriques-Babbage).

**Route C — analytic point evaluation.** The 13 forms evaluated at five
points tau in H, canonical image formed, all 55 quadrics checked.
Worst **relative** residual (cancellation measured against the size of
the terms being cancelled, not against 1): **2.34e-41** at dps 40 — the
working-precision floor. A decoy random non-kernel quadric with the same
coefficient range returns 0.132, i.e. O(1). The test is not trivially
satisfied.

## Self-certifying exclusions

- **Hyperelliptic would have forced dim I_2 = C(12,2) = 66**, since the
  canonical image would be the rational normal curve of degree 12.
  We got 55. X0(143) is not hyperelliptic — consistent with Ogg's list,
  now independently confirmed rather than cited.
- **Trigonal or plane quintic would have forced rank(V . I_2) < 395**,
  with the quadrics cutting out a scroll or a Veronese surface instead
  of the curve. We got exactly 395.

## Atkin-Lehner grading of the ideal

W_11 and W_13 are genuine automorphisms of X0(143), so they act on I_2.
(Hecke operators do not — this is why the ideal is graded by
Atkin-Lehner and not by the Hecke algebra.) I_2 is verified stable under
both, and splits:

| sector | dim Sym^2 | dim I_2 | dim H^0(2K) |
|---|---|---|---|
| (+,+) | 38 | 26 | 12 |
| (+,-) | 11 | 5 | 6 |
| (-,+) | 11 | 4 | 7 |
| (-,-) | 31 | 20 | 11 |
| **total** | **91** | **55** | **36** |

E2: computed twice — once by isotypic-projector rank, once by the
character formula on traces restricted to I_2 (tr W_11 = 7, tr W_13 = 5,
tr W_143 = 37). The two routes share no step and agree on all four
sectors.

Quotient genera fall out of the same traces: X0(143)/W_11 has genus 7,
/W_13 genus 6, /W_143 genus 2, and the full quotient X0(143)* has
**genus 1**.

## Galois block x Atkin-Lehner sector table

charpoly(T_2 | S_2) = z (z+2)^2 (z^4-3z^3-z^2+5z+1)
(z^6-10z^4+2z^3+24z^2-7z-12), degrees 1+2+4+6 = 13 — the S_2 shadow of
the ledger's H_1 blocks [2,4,8,12].

| block | dim | (+,+) | (+,-) | (-,+) | (-,-) |
|---|---|---|---|---|---|
| f1 = 143a1 | 1 | 1 | 0 | 0 | 0 |
| oldspace ghost (11a1) | 2 | 0 | 0 | 1 | 1 |
| f2 orbit | 4 | 0 | 0 | 4 | 0 |
| f3 orbit | 6 | 0 | 6 | 0 | 0 |
| **total** | **13** | **1** | **6** | **5** | **1** |

Every Galois orbit is Atkin-Lehner **pure**. This is structurally forced,
not accidental: for sigma in Gal(Qbar/Q), f^sigma | W_Q = sigma(eps_Q)
f^sigma = eps_Q f^sigma because eps_Q = +/-1 lies in Q. A Galois orbit
cannot split across sign combinations.

Removing the oldspace leaves the newform sectors
{(+,+):1, (-,+):4, (+,-):6, (-,-):0}, which **reproduces the package
constant `x0_143.AL_DECOMPOSITION` exactly**, and identifies the extra
(-,+) and (-,-) lines in the full space as the two oldspace directions.
W_143 = +1 on f1 alone, so f1 is the unique newform with root number
eps = -1 — the "electron is unique" claim, re-derived here from
Atkin-Lehner data rather than from L-values.

## Ledger reconciliation

- `FIELD_DISCRIMINANTS`: disc = 1957 (f2) and 194616205 = 5 . 7 . 5560463
  (f3). Both reproduced.
- `FIELD_POLY_F3` equals h6(-z) exactly, as recorded.
- `FIELD_POLY_F2` is **polredabs(charpoly(T_2|f2))**, i.e. z^4-4z^2-z+1,
  not charpoly(T_2|f2) = z^4-3z^3-z^2+5z+1 itself. Same field
  (both polredabs to z^4-4z^2-z+1), same discriminant. The [Mem] shorthand
  "FIELD_POLY_F2 = g4" is imprecise; the stored value is the reduced
  generator. Recording the distinction, not filing a correction — the
  field is right and no downstream number moves.

## Correction filed

**CC-08 (PROPOSED)** — Arithmetica Generale Pr 3.7.5 states that orbit
Omega2 "splits evenly among the four sign combinations" and that Omega3
"has two pairs of conjugate newforms with eigenvalues (+1,-1) and
(-1,+1)". Both clauses are false, and not merely empirically: Atkin-Lehner
eigenvalues are constant on Galois orbits (proof above), so no orbit can
split. Computed truth: Omega2 is uniformly (-1,+1); Omega3 is uniformly
(+1,-1). Omega1 = (+1,+1) as stated. Note that the package constant
`AL_DECOMPOSITION` was already correct — the error is prose-side only,
and the parenthetical "[verified via LMFDB Atkin-Lehner data]" in
Pr 3.7.5 is not supported by what LMFDB actually shows.

## Not claimed

No physics reading of the quadrics. The Schottky problem is untouched —
tau is a Jacobian point by construction, not by certification. Nothing
here is stated about the ideal in a Hecke eigenbasis over the number
fields; the computation is entirely over Q.

## Next moves

1. Re-express the 55 quadrics in an Atkin-Lehner-adapted basis, where the
   26/5/4/20 grading becomes manifest and the quadrics get sparse. This
   is the presentation that would make the ideal readable.
2. `mtft.canonical` — ship the quadric basis and the q-expansions as
   frozen certified data, with P1-P8 as gates, following the
   `jacobian_order` CSV pattern.
3. `mtft.periods` — the higher-value port. tau in H_13 still lives only
   in PARI/GP output.

## Files

- `X0_143_I2_quadric_basis.txt` — 91 x 55 integer matrix, the ideal.
- `X0_143_S2_qexpansions.txt` — the basis the coordinates refer to,
  q^0..q^140. Without this the quadric matrix is meaningless.
- `X0_143_I2_sample_quadrics.txt` — Q1, Q2, Q3 written out.
- `X0_143_canonical_ideal_cert.json` — machine-readable certificate.
- `PREREG_canonical_ideal.md` — pre-registration, filed first.
