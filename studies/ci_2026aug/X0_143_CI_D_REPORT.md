# CI-D: P located. The descent is fully explicit.

**Session 2026-08-16. Pre-registered in `PREREG_CI_D.md`. All four
predictions confirmed, and the result is verified by a second route
sharing no step with the first.**

---

## Headline

**L_-- = O_E(P) with P = (4, -7), a generator of the Mordell-Weil group
of 143a1.**

E = X0(143)* is the elliptic curve 143a1: y^2 + y = x^3 - x^2 - x - 2,
rank 1, trivial torsion, so E(Q) = Z. The degree-1 line bundle carrying
the level-11 ghost line is the one attached to a generator of that group.

## Route 1 — analytic (modular parametrisation of the CM points)

- **D2 confirmed.** The modular degree of 143a1 is **4**, exactly the
  degree of X0(143) -> X0(143)*. So the AL quotient map *is* the modular
  parametrisation, and X0(143)* = 143a1 on the nose, not merely up to
  isogeny.
- The W13-fixed points solve 143c t^2 - 26a t - b = 0 with a + d = 0,
  giving t = (13a +- sqrt(-13))/(143c) — discriminant **-52**, as the
  class-number route said independently.
- Evaluating z(t) = sum a_n q^n / n at 90-digit precision and mapping
  through the period lattice, the four fixed points land on **exactly two
  points** (D3 confirmed), from two different matrix families (c = 1 and
  c = 2) that give the same pair:

      Q1 = ( 2i, -3 + 2i)        Q2 = (-2i, -3 - 2i)

  Both lie in E(Q(i)) and are exactly on the curve. Landing on algebraic
  points with coordinates in Z[i] to 90 digits is itself the certificate
  that the normalisation and Manin constant are right — nothing was tuned.
- **Q1 + Q2 = (2, 0) = -2G** with G = (4, 6). So 2P = -2G, hence
  P + G is 2-torsion. E(Q) has trivial torsion and the 2-division
  polynomial 4x^3 - 4x^2 - 4x - 7 is irreducible over Q, so there is **no
  rational 2-torsion** and P is the **unique** rational square root:

      **P = -G = (4, -7)**,   2P = Q1 + Q2 verified exactly.

- **D1 confirmed.** h(P) = 0.2400338318285875365... equals the canonical
  height of the Heegner generator returned by `ellheegner`, so P is a
  generator of E(Q) = Z, not a proper multiple.

## Route 2 — arithmetic (point counting on the reconstructed genus-2 curve)

Independent of everything above. Build the cover from P and nothing else:

- The line through Q1, Q2 is **l1 = y - x + 3**, meeting E again at
  R = -(Q1 + Q2) = (2, -1).
- The tangent at P has slope -3, giving **l2 = y + 3x - 5**, and it meets
  E again at the same R = (2, -1).
- So h = l1 / l2 has div(h) = Q1 + Q2 - 2P exactly, and

      **C2 :  w^2 = (y - x + 3)/(y + 3x - 5)   over   y^2 + y = x^3 - x^2 - x - 2**

  is the double cover with L_-- = O(P), i.e. X0(143)/W143.

Counting C2(F_p) directly, with the three special fibres handled exactly
(h(O) = 1, h(R) = 2, and the leading coefficient at the double pole
h_0(P) = 52, from s(11 - 6x) + s^2 = (x-4)^2(x-2)):

| | |
|---|---|
| primes tested, 3 to 149 (excluding 11, 13) | 32 |
| **a_p(C2) = a_p(143a1) + a_p(11a1)** | **32 / 32** |
| mismatches | **0** |

**Jac(X0(143)/W143) ~ 143a1 x 11a1**, confirmed prime by prime. This
single test simultaneously validates Q1 and Q2, validates P (the cover
was built from P's tangent line), validates that the ghost really is
11a1, and shows the scalar normalisation of h needs no quadratic twist.

Route 1 is a complex-analytic evaluation of CM points through a period
lattice; Route 2 is mod-p counting on an explicit affine curve. They
share no step.

## The descent picture, complete

| object | value |
|---|---|
| E = X0(143)* | **143a1**, y^2+y = x^3-x^2-x-2, rank 1, trivial torsion |
| origin O | common image of the four cusps (a single free AL orbit) |
| deg (L++, L+-, L-+, L--) | **(0, 6, 5, 1)**, two derivations |
| L++ | O_E |
| **L--** | **O_E(P), P = (4,-7), a generator of E(Q)** |
| section of L-- | **g(tau) - 13 g(13 tau)**, g = 11a1 = eta(t)^2 eta(11t)^2 |
| partner section | e8 = g(tau) + 13 g(13 tau), sector (-,+) |
| W11 | acts **freely** (0 fixed points) |
| D_W13 | **Q1 + Q2**, the disc -52 CM points, x = +-2i |
| D_W143 | degree 10, the disc -143 and -572 CM points |
| C2 = X0(143)/W143 | w^2 = (y-x+3)/(y+3x-5), Jac ~ 143a1 x 11a1 |

Every quantity in the four sessions of this arc — the 55 quadrics, the
26/5/4/20 grading, the 12/6/7/11 grading of H^0(2K), the projection
table, the CI-A relation, the CI-B ghost-only 3-space — now traces back
to this one table.

## Status

| item | status |
|---|---|
| modular degree 143a1 = 4 | EXACT |
| Q1, Q2 = (+-2i, -3 +- 2i) | Cert (90-digit landing on exact points) |
| Q1 + Q2 = -2G | EXACT |
| P = (4,-7), unique rational square root | Pr (trivial torsion + irreducible 2-division poly) |
| P generates E(Q) | Cert (height equality with Heegner generator) |
| Jac(C2) ~ 143a1 x 11a1 | Cert (32/32 primes) |

No physics reading is claimed and none follows. The pre-registration for
this arc excludes it.

## Remaining open from the arc

1. **The disc^2 observation** from CI-A (a = 7^2 · 13 · 1957^2) — still
   flagged, not claimed. Cross-level falsification as registered: run the
   same invariant at other levels whose X0(N) has a 1-dimensional (+,+)
   sector and a Galois orbit filling another sector.
2. **CC-08 and CC-09 adoption**, plus the audit of anything in the corpus
   keyed to h(-143) = 7 or to a degree-7 CM minimal polynomial.
3. **`mtft.canonical`** — ship the adapted basis, the adapted
   q-expansions, the four class matrices and this descent table as frozen
   certified data, with the ten bundle-rank tests and the projection table
   as gates.
4. **`mtft.periods`** — still the larger gap. tau in H_13 lives only in
   PARI/GP output.
