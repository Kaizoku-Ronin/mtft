# CS-01 — a canonical prescription for the curve-to-torus area ratio? (2026-10-05)

Module `mtft.matter.canonical_shape`; tests `tests/test_cs01_canonical_shape.py`.

**Pre-registration (fixed before evaluation).** Target: V_cb's window A_X/A_E in [18.1, 36.4] (SP-2). A hit needs (i) a definition that never refers to the target, (ii) the
natural normalisation within its family, (iii) its whole convention range inside the window. Every candidate is logged. Chance level: a log-uniform prior on [1, 1000] puts
0.10 of the range in the window, so with four families about 0.35 chance of at least one accidental hit.

| Candidate | A_X/A_E | Verdict |
|---|---|---|
| Pullback of the flat torus through phi (deg 4); Petersson norm | 4 | miss |
| Hyperbolic area 56 pi over the Neron lattice area 3.9709 | 44.30 (22.15 with the factor-2 convention) | straddles: not a hit |
| Arakelov metric against the flat torus, both area 1 | 1 | miss |
| Dimensionless invariants (genus 13, abs chi 24, 2g 26, area/2 pi 28, index 168, cusp widths 1, 11, 13, 143) | 24, 26, 28 in the window | no derivation: numerology |

**Verdict (STRUCTURAL).** No canonical prescription survives. The near-misses are what a window this wide collects from a genus-13 curve. The shape stays a parameter measured
by V_cb. Open: a physical mechanism (flux quantisation of the Kaehler class, a non-perturbative potential) that would turn one of these numbers into an area ratio.
