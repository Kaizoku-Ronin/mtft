# CW-01 — one-loop lifting of the Higgs direction: failed computation, and the conditional protection statement
(computed 2026-09-15; register written 2026-09-18 — the v0.31.4 audit correctly noted that no CW-01 derivation had been
filed; this is the record.  Supersedes the chat-only statement.  CC-28 below.)

## What was computed
On the h = 0.2 and h = 0.3 magnetic meshes (surface.magnetic), the fermion-tower contribution to the one-loop Higgs mass
matrix on the 18-dimensional Higgs space: Q_kl = - sum_{n,m} h^k_{nm} conj(h^l_{nm}) f(M_n^2, M_m^2), with the triple
overlaps h^k_{nm} = int psi_n psi_m conj(psi^H_k) over the family tower, M_n^2 = lambda_n - lambda_0, and f the finite
part of the bubble ((a log a - b log b)/(a - b)); zero-mode x zero-mode pairs excluded (quartic order).

## Why it failed (lessons)
1. The claimed completeness relation sum_{nm} |h^k_{nm}|^2 -> delta_kl is FALSE: the double sum is [delta(x-y)]^2, the
   6D quartic UV divergence.  Measured: 0.10 (25 modes), 0.52 (100), 1.12 (200), 2.38 (400) — linear in the mode count.
2. The "spread" of the form grew with the cutoff (0.10 -> 0.25) and its preferred direction stabilised only where the
   lattice artefact (the unweighted nodal sum of |psi_H|^2, a mesh-density effect) dominates.  Extracting the physical
   direction dependence needs a continuum-limit renormalisation in h at fixed cutoff — not done.
3. A "parameter-free finite" quantity must be checked for convergence before it is called finite.

## The protection statement, corrected (CC-28)
Original (chat, 2026-09-15): "the flux compactification preserves 4D N=1 SUSY, so the Higgs directions are exact moduli
and no loop selects one."  Audit (H-04): the M1 specification declares no supersymmetric parent; an N=1 label alone
does not imply a flat Higgs branch.  Corrected statement: IF the parent is a supersymmetric 6D gauge theory whose
compactification satisfies the Killing-spinor/twist conditions (of which HYM is the D-flatness part), AND the Higgs
branch is F-flat (no mu-term: H_u H_d is forbidden by the extra U(1)s unless axion-dressed, AXG-01), THEN the Higgs
directions are perturbatively protected and the fermion-only loop above is one half of a cancelling supertrace.
None of the antecedents is established for M1 as specified.  Establishing them (or choosing a non-supersymmetric
parent, as C3X does) is the prerequisite for any vacuum-selection computation.
Status: computation FAILED (documented); protection statement CONDITIONAL.
