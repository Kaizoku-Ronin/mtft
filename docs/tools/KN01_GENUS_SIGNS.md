# KN-01 — genus characters and the windings of the systole of X₀(143) round the cusps (2026-11-19)

Module `mtft.lines.genus`; tests `tests/test_lines_genus.py`. Ported from the project record KN-1 (19 November 2026: is the
horoforce a knotting force?), filed under `claude/experiments/` in the MTFT Python project.

The answer KN-1 recorded uses what LN-01 already ships: a closed line's linking number with the trefoil is Rademacher's Ψ
(`mtft.lines.rademacher.Psi`, Ghys), its cyclic word is `mtft.lines.lorenz.lorenz_word`, and its windings round the cusps of X₀(143)
are `mtft.lines.rademacher.windings`. What is new here is the law behind one unforeseen pattern.

For a trace-4 element A of Γ₀(N) (the form Q_A has discriminant 12 = (−3)(−4), one narrow class per genus) and a prime p | N,
Dirichlet composition moves the class of Q_A by a prime ideal of norm p, so

    Ψ(A_p) = χ₋₃(p) Ψ(A),     A_p = [[a, pb], [c/p, d]],

and the winding of a unit ∏ η(δτ)^{r_δ} along the line is (Ψ(A)/24) Σ_δ r_δ χ₋₃(δ). At 143 the signs are (1, −1, 1, −1), so the
systole winds Ψ = ±1 times round the divisor of (η(τ)/η(11τ))¹² and not round the other two units' divisors: 11 ≡ 2 and 13 ≡ 1 (mod 3).

Status: EXACT, two routes (Dedekind sums against genus theory). KN-1 compared them on 1193 elements at the twenty primes p < 200 with
p ≡ ±1 (mod 12) and on the 62 found at 143; the gates rerun eleven primes and 143, and check that no trace-4 element exists where 3 is
not a square. For longer lines Ψ is not a class character, and no such law is claimed.
