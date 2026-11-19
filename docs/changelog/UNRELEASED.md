# Unreleased — after v1.1.0

Everything after v1.0.0 shipped in v1.1.0 ([release notes](CHANGELOG_v110.md)). Since then:

## Wave 7 (2026-11-19): twists and knots

- `mtft.matter.twist_costs` (TC-01): the exact Hessian at the Pati–Salam pattern for any winding weight, the twist table of the spin
  structures of 143a1 (tree Yukawas, V = 0, M_A = 0 and the twelve's flatness go together on all nine pairs), the continuous twist's
  KK spectrum, and `hold_pull`, the pull of a one-loop energy on τ against the γ-flux's well. From records FF-12 and FF-13.
- `mtft.lines.genus` (KN-01): genus characters and the trace-4 law Ψ(A_p) = χ₋₃(p)Ψ(A), which is why the systole of X₀(143) winds
  round the divisor of (η(τ)/η(11τ))¹² and not round the other cusps' units. From record KN-1.
- Ledger entries `docs/tools/TC01_TWIST_COSTS.md` and `docs/tools/KN01_GENUS_SIGNS.md` (kept off the puzzle map, which waits for results); tests `test_twist_costs.py`, `test_lines_genus.py`.
