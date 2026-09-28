# The parent-action ledger — what any parent of MTFT on X₀(143) × 143a1 must satisfy, and the tools to enforce it (28 September 2026)

This is the "Sudoku" for the parent action: every constraint proved so far, the theorem it comes from, and the package gate that checks it.
A candidate parent is a finite data set — dimension, gauge group, fermion representations, spinor twist, extra terms — and the gates decide.

## A. Constraints already proved (exact)

| # | Requirement | Theorem / source | Gate (package) | Status of the three-stack 8D gauge theory |
|---|---|---|---|---|
| 1 | An odd number of chiral families | chirality parity (notes 12): adjoint matter on the spin surface has even net chirality, zero for the spin twist | `adjoint_net_chirality`, `chirality_parity_theorem` | fails: S1/S2 families are vector-like pairs |
| 2 | A Yukawa between chiral families | one-ray lemma (notes 13a): split polystable vacua have none | `three_stack_trilemma`, `block_census` | fails for split vacua |
| 3 | u^c and d^c in different summands | hypercharge lemma (notes 13b) | (to build: `hypercharge_assignment`) | forces ≥ 2 singlet classes |
| 4 | SU(3)×SU(2)×U(1)_Y unbroken at M_KK | colour slope condition + trilemma (notes 11) | `colour_slope_condition`, `three_stack_trilemma` | fails for all four VI.3 solutions |
| 5 | No destabilising twisted colour line | S1 theorem (notes 10), rank-2 lemma (notes 13c) | `s1_iterated_extension_instability`, `rank2_net_index_lemma` | fails for rank ≤ 2 recombinations |
| 6 | Tachyon-free background | census + Hodge index | `block_census`, `census_scan` | fails at every r for split backgrounds |
| 7 | Light Higgs at leading order | slope formula 2πμ/(A_X A_E) (VI.3 remark) | `higgs_slope_mass` | needs μ(L_H) = 0 |
| 8 | Strong CP after the m_u lift | CP structure (notes 9): torus CP-conserving, curve breaks CP by the flux choice | `cm_point_cp_structure`, `torus_cp_test` | open: lift must be phase-aligned or an axion exists |
| 9 | 6D tensor/anomaly integrality (for 6D parents) | CC-33 | `tensor_integrality`, `gaugino_p2_integrality` | none of M1's assignments survives |
| 10 | Kinetic normalisation of hypercharge | III.7 (CC-26) | `cc26_normalisation` | done |

## B. The three candidate classes that survive A.1 (each needs its own gates)

**B1. A chiral 8D parent** — fermions in complex (non-self-conjugate) representations, so one block per pair carries the families.
Then VI.1's index ab is the family number, A.2–A.7 apply with rank ≥ 3 summands required, and the 8D gauge anomaly (a ten-form I₁₀,
with the quintic invariant tr F⁵ present for U(N), N ≥ 5, and mixed gauge–gravitational terms) must cancel or be Green–Schwarz-compensated.
Tools to build: `anomaly8d.I10(field_content)` (Chern-root expansion, exact rationals, in the style of T25/CC-33), a field-content
enumerator with the parity and index gates, and the rank ≥ 3 vacuum search.

**B2. Matter on curves** — families localised on a curve C ⊂ S (X × {e₀}, the graph of the modular parametrisation, or {x} × E), where
chirality is a degree on a curve and can be 3 (the Chapter III–IV picture, embedded). Needs a 6D defect in the parent and the 6D
anomaly/tensor integrality of A.9 on the defect. Tools: `surface.curves` (classes and self-intersections in NS = U ⊕ ⟨−2δ⟩, δ the modular
degree), a localised-matter index calculator, the bulk–defect Yukawa rule (the analogue of VI.2 with one curve leg).

**B3. A non-spin surface** — blow up X₀(143) × 143a1 at an arithmetic point (a CM point × e₀): the exceptional curve makes NS odd, so
c₁(F)² and hence the net chirality can be odd with adjoint matter, keeping the 8D gauge theory. Tools: `surface.blowup` (the lattice
U ⊕ ⟨−2δ⟩ ⊕ ⟨−1⟩, spin^c structures, the parity gate on it), then the census and slope gates on the new lattice.

## C. Cosmology, in the same discipline
- CC-34 (`CC34_DARK_ENERGY_FORMULA_AUDIT.md`): the dark-energy formula is DIAGNOSTIC; "Λ constant" is the falsification watch against DESI.
- Tool: `cosmology.vacuum_energy(r)` — the flux energy of a background (topological part plus the ‖ΛF‖² term from the census slopes) and
  the resulting potential for the shape modulus r and the volume; decides constant Λ versus a rolling modulus once a parent is fixed.
- Dark-matter candidates to register as open items: the singlet stack's 14 Wilson-line moduli; a W₁₃-protected lightest odd state.

## D. Repository upgrades that serve all of the above
1. `docs/SM/PARENT_ACTION_REQUIREMENTS.md` — this ledger, kept append-only; `tests/test_parent_ledger.py` checks that every listed gate
   exists and returns its recorded verdict (the ledger cannot drift from the code).
2. `mtft.parent` — a candidate-parent data class plus `gates(candidate)` returning EXACT verdicts per requirement; a `mtft parent` CLI
   subcommand and Legend entries.
3. `anomaly8d`, `surface.curves`, `surface.blowup`, `cosmology` as above.
4. An "external checks" page: the exact lemmas an arithmetic geometer can verify independently (gonality from the 𝔽₄ count, the Petri
   gate on the W₁₃ quotient, the class polynomial of the CM points, the closed-form theta factors) — the validation contract.
5. Viewer inputs: the density atlas, the Lorentzian NS picture and the particle gallery as exported grids with a small loader.

Release shape: v0.33.1 (eight commits: atlas, review corrections, census, S1 wall, S1 theorem, trilemma, parity, this ledger with CC-34) now; v0.34.0
when the `parent` gates and at least one of B1–B3's toolchains land.
