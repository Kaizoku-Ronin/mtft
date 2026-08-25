# mtft v0.20.0 — the origami tier

## New: `mtft.origami` (subpackage)

Dimers, origami / t-embeddings, and the **observable insertion calculus**,
built on P. Galashin, *Amplituhedra and Origami, I: Tree Level*
(arXiv:2410.09574v2).  Three reusable layers plus two certified instances.

### `mtft.origami.dimer`
`DimerGraph` — any weighted planar bipartite graph in a disk.  Almost-perfect
matchings, boundary measurement (1.1), full Plücker vector,
**Grassmann–Plücker certification** (proves the vector really is a point of
Gr(k,n) rather than bookkeeping), cyclic-symmetry test, coarse-graining by
sufficient statistics, and the **`ensemble_conservation` gate** — aggregating
micro-states into equal-statistic classes must preserve Z exactly.  That gate
exists because an external audit caught a dropped matching; it is now
structurally impossible to repeat that error silently.

### `mtft.origami.insertion`
The one formula behind the Tano moments, the Dirichlet log-n metric, dimer
edge occupations, the Ising models, and differentiation under the integral:

    d^k/dλ^k log Z = κ_k(A)

`cumulants`, `fisher_metric`, `cubic_tensor` (Amari), `brioschi_curvature`,
`cumulant_curvature` (independent E2 route, documented sign convention), and
the **`path_independence` gate**: dψ = Σ⟨N_e⟩dθ_e is an exact differential, so
line integrals depend only on endpoints and every loop integrates to zero —
the statistical-mechanical form of a conservative field.  All in
multiplicative coordinates with D = X d/dX, which keeps everything rational.

### `mtft.origami.perfect`
Galashin (9.2)–(9.3): the Θ involution, `t_coefficients`, winding and
sign-flip validity, the perfect-system solver `solve_perfect_branches`,
`cyclic_matrix` / `orbit_structure` for branch orbits, and
`equivariant_kasteleyn_factor` for character-block factorization.

### `mtft.origami.instances`
- **(2,4)** Galashin Example 2.18 with an explicit Kenyon–Smirnov t-embedding,
  Lemma 1.10 verification, and the exact channel identities
  **S(1,3) = pr, S(2,4) = qs** — the two terms of the Gr(2,4) Plücker exchange
  relation are individually the two Mandelstam channels.
- **(3,6)** the C₃-symmetric hexagonal prism (Figure 14 / Example 9.21):
  91 APMs, Z = 280, top cell, plus `PRISM_C` (certified rational
  representative) and `PRISM_LAMBDA0` (the exact C₃-fixed perfect branch).

### `mtft.origami.gates`
Fourteen gates, runnable end to end:
`from mtft.origami.gates import run_all; run_all()`

Certified results now carried by the package:
- section A is the full trinomial ⇒ **K ≡ 1/4 exactly**, coefficients symbolic;
- section B closed forms, both curvature routes agreeing with residual 0:
  `det g = XY(XY+cX+cY+c)/(c+X+Y+XY)³`,
  `K = −(c−1)(c−X²)(c−Y²)/(4(XY+cX+cY+c)²)`;
- **flatness by symmetry**: g even and the cubic tensor odd under the
  complement involution, so K = 0 at the fixed point follows from symmetry
  alone, with ρ = tanh J and g* = ¼[[1,ρ],[ρ,1]];
- strict bounds −(c−1)/4 < K < (c−1)/(4c) from positive-coefficient brackets;
- the Fisher quadric z₀₀z₁₁ = √c·z₁₀z₀₁ in S³(2), the c = 1 case being the
  Clifford torus (the geometric meaning of the flat/product anchor);
- **Δ₀₂₄(w) = (2w+3)(w²+3)** — trivial character × Eisenstein norm, product
  forced by C₃, with w = 2 the unique positive balance point giving
  49 = 7_trivial × 7_Eisenstein (w = 0 gives the ramified 9);
- four perfect branches organized **1 + 3** under C₃, verified as a group
  action, with the fixed branch known exactly.

## New: `mtft.hardy_ramanujan`
An orthodox end-to-end benchmark that claims nothing new and exercises the
whole stack: modular form → partition function → cumulant geometry → saddle →
integer combinatorics.  `psi_direct` vs `psi_modular` is an E2 pair agreeing to
~1e-40 (η modularity as an **identity**, not an asymptotic); `saddle_partition`
inverts on the exact ψ and beats the closed form ~3×; `prefactor_identity`
certifies symbolically that 1/(4√3·n) is (modular half-log) × (Gaussian
determinant).

## Honest negatives carried forward
- Positroid boundaries are **nonsingular** in Fisher curvature and lie at
  finite Fisher distance; only singularity *detection* fails.  The naive
  "Fisher sees positroid degeneration" conjecture is dead at (2,4).
- The Q(L)/C_Eisenstein conjecture from v0.19 remains falsified for f3 on
  support grounds.

## Version
Three-way guard: `pyproject.toml`, `__init__.py`, `CITATION.cff` all 0.20.0.
