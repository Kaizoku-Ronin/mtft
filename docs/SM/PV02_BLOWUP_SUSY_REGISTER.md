# PV-02 — supersymmetry read arithmetically, and the blow-up (B3) with the supersymmetric twist (29 September 2026)

Companion to PV-01.  EXACT unless marked; the finite search is exact within its boxes.

## 1. Supersymmetry in the existing code, named
- Every index in the ledger is a Witten index of Q = ∂̄ + ∂̄*: on X0(143) the de Rham version has 2 bosonic and 26 fermionic zero modes, index χ = −24.
- Stability (Donaldson–Uhlenbeck–Yau) is D-flatness; the slope of a block is its Fayet–Iliopoulos term; `higgs_slope_mass` is the FI mass.
- For one block (a, b) with field strengths x, y along the two factors: energy ∝ rms², charge ∝ gm², so rms ≥ gm is the BPS bound; equality
  (all five means equal) is the anti-self-dual, zero-slope, supersymmetric point for opposite-sign blocks; (x − y)² = 2(rms² − gm²) is the D-term.
- `vacuum_energy`: E_split − E_HYM = 4π² Σ nᵢ(μᵢ − μ̄)²/r, the slope variance — the D-term energy of a split background
  (S1: (41r² − 14r + 89)/6r, S2: 5(13r² − 10r + 13)/6r, both positive: supersymmetry broken at every shape).  `tests/test_pv02.py`.
- Arithmetic: the Möbius function is (−1)^F (Spector 1990); Σ_{d|n} μ(d) = [n = 1] is level-by-level cancellation; PNT ⟺ Σμ(n)/n = 0;
  RH ⟺ Mertens M(x) = O(x^{1/2+ε}).  The alphabet: Σ wₙ n^{−s} = −ζ(s)ζ′(s+1), and μ ∗ w = (log n)/n.  Continuation to s = −1:
  −ζ(−1)ζ′(0) = −(log 2π)/24 (exact identity; meaning open).

## 2. Eleven, thirteen, twelve
Twin-prime levels 6k ∓ 1: gm² = 36k² − 1, am² = 36k², rms² = 36k² + 1; genus of X0(pq) is k(3k + 1) − 1; g − 1 = 6k forces (3k + 1)(k − 2) = 0.
So 12 = (11 + 13)/2 = g − 1 = deg K^{1/2} of X0(143) — the spin twist (12, 0) — is unique to 143 among twin-prime levels (test).

## 3. What a supersymmetric parent implies (exact)
- 8D supersymmetry (16 supercharges) has only vector multiplets: bulk matter is adjoint; PV-01's U(8) class with Λ²(8) + 8·8 is non-supersymmetric.
- On a curved surface supersymmetry needs the topological twist R = O: adjoint blocks then have net chirality −c₁(L)·K = −24b on the product
  (χ(X0(143)) × torus flux; `adjoint_net_chirality((a, b), (0, 0))`), never 3.  Families need curves (B2) or a blow-up (B3).

## 4. `mtft.surface.blowup` (B3)
- Lattice U ⊕ ⟨−8⟩ ⊕ ⟨−1⟩ (modular degree 4, exceptional curve), classes (a, b, c, e), K̃ = 24f₁ + E, K̃² = −1, χ(O) = 0, not spin.
- SUSY twist: net(L) = e − 24b, linear: for adjoint blocks a coboundary φᵢ − φⱼ with φ(d) = e − 24b.  Consequence: with hypercharge from the
  stack U(1)s the exact Standard Model needs n_A = n_B + 1 singlet stacks — minimal U(8) = U(3) × U(2) × U(1)²_A × U(1)_B, with
  φ_c − φ_L = 3, φ_B − φ_c = 3, φ_A1 + φ_A2 − 2φ_c = 3, and at least six vector-like doublet pairs in every solution.
- Kähler class A_X f₁ + A_E f₂ − εE ample iff 0 < ε < min(A_X, A_E); slope μ = aA_E + bA_X + eε.
- The block (−1, 0, 3): index 3, slope zero at ε = A_E/3 — a D-flat three-family block.
- `d_flat_sm_search(n_box=8)`: 116 D-flat (all slopes equal, ample) exact-SM assignments in the box; fewest vector-like extras 31
  (non-collinear, exceptional fluxes up to 29); the single-line-bundle family v = (1, 0, −3) with charges (0, 1, −1, 2, −3) is traceless
  (an SU(8) flux), D-flat on ε = A_E/3 for every r > 1/3, with 6 u^c, 6 e^c, 9 doublet pairs and 15 sterile singlets.  The floor of 9 extras
  (φ_A1 ∈ {0, 3}) is not reached in the box.
- Open: h¹ and h⁰(K̃ ⊗ L) on the blow-up (the Yukawa rule under the twist), the choice of blow-up point (a CM point × e₀), global anomalies,
  the Green–Schwarz structure of the extra U(1)s, and every scale.

Visual companion: `viz/arithmetic_susy_picture.html` (five interactive pictures; published from this session).
