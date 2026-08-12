# X₀(143): Period Matrix in the Manin Basis and the Genus-13 Riemann Matrix — v6

**Session:** 2026-08-11 · **Engine:** PARI/GP 2.15.4 (container) · **Precision:** realprecision 50, tolerance 10⁻⁴⁰ · **Tally: 19 PASS / 1 FAIL**

**Artifacts:** `mtft_period_matrix_manin_v6.gp` (self-contained, Windows: `cmd /c gp -q mtft_period_matrix_manin_v6.gp`) · `period_matrix_manin_v6.txt` (full 50-digit data) · `x0143_period_data_v6.json` (machine-readable) · `run_v6.log`

This closes **Paper 32, Open Problem 3** in full normalized form: not only the period matrix over a Manin basis, but the Riemann matrix τ over a certified symplectic basis, with the exact integer change of basis recorded so every cycle is reconstructible.

## 1. Objects built

**Ω (13 × 29)** — rows are the thirteen holomorphic differentials [f₁ | f₂ embeddings 1–4 | f₃ embeddings 1–6 | 11a(q) | 11a(q¹³)]; columns are the 29 unimodular `mspathgens` paths, which in weight 2 form a *free* ℚ-basis of H₁(X₀(143), cusps; ℚ) because the single Manin relation has all group-ring degrees zero. Convention (Df): entry = ∫_path f dz, no 2πi.

**K (29 × 26), Ω_cusp = Ω·K** — integer kernel of the rank-3 boundary map to the four cusps {1, 11, 13, 143}.

**Q, E (26 × 26)** — `mspetersson` returns the intersection product on the 29-dimensional symbol space as an integral antisymmetric matrix of rank exactly 26. Pulled back along the dual lift L = K(KᵀK)⁻¹, this gives Q = LᵀPL on cuspidal homology: integral, antisymmetric, det 1. The intersection form is E = Q⁻¹.

**S (26 × 26) ∈ GL₂₆(ℤ)** — from a from-scratch skew-Smith reduction over ℤ (minimal-pivot search, Euclid divisibility loop, clear, sign-normalize, recurse), giving SᵀES = J = [[0, I₁₃], [−I₁₃, 0]] exactly.

**τ (13 × 13)** — Ω_cusp·S = (A | B), τ = A⁻¹B. Symmetric with Im τ positive definite: an explicit point of Siegel upper half-space 𝓗₁₃. τ₁₁ = −0.54625994325611195984 + 1.13059702944511923787 i.

**Auxiliary:** T₂, T₃ on homology from Merel cosets + `mspathlog` (no q-expansions), in both the Manin and symplectic bases; the complete P¹(ℤ/143) table — all 168 Manin symbols with (c:d), a unimodular path, and 29-generator coordinates.

## 2. Certificate ledger

| # | Statement | Status | Margin |
|---|---|---|---|
| C1 | Weight-2 Manin relation trivial; 29 paths free basis | Cert (exact) | integer |
| C2 | Boundary rank 3; cuspidal kernel dim 26 | Cert (exact) | integer |
| C3 | All 13 rows: Manin-route L-value = `lfun` L-value | Cert (numeric, E2) | 2.4×10⁻⁵⁸ |
| C4 | diag(a_p)·Ω = Ω·T_p, p = 2, 3, Merel-coset T_p | Cert (numeric, E2) | 4.4×10⁻⁵⁶, 5.4×10⁻⁵⁶ |
| C5 | charpoly(T₂\|cusp) = x²(x+2)⁴g₄²h₆² | Cert (exact, third route) | integer |
| C5b | charpoly(T₃\|cusp) = (x+1)⁶g₄₃²h₆₃² | Cert (exact) | integer |
| C6 | per₁₁ vs Paper 33 archive; \|Re λ₁\| = ½; \|Im λ₁\| anchor | Cert (numeric) | 4.7×10⁻⁵¹ |
| C7 | Gram rank 13 | Cert (numeric) | eig ∈ [7.3×10⁻³, 17.9] |
| C8 | Direct `mfsymboleval` = Ω × pathlog coordinates | Cert (numeric) | 0 |
| C9a | f₁ period lattice = Λ(143a1)/(2πi), Manin constant 1 | Cert (numeric, E2 — AGM) | 1.3×10⁻⁵⁷ |
| C9b | Sign discriminator: exactly one candidate ∈ Λ | **FAIL — honest negative** | discriminator degenerate |
| R8 | `mspetersson` radical dim exactly 3 = 29 − 26 | Cert (exact) | integer |
| R1 | E integral, antisymmetric, det 1 | Cert (exact) | integer |
| **R2** | **Ω_cusp·Q·Ω_cuspᵀ = 0 (first Riemann bilinear relation)** | **Cert (numeric, E2)** | **2.7×10⁻⁵⁶** |
| R3 | SᵀES = J exactly, \|det S\| = 1 | Cert (exact) | integer |
| R4 | τ symmetric | Cert (numeric) | 7.0×10⁻⁵⁷ |
| R5 | Im τ positive definite | Cert (numeric) | eig ∈ [0.0928480858, 4.1143928656] |
| R6 | Hecke integral + equivariant in symplectic basis | Cert (numeric) | 6.8×10⁻⁵⁶, 6.5×10⁻⁵⁶ |
| R6b | charpoly invariant under symplectic conjugation | Cert (exact) | integer |
| R7 | Genus-1 replay at level 11: j(τ₁₁) = j(11a1) = −122023936/161051 | Cert (numeric, E2) | 2.7×10⁻⁵⁵ |

R1–R8 were pre-registered in the script and printed before any of them was computed.

**The load-bearing E2 pairs share no computational steps.** R2 is the strongest: a combinatorially computed *integer* pairing annihilates a transcendental period matrix obtained by a wholly disjoint analytic route. C4 couples Merel-coset homology to embedded eigenvalues; C3 couples modular-symbol evaluation to independent `lfun` machinery; C9a couples the modular-symbol period lattice to the AGM lattice of 143a1; R7 validates every stage of the symplectic pipeline end-to-end where the answer is known in closed form.

**Orientation** came out E = +Q⁻¹, fixed by positivity of Im τ and corroborated by the level-11 replay. Recorded as Df.

## 3. Open flag (unchanged from v5): the [∞, 2/77] sign

[∞, 1/11] agrees with the Paper 33 v2 archive to 10⁻⁵¹, but v5/v6 (PARI 2.15.4) evaluates [∞, 2/77] with the opposite overall sign to the archive (PARI 2.17): λ₁ = −½ − 1.0232745926964612…i here versus +½ + 1.0233i there. Since ∫ f dz is path-independent, one is wrong; magnitudes are unaffected either way.

The pre-registered adjudicator C9b is **filed as an honest negative on the discriminator design**: the data shows per₁₁ = b₁ and per₂₇₇ = −b₁ + b₂, so *both* sign candidates for the difference are lattice points and the test has a structural kernel. Adjudication requires a route outside PARI's `mfsymboleval`: Sage/Cremona modular symbols, Magma, or direct contour quadrature of the q-expansion. **Assigned to the Kimi queue. No CC entry is filed pending that result**; if the v6 sign survives, the Paper 33 correction is sign-only.

## 4. Engine notes (for the corpus)

Three PARI gotchas cost debugging cycles and are logged: `M[,perm]` will not accept a vector of column indices — use `matrix(n,n,i,j,S[i,perm[j]])`; `qfjacobi` rejects 1×1 input, so genus-1 replays need a special case; and `mfsymboleval` on a degree-1 orbit returns a degree-0 t_POL that must pass through the `myval` coercion or it silently propagates as a t_SER and poisons downstream arithmetic. Also: `mspetersson(M)` with no symbol arguments returns the entire intersection matrix in ~1 ms — no need to loop over basis pairs.

## 5. Downstream

The symplectic basis makes three things immediately computable: the Hecke action on τ itself; the isogeny decomposition of J₀(143) into factors matching the [1, 4, 6] Galois-orbit dimensions, read off the invariant subspaces of T₂ˢ in the symplectic basis; and theta-function/Schottky work at the τ point. The 168-symbol table plus S also permits the Atkin–Lehner involutions to be written as explicit integer symplectic matrices, which is the natural route into the (−,−) sector question.
