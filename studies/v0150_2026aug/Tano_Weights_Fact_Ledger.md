# Tano Weight Fact Ledger
**w_n = Σ_{d|n} (log d)/d — every proved fact, certified computation, and standing connection in the MTFT corpus**

Compiled August 11, 2026. Status tags follow house convention: **Df** definition · **Pp** proposition · **Pr** proved · **Th** theorem · **Cert** machine certificate · **Measured** reproducible computation · **Heur** interpretation · **Conj/Open** program target · **EXACT / CERTIFIED(tol) / PHENO** exactness classes · **[Ext]** external artifact · **SUPERSEDED** retracted from corpus.

**Sourcing discipline.** Everything below traces to a mounted artifact this session (papers 1–37, Arithmetica Generale, the two dictionaries, the Three Ensembles draft, the coalescence/audit reports, `arithmetic_wick.py`, the mtft skill's v0.11.4 certified-anchor table) — except items tagged **[Mem]**, which come from the v0.12.x–v0.14.0 session record; the transcript mount is empty this session, so **[Mem] items must be re-verified against the repo/changelogs before citation**. Two corrections are newly flagged by this compilation (§3.4, §12.1) with in-session machine verification.

---

## 1. Definition and elementary structure

**1.1 (Df — the weight).** w_n := Σ_{d|n} (log d)/d, the normalized divisor-log operator. w_1 = 0. First ten values (corrected — see §12.1):

| n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| w_n | 0 | 0.346574 | 0.366204 | 0.693147 | 0.321888 | **1.011404** | 0.277987 | **0.953077** | **0.610340** | **0.898720** |

Bold entries correct the Chapter 5 table (§12.1). w_4 = log 2 exactly.

**1.2 (Pr — convolution identity).** n·w_n = f(n) = (σ ∗ Λ)(n), where f(n) = Σ_{k|n} (n/k) log k, σ the divisor-sum, Λ von Mangoldt. [Master Paper Th; Three Ensembles Df 1.1; audit-verified at n = 4, 6; re-verified this session at n = 6, 8, 9, 10 to 10⁻¹².]

**1.3 (Pr — spectral-derivative representation).** w_n = σ′₋₁(n) = ∂_α σ_α(n)|_{α=−1}: the weight is the α-derivative of the divisor-power family at the harmonic point. [Paper 1 Prop 9.2.]

**1.4 (Pr — prime values and UV regulation).** w_p = (log p)/p → 0; the prime skeleton {(log p)/p at primes, 0 else} is the EXTRACT of the full weight sequence. Σ_{p≤x} w_p ~ log x (Mertens), giving logarithmic UV control: high modes are suppressed by the arithmetic itself, a regulator absent from the Wilson action. [Paper 24 §4.1; `arithmetic_wick.py`.]

---

## 2. Generating functions and the CC-02 identity chain

**2.1 (Th — Lambert series).** W(q) = Σ_{n≥1} w_n qⁿ = Σ_{d≥1} (log d)/d · q^d/(1−q^d). This is the exact q-series feeding the modular holonomy potential. [Paper 1 Th 10.1.]

**2.2 (Th — master Dirichlet identity).** F(s) := Σ f(n) n^{−s} = −ζ(s−1)·ζ′(s), Re s > 2. Proof: series of σ is ζ(s)ζ(s−1), of Λ is −ζ′/ζ; convolve. [Master Paper Th 6.1; AG Pr 0.5.2; audit ✓.]

**2.3 (Pr/Cert — the shift chain; CC-02).** One meromorphic object, three shifts:
- Σ w_n n^{−s} = **−ζ(s)·ζ′(s+1)** = F(s+1)  (the weight series itself)
- Σ n·w_n n^{−s} = −ζ(s−1)·ζ′(s) = F(s)
- Σ n²·w_n n^{−s} = −ζ(s−2)·ζ′(s−1) = F(s−1)  (the stiffness weighting)

The E = n (Laplace) vs E = log n (Dirichlet) conventions differ by exactly one unit shift ("shift discipline," Three Ensembles Pr 3.1). The coalescence report states the correction explicitly: "the series of w_n itself is −ζ(s)ζ′(s+1)." Known corpus offender: the Mathematical Dictionary line "F(s) = −ζ(s−1)ζ′(s) = Σ n²w_n n^{−s}" — that RHS is F(s−1), not F(s). CC-02 filing with three independent certification routes: **[Mem]**. **Fresh route this session:** Σ_{n≤2×10⁵} w_n n^{−3} agrees with −ζ(3)ζ′(4) to 1.2×10⁻¹¹ (truncation-tail scale). **Cert.**

**2.4 (Measured — weighted theta modularity no-go).** `weighted_theta_cusp_fit`: fitted cusp constant A = −0.24963 against predicted −1/4, C ≈ ζ′(0); filed as the modularity no-go for the weighted theta of w_n. [v0.7.1 ship record, audit §J.]

---

## 3. Mean values and Tauberian laws

**3.1 (Th — vacuum torque flux).** With S(N) = Σ_{n≤N} n·w_n = Σ_{n≤N} f(n):
T∞ := lim S(N)/N² = **−ζ′(2)/2 = 0.46877412715792188…**
Tauberian from the simple pole of F(s) at s = 2 with residue −ζ′(2). [Master Paper Th 7.2. In-session partial sum at N = 2×10⁵: 0.4687467 ✓.]

**3.2 (Cert/EXACT — Cesàro mean).** lim (1/N) Σ_{n≤N} w_n = **−ζ′(2) = 0.9375482543158438** (same constant as the marked-primon-gas cold α, EXACT closed form; audit item "Cesàro average of w_n" ✓). [In-session at N = 2×10⁵: 0.9373063, residual at the O((log N)²/N) scale ✓.]

**3.3 (Cert — Paper 33 Tauberian law).** S(y) ~ c/y² with c → −ζ′(2)/(2π)² is a genuine Tauberian law. [Audit ✓.]

**3.4 (FLAG — normalization tension, this compilation).** Two corpus statements conflate 3.1 and 3.2: AG §0.5 writes "T∞ = lim (1/N) Σ w_n = −ζ′(2)/2" (correct value, wrong defining sum — the (1/N)-normalized limit is −ζ′(2)); the Mathematical Dictionary writes "torque normalization T∞ = −ζ′(2) ≈ 0.9375" (correct Cesàro value attached to the torque symbol). Correct pairing is exactly 3.1 + 3.2. **Candidate CC entry** if not already inside CC-02's scope — check the CC ledger before filing.

---

## 4. The two-species Tano weight lattice

**4.1 (Df — lattice; companion Df 0).** For (a,b) ∈ ℤ²:
- *skeleton*: w_n^{(a,b)} = Σ_{de=n} d^a e^b Λ(e), series D_{a,b}(s) = −ζ(s−a)·(ζ′/ζ)(s−b) — poles at s = b+ρ for every nontrivial zero: **zero-tuned**;
- *bulk*: w̄_n^{(a,b)} = Σ_{de=n} d^a e^b log e, series D̄_{a,b}(s) = −ζ(s−a)·ζ′(s−b) — entire off the real points: **zero-blind**.

w_n is the **bulk (0,−1) point**. [Three Ensembles Df 1.2.]

**4.2 (Pr — coalescence).** Three previously separate corpus objects are lattice points: skeleton (2,1) generates the Th-1 RH diagnostic; bulk (2,1) is the bulk stiffness series −ζ(s−2)ζ′(s−1); bulk (0,−1) on the real axis is the Dirichlet partition function Z_D. [Three Ensembles Pr 1.3.]

**4.3 (Working rule).** New ensembles are built by choosing (a,b); **check the pole structure before interpreting asymptotics** (skeleton species encode RH in their poles; bulk species do not). [Skill recipes.]

---

## 5. Laplace ensemble — stiffness and mass gap

**5.1 (Df).** T(y) = Σ_{n≥2} w_n e^{−2πyn}; gauge-filtered stiffness μ_N(y) = Σ n² w_n e^{−2πyn}(1 − cos 2πn/N) (dictionary form carries (N/2)·min_m over m). Skeleton stiffness μ^Λ uses the (2,1) skeleton species.

**5.2 (Pr — filtered moment identity, EXACT).** μ_N(y) = (1/4π²)·[T″(y) − Re T″(y − i/N)]. The gauge filter is a second difference of one analytic object in the imaginary direction; the mass gap is the difference between free-energy curvature at real and at ℤ_N-center-twisted coupling. Machine: rel_diff = 2.5×10⁻¹⁶ at (y_c, N=3). [Three Ensembles Pr 2.2; anchor ✓ on 0.11.4.]

**5.3 (Th — unconditional positivity).** μ_N(y) > 0 for all N ≥ 2 and all y > 0 — a nonperturbative arithmetic lower bound valid at **all** couplings. [Paper 24 Th 4.4.]

**5.4 (Pr — gauge filter / sieve identities).** Σ_{m=1}^{N−1}(1 − cos 2πnm/N) = N·1_{N∤n} (character orthogonality); even-N universality: min_m at m = N/2 collapses every even SU(N) to the odd-integer sieve 1−(−1)ⁿ; SU(p) suppresses exactly the p-th Euler factor of ζ; the mode average ⟨μ_N⟩(N−1)/N → μ_∞ reconstructs the full sum — Eratosthenes in disguise. [Dictionary; Paper 24 §4.3.]

**5.5 (Measured — head domination).** In the explicit-formula zero sum Σ_ρ c_ρ X^{−(ρ+1)} with c_ρ = −Γ(ρ+1)ζ(ρ−1), Stirling gives |c_ρ| ≪ γ^{β+1/2}e^{−πγ/2}. At X = 2πy_c the first zero pair carries **0.999946** of the 12-zero magnitude sum: the Laplace kernel is effectively a single-zero probe. [Three Ensembles Measured 2.4.]

**5.6 (Cert — canonical depths).** y_c = **0.18174** (`CriticalDepths.y_conf`, canon-test enforced). Drifting paper/module values to canonicalize: 0.1821304, 0.1812561, ln π/2π = 0.1821894, and `arithmetic_wick.py`'s 0.1812. Skeleton zeros y_s1 = 0.1236, y_s2 = 0.2106 (module constants).

**5.7 (Th — the corrected RH equivalence).** RH ⟺ limsup_{y→0⁺} |𝒟(y)| < ∞, where 𝒟 = (κ^Λ − κ_Main)·(2πy)^{−3/2} is the normalized curvature oscillation of the skeleton stiffness. Proof direction L1 complete; L2 complete with routine gap ledger; adversarial off-line conspiracy scans (slope = 1/2 − max β, crossover verified to y = 10⁻⁹⁵) back the mechanism. [Companion Th 1, two-engine validated; Three Ensembles Th 2.3.]

**5.8 (SUPERSEDED).** The four κ ≥ 0 propositions — AG Pr 0.6.2, Pr 5.5.3, Pr 12.4.1, the RH-MTFT box, and the Dictionary's "RH ⇔ κ(y) ≥ 0" line — are **retracted**, superseded by 5.7. Do not cite the κ ≥ 0 form.

---

## 6. Dirichlet ensemble — the arithmetic Wick bridge and Speiser

**6.1 (Pr — closed form, EXACT).** Z_D(β) = Σ_{n≥2} w_n n^{−β} = **−ζ(β)·ζ′(β+1)**, β > 1 — the arithmetic Wick rotation: Laplace ↔ Dirichlet via Mellin, e^{−2πyn} ↔ n^{−β}, E = n ↔ E = log n, y ↔ β. [Three Ensembles Pr 3.1; `arithmetic_wick.py`; audit "exact ✓".]

**6.2 (Pr — curvature split).** g_D(β) = Var_β(log n) = ∂²_β log ζ(β) + ∂²_β log(−ζ′(β+1)): a von Mangoldt piece Σ Λ(n)(log n)n^{−β} > 0 plus the logarithmic curvature of ζ′. Anchor: g_D(3) = 0.3351038786441419, ζ′ share **48.5889 %** (EXACT). [Pr 3.2; Kimi-verified independently at β = 2.5, 3, 4.]

**6.3 (Pr/Cert — Hadamard over ζ′ zeros).** (s−1)²ζ′(s) is entire of order 1 (genus ≤ 1), whence ∂²_s log(−ζ′(s)) = 2/(s−1)² − Σ_{ρ′}(s−ρ′)^{−2}, sum over all ζ′ zeros. Pole coefficient **exactly 2**, confirmed Hadamard-free. CERTIFIED |residual| < 10⁻⁵ on s ∈ [3,10]; catastrophically conditioned outside — a usage rule, not a doubt. Do not quote the per-point s = 3 residual (4.0×10⁻⁹) as a global bound. [Pr 3.3; audit §I.]

**6.4 (Pr — Speiser intrinsic; census).** Speiser 1935: RH ⟺ ζ′ ≠ 0 in 0 < Re s < 1/2 — so Z_D's own second factor carries an RH-equivalent analyticity statement. Census: 19 ζ′ zeros to γ < 100, all Re ρ′ > 1/2 (min 0.7806), refined to |ζ′(ρ′)| < 10⁻²⁹; negative-axis zeros one per interval (−2n−2, −2n), first at −2.7172628292, via the exact H₂ functional-equation solver (no negative-axis ζ evaluation). Berndt: N′(T) = (T/2π)log(T/4πe) + O(log T) — the ζ′ ensemble is sparser than ζ's by (T/2π)log 2. [Pr 3.4; v0.7.1 API; skill anchors.]

---

## 7. Critical ensemble — Li coefficients

**7.1 (Df + criterion).** λ_n = (1/(n−1)!)(d/ds)ⁿ[s^{n−1} log ξ(s)]|_{s=1}. Li 1997: RH ⟺ λ_n ≥ 0 **for all n**. Bombieri–Lagarias caveat rides every report: any finite prefix of positivity is logically empty. [Three Ensembles Df 4.1.]

**7.2 (Pr/Cert — values).** λ₁ = 1 + γ/2 − ½ log 4π = **0.02309570896612103** (EXACT; three independent methods, errors 7×10⁻³² / 1.4×10⁻³⁷; fourth FFT route to the float64 floor). Certified: λ₂ = 0.09234573522804667, λ₃ = 0.20763892055432480, λ₄ = 0.36879047949224164, λ₅ = 0.57554271446117745, λ₁₂ = 3.2632553206246199; Keiper (1992) agreement to all comparable digits. [Audit §K.]

**7.3 (Pr — monotone lower bounds).** λ_n = Σ_ρ [1 − (1−1/ρ)ⁿ] paired; on the line the pair term is 4 sin²(n·arctan(1/2t)) ≥ 0, so truncations are monotone lower bounds. [Pr 4.3.]

**7.4 (Measured — density domination).** The λ tail is density-dominated: zero ordinates hug the Riemann–von Mangoldt density (direct ratio 0.99999); S(T) fluctuations contribute ~0.1 % at γ ≈ 236; bulk fluctuations cancel to a uniform 1.5×10⁻⁴. The density-first comparison program (kernel inequalities against the smooth density, then S(T) control) is the stated Conj/Open target. [Measured 4.5; §5 of the draft.]

---

## 8. The Euler shift deformation (Paper 20 — all Pr within analytic number theory)

**8.1** f^(e)(n) = Σ_{k|n}(n/k) log(ek) = f(n) + σ(n); equivalently f^(e) = σ ∗ (Λ + 1) (the deformed von Mangoldt Λ+1: constant background field of weight 1 on top of the prime spikes). Weights decompose **w^(e)_n = w_n + σ₋₁(n)** — prime-sensitive logarithmic layer plus prime-blind harmonic layer. Electron mode: w^(e)_1 = 1.

**8.2** Dirichlet series: Σ w^(e)_n n^{−s} = ζ(s)·[ζ(s+1) − ζ′(s+1)] (consistent with the CC-02 chain plus Σ σ₋₁(n)n^{−s} = ζ(s)ζ(s+1)).

**8.3** Lambert series of the correction: Σ σ₋₁(n)qⁿ = −log η(τ) + πiτ/12 — the shift connects the weights to the Dedekind eta function and the partition function.

**8.4** Spectral reading: the shift is the covariant derivative e·∂_α(e^α σ_α)|_{α=−1} along the Eisenstein α-line; continuous flow α ∈ [−1, 0].

**8.5** Endpoint identity: w^(0)_n = ½ d(n) log n at α = 0.

**8.6** RH link: growth of σ₋₁(n) = σ(n)/n is governed by Robin's inequality, making RH equivalent to a boundedness condition on the shifted weights (Euler-shift vacuum stability).

---

## 9. Weights in the geometry / information layer

**9.1 (Df).** Arithmetic Hosotani potential V_Hos(θ) = −Σ_n w_n log|1 − U_n e^{iθ}|; the mass gap is the smallest eigenvalue of CURVE(V_Hos) at the equilibrium θ₀. [AG Df 0.6.2.]

**9.2 (Df).** Torsion field T_N(q, β) = Σ_{n=1}^N w_n · ∂² log Z_n/∂β∂q — CURVE of each log-assembly, reassembled with the weights. [AG Df 0.6.1.]

**9.3 (Df/Pp).** Fisher–Rao metrics of the two Gibbs families: g_L(y) = Var_y(n), g_D(β) = Var_β(log n); their critical points are the pipeline's fixed points in each picture. Vacuum components (AG Pr 6.1.3, 6.4.3): g_ββ = Var[log w_n], κ_∞(1) = H₂(w) − T∞² + T∞ ≥ 0 automatic at y = 1. **Caution:** these AG formulas sit next to the §3.4 normalization slip — re-derive the constants when reusing. H(w) has no elementary closed form (AG Nb).

**9.4 (Measured/PHENO).** Shannon information I(N,y) = −Σ p_n ln p_n of the weight distribution discriminates number classes. [Dictionary.]

**9.5 (Cert/EXACT).** Marked primon gas on the same arithmetic: cold-gas α = −ζ′(2) = 0.9375482543158438, B = e^{−ζ″(2)} = 0.13679384954115162 (closed forms); KMS modular-flow check; spectral edge softness. [Skill anchors, 0.11.4.]

**9.6 (Pp).** F(s) gives the lattice sum a meromorphic continuation; RG flow via Mellin; ζ zeros enter the stiffness as oscillatory corrections through the explicit formula — the analytic structure absent from the Wilson action. [Paper 24 §4.2.]

**9.7 (Connection).** The falsifiability engine's 23 pre-registered zero-parameter predictions route through w_n and X₀(143) (`mtft.falsify.honest_report()`).

---

## 10. Dynamical-units results touching the weights

**10.1 (Cert — du02 free-level obstruction).** On the shared 26-dimensional harmonic stage of X₀(143), the graph clock acts as **zero** while the Hecke clock runs all 13 lines doubled — so no free exchange rate χ_H/χ_g exists on homology; an internal rate requires an interaction lifting the harmonic degeneracy. du03 (the dispersion relation of the 26 modes against the 13 Hecke lines) is the open frontier. [Skill, artifact-backed.]

**10.2 [Mem].** Two-clock ledger invariants (du01), parity selection rule certified, systole ℓ_sys = 2 arccosh(2), Lamb-type effective Hamiltonian — session record within v0.11.4; re-verify against `references/particle-box.md` and the du study artifacts.

---

## 11. Recent session layer (v0.12.x–v0.14.0) — all [Mem]; re-verify against repo before citing

**11.1** CC-02 filing: the §2.3 chain certified by **three independent routes** (this compilation adds an in-session fourth at s = 3).

**11.2** `curvature.py` (v0.14.0): Tano-manifold **sign-changing Gaussian curvature** certified; zero-crossing at β₀ = 8.8565170425; theorem that the fourth cumulant κ₄ contributes **exactly zero** to the Gaussian curvature of any exponential family.

**11.3** `moments.py`, `hecke.py`, `eisenstein.py` (v0.14.0): weight-moment machinery plus Hecke/Eisenstein companions, 29 tests across the four modules.

**11.4** `weil.py` (staged, **gated pending Kimi audit** — four tasks): Gabor-compressed Weil Hermitian form; W1 pre-registered study complete. **Honest negative WI-N1 on file:** no critical-line certificate is achievable for the X₀(143) newform L-functions via Weil compression (degree-2 Montgomery–Vaughan length wall). Follow-ons W2 (degree-1 Dirichlet characters mod 11/13/143, tying to `eisenstein.py`) and W3 (Weil form of ξ′ against the Speiser census) are queued.

**11.5** CC-03: `finite_atom_curvature` made adaptive in precision.

---

## 12. Corrections ledger for the weights (append-only)

**12.1 (NEW FLAG, this compilation — Chapter 5 weight table).** MTFT_Chapter5's worked example and Table 1 misstate four values: it gives w_6 ≈ 1.060, w_8 ≈ 1.040, w_9 ≈ 0.732, w_10 ≈ 1.036. Correct values: **w_6 = 1.011404, w_8 = 0.953077, w_9 = 0.610340, w_10 = 0.898720.** Verified this session by two independent routes agreeing to 10⁻¹²: (i) direct divisor sieve to N = 2×10⁵; (ii) closed forms via f = σ∗Λ — w_6 = (4log2+3log3)/6, w_8 = 11log2/8, w_9 = 5log3/9, w_10 = (6log2+3log5)/10. Note the Chapter 5 prose lists the correct summands for w_6 but sums them wrong. Prime and prime-power entries n = 2,3,4,5,7 in the table are correct. **Candidate CC-04** (or next free number) — file per append-only protocol.

**12.2 (FLAG — §3.4).** T∞ / Cesàro normalization tension across AG §0.5 and the Mathematical Dictionary; correct statements are 3.1 and 3.2.

**12.3 (Recorded — §2.3).** Dictionary "Σ n²w_n n^{−s}" shift error; corrected chain in §2.3 (CC-02 scope).

**12.4 (Recorded — §5.8).** κ ≥ 0 family superseded by the boundedness equivalence Th 5.7.

**12.5 (Context).** CC-01 (Paper 26 §8 Tano Mass Formula) is retracted — a Hecke-eigenvalue claim, not a weight statement, listed because it consumed weight-adjacent structure.

---

## 13. Verification stamp (this session)

Direct sieve (N = 2×10⁵) vs f(n)/n closed forms: agreement 10⁻¹² at n = 6, 8, 9, 10. Cesàro partial 0.9373063 → −ζ′(2) = 0.9375482543158438; torque partial 0.4687467 → −ζ′(2)/2 = 0.46877412715792188; both residuals at their predicted convergence scales. CC-02 identity at s = 3: Σ w_n n^{−3} vs −ζ(3)ζ′(4) agree to 1.2×10⁻¹¹ (truncation-tail scale). mpmath dps = 30.
