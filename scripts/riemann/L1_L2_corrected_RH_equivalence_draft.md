# The Corrected RH Equivalence: Boundedness of the Normalized Curvature Oscillation

**Status.** Theorem draft for corpus integration. Supersedes Pr 12.4.1, Pr 5.5.3, Pr 0.6.2, and the RH-MTFT box (all false as stated: κ^Λ(y) < 0 unconditionally for all y > 0, verified June 2026). Numerical anchors from the July 2026 sessions are collected in Appendix A. Proof direction L1 is written in full; L2 is written as a complete argument with a gap ledger (§6) listing the steps that remain routine-but-unwritten.

---

## 1. Setup

**Df 0 (Tano weights).** For (a,b) ∈ ℤ², define the two-species weight lattice

  skeleton: w_n^{(a,b)} = Σ_{de=n} d^a e^b Λ(e),  Dirichlet series D_{a,b}(s) = −ζ(s−a)·(ζ′/ζ)(s−b);
  bulk:  w̄_n^{(a,b)} = Σ_{de=n} d^a e^b log e,  Dirichlet series D̄_{a,b}(s) = −ζ(s−a)·ζ′(s−b).

Skeleton series have poles at s = b+ρ for every nontrivial zero ρ (zero-tuned); bulk series are entire off the real points (zero-blind). The (0,0) skeleton edge is the classical identity Σ_{d|n}Λ(d) = log n. These are distinct from Selberg's generalized von Mangoldt functions Λ_k = μ ∗ log^k (Möbius-side convolution). Coalescence: three existing corpus objects are lattice points — the skeleton (2,1) weight is the parent of the present RH diagnostic (Df 1 below); the bulk (2,1) weight is the June-2026 bulk stiffness −ζ(s−2)ζ′(s−1); and the bulk (0,−1) series evaluated on the real axis is the Arithmetica Wick Dirichlet-ensemble partition function Z_D(β) = −ζ(β)ζ′(β+1).

**Df 1 (Skeleton stiffness).** For y > 0 set X = 2πy and

  μ^Λ(y) = Σ_{n≥1} w_n e^{−nX},  w_n = Σ_{dm=n} d² · m Λ(m),

Λ the von Mangoldt function. The Mellin parent is F(s) = ζ(s−2) G(s−1), G(w) = −ζ′(w)/ζ(w), absolutely convergent for Re(s) > 3.

**Pr 1 (Explicit formula).** For every y > 0,

  μ^Λ(y) = 2G(2) X^{−3} + ζ(0) X^{−2} + Σ_ρ m_ρ c_ρ X^{−(ρ+1)} + T(X),

where the sum runs over nontrivial zeros ρ = β+iγ of ζ (m_ρ the multiplicity),

  c_ρ = −Γ(ρ+1) ζ(ρ−1),

and T(X) is the trivial tower: contributions of the double poles of Γ(s) X^{−s} F(s) at s = 1−2k, k ≥ 1, each of the form X^{2k−1}(a_k − b_k ln X); the poles of Γ at s = 0, −2, −4, … are cancelled by ζ(−2) = ζ(−4) = ⋯ = 0. The zero sum converges absolutely and locally uniformly on (0, ∞): by Stirling |Γ(ρ+1)| ≪ γ^{β+1/2} e^{−πγ/2} and |ζ(ρ−1)| ≪ γ^{3/2+ε} in the strip, so Σ_ρ |c_ρ| (1+γ²) < ∞.

*Proof sketch.* Perron/Mellin inversion of Γ(s) X^{−s} F(s) along Re(s) = 4 and contour shift; the vertical decay of Γ(s) justifies the shift. Residues: s = 3 (pole of ζ(s−2), residue 1) gives 2G(2)X^{−3}; s = 2 (pole of G(s−1), residue +1) gives ζ(0)X^{−2}; s = ρ+1 (G(s−1) ~ −m_ρ/(s−1−ρ)) gives m_ρ c_ρ X^{−(ρ+1)}; s = 1−2k the trivial tower. Full tail estimates: gap G1, §6. ∎

*Numerical seal:* direct sum vs. formula at y = 0.02, N = 700, K = 15, four trivial poles: relative residual 2.7×10⁻¹⁸ (Appendix A.1).

**Df 2 (Curvatures and diagnostic).** With D := d/d ln X = d/d ln y:

  M(X) := 2G(2)X^{−3} + ζ(0)X^{−2} + T(X)  (zero-free part),
  Z(X) := Σ_ρ m_ρ c_ρ X^{−(ρ+1)}     (zero part),
  κ^Λ(y) := D² log μ^Λ,  κ_Main(y) := D² log M,
  Δκ(y) := κ^Λ − κ_Main,  𝒟(y) := Δκ(y) · X^{−3/2}.

**Df 3 (Stable decomposition).** Δκ = D[P/Q] with

  P = (DZ)·M − Z·(DM),  DP = (D²Z)·M − Z·(D²M),  Q = M(M+Z),

an exact identity (the cross terms DZ·DM cancel). All subtractions in P, DP, and DP·Q − P·DQ are O(1)-relative; this is simultaneously the numerically stable evaluation and the analytically convenient normal form.

**Df 4 (Admissible multiset; the functional 𝒟[Z]).** Let 𝒵 ⊂ {σ+iτ : 0 < σ < 1, |τ| ≥ τ₀ > 0} be a multiset, closed under conjugation, with counting function N_𝒵(T) ≪ T^A for some A. Define Z(X) = Σ_{ρ∈𝒵} c_ρ X^{−(ρ+1)} with the same coefficient formula c_ρ = −Γ(ρ+1)ζ(ρ−1), and 𝒟[𝒵](y) by Df 2 with this Z. For 𝒵 = the zeta zeros this recovers 𝒟. (The coefficient formula is part of the definition of the functional; c_ρ ≠ 0 for all admissible ρ by Pr 2.)

**Pr 2 (Nonvanishing of coefficients).** For 0 < β < 1: c_ρ ≠ 0 and (2−ρ) ≠ 0.

*Proof.* Γ never vanishes. ζ(ρ−1) has Re(ρ−1) ∈ (−1, 0); ζ(s) ≠ 0 on −1 < Re(s) ≤ 0: nontrivial zeros lie in 0 < Re < 1 and Re = 0 is excluded by the functional equation together with the nonvanishing on Re = 1 (Hadamard–de la Vallée Poussin); trivial zeros are at −2, −4, …; ζ(0) = −1/2. Finally 2−ρ has Re ≥ 1. ∎

---

## 2. The theorem

**Th 1 (Corrected RH equivalence).**

  RH ⟺ limsup_{y→0⁺} |𝒟(y)| < ∞.

More precisely, with Θ := sup{Re ρ}:

(a) If Θ = 1/2 (RH), then |Δκ(y)| ≤ C X^{3/2} for 0 < y ≤ y₀; hence |𝒟| ≤ C.
(b) If Θ > 1/2, then for every β₀ ∈ (1/2, Θ] realized (or approached) by zeros and every ε > 0: limsup_{y→0⁺} |Δκ(y)| X^{−(2−β₀)−ε} = ∞; hence limsup |𝒟(y)| = ∞, with divergence rate ≥ X^{(1/2−β₀)+ε}.

The envelope exponent 3/2 in (a) and the divergence rate β₀ − 1/2 in (b) are both confirmed numerically: measured slope 1.504 vs. 3/2 (true zeros); measured divergence slopes −0.237, −0.394 vs. predicted −0.25, −0.4 for β₀ = 0.75, 0.9, with β₀ = 0.6 in the pre-asymptotic mixture regime (Appendix A.2–A.3, A.7).

---

## 3. Lemma L1 (RH ⟹ bounded): full proof

Assume Θ = 1/2. Write r := Z/M. On (0, y₀]:

1. M(X) = 2G(2)X^{−3}(1 + u(X)) with u(X) = (ζ(0)/2G(2))X + O(X⁴ ln X) from the trivial tower; in particular M > 0 and |u| ≤ 1/2 for y₀ small.
2. |Z(X)| ≤ X^{−5/2} Σ_ρ m_ρ|c_ρ| =: C₁ X^{−5/2} (each |X^{−(ρ+1)}| = X^{−3/2−1} — note Re(ρ+1) = 3/2, so modulus X^{−3/2}; combined with the density bound the sum converges absolutely by Pr 1). Hence |r| ≤ C₂ X^{3/2} → 0.
3. Termwise differentiation is justified by absolute convergence of Σ m_ρ|c_ρ|(1+γ²): |DZ| ≤ C₃X^{−3/2}·sup-weights, |D²Z| ≤ C₄X^{−3/2}, and similarly Dr, D²r = O(X^{3/2}) with constants Σ m_ρ|c_ρ||2−ρ|², Σ m_ρ|c_ρ||2−ρ| absorbed (the shift from exponent −(ρ+1) to 2−ρ comes from dividing by the leading X^{−3}; the correction from u contributes exponents larger by ≥ 1).
4. Δκ = D² log(1+r) = D²r − D(r Dr/(1+r)) and with |r| ≤ 1/2: |Δκ| ≤ |D²r| + 4(|r||D²r| + |Dr|²) ≤ C X^{3/2}.

Therefore |𝒟| = |Δκ|X^{−3/2} ≤ C. ∎

(The same argument gives L1 for any admissible multiset 𝒵 contained in the critical line: 𝒟[𝒵] bounded.)

---

## 4. Lemma L2 (a zero off the line ⟹ unbounded): the abscissa argument

Assume some zero ρ₀ = β₀ + iγ₀ with β₀ > 1/2. (If Θ is not attained, fix any zero with β₀ ∈ (1/2, Θ); the argument is uniform.)

**Step 1 (Normal form).** From Df 3, expanding P/Q = D log(1+r) and using the bounds of L1-steps 1–3 with 3/2 replaced by 2−Θ:

  Δκ(X) = Σ_ρ m_ρ a_ρ X^{2−ρ} + R(X),  a_ρ = c_ρ (2−ρ)² / (2G(2)),

where the leading sum converges absolutely and R collects (i) corrections from u (exponents shifted right by ≥ 1), (ii) quadratic-and-higher terms in r (exponent real parts ≥ 2(2−Θ), i.e. shifted right by ≥ 2−Θ ≥ 1). Hence there is δ₀ ≥ 1 with

  |R(X)| ≤ C X^{(2−Θ)+δ₀}  (0 < X ≤ X₀).  (4.1)

**Step 2 (Laplace transform).** Put v = ln(1/X), f(v) := Δκ(e^{−v}), and for Re(s) > Θ−2 define

  ℱ(s) := ∫_{v₀}^∞ f(v) e^{−sv} dv.

The integral converges there because |f| ≤ C e^{−(2−Θ)v}. Termwise (dominated convergence via Γ-decay),

  ℱ(s) = Σ_ρ m_ρ a_ρ e^{−(2−ρ)(v₀)+…}/(s+(2−ρ)) + ℛ(s)
     = Σ_ρ m_ρ a_ρ E_ρ(s)/(s − (ρ−2)) + ℛ(s),

with E_ρ entire and nonvanishing at ρ−2, and ℛ(s) analytic in Re(s) > Θ−2−δ₀ by (4.1). The zero-sum part is meromorphic in Re(s) > Θ−2−δ₀ with poles exactly at s = ρ−2; these are isolated (zeros are discrete) and the pole at s = ρ₀−2 has residue m_{ρ₀} a_{ρ₀} E_{ρ₀}(ρ₀−2) ≠ 0 by Pr 2. Hence **ℱ has a genuine singularity at s = ρ₀ − 2, with Re(ρ₀−2) = β₀−2.**

**Step 3 (Abscissa contradiction).** Suppose, for contradiction, that |Δκ(X)| ≤ C_ε X^{(2−β₀)+ε} as X → 0⁺ for some ε ∈ (0, δ₀). Then |f(v)| ≤ C_ε e^{−((2−β₀)+ε)v}, so ℱ(s) converges — and is therefore **analytic** — in the half-plane Re(s) > β₀−2−ε, which contains the point ρ₀−2. This contradicts Step 2. Hence for every ε > 0,

  limsup_{X→0⁺} |Δκ(X)| · X^{−(2−β₀)−ε} = ∞.

**Step 4 (Conclusion).** Since 2−β₀ < 3/2,

  limsup |𝒟(y)| = limsup |Δκ| X^{−3/2} ≥ limsup |Δκ| X^{−(2−β₀)−ε} · X^{(1/2−β₀)+ε} = ∞

along the witnessing sequence, for ε < β₀ − 1/2. ∎

**Nb 1 (No conspiracy).** Cancellation among several off-line zeros is impossible: the poles s = ρ−2 sit at distinct isolated points, so no combination of other terms can remove the singularity at ρ₀−2. Verified numerically with two and four simultaneous off-line quadruplets (slope unchanged from the single-quadruplet value) and with an adversarial coefficient sign flip (still divergent) — Appendix A.4.

**Nb 2 (Effectivity and crossovers).** The limsup statement is not effective: when two off-line zeros at (β, γ) and (β′, γ′) with β′ > β, γ′ > γ compete, the larger β′ dominates only below the crossover scale X* given exactly by the coefficient balance ln(1/X*) = ln|a_ρ/a_{ρ′}| / (β′−β), whose leading (Γ-decay) part is π(γ′−γ)/(2(β′−β)); the Stirling polynomial prefactors in a_ρ = c_ρ(2−ρ)²/(2G(2)) shift the crossover materially and must be kept. For (0.75, γ₁) vs (0.9, γ₃): exact balance gives ln(1/X*) = 95.37 nats = 41.42 decades, i.e. X* ≈ 10⁻⁴¹·⁴, versus 49.46 decades from the bare-Γ formula — an 8-decade discrepancy. The measured sliding-window slope crosses the midpoint at ≈ 10⁻⁴², matching the exact balance (Appendix A.5, A.7). Any quantitative version of Th 1 must carry these exponentially small, prefactor-corrected scales.

---

## 5. Corollary: the de Bruijn–Newman dictionary

Let H_t(x) = ∫₀^∞ e^{tu²} Φ(u) cos(xu) du be the Rodgers–Tao flow, H₀(x) = ξ(1/2 + ix/2)/8, and Λ the de Bruijn–Newman constant: H_t has only real zeros ⟺ t ≥ Λ; RH ⟺ Λ ≤ 0; and Λ ≥ 0 (Rodgers–Tao 2018).

**Cor 1.** Let 𝒵(t) := {1/2 + i x_j(t)/2 : H_t(x_j(t)) = 0} interpreted as an admissible multiset (real zeros x_j contribute on-line points; a non-real zero x = a+ib of H_t contributes the off-line point 1/2 + b/2 + ia/2 and its partners). Then, granting L1/L2 at the multiset level (Df 4) and the density bounds for zeros of H_t:

  Λ = inf{ t : limsup_{y→0⁺} |𝒟[𝒵(t)](y)| < ∞ },

and in particular RH ⟺ limsup |𝒟[𝒵(0)]| < ∞ ⟺ Λ = 0 (using Λ ≥ 0). The corrected equivalence is thereby anchored to the named, studied constant: the curvature diagnostic is a boundedness coordinate for the de Bruijn–Newman flow.

*Numerical anchors:* H₀ = ξ(1/2+ix/2)/8 sealed at 10⁻¹³; zeros x_j(0) = 2γ_j to 4×10⁻¹¹; drift dx_j/dt matches the repulsion ODE ẋ_j = Σ_{k≠j} 2/(x_j−x_k) (K = 40 + density tail) at ratios 1.000, 1.000, 1.001; amplitude along the flow (1/A)dA/dt|₀ = +0.340 = (dγ₁/dt)(d ln A/dγ₁) = (−0.290)(−1.17), with the analytic single-zero amplitude matching dense RMS to ≤ 0.5% (Appendix A.6).

**Nb 3.** For t < Λ the non-real zeros of H_t occur at heights far beyond numerical reach (Lehmer-pair pinching); the synthetic off-line multisets of §4 are the honest laboratory for that regime, legitimate because 𝒟[·] depends only on the zero multiset.

---

## 6. Gap ledger (to full rigor)

G1. Pr 1 contour shift: write the T → ∞ tail estimate (vertical Γ-decay makes this standard) and the convergence/asymptotic status of the trivial tower T(X) on (0, X₀].
G2. L1 step 3 / L2 step 1: write the explicit constants in the remainder bound (4.1), including the u-corrections and the quadratic-in-r terms, and verify δ₀ ≥ 1 uniformly for Θ ∈ (1/2, 1).
G3. L2 step 2: the interchange (termwise Laplace transform) with the stated dominated-convergence bound; meromorphy of the zero-sum transform on Re(s) > Θ−2−δ₀ away from the pole set.
G4. Cor 1: (i) admissibility (density bound) for zeros of H_t — available from Polymath 15 literature; (ii) the multiset extension of Pr 1 is definitional (Df 4), but the statement "non-real H_t zero ⟹ off-line point of 𝒵(t) with Re ≠ 1/2" needs the b ≠ 0 ⟹ |Re − 1/2| = |b|/2 > 0 bookkeeping written out.
G5. Multiplicities: L2 assumes the residue m_ρ a_ρ ≠ 0; m_ρ ≥ 1 integer and a_ρ ≠ 0, so this is automatic — record it.

None of G1–G5 appears deep; G2 is the longest write-out.

## 7. Corpus edits entailed

- Pr 12.4.1, Pr 5.5.3, Pr 0.6.2, RH-MTFT box: replace "κ_∞(y) ≥ 0 ⟺ RH" by Th 1. κ^Λ(y) is strictly negative for all y > 0 unconditionally; the RH-sensitive object is the boundedness of the normalized oscillation, not the sign of the curvature.
- §12.7(c) (modular approach): unchanged as a proof strategy for Th 1(a)-type bounds via SL(2,ℤ) constraints on W(q).
- Add Cor 1 as the bridge proposition to the de Bruijn–Newman constant; cite Rodgers–Tao (2018) for Λ ≥ 0 and Polymath 15 for effective bounds.

---

## Appendix A: numerical anchors (July 2026 sessions; mpmath, dps 30–60; scripts in repo)

A.1 Explicit-formula seal: y* = 0.02, relative residual 2.68×10⁻¹⁸; zero sum resolved at −2.36×10⁻⁷ against that floor.
A.2 True zeros: slope of log|Δκ| vs log y = 1.504 (pred. 3/2) over y ∈ [10⁻⁷, 10⁻²]; max|𝒟| = 3.6×10⁻⁶.
A.3 Single off-line quadruplets: slopes of log|𝒟| = −0.067 / −0.237 / −0.394 for β₀ = 0.6 / 0.75 / 0.9 (pred. −0.1 / −0.25 / −0.4). The β₀ = 0.6 window [10⁻⁷,10⁻²] is not yet asymptotic for so slow a divergence (mixture with the on-line background); β₀ = 0.75, 0.9 land within 0.013 and 0.007.
A.4 Conspiracy tests: 2× and 4× quadruplets at β = 0.75: slope unchanged from single (no cancellation); adversarial sign flip: still divergent.
A.5 Crossover: (0.75, γ₁) + (0.9, γ₃): corrected-estimator slopes −0.25027 on [10⁻³⁰, 10⁻²⁰] and −0.40052 on [10⁻⁶⁴, 10⁻⁵⁶] (pred. −0.25, −0.40); sliding-window midpoint at ≈ 10⁻⁴², matching the exact coefficient-balance crossover 10⁻⁴¹·⁴ (Nb 2). Computed with the stable normal form Δκ = D[P/Q] (Df 3); the naive two-κ subtraction loses ~|log₁₀ X| digits and fails below 10⁻²⁰.
A.6 dBN anchor: H₀ vs ξ relative error ≤ 1.4×10⁻¹²; x_j(0)/2 − γ_j ≤ 4×10⁻¹¹ (j = 1); drift ratios 1.000/1.000/1.001; A(t) = 2.364, 2.445, 2.487, 2.529, 2.620 (×10⁻⁶) at t = −0.15, −0.05, 0, +0.05, +0.15.
A.7 Estimator note: bin-RMS slope fits must drop terminal bins containing few samples; a grid endpoint on a bin boundary creates a one-sample bin whose phase-lottery value has maximal regression leverage (observed bias up to 0.035 in slope). The artifact can masquerade as systematic when the window stride is near-resonant with the oscillation (stride 6 decades vs γ₃: 6γ₃ln10/2π = 54.99 cycles, within 1.4% of integer). All A.2–A.5 values above use the corrected estimator (bins with ≥10 samples); pre-correction values (1.479, −0.258, −0.413, −0.2534, −0.4068…−0.4353) are superseded.
