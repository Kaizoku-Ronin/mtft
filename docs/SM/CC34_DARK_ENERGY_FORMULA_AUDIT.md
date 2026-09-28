# CC-34 — the dark-energy formula ρ_Λ/M_P⁴ = δ⁻⁶ e^{−2α⁻¹} (Paper 3 §3, Master Paper Thm 33.1, final paper §38) — DIAGNOSTIC

## What the papers state
ρ_Λ/M_P⁴ = δ⁻⁶ exp(−2α⁻¹) with the Feigenbaum constant δ = 4.6692… and α⁻¹ = 137.036…, evaluated as 1.296 × 10⁻¹²³ and reported as matching
observation "within experimental uncertainty" and "without fine-tuning"; corollaries: Λ is exactly constant, dark matter is not a particle.

## Findings (recomputed in `scripts/cosmology/cc34_cosmology_audit.py` and `tests/test_parent_ledger.py`; every number reproducible)
1. **Arithmetic slip.** exp(−2 × 137.036) = exp(−274.072) = 9.37 × 10⁻¹²⁰, not the 1.342 × 10⁻¹¹⁹ used in the papers. The formula's value is
   **9.05 × 10⁻¹²⁴**, not 1.296 × 10⁻¹²³.
2. **Convention dependence of the target.** Planck 2018 gives ρ_Λ = (2.24 meV)⁴ = 2.52 × 10⁻⁴⁷ GeV⁴. Against the ordinary Planck mass this
   is 1.13 × 10⁻¹²³ (the corrected formula is 20% low); against the reduced Planck mass, which the gravitational sector otherwise uses,
   it is 7.15 × 10⁻¹²¹ (the formula is 790 times too small). A derivation must state which mass appears and why.
3. **Fragility.** Two structural choices are unexplained by any action: the power 6 of δ and the factor 2 in the exponent. Replacing 2 by
   2.01 changes the result by a factor 4; replacing δ⁶ by δ⁵ by a factor 4.7; the Hubble tension alone moves the target by 17%.
   Under the E2 rule a one-number coincidence with this sensitivity is DIAGNOSTIC, not CERTIFIED, pending a route from an action.
4. **The corollary "Λ exactly constant" is now under observational pressure.** DESI's BAO results, combined with supernova and CMB data,
   give strong hints that dark energy evolves; the cause (new physics or systematics) is unsettled. This is the programme's most
   falsifiable cosmological statement and should be tracked as such rather than restated.
5. **Dark matter as information curvature** (G_eff/G_N = 1 + η|R_F|) has η unfixed (acknowledged in the compilation as open), and
   "lensing without clumping" must confront the Bullet Cluster and the CMB acoustic peaks before it can be called a prediction.

## Proposed register text
"The formula ρ_Λ/M_P⁴ = δ⁻⁶e^{−2α⁻¹} evaluates to 9.05 × 10⁻¹²⁴ (an arithmetic slip in the papers gave 1.296 × 10⁻¹²³). Its agreement with
observation holds only for the non-reduced Planck mass, at the 20% level, and rests on two undetermined structural choices; it is retagged
DIAGNOSTIC. 'Λ exactly constant' is kept as a falsifiable prediction with DESI as the watch. The information-curvature dark-matter
proposal is retagged as a proposal with one free constant. None of this touches the arithmetic layer (Parts A–E)."

## What would upgrade it
A vacuum energy computed from the parent action on X₀(143) × 143a1: the flux energy (topological plus the ‖ΛF‖² term, which the census
already expresses through the slopes as a function of r) together with the moduli potential. That calculation decides between a constant Λ
and a rolling shape modulus, and it is where the Feigenbaum and α⁻¹ inputs would have to reappear if the formula is to be more than numerology.
