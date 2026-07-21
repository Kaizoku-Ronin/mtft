# Anatomy of the Degree-7 Jacobian Conjecture Counterexample
### (Alpöge / Claude Fable 5, July 20 2026 — reverse-engineered July 20 2026)

**The map.** F(x,y,z) = (F₁,F₂,F₃), F₁ = (1+xy)³z + y²(1+xy)(4+3xy), F₂ = y + 3x(1+xy)²z + 3xy²(4+3xy), F₃ = 2x − 3x²y − x³z. det DF ≡ −2 (sympy exact, PARI exact, numeric spot-checks). Component degrees (7, 6, 4).

**Verified facts** (all exact arithmetic; scripts jc_check.py, jc_reverse.py, jc_cubic.py, jc_verify_mechanism.py):

1. **Degree 3.** Every tested fiber has 3 points with multiplicity; étale ⟹ reduced. Geometric degree d = 3.
2. **ℤ/2-equivariance.** F∘ι = ι′∘F for ι(x,y,z) = (−x,−y,z), ι′(X,Y,Z) = (X,−Y,−Z). F₁ even; F₂, F₃ odd.
3. **The fiber cubic.** F₁ and F₃ are linear in z; eliminating z then y yields the tautological identity P(x; F₁,F₂,F₃) ≡ 0 for the **depressed** cubic
   P(T; A,B,C) = p₃T³ + p₁T + p₀,
   p₃ = 27A²C² − 18ABC + B³C + 16A − B²,  p₁ = 4 − 3BC,  p₀ = −2C.
   The three fiber x-coordinates are the roots of P; each root extends to exactly one fiber point (y, z rational in x and the target). No T² term ⟹ fiber x-coordinates sum to 0, matching the equivariance. Sanity: at (−1/4,0,0), P ∝ T(T−1)(T+1) — the observed collision fiber {0, ±1}.
4. **Escape wall = {p₃ = 0}.** There the cubic degenerates directly to linear (T² already absent): the two ι-paired sheets **escape to infinity together**; fiber count drops 3 → 1. Verified at three wall points (all fiber = 1). This is Jelonek's non-properness hypersurface A(F).
5. **Missed curve = {p₃ = 0} ∩ {p₁ = 0}.** Eliminating: 144A² − 24AB² + B⁴ = (12A − B²)², so the locus is the smooth curve {A = B²/12, C = 4/(3B)}. There P ≡ p₀ ≠ 0: **no roots, empty fiber.** Verified empty at (1/3, 2, 2/3) and (1/12, 1, 4/3). Hence Image(F) = ℂ³ ∖ {curve} — open, complement of codimension 2. F is étale, non-proper, **non-surjective**, degree 3.
6. **Monodromy S₃, one organizing hypersurface.** disc_T(P) = −4·S²·p₃ with S = 27AC² − 9BC + 8. Not a square ⟹ Galois group of the generic fiber = S₃. The quadratic resolvent ℚ(A,B,C)(√(−p₃)) is branched along the *same* hypersurface {p₃=0} that is the escape wall and contains the missed curve. S∘F factors into two irreducibles (deg 5 × deg 8), consistent with the distinguished-sheet splitting inside the S₃-closure.
7. **Invariant spine.** The x- and z-axes swap with doubling: F(a,0,0) = (0,0,2a), F(0,0,c) = (c,0,0); F² doubles along each axis. The collision target (−1/4,0,0) sits on this spine.

**The mechanism, in one paragraph.** Over W := ℂ³ ∖ {p₃=0}, F restricts to a *proper* 3-sheeted étale covering space, classified by a transitive S₃ representation of π₁(W). The classical intuition for JC — "ℂⁿ is simply connected, so étale self-maps should be injective" — silently assumes étale maps are covers; they are only covers where proper. This construction realizes a genuinely nontrivial finite cover of a hypersurface complement and extends it to a polynomial endomorphism of all of ℂ³ by arranging the only compatible escape pattern for a depressed cubic: the two conjugate sheets leave together across the wall, and on the codimension-2 stratum where even the last sheet leaves, the image simply omits a curve. Non-injectivity is then forced globally: monodromy permutes the sheets, and any S₃-orbit over a ι′-fixed target produces the observed rational triple collision. The naive prototype — the tautological root-cover (r,p) ↦ (p, −r³−pr) of the depressed cubic family — has non-constant Jacobian −(3r²+p); the achievement of the third variable z is exactly the room needed to re-coordinate the family until the Jacobian is the constant −2. Why n=3 and not n=2 is now a sharp question: whether a 2-variable polynomial family of cubics (or any d ≥ 2 family) admits the constant-Jacobian tuning. That is the surviving open case.

**Consequences recorded.** JC false for all n ≥ 3 (pad by identity); n = 2 open; Dixmier conjecture false for A₃ (via DCₙ ⟹ JCₙ); cubic-homogeneous counterexamples exist in some dimension (BCW/Drużkowski); Ax–Grothendieck untouched.

**Flag-and-dismiss.** Bézout product of component degrees 7·6·4 = 168 = [PSL₂(ℤ) : Γ₀(143)]: coincidence class AG-D5; no mechanism; not corpus material.
