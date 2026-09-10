# mtft v0.27.4 — two-prime old sector and the Hermitian-part selection rule; docs fixes (2026-09-10)

`surface.oldsector` (EXACT, sympy): local blocks B_p, Petersson Gram G_p, adjoints, Hermitian/skew
parts; the level-N0·pq tensor sector with U_p, U_q, H_add, H_prod, H_pq = Herm(U_p U_q);
connected part by partial traces in the G-orthonormal frame; Schmidt rank; and
`hermitian_hecke_selection`: within H_α = A_p⊗A_q + α R_p⊗R_q the Hermitian parts of the Hecke
algebra Q[U_p,U_q] form a 4-dimensional space containing exactly one interaction direction,
α = 1 (H_pq).  The interaction term is the product of the skew parts of the U's, which vanish at
good primes.  Reproduces Astra's composite-correspondence study for 11a1 at (13, 17): connected
spectrum {−108,−12,12,108}, raw charpoly (x²−92x−1584)(x²+100x−1584), Schmidt ranks 0/1/2.
Status of the rule: a stated arithmetic-native principle, canonical and basis-independent, not
derived from the curve; it is the candidate for the disc lab's "independently motivated rule".
Docs: intertwiner header no longer says the periods frame has no integer map; provenance string
records the Hodge orientation minus sign.  Base: live 0.27.3 tarball.  Tests: 2.  Pin four-way.
