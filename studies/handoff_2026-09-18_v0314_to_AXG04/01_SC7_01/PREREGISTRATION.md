# SC7-01 — spin-circle extension of the supplied MTFT model

Recorded 2026-09-16 before running this experiment's code. This is a bounded mathematical and effective-field-theory test, not a construction of an ultraviolet-complete parent.

## Fixed inputs

- Supplied MTFT 0.31.4 archive, SHA-256 2e6fcd0fde3b1601e3cf58aaec338790aa760034e548ecded3086bfa2a110db5.
- Compact curve X = X0(143), genus 13, and S0 = O(6(0)+6(1/11)), S0^2 = K_X.
- P = S(S0), the unit-circle bundle, with fiber orientation chosen so Euler class = c1(S0) = 12.
- Same independent line-bundle stack data as M1: ranks (3,2,1,1,1), degrees (0,-3,3,3,0).
- Geometric Q8 lift from the preceding audit; this is a symmetry of the underlying spin geometry, not necessarily of the flux vacuum.

## Model being tested

The 7D gauge sector lives on M^(1,3) x P. Use the standard partial topological twist of 7D supersymmetric Yang–Mills, with smooth commuting stack backgrounds and no boundary, defect, singular source, or new projection. In an anti-Hermitian convention its internal BPS equations are F_A - phi wedge phi = 0, d_A phi = 0, d_A^* phi = 0. For commuting diagonal fields the gauge connection must be flat. Massless modes are computed by the corresponding flat complex-connection cohomology, as in the cited 7D Higgs-bundle framework.

This convention differs from simply transplanting the original 6D Weyl index. A failure here excludes this specified extension, not every possible 7D parent.

## Gates and controls

1. **Topology:** independently recover integral homology using the circle-bundle Gysin sequence and abelianization of its fundamental-group presentation. Track how all M1 bundle degrees reduce modulo 12.
2. **BPS connection:** distinguish the naive pulled-back curved connection from a flat connection on the same torsion bundle. Derive the fiber holonomy of the flat extension, with explicit orientation convention.
3. **Zero modes:** calculate the exact twisted cellular cohomology using Fox derivatives of the Seifert presentation. Check it independently using the fiber-cohomology spectral sequence. A nontrivial rank-one fiber holonomy predicts acyclicity; this is an analytic expectation, not a blind numerical discovery.
4. **Controls:** the trivial character must give Betti numbers (1,26,26,1). A generic nontrivial base character with trivial fiber holonomy must give (0,24,24,0). Scan every torsion residue and all actual charged M1 sectors. Recovering the model requires the appropriate three chiral families, not merely a total of three states.
5. **Symmetry:** retain exact Q8 square-map compatibility; do not infer that W11 survives the odd-degree M1 flux assignment. Check what ordinary Fourier sections on P retain from S0, while distinguishing them from 7D physical zero modes. A full Hecke lift is not presumed.
6. **Gravity:** for a separately assumed two-derivative unwarped 7D Einstein action, compute the connection metric's scalar curvature directly and compare it with the submersion formula. Test the uniform radius direction in 4D Einstein frame with nonnegative bulk vacuum energy and ordinary nonnegative flux energy. No added negative-tension sources, quantum corrections, or fitted stabilizing terms.

## Stopping and reporting rules

- Exact algebra uses rational/Gaussian-rational arithmetic, Smith normal form, and polynomial identities where possible.
- Acyclicity of the physical charged sectors is a decisive failure of this minimal matter construction. Do not add defects or tune boundary conditions after seeing it and then report the original gate as passing.
- Symbolic identities under the declared model assumptions are EXACT. No numerical mesh result, absolute physical mass, or supersymmetric gravitational embedding is certified.
- Preserve the supplied source. Deliver a separate report, runnable code, result ledger, input hashes, and reproducible bundle.
