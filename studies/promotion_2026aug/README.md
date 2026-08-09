# promotion_2026aug — provenance studies for v0.14.0

Two study chains plus two standalone studies: 9 studies, 76 gates
total, all green. Each script writes its JSON ledger alongside; the
ledgers here are the session records (byte-preserved, including the
falsified pre-registrations they carry).

## X0(143) skeleton chain (4 studies, 33 gates)

1. `x0143_graph_uncertainty.py` (7 gates) — the ancestry commutator
   [D,L] = A(d_v − d_u) evaluated on the house graph: where the
   noncommutativity lives. Self-loop edge identified.
2. `x0143_ribbon_embedding.py` (10 gates) — R-cyclic orderings make the
   dual skeleton a ribbon graph; surface embedding, genus, spanning
   trees, the real involution's fixed triangles.
3. `x0143_hecke_particles.py` (9 gates) — Manin symbols on the same
   P^1(Z/143) flag set; charpoly(T_2|H_1) = x^2 (x+2)^4 g4^2 h6^2;
   the four particle blocks. Promoted to `mtft.hecke`.
4. `eisenstein_congruences.py` (7 gates) — the congruence moduli of the
   four blocks (1, 5^4, 7^2, 12^2), bad-prime U-operators, the
   Sturm certificate at bound 28. Promoted to `mtft.eisenstein`.

## Tano ensemble chain (3 studies, 28 gates)

5. `w2_susceptibility.py` (9 gates) — <w^2> closed form, chi_w,
   Cov(log n, w) = zeta''(beta+1), cold constants. Promoted to
   `mtft.moments`.
6. `w3_cumulants.py` (8 gates) — the triple-Euler U engine,
   Amari-Chentsov tensor, skewness census. Promoted to `mtft.moments`.
7. `curvature_tano_manifold.py` (11 gates) — Brioschi curvature of the
   (beta, lambda) manifold: cancellation theorem, Hagedorn slope,
   summit, flat temperature, cold dive 6/5. Promoted to
   `mtft.curvature`. Contains the falsified 3-atom pre-registration
   (G7b) — preserved per protocol.

## Standalone studies (2 studies, 15 gates)

8. `triangular_layers.py` (7 gates) — honest negative: the triangular
   Faulhaber layers are a re-encoding of the Bernoulli tail
   (a_{p,3} = −4 a_{p,2} exact; a layer-4 closed form; higher layers
   flag only false positives; layer-4 hit rate inside the decoy band).
   No module; the finding is registered in the legend.
9. `m7_graph_channel.py` (8 gates) — the canonical zero-parameter
   graph-side coupling V (triangle cusp-width potential compressed to
   homology through the harmonic embedding): NOT in the Hecke algebra
   ([V,T_p] nonzero in all 676 entries for p = 2, 3, 5), genuinely
   flavor-changing (all 12 off-block amplitudes nonzero), exactly
   rational and G-self-adjoint — the mixing is a real orthogonal
   rotation, CP-EVEN; iota*-even under the star involution; the
   degree-potential null control is trivial exactly as required
   (V = 3I); the graph-distance second potential reproduces the
   channel (non-Hecke, flavor-changing, still real). Study-only
   landing (M6/M7 physics arc); no module.

## Auditor notes

- F1 (stale record): the curvature ledger's gate G7c records
  K_12346 = 1.3549364825; the correct limit is 1.3549368866023 (the
  ledger's top-level key and `mtft.curvature` agree).
- F2 (resolved in the second drop): `eisenstein_congruences.py` and
  `triangular_layers.py` landed; both re-ran ledgers identical to the
  uploaded records (7/7 gates each) in the auditor's environment.
- `m7_graph_channel.py` landed with a one-line portability patch,
  disclosed here: `from mtftpkg import hecke as H` →
  `from mtft import hecke as H` (v0.14.0 exposes `mtft.hecke`
  publicly; `mtftpkg` was the author's local alias). No other byte
  touched.
- `m7_graph_channel_ledger.json`: the ledger uploaded with the script
  predated the script's P7a/P7b revision (its failed degree-function
  P7 is exactly the null control P7a in the current script). The
  landed ledger was regenerated from the author's unmodified script —
  8/8 gates pass, 247.5 s; the pre-revision record is preserved
  verbatim in the v0.14.0 audit report. Gate P7b (graph-distance
  potential, uncertified at drop time) is now certified on two
  independent routes: the author's engine and the auditor's exact
  sympy route agree to all printed digits (off-block norm 16.772862,
  676/676 commutator entries vs T_2, G-self-adjoint).
- Gate-total correction: an earlier draft of this README counted the
  Tano chain as 30 gates / 56 total; the ledgers hold 28 / 54 for the
  two chains (76 across all nine studies). Counts above are from the
  ledgers themselves.
