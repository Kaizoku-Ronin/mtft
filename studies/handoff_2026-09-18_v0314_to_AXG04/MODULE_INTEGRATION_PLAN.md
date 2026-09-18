# Proposed module and tool integration

These are **proposed APIs**, not functions already present in v0.31.4. The source archive and studies are authoritative about existing behavior. Claude should adapt names to the repository's current conventions after inspecting the target branch; do not assume its working tree still equals this supplied candidate.

## 1. Preserve the experimental boundaries

The baseline already has `mtft.surface.arithspin.al_lift_on_S0`, `smflux.abelian_anomaly_polynomial`, `smflux.mixed_nonabelian_anomalies`, `smflux.majorana_obstruction`, and `hym.normalise_yukawa`. Extend these narrowly where appropriate. The new parent/action tools deserve an explicitly experimental namespace (for example `mtft.research.parents`) until their scope and API stabilize. The example namespace is a proposal; it is absent from the supplied source.

Archive these studies before refactoring them. Most audit scripts execute calculations at module scope and write JSON beside themselves; importing them directly into the public package would create unwanted work and file writes. Extract pure functions, keep the studies as regression references, and provide explicit CLI/report entry points.

## 2. Result and model records

Every public research result should carry enough information to reconstruct what was checked:

```text
study_id, model_id, schema_version
mathematical_status: EXACT | DIAGNOSTIC | OPEN
model_gate: PASS | FAIL | NOT_TESTED | CONDITIONAL
assumptions, conventions, source_hashes
inputs, outputs, witnesses, limitations
supersedes (only where applicable)
```

This separates successful computation from a failed proposed spectrum. Use exact integer/rational JSON encodings, or strings for algebraic expressions with a documented symbol table. Do not serialize exact rationals as approximate floats. Matrix records need row/column bases, action side, homology/cohomology convention, orientation, kinetic metric and transformation law when relevant.

Model records must distinguish:

- global gauge group and any quotient, representations and multiplicities;
- physical 6D chirality from presentation as left-handed 4D fields;
- flux/cocharacter lattice, periods and charge normalization;
- internal manifold, compact/cusped choice, bundle and actual kinetic/mode operator;
- parent field multiplicity from internal zero-mode multiplicity;
- supplied couplings/potentials from computed consequences;
- local cancellation from integral and global consistency checks.

`C3X` is a candidate identifier. It is not a replacement definition for M1 or a new release version.

## 3. Implementation map

| ID | Proposed home / functions | Source to extract | Required acceptance gates |
|---|---|---|---|
| INT-01 | `surface.arithspin`: retain `al_lift_on_S0`; add `theta_compatible_al_lift` and `check_theta_square_map` | Initial audit `eta_transform`, `symmetric_square`, `group_stats`; teaching matrix checks | D8 section squares; Q8 compatible squares; both commutators −I; rank-three product map; exact eta and frozen-basis agreement; raw B0 square −I/13 |
| INT-02 | `surface.spin_circle`: `circle_bundle_topology`, `fiber_holonomy`, `rank_one_cohomology`, `spin_circle_metric` | SC7 `seifert_relators`, `fox_row`, `all_cohomology`, `curvature_check` | e=12 torsion; degree ±3 and ±6 acyclicity; trivial/nontrivial base controls; duality uses inverse character; physical gate FAIL retained |
| INT-03 | `surface.hopf_geometry`: `half_form_pencil`, `hopf_connection`, `pullback_line_degrees` | HOPF `qm`, curvature and Chern calculations | Degree 12; anti-self-duality in declared orientation; c2=1 on S4; pulled-back line degrees ±12; arithmetic left action distinguished from right fiber action |
| INT-04 | `hecke` or dedicated surface helpers: `crt_projective_line`, `coprime_block_projectors`, `quaternionic_multiplicity_gate` | HOPF Manin reconstruction, `evaluate`, CRT permutations, rational commutator system | 168 points; cusp widths 1,11,13,143; projectors ranks 12/14, idempotence and complementarity; original 8D quaternionic obstruction and doubled 16D control |
| INT-05 | `research.anomalies`: `chiral_anomaly_polynomial`, `pushforward_curve_flux`, `factorized_descent` | AXG-01 anomaly audit; AXG-02 parent anomaly; AXG-04 independent Cartan reconstruction | I6 cubic/linear/non-Abelian terms with conventions; exact I8→I6; irreducible p2 and f c3 witnesses; zero-index chirality flip changes I8 but preserves pushforward |
| INT-06 | `research.charge_lattices`: `shift_kernel`, `smith_remnant`, `integer_dressing`, `canonical_vector_masses` | AXG-01 anomaly/operator audits; AXG-03 tensor flux audit | K kernel p,Y,B−L; SNF(1,1); explicit dressings/forbidden Majorana; kinetic-normalized spectrum; reject rational-only dressing where integer periods required |
| INT-07 | `research.tensor_gs`: `native_scalar_flux_gate`, `flux_transgression`, `tensor_factorization`, `integral_lattice_gate` | AXG-02 flux/scalar; AXG-03 tensor factor/flux; AXG-04 charged action | Km=(12,12); kΔm=0; integral 3k1; inertia(2,1) witness; C3X Ω and characteristic vector; one-parent fractional norms FAIL; no implicit change of axion period |
| INT-08 | `research.discrete_anomalies`: `spin_zn_fermion_test`, `operator_residue_table` | AXG-03 discrete/operator; AXG-04 charged action | Z3 S1=S3=24 fails cubic mod9; isolated repair ledger includes 12 new modes; C3X Z12 tests; explicit fermion-only vs full-theory scope |
| INT-09 | `research.mode_operators`: `spin_dirac_index`, `degree_one_purity_certificate`, `bochner_bound`, `six_d_scalar_yukawa_gate` | AXG-03 Higgs parent; AXG-04 Higgs action | Index distinct from h0/h1; O(P)→(2,1) versus O(P+Q−R)→(1,0); degree−6 positive scalar bound; same/opposite chirality bilinear gates |
| INT-10 | `research.parents`: immutable M1-parent controls and `c3x_model`; `validate_spectrum`, `action_data`, `reduce_yukawas` | AXG-02 spectrum; AXG-04 chiral/charged/Higgs audits | Six physical C3X reps ×3, one Higgs; signed dimension zero; correct conjugation in 4D; supplied arbitrary 3×3 Yukawas; direct-product global form required |
| INT-11 | `research.compactification`: `einstein_frame_potential`, `product_background`, `radion_hessian`, `remaining_modes` | Initial/SC7 gravity; AXG-01/02 vacuum; AXG-04 gravity background | 6D exponents 4,6,2 and 7D 5,7,3; nonnegative-energy runaways; negative-potential AdS controls; canonical radion; λ4<0 and ell4²<3R² for C3X; untested KK modes marked |
| INT-12 | Experimental topology study, not a generic solver: `c3x_spin_bordism_certificate` | AXG-04 `global_bordism_audit.py` | Explicit Sq2 matrices over F2; dimensions 5,9,16; ranks 2,7; ker=image; stated ordinary Spin-BG category; no automatic promotion to global GS completion |
| INT-13 | Parent comparison controls: `split_bundle_morse_index`, `flux_centralizer` | AXG-04 larger-parent audit/notes | U8 mirror indices, U4×U2×U2 enhancement, complex Morse count 312 with 180 colored directions; do not equate a fixed-YM negative mode with a completed coupled-AdS stability result |
| INT-14 | Documentation and release registers | All reports; teaching manuscript and figures | Append correction links; retain retracted values with provenance; update stale docstrings; ship data, licenses, tests and reproducible CLI examples |

## 4. Suggested sequence of work

**First: preserve and correct.** Import the frozen studies and data. Add the theta-square compatibility distinction. Record the actual archived baseline hash and inspect the target repository for any intervening fixes. Keep the corrected CC-26 contraction and diagnostic classification. Confirm the two mixed-anomaly helpers' factor-of-two convention before combining outputs.

**Second: promote general exact tools.** Implement INT-02–09 as small pure functions with clear inputs. Maintain exact arithmetic for ranks, Smith forms, divisors, anomalies and Clifford identities. Generic routines must return unsupported when a proposed gauge quotient, bundle or mode operator lies outside the implemented assumptions.

**Third: compose named candidates.** Build M1 controls and C3X from explicit immutable inputs. Connect field inventory → anomaly polynomial → charge-lattice tests → mode count → action reduction → vacuum diagnostics. This should expose which failure changes when one ingredient changes. A “viable” boolean is too coarse to represent the current results.

**Fourth: release integration.** Add Legend/registry entries where appropriate, examples and separate exact/diagnostic tests. Ensure imports have no computation or writing side effects. Verify a fresh source distribution contains every datum needed by its gates. Run the repository's required full suite after refactoring; historical check counts are not a substitute.

## 5. Regression strategy

Tests should distinguish mathematical possibilities that earlier reasoning conflated:

| Pair of cases | Regression purpose |
|---|---|
| D8 section action / Q8 theta-compatible action | Preserve the square-map criterion rather than merely checking a group name |
| Torsion flat bundle with nontrivial fiber holonomy / trivial fiber holonomy | Recover acyclicity and nonzero controls |
| Original quartic block / multiplicity-doubled block | Localize the quaternionic obstruction to the actual representation |
| Identity kinetic metric / stack kinetic metric | Prevent raw charge eigenvalues from being presented as physical masses |
| Rational row-span membership / integer dressing | Preserve axion periodicity and finite remnants |
| 4D pushforward agreement / unequal 6D polynomials | Preserve parent sensitivity of zero-index sectors |
| Elementary scalar Laplacian / Dolbeault H1 / gauge Hessian | Prevent a cohomology count from becoming a claim about a different particle operator |
| Index-one O(P) / pure O(P+Q−R) | Verify mirror removal, not only the index |
| Local GS factorization / integral lattice / global partition function | Keep three different consistency gates |
| Stable radion truncation / complete KK stability | Avoid overstating the background control |

Use independent constructions where supplied: direct charges versus polynomial factors, eta transformations versus frozen matrices, Fox derivatives versus Gysin/fiber cohomology, and Cartan expansion versus characteristic-class factorization. Do not replace these by tests that simply assert a hardcoded result returned by the same implementation.

## 6. Machine-readable backlog

`integration_backlog.json` mirrors these task IDs with dependencies, source studies and completion status `PROPOSED`. It is an intake aid, not a declaration that code has been merged or that open research gates have been solved.
