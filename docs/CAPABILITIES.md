# Capability index

This index covers every public Python module in `src/mtft`, including the
early arithmetic/phenomenology tools and current geometry/research engines.
Descriptions come from module docstrings; they describe software scope and
do not independently certify mathematical or physical claims. Check the
[SM correction registers](SM/) and each engine's assumptions before use.

Names in the last column are a few defined entry points, not an exhaustive
API listing. Follow the source link for signatures, exports and limitations.
Use `help(module)` or `python -m mtft.legend search TERM` for more detail.

Regenerate with `python scripts/build_capability_index.py`; use `--check` in
reviews to detect modules missing from this index.

## Core modules

| Module | Scope | Entry points |
|---|---|---|
| [`mtft`](../src/mtft/__init__.py) | MTFT — Modular Time Field Theory | Package exports / data |
| [`mtft.__main__`](../src/mtft/__main__.py) | MTFT Command-Line Interface | `cmd_verify`, `cmd_report`, `cmd_tower`, `cmd_screen` (+3 more in source) |
| [`mtft.al_morphology`](../src/mtft/al_morphology.py) | mtft.al_morphology — what Atkin-Lehner does to shape (v0.24.0). | `al_traces`, `morphology`, `cusp_torsor` |
| [`mtft.arithmetic`](../src/mtft/arithmetic.py) | MTFT Arithmetic Weights and Holonomy Stiffness | `weight`, `weight_array`, `weight_euler`, `damped_weight` (+10 more in source) |
| [`mtft.arithmetic_machine`](../src/mtft/arithmetic_machine.py) | Arithmetic Machine: Computation as a Five-Primitive Object | `Primitive`, `PrimitiveLevel`, `PrimitiveDecomposition`, `decompose_turing_machine` (+29 more in source) |
| [`mtft.arithmetic_wick`](../src/mtft/arithmetic_wick.py) | Arithmetic Wick Rotation: Two Ensembles on the Same Weights | `sieve_primes`, `compute_weights`, `compute_skeleton_weights`, `DirichletEnsemble` (+14 more in source) |
| [`mtft.burning_ship`](../src/mtft/burning_ship.py) | Burning Ship Fractal: The Fermion Vacuum | `burning_ship_iterate`, `burning_ship_array`, `sector`, `jacobian` (+4 more in source) |
| [`mtft.busy_beaver`](../src/mtft/busy_beaver.py) | Arithmetic Busy Beavers: Modular Constraints on Uncomputability | `HeckeSign`, `hecke_sign`, `hecke_sign_pattern`, `hecke_constraint_density` (+29 more in source) |
| [`mtft.chain`](../src/mtft/chain.py) | mtft.chain — the internal (rung-4) model, de-duplicated | `rho`, `Internal`, `internal`, `gap` (+4 more in source) |
| [`mtft.codifferent`](../src/mtft/codifferent.py) | mtft.codifferent — the Canonical Codifferent Theorem at X0(143) (v0.19.0). | `field_trace_table`, `eigen_an`, `gamma_table`, `verify_orbit` (+1 more in source) |
| [`mtft.combinatorial`](../src/mtft/combinatorial.py) | Combinatorial Ancestry — the 2024 lineage as working tools | `bernoulli_plus`, `faulhaber_coeffs`, `power_sum`, `power_sum_backward` (+43 more in source) |
| [`mtft.constants`](../src/mtft/constants.py) | MTFT Constants: The Arithmetic Alphabet | `PhysicalConstants`, `CriticalDepths`, `QuarkMasses`, `LeptonMasses` (+1 more in source) |
| [`mtft.coset_reps`](../src/mtft/coset_reps.py) | mtft.coset_reps — PSL(2,p) accounting for the X0(143) coset layer. | `order_psl2`, `sum_squares_check`, `torus_char_count`, `stage_decomposition` |
| [`mtft.cosmology`](../src/mtft/cosmology.py) | MTFT cosmology: modified Friedmann equations with τ-field dark sector. | `FriedmannMTFT` |
| [`mtft.coupled`](../src/mtft/coupled.py) | mtft.coupled — the spatial sector (rungs 5 and 5b) | `Measure`, `kesten`, `H_of`, `band` (+5 more in source) |
| [`mtft.critical_ensemble`](../src/mtft/critical_ensemble.py) | Critical Ensemble: Li Coefficients as the Third Curvature Family | `lambda_1_closed_form`, `logxi_taylor`, `li_lambda`, `li_lambda_batch` (+8 more in source) |
| [`mtft.curvature`](../src/mtft/curvature.py) | Curvature of the Tano statistical manifold | `brioschi`, `gaussian_family_curvature`, `metric_components`, `gaussian_curvature` (+4 more in source) |
| [`mtft.cuspidal`](../src/mtft/cuspidal.py) | mtft.cuspidal — cuspidal subgroups and Eisenstein torsion (v0.24.0). | `charpoly`, `eisenstein_subspace`, `cuspidal_group`, `eisenstein_kernel_mod2` (+3 more in source) |
| [`mtft.dark_sector`](../src/mtft/dark_sector.py) | τ-field dark sector: vortex halos, flat rotation curves, Tully-Fisher. | `TauVortexHalo`, `rotation_curve`, `tully_fisher`, `rotation_curve_kpc` (+1 more in source) |
| [`mtft.decay`](../src/mtft/decay.py) | Radioactive Decay in Modular Time | `ModularDecay`, `phase_evolution`, `u238_example`, `neutron_beta_decay` |
| [`mtft.dimensional_bridge`](../src/mtft/dimensional_bridge.py) | The Dimensional Bridge: Electron Mass from η(τ) | `charge_from_feigenbaum`, `alpha_F_from_impedance`, `feigenbaum_product_lock`, `modular_impedance` (+2 more in source) |
| [`mtft.divisor_gas`](../src/mtft/divisor_gas.py) | divisor_gas.py — the weights as the energy of the primon gas (PG-01, 2026-10-12). | `divisor_energy`, `divisor_partition`, `divisor_gas`, `weight_is_energy` (+17 more in source) |
| [`mtft.eisenstein`](../src/mtft/eisenstein.py) | Eisenstein congruences — the congruence primes of X_0(143) | `sturm_bound`, `hecke_on_block`, `eisenstein_modulus`, `congruence_census` (+1 more in source) |
| [`mtft.ep`](../src/mtft/ep.py) | mtft.ep — exceptional points: extraction, winding, census, staircase | `levels_of`, `gsq`, `closest_pair`, `newton` (+7 more in source) |
| [`mtft.estimator_standards`](../src/mtft/estimator_standards.py) | estimator_standards.py — A.7 discipline for log-log slope fits (mtft repo). | `binned_log_slope`, `stride_resonance_check`, `recommended_samples_per_decade` |
| [`mtft.exception_spectrum`](../src/mtft/exception_spectrum.py) | mtft.exception_spectrum — the Exception-Spacing Curvature Law (v0.19.0). | `K_from_AB`, `K_atoms`, `K_marked_set`, `two_exception_C` (+5 more in source) |
| [`mtft.expansion`](../src/mtft/expansion.py) | mtft.expansion — the remainder hierarchy, and safe extraction | `incident_couplings`, `channels_A`, `channels_C`, `channels_C3` (+12 more in source) |
| [`mtft.falsify`](../src/mtft/falsify.py) | MTFT Falsifiability Engine | `Prediction`, `CouplingShift`, `coupling_shift`, `coupling_shift_table` (+6 more in source) |
| [`mtft.forms`](../src/mtft/forms.py) | Modular forms relevant to MTFT. | `nome`, `nome_array`, `dedekind_eta`, `dedekind_eta_array` (+7 more in source) |
| [`mtft.gl2_peel`](../src/mtft/gl2_peel.py) | mtft.gl2_peel — the GL(2) peel of f1 = 143a1: BSD rank from the skeleton. | `conductor_certificate`, `ap_point_count`, `lamf_sieve`, `S_direct` (+3 more in source) |
| [`mtft.hardy_ramanujan`](../src/mtft/hardy_ramanujan.py) | mtft.hardy_ramanujan — an orthodox end-to-end benchmark (v0.20.0). | `psi_direct`, `psi_modular`, `modularity_residual`, `layers` (+3 more in source) |
| [`mtft.hecke`](../src/mtft/hecke.py) | Manin symbols and Hecke blocks of X_0(143) | `merel`, `model`, `hecke_matrix`, `cuspidal_hecke` (+7 more in source) |
| [`mtft.hodge_polarization`](../src/mtft/hodge_polarization.py) | mtft.hodge_polarization — pinned polarization data for H1(X0(143)). | Package exports / data |
| [`mtft.hosotani`](../src/mtft/hosotani.py) | Hosotani mechanism for electroweak symmetry breaking in MTFT. | `HosotaniPotential`, `HosotaniMTFT` |
| [`mtft.info_geometry`](../src/mtft/info_geometry.py) | Information geometry layer: the Counting → Geometry → Dynamics bridge. | `logistic_iterate`, `lyapunov_exponent`, `fisher_rao_metric`, `ricci_scalar_logistic` (+2 more in source) |
| [`mtft.integral_lattice`](../src/mtft/integral_lattice.py) | mtft.lattice — exact integer-lattice toolkit (v0.19.0). | `InexactInputError`, `clear_denominators`, `kernel_modp`, `rank_modp` (+11 more in source) |
| [`mtft.jacobian`](../src/mtft/jacobian.py) | The 3×3 Jacobian Stiffness Engine on J₀(143)   (Paper 30) | `JacobianStiffness` |
| [`mtft.jc_counterexample`](../src/mtft/jc_counterexample.py) | The Jacobian Conjecture Counterexample — Machine Certificate | `P_const`, `P_var`, `P_add`, `P_neg` (+19 more in source) |
| [`mtft.kakeya`](../src/mtft/kakeya.py) | mtft.kakeya — Arf parity, direction sets, and finite Kakeya geometry. | `direction_set`, `affine_frame`, `radical_of_parity`, `arf_direction_theorem` (+6 more in source) |
| [`mtft.koide`](../src/mtft/koide.py) | Geometric Koide Theorem | `koide_ratio`, `koide_leptons`, `koide_up_quarks`, `koide_down_quarks` (+10 more in source) |
| [`mtft.lattice`](../src/mtft/lattice.py) | MTFT Lattice Gauge Theory | `random_su_n`, `random_su_n_defective_v0261`, `su_n_identity`, `su_n_center` (+10 more in source) |
| [`mtft.lchannels`](../src/mtft/lchannels.py) | mtft.lchannels — SU(p) gauge filter as Dirichlet L-function channels. | `primitive_root`, `char_table`, `gauss_sum`, `even_js` (+6 more in source) |
| [`mtft.ledger`](../src/mtft/ledger.py) | mtft.ledger — every certified constant of the MTFT program, as data. | `Entry`, `Family`, `rho`, `gap` (+12 more in source) |
| [`mtft.ledger_peel`](../src/mtft/ledger_peel.py) | mtft.ledger_peel — v0.12.0 addendum entries (peel wave), Entry-schema. | `verify` |
| [`mtft.legend`](../src/mtft/legend.py) | The Legend — a map key to the arithmetic territory | `LegendEntry`, `legend_map`, `card`, `trace` (+4 more in source) |
| [`mtft.levels`](../src/mtft/levels.py) | mtft.levels — level-generic modular symbols for squarefree N (v0.24.0). | `UnsupportedLevelError`, `level_data`, `genus`, `is_supported` (+8 more in source) |
| [`mtft.lhcb_analysis`](../src/mtft/lhcb_analysis.py) | MTFT × LHCb Open Data Analysis Bridge | `LHCbNtuple`, `combine_ntuples`, `setup_check` |
| [`mtft.liealg`](../src/mtft/liealg.py) | Numerical Lie-algebra fingerprinting for subalgebras of u(n). | `LieGateAmbiguous`, `inner`, `antiherm`, `vec_u` (+16 more in source) |
| [`mtft.marked_gap`](../src/mtft/marked_gap.py) | mtft.marked_gap — Rung-4 mass gap in two-temperature ensemble language. | `Lambda_of`, `eps_level`, `delta_gap`, `predicted_spectrum` (+1 more in source) |
| [`mtft.marked_gas`](../src/mtft/marked_gas.py) | marked_gas.py — The Marked Primon Gas: a KMS construction for the | `Certified`, `z1`, `z2`, `zD_certified_interval` (+12 more in source) |
| [`mtft.modular`](../src/mtft/modular.py) | Modular time field and SL(2,ℤ) geometry on the upper half-plane ℍ. | `sl2z_transform`, `reduce_to_fundamental`, `hyperbolic_distance`, `poincare_metric` (+1 more in source) |
| [`mtft.modular_curve`](../src/mtft/modular_curve.py) | mtft.modular_curve | `CuspClass`, `VortexConfig`, `HeckeSpectrum`, `HomologyData` (+2 more in source) |
| [`mtft.moments`](../src/mtft/moments.py) | Tano weight moments — closed forms for the arithmetic ensemble | `sieve_primes`, `primesum_logk`, `rep_mul`, `S_r` (+18 more in source) |
| [`mtft.monster_hash`](../src/mtft/monster_hash.py) | MonsterHash — SL(2,Z)-Sponge Hash Function | `MonsterHash`, `compare_hashes` |
| [`mtft.music`](../src/mtft/music.py) | MTFT Music Module — Arithmetic Sonification Engine | `VacuumSonifier`, `ModularScale`, `Note`, `Phrase` (+4 more in source) |
| [`mtft.parent`](../src/mtft/parent.py) | mtft.parent — candidate parent actions as data, and the parent-action ledger's gates as verdicts (PV-01, 2026-09-29). | `ParentCandidate`, `gates`, `ledger_table` |
| [`mtft.particles`](../src/mtft/particles.py) | Standard Model particle database with MTFT modular-time embeddings. | `ParticleType`, `Particle`, `StandardModel` |
| [`mtft.peel`](../src/mtft/peel.py) | mtft.peel — Mellin peel engine for the bulk and skeleton stiffness. | `w_sieve`, `lambda_sieve`, `F_bulk`, `mu_bulk_direct` (+5 more in source) |
| [`mtft.quadratic_forms`](../src/mtft/quadratic_forms.py) | mtft.quadratic_forms — the Gauss-Legendre three-squares layer (v0.19.0). | `v2`, `forbidden`, `forbidden_projector`, `r3_array` (+3 more in source) |
| [`mtft.quantum`](../src/mtft/quantum.py) | MTFT Quantum Computing Primitives | `gell_mann_matrices`, `HolonomyGate`, `holonomy_gate_set`, `TopologicalQudit` (+6 more in source) |
| [`mtft.riemann`](../src/mtft/riemann.py) | MTFT × Riemann: Explicit Formula and ζ-Zero Connection | `gamma_suppression`, `gamma_suppression_table`, `zero_contribution`, `zero_sum` (+24 more in source) |
| [`mtft.tano_metric`](../src/mtft/tano_metric.py) | MTFT Materials Science Layer: The Tano Metric | `Element`, `get_element`, `tano_contrast`, `geometry_index` (+5 more in source) |
| [`mtft.thetafun`](../src/mtft/thetafun.py) | Numerical genus-g theta functions with characteristics for X0(143). | `char_to_ab`, `lll_reduce_gram`, `siegel_ready`, `reduce_char` (+7 more in source) |
| [`mtft.tower`](../src/mtft/tower.py) | MTFT Multi-N Tower: Yang-Mills Confinement Landscape | `tower_stiffness`, `even_n_universality`, `confinement_boundary`, `phase_transition_scaling` (+5 more in source) |
| [`mtft.verify`](../src/mtft/verify.py) | MTFT Verification Scorecard | `Prediction`, `gauge_sector`, `higgs_sector`, `mass_ratios` (+5 more in source) |
| [`mtft.viz`](../src/mtft/viz.py) | MTFT Visualization Helpers | `stiffness_landscape`, `hosotani_potential_plot`, `rotation_curve_plot`, `koide_manifold_3d` (+3 more in source) |
| [`mtft.weil`](../src/mtft/weil.py) | mtft.weil -- Gabor-compressed Weil explicit-formula form (CANDIDATE, W1). | `Window`, `gabor`, `prime_powers`, `nu_parts` (+8 more in source) |
| [`mtft.x0_143`](../src/mtft/x0_143.py) | The Modular Curve X₀(143): Arithmetic Stage of MTFT | `EllipticCurve143a1`, `hecke_polynomial_f2_T2`, `hecke_polynomial_f3_T2`, `hecke_polynomial_f2_T3` (+12 more in source) |

## boundary

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.boundary`](../src/mtft/boundary/__init__.py) | mtft.boundary — the cusp/elliptic boundary layer of the canonical ring. | `data_path`, `s4_qexpansions`, `cusp_functionals`, `operator` (+3 more in source) |
| [`mtft.boundary.gates`](../src/mtft/boundary/gates.py) | Gates for `mtft.boundary`. | `gate_census_consistency`, `gate_h0_2k`, `gate_cusp_functionals`, `gate_leakage` (+3 more in source) |

## canonical

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.canonical`](../src/mtft/canonical/__init__.py) | mtft.canonical — the canonical ideal of X0(143) and the Atkin-Lehner descent. | `data_path`, `s2_qexpansions`, `adapted_basis`, `adapted_qexpansions` (+6 more in source) |
| [`mtft.canonical.gates`](../src/mtft/canonical/gates.py) | Gates for `mtft.canonical`. | `gate_petri`, `gate_generation`, `gate_sector_grading`, `gate_bundles` (+7 more in source) |
| [`mtft.canonical.integral`](../src/mtft/canonical/integral.py) | mtft.canonical.integral — integral models, saturation, and the 2026-08-24 arc. | `adapted_matrix`, `saturated_qexpansions`, `sector_columns`, `count_points_modp` (+6 more in source) |
| [`mtft.canonical.integral_gates`](../src/mtft/canonical/integral_gates.py) | Gates for `mtft.canonical.integral` — the 2026-08-24 arc, recomputed live. | `gate_frame`, `gate_integral_model`, `gate_counts_mod2`, `gate_counts_mod3` (+4 more in source) |

## crypto

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.crypto`](../src/mtft/crypto/__init__.py) | MTFT Cryptographic Primitives | `ArithmeticHash`, `BurningShipPRNG`, `sl2z_power`, `ModularKeyExchange` (+3 more in source) |
| [`mtft.crypto.jacobian_order`](../src/mtft/crypto/jacobian_order.py) | mtft.crypto.jacobian_order | `JacobianOrder`, `JacobianOrderTable`, `total_jacobian_order`, `hasse_weil_bounds` (+1 more in source) |

## gprun

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.gprun`](../src/mtft/gprun/__init__.py) | mtft.gprun — a local job runner for long PARI/GP computations. | `find_gp`, `Job`, `serve`, `main` |
| [`mtft.gprun.__main__`](../src/mtft/gprun/__main__.py) | See module source for its public interface. | Package exports / data |

## homology

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.homology`](../src/mtft/homology/__init__.py) | Canonical integral homology of X0(143). | `standard_J`, `int_inverse`, `is_symplectic`, `is_anti_symplectic` (+4 more in source) |

## interactions

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.interactions`](../src/mtft/interactions/__init__.py) | SM interaction reference and MTFT evidence ledger (no amplitude evaluation). | Package exports / data |
| [`mtft.interactions.__main__`](../src/mtft/interactions/__main__.py) | Command-line interface to the SM reference catalog and mapping ledger. | `main` |
| [`mtft.interactions.catalog`](../src/mtft/interactions/catalog.py) | Reference catalog queries and evidence-ledger validation. | `load_catalog`, `catalog_digest`, `validate_catalog`, `select_vertices` (+3 more in source) |
| [`mtft.interactions.render`](../src/mtft/interactions/render.py) | Portable, offline HTML and deterministic SVG for the interaction atlas. | `vertex_svg`, `render_html` |
| [`mtft.interactions.ufo`](../src/mtft/interactions/ufo.py) | Read a conservative, declarative subset of tree-level UFO without executing it. | `UFOError`, `import_ufo` |

## matter

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.matter`](../src/mtft/matter/__init__.py) | Matter-01: explicitly labelled reference quantum-field calculations. | Package exports / data |
| [`mtft.matter.__main__`](../src/mtft/matter/__main__.py) | Matter command line: Matter-01's measured-input QED benchmark, and the matter-programme toolkit (presets, background, ledger, epsilon, puzzle, levels, timerev, prereg, forces, atlas). | `main` |
| [`mtft.matter.arakelov_kernel`](../src/mtft/matter/arakelov_kernel.py) | AG-01 (2026-10-07): the canonical Green's-function kernel of X0(143) — the pre-registered test of FA-01, run. | `multiplier_orders`, `kernel`, `spread`, `ratio_to_flat` (+5 more in source) |
| [`mtft.matter.background`](../src/mtft/matter/background.py) | mtft.matter.background — one object for a Pati–Salam background on blow-ups of X0(143) x 143a1 (toolkit, 2026-10-02). | `Background`, `preset`, `ps_symmetric_gamma_tradeoff` |
| [`mtft.matter.bl_wilson_line`](../src/mtft/matter/bl_wilson_line.py) | BL-01 (2026-10-06): a B–L Wilson line on 143a1 — the right shape, the wrong size, and an unstable vacuum. | `tau_143a1`, `half_period`, `character`, `image_sum` (+9 more in source) |
| [`mtft.matter.bridged_top`](../src/mtft/matter/bridged_top.py) | YB-01 (2026-10-08): the bridged top — the second redesign fails its pre-registered test. | `tail_value`, `jet_factor`, `bridged_ratio`, `y_bridged` (+1 more in source) |
| [`mtft.matter.cabibbo`](../src/mtft/matter/cabibbo.py) | CA-01 (2026-10-05): the Cabibbo angle from the gamma-flux, and two measurements of the curve's size. | `r_for_cabibbo`, `r_from_vcb`, `gap`, `c_bracket` (+2 more in source) |
| [`mtft.matter.canonical_shape`](../src/mtft/matter/canonical_shape.py) | CS-01 (2026-10-05): is there a canonical prescription for the curve-to-torus area ratio?  Pre-registered test. | `chance_per_candidate`, `candidates`, `invariants`, `hit` (+1 more in source) |
| [`mtft.matter.charged_modes`](../src/mtft/matter/charged_modes.py) | Matter-03 charged holomorphic Higgs-deformation complex. | `split_charged_census`, `paired_charged_report`, `certified_paired_charged_report`, `main` |
| [`mtft.matter.ckm_cusps`](../src/mtft/matter/ckm_cusps.py) | CK-01 (2026-10-02): quark mixing and the cusps of X0(143). | `cusp_distances`, `mixing_bounds`, `split_assignments`, `texture` (+5 more in source) |
| [`mtft.matter.collinear`](../src/mtft/matter/collinear.py) | Matter-02 exact cohomology and a holomorphic Yukawa rank witness. | `Cohomology`, `cohomology`, `cusp_jet_certificate`, `YukawaWitness` (+3 more in source) |
| [`mtft.matter.condensate_dressing`](../src/mtft/matter/condensate_dressing.py) | CD-01 (2026-10-05): the condensate dresses u^c and e^c but never d^c — and that carries the up–down asymmetry. | `name`, `one_step`, `dressed`, `undressed` (+2 more in source) |
| [`mtft.matter.condensate_insertion`](../src/mtft/matter/condensate_insertion.py) | CI-01 (2026-10-06): the condensate insertion — the lepton pattern lives in the tails, where the condensate cannot reach. | `r_B`, `r_A`, `product_relation`, `tail_shortfall` (+5 more in source) |
| [`mtft.matter.condensate_phases`](../src/mtft/matter/condensate_phases.py) | CP-01 (2026-10-03): the condensate's jet alignment — spin-1 phases, the tree-level potential, and a Georgi–Jarlskog structure in the texture. | `sphere_moment`, `poly`, `norm2`, `quartic` (+19 more in source) |
| [`mtft.matter.conjugate_orientation`](../src/mtft/matter/conjugate_orientation.py) | CO-01 (2026-10-05): is there a conjugate condensate?  No — and the one thing that tells right from left. | `d_term_potential`, `d_term_minimum`, `lifted_minimum`, `scalar_masses` (+2 more in source) |
| [`mtft.matter.cusp_families`](../src/mtft/matter/cusp_families.py) | AL-01 (2026-10-02): blowing up cusps of X0(143) x 143a1 — what Atkin–Lehner symmetry allows, and the one-family-per-cusp background. | `dot`, `K`, `chi`, `jets` (+8 more in source) |
| [`mtft.matter.double_duty`](../src/mtft/matter/double_duty.py) | DD-01 (2026-10-05): making the shape do double duty — a parity selection rule for the CKM matrix. | `parities`, `allowed`, `decoupled_family`, `required_coefficient` (+3 more in source) |
| [`mtft.matter.down_parity`](../src/mtft/matter/down_parity.py) | DS-01 (2026-10-07): the down sector is parity-symmetric — and V_ub then has a floor the data sit below. | `touch`, `touched`, `angles_2x2`, `c_floor` (+4 more in source) |
| [`mtft.matter.dressing_profiles`](../src/mtft/matter/dressing_profiles.py) | DP-01 (2026-10-05): how strong the dressing must be, and why one dressing cannot serve both quarks and leptons. | `weight`, `masses_for_z`, `z_leptons`, `conflict` (+2 more in source) |
| [`mtft.matter.epsilon`](../src/mtft/matter/epsilon.py) | EPS-01 (2026-10-02): what sets epsilon, the family-hierarchy parameter of CK-01. | `allowed_entries`, `d_charges_from_condensate`, `texture`, `texture_properties` (+8 more in source) |
| [`mtft.matter.flavour_audit`](../src/mtft/matter/flavour_audit.py) | FA-01 (2026-10-06): the flavour audit — one shared artefact found, four independent failures confirmed, and the placement fixed by arithmetic. | `vcb_images`, `vcb_2d`, `sx_for`, `r_images` (+6 more in source) |
| [`mtft.matter.flux_scan`](../src/mtft/matter/flux_scan.py) | FS-01 (2026-10-02): scanning SU(2)_R-breaking fluxes at the single blown-up cusp — the toolkit's stack-level extension. | `StackFlux`, `scan`, `best`, `report` (+3 more in source) |
| [`mtft.matter.flux_tails`](../src/mtft/matter/flux_tails.py) | FL-01 (2026-10-05): which tail sets which angle, what the data demand, and why curve flux cannot supply it. | `required_angles`, `required_tails`, `flux_ratio`, `best_flux_ratio` (+1 more in source) |
| [`mtft.matter.forces`](../src/mtft/matter/forces.py) | GF-01 (2026-10-11): forces first — the eight dots of the parent U(8), the lines between them, and what the forces fix before a particle is chosen. | `e`, `vec`, `line`, `field_vectors` (+51 more in source) |
| [`mtft.matter.fricke_home`](../src/mtft/matter/fricke_home.py) | FH-01 (2026-10-07): the top at the Fricke point — parity holds, four r measurements meet, and V_ub tests the parity breaking, not the home. | `tau_143a1`, `shortest_dual`, `separation`, `torus_suppression` (+4 more in source) |
| [`mtft.matter.gauge_couplings`](../src/mtft/matter/gauge_couplings.py) | GC-01 (2026-10-05): what the gauge couplings fix — not the curve's shape, but the Pati–Salam scale, and they disagree with b–tau. | `ratio`, `matching_scale`, `alpha23_scale`, `design_ratios` (+4 more in source) |
| [`mtft.matter.gaugino_towers`](../src/mtft/matter/gaugino_towers.py) | GT-01 (2026-10-05): the two gaugino towers — the semi-Dirac-like split exists, with the right shape and the wrong strength. | `share`, `tau_at`, `lepton_bound`, `reversal_requirements` (+2 more in source) |
| [`mtft.matter.green_schwarz`](../src/mtft/matter/green_schwarz.py) | GS-01 (2026-10-05): Green–Schwarz relaxation, the slepton direction, and what protects the vector-like pairs. | `d_potential`, `minimum_fixed_xi`, `minimum_dynamical_xi`, `su4_shift` (+3 more in source) |
| [`mtft.matter.higgs_profile`](../src/mtft/matter/higgs_profile.py) | HP-01 (2026-10-04): the light Higgs of the ratio-3 background — where it lives, what type it is, and why it cannot carry the Yukawas. | `localisation`, `phi_taylor`, `chi_entry`, `superpotential_allowed` (+2 more in source) |
| [`mtft.matter.higgs_vacuum`](../src/mtft/matter/higgs_vacuum.py) | Matter-03 exact inputs to a conditional classical Higgs-branch theorem. | `single_arrow_gate`, `HiggsPairWitness`, `paired_higgs_witness`, `classical_gauge_certificate` (+2 more in source) |
| [`mtft.matter.jet_vertex`](../src/mtft/matter/jet_vertex.py) | HJ-01 (2026-10-04): the family vertex of the jet classes — selection rule, locality, the round blow-up, and where X_0(143) enters. | `curve_charge`, `round_texture`, `singular_ratios`, `newform_coefficients` (+4 more in source) |
| [`mtft.matter.kinetics`](../src/mtft/matter/kinetics.py) | Matter-03 field-coordinate normalization with supplied physical inputs. | `canonical_frame`, `canonical_yukawa`, `canonical_mass_matrix`, `mass_singular_values` |
| [`mtft.matter.lattice_breaking`](../src/mtft/matter/lattice_breaking.py) | LB-01 (2026-10-04): how the lattice of 143a1 breaks the blown-up point's rotation, and whether jet families can carry an O(1) Yukawa. | `lattice_data`, `blowup_radius`, `breaking`, `parity_channels` (+3 more in source) |
| [`mtft.matter.ledger`](../src/mtft/matter/ledger.py) | mtft.matter.ledger — machine-readable ledger of the matter programme from Matter-04 to EPS-01 (toolkit, 2026-10-02). | `entries`, `to_markdown`, `check` |
| [`mtft.matter.leptons`](../src/mtft/matter/leptons.py) | LEP-01 (2026-10-02): the charged-lepton texture — the first place the construction touches an atom. | `koide`, `gj_at_matching`, `vertex_charges`, `su4_identity_holds` (+7 more in source) |
| [`mtft.matter.lever`](../src/mtft/matter/lever.py) | LV-01 (2026-10-07): the lever at the jet point — what V_ub measures, and a correction to FH-4. | `y_down`, `down_only`, `formulas`, `opposite_arms` (+2 more in source) |
| [`mtft.matter.messenger_yukawa`](../src/mtft/matter/messenger_yukawa.py) | PB-01 (2026-10-02): messenger-exchange Yukawas for the jet families of the Pati–Salam branch. | `su_invariants`, `holomorphic_orders`, `d_term_vacuum`, `d_term_vacuum_scan` (+3 more in source) |
| [`mtft.matter.moduli_casimir`](../src/mtft/matter/moduli_casimir.py) | MS-01 (2026-10-05): can the Scherk–Schwarz Casimir energy fix the curve's size?  Only if one sign comes out right — and a first count says it doesn't. | `odd_sums`, `minimum_ratio`, `a1_counts`, `dolbeault_a1` (+4 more in source) |
| [`mtft.matter.neutrino_gate`](../src/mtft/matter/neutrino_gate.py) | NU-01 (2026-10-02): the neutrino gate of the Pati–Salam branch — what neutral fermions the construction provides, where they live, which masses | `stack_class`, `entry`, `fi_check`, `localised_extras` (+6 more in source) |
| [`mtft.matter.neutrino_mu`](../src/mtft/matter/neutrino_mu.py) | NU-02 (2026-10-02): the size of the lepton-number-violating mass mu — the arithmetic geometry of the localisation points and the | `newform_143a`, `ModularParametrisation`, `act`, `jpair` (+5 more in source) |
| [`mtft.matter.normalization`](../src/mtft/matter/normalization.py) | Conditional U(1) dimensional reduction, not a selected MTFT vacuum. | `canonical_u1_vertex`, `four_dimensional_coupling`, `electroweak_matching`, `breathing_gauge_kinetics` (+1 more in source) |
| [`mtft.matter.odd_modes`](../src/mtft/matter/odd_modes.py) | OM-01 (2026-10-08): the forced-odd-mode count — the parity audit's verdict made topological. | `split_by_poles`, `lefschetz`, `first_forced_odd_degree`, `b_stack_blocks` (+2 more in source) |
| [`mtft.matter.one_loop_adjoints`](../src/mtft/matter/one_loop_adjoints.py) | Matter-05 (2026-10-02): one-loop masses of the 27 bosonic adjoints per gauge generator under the Scherk–Schwarz twist of Matter-04 — | `tau_143a1`, `lattice_sums`, `spin_structure_table`, `mass_coefficients` (+6 more in source) |
| [`mtft.matter.oracle`](../src/mtft/matter/oracle.py) | Independent massive Born e- e+ -> gamma gamma reference formulas. | `mandelstam`, `msq_mandelstam`, `msq_cm`, `dsigma_domega` (+5 more in source) |
| [`mtft.matter.overlaps`](../src/mtft/matter/overlaps.py) | XQ-01 (2026-10-03): the overlaps X_qq' — family normalisation, the curve factor from X_0(143), and what they do to the texture. | `cusp_bergman_jets`, `correlations`, `antidiagonal_weights`, `texture_masses` (+3 more in source) |
| [`mtft.matter.parity_audit`](../src/mtft/matter/parity_audit.py) | PAR-01 (2026-10-08): the parity audit (filed as PA-01, which is Path A's label; relabelled by CC-62; claims PA-1 to PA-4) — every source odd under the torus involution, what it reaches at first order, and the verdict on the two-point design. | `down_sector_odd_sources`, `odd_moduli`, `floor_vub`, `loophole` (+1 more in source) |
| [`mtft.matter.path_a`](../src/mtft/matter/path_a.py) | PA-01 (2026-10-02): Path A — the matter-untwisted completion of YK-01: its vector-like states, what could give them mass, and whether the background | `content`, `mass_sources`, `T_charges`, `potential_profile` (+3 more in source) |
| [`mtft.matter.pati_salam`](../src/mtft/matter/pati_salam.py) | Pati–Salam first look (2026-10-02): the minimal fermion flux of Matter-04 regrouped as SU(4) x SU(2)_L x SU(2)_R.  Index counts only. | `group_phi`, `block_counts`, `right_handed_assignment`, `first_look` |
| [`mtft.matter.prereg`](../src/mtft/matter/prereg.py) | PR-01 (2026-10-10): sealed pre-registrations and the trials counter. | `seal_hash`, `verify`, `evaluate`, `trials` (+3 more in source) |
| [`mtft.matter.ps_breaking`](../src/mtft/matter/ps_breaking.py) | PS-01 (2026-10-02): the Pati–Salam branch of the minimal Scherk–Schwarz flux — how it breaks, what is anomalous, and the first comparison | `stack_hypercharge`, `block_hypercharge`, `branching`, `anomalies` (+8 more in source) |
| [`mtft.matter.puzzle`](../src/mtft/matter/puzzle.py) | mtft.matter.puzzle — the programme map: every piece of the Standard-Model puzzle on X0(143) x 143a1, its status, and its own puzzles. | `Sub`, `Piece`, `ordered`, `pieces` (+17 more in source) |
| [`mtft.matter.qed`](../src/mtft/matter/qed.py) | Matter-01: a massive tree-level QED reference process. | `minkowski_dot`, `gamma_matrices`, `slash`, `u_spinor` (+12 more in source) |
| [`mtft.matter.quark_lepton_split`](../src/mtft/matter/quark_lepton_split.py) | QL-01 (2026-10-06): the quark–lepton split — the b–tau bound closes the dressing route, and names what the vertex must supply. | `z_e`, `first_over_third`, `dressing_bound`, `tau_over_b` (+5 more in source) |
| [`mtft.matter.ratio_three`](../src/mtft/matter/ratio_three.py) | RT-01 (2026-10-02): the ratio-3 working background of FS-01 — does eps_d/eps_u = 3 survive its order-one factors, and what do its neutrinos do? | `leading_ratio`, `eps_sector`, `corrected_ratio_window`, `neutral_matrix` (+3 more in source) |
| [`mtft.matter.report`](../src/mtft/matter/report.py) | Deterministic numerical diagnostics and a browser-ready QED reference report. | `numerical_gates`, `benchmark`, `benchmark_report` |
| [`mtft.matter.scales`](../src/mtft/matter/scales.py) | SC-01 (2026-10-17): the scale ladder of the twisted compactification — what is fixed, by what, and what the measured background implies. | `tau`, `lattice_sums`, `ladder`, `g4_squared` (+3 more in source) |
| [`mtft.matter.special_points`](../src/mtft/matter/special_points.py) | SP-01 (2026-10-05): where the top can live — distances on X_0(143) in its canonical metric. | `density_setup`, `arakelov_density`, `path_length`, `curve_area_for_Vcb` (+2 more in source) |
| [`mtft.matter.spin_twists`](../src/mtft/matter/spin_twists.py) | Matter-04 (2026-10-02): the adjoint question answered through spinors — spin structures on 143a1, a Scherk–Schwarz twist of the | `two_division_polynomial`, `elliptic_spin_structures`, `flat_cohomology_curve`, `kunneth` (+11 more in source) |
| [`mtft.matter.star_quotients`](../src/mtft/matter/star_quotients.py) | SQ-01 (2026-10-08): the star quotients X0(pq)* for supersingular p < q — is X0(143) special? | `kronecker`, `conductor`, `eichler`, `class_number` (+23 more in source) |
| [`mtft.matter.third_observable`](../src/mtft/matter/third_observable.py) | TO-01 (2026-10-05): a pre-registered third observable — and why V_cb never measured the curve's size. | `sector_ratios`, `vcb_range`, `direct_entry_bound`, `passes` (+2 more in source) |
| [`mtft.matter.thresholds`](../src/mtft/matter/thresholds.py) | TH-01 (2026-10-05): can thresholds reconcile the gauge couplings with the b–tau scale? | `pair_b`, `slopes`, `N_D`, `required_nl` (+4 more in source) |
| [`mtft.matter.tilt`](../src/mtft/matter/tilt.py) | TL-01 (2026-10-07): the gamma tilt is a rotation — it cancels between up and down quarks, and where it survives, V_ub measures the down quarks' share of V_cb. | `sector`, `ckm`, `left_angle23`, `block_tilts` (+4 more in source) |
| [`mtft.matter.time_reversal`](../src/mtft/matter/time_reversal.py) | TR-01 (2026-10-10): time reversal in the programme — the two mirror classes of the modular tiling, the period lattices of the thirteen | `is_reflection`, `reflect`, `fixed_form`, `content` (+23 more in source) |
| [`mtft.matter.top_family`](../src/mtft/matter/top_family.py) | TF-01 (2026-10-05): can any jet configuration host the top?  No — and what to build instead. | `min_breaking_units`, `two_fold_anisotropy`, `flux_family_design`, `report` |
| [`mtft.matter.top_point`](../src/mtft/matter/top_point.py) | TP-01 (2026-10-05): the flux family fails, the constant-jet top survives — a two-point design. | `flux_family_chirality`, `global_cost`, `two_point_families`, `higgs_count` (+2 more in source) |
| [`mtft.matter.two_point_audit`](../src/mtft/matter/two_point_audit.py) | TA-01 (2026-10-05): auditing the two-point design — the top's blow-up is too big to fit beside the jets on one torus. | `fibre_geometry`, `radius`, `blowups`, `tail` (+3 more in source) |
| [`mtft.matter.two_scales`](../src/mtft/matter/two_scales.py) | GC-02 (2026-10-17): two scales, not one — the flux-split couplings matched at the curve's compactification scale M_C, above the | `ps_b`, `ps_inverse_couplings`, `two_scale_ratio`, `compactification_scale` (+1 more in source) |
| [`mtft.matter.vacuum_atlas`](../src/mtft/matter/vacuum_atlas.py) | PG-01 (2026-10-12): the vacuum as a dynamical system — the critical-point atlas of the potentials the programme has, the shape modulus's cosh well, | `classify`, `gradient_flow_eigs`, `shape_well`, `volume_runaway` (+15 more in source) |
| [`mtft.matter.vacuum_energy`](../src/mtft/matter/vacuum_energy.py) | VE-01 (2026-10-17): the vacuum energy of the twisted compactification — the programme's own terms, their signs, and the one point | `tau`, `odd_windings`, `lattice_sums`, `trace_holonomy` (+18 more in source) |
| [`mtft.matter.vacuum_stability`](../src/mtft/matter/vacuum_stability.py) | VS-01 (2026-10-05): the condensate's feedback on the (4,2,1) scalars — a D-flat Q block breaks colour. | `su4_pair`, `scalar_masses_4L`, `stable`, `stability_bound` (+3 more in source) |
| [`mtft.matter.vector_like_masses`](../src/mtft/matter/vector_like_masses.py) | VM-01 (2026-10-05): what gives the vector-like pairs their masses — not the Pati–Salam condensate, at any order, in our vacuum. | `su4_invariant_possible`, `operators`, `vacuum_sources`, `tree_level_mass_possible` (+2 more in source) |
| [`mtft.matter.yukawa_audit`](../src/mtft/matter/yukawa_audit.py) | YA-01 (2026-10-04): a character audit of the Yukawa amplitude for jet families — why the minimal background does not rescue the leading coupling, and | `higgs_character`, `holomorphic_amplitude_character`, `forbidden_holomorphic`, `curve_allowed` (+5 more in source) |
| [`mtft.matter.yukawa_jets`](../src/mtft/matter/yukawa_jets.py) | YK-01 (2026-10-02): Yukawa couplings of exceptional-jet families — a vanishing theorem, the only tree-level way out, and its price. | `twist_table`, `tree_yukawa_possible`, `classify`, `sm_pattern` (+6 more in source) |

## origami

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.origami`](../src/mtft/origami/__init__.py) | mtft.origami — dimers, origami/t-embeddings, and the observable insertion calculus. | Package exports / data |
| [`mtft.origami.dimer`](../src/mtft/origami/dimer.py) | mtft.origami.dimer — weighted planar bipartite graphs and their dimer ensembles. | `DimerGraph`, `ensemble_conservation` |
| [`mtft.origami.gates`](../src/mtft/origami/gates.py) | mtft.origami.gates — the v0.20.0 gate battery, runnable end to end. | `gate_boundary_measurement_24`, `gate_mandelstam_24`, `gate_ensemble_conservation_24`, `gate_simplex_curvature` (+11 more in source) |
| [`mtft.origami.insertion`](../src/mtft/origami/insertion.py) | mtft.origami.insertion — the observable insertion calculus. | `D_log`, `cumulants`, `fisher_metric`, `cubic_tensor` (+4 more in source) |
| [`mtft.origami.instances`](../src/mtft/origami/instances.py) | mtft.origami.instances — the two certified instances. | `galashin_24`, `t_embedding_24`, `mandelstams_24`, `closed_curvature_B` (+1 more in source) |
| [`mtft.origami.perfect`](../src/mtft/origami/perfect.py) | mtft.origami.perfect — perfect t-embeddings and their branch structure. | `bracket`, `t_coefficients`, `Theta`, `winding` (+7 more in source) |

## periods

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.periods`](../src/mtft/periods/__init__.py) | mtft.periods — genuine period/Hodge geometry of X0(143). | Package exports / data |
| [`mtft.periods.__main__`](../src/mtft/periods/__main__.py) | See module source for its public interface. | `main` |
| [`mtft.periods.bridge`](../src/mtft/periods/bridge.py) | Exact bridge between the v6 period-Manin basis and ``mtft.hecke``. | `relative_basis_change`, `cuspidal_basis_change`, `hecke_to_symplectic_change` |
| [`mtft.periods.channels`](../src/mtft/periods/channels.py) | Bergman harmonic channels on X_0(143). | `bergman_bilinear`, `bergman_channel`, `channel_density`, `mode_crossover` |
| [`mtft.periods.core`](../src/mtft/periods/core.py) | Core period and Hodge geometry for X_0(143). | `data_path`, `period_record`, `intersection_inverse`, `intersection_form` (+11 more in source) |
| [`mtft.periods.forms`](../src/mtft/periods/forms.py) | Native q-expansion evaluation and Bergman density for X_0(143). | `raw_qexpansions`, `q_tail_bound`, `raw_form_values`, `normalized_form_values` (+1 more in source) |
| [`mtft.periods.gates`](../src/mtft/periods/gates.py) | Call-time gates for ``mtft.periods``.  Fast tier is no-PARI. | `gate_integral_symplectic`, `gate_period_reconstruction`, `gate_basis_bridge`, `gate_hodge_bridge` (+7 more in source) |
| [`mtft.periods.hamiltonian`](../src/mtft/periods/hamiltonian.py) | Quadratic-Hamiltonian layer on the X_0(143) Hodge stage. | `hodge_adjoint`, `hermitian_split`, `hamiltonian_split`, `channel_report` (+4 more in source) |
| [`mtft.periods.involutions`](../src/mtft/periods/involutions.py) | Involutions on the promoted X_0(143) homology stage. | `al_matrix`, `transported_intersection`, `hecke_blocks`, `al_signs` (+5 more in source) |
| [`mtft.periods.oldtorus`](../src/mtft/periods/oldtorus.py) | The level-11 oldspace of X_0(143) as an arithmetic abelian surface. | `old_lattice`, `intersection_on_old`, `polarization_type`, `j_arith` (+5 more in source) |
| [`mtft.periods.physics`](../src/mtft/periods/physics.py) | Physics-facing, but epistemically separated, tools built from X_0(143). | `hodge_structure_hecke`, `hodge_metric_hecke`, `graph_coupling`, `complex_linear_decomposition` (+3 more in source) |

## research

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.research`](../src/mtft/research/__init__.py) | Experimental research namespace (Wave B/C of the v0.32.0 plan): parent theories, anomaly polynomials, charge lattices, | Package exports / data |
| [`mtft.research.anomalies`](../src/mtft/research/anomalies.py) | INT-05 (AXG-01 §3): 4D anomaly polynomial of a stack model with line fluxes, and its factorization through the shift matrix. | `stack_model`, `anomaly_polynomials`, `directional`, `polarised_trace` (+3 more in source) |
| [`mtft.research.anomaly8d`](../src/mtft/research/anomaly8d.py) | PV-01 (2026-09-29): the ten-form anomaly polynomial of chiral eight-dimensional field content, exactly — the gate the | `Ch`, `s`, `rank`, `fund` (+24 more in source) |
| [`mtft.research.bordism`](../src/mtft/research/bordism.py) | INT-12 (AXG-04 §5): the C3X ordinary spin-bordism certificate, Omega_7^Spin(B(SU3 x SU2 x U1^2)) = 0. | `c3x_spin_bordism_certificate` |
| [`mtft.research.charge_lattices`](../src/mtft/research/charge_lattices.py) | INT-06 (AXG-01 §§4, 6; AXG-03 §4): shift kernel, Smith remnants, integer dressings and kinetic-normalised vector masses. | `shift_kernel`, `smith_remnant`, `integer_dressing`, `canonical_vector_masses` |
| [`mtft.research.chirality`](../src/mtft/research/chirality.py) | R2C-01 (Astra, 2026-09-18; verified independently 2026-09-19): 6D chirality assignments for M1's ten oriented bifundamentals. | `sector_degrees`, `clifford_6d`, `bilinear_selection`, `signed_index` (+3 more in source) |
| [`mtft.research.compactification`](../src/mtft/research/compactification.py) | INT-11 (initial audit §7, SC7-01 §8, AXG-01/02 vacuum, AXG-04 §8): Einstein-frame radius potentials and the unwarped | `einstein_frame_potential`, `product_background`, `c3x_ads_control`, `planck_reduction` (+1 more in source) |
| [`mtft.research.discrete_anomalies`](../src/mtft/research/discrete_anomalies.py) | INT-08 (AXG-03 §5): Spin x Z_n fermion anomaly test (Hsieh) on a net chiral charge ledger. | `m1_ledger`, `spin_zn_fermion_test`, `m1_z3_generator_certificate` |
| [`mtft.research.gravitational_anomaly`](../src/mtft/research/gravitational_anomaly.py) | R2C-03 (2026-09-19): tensor integrality for the p2 term, and the colour-cubic (2-form x 6-form) obstruction of M1's declared content. | `tensor_coefficient_witnesses`, `tensor_integrality`, `m1_tensor_survivors`, `colour_cubic_ledger` (+4 more in source) |
| [`mtft.research.involutions`](../src/mtft/research/involutions.py) | Oloid-to-Hodge handoff (Astra, 2026-09-19), exact algebra: one rational involution in three settings. | `oloid_involution`, `triangle_number`, `conjugating_coordinate`, `involution_checks` (+2 more in source) |
| [`mtft.research.mode_operators`](../src/mtft/research/mode_operators.py) | INT-09 (AXG-03 §§6–7, AXG-04 §3): index versus cohomology, purity certificates, the elementary-scalar Bochner bound, | `spin_dirac_index`, `s0_twist_cohomology`, `purity_certificate`, `bochner_bound` (+1 more in source) |
| [`mtft.research.parents`](../src/mtft/research/parents.py) | INT-10: immutable model records — M1, the rejected parent controls, and C3X — with the conventions each carries. | Package exports / data |
| [`mtft.research.pipeline`](../src/mtft/research/pipeline.py) | Composed gate battery (Wave C): field inventory -> anomaly polynomial -> charge lattice -> discrete anomaly -> flux/scalar gate | `m1_gate_report`, `c3x_gate_report`, `route_2_requirements` |
| [`mtft.research.product_surface`](../src/mtft/research/product_surface.py) | R2C-07 (2026-09-19): the product surface S = X0(143) x E (E = 143a1) — index and slope decoupled (exact arithmetic). | `curve_cohomology`, `torus_cohomology`, `product_block`, `triangle` (+28 more in source) |
| [`mtft.research.reflection_positivity`](../src/mtft/research/reflection_positivity.py) | CC-35 (2026-09-29): reflection positivity of the MTFT lattice action — the Polyakov-loop term is NOT covered by the Osterwalder–Seiler | `a_n`, `hook_expansion`, `reduce_to_su`, `conjugate_rep` (+10 more in source) |
| [`mtft.research.repaired_action`](../src/mtft/research/repaired_action.py) | CC-69 (2026-10-17): the character repair of CC-35 adopted — the MTFT plaquette action with the arithmetic tower on Sym^n characters, | `a_n`, `sym_dim`, `WeylGroupIntegrator`, `integrator` (+11 more in source) |
| [`mtft.research.tensor_gs`](../src/mtft/research/tensor_gs.py) | INT-07 (AXG-02 §6, AXG-03 §4, AXG-04 §5): flux/scalar gates, tensor transgression, factorization and integral lattices. | `native_scalar_flux_gate`, `flux_transgression`, `factorization_check`, `c3x_anomaly_polynomial` (+1 more in source) |
| [`mtft.research.theta_torus`](../src/mtft/research/theta_torus.py) | R2C-09 (2026-09-19): the torus factor of the surface Yukawa — theta functions on 143a1. | `tau_143a1`, `theta_basis`, `inner`, `torus_factor` (+7 more in source) |
| [`mtft.research.unified_parent`](../src/mtft/research/unified_parent.py) | R2C-04 (2026-09-19): where an MTFT-native parent can live — exact representation-theoretic facts. | `parent_admissible`, `complex_block_families`, `e7_three_27_families`, `e6_three_16_families` (+3 more in source) |
| [`mtft.research.vacuum_energy`](../src/mtft/research/vacuum_energy.py) | PV-01 (2026-09-29): the flux (Yang–Mills) energy of a split background on S = X0(143) x 143a1 as an exact function of the shape modulus | `intersection`, `line_bundle_energy`, `split_topology`, `flux_energy` (+3 more in source) |
| [`mtft.research.vector_higgs`](../src/mtft/research/vector_higgs.py) | R2C-02 (2026-09-19): the internal-vector Higgs of M1 — mode operator, tachyon mass, and the Yukawa as the gauge interaction. | `tachyon_mass2`, `mass_operator_spectrum`, `m1_vector_higgs_spectrum`, `yukawa_gauge_ratio` (+5 more in source) |
| [`mtft.research.yukawa_triangle`](../src/mtft/research/yukawa_triangle.py) | R2C-06 (2026-09-19): the telescoping theorem for gauge-vertex Yukawas, and the light-Higgs no-go of the class on a curve. | `triangle_closure`, `forced_higgs_hypercharge`, `m1_triangles`, `e7_triangle_test` (+1 more in source) |

## surface

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.surface`](../src/mtft/surface/__init__.py) | mtft.surface — the Modular Surface Laboratory as an mtft subpackage (v0.27.0). | `report`, `all_pass` |
| [`mtft.surface.arithspin`](../src/mtft/surface/arithspin.py) | mtft.surface.arithspin — arithmetic spin structures, cuspidal group, CM fixed points, three-family purity | `divisors`, `eta_quotient_lattice`, `in_lattice`, `cuspidal_group` (+15 more in source) |
| [`mtft.surface.bimodule`](../src/mtft/surface/bimodule.py) | mtft.surface.bimodule — doubled-space (real spectral triple) census on H_1(X0(N), R). | `orthonormal_frame`, `sector_dimensions`, `tensor_sector_dimensions`, `adjoint_identity` (+10 more in source) |
| [`mtft.surface.blowup`](../src/mtft/surface/blowup.py) | mtft.surface.blowup — B3 of the parent-action ledger (PV-02, 2026-09-29): the blow-up of X0(143) x 143a1 at an arithmetic point, its | `dot`, `self_intersection`, `is_spin`, `net_chirality` (+11 more in source) |
| [`mtft.surface.condensation`](../src/mtft/surface/condensation.py) | mtft.surface.condensation — tachyon condensation as bundle extension (COND-01, v0.29.0). | `tachyon_mass2`, `condensation_energy`, `extension_cohomology` |
| [`mtft.surface.crt_dessin`](../src/mtft/surface/crt_dessin.py) | INT-04a (HOPF-02 §§1–3): P^1(Z/143) = P^1(F_11) x P^1(F_13) by CRT and the canonical dessin of X0(143). | `crt`, `canon_prime`, `crt_projective_line`, `act` (+4 more in source) |
| [`mtft.surface.cycles`](../src/mtft/surface/cycles.py) | mtft.surface.cycles — EXACT layer: deterministic integral basis of H_1(X0(N), Z). | `bareiss_det`, `CycleBasis`, `tree_cotree` |
| [`mtft.surface.dynamics`](../src/mtft/surface/dynamics.py) | mtft.surface.dynamics — linear Hamiltonian dynamics on H_1(X0(143), R) with genericity controls. | `AmbiguousClosure`, `Stage`, `hecke_block_projectors`, `stage_from_frozen` (+7 more in source) |
| [`mtft.surface.frozen`](../src/mtft/surface/frozen.py) | mtft.surface.frozen — frozen certified data for X0(143), re-verified at call time without PARI/GP. | `x0143`, `verify_gates`, `periods_frame_gates_from`, `canonical_frame_gates_from` (+2 more in source) |
| [`mtft.surface.gauge`](../src/mtft/surface/gauge.py) | mtft.surface.gauge — gauge theory ON the modular surface (reading A). | `closed_area`, `flux_action`, `flux_partition_sum`, `flat_connection_torus_dimension` (+12 more in source) |
| [`mtft.surface.hodge`](../src/mtft/surface/hodge.py) | mtft.surface.hodge — discrete Hodge layer on the Manin complex. | `unweighted_hodge`, `Mesh`, `base_mesh`, `refine` (+8 more in source) |
| [`mtft.surface.hodge_blocks`](../src/mtft/surface/hodge_blocks.py) | INT-04b (HOPF-02 §§2, 7): coprime block projectors of T2 on H_1(X0(143)) and the quaternionic-multiplicity gate. | `coprime_block_projectors`, `quaternionic_multiplicity_gate` |
| [`mtft.surface.hodge_structure`](../src/mtft/surface/hodge_structure.py) | mtft.surface.hodge_structure — the true Hodge structure on H_1(X0(N), R) in the cycle basis. | `HodgeStructure`, `from_periods`, `gates_pass`, `family_distances` (+4 more in source) |
| [`mtft.surface.hopf_geometry`](../src/mtft/surface/hopf_geometry.py) | INT-03 (HOPF-02): the half-form pencil of X0(143) and the pulled-back quaternionic Hopf connection (exact bookkeeping). | `half_form_pencil`, `hopf_connection_pullback`, `pullback_line_degrees`, `w13_fixed_points_under_pencil` |
| [`mtft.surface.hym`](../src/mtft/surface/hym.py) | mtft.surface.hym — the compact uniformising metric of X0(143), Green's functions and Hermitian–Yang–Mills | `assemble`, `solve_liouville`, `locate`, `solve_neumann` (+17 more in source) |
| [`mtft.surface.intertwiner`](../src/mtft/surface/intertwiner.py) | mtft.surface.intertwiner — the exact map between the surface cycle frame and the canonical | `sector_j`, `canonical_frame_gates`, `cross_frame_hodge_check`, `periods_frame_map` |
| [`mtft.surface.ising`](../src/mtft/surface/ising.py) | mtft.surface.ising — the Ising model on the dual Manin graph: a first QFT on X0(N). | `gf2_rref`, `gf2_solve`, `gf2_nullspace`, `gf2_complement_reps` (+11 more in source) |
| [`mtft.surface.magnetic`](../src/mtft/surface/magnetic.py) | KK-TOWER-01: magnetic Bochner spectra of line bundles O(D) on X0(143) in the compact (Liouville) metric. | `MagneticMesh`, `m1_towers`, `m1_fem_yukawa`, `ratio_distribution` |
| [`mtft.surface.manin`](../src/mtft/surface/manin.py) | mtft.surface.manin — EXACT layer: the Manin complex of X0(N). | `factorize`, `divisors`, `euler_phi`, `kronecker` (+9 more in source) |
| [`mtft.surface.marked`](../src/mtft/surface/marked.py) | mtft.surface.marked — marked two-factor geometry: decomposition, degeneracy maps, readouts (EXACT). | `decomposition`, `degeneracy_maps`, `readout_coefficient`, `cross_map_readout` (+3 more in source) |
| [`mtft.surface.oldsector`](../src/mtft/surface/oldsector.py) | mtft.surface.oldsector — exact two-prime oldform sectors and the Hermitian-part selection rule. | `local_block`, `two_prime_sector`, `connected_part`, `schmidt_rank` (+2 more in source) |
| [`mtft.surface.petersson`](../src/mtft/surface/petersson.py) | mtft.surface.petersson — Petersson norms by hyperbolic quadrature on the glued Manin mesh (PET-01, v0.29.0). | `quadrature`, `reduce_points`, `ext`, `eval_forms` (+1 more in source) |
| [`mtft.surface.rrspace`](../src/mtft/surface/rrspace.py) | mtft.surface.rrspace — Riemann–Roch spaces on X0(143) with poles at CM points (SM-02 / step 5, v0.30.1). | `w13_fixed_representatives`, `w143_fixed_representatives`, `eigenforms_at`, `cm_classes` (+24 more in source) |
| [`mtft.surface.smflux`](../src/mtft/surface/smflux.py) | mtft.surface.smflux — multi-stack flux models on X0(143): search, anomaly ledger, Higgs sectors, purity | `conj`, `canon`, `pair_rep`, `search` (+7 more in source) |
| [`mtft.surface.spectral`](../src/mtft/surface/spectral.py) | mtft.surface.spectral — hyperbolic Laplace spectrum of X0(N) by cusp-truncated FEM (SPEC-01, v0.28.2). | `in_domain`, `face_mesh`, `surface_spectrum`, `assemble_local` (+8 more in source) |
| [`mtft.surface.spin_circle`](../src/mtft/surface/spin_circle.py) | INT-02 (SC7-01, Astra 2026-09-16): the spin circle bundle P = S(S0) -> X0(143) and rank-one local systems on it. | `commutator_word`, `seifert_relators`, `fox_row`, `low_cohomology` (+6 more in source) |
| [`mtft.surface.transport`](../src/mtft/surface/transport.py) | mtft.surface.transport — EXACT Hecke / Atkin-Lehner transport onto the cycle lattice. | `egcd`, `lift_dart`, `cusp_class`, `atkin_lehner_matrix` (+4 more in source) |
| [`mtft.surface.yukawa`](../src/mtft/surface/yukawa.py) | mtft.surface.yukawa — exact Yukawa tensors on X0(143) as canonical-ring multiplication (YUK-01, v0.29.0). | `load_basis`, `al_eigenbasis`, `conv`, `rank_mod` (+2 more in source) |

## thetachar

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.thetachar`](../src/mtft/thetachar/__init__.py) | Theta characteristics / spin structures mod 2 for X0(143). | `gf2_rank`, `gf2_inv`, `gf2_solve_affine`, `standard_J2` (+5 more in source) |
