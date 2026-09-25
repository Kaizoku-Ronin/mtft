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
| [`mtft.__main__`](../src/mtft/__main__.py) | MTFT Command-Line Interface | `cmd_verify`, `cmd_report`, `cmd_tower`, `cmd_screen` (+2 more in source) |
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
| [`mtft.research.product_surface`](../src/mtft/research/product_surface.py) | R2C-07 (2026-09-19): the product surface S = X0(143) x E (E = 143a1) — index and slope decoupled (exact arithmetic). | `curve_cohomology`, `torus_cohomology`, `product_block`, `triangle` (+13 more in source) |
| [`mtft.research.tensor_gs`](../src/mtft/research/tensor_gs.py) | INT-07 (AXG-02 §6, AXG-03 §4, AXG-04 §5): flux/scalar gates, tensor transgression, factorization and integral lattices. | `native_scalar_flux_gate`, `flux_transgression`, `factorization_check`, `c3x_anomaly_polynomial` (+1 more in source) |
| [`mtft.research.theta_torus`](../src/mtft/research/theta_torus.py) | R2C-09 (2026-09-19): the torus factor of the surface Yukawa — theta functions on 143a1. | `tau_143a1`, `theta_basis`, `inner`, `torus_factor` (+6 more in source) |
| [`mtft.research.unified_parent`](../src/mtft/research/unified_parent.py) | R2C-04 (2026-09-19): where an MTFT-native parent can live — exact representation-theoretic facts. | `parent_admissible`, `complex_block_families`, `e7_three_27_families`, `e6_three_16_families` (+3 more in source) |
| [`mtft.research.vector_higgs`](../src/mtft/research/vector_higgs.py) | R2C-02 (2026-09-19): the internal-vector Higgs of M1 — mode operator, tachyon mass, and the Yukawa as the gauge interaction. | `tachyon_mass2`, `mass_operator_spectrum`, `m1_vector_higgs_spectrum`, `yukawa_gauge_ratio` (+5 more in source) |
| [`mtft.research.yukawa_triangle`](../src/mtft/research/yukawa_triangle.py) | R2C-06 (2026-09-19): the telescoping theorem for gauge-vertex Yukawas, and the light-Higgs no-go of the class on a curve. | `triangle_closure`, `forced_higgs_hypercharge`, `m1_triangles`, `e7_triangle_test` (+1 more in source) |

## surface

| Module | Scope | Entry points |
|---|---|---|
| [`mtft.surface`](../src/mtft/surface/__init__.py) | mtft.surface — the Modular Surface Laboratory as an mtft subpackage (v0.27.0). | `report`, `all_pass` |
| [`mtft.surface.arithspin`](../src/mtft/surface/arithspin.py) | mtft.surface.arithspin — arithmetic spin structures, cuspidal group, CM fixed points, three-family purity | `divisors`, `eta_quotient_lattice`, `in_lattice`, `cuspidal_group` (+15 more in source) |
| [`mtft.surface.bimodule`](../src/mtft/surface/bimodule.py) | mtft.surface.bimodule — doubled-space (real spectral triple) census on H_1(X0(N), R). | `orthonormal_frame`, `sector_dimensions`, `tensor_sector_dimensions`, `adjoint_identity` (+10 more in source) |
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
