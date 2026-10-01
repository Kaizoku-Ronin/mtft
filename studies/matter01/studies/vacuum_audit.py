"""Matter-01 B3 vacuum-to-electromagnetism audit.

Rechecks one frozen, finite-box blow-up witness; does not repeat a full search.
Numerical parameters of the separate QED reference engine are not MTFT outputs.
Run: python studies/vacuum_audit.py [--repo path/to/mtft] [--output path.json]
Requires sympy, as does mtft.surface.blowup. Uses exact integer/rational arithmetic.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path


ORDER = ("c", "L", "A1", "A2", "B")
RANKS = {"c": 3, "L": 2, "A1": 1, "A2": 1, "B": 1}
STACKS = {
    "c": (0, 0, 0), "L": (-3, 1, 21),
    "A1": (-1, 0, 8), "A2": (4, -1, -29), "B": (3, -1, -21),
}
NORMAL = (8, 3, 1)  # (A_E, A_X, epsilon), arbitrary common area unit
HYPERCHARGE = {"c": F(1, 6), "L": F(0), "A1": F(-1, 2),
               "A2": F(-1, 2), "B": F(1, 2)}


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def jsonable(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [jsonable(v) for v in x]
    if isinstance(x, (str, int, float, bool)) or x is None:
        return x
    return str(x)


def audit(repo: Path):
    B = load_module(repo / "src/mtft/surface/blowup.py", "matter01_blowup")
    A8 = load_module(repo / "src/mtft/research/anomaly8d.py", "matter01_anomaly8d")
    checks = []

    def check(name, passed, witness):
        checks.append({"name": name, "pass": bool(passed), "status": "EXACT",
                       "witness": jsonable(witness)})

    def independent_phi(d):
        return d[2] - 24 * d[1]

    def independent_dot(d, q):
        return d[0] * q[1] + q[0] * d[1] - d[2] * q[2]

    def difference(i, j):
        return tuple(a - b for a, b in zip(STACKS[i], STACKS[j]))

    AE, AX, eps = NORMAL
    phis = {k: independent_phi(d) for k, d in STACKS.items()}
    slopes = {k: d[0] * AE + d[1] * AX + d[2] * eps for k, d in STACKS.items()}
    check("ampleness", 0 < eps < min(AE, AX) and B.ample(AX, AE, eps), NORMAL)
    check("all_stack_slopes_zero", all(v == 0 for v in slopes.values()), slopes)
    for name, d in STACKS.items():
        check(f"phi_{name}_vs_release", phis[name] == B.phi(d) == B.net_chirality(d), phis[name])
        check(f"slope_{name}_vs_release", slopes[name] == B.slope(d, AX, AE, eps), slopes[name])
    check("five_distinct_stack_classes", len(set(STACKS.values())) == 5, STACKS)
    check("noncollinear_witness", not B._collinear(list(STACKS.values())), STACKS)
    trace_flux = tuple(sum(RANKS[k] * STACKS[k][j] for k in ORDER) for j in range(3))
    check("rank_weighted_flux_traceless", trace_flux == (0, 0, 0), trace_flux)

    block_counts = {"Q": [phis["c"] - phis["L"]],
                    "u^c": [phis[a] - phis["c"] for a in ("A1", "A2")],
                    "d^c": [phis["B"] - phis["c"]],
                    "L": [phis[a] - phis["L"] for a in ("A1", "A2")] +
                         [-(phis["B"] - phis["L"])],
                    "e^c": [phis["B"] - phis[a] for a in ("A1", "A2")]}
    net = {k: sum(v) for k, v in block_counts.items()}
    pairs = {k: (sum(abs(n) for n in v) - abs(sum(v))) // 2 for k, v in block_counts.items()}
    sterile_lower_bound = abs(phis["A1"] - phis["A2"])
    release = B.sm_phi_conditions(phis)
    check("independent_net_content_vs_release", net == release["net"], net)
    check("three_net_SM_families", all(n == 3 for n in net.values()), net)
    lower_bounds = {"u^c_pairs": pairs["u^c"], "e^c_pairs": pairs["e^c"],
                    "doublet_pairs": pairs["L"], "neutral_singlets": sterile_lower_bound}
    expected = {"u^c_pairs": 5, "e^c_pairs": 5, "doublet_pairs": 8, "neutral_singlets": 13}
    check("index_forced_lower_bounds", lower_bounds == expected, lower_bounds)
    check("archived_extras_score", sum(lower_bounds.values()) == 31, sum(lower_bounds.values()))

    # RR provides useful bounds without claiming a computation of h0/h1/h2.
    # chi = h0-h1+h2 and hq>=0 imply h1 >= max(0,-chi).
    # Classes include no gamma component; intersection here is independent of B.chi.
    rr_bounds = []
    for i, left in enumerate(ORDER):
        for right in ORDER[i + 1:]:
            d = difference(left, right)
            dual = tuple(-v for v in d)
            self2 = independent_dot(d, d)
            kd = 24 * d[1] - d[2]
            ci, cid = F(self2 - kd, 2), F(self2 + kd, 2)
            bound, bound_dual = max(0, -ci), max(0, -cid)
            row = {"orientation": f"{left}/{right}", "class_(a,b,e)": d,
                   "chi_L": ci, "chi_dual": cid, "net_chi_difference": ci-cid,
                   "h1_L_lower_bound": bound, "h1_dual_lower_bound": bound_dual,
                   "status": "EXACT lower bounds for cohomology; not full dimensions or physical masses"}
            rr_bounds.append(row)
            check(f"RR_bound_{left}_{right}", ci == B.chi(d) and cid == B.chi(dual)
                  and ci-cid == independent_phi(d), row)

    # Traceless U(8) representative of the root-charge hypercharge assignment.
    trY = sum(RANKS[k] * HYPERCHARGE[k] for k in ORDER)
    trY2 = sum(RANKS[k] * HYPERCHARGE[k] ** 2 for k in ORDER)
    trT32 = F(1, 2)
    check("hypercharge_traceless", trY == 0, trY)
    check("hypercharge_generator_norm", trY2 == F(5, 6), trY2)
    shift = B.sp.Symbol("s", real=True)
    shift_norm = B.sp.expand(sum(RANKS[k] * (B.sp.Rational(HYPERCHARGE[k].numerator,
                    HYPERCHARGE[k].denominator) + shift) ** 2 for k in ORDER))
    check("common_center_hypercharge_lift_norm", shift_norm == B.sp.Rational(5,6)+8*shift**2,
          str(shift_norm))
    charges = {"Q": HYPERCHARGE["c"] - HYPERCHARGE["L"],
               "u^c": HYPERCHARGE["A1"] - HYPERCHARGE["c"],
               "d^c": HYPERCHARGE["B"] - HYPERCHARGE["c"],
               "L": HYPERCHARGE["A1"] - HYPERCHARGE["L"],
               "e^c": HYPERCHARGE["B"] - HYPERCHARGE["A1"]}
    check("SM_hypercharge_root_assignments", charges == {"Q": F(1, 6), "u^c": F(-2, 3),
          "d^c": F(1, 3), "L": F(-1, 2), "e^c": F(1)}, charges)
    charges_after_ew = {"L_doublet": [F(1, 2) + charges["L"], F(-1, 2) + charges["L"]],
                       "e^c": charges["e^c"],
                       "H_u": [F(1, 2) + F(1, 2), F(-1, 2) + F(1, 2)],
                       "H_d": [F(1, 2) - F(1, 2), F(-1, 2) - F(1, 2)]}
    check("electric_charge_and_neutral_Higgs_components", charges_after_ew == {
          "L_doublet": [0, -1], "e^c": 1, "H_u": [1, 0], "H_d": [0, -1]}, charges_after_ew)
    triangle_checks = []
    for a in ("A1", "A2"):
        ds = [difference(a, "L"), difference("B", a), difference("L", "B")]
        closure = tuple(sum(d[j] for d in ds) for j in range(3))
        triangle_checks.append({"stack": a, "classes": ds, "sum": closure,
                                "Y_sum": charges["L"] + charges["e^c"] - F(1, 2)})
    check("charged_lepton_Yukawa_topological_triangle_only", all(t["sum"] == (0, 0, 0)
          and t["Y_sum"] == 0 for t in triangle_checks), triangle_checks)

    # These two topological projections are NOT a specified axion mass matrix.
    trYF = tuple(sum(RANKS[k] * HYPERCHARGE[k] * STACKS[k][j] for k in ORDER) for j in range(3))
    trYF2 = sum(RANKS[k] * HYPERCHARGE[k] * independent_dot(STACKS[k], STACKS[k]) for k in ORDER)
    check("linear_hypercharge_flux_projection", trYF == (0, 0, 0), trYF)
    check("quadratic_hypercharge_flux_projection", trYF2 == 233, trYF2)

    # Necessary anomalies of the NET low-energy spectrum, not parent completion.
    rows = [("Q", 3, 2, F(1, 6), 1), ("u^c", 3, 1, F(-2, 3), -1),
            ("d^c", 3, 1, F(1, 3), -1), ("L", 1, 2, F(-1, 2), 0),
            ("e^c", 1, 1, F(1), 0)]
    anomalies = {
        "gravitational_U1Y": sum(net[k] * nc * nw * y for k, nc, nw, y, _ in rows),
        "U1Y_cubed": sum(net[k] * nc * nw * y ** 3 for k, nc, nw, y, _ in rows),
        "SU3_squared_U1Y": sum(net[k] * nw * F(1, 2) * y for k, nc, nw, y, _ in rows if nc == 3),
        "SU2_squared_U1Y": sum(net[k] * nc * F(1, 2) * y for k, nc, nw, y, _ in rows if nw == 2),
        "SU3_cubed": sum(net[k] * nw * cubic for k, nc, nw, y, cubic in rows),
        "net_SU2_doublet_multiplicity": sum(net[k] * nc for k, nc, nw, y, _ in rows if nw == 2),
    }
    check("4d_net_perturbative_anomaly_cancellation", all(v == 0 for k, v in anomalies.items()
          if k != "net_SU2_doublet_multiplicity"), anomalies)
    check("4d_net_SU2_Witten_parity", anomalies["net_SU2_doublet_multiplicity"] % 2 == 0,
          anomalies["net_SU2_doublet_multiplicity"])
    I10 = A8.I10([(A8.adjoint("A", 8), 1)])
    check("8d_adjoint_local_I10_zero", I10 == 0, str(I10))
    conditional_gY2_over_g2 = trT32 / trY2
    conditional_sin2 = conditional_gY2_over_g2 / (1 + conditional_gY2_over_g2)
    check("conditional_universal_YM_coupling_normalisation", conditional_gY2_over_g2 == F(3, 5)
          and conditional_sin2 == F(3, 8), {"gY_squared_over_g2_squared": conditional_gY2_over_g2,
          "sin2_theta_at_matching": conditional_sin2})

    # A comparison control exposes why minimizing the index score is not enough.
    fam = B.collinear_family()
    alt_stacks = {k: tuple(fam["stacks"][k]) for k in ORDER}
    alt_phis = {k: independent_phi(v) for k, v in alt_stacks.items()}
    alt_rec = B.sm_phi_conditions(alt_phis)
    alt_rr = []
    alt_ok = True
    for i, left in enumerate(ORDER):
        for right in ORDER[i+1:]:
            d = tuple(a-b for a, b in zip(alt_stacks[left], alt_stacks[right]))
            sq, kd = independent_dot(d, d), 24*d[1]-d[2]
            ci, cid = F(sq-kd, 2), F(sq+kd, 2)
            alt_ok &= ci == B.chi(d) and cid == B.chi(tuple(-x for x in d))
            alt_rr.append({"orientation": f"{left}/{right}", "class_(a,b,e)": d,
                           "chi_L": ci, "chi_dual": cid,
                           "h1_L_lower_bound": max(0,-ci), "h1_dual_lower_bound": max(0,-cid)})
    alt_bound_sum = sum(r["h1_L_lower_bound"]+r["h1_dual_lower_bound"] for r in alt_rr)
    primary_bound_sum = sum(r["h1_L_lower_bound"]+r["h1_dual_lower_bound"] for r in rr_bounds)
    check("collinear_control_net_families", all(v == 3 for v in alt_rec["net"].values()), alt_rec["net"])
    check("collinear_control_Dflat_at_3_3_1", all(B.slope(d, 3, 3, 1) == 0 for d in alt_stacks.values())
          and B.ample(3,3,1), alt_stacks)
    check("collinear_control_RR_bounds", alt_ok, alt_rr)
    check("collinear_control_geometric_count_proxy", alt_bound_sum == 666 and primary_bound_sum == 8588,
          {"collinear_H1_bound_sum": alt_bound_sum, "primary_H1_bound_sum": primary_bound_sum})

    paths = ("src/mtft/surface/blowup.py", "src/mtft/research/anomaly8d.py", "src/mtft/parent.py",
             "src/mtft/legend.py", "docs/SM/PV02_BLOWUP_SUSY_REGISTER.md",
             "docs/SM/PARENT_ACTION_REQUIREMENTS.md", "docs/SM/CC26_NORMALISATION_RETRACTION.md")
    source_hashes = {p: hashlib.sha256((repo / p).read_bytes()).hexdigest() for p in paths}
    open_gates = [
        {"gate": "arithmetic_blowup_point_and_holomorphic_bundles", "status": "OPEN", "reason": "An actual point, line bundles including Picard data, and a compatible holomorphic metric are not fixed by NS classes."},
        {"gate": "full_cohomology_and_particle_spectrum", "status": "OPEN", "reason": "No h0, h1, h2 or normalized zero-mode bases on this blow-up have been computed."},
        {"gate": "nonzero_normalized_electron_Yukawa_and_mass", "status": "OPEN", "reason": "Triangle class/charge closure is necessary but does not compute overlaps, kinetic Grams, a selected VEV, or the mass matrix."},
        {"gate": "vectorlike_exotic_lifting", "status": "OPEN", "reason": "Index-forced pairs have not been assigned a full-rank gauge-invariant mass matrix; further pairs may exist."},
        {"gate": "hypercharge_Stueckelberg_kernel", "status": "OPEN", "reason": "No axion/GS charge matrix is specified. Tr(Yc1)=0 does not establish every relevant coupling vanishes; Tr(Yc1^2)=233 is nonzero."},
        {"gate": "physical_gauge_kinetic_matrix_and_mixing", "status": "OPEN", "reason": "Tr(Y^2)=5/6 is a generator norm. The universal parent kinetic action, volume, thresholds, extra-U1 mixing and RG matching are not derived."},
        {"gate": "selected_electroweak_vacuum_and_massless_photon", "status": "OPEN", "reason": "Neutral Higgs components exist algebraically; an actual stable VEV and all axion/gauge mass kernels remain to be demonstrated."},
        {"gate": "global_parent_consistency", "status": "OPEN", "reason": "Local I10 cancellation does not supply global anomalies, charge/flux quantisation and a UV/global completion of the non-spin twisted parent."},
        {"gate": "complete_stability_and_moduli", "status": "OPEN", "reason": "The exact split-slope/D-flat check is conditional internal HYM data, not a computed full fluctuation Hessian or moduli-stabilized spacetime vacuum."},
        {"gate": "dimensionful_scales", "status": "OPEN", "reason": "Kähler ratios are dimensionless; compactification, gauge and Higgs/electron scales are not determined."},
        {"gate": "B3_C3X_common_parent", "status": "OPEN", "reason": "B3 is the 8D blow-up adjoint/twisted branch. ADS-03 evolves a separate 6D C3X-derived classical scalar benchmark; there is no common reduction yet."},
    ]
    return jsonable({
        "study": "Matter-01 selected B3 witness and U(1) bridge", "version": 1,
        "source_repository": "mtft 0.34.0", "source_hashes": source_hashes,
        "scope": "Exact index, slope, charge and necessary-anomaly audit. No claim of a complete physical vacuum or derived QED.",
        "selection": {"method": "Frozen first minimum-score witness of the archived restricted box, independently rechecked without a new full scan.",
                      "archived_n_solutions": 116, "archived_min_extras_score": 31,
                      "archived_boxes": {"n_box": 8, "box_a": 4, "box_b": 1, "box_e": 30},
                      "archived_counts_status": "EXACT in recorded finite box; not re-enumerated in Matter-01; provenance PV02 and v0340_review parent audit",
                      "normal_slice_selection_check": "Four solutions at normal (8,3,1); minimum31 reproduced during this audit, not a global minimum claim."},
        "witness": {"gauge_parent": "conditional topologically twisted 8D U(8) adjoint gauge parent",
                    "stack_order": ORDER, "ranks": RANKS, "classes_(a,b,e)": STACKS,
                    "kahler_normal_(A_E,A_X,epsilon)": NORMAL,
                    "r_A_X_over_A_E": F(AX, AE), "epsilon_over_A_E": F(eps, AE),
                    "phi": phis, "slopes": slopes, "rank_weighted_c1": trace_flux,
                    "net_SM_indices": net, "signed_index_blocks": block_counts,
                    "index_forced_lower_bounds": lower_bounds,
                    "RR_cohomology_lower_bounds": rr_bounds,
                    "RR_warning": "Negative Euler characteristics force much larger H1 dimensions than the abs(net)-based extras score suggests. E.g. c/L has H1>=222 and dual>=225; L/B has H1>=897 and dual>=891. Their lifting and physical interpretation must be handled by a complete reduction.",
                    "extras_score_warning": "31 adds representation-pair counts to neutral-singlet multiplicity; it is neither a particle count nor the full massless spectrum."},
    "U1_bridge": {"convention": "Q=T3+Y; L has Y=-1/2; e^c is left-handed anti-electron with Y=+1",
                      "stack_Y": HYPERCHARGE, "root_Y": charges, "Tr_Y_squared": trY2,
                      "common_center_lift": "Adding s to all five stack hypercharges leaves every adjoint root charge unchanged; Tr(Y(s)^2)=5/6+8s^2. The displayed s=0 representative is traceless. A full GS/kinetic sector must decide the physical lift.",
                      "electric_charges_after_EW": charges_after_ew,
                      "lepton_Yukawa_triangles": triangle_checks,
                      "Tr_Y_c1": trYF, "integral_Tr_Y_c1_squared": trYF2,
                      "conditional_universal_YM": {"assumption": "One universal U8 parent kinetic term, no threshold or kinetic mixing correction, generator tr(TaTb)=delta_ab/2.",
                          "gY_squared_over_g2_squared": conditional_gY2_over_g2,
                          "sin2_theta_matching": conditional_sin2,
                          "status": "EXACT conditional group normalization; not a measured low-energy coupling prediction"},
                      "photon_conditions": ["All axion mass rows annihilate the hypercharge direction before EWSB.",
                          "Selected Higgs VEV is neutral under Q=T3+Y.", "Positive kinetic matrix is canonically normalized including extra U1 mixing.",
                          "The combined gauge mass matrix has the electromagnetic null direction."],
                      "standard_no_extra_mixing_relation": "1/e^2 = 1/g2^2 + 1/gY^2", "derived_electron_mass": None,
                      "derived_alpha_EM": None, "physical_QED_parameters": "INPUT in the separate reference benchmark"},
        "anomalies": {"4d_net": anomalies, "8d_local_adjoint_I10": str(I10),
                      "status": "EXACT necessary local/net checks; globally consistent parent and extra-U1 sectors OPEN"},
        "collinear_comparison_control": {"classes_(a,b,e)": alt_stacks, "phi": alt_phis,
                  "ample_Dflat_locus": "epsilon/A_E=1/3, A_X/A_E>1/3", "net_SM_indices": alt_rec["net"],
                  "index_forced_lower_bounds": alt_rec["vector_like_pairs"], "index_extras_score": 36,
                  "RR_cohomology_lower_bounds": alt_rr,
                  "unweighted_H1_lower_bound_sum_all_Hom_orientations": alt_bound_sum,
                  "primary_unweighted_H1_lower_bound_sum": primary_bound_sum,
                  "status": "EXACT comparison; neither sum is a particle count, and this control is not established viable or optimal.",
                  "recommendation": "Use this lower-degree one-line-bundle family as a cohomology control before preferring the minimum-net-index-extras witness."},
        "checks": checks, "summary": {"n_exact_checks": len(checks), "n_failed": sum(not c["pass"] for c in checks),
                  "complete_vacuum_established": False, "MTFT_derived_electron_photon_sector": False},
        "open_gates": open_gates,
        "physics_references": [
            {"url": "https://pdg.lbl.gov/2025/reviews/standard_model_and_related.html", "use": "Electroweak charge and coupling conventions."},
            {"url": "https://cds.cern.ch/record/819632/files/0502010.pdf", "use": "A. Pich, Standard Model of Electroweak Interactions: Q=T3+Y and e=g sin(theta)=g' cos(theta)."},
            {"url": "https://arxiv.org/abs/hep-th/0304245", "use": "Ruegg/Ruiz-Altaba, The Stueckelberg Field; a gauge symmetry alone does not exclude a Stueckelberg photon mass."},
            {"url": "https://arxiv.org/abs/1411.0034", "use": "Heterotic example of model-dependent axion/hypercharge constraints; not assumed to be the MTFT parent."}],
    })


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    default_repo = root / "repo_candidate" if (root / "repo_candidate").exists() else root.parents[1]
    parser.add_argument("--repo", type=Path, default=default_repo)
    parser.add_argument("--output", type=Path, default=root / "results/vacuum_audit.json")
    args = parser.parse_args()
    result = audit(args.repo.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["summary"], indent=2))
    if result["summary"]["n_failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
