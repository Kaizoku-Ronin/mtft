"""Make portable display evidence from the Matter-01 archived calculations.

The browser formula is checked against the independent Python angular oracle.
The displayed residuals are finite-grid numerical diagnostics, not global bounds.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def compose(results: Path, engine_file: Path | None = None) -> dict:
    oracle = json.loads((results / "qed_oracle_checks.json").read_text())
    vacuum = json.loads((results / "vacuum_audit.json").read_text())
    reduction = json.loads((results / "u1_reduction_checks.json").read_text())
    inputs = oracle["inputs"]
    rows = []
    engine = json.loads(engine_file.read_text()) if engine_file else None
    if engine:
        gates = engine["verification"]["gates"]
        rows.append({"label": "Full-sum Ward check · normalized residual", "value": gates["ward_full_sum"]["measured_error"]})
        rows.append({"label": "Polarization gauge shift · normalized residual", "value": gates["gauge_shift_invariance"]["measured_error"]})
        rows.append({"label": "Single diagram in Ward test · nonzero amplitude", "value": gates["ward_requires_two_diagrams"]["individual_diagram_amplitude"]})
    independent = next(g for g in oracle["gates"] if g["name"] == "explicit spinor Msq versus independent oracle")
    rows.append({"label": "Independent |𝓜|² comparison · max relative error", "value": independent["maximum_relative_error"]})
    integration = next(g for g in oracle["gates"] if g["name"] == "engine integrated event total versus oracle")
    rows.append({"label": "Integrated event total · max relative error", "value": integration["maximum_relative_error"]})
    notes = "Ward and gauge-shift residuals are archived diagnostics at √s = 10 MeV, cos θ = 0.37, φ = 0.29 over all 16 spin/polarization choices; they do not change with the sliders. The independent |𝓜|² comparison covers 126 configurations. These are finite-sample numerical checks, not a rigorous global bound."
    if not engine:
        notes += " The engine Ward evidence has not yet been embedded in this internal build."
    return {
        "study": "Matter-01 U(1) lab display evidence",
        "constants": {
            "alpha": 1 / inputs["alpha_inverse"],
            "electron_mass_mev": 1000 * inputs["electron_mass_gev"],
            "mev_minus2_to_nb": 1e6 * inputs["gev2_to_nb"],
            "provenance": inputs["provenance"],
        },
        "angular_oracle": oracle["angular_grid"],
        "total_oracle": oracle["total_grid"],
        "verification": {"rows": rows, "notes": notes},
        "spectrum_status": [
            {"item": "Audited B3 flux-class witness", "status": "EXACT conditional", "supported": True, "detail": "Zero stack slopes and net indices of three for Q, uᶜ, dᶜ, L, and eᶜ. The marked point and holomorphic representatives remain open; net indices are not the complete massless spectrum."},
            {"item": "Electric-charge embedding", "status": "EXACT conditional", "supported": True, "detail": "Q = T₃ + Y assigns neutral/neutrino and charge −1/electron components to the lepton doublet. A charge assignment does not establish a massless physical photon."},
            {"item": "Extra matter and Yukawas", "status": "OPEN", "supported": False, "detail": "Actual cohomology, nonzero overlap tensors, extra-matter masses, and the normalized electron mass remain to be calculated."},
            {"item": "Massless photon in the MTFT vacuum", "status": "OPEN", "supported": False, "detail": "Axion/Stückelberg terms, extra-U(1) mixing, and the Higgs mass matrix must retain the electromagnetic null direction."},
            {"item": "Electron mass and α in this lab", "status": "INPUT", "supported": False, "detail": "Measured PDG inputs. The reference electron/positron/photon fields are defined by the QED action; their spectrum is not derived from MTFT here."},
            {"item": "Gauge–breathing reduction", "status": "EXACT conditional", "supported": True, "detail": "C3X has Z(χ) = exp(χ/Mₚ). This conditional reduction and the B3 vacuum audit are separate constructions; no common parent or coupled gauge evolution is established."},
        ],
        "provenance": {
            "oracle_source": "results/qed_oracle_checks.json",
            "vacuum_source": "results/vacuum_audit.json",
            "reduction_source": "results/u1_reduction_checks.json",
            "engine_source": "results/" + engine_file.name if engine_file else None,
            "witness_stack_order": vacuum["witness"]["stack_order"],
            "reduction_open": reduction["open"],
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path(__file__).resolve().parents[1] / "results")
    parser.add_argument("--engine", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "lab_evidence.json")
    args = parser.parse_args()
    evidence = compose(args.results, args.engine)
    args.output.write_text(json.dumps(evidence, indent=2) + "\n")
    print(args.output)
