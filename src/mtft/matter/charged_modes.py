"""Matter-03 charged holomorphic Higgs-deformation complex.

Counts refer to the fixed Matter-02 blow-up representative and BHV bulk
dictionary. They are cohomology/hypercohomology dimensions, not observed
particle numbers or physical masses. The paired neutral Higgs scales
are assumed nonzero; no a-field expectation value is present.
"""
from __future__ import annotations
from fractions import Fraction

from .collinear import cohomology

_COMPONENTS = (("c1", 0, Fraction(1, 6)),
               ("c2", 0, Fraction(1, 6)),
               ("c3", 0, Fraction(1, 6)),
               ("L+", 1, Fraction(1, 2)),
               ("L-", 1, Fraction(-1, 2)),
               ("A1", 2, Fraction(-1, 2)),
               ("A2", -3, Fraction(-1, 2)),
               ("B", -1, Fraction(1, 2)))


def split_charged_census() -> dict:
    """Complete |Q|=1 split a/Phi/gauge census, including weak q0 roots."""
    rows = {"Qminus": [], "Qplus": []}
    for i, ni, qi in _COMPONENTS:
        for j, nj, qj in _COMPONENTS:
            charge = qi - qj
            if charge not in (-1, 1):
                continue
            cc = cohomology(ni - nj)
            key = "Qminus" if charge == -1 else "Qplus"
            rows[key].append({"target": i, "source": j, "exponent": ni - nj,
                              "electric_charge": int(charge), "a_dim": cc.h1,
                              "phi_dim": cc.phi_dim, "gauge_dim": cc.h0})
    totals = {key: {name: sum(row[name] for row in table)
                    for name in ("a_dim", "phi_dim", "gauge_dim")}
              for key, table in rows.items()}
    for total in totals.values():
        total["chiral_dim"] = total["a_dim"] + total["phi_dim"]
    return {"status": "EXACT CONDITIONAL SPLIT COHOMOLOGY",
            "rows": rows, "totals": totals,
            "physical_particle_counts": False,
            "note": "q0 gauge generators are separate from a and Phi chiral classes."}


def _six_blocks() -> list[dict]:
    # Negative-charge root orientation i<-L+, i<-B for three spectators.
    blocks = []
    for spectator, exponent in (("L-", 1), ("A1", 2), ("A2", -3)):
        for entry, shift in (("down", 2), ("up", -2)):
            q = exponent - 1 if entry == "down" else exponent + 1
            partner = -q - shift
            shape = (cohomology(q).h1, cohomology(partner).h1)
            blocks.append({"spectator": spectator, "higgs_entry": entry,
                           "negative_a_exponent": q,
                           "positive_a_exponent": partner,
                           "higgs_exponent": shift,
                           "shape": shape, "rank_capacity": min(shape)})
    return blocks


def paired_charged_report() -> dict:
    """A bound valid for any nonzero neutral down/up pair in one weak direction.

    The aa map has six independent root blocks. Their capacities bound
    its rank by118. Each charged weak gauge generator removes one Phi
    class. Therefore the two-term holomorphic complex has at least160
    degree-one classes per charge sign. This is a conditional geometric
    obstruction to light-family selection in the restricted parent.
    """
    census = split_charged_census()
    blocks = _six_blocks()
    capacity = sum(row["rank_capacity"] for row in blocks)
    lower = {}
    for sign, total in census["totals"].items():
        a_kernel_min = total["a_dim"] - capacity
        phi_coker = total["phi_dim"] - 1
        lower[sign] = {"a_kernel_lower_bound": a_kernel_min,
                       "phi_cokernel_dimension": phi_coker,
                       "hyperH1_lower_bound": a_kernel_min + phi_coker}
    return {"status": "EXACT CONDITIONAL LOWER BOUND",
            "split_census": census, "six_a_mass_blocks": blocks,
            "a_map_rank_upper_bound": capacity,
            "charged_gauge_map_rank": {"Qminus": 1, "Qplus": 1},
            "paired_lower_bounds": lower,
            "assumptions": ["Fixed Matter-02 cusp representative and BHV bulk dictionary",
                            "Nonzero Phi entries L+<-B and B<-L+; no a expectation value",
                            "End(V) --adPhi--> K_S tensor End(V) deformation complex"],
            "physical_mass_spectrum": False,
            "open": ["Hermitian BPS realization", "physical kinetic metrics",
                     "additional lifting interactions", "observed light-family selection"]}


def certified_paired_charged_report() -> dict:
    """Exact degree-one counts for canonical orders down3/up2, not masses.

    Projecting S to C gives a Leray curve subspace and elliptic/torsion
    quotient. The down quotient map is injective; the up curve map is
    surjective. These conditions control off-diagonal extension terms
    and give six aa ranks2,2,5,3,17,12, total41. Independent proof review
    is archived in CHARGED_COMPLEX_DERIVATION.md.

    Both Higgs scales are nonzero and the a background is zero. This
    profile is a chosen representative, not a derived electroweak scale.
    """
    result = paired_charged_report()
    ranks = {("L-", "down"): 2, ("L-", "up"): 2,
             ("A1", "down"): 5, ("A1", "up"): 3,
             ("A2", "down"): 17, ("A2", "up"): 12}
    for block in result["six_a_mass_blocks"]:
        block["certified_rank"] = ranks[block["spectator"], block["higgs_entry"]]
    total_rank = sum(ranks.values())
    result["status"] = "EXACT CONDITIONAL HOLOMORPHIC deformation complex"
    result["canonical_vanishing_orders"] = {"down": 3, "up": 2}
    result["a_map_rank"] = total_rank
    result["paired_hyperH1"] = {
        sign: {"a_kernel_dimension": total["a_dim"] - total_rank,
               "phi_cokernel_dimension": total["phi_dim"] - 1,
               "hyperH1_dimension": total["a_dim"] - total_rank + total["phi_dim"] - 1}
        for sign, total in result["split_census"]["totals"].items()}
    result["assumptions"].append("Chosen canonical vanishing orders down3/up2 at the ordinary cusp")
    result["hermitian_bps_existence"] = "Established separately for the stated stable paired branch under its parent assumptions."
    result["open"] = ["actual Hermitian metric and physical kinetic integrals",
                      "additional lifting interactions", "observed light-family selection"]
    return result


def main():
    """Print or save the certified conditional charged-complex report."""
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    payload = json.dumps(certified_paired_charged_report(), indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding='utf-8')
        print('Conditional paired charged-complex count: 237 per sign; physical masses OPEN.')
    else:
        print(payload, end='')


if __name__ == '__main__':
    main()
