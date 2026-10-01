"""Matter-03 exact inputs to a conditional classical Higgs-branch theorem.

S=Bl_(infinity,O)(X0(143) x 143a1), L=pi*O(P) tensor O(-3D).
These computations concern the canonical-twisted, defect-free BHV bulk
parent. A surface Hitchin--Kobayashi theorem supplies metric existence;
this module does not compute that metric or select a Higgs/mass scale.
All arithmetic in this module uses the Python standard library.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path

from .collinear import (
    _json_safe, _load_inputs, _rref, cohomology, cusp_jet_certificate,
    yukawa_witness,
)


HD_VECTOR = tuple(map(F, (
    "8235/34426", "0", "-603/17213", "171/9836", "81/9836",
    "593/68852", "0", "0", "0", "0", "0", "0", "0",
)))
HU_VECTOR = tuple(map(F, (
    "27/3083", "0", "-37/6166", "25/6166", "139/6166",
    "0", "0", "0", "0", "0", "0", "0", "0",
)))
POWERS = (0, 0, 0, 1, 1, 2, -3, -1)
Y = (F(1, 6), F(1, 6), F(1, 6), F(0), F(0), F(-1, 2), F(-1, 2), F(1, 2))
T3 = (F(0), F(0), F(0), F(1, 2), F(-1, 2), F(0), F(0), F(0))
Q = tuple(y+t for y, t in zip(Y, T3))
PRIMARY_REFERENCES = (
    {"url": "https://arxiv.org/pdf/1312.2673",
     "locator": "Tanaka, Definitions 1.1-1.2 and Theorem 1.3: End(V) tensor K_S pairs on surfaces"},
    {"url": "https://arxiv.org/pdf/0802.3391",
     "locator": "Beasley-Heckman-Vafa, equations 3.28 and 3.32-3.34: classical bulk F/D equations"},
)


def _rational(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise TypeError(f"{name} must be an exact integer or Fraction")
    return F(value)


def single_arrow_gate(A_E, A_X, epsilon) -> dict:
    """Exact necessary slope gate, including the full spectator constraints.

    On the wall, equality for the invariant weak subbundle excludes every
    nonzero single-Hd solution with any smooth positive Hermitian metric.
    Off the wall, the invariant L^2 and L^-3 spectators already obstruct
    semistability. Changing this conclusion requires changing the parent.
    """
    ae = _rational(A_E, "A_E")
    ax = _rational(A_X, "A_X")
    eps = _rational(epsilon, "epsilon")
    if not (0 < eps < min(ae, ax)):
        raise ValueError("The stated product blow-up polarization must be ample")
    x = ae-3*eps
    spectator_gate = 2*x <= 0 and -3*x <= 0
    return {"mu_L": x, "mu_V": F(0), "A1_slope": 2*x,
            "A2_slope": -3*x, "spectator_necessary_gate": spectator_gate,
            "on_wall": x == 0, "nonzero_single_Hd_D_flat": False,
            "obstruction": "equal-slope invariant weak subbundle" if x == 0
            else "invariant spectator has positive slope",
            "metric_scope": "every smooth positive Hermitian metric"}


@dataclass(frozen=True)
class HiggsPairWitness:
    Hd_basis_coefficients: tuple[F, ...]
    Hu_basis_coefficients: tuple[F, ...]
    Hd_q_coefficients: tuple[F, ...]
    Hu_q_coefficients: tuple[F, ...]
    Hd_curve_order: int
    Hu_curve_order: int
    product_curve_order: int
    product_exceptional_valuation: int
    product_leading_coefficient: F
    Hd_space_dimension: int
    Hu_space_dimension: int
    spectral_product_is_square: bool
    active_rank2_Higgs_pair_stable: bool
    full_rank8_polystable_on_wall: bool
    Hermitian_metric_exists_by_surface_HK: bool
    Hermitian_metric_computed: bool = False
    selected_scale: bool = False
    physical_mass_values_derived: bool = False


@lru_cache(maxsize=1)
def paired_higgs_witness() -> HiggsPairWitness:
    """Freeze a spectral product of odd exceptional valuation seven.

    The sections use the Matter-02 Hd with a4=1,a5=0 and reverse Hu with
    a3=1,a4=0. Since f(z) dz/z has curve order n-1, their orders are3
    and2. The blow-up top-form Jacobian contributes one exceptional zero
    per form, so the canonical-square product has order3+2+2=7.
    A square cannot have odd divisorial valuation; the rank2 spectral
    polynomial is therefore irreducible over the surface function field.
    """
    cusp_jet_certificate()
    inputs, _ = _load_inputs()
    f = inputs["X0_143_S2_qexpansions.txt"]
    down = tuple(sum(F(row[i])*HD_VECTOR[i] for i in range(13)) for row in f)
    up = tuple(sum(F(row[i])*HU_VECTOR[i] for i in range(13)) for row in f)
    if (HD_VECTOR != yukawa_witness().canonical_basis_coefficients
            or down[1:6] != (0, 0, 0, 1, 0)
            or up[1:5] != (0, 0, 1, 0)):
        raise ValueError("Frozen paired Higgs sections failed their exact cusp constraints")
    down_order = next(i-1 for i, value in enumerate(down) if value)
    up_order = next(i-1 for i, value in enumerate(up) if value)
    product = tuple(sum(down[i]*up[n-i] for i in range(n+1)) for n in range(20))
    q_order = next(i for i, value in enumerate(product) if value)
    curve_order = q_order-2
    exceptional_order = curve_order+2
    if (down_order, up_order, curve_order, exceptional_order, product[q_order]) != (3, 2, 5, 7, 1):
        raise ValueError("Canonical product valuation certificate failed")
    if cohomology(2).phi_dim != 10 or cohomology(-2).phi_dim != 11:
        raise ValueError("Pair mode spaces disagree with the exact cohomology census")
    return HiggsPairWitness(HD_VECTOR, HU_VECTOR, down, up, down_order, up_order,
                            curve_order, exceptional_order, product[q_order], 10, 11,
                            False, True, True, True)


def _unit_commutator(i, j, k, l):
    entries = {}
    if j == k:
        entries[(i, l)] = entries.get((i, l), 0)+1
    if l == i:
        entries[(k, j)] = entries.get((k, j), 0)-1
    return {key: value for key, value in entries.items() if value}


def classical_gauge_certificate() -> dict:
    """Constant compact centralizer and neutral root constraints only.

    The result does not include axion/Stueckelberg terms, and does not
    prove that Q is the unique massless Abelian generator of a completion.
    """
    allowed = tuple((i, j) for i in range(8) for j in range(8) if POWERS[i] == POWERS[j])
    roots = ((3, 7), (7, 3))
    rows = [[_unit_commutator(i, j, k, l).get((a, b), 0) for i, j in allowed]
            for k, l in roots for a in range(8) for b in range(8)]
    rank = len(_rref(rows)[1])
    if (len(allowed), rank, sum(Q), Q[3]-Q[7]) != (16, 3, 0, 0):
        raise ValueError("Classical paired-root gauge constraints failed")
    return {"U8_centralizer": "U(3) x U(1)^4", "dimension": 16-rank,
            "SU8_centralizer": "S(U(3) x U(1)^4)", "dimension_SU8": 15-rank,
            "ordered_indices": ("c1", "c2", "c3", "w1", "w2", "A1", "A2", "B"),
            "Q": Q, "Q_commutes_with_both_arrows": True,
            "unique_photon": False, "axionic_photon_kernel": "OPEN"}


def report() -> dict:
    """JSON-safe exact inputs and explicit theorem/physics boundaries."""
    pair = paired_higgs_witness()
    return _json_safe({
        "milestone": "Matter-03", "status": "EXACT CONDITIONAL CLASSICAL EXISTENCE",
        "parent": "defect-free classical BHV bulk canonical twist on the fixed Matter-02 surface",
        "single_Hd_gate": single_arrow_gate(3, 3, 1),
        "paired_higgs": asdict(pair), "classical_gauge": classical_gauge_certificate(),
        "branch": {"common_scaling_parameter": "t in C*", "t0_split_pair_polystable": True,
                   "selected_scale": False, "Higgs_modulus_stabilized": False},
        "existence_method": "irreducible generic rank2 spectral polynomial, wall polystability, surface HK theorem",
        "primary_references": PRIMARY_REFERENCES,
        "assumptions": ["Archived columns are the full S2(Gamma0(143)) cusp basis.",
                        "The specified smooth compact projective holomorphic surface and line bundle, without flat twists.",
                        "An ample slope-zero wall polarization, such as A_E=A_X=3, epsilon=1.",
                        "The canonical-twisted surface Hitchin-Kobayashi correspondence.",
                        "No added defects, soft terms or quantum/axion potential."],
        "open": ["Higgs modulus and light-family selection", "other charged and neutral deformed sectors",
                 "numeric Hermitian D-term metric", "positive L2 kinetic normalization",
                 "physical masses and axionic photon kernel"],
        "input_hashes": cusp_jet_certificate()["input_hashes"],
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = json.dumps(report(), indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding='utf-8')
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
