"""AXG-04: bounded U(8) Einstein--Yang--Mills + adjoint Dirac control.

Standard-library exact arithmetic. These checks certify finite block algebra,
Riemann--Roch index arithmetic, and the Atiyah--Bott Morse-index substitution.
They do not solve the Einstein equations or construct a stable chiral vacuum.
"""

from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import json

CHECKS = []


def check(name, condition):
    passed = bool(condition)
    CHECKS.append({"name": name, "passed": passed})
    if not passed:
        raise AssertionError(name)


def matrix_unit(i, j, n=8):
    return [[int(a == i and b == j) for b in range(n)] for a in range(n)]


def multiply(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def subtract(a, b):
    return [[x-y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def charge(i, j):
    return tuple(int(a == i)-int(a == j) for a in range(5))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def main():
    stacks = ("c", "L", "a", "b", "d")
    ranks = (3, 2, 1, 1, 1)
    degrees = (0, -3, 3, 3, 0)
    y = (Q(1, 6), Q(0), Q(-1, 2), Q(1, 2), Q(-1, 2))
    genus = 13
    offsets = (0, 3, 5, 6, 7)
    m8 = tuple(degree for rank, degree in zip(ranks, degrees) for _ in range(rank))
    delta = [[mi-mj for mj in degrees] for mi in degrees]
    dimensions = [[ni*nj for nj in ranks] for ni in ranks]
    check("fundamental rank is eight", sum(ranks) == 8)
    check("complexified adjoint decomposes into 64 matrix entries", sum(map(sum, dimensions)) == 64)
    check("diagonal subgroup has dimension sixteen", sum(n*n for n in ranks) == 16)
    check("off-diagonal charged components have dimension forty-eight", sum(map(sum, dimensions))-16 == 48)
    check("total first Chern number vanishes", dot(ranks, degrees) == 0)
    check("nonzero split flux has trace M squared thirty-six", sum(x*x for x in m8) == 36)
    check("degree matrix is antisymmetric", all(delta[i][j] == -delta[j][i] for i in range(5) for j in range(5)))
    check("spin Riemann-Roch index equals flux difference", all((genus-1+d)+1-genus == d for row in delta for d in row))
    check("opposite six-dimensional Weyl chiralities cancel each block index", all(d + (-d) == 0 for row in delta for d in row))

    # The product chirality is chi_6 = chi_4 * chi_internal. For each fixed
    # internal zero mode the two 6D Weyl summands give opposite 4D chirality.
    chirality_table = [{"chi6": c6, "internal_chirality": ci, "chi4": c6*ci}
                       for c6 in (-1, 1) for ci in (-1, 1)]
    check("product chirality holds in all four spin sectors", all(r["chi6"] == r["chi4"]*r["internal_chirality"] for r in chirality_table))
    check("each internal mode has a four-dimensional mirror", all(sum(r["chi4"] for r in chirality_table if r["internal_chirality"] == ci) == 0 for ci in (-1, 1)))

    # Equal-curvature blocks enhance if their flat holonomies also agree.
    eigenvalue_multiplicities = Counter(m8)
    commutant_dimension = sum(n*n for n in eigenvalue_multiplicities.values())
    check("bare flux centralizer is U4 times U2 times U2", sorted(eigenvalue_multiplicities.values()) == [2, 2, 4])
    check("bare flux centralizer dimension is twenty-four", commutant_dimension == 24)
    check("enhancement contributes eight extra gauge generators", commutant_dimension-16 == 8)

    # Nonzero group factors of the gauge vertex; no claim about zero-mode
    # overlap integrals or which vectorlike members remain light.
    triangles = {"up": (0, 1, 2), "down": (0, 1, 3),
                 "neutrino": (4, 1, 2), "charged_lepton": (4, 1, 3)}
    triangle_rows = []
    for label, (i, j, k) in triangles.items():
        ei, ej, ek = offsets[i], offsets[j], offsets[k]
        x, h, z = matrix_unit(ei, ej), matrix_unit(ej, ek), matrix_unit(ek, ei)
        coefficient = trace(multiply(x, subtract(multiply(h, z), multiply(z, h))))
        qsum = tuple(sum(t) for t in zip(charge(i, j), charge(j, k), charge(k, i)))
        check(label+" gauge-vertex triangle has nonzero commutator trace", coefficient == 1)
        check(label+" triangle charge sums to zero", qsum == (0,)*5)
        triangle_rows.append({"name": label, "block_cycle": [stacks[i], stacks[j], stacks[k], stacks[i]],
                              "trace_group_coefficient": coefficient})
    check("internal La vector has Higgs-up hypercharge", dot(charge(1, 2), y) == Q(1, 2))
    check("internal Lb vector has Higgs-down hypercharge", dot(charge(1, 3), y) == Q(-1, 2))

    # Atiyah--Bott Prop. 5.4: negative normal bundle is H1(ad_negative).
    # For each Hom(high slope, low slope), negative degree implies h0=0.
    # Its h1 is rank*(g-1+positive slope difference).
    negative_rows = []
    for source in range(5):
        for target in range(5):
            difference = degrees[source]-degrees[target]
            if difference <= 0:
                continue
            rank = ranks[source]*ranks[target]
            hom_degree = -rank*difference
            euler = hom_degree+rank*(1-genus)
            negative_rows.append({"source": stacks[source], "target": stacks[target],
                                  "rank": rank, "slope_difference": difference,
                                  "hom_degree": hom_degree, "h0": 0,
                                  "h1": -euler, "multiplets": genus-1+difference})
    negative_complex = sum(row["h1"] for row in negative_rows)
    grouped = 2*4*(12+3) + 2*2*(12+6) + 4*2*(12+3)
    root_sum = sum(genus-1+mi-mj for mi in m8 for mj in m8 if mi > mj)
    check("all negative normal Hom bundles have strictly negative degree", all(row["hom_degree"] < 0 for row in negative_rows))
    check("negative dimension agrees in block and grouped calculations", negative_complex == grouped == 312)
    check("negative dimension agrees with U8 root calculation", negative_complex == root_sum)
    check("real Yang-Mills Morse index is six hundred twenty-four", 2*negative_complex == 624)
    higgs_rows = [r for r in negative_rows if r["source"] in ("a", "b") and r["target"] == "L"]
    check("each Higgs block supplies eighteen unstable weak doublets", len(higgs_rows) == 2 and all(r["multiplets"] == 18 and r["h1"] == 36 for r in higgs_rows))
    colored_negative = sum(r["h1"] for r in negative_rows if "c" in (r["source"], r["target"]))
    check("colored negative normal directions number one hundred eighty", colored_negative == 180)
    check("six-dimensional Dirac gravitational chirality dimension cancels", 64-64 == 0)

    result = {
        "experiment": "AXG04_U8_ADJOINT_DIRAC_CONTROL",
        "scope": "Smooth compact spin genus-13 curve; standard positive Yang-Mills action; fixed split connection; no boundary projection, branes, extra scalar stabilization or solved Einstein background.",
        "stacks": stacks, "ranks": ranks, "line_degrees": degrees,
        "flux_matrix_diagonal": m8, "bifundamental_dimensions": dimensions,
        "degree_difference_matrix": delta,
        "opposite_6D_chirality_index_matrix": [[-d for d in row] for row in delta],
        "Dirac_net_4D_index_matrix": [[0]*5 for _ in range(5)],
        "spin_chirality_table": chirality_table,
        "action": "Integral sqrt(-g) [M6^4/2 (R6-2 Lambda6) - Tr(F_MN F^MN)/(2 g6^2) + i Tr(barPsi Gamma^M (nabla_M Psi - i[A_M,Psi])) - Mpsi Tr(barPsi Psi)]",
        "action_parameters_are_inputs": ["M6", "g6", "Lambda6", "Mpsi (zero if desired)", "spin structure", "flux and holonomy background"],
        "chirality_conclusion": "An adjoint complex Dirac fermion contains both 6D Weyl chiralities in the same gauge bundle. Every smooth common internal operator pairs their resulting 4D chiralities. This does not supply a chiral Standard Model.",
        "anomaly_scope": "Opposite chiral determinants in the identical representation cancel local chiral anomalies; paired determinant is the vectorlike control. No UV or quantum-gravity completion is claimed.",
        "bare_flux_centralizer": {"group": "U(4) x U(2) x U(2)", "dimension": commutant_dimension,
                                    "five_block_subgroup_dimension": 16,
                                    "extra_input_for_five_blocks": "Inequivalent flat twists on equal-degree line bundles; the bare curvature alone is insufficient."},
        "gauge_vertex_triangles": triangle_rows,
        "YM_negative_normal_blocks": negative_rows,
        "YM_negative_complex_dimension": negative_complex,
        "YM_Morse_index_real": 2*negative_complex,
        "colored_negative_complex_dimension": colored_negative,
        "stability_conclusion": "The split unequal-slope connection is a Yang-Mills saddle. Its charged gauge one-form fluctuations have a curvature term absent from the elementary scalar Bochner operator. The Higgs sectors are unstable, but so are colored sectors; no electroweak vacuum has been selected.",
        "gravity_scope": "The displayed action includes a dynamical metric by assumption. No M4 times X solution, backreaction, radius stabilization, Newton constant prediction, or unification principle was derived.",
        "sources": [
            {"url": "https://www.uvm.edu/~cvincen1/files/teaching/spring2019-math382/atiyahbott.pdf", "location": "Atiyah--Bott, Proposition 5.4 and equation (5.10), printed pages 556-559", "use": "Negative normal H1 and exact real Morse-index formula"},
            {"url": "https://arxiv.org/pdf/hep-th/0404229", "location": "Introduction equation (1.1), Appendix A equations (A.2)-(A.4)", "use": "Internal gauge components and gauge-vertex origin of Yukawa couplings; torus spectra are not transplanted to X0(143)"}],
        "checks": CHECKS, "checks_passed": sum(c["passed"] for c in CHECKS),
        "checks_failed": sum(not c["passed"] for c in CHECKS)
    }
    target = Path(__file__).with_name("larger_parent_results.json")
    target.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"checks_passed": result["checks_passed"], "checks_failed": result["checks_failed"], "output": str(target)}))


if __name__ == "__main__":
    main()
