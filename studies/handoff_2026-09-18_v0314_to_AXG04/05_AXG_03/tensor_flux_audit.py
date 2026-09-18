"""AXG-03: normalized tensor transgression and the factor-three lattice gate.

Python standard library only. This verifies the explicitly stated conditional
algebra. It does not construct a six-dimensional anomaly-free tensor theory.

Conventions: f_a=F_a/(2*pi), integral axion coordinates have period 1,
X4_abelian=(1/2) b_ab f_a f_b with symmetric b. Then the flux pushforward
is (b*m)_a f_a. The physical 2*pi-periodic axion has delta theta=-q_a Lambda_a.

Global scope: the coefficient pairing is integral on the actual cocharacter
lattice into the primitive string/axion lattice. A rational local solution
does not establish that scope. See report for self-duality and BF modes.
"""

from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
from pathlib import Path
import json


CHECKS = []


def check(name, condition):
    ok = bool(condition)
    CHECKS.append({"name": name, "passed": ok})
    if not ok:
        raise AssertionError(name)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def matvec(a, b):
    return tuple(dot(row, b) for row in a)


def scale(c, a):
    return tuple(c*x for x in a)


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def zero_matrix(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def smith_two_rows(rows):
    """Nonzero Smith invariants from determinantal divisors, for rank 2."""
    entries = [abs(int(x)) for row in rows for x in row]
    d1 = reduce(gcd, entries)
    minors = [abs(rows[0][i]*rows[1][j]-rows[0][j]*rows[1][i])
              for i, j in combinations(range(len(rows[0])), 2)]
    determinant_divisor = reduce(gcd, [int(x) for x in minors])
    if not d1 or not determinant_divisor:
        raise ValueError("Expected rank two")
    return [d1, determinant_divisor//d1]


def encode(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, (tuple, list)):
        return [encode(x) for x in obj]
    if isinstance(obj, dict):
        return {k: encode(v) for k, v in obj.items()}
    return obj


def main():
    m = (0, -3, 3, 3, 0)
    m0 = (0, -1, 1, 1, 0)
    k1 = (0, -2, 1, 1, 0)
    k2 = (3, -4, 0, 0, 1)
    kd = (3, -2, -1, -1, 1)
    p = (1, 1, 1, 1, 1)
    Y = (F(1, 6), 0, F(-1, 2), F(1, 2), F(-1, 2))
    BL = (F(1, 3), 0, 0, 0, -1)
    check("flux has cocharacter factor three", m == scale(3, m0))
    check("primitive native direction", kd == add(k2, scale(-1, k1)))
    check("native direction annihilates flux", dot(kd, m) == 0)
    check("primitive k1 fails native-scalar degree test", dot(k1, m) == 12)

    # X4 = f_L^2 + f_a f_b in the central restriction. Globally it is
    # c2(E_L) + c1(L_a)c1(L_b), an ordinary integral characteristic class.
    b = zero_matrix(5)
    b[1][1] = F(2)
    b[2][3] = b[3][2] = F(1)
    q = matvec(b, m)
    check("integer symmetric tensor pairing", all(b[i][j] == b[j][i]
          and b[i][j].denominator == 1 for i in range(5) for j in range(5)))
    check("even diagonal makes integral central quadratic class",
          all(int(b[i][i]) % 2 == 0 for i in range(5)))
    check("primitive-flux control gives k1", matvec(b, m0) == k1)
    check("actual-flux control gives 3 k1", q == scale(3, k1))
    check("tensor row need not annihilate flux", dot(q, m) == 36)

    # Verify the polarization normalization, rather than silently adopting
    # an ambiguous K=bm convention. X=(1/2) f^T b f; coefficient of w in
    # X(f+m*w), w^2=0, is the sum of the two symmetric cross terms.
    cross = tuple(sum((b[i][j]*m[j]+b[j][i]*m[j])/2
                      for j in range(5)) for i in range(5))
    check("quadratic polarization includes both flux cross terms", cross == q)
    check("naive half-descent row differs from full transgression",
          scale(F(1, 2), q) != q)

    # c2 of a rank-two bundle has roots x1,x2. Twisting by a line of
    # degree -3 gives c2' = c2 - 3(x1+x2) w. The a,b product adds
    # 3(f_a+f_b)w. In center coordinates x1+x2=2*f_L.
    determinant_row = (0, -3, 3, 3, 0)
    determinant_to_center = (3, 2, 1, 1, 1)
    recovered_center = tuple(determinant_row[i]*determinant_to_center[i]
                             for i in range(5))
    check("integral U(2) characteristic class gives same center row",
          recovered_center == q)
    native_det_row = (1, -1, -1, -1, 1)
    check("native row is an actual determinant character",
          tuple(native_det_row[i]*determinant_to_center[i]
                for i in range(5)) == kd)

    # There is no local linear-algebra obstruction: dividing b by three
    # yields q=k1. The new denominators are exactly the missing global data.
    b_rat = [[x/3 for x in row] for row in b]
    check("rational local control reaches primitive k1", matvec(b_rat, m) == k1)
    check("rational control fails primitive integer coefficient lattice",
          any(x.denominator == 3 for row in b_rat for x in row))
    check("rational class has a one-third product period on S2 x S2",
          b_rat[2][3] == F(1, 3))

    # Symbolically, for every integral symmetric b, (bm)_i is
    # 3*(-b_iL+b_ia+b_ib). Testing its coefficient vector is an identity,
    # not a bounded search over b values.
    check("universal transgression coefficient identity", m == scale(3, m0))
    check("target is nonzero modulo three", any(x % 3 for x in k1))
    check("adding any integral native direction cannot repair divisibility",
          all(any((k1[i]-n*kd[i]) % 3 for i in range(5))
              for n in range(3)))
    original_snf = smith_two_rows([kd, k1])
    changed_snf = smith_two_rows([kd, q])
    actual_group_snf = smith_two_rows([native_det_row, determinant_row])
    check("original periodic map is primitive", original_snf == [1, 1])
    check("native plus integral tensor has index three", changed_snf == [1, 3])
    check("unitary-center quotient does not remove index three",
          actual_group_snf == [1, 3])
    check("factor three cannot be removed by integral row operations",
          original_snf != changed_snf)
    check("same common-phase real kernel", dot(kd, p) == dot(q, p) == 0)
    check("same hypercharge real kernel", dot(kd, Y) == dot(q, Y) == 0)
    check("same B-L real kernel", dot(kd, BL) == dot(q, BL) == 0)

    # The explicit integral control is not complete 6D GS cancellation.
    # No anomaly polynomial, tensor signature, action, or scalar potential
    # is inferred from its existence.
    result = {
        "study": "AXG-03 tensor flux transgression",
        "scope": "Conditional exact lattice calculation, not a completed tensor parent",
        "normalization": {
            "f": "F/(2*pi)", "normalized_axion_period": 1,
            "physical_axion_period": "2*pi",
            "abelian_class": "X4=(1/2)*b_ab*f_a*f_b",
            "induced_shift": "q=b*m; delta(theta)=-q_a*Lambda_a",
            "patching": "Includes both symmetric flux cross terms; naive local half-descent is insufficient",
        },
        "m": m, "m0": m0, "k1": k1, "k_delta": kd,
        "integral_control_class": "c2(E_L)+c1(L_a)*c1(L_b)",
        "integral_control_pairing": b,
        "integral_control_shift": q,
        "rational_local_control_pairing": b_rat,
        "smith_invariants_original": original_snf,
        "smith_invariants_native_plus_integral_tensor": changed_snf,
        "smith_invariants_actual_unitary_determinant_basis": actual_group_snf,
        "global_gate": "If the transgression pairing is integral on the actual primitive cocharacter and axion/string lattices, m=3*m0 forces q in 3 times the axion charge lattice. Primitive k1 is excluded in that class.",
        "possible_changes_requiring_new_data": [
            "Different primitive flux or gauge global form and allowed bundles",
            "Fractional tensor/axion lattice with specified discrete or relative sector",
            "Localized defects or additional non-flux shift sources",
            "Additional native axions or additional Higgs condensates",
        ],
        "tensor_modes_before_flux_interactions": {
            "genus": 13,
            "nonchiral": {"independent_scalars": 2, "harmonic_one_form_vectors": 26},
            "self_dual": {"independent_scalars": 1, "harmonic_one_form_vectors": 13},
            "note": "Self-duality relates internal-integral and dual-external scalar; flux can mix/lift these kinematic modes.",
        },
        "sources": [
            {"url": "https://arxiv.org/pdf/1506.05771", "location": "equations (55)-(60), (73)-(76)", "use": "Explicit factor-two flux patching and electric/dual axions"},
            {"url": "https://arxiv.org/pdf/1711.04777", "location": "equations (2.27), (3.32)", "use": "String charge quantization and actual cocharacter lattice"},
        ],
        "checks": CHECKS,
        "checks_passed": sum(c["passed"] for c in CHECKS),
        "checks_failed": sum(not c["passed"] for c in CHECKS),
    }
    target = Path(__file__).with_name("tensor_flux_results.json")
    target.write_text(json.dumps(encode(result), indent=2)+"\n")
    print(json.dumps({"checks_passed": result["checks_passed"],
                      "checks_failed": result["checks_failed"],
                      "output": str(target)}))


if __name__ == "__main__":
    main()
