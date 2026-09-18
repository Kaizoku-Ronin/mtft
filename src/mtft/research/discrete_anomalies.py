"""INT-08 (AXG-03 §5): Spin x Z_n fermion anomaly test (Hsieh) on a net chiral charge ledger.

For a Z_n generator with integer charges q on Weyl fermions: S1 = sum q, S3 = sum q^3 (with multiplicities).  Conditions:
(n^2 + 3n + 2) S3 = 0 mod 6n (cubic) and 2 S1 = 0 mod n (linear, gravitational).  M1's candidate Z3: the residual gauge transformation
exp(2 pi i g/3) with g = (0, 0, 1, 0, 1) fixes both axion circles (K_hat g = (0, 3)) and lies in the native-axion kernel
(k_Delta g = 0); its fermion ledger has S1 = S3 = 24: cubic residue 12 mod 18 (S3 = 6 mod 9) -> the fermion-only test FAILS; the full
GS/topological contribution was not constructed, so this is not a theorem about every completion (H-18)."""
import sympy as sp

def m1_ledger(generator=(0, 0, 1, 0, 1), ranks=(3, 2, 1, 1, 1), degrees=(0, -3, 3, 3, 0), names=("c", "L", "a", "b", "d")):
    """Net chiral bifundamentals with orientation sign(m_i - m_j), multiplicity |m_i - m_j| N_i N_j, integer charge
    orientation * (g_i - g_j)."""
    out = []
    for i in range(5):
        for j in range(i + 1, 5):
            index = degrees[i] - degrees[j]
            if index == 0: continue
            o = 1 if index > 0 else -1; q = o * (generator[i] - generator[j]); count = abs(index) * ranks[i] * ranks[j]
            out.append({"sector": names[i] + names[j], "orientation": o, "multiplicity": abs(index), "dimension": ranks[i] * ranks[j], "charge": q, "count": count})
    return out

def spin_zn_fermion_test(ledger, n):
    S1 = sum(e["count"] * e["charge"] for e in ledger); S3 = sum(e["count"] * e["charge"] ** 3 for e in ledger)
    cubic = ((n * n + 3 * n + 2) * S3) % (6 * n); linear = (2 * S1) % n
    return {"n": n, "S1": S1, "S3": S3, "cubic_residue_mod_6n": cubic, "linear_residue_mod_n": linear, "passes": cubic == 0 and linear == 0,
            "phase": sp.frac(sp.Rational(n * n + 3 * n + 2, 6 * n) * S3), "scope": "fermion-only; topological/GS contributions not included"}

def m1_z3_generator_certificate(K=None):
    """g = (0, 0, 1, 0, 1) generates the Z3 remnant: K_hat g = (0, 3) (both axion circles fixed mod the period), k_Delta g = 0."""
    K = K if K is not None else sp.Matrix([[0, -2, 1, 1, 0], [3, -4, 0, 0, 1]]); k1 = K.row(0); kd = K.row(1) - k1; Kh = sp.Matrix.vstack(kd, 3 * k1); g = sp.Matrix([0, 0, 1, 0, 1])
    return {"K_hat": Kh, "K_hat_g": list(Kh * g), "k_delta_g": (kd * g)[0], "generator": list(g), "certifies_Z3": list(Kh * g) == [0, 3] and (kd * g)[0] == 0}
