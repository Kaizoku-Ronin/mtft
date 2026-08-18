# CC-08 (append-only) — AG Pr 3.7.5: Atkin–Lehner does not split Galois orbits

**Filed:** 2026-08-18, v0.16.0 integration, first raised as a **proposal** in the
CI wave (`X0_143_CI_C_REPORT.md` lineage; wave files byte-preserved).
**Class:** [EXACT] — structural argument plus the auditor's independent exact
recomputation of the full W_Q action (`ci_verify_kimi.py`, `ci_verify_kimi.json`).

## The claim being corrected

Arithmetica Generale **Pr 3.7.5** asserts that the Atkin–Lehner involutions
W_11, W_13 split the newform Galois orbits of S_2(Gamma_0(143)): that within
the quartic orbit f2 (resp. the sextic orbit f3), different embeddings carry
different W_Q eigenvalues.

## Why it is false

For a newform f with trivial character at level N with p || N, the W_p
eigenvalue is -a_p(f), and a_p(f) = ±1 exactly. Since ±1 is rational, the
eigenvalue is fixed by every embedding of the Hecke field: **W_Q eigenvalues
are constant on Galois orbits.** No orbit split is possible.

## Corrected sector table (auditor exact replay)

W_Q built from scratch on the 26-dim cuspidal homology of the `mtft.hecke`
Manin model via the Cremona endpoint route (flag → primitive SL2(Z) lift →
Mobius action on the endpoints {b/d, a/c} → continued-fraction reconversion),
calibrated by reproducing the model's star involution exactly from the
endpoint action of J = diag(-1, 1). Results:

    W_11^2 = W_13^2 = W_143^2 = I,   W_11 W_13 = W_143 (both orders),
    [W_Q, T_p] = 0 (p = 2,3,5,7),    [W_Q, iota*] = 0
    traces:            W_11: 2,  W_13: -2,  W_143: -18
    eigenspace dims:   W_11 (14,12),  W_13 (12,14),  W_143 (4,22)
    quotient genera:   X/W_11 = 7,  X/W_13 = 6,  X/W_143 = 2,  X* = 1

Block purity (homology dims; halve for S_2):

    ell  (143a1, dim 2):  (+,+) : 2
    old  (11a1 ghost, 4): (-,+) : 2,  (-,-) : 2
    q4   (f2 orbit,   8): (-,+) : 8   — uniform on the whole orbit
    q6   (f3 orbit,  12): (+,-) : 12  — uniform on the whole orbit

S_2 joint sectors: (+,+): 1, (-,+): 5, (+,-): 6, (-,-): 1  (total 13).

## Downstream

Any corpus statement keyed to an AL sign varying inside a Galois orbit is
void. The CI wave's own sector table (used for the adapted basis and the
descent structure) is the corrected one and is unaffected.
