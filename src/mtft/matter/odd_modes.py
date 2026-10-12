"""OM-01 (2026-10-08): the forced-odd-mode count — the parity audit's verdict made topological.

Artebani–Quispe's lesson (an odd count of fixed points forces something to exist), transplanted: the torus flip w -> -w on 143a1 fixes the four 2-torsion
points, and the holomorphic Lefschetz formula fixes how many odd modes a symmetric line bundle must carry.

1. The count (EXACT).  For the symmetric bundle O(N O) with the natural lift, sections are functions with poles of order <= N at O: 1, wp, wp^2, ... are even and
   wp', wp wp', ... are odd.  So #even = 1 + floor(N/2), #odd = floor((N - 1)/2), and #even - #odd = ((-1)^N + 3)/2 — the Lefschetz sum (1/2) sum over the four fixed
   points of the fibre signs ((-1)^N at O, +1 at the three others).  Odd modes are forced from N = 3 up.
2. The design has no torus degree anywhere (EXACT).  One unit of torus flux on a block makes 24 chiral families (TP-01), so every block has fibre degree 0; the gamma
   charges are Poincare-type, also fibre degree 0 (EPS-01).  Flat bundles' zero modes are constant along the torus — even.  The only modes topology forces to be
   odd are the w-jet families themselves (one per block, from the exceptional charge at the jet point): states, not sources.
3. Verdict (EXACT).  No field that could carry an odd vev or an odd coupling in the down sector — the down Higgs, any B-stack scalar, the vector-like pairs — has
   an odd mode at all.  Giving a B-stack block the three units of torus flux an odd mode needs would cost 72 chiral families.  The parity audit's verdict (PA-1)
   is a topological statement, not an inventory."""

BLOCKS = {"Q": ("c", "L"), "L_doublet": ("A1", "L"), "u^c": ("c", "A2"), "d^c": ("c", "B"), "nu^c": ("A1", "A2"), "e^c": ("A1", "B"),
          "H_u": ("L", "A2"), "H_d": ("L", "B")}
FIBRE_DEGREE = {k: 0 for k in BLOCKS}            # TP-01: b = 0 throughout; gamma charges are fibre-degree 0
FAMILIES_PER_TORUS_UNIT = 24


def split_by_poles(N):
    """Even and odd sections of O(N O) under w -> -w, from pole orders (the natural lift)."""
    if N <= 0:
        return (1, 0) if N == 0 else (0, 0)
    orders = [0] + list(range(2, N + 1))
    return sum(1 for k in orders if k % 2 == 0), sum(1 for k in orders if k % 2 == 1)

def lefschetz(N):
    """(1/2) sum over the four fixed points of the fibre signs: (-1)^N at O, +1 at the three 2-torsion points."""
    return ((-1) ** N + 3) // 2

def first_forced_odd_degree():
    return min(N for N in range(1, 20) if split_by_poles(N)[1] > 0)

def b_stack_blocks():
    return [k for k, v in BLOCKS.items() if "B" in v]

def odd_modes_on_b_stack():
    """Blocks on the B stack whose torus profile can be odd: needs |fibre degree| >= 3."""
    return [k for k in b_stack_blocks() if abs(FIBRE_DEGREE[k]) >= first_forced_odd_degree()]

def family_cost(b):
    return FAMILIES_PER_TORUS_UNIT * abs(b)
