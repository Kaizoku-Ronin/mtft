"""
Computation as an Arithmetic Object
====================================

The v0.7.0 computation tier: a Turing machine decomposed into the
five primitives of the Arithmetica Generale, the Hecke-constrained
Busy Beaver hierarchy on the genus-13 tape window, computational
stiffness (how much the X_0(143) constraint compresses the search
space), the Faulhaber correction structure, and the arithmetic Wick
rotation between the mass-gap (Laplace) and zeta (Dirichlet) pictures.

  pip install mtft
  python examples/14_computation.py
"""

from mtft.arithmetic_machine import (
    decompose_turing_machine,
    level_hierarchy,
    search_space_compression,
)
from mtft.arithmetic_wick import wick_at_critical_depths
from mtft.busy_beaver import (
    bb_genus,
    faulhaber_decompose,
    hecke_constraint_density,
    verify_telescoping,
)

# ── 1. A Turing machine IS a five-primitive object ───────────
print("THE FIVE-PRIMITIVE DECOMPOSITION")
print("=" * 55)
decomp = decompose_turing_machine()
print(decomp)
print()
print("Level hierarchy of computational formalisms:")
for level, names in level_hierarchy().items():
    print(f"  {level.name:>10s}: {', '.join(names)}")
print()

# ── 2. The Hecke sign oracle ─────────────────────────────────
print("HECKE CONSTRAINT DENSITY (200 traces)")
print("=" * 55)
d = hecke_constraint_density()
print(f"  fermionic (trace < 0, write 1): {d['fermionic']:.1%}")
print(f"  bosonic   (trace > 0, write 0): {d['bosonic']:.1%}")
print(f"  free      (trace = 0, either):  {d['free']:.1%}")
print()

# ── 3. The Busy Beaver hierarchy ─────────────────────────────
print("BUSY BEAVER: HECKE-CONSTRAINED vs UNRESTRICTED")
print("=" * 55)
print(f"  {'n':>3s} {'BB_Hecke':>9s} {'BB_unc':>7s} {'machines':>9s}")
for n in (1, 2):
    rh = bb_genus(n)
    ru = bb_genus(n, hecke_constrained=False)
    print(f"  {n:>3d} {rh.bb_value:>9d} {ru.bb_value:>7d} {rh.total_machines:>9,d}")
print("  (BB_Hecke(3) = 9 — exact, ~1 min of enumeration; try")
print("   `python -m mtft.busy_beaver bb 3` yourself)")
print()

# ── 4. Search-space compression ──────────────────────────────
print("SEARCH-SPACE COMPRESSION FROM THE CONSTRAINT")
print("=" * 55)
for n in (1, 2, 3):
    print(f"  n={n}: constrained space is 2^{search_space_compression(n):.2f}"
          f" times smaller")
print()

# ── 5. Faulhaber correction structure ────────────────────────
print("FAULHABER DECOMPOSITION (mass-gap analogue)")
print("=" * 55)
fd = faulhaber_decompose({1: bb_genus(1), 2: bb_genus(2)})
print(f"  naive dominant D = 2^13 = {fd.naive_dominant}")
print(f"  corrections R(n) = {fd.corrections}")
print(f"  degree bound     = {fd.degree_bound}  (canonical degree)")
print(f"  telescoping identity holds: {verify_telescoping(fd)}")
print()

# ── 6. The arithmetic Wick rotation ──────────────────────────
print("WICK ROTATION AT THE CRITICAL DEPTHS")
print("=" * 55)
print(f"  {'depth':>6s} {'y':>8s} {'g_D/g_L':>10s} {'S_D - S_L':>10s}")
for label, wr in wick_at_critical_depths(n_max=500).items():
    print(f"  {label:>6s} {wr.y:>8.5f} {wr.curvature_ratio:>10.4f}"
          f" {wr.entropy_difference:>10.4f}")
print()
print("Same weights w_n, two assemblies: e^(-2*pi*y*n) gives the")
print("mass-gap picture, n^(-beta) gives the zeta picture.")
