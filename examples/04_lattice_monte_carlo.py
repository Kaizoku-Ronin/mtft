"""
Lattice Monte Carlo: MTFT Gauge Theory in Action
=================================================

Run an actual SU(3) lattice simulation with the MTFT action:

  S_MTFT = S_Wilson + S_arithmetic

where S_arithmetic adds arithmetic-weighted Polyakov loop terms.
Watch the plaquette thermalize, measure Wilson loops, and extract
the string tension via the Creutz ratio.

  pip install mtft
  python examples/04_lattice_monte_carlo.py
"""

from mtft.lattice import (
    LatticeConfig, MTFTAction,
    metropolis_sweep, avg_plaquette, avg_polyakov,
    wilson_loop, creutz_ratio,
)

# ── 1. Set up a small lattice ────────────────────────────────
print("SU(3) LATTICE GAUGE THEORY WITH MTFT ACTION")
print("=" * 55)
cfg = LatticeConfig(N=3, L=4)          # SU(3) on a 4^4 lattice
action = MTFTAction(beta=6.0, kappa=1.0, y=0.18)

print(f"  Gauge group: SU({cfg.N})")
print(f"  Lattice:     {cfg.L}^4 = {cfg.L**4} sites")
print(f"  Wilson beta: {action.beta}")
print(f"  MTFT kappa:  {action.kappa}")
print(f"  Depth y:     {action.y}")
print()

# ── 2. Cold start and thermalization ─────────────────────────
print("THERMALIZATION (cold start -> equilibrium)")
print("-" * 55)
print(f"  {'Sweep':>5s}  {'<plaq>':>8s}  {'|<Poly>|':>10s}  {'Acceptance':>10s}")
print(f"  {'-'*38}")
print(f"  {'cold':>5s}  {avg_plaquette(cfg):8.4f}  {abs(avg_polyakov(cfg)):10.4f}  {'---':>10s}")

for sweep in range(1, 21):
    acc = metropolis_sweep(cfg, action, n_hits=5, epsilon=0.2)
    if sweep <= 5 or sweep % 5 == 0:
        plaq = avg_plaquette(cfg)
        poly = abs(avg_polyakov(cfg))
        print(f"  {sweep:5d}  {plaq:8.4f}  {poly:10.4f}  {acc:10.1%}")

# ── 3. Wilson loops and string tension ───────────────────────
print(f"\nWILSON LOOPS W(R, T)")
print("-" * 55)
print(f"  {'R':>3s}  {'T':>3s}  {'Re Tr W(R,T)':>14s}")
for R in range(1, 3):
    for T in range(1, 3):
        W = wilson_loop(cfg, (0,0,0,0), R, T)
        tr = float(W.trace().real / cfg.N)
        print(f"  {R:3d}  {T:3d}  {tr:14.4f}")

# ── 4. Creutz ratio (string tension proxy) ───────────────────
print(f"\nCREUTZ RATIO chi(R, T) = -ln[W(R,T)W(R-1,T-1) / W(R,T-1)W(R-1,T)]")
print("-" * 55)
for R in range(2, 4):
    for T in range(2, 4):
        chi = creutz_ratio(cfg, R, T)
        print(f"  chi({R},{T}) = {chi:.4f}")

print(f"\n  The Creutz ratio approaches the string tension sigma*a^2")
print(f"  as R,T -> infinity.  MTFT predicts sigma is controlled by")
print(f"  the arithmetic stiffness mu_3(y).")
