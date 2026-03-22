"""
Full Verification Scorecard
============================

Every MTFT prediction vs PDG/CODATA experiment, with formula,
value, error, and paper reference.

This is the referee-ready comparison. Run it and see the theory
measured against nature.

  pip install mtft
  python examples/10_full_scorecard.py
"""

from mtft.verify import full_report, print_report

# ── 1. The pretty-printed report ─────────────────────────────
print_report()

# ── 2. Programmatic access ───────────────────────────────────
print("\n\nPROGRAMMATIC ACCESS")
print("=" * 55)
report = full_report()

# Count by accuracy tier
tier_1 = [p for p in report if p.error_percent < 0.1]
tier_2 = [p for p in report if 0.1 <= p.error_percent < 1.0]
tier_3 = [p for p in report if 1.0 <= p.error_percent < 5.0]
tier_4 = [p for p in report if p.error_percent >= 5.0]

print(f"  Tier 1 (< 0.1%):  {len(tier_1)} predictions")
for p in tier_1:
    print(f"    {p.name}: {p.error_percent:.4f}%")

print(f"\n  Tier 2 (0.1-1%):  {len(tier_2)} predictions")
for p in tier_2:
    print(f"    {p.name}: {p.error_percent:.4f}%")

print(f"\n  Tier 3 (1-5%):    {len(tier_3)} predictions")
for p in tier_3:
    print(f"    {p.name}: {p.error_percent:.2f}%")

print(f"\n  Tier 4 (> 5%):    {len(tier_4)} predictions")
for p in tier_4:
    print(f"    {p.name}: {p.error_percent:.1f}%  (noted conjecture)")

# ── 3. By sector ─────────────────────────────────────────────
print(f"\nBY SECTOR")
print("-" * 55)
sectors = {}
for p in report:
    sectors.setdefault(p.sector, []).append(p)

for sector, preds in sectors.items():
    avg_err = sum(p.error_percent for p in preds) / len(preds)
    best = min(preds, key=lambda p: p.error_percent)
    print(f"  {sector:15s}  {len(preds)} predictions  "
          f"avg err = {avg_err:.3f}%  best = {best.name} ({best.error_percent:.4f}%)")

# ── 4. Export to CSV ─────────────────────────────────────────
print(f"\nCSV EXPORT (paste into a spreadsheet)")
print("-" * 55)
print("name,formula,mtft_value,experimental,error_pct,paper,sector")
for p in report:
    f = p.formula.replace(",", ";")
    print(f"{p.name},{f},{p.mtft_value},{p.experimental},{p.error_percent:.4f},{p.paper},{p.sector}")
