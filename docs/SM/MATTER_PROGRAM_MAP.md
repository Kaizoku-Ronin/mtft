# The matter programme — map and toolkit (2026-10-02)

> Superseded as the programme map by `MATTER_PUZZLE_MAP.md` (generated from `mtft.matter.puzzle`, 2026-10-08). Kept for the toolkit usage notes below; the chain and open gates here are the state of 2026-10-02.

## The chain
Matter-04 (spinor twist removes adjoint fermions; jet families) -> Matter-05 (one-loop adjoint masses) -> PS-01 (Pati–Salam branch, breaking, anomalies, b–tau)
-> NU-01 / NU-02 (neutral inventory, localisation, inverse-seesaw requirements, modular parametrisation from scratch) -> AL-01 (Atkin–Lehner obstruction; one
family per cusp) -> YK-01 (jet–jet Yukawas vanish; CC-36) -> PA-01 (matter untwisting is a one-loop maximum; messengers) -> PB-01 (holomorphic route closed;
D-term mixing; cusp kernel) -> CK-01 (1/143 bound; CC-37; family triplet) -> EPS-01 (what sets epsilon).

**Canonical background:** one blown-up cusp (Matter-04), Scherk–Schwarz twist by an even spin structure (complex pair preferred), Pati–Salam broken by the
cheap gamma-flux (global Higgs kept). **Open gates:** the heavy bidoublet spectrum and the condensate's alignment (fix eps_u, eps_d); an SU(2)_R-breaking gamma-flux
(TK-1, EPS-01 item 5); Green–Schwarz axions for the stack U(1)s; the Dirac pairing M of the inverse seesaw; the Georgi–Jarlskog relations.

## Using the toolkit
```python
from mtft.matter.background import PRESETS, Background
print(PRESETS["matter04"].text())               # indices, zero modes, twists, Yukawa verdict, stability, messengers, CKM admissibility, families
b = Background("mine", ("inf",), {"4": (0,), "L": (-3,), "R": (3,)}, gamma={"4": 1, "L": -1, "R": -1})
b.analyse()                                     # the same checks for any configuration
from mtft.matter import ledger; ledger.entries(status="active"); ledger.check()
```
Command line: `python -m mtft matter presets | background NAME [--json] | ledger [--status S] [--markdown] | epsilon` (also `python -m mtft.matter ...`).
The ledger (`docs/SM/MATTER_LEDGER.md`, generated) lists every claim with its class, status, module and test; retractions cite their correction.
