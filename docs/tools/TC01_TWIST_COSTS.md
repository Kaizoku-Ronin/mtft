# TC-01 — what a twist buys and costs at the Pati–Salam point (2026-11-19)

Module `mtft.matter.twist_costs`; tests `tests/test_twist_costs.py`. Ported from the project records FF-12 (the trace test) and FF-13
(the twist table), 19 November 2026, whose plans, scripts and numbers are filed under `claude/experiments/` in the MTFT Python project.

What it is: VE-01's one-loop Wilson-line energy V = −Σ_w W(w)|Tr U^w|² with any non-negative weight W on the windings of 143a1, and
its exact Hessian at the Pati–Salam pattern on an offset h ∈ E[2]∖O:

    M_A = Σ_{χ_h(w)=+1} W w wᵀ,   M_B = Σ_{χ_h(w)=−1} W w wᵀ,   Hess V = −8π² M_B ⊗ s sᵀ + 8π² M_A ⊗ (8·1 − J).

T_S has curvature 64π² eig(M_A − M_B), the twelve 64π² eig(M_A) (six times each). Route 2 is the direct sum with Richardson second
differences, compared on the Hessian's own scale (agreement below 10⁻⁹; a per-direction relative error is meaningless where a
curvature nearly cancels, FF-13's recorded fault).

The twist table (EXACT at one loop). For P ∈ E[2] (a change of spin structure, MATTER-04), tree Yukawas ⟺ P = h ⟺ V = 0 ⟺ M_A = 0
⟺ the twelve flat, on all nine pairs. Among the spin structures of 143a1 the pin of the twelve and the tree Yukawas exclude each
other; for P ≠ h the energy (−1.8 to −6.2 M_KK⁴) also pulls τ off 143a1 with sinh d ≥ 620 (ESTIMATE, `hold_pull`).

The continuous twist P = (0, (1 − δ)/2) (a model setting: it needs a flat R-symmetry bundle L with L² ≠ O, outside the Calabi–Yau
case): quark–lepton fermions δ·M_SS, Φ-type Higgs bosons 2δ·M_SS, Φ in offset-h blocks (1 − 2δ)·M_SS (EXACT, closed form against a
brute-force minimum to 10⁻¹²). Its pin reproduces FF-12's even-winding pin.

Gates: the anchors of FF-13 (T_S, the twelve, V on three pairs), FF-12's pin at δ = 0.05, the jet vacuum's pull 8.73–10.63 M_KK⁴
(FF-4), the continuous twist's hold (max sinh d = 0.011 at δ = 0.05), and every raise path.
