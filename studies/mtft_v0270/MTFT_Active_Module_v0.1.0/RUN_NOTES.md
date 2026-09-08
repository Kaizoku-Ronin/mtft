# Execution notes

The first execution completed the baseline module and hull calculations but
stopped at the secondary Lie-closure comparison for frame seed 2026. MTFT raised
`LieGateAmbiguous`: one residual 1.2076478544427825e-7 lay in its registered
(1e-7,1e-5] ambiguity band. No package gate was loosened or source changed.

The runner was amended to record that ambiguity and continue the independently
defined primitive-kernel and arithmetic measurements. Failure of this numerical
closure comparison does not itself establish instability of the independently
computed projector. The frame/precision comparisons and independent localized
support reconstruction adjudicate that separate question. Failed closure runs
remain explicit in the result ledger.

The completed run also records ambiguity at seed 31415, with residual
1.5844900252434873e-7. Both period precisions at seed 143 pass the package gate.

After resolved primary coupling was observed, an exploratory check compared
the exact commuting full-space pairs (T2,T3) and (T2,W11) with their compressed
active-space actions. The motivation was to test whether apparent new
noncommutativity could be explained by discarded complement coupling. This
was not a preregistered primary observable. Exact commutation and the numerical
compression identity are recorded in `independent_coupling_check_results.json`.

The rank-four support explanation was discovered while independently
reconstructing the projector, before the coupling output was read. Its general
rank bound is proved separately; saturation and equality with the active
projector remain numerical, not promoted to exact by this discovery.
