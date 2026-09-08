# Execution notes and explanatory follow-up

The primary protocol was frozen before reading the new overlap and coupling
observables. No thresholds, operators, arithmetic blocks, triangle choices,
or local projectors were optimized in this run.

After observing W143's coupling rank four, an explanatory follow-up verified
its exact positive-sign projector rank and derived the general identity
QWP=2Q R_+ P. The rank-four upper bound and equal-sign q4/q6 cancellation are
recorded in `exact_blocks_results.json` and `W143_COUPLING_LEMMA.md`.
These explanations are not additional preregistered endpoints. Saturation of
the bound by the selected numerical projector remains DIAGNOSTIC.

The independent geometry check was extended after the primary computation
to verify its additional singular spectra and signed Gram tables through
10x16 coefficient matrices rather than the primary 26x26 maps. Agreement
was assessed against all three registered coordinate seeds. All numerical
results use double precision and shared upstream period/arithmetic inputs.

The source release was left unchanged. Known Lie-closure ambiguities from
the prior study remain recorded there; no new result requires that closure.
