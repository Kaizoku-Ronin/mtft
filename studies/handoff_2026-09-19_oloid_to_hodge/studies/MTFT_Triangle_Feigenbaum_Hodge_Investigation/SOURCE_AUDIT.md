# Source snapshot and audit scope

Reviewed source archive: `mtft-0.32.0.tar.gz`.

SHA-256: `46ffc563a0cd81a7a82d975497bf64d82f8471eb424409e3d154532b5d8d6069`.

The audit concerns this source distribution, not an assertion that every later
repository commit or external project document was inspected. The archive is not
redistributed in this bundle.

| Module beneath `src/mtft/` | Finding relevant to the proposed dynamics |
|---|---|
| `combinatorial.py` | `sigma_reflect`, `sigma_split`, `to_T_basis`, and `odd_sector_factor` implement the triangular reflection and even/odd decomposition. |
| `quadratic_forms.py` | Three-square obstruction and triangular/odd-square arithmetic; an invariant or projector does not itself specify time evolution. |
| `info_geometry.py` | `logistic_iterate` accepts the parameter r; Feigenbaum constants and the accumulation parameter are stored literals. |
| `surface/dynamics.py` | Constructs linear Hamiltonian evolution from supplied quadratic energy data. |
| `periods/hamiltonian.py` | Specified graph-derived Hermitian/quadratic Hamiltonians and symplectic frequencies; no triangle-to-logistic parameter selection found. |
| `surface/oldsector.py` | Conditional parameter selection inside a different operator family, under an imposed Hermitian Hecke selection principle. |
| `surface/marked.py` | Distinct objective functions select distinct endpoints; they do not identify the logistic parameter. |
| `research/mode_operators/` | Cohomological and operator certificate machinery; no selection of the proposed quadratic dynamics found in the inspected material. |

The reviewed `info_geometry.py` also contains a legacy nonzero expression labeled
as one-dimensional Fisher Ricci curvature. Intrinsic curvature in dimension one
vanishes. That expression was excluded from the investigation; this audit does
not establish whether a correction has already been recorded elsewhere.

No full MTFT test suite was run. Verification reported in this bundle belongs to
the standalone mathematical calculations, whose assumptions and scope are explicit
in the report and programs.
