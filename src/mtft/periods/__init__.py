"""mtft.periods — genuine period/Hodge geometry of X0(143). (v0.21 dev)

NOTE: this __init__ and __main__ are verification scaffolding written by
Claude to run Sol's five modules; replace with Sol's originals before release.
"""
from .core import (GENUS, LEVEL, data_path, period_record, intersection_form,
                   intersection_inverse, symplectic_change, symplectic_form,
                   omega_cusp, omega_symplectic, legacy_omega_symplectic,
                   frozen_riemann_matrix, riemann_matrix, normalized_periods,
                   hodge_complex_structure, hodge_metric, charge_energy)
from .bridge import (relative_basis_change, cuspidal_basis_change,
                     hecke_to_symplectic_change)
from .forms import (raw_qexpansions, q_tail_bound, raw_form_values,
                    normalized_form_values, bergman_density)
from .physics import (hodge_structure_hecke, hodge_metric_hecke, graph_coupling,
                      complex_linear_decomposition, metric_hs_norm,
                      cp_channel_report, finite_charge_partition)
from . import gates
