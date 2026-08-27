"""mtft.periods — genuine period/Hodge geometry of X0(143).

v0.22 candidate __init__: replaces the v0.21.0 scaffolding version, which
shipped without ``__all__`` (compatibility_audit: period_declares___all__
false).  The export surface below is exactly the 36 names observed live in
v0.21.0, so nothing that worked before breaks; the change is that the public
API is now declared rather than implied.

Legend registration for these names remains an open item (the curated
registry still stops at v0.15-era entries plus c9b_period_sign_adjudication).
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
from .involutions import (al_matrix, al_signs, sector_census,
                          route2_fixed_intersections, oldspace_projector,
                          star_symplectic, star_charge_orbit)
from .oldtorus import (old_lattice, polarization_type, l9_index,
                       principal_form, product_charpoly, entropy)
from .hamiltonian import (hodge_adjoint, hermitian_split, hamiltonian_split,
                          channel_report, pairing_stability,
                          symplectic_frequencies, oldspace_routing,
                          hecke_block_routing)
from .channels import (bergman_bilinear, bergman_channel, channel_density,
                       mode_crossover)
from . import (bridge, channels, core, forms, gates, hamiltonian,
               involutions, oldtorus, physics)

__all__ = [
    # constants and data access
    "GENUS", "LEVEL", "data_path", "period_record",
    # exact integral structures
    "intersection_form", "intersection_inverse", "symplectic_change",
    "symplectic_form",
    # periods and Hodge geometry
    "omega_cusp", "omega_symplectic", "legacy_omega_symplectic",
    "frozen_riemann_matrix", "riemann_matrix", "normalized_periods",
    "hodge_complex_structure", "hodge_metric", "charge_energy",
    # exact basis bridge
    "relative_basis_change", "cuspidal_basis_change",
    "hecke_to_symplectic_change",
    # q-expansions and Bergman density
    "raw_qexpansions", "q_tail_bound", "raw_form_values",
    "normalized_form_values", "bergman_density",
    # physics-facing diagnostics
    "hodge_structure_hecke", "hodge_metric_hecke", "graph_coupling",
    "complex_linear_decomposition", "metric_hs_norm", "cp_channel_report",
    "finite_charge_partition",
    # involutions and sectors (v0.22)
    "al_matrix", "al_signs", "sector_census", "route2_fixed_intersections",
    "oldspace_projector", "star_symplectic", "star_charge_orbit",
    # oldspace abelian surface (v0.22)
    "old_lattice", "polarization_type", "l9_index", "principal_form",
    "product_charpoly", "entropy",
    # quadratic Hamiltonian layer (v0.22)
    "hodge_adjoint", "hermitian_split", "hamiltonian_split",
    "channel_report", "pairing_stability", "symplectic_frequencies",
    "oldspace_routing", "hecke_block_routing",
    # Bergman harmonic channels (v0.22)
    "bergman_bilinear", "bergman_channel", "channel_density",
    "mode_crossover",
    # submodules
    "bridge", "channels", "core", "forms", "gates", "hamiltonian",
    "involutions", "oldtorus", "physics",
]
