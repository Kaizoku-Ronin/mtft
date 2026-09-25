"""SM interaction reference and MTFT evidence ledger (no amplitude evaluation).

The bundled UFO snapshot is external reference data, not an MTFT derivation.
Use ``python -m mtft.interactions --help`` for the offline atlas and exports.
"""

from .catalog import (
    SECTORS,
    STATUSES,
    catalog_digest,
    load_catalog,
    new_ledger,
    select_vertices,
    summary,
    validate_catalog,
    validate_ledger,
)

__all__ = [
    "SECTORS", "STATUSES", "catalog_digest", "load_catalog", "new_ledger",
    "select_vertices", "summary", "validate_catalog", "validate_ledger",
]
