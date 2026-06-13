"""
goSPL analysis tools (post-processing of model output).

Currently provides :mod:`gospl.analyse.provenance` — sediment source-to-sink
provenance attribution and a copper-fertility layer. See
``docs/DESIGN_PROVENANCE.md``.
"""

from .provenance import (  # noqa: F401
    GosplOutput,
    ProvenanceTracker,
    build_neighbours,
    downhill_edges,
    steepest_receivers,
    write_xdmf,
)
