"""
goSPL analysis tools (post-processing of model output).

- :mod:`gospl.analyse.provenance` — sediment source-to-sink provenance
  attribution and a copper-fertility layer (``docs/DESIGN_PROVENANCE.md``).
- :mod:`gospl.analyse.stratamesh` — build a 3-D wedge (triangular-prism) volume
  of the recorded stratigraphy as an XDMF for ParaView (lithology or provenance).
- :mod:`gospl.analyse.gridexport` — rasterise a surface to a CF-NetCDF grid for
  PyGMT / ArcGIS (fields + drainage basins + chi) and extract/plot per-basin
  river longitudinal profiles.
- :mod:`gospl.analyse.stratasection` — publication stratigraphic cross-sections,
  horizontal (depth) slices, synthetic wells and Wheeler diagrams (matplotlib).

Runnable as commands — ``gospl-provenance`` / ``gospl-strata-volume`` /
``gospl-grid`` / ``gospl-section`` (after install) — or via
``python -m gospl.analyse.<tool>``. See :doc:`the running guide </user_guide/running>`.
"""

from . import gridexport  # noqa: F401  (submodule; also runnable via -m)
from . import stratamesh  # noqa: F401  (submodule; also runnable via -m)
from . import stratasection  # noqa: F401  (submodule; also runnable via -m)
from .gridexport import (  # noqa: F401
    basin_rivers,
    grid_export,
    plot_basin_map,
    plot_long_profile,
    to_netcdf,
)
from .provenance import (  # noqa: F401
    GosplOutput,
    ProvenanceTracker,
    build_neighbours,
    downhill_edges,
    steepest_receivers,
    write_xdmf,
)
from .stratamesh import (  # noqa: F401
    build_partition,
    detect_first_layer,
    discover_steps,
)
from .stratasection import (  # noqa: F401
    cross_section,
    horizontal_slice,
    load_strata,
    synthetic_well,
    wheeler,
)
