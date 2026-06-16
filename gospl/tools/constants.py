"""
Shared numeric constants used as sentinels, thresholds, and floors
throughout goSPL.

Every constant defined here is documented in AGENTS.md under the
"Magic numbers" section. The rule from AGENTS.md is:

    These literal values MUST NOT be reintroduced inline. New uses
    should reference this module.

When you add a constant here:

  1. Add or update its entry in AGENTS.md > Magic numbers so the audit
     list and the code stay in sync.
  2. Replace every inline literal of the same VALUE that ALSO matches
     the SAME SEMANTIC ROLE. Inline literals that share a value but
     mean something different (e.g. a solver tolerance that happens to
     equal a deposit floor) MUST NOT be replaced — leave them with a
     `# TODO-REFACTOR: value matches X but distinct role` comment so
     the next reader knows the coincidence is intentional.
  3. Run `pytest tests/` to confirm no regression.
"""

# Pre-fill value before MPI.COMM_WORLD.Allreduce(MAX) for nodes that
# are not owned by the local partition, AND boundary marker before
# fitedges interpolation. Plausible physical range is bounded by ~9 km
# (Earth's tallest point), so a value of -1e8 m is unambiguously
# recognisable as "no data" after the reduction. Used by:
#   - mesher/tectonics.py (Allreduce pre-fill + idBorders marker)
#   - mesher/unstructuredmesh.py (idBorders marker before fitedges)
#   - tools/addprocess.py (Allreduce pre-fill in applyFlexure / cptOrography)
#   - flow/pitfilling.py (sentinel value in _performFilling graph)
# See AGENTS.md > Magic numbers.
MISSING_DATA_SENTINEL = -1.0e8

# Same role as MISSING_DATA_SENTINEL but for fields whose plausible
# magnitude can exceed 1e8 over a single step (cumulative erosion/
# deposition, flexural isostasy deflection, soil column thickness).
# Used in mesher/tectonics.py::_advectPlates as the Allreduce(MAX)
# pre-fill for the gED / gFI / gSL global arrays. Two orders of
# magnitude below MISSING_DATA_SENTINEL so the two sentinels never
# collide.
MISSING_LARGE_SENTINEL = -1.0e10

# Minimum flow accumulation / discharge / sediment-load value written
# to the on-disk HDF5 outputs. Anything at or below this floor is
# clamped to this floor (avoids -inf in log10 visualisation and
# absorbs the trickle of numerical noise from KSP residuals). Used by:
#   - tools/outmesh.py (FA, fillFA, iceH, iceUb, iceMelt, iceAbr, sedLoad output fields)
# See AGENTS.md > Magic numbers.
DISCHARGE_FLOOR = 1.0e-8

# Drop sub-millimetre marine sediment deposits as numerical noise.
# Used by sed/seaplex.py::seaChange right before the diffusion step.
# Prevents the round-off accumulation of vanishingly thin deposits
# across many timesteps from drifting the mass balance. See
# AGENTS.md > Magic numbers.
DEPOSIT_FLOOR = 1.0e-3

# Soil thickness below which the cell is treated as bedrock-exposed in
# the soil-aware SPL kernels. Used by eroder/soilSPL.py both before
# the SNES residual computation and after the soil-thickness update.
# See AGENTS.md > Magic numbers.
BEDROCK_EXPOSED = 1.0e-1

# Sentinel thickness for the layer-0 "infinite bedrock" reservoir
# when no initial stratigraphy file is supplied. The 1e6 m sentinel
# is added and subtracted around the per-node cumsum in
# sed/stratplex.py::erodeStrat so that the bedrock reservoir cannot
# ever be exhausted by a single erosion step, while the sentinel
# offset cancels in the cumulative-thickness arithmetic and does not
# contaminate eroded volumes. See AGENTS.md > Magic numbers.
BEDROCK_SENTINEL = 1.0e6

# Sentinel passed to the Fortran multi-flow-direction kernels
# (`mfdreceivers`, `mfdrcvrs`) as the "this neighbour is no-data,
# skip it for flow routing" floor. Also used to mask deep-ocean and
# open-boundary cells before the kernel call so they are not
# considered as downstream receivers.
#
# DISTINCT from MISSING_DATA_SENTINEL: this one means "ignore for
# flow routing", not "pre-fill before Allreduce(MAX)". The two values
# are kept far enough apart (-1e6 vs -1e8) that a mistake would be
# obvious in any diagnostic dump.
#
# Used by:
#   - flow/iceplex.py (border mask + mfdreceivers no-data floor)
#   - sed/seaplex.py (border + deep-ocean masks + mfdrcvrs no-data floor)
# See AGENTS.md > Magic numbers.
BOUNDARY_FLOW_SENTINEL = -1.0e6

# Upper-bound clamp for the global spillover graph in
# flow/pitfilling.py::_performFilling. Entries above this value are
# treated as outliers (stale initialiser values or numerical blow-up)
# and rewritten to MISSING_DATA_SENTINEL. The value 1e7 m is two
# orders above Earth's ~9 km elevation cap, so no real physical
# elevation can ever exceed it. Paired with MISSING_DATA_SENTINEL —
# see the sentinel-filter comment in pitfilling.py for the full
# rationale.
GRAPH_OUTLIER_CAP = 1.0e7

# Dimensionless diffusion numbers for the marine-morphology smoothing used
# ONLY to derive coherent marine flow directions in sed/seaplex.py::_matOcean
# (_hillSlope(smooth=2)). This smoothing never alters elevation — it produces a
# de-noised bathymetry so the multiple-flow-direction kernel routes marine
# sediment down the large-scale seafloor slope instead of trapping it in small
# local minima.
#
# The smoothing solves one implicit diffusion step (I - L) h_s = h with a
# per-node coefficient Kd = N * cell_area (m^2). Because the FV Laplacian
# stencil scales as Kd / Δx^2, folding the cell area into Kd makes N a
# dimensionless smoothing strength that is BOTH timestep-independent AND
# mesh-resolution-independent — replacing the former hard-coded Kd = Cd*dt with
# Cd in {1e5 (land), 5e6 (sea)} m^2/yr, which scaled with dt and with 1/Δx^2 so
# the effective smoothing silently varied with the timestep and the mesh.
#
# The seafloor (SEA) is smoothed more heavily than emergent land (5x here):
# routing only needs coherent bathymetric directions, while coastal land just
# needs to stay continuous for the pit-fill + flow-direction pass. These are
# moderate de-noising strengths (v ~ N per neighbour) — enough to make the
# seafloor monotone toward sinks without flattening the large-scale slope;
# increase for stronger smoothing. NB the marine routing (and hence the
# downstream conservation residual on a coarse mesh) is sensitive to this, so
# retune against the dual/marine conservation tests if you change it. See
# AGENTS.md > Magic numbers.
MARINE_SMOOTH_N_LAND = 1.0
MARINE_SMOOTH_N_SEA = 5.0
