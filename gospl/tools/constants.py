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

# Round-off floor of deposits (1 nm): the routed marine deposit is set to zero
# below it right before the marine diffusion (sed/seaplex.py::seaChange), and
# deposeStrat records every deposit above it (it skipped those under 0.1 mm). It
# was 1 mm until 2026-11-25: real sub-millimetre deposits were discarded
# without putting the volume back (4e-5 of the sediment moved on the closed
# test sphere, ~5 % of thinly spread coastal-erosion material). See AGENTS.md
# > Magic numbers.
DEPOSIT_FLOOR = 1.0e-9

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

# Ice thickness (m) above which a land cell is treated as ice-covered for the
# subaerial soil gate: pedogenic soil production is suppressed and the existing
# regolith is held FROZEN INERT (preserved, not zeroed — unlike the subaqueous
# marine/lake case) under ice. Matches the ice-presence threshold used inline in
# flow/iceplex.py (_glacialMeltwater / _routeTill), which carries a distinct role
# there. Used by eroder/soilSPL.py::_iceFrozenMask. See AGENTS.md > Magic numbers.
ICE_COVER_MIN = 1.0e-2

# Relative floor (fraction of the domain's maximum smoothed ice discharge) below
# which the smoothed discharge in IceMesh._iceFlowMFD is treated as zero before
# the Bahr thickness H = icewe*icewf*Q^0.3. The smoothing is one implicit
# diffusion solve, which gives the discharge an exponentially decaying tail over
# every cell, and the power 0.3 inflates it (Q ~ 1e-3 m3/yr becomes 10 cm of
# ice); below ~100x the flow solver's relative tolerance (rtol 1e-8) the tail is
# not resolved by the solve, so whether a fringe cell got ice depended on solver
# noise and hence on the partition. Every cell this removes lies in the pure
# smoothing tail (zero routed discharge); on goSPL-examples glacial_erosion it
# removes 0.7% of the ice volume, all fringe thinner than ~2 m. See AGENTS.md >
# Magic numbers.
ICE_DISCHARGE_REL_FLOOR = 1.0e-6

# Mean cell size (m) above which goSPL warns that gravity-current slope erosion
# (`gravity: erosion_length`) is not meaningful. Submarine channels are a few
# hundred metres to a few kilometres wide and need several cells across to be
# resolved; at 615 m cells (the regional example) the erosion only roughened
# the toe of the slope. An estimate, not a measured threshold: no finer mesh has
# been tested (2026-11-23).
GRAVITY_EROSION_MAX_CELL = 200.0
# Tolerance (m) between an elevation and a depression-filled surface (lFill,
# waterFilled). The fill leaves 1-ulp noise between them, and that noise
# depends on the partition and the platform. A cell within this distance of the
# filled surface is AT it: routing reads lFill there (FAMesh._routingSurface);
# during an overflow pass it is not under water (_distributeDownstream, the
# original strict `<`), and at the end of flowAccumulation it belongs to the
# lake (the original inclusive `<=`). Without the tolerance a pit's spill
# node, whose elevation IS the water level, kept or lost its whole discharge
# on a 1-ulp coin flip that differed with the rank count (2026-11-24).
FILL_LEVEL_TOL = 1.0e-9


# Growth (m) per sub-step of the carbonate production integration
# (carbplex._carbSubsteps): the step is cut into ceil(sum(g_f) * dt / this)
# second-order sub-steps, from the inputs only, so every rank uses the same
# number. 0.25 m matches the exact Bosscher & Schlager ODE to 1e-4 relative
# at 1 and 10 kyr steps (benchmarks/test_carbonate.py).
CARB_SUBSTEP_GROWTH = 0.25

# Fraction of a node's size-group volume below which the routed lithology mix
# of that group is not used to split it (seaplex._marineLitShare): the node
# starts from the domain mix instead. A deep node fed only by resuspended
# pelagic rain must put the depth-biased coarse part of its deposit in the
# coarse group, but its routed coarse mix was ~1e-6 of that volume, the
# difference of two nearly equal columns: diffusion round-off then set the
# carb_sand share (0.01 vs 0.15 at np=1 vs np=3, 2026-11-25). 1e-2 left 4e-7 m
# of difference (elevation: 1.2e-7), 1e-3 left 1.5e-6 m.
LITH_PRIOR_REL_FLOOR = 1.0e-2

# Minimum discharge rate (m/yr) on the diagonal of the groundwater solute
# transport operator (gwplex._updateSolute), added to the export too so the
# budget closes. Under the sea or a lake on a flat bed the head is the surface
# and neighbours share it: the row was empty or round-off (1e-19), the direct
# LU failed and every solute field became NaN (a bathymetry clipped at a
# constant depth, 2026-11-25). Far below any recharge: well-posed rows are
# unchanged.
SOLUTE_DISCHARGE_FLOOR = 1.0e-6
