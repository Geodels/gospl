"""
Elevation-spike probe: a verbose-only diagnostic that attributes a runaway
local deposit to the stage that produced it.

A goSPL timestep raises elevation in several places (fluvial deposition,
continental cascade closed-sink deposit, pit/lake infill, marine clinoform
routing, marine diffusion, hillslope creep). Most of those are bounded by a
physical envelope (a pit's spill rim, the clinoform surface, the pre-deposition
global maximum), but a handful are deliberately unbounded so that mass is
conserved when the routing has nowhere else to put the sediment. When one of
those fires on a degenerate cell it produces a km-scale needle that is only
visible several output steps later, by which time the stage that caused it is
impossible to identify from the output alone.

``probeZ`` records the global elevation extremes at a named point in the step
and prints how far the maximum moved since the previous probe, so the offending
stage names itself in the log. ``reportVolume`` reports the volume and the peak
thickness of one deposit increment, which is what distinguishes "this stage
placed a lot of sediment" from "this stage placed sediment badly".

SCOPE — this is a **deliberately partial** debugging facility, not a complete
audit. What is instrumented:

* the per-stage elevation maximum around SPL, the continental cascade, the
  marine path, hillslope creep and flexure (``model.runProcesses``);
* every deposition site that is **deliberately unbounded**, i.e. the ones that
  can lift a cell arbitrarily far because they exist to conserve mass when the
  routing has nowhere else to put the sediment: the closed-sink deposit
  (``sedplex._distributeSediment``), the force-deposit at marine terminal sinks
  and the residual drained at cascade exit (``seaplex._distOcean``);
* the two known amplifiers/conditioners that turn a normal deposit into a
  needle: the marine-diffusion mass rescale (``hillslope._diffuseOcean``) and
  the ``fDep`` 0.99-cap population (``SPL.erodepSPL``);
* how much the pit-fill ``gmax`` guard had to clip (``sedplex._updateSinks``),
  which also shows when that guard has been made useless by an earlier
  unbounded stage having already raised ``gmax``.

What is NOT instrumented (add a probe when you need it): ice/glacial till, soil
production and creep, groundwater/duricrust, tectonics and horizontal
advection, and the ``nlSPL``/``soilSPL`` erosion flavours.

**Every call is COLLECTIVE** (``Vec.max``/``Vec.min`` and the ``allreduce`` in
``reportVolume`` wrap ``MPI_Allreduce``), so it MUST be reached by every rank.
The only gate is ``self.verbose``, a config scalar that is identical on every
rank, which keeps that safe — see AGENTS.md > MPI contract > "THE #1 parallel
deadlock class". Do NOT put a call behind a rank-local condition, and do NOT
gate it on ``MPIrank == 0``.

Cost, verbose runs only: ~6 ``Vec.max``/``Vec.min`` plus ~8 ``allreduce`` per
timestep. Negligible against a real step (tens of seconds on a multi-million
node mesh), but it is not free — that is why it is verbose-gated rather than
always on.
"""

from mpi4py import MPI

MPIrank = MPI.COMM_WORLD.Get_rank()


def probeZ(obj, tag):
    """
    Report the global elevation range and the change in the global maximum
    since the previous probe.

    :arg obj: the goSPL ``Model`` (anything exposing ``hGlobal`` and ``verbose``)
    :arg tag: short label naming the stage that just ran
    """

    # `verbose` is a YAML scalar -> identical on every rank -> safe gate for
    # the collective reductions below.
    if not getattr(obj, "verbose", False):
        return

    imax, zmax = obj.hGlobal.max()
    imin, zmin = obj.hGlobal.min()

    prev = getattr(obj, "_zprobe_last", None)
    obj._zprobe_last = zmax
    obj._zprobe_node = imax

    if MPIrank != 0:
        return

    if prev is None:
        print(
            "[zprobe] %-22s zmax %12.2f (node %d)  zmin %12.2f"
            % (tag, zmax, imax, zmin),
            flush=True,
        )
    else:
        dz = zmax - prev
        flag = "  <== RAISED THE GLOBAL MAXIMUM" if dz > 1.0 else ""
        print(
            "[zprobe] %-22s zmax %12.2f (node %d)  dzmax %+12.2f%s"
            % (tag, zmax, imax, dz, flag),
            flush=True,
        )

    return


def resetZ(obj):
    """
    Drop the running reference so the first probe of a step prints an absolute
    value rather than a delta against the previous step.

    :arg obj: the goSPL ``Model``
    """

    obj._zprobe_last = None

    return


def reportVolume(obj, tag, thickness, note=""):
    """
    Report the volume and the peak thickness of a deposit increment.

    Used at the deliberately unbounded deposition sites so the log shows both
    how much sediment they placed and how thick the worst cell got.

    **Collective** (two ``Allreduce``), gated only on ``obj.verbose``.

    :arg obj: the goSPL ``Model``
    :arg tag: short label naming the deposition site
    :arg thickness: local per-node deposit thickness (m), length ``lpoints``
    :arg note: optional extra text appended to the line
    """

    if not getattr(obj, "verbose", False):
        return

    owned = obj.inIDs == 1
    vol = MPI.COMM_WORLD.allreduce(
        float((thickness * obj.larea)[owned].sum()), op=MPI.SUM
    )
    thmax = MPI.COMM_WORLD.allreduce(
        float(thickness[owned].max()) if owned.any() else 0.0, op=MPI.MAX
    )

    if MPIrank == 0 and (vol != 0.0 or thmax != 0.0):
        print(
            "[zprobe] %-22s volume %11.4e m3  max thickness %11.2f m %s"
            % (tag, vol, thmax, note),
            flush=True,
        )

    return
