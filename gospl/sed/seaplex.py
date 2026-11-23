import os

import gc
import sys
import vtk
vtk.vtkObject.GlobalWarningDisplayOff()

import warnings
import petsc4py
from gospl.tools.petscgc import safe_garbage_cleanup
from gospl.tools.zprobe import reportVolume
import numpy as np
from scipy import spatial
from scipy.special import ndtr

from mpi4py import MPI
from time import process_time
from vtk.util import numpy_support  # type: ignore

from gospl.tools.constants import (
    BOUNDARY_FLOW_SENTINEL,
    DEPOSIT_FLOOR,
    GRAVITY_EROSION_MAX_CELL,
    MISSING_DATA_SENTINEL,
)

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import mfdrcvrs
    from gospl._fortran import epsfill
    from gospl._fortran import cellslope

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()


class SEAMesh(object):
    """
    This class encapsulates all the functions related to sediment transport and deposition in the marine environment for **river delivered sediments**.

    .. note::
        For an overview of solution to nonlinear ODE and PDE problems, one might find the online book from `Langtangen (2016) <http://hplgit.github.io/num-methods-for-PDEs/doc/pub/nonlin/pdf/nonlin-4print-A4-2up.pdf>`_ relevant.
    """

    def __init__(self, *args, **kwargs):
        """
        Initialisation of the `SEAMesh` class.
        """

        self.coastDist = None
        self.zMat = self._matrix_build_diag(np.zeros(self.lpoints))
        # Marine gravity currents (opt-in `gravity:` block): cumulative deposit
        # thickness, numpy (no PETSc object, nothing to destroy).
        if getattr(self, "gravOn", False):
            self.gravDep = np.zeros(self.lpoints, dtype=np.float64)
            if self.gravErosionLength is not None:
                # Slope erosion needs channels resolved by the mesh: warn on
                # meshes coarser than GRAVITY_EROSION_MAX_CELL (collective).
                own = self.inIDs == 1
                tot = MPI.COMM_WORLD.allreduce(
                    np.array([np.sqrt(self.larea[own]).sum(), float(own.sum())]),
                    op=MPI.SUM)
                dx = tot[0] / max(tot[1], 1.0)
                if MPIrank == 0 and dx > GRAVITY_EROSION_MAX_CELL:
                    print(
                        "[gravity] WARNING: slope erosion (erosion_length) on a "
                        "mesh of %0.0f m cells; it needs channels resolved by the "
                        "mesh (cells of %0.0f m or less) and otherwise mostly "
                        "roughens the toe of the slope." % (dx, GRAVITY_EROSION_MAX_CELL),
                        flush=True,
                    )

        return

    def _gravFloodFraction(self, sedRate, disch):
        r"""
        Fraction of the long-term river load carried in floods above the
        hyperpycnal threshold (``gravity: flood``), per node.

        Daily discharge ``X = Q/Qbar`` is log-normal with log-std ``sigma``
        (mean 1) and the load follows a rating curve ``Qs ~ Q^b``, so the
        concentration is ``C(X) = Cbar X^(b-1) / E[X^b]`` with the mean-flux
        concentration ``Cbar = rho_s Qs / Q``. The fraction of the load carried
        above ``C_c`` is

        .. math::

            f_h = \Phi\left(\frac{\mu + b\sigma^2 - \ln X_0}{\sigma}\right),
            \qquad X_0 = \left(\frac{C_c E[X^b]}{\bar C}\right)^{1/(b-1)},
            \qquad \mu = -\sigma^2/2

        (``docs/DESIGN_GRAVITY_CURRENTS.md`` §2.3). Zero where the load or the
        discharge is zero. Rank-local.

        :arg sedRate: river sediment flux at the node (m3/yr)
        :arg disch: water discharge at the node (m3/yr)

        :return: fraction in [0, 1] (local array)
        """
        b = self.gravRating
        sig = self.gravSigma
        mu = -0.5 * sig * sig
        frac = np.zeros(self.lpoints, dtype=np.float64)
        ok = (sedRate > 0.0) & (disch > 0.0)
        if not ok.any():
            return frac
        cbar = self.gravRhoS * sedRate[ok] / disch[ok]
        lnEXb = b * mu + 0.5 * b * b * sig * sig
        lnX0 = (np.log(self.gravCcrit) + lnEXb - np.log(cbar)) / (b - 1.0)
        frac[ok] = ndtr((mu + b * sig * sig - lnX0) / sig)
        return frac

    def _gravNetwork(self, hl):
        r"""
        Routing network of the gravity currents (phase G3b): sets ``self._gW``
        (the transposed receiver-weight matrix ``W^T``) and ``self._gSink`` (the
        terminal sinks of that network).

        ``gravity: network: smoothed`` reuses the clinoform cascade's matrix
        ``self.dMat1`` (smoothed bathymetry, up to 12 receivers, made to SPREAD
        river sediment). The default ``raw`` builds a separate matrix on the
        eps-filled ACTUAL bathymetry with ``gravity: receivers`` receivers and a
        flow-partition exponent ``gravity: flow_exp`` (``_gravReceivers``): on steep
        slopes the steepest path takes most of the flux (channels), on gentle
        floors the gradients are alike and the flux spreads (lobes). Same
        construction as ``_matOcean`` (serial ``epsfill`` on rank 0 from the
        owned nodes, single-pass CSR assembly); release with ``_gravRelease``.

        Collective (min reduction, gather, bcast, matrix assembly).

        :arg hl: pre-deposition local elevation
        """
        if self.gravNetwork == "smoothed":
            self._gW, self._gSink, self._gWown = self.dMat1, self.is_sink_local, False
            self._gWF, self._gSinkF = None, None
            return
        minh = max(self.hGlobal.min()[1] + 0.1, self.sealevel + self.oFill)
        surf = self._gravSmooth(hl, self.gravNetworkSmooth)
        far = None
        if self.flatModel:
            surf[self.outletIDs] = BOUNDARY_FLOW_SENTINEL
        else:
            far = self.coastDist > self.offshore
            surf[far] = BOUNDARY_FLOW_SENTINEL
        fillz = self._gatherGlobalOnRoot(surf)
        if MPIrank == 0:
            fillz = epsfill(minh, fillz)
        fillz = self._scatterFromRoot(fillz)
        if far is not None:
            fillz[far] = surf[far]
        self._gW, self._gSink = self._gravAssemble(
            fillz, self.gravReceivers, self.gravFlowExp)
        self._gWown = True
        # dual lithology with its own fine network (spreading lobes)
        self._gWF, self._gSinkF = None, None
        if self.stratLith and (self.gravReceiversFine, self.gravFlowExpFine) != (
                self.gravReceivers, self.gravFlowExp):
            self._gWF, self._gSinkF = self._gravAssemble(
                fillz, self.gravReceiversFine, self.gravFlowExpFine)
        return

    def _gravAssemble(self, fillz, nrcv, fexp):
        """Assemble one gravity network matrix ``W^T`` on the filled surface
        ``fillz`` with ``nrcv`` receivers and weights ``S^fexp``; returns the
        matrix and the sink mask. Collective (assembly)."""
        rcv, wght = self._gravReceivers(fillz, nrcv, fexp)
        nodes = np.arange(self.lpoints)
        sink = wght.sum(axis=1) <= 0.0
        if self.flatModel:
            rcv[self.outletIDs, :] = np.tile(self.outletIDs, (12, 1)).T
            wght[self.outletIDs, :] = 0.0
        if self.flatModel:
            sink[self.outletIDs] = False

        IntType = petsc4py.PETSc.IntType
        rcvFlat = rcv.astype(IntType)
        data = wght.copy()
        data[rcvFlat == nodes[:, None]] = 0.0
        indptr = np.arange(0, self.lpoints * 12 + 1, 12, dtype=IntType)
        off = self._matrix_build(nnz=(12, 12))
        off.assemblyBegin()
        off.setValuesLocalCSR(
            indptr, rcvFlat.ravel(), data.ravel(),
            addv=petsc4py.PETSc.InsertMode.ADD_VALUES,
        )
        off.assemblyEnd()
        W = self.zMat.copy()
        W.axpy(1.0, off)
        off.destroy()
        W.transpose()
        return W, sink

    def _gravReceivers(self, fillz, nrcv, fexp):
        r"""
        Receivers and weights of the gravity network on the filled surface
        ``fillz``: the ``gravity: receivers`` steepest strictly-downhill
        neighbours of each node, weighted by :math:`\mathrm{S^{p}}` (``flow_exp``)
        and normalised; a node with no lower neighbour is a sink (no weight).

        Not ``mfdrcvrs``: that kernel gives every selected receiver the SAME
        weight ``1/k`` (so the flow-partition exponent had no effect and the
        flux could not concentrate), does not renormalise when it keeps fewer
        receivers than it found (flux lost), and its ``k > receivers`` branch
        left nodes with lower neighbours as sinks. Exact slope ties break on
        the partition-invariant input id, as in ``mfdrcvrs``. Ghost rows are
        computed on a truncated neighbourhood but dropped at assembly.

        :arg fillz: eps-filled local elevation (partition invariant)
        :arg nrcv: number of receivers
        :arg fexp: flow-partition exponent (weights ``S^fexp``)
        :return: ``rcv`` (lpoints x 12, local ids, self-padded) and ``wght``
        """
        ngb = self.FVmesh_ngbID
        nodes = np.arange(self.lpoints)
        valid = ngb >= 0
        nb = np.where(valid, ngb, nodes[:, None])
        dist = np.linalg.norm(self.lcoords[nb] - self.lcoords[:, None, :], axis=2)
        drop = fillz[:, None] - fillz[nb]
        down = valid & (drop > 0.0) & (dist > 0.0)
        slp = np.where(down, drop / np.where(dist > 0.0, dist, 1.0), 0.0)
        # deterministic exact-tie break (lower input id wins), as mfdrcvrs
        key = np.where(down, slp * (1.0 - 1.0e-15 * self.gid[nb]), -1.0)
        order = np.argsort(-key, axis=1, kind="stable")
        keep = np.zeros_like(down)
        nr = min(nrcv, ngb.shape[1])
        np.put_along_axis(keep, order[:, :nr], True, axis=1)
        keep &= down
        w = np.where(keep, slp ** fexp, 0.0)
        wsum = w.sum(axis=1)
        w = np.divide(w, wsum[:, None], out=np.zeros_like(w), where=wsum[:, None] > 0.0)
        rcv = np.where(keep, nb, nodes[:, None])
        if rcv.shape[1] < 12:
            pad = 12 - rcv.shape[1]
            rcv = np.hstack([rcv, np.repeat(nodes[:, None], pad, axis=1)])
            w = np.hstack([w, np.zeros((self.lpoints, pad))])
        return rcv, w

    def _gravRelease(self):
        """Destroy the gravity network matrix if ``_gravNetwork`` built one
        (collective)."""
        if getattr(self, "_gWown", False) and self._gW is not None:
            self._gW.destroy()
        if getattr(self, "_gWF", None) is not None:
            self._gWF.destroy()
        self._gW, self._gWown, self._gWF = None, False, None
        return

    def _gravSmooth(self, arr, n):
        """``n`` passes of area-weighted averaging of a local array over each
        node and its neighbours, owner -> ghost sync after every pass (so the
        owned values do not depend on the partition). ``n = 0`` returns a copy.
        Collective when ``n > 0`` (YAML scalar: same on every rank)."""
        out = np.array(arr, dtype=np.float64, copy=True)
        if n <= 0:
            return out
        ngb = self.FVmesh_ngbID
        valid = ngb >= 0
        nb = np.where(valid, ngb, 0)
        a = self.larea
        wn = np.where(valid, a[nb], 0.0)
        den = a + wn.sum(axis=1)
        out = self._gravSync(out)
        for _ in range(n):
            out = self._gravSync((a * out + (wn * out[nb]).sum(axis=1)) / den)
        return out

    def _gravSpreadPasses(self):
        """Number of ``_gravSpread`` passes for the spreading length
        ``deposit_spread`` (m): each pass spreads a deposit by a variance of
        about half a squared cell size, so ``n = ceil(2 (L / dx)^2)`` with
        ``dx`` the mean cell size of the owned nodes (global reduction, the
        same on every rank), and the same length smooths the same width at any
        resolution. Collective."""
        own = self.inIDs == 1
        tot = MPI.COMM_WORLD.allreduce(
            np.array([np.sqrt(self.larea[own]).sum(), float(own.sum())]), op=MPI.SUM)
        dx = tot[0] / max(tot[1], 1.0)
        # default (key absent): two mean cell sizes
        length = 2.0 * dx if self.gravDepositSpread is None else self.gravDepositSpread
        return int(np.ceil(2.0 * (length / dx) ** 2))

    def _gravSpread(self, vol, n, hl, exclude=None):
        r"""
        Spread gravity deposits (m3 per node, ``(lpoints,)`` or
        ``(lpoints, k)``) laterally by ``n`` passes of a conservative
        diffusion of their thickness ``t = v / a``:

        .. math::

            v_i' = v_i + \sum_j w_{ij} (t_j - t_i), \qquad
            w_{ij} = \frac{1}{2} \frac{\min(a_i, a_j)}{\max(n_i, n_j)}

        (``n_i`` the number of neighbours). The exchange is antisymmetric, so
        volume is conserved exactly, a uniform thickness is left unchanged on
        any mesh, and ``sum_j w_ij <= a_i / 2`` keeps every pass positive.
        Only cells deeper than the bypass depth (not open outlets, not the
        eroding or failing cells ``exclude``) take
        part. The same operator for every column, so class columns still sum
        to column 0. Gives the belts of a concentrated network a width of a
        few cells. An earlier donor-share operator (``v_j`` split in
        proportion to the recipients' areas) conserved volume but not a
        uniform thickness: odd-sized cells gained or lost a little every
        step, and chains of 2-5 m dots grew along the mesh seams.

        Ghost thicknesses, areas and degrees are synced from their owners, so
        owned values are partition independent. Collective when ``n > 0``.
        """
        out = np.array(vol, dtype=np.float64, copy=True)
        if n <= 0:
            return out
        two = out.ndim == 2
        if not two:
            out = out[:, None]
        ngb = self.FVmesh_ngbID
        valid = ngb >= 0
        nb = np.where(valid, ngb, 0)
        ok = (self.sealevel - hl) > self.gravBypassDepth
        ok[self.outletIDs] = False
        if exclude is not None:
            # eroding / failing cells receive no deposit in the step they erode
            ok &= ~exclude
        okf = self._gravSync(ok.astype(np.float64)) > 0.5
        a = self._gravSync(self.larea.copy())
        deg = self._gravSync(valid.sum(axis=1).astype(np.float64))
        w = 0.5 * np.minimum(a[:, None], a[nb]) / np.maximum(deg[:, None], deg[nb])
        w = np.where(valid & okf[:, None] & okf[nb], w, 0.0)
        for _ in range(n):
            for c in range(out.shape[1]):
                t = self._gravSync(out[:, c]) / a
                out[:, c] = out[:, c] + (w * (t[nb] - t[:, None])).sum(axis=1)
        # ghost rows were computed on truncated neighbourhoods: take the owners'
        for c in range(out.shape[1]):
            out[:, c] = self._gravSync(out[:, c])
        return out if two else out[:, 0]

    def _gravSlope(self, hl):
        """Seafloor gradient used by the gravity tapers and erosion:
        ``cellslope`` averaged ``slope_smooth`` times (``_gravSmooth``), so the
        deposition and erosion follow the slope at the scale of a few cells,
        not the cell-to-cell noise. Collective when smoothing."""
        # synced first: a ghost's cellslope uses its truncated neighbourhood
        return self._gravSmooth(self._gravSync(cellslope(self.lpoints, hl)),
                                self.gravSlopeSmooth)

    def _gravSync(self, arr):
        """Owner -> ghost copy of a local array (collective; scratch
        ``self.tmpL``, ``self.tmp1``)."""
        self.tmpL.setArray(arr)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        return self.tmpL.getArray().copy()

    def _gravLambda(self, runout):
        r"""
        Deposited fraction of the gravity flux at each node,
        ``Lambda = 1 - exp(-sqrt(cell area)/runout)``, tapered where the flow
        bypasses: by ``clip((d - d_b)/d_b, 0, 1)`` with ``d`` the water depth and
        ``d_b = bypass_depth`` (wave base: no deposition on the shelf, full
        deposition from twice that depth) and, when set, by
        ``1/(1 + (S/S_b)^4)`` with ``S`` the gradient and ``S_b = bypass_slope``
        (one half at ``S_b``: the flow crosses the foresets and deposits at
        their toe). Smooth tapers, not switches, so that noisy gradients do not
        make neighbouring cells alternate between deposit and bypass. Terminal
        sinks of the gravity network (``self._gSink``) keep what reaches them
        (1) and open outlets export it (0). Synced owner -> ghost (the gradient
        of a ghost node, computed on its truncated neighbourhood, is not the
        owner's).

        Scratch Vecs: ``self.tmpL``, ``self.tmp1``. Collective.

        :arg runout: run-out length (m)

        :return: local array of deposited fractions in [0, 1]
        """
        lam = 1.0 - np.exp(-np.sqrt(self.larea) / runout)
        hl = self.hLocal.getArray().copy()
        # Smooth tapers, not on/off switches: a hard threshold on a gradient
        # that is noisy from cell to cell (or on the depth) made neighbouring
        # cells alternate between deposit and bypass, and the deposits then
        # flattened their own cells, giving a speckled, lumpy deposit.
        if self.gravBypassDepth > 0.0:
            depth = self.sealevel - hl
            lam *= np.clip((depth - self.gravBypassDepth) / self.gravBypassDepth, 0.0, 1.0)
        if self.gravBypassSlope is not None:
            ratio = self._gravSlope(hl) / self.gravBypassSlope
            lam *= 1.0 / (1.0 + ratio ** 4)
        lam[self._gSink] = 1.0
        lam[self.outletIDs] = 0.0
        return self._gravSync(lam)

    def _gravEps(self, hl):
        r"""
        Slope-erosion coefficient of the gravity currents (phase G3b): the
        fraction of the flux picked up at each node,

        .. math::

            \varepsilon_i = \frac{\ell_i}{L_e}\,\frac{(S_i/S_e)^4}{1 + (S_i/S_e)^4}

        with the cell size ``l_i``, ``erosion_length`` ``L_e`` and
        ``erosion_slope`` ``S_e`` (0 on flat ground, one half at ``S_e``), times
        the depth taper of ``_gravLambda`` (no erosion on the shelf above the
        bypass depth). Zero at the network's sinks and at open outlets. Synced
        owner -> ghost. Collective.

        The gradient is ``_gravSlope`` (smoothed ``slope_smooth`` times).

        :arg hl: pre-deposition local elevation

        :return: local array of erosion coefficients (>= 0)
        """
        slope = self._gravSlope(hl)
        ratio = (slope / self.gravErosionSlope) ** 4
        eps = np.sqrt(self.larea) / self.gravErosionLength * ratio / (1.0 + ratio)
        # eps <= 1: a current at most doubles its load in one cell. With an
        # erosion_length shorter than the cell size, eps reached ~18 per cell
        # on a 360 km global mesh, the flux grew as prod(1 + eps) along a path
        # (1e22) and every erosion solve failed and was discarded.
        eps = np.minimum(eps, 1.0)
        if self.gravBypassDepth > 0.0:
            depth = self.sealevel - hl
            eps *= np.clip((depth - self.gravBypassDepth) / self.gravBypassDepth, 0.0, 1.0)
        eps[self._gSink] = 0.0
        eps[self.outletIDs] = 0.0
        return self._gravSync(eps)

    def _gravSolve(self, factor, col):
        r"""
        One routed solve ``(I - W^T diag(factor)) q = col`` on the gravity
        network, ``factor = 1 - Lambda`` (+ ``eps`` with slope erosion). Shared
        routing KSP (``_solve_KSP(route=True)``); returns the local flux and the
        volume ponded by a failed solve (or None). Scratch ``self.tmp``,
        ``self.tmp1``, ``self.tmpL``. Collective.

        Two guards keep the gravity solves from changing the river routing:
        a supply that is zero on every rank returns ``q = 0`` without a solve
        (the routing solver took a zero right-hand side for a failed solve and
        ran its exact rescue), and the sticky exact-rescue switch of the flow
        solver (``_ida_exact``) is restored after the solve. Before these,
        ``fraction: 0`` (or a step with no supply) turned the switch on and
        every later river/sediment routing solve of the run used the exact
        path: ulp-level changes that broke bitwise identity on Linux.
        """
        own = self.inIDs == 1
        cmax = MPI.COMM_WORLD.allreduce(
            float(np.abs(col[own]).max()) if own.any() else 0.0, op=MPI.MAX)
        if cmax == 0.0:
            return np.zeros(self.lpoints, dtype=np.float64), None
        exact = getattr(self, "_ida_exact", None)
        self.tmpL.setArray(factor)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        M = self._gW.copy()
        M.diagonalScale(R=self.tmp1)
        M.scale(-1.0)
        M.shift(1.0)
        self.tmpL.setArray(col)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self._solve_KSP(False, M, self.tmp, self.tmp1, route=True)
        M.destroy()
        if exact is not None:
            self._ida_exact = exact
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        q = np.maximum(self.tmpL.getArray().copy(), 0.0)
        return q, self._takeRoutingPond()

    def _gravFailure(self, hl):
        r"""
        Slope-failure trigger of the gravity currents (phase G3). On marine
        cells below the bypass depth whose gradient :math:`S` (``_gravSlope``,
        smoothed ``slope_smooth`` times) exceeds ``failure_slope``
        :math:`S_f`, the sediment fails down to the critical slope:

        .. math::

            \Delta z_i = T_d(d_i)\, \min\left(\frac{(S_i - S_f)\,\ell_i}{2},\;
            z^{max}_f,\; D_i\right)

        with :math:`\ell_i` the cell size, :math:`z^{max}_f` ``max_failure``
        (m per step), :math:`D_i` the sediment deposited during the run
        (positive cumulative erosion-deposition, so bedrock never fails) and
        :math:`T_d` the depth taper of the deposition. Lowering a cell by
        :math:`(S - S_f)\ell/2` takes about :math:`S - S_f` off its gradient.
        The failed volume becomes a gravity-current source at that cell, with
        the composition removed from the stratigraphy, exactly like the slope
        erosion. Zero at the network's sinks and at open outlets. Collective
        (slope smoothing).

        :arg hl: pre-deposition local elevation

        :return: failed bulk volume (m3, local) and the failing-node mask
        """
        S = self._gravSlope(hl)
        dz = 0.5 * np.maximum(S - self.gravFailureSlope, 0.0) * np.sqrt(self.larea)
        avail = np.maximum(self.cumEDLocal.getArray(), 0.0)
        dz = np.minimum(np.minimum(dz, self.gravMaxFailure), avail)
        depth = self.sealevel - hl
        db = self.gravBypassDepth
        if db > 0.0:
            dz *= np.clip((depth - db) / db, 0.0, 1.0)
        else:
            dz[depth <= 0.0] = 0.0
        dz[self._gSink] = 0.0
        dz[self.outletIDs] = 0.0
        dz = self._gravSync(dz)
        vol = dz * self.larea
        return vol, vol > 0.0

    def _gravErosion(self, supply, hl):
        r"""
        Slope erosion by the gravity currents (phase G3b,
        ``docs/DESIGN_GRAVITY_CURRENTS.md`` 6c). A node ERODES where the
        pick-up coefficient exceeds the deposited fraction (``eps > Lambda`` of
        the coarse or single fraction) and then receives no deposit
        (each node erodes or deposits in a step, never both, so the
        stratigraphy is not counted twice); elsewhere it only deposits. The
        total flux with pick-up solves

        .. math::

            (I - W^T \mathrm{diag}(1 - \Lambda + \varepsilon))\, Q = S

        and the eroded bulk volume is ``e = eps Q``, at most
        ``max_erosion x cell area`` (``eps`` is lowered where exceeded and the
        solve repeated, at most 20 passes; this also guards against a runaway).
        The eroded volume then enters the deposition routing as an explicit
        source, so conservation is exact whatever the coupling error.

        With ``erosion_exp: a > 1`` the pick-up grows with the flux,
        ``eps = eps0 (Q / Q_e)^(a-1)`` (``erosion_flux`` ``Q_e`` in m3/yr), so
        the strongest path erodes most and captures the flow (channel
        selection); ``Q`` is the flux of the linear solve, evaluated once (a
        fixed point within the step diverges), and ``eps <= 1``.

        Scratch ``self.tmp``, ``self.tmp1``, ``self.tmpL``. Collective.

        :arg supply: total gravity supply (m3 over the step, local array)
        :arg hl: pre-deposition local elevation

        :return: eroded bulk volume (m3, local) and the eroding-node mask
        """
        eps0 = self._gravEps(hl)
        lam0 = self._gravLambda(self.gravRunout)
        ero = eps0 > lam0
        eps = np.where(ero, eps0, 0.0)
        lam = np.where(ero, 0.0, lam0)
        Q, _ = self._gravSolve(1.0 - lam + eps, supply)
        if self.gravErosionExp != 1.0:
            # Non-linear flux dependence eps = eps0 (Q / Q_e)^(a-1), with Q the
            # flux of the LINEAR solve above, evaluated once. Iterating on Q
            # within a step diverges: the pick-up multiplies the flux along a
            # path by prod(1 + eps) and eps grows with it (Q reached 1e10 Q_e
            # in 4 passes on the regional margin, then the solve failed). The
            # channel feedback acts between steps, through the topography.
            # eps <= 1: a current at most doubles its load in one cell.
            qe = self.gravErosionFlux * self.dt
            epsk = eps0 * (np.maximum(Q, 0.0) / qe) ** (self.gravErosionExp - 1.0)
            epsk = np.minimum(epsk, 1.0)
            ero = epsk > lam0
            eps = self._gravSync(np.where(ero, epsk, 0.0))
            lam = np.where(ero, 0.0, lam0)
            Q, _ = self._gravSolve(1.0 - lam + eps, supply)
        cap = self.gravMaxErosion * self.larea
        owned = self.inIDs == 1
        for _ in range(20):
            over = ero & (eps * Q > cap * (1.0 + 1.0e-9)) & (Q > 0.0)
            more = MPI.COMM_WORLD.allreduce(bool((over & owned).any()), op=MPI.LOR)
            if not more:
                break
            eps = eps.copy()
            eps[over] = cap[over] / Q[over]
            eps = self._gravSync(eps)
            Q, _ = self._gravSolve(1.0 - lam + eps, supply)
        e = eps * Q
        e[~ero] = 0.0
        return e, ero

    def _gravRoute(self, supply, runout, room=None, nodep=None):
        r"""
        Route gravity-current supplies down the gravity network with a
        run-out length, ``(I - W^T diag(1 - Lambda)) q = s``, deposit
        ``Lambda q`` (``docs/DESIGN_GRAVITY_CURRENTS.md`` §2.1-2.2).

        ``W^T`` is the gravity network matrix ``self._gW`` (``_gravNetwork``).
        ``Lambda`` comes from ``_gravLambda`` (run-out length, bypass on the
        shelf and, optionally, on steep slopes, 1 at terminal sinks, 0 at open
        outlets) and is 0 at the eroding nodes ``nodep`` (phase G3b).

        Accommodation: with ``room`` (m3 per node), a node cannot receive more
        than its room. The first column (the total supply) is solved, every
        node that overfills gets ``Lambda = room / q`` (it deposits exactly its
        room and passes the rest on) and the solve is repeated until nothing
        overfills (``Lambda`` only decreases, so a few passes suffice; at most
        20). Terminal sinks keep everything whatever their room (mass must stay
        in the domain). The remaining columns (provenance classes) are then
        solved with the final ``Lambda``, so they sum to the total.

        Scratch Vecs: overwrites ``self.tmp``, ``self.tmp1`` and ``self.tmpL``
        (all free at this point of ``seaChange``). Collective.

        :arg supply: gravity supply per node (m3 over the step), local array of
            shape ``(lpoints,)`` or ``(lpoints, k)`` (column 0 = the total)
        :arg runout: run-out length (m)
        :arg room: accommodation per node (m3, local array) or None
        :arg nodep: eroding-node mask (no deposit there) or None

        :return: deposited and exported volumes (m3, local arrays of the same
            shape as ``supply``; export is non-zero at outlets only)
        """
        single = np.ndim(supply) == 1
        sup = np.asarray(supply, dtype=np.float64)
        if single:
            sup = sup[:, None]
        lam = self._gravLambda(runout)
        if nodep is not None:
            lam[nodep] = 0.0
        owned = self.inIDs == 1

        q, pond = self._gravSolve(1.0 - lam, sup[:, 0])
        if room is not None:
            cap = np.maximum(room, 0.0)
            capped = ~self._gSink
            for _ in range(20):
                over = capped & (lam * q > cap * (1.0 + 1.0e-9) + 1.0e-9) & (q > 0.0)
                more = MPI.COMM_WORLD.allreduce(bool((over & owned).any()), op=MPI.LOR)
                if not more:
                    break
                lam = lam.copy()
                lam[over] = cap[over] / q[over]
                # sync owner -> ghost so every rank scales the same columns
                lam = self._gravSync(lam)
                q, pond = self._gravSolve(1.0 - lam, sup[:, 0])

        dep = np.zeros_like(sup)
        exported = np.zeros_like(sup)
        for k in range(sup.shape[1]):
            if k > 0:
                q, pond = self._gravSolve(1.0 - lam, sup[:, k])
            dep[:, k] = lam * q
            exported[self.outletIDs, k] = q[self.outletIDs]
            if pond is not None:
                dep[:, k] += pond
        if single:
            return dep[:, 0], exported[:, 0]
        return dep, exported

    def _globalCoastsTree(self, coastXYZ, k_neighbors=1):
        """
        This function takes all local coastline points and computes locally the distance of all marine points to the coastline.

        :arg coastXYZ: local coastline coordinates
        :arg k_neighbors: number of nodes to use when querying the kd-tree
        """

        label_offset = np.zeros(MPIsize + 1, dtype=int)
        label_offset[MPIrank + 1] = len(coastXYZ)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, label_offset, op=MPI.MAX)
        offset = np.cumsum(label_offset)[:-1]
        totpts = np.sum(label_offset)

        # Edge case: no rank has any coastline points (all-continental run or
        # all marine), so the cKDTree below would raise on an empty input.
        if totpts == 0:
            return

        coastpts = np.zeros(((totpts) * 3,), dtype=np.float64)
        MPI.COMM_WORLD.Allgatherv(
            [coastXYZ, MPI.DOUBLE],
            [coastpts, label_offset[1:] * 3, offset * 3, MPI.DOUBLE],
        )

        # Get coastlines points globally
        tree = spatial.cKDTree(coastpts.reshape((totpts, 3)), leafsize=10)
        self.coastDist[self.seaID], _ = tree.query(
            self.lcoords[self.seaID, :], k=k_neighbors
        )

        if self.memclear:
            del tree, coastpts, offset, totpts, label_offset
            gc.collect()

        return

    def _distanceCoasts(self, data, k_neighbors=1):
        """
        This function computes for every marine vertices the distance to the closest coastline. It calls the private function:

        - _globalCoastsTree

        .. important::

            The calculation takes advantage of the `vtkContourFilter` function from VTK library which is performed on the **global** VTK mesh. Once the coastlines have been extracted, the distances are obtained by querying a kd-tree (initialised with the coastal nodes) for marine vertices contained within each partition.

        :arg data: local elevation Numpy Array
        :arg k_neighbors: number of nodes to use when querying the kd-tree
        """

        t0 = process_time()

        self.coastDist = np.zeros(self.lpoints)
        pointData = self.vtkMesh.GetPointData()
        array = numpy_support.numpy_to_vtk(data, deep=1)
        array.SetName("elev")
        pointData.AddArray(array)
        self.vtkMesh.SetFieldData(pointData)

        cf = vtk.vtkContourFilter()
        cf.SetInputData(self.vtkMesh)
        cf.SetValue(0, self.sealevel)
        cf.SetInputArrayToProcess(0, 0, 0, 0, "elev")
        cf.GenerateTrianglesOff()
        cf.Update()
        if cf.GetOutput().GetPoints() is not None:
            coastXYZ = numpy_support.vtk_to_numpy(cf.GetOutput().GetPoints().GetData())
        else:
            coastXYZ = np.zeros((0, 3), dtype=np.float64)

        # Get coastlines points globally
        self._globalCoastsTree(coastXYZ, k_neighbors)

        if self.memclear:
            del array, pointData, cf, coastXYZ
            gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Construct Distance to Coast (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _matOcean(self):
        """
        This function builds from neighbouring slopes the downstream directions in the marine environment. It calls a Fortran subroutine that locally computes for each vertice:

        - the indices of receivers (downstream) nodes depending on the desired number of flow directions (SFD to MFD).
        - the distances to the receivers based on mesh resolution.
        - the associated weights calculated based on the number of receivers and proportional to the slope.

        From these downstream directions, a local ocean downstream matrix is computed.
        """

        # Define multiple flow directions for filled + eps elevations.
        #
        # `oFill` is a depth BELOW SEA LEVEL, so the cut-off handed to `epsfill`
        # is `sealevel + oFill` — exactly as `pitfilling.fillElevation` and
        # `fillIceElevation` compute it. This site used to compare the raw
        # `self.oFill` as an ABSOLUTE elevation, which is the same number only
        # when the sea level is 0 (every Earth example, fixture and benchmark)
        # and diverges by `sealevel` otherwise, so one YAML key meant two
        # different things in two places.
        #
        # It matters because `epsfill` only eps-fills the cells AT OR ABOVE its
        # cut-off (every cell below it is pre-flagged in the seeding loop and
        # never revisited). Putting the cut-off above sea level therefore skips
        # the entire marine domain: every closed bathymetric pocket on the shelf
        # survives into the routing surface, and sediment gets trapped at the
        # coast instead of routing basinward. On a run with a deep datum
        # (`sea: position: -2200` with `oFill: -1500`) the absolute reading gave
        # a cut-off of -1500 m, i.e. 700 m ABOVE sea level.
        hl = self.hLocal.getArray().copy()
        minh = self.hGlobal.min()[1] + 0.1
        minh = max(minh, self.sealevel + self.oFill)
        # Recorded so the contract is assertable from a test (and visible when
        # debugging a marine-routing problem); nothing reads it in the model.
        self._oceanFillLevel = minh
        if not self.flatModel:
            # Only consider filleps in the first kms offshore
            hsmth = self._hillSlope(smooth=2)
            hsmth[self.coastDist > self.offshore] = BOUNDARY_FLOW_SENTINEL
            # The cut-off must come from the surface that is FILLED (the
            # smoothed one), not the raw elevation: smoothing raises the
            # deepest cells, so `hGlobal.min() + 0.1` sat below every smoothed
            # value, no cell seeded `epsfill`, and it silently returned the
            # surface unchanged. The marine routing then ran on unfilled
            # bathymetry with hundreds of terminal sinks, and river loads were
            # force-deposited km-thick in the first shallow pit they reached
            # (fixed 2026-11-12). The deepest smoothed cell now always seeds.
            owned = self.inIDs == 1
            smin = MPI.COMM_WORLD.allreduce(
                float(hsmth[owned].min()) if owned.any() else np.inf, op=MPI.MIN
            )
            minh = max(smin + 0.1, self.sealevel + self.oFill)
            self._oceanFillLevel = minh
        else:
            hsmth = hl.copy()
            hsmth[self.outletIDs] = BOUNDARY_FLOW_SENTINEL

        # The filled + eps is done on the global grid, serially on rank 0, so
        # assemble it there from the owned nodes only (hsmth is ghost-synced
        # via _hillSlope's globalToLocal, so owned values are exact).
        fillz = self._gatherGlobalOnRoot(hsmth)
        if MPIrank == 0:
            fillz = epsfill(minh, fillz)
        # Send each processor its local rows (not the whole mesh)
        fillz = self._scatterFromRoot(fillz)
        if not self.flatModel:
            fillz[self.coastDist > self.offshore] = hl[self.coastDist > self.offshore]
        # Pass `self.gid` for deterministic exact-tie-break on slope —
        # see fortran/functions.F90:mfdrcvrs.
        rcv, _, wght = mfdrcvrs(
            12, self.flowExp, fillz, BOUNDARY_FLOW_SENTINEL, self.gid,
        )

        # Set borders nodes
        if self.flatModel:
            rcv[self.outletIDs, :] = np.tile(self.outletIDs, (12, 1)).T
            wght[self.outletIDs, :] = 0.0

        # Detect terminal sinks: cells whose every flow direction points
        # back to themselves. After the dMat1 transpose below the column
        # at each such cell is identically zero, so dMat1.mult would
        # annihilate any sediment that lands there. _distOcean uses this
        # mask to force-deposit at sinks instead of routing through dMat1.
        # On flat models the open-boundary cells were artificially marked
        # self-draining just above; exclude them here so the existing
        # vdep[idBorders]=0 rule keeps handling outflux at the edge.
        nodes_idx = np.arange(self.lpoints)
        self.is_sink_local = np.all(
            rcv.astype(int) == nodes_idx[:, None], axis=1
        )
        if self.flatModel:
            self.is_sink_local[self.outletIDs] = False

        # Define downstream matrix based on filled + dir elevations
        self.dMat1 = self.zMat.copy()
        if not self.flatModel and self.Gmar > 0.:
            self.dMat2 = self.iMat.copy()
        nodes = np.arange(self.lpoints, dtype=petsc4py.PETSc.IntType)

        # Assemble all 12 flow directions in a single CSR pass instead of 12
        # build+axpy iterations. ADD_VALUES makes receivers shared by several
        # directions sum exactly as the per-direction axpy loop did, so the
        # assembled matrix is identical; this just replaces 12 matrix builds +
        # 12 sparse axpys (the dominant cost of this step) with one build and
        # one (or two) axpy.
        rcvFlat = rcv.astype(petsc4py.PETSc.IntType)
        data = wght.copy()
        data[rcvFlat == nodes[:, None]] = 0.0
        indptr = np.arange(0, self.lpoints * 12 + 1, 12, dtype=petsc4py.PETSc.IntType)
        offMat = self._matrix_build(nnz=(12, 12))
        offMat.assemblyBegin()
        offMat.setValuesLocalCSR(
            indptr,
            rcvFlat.ravel(),
            data.ravel(),
            addv=petsc4py.PETSc.InsertMode.ADD_VALUES,
        )
        offMat.assemblyEnd()
        self.dMat1.axpy(1.0, offMat)
        if not self.flatModel and self.Gmar > 0.:
            self.dMat2.axpy(-1.0, offMat)
        offMat.destroy()

        if self.memclear:
            del data, indptr, nodes
            del hl, fillz, rcv, wght
            gc.collect()

        # Store flow direction matrix
        self.dMat1.transpose()
        if not self.flatModel and self.Gmar > 0.:
            self.dMat2.transpose()

        return

    def _depMarineSystem(self, sedflux):
        r"""
        Setup matrix for the marine sediment deposition.

        The upstream incoming sediment flux is obtained from the total river sediment flux :math:`\mathrm{Q_{t_i}}` where:

        .. math::

            \mathrm{Q_{t_i}^{t+\Delta t} - \sum_{ups} w_{i,j} Q_{t_u}^{t+\Delta t}}= \mathrm{(\eta_i^{t} - \eta_i^{t+\Delta t}) \frac{\Delta t}{\Omega_i}}

        which gives:

        .. math::

            \mathrm{Q_{s_i}} = \mathrm{Q_{t_i}} - \mathrm{(\eta_i^{t} - \eta_i^{t+\Delta t}) \frac{\Delta t}{\Omega_i}}

        And the evolution of marine elevation is based on incoming sediment flux:

        .. math::

            \mathrm{\frac{\eta_i^{t+\Delta t}-\eta_i^t}{\Delta t}} = \mathrm{G{_m} Q_{s_i} / \Omega_i}

        This system of coupled equations is solved implicitly using PETSc by assembling the matrix and vectors using the nested submatrix and subvectors and by using the ``fieldsplit`` preconditioner combining two separate preconditioners for the collections of variables.

        :arg sedflux: incoming marine sediment volumes

        :return: volDep (the deposited volume of the distributed sediments)
        """

        hl = self.hLocal.getArray()
        fDepm = np.full(self.lpoints, self.Gmar)
        fDepm[fDepm > 0.99] = 0.99
        fDepm[hl > self.sealevel] = 0.

        # Define submatrices
        A00 = self._matrix_build_diag(-fDepm)
        A00.axpy(1.0, self.iMat)
        A01 = self._matrix_build_diag(-fDepm * self.dt / self.larea)
        A10 = self._matrix_build_diag(self.larea / self.dt)

        # Assemble the matrix for the coupled system
        mats = [[A00, A01], [A10, self.dMat2]]
        sysMat = petsc4py.PETSc.Mat().createNest(mats=mats, comm=MPIcomm)
        sysMat.assemblyBegin()
        sysMat.assemblyEnd()

        # Clean up. 
        A00.destroy()
        A01.destroy()
        A10.destroy()
        self.dMat2.destroy()

        # Create nested vectors
        self.tmpL.setArray(1. - fDepm)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.tmp.pointwiseMult(self.tmp, self.hGlobal)

        self.tmpL.setArray(sedflux / self.dt)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        self.h.pointwiseMult(self.hGlobal, self.areaGlobal)
        self.h.scale(1. / self.dt)
        self.tmp1.axpy(1., self.h)

        rhs_vec = petsc4py.PETSc.Vec().createNest([self.tmp, self.tmp1], comm=MPIcomm)
        rhs_vec.setUp()
        hq_vec = rhs_vec.duplicate()

        # Define solver and precondition conditions
        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        # Options prefix so the coupled marine-deposition solve can be tuned from
        # the command line without editing code, e.g.
        # `-marine_dep_pc_fieldsplit_type additive`.
        ksp.setOptionsPrefix("marine_dep_")
        ksp.setType(petsc4py.PETSc.KSP.Type.TFQMR)
        ksp.setOperators(sysMat)
        # Bound the iteration count. The Schur fieldsplit converges in a handful
        # of outer iterations on a well-conditioned step, but with no cap a
        # badly-conditioned step let TFQMR grind to PETSc's default 10000 -- one
        # such solve was the P=48 `sea`-phase blow-up (760 s) in the rank sweep.
        # Capping low makes a stalling step fail FAST; divergence is already
        # logged-and-continued below, so deposition that step is at worst
        # approximate rather than catastrophically slow.
        ksp.setTolerances(rtol=self.rtol, max_it=500)

        pc = ksp.getPC()
        pc.setType("fieldsplit")
        nested_IS = sysMat.getNestISs()
        pc.setFieldSplitIS(('h', nested_IS[0][0]), ('q', nested_IS[0][1]))

        # Schur-complement fieldsplit. This marine-deposition system has the same
        # two-field (elevation h / sediment-flux q) coupled structure as the
        # detachment-limited SPL solve (see eroder/SPL._coupledEDSystem): the h
        # and q blocks are near-triangular M-matrices the per-block ILU solves
        # almost exactly, so a Schur factorisation captures the full
        # deposition<->routing coupling and converges in a handful of outer
        # iterations -- versus ~100+ for the default additive (block-diagonal)
        # split, which only sees part of the coupling. This changes only the
        # preconditioner, so the solution (to `rtol`) is unchanged. Defaults are
        # injected into the options DB guarded by `hasName`, so PETSC_OPTIONS /
        # `-marine_dep_*` still override.
        opts = petsc4py.PETSc.Options()
        for key, val in (
            ("marine_dep_pc_fieldsplit_type", "schur"),
            ("marine_dep_pc_fieldsplit_schur_fact_type", "full"),
            ("marine_dep_pc_fieldsplit_schur_precondition", "selfp"),
            # Block solves. They MUST go through the options database: the
            # `subksps[i].setType("preonly")` calls below run before
            # `setFromOptions` switches the split to Schur, so they never took
            # effect and both blocks ran GMRES (rtol 1e-5, up to 10,000 its)
            # inside every outer iteration. On the routing (q) block, with long
            # drainage chains, that cost grew step after step (~300 s -> ~5,000
            # s per step on a 23.6M-node mesh, for the river SPL solve; fixed
            # 2026-11-12). h: one ILU application per rank block (GASM).
            ("marine_dep_fieldsplit_h_ksp_type", "preonly"),
            ("marine_dep_fieldsplit_h_pc_type", "gasm"),
            ("marine_dep_fieldsplit_q_ksp_type", "preonly"),
            # Marine routing uses up to 12 flow directions, so exact LU fills in
            # heavily (27 s symbolic on a 655k-node sphere) and GASM's set-up
            # is as slow (35 s); block-Jacobi ILU solves it in 0.29 s, 6 outer
            # iterations. (The river SPL q block keeps exact LU: 5 directions.)
            ("marine_dep_fieldsplit_q_pc_type", "bjacobi"),
            ("marine_dep_fieldsplit_q_sub_pc_type", "ilu"),
            ("marine_dep_fieldsplit_q_sub_pc_factor_shift_type", "nonzero"),
        ):
            if not opts.hasName(key):
                opts.setValue(key, val)

        subksps = pc.getFieldSplitSubKSP()
        subksps[0].setType("preonly")
        subksps[0].getPC().setType("gasm")
        subksps[1].setType("preonly")
        subksps[1].getPC().setType("gasm")

        # Apply the Schur defaults (and any `-marine_dep_*` overrides).
        ksp.setFromOptions()

        ksp.solve(rhs_vec, hq_vec)
        r = ksp.getConvergedReason()
        # Note: on KSP divergence we log and continue. Marine deposition this
        # step may be inaccurate; operators should monitor this warning.
        if r < 0:
            KSPReasons = self._make_reasons(petsc4py.PETSc.KSP.ConvergedReason())
            if MPIrank == 0:
                print(
                    "Linear solver for marine deposition failed to converge after iterations",
                    ksp.getIterationNumber(),
                    flush=True,
                )
                print("with reason: ", KSPReasons[r], flush=True)
        else:
            if MPIrank == 0 and self.verbose:
                print(
                    "Linear solver for marine deposition converge after %d iterations"
                    % ksp.getIterationNumber(),
                    flush=True,
                )

        # Update the solution
        self.newH = hq_vec.getSubVector(nested_IS[0][0])

        # Clean up
        subksps[0].destroy()
        subksps[1].destroy()
        nested_IS[0][0].destroy()
        nested_IS[1][0].destroy()
        nested_IS[0][1].destroy()
        nested_IS[1][1].destroy()
        pc.destroy()
        ksp.destroy()
        sysMat.destroy()
        hq_vec.destroy()
        rhs_vec.destroy()

        # Get the marine deposition volume
        self.tmp.waxpy(-1.0, self.hGlobal, self.newH)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        volDep = self.tmpL.getArray().copy() * self.larea
        volDep[volDep < 0] = 0.
        safe_garbage_cleanup()

        return volDep

    def _distOcean(self, sedflux, provFlux=None):
        """
        Based on the incoming marine volumes of sediment and maximum clinoforms slope we distribute locally sediments downslope.

        In-model provenance (B2b-marine, lockstep): when ``provFlux`` is given
        (the per-class incoming marine flux, shape ``(lpoints, provNb)`` with
        ``Σ_c provFlux == sedflux``), each class's sub-flux is routed through the
        **same** ``dMat1`` cascade in lockstep with the total and deposited in
        the **same proportion** at every node (capacity ``marVol`` is shared, set
        by the total flux; provenance is a passive label, so the deposited and
        continuing fractions of a node's arriving flux carry its arriving
        composition unchanged). This yields a spatially-resolved per-node
        deposited composition ``self._marDepProv`` (shape ``(lpoints, provNb)``,
        ``Σ_c == vdep`` exactly) — the source that actually drained to each
        offshore cell — instead of the old domain-wide average. This is the
        marine analogue of the continental ``_moveDownstream`` lockstep.

        :arg sedflux: incoming marine sediment volumes
        :arg provFlux: optional per-class incoming marine flux (lpoints, provNb)

        :return: vdep (the deposited volume of the distributed sediments)
        """

        marVol = self.maxDepQs.copy()
        sinkVol = sedflux.copy()
        vdep = np.zeros(self.lpoints, dtype=float)
        # Diagnostic accumulator for the two UNBOUNDED escapes in this routine
        # (the clinoform capacity `marVol` bounds everything else): the
        # force-deposit at terminal sinks and the residual drained at loop exit.
        # Both can pile an arbitrary thickness onto a single degenerate cell.
        sinkForced = np.zeros(self.lpoints, dtype=float)

        # Provenance lockstep state: one routed sub-flux per source class,
        # carried in parallel with the total `sinkVol`/`self.tmp`. `pvCur[c]`
        # holds class c's "flux to route this pass" (the analogue of self.tmp);
        # `vdepP[:, c]` accumulates its deposited volume (Σ_c vdepP == vdep).
        prov = provFlux is not None
        if prov:
            nC = self.provNb
            sinkVolP = provFlux.copy()
            vdepP = np.zeros((self.lpoints, nC), dtype=float)
            pvCur = [self.tmp.duplicate() for _ in range(nC)]
            pvNxt = self.tmp.duplicate()

        # Drain initial sediment at terminal sinks BEFORE the routing loop.
        # dMat1 has zero columns at sink cells (rcv==self for every flow
        # direction, see _matOcean), so the very first dMat1.mult below
        # would annihilate any sedflux that landed directly at a sink. On
        # global sphere fixtures this is the dominant mass-loss path:
        # vSed naturally accumulates at deep ocean basins (sinks are
        # sediment attractors), so a large fraction of the marine supply
        # arrives at sinks in the very first iteration. Without this
        # pre-drain ~25% of total redistributed sediment vanishes per run
        # (regression: test_mass_conservation; AGENTS.md Known bugs entry
        # dated 2026-06). The in-loop force-deposit below additionally
        # catches sediment that arrives at sinks from upstream during
        # later iterations.
        if self.is_sink_local.any():
            sink_mass = sinkVol[self.is_sink_local]
            if (sink_mass != 0).any():
                vdep[self.is_sink_local] += sink_mass
                sinkForced[self.is_sink_local] += sink_mass
                sinkVol[self.is_sink_local] = 0.0
                if prov:
                    vdepP[self.is_sink_local, :] += sinkVolP[self.is_sink_local, :]
                    sinkVolP[self.is_sink_local, :] = 0.0

        self.tmpL.setArray(sinkVol)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        if prov:
            for c in range(nC):
                self.tmpL.setArray(sinkVolP[:, c])
                self.dm.localToGlobal(self.tmpL, pvCur[c])

        # Convergence threshold: relative to initial input volume (1e-6 of
        # total) with a fixed floor of 1 m^3. Scales with input magnitude so
        # continental-scale runs converge in far fewer iterations than the
        # previous fixed `> 1.0 m^3` rule, while small grids still terminate
        # cleanly. A max-iteration cap protects against pathological cases.
        sumExcess = self.tmp.sum()
        initial_total = sumExcess
        threshold = max(1.0, 1.0e-6 * initial_total)
        max_iters = 5000

        step = 0
        while sumExcess > threshold:
            if step >= max_iters:
                if MPIrank == 0:
                    print(
                        "  --- Marine routing did not converge after %d "
                        "iterations (residual %.3e m^3); continuing."
                        % (max_iters, sumExcess),
                        flush=True,
                    )
                break

            # Move to downstream nodes
            self.dMat1.mult(self.tmp, self.tmp1)
            self.dm.globalToLocal(self.tmp1, self.tmpL)

            # In case there is too much sediment coming in
            sinkVol = self.tmpL.getArray().copy()

            # Provenance: route each class's flux through the SAME dMat1 (linear,
            # so Σ_c stays == the total) and snapshot the per-node arriving total
            # `arr` before the deposition block consumes `sinkVol`.
            if prov:
                arr = sinkVol.copy()
                sinkVolP = np.zeros((self.lpoints, nC), dtype=float)
                for c in range(nC):
                    self.dMat1.mult(pvCur[c], pvNxt)
                    self.dm.globalToLocal(pvNxt, self.tmpL)
                    sinkVolP[:, c] = self.tmpL.getArray()

            excess = sinkVol >= marVol
            sinkVol[excess] -= marVol[excess]
            vdep[excess] += marVol[excess]
            marVol[excess] = 0.0

            # In case there is some room to deposit sediments
            noexcess = np.invert(excess)
            marVol[noexcess] -= sinkVol[noexcess]
            vdep[noexcess] += sinkVol[noexcess]
            vdep[self.outletIDs] = 0.0
            sinkVol[noexcess] = 0.0
            sinkVol[self.outletIDs] = 0.0

            # Force-deposit at terminal sinks. Their dMat1 column is zero
            # (rcv==self for all directions, see _matOcean), so the next
            # dMat1.mult would annihilate any remaining sinkVol at these
            # cells once marVol is exhausted. Without this guard ~25% of
            # marine sediment vanished at saturated basin floors on global
            # sphere runs (regression: test_mass_conservation; see
            # AGENTS.md Known bugs entry dated 2026-06). The overflow
            # accumulates above clinoH at the sink; lateral spreading is
            # subsequently handled by _diffuseOcean in seaChange.
            sink_mask = self.is_sink_local
            if sink_mask.any():
                vdep[sink_mask] += sinkVol[sink_mask]
                sinkForced[sink_mask] += sinkVol[sink_mask]
                sinkVol[sink_mask] = 0.0

            # Provenance: split each node's arriving per-class flux into the
            # deposited part and the continuing part in the SAME proportion the
            # total split (a passive label has no settling bias). `cont` is the
            # post-block continuing total; `dep_inc` what deposited this pass.
            # Outlet cells discard their inflow (outflux leaves the domain), so
            # force their deposit increment to zero — matching vdep[outlet]=0.
            # Only split where the arriving flux is non-negligible: below the
            # tiny floor `arr` is denormal round-off (1/arr would overflow to
            # inf, then 0*inf -> nan), and the per-class flux there (Σ_c == arr)
            # is equally negligible, so dropping it is mass-neutral for the label.
            if prov:
                cont = sinkVol
                dep_inc = arr - cont
                dep_inc[self.outletIDs] = 0.0
                safe = arr > 1.0e-12
                frac_dep = np.zeros(self.lpoints)
                frac_cont = np.zeros(self.lpoints)
                np.divide(dep_inc, arr, out=frac_dep, where=safe)
                np.divide(cont, arr, out=frac_cont, where=safe)
                vdepP += sinkVolP * frac_dep[:, None]
                sinkVolP *= frac_cont[:, None]

            # Find where excess and sink
            self.tmpL.setArray(sinkVol)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            if prov:
                for c in range(nC):
                    self.tmpL.setArray(sinkVolP[:, c])
                    self.dm.localToGlobal(self.tmpL, pvCur[c])

            # One Allreduce per iteration that doubles as the loop condition
            # for the next pass and the print value for this pass.
            sumExcess = self.tmp.sum()
            if MPIrank == 0 and self.verbose:
                if step % 100 == 0:
                    print(
                        "  --- Marine excess (sum in km3) %0.05f | iter %d"
                        % (sumExcess * 1.0e-9, step),
                        flush=True
                    )

            step += 1

        # Diagnostic (verbose): marine-routing cascade iteration count. Each
        # pass is a cheap dMat1.mult (not a KSP solve); a high count is the
        # marine routing taking many spillover rounds, not a solver issue.
        if MPIrank == 0 and self.verbose and step > 20:
            print(
                "[sea] marine routing cascade ran %d iterations" % step,
                flush=True,
            )

        # Drain any residual sediment still in self.tmp at loop exit
        # (because the convergence threshold is non-zero, or because the
        # max_iters cap was hit). The in-loop force-deposit at sinks
        # should already drain almost everything; this is the safety net
        # that converts the few remaining cubic metres to a local deposit
        # rather than silently dropping them.
        self.dm.globalToLocal(self.tmp, self.tmpL)
        residual = self.tmpL.getArray().copy()
        if residual.any():
            vdep += residual
            vdep[self.outletIDs] = 0.0

        # Diagnostic: report the two unbounded escapes as THICKNESSES so a
        # needle is obvious (marVol-bounded deposition can never lift a cell
        # above the clinoform surface; these two can).
        reportVolume(self, "marine sink force-dep", sinkForced / self.larea,
                     note="(%d terminal sink cells)"
                     % int(self.is_sink_local.sum()))
        reportVolume(self, "marine residual drain", residual / self.larea,
                     note="(cascade exited after %d iters)" % step)

        # Provenance: drain the per-class residual in lockstep (Σ_c == residual),
        # store the per-node deposited composition for _marineProvFraction, and
        # release the per-class scratch vectors.
        if prov:
            for c in range(nC):
                self.dm.globalToLocal(pvCur[c], self.tmpL)
                vdepP[:, c] += self.tmpL.getArray()
                pvCur[c].destroy()
            vdepP[self.outletIDs, :] = 0.0
            pvNxt.destroy()
            self._marDepProv = vdepP

        if self.memclear:
            del marVol, sinkVol
            gc.collect()

        return vdep

    def _gravFineFraction(self, hl, sedFlux, fineFlux):
        """
        Dual lithology with gravity currents: fine fraction of the combined
        marine deposit. ``_marineFineFraction`` sets the fine fraction of the
        cascade + diffusion deposit only (from the load left to it), and each
        node's fraction is then the volume-weighted blend of that deposit and
        the gravity coarse/fine deposits (``self._gravThick``), so the fine
        volume of both pathways is conserved exactly. Rank-local apart from the
        reductions inside ``_marineFineFraction``.

        :arg hl: pre-deposition local elevation
        :arg sedFlux: marine sediment volume left to the cascade (m^3)
        :arg fineFlux: fine part of ``sedFlux`` (m^3)
        """
        casc, gC, gF = self._gravThick
        self._marineFineFraction(hl, sedFlux, fineFlux, mdep=casc)
        grav = (gC + gF) > 0.0
        if grav.any():
            ffc = np.where(casc > 0.0, self.depoFineFrac, 0.0)
            tot = casc + gC + gF
            self.depoFineFrac[grav] = (ffc[grav] * casc[grav] + gF[grav]) / tot[grav]
        return

    def _gravProvFraction(self):
        """
        Provenance with gravity currents: source composition of the combined
        marine deposit. ``_marineProvFraction`` has set the composition of the
        cascade + diffusion deposit; at nodes that also received a gravity
        deposit, the composition becomes the volume-weighted blend of that
        deposit and the per-class gravity deposit thicknesses
        (``self._gravProvThick``), so each class volume is conserved.
        ``deposeStrat`` then normalises the composition of every depositing
        node. Rank-local.
        """
        casc, gC, gF = self._gravThick
        gp = self._gravProvThick
        grav = (gC + gF) > 0.0
        if gp is None or not grav.any():
            return
        tot = casc + gC + gF
        self.depoProvFrac[grav, :] = (
            self.depoProvFrac[grav, :] * casc[grav, None] + gp[grav, :]
        ) / tot[grav, None]
        return

    def _marineFineFraction(self, hl, sedFlux, fineFlux, mdep=None):
        """
        Dual-lithology (3c): set the per-node fine fraction of the marine
        deposit so fine concentrates in deep / distal water and coarse stays
        proximal (near the coast), conserving the marine fine volume.

        Composition only — the marine deposit geometry (already in ``self.tmp``
        after ``_distOcean``/``_depMarineSystem``/``_diffuseOcean``) is NOT
        modified, so elevations/flow are untouched. This is the subaqueous
        analogue of ``sedplex._pitFineFraction``: the lower-Gmar "fines deposit
        less readily / travel farther" behaviour is expressed as fine biased
        toward the deeper deposits rather than re-routed.

        Mechanism:
          - Domain marine fine fraction ``ff_mar = Σ(fineFlux) / Σ(sedFlux)`` —
            the fine that actually reached the sea after the continental
            cascade (``fineFlux`` is the post-cascade fine input, ``sedFlux``
            the total input, both zero outside the marine sinks). With the
            fine-enriched overspill in the cascade, this is genuinely enriched
            relative to the upstream-eroded composition.
          - Water depth ``d = max(sealevel − hl, 0)`` is the distal/deep proxy.
          - Fine fraction biased ∝ d, renormalised to the deposit-weighted mean
            depth, so ``Σ(mdep·larea·ffrac) == ff_mar·Σ(mdep·larea)``.

        :arg hl: pre-deposition local elevation
        :arg sedFlux: marine input sediment volume per node (m^3)
        :arg fineFlux: marine input fine sediment volume per node (m^3)
        :arg mdep: marine deposit thickness to use (local array); default the
            deposit held in ``self.tmp`` (the gravity-current path passes the
            cascade part only)
        """
        if mdep is None:
            self.dm.globalToLocal(self.tmp, self.tmpL)
            mdep = self.tmpL.getArray().copy()
        owned = self.inIDs == 1

        # Fine fraction of sediment reaching the marine domain (post-cascade).
        fineNum = float(np.sum(fineFlux[owned]))
        den = float(np.sum(sedFlux[owned]))
        fineNum = MPI.COMM_WORLD.allreduce(fineNum, op=MPI.SUM)
        den = MPI.COMM_WORLD.allreduce(den, op=MPI.SUM)
        ff_mar = fineNum / den if den > 0.0 else 0.0

        # Deposit-weighted mean water depth.
        depth = np.maximum(self.sealevel - hl, 0.0)
        mvol = mdep * self.larea
        dw = MPI.COMM_WORLD.allreduce(float(np.sum(mvol[owned])), op=MPI.SUM)
        dwd = MPI.COMM_WORLD.allreduce(float(np.sum((mvol * depth)[owned])), op=MPI.SUM)
        depthbar = dwd / dw if dw > 0.0 else 0.0

        if depthbar > 0.0:
            marine = mdep > 0.0
            self.depoFineFrac[marine] = self._conservingFineFraction(
                ff_mar, depth / max(depthbar, 1.0e-30), mvol, owned, dw
            )[marine]

        return

    def _conservingFineFraction(self, ff_mar, w, mvol, owned, dw):
        """
        Per-node fine fraction ``clip(s * w, 0, 1)`` with the scale ``s`` chosen
        so the deposited fine volume is exactly ``ff_mar * Σ mvol``.

        The plain ``clip(ff_mar * w, 0, 1)`` conserves only while no node needs
        a fraction above 1: deep deposits (large ``w``) were clipped and the
        clipped fine was lost (0.97% of the fine budget on `minimal_dual` once
        the marine eps-fill routed sediment into deep water; fixed 2026-11-12).
        ``Σ mvol * clip(s * w, 0, 1)`` is monotone in ``s``, so ``s`` is found by
        bisection. If even ``clip = 1`` wherever ``w > 0`` cannot absorb the
        fine, the remainder goes uniformly to the ``w == 0`` deposits.

        Collective: one ``allreduce`` per bisection step; every rank runs the
        same number of steps (global sums and a fixed iteration count).

        :arg ff_mar: domain fine fraction of the marine input
        :arg w: per-node depth weight (``depth / depth_bar``)
        :arg mvol: per-node deposited volume (m^3)
        :arg owned: owned-node mask
        :arg dw: ``Σ mvol`` over owned nodes of all ranks
        :return: per-node fine fraction (local array)
        """

        comm = MPI.COMM_WORLD
        target = ff_mar * dw
        if target <= 0.0:
            # No fine reached the sea (e.g. an all-coarse dual run): exactly 0,
            # keeping all-coarse dual bit-identical to single-fraction.
            return np.zeros_like(w)

        def total(s):
            f = np.clip(s * w, 0.0, 1.0)
            return comm.allreduce(float(np.sum((mvol * f)[owned])), op=MPI.SUM)

        pos = w > 0.0
        cap = comm.allreduce(float(np.sum(mvol[owned & pos])), op=MPI.SUM)
        if target >= cap:
            # Every w > 0 deposit is all-fine; the rest fills the w == 0 ones.
            f = np.where(pos, 1.0, 0.0)
            rest = dw - cap
            if rest > 0.0:
                f[~pos] = min(1.0, (target - cap) / rest)
            return f

        lo, hi = 0.0, max(ff_mar, 1.0e-30)
        while total(hi) < target:
            hi *= 2.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if total(mid) < target:
                lo = mid
            else:
                hi = mid
        return np.clip(hi * w, 0.0, 1.0)

    def _marineProvFraction(self, sedFlux):
        """
        In-model provenance (B2b, marine, lockstep): set each marine node's
        deposit source composition to the **per-node mix that was actually
        routed there** by ``_distOcean`` — ``depoProvFrac[node, c] =
        _marDepProv[node, c] / Σ_c _marDepProv[node, c]`` — so a source confined
        to one margin's rivers stays confined to its offshore depocenter instead
        of being smeared across the whole ocean (the previous domain-uniform
        average, see DESIGN_PROVENANCE.md §6). Provenance is a passive label, so
        deposition does not sort by class — the lockstep split is exactly
        proportional, hence ``Σ_c depoProvFrac == 1`` and ``stratP`` still
        partitions ``stratH``.

        Fallback: ``_diffuseOcean`` may give a marine node a thin deposit by
        lateral diffusion from neighbours after the routing pass, so it has
        deposit thickness but no routed composition (``Σ_c _marDepProv == 0``
        there). Those few nodes get the domain-wide basin-delivered mix
        ``Σ(provFrac·sedFlux)/Σ sedFlux`` — the old uniform value — which keeps
        them sensible and conserving without affecting the spatially-resolved
        bulk.

        :arg sedFlux: marine input sediment volume per node (m^3).
        """
        # Prefer the post-diffusion per-class composition (_marDiffProv): the
        # class thicknesses spread by the SAME marine diffusion as the total, so
        # the diffused far field carries the locally-mixed source rather than the
        # domain average. Fall back to the pre-diffusion routed deposit if the
        # diffusion tracer was not run (e.g. diffusion disabled).
        mdp = getattr(self, "_marDiffProv", None)
        if mdp is None:
            mdp = getattr(self, "_marDepProv", None)
        if mdp is None:
            return

        # Marine deposit mask: post-diffusion deposited thickness in self.tmp.
        self.dm.globalToLocal(self.tmp, self.tmpL)
        marine = self.tmpL.getArray() > 0.0

        # Per-node routed composition (spatially resolved).
        denom = mdp.sum(axis=1)
        routed = denom > 0.0
        comp = np.divide(
            mdp, denom[:, None], out=np.zeros_like(mdp), where=denom[:, None] > 0.0
        )
        sel = marine & routed
        self.depoProvFrac[sel, :] = comp[sel, :]

        # Domain-average fallback for diffusion-only marine nodes (no routed
        # composition). Computed collectively (every rank participates, so the
        # Allreduce is safe even when only some ranks own such nodes).
        owned = self.inIDs == 1
        den = MPI.COMM_WORLD.allreduce(float(np.sum(sedFlux[owned])), op=MPI.SUM)
        if den > 0.0:
            ffmar = np.zeros(self.provNb)
            for c in range(self.provNb):
                num = MPI.COMM_WORLD.allreduce(
                    float(np.sum((self.provFrac[:, c] * sedFlux)[owned])), op=MPI.SUM
                )
                ffmar[c] = num / den
            fallback = marine & ~routed
            if fallback.any():
                self.depoProvFrac[fallback, :] = ffmar
        return

    def seaChange(self):
        """
        This function is the main entry point to perform marine river-induced deposition. It calls the private functions:

        - _distanceCoasts
        - _matOcean
        - _diffuseOcean
        """

        t0 = process_time()

        # Set all nodes below sea-level as sinks
        self.sinkIDs = self.lFill <= self.sealevel

        # Define coastal distance for marine points
        hl = self.hLocal.getArray().copy()
        if self.clinSlp > 0.0:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self._distanceCoasts(hl)
            # From the distance to coastline define the upper limit of the shelf to ensure a maximum slope angle
            self.clinoH = self.sealevel - 1.0e-3 - self.coastDist * self.clinSlp  # TODO-REFACTOR: value matches DEPOSIT_FLOOR but distinct role (clinoH 1mm offset below sealevel); do not replace
        else:
            self.clinoH = np.full(self.lpoints, self.sealevel - 1.0e-3, dtype=np.float64)  # TODO-REFACTOR: value matches DEPOSIT_FLOOR but distinct role (clinoH 1mm offset); do not replace
        self.clinoH[hl >= self.sealevel] = hl[hl >= self.sealevel]
        self.clinoH[self.clinoH < hl] = hl[self.clinoH < hl]
        self.maxDepQs = (self.clinoH - hl) * self.larea

        # Get the volumetric marine sediment (m3) to distribute during the time step.
        self.vSedLocal.copy(result=self.QsL)
        sedFlux = self.QsL.getArray().copy() * self.dt
        sedFlux[np.invert(self.sinkIDs)] = 0.0
        sedFlux[sedFlux < 0] = 0.0

        # Gravity currents (opt-in): divert a fraction of the arriving load to
        # the routed gravity pathway; the clinoform cascade gets the rest.
        gravOn = getattr(self, "gravOn", False)
        if gravOn:
            if self.gravFraction is not None:
                fgrav = np.full(self.lpoints, self.gravFraction, dtype=np.float64)
            else:
                fgrav = self._gravFloodFraction(
                    self.QsL.getArray().copy(), self.FAL.getArray().copy()
                )
            fgrav[np.invert(self.sinkIDs)] = 0.0
            gravFine = None
            dual = self.stratNb > 0 and self.stratLith
            if dual and self.gravCoarseFactor != 1.0:
                # Gravity currents carry mostly the suspended (fine) load: the
                # coarse part of the load is diverted with coarse_factor times
                # the fraction (dual lithology only).
                fineAll = self.vSedFLocal.getArray().copy() * self.dt
                fineAll[np.invert(self.sinkIDs)] = 0.0
                fineAll = np.clip(fineAll, 0.0, np.maximum(sedFlux, 0.0))
                gravFine = fgrav * fineAll
                gravSed = gravFine + self.gravCoarseFactor * fgrav * (sedFlux - fineAll)
                sedFlux = sedFlux - gravSed
            else:
                gravSed = fgrav * sedFlux
                sedFlux = sedFlux - gravSed
                if dual:
                    fineAll = self.vSedFLocal.getArray().copy() * self.dt
                    fineAll[np.invert(self.sinkIDs)] = 0.0
                    fineAll[fineAll < 0] = 0.0
                    gravFine = np.minimum(fgrav * fineAll, gravSed)

        # Downstream direction matrix for ocean distribution
        self._matOcean()
        # In-model provenance (B2b-marine): the per-class incoming marine flux
        # is `provFrac · sedFlux` (Σ_c == sedFlux, nonzero only at the river-fed
        # marine nodes), routed in lockstep with the total through _distOcean so
        # each offshore deposit carries the source that drained to it.
        provOn = getattr(self, "provOn", False)
        provFlux = (self.provFrac * sedFlux[:, None]) if provOn else None
        marDep = self._distOcean(sedFlux, provFlux=provFlux)
        reportVolume(self, "marine routed deposit", marDep / self.larea)
        if not self.flatModel and self.Gmar > 0.:
            vdep = self._depMarineSystem(marDep)
            if provOn:
                # The coupled deposition solve re-distributes the routed volume;
                # carry the pass-1 per-node composition onto it (passive label),
                # so the second routing pass keeps tracking provenance.
                mdp = self._marDepProv
                denom = mdp.sum(axis=1)
                comp = np.divide(
                    mdp, denom[:, None],
                    out=np.zeros_like(mdp), where=denom[:, None] > 0.0,
                )
                provFlux = vdep[:, None] * comp
            marDep = self._distOcean(vdep, provFlux=provFlux)
        if gravOn:
            # Gravity network (own matrix on the actual bathymetry by default).
            self._gravNetwork(hl)

            # Slope erosion (G3b, optional): bulk volume picked up on steep
            # ground, removed from the stratigraphy here (erodeStrat gives its
            # deposit-equivalent coarse/fine/class split; the river step's
            # thCoarse/thFine/provEro are restored around the call) and added
            # to the routed supply as an explicit source.
            eroMask = None
            eroBulk = np.zeros(self.lpoints, dtype=np.float64)
            failBulk = np.zeros(self.lpoints, dtype=np.float64)
            eroC = np.zeros(self.lpoints, dtype=np.float64)
            eroF = np.zeros(self.lpoints, dtype=np.float64)
            eroP = None
            if self.gravErosionLength is not None or self.gravFailureSlope is not None:
                eroMask = np.zeros(self.lpoints, dtype=bool)
                if self.gravErosionLength is not None:
                    eroBulk, eroMask = self._gravErosion(gravSed, hl)
                # Slope failure (G3): failed sediment joins the currents at
                # its cell, through the same stratigraphy/routing path.
                if self.gravFailureSlope is not None:
                    failBulk, failMask = self._gravFailure(hl)
                    eroBulk = eroBulk + failBulk
                    eroMask = eroMask | failMask
                if self.stratNb > 0:
                    saved = {k: getattr(self, k).copy()
                             for k in ("thCoarse", "thFine", "provEro")
                             if getattr(self, k, None) is not None}
                    self.tmpL.setArray(-eroBulk / self.larea)
                    self.dm.localToGlobal(self.tmpL, self.tmp)
                    self.erodeStrat()
                    eroC = self.thCoarse * self.dt * self.larea
                    if self.stratLith:
                        eroF = self.thFine * self.dt * self.larea
                    if provOn:
                        eroP = self.provEro * self.dt * self.larea[:, None]
                    for k, v in saved.items():
                        setattr(self, k, v)
                else:
                    eroC = eroBulk.copy()
            eroTot = eroC + eroF
            # share of the removed volume that failed (budget reporting)
            failTot = eroTot * np.divide(failBulk, eroBulk, out=np.zeros_like(eroBulk),
                                         where=eroBulk > 0.0)

            # Column 0: the supply; with provenance, one more column per source
            # class (the arriving composition provFrac, which a passive label
            # shares between the coarse and fine parts, plus the eroded classes
            # split by the eroded coarse/fine share). Same operator, so the
            # class deposits sum to the total.
            comp = self.provFrac if provOn else None

            def _cols(v, eroPart):
                if comp is None:
                    return v[:, None]
                cls = comp * (v - eroPart)[:, None]
                if eroP is not None:
                    cls = cls + eroP * np.divide(
                        eroPart, eroTot, out=np.zeros_like(eroTot), where=eroTot > 0
                    )[:, None]
                return np.column_stack([v, cls])

            # Accommodation: a gravity deposit cannot build a cell above the
            # bypass depth (the fine fraction gets what the coarse one leaves).
            room = None
            if self.gravBypassDepth > 0.0:
                room = np.maximum(
                    (self.sealevel - self.gravBypassDepth - hl) * self.larea - marDep,
                    0.0,
                )
            coarse = (gravSed - gravFine if gravFine is not None else gravSed) + eroC
            dC, eC = self._gravRoute(_cols(coarse, eroC), self.gravRunout,
                                     room=room, nodep=eroMask)
            if gravFine is not None:
                roomF = None if room is None else np.maximum(room - dC[:, 0], 0.0)
                fine = gravFine + eroF
                if self._gWF is not None:
                    # the fine fraction on its own (spreading) network
                    self._gW, self._gWF = self._gWF, self._gW
                    self._gSink, self._gSinkF = self._gSinkF, self._gSink
                dF, eF = self._gravRoute(_cols(fine, eroF), self.gravRunoutFine,
                                         room=roomF, nodep=eroMask)
                if self._gWF is not None:
                    self._gW, self._gWF = self._gWF, self._gW
                    self._gSink, self._gSinkF = self._gSinkF, self._gSink
            else:
                dF, eF = np.zeros_like(dC), np.zeros_like(eC)
            if self.gravDepositSpread is None or self.gravDepositSpread > 0.0:
                npass = self._gravSpreadPasses()
                dC = self._gravSpread(dC, npass, hl, exclude=eroMask)
                if gravFine is not None:
                    dF = self._gravSpread(dF, npass, hl, exclude=eroMask)
            gDepC, gExpC, gDepF, gExpF = dC[:, 0], eC[:, 0], dF[:, 0], eF[:, 0]
            gravProvDep = (dC[:, 1:] + dF[:, 1:]) if comp is not None else None
            self._gravRelease()
        self.dMat1.destroy()

        # Diffuse downstream
        dh = np.divide(marDep, self.larea, out=np.zeros_like(self.larea), where=self.larea != 0)
        # Drop sub-millimeter deposits: numerical cleanup to avoid accumulating
        # round-off thicknesses across many timesteps.
        dh[dh < DEPOSIT_FLOOR] = 0.
        self.tmpL.setArray(dh)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        dh = self.tmpL.getArray().copy()

        # Provenance: hand the per-class deposit thickness (composition of the
        # routed deposit × the floored total dh, so Σ_c == dh) to _diffuseOcean
        # so each class is spread by the SAME marine diffusion as the total and
        # the source label follows the sediment into the diffused far field.
        provThick = None
        if provOn:
            mdp = self._marDepProv
            denom = mdp.sum(axis=1)
            comp = np.divide(
                mdp, denom[:, None], out=np.zeros_like(mdp), where=denom[:, None] > 0.0
            )
            provThick = comp * dh[:, None]
        self._diffuseOcean(dh, provThick=provThick)

        if gravOn:
            # Add the gravity deposit (not diffused this step) to the marine
            # deposit held in self.tmp, keep the cascade part for the fine
            # fraction, and check the budget (supplied = deposited + exported).
            self.dm.globalToLocal(self.tmp, self.tmpL)
            gravCasc = self.tmpL.getArray().copy()
            gC = gDepC / self.larea
            gF = gDepF / self.larea
            self._gravThick = (gravCasc, gC, gF)
            self._gravProvThick = (
                gravProvDep / self.larea[:, None] if gravProvDep is not None else None
            )
            eTh = eroBulk / self.larea
            # The stratigraphy gets the DEPOSITS only (cascade + gravity): the
            # slope erosion was already removed from it by erodeStrat. Elevation
            # and erosion/deposition rates get the net change.
            self._gravDepoOnly = gravCasc + gC + gF
            self.tmpL.setArray(gravCasc + gC + gF - eTh)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            # net gravity change: deposit (+) or slope erosion (-)
            self.gravDep += gC + gF - eTh
            owned = self.inIDs == 1
            sums = MPI.COMM_WORLD.allreduce(
                np.array([gravSed[owned].sum(),
                          (gDepC + gDepF)[owned].sum(),
                          (gExpC + gExpF)[owned].sum(),
                          eroTot[owned].sum(),
                          failTot[owned].sum()]),
                op=MPI.SUM,
            )
            # (supplied, deposited, exported, eroded incl. failures, failed)
            self.gravBudget = sums
            if MPIrank == 0 and self.verbose:
                err = (sums[0] + sums[3] - sums[1] - sums[2]) / max(sums[0] + sums[3], 1.0)
                print(
                    "[gravity] supplied %0.4g km3, eroded %0.4g km3 (of which "
                    "failed %0.4g km3), deposited %0.4g km3, exported %0.4g km3 "
                    "(closure %0.1e)"
                    % (sums[0] * 1e-9, sums[3] * 1e-9, sums[4] * 1e-9,
                       sums[1] * 1e-9, sums[2] * 1e-9, err),
                    flush=True,
                )

        # Update cumulative erosion and deposition as well as elevation
        self.dm.globalToLocal(self.tmp, self.tmpL)
        reportVolume(self, "marine diffused depo", self.tmpL.getArray().copy())
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update soil thickness (both modes — updateSoilThickness applies the
        # subaqueous/ice gates; in regolith mode it skips the deposition increment
        # since marine deposition is stratigraphy, not soil). Marine soil is 0
        # anyway via the subaerial gate. DESIGN_SOIL_REGOLITH.md §5.
        if self.cptSoil:
            self.updateSoilThickness()

        # Update erosion/deposition rates
        self.dm.globalToLocal(self.tmp, self.tmpL)
        add_rate = self.tmpL.getArray() / self.dt
        self.tmpL.setArray(add_rate)
        self.EbLocal.axpy(1.0, self.tmpL)

        # Update stratigraphic layer parameters
        if self.stratNb > 0:
            if gravOn:
                # deposits only: the gravity slope erosion is already out of
                # the stratigraphy (erodeStrat in the gravity step)
                self.tmpL.setArray(self._gravDepoOnly)
                self.dm.localToGlobal(self.tmpL, self.tmp)
            if self.stratLith:
                # Post-cascade fine reaching the sea (same sink mask as
                # sedFlux), carried by the continental fine-enriched overspill.
                fineFlux = self.vSedFLocal.getArray().copy() * self.dt
                fineFlux[np.invert(self.sinkIDs)] = 0.0
                fineFlux[fineFlux < 0] = 0.0
                if gravOn:
                    self._gravFineFraction(hl, sedFlux, fineFlux - gravFine)
                else:
                    self._marineFineFraction(hl, sedFlux, fineFlux)
            if getattr(self, "provOn", False):
                self._marineProvFraction(sedFlux)
                if gravOn:
                    self._gravProvFraction()
            self.deposeStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Distribute River Sediments in the Ocean (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        if self.memclear:
            del marDep, dh, add_rate, hl, sedFlux
            gc.collect()

        return
