import os

import gc
import sys
import vtk
vtk.vtkObject.GlobalWarningDisplayOff()

import warnings
import petsc4py
from gospl.tools.petscgc import safe_garbage_cleanup
import numpy as np
from scipy import spatial

from mpi4py import MPI
from time import process_time
from vtk.util import numpy_support  # type: ignore

from gospl.tools.constants import (
    BOUNDARY_FLOW_SENTINEL,
    DEPOSIT_FLOOR,
    MISSING_DATA_SENTINEL,
)

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import mfdrcvrs
    from gospl._fortran import epsfill

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

        return

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

        # Define multiple flow directions for filled + eps elevations
        hl = self.hLocal.getArray().copy()
        minh = self.hGlobal.min()[1] + 0.1
        if not self.flatModel:
            # Only consider filleps in the first kms offshore
            minh = max(minh, self.oFill)
            hsmth = self._hillSlope(smooth=2)
            hsmth[self.coastDist > self.offshore] = BOUNDARY_FLOW_SENTINEL
        else:
            minh = max(minh, self.oFill)
            hsmth = hl.copy()
            hsmth[self.outletIDs] = BOUNDARY_FLOW_SENTINEL

        # The filled + eps is done on the global grid, serially on rank 0, so
        # assemble it there from the owned nodes only (hsmth is ghost-synced
        # via _hillSlope's globalToLocal, so owned values are exact).
        fillz = self._gatherGlobalOnRoot(hsmth)
        if MPIrank == 0:
            fillz = epsfill(minh, fillz)
        # Send elevation + eps to other processors
        fillEPS = MPI.COMM_WORLD.bcast(fillz, root=0)
        fillz = fillEPS[self.locIDs]
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
            del hl, fillz, fillEPS, rcv, wght
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
                        % (sumExcess * 10.e-9, step),
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

    def _marineFineFraction(self, hl, sedFlux, fineFlux):
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
        """
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
            ffrac = np.clip(ff_mar * depth / max(depthbar, 1.0e-30), 0.0, 1.0)
            marine = mdep > 0.0
            self.depoFineFrac[marine] = ffrac[marine]

        return

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

        # Downstream direction matrix for ocean distribution
        self._matOcean()
        # In-model provenance (B2b-marine): the per-class incoming marine flux
        # is `provFrac · sedFlux` (Σ_c == sedFlux, nonzero only at the river-fed
        # marine nodes), routed in lockstep with the total through _distOcean so
        # each offshore deposit carries the source that drained to it.
        provOn = getattr(self, "provOn", False)
        provFlux = (self.provFrac * sedFlux[:, None]) if provOn else None
        marDep = self._distOcean(sedFlux, provFlux=provFlux)
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

        # Update cumulative erosion and deposition as well as elevation
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update soil thickness
        if self.cptSoil:
            self.updateSoilThickness()

        # Update erosion/deposition rates
        self.dm.globalToLocal(self.tmp, self.tmpL)
        add_rate = self.tmpL.getArray() / self.dt
        self.tmpL.setArray(add_rate)
        self.EbLocal.axpy(1.0, self.tmpL)

        # Update stratigraphic layer parameters
        if self.stratNb > 0:
            if self.stratLith:
                # Post-cascade fine reaching the sea (same sink mask as
                # sedFlux), carried by the continental fine-enriched overspill.
                fineFlux = self.vSedFLocal.getArray().copy() * self.dt
                fineFlux[np.invert(self.sinkIDs)] = 0.0
                fineFlux[fineFlux < 0] = 0.0
                self._marineFineFraction(hl, sedFlux, fineFlux)
            if getattr(self, "provOn", False):
                self._marineProvFraction(sedFlux)
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
