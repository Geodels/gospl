import os
import sys
import petsc4py
import numpy as np

from mpi4py import MPI
from scipy import spatial
from time import process_time

from gospl.tools.constants import MISSING_DATA_SENTINEL, MISSING_LARGE_SENTINEL

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import adveciioe
    from gospl._fortran import adveciioe2
    from gospl._fortran import advecupwind
    from gospl._fortran import fitedges
    from gospl._fortran import getrange
    from gospl._fortran import getfacevelocity

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()


class Tectonics(object):

    """
    This class defines how 2D and spherical mesh surfaces are changing given a set of horizontal and vertical rates.

    .. note::

        Three advection approaches are proposed, a standard upwind scheme, an Inflow-Implicit/Outflow-Explicit scheme and a semi-lagrangian approach based on nearest neighbour search based on a kdTree-search where interpolation is weighted by distance.
    """

    def __init__(self):
        """
        The initialisation of the `Tectonics` class.
        """

        self.hdisp = None
        self.tecNb = -1
        self.paleoZ = None
        self.plateStep = False

        self.fiso = self.hGlobal.duplicate()

        # Cached horizontal-advection operator (FV upwind / IIOE). The matrix
        # coefficients depend only on the face velocities (refreshed when the
        # tectonic interval advances) and on `self.dt`, so within an interval at
        # fixed dt the operator is constant: build it once and reuse it, along
        # with its KSP (so the block-Jacobi PCSetUp is factorised once too).
        self._advMatLeft = None
        self._advMatRight = None
        self._advKSP = None
        self._advNbOut = None
        self._advDt = None
        self._advRebuild = True

        return

    def getTectonics(self):
        """
        Finds the current tectonic regimes (horizontal and vertical) for the considered time interval.

        .. note::
            The approach here does not allow for mesh refinement in zones of convergence and thus can be limiting.

            Yet using a fixed mesh has one main advantage: the mesh and Finite Volume discretisation do not have to be rebuilt each time the mesh is advected.
        """

        if self.tecdata is None:
            return

        nb = self.tecNb
        if nb < len(self.tecdata) - 1:
            if self.tecdata.at[nb + 1, "start"] < self.tNow + self.dt:
                nb += 1

        if nb > self.tecNb or nb == -1:
            if nb == -1:
                nb = 0

            self.tecNb = nb
            if nb < len(self.tecdata.index) - 1:
                timer = self.tecdata.at[nb, "end"] - self.tecdata.at[nb, "start"]
            else:
                timer = self.tEnd - self.tecdata.at[nb, "start"]

            # Horizontal displacements
            # was iloc[nb, -1]; named access is safer — see AGENTS.md
            # Forcing DataFrame layout contract
            if self.tecdata.at[nb, "hMap"] != "empty":
                fname = self.tecdata.at[nb, "hMap"][0] + ".npz"
                mdata = np.load(fname)
                key = self.tecdata.at[nb, "hMap"][1]
                self.hdisp = mdata[key][self.locIDs, :]
                # In case of advection based on interpolation from plate position
                if self.advscheme == 0:
                    self._readAdvectionData(mdata[key], timer)
                    self.plateStep = True
                    self.plateTimer = self.tecdata.at[nb, "end"]
                else:
                    # Get the velocity from the input file.
                    nodeVel = np.zeros((self.lpoints, 3))
                    if self.cyclicBC:
                        # Cyclic 2D model: the mesh is a cylinder. The user keeps
                        # supplying the displacement in the flat (vx, vy) frame
                        # (as for a planar model); map it onto the cylinder
                        # tangent so the periodic-axis component advects AROUND
                        # the seam instead of pushing radially through it.
                        nodeVel = self._cylinderVelocity(self.hdisp)
                    elif self.flatModel:
                        # Planar 2D model: the displacement lives in the (x, y)
                        # plane (the z component is zero).
                        nodeVel[:, :2] = self.hdisp[:, :2]
                    else:
                        nodeVel = self.hdisp.copy()
                    # Store velocity on voronoi edges and dot product based on
                    # vecocity vector and face normals. The face velocities just
                    # changed, so the cached advection operator must be rebuilt.
                    getfacevelocity(self.lpoints, nodeVel)
                    self._advRebuild = True
            else:
                self.hdisp = None

            # Vertical displacements
            if self.tecdata.at[nb, "tMap"] != "empty":
                fname = self.tecdata.at[nb, "tMap"][0] + ".npz"
                mdata = np.load(fname)
                key = self.tecdata.at[nb, "tMap"][1]
                self.upsub = mdata[key][self.locIDs]
            else:
                self.upsub = None

            # Paleo-elevation fitting
            self.paleoZ = None
            if self.tecdata.at[nb, "zMap"] != "empty":
                fname = self.tecdata.at[nb, "zMap"][0] + ".npz"
                mdata = np.load(fname)
                key = self.tecdata.at[nb, "zMap"][1]
                if len(self.tecdata.at[nb, "zMap"]) == 3:
                    self.paleoZ = mdata[key][self.locIDs]

            del mdata

        # Perform advection based on the flow-Implicit/Outflow-Explicit
        # or a `first-order upwind implicitly scheme.
        if self.hdisp is not None and self.advscheme > 0:
            self._varAdvector()

        return

    def _cylinderVelocity(self, hdisp):
        """
        Map a flat ``(vx, vy)`` horizontal displacement onto the tangent of the
        cylinder used for a cyclic (periodic) 2D model.

        A cyclic mesh is a cylinder: the periodic flat axis is wrapped into a
        circle lying in the plane that contains the synthetic ``z`` dimension,
        while the other flat axis stays the (straight) cylinder axis. A flat
        velocity in the periodic direction must therefore become a velocity
        *around* the cylinder — tangent to that circle (the local ``θ̂``) — not a
        Cartesian translation, which would push radially through the surface and
        be discarded by the face-normal dot product. The non-periodic component
        is unchanged (it is already tangent, along the cylinder axis).
        """

        coords = self.lcoords
        vel = np.zeros((self.lpoints, 3))
        if self.east == 2 or self.west == 2:
            # Periodic in x -> circle in (x, z), cylinder axis = y.
            rad = np.sqrt(coords[:, 0] ** 2 + coords[:, 2] ** 2)
            rad[rad == 0.0] = 1.0
            cos, sin = coords[:, 0] / rad, coords[:, 2] / rad
            vel[:, 0] = -hdisp[:, 0] * sin   # around-seam component (θ̂)
            vel[:, 2] = hdisp[:, 0] * cos
            vel[:, 1] = hdisp[:, 1]          # axial component (unchanged)
        else:
            # Periodic in y -> circle in (y, z), cylinder axis = x.
            rad = np.sqrt(coords[:, 1] ** 2 + coords[:, 2] ** 2)
            rad[rad == 0.0] = 1.0
            cos, sin = coords[:, 1] / rad, coords[:, 2] / rad
            vel[:, 1] = -hdisp[:, 1] * sin   # around-seam component (θ̂)
            vel[:, 2] = hdisp[:, 1] * cos
            vel[:, 0] = hdisp[:, 0]          # axial component (unchanged)
        return vel

    def _buildAdvecMat(self, iioe, lCoeffs, rCoeffs=None):
        """
        Build the advection operator(s) from the per-cell FV coefficients
        ``(lpoints, 1+maxnb)`` (col 0 = diagonal, cols ``1..maxnb`` = neighbour
        entries for ``self.FVmesh_ngbID``) in a single-pass CSR assembly
        (``_assembleDiffMatCSR``) rather than the old ``maxnb``-matrix ``axpy``
        loop.

        :arg iioe: True for the IIOE scheme (build both left/right operators);
            False for the implicit upwind scheme (left operator only).
        :arg lCoeffs: left (implicit) coefficients.
        :arg rCoeffs: right (explicit) coefficients (IIOE only).

        :return: ``(advMat_left, advMat_right)`` if ``iioe`` else ``advMat_left``.
        """

        if self.flatModel:
            # Pin the (non-periodic) domain edges with a Dirichlet row. A cyclic
            # seam stays free so material advects across it (advectBorders
            # excludes the seam nodes).
            lCoeffs[self.advectBorders, 1:] = 0.0
            lCoeffs[self.advectBorders, 0] = 1.0
            if iioe:
                rCoeffs[self.advectBorders, 1:] = 0.0
                rCoeffs[self.advectBorders, 0] = 1.0

        advMat_left = self._assembleDiffMatCSR(lCoeffs)
        if iioe:
            advMat_right = self._assembleDiffMatCSR(rCoeffs)
            return advMat_left, advMat_right
        return advMat_left

    def _advSolve(self, rhs, sol):
        """
        Solve one advection system with the cached operator + KSP. ``sol`` must
        carry the warm-start guess on entry (the field's pre-advection value),
        which is close to the solution since advection moves it only slightly.
        Falls back to the fgmres/asm solver if the cached KSP diverges.
        """

        ksp = self._advKSP
        ksp.solve(rhs, sol)
        if ksp.getConvergedReason() < 0:
            sol = self._solve_KSP2(self._advMatLeft, rhs, sol)

        return sol

    def _advectorIIOE2(self, gvec, vL, newv, vmin, vmax, nbOut):
        """
        Perform the advection for the Inflow-Implicit/Outflow-Explicit Scheme 2:
        an anti-diffusion correction applied where the Scheme-1 result
        (``newv``) over/undershoots the local neighbourhood range
        ``[vmin, vmax]``. The correction operator depends on the per-field range,
        so it is built fresh each call (not cached) and writes the corrected
        result back into ``self.tmp``.

        :arg gvec: the field's pre-advection global Vec (the explicit operator
            multiplies THIS field — previously hard-wired to ``self.hGlobal``,
            which corrupted the correction for cumED/flexure/soil).
        :arg vL: the field's pre-advection local array.
        :arg newv: the Scheme-1 advected local array.
        :arg vmin, vmax: local neighbourhood min/max of the pre-advection field.
        :arg nbOut: outflow-neighbour count from ``adveciioe``.
        """

        diffmax = newv - vmax
        diffmax[diffmax < 0] = 0.
        diffmin = newv - vmin
        diffmin[diffmin > 0] = 0.
        diff = np.abs(diffmax) + np.abs(diffmin)
        # Sum the overshoot into `tmp1` (NOT `tmp`): `tmp` still holds the
        # Scheme-1 advected field, which is the result to keep when there is no
        # overshoot (excess == 0). Clobbering it here would leave the field as
        # the all-zero `diff` and zero the elevation on every no-overshoot step.
        self.tmpL.setArray(diff)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        excess = self.tmp1.sum()

        if excess > 0.:
            lCoeffs, rCoeffs = adveciioe2(self.lpoints, self.dt, nbOut, vL, vmin, vmax)
            advMat_left2, advMat_right2 = self._buildAdvecMat(True, lCoeffs, rCoeffs)
            advMat_right2.mult(gvec, self.tmp1)
            self._solve_KSP(True, advMat_left2, self.tmp1, self.tmp)
            advMat_left2.destroy()
            advMat_right2.destroy()

        return

    def _varAdvector(self):
        """
        Perform the advection of elevation, erosion (and flexure and soil production if defined) by solving the advection equation based on the finite volume space discretization and a semi-implicit temporal discretization. 
        
        The approach relies on an Inflow-Implicit/Outflow-Explicit scheme following the work from `Mikula & Ohlberger, 2014 <https://www.math.sk/mikula/mo-FVCA6.pdf>`_ or a `first-order upwind implicitly scheme <https://www.sciencedirect.com/science/article/pii/S0168927414001032>`_ depending on the user configuration.

        .. note::

           Here, the outflow from a cell is treated explicitly while inflow is treated implicitly.

            Since the matrix of the system is determined by the inflow fluxes it is a M-matrix yielding favourable solvability and stability properties.

        .. important::

            The method allows for large time steps without losing too much stability and precision.

            It is formally second order accurate in space and time for 1D advection problems with variable velocity. In addition, numerical experiments indicate its second order accuracy for smooth solutions in general.

        .. note::

            Velocity at the face is taken to be the linear interpolation for each vertex (in a vertex-centered discretisation the dual of the delaunay triangulation (i.e. the voronoi mesh has its edges intersecting the middle of the nodes edges)).

            Similarly, we consider that the advected variables at the face are defined by linear interpolation from each connected vertex.
        """

        t0 = process_time()
        iioe = self.advscheme != 1

        # Build (or reuse) the cached advection operator + KSP. The coefficients
        # depend only on the face velocities (refreshed on tectonic-interval
        # change, which sets `_advRebuild`) and on `self.dt`, so within an
        # interval at fixed dt the operator is constant: assemble it once and
        # reuse it (and its block-Jacobi PCSetUp) across every step.
        if self._advRebuild or self._advDt != self.dt:
            if self._advMatLeft is not None:
                self._advMatLeft.destroy()
            if self._advMatRight is not None:
                self._advMatRight.destroy()
                self._advMatRight = None
            if iioe:
                self._advNbOut, lCoeffs, rCoeffs = adveciioe(self.lpoints, self.dt)
                self._advMatLeft, self._advMatRight = self._buildAdvecMat(
                    True, lCoeffs, rCoeffs
                )
            else:
                lCoeffs = advecupwind(self.lpoints, self.dt)
                self._advMatLeft = self._buildAdvecMat(False, lCoeffs)
            if self._advKSP is None:
                self._advKSP = self._makeDiffusionKSP("advect_")
            self._advKSP.setOperators(self._advMatLeft, self._advMatLeft)
            self._advDt = self.dt
            self._advRebuild = False

        nbOut = self._advNbOut if iioe else None

        # Each field is advected the same way; only the boundary-edge treatment
        # differs (elevation/cumED are re-extrapolated, flexure/soil are not).
        self._advectField(self.hGlobal, iioe, nbOut)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        self._resetEdges(self.hLocal, self.hGlobal, MISSING_DATA_SENTINEL)

        # cumED uses the large-magnitude sentinel (consistent with
        # `_advectPlates`): cumulative erosion/deposition can exceed the
        # MISSING_DATA_SENTINEL magnitude, so a smaller marker could collide.
        self._advectField(self.cumED, iioe, nbOut)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self._resetEdges(self.cumEDLocal, self.cumED, MISSING_LARGE_SENTINEL)

        # Flexural isostasy (carried in the numpy `localFlex`; no edge reset).
        if self.flexOn:
            self.tmpL.setArray(self.localFlex)
            self.dm.localToGlobal(self.tmpL, self.fiso)
            self._advectField(self.fiso, iioe, nbOut)
            self.dm.globalToLocal(self.fiso, self.tmpL)
            self.localFlex = self.tmpL.getArray().copy()

        # Soil thickness (no edge reset; clamped to a physical range after).
        if self.cptSoil:
            self._advectField(self.Gsoil, iioe, nbOut)
            self.dm.globalToLocal(self.Gsoil, self.Lsoil)
            lsoil = self.Lsoil.getArray().copy()
            np.clip(lsoil, 0.0, self.soil_transition, out=lsoil)
            self.Lsoil.setArray(lsoil)
            self.dm.localToGlobal(self.Lsoil, self.Gsoil)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Advection Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _advectField(self, gvec, iioe, nbOut):
        """
        Advect one global field Vec **in place** with the cached advection
        operator (and, for IIOE Scheme 2, the anti-diffusion correction).

        The field's current value is used as the warm-start guess (advection
        moves it only slightly, so this converges in very few iterations).

        Scratch used: ``self.tmp`` (solution), ``self.tmp1`` (RHS), ``self.tmpL``
        (local views for the Scheme-2 range/correction).

        :arg gvec: global field Vec (overwritten with the advected field).
        :arg iioe: True for the IIOE schemes, False for implicit upwind.
        :arg nbOut: outflow-neighbour count (IIOE Scheme 2 only; else None).
        """

        if iioe and self.advscheme == 3:
            # Neighbourhood min/max of the pre-advection field for Scheme 2.
            self.dm.globalToLocal(gvec, self.tmpL)
            fL = self.tmpL.getArray().copy()
            fmin, fmax = getrange(self.lpoints, fL)

        gvec.copy(result=self.tmp)            # warm-start guess
        if iioe:
            self._advMatRight.mult(gvec, self.tmp1)
            self._advSolve(self.tmp1, self.tmp)
            if self.advscheme == 3:
                self.dm.globalToLocal(self.tmp, self.tmpL)
                self._advectorIIOE2(gvec, fL, self.tmpL.getArray(), fmin, fmax, nbOut)
        else:
            self._advSolve(gvec, self.tmp)

        self.tmp.copy(result=gvec)

        return

    def _resetEdges(self, lvec, gvec, sentinel):
        """
        Flat-model only: mark the (non-periodic) domain-edge nodes with
        ``sentinel`` and re-extrapolate them from their interior neighbours with
        ``fitedges``, then sync the global Vec. No-op on a global/sphere mesh.

        :arg lvec: local Vec of the field (read + rewritten).
        :arg gvec: matching global Vec (rewritten from ``lvec``).
        :arg sentinel: edge marker (``fitedges`` treats any value ``<= -1e7`` as
            missing, so both data/large sentinels are recognised).
        """

        if not self.flatModel:
            return

        arr = lvec.getArray().copy()
        arr[self.advectBorders] = sentinel
        arr = fitedges(arr)
        lvec.setArray(arr)
        self.dm.localToGlobal(lvec, gvec)

        return

    def _readAdvectionData(self, hdisp, timer):
        """
        From a tectonic input file, this function reads the horizontal displacements information, containing the displacement rates along each axis.

        .. note::
            A cKDTree is built with the advected coordinates and used to interpolate the mesh variables on the initial local mesh position.

            The interpolation is based on a weighting distance function accounting for the 3 closest advected vertices.

        :arg hdisp: displacement rates in 3D
        :arg timer: displacement time step in years
        """

        t0 = process_time()
        XYZ = self.mCoords + hdisp * timer

        # Build a tree with the advected nodes
        tree = spatial.cKDTree(XYZ, leafsize=10)

        # Query the distances from unstructured nodes
        distances, self.tec_IDs = tree.query(self.mCoords, k=3)

        # Inverse weighting distance...
        self.tec_weights = np.divide(
            1.0, distances, out=np.zeros_like(distances), where=distances != 0
        )
        self.tec_onIDs = np.where(distances[:, 0] == 0)[0]
        self.tec_sumw = np.sum(self.tec_weights, axis=1)

        del tree, distances

        if MPIrank == 0 and self.verbose:
            print(
                "Read advection data (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return

    def _advectPlates(self):
        """
        Advects surface information based on plate evolution and performs interpolation.

        .. important::

            The interpolated values are the elevation,  cumulative erosion deposition, flexural response and soil production. The interpolation is done on the spherical mesh on a single processor. In case the stratigraphic information is also recorded then this is also interpolated.
        """

        # Send local elevation globally
        t0 = process_time()
        hl = self.hLocal.getArray().copy()
        gZ = np.zeros(self.mpoints, dtype=np.float64) + MISSING_DATA_SENTINEL
        gZ[self.locIDs] = hl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)

        # Send local erosion deposition globally
        edl = self.cumEDLocal.getArray().copy()
        gED = np.zeros(self.mpoints, dtype=np.float64) + MISSING_LARGE_SENTINEL
        gED[self.locIDs] = edl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gED, op=MPI.MAX)

        # Send local flexural isostasy globally
        if self.flexOn:
            gFI = np.zeros(self.mpoints, dtype=np.float64) + MISSING_LARGE_SENTINEL
            gFI[self.locIDs] = self.localFlex
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gFI, op=MPI.MAX)

        # Send local soil thickness globally
        if self.cptSoil:
            lsoil = self.Lsoil.getArray().copy()
            gSL = np.zeros(self.mpoints, dtype=np.float64) + MISSING_LARGE_SENTINEL
            gSL[self.locIDs] = lsoil
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gSL, op=MPI.MAX)

        if MPIrank == 0 and self.verbose:
            print(
                "Transfer local elevation, erosion deposition, soil and flexural information globally (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        # Perform interpolation
        t0 = process_time()

        tmp = np.sum(self.tec_weights * gZ[self.tec_IDs], axis=1)
        nelev = np.divide(tmp, self.tec_sumw, out=np.zeros_like(self.tec_sumw), where=self.tec_sumw != 0)
        tmp = np.sum(self.tec_weights * gED[self.tec_IDs], axis=1)
        nerodep = np.divide(tmp, self.tec_sumw, out=np.zeros_like(self.tec_sumw), where=self.tec_sumw != 0)
        if self.flexOn:
            tmp = np.sum(self.tec_weights * gFI[self.tec_IDs], axis=1)
            nfi = np.divide(tmp, self.tec_sumw, out=np.zeros_like(self.tec_sumw), where=self.tec_sumw != 0)
        if self.cptSoil:
            tmp = np.sum(self.tec_weights * gSL[self.tec_IDs], axis=1)
            nsoil = np.divide(tmp, self.tec_sumw, out=np.zeros_like(self.tec_sumw), where=self.tec_sumw != 0)

        if len(self.tec_onIDs) > 0:
            nelev[self.tec_onIDs] = gZ[self.tec_IDs[self.tec_onIDs, 0]]
            nerodep[self.tec_onIDs] = gED[self.tec_IDs[self.tec_onIDs, 0]]
            if self.flexOn:
                nfi[self.tec_onIDs] = gFI[self.tec_IDs[self.tec_onIDs, 0]]
            if self.cptSoil:
                nsoil[self.tec_onIDs] = gSL[self.tec_IDs[self.tec_onIDs, 0]]

        if MPIrank == 0 and self.verbose:
            print(
                "Define local elevation and erosion/deposition after advection (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        self.hGlobal.setArray(nelev[self.glbIDs])
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        self.cumED.setArray(nerodep[self.glbIDs])
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)

        if self.flexOn:
            self.localFlex = nfi[self.locIDs].copy()

        if self.cptSoil:
            nsoil[nsoil < 0.] = 0.
            nsoil[nsoil > self.soil_transition] = self.soil_transition
            self.Gsoil.setArray(nsoil[self.glbIDs])
            self.dm.globalToLocal(self.Gsoil, self.Lsoil)

        # Update stratigraphic record
        if self.stratNb > 0 and self.stratStep > 0:
            self._advectStrati()

        return

    def _findPts2Reduce(self):
        """
        Finds the nodes part of another partition that are required for interpolation.

        :return: Indices of needed nodes present in other partitions
        """

        # Global neighbours for each local partition
        lgNghbs = self.idNbghs[self.locIDs, :].flatten()

        # Global IDs required locally but part of another partition
        gids = lgNghbs[~np.isin(lgNghbs, self.locIDs)]

        # Get all points that have to be transferred
        vals = np.zeros(self.mpoints, dtype=int)
        vals[gids] = 1
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, vals, op=MPI.MAX)

        return np.where(vals > 0)[0]

    def _advectStrati(self):
        """
        From advected stratigraphic information, retrieve points that are part of another partition and perform the interpolation.
        """

        # Get lobal point IDs that will need to be transferred globally
        t0 = process_time()
        redIDs = self._findPts2Reduce()
        if MPIrank == 0 and self.verbose:
            print(
                "Send stratigraphy information to local partition (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        # Transfer the variables accordingly and perform operation
        t0 = process_time()
        nghs = self.tec_IDs[self.locIDs]
        wgts = self.tec_weights[self.locIDs, :, None]
        sumw = self.tec_sumw[self.locIDs, :, None]
        # onID = self.tec_onIDs[self.locIDs]
        # tec_onIDs holds global-mesh indices where the advected node
        # coincides exactly. Remap to local-mesh indices and keep only
        # entries this rank owns.
        g2l = np.full(self.mpoints, -1, dtype=int)
        g2l[self.locIDs] = np.arange(len(self.locIDs))
        onID = g2l[self.tec_onIDs]
        onID = onID[onID >= 0]

        self.stratH[:, : self.stratStep] = self._updateStratInfo(
            redIDs, self.stratH[:, : self.stratStep], sumw, onID, wgts, nghs
        )
        self.stratZ[:, : self.stratStep] = self._updateStratInfo(
            redIDs, self.stratZ[:, : self.stratStep], sumw, onID, wgts, nghs
        )
        self.phiS[:, : self.stratStep] = self._updateStratInfo(
            redIDs, self.phiS[:, : self.stratStep], sumw, onID, wgts, nghs
        )
        if MPIrank == 0 and self.verbose:
            print(
                "Define local stratigraphy variables after advection (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )
        del nghs, wgts, sumw, onID

        return

    def _updateStratInfo(self, reduceIDs, variable, sumw, onIDs, weights, nghs):
        """
        Update stratigraphic variables based on interpolated values.

        :arg reduceIDs: indices of nodes information that need to be available globally
        :arg variable: stratigraphic variable to interpolate
        :arg sumw: sum of weights for interpolation
        :arg onIDs: nodes that remain fixed
        :arg weights: weights for interpolation
        :arg nghs: neighbours for interpolation

        :return: nvar interpolated stratigraphic variable
        """

        # Reduce variable so that all values that need to be read from another partition can be
        vals = np.zeros((self.mpoints, self.stratStep), dtype=np.float64) + MISSING_DATA_SENTINEL
        vals[self.locIDs, :] = variable
        redVals = vals[reduceIDs, :]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, redVals, op=MPI.MAX)
        vals[reduceIDs, :] = redVals

        # Inverse weighting distance...
        tmp = np.sum(weights * vals[nghs, :], axis=1)
        nvar = np.divide(tmp, sumw, where=sumw != 0, out=np.zeros_like(tmp))
        if len(onIDs) > 0:
            nvar[onIDs, :] = vals[nghs[onIDs, 0], :]

        return nvar

    def updatePaleoZ(self):
        """
        Update surface information based on paleo-reconstruction.

        :arg paleoZ: paleo-elevation data
        :arg minZ: elevation threshold used for fitting
        """

        if self.plateStep:
            if self.plateTimer == self.tNow:
                self._advectPlates()

        if self.paleoZ is None:
            return

        t0 = process_time()

        # Send local elevation globally
        hl = self.hLocal.getArray().copy()
        ids = np.where(hl <= self.sealevel)[0]
        hl[ids] = self.paleoZ[ids]

        # Fit simulated elevations to paleoelevation ones
        self.hLocal.setArray(hl)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        if MPIrank == 0 and self.verbose:
            print(
                "Update model based on paleo-elevations (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return
