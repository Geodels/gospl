import os
import gc
import sys
import scipy
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

from gospl.tools.constants import BOUNDARY_FLOW_SENTINEL, MISSING_DATA_SENTINEL

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import mfdreceivers
    from gospl._fortran import epsfill

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class IceMesh(object):
    """
    This class calculates **ice flow acculation** based on a multiple flow direction paths (MFD).

    Useful links for improvements:
    - `Braun et al. (1999) <https://doi.org/10.3189/172756499781821797>`_
    - `Herman & Braun (2008) <https://doi.org/10.1029/2007JF000807>`_
    - `Tomkin (2009) <https://doi.org/10.1016/j.geomorph.2008.04.021>`_
    - `Egholm (2012) <https://doi.org/10.1016/j.geomorph.2011.12.019>`_
    - `Deal & Prasicek (2020) <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL089263>`_
    - `Hergarten (2021) <http://dx.doi.org/10.5194/esurf-2021-1>`_
    - `Liebl et al. (2023) <https://gmd.copernicus.org/articles/16/1315/2023/>`_
    """

    def __init__(self, *args, **kwargs):
        """
        Initialisation of the `IceMesh` class.

        This method initializes the ice-related fields (ice height, flow accumulation, and flexural response)
        based on the configuration flags `iceOn` and `flexOn`. Memory is allocated for these fields.
        """

        if self.iceOn:
            self.iceHL = self.hLocal.duplicate()
            self.iceFAG = self.hGlobal.duplicate()
            self.iceFAL = self.hLocal.duplicate()
            # Local meltwater field captured at the end of iceAccumulation
            # and re-injected into the river FA source term in flowplex so
            # glacial discharge is not lost downstream.
            self.iceMeltL = self.hLocal.duplicate()
            if self.flexOn:
                self.iceFlex = self.hLocal.duplicate()
            self.iceHL.set(0.0)
            self.iceMeltL.set(0.0)

        return

    def _matrixIceFlow(self, dir_ice=1):
        """
        Compute Flow Direction Matrix for Ice Flow.

        This function calculates the flow direction matrix for ice flow using a multiple flow direction (MFD) algorithm.
        It fills in elevation data (using `epsfill`) and constructs the matrix for flow direction and weighting.

        Parameters:
        -----------
        dir_ice : int, optional
            Number of flow directions to consider. Defaults to 1.
        """

        # The filled + eps is done on the global grid!
        hl = self.hLocal.getArray().copy()
        minh = self.hGlobal.min()[1] + 0.1
        minh = max(minh, self.sealevel)
        if self.flatModel:
            hl[self.idBorders] = BOUNDARY_FLOW_SENTINEL
        fillz = np.zeros(self.mpoints, dtype=np.float64) + MISSING_DATA_SENTINEL
        fillz[self.locIDs] = hl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, fillz, op=MPI.MAX)
        if MPIrank == 0:
            fillz = epsfill(minh, fillz)
        # Broadcast filled elevation data across processors
        fillEPS = MPI.COMM_WORLD.bcast(fillz, root=0)
        fillz = fillEPS[self.locIDs]

        # Calculate receivers and weights for the flow direction matrix
        rcv, _, wght = mfdreceivers(dir_ice, 1.0, fillz, BOUNDARY_FLOW_SENTINEL)

        # Handle borders for flat models
        if self.flatModel:
            rcv[self.idBorders, :] = np.tile(self.idBorders, (dir_ice, 1)).T
            wght[self.idBorders, :] = 0.0

        self.iceRcv = rcv.copy()
        self.iceWght = wght.copy()

        # Create and assemble the flow direction matrix
        self.iceMat = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]

        for k in range(dir_ice):
            # Flow direction matrix for a specific direction
            tmpMat = self._matrix_build()
            data = wght[:, k].copy()
            data[rcv[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                rcv[:, k].astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()

            # Accumulate weights for each direction
            self.iceMat.axpy(-1.0, tmpMat)
            tmpMat.destroy()

        if self.memclear:
            del data, indptr, nodes
            del hl, fillz, fillEPS, rcv, wght
            gc.collect()

        # Transpose the flow matrix for further computations
        self.iceMat.transpose()

        return

    def iceAccumulation(self):
        """
        Main Ice Accumulation Calculation.

        This method calculates ice accumulation based on the elevation, equilibrium-line altitude (ELA), and terminus positions. It integrates the flow direction matrix and computes ice flow across the landscape.

        The method also accounts for flexural responses if enabled.
        """

        ti = process_time()

        # Copy ice height to flexural parameter if flexure modeling is active
        if self.flexOn:
            self.iceHL.copy(result=self.iceFlex)

        # Get dynamic properties such as ice cap altitude and equilibrium-line altitude
        iceH = self.iceH(self.tNow)  # Ice Cap Altitude
        elaH = self.elaH(self.tNow)  # Equilibrium-Line Altitude
        iceT = self.iceT(self.tNow)  # Glacier terminus

        # If maximum elevation is below ELA, no ice accumulation occurs
        max_elev = self.hGlobal.max()[1]
        if max_elev < elaH:
            self.iceHL.set(0.)
            self.iceFAL.set(0.)
            self.iceFAG.set(0.)
            self.iceMeltL.set(0.)
            if self.flexOn:
                self.iceFlex.set(0.)
            return

        # Degenerate config: ice-cap altitude must lie above the ELA, otherwise
        # the (hl - elaH) / (iceH - elaH) ratio below produces NaN/Inf.
        if iceH <= elaH:
            self.iceHL.set(0.)
            self.iceFAL.set(0.)
            self.iceFAG.set(0.)
            self.iceMeltL.set(0.)
            if self.flexOn:
                self.iceFlex.set(0.)
            return

        # Compute the flow direction matrix for ice
        self._matrixIceFlow(self.iceDir)

        # Ice accumulation calculation
        hl = self.hLocal.getArray().copy()
        rainA = self.bL.getArray().copy()
        tmp = (hl - elaH) / (iceH - elaH)
        tmp[tmp > 1.] = 1.0
        tmp[tmp < 0.] *= self.meltfac

        # Calculate accumulation rates
        iceA = np.multiply(rainA, tmp)
        self.tmpL.setArray(iceA)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self._solve_KSP(True, self.iceMat, self.tmp, self.iceFAG)
        self.dm.globalToLocal(self.iceFAG, self.iceFAL)
        ice_clamped = self.iceFAL.getArray()
        ice_clamped[ice_clamped < 0.] = 0.0
        ice_clamped[hl < iceT] = 0.0
        self.iceFAL.setArray(ice_clamped)
        self.dm.localToGlobal(self.iceFAL, self.iceFAG)

        # Capture glacial meltwater for the river FA. Below ELA the local
        # ablation rate is `|hl - elaH|/(iceH - elaH)` of the precipitation;
        # we gate by ice presence (clamped, pre-smoothing iceFAL > 0) so we
        # only release water where ice actually reached. `meltfac` is
        # deliberately NOT applied here because it is a numerical
        # sink-amplifier inside the implicit ice solver, not a physical
        # melt multiplier. This is consumed by flowplex.flowAccumulation
        # and added back to `rainA` before the river FA is solved.
        melt_local = np.zeros(self.lpoints, dtype=np.float64)
        below_ela = hl < elaH
        if below_ela.any():
            ablation = np.zeros(self.lpoints, dtype=np.float64)
            ablation[below_ela] = np.clip(
                -(hl[below_ela] - elaH) / (iceH - elaH), 0.0, 1.0
            )
            has_ice = ice_clamped > 1.0e-8  # TODO-REFACTOR: value matches DISCHARGE_FLOOR but distinct role (ice-presence threshold for melt capture); do not replace
            melt_local = ablation * rainA * has_ice
        self.iceMeltL.setArray(melt_local)

        # Smooth and diffuse glacier flow accumulation.
        # NB: the linear-diffusion smoothing below is NOT mass-conservative
        # (it spreads the field for visualisation and erosion-driver
        # robustness). The output `iceFA` therefore differs slightly from
        # the strict accumulated source. The meltwater field above was
        # captured BEFORE smoothing, so the river re-injection is based on
        # the true (un-smeared) ice presence.
        self.iceFAG.copy(result=self.tmp1)
        smthIce = self._hillSlope(smooth=1)
        smthIce[smthIce < 0.] = 0.
        self.iceFAL.setArray(smthIce)
        self.dm.localToGlobal(self.iceFAL, self.iceFAG)

        # Ice thickness from Bahr-style width-area scaling. Computed
        # whenever ice is on (used by both flexure and the `iceH` output);
        # the flex copy stays gated behind `flexOn`.
        tmp = self.icewe * self.icewf * smthIce**0.3
        tmp[tmp < 1.e-1] = 0.  # TODO-REFACTOR: value matches BEDROCK_EXPOSED but distinct role (minimum ice thickness for visualization); do not replace
        self.tmpL.setArray(tmp)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        tmp = self._hillSlope(smooth=1)
        self.iceHL.setArray(tmp)

        if self.flexOn:
            # If simulation starts then set the ice flex variable to the initial glacier thickness
            if self.tNow == self.tStart:
                self.iceHL.copy(result=self.iceFlex)

        if MPIrank == 0 and self.verbose:
            print(
                "Glaciers Accumulation (%0.02f seconds)" % (process_time() - ti),
                flush=True,
            )

        return
