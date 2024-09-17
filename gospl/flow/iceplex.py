import os
import gc
import sys
import scipy
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import mfdreceivers

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class IceMesh(object):
    """
    This class calculates **ice flow acculation** based on a multiple flow direction paths (MFD) performed on a smooth landscape elevation mesh.

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
        """

        if self.iceOn:
            self.iceHL = self.hLocal.duplicate()
            self.iceFAG = self.hGlobal.duplicate()
            self.iceFAL = self.hLocal.duplicate()
            self.iceFlex = self.hLocal.duplicate()
            self.iceHL.set(0.0)

        return

    def _buildIceDirection(self, h):
        """
        Compute the ice flow direction matrix.
        """

        # Define multiple flow directions for unfilled elevation
        ice_rcvID, _, ice_wghtVal = mfdreceivers(
            8, 1., h, self.sealevel
        )

        ice_rcvID[self.ghostIDs, :] = -1
        ice_wghtVal[self.ghostIDs, :] = 0

        # Set borders nodes
        if self.flatModel:
            ice_rcvID[self.idBorders, :] = np.tile(self.idBorders, (8, 1)).T
            ice_wghtVal[self.idBorders, :] = 0.0

        iceMat = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]

        for k in range(0, 8):
            # Flow direction matrix for a specific direction
            tmpMat = self._matrix_build()
            data = -ice_wghtVal[:, k].copy()
            data[ice_rcvID[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                ice_rcvID[:, k].astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            # Add the weights from each direction
            iceMat.axpy(1.0, tmpMat)
            tmpMat.destroy()

        if self.memclear:
            del data, indptr, nodes
            gc.collect()

        petsc4py.PETSc.garbage_cleanup()

        return iceMat.transpose()

    def iceAccumulation(self):
        """
        This function is the **main entry point** for ice accumulation computation.
        """

        ti = process_time()

        if self.flexOn:
            self.iceHL.copy(result=self.iceFlex)

        # Define a smoothed surface to compute ice flow
        hl = self.hLocal.getArray()

        # Get amount of ice
        rainA = self.bL.getArray()
        elaH = self.elaH(self.tNow)
        iceH = self.iceH(self.tNow)
        tmp = (hl - elaH) / (iceH - elaH)
        tmp[tmp > 1.] = 1.0
        iceA = np.multiply(rainA, tmp)

        self.hGlobal.copy(result=self.tmp1)
        smthH = self._hillSlope(smooth=1)

        # Get ice flow directions
        iceMat = self._buildIceDirection(smthH)

        # Solve ice accumulation
        self.tmpL.setArray(iceA)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self._solve_KSP(True, iceMat, self.tmp, self.iceFAG)
        self.dm.globalToLocal(self.iceFAG, self.iceFAL)
        iceMat.destroy()

        # Get corresponding ice thicknesses
        tmp = self.iceFAL.getArray()
        tmp[tmp < 0] = 0.
        self.iceFAL.setArray(tmp)
        self.dm.localToGlobal(self.iceFAL, self.iceFAG)
        tmp = self.iceFAL.getArray().copy()
        iceH = self.scaleIce * tmp**0.3
        iceH[tmp <= 1.e-3] = 0.

        # Smooth ice thicknesses
        self.tmpL.setArray(iceH + hl)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        smthIce = self._hillSlope(smooth=1)
        self.hGlobal.copy(result=self.tmp1)
        smthH = self._hillSlope(smooth=1)
        tmp = smthIce - smthH
        tmp[tmp < 0.1] = 0.
        self.iceHL.setArray(tmp)

        # Update ice accumulation
        smthIce = (tmp / self.scaleIce)**(1. / 0.3)
        self.iceFAL.setArray(smthIce)
        self.dm.localToGlobal(self.iceFAL, self.iceFAG)
        petsc4py.PETSc.garbage_cleanup()

        # If simulation starts then set the ice flex variable to the initial glacier thickness
        if self.tNow == self.tStart:
            self.iceHL.copy(result=self.iceFlex)

        if MPIrank == 0 and self.verbose:
            print(
                "Glaciers Accumulation (%0.02f seconds)" % (process_time() - ti),
                flush=True,
            )

        return
