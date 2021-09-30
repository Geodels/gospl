import sys
import petsc4py
import numpy as np

from mpi4py import MPI
from scipy import spatial
from time import process_time

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class EarthPlate(object):

    """
    This class defines how spherical mesh surface is changing given a plate reconstruction model.
    """

    def __init__(self):
        """
        The initialisation of `EarthPlate` class.
        """

        self.plateMov = 0

        return

    def advectPlates(self):
        """
        Advect surface and stratigraphic information based on plate evolution.
        """

        if self.platedata is None:
            return

        t0 = process_time()
        # Check if we reached a plate advection time step
        nb = self.plateMov
        if nb < len(self.platedata):
            if self.platedata.iloc[nb, 0] < self.tNow + self.dt:
                nb += 1
        if nb == self.plateMov:
            return
        if self.platedata.iloc[self.plateMov, 1] == "empty":
            self.plateMov = nb
            return

        # dt = self.platedata.iloc[nb, 0] - self.platedata.iloc[self.plateMov, 0]
        # Plate advection starting time is define as an integer in Myrs
        # and follows gPlates convention positive when moving back in time
        # self.advecttime = -int(self.platedata.iloc[self.plateMov, 0] / 1.0e6)
        # Plate advection time step is defined as an integer in Myrs
        # self.advectdt = int(dt / 1.0e6)

        # From a plate input file read the plate advection information for each point of the mesh.
        fplate = self.platedata.iloc[self.plateMov, 1]
        mdata = np.load(fplate)
        # self.plateIds = mdata["iplate"]
        self.isCluster = mdata["clust"]
        self.clustNgbhs = mdata["cngbh"]
        self.distNbghs = mdata["dngbh"]
        self.idNbghs = mdata["ingbh"]

        # Send local elevation globally
        t0 = process_time()
        hl = self.hLocal.getArray().copy()
        gZ = np.zeros(self.mpoints, dtype=np.float64) - 1.0e8
        gZ[self.locIDs] = hl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)

        # Send local erosion deposition globally
        edl = self.cumEDLocal.getArray().copy()
        gED = np.zeros(self.mpoints, dtype=np.float64) - 1.0e10
        gED[self.locIDs] = edl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gED, op=MPI.MAX)

        if MPIrank == 0 and self.verbose:
            print(
                "Transfer local elevation/erosion globally (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        # For clustered points get heights of nearest neighbours
        t0 = process_time()
        idCluster = self.isCluster > 0
        heightsInCluster = gZ[idCluster]
        neighbourHeights = heightsInCluster[self.clustNgbhs]
        erodepInCluster = gED[idCluster]
        neighbourErodep = erodepInCluster[self.clustNgbhs]
        # Set new heights to the maximum height of nearest neighbours
        gZ[idCluster] = np.max(neighbourHeights, axis=1)
        # Set new erosion deposition to the maximum erosion of nearest neighbours
        gED[idCluster] = np.min(neighbourErodep, axis=1)

        # Update elevation and erosion/deposition
        if self.interp == 1:
            nelev = gZ[self.idNbghs]
            nerodep = gED[self.idNbghs]
        else:
            # Inverse weighting distance...
            weights = np.divide(
                1.0,
                self.distNbghs,
                out=np.zeros_like(self.distNbghs),
                where=self.distNbghs != 0,
            )
            temp = np.sum(weights, axis=1)
            onIDs = np.where(self.distNbghs[:, 0] == 0)[0]
            tmp = np.sum(weights * gZ[self.idNbghs], axis=1)
            nelev = np.divide(tmp, temp, out=np.zeros_like(temp), where=temp != 0)
            tmp = np.sum(weights * gED[self.idNbghs], axis=1)
            nerodep = np.divide(tmp, temp, out=np.zeros_like(temp), where=temp != 0)
            if len(onIDs) > 0:
                nelev[onIDs] = gZ[self.idNbghs[onIDs, 0]]
                nerodep[onIDs] = gED[self.idNbghs[onIDs, 0]]

        if MPIrank == 0 and self.verbose:
            print(
                "Define local elevation and erosion/deposition after advection (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        # From a plate input file read the associated uplift/subsidence for each point of the mesh.
        fplate = self.platedata.iloc[self.plateMov, 2]
        if fplate != "empty":
            mdata = np.load(fplate)
            nelev += mdata["z"]

        self.plateMov = nb
        if mdata is not None:
            del mdata

        self.hGlobal.setArray(nelev[self.glbIDs])
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        self.cumED.setArray(nerodep[self.glbIDs])
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)

        return
