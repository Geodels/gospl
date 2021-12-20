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
        if self.platedata is not None:
            self.tecL = self.hLocal.duplicate()

        return

    def _readAdvectionData(self):
        """
        From a plate input file read the plate advection information, elevation and subsidence/uplift for each point of the mesh.
        """

        fplate = self.platedata.iloc[self.plateMov, 1]
        if fplate != "empty":
            mdata = np.load(fplate)
            self.isCluster = mdata["clust"]
            self.clustNgbhs = mdata["cngbh"]
            self.distNbghs = mdata["dngbh"]
            self.idNbghs = mdata["ingbh"]
            del mdata

        tplate = "empty"
        tplate = self.platedata.iloc[self.plateMov, 2]
        if tplate != "empty":
            if self.plateMov < len(self.platedata) - 1:
                dstep = (
                    self.platedata.iloc[self.plateMov + 1, 0]
                    - self.platedata.iloc[self.plateMov, 0]
                )
            else:
                dstep = self.tEnd - self.platedata.iloc[self.plateMov, 0]
            mdata = np.load(tplate)
            self.uplift = mdata["t"][self.locIDs] / dstep
            del mdata

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

        # Get relevant files information
        self._readAdvectionData()

        # For clustered points get heights of nearest neighbours
        t0 = process_time()
        idCluster = self.isCluster > 0
        tmp = gZ[idCluster]
        tmpngb = tmp[self.clustNgbhs]
        # Set new heights to the maximum height of nearest neighbours
        gZ[idCluster] = np.max(tmpngb, axis=1)

        # For clustered points get erosion/deposition of nearest neighbours
        tmp = gED[idCluster]
        tmpngb = tmp[self.clustNgbhs]
        # Set new erosion deposition to the maximum erosion of nearest neighbours
        gED[idCluster] = np.min(tmpngb, axis=1)

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
        self.plateMov = nb

        self.hGlobal.setArray(nelev[self.glbIDs])
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        self.cumED.setArray(nerodep[self.glbIDs])
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)

        # Update stratigraphic record
        if self.stratNb > 0 and self.stratStep > 0:
            self._advectStrata(weights, onIDs)

        return

    def _updateStrataVars(self, reduceIDs, variable, sumw, onIDs, weights, nghs):

        # Reduce variable so that all values that needs to be read from another partition can be
        vals = np.zeros((self.mpoints, self.stratStep), dtype=np.float64) - 1.0e8
        vals[self.locIDs, :] = variable
        redVals = vals[reduceIDs, :]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, redVals, op=MPI.MAX)
        vals[reduceIDs, :] = redVals

        if self.interp == 1:
            nvar = vals[nghs, :]
        else:
            # Inverse weighting distance...
            nghs = self.idNbghs[self.locIDs]
            tmp = np.sum(weights * vals[nghs, :], axis=1)
            nvar = np.divide(tmp, sumw, where=sumw != 0, out=np.zeros_like(tmp))
            if len(onIDs) > 0:
                nvar[onIDs, :] = vals[nghs[onIDs, 0], :]

        return nvar

    def _findPts2Reduce(self):

        # Global neighbours for each local partition
        lgNghbs = self.idNbghs[self.locIDs, :].flatten()

        # Global IDs required locally but part of another partition
        gids = lgNghbs[~np.in1d(lgNghbs, self.locIDs)]

        # Get all points that have to be transferred
        vals = np.zeros(self.mpoints, dtype=int)
        vals[gids] = 1
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, vals, op=MPI.MAX)

        return np.where(vals > 0)[0]

    def _advectStrata(self, weights, onIDs):

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
        nghs = self.idNbghs[self.locIDs]
        wgts = weights[self.locIDs, :, None]
        sumw = np.sum(weights, axis=1)[self.locIDs, None]
        onID = np.where(self.distNbghs[self.locIDs, 0] == 0)[0]

        self.stratH[:, : self.stratStep] = self._updateStrataVars(
            redIDs, self.stratH[:, : self.stratStep], sumw, onID, wgts, nghs
        )
        self.stratZ[:, : self.stratStep] = self._updateStrataVars(
            redIDs, self.stratZ[:, : self.stratStep], sumw, onID, wgts, nghs
        )
        self.phiS[:, : self.stratStep] = self._updateStrataVars(
            redIDs, self.phiS[:, : self.stratStep], sumw, onID, wgts, nghs
        )
        if MPIrank == 0 and self.verbose:
            print(
                "Define local stratigraphy variables after advection (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return

    def forcePaleoElev(self):
        """
        Advect surface and stratigraphic information based on plate evolution.
        """

        if self.platedata is None:
            return

        if self.plateMov <= 0:
            return

        t0 = process_time()
        # Check if we reached a plate advection time step
        nb = self.plateMov

        if nb < len(self.platedata):
            if self.platedata.iloc[nb, 0] < self.tNow + self.dt:
                nb += 1
        if nb == self.plateMov:
            return

        if (self.platedata.iloc[self.plateMov, 1] == "empty") & (
            self.platedata.iloc[self.plateMov, 2] == "empty"
        ):
            return

        # Get relevant file information
        tplate = "empty"
        if self.plateMov > 0:
            tplate = self.platedata.iloc[self.plateMov - 1, 2]

        if tplate == "empty":
            return

        # Get the tectonic forcing required to fit with paleo-elevation model
        mdata = np.load(tplate)

        # If we want to fit all the paleo-elevation model
        if "z" in list(mdata.keys()):
            # Send local elevation globally
            t0 = process_time()
            hl = self.hLocal.getArray().copy()
            gZ = np.zeros(self.mpoints, dtype=np.float64) - 1.0e8
            gZ[self.locIDs] = hl
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)
            tec = mdata["z"][self.locIDs] - gZ[self.locIDs]
            self.tecL.setArray(tec)
            self.hGlobal.setArray(mdata["z"][self.glbIDs])
            self.dm.globalToLocal(self.hGlobal, self.hLocal)
            del mdata

            if MPIrank == 0 and self.verbose:
                print(
                    "Force paleo-elevation (%0.02f seconds)" % (process_time() - t0),
                    flush=True,
                )

        # If we want to fit the paleo-bathymetry model
        # if "zm" in list(mdata.keys()):
        #     # Send local bathymetry globally
        #     t0 = process_time()
        #     hl = self.hLocal.getArray().copy()
        #     gZ = np.zeros(self.mpoints, dtype=np.float64) - 1.0e8
        #     gZ[self.locIDs] = hl
        #     MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)
        #     # Only assign tectonic in the marine domain
        #     tec = mdata["zm"][self.locIDs] - gZ[self.locIDs]
        #     tec[gZ > self.sealevel] = 0.0
        #     self.tecL.setArray(tec)
        #     # Only force bathymentry values from paleo-elevation dataset
        #     gZ[gZ <= self.sealevel] = mdata["zm"][gZ <= self.sealevel]
        #     self.hGlobal.setArray(gZ[self.glbIDs])
        #     self.dm.globalToLocal(self.hGlobal, self.hLocal)
        #     del mdata
        #     if MPIrank == 0 and self.verbose:
        #         print(
        #             "Force paleo-bathymetry (%0.02f seconds)" % (process_time() - t0),
        #             flush=True,
        #         )
        return
