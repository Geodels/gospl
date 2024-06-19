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

        self.paleoZ = None
        self.plateMov = 0
        self.paleoMov = 0
        self.paleoTime = None
        if self.platedata is not None:
            self.tecL = self.hLocal.duplicate()

        return

    def _readAdvectionData(self):
        """
        From a plate input file read the plate advection information, elevation and subsidence/uplift for each point of the mesh.
        """

        fplate = self.platedata.iloc[self.plateMov, 1]
        mdata = np.load(fplate)
        if "clust" in list(mdata.keys()):
            self.isCluster = mdata["clust"]
        else:
            self.isCluster = None
        if "cngbh" in list(mdata.keys()):
            self.clustNgbhs = mdata["cngbh"]
        else:
            self.clustNgbhs = None
        self.distNbghs = mdata["dngbh"]
        self.idNbghs = mdata["ingbh"]
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
            if self.platedata.iloc[nb, 0] == self.tNow:
                nb += 1
        if nb == self.plateMov:
            return

        # Check relevant file information
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

        if self.gflexOn:
            # Send local erosion deposition from flexure globally
            ced = self.cumEDFlex.getArray().copy()
            gcED = np.zeros(self.mpoints, dtype=np.float64) - 1.0e10
            gcED[self.locIDs] = ced
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gcED, op=MPI.MAX)
            
            # Send local flexural isostasy globally
            cfl = self.cumFlexL.getArray().copy()
            gcFL = np.zeros(self.mpoints, dtype=np.float64) - 1.0e10
            gcFL[self.locIDs] = cfl
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gcFL, op=MPI.MAX)

        if MPIrank == 0 and self.verbose:
            print(
                "Transfer local elevation/erosion globally (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        # Get relevant files information
        t0 = process_time()
        self._readAdvectionData()
        if MPIrank == 0 and self.verbose:
            print(
                "Read advection data (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )
        t0 = process_time()
        if self.clustNgbhs is not None:
            # For clustered points get heights of nearest neighbours
            idCluster = self.isCluster > 0
            tmp = gZ[idCluster]
            tmpngb = tmp[self.clustNgbhs]
            # Set new heights to the maximum height of nearest neighbours
            gZ[idCluster] = np.max(tmpngb, axis=1)

            # For clustered points get erosion/deposition of nearest neighbours
            tmp = gED[idCluster]
            tmpngb = tmp[self.clustNgbhs]
            # Set new erosion deposition to the maximum values of nearest neighbours
            gED[idCluster] = np.min(tmpngb, axis=1)

            if self.gflexOn:
                tmp = gcED[idCluster]
                tmpngb = tmp[self.clustNgbhs]
                gcED[idCluster] = np.min(tmpngb, axis=1)

                tmp = gcFL[idCluster]
                tmpngb = tmp[self.clustNgbhs]
                gcFL[idCluster] = np.min(tmpngb, axis=1)

        # Update elevation and erosion/deposition
        if self.interp == 1:
            nelev = gZ[self.idNbghs]
            nerodep = gED[self.idNbghs]
            if self.gflexOn:
               nf_ed = gcED[self.idNbghs]
               nflex = gcFL[self.idNbghs]
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
            if self.gflexOn:
                tmp = np.sum(weights * gcED[self.idNbghs], axis=1)
                nf_ed = np.divide(tmp, temp, out=np.zeros_like(temp), where=temp != 0)
                tmp = np.sum(weights * gcFL[self.idNbghs], axis=1)
                nflex = np.divide(tmp, temp, out=np.zeros_like(temp), where=temp != 0)
            if len(onIDs) > 0:
                nelev[onIDs] = gZ[self.idNbghs[onIDs, 0]]
                nerodep[onIDs] = gED[self.idNbghs[onIDs, 0]]
                if self.gflexOn:
                    nf_ed[onIDs] = gcED[self.idNbghs[onIDs, 0]]
                    nflex[onIDs] = gcFL[self.idNbghs[onIDs, 0]]

        # Remove marine erosion
        if self.fitMarine:
            seanIDs = np.where(np.logical_and(nelev < self.sealevel, nerodep < 0))[0]
            nerodep[seanIDs] = 0.0

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

        if self.gflexOn:
            self.cumEDFlex.setArray(nf_ed[self.locIDs])
            self.cumFlexL.setArray(nflex[self.locIDs])

        # Update stratigraphic record
        if self.stratNb > 0 and self.stratStep > 0:
            self._advectStrata(weights, onIDs)

        # Get the tectonic forcing from the paleo-reconstruction data
        self._getPaleoInfo()
        # If paleo-elevations are provided update elevations
        self._updatePaleoElev()

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

    def _updatePaleoElev(self):
        """
        Update surface information based on paleo-reconstruction.
        """

        if self.paleoZ is not None:

            if self.paleoTime is not None:
                if self.paleoTime > self.tNow:
                    return

            t0 = process_time()

            # Send local elevation globally
            hl = self.hLocal.getArray().copy()
            gZ = np.zeros(self.mpoints, dtype=np.float64) - 1.0e8
            gZ[self.locIDs] = hl
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)

            # Compute tectonic missmatch (will be outputed) in metres
            tec = self.paleoZ - gZ[self.locIDs]
            self.tecL.setArray(tec)

            # Only fit marine elevation
            if self.fitMarine:
                ids = np.where(gZ[self.locIDs] >= self.sealevel)[0]
                self.paleoZ[ids] = gZ[self.locIDs][ids]

            # Fit simulated elevations to paleoelevation ones
            self.hLocal.setArray(self.paleoZ)
            self.dm.localToGlobal(self.hLocal, self.hGlobal)

            if MPIrank == 0 and self.verbose:
                print(
                    "Update model based on paleo-elevations (%0.02f seconds)"
                    % (process_time() - t0),
                    flush=True,
                )

        return

    def _getPaleoInfo(self):
        """
        Retreive paleoelevation information from inputs (vertical movements and elevations).
        """

        if self.platedata is None:
            return

        t0 = process_time()

        # Check if we reached a plate advection time step
        nb = self.paleoMov

        if nb < len(self.platedata):
            if self.platedata.iloc[nb, 0] == self.tNow:
                nb += 1
        if nb == self.paleoMov:
            return

        # Get relevant file information
        if self.platedata.iloc[self.paleoMov, 2] == "empty":
            self.paleoMov = nb
            return

        # Read and load the file
        tplate = self.platedata.iloc[self.paleoMov, 2]
        mdata = np.load(tplate)
        # If tectonic conditions exist
        if "t" in list(mdata.keys()):
            self.uplift = mdata["t"][self.locIDs]
        else:
            self.uplift = None

        # If paleo-elevations information are provided
        if "z" in list(mdata.keys()):
            self.paleoZ = mdata["z"][self.locIDs].copy()
            if nb >= len(self.platedata):
                self.paleoTime = self.tEnd + self.dt
            else:
                self.paleoTime = self.platedata.iloc[nb, 0]
            if MPIrank == 0 and self.verbose:
                print(
                    "Store paleo-elevation values (%0.02f seconds)"
                    % (process_time() - t0),
                    flush=True,
                )

        del mdata

        # Update next paleoelevation dataset parameters
        self.paleoMov = nb

        return
