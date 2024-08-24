import os
import gc
import sys
import petsc4py
import numpy as np

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import strataonesed

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class STRAMesh(object):
    """
    This class encapsulates all the functions related to underlying stratigraphic information. Sediment compaction in stratigraphic layers geometry and properties change are also considered.

    """

    def __init__(self):
        """
        The initialisation of `STRAMesh` class related to stratigraphic informations.
        """

        self.stratH = None
        self.stratZ = None
        self.phiS = None

        return

    def readStratLayers(self):
        """
        When stratigraphic layers are turned on, this function reads any initial stratigraphic layers provided by within the YAML input file (key: `npstrata`).

        The following variables will be read from the file:

        - thickness of each stratigrapic layer `strataH` accounting for both erosion & deposition events.
        - elevation at time of deposition, considered to be to the current elevation for the top stratigraphic layer `strataZ`.
        - porosity of coarse sediment `phiS` in each stratigraphic layer computed at center of each layer.
        """

        if self.strataFile is not None:
            fileData = np.load(self.strataFile)
            stratVal = fileData["strataH"]
            self.initLay = stratVal.shape[1]
            self.stratNb += self.initLay

            # Create stratigraphic arrays
            self.stratH = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratH[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            stratVal = fileData["strataZ"]
            self.stratZ = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratZ[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            stratVal = fileData["phiS"]
            self.phiS = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.phiS[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            if self.memclear:
                del fileData, stratVal
                gc.collect()
        else:
            self.stratH = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.phiS = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratZ = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratH[:, 0] = 1.0e6
            self.phiS[:, 0] = self.phi0s

        return

    def deposeStrat(self):
        """
        Add deposition on top of an existing stratigraphic layer. The following variables will be recorded:

        - thickness of each stratigrapic layer `stratH` accounting for both erosion & deposition events.
        - porosity of sediment `phiS` in each stratigraphic layer computed at center of each layer.
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        depo = self.tmpL.getArray().copy()
        depo[depo < 1.0e-4] = 0.0
        self.stratH[:, self.stratStep] += depo
        ids = np.where(depo > 0)[0]
        self.phiS[ids, self.stratStep] = self.phi0s

        # Cleaning arrays
        if self.memclear:
            del depo, ids
            gc.collect()

        return

    def erodeStrat(self):
        """
        This function removes eroded sediment thicknesses from the stratigraphic pile. The function takes into account the porosity values of considered lithologies in each eroded stratigraphic layers.

        It follows the following assumptions:

        - Eroded thicknesses from stream power law and hillslope diffusion are considered to encompass both the solid and void phase.
        - Only the solid phase will be moved dowstream by surface processes.
        - The corresponding deposit thicknesses for those freshly eroded sediments correspond to uncompacted thicknesses based on the porosity at surface given from the input file.
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        ero = self.tmpL.getArray().copy()
        ero[ero > 0] = 0.0

        # Nodes experiencing erosion
        nids = np.where(ero < 0)[0]
        if len(nids) == 0:
            self.thCoarse = np.zeros(self.lpoints)
            return

        # Cumulative thickness for each node
        self.stratH[nids, 0] += 1.0e6
        cumThick = np.cumsum(self.stratH[nids, self.stratStep :: -1], axis=1)[:, ::-1]
        boolMask = cumThick < -ero[nids].reshape((len(nids), 1))
        mask = boolMask.astype(int)

        thickS = self.stratH[nids, 0 : self.stratStep + 1]
        thCoarse = thickS * (1.0 - self.phiS[nids, 0 : self.stratStep + 1])
        thCoarse = np.sum((thCoarse * mask), axis=1)

        # Clear all stratigraphy points which are eroded
        cumThick[boolMask] = 0.0
        tmp = self.stratH[nids, : self.stratStep + 1]
        tmp[boolMask] = 0
        self.stratH[nids, : self.stratStep + 1] = tmp

        # Erode remaining stratal layers
        # Get non-zero top layer number
        eroLayNb = np.bincount(np.nonzero(cumThick)[0]) - 1
        eroVal = cumThick[np.arange(len(nids)), eroLayNb] + ero[nids]

        self.thCoarse = np.zeros(self.lpoints)
        # From sand thickness extract the solid phase that is eroded from this last layer
        tmp = self.stratH[nids, eroLayNb] - eroVal
        tmp[tmp < 1.0e-8] = 0.0
        # Define the uncompacted sand thickness that will be deposited dowstream
        thCoarse += tmp * (1.0 - self.phiS[nids, eroLayNb])
        self.thCoarse[nids] = thCoarse / (1.0 - self.phi0s)
        self.thCoarse[self.thCoarse < 0.0] = 0.0

        # Update thickness of top stratigraphic layer
        self.stratH[nids, eroLayNb] = eroVal
        self.stratH[nids, 0] -= 1.0e6
        self.stratH[self.stratH < 0] = 0.0
        self.phiS[self.stratH < 0] = 0.0
        self.thCoarse /= self.dt

        return

    def elevStrat(self):
        """
        This function updates the current stratigraphic layer elevation.
        """

        self.stratZ[:, self.stratStep] = self.hLocal.getArray()

        return

    def _depthPorosity(self, depth):
        """
        This function uses the depth-porosity relationships to compute the porosities for each lithology and then the solid phase to get each layer thickness changes due to compaction.

        .. note::

            We assume that porosity cannot increase after unloading.

        :arg depth: depth below basement for each sedimentary layer

        :return: newH updated sedimentary layer thicknesses after compaction
        """

        # Depth-porosity functions
        phiS = self.phi0s * np.exp(depth / self.z0s)
        phiS = np.minimum(phiS, self.phiS[:, : self.stratStep + 1])

        # Compute the solid phase in each layers
        tmpS = self.stratH[:, : self.stratStep + 1]
        tmpS *= 1.0 - self.phiS[:, : self.stratStep + 1]
        solidPhase = tmpS

        # Get new layer thickness after porosity change
        tot = 1.0 - phiS[:, : self.stratStep + 1]

        ids = np.where(tot > 0.0)
        newH = np.zeros(tot.shape)
        newH[ids] = solidPhase[ids] / tot[ids]
        newH[newH <= 0] = 0.0
        phiS[newH <= 0] = 0.0

        # Update porosities in each sedimentary layer
        self.phiS[:, : self.stratStep + 1] = phiS

        if self.memclear:
            del phiS, solidPhase
            del ids, tmpS, tot
            gc.collect()

        return newH

    def getCompaction(self):
        """
        This function computes the changes in sedimentary layers porosity and thicknesses due to compaction.

        .. note::

            We assume simple depth-porosiy relationships for each sediment type available in each layers.
        """

        t0 = process_time()
        topZ = np.vstack(self.hLocal.getArray())
        totH = np.sum(self.stratH[:, : self.stratStep + 1], axis=1)

        # Height of the sediment column above the center of each layer is given by
        cumZ = -np.cumsum(self.stratH[:, self.stratStep :: -1], axis=1) + topZ
        elev = np.append(topZ, cumZ[:, :-1], axis=1)
        zlay = np.fliplr(elev - np.fliplr(self.stratH[:, : self.stratStep + 1] / 2.0))

        # Compute lithologies porosities for each depositional layers
        # Get depth below basement
        depth = zlay - topZ

        # Now using depth-porosity relationships we compute the porosities
        newH = self._depthPorosity(depth)

        # Get the total thickness changes induced by compaction and
        # update the elevation accordingly
        dz = totH - np.sum(newH, axis=1)
        dz[dz <= 0] = 0.0
        self.hLocal.setArray(topZ.flatten() - dz.flatten())
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        # Update each layer thicknesses
        self.stratH[:, : self.stratStep + 1] = newH
        if self.memclear:
            del dz, newH, totH, topZ
            del depth, zlay, cumZ, elev
            gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Lithology Porosity Values (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return

    def stratalRecord(self, indices, weights, onIDs):
        """
        Once the interpolation has been performed, the following function updates the stratigraphic records based on the advected mesh.

        The function relies on fortran subroutines strataonesed.

        :arg indices: indices of the closest nodes used for interpolation
        :arg weights: weights based on the distances to closest nodes
        :arg onIDs: index of nodes remaining at the same position.

        """

        # Get local stratal dataset after displacements
        loc_stratH = self.stratH[:, : self.stratStep]
        loc_stratZ = self.stratZ[:, : self.stratStep]
        loc_phiS = self.phiS[:, : self.stratStep]
        nstratH, nstratZ, nphiS = strataonesed(
            self.lpoints,
            self.stratStep,
            indices,
            weights,
            loc_stratH,
            loc_stratZ,
            loc_phiS,
        )

        if len(onIDs) > 0:
            nstratZ[onIDs, :] = loc_stratZ[indices[onIDs, 0], :]
            nstratH[onIDs, :] = loc_stratH[indices[onIDs, 0], :]
            nphiS[onIDs, :] = loc_phiS[indices[onIDs, 0], :]

        # Updates stratigraphic records after mesh advection on the edges of each partition
        # to ensure that all stratigraphic information on the adjacent nodes of the neighbouring
        # partition are equals on all processors sharing a common number of nodes.
        for k in range(self.stratStep):
            self.tmp.setArray(nstratZ[:, k])
            self.dm.globalToLocal(self.tmp, self.tmpL)
            self.stratZ[:, k] = self.tmpL.getArray().copy()

            self.tmp.setArray(nstratH[:, k])
            self.dm.globalToLocal(self.tmp, self.tmpL)
            self.stratH[:, k] = self.tmpL.getArray().copy()
            self.tmp.setArray(nphiS[:, k])
            self.dm.globalToLocal(self.tmp, self.tmpL)
            self.phiS[:, k] = self.tmpL.getArray().copy()

        return
