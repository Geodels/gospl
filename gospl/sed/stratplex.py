import os
import gc
import sys
import petsc4py
import numpy as np

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import stratasimple
    from gospl._fortran import stratabuild
    from gospl._fortran import stratabuildcarb

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class STRAMesh(object):
    """
    This class encapsulates all the functions related to underlying stratigraphic information. As mentionned previously, `gospl` has the ability to track different types of clastic sediment size and one type of carbonate (still under development). Sediment compaction in stratigraphic layers geometry and properties change are also considered.

    """

    def __init__(self):
        """
        The initialisation of `STRAMesh` class related to stratigraphic informations.
        """

        self.stratH = None
        self.stratF = None
        self.stratW = None
        self.stratZ = None
        self.stratC = None

        self.phiS = None
        self.phiF = None
        self.phiC = None

        return

    def readStratLayers(self):
        """
        When stratigraphic layers are turned on, this function reads any initial stratigraphic layers provided by within the YAML input file (key: `npstrata`).

        The following variables will be read from the file:

        - thickness of each stratigrapic layer `strataH` accounting for both erosion & deposition events.
        - proportion of fine sediment `strataF` contains in each stratigraphic layer.
        - proportion of weathered sediment `strataW` contains in each stratigraphic layer.
        - elevation at time of deposition, considered to be to the current elevation for the top stratigraphic layer `strataZ`.
        - porosity of coarse sediment `phiS` in each stratigraphic layer computed at center of each layer.
        - porosity of fine sediment `phiF` in each stratigraphic layer computed at center of each layer.
        - porosity of weathered sediment `phiW` in each stratigraphic layer computed at center of each layer.
        - proportion of carbonate sediment `strataC` contains in each stratigraphic layer if the carbonate module is turned on.
        - porosity of carbonate sediment `phiC` in each stratigraphic layer computed at center of each layer when the carbonate module is turned on.

        """

        if self.strataFile is not None:
            fileData = np.load(self.strataFile)
            stratVal = fileData["strataH"]
            self.initLay = stratVal.shape[1]
            self.stratNb += self.initLay

            # Create stratigraphic arrays
            self.stratH = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratH[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            stratVal = fileData["strataF"]
            self.stratF = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratF[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            stratVal = fileData["strataW"]
            self.stratW = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratW[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            stratVal = fileData["strataZ"]
            self.stratZ = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratZ[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            stratVal = fileData["phiS"]
            self.phiS = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.phiS[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            stratVal = fileData["phiF"]
            self.phiF = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.phiF[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            stratVal = fileData["phiW"]
            self.phiW = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.phiW[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            if self.carbOn:
                stratVal = fileData["strataC"]
                self.stratC = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
                self.stratC[:, 0 : self.initLay] = stratVal[
                    self.locIDs, 0 : self.initLay
                ]

                stratVal = fileData["phiC"]
                self.phiC = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
                self.phiC[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            if self.memclear:
                del fileData, stratVal
                gc.collect()
        else:
            self.stratH = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.phiS = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratZ = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)

        return

    def deposeStrat(self, stype):
        """
        Add deposition on top of an existing stratigraphic layer. The following variables will be recorded:

        - thickness of each stratigrapic layer `stratH` accounting for both erosion & deposition events.
        - proportion of fine sediment `stratF` contains in each stratigraphic layer.
        - proportion of weathered sediment `stratW` contains in each stratigraphic layer.
        - porosity of coarse sediment `phiS` in each stratigraphic layer computed at center of each layer.
        - porosity of fine sediment `phiF` in each stratigraphic layer computed at center of each layer.
        - porosity of weathered sediment `phiW` in each stratigraphic layer computed at center of each layer.
        - proportion of carbonate sediment `stratC` contains in each stratigraphic layer if the carbonate module is turned on.
        - porosity of carbonate sediment `phiC` in each stratigraphic layer computed at center of each layer when the carbonate module is turned on.

        :arg stype: sediment type (integer)
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        depo = self.tmpL.getArray().copy()
        depo[depo < 1.0e-4] = 0.0
        if self.stratF is not None:
            fineH = self.stratH[:, self.stratStep] * self.stratF[:, self.stratStep]
        if self.stratW is not None:
            clayH = self.stratH[:, self.stratStep] * self.stratW[:, self.stratStep]
        if self.carbOn:
            carbH = self.stratH[:, self.stratStep] * self.stratC[:, self.stratStep]
        self.stratH[:, self.stratStep] += depo
        ids = np.where(depo > 0)[0]

        if stype == 0:
            self.phiS[ids, self.stratStep] = self.phi0s
        elif stype == 1:
            fineH[ids] += depo[ids]
            self.phiF[ids, self.stratStep] = self.phi0f
        elif stype == 2:
            carbH[ids] += depo[ids]
            self.phiC[ids, self.stratStep] = self.phi0c
        elif stype == 3:
            clayH[ids] += depo[ids]
            self.phiW[ids, self.stratStep] = self.phi0w

        if self.stratF is not None:
            self.stratF[ids, self.stratStep] = (
                fineH[ids] / self.stratH[ids, self.stratStep]
            )
            if self.memclear:
                del fineH
        if self.stratW is not None:
            self.stratW[ids, self.stratStep] = (
                clayH[ids] / self.stratH[ids, self.stratStep]
            )
            if self.memclear:
                del clayH
        if self.carbOn:
            self.stratC[ids, self.stratStep] = (
                carbH[ids] / self.stratH[ids, self.stratStep]
            )
            if self.memclear:
                del carbH

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
            return

        # Cumulative thickness for each node
        self.stratH[nids, 0] += 1.0e6
        cumThick = np.cumsum(self.stratH[nids, self.stratStep :: -1], axis=1)[:, ::-1]
        boolMask = cumThick < -ero[nids].reshape((len(nids), 1))
        mask = boolMask.astype(int)

        if self.stratF is not None:
            # Get fine sediment eroded from river incision
            thickF = (
                self.stratH[nids, 0 : self.stratStep + 1]
                * self.stratF[nids, 0 : self.stratStep + 1]
            )
            # From fine thickness extract the solid phase that is eroded
            thFine = thickF * (1.0 - self.phiF[nids, 0 : self.stratStep + 1])
            thFine = np.sum((thFine * mask), axis=1)

        if self.stratW is not None:
            # Get weathered sediment eroded from river incision
            thickW = (
                self.stratH[nids, 0 : self.stratStep + 1]
                * self.stratW[nids, 0 : self.stratStep + 1]
            )
            # From weathered thickness extract the solid phase that is eroded
            thClay = thickW * (1.0 - self.phiW[nids, 0 : self.stratStep + 1])
            thClay = np.sum((thClay * mask), axis=1)

        # Get carbonate sediment eroded from river incision
        if self.carbOn:
            thickC = (
                self.stratH[nids, 0 : self.stratStep + 1]
                * self.stratC[nids, 0 : self.stratStep + 1]
            )
            # From carbonate thickness extract the solid phase that is eroded
            thCarb = thickC * (1.0 - self.phiC[nids, 0 : self.stratStep + 1])
            thCarb = np.sum((thCarb * mask), axis=1)
            # From sand thickness extract the solid phase that is eroded
            thickS = (
                self.stratH[nids, 0 : self.stratStep + 1] - thickF - thickC - thickW
            )
            thCoarse = thickS * (1.0 - self.phiS[nids, 0 : self.stratStep + 1])
            thCoarse = np.sum((thCoarse * mask), axis=1)
        else:
            # From sand thickness extract the solid phase that is eroded
            if self.stratF is not None:
                thickS = self.stratH[nids, 0 : self.stratStep + 1] - thickF - thickW
            else:
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

        if self.stratF is not None:
            # Get thickness of each sediment type eroded in the remaining layer
            self.thFine = np.zeros(self.lpoints)
            # From fine thickness extract the solid phase that is eroded from this last layer
            tmp = (self.stratH[nids, eroLayNb] - eroVal) * self.stratF[nids, eroLayNb]
            tmp[tmp < 1.0e-8] = 0.0
            thFine += tmp * (1.0 - self.phiF[nids, eroLayNb])
            # Define the uncompacted fine thickness that will be deposited dowstream
            self.thFine[nids] = thFine / (1.0 - self.phi0f)
            self.thFine[self.thFine < 0.0] = 0.0

        if self.stratW is not None:
            # Get thickness of each sediment type eroded in the remaining layer
            self.thClay = np.zeros(self.lpoints)
            # From weathered thickness extract the solid phase that is eroded from this last layer
            tmp = (self.stratH[nids, eroLayNb] - eroVal) * self.stratW[nids, eroLayNb]
            tmp[tmp < 1.0e-8] = 0.0
            thClay += tmp * (1.0 - self.phiW[nids, eroLayNb])
            # Define the uncompacted weathered thickness that will be deposited dowstream
            self.thClay[nids] = thClay / (1.0 - self.phi0w)
            self.thClay[self.thClay < 0.0] = 0.0

        self.thCoarse = np.zeros(self.lpoints)
        if self.carbOn:
            # From carb thickness extract the solid phase that is eroded from this last layer
            self.thCarb = np.zeros(self.lpoints)
            tmp = (self.stratH[nids, eroLayNb] - eroVal) * self.stratC[nids, eroLayNb]
            tmp[tmp < 1.0e-8] = 0.0
            thCarb += tmp * (1.0 - self.phiC[nids, eroLayNb])
            # Define the uncompacted carbonate thickness that will be deposited dowstream
            self.thCarb[nids] = thCarb / (1.0 - self.phi0c)
            self.thCarb[self.thCarb < 0.0] = 0.0
            # From sand thickness extract the solid phase that is eroded from this last layer
            tmp = self.stratH[nids, eroLayNb] - eroVal
            tmp *= (
                1.0
                - self.stratC[nids, eroLayNb]
                - self.stratW[nids, eroLayNb]
                - self.stratF[nids, eroLayNb]
            )
            # Define the uncompacted sand thickness that will be deposited dowstream
            thCoarse += tmp * (1.0 - self.phiS[nids, eroLayNb])
            self.thCoarse[nids] = thCoarse / (1.0 - self.phi0s)
            self.thCoarse[self.thCoarse < 0.0] = 0.0
        else:
            # From sand thickness extract the solid phase that is eroded from this last layer
            tmp = self.stratH[nids, eroLayNb] - eroVal
            if self.stratF is not None:
                tmp *= 1.0 - self.stratF[nids, eroLayNb] - self.stratW[nids, eroLayNb]
            tmp[tmp < 1.0e-8] = 0.0
            # Define the uncompacted sand thickness that will be deposited dowstream
            thCoarse += tmp * (1.0 - self.phiS[nids, eroLayNb])
            self.thCoarse[nids] = thCoarse / (1.0 - self.phi0s)
            self.thCoarse[self.thCoarse < 0.0] = 0.0

        # Update thickness of top stratigraphic layer
        self.stratH[nids, eroLayNb] = eroVal
        self.stratH[nids, 0] -= 1.0e6
        self.stratH[self.stratH <= 0] = 0.0
        self.phiS[self.stratH <= 0] = 0.0
        if self.stratF is not None:
            self.stratF[self.stratH <= 0] = 0.0
            self.phiF[self.stratH <= 0] = 0.0
        if self.stratW is not None:
            self.stratW[self.stratH <= 0] = 0.0
            self.phiW[self.stratH <= 0] = 0.0
        if self.carbOn:
            self.stratC[self.stratH <= 0] = 0.0
            self.phiC[self.stratH <= 0] = 0.0

        self.thCoarse /= self.dt
        if self.stratF is not None:
            self.thFine /= self.dt
        if self.stratW is not None:
            self.thClay /= self.dt
        if self.carbOn:
            self.thCarb /= self.dt

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
        if self.stratF is not None:
            phiF = self.phi0f * np.exp(depth / self.z0f)
            phiF = np.minimum(phiF, self.phiF[:, : self.stratStep + 1])
            phiW = self.phi0w * np.exp(depth / self.z0w)
            phiW = np.minimum(phiW, self.phiW[:, : self.stratStep + 1])
        if self.carbOn:
            phiC = self.phi0c * np.exp(depth / self.z0c)
            phiC = np.minimum(phiC, self.phiC[:, : self.stratStep + 1])

        # Compute the solid phase in each layers
        if self.stratF is not None:
            tmpF = (
                self.stratH[:, : self.stratStep + 1]
                * self.stratF[:, : self.stratStep + 1]
            )
            tmpF *= 1.0 - self.phiF[:, : self.stratStep + 1]
            tmpW = (
                self.stratH[:, : self.stratStep + 1]
                * self.stratW[:, : self.stratStep + 1]
            )
            tmpW *= 1.0 - self.phiW[:, : self.stratStep + 1]

        if self.carbOn:
            tmpC = (
                self.stratH[:, : self.stratStep + 1]
                * self.stratC[:, : self.stratStep + 1]
            )
            tmpC *= 1.0 - self.phiC[:, : self.stratStep + 1]
            tmpS = (
                self.stratC[:, : self.stratStep + 1]
                + self.stratF[:, : self.stratStep + 1]
                + self.stratW[:, : self.stratStep + 1]
            )
            tmpS = self.stratH[:, : self.stratStep + 1] * (1.0 - tmpS)
            tmpS *= 1.0 - self.phiS[:, : self.stratStep + 1]
            solidPhase = tmpC + tmpS + tmpF + tmpW
        else:
            if self.stratF is not None:
                tmpS = self.stratH[:, : self.stratStep + 1] * (
                    1.0
                    - self.stratF[:, : self.stratStep + 1]
                    - self.stratW[:, : self.stratStep + 1]
                )
                tmpS *= 1.0 - self.phiS[:, : self.stratStep + 1]
                solidPhase = tmpS + tmpF + tmpW
            else:
                tmpS = self.stratH[:, : self.stratStep + 1]
                tmpS *= 1.0 - self.phiS[:, : self.stratStep + 1]
                solidPhase = tmpS

        # Get new layer thickness after porosity change
        if self.stratF is not None:
            tmpF = self.stratF[:, : self.stratStep + 1] * (
                1.0 - phiF[:, : self.stratStep + 1]
            )
            tmpW = self.stratW[:, : self.stratStep + 1] * (
                1.0 - phiW[:, : self.stratStep + 1]
            )
        if self.carbOn:
            tmpC = self.stratC[:, : self.stratStep + 1] * (
                1.0 - phiC[:, : self.stratStep + 1]
            )
            tmpS = (
                1.0
                - self.stratF[:, : self.stratStep + 1]
                - self.stratC[:, : self.stratStep + 1]
                - self.stratW[:, : self.stratStep + 1]
            )
            tmpS *= 1.0 - phiS[:, : self.stratStep + 1]
            tot = tmpS + tmpC + tmpF + tmpW
        else:
            if self.stratF is not None:
                tmpS = (
                    1.0
                    - self.stratF[:, : self.stratStep + 1]
                    - self.stratW[:, : self.stratStep + 1]
                ) * (1.0 - phiS[:, : self.stratStep + 1])
                tot = tmpS + tmpF + tmpW
            else:
                tot = 1.0 - phiS[:, : self.stratStep + 1]

        ids = np.where(tot > 0.0)
        newH = np.zeros(tot.shape)
        newH[ids] = solidPhase[ids] / tot[ids]
        newH[newH <= 0] = 0.0
        phiS[newH <= 0] = 0.0
        if self.stratF is not None:
            phiF[newH <= 0] = 0.0
            phiW[newH <= 0] = 0.0
        if self.carbOn:
            phiC[newH <= 0] = 0.0

        # Update porosities in each sedimentary layer
        self.phiS[:, : self.stratStep + 1] = phiS
        if self.stratF is not None:
            self.phiF[:, : self.stratStep + 1] = phiF
            self.phiW[:, : self.stratStep + 1] = phiW
        if self.carbOn:
            self.phiC[:, : self.stratStep + 1] = phiC

        if self.memclear:
            del phiS, solidPhase
            del ids, tmpS, tot
            if self.stratF is not None:
                del tmpF, phiF, tmpW, phiW
            if self.carbOn:
                del phiC, tmpC
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

    def updateDispStrata(self):
        """
        Gets the stratigraphic records relevant to each partition after mesh advection.

        The functions returns local stratigraphic layer information:

        - thickness of each stratigrapic layer `loc_stratH` accounting for both erosion & deposition events.
        - proportion of fine sediment `loc_stratF` contains in each stratigraphic layer.
        - proportion of weathered sediment `loc_stratW` contains in each stratigraphic layer.
        - elevation at time of deposition, considered to be to the current elevation for the top stratigraphic layer `loc_stratZ`.
        - porosity of coarse sediment `loc_phiS` in each stratigraphic layer computed at center of each layer.
        - porosity of fine sediment `loc_phiF` in each stratigraphic layer computed at center of each layer.
        - porosity of waethered sediment `loc_phiW` in each stratigraphic layer computed at center of each layer.
        - proportion of carbonate sediment `loc_strataC` contains in each stratigraphic layer if the carbonate module is turned on.
        - porosity of carbonate sediment `loc_phiC` in each stratigraphic layer computed at center of each layer when the carbonate module is turned on.

        .. note::

            In `gospl`, the stratigraphic layers are only defined locally. For interpolation on the edges of each partition it is important to ensure that all stratigraphic information on the adjacent nodes of the neighbouring partition are accessible. This is done by applying MPI `Allreduce` operation to the nodes parts of the overlaid ('shadow') zone.

        :return: loc_stratH, loc_stratZ, loc_stratF, loc_stratW, loc_stratC, loc_phiS, loc_phiF, loc_phiW, loc_phiC
        """

        stratH = self.stratH[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp = np.full((self.shadowgNb, self.stratStep), -1.0e8, dtype=np.float64)
            temp[self.gshadinIDs, :] = stratH[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            stratH[self.shadowAlls, :] = temp
            loc_stratH = stratH[self.locIDs, :]
        else:
            loc_stratH = stratH[self.locIDs, :]

        stratZ = self.stratZ[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = stratZ[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            stratZ[self.shadowAlls, :] = temp
            loc_stratZ = stratZ[self.locIDs, :]
        else:
            loc_stratZ = stratZ[self.locIDs, :]

        phiS = self.phiS[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = phiS[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            phiS[self.shadowAlls, :] = temp
            loc_phiS = phiS[self.locIDs, :]
        else:
            loc_phiS = phiS[self.locIDs, :]

        if self.stratF is not None:
            stratF = self.stratF[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = stratF[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                stratF[self.shadowAlls, :] = temp
                loc_stratF = stratF[self.locIDs, :]
            else:
                loc_stratF = stratF[self.locIDs, :]

            phiF = self.phiF[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = phiF[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                phiF[self.shadowAlls, :] = temp
                loc_phiF = phiF[self.locIDs, :]
            else:
                loc_phiF = phiF[self.locIDs, :]
        else:
            loc_stratF = None
            loc_phiF = None

        if self.stratW is not None:
            stratW = self.stratW[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = stratW[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                stratW[self.shadowAlls, :] = temp
                loc_stratW = stratW[self.locIDs, :]
            else:
                loc_stratW = stratW[self.locIDs, :]

            phiW = self.phiW[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = phiW[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                phiW[self.shadowAlls, :] = temp
                loc_phiW = phiW[self.locIDs, :]
            else:
                loc_phiW = phiW[self.locIDs, :]
        else:
            loc_stratW = None
            loc_phiW = None

        if self.carbOn:
            stratC = self.stratC[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = stratZ[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                stratC[self.shadowAlls, :] = temp
                loc_stratC = stratC[self.locIDs, :]
            else:
                loc_stratC = stratC[self.locIDs, :]

            phiC = self.phiC[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = phiC[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                phiC[self.shadowAlls, :] = temp
                loc_phiC = phiC[self.locIDs, :]
            else:
                loc_phiC = phiC[self.locIDs, :]

            return (
                loc_stratH,
                loc_stratZ,
                loc_stratF,
                loc_stratW,
                loc_stratC,
                loc_phiS,
                loc_phiF,
                loc_phiW,
                loc_phiC,
            )
        else:
            return (
                loc_stratH,
                loc_stratZ,
                loc_stratF,
                loc_stratW,
                loc_phiS,
                loc_phiF,
                loc_phiW,
            )

    def stratalRecord(
        self,
        indices,
        weights,
        loc_stratH,
        loc_stratZ,
        loc_stratF,
        loc_stratW,
        loc_stratC,
        loc_phiS,
        loc_phiF,
        loc_phiW,
        loc_phiC,
    ):
        """
        Once the interpolation has been performed, the following function updates the stratigraphic records based on the advected mesh.

        The function relies on 3 fortran subroutines (for loop performance purposes):

        1. stratasimple
        2. stratabuild
        3. stratabuildcarb

        :arg indices: indices of the closest nodes used for interpolation
        :arg weights: weights based on the distances to closest nodes
        :arg loc_stratH: thickness of each stratigrapic layer accounting for both erosion & deposition events.
        :arg loc_stratF: proportion of fine sediment contains in each stratigraphic layer.
        :arg loc_stratW: proportion of weathered sediment contains in each stratigraphic layer.
        :arg loc_stratZ: elevation at time of deposition, considered to be to the current elevation for the top stratigraphic layer.
        :arg loc_phiS: porosity of coarse sediment in each stratigraphic layer computed at center of each layer.
        :arg loc_phiF: porosity of fine sediment in each stratigraphic layer computed at center of each layer.
        :arg loc_phiW: porosity of weathered sediment in each stratigraphic layer computed at center of each layer.
        :arg loc_strataC: proportion of carbonate sediment contains in each stratigraphic layer if the carbonate module is turned on.
        :arg loc_phiC: porosity of carbonate sediment in each stratigraphic layer computed at center of each layer when the carbonate module is turned on.

        """

        if self.carbOn:
            (
                self.stratH[:, : self.stratStep],
                self.stratZ[:, : self.stratStep],
                self.stratF[:, : self.stratStep],
                self.stratW[:, : self.stratStep],
                self.stratC[:, : self.stratStep],
                self.phiS[:, : self.stratStep],
                self.phiF[:, : self.stratStep],
                self.phiW[:, : self.stratStep],
                self.phiC[:, : self.stratStep],
            ) = stratabuildcarb(
                self.lpoints,
                self.stratStep,
                indices,
                weights,
                loc_stratH,
                loc_stratZ,
                loc_stratF,
                loc_stratW,
                loc_stratC,
                loc_phiS,
                loc_phiF,
                loc_phiW,
                loc_phiC,
            )
        elif self.stratF is not None:
            (
                self.stratH[:, : self.stratStep],
                self.stratZ[:, : self.stratStep],
                self.stratF[:, : self.stratStep],
                self.stratW[:, : self.stratStep],
                self.phiS[:, : self.stratStep],
                self.phiF[:, : self.stratStep],
                self.phiW[:, : self.stratStep],
            ) = stratabuild(
                self.lpoints,
                self.stratStep,
                indices,
                weights,
                loc_stratH,
                loc_stratZ,
                loc_stratF,
                loc_stratW,
                loc_phiS,
                loc_phiF,
                loc_phiW,
            )
        else:
            (
                self.stratH[:, : self.stratStep],
                self.stratZ[:, : self.stratStep],
                self.phiS[:, : self.stratStep],
            ) = stratasimple(
                self.lpoints,
                self.stratStep,
                indices,
                weights,
                loc_stratH,
                loc_stratZ,
                loc_phiS,
            )

        return

    def localStrat(self):
        """
        Updates stratigraphic records after mesh advection on the edges of each partition to ensure that all stratigraphic information on the adjacent nodes of the neighbouring partition are equals on all processors sharing a common number of nodes.
        """

        stratH = self.stratH[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp = np.full((self.shadowgNb, self.stratStep), -1.0e8, dtype=np.float64)
            temp[self.gshadinIDs, :] = stratH[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            stratH[self.shadowAlls, :] = temp
            self.stratH[:, : self.stratStep] = stratH[self.locIDs, :]
        else:
            self.stratH[:, : self.stratStep] = stratH[self.locIDs, :]
        stratZ = self.stratZ[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = stratZ[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            stratZ[self.shadowAlls, :] = temp
            self.stratZ[:, : self.stratStep] = stratZ[self.locIDs, :]
        else:
            self.stratZ[:, : self.stratStep] = stratZ[self.locIDs, :]
        phiS = self.phiS[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = phiS[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            phiS[self.shadowAlls, :] = temp
            self.phiS[:, : self.stratStep] = phiS[self.locIDs, :]
        else:
            self.phiS[:, : self.stratStep] = phiS[self.locIDs, :]
        if self.memclear:
            del stratH, stratZ, phiS

        if self.stratF is not None:
            stratF = self.stratF[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = stratF[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                stratF[self.shadowAlls, :] = temp
                self.stratF[:, : self.stratStep] = stratF[self.locIDs, :]
            else:
                self.stratF[:, : self.stratStep] = stratF[self.locIDs, :]
            phiF = self.phiF[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = phiF[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                phiF[self.shadowAlls, :] = temp
                self.phiF[:, : self.stratStep] = phiF[self.locIDs, :]
            else:
                self.phiF[:, : self.stratStep] = phiF[self.locIDs, :]
            if self.memclear:
                del stratF, phiF

        if self.stratW is not None:
            stratW = self.stratW[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = stratW[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                stratW[self.shadowAlls, :] = temp
                self.stratW[:, : self.stratStep] = stratW[self.locIDs, :]
            else:
                self.stratW[:, : self.stratStep] = stratW[self.locIDs, :]
            phiW = self.phiW[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = phiW[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                phiW[self.shadowAlls, :] = temp
                self.phiW[:, : self.stratStep] = phiW[self.locIDs, :]
            else:
                self.phiW[:, : self.stratStep] = phiW[self.locIDs, :]
            if self.memclear:
                del stratW, phiW

        if self.carbOn:
            stratC = self.stratC[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = stratC[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                stratC[self.shadowAlls, :] = temp
                self.stratC[:, : self.stratStep] = stratC[self.locIDs, :]
            else:
                self.stratC[:, : self.stratStep] = stratC[self.locIDs, :]
            phiC = self.phiC[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = phiC[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                phiC[self.shadowAlls, :] = temp
                self.phiC[:, : self.stratStep] = phiC[self.locIDs, :]
            else:
                self.phiC[:, : self.stratStep] = phiC[self.locIDs, :]
            if self.memclear:
                del stratC, phiC

        return
