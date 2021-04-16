import os
import gc
import sys
import h5py
import glob
import shutil
import petsc4py
import numpy as np
from scipy import spatial
from scipy import ndimage

from mpi4py import MPI
from time import process_time

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class WriteMesh(object):
    """
    Class for writing model outputs. The outputs are written as hdf5 files for each
    mesh partition. The outputs structure is composed of library

    .. note::

        The model outputs are all located in an output folder (`dir` key in the inputfile documentation) and consist of a time series file named `gospl.xdmf` and 2 other folders (`h5` and `xmf`).

    The `XDMF` file is the main entry point for visualising the output and should be sufficient for most users. This file can easely be opened within `Paraview <https://www.paraview.org/download/>`_.
    """

    def __init__(self):
        """
        Initialise model outputs paramters.
        """

        self.step = 0
        self.stratStep = 0
        self.file = "gospl"
        if MPIrank == 0:
            self._createOutputDir()

        # Sync the chosen output dir to all processors
        self.outputDir = MPIcomm.bcast(self.outputDir, root=0)

        # In case of a restarting simulation, get the corresponding time
        # values according to the restarting step
        if self.rStep > 0:
            self.step = self.rStep + 1
            if self.strat > 0:
                n = int(self.tout / self.strat)
                self.stratStep = self.rStep * n
            self.tNow = self.tStart + self.rStep * self.tout
            self.rStart = self.tStart + self.rStep * self.tout
            self.saveTime = self.tNow + self.tout
            self.saveStrat = self.tEnd + self.tout
            if self.strat > 0:
                self.saveStrat = self.tNow + self.strat
        else:
            self.rStart = self.tStart

        if self.strataFile is not None:
            self.stratStep = self.initLay

        # If series of paleomaps are used to force the model through time.
        if self.forceStep >= 0:

            res = 0.5
            nx = int(360.0 / res)
            ny = int(180.0 / res)
            lon = np.linspace(0.0, 360.0, nx)
            lat = np.linspace(0, 180, ny)

            lon, lat = np.meshgrid(lon, lat)
            latlonreg = np.dstack([lat.flatten(), lon.flatten()])[0]
            self.regshape = lon.shape
            outTree = spatial.cKDTree(self.lLatLon, leafsize=10)
            distances, self.forceIDs = outTree.query(latlonreg, k=3)
            self.forceWeights = 1.0 / distances ** 2
            self.forceTop = 1.0 / np.sum(self.forceWeights, axis=1)
            self.forceonIDs = np.where(distances[:, 0] == 0)[0]

            outTree = spatial.cKDTree(latlonreg, leafsize=10)
            distances, self.regIDs = outTree.query(self.lLatLon, k=3)
            self.regWeights = 1.0 / distances ** 2
            self.regTop = 1.0 / np.sum(self.regWeights, axis=1)
            self.regonIDs = np.where(distances[:, 0] == 0)[0]

            del distances, outTree, lon, lat, latlonreg
            gc.collect()

        # Petsc vectors
        self.upsG = self.hGlobal.duplicate()
        self.upsL = self.hLocal.duplicate()

        return

    def visModel(self):
        """
        Main function to write model outputs on disk.
        """

        # Output time step for first step
        if self.saveTime == self.tStart:
            self._outputMesh()
            self.saveTime += self.tout

        # Output time step after start time
        elif self.tNow >= self.saveTime:
            self._outputMesh()
            self.saveTime += self.tout

            # Forcing with backward model
            if self.forceStep >= 0 and self.newForcing:
                self._forcePaleo()
                self.steppaleo += 1

            if self.steppaleo == 1:
                self.steppaleo = 0
                self.newForcing = False
                self.forceStep += 1
            else:
                self.newForcing = True

        return

    def _createOutputDir(self):
        """
        Create a directory to store outputs. By default the folder will be called
        `output`. If a folder name is specified in the YAML input file, this name
        will be used.

        .. note::
            The input option `makedir` gives the ability to delete any existing output folder with the same name (if set to `False`) or to create a new folder with the given dir name plus a number at the end (*e.g.* `outputDir_XX` if set to `True` with `XX` the run number). It prevents overwriting on top of previous runs.

        """

        # Get output directory
        if self.outputDir is not None:
            self.outputDir = os.getcwd() + "/" + self.outputDir
        else:
            self.outputDir = os.getcwd() + "/output"

        if self.rStep == 0:
            if self.makedir:
                if os.path.exists(self.outputDir):
                    self.outputDir += "_" + str(
                        len(glob.glob(self.outputDir + str("*"))) - 1
                    )
            else:
                if os.path.exists(self.outputDir):
                    shutil.rmtree(self.outputDir, ignore_errors=True)

            os.makedirs(self.outputDir)
            os.makedirs(self.outputDir + "/h5")
            os.makedirs(self.outputDir + "/xmf")

        return

    def _outputStrat(self):
        """
        Saves mesh local stratigraphic information stored in the DMPlex
        as HDF5 file. The following variables will be recorded:

        - elevation at time of deposition, considered to be to the current elevation
          for the top stratigraphic layer `stratZ`.
        - thickness of each stratigrapic layer `stratH` accounting for both
          erosion & deposition events.
        - proportion of fine sediment `stratF` contains in each stratigraphic layer.
        - proportion of weathered sediment `stratW` contains in each stratigraphic layer.
        - porosity of coarse sediment `phiS` in each stratigraphic layer computed at
          center of each layer.
        - porosity of fine sediment `phiF` in each stratigraphic layer computed at
          center of each layer.
        - porosity of weathered sediment `phiW` in each stratigraphic layer computed at
          center of each layer.
        - proportion of carbonate sediment `stratC` contains in each stratigraphic layer
          if the carbonate module is turned on.
        - porosity of carbonate sediment `phiC` in each stratigraphic layer computed at
          center of each layer when the carbonate module is turned on.

        .. important::

            It is worth mentioning that the stratigraphic architecture is only outputed as
            HDF5 files and does not record the XMF and XDMF files. A set of post-processing
            scripts are then required to extract the informations and visualise the stratigraphic
            records of any specific simulations.
        """

        t = process_time()
        h5file = (
            self.outputDir
            + "/h5/stratal."
            + str(self.step)
            + ".p"
            + str(MPIrank)
            + ".h5"
        )
        with h5py.File(h5file, "w") as f:

            # Write stratal layers elevations per layers
            f.create_dataset(
                "stratZ",
                shape=(self.lpoints, self.stratStep),
                dtype="float32",
                compression="gzip",
            )
            f["stratZ"][:, : self.stratStep] = self.stratZ[:, : self.stratStep]

            # Write stratal layers thicknesses per layers
            f.create_dataset(
                "stratH",
                shape=(self.lpoints, self.stratStep),
                dtype="float64",
                compression="gzip",
            )
            f["stratH"][:, : self.stratStep] = self.stratH[:, : self.stratStep]

            # Write porosity values for coarse sediments
            f.create_dataset(
                "phiS",
                shape=(self.lpoints, self.stratStep),
                dtype="float64",
                compression="gzip",
            )
            f["phiS"][:, : self.stratStep] = self.phiS[:, : self.stratStep]

            if self.stratF is not None:
                # Write stratal layers fine percentage per layers
                f.create_dataset(
                    "stratF",
                    shape=(self.lpoints, self.stratStep),
                    dtype="float64",
                    compression="gzip",
                )
                f["stratF"][:, : self.stratStep] = self.stratF[:, : self.stratStep]

                # Write porosity values for fine sediments
                f.create_dataset(
                    "phiF",
                    shape=(self.lpoints, self.stratStep),
                    dtype="float64",
                    compression="gzip",
                )
                f["phiF"][:, : self.stratStep] = self.phiF[:, : self.stratStep]

            if self.stratW is not None:
                # Write stratal layers weathered percentage per layers
                f.create_dataset(
                    "stratW",
                    shape=(self.lpoints, self.stratStep),
                    dtype="float64",
                    compression="gzip",
                )
                f["stratW"][:, : self.stratStep] = self.stratW[:, : self.stratStep]

                # Write porosity values for weathered sediments
                f.create_dataset(
                    "phiW",
                    shape=(self.lpoints, self.stratStep),
                    dtype="float64",
                    compression="gzip",
                )
                f["phiW"][:, : self.stratStep] = self.phiW[:, : self.stratStep]

            # Write stratal layers carbonate percentage per layers
            if self.carbOn:
                f.create_dataset(
                    "stratC",
                    shape=(self.lpoints, self.stratStep),
                    dtype="float64",
                    compression="gzip",
                )
                f["stratC"][:, : self.stratStep] = self.stratC[:, : self.stratStep]

                # Write porosity values for carbonate sediments
                f.create_dataset(
                    "phiC",
                    shape=(self.lpoints, self.stratStep),
                    dtype="float64",
                    compression="gzip",
                )
                f["phiC"][:, : self.stratStep] = self.phiC[:, : self.stratStep]

        MPIcomm.Barrier()

        if MPIrank == 0 and self.verbose:
            print(
                "Creating stratal outputfile \
                  (%0.02f seconds)"
                % (process_time() - t),
                flush=True,
            )

        return

    def _outputMesh(self):
        """
        Saves mesh local information stored in the DMPlex to HDF5 file. If the file already
        exists, it will be overwritten. Mesh characteristics are recorded for each partition.
        The following variables will be available:

        - surface elevation `elev`.
        - cumulative erosion & deposition values `erodep`.
        - flow accumulation `flowAcc` before pit filling.
        - flow accumulation `fillAcc` for depressionless surface.
        - river sediment load `sedLoad`.
        - fine sediment load `sedLoadf` when dual lithologies are accounted for.
        - carbonate sediment load `sedLoadc` when carbonate module is turned on.
        - uplift subsidence values if vertical tectonic forcing is considered `uplift`.
        - horizontal displacement values when considered `hdisp`.
        - precipitation maps based on forcing conditions `rain`.

        """

        t = process_time()
        if self.step == 0:
            topology = self.outputDir + "/h5/topology.p" + str(MPIrank) + ".h5"
            with h5py.File(topology, "w") as f:
                f.create_dataset(
                    "coords",
                    shape=(self.lpoints, 3),
                    dtype="float32",
                    compression="gzip",
                )
                f["coords"][:, :] = self.lcoords
                f.create_dataset(
                    "cells",
                    shape=(len(self.lcells[:, 0]), 3),
                    dtype="int32",
                    compression="gzip",
                )
                f["cells"][:, :] = self.lcells + 1
            self.elems = MPIcomm.gather(len(self.lcells[:, 0]), root=0)
            self.nodes = MPIcomm.gather(self.lpoints, root=0)

        h5file = (
            self.outputDir
            + "/h5/"
            + self.file
            + "."
            + str(self.step)
            + ".p"
            + str(MPIrank)
            + ".h5"
        )
        with h5py.File(h5file, "w") as f:
            f.create_dataset(
                "elev",
                shape=(self.lpoints, 1),
                dtype="float32",
                compression="gzip",
            )
            f["elev"][:, 0] = self.hLocal.getArray()
            f.create_dataset(
                "erodep",
                shape=(self.lpoints, 1),
                dtype="float32",
                compression="gzip",
            )
            f["erodep"][:, 0] = self.cumEDLocal.getArray()
            f.create_dataset(
                "flowAcc",
                shape=(self.lpoints, 1),
                dtype="float32",
                compression="gzip",
            )
            data = self.FAL.getArray().copy()
            data[data <= 0.0] = 1.0
            f["flowAcc"][:, 0] = data
            f.create_dataset(
                "fillAcc",
                shape=(self.lpoints, 1),
                dtype="float32",
                compression="gzip",
            )
            data = self.fillFAL.getArray().copy()
            data[data <= 0.0] = 1.0
            f["fillAcc"][:, 0] = data
            f.create_dataset(
                "sedLoad",
                shape=(self.lpoints, 1),
                dtype="float32",
                compression="gzip",
            )
            f["sedLoad"][:, 0] = self.vSedLocal.getArray().copy()
            if self.stratNb > 0:
                if self.stratF is not None:
                    f.create_dataset(
                        "sedLoadf",
                        shape=(self.lpoints, 1),
                        dtype="float32",
                        compression="gzip",
                    )
                    f["sedLoadf"][:, 0] = self.vSedfLocal.getArray().copy()
                if self.stratW is not None:
                    f.create_dataset(
                        "sedLoadw",
                        shape=(self.lpoints, 1),
                        dtype="float32",
                        compression="gzip",
                    )
                    f["sedLoadw"][:, 0] = self.vSedwLocal.getArray().copy()
                if self.carbOn:
                    f.create_dataset(
                        "sedLoadc",
                        shape=(self.lpoints, 1),
                        dtype="float32",
                        compression="gzip",
                    )
                    f["sedLoadc"][:, 0] = self.vSedcLocal.getArray().copy()
            if self.uplift is not None:
                f.create_dataset(
                    "uplift",
                    shape=(self.lpoints, 1),
                    dtype="float32",
                    compression="gzip",
                )
                f["uplift"][:, 0] = self.uplift
            if self.hdisp is not None:
                f.create_dataset(
                    "hdisp",
                    shape=(self.lpoints, 3),
                    dtype="float32",
                    compression="gzip",
                )
                f["hdisp"][:, :] = self.hdisp
            if self.rainVal is not None:
                f.create_dataset(
                    "rain",
                    shape=(self.lpoints, 1),
                    dtype="float32",
                    compression="gzip",
                )
                f["rain"][:, 0] = self.rainVal

            del data

        if self.stratNb > 0 and self.stratStep > 0:
            self._outputStrat()

        if MPIrank == 0:
            self._save_DMPlex_XMF()
            self._save_XDMF()

        MPIcomm.Barrier()
        if MPIrank == 0 and self.verbose:
            print(
                "Creating outputfile (%0.02f seconds)" % (process_time() - t),
                flush=True,
            )

        if MPIrank == 0:
            print("+++ Output Simulation Time: %0.02f years" % (self.tNow), flush=True)

        self.step += 1
        gc.collect()

        return

    def readData(self):
        """
        When a simulation restarts, variables from previous HDF5 output files are read and assigned
        to the restarting run.

        The following variables are used:

        - surface elevation `elev`.
        - cumulative erosion & deposition values `erodep`.
        - flow accumulation `flowAcc` before pit filling.
        - flow accumulation `fillAcc` for depressionless surface.
        - river sediment load `sedLoad`.
        - fine sediment load `sedLoadf` when dual lithologies are accounted for.
        - weathered sediment load `sedLoadw` when dual lithologies are accounted for.
        - carbonate sediment load `sedLoadc` when carbonate module is turned on.

        """

        h5file = (
            self.outputDir
            + "/h5/"
            + self.file
            + "."
            + str(self.step - 1)
            + ".p"
            + str(MPIrank)
            + ".h5"
        )
        if os.path.exists(h5file):
            hf = h5py.File(h5file, "r")
        else:
            raise ValueError("Restart file is missing...")
        self.hLocal.setArray(np.array(hf["/elev"])[:, 0])
        self.dm.localToGlobal(self.hLocal, self.hGlobal)
        self.cumEDLocal.setArray(np.array(hf["/erodep"])[:, 0])
        self.dm.localToGlobal(self.cumEDLocal, self.cumED)
        self.vSedLocal.setArray(np.array(hf["/sedLoad"])[:, 0])
        self.dm.localToGlobal(self.vSedLocal, self.vSed)
        self.FAL.setArray(np.array(hf["/flowAcc"])[:, 0])
        self.dm.localToGlobal(self.FAL, self.FAG)
        self.fillFAL.setArray(np.array(hf["/fillAcc"])[:, 0])
        self.dm.localToGlobal(self.fillFAL, self.fillFAG)
        if self.stratNb > 0:
            if self.stratF is not None:
                self.vSedfLocal.setArray(np.array(hf["/sedLoadf"])[:, 0])
                self.dm.localToGlobal(self.vSedfLocal, self.vSedf)
            if self.stratW is not None:
                self.vSedwLocal.setArray(np.array(hf["/sedLoadw"])[:, 0])
                self.dm.localToGlobal(self.vSedwLocal, self.vSedw)
            if self.carbOn:
                self.vSedcLocal.setArray(np.array(hf["/sedLoadc"])[:, 0])
                self.dm.localToGlobal(self.vSedcLocal, self.vSedc)
        self.elems = MPIcomm.gather(len(self.lcells[:, 0]), root=0)
        self.nodes = MPIcomm.gather(len(self.lcoords[:, 0]), root=0)

        hf.close()

        if self.stratNb > 0 and self.stratStep > 0:
            h5file = (
                self.outputDir
                + "/h5/stratal."
                + str(self.step - 1)
                + ".p"
                + str(MPIrank)
                + ".h5"
            )
            if os.path.exists(h5file):
                hf = h5py.File(h5file, "r")
            else:
                raise ValueError("Restart file is missing...")

            self.stratZ.fill(0.0)
            self.stratZ[:, : self.stratStep] = np.array(hf["/stratZ"])
            self.stratH.fill(0.0)
            self.stratH[:, : self.stratStep] = np.array(hf["/stratH"])
            self.phiS.fill(0.0)
            self.phiS[:, : self.stratStep] = np.array(hf["/phiS"])

            if self.stratF is not None:
                self.stratF.fill(0.0)
                self.stratF[:, : self.stratStep] = np.array(hf["/stratF"])
                self.phiF.fill(0.0)
                self.phiF[:, : self.stratStep] = np.array(hf["/phiF"])
            if self.stratW is not None:
                self.stratW.fill(0.0)
                self.stratW[:, : self.stratStep] = np.array(hf["/stratW"])
                self.phiW.fill(0.0)
                self.phiW[:, : self.stratStep] = np.array(hf["/phiW"])
            if self.carbOn:
                self.stratC.fill(0.0)
                self.stratC[:, : self.stratStep] = np.array(hf["/stratC"])
                self.phiC.fill(0.0)
                self.phiC[:, : self.stratStep] = np.array(hf["/phiC"])

            hf.close()

        return

    def _forcePaleo(self):
        """
        Forcing the model based on backward simulation results. This function computes
        difference between backward model elevation and current one.

        .. note::

            The backward elevation is often computed based on paleo-elevation maps at specific
            time interval. During the pre-processing phase a series of backward elevations are
            then created based on these paleo-elevation maps. These backward elevations are then
            used to force the forward model over time. This is done by computing the differences
            at any given time step between the obtained forward elevations and the predicted backward
            ones. These differences are then used to define an uplift map that matches the differences
            modulated by an imposed factor (`alpha`). Once the uplift map has been computed the model
            is then rerun from the previous time step to the current one.

        """

        # Store the vertical displacement in the uplift variable
        self._upliftBF(self.alpha[self.forceStep], self.stepb[self.forceStep])

        # Update model parameters to the former tectonic time step
        self.tNow -= self.tout
        self.step -= 1
        n = int(self.tout / self.tecStep)
        self.tecNb -= n
        if self.strat > 0:
            n = int(self.tout / self.strat)
            self.stratStep = (self.step - 1) * n
        self.saveTime = self.tNow + self.tout

        if self.saveStrat <= self.tEnd + self.strat:
            self.saveStrat = self.tNow + self.strat

        # Getting PETSc vectors values
        if self.step == 0:
            loadData = np.load(self.meshFile)
            gZ = loadData["z"]
            self.hLocal.setArray(gZ[self.locIDs])
            self.dm.localToGlobal(self.hLocal, self.hGlobal)
            self.vSed.set(0.0)
            self.vSedLocal.set(0.0)
            if self.stratNb > 0:
                if self.stratF is not None:
                    self.vSedf.set(0.0)
                    self.vSedfLocal.set(0.0)
                if self.stratW is not None:
                    self.vSedw.set(0.0)
                    self.vSedwLocal.set(0.0)
                if self.carbOn:
                    self.vSedc.set(0.0)
                    self.vSedcLocal.set(0.0)
            self.cumED.set(0.0)
            self.cumEDLocal.set(0.0)
        else:
            self.readData()

        return

    def _upliftBF(self, alpha, step):
        """
        Reads locally the backward model elevation at any given step.
        Compute elevation difference between backward and forward models.

        :arg alpha: fitting coefficient for backward uplift differences
        :arg step: backward model time step to use
        """

        h5file = self.forceDir + "/h5/gospl." + str(step) + ".p" + str(MPIrank) + ".h5"
        if os.path.exists(h5file):
            hf = h5py.File(h5file, "r")
        else:
            print("Backward file: {} is missing!".format(h5file), flush=True)
            raise ValueError("Backward file is missing...")

        # Get the difference between backward and forward model on TIN
        ldiff = np.array(hf["/elev"])[:, 0] - self.hLocal.getArray()

        # erosional features using a low-pass filter
        gdiff = ldiff[self.lgIDs]
        if MPIsize > 1:
            temp = np.full(self.shadowgNb, -1.0e8, dtype=np.float64)
            temp[self.gshadinIDs] = gdiff[self.gshadowIDs]
            temp[self.gshadoutIDs] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            gdiff[self.shadowAlls] = temp
        ldiff = gdiff[self.locIDs]

        # self.uplift = ldiff
        # self.uplift *= alpha / self.tecStep
        # hf.close()
        # return

        # Map it on a regularly spaced mesh (lon/lat)
        diffreg = (
            np.sum(self.forceWeights * ldiff[self.forceIDs], axis=1) * self.forceTop
        )
        if len(self.forceonIDs) > 0:
            diffreg[self.forceonIDs] = ldiff[self.forceIDs[self.forceonIDs, 0]]

        # Apply a low-pass filter to the regular array
        # median = False
        kernel_size = 3
        # if median:
        #     filter = ndimage.median_filter(
        #         diffreg.reshape(self.regshape), size=kernel_size, mode="wrap"
        #     )
        # else:
        filter = ndimage.gaussian_filter(
            diffreg.reshape(self.regshape), sigma=kernel_size, mode="wrap"
        )
        filter = filter.flatten()

        # Map it back to the spherical mesh
        lfilter = np.sum(self.regWeights * filter[self.regIDs], axis=1) * self.regTop
        if len(self.regonIDs) > 0:
            lfilter[self.regonIDs] = filter[self.regIDs[self.regonIDs, 0]]

        # From local to global
        gfilter = lfilter[self.lgIDs]
        if MPIsize > 1:
            temp = np.full(self.shadowgNb, -1.0e8, dtype=np.float64)
            temp[self.gshadinIDs] = gfilter[self.gshadowIDs]
            temp[self.gshadoutIDs] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            gfilter[self.shadowAlls] = temp

        # Specify new uplift value matching backward elevation removing the
        # erosional features using a low-pass filter
        self.uplift = gfilter[self.locIDs]
        self.uplift[self.seaID] = ldiff[self.seaID]

        self.uplift *= alpha / self.tecStep
        hf.close()

        del diffreg, filter, lfilter, gfilter, ldiff, hf
        if MPIsize > 1:
            del temp
        gc.collect()

        return

    def _save_DMPlex_XMF(self):
        """
        Saves mesh local information stored in the HDF5 to XmF file. The XMF files are XML schema
        explaining how to read `gospl` data files.

        The XmF file is written by a single processor (rank 0) and contains each partition HDF5 files
        in blocks. The variables described for the HDF5 file (function `_outputMesh` above) are all
        accessible from this file.

        """

        xmf_file = self.outputDir + "/xmf/" + self.file + str(self.step) + ".xmf"

        f = open(xmf_file, "w")

        # Header for xml file
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        f.write(" <Domain>\n")
        f.write('    <Grid GridType="Collection" CollectionType="Spatial">\n')
        f.write('      <Time Type="Single" Value="%0.02f"/>\n' % self.saveTime)

        for p in range(MPIsize):
            pfile = (
                "h5/" + str(self.file) + "." + str(self.step) + ".p" + str(p) + ".h5"
            )
            tfile = "h5/topology.p" + str(p) + ".h5"
            f.write('      <Grid Name="Block.%s">\n' % (str(p)))
            f.write(
                '         <Topology Type="Triangle" NumberOfElements="%d" BaseOffset="1">\n'
                % self.elems[p]
            )
            f.write('          <DataItem Format="HDF" DataType="Int" ')
            f.write('Dimensions="%d 3">%s:/cells</DataItem>\n' % (self.elems[p], tfile))
            f.write("         </Topology>\n")

            f.write('         <Geometry Type="XYZ">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write(
                'Dimensions="%d 3">%s:/coords</DataItem>\n' % (self.nodes[p], tfile)
            )
            f.write("         </Geometry>\n")

            f.write('         <Attribute Type="Scalar" Center="Node" Name="Z">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write('Dimensions="%d 1">%s:/elev</DataItem>\n' % (self.nodes[p], pfile))
            f.write("         </Attribute>\n")

            f.write('         <Attribute Type="Scalar" Center="Node" Name="ED">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write(
                'Dimensions="%d 1">%s:/erodep</DataItem>\n' % (self.nodes[p], pfile)
            )
            f.write("         </Attribute>\n")

            f.write('         <Attribute Type="Scalar" Center="Node" Name="FA">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write(
                'Dimensions="%d 1">%s:/flowAcc</DataItem>\n' % (self.nodes[p], pfile)
            )
            f.write("         </Attribute>\n")

            f.write('         <Attribute Type="Scalar" Center="Node" Name="fillFA">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write(
                'Dimensions="%d 1">%s:/fillAcc</DataItem>\n' % (self.nodes[p], pfile)
            )
            f.write("         </Attribute>\n")

            f.write('         <Attribute Type="Scalar" Center="Node" Name="SL">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write(
                'Dimensions="%d 1">%s:/sedLoad</DataItem>\n' % (self.nodes[p], pfile)
            )
            f.write("         </Attribute>\n")
            if self.stratNb > 0:
                if self.stratF is not None:
                    f.write(
                        '         <Attribute Type="Scalar" Center="Node" Name="SLf">\n'
                    )
                    f.write(
                        '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
                    )
                    f.write(
                        'Dimensions="%d 1">%s:/sedLoadf</DataItem>\n'
                        % (self.nodes[p], pfile)
                    )
                    f.write("         </Attribute>\n")
                if self.stratW is not None:
                    f.write(
                        '         <Attribute Type="Scalar" Center="Node" Name="SLw">\n'
                    )
                    f.write(
                        '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
                    )
                    f.write(
                        'Dimensions="%d 1">%s:/sedLoadw</DataItem>\n'
                        % (self.nodes[p], pfile)
                    )
                    f.write("         </Attribute>\n")
                if self.carbOn:
                    f.write(
                        '         <Attribute Type="Scalar" Center="Node" Name="SLc">\n'
                    )
                    f.write(
                        '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
                    )
                    f.write(
                        'Dimensions="%d 1">%s:/sedLoadc</DataItem>\n'
                        % (self.nodes[p], pfile)
                    )
                    f.write("         </Attribute>\n")

            if self.hdisp is not None:
                f.write(
                    '         <Attribute Type="Vector" Center="Node" Name="hTec">\n'
                )
                f.write(
                    '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
                )
                f.write(
                    'Dimensions="%d 3">%s:/hdisp</DataItem>\n' % (self.nodes[p], pfile)
                )
                f.write("         </Attribute>\n")
            if self.uplift is not None:
                f.write(
                    '         <Attribute Type="Scalar" Center="Node" Name="vTec">\n'
                )
                f.write(
                    '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
                )
                f.write(
                    'Dimensions="%d 1">%s:/uplift</DataItem>\n' % (self.nodes[p], pfile)
                )
                f.write("         </Attribute>\n")
            if self.rainVal is not None:
                f.write(
                    '         <Attribute Type="Scalar" Center="Node" Name="Rain">\n'
                )
                f.write(
                    '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
                )
                f.write(
                    'Dimensions="%d 1">%s:/rain</DataItem>\n' % (self.nodes[p], pfile)
                )
                f.write("         </Attribute>\n")

            f.write('         <Attribute Type="Scalar" Center="Node" Name="sea">\n')
            f.write(
                '          <DataItem ItemType="Function" Function="$0 * 0.00000000000000001 + %f" Dimensions="%d 1">\n'
                % (self.sealevel, self.nodes[p])
            )
            f.write(
                '           <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write(
                'Dimensions="%d 1">%s:/erodep</DataItem>\n' % (self.nodes[p], pfile)
            )
            f.write("          </DataItem>\n")
            f.write("         </Attribute>\n")

            f.write("      </Grid>\n")

        f.write("    </Grid>\n")
        f.write(" </Domain>\n")
        f.write("</Xdmf>\n")
        f.close()

        return

    def _save_XDMF(self):
        """
        This function writes the XDmF file which is calling the XmF files above. The XDmF file
        represents the *time series* of the model outputs and can be directly loaded and visualised
        with `Paraview <https://www.paraview.org/download/>`_.

        .. note::

            For a brief overview of the approach used to record `gospl` outputs, user can read
            this `visit documentation <https://www.visitusers.org/index.php?title=Using_XDMF_to_read_HDF5>`_
        """

        xdmf_file = self.outputDir + "/" + self.file + ".xdmf"
        f = open(xdmf_file, "w")

        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        f.write(" <Domain>\n")
        f.write('    <Grid GridType="Collection" CollectionType="Temporal">\n')

        for s in range(self.step + 1):
            xmf_file = "xmf/" + self.file + str(s) + ".xmf"
            f.write(
                '      <xi:include href="%s" xpointer="xpointer(//Xdmf/Domain/Grid)"/>\n'
                % xmf_file
            )

        f.write("    </Grid>\n")
        f.write(" </Domain>\n")
        f.write("</Xdmf>\n")
        f.close()

        return
