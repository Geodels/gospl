import os
import gc
import sys
import h5py
import glob
import shutil
import petsc4py
import numpy as np

from mpi4py import MPI
from time import process_time

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class WriteMesh(object):
    """
    Class for writing model outputs. The outputs are written as hdf5 files for each mesh partition.

    .. note::

        The model outputs are all located in an output folder (`dir` key in the inputfile documentation) and consist of a time series file named `gospl.xdmf` and 2 other folders (`h5` and `xmf`).

    The `XDMF` file is the main entry point for visualising the output and should be sufficient for most users. This file can easely be opened within `Paraview <https://www.paraview.org/download/>`_.
    """

    def __init__(self):
        """
        Initialise model outputs parameters.
        """

        self.step = 0
        self.stratStep = 1
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
                self.stratStep = self.rStep * n + 1
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

        return

    def _createOutputDir(self):
        """
        Create a directory to store outputs. By default the folder will be called `output`. If a folder name is specified in the YAML input file, this name will be used.

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

        shutil.copy(self.finput, self.outputDir)

        return

    def _outputStrat(self):
        """
        Saves mesh local stratigraphic information stored in the DMPlex as HDF5 file. The following variables will be recorded:

        - elevation at time of deposition, considered to be to the current elevation for the top stratigraphic layer `stratZ`.
        - thickness of each stratigrapic layer `stratH` accounting for both erosion & deposition events.
        - porosity of sediment `phiS` in each stratigraphic layer computed at center of each layer.

        .. important::

            It is worth mentioning that the stratigraphic architecture is only outputed as HDF5 files and does not record the XMF and XDMF files. A set of post-processing scripts are then required to extract the informations and visualise the stratigraphic records of any specific simulations.
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
                shape=(self.lpoints, self.stratStep + 1),
                dtype="float32",
                compression="gzip",
            )
            f["stratZ"][:, : self.stratStep + 1] = self.stratZ[:, : self.stratStep + 1]

            # Write stratal layers thicknesses per layers
            f.create_dataset(
                "stratH",
                shape=(self.lpoints, self.stratStep + 1),
                dtype="float64",
                compression="gzip",
            )
            f["stratH"][:, : self.stratStep + 1] = self.stratH[:, : self.stratStep + 1]

            # Write porosity values for coarse sediments
            f.create_dataset(
                "phiS",
                shape=(self.lpoints, self.stratStep + 1),
                dtype="float64",
                compression="gzip",
            )
            f["phiS"][:, : self.stratStep + 1] = self.phiS[:, : self.stratStep + 1]

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
        Saves mesh local information stored in the DMPlex to HDF5 file. If the file already exists, it will be overwritten. Mesh characteristics are recorded for each partition. The following variables will be available:

        - surface elevation `elev`.
        - cumulative erosion & deposition values `erodep`.
        - erosion & deposition rate values `EDrate` for the considered time step.
        - flow accumulation `FA`.
        - flow accumulation `fillFA` considering pit filling.
        - river sediment load `sedLoad`.
        - uplift subsidence values if vertical tectonic forcing is considered `uplift`.
        - flexural isostasy rebound `flexIso` if flexure is considered.
        - precipitation maps based on forcing conditions `rain` (could also correspond to the orographic rain if the functionality is turned on).

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
                "EDrate",
                shape=(self.lpoints, 1),
                dtype="float32",
                compression="gzip",
            )
            data = self.EbLocal.getArray().copy()
            f["EDrate"][:, 0] = data
            if not self.fast:
                f.create_dataset(
                    "waterFill",
                    shape=(self.lpoints, 1),
                    dtype="float32",
                    compression="gzip",
                )
                f["waterFill"][:, 0] = self.waterFilled
            f.create_dataset(
                "fillFA",
                shape=(self.lpoints, 1),
                dtype="float32",
                compression="gzip",
            )
            data = self.fillFAL.getArray().copy()
            data[data <= 1.0e-8] = 1.0e-8
            if not self.fast:
                data[self.seaID] = 1.0
            f["fillFA"][:, 0] = data
            f.create_dataset(
                "FA",
                shape=(self.lpoints, 1),
                dtype="float32",
                compression="gzip",
            )
            data = self.FAL.getArray().copy()
            data[data <= 1.0e-8] = 1.0e-8
            if not self.fast:
                data[self.seaID] = 1.0
            f["FA"][:, 0] = data

            if self.iceOn:
                f.create_dataset(
                    "iceFA",
                    shape=(self.lpoints, 1),
                    dtype="float32",
                    compression="gzip",
                )
                data = self.iceFAL.getArray().copy()
                data[data <= 1.0e-8] = 1.0e-8
                if not self.fast:
                    data[self.seaID] = 1.0
                f["iceFA"][:, 0] = data
            if self.flexOn:
                f.create_dataset(
                    "flexIso",
                    shape=(self.lpoints, 1),
                    dtype="float32",
                    compression="gzip",
                )
                f["flexIso"][:, 0] = self.localFlex
            f.create_dataset(
                "sedLoad",
                shape=(self.lpoints, 1),
                dtype="float32",
                compression="gzip",
            )
            data = self.vSedLocal.getArray().copy()
            data[data <= 1.0e-8] = 1.0e-8
            f["sedLoad"][:, 0] = data
            if self.upsub is not None:
                f.create_dataset(
                    "uplift",
                    shape=(self.lpoints, 1),
                    dtype="float32",
                    compression="gzip",
                )
                f["uplift"][:, 0] = self.upsub
            if self.rainVal is not None:
                f.create_dataset(
                    "rain",
                    shape=(self.lpoints, 1),
                    dtype="float32",
                    compression="gzip",
                )
                f["rain"][:, 0] = self.rainVal
            if self.memclear:
                del data
                gc.collect()

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

        return

    def readData(self):
        """
        When a simulation restarts, variables from previous HDF5 output files are read and assigned to the restarting run.

        The following variables are used:

        - surface elevation `elev`.
        - cumulative erosion & deposition values `erodep`.
        - erosion & deposition values `EDrate` for the considered time step.
        - flow accumulation `fillFA` considering pit filling.
        - river sediment load `sedLoad`.
        - flexural isostasy induced tectonics `flexIso`.

        .. note::

            If stratigraphy is turned on, the function also reads underlying stratigraphic information.

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
        self.EbLocal.setArray(np.array(hf["/EDrate"])[:, 0])
        self.dm.localToGlobal(self.EbLocal, self.Eb)
        self.FAL.setArray(np.array(hf["/FA"])[:, 0])
        self.fillFAL.setArray(np.array(hf["/fillFA"])[:, 0])
        self.dm.localToGlobal(self.FAL, self.FAG)
        self.elems = MPIcomm.gather(len(self.lcells[:, 0]), root=0)
        self.nodes = MPIcomm.gather(len(self.lcoords[:, 0]), root=0)

        if self.flexOn:
            self.localFlex = np.array(hf["/flexIso"])[:, 0]
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

            hf.close()

        return

    def _save_DMPlex_XMF(self):
        """
        Saves mesh local information stored in the HDF5 to XmF file. The XMF files are XML schema explaining how to read `gospl` data files.

        The XmF file is written by a single processor (rank 0) and contains each partition HDF5 files in blocks. The variables described for the HDF5 file (function `_outputMesh` above) are all accessible from this file.

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
            if not self.fast:
                f.write(
                    '         <Attribute Type="Scalar" Center="Node" Name="waterFill">\n'
                )
                f.write(
                    '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
                )
                f.write(
                    'Dimensions="%d 1">%s:/waterFill</DataItem>\n'
                    % (self.nodes[p], pfile)
                )
                f.write("         </Attribute>\n")

            f.write('         <Attribute Type="Scalar" Center="Node" Name="EDrate">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write(
                'Dimensions="%d 1">%s:/EDrate</DataItem>\n' % (self.nodes[p], pfile)
            )
            f.write("         </Attribute>\n")

            f.write('         <Attribute Type="Scalar" Center="Node" Name="FA">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write(
                'Dimensions="%d 1">%s:/FA</DataItem>\n' % (self.nodes[p], pfile)
            )
            f.write("         </Attribute>\n")

            f.write('         <Attribute Type="Scalar" Center="Node" Name="fillFA">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write(
                'Dimensions="%d 1">%s:/fillFA</DataItem>\n' % (self.nodes[p], pfile)
            )
            f.write("         </Attribute>\n")

            if self.iceOn:
                f.write('         <Attribute Type="Scalar" Center="Node" Name="iceFA">\n')
                f.write(
                    '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
                )
                f.write(
                    'Dimensions="%d 1">%s:/iceFA</DataItem>\n' % (self.nodes[p], pfile)
                )
                f.write("         </Attribute>\n")

            if self.flexOn:
                f.write('         <Attribute Type="Scalar" Center="Node" Name="flexIso">\n')
                f.write(
                    '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
                )
                f.write(
                    'Dimensions="%d 1">%s:/flexIso</DataItem>\n' % (self.nodes[p], pfile)
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

            if self.upsub is not None:
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
        This function writes the XDmF file which is calling the XmF files above. The XDmF file represents the *time series* of the model outputs and can be directly loaded and visualised with `Paraview <https://www.paraview.org/download/>`_.

        .. note::

            For a brief overview of the approach used to record `gospl` outputs, user can read this `visit documentation <https://www.visitusers.org/index.php?title=Using_XDMF_to_read_HDF5>`_
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
