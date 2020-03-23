import os
import gc
import sys
import meshplex
import petsc4py
import numpy as np
import pandas as pd

from time import clock
from mpi4py import MPI
from scipy import spatial
from petsc4py import PETSc

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import defineTIN
    from gospl._fortran import ngbGlob

petsc4py.init(sys.argv)
MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIsize = PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class UnstMesh(object):
    """
    Creating a distributed DMPlex and global vector based on imported triangulated mesh.
    """

    def __init__(self):
        """
        Build the DMPlex.
        """

        self.hdisp = None
        self.uplift = None
        self.rainVal = None

        self._buildMesh()

        return

    def _meshfrom_cell_list(self, dim, cells, coords):
        """
        Create a DMPlex from a list of cells and coords.

        :arg dim: The topological dimension of the mesh
        :arg cells: The vertices of each cell
        :arg coords: The coordinates of each vertex
        :arg comm: communicator to build the mesh on
        """

        if MPIrank == 0:
            cells = np.asarray(cells, dtype=np.int32)
            coords = np.asarray(coords, dtype=np.double)
            MPIcomm.bcast(cells.shape, root=0)
            MPIcomm.bcast(coords.shape, root=0)
            # Provide the actual data on rank 0.
            self.dm = PETSc.DMPlex().createFromCellList(
                dim, cells, coords, comm=PETSc.COMM_WORLD
            )
            del cells, coords
        else:
            cell_shape = list(MPIcomm.bcast(None, root=0))
            coord_shape = list(MPIcomm.bcast(None, root=0))
            cell_shape[0] = 0
            coord_shape[0] = 0
            self.dm = PETSc.DMPlex().createFromCellList(
                dim,
                np.zeros(cell_shape, dtype=np.int32),
                np.zeros(coord_shape, dtype=np.double),
                comm=PETSc.COMM_WORLD,
            )
        return

    def _meshStructure(self):
        """
        Define the mesh structure and the associated voronoi
        """

        # Create mesh structure with meshplex
        t0 = clock()
        Tmesh = meshplex.MeshTri(self.lcoords, self.lcells)
        self.area = np.abs(Tmesh.control_volumes)
        self.area[np.isnan(self.area)] = 1.0

        # Voronoi and simplices declaration
        Tmesh.create_edges()
        cc = Tmesh.cell_circumcenters
        edges_nodes = Tmesh.edges["nodes"]
        cells_nodes = Tmesh.cells["nodes"]
        cells_edges = Tmesh.cells["edges"]

        # Finite volume discretisation
        self.FVmesh_ngbID = defineTIN(
            self.lcoords, cells_nodes, cells_edges, edges_nodes, self.area, cc.T
        )
        del Tmesh, edges_nodes, cells_nodes, cells_edges, cc
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print("FV discretisation (%0.02f seconds)" % (clock() - t0), flush=True)

        return

    def reInitialiseModel(self):
        """
        Reinitialise gospl model for paleo-fitting experiments.
        """

        t0step = clock()

        # Restart time
        self.tNow = self.tStart
        self.step = 0
        self.stratStep = 0
        self.rStart = self.tStart
        self.saveTime = self.tNow
        if self.strat > 0:
            self.saveStrat = self.tNow + self.strat
        else:
            self.saveStrat = self.tEnd + self.tout

        # Forcing functions
        self.rainNb = -1
        self.tecNb = -1
        self.flexNb = -1

        # Getting PETSc vectors values
        loadData = np.load(self.meshFile)
        gZ = loadData["z"]
        self.hLocal.setArray(gZ[self.glIDs])
        self.dm.localToGlobal(self.hLocal, self.hGlobal)
        self.vSed.set(0.0)
        self.vSedLocal.set(0.0)
        self.cumED.set(0.0)
        self.cumEDLocal.set(0.0)

        # Update external forces
        self.applyForces()
        self.applyTectonics()

        del gZ, loadData
        gc.collect()

        if MPIrank == 0:
            print(
                "--- Reinitialise Phase \
                  (%0.02f seconds)\n+++"
                % (clock() - t0step),
                flush=True,
            )

        return

    def _buildMesh(self):
        """
        Construct a PETSC mesh and distribute it amongst the different
        processors.
        """

        # Read mesh attributes from file
        t0 = clock()
        loadData = np.load(self.meshFile)
        self.gCoords = loadData["v"]
        self.gpoints = len(self.gCoords)
        gZ = loadData["z"]
        ngbGlob(self.gpoints, loadData["n"])
        if MPIrank == 0 and self.verbose:
            print(
                "Reading mesh information (%0.02f seconds)" % (clock() - t0), flush=True
            )

        # Create DMPlex
        self._meshfrom_cell_list(2, loadData["c"], self.gCoords)
        del loadData
        gc.collect()
        if MPIrank == 0 and self.verbose:
            print("Create DMPlex (%0.02f seconds)" % (clock() - t0), flush=True)

        # Define one DoF on the nodes
        t0 = clock()
        self.dm.setNumFields(1)
        origSect = self.dm.createSection(1, [1, 0, 0])
        origSect.setFieldName(0, "points")
        origSect.setUp()
        self.dm.setDefaultSection(origSect)
        origVec = self.dm.createGlobalVector()
        if MPIrank == 0 and self.verbose:
            print(
                "Define one DoF on the nodes (%0.02f seconds)" % (clock() - t0),
                flush=True,
            )

        # Distribute to other processors if any
        t0 = clock()
        if MPIsize > 1:
            partitioner = self.dm.getPartitioner()
            partitioner.setType(partitioner.Type.PARMETIS)
            partitioner.setFromOptions()
            sf = self.dm.distribute(overlap=1)
            newSect, newVec = self.dm.distributeField(sf, origSect, origVec)
            self.dm.setDefaultSection(newSect)
            newSect.destroy()
            newVec.destroy()
            sf.destroy()
        MPIcomm.Barrier()
        origVec.destroy()
        origSect.destroy()
        if MPIrank == 0 and self.verbose:
            print("Distribute DMPlex (%0.02f seconds)" % (clock() - t0), flush=True)

        # Define local vertex & cells
        t0 = clock()
        self.lcoords = self.dm.getCoordinatesLocal().array.reshape(-1, 3)
        self.npoints = self.lcoords.shape[0]
        cStart, cEnd = self.dm.getHeightStratum(0)
        self.lcells = np.zeros((cEnd - cStart, 3), dtype=PETSc.IntType)
        point_closure = None
        for c in range(cStart, cEnd):
            point_closure = self.dm.getTransitiveClosure(c)[0]
            self.lcells[c, :] = point_closure[-3:] - cEnd
        if point_closure is not None:
            del point_closure
        gc.collect()
        if MPIrank == 0 and self.verbose:
            print("Mesh coords/cells (%0.02f seconds)" % (clock() - t0), flush=True)

        # Define local vertex & cells
        t = clock()
        cStart, cEnd = self.dm.getHeightStratum(0)
        # Dealing with triangular cells only
        self.lcells = np.zeros((cEnd - cStart, 3), dtype=PETSc.IntType)
        for c in range(cStart, cEnd):
            point_closure = self.dm.getTransitiveClosure(c)[0]
            self.lcells[c, :] = point_closure[-3:] - cEnd
        del point_closure
        if MPIrank == 0 and self.verbose:
            print("Defining local DMPlex (%0.02f seconds)" % (clock() - t), flush=True)

        # From global values to local ones...
        tree = spatial.cKDTree(self.gCoords, leafsize=10)
        distances, self.glIDs = tree.query(self.lcoords, k=1)
        nelev = gZ[self.glIDs]

        # From local to global values...
        self.lgIDs = -np.ones(self.gpoints, dtype=int)
        self.lgIDs[self.glIDs] = np.arange(self.npoints)
        self.outIDs = np.where(self.lgIDs < 0)[0]
        self.lgIDs[self.lgIDs < 0] = 0
        del tree, distances

        # Local/Global mapping
        self.lgmap_row = self.dm.getLGMap()
        l2g = self.lgmap_row.indices.copy()
        offproc = l2g < 0
        l2g[offproc] = -(l2g[offproc] + 1)
        self.lgmap_col = PETSc.LGMap().create(l2g, comm=PETSc.COMM_WORLD)

        # Vertex part of an unique partition
        vIS = self.dm.getVertexNumbering()
        self.inIDs = np.zeros(self.npoints, dtype=int)
        self.inIDs[vIS.indices >= 0] = 1
        vIS.destroy()

        # Local/Global vectors
        self.hGlobal = self.dm.createGlobalVector()
        self.hLocal = self.dm.createLocalVector()
        self.sizes = self.hGlobal.getSizes(), self.hGlobal.getSizes()

        self.hLocal.setArray(nelev)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)
        del l2g, offproc, nelev, gZ
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print("Local/global mapping (%0.02f seconds)" % (clock() - t0), flush=True)

        # Create mesh structure with meshplex
        self._meshStructure()

        self.cumED = self.hGlobal.duplicate()
        self.cumED.set(0.0)
        self.cumEDLocal = self.hLocal.duplicate()
        self.cumEDLocal.set(0.0)

        self.sealevel = self.seafunction(self.tNow + self.dt)

        areaLocal = self.hLocal.duplicate()
        self.areaGlobal = self.hGlobal.duplicate()
        areaLocal.setArray(self.area)
        self.dm.localToGlobal(areaLocal, self.areaGlobal)
        areaLocal.destroy()

        # Forcing event number
        self.bG = self.hGlobal.duplicate()
        self.bL = self.hLocal.duplicate()

        self.rainNb = -1
        self.tecNb = -1
        self.flexNb = -1

        return

    def initExtForce(self):
        """
        Initialise external forces.
        """

        self.applyForces()
        if self.backward:
            self.applyTectonics()

        return

    def applyForces(self):
        """
        Find the different values for climatic and sea-level forces that will
        be applied to the considered time interval
        """

        t0 = clock()
        # Sea level
        self.sealevel = self.seafunction(self.tNow + self.dt)

        # Climate
        self._updateRain()

        if MPIrank == 0 and self.verbose:
            print(
                "Update Climatic Forces (%0.02f seconds)" % (clock() - t0), flush=True
            )

        return

    def applyTectonics(self):
        """
        Find the different values for tectonic forces that will be applied to
        the considered time interval
        """

        t0 = clock()
        self._updateTectonic()

        if MPIrank == 0 and self.verbose:
            print(
                "Update Tectonic Forces (%0.02f seconds)" % (clock() - t0), flush=True
            )

        return

    def updatePaleomap(self):
        """
        Force model to match paleomaps at given time interval
        """

        for k in range(self.paleoNb):
            if self.tNow == self.paleodata.iloc[k, 0]:
                loadData = np.load(self.paleodata.iloc[k, 1])
                gZ = loadData["z"]
                self.hLocal.setArray(gZ[self.glIDs])
                self.dm.localToGlobal(self.hLocal, self.hGlobal)
                del loadData, gZ
                gc.collect()

        return

    def _updateRain(self):
        """
        Find current rain values for the considered time interval.
        """

        nb = self.rainNb
        if nb < len(self.raindata) - 1:
            if self.raindata.iloc[nb + 1, 0] <= self.tNow + self.dt:
                nb += 1

        if nb > self.rainNb or nb == -1:
            if nb == -1:
                nb = 0

            self.rainNb = nb
            if pd.isnull(self.raindata["rUni"][nb]):
                loadData = np.load(self.raindata.iloc[nb, 2])
                rainArea = loadData[self.raindata.iloc[nb, 3]]
                del loadData
            else:
                rainArea = np.full(self.gpoints, self.raindata.iloc[nb, 1])
            self.rainArea = rainArea[self.glIDs] * self.area

        self.rainVal = self.rainArea / self.area
        localZ = self.hLocal.getArray()
        rainArea = self.rainArea.copy()
        rainArea[localZ < self.sealevel] = 0.0

        self.bL.setArray(rainArea)
        self.dm.localToGlobal(self.bL, self.bG, 1)

        del rainArea, localZ
        gc.collect()

        return

    def _updateTectonic(self):
        """
        Find the current tectonic rates for the considered time interval.
        """

        if self.tecdata is None:
            self.tectonic = None
            return

        nb = self.tecNb
        if nb < len(self.tecdata) - 1:
            if self.tecdata.iloc[nb + 1, 0] < self.tNow + self.dt:
                nb += 1

        if nb > self.tecNb or nb == -1:
            if nb == -1:
                nb = 0

            self.tecNb = nb
            self.upsubs = False
            if nb < len(self.tecdata.index) - 1:
                timer = self.tecdata.iloc[nb + 1, 0] - self.tecdata.iloc[nb, 0]
            else:
                timer = self.tEnd - self.tecdata.iloc[nb, 0]

            mdata = None
            if self.tecdata.iloc[nb, 1] != "empty":
                mdata = np.load(self.tecdata.iloc[nb, 1])
                self.hdisp = mdata["xyz"][self.glIDs, :]
                self._meshAdvectorSphere(mdata["xyz"], timer)

            if self.tecdata.iloc[nb, 2] != "empty":
                mdata = np.load(self.tecdata.iloc[nb, 2])
                self._meshUpliftSubsidence(mdata["z"])
                self.upsubs = True

            if mdata is not None:
                del mdata

        elif self.upsubs and self.tNow + self.dt < self.tEnd:
            tmp = self.hLocal.getArray().copy()
            self.hLocal.setArray(tmp + self.uplift * self.dt)
            self.dm.localToGlobal(self.hLocal, self.hGlobal)
            del tmp
            gc.collect()

        elif self.forceStep >= 0 and not self.newForcing:
            self._forceUpliftSubsidence()

        return

    def _meshUpliftSubsidence(self, tectonic):
        """
        Apply vertical displacements based on tectonic rates.

        :arg tectonic: local tectonic rates
        """

        # Define vertical displacements
        tmp = self.hLocal.getArray().copy()
        self.uplift = tectonic[self.glIDs]
        self.hLocal.setArray(tmp + self.uplift * self.dt)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)
        del tmp
        gc.collect()

        return

    def _forceUpliftSubsidence(self):
        """
        Perform tectonic uplift at given tectonic interval.
        """

        # Define vertical displacements
        tmp = self.hLocal.getArray().copy()
        self.hLocal.setArray(tmp + self.uplift * self.dt)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)
        del tmp
        gc.collect()

        return

    def _meshAdvectorSphere(self, tectonic, timer):
        """
        Advect spherical mesh horizontally and interpolate mesh information.

        :arg tectonic: local tectonic rates
        :arg timer: tectonic time step in years
        """

        # Move coordinates
        XYZ = self.gCoords + tectonic * timer

        # Elevation
        tmp = self.hLocal.getArray().copy()
        elev = np.zeros(self.gpoints)
        elev = tmp[self.lgIDs]
        elev[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, elev, op=MPI.MAX)

        # Erosion/deposition
        tmp = self.cumEDLocal.getArray().copy()
        erodep = np.zeros(self.gpoints)
        erodep = tmp[self.lgIDs]
        erodep[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, erodep, op=MPI.MAX)

        # Build kd-tree
        kk = 3
        tree = spatial.cKDTree(XYZ, leafsize=10)
        distances, indices = tree.query(self.lcoords, k=kk)

        # Inverse weighting distance...
        if kk == 1:
            nelev = elev[indices]
            nerodep = erodep[indices]
        else:
            weights = np.divide(
                1.0, distances, out=np.zeros_like(distances), where=distances != 0
            )
            onIDs = np.where(distances[:, 0] == 0)[0]
            tmp = np.sum(weights * elev[indices], axis=1)
            tmp2 = np.sum(weights, axis=1)
            nelev = np.divide(tmp, tmp2, out=np.zeros_like(tmp2), where=tmp2 != 0)
            tmp = np.sum(weights * erodep[indices], axis=1)
            nerodep = np.divide(tmp, tmp2, out=np.zeros_like(tmp2), where=tmp2 != 0)
            if len(onIDs) > 0:
                nelev[onIDs] = elev[indices[onIDs, 0]]
                nerodep[onIDs] = erodep[indices[onIDs, 0]]
            del weights, tmp, tmp2

        self.hLocal.setArray(nelev)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

        self.cumEDLocal.setArray(nerodep)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        del XYZ, elev, nelev, erodep, nerodep
        del tree, distances, indices
        gc.collect()

        return

    def destroy_DMPlex(self):
        """
        Destroy PETSc DMPlex objects and associated Petsc local/global
        Vectors and Matrices.
        """

        t0 = clock()

        self.hLocal.destroy()
        self.hGlobal.destroy()
        self.FAG.destroy()
        self.FAL.destroy()
        self.FillG.destroy()
        self.FillL.destroy()
        self.cumED.destroy()
        self.cumEDLocal.destroy()
        self.vSed.destroy()
        self.vSedLocal.destroy()
        self.areaGlobal.destroy()
        self.bG.destroy()
        self.bL.destroy()
        self.hOld.destroy()
        self.hOldLocal.destroy()
        self.stepED.destroy()
        self.tmpL.destroy()
        self.tmp.destroy()
        self.Eb.destroy()
        self.EbLocal.destroy()
        self.vGlob.destroy()

        self.iMat.destroy()
        if not self.fast:
            self.wMat.destroy()
            self.wMat0.destroy()
        self.Diff.destroy()
        self.lgmap_col.destroy()
        self.lgmap_row.destroy()
        self.dm.destroy()

        del self.lcoords
        del self.lcells
        del self.inIDs

        if not self.fast:
            del self.distRcv
            del self.wghtVal
            del self.wghtVal0
            del self.rcvID
            del self.rcvID0
            del self.slpRcv

        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Cleaning Model Dataset (%0.02f seconds)" % (clock() - t0), flush=True
            )

        if self.showlog:
            self.log.view()

        if MPIrank == 0:
            print(
                "\n+++\n+++ Total run time (%0.02f seconds)\n+++"
                % (clock() - self.modelRunTime),
                flush=True,
            )

        return
