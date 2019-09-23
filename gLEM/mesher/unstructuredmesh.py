import gc
import numpy as np
import pandas as pd
from mpi4py import MPI

import sys,petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc
from time import clock

from gLEM._fortran import defineTIN
from gLEM._fortran import ngbGlob

import meshplex
from pykdtree.kdtree import KDTree

MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIsize = PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD

try: range = xrange
except: pass

class UnstMesh(object):
    """
    Creating a distributed DMPlex and global vector from it based on triangulated mesh
    """
    def __init__(self):

        self.icosahedralSphereMesh()

        return

    def meshfrom_cell_list(self, dim, cells, coords):
        """
        Create a DMPlex from a list of cells and coords.

        :arg dim: The topological dimension of the mesh
        :arg cells: The vertices of each cell
        :arg coords: The coordinates of each vertex
        :arg comm: communicator to build the mesh on.
        """

        if MPIrank == 0:
            cells = np.asarray(cells, dtype=np.int32)
            coords = np.asarray(coords, dtype=np.double)
            MPIcomm.bcast(cells.shape, root=0)
            MPIcomm.bcast(coords.shape, root=0)
            # Provide the actual data on rank 0.
            self.dm = PETSc.DMPlex().createFromCellList(dim, cells, coords, comm=PETSc.COMM_WORLD)
            del cells, coords
        else:
            cell_shape = list(MPIcomm.bcast(None, root=0))
            coord_shape = list(MPIcomm.bcast(None, root=0))
            cell_shape[0] = 0
            coord_shape[0] = 0
            self.dm = PETSc.DMPlex().createFromCellList(dim,
                                        np.zeros(cell_shape, dtype=np.int32),
                                        np.zeros(coord_shape, dtype=np.double),
                                        comm=PETSc.COMM_WORLD)

        return

    def readDataMesh(self):

        t0 = clock()
        loadData = np.load(self.meshFile)
        self.gCoords = loadData['v']
        self.gpoints = len(self.gCoords)
        gZ = loadData['z']
        ngbGlob(self.gpoints,loadData['n'])
        del loadData
        # self.gCoords = self.swarmXYZ.copy()

        # From global values to local ones...
        tree = KDTree(self.gCoords,leafsize=10)
        distances, self.glIDs = tree.query(self.lcoords, k=1)
        nelev = gZ[self.glIDs] #self.swarmZ[self.glIDs]

        # From local to global values...
        tree = KDTree(self.lcoords,leafsize=10)
        distances, self.lgIDs = tree.query(self.gCoords, k=1)
        self.outIDs = np.where(distances>0.1)[0]
        del tree, distances

        # Local/Global mapping
        self.lgmap_row = self.dm.getLGMap()
        l2g = self.lgmap_row.indices.copy()
        offproc = l2g < 0
        l2g[offproc] = -(l2g[offproc] + 1)
        self.lgmap_col = PETSc.LGMap().create(l2g, comm=PETSc.COMM_WORLD)

        # Vertex part of an unique partition
        vIS = self.dm.getVertexNumbering()
        self.inIDs = np.zeros(self.npoints,dtype=int)
        self.inIDs[vIS.indices>=0] = 1
        vIS.destroy()

        # Local/Global vectors
        self.hGlobal = self.dm.createGlobalVector()
        self.hLocal = self.dm.createLocalVector()
        self.sizes = self.hGlobal.getSizes(), self.hGlobal.getSizes()

        self.hLocal.setArray(nelev)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)
        del l2g, offproc, nelev #, gZ
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print('local/global mapping (%0.02f seconds)'% (clock() - t0))

        return

    def icosahedralSphereMesh(self):
        """Generate an icosahedral approximation to the surface of the
        sphere.
        """

        t0 = clock()
        if self.reflevel < 0 or self.reflevel % 1:
            raise RuntimeError("Number of refinements must be a non-negative integer")

        from math import sqrt
        phi = (1 + sqrt(5)) / 2

        # vertices of an icosahedron with an edge length of 2
        vertices = np.array([[-1, phi, 0],
                             [1, phi, 0],
                             [-1, -phi, 0],
                             [1, -phi, 0],
                             [0, -1, phi],
                             [0, 1, phi],
                             [0, -1, -phi],
                             [0, 1, -phi],
                             [phi, 0, -1],
                             [phi, 0, 1],
                             [-phi, 0, -1],
                             [-phi, 0, 1]],
                            dtype=np.double)

        # faces of the base icosahedron
        faces = np.array([[0, 11, 5],
                          [0, 5, 1],
                          [0, 1, 7],
                          [0, 7, 10],
                          [0, 10, 11],
                          [1, 5, 9],
                          [5, 11, 4],
                          [11, 10, 2],
                          [10, 7, 6],
                          [7, 1, 8],
                          [3, 9, 4],
                          [3, 4, 2],
                          [3, 2, 6],
                          [3, 6, 8],
                          [3, 8, 9],
                          [4, 9, 5],
                          [2, 4, 11],
                          [6, 2, 10],
                          [8, 6, 7],
                          [9, 8, 1]], dtype=np.int32)

        self.meshfrom_cell_list(2, faces, vertices)
        del faces, vertices

        # Distribute to other processors if any
        if MPIsize > 1:
            t0 = clock()
            sf = self.dm.distribute(overlap=0)

        for i in range(self.reflevel):
            self.dm.setRefinementUniform(True)
            self.dm = self.dm.refine()

        if MPIsize > 1:
            self.dm.distributeOverlap(1)
            if MPIrank == 0 and self.verbose:
                print('create mesh (%0.02f seconds)'% (clock() - t0))

        # Define one DoF on the nodes
        t0 = clock()
        self.dm.setNumFields(1)
        origSect = self.dm.createSection(1, [1,0,0])
        origSect.setFieldName(0, "points")
        origSect.setUp()
        self.dm.setDefaultSection(origSect)
        origSect.destroy()

        self.lcoords = self.dm.getCoordinatesLocal().array.reshape(-1, 3)
        scale = (self.radius / np.linalg.norm(self.lcoords, axis=1)).reshape(-1, 1)
        self.lcoords *= scale
        self.npoints = self.lcoords.shape[0]

        # Define local vertex & cells
        cStart, cEnd = self.dm.getHeightStratum(0)
        self.lcells = np.zeros((cEnd-cStart,3), dtype=PETSc.IntType)
        for c in range(cStart, cEnd):
            point_closure = self.dm.getTransitiveClosure(c)[0]
            self.lcells[c,:] = point_closure[-3:]-cEnd
        del point_closure
        gc.collect()
        if MPIrank == 0 and self.verbose:
            print('mesh coords/cells (%0.02f seconds)'% (clock() - t0))

        # Interpolate values from input using inverse weighting distance
        self.readDataMesh()

        # Create mesh structure with meshplex
        t0 = clock()
        Tmesh = meshplex.MeshTri(self.lcoords,self.lcells)
        self.area = np.abs(Tmesh.control_volumes)
        self.area[np.isnan(self.area)] = 1.

        # Voronoi and simplices declaration
        Tmesh.create_edges()
        cc = Tmesh.cell_circumcenters
        edges_nodes = Tmesh.edges['nodes']
        cells_nodes = Tmesh.cells['nodes']
        cells_edges = Tmesh.cells['edges']
        # Finite volume discretisation
        self.FVmesh_ngbID = defineTIN(self.lcoords, cells_nodes, cells_edges,
                                      edges_nodes, self.area, cc.T)
        del Tmesh, edges_nodes, cells_nodes, cells_edges, cc
        gc.collect()
        if MPIrank == 0 and self.verbose:
            print('fv discretisation (%0.02f seconds)'% (clock() - t0))

        self.cumED = self.hGlobal.duplicate()
        self.cumED.set(0.0)
        self.cumEDLocal = self.hLocal.duplicate()
        self.cumEDLocal.set(0.0)

        self.sealevel = self.seafunction(self.tNow+self.dt)

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

        return

    def applyForces(self):
        """
        Find the different values for climatic and tectonic forces that will be applied to the
        considered time interval
        """

        t0 = clock()
        # Sea level
        self.sealevel = self.seafunction(self.tNow+self.dt)

        # Climate
        self._updateRain()

        # Tectonic
        self._updateTectonic()

        if MPIrank == 0 and self.verbose:
            print('Update External Forces (%0.02f seconds)'% (clock() - t0))

        return

    def updatePaleomap(self):
        """
        Find the paleomap for the considered time interval
        """

        if self.paleodata is None :
            return

        for k in range(self.paleoNb):
            if self.tNow == self.paleodata.iloc[k,0]:
                loadData = np.load(self.paleodata.iloc[k,1])
                gZ = loadData['z']
                self.hLocal.setArray(gZ[self.glIDs])
                self.dm.localToGlobal(self.hLocal, self.hGlobal)
                del loadData, gZ
                gc.collect()

        return

    def _updateRain(self):
        """
        Find the current rain values for the considered time interval
        """

        nb = self.rainNb
        if nb < len(self.raindata)-1 :
            if self.raindata.iloc[nb+1,0] <= self.tNow+self.dt :
                nb += 1

        if nb > self.rainNb or nb == -1:
            if nb == -1:
                nb = 0

            self.rainNb = nb
            if pd.isnull(self.raindata['rUni'][nb]):
                loadData = np.load(self.raindata.iloc[nb,2])
                rainArea = loadData[self.raindata.iloc[nb,3]]
                del loadData
            else:
                rainArea = np.full(self.gpoints,self.raindata.iloc[nb,1])
            self.rainArea = rainArea[self.glIDs]*self.area

        localZ = self.hLocal.getArray()
        rainArea = self.rainArea.copy()
        rainArea[localZ<self.sealevel] = 0.

        self.bL.setArray(rainArea)
        self.dm.localToGlobal(self.bL, self.bG, 1)

        del rainArea, localZ
        gc.collect()

        return

    def _updateTectonic(self):
        """
        Find the current tectonic rates for the considered time interval
        """

        if self.tecdata is None :
            self.tectonic = None
            return

        nb = self.tecNb
        if nb < len(self.tecdata)-1 :
            if self.tecdata.iloc[nb+1,0] <= self.tNow+self.dt :
                nb += 1

        if nb > self.tecNb or nb == -1:
            if nb == -1:
                nb = 0

            self.tecNb = nb
            mdata = np.load(self.tecdata.iloc[nb,1])

            if nb < len(self.tecdata.index)-1:
                timer = self.tecdata.iloc[nb+1,0]-self.tecdata.iloc[nb,0]
            else:
                timer = self.tEnd-self.tecdata.iloc[nb,0]

            self._meshAdvectorSphere(mdata['xyz'], timer)
            del mdata

        return

    def _meshAdvectorSphere(self, tectonic, timer):
        """
        Advect spherical mesh horizontally and interpolate mesh information
        """

        # Move coordinates
        XYZ = self.gCoords + tectonic*timer

        # Elevation
        tmp = self.hLocal.getArray().copy()
        elev = np.zeros(self.gpoints)
        elev = tmp[self.lgIDs]
        elev[self.outIDs] = -1.e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, elev, op=MPI.MAX)

        # Erosion/deposition
        tmp = self.cumEDLocal.getArray().copy()
        erodep = np.zeros(self.gpoints)
        erodep = tmp[self.lgIDs]
        erodep[self.outIDs] = -1.e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, erodep, op=MPI.MAX)

        # Build kd-tree
        kk = 3
        tree = KDTree(XYZ,leafsize=10)
        distances, indices = tree.query(self.lcoords, k=kk)

        # Inverse weighting distance...
        if kk == 1:
            nelev = elev[indices]
            nerodep = erodep[indices]
        else:
            weights = 1.0 / distances**2
            onIDs = np.where(distances[:,0] == 0)[0]
            nelev = np.sum(weights*elev[indices],axis=1)/np.sum(weights, axis=1)
            nerodep = np.sum(weights*erodep[indices],axis=1)/np.sum(weights, axis=1)
            if len(onIDs)>0:
                nelev[onIDs] = elev[indices[onIDs,0]]
                nerodep[onIDs] = erodep[indices[onIDs,0]]
            del weights

        self.hLocal.setArray(nelev)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

        self.cumEDLocal.setArray(nerodep)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        del XYZ,  elev, nelev, erodep, nerodep
        del tree, distances, indices
        gc.collect()

        return

    def destroy_DMPlex(self):
        """
        Destroy PETSc DMPlex objects and associated Petsc local/global Vectors and Matrices.
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

        self.iMat.destroy()
        self.wMat.destroy()
        self.wMat0.destroy()
        self.lgmap_col.destroy()
        self.lgmap_row.destroy()
        self.dm.destroy()

        del self.lcoords
        del self.lcells
        del self.inIDs
        del self.distRcv
        del self.wghtVal
        del self.rcvID
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print('Cleaning Model Dataset (%0.02f seconds)'% (clock() - t0))

        return
