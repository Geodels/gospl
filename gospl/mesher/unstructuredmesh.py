import os
os.environ["VTK_USE_VISKORES"] = "0"
import gc
import sys
import vtk
vtk.vtkObject.GlobalWarningDisplayOff()

import warnings
import petsc4py
import numpy as np
import pandas as pd

from mpi4py import MPI
from scipy import spatial
from time import process_time

from vtk.util import numpy_support  # type: ignore

from gospl.tools.constants import MISSING_DATA_SENTINEL

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import globalngbhs
    from gospl._fortran import definetin
    from gospl._fortran import fitedges
    from gospl._fortran import updatearea

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class UnstMesh(object):
    """
    This class defines the 2D or spherical mesh characteristics and builds a PETSc DMPlex that encapsulates this unstructured mesh, with interfaces for both topology and geometry. The PETSc DMPlex is used for parallel redistribution and for load balancing.

    .. note::

        goSPL is built around a **Finite-Volume** method (FVM) for representing and evaluating  partial differential equations. It requires the definition of several mesh variables such as:

            - the number of neighbours surrounding every node,
            - the cell area defined using Voronoi area,
            - the length of the edge connecting every nodes, and
            - the length of the Voronoi faces shared by each node with his neighbours.

    In addition to mesh defintions, the class declares several functions related to forcing conditions (*e.g.* paleo-precipitation maps, tectonic (vertical and horizontal) displacements, stratigraphic layers...). These functions are defined within the `UnstMesh` class as they rely heavily on the mesh structure.

    .. important::

        The grid (2D or spherical) requires locally-orthogonal Voronoi/Delaunay staggering, or an unstructured C-grid type numerical formulation as described in `Engwirda 2017 <https://arxiv.org/pdf/1611.08996>`_

    Finally a function to clean all PETSc variables is defined and called at the end of a simulation.
    """

    def __init__(self):
        """
        The initialisation of `UnstMesh` class calls the private function **_buildMesh**.
        """

        self.upsub = None
        self.rainVal = None
        self.evapVal = None
        self.evapLoss = 0.0
        self.sedfacVal = None
        self.memclear = False
        self.southPts = None

        # Define the mesh variables and build PETSc DMPLEX.
        self._buildMesh()

        return

    def _meshfrom_cell_list(self, dim, cells, coords):
        """
        Creates a DMPlex from a list of cells and coordinates.

        .. note::

            PETSc DMPlex requires to be initialised on one processor before load balancing.

        :arg dim: topological dimension of the mesh
        :arg cells: vertices of each cell
        :arg coords: coordinates of each vertex
        """

        if MPIrank == 0:
            cells = np.asarray(cells, dtype=np.int32)
            coords = np.asarray(coords, dtype=np.float64)
            MPIcomm.bcast(cells.shape, root=0)
            MPIcomm.bcast(coords.shape, root=0)
            # Provide the actual data on rank 0.
            self.dm = petsc4py.PETSc.DMPlex().createFromCellList(
                dim, cells, coords, comm=petsc4py.PETSc.COMM_WORLD
            )
            del cells, coords
        else:
            cell_shape = list(MPIcomm.bcast(None, root=0))
            coord_shape = list(MPIcomm.bcast(None, root=0))
            cell_shape[0] = 0
            coord_shape[0] = 0
            self.dm = petsc4py.PETSc.DMPlex().createFromCellList(
                dim,
                np.zeros(cell_shape, dtype=np.int32),
                np.zeros(coord_shape, dtype=np.float64),
                comm=petsc4py.PETSc.COMM_WORLD,
            )

        return

    def _meshStructure(self):
        """
        Defines the mesh structure and the associated voronoi parameters used in the Finite Volume method.

        .. important::
            The mesh structure is built locally on a single partition of the global mesh.

        Once the voronoi definitions have been obtained a call to the fortran subroutine `definetin` is performed to order each node and the dual mesh components, it records:

        - all cells surrounding a given vertice,
        - all edges connected to a given vertice,
        - the triangulation edge lengths,
        - the voronoi edge lengths.
        """

        # Create mesh structure and voronoi parameters used for
        # Centroidal Voronoi Tessellation and Spherical Centroidal Voronoi
        # Tessellation
        t0 = process_time()
        Tmesh = self.initVoronoi(self.lcoords, self.lcells)
        larea = np.abs(self.control_volumes)
        larea[np.isnan(larea)] = 1.0

        # Voronoi and simplices declaration
        self.create_edges()
        cc = self.cell_circumcenters
        if not self.flatModel:
            # Ensure voronoi points are properly set on the sphere
            radius = np.linalg.norm(self.lcoords[0])
            cc = cc * (radius / np.linalg.norm(cc, axis=1)).reshape((len(cc), 1))

        edges_nodes = self.edges["nodes"]
        cells_nodes = self.cells["nodes"]
        cells_edges = self.cells["edges"]

        # Finite volume discretisation
        self.FVmesh_ngbID, self.larea = definetin(
            self.lcoords,
            cells_nodes,
            cells_edges,
            edges_nodes,
            cc.T,
        )
        self.larea[np.isnan(self.larea)] = 1.0
        issues = np.zeros(1)
        issues[0] = np.max(np.abs(larea - self.larea))
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, issues, op=MPI.MAX)
        if MPIrank == 0 and issues[0] > 1.e2 and self.flatModel:
            print(
                "\n--------------\n"
                "Warning:\n"
                "Some issues have been encountered in the Finite Volume declaration.\n"
                "This is likely due to your initial grid discretisation. Use some of\n"
                "the meshing approaches proposed in the pre-processing workflow.\n"
                "Your grid needs to be a Delaunay with optimal voronoi (C-grid).\n"
                "--------------\n",
                flush=True
            )
        if issues[0] > 1.e2 and self.flatModel:
            self.larea = larea.copy()
            updatearea(larea)

        self.maxarea = np.zeros(1, dtype=np.float64)
        self.maxarea[0] = self.larea.max()
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, self.maxarea, op=MPI.MAX)

        del Tmesh, edges_nodes, cells_nodes, cells_edges, cc, larea
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "FV discretisation (%0.02f seconds)" % (process_time() - t0), flush=True
            )

        return

    def _generateVTKmesh(self, points, cells):
        """
        A global VTK mesh is generated to compute the distance between mesh vertices and coastlines position.

        .. note::
            The distance to the coastline for every marine vertices is used to define a maximum shelf slope during deposition.

            The coastline contours are efficiently obtained from VTK contouring function.
        """

        self.vtkMesh = vtk.vtkUnstructuredGrid()

        # Define mesh vertices
        vtk_points = vtk.vtkPoints()
        vtk_array = numpy_support.numpy_to_vtk(points, deep=True)
        vtk_points.SetData(vtk_array)
        self.vtkMesh.SetPoints(vtk_points)

        # Define mesh cells
        cell_types = []
        cell_offsets = []
        cell_connectivity = []
        len_array = 0
        numcells, num_local_nodes = cells.shape
        cell_types.append(np.empty(numcells, dtype=np.ubyte))
        cell_types[-1].fill(vtk.VTK_TRIANGLE)
        cell_offsets.append(
            np.arange(
                len_array,
                len_array + numcells * (num_local_nodes + 1),
                num_local_nodes + 1,
                dtype=np.int64,
            )
        )
        cell_connectivity.append(
            np.c_[
                num_local_nodes * np.ones(numcells, dtype=cells.dtype), cells
            ].flatten()
        )
        len_array += len(cell_connectivity[-1])
        cell_types = np.concatenate(cell_types)
        cell_offsets = np.concatenate(cell_offsets)
        cell_connectivity = np.concatenate(cell_connectivity)

        # Connectivity
        connectivity = numpy_support.numpy_to_vtkIdTypeArray(
            cell_connectivity.astype(np.int64), deep=1
        )
        cell_array = vtk.vtkCellArray()
        cell_array.SetCells(len(cell_types), connectivity)

        self.vtkMesh.SetCells(
            numpy_support.numpy_to_vtk(
                cell_types, deep=1, array_type=vtk.vtkUnsignedCharArray().GetDataType()
            ),
            numpy_support.numpy_to_vtk(
                cell_offsets, deep=1, array_type=vtk.vtkIdTypeArray().GetDataType()
            ),
            cell_array,
        )

        # Cleaning function parameters here...
        del vtk_points, vtk_array, connectivity, numcells, num_local_nodes
        del cell_array, cell_connectivity, cell_offsets, cell_types
        gc.collect()

        return

    def _buildMesh(self):
        """
        This function is at the core of the `UnstMesh` class. It encapsulates both mesh construction (triangulation and voronoi representation for the Finite Volume discretisation), PETSc DMPlex distribution and several PETSc vectors allocation.

        The function relies on several private functions from the class:

        - _generateVTKmesh
        - _meshfrom_cell_list
        - _meshStructure
        - _readErosionDeposition
        - readStratLayers

        .. note::

            It is worth mentionning that partitioning and field distribution from global to local PETSc DMPlex takes a lot of time for large mesh.
        """

        # Read mesh attributes from file
        t0 = process_time()
        loadData = np.load(self.meshFile)
        self.mCoords = loadData[self.infoCoords]
        self.mpoints = len(self.mCoords)
        gZ = loadData[self.infoElev]
        mCells = loadData[self.infoCells].astype(int)

        # Get global mesh vertex neighbors
        if MPIrank == 0:
            globalngbhs(self.mpoints, mCells)
        mCells = None
        self.vtkMesh = None
        self.flatModel = False
        if MPIrank == 0 and self.verbose:
            print(
                "Reading mesh information (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Create DMPlex
        t0 = process_time()
        self._meshfrom_cell_list(2, loadData[self.infoCells], self.mCoords)
        del loadData
        gc.collect()
        if MPIrank == 0 and self.verbose:
            print("Create DMPlex (%0.02f seconds)" % (process_time() - t0), flush=True)

        # Define one degree of freedom on the nodes
        t0 = process_time()
        self.dm.setNumFields(1)
        origSect = self.dm.createSection(1, [1, 0, 0])
        origSect.setFieldName(0, "points")
        origSect.setUp()
        self.dm.setDefaultSection(origSect)
        origVec = self.dm.createGlobalVector()
        if MPIrank == 0 and self.verbose:
            print(
                "Define one DoF on the nodes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Distribute to other processors if any
        t0 = process_time()
        if MPIsize > 1:
            partitioner = self.dm.getPartitioner()
            partitioner.setType(partitioner.Type.PARMETIS)
            partitioner.setFromOptions()
            sf = self.dm.distribute(overlap=self.overlap)
            newSect, newVec = self.dm.distributeField(sf, origSect, origVec)
            self.dm.setDefaultSection(newSect)
            newSect.destroy()
            newVec.destroy()
            sf.destroy()
        MPIcomm.Barrier()
        origVec.destroy()
        origSect.destroy()
        if MPIrank == 0 and self.verbose:
            print(
                "Distribute DMPlex (%0.02f seconds)" % (process_time() - t0), flush=True
            )

        # Define local vertex & cells
        t0 = process_time()
        self.gcoords = self.dm.getCoordinates().array.reshape(-1, 3)
        self.lcoords = self.dm.getCoordinatesLocal().array.reshape(-1, 3)
        self.gpoints = self.gcoords.shape[0]
        self.lpoints = self.lcoords.shape[0]

        cStart, cEnd = self.dm.getHeightStratum(0)
        nlocal_cells = cEnd - cStart
        self.lcells = np.zeros((nlocal_cells, 3), dtype=petsc4py.PETSc.IntType)

        # For an uninterpolated 2D DMPlex, `getCone(c)` is the same 3
        # vertices as `getTransitiveClosure(c)[0][-3:]` and is ~5x cheaper
        # per call (no closure expansion). We probe the first cell to
        # confirm both calls agree (same values in the same order, which
        # preserves cell orientation downstream); if they don't, fall
        # back to the original closure walk.
        use_cone = False
        if nlocal_cells > 0:
            closure0 = self.dm.getTransitiveClosure(cStart)[0][-3:]
            cone0 = self.dm.getCone(cStart)
            if len(cone0) == 3 and np.array_equal(cone0, closure0):
                use_cone = True

        if use_cone:
            for c in range(cStart, cEnd):
                self.lcells[c - cStart, :] = self.dm.getCone(c) - cEnd
        else:
            point_closure = None
            for c in range(cStart, cEnd):
                point_closure = self.dm.getTransitiveClosure(c)[0]
                self.lcells[c - cStart, :] = point_closure[-3:] - cEnd
            if point_closure is not None:
                del point_closure
        gc.collect()
        if MPIrank == 0 and self.verbose:
            print(
                "Defining local DMPlex (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Build local VTK mesh
        if self.clinSlp > 0.0:
            t0 = process_time()
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self._generateVTKmesh(self.lcoords, self.lcells)
            if MPIrank == 0 and self.verbose:
                print(
                    "Generate VTK mesh (%0.02f seconds)" % (process_time() - t0),
                    flush=True,
                )

        # From mesh values to local and global ones...
        t0 = process_time()
        # Larger leafsize + skipping the rebalance/compact passes cuts the
        # cKDTree build by ~2-3x for million-vertex meshes; the per-query
        # cost grows only marginally because the queries here are k=1
        # nearest-neighbour on an essentially identical point set.
        tree = spatial.cKDTree(
            self.mCoords,
            leafsize=32,
            balanced_tree=False,
            compact_nodes=False,
        )
        distances, self.locIDs = tree.query(self.lcoords, k=1)
        distances, self.glbIDs = tree.query(self.gcoords, k=1)

        # Local/Global mapping
        self.lgmap_row = self.dm.getLGMap()
        l2g = self.lgmap_row.indices.copy()
        offproc = l2g < 0
        l2g[offproc] = -(l2g[offproc] + 1)
        self.lgmap_col = petsc4py.PETSc.LGMap().create(
            l2g, comm=petsc4py.PETSc.COMM_WORLD
        )
        # Per-local-node global vertex ID, consistent across MPI
        # decompositions (PETSc DMPlex guarantees a node's global ID is
        # the same on every rank that has it as either owned or ghost).
        # Used by `mfdreceivers` / `mfdrcvrs` to break exact slope ties
        # deterministically — without this, quicksort can pick a different
        # receiver on different ranks for the same global node and
        # cascade into divergent drainage statistics under partitioning.
        # int32 because f2py declares the Fortran arg as default integer.
        self.gid = l2g.astype(np.int32)

        # Vertex part of an unique partition
        vIS = self.dm.getVertexNumbering()
        self.inIDs = np.zeros(self.lpoints, dtype=int)
        self.inIDs[vIS.indices >= 0] = 1
        # glIDs are the indices part of a local partition which are not ghosts
        self.glIDs = np.where(self.inIDs == 1)[0]
        # ghostIDs are the shadow nodes indices on each partition
        self.ghostIDs = np.where(self.inIDs == 0)[0]

        # Local mesh boundary points in 2D model
        localBound = self._get_boundary()
        idLocal = np.where(vIS.indices >= 0)[0]
        self.idBorders = np.where(np.isin(idLocal, localBound))[0]
        nib = np.zeros(1, dtype=np.int64)
        nib[0] = len(self.idBorders)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, nib, op=MPI.MAX)
        if nib[0] > 0:
            self.flatModel = True
            self.south = int(self.boundCond[0])
            self.east = int(self.boundCond[1])
            self.north = int(self.boundCond[2])
            self.west = int(self.boundCond[3])
            xmin = self.mCoords[:, 0].min()
            xmax = self.mCoords[:, 0].max()
            ymin = self.mCoords[:, 1].min()
            ymax = self.mCoords[:, 1].max()
            # Use np.isclose so boundary detection survives any future
            # coordinate transformation that introduces float drift.
            tol = 1.0e-9
            self.southPts = np.where(np.isclose(self.lcoords[:, 1], ymin, atol=tol))[0]
            self.northPts = np.where(np.isclose(self.lcoords[:, 1], ymax, atol=tol))[0]
            self.eastPts = np.where(np.isclose(self.lcoords[:, 0], xmax, atol=tol))[0]
            self.westPts = np.where(np.isclose(self.lcoords[:, 0], xmin, atol=tol))[0]

        del idLocal
        vIS.destroy()

        # Local/Global vectors
        self.hGlobal = self.dm.createGlobalVector()
        self.hLocal = self.dm.createLocalVector()
        self.sizes = self.hGlobal.getSizes(), self.hGlobal.getSizes()

        self.hGlobal.setArray(gZ[self.glbIDs])
        self.hLocal.setArray(gZ[self.locIDs])

        if MPIrank == 0 and self.verbose:
            print(
                "Local/global mapping (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Create mesh structure
        self._meshStructure()

        # Get local mesh borders not included in the shadow regions for parallel pit filling
        masknodes = ~np.isin(self.lcells, self.ghostIDs)
        tmp2 = np.sum(masknodes.astype(int), axis=1)
        out = np.where(np.logical_and(tmp2 > 0, tmp2 < 3))[0]
        ptscells = self.lcells[out, :].flatten()
        self.idLBounds = np.setdiff1d(ptscells, self.ghostIDs)

        # Define cumulative erosion deposition arrays
        self._readErosionDeposition()

        self.sealevel = self.seafunction(self.tNow)
        self.areaGlobal = self.hGlobal.duplicate()
        self.areaLocal = self.hLocal.duplicate()
        self.areaLocal.setArray(self.larea)
        self.dm.localToGlobal(self.areaLocal, self.areaGlobal)
        self.dm.globalToLocal(self.areaGlobal, self.areaLocal)
        self.larea = self.areaLocal.getArray().copy()

        # Forcing event number
        self.bG = self.hGlobal.duplicate()
        self.bL = self.hLocal.duplicate()
        self.rainNb = -1
        self.evapNb = -1
        self.flexNb = -1
        self.teNb = -1
        self.sedfactNb = -1

        del tree, distances, tmp2
        del l2g, offproc, gZ, out, ptscells
        gc.collect()

        # Build stratigraphic data if any
        if self.stratNb > 0:
            self.readStratLayers()

        return

    def _set_DMPlex_boundary_points(self, label):
        """
        In case of a 2D mesh, this function finds the points that join the edges that have been marked as "boundary" faces in the directed acyclic graph (DAG) then sets them as boundaries.
        """

        self.dm.createLabel(label)
        self.dm.markBoundaryFaces(label)

        # pStart, pEnd = self.dm.getDepthStratum(0)  # points
        eStart, eEnd = self.dm.getDepthStratum(1)  # edges
        edgeIS = self.dm.getStratumIS(label, 1)

        if edgeIS and eEnd - eStart > 0:
            edge_mask = np.logical_and(edgeIS.indices >= eStart, edgeIS.indices < eEnd)
            boundary_edges = edgeIS.indices[edge_mask]

            # Query the DAG  (directed acyclic graph) for points that join an edge
            for edge in boundary_edges:
                vertices = self.dm.getCone(edge)
                # mark the boundary points
                for vertex in vertices:
                    self.dm.setLabelValue(label, vertex, 1)

        edgeIS.destroy()

        return

    def _get_boundary(self, label="boundary"):
        """
        In case of a 2D mesh, this function finds the nodes on the boundary from the DM.
        """

        label = "boundary"
        self._set_DMPlex_boundary_points(label)

        pStart, pEnd = self.dm.getDepthStratum(0)

        labels = []
        for i in range(self.dm.getNumLabels()):
            labels.append(self.dm.getLabelName(i))

        if label not in labels:
            raise ValueError("There is no {} label in the DM".format(label))

        stratSize = self.dm.getStratumSize(label, 1)
        if stratSize > 0:
            labelIS = self.dm.getStratumIS(label, 1)
            pt_range = np.logical_and(labelIS.indices >= pStart, labelIS.indices < pEnd)
            indices = labelIS.indices[pt_range] - pStart
            labelIS.destroy()
        else:
            indices = np.zeros((0,), dtype=np.int32)

        return indices

    def _readErosionDeposition(self):
        """
        Reads existing cumulative erosion depostion from a previous experiment as defined in the YAML input file following the `nperodep` key.

        This functionality can be used when restarting from a previous simulation in which the mesh has been modified either to account for horizontal advection or to refine/coarsen a specific region during a given time period.
        """

        # Build PETSc vectors
        self.cumED = self.hGlobal.duplicate()
        self.cumED.set(0.0)
        self.cumEDLocal = self.hLocal.duplicate()
        self.cumEDLocal.set(0.0)

        # Read mesh value from file
        if self.dataFile is not None:
            fileData = np.load(self.dataFile)
            gED = fileData["ed"]
            del fileData

            self.cumEDLocal.setArray(gED[self.locIDs])
            self.cumED.setArray(gED[self.glbIDs])
            del gED
            gc.collect()

        return

    def applyForces(self):
        """
        Finds the different values for climatic, tectonic and sea-level forcing that will be applied at any given time interval during the simulation.
        """

        t0 = process_time()
        # Sea level
        if self.tNow == self.tStart:
            self.sealevel = self.seafunction(self.tNow)
            self.oldsealevel = self.sealevel.copy()
        else:
            self.oldsealevel = self.sealevel.copy()
            self.sealevel = self.seafunction(self.tNow + self.dt)

        # Climate information
        if self.oroOn:
            self.cptOrography()
        else:
            self._updateRain()

        # Evaporation forcing (independent of orographic vs uniform rain;
        # populates self.evapVal which the flow solver consumes downstream).
        if self.evapdata is not None:
            self._updateEvap()

        # Erodibility factor information
        if self.sedfacdata is not None:
            self._updateEroFactor()

        # Glacier-geometry maps (per-vertex / time-series hela/hice/hterm).
        if self.iceOn and self._iceTimeSeries is not None:
            self._updateIce()

        if MPIrank == 0 and self.verbose:
            print(
                "Update Climatic Forces (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Tectonic forcing
        self.applyTectonics()

        # Assign mesh boundaries
        if self.flatModel:
            tmp = self.hLocal.getArray().copy()
            if self.south == 0 and len(self.southPts) > 0:
                # tmp[self.southPts] = getbc(len(self.southPts), tmp, self.southPts)
                tmp[self.southPts] = MISSING_DATA_SENTINEL
                tmp = fitedges(tmp)
            if self.north == 0 and len(self.northPts) > 0:
                # tmp[self.northPts] = getbc(len(self.northPts), tmp, self.northPts)
                tmp[self.northPts] = MISSING_DATA_SENTINEL
                tmp = fitedges(tmp)
            if self.east == 0 and len(self.eastPts) > 0:
                # tmp[self.eastPts] = getbc(len(self.eastPts), tmp, self.eastPts)
                tmp[self.eastPts] = MISSING_DATA_SENTINEL
                tmp = fitedges(tmp)
            if self.west == 0 and len(self.westPts) > 0:
                # tmp[self.westPts] = getbc(len(self.westPts), tmp, self.westPts)
                tmp[self.westPts] = MISSING_DATA_SENTINEL
                tmp = fitedges(tmp)
            self.hLocal.setArray(tmp)
            self.dm.localToGlobal(self.hLocal, self.hGlobal)

        return

    def applyTectonics(self):
        """
        Finds the different values for tectonic forces that will be applied at any given time interval during the simulation.
        """

        t0 = process_time()

        if self.upsub is not None:
            # Define vertical displacements
            tmp = self.hLocal.getArray().copy()
            self.hLocal.setArray(tmp + self.upsub * self.dt)
            self.dm.localToGlobal(self.hLocal, self.hGlobal)

        if MPIrank == 0 and self.verbose:
            print(
                "Update Tectonic Forces (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _updateRain(self):
        """
        Finds current rain values for the considered time interval and computes the **volume** of water available for runoff on each vertex.

        .. note::

            It is worth noting that the precipitation maps are considered as runoff water. If one wants to account for evaporation and infiltration you will need to modify the precipitation maps accordingly as a pre-processing step.

        """

        nb = self.rainNb
        if nb < len(self.raindata) - 1:
            if self.raindata.at[nb + 1, "start"] <= self.tNow:  # + self.dt:
                nb += 1

        if nb > self.rainNb or nb == -1:
            if nb == -1:
                nb = 0

            self.rainNb = nb
            if pd.isnull(self.raindata["rUni"][nb]):
                if pd.isnull(self.raindata["rzA"][nb]):
                    loadData = np.load(self.raindata.at[nb, "rMap"])
                    rainVal = loadData[self.raindata.at[nb, "rKey"]]
                    self.rainA = None
                    self.rainB = None
                    del loadData
                else:
                    self.rainA = self.raindata.at[nb, "rzA"]
                    self.rainB = self.raindata.at[nb, "rzB"]
            else:
                self.rainA = None
                self.rainB = None
                rainVal = np.full(self.mpoints, self.raindata.at[nb, "rUni"])
            
            if self.rainA is None:
                rainVal[rainVal < 0] = 0.0
                self.rainMesh = rainVal

        if self.rainA is None:
            self.rainVal = self.rainMesh[self.locIDs]
        else:
            tmp = self.hLocal.getArray().copy()
            self.rainVal = tmp  * self.rainA + self.rainB
            self.rainVal[self.rainVal < 0] = 0.0
        self.bL.setArray(self.rainVal * self.larea)
        self.dm.localToGlobal(self.bL, self.bG)

        return

    def _updateIce(self):
        """
        Refresh the glacier-geometry fields (terminus / ELA / ice-cap altitude)
        for the current time from the ``_iceTimeSeries`` built in the input
        parser — the ice analogue of ``_updateRain``.

        Each interval supplies, per field, either a uniform scalar or a
        per-vertex ``[file, key]`` map; both are materialised here to full-mesh
        arrays (``elaMesh`` / ``iceMesh`` / ``termMesh``) which ``iceAccumulation``
        indexes by ``locIDs``. Maps are (re)loaded only when the active interval
        changes (step changes, like the precipitation maps).
        """
        ts = self._iceTimeSeries

        # Active interval: the latest event whose start time is at or before now.
        nb = 0
        for k in range(len(ts)):
            if ts[k]["start"] <= self.tNow:
                nb = k
            else:
                break
        if nb == self._iceSeriesIdx:
            return
        self._iceSeriesIdx = nb

        iv = ts[nb]
        self.elaMesh = self._resolveIceField(iv["hela"])
        self.iceMesh = self._resolveIceField(iv["hice"])
        self.termMesh = self._resolveIceField(iv["hterm"])

        return

    def _resolveIceField(self, field):
        """
        Materialise one glacier-geometry field (a ``(scalar, map_spec)`` pair,
        exactly one non-None) to a full-mesh array.
        """
        scalar, spec = field
        if spec is not None:
            return self._loadIceMap(spec, "ice map")
        return np.full(self.mpoints, scalar, dtype=np.float64)

    def _updateEvap(self):
        """
        Finds current evaporation values for the considered time interval
        and stores them as a per-node m/yr numpy array on `self.evapVal`.

        Mirrors `_updateRain` but only supports two source types:
        a `eUniform` scalar (m/yr) or an `eMap`+`eKey` pair pointing to a
        per-node array in an .npz file. Elevation-banded evaporation is
        intentionally not supported in v1 (see DESIGN_EVAPORATION.md D4).

        .. note::

            Evaporation is treated as a within-step sink at two points in
            the flow solver: (i) subtracted from rainfall before the IDA
            solve (channel runoff), and (ii) subtracted from per-pit inflow
            inside `_distributeDownstream` before the spillover decision
            (lake-surface evap). Lakes that cannot sustain themselves
            against the local evap budget simply do not form. The same
            cell-area assumption is used on land and over lake surfaces;
            users who need a stronger lake-surface rate should encode the
            spatial pattern directly in `eMap`.
        """

        nb = self.evapNb
        if nb < len(self.evapdata) - 1:
            if self.evapdata.at[nb + 1, "start"] <= self.tNow:
                nb += 1

        if nb > self.evapNb or nb == -1:
            if nb == -1:
                nb = 0

            self.evapNb = nb
            if pd.isnull(self.evapdata["eUni"][nb]):
                loadData = np.load(self.evapdata.at[nb, "eMap"])
                evapVal = loadData[self.evapdata.at[nb, "eKey"]]
                del loadData
            else:
                evapVal = np.full(self.mpoints, self.evapdata.at[nb, "eUni"])

            evapVal[evapVal < 0] = 0.0
            self.evapMesh = evapVal

        self.evapVal = self.evapMesh[self.locIDs]

        return

    def _updateEroFactor(self):
        """
        Finds current erodibility factor values for the considered time interval.

        .. note::

            It is worth noting that the erodibility factor is an indice representing different lithological classes (see `Moosdorf et al., 2018 <https://www.sciencedirect.com/science/article/pii/S0143622817306859>`_).

        """

        nb = self.sedfactNb
        if nb < len(self.sedfacdata) - 1:
            if self.sedfacdata.at[nb + 1, "start"] <= self.tNow:  # + self.dt:
                nb += 1

        if nb > self.sedfactNb or nb == -1:
            if nb == -1:
                nb = 0

            self.sedfactNb = nb
            if pd.isnull(self.sedfacdata["sUni"][nb]):
                loadData = np.load(self.sedfacdata.at[nb, "sMap"])
                sedfacVal = loadData[self.sedfacdata.at[nb, "sKey"]]
                del loadData
            else:
                sedfacVal = np.full(self.mpoints, self.sedfacdata.at[nb, "sUni"])
            sedfacVal[sedfacVal < 0.1] = 0.1
            self.sedFacMesh = sedfacVal

        self.sedfacVal = self.sedFacMesh[self.locIDs]

        return

    def destroy_DMPlex(self):
        """
        Destroys PETSc DMPlex objects and associated PETSc local/global Vectors and Matrices at the end of the simulation.
        """

        t0 = process_time()

        self.hLocal.destroy()
        self.hGlobal.destroy()
        self.hOldFlex.destroy()
        self.h.destroy()
        self.hl.destroy()
        self.dh.destroy()
        self.FAG.destroy()
        self.FAL.destroy()
        self.fillFAL.destroy()
        self.cumED.destroy()
        self.cumEDLocal.destroy()
        self.vSed.destroy()
        self.vSedLocal.destroy()
        self.vSedF.destroy()
        self.vSedFLocal.destroy()
        if getattr(self, "provOn", False):
            for v in self.vSedP:
                v.destroy()
            self.vSedPLocal.destroy()
        self.areaGlobal.destroy()
        self.bG.destroy()
        self.bL.destroy()
        self.hOld.destroy()
        self.hOldLocal.destroy()
        self.Qs.destroy()
        self.newH.destroy()
        self.tmpL.destroy()
        self.tmp.destroy()
        self.QsL.destroy()
        self.nQs.destroy()
        self.tmp1.destroy()
        self.stepED.destroy()
        self.Eb.destroy()
        self.EbLocal.destroy()
        if self.iceOn:
            self.iceHL.destroy()
            self.iceMeltL.destroy()
            self.iceUbL.destroy()
            self.iceAbrL.destroy()
            if self.flexOn:
                self.iceFlex.destroy()

        self.iMat.destroy()
        self.lgmap_col.destroy()
        self.lgmap_row.destroy()
        self.dm.destroy()
        self.zMat.destroy()
        self.mat.destroy()

        # Cached KSP/SNES/TS helpers (created lazily on first use)
        for name in (
            "_ksp_main", "_ksp_fallback",
            "_snes_ed", "_snes_ed_f", "_snes_ed_x",
            "_snes_nl", "_snes_nl_f", "_snes_nl_x", "_snes_nl_J",
            "_snes_soil", "_snes_soil_f", "_snes_soil_x",
            "_snes_hill", "_snes_hill_f", "_snes_hill_x",
            "_snes_ice", "_snes_ice_f", "_snes_ice_x",
            "_ts_marine", "_ts_marine_x",
            "_ts_soil", "_ts_soil_x", "_ts_soil_f",
        ):
            obj = getattr(self, name, None)
            if obj is not None:
                obj.destroy()

        petsc4py.PETSc.garbage_cleanup()

        del self.lcoords, self.lcells, self.inIDs

        del self.stratH, self.stratZ, self.phiS

        if not self.fast:
            del self.distRcv, self.wghtVal, self.rcvID

        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Cleaning Model Dataset (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        if self.showlog:
            self.log.view()

        if MPIrank == 0:
            print(
                "\n+++\n+++ Total run time (%0.02f seconds)\n+++"
                % (process_time() - self.modelRunTime),
                flush=True,
            )

        return
