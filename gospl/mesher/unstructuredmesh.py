import os
import gc
import sys
import vtk

import warnings
import meshplex
import petsc4py
import numpy as np
import pandas as pd

from mpi4py import MPI
from scipy import spatial
from time import process_time

from vtk.util import numpy_support

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import definetin
    from gospl._fortran import ngbglob

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class UnstMesh(object):
    """
    This class defines the spherical mesh characteristics and builds a PETSc DMPlex that encapsulates this unstructured mesh, with interfaces for both topology and geometry. The PETSc DMPlex is used for parallel redistribution for load balancing.

    .. note::

        `gospl` is built around a **Finite-Volume** method (FVM) for representing and evaluating  partial differential equations. It requires the definition of several mesh variables such as:

            - the number of neighbours surrounding every node,
            - the cell area defined using  Voronoi area,
            - the length of the edge connecting every nodes, and
            - the length of the Voronoi faces shared by each node with his neighbours.

    In addition to mesh defintions, the class declares several functions related to forcing conditions (*e.g.* paleo-precipitation maps, tectonic (vertical and horizontal) displacements, stratigraphic layers...). These functions are defined within the `UnstMesh` class as they rely heavely on the mesh structure.

    Finally a function to clean all PETSc variables is defined and called at the end of a simulation.
    """

    def __init__(self):
        """
        The initialisation of `UnstMesh` class calls the private function **_buildMesh**.
        """

        self.hdisp = None
        self.uplift = None
        self.rainVal = None
        self.memclear = True
        # Let us define the mesh variables and build PETSc DMPLEX.
        self._buildMesh()

        return

    def _meshfrom_cell_list(self, dim, cells, coords):
        """
        Creates a DMPlex from a list of cells and coordinates.

        .. note::

            As far as I am aware, PETSc DMPlex requires to be initialised on one processor before load balancing.

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
        Defines the mesh structure and the associated voronoi parameter used in the Finite Volume method.

        .. important::
            The mesh structure is built locally on a single partition of the global mesh.

        This function uses the `meshplex` `library <https://meshplex.readthedocs.io>`_ to compute from the list of coordinates and cells the volume of each voronoi and their respective characteristics.

        Once the voronoi definitions have been obtained a call to the fortran subroutine `definetin` is performed to order each node and the dual mesh components, it records:

        - all cells surrounding a given vertice,
        - all edges connected to a given vertice,
        - the triangulation edge lengths,
        - the voronoi edge lengths.

        """

        # Create mesh structure with meshplex
        t0 = process_time()
        Tmesh = meshplex.MeshTri(self.mCoords, self.mCells)
        self.marea = np.abs(Tmesh.control_volumes)
        self.marea[np.isnan(self.marea)] = 1.0
        self.larea = self.marea[self.locIDs]
        self.garea = self.marea[self.glbIDs]

        # Voronoi and simplices declaration
        Tmesh.create_edges()
        cc = Tmesh.cell_circumcenters
        if meshplex.__version__ >= "0.16.0":
            edges_nodes = Tmesh.edges["points"]
            cells_nodes = Tmesh.cells("points")
            cells_edges = Tmesh.cells("edges")
        elif meshplex.__version__ >= "0.14.0":
            edges_nodes = Tmesh.edges["points"]
            cells_nodes = Tmesh.cells["points"]
            cells_edges = Tmesh.cells["edges"]
        else:
            edges_nodes = Tmesh.edges["nodes"]
            cells_nodes = Tmesh.cells["nodes"]
            cells_edges = Tmesh.cells["edges"]

        # Finite volume discretisation
        self.FVmesh_ngbID, self.edgeMax = definetin(
            self.lpoints,
            self.mCoords,
            self.lgIDs,
            cells_nodes,
            cells_edges,
            edges_nodes,
            self.marea,
            cc.T,
        )

        del Tmesh, edges_nodes, cells_nodes, cells_edges, cc
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "FV discretisation (%0.02f seconds)" % (process_time() - t0), flush=True
            )

        return

    def reInitialiseElev(self):
        """
        This functions reinitialises `gospl` elevation if one wants to restart a simulation that does not necessitate to rebuild the mesh structure. It can be used when running a simulation from the Jupyter environment.

        .. note:
            The elevation is read directly from the mesh elevation file defined in the YAML input file from the key: **npdata**.

        """

        t0step = process_time()

        # Restart time
        self.tNow = self.tStart
        self.step = 0
        self.stratStep = 0
        self.rStart = self.tStart
        self.saveTime = self.tNow
        if self.strat > 0:
            self.saveStrat = self.tNow + self.strat
        else:
            self.saveStrat = self.tEnd + 2.0 * self.tout

        # Forcing functions
        self.rainNb = -1
        self.tecNb = -1
        self.flexNb = -1

        # Getting PETSc vectors values
        loadData = np.load(self.meshFile)
        gZ = loadData["z"]
        self.hLocal.setArray(gZ[self.locIDs])
        self.hGlobal.setArray(gZ[self.glbIDs])
        # self.dm.localToGlobal(self.hLocal, self.hGlobal)
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

        # Update external forces
        self.applyForces()
        self.applyTectonics()

        del gZ, loadData
        gc.collect()

        if MPIrank == 0:
            print(
                "--- Reinitialise Phase \
                  (%0.02f seconds)\n+++"
                % (process_time() - t0step),
                flush=True,
            )

        return

    def _xyz2lonlat(self):
        """
        Converts local x,y,z representation of cartesian coordinates from the spherical triangulation to latitude, longitude in degrees.

        The latitudinal and longitudinal extend are between [0,180] and [0,360] respectively. Lon/lat coordinates are used when forcing the simulation with paleo-topoography maps.
        """

        self.lLatLon = np.zeros((self.lpoints, 2))
        self.lLatLon[:, 0] = np.arcsin(self.lcoords[:, 2] / self.radius)
        self.lLatLon[:, 1] = np.arctan2(self.lcoords[:, 1], self.lcoords[:, 0])
        self.lLatLon[:, 0] = np.mod(np.degrees(self.lLatLon[:, 0]) + 90, 180.0)
        self.lLatLon[:, 1] = np.mod(np.degrees(self.lLatLon[:, 1]) + 180.0, 360.0)

        return

    def _generateVTKmesh(self, points, cells):
        """
        A global VTK mesh is generated to compute the distance between mesh vertices and coastlines position.

        The distance to the coastline for every marine vertices is used to define a maximum shelf slope during deposition. The coastline contours are efficiently obtained from VTK contouring function. This function is performed on a VTK mesh which is built in this function.
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
        This function is at the core of the `UnstMesh` class. It encapsulates both spherical mesh construction (triangulation and voronoi representation for the Finite Volume discretisation), PETSc DMPlex distribution and several PETSc vectors allocation.

        The function relies on several private functions from the class:

        - _generateVTKmesh
        - _meshfrom_cell_list
        - _meshStructure
        - _readErosionDeposition
        - _xyz2lonlat
        - readStratLayers

        .. note::

            It is worth mentionning that partitioning and field distribution from global to local PETSc DMPlex takes a lot of time for large mesh and there might be some quicker way in PETSc to perform this step that I am unaware of...

        """

        # Read mesh attributes from file
        t0 = process_time()
        loadData = np.load(self.meshFile)
        self.mCoords = loadData["v"]
        self.mpoints = len(self.mCoords)
        gZ = loadData["z"]
        self.mCells = loadData["c"].astype(int)
        self.vtkMesh = None
        self.flatModel = False
        if self.mCoords[:, 2].max() == 0.0:
            self.flatModel = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self._generateVTKmesh(self.mCoords, self.mCells)

        # Store global neighbouring on process rank 0
        if MPIrank == 0:
            ngbglob(self.mpoints, loadData["n"])
        if MPIrank == 0 and self.verbose:
            print(
                "Reading mesh information (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Create DMPlex
        self._meshfrom_cell_list(2, loadData["c"], self.mCoords)
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
        self.lcells = np.zeros((cEnd - cStart, 3), dtype=petsc4py.PETSc.IntType)
        point_closure = None
        for c in range(cStart, cEnd):
            point_closure = self.dm.getTransitiveClosure(c)[0]
            self.lcells[c, :] = point_closure[-3:] - cEnd
        if point_closure is not None:
            del point_closure
        gc.collect()
        if MPIrank == 0 and self.verbose:
            print(
                "Defining local DMPlex (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # From mesh values to local and global ones...
        t0 = process_time()
        tree = spatial.cKDTree(self.mCoords, leafsize=10)
        distances, self.locIDs = tree.query(self.lcoords, k=1)
        distances, self.glbIDs = tree.query(self.gcoords, k=1)

        # From local to global values...
        self.lgIDs = -np.ones(self.mpoints, dtype=int)
        self.lgIDs[self.locIDs] = np.arange(self.lpoints)
        self.outIDs = np.where(self.lgIDs < 0)[0]

        # Local/Global mapping
        self.lgmap_row = self.dm.getLGMap()
        l2g = self.lgmap_row.indices.copy()
        offproc = l2g < 0
        l2g[offproc] = -(l2g[offproc] + 1)
        self.lgmap_col = petsc4py.PETSc.LGMap().create(
            l2g, comm=petsc4py.PETSc.COMM_WORLD
        )

        # Vertex part of an unique partition
        vIS = self.dm.getVertexNumbering()
        self.inIDs = np.zeros(self.lpoints, dtype=int)
        self.inIDs[vIS.indices >= 0] = 1
        self.lIDs = np.where(self.inIDs == 1)[0]

        if self.flatModel:
            # Local mesh boundary points in 2D model
            localBound = self._get_boundary()
            idLocal = np.where(vIS.indices >= 0)[0]
            self.idBorders = np.where(np.isin(idLocal, localBound))[0]
            self.glBorders = np.zeros(self.mpoints, dtype=int)
            idLocal = self.glBorders.copy()
            localBound = np.zeros(self.lpoints, dtype=int)
            localBound[self.idBorders] = 1
            idLocal[self.locIDs] = localBound
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, idLocal, op=MPI.MAX)
            self.glBorders = np.where(idLocal == 1)[0]
            del idLocal, localBound
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

        # Create mesh structure with meshplex
        self._meshStructure()

        # Find global shadow nodes indices
        self.shadow_lOuts = np.where(self.inIDs == 0)[0]
        distances, indices = tree.query(self.lcoords[self.shadow_lOuts, :], k=1)
        tmp = -np.ones(self.mpoints, dtype=int)
        tmp[indices] = 1.0
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, tmp, op=MPI.MAX)
        self.shadowAlls = np.unique(np.where(tmp > 0)[0])

        # Find local shadow nodes global indices
        tree = spatial.cKDTree(self.mCoords[self.shadowAlls, :], leafsize=10)
        distances, indices = tree.query(self.lcoords, k=1)
        if MPIsize > 1:
            self.gshadowIDs = np.unique(self.shadowAlls[indices])
            tmp.fill(0.0)
            tmp[self.gshadowIDs] = 1.0
            ltmp = tmp[self.locIDs]
            self.lshadowIDs = np.unique(np.where(ltmp > 0)[0])
            distances, indices = tree.query(self.lcoords[self.shadow_lOuts], k=1)
            self.shadow_gOuts = np.unique(self.shadowAlls[indices])

            # Get shadow zones sorted indices
            self.shadowgNb = len(self.shadowAlls)
            self.gshadinIDs = np.where(np.in1d(self.shadowAlls, self.gshadowIDs))[0]
            self.gshadoutIDs = np.where(np.in1d(self.shadowAlls, self.shadow_gOuts))[0]
            tmp = np.concatenate((self.shadow_lOuts, np.arange(self.lpoints)))
            cbin = np.bincount(tmp)
            self.lshadinIDs = np.where(cbin == 1)[0]
            del ltmp, cbin

        # Define cumulative erosion deposition arrays
        self._readErosionDeposition()

        self.sealevel = self.seafunction(self.tNow)
        self.areaGlobal = self.hGlobal.duplicate()
        self.areaGlobal.setArray(self.garea)
        self.minArea = self.areaGlobal.min()

        # Forcing event number
        self.bG = self.hGlobal.duplicate()
        self.bL = self.hLocal.duplicate()

        self.rainNb = -1
        self.tecNb = -1
        self.flexNb = -1

        del tree, distances, indices, tmp
        del l2g, offproc, gZ
        gc.collect()

        # Map longitude/latitude coordinates
        self._xyz2lonlat()

        # Build stratigraphic data if any
        if self.stratNb > 0:
            self.readStratLayers()

        return

    def _set_DMPlex_boundary_points(self, label):
        """
        In case of a flat mesh (non global), this function finds the points that join the edges that have been marked as "boundary" faces in the DAG then sets them as boundaries.
        """

        self.dm.createLabel(label)
        self.dm.markBoundaryFaces(label)

        pStart, pEnd = self.dm.getDepthStratum(0)  # points
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
            del vertices

        edgeIS.destroy()

        return

    def _get_boundary(self, label="boundary"):
        """
        In case of a flat mesh (non global), this function finds the nodes on the boundary from the DM.
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
            indices = np.zeros((0,), dtype=np.int)

        return indices

    def _readErosionDeposition(self):
        """
        Reads existing cumulative erosion depostion from a previous experiment if any as defined in the YAML input file following the  `nperodep` key.

        This functionality can be used when restarting from a previous simulation in which the sperical mesh has been modified either to account for horizontal advection or to refine/coarsen a specific region during a given time period.
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
            # self.dm.localToGlobal(self.cumEDLocal, self.cumED)
            del gED
            gc.collect()

        return

    def initExtForce(self):
        """
        Initialises the forcing conditions: sea level condition, rainfall maps and tectonics.
        """

        self.applyForces()
        if self.backward:
            self.applyTectonics()

        return

    def applyForces(self):
        """
        Finds the different values for climatic and sea-level forcing that will be applied at any given time interval during the simulation.
        """

        t0 = process_time()
        # Sea level
        if self.tNow == self.tStart:
            self.sealevel = self.seafunction(self.tNow)
        else:
            self.sealevel = self.seafunction(self.tNow + self.dt)

        # Climate information
        self._updateRain()

        if MPIrank == 0 and self.verbose:
            print(
                "Update Climatic Forces (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def applyTectonics(self):
        """
        Finds the different values for tectonic forces that will be applied at any given time interval during the simulation.
        """

        t0 = process_time()
        self._updateTectonics()

        if MPIrank == 0 and self.verbose:
            print(
                "Update Tectonic Forces (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def updatePaleomap(self):
        """
        Forces model to match provided paleomaps at specific time interval during the simulation.

        This function is only used when the `paleomap` key is defined in the YAML input file. It will read the paleotopography map and assign the paleo elevation on the spherical mesh.
        """

        for k in range(self.paleoNb):
            if self.tNow == self.paleodata.iloc[k, 0]:
                loadData = np.load(self.paleodata.iloc[k, 1])
                gZ = loadData["z"]
                self.hLocal.setArray(gZ[self.locIDs])
                self.hGlobal.setArray(gZ[self.glbIDs])
                del loadData, gZ
                gc.collect()

        return

    def _updateRain(self):
        """
        Finds current rain values for the considered time interval and computes the **volume** of water available for runoff on each vertex.

        .. note::

            It is worth noting that the precipitation maps are considered as runoff water. If one wants to account for evaporation and infiltration you will need to modify the precipitation maps accordingly as a pre-processing step.

        """

        nb = self.rainNb
        if nb < len(self.raindata) - 1:
            if self.raindata.iloc[nb + 1, 0] <= self.tNow:  # + self.dt:
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
                rainArea = np.full(self.mpoints, self.raindata.iloc[nb, 1])

            self.rainArea = rainArea * self.marea

        self.rainVal = self.rainArea[self.locIDs] / self.marea[self.locIDs]
        tmpZ = self.hLocal.getArray()
        rainArea = self.rainArea[self.locIDs]
        rainArea[tmpZ < self.sealevel] = 0.0
        self.bL.setArray(rainArea)

        tmpZ = self.hGlobal.getArray()
        rainArea = self.rainArea[self.glbIDs]
        rainArea[tmpZ < self.sealevel] = 0.0
        self.bG.setArray(rainArea)

        if self.memclear:
            del rainArea, tmpZ
            gc.collect()

        return

    def _updateTectonics(self):
        """
        Finds the current tectonic regimes (horizontal and vertical) for the considered time interval.

        For horizontal displacements, the mesh variables will have to be first advected over the grid and then reinterpolated on the initial mesh coordinates. The approach here does not allow for mesh refinement in zones of convergence and thus can be limiting but using a fixed mesh has one main advantage: the mesh and Finite Volume discretisation do not have to be rebuilt each time the mesh is advected.

        This function calls the following 2 private functions:

        - _meshAdvectorSphere
        - _meshUpliftSubsidence

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
                self.hdisp = mdata["xyz"][self.locIDs, :]
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

        if self.forceStep >= 0 and not self.newForcing:
            self._meshUpliftSubsidence(None)

        return

    def _meshUpliftSubsidence(self, tectonic):
        """
        Applies vertical displacements based on tectonic rates.

        :arg tectonic: local tectonic rates
        """

        # Define vertical displacements
        tmp = self.hLocal.getArray().copy()
        if tectonic is not None:
            self.uplift = tectonic[self.locIDs]
        self.hLocal.setArray(tmp + self.uplift * self.dt)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)
        del tmp
        gc.collect()

        return

    def _meshAdvectorSphere(self, tectonic, timer):
        """
        Advects the spherical mesh horizontally and interpolates mesh information.

        The advection proceeds in each partition seprately in the following way:

        1. based on the horizontal displacement velocities, the mesh coordinates and associated variables (cumulative erosion deposition and stratigraphic layers composition) are moved.
        2. a kdtree is built with the advected coordinates and used to interpolate the mesh variables on the initial local mesh position. The interpolation is based on a weighting distance function accounting for the 3 closest advected vertices.
        3. interpolated variables on the initial mesh coordinates are then stored in PETSc vectors and class parameters are updated accordingly.


        :arg tectonic: local tectonic rates in 3D
        :arg timer: tectonic time step in years
        """

        t1 = process_time()

        # Move local coordinates
        XYZ = self.lcoords + tectonic[self.locIDs, :] * timer

        # Update local elevation
        tmp = self.hLocal.getArray().copy()
        elev = tmp[self.lgIDs]
        if MPIsize > 1:
            temp = np.full(self.shadowgNb, -1.0e8, dtype=np.float64)
            temp[self.gshadinIDs] = elev[self.gshadowIDs]
            temp[self.gshadoutIDs] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            elev[self.shadowAlls] = temp
            loc_elev = elev[self.locIDs]
        else:
            loc_elev = elev[self.locIDs]

        # Update local erosion/deposition
        tmp = self.cumEDLocal.getArray().copy()
        erodep = tmp[self.lgIDs]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs] = erodep[self.gshadowIDs]
            temp[self.gshadoutIDs] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            erodep[self.shadowAlls] = temp
            loc_erodep = erodep[self.locIDs]
        else:
            loc_erodep = erodep[self.locIDs]

        # Update local stratal dataset after displacements
        if self.stratNb > 0 and self.stratStep > 0:
            if self.carbOn:
                (
                    loc_stratH,
                    loc_stratZ,
                    loc_stratF,
                    loc_stratW,
                    loc_stratC,
                    loc_phiS,
                    loc_phiF,
                    loc_phiW,
                    loc_phiC,
                ) = self.updateDispStrata()
            else:
                (
                    loc_stratH,
                    loc_stratZ,
                    loc_stratF,
                    loc_stratW,
                    loc_phiS,
                    loc_phiF,
                    loc_phiW,
                ) = self.updateDispStrata()

        # Build and query local kd-tree
        tree = spatial.cKDTree(XYZ, leafsize=10)
        distances, indices = tree.query(self.lcoords, k=3)

        # Inverse weighting distance...
        if self.interp == 1:
            nelev = loc_elev[indices[:, 0]]

        weights = np.divide(
            1.0, distances, out=np.zeros_like(distances), where=distances != 0
        )
        onIDs = np.where(distances[:, 0] == 0)[0]
        temp = np.sum(weights, axis=1)

        # Update elevation
        if self.interp > 1:
            tmp = np.sum(weights * loc_elev[indices], axis=1)
            nelev = np.divide(tmp, temp, out=np.zeros_like(temp), where=temp != 0)

        # Update erosion deposition
        tmp = np.sum(weights * loc_erodep[indices], axis=1)
        nerodep = np.divide(tmp, temp, out=np.zeros_like(temp), where=temp != 0)

        # Update stratigraphic record
        if self.stratNb > 0 and self.stratStep > 0:
            if self.carbOn:
                self.stratalRecord(
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
            else:
                self.stratalRecord(
                    indices,
                    weights,
                    loc_stratH,
                    loc_stratZ,
                    loc_stratF,
                    loc_stratW,
                    None,
                    loc_phiS,
                    loc_phiF,
                    loc_phiW,
                    None,
                )

        if len(onIDs) > 0:
            nerodep[onIDs] = loc_erodep[indices[onIDs, 0]]
            if self.interp > 1:
                nelev[onIDs] = loc_elev[indices[onIDs, 0]]
            if self.stratNb > 0 and self.stratStep > 0:
                self.stratZ[onIDs, : self.stratStep] = loc_stratZ[indices[onIDs, 0], :]
                self.stratH[onIDs, : self.stratStep] = loc_stratH[indices[onIDs, 0], :]
                self.phiS[onIDs, : self.stratStep] = loc_phiS[indices[onIDs, 0], :]
                if self.stratF is not None:
                    self.stratF[onIDs, : self.stratStep] = loc_stratF[
                        indices[onIDs, 0], :
                    ]
                    self.phiF[onIDs, : self.stratStep] = loc_phiW[indices[onIDs, 0], :]
                if self.stratW is not None:
                    self.stratW[onIDs, : self.stratStep] = loc_stratW[
                        indices[onIDs, 0], :
                    ]
                    self.phiW[onIDs, : self.stratStep] = loc_phiF[indices[onIDs, 0], :]
                if self.carbOn:
                    self.stratC[onIDs, : self.stratStep] = loc_stratC[
                        indices[onIDs, 0], :
                    ]
                    self.phiC[onIDs, : self.stratStep] = loc_phiC[indices[onIDs, 0], :]

        self.hLocal.setArray(nelev)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        self.cumEDLocal.setArray(nerodep)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED)

        # Update local stratal dataset
        if self.stratNb > 0 and self.stratStep > 0:
            self.localStrat()
            del loc_stratH, loc_stratZ, loc_stratF, loc_stratW
            del loc_phiS, loc_phiF, loc_phiW
            if self.carbOn:
                del loc_stratC, loc_phiC

        if MPIrank == 0 and self.verbose:
            print(
                "Advect Mesh Information Horizontally (%0.02f seconds)"
                % (process_time() - t1),
                flush=True,
            )

        del XYZ, elev, nelev, erodep, nerodep
        del loc_elev, loc_erodep, onIDs
        del weights, tmp, temp
        del tree, distances, indices
        gc.collect()

        return

    def destroy_DMPlex(self):
        """
        Destroys PETSc DMPlex objects and associated PETSc local/global Vectors and Matrices at the end of the simulation.
        """

        t0 = process_time()

        self.hLocal.destroy()
        self.hGlobal.destroy()
        self.FAG.destroy()
        self.FAL.destroy()
        self.fillFAL.destroy()
        self.cumED.destroy()
        self.cumEDLocal.destroy()
        self.vSed.destroy()
        self.vSedLocal.destroy()
        if self.stratNb > 0:
            if self.stratF is not None:
                self.vSedf.destroy()
                self.vSedfLocal.destroy()
            if self.stratW is not None:
                self.vSedw.destroy()
                self.vSedwLocal.destroy()
            if self.carbOn:
                self.vSedc.destroy()
                self.vSedcLocal.destroy()
        self.areaGlobal.destroy()
        self.bG.destroy()
        self.bL.destroy()
        self.hOld.destroy()
        self.hOldLocal.destroy()
        self.Qs.destroy()
        self.tmpL.destroy()
        self.tmp.destroy()
        self.tmp1.destroy()
        self.stepED.destroy()
        self.Eb.destroy()
        self.EbLocal.destroy()
        self.upsG.destroy()
        self.upsL.destroy()

        self.iMat.destroy()
        if not self.fast:
            self.fMat.destroy()
        self.lgmap_col.destroy()
        self.lgmap_row.destroy()
        self.dm.destroy()

        del self.lcoords, self.lcells, self.inIDs

        del self.stratH, self.stratZ, self.phiS

        if self.stratF is not None:
            del self.stratF, self.phiF

        if self.stratW is not None:
            del self.stratW, self.phiW

        if self.carbOn:
            del self.stratC, self.phiC

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
