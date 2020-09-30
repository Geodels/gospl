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
    from gospl._fortran import defineTIN
    from gospl._fortran import ngbGlob
    from gospl._fortran import strataBuild
    from gospl._fortran import strataBuildCarb

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class UnstMesh(object):
    """
    This class defines the spherical mesh characteristics and builds a PETSc DMPlex that
    encapsulates this unstructured mesh, with interfaces for both topology and geometry.
    The PETSc DMPlex is used for parallel redistribution for load balancing.

    .. note::

        `gospl` is built around a **Finite-Volume** method (FVM) for representing and evaluating
        partial differential equations. It requires the definition of several mesh variables such as:

            - the number of neighbours surrounding every node,
            - the cell area defined using  Voronoi area,
            - the length of the edge connecting every nodes, and
            - the length of the Voronoi faces shared by each node with his neighbours.

    In addition to mesh defintions, the class declares several functions related to forcing conditions
    (*e.g.* paleo-precipitation maps, tectonic (vertical and horizontal) displacements, stratigraphic
    layers...). These functions are defined within the `UnstMesh` class as they rely heavely on the mesh
    structure.

    Finally a function to clean all PETSc variables is defined and called at the end of a simulation.
    """

    def __init__(self):
        """
        The initialisation of `UnstMesh` class calls the private function **_buildMesh**.
        """

        self.hdisp = None
        self.uplift = None
        self.rainVal = None
        self.stratH = None
        self.stratF = None
        self.stratZ = None
        self.stratC = None
        self.phiS = None
        self.phiF = None
        self.phiC = None

        # Let us define the mesh variables and build PETSc DMPLEX.
        self._buildMesh()

        return

    def _meshfrom_cell_list(self, dim, cells, coords):
        """
        Creates a DMPlex from a list of cells and coordinates.

        .. note::

            As far as I am aware, PETSc DMPlex requires to be initialised on one
            processor before load balancing.

        :arg dim: topological dimension of the mesh
        :arg cells: vertices of each cell
        :arg coords: coordinates of each vertex
        """

        if MPIrank == 0:
            cells = np.asarray(cells, dtype=np.int32)
            coords = np.asarray(coords, dtype=np.double)
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
                np.zeros(coord_shape, dtype=np.double),
                comm=petsc4py.PETSc.COMM_WORLD,
            )
        return

    def _meshStructure(self):
        """
        Defines the mesh structure and the associated voronoi parameter used in the
        Finite Volume method.

        .. important::
            The mesh structure is built locally on a single partition of the global mesh.

        This function uses the `meshplex` `library <https://meshplex.readthedocs.io>`_ to compute
        from the list of coordinates and cells the volume of each voronoi and their respective
        characteristics.

        Once the voronoi definitions have been obtained a call to the fortran subroutine `defineTIN`
        is performed to order each node and the dual mesh components, it records:

        - all cells surrounding a given vertice,
        - all edges connected to a given vertice,
        - the triangulation edge lengths,
        - the voronoi edge lengths.

        """

        # Create mesh structure with meshplex
        t0 = process_time()
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
        edgeMax = np.zeros(1, dtype=np.float64)
        self.FVmesh_ngbID, edgeMax[0] = defineTIN(
            self.lcoords, cells_nodes, cells_edges, edges_nodes, self.area, cc.T
        )
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, edgeMax, op=MPI.MAX)
        self.edgeMax = edgeMax[0]

        self.garea = np.zeros(self.gpoints)
        self.garea = self.area[self.lgIDs]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, self.garea, op=MPI.MAX)

        del Tmesh, edges_nodes, cells_nodes, cells_edges, cc
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "FV discretisation (%0.02f seconds)" % (process_time() - t0), flush=True
            )

        return

    def reInitialiseElev(self):
        """
        This functions reinitialises `gospl` elevation if one wants to restart a simulation that does
        not necessitate to rebuild the mesh structure. It can be used when running a simulation
        from the Jupyter environment.

        .. note:
            The elevation is read directly from the mesh elevation file defined in the YAML
            input file from the key: **npdata**.

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
        self.hLocal.setArray(gZ[self.glIDs])
        self.dm.localToGlobal(self.hLocal, self.hGlobal)
        self.vSed.set(0.0)
        self.vSedLocal.set(0.0)
        if self.stratNb > 0:
            self.vSedf.set(0.0)
            self.vSedfLocal.set(0.0)
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
        Converts local x,y,z representation of cartesian coordinates from the spherical
        triangulation to latitude, longitude in degrees.

        The latitudinal and longitudinal extend are between [0,180] and [0,360] respectively.
        Lon/lat coordinates are used when forcing the simulation with paleo-topoography maps.
        """

        # self.gLatLon = np.zeros((self.gpoints, 2))
        # self.gLatLon[:, 0] = np.arcsin(self.gCoords[:, 2] / self.radius)
        # self.gLatLon[:, 1] = np.arctan2(self.gCoords[:, 1], self.gCoords[:, 0])
        # self.gLatLon[:, 0] = np.mod(np.degrees(self.gLatLon[:, 0]) + 90, 180.0)
        # self.gLatLon[:, 1] = np.mod(np.degrees(self.gLatLon[:, 1]) + 180.0, 360.0)

        self.lLatLon = np.zeros((self.npoints, 2))
        self.lLatLon[:, 0] = np.arcsin(self.lcoords[:, 2] / self.radius)
        self.lLatLon[:, 1] = np.arctan2(self.lcoords[:, 1], self.lcoords[:, 0])
        self.lLatLon[:, 0] = np.mod(np.degrees(self.lLatLon[:, 0]) + 90, 180.0)
        self.lLatLon[:, 1] = np.mod(np.degrees(self.lLatLon[:, 1]) + 180.0, 360.0)

        return

    def _generateVTKmesh(self, points, cells):
        """
        A global VTK mesh is generated when the `shelfslope` key is set to True in the YAML
        input file.

        In this case, the distance to the coastline for every marine vertices is used to define
        a maximum shelf slope during deposition. The coastline contours are efficiently obtained
        from VTK contouring function. This function is performed on a VTK mesh which is built in
        here.
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
        This function is at the core of the `UnstMesh` class. It encapsulates both spherical mesh
        construction (triangulation and voronoi representation for the Finite Volume discretisation),
        PETSc DMPlex distribution and several PETSc vectors allocation.

        The function relies on several private functions from the class:

        - _generateVTKmesh
        - _meshfrom_cell_list
        - _meshStructure
        - _readErosionDeposition
        - _xyz2lonlat
        - _readStratLayers

        .. note::

            It is worth mentionning that partitioning and field distribution from global to local
            PETSc DMPlex takes a lot of time for large mesh and there might be some quicker way
            in PETSc to perform this step that I am unaware of...

        """

        # Read mesh attributes from file
        t0 = process_time()
        loadData = np.load(self.meshFile)
        self.gCoords = loadData["v"]
        self.gpoints = len(self.gCoords)
        gZ = loadData["z"]

        self.vtkMesh = None
        if self.shelfslope:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self._generateVTKmesh(self.gCoords, loadData["c"])

        # Store global neighbouring on process rank 0
        if MPIrank == 0:
            ngbGlob(self.gpoints, loadData["n"])
        if MPIrank == 0 and self.verbose:
            print(
                "Reading mesh information (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Create DMPlex
        self._meshfrom_cell_list(2, loadData["c"], self.gCoords)
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
        self.lcoords = self.dm.getCoordinatesLocal().array.reshape(-1, 3)
        self.npoints = self.lcoords.shape[0]
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

        # From global values to local ones...
        t0 = process_time()
        tree = spatial.cKDTree(self.gCoords, leafsize=10)
        distances, self.glIDs = tree.query(self.lcoords, k=1)
        nelev = gZ[self.glIDs]

        # From local to global values...
        self.lgIDs = -np.ones(self.gpoints, dtype=int)
        self.lgIDs[self.glIDs] = np.arange(self.npoints)
        self.outIDs = np.where(self.lgIDs < 0)[0]
        self.lgIDs[self.lgIDs < 0] = 0

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
        self.inIDs = np.zeros(self.npoints, dtype=int)
        self.inIDs[vIS.indices >= 0] = 1
        self.lIDs = np.where(self.inIDs == 1)[0]
        vIS.destroy()

        # Local/Global vectors
        self.hGlobal = self.dm.createGlobalVector()
        self.hLocal = self.dm.createLocalVector()
        self.sizes = self.hGlobal.getSizes(), self.hGlobal.getSizes()

        self.hLocal.setArray(nelev)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

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
        tmp = -np.ones(self.gpoints, dtype=int)
        tmp[indices] = 1.0
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, tmp, op=MPI.MAX)
        self.shadowAlls = np.unique(np.where(tmp > 0)[0])

        # Find local shadow nodes global indices
        tree = spatial.cKDTree(self.gCoords[self.shadowAlls, :], leafsize=10)
        distances, indices = tree.query(self.lcoords, k=1)
        if MPIsize > 1:
            self.gshadowIDs = np.unique(self.shadowAlls[indices])
            tmp.fill(0.0)
            tmp[self.gshadowIDs] = 1.0
            ltmp = tmp[self.glIDs]
            self.lshadowIDs = np.unique(np.where(ltmp > 0)[0])
            distances, indices = tree.query(self.lcoords[self.shadow_lOuts], k=1)
            self.shadow_gOuts = np.unique(self.shadowAlls[indices])

            # Get shadow zones sorted indices
            self.shadowgNb = len(self.shadowAlls)
            self.gshadinIDs = np.where(np.in1d(self.shadowAlls, self.gshadowIDs))[0]
            self.gshadoutIDs = np.where(np.in1d(self.shadowAlls, self.shadow_gOuts))[0]
            tmp = np.concatenate((self.shadow_lOuts, np.arange(self.npoints)))
            cbin = np.bincount(tmp)
            self.lshadinIDs = np.where(cbin == 1)[0]
            del ltmp, cbin

        # Define cumulative erosion deposition arrays
        self._readErosionDeposition()

        self.sealevel = self.seafunction(self.tNow)

        areaLocal = self.hLocal.duplicate()
        self.areaGlobal = self.hGlobal.duplicate()
        areaLocal.setArray(self.area)
        self.dm.localToGlobal(areaLocal, self.areaGlobal)
        areaLocal.destroy()
        self.minArea = self.areaGlobal.min()

        # Forcing event number
        self.bG = self.hGlobal.duplicate()
        self.bL = self.hLocal.duplicate()

        self.rainNb = -1
        self.tecNb = -1
        self.flexNb = -1

        del tree, distances, indices, tmp
        del l2g, offproc, nelev, gZ
        gc.collect()

        # Map longitude/latitude coordinates
        self._xyz2lonlat()

        # Build stratigraphic data if any
        if self.stratNb > 0:
            self._readStratLayers()

        return

    def _readErosionDeposition(self):
        """
        Reads existing cumulative erosion depostion from a previous experiment if any
        as defined in the YAML input file following the  `nperodep` key.

        This functionality can be used when restarting from a previous simulation in which
        the sperical mesh has been modified either to account for horizontal advection or
        to refine/coarsen a specific region during a given time period.
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
            nED = gED[self.glIDs]

            self.cumEDLocal.setArray(nED)
            self.dm.localToGlobal(self.cumEDLocal, self.cumED)
            del nED, gED
            gc.collect()

        return

    def _readStratLayers(self):
        """
        When stratigraphic layers are turned on, this function reads any initial stratigraphic
        layers provided by within the YAML input file (key: `npstrata`).

        The following variables will be read from the file:

        - thickness of each stratigrapic layer `strataH` accounting for both
          erosion & deposition events.
        - proportion of fine sediment `strataF` contains in each stratigraphic layer.
        - elevation at time of deposition, considered to be to the current elevation
          for the top stratigraphic layer `strataZ`.
        - porosity of coarse sediment `phiS` in each stratigraphic layer computed at
          center of each layer.
        - porosity of fine sediment `phiF` in each stratigraphic layer computed at
          center of each layer.
        - proportion of carbonate sediment `strataC` contains in each stratigraphic layer
          if the carbonate module is turned on.
        - porosity of carbonate sediment `phiC` in each stratigraphic layer computed at
          center of each layer when the carbonate module is turned on.

        """

        if self.strataFile is not None:
            fileData = np.load(self.strataFile)
            stratVal = fileData["strataH"]
            self.initLay = stratVal.shape[1]
            self.stratNb += self.initLay

            # Create stratigraphic arrays
            self.stratH = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
            self.stratH[:, 0 : self.initLay] = stratVal[self.glIDs, 0 : self.initLay]

            stratVal = fileData["strataF"]
            self.stratF = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
            self.stratF[:, 0 : self.initLay] = stratVal[self.glIDs, 0 : self.initLay]

            stratVal = fileData["strataZ"]
            self.stratZ = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
            self.stratZ[:, 0 : self.initLay] = stratVal[self.glIDs, 0 : self.initLay]

            stratVal = fileData["phiS"]
            self.phiS = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
            self.phiS[:, 0 : self.initLay] = stratVal[self.glIDs, 0 : self.initLay]

            stratVal = fileData["phiF"]
            self.phiF = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
            self.phiF[:, 0 : self.initLay] = stratVal[self.glIDs, 0 : self.initLay]

            if self.carbOn:
                stratVal = fileData["strataC"]
                self.stratC = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
                self.stratC[:, 0 : self.initLay] = stratVal[
                    self.glIDs, 0 : self.initLay
                ]

                stratVal = fileData["phiC"]
                self.phiC = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
                self.phiC[:, 0 : self.initLay] = stratVal[self.glIDs, 0 : self.initLay]

            del fileData, stratVal
            gc.collect()
        else:
            self.stratH = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
            self.stratF = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
            self.stratZ = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
            self.phiS = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
            self.phiF = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
            if self.carbOn:
                self.stratC = np.zeros((self.npoints, self.stratNb), dtype=np.float64)
                self.phiC = np.zeros((self.npoints, self.stratNb), dtype=np.float64)

        return

    def initExtForce(self):
        """
        Initialises the forcing conditions: sea level condition, rainfall maps and
        tectonics.
        """

        self.applyForces()
        if self.backward:
            self.applyTectonics()

        return

    def applyForces(self):
        """
        Finds the different values for climatic and sea-level forcing that will
        be applied at any given time interval during the simulation.
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
        Finds the different values for tectonic forces that will be applied at any given time
        interval during the simulation.
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

        This function is only used when the `paleomap` key is defined in the YAML input file. It will
        read the paleotopography map and assign the paleo elevation on the spherical mesh.
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
        Finds current rain values for the considered time interval and computes the **volume** of water
        available for runoff on each vertex.

        .. note::

            It is worth noting that the precipitation maps are considered as runoff water. If one wants to
            account for evaporation and infiltration you will need to modify the precipitation maps
            accordingly as a pre-processing step.

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

    def _updateTectonics(self):
        """
        Finds the current tectonic regimes (horizontal and vertical) for the considered time interval.

        For horizontal displacements, the mesh variables will have to be first advected over the grid and then reinterpolated on the initial mesh coordinates. The approach here does not allow for mesh
        refinement in zones of convergence and thus can be limiting but using a fixed mesh has one main
        advantage: the mesh and Finite Volume discretisation do not have to be rebuilt each time the mesh
        is advected.

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
            self.uplift = tectonic[self.glIDs]
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
        2. a kdtree is built with the advected coordinates and used to interpolate the mesh variables
        on the initial local mesh position. The interpolation is based on a weighting distance function accounting for the 3 closest advected vertices.
        3. interpolated variables on the initial mesh coordinates are then stored in PETSc vectors and class parameters are updated accordingly.

        :arg tectonic: local tectonic rates in 3D
        :arg timer: tectonic time step in years
        """

        t1 = process_time()

        # Move local coordinates
        XYZ = self.lcoords + tectonic[self.glIDs, :] * timer

        # Update local elevation
        tmp = self.hLocal.getArray().copy()
        elev = tmp[self.lgIDs]
        if MPIsize > 1:
            temp = np.full(self.shadowgNb, -1.0e8, dtype=np.float64)
            temp[self.gshadinIDs] = elev[self.gshadowIDs]
            temp[self.gshadoutIDs] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            elev[self.shadowAlls] = temp
            loc_elev = elev[self.glIDs]
        else:
            loc_elev = elev[self.glIDs]

        # Update local erosion/deposition
        tmp = self.cumEDLocal.getArray().copy()
        erodep = tmp[self.lgIDs]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs] = erodep[self.gshadowIDs]
            temp[self.gshadoutIDs] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            erodep[self.shadowAlls] = temp
            loc_erodep = erodep[self.glIDs]
        else:
            loc_erodep = erodep[self.glIDs]

        # Update local stratal dataset after displacements
        if self.stratNb > 0 and self.stratStep > 0:
            if self.carbOn:
                (
                    loc_stratH,
                    loc_stratZ,
                    loc_stratF,
                    loc_stratC,
                    loc_phiS,
                    loc_phiF,
                    loc_phiC,
                ) = self._updateDispStrata()
            else:
                (
                    loc_stratH,
                    loc_stratZ,
                    loc_stratF,
                    loc_phiS,
                    loc_phiF,
                ) = self._updateDispStrata()

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
                self._stratalRecord(
                    indices,
                    weights,
                    loc_stratH,
                    loc_stratZ,
                    loc_stratF,
                    loc_stratC,
                    loc_phiS,
                    loc_phiF,
                    loc_phiC,
                )
            else:
                self._stratalRecord(
                    indices,
                    weights,
                    loc_stratH,
                    loc_stratZ,
                    loc_stratF,
                    None,
                    loc_phiS,
                    loc_phiF,
                    None,
                )

        if len(onIDs) > 0:
            nerodep[onIDs] = loc_erodep[indices[onIDs, 0]]
            if self.interp > 1:
                nelev[onIDs] = loc_elev[indices[onIDs, 0]]
            if self.stratNb > 0 and self.stratStep > 0:
                self.stratZ[onIDs, : self.stratStep] = loc_stratZ[indices[onIDs, 0], :]
                self.stratH[onIDs, : self.stratStep] = loc_stratH[indices[onIDs, 0], :]
                self.stratF[onIDs, : self.stratStep] = loc_stratF[indices[onIDs, 0], :]
                self.phiS[onIDs, : self.stratStep] = loc_phiS[indices[onIDs, 0], :]
                self.phiF[onIDs, : self.stratStep] = loc_phiF[indices[onIDs, 0], :]
                if self.carbOn:
                    self.stratC[onIDs, : self.stratStep] = loc_stratC[
                        indices[onIDs, 0], :
                    ]
                    self.phiC[onIDs, : self.stratStep] = loc_phiC[indices[onIDs, 0], :]

        self.hLocal.setArray(nelev)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

        self.cumEDLocal.setArray(nerodep)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        # Update local stratal dataset
        if self.stratNb > 0 and self.stratStep > 0:
            self._localStrat()
            del loc_stratH, loc_stratZ, loc_stratF, loc_phiS, loc_phiF
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

    def _updateDispStrata(self):
        """
        Gets the stratigraphic records relevant to each partition after mesh advection.

        The functions returns local stratigraphic layer information:

        - thickness of each stratigrapic layer `loc_stratH` accounting for both
          erosion & deposition events.
        - proportion of fine sediment `loc_stratF` contains in each stratigraphic layer.
        - elevation at time of deposition, considered to be to the current elevation
          for the top stratigraphic layer `loc_stratZ`.
        - porosity of coarse sediment `loc_phiS` in each stratigraphic layer computed at
          center of each layer.
        - porosity of fine sediment `loc_phiF` in each stratigraphic layer computed at
          center of each layer.
        - proportion of carbonate sediment `loc_strataC` contains in each stratigraphic layer
          if the carbonate module is turned on.
        - porosity of carbonate sediment `loc_phiC` in each stratigraphic layer computed at
          center of each layer when the carbonate module is turned on.

        .. note::

            In `gospl`, the stratigraphic layers are only defined locally. For interpolation on the
            edges of each partition it is important to ensure that all stratigraphic information on
            the adjacent nodes of the neighbouring partition are accessible. This is done by applying
            MPI `Allreduce` operation to the nodes parts of the overlaid ('shadow') zone.

        :return: loc_stratH, loc_stratZ, loc_stratF, loc_stratC, loc_phiS, loc_phiF, loc_phiC
        """

        stratH = self.stratH[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp = np.full((self.shadowgNb, self.stratStep), -1.0e8, dtype=np.float64)
            temp[self.gshadinIDs, :] = stratH[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            stratH[self.shadowAlls, :] = temp
            loc_stratH = stratH[self.glIDs, :]
        else:
            loc_stratH = stratH[self.glIDs, :]

        stratZ = self.stratZ[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = stratZ[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            stratZ[self.shadowAlls, :] = temp
            loc_stratZ = stratZ[self.glIDs, :]
        else:
            loc_stratZ = stratZ[self.glIDs, :]

        stratF = self.stratF[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = stratF[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            stratF[self.shadowAlls, :] = temp
            loc_stratF = stratF[self.glIDs, :]
        else:
            loc_stratF = stratF[self.glIDs, :]

        phiS = self.phiS[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = phiS[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            phiS[self.shadowAlls, :] = temp
            loc_phiS = phiS[self.glIDs, :]
        else:
            loc_phiS = phiS[self.glIDs, :]

        phiF = self.phiF[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = phiF[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            phiF[self.shadowAlls, :] = temp
            loc_phiF = phiF[self.glIDs, :]
        else:
            loc_phiF = phiF[self.glIDs, :]

        if self.carbOn:
            stratC = self.stratC[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = stratZ[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                stratC[self.shadowAlls, :] = temp
                loc_stratC = stratC[self.glIDs, :]
            else:
                loc_stratC = stratC[self.glIDs, :]

            phiC = self.phiC[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = phiC[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                phiC[self.shadowAlls, :] = temp
                loc_phiC = phiC[self.glIDs, :]
            else:
                loc_phiC = phiC[self.glIDs, :]

            return (
                loc_stratH,
                loc_stratZ,
                loc_stratF,
                loc_stratC,
                loc_phiS,
                loc_phiF,
                loc_phiC,
            )
        else:
            return loc_stratH, loc_stratZ, loc_stratF, loc_phiS, loc_phiF

    def _stratalRecord(
        self,
        indices,
        weights,
        loc_stratH,
        loc_stratZ,
        loc_stratF,
        loc_stratC,
        loc_phiS,
        loc_phiF,
        loc_phiC,
    ):
        """
        Once the interpolation has been performed, the following function updates the
        stratigraphic records based on the advected mesh.

        The function relies on 2 fortran subroutines (for loop performance purposes):

        1. strataBuild
        2. strataBuildCarb

        :arg indices: indices of the closest nodes used for interpolation
        :arg weights: weights based on the distances to closest nodes
        :arg loc_stratH: thickness of each stratigrapic layer accounting for both erosion & deposition events.
        :arg loc_stratH: proportion of fine sediment contains in each stratigraphic layer.
        :arg loc_stratZ: elevation at time of deposition, considered to be to the current elevation for the top stratigraphic layer.
        :arg loc_phiS: porosity of coarse sediment in each stratigraphic layer computed at center of each layer.
        :arg loc_phiF: porosity of fine sediment in each stratigraphic layer computed at center of each layer.
        :arg loc_strataC: proportion of carbonate sediment contains in each stratigraphic layer if the carbonate module is turned on.
        :arg loc_phiC: porosity of carbonate sediment in each stratigraphic layer computed at center of each layer when the carbonate module is turned on.

        """

        if self.carbOn:
            (
                self.stratH[:, : self.stratStep],
                self.stratZ[:, : self.stratStep],
                self.stratF[:, : self.stratStep],
                self.stratC[:, : self.stratStep],
                self.phiS[:, : self.stratStep],
                self.phiF[:, : self.stratStep],
                self.phiC[:, : self.stratStep],
            ) = strataBuildCarb(
                self.npoints,
                self.stratStep,
                indices,
                weights,
                loc_stratH,
                loc_stratZ,
                loc_stratF,
                loc_stratC,
                loc_phiS,
                loc_phiF,
                loc_phiC,
            )
        else:
            (
                self.stratH[:, : self.stratStep],
                self.stratZ[:, : self.stratStep],
                self.stratF[:, : self.stratStep],
                self.phiS[:, : self.stratStep],
                self.phiF[:, : self.stratStep],
            ) = strataBuild(
                self.npoints,
                self.stratStep,
                indices,
                weights,
                loc_stratH,
                loc_stratZ,
                loc_stratF,
                loc_phiS,
                loc_phiF,
            )

        return

    def _localStrat(self):
        """
        Updates stratigraphic records after mesh advection on the edges of each partition to ensure that
        all stratigraphic information on the adjacent nodes of the neighbouring partition are equals on
        all processors sharing a common number of nodes.
        """

        stratH = self.stratH[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp = np.full((self.shadowgNb, self.stratStep), -1.0e8, dtype=np.float64)
            temp[self.gshadinIDs, :] = stratH[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            stratH[self.shadowAlls, :] = temp
            self.stratH[:, : self.stratStep] = stratH[self.glIDs, :]
        else:
            self.stratH[:, : self.stratStep] = stratH[self.glIDs, :]

        stratZ = self.stratZ[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = stratZ[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            stratZ[self.shadowAlls, :] = temp
            self.stratZ[:, : self.stratStep] = stratZ[self.glIDs, :]
        else:
            self.stratZ[:, : self.stratStep] = stratZ[self.glIDs, :]

        stratF = self.stratF[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = stratF[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            stratF[self.shadowAlls, :] = temp
            self.stratF[:, : self.stratStep] = stratF[self.glIDs, :]
        else:
            self.stratF[:, : self.stratStep] = stratF[self.glIDs, :]

        phiS = self.phiS[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = phiS[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            phiS[self.shadowAlls, :] = temp
            self.phiS[:, : self.stratStep] = phiS[self.glIDs, :]
        else:
            self.phiS[:, : self.stratStep] = phiS[self.glIDs, :]

        phiF = self.phiF[self.lgIDs, : self.stratStep]
        if MPIsize > 1:
            temp.fill(-1.0e8)
            temp[self.gshadinIDs, :] = phiF[self.gshadowIDs, :]
            temp[self.gshadoutIDs, :] = -1.0e8
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
            phiF[self.shadowAlls, :] = temp
            self.phiF[:, : self.stratStep] = phiF[self.glIDs, :]
        else:
            self.phiF[:, : self.stratStep] = phiF[self.glIDs, :]
        del stratH, stratZ, stratF, phiS, phiF

        if self.carbOn:
            stratC = self.stratC[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = stratC[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                stratC[self.shadowAlls, :] = temp
                self.stratC[:, : self.stratStep] = stratC[self.glIDs, :]
            else:
                self.stratC[:, : self.stratStep] = stratC[self.glIDs, :]

            phiC = self.phiC[self.lgIDs, : self.stratStep]
            if MPIsize > 1:
                temp.fill(-1.0e8)
                temp[self.gshadinIDs, :] = phiC[self.gshadowIDs, :]
                temp[self.gshadoutIDs, :] = -1.0e8
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, temp, op=MPI.MAX)
                phiC[self.shadowAlls, :] = temp
                self.phiC[:, : self.stratStep] = phiC[self.glIDs, :]
            else:
                self.phiC[:, : self.stratStep] = phiC[self.glIDs, :]
            del stratC, phiC

        return

    def destroy_DMPlex(self):
        """
        Destroys PETSc DMPlex objects and associated PETSc local/global
        Vectors and Matrices at the end of the simulation.
        """

        t0 = process_time()

        self.hLocal.destroy()
        self.hGlobal.destroy()
        self.FAG.destroy()
        self.FAL.destroy()
        self.fillFAG.destroy()
        self.fillFAL.destroy()
        self.cumED.destroy()
        self.cumEDLocal.destroy()
        self.vSed.destroy()
        self.vSedLocal.destroy()
        if self.stratNb > 0:
            self.vSedf.destroy()
            self.vSedfLocal.destroy()
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
        self.stepED.destroy()
        self.Eb.destroy()
        self.EbLocal.destroy()
        self.upsG.destroy()
        self.upsL.destroy()

        self.iMat.destroy()
        if not self.fast:
            self.wMat.destroy()
            self.fillMat.destroy()
            if self.Cda > 0.0 or self.Cdm > 0.0:
                self.Diff.destroy()
        self.lgmap_col.destroy()
        self.lgmap_row.destroy()
        self.dm.destroy()

        del self.lcoords
        del self.lcells
        del self.inIDs
        del self.stratH
        del self.stratZ
        del self.stratF
        del self.phiS
        del self.phiF
        if self.carbOn:
            del self.stratC
            del self.phiC

        if not self.fast:
            del self.distRcv
            del self.wghtVal
            del self.rcvID

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
