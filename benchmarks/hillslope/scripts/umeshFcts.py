import os
import gc
import numpy as np
import xarray as xr
import uxarray as uxr
from scipy import interpolate
from scipy import spatial
from scipy.interpolate import griddata
import pyinterp

import vtk
from vtk.util import numpy_support

import jigsawpy
from jigsawpy.savejig import savejig
from mpas_tools.logging import check_call

import xesmf as xe

from mpas_tools.mesh.interpolation import interp_bilin
from mpas_tools.cime.constants import constants
from mpas_tools.mesh.creation import build_spherical_mesh
from mpas_tools.mesh.creation.jigsaw_to_netcdf import jigsaw_to_netcdf
from mpas_tools.viz.paraview_extractor import extract_vtk


def buildGlobalMeshSimple(widthCell, foldername):
    '''
    This function builds a global UGRID mesh using UXarray.

    .. note::

       UGRID conventions are a foundation to represent Unstructured Grids.
       These conventions are intended to describe how these grids should be
       stored  within a NetCDF file, with a particular focus on environmental
       and geoscience applications.

    By default the constructed UGRID mesh will be name mesh_XXkm.nc where
    ``XX`` is the provided widthCell.

    This function uses the library jigsaw to build the mesh.

    :arg widthCell: the width of the voronoi mesh
    :arg foldername: the name of the folder used to store the mesh
    '''

    if not os.path.exists(foldername):
        os.makedirs(foldername)

    cellWidth, lon, lat = cellWidthVsLatLon(widthCell)
    earthRadius = constants['SHR_CONST_REARTH']
    ufile = foldername + '/mesh_'+str(widthCell)+'km.nc'

    if not os.path.exists(ufile):
        build_spherical_mesh(cellWidth, lon, lat,
                             earth_radius=earthRadius,
                             out_filename=ufile,
                             plot_cellWidth=False)

    return


def refineGlobalMesh(widthCell, lon, lat, foldername):
    '''
    This function builds a refined UGRID mesh using UXarray and variable cell widths.

    .. note::

       UGRID conventions are a foundation to represent Unstructured Grids.
       These conventions are intended to describe how these grids should be
       stored  within a NetCDF file, with a particular focus on environmental
       and geoscience applications.

    By default the constructed UGRID mesh will be name mesh_XXkm.nc where
    ``XX`` is the provided widthCell.

    This function uses the library jigsaw to build the mesh.

    :arg widthCell: the width of the voronoi mesh
    :arg lat: latitudinal extent of the data grid (nc file)
    :arg lon: longitudinal extent of the data grid (nc file)
    :arg foldername: the name of the folder used to store the mesh
    '''

    if not os.path.exists(foldername):
        os.makedirs(foldername)

    earthRadius = constants['SHR_CONST_REARTH']
    ufile = foldername + '/mesh_refine.nc'

    if not os.path.exists(ufile):
        build_spherical_mesh(widthCell, lon, lat,
                             earth_radius=earthRadius,
                             out_filename=ufile,
                             plot_cellWidth=False)

    return


def xyz2lonlat(X, Y, Z):

    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    xs = np.array(X)
    ys = np.array(Y)
    zs = np.array(Z/r)

    lons = np.arctan2(ys, xs)
    lats = np.arcsin(zs)

    # Convert spherical mesh longitudes and latitudes to degrees
    lonlat = np.empty((len(X), 2))
    lonlat[:, 0] = np.mod(np.degrees(lons) + 180.0, 360.0) - 180.0
    lonlat[:, 1] = np.mod(np.degrees(lats) + 90, 180.0) - 90.0

    return lonlat

def planarMesh(nds, outfolder, fvtk=None, fumpas=True, voro=True):

    opts = jigsawpy.jigsaw_jig_t()
    opts.geom_file = outfolder+'/mesh2D.msh'
    opts.jcfg_file = outfolder+'/mesh2D.jig'
    opts.mesh_file = outfolder+'/mesh2D-MESH.msh'
    opts.hfun_file = outfolder+'/mesh2D-HFUN.msh'

    newx = nds.x.values
    newy = nds.y.values
    hmat = jigsawpy.jigsaw_msh_t()
    hmat.mshID = 'EUCLIDEAN-GRID'
    hmat.xgrid = newx
    hmat.ygrid = newy
    hmat.value = nds.cellwidth.values
    jigsawpy.savemsh(opts.hfun_file, hmat)

    typeEDGE2 = {'names': ['index', 'IDtag'], 'formats': [('<i4', (2,)), '<i4'], 'offsets': [0, 8], 'itemsize': 12, 'aligned': True}
    typeVERT2 = {'names': ['coord', 'IDtag'], 'formats': [('<f8', (2,)), '<i4'], 'offsets': [0, 16], 'itemsize': 24, 'aligned': True}

    geom_points = np.array([
        ((newx[0], newy[0]),0),          # outer square
        ((newx[-1], newy[0]),0),
        ((newx[-1], newy[-1]),0),
        ((newx[0], newy[-1]),0)],
        dtype=typeVERT2)

    geom_edges = np.array([
        ((0, 1), 0),          # outer square
        ((1, 2), 0),
        ((2, 3), 0),
        ((3, 0), 0)],
        dtype=typeEDGE2)

    geom = jigsawpy.jigsaw_msh_t()
    geom.mshID = 'EUCLIDEAN-MESH'
    geom.vert2 = geom_points
    geom.edge2 = geom_edges
    jigsawpy.savemsh(opts.geom_file, geom)

    opts.hfun_scal = 'absolute'
    opts.hfun_hmax = float("inf")
    opts.hfun_hmin = 0.0
    opts.mesh_dims = +2  # 2-dim. simplexes
    opts.optm_qlim = 0.9375
    opts.verbosity = +1

    savejig(opts.jcfg_file, opts)
    check_call(['jigsaw', opts.jcfg_file], logger=None)

    jigsawpy.loadmsh(opts.mesh_file, geom)
    if fvtk is not None:
        jigsawpy.savevtk(outfolder+'/'+fvtk, geom)

    if fumpas:
        jigsaw_to_netcdf(msh_filename=opts.mesh_file,
                         output_name=outfolder+'/mesh2D.nc', on_sphere=False)
        args = ['MpasMeshConverter.x', outfolder+'/mesh2D.nc', 
                outfolder+'/base2D.nc']
        check_call(args=args)

    if fumpas and voro:
        extract_vtk(
            filename_pattern=outfolder+'/base2D.nc',
            variable_list='areaCell',
            dimension_list=['maxEdges=', 'nVertLevels=', 'nParticles='], 
            mesh_filename=outfolder+'/base2D.nc',
            out_dir=outfolder,
            ignore_time=True,
            lonlat=False,
            xtime='none'
            )

    return


def cellWidthVsLatLon(width=50):
    """
    Create cell width array for this mesh on a regular latitude-longitude grid.
    Returns
    -------
    cellWidth : ndarray
        m x n array of cell width in km
    lon : ndarray
        longitude in degrees (length n and between -180 and 180)
    lat : ndarray
        longitude in degrees (length m and between -90 and 90)
    """
    dlat = 1
    dlon = 1
    # 20 km | 10 km | 5 km
    constantCellWidth = width

    nlat = int(180/dlat) + 1
    nlon = int(360/dlon) + 1

    lat = np.linspace(-90., 90., nlat)
    lon = np.linspace(-180., 180., nlon)

    cellWidth = constantCellWidth * np.ones((lat.size, lon.size))

    return cellWidth, lon, lat


def cellWidthVsLatLonFunc_simple(ncgrid):
    """
    Create cell width array for this mesh on a regular latitude-longitude grid.
    Returns
    -------
    cellWidth : ndarray
        m x n array of cell width in km
    lon : ndarray
        longitude in degrees (length n and between -180 and 180)
    lat : ndarray
        longitude in degrees (length m and between -90 and 90)
    """

    ds = ncgrid[['h']]
    lat = ds.lat.values
    lon = ds.lon.values
    ds['cellwidth'] = (['lat', 'lon'], 50*np.ones((lat.size, lon.size)))
    ds['cellwidth'] = ds['cellwidth'].where(ds.h < 0, 10)
    ds['cellwidth'] = ds['cellwidth'].where((ds.h < -500) | (ds.h >= 0), 15)
    ds['cellwidth'] = ds['cellwidth'].where((ds.h < -1000) | (ds.h >= -500), 20)
    ds['cellwidth'] = ds['cellwidth'].where((ds.h < -2000) | (ds.h >= -1000), 25)
    ds['cellwidth'] = ds['cellwidth'].where((ds.h < -3000) | (ds.h >= -2000), 30)
    ds['cellwidth'] = ds['cellwidth'].where((ds.h < -5000) | (ds.h >= -3000), 35)
    ds['cellwidth'] = ds['cellwidth'].where((ds.h < -8000) | (ds.h >= -5000), 40)
    ds['cellwidth'] = ds['cellwidth'].where((ds.h < -14000) | (ds.h >= -8000), 45)

    cellWidth = ds['cellwidth'].values

    return ds, cellWidth, lon, lat


def cellWidthVsLatLonFuncPower(ncgrid, width=[10, 50], power=1.2,
                               maxdepth=10000):
    """
    Create cell width array for this mesh on a regular latitude-longitude grid.
    Returns
    -------
    cellWidth : ndarray
        m x n array of cell width in km
    lon : ndarray
        longitude in degrees (length n and between -180 and 180)
    lat : ndarray
        longitude in degrees (length m and between -90 and 90)
    """

    ds = ncgrid[['h']]
    lat = ds.lat.values
    lon = ds.lon.values

    # Normalise
    val = -ncgrid.h.copy()
    val = val.where(val > 0, 0)
    val = val.where(val < maxdepth, maxdepth)
    val = val/maxdepth
    # Power fitting
    xp = np.linspace(0, 1, 100)
    fpow = interpolate.interp1d(xp, (width[1]-width[0])*xp**power+width[0])
    interv = fpow(val.values.flatten())
    ds['cellwidth'] = (['lat', 'lon'], interv.reshape(val.shape))
    cellWidth = ds['cellwidth'].values

    return ds, cellWidth, lon, lat


def cellWidthVsLatLonFuncDist(ncgrid, width=[10, 50], maxdist=10000,
                              latlon=True):
    """
    Create cell width array for this mesh on a regular latitude-longitude grid.
    Returns
    -------
    cellWidth : ndarray
        m x n array of cell width in km
    lon : ndarray
        longitude in degrees (length n and between -180 and 180)
    lat : ndarray
        longitude in degrees (length m and between -90 and 90)
    """

    # Normalise
    val = ncgrid.coast.copy()
    val = val.where(val < maxdist, maxdist)
    val = val/maxdist
    # Power fitting
    xp = np.linspace(0, 1, 100)
    fpow = interpolate.interp1d(xp, (width[1]-width[0])*xp+width[0])
    interv = fpow(val.values.flatten())

    if latlon:
        ncgrid['cellwidth'] = (['lat', 'lon'], interv.reshape(val.shape))
    else:
        ncgrid['cellwidth'] = (['y', 'x'], interv.reshape(val.shape))
    cellWidth = ncgrid['cellwidth'].values

    return ncgrid, cellWidth


def inter2UGRID(ncgrid, ugrid, nfolder, varname, type='face', latlon=True):
    '''
    This function interpolates from a regular netcdf grid to a UGRID mesh.

    :arg ncgrid: the regular grid containing the variables to interpolate onto the UGRID
    :arg ugrid: the UGRID where the variables would be interpolated onto
    :arg nfolder: the name of the folder used to variables
    :arg varname: the name of the UGRID variables file
    :arg type: are the variables interpolated on the nodes or the face of the UGRID
    :arg latlon: set to True if running a global model
    '''

    if not os.path.exists(nfolder):
        os.makedirs(nfolder)
    
    # Pick mesh target coordinates
    if type == 'face':
        if latlon:
            meshLon = ugrid['lonCell'].values
            meshLat = ugrid['latCell'].values
            dim_name = 'nCells'
        else:
            meshLon = ugrid['xCell'].values
            meshLat = ugrid['yCell'].values
            dim_name = 'nCells'
    elif type == 'node':
        if latlon:
            meshLon = ugrid['lonVertex'].values
            meshLat = ugrid['latVertex'].values
            dim_name = 'nVertices'
        else:
            meshLon = ugrid['xVertex'].values
            meshLat = ugrid['yVertex'].values
            dim_name = 'nVertices'
    else:
        print('Function only allows 2 types: face or node')
        return
    if latlon:
        dlon = ncgrid.lon.values
        dlat = ncgrid.lat.values
    else:
        dlon = ncgrid.x.values
        dlat = ncgrid.y.values

    # Interpolate each variable and collect into a plain xarray Dataset
    data_vars = {}
    for key in ncgrid.data_vars:
        vData = interp_bilin(dlon, dlat, ncgrid[key].values, meshLon, meshLat)
        data_vars[key] = xr.DataArray(vData, dims=[dim_name], name=key)
    out_ds = xr.Dataset(data_vars)
    out_ds.to_netcdf(os.path.join(nfolder, varname + '.nc'))

    return


def performInterp(coords, ncoords, data, mthd='IDW'):

    mesh = pyinterp.RTree()
    mesh.packing(ncoords, data)

    if mthd == 'IDW':
        val, neighbors = mesh.inverse_distance_weighting(
            coords,
            within=False,  # Extrapolation is forbidden
            k=11,  # We are looking for at most 11 neighbors
            num_threads=0)
    elif mthd == 'RBF':
        val, neighbors = mesh.radial_basis_function(
            coords,
            within=False,  # Extrapolation is forbidden
            k=11,  # We are looking for at most 11 neighbors
            rbf='linear',
            smooth=1e-4,
            num_threads=0)
    elif mthd == 'KRG':
        val, neighbors = mesh.universal_kriging(
            coords,
            within=False,  # Extrapolation is forbidden
            k=11,
            covariance='matern_12',
            alpha=100_000,
            num_threads=0)
    else:
        print('Specified method is not supported')

    return val


def mvNodes(x, y, z, vx, vy, vz, dt):

    nx = x + vx * dt
    ny = y + vy * dt
    nz = z + vz * dt

    # Convert spherical coordinates to lon/lat coordinates
    mvlonlat = xyz2lonlat(nx, ny, nz)

    return np.vstack((mvlonlat[:, 0], mvlonlat[:, 1])).T


def get_Tectonic(ufile, data_file, vkeys, zkeys, dkey, dt, mthd='IDW'):

    # Open the UGRID file
    dual_mesh = uxr.open_dataset(ufile, *data_file, use_dual=True)
    coords = np.vstack((dual_mesh.uxgrid.node_lon, dual_mesh.uxgrid.node_lat)).T

    # Move nodes according to displacement
    ncoords = mvNodes(dual_mesh.uxgrid.node_x, dual_mesh.uxgrid.node_y,
                      dual_mesh.uxgrid.node_z, dual_mesh[vkeys[0]],
                      dual_mesh[vkeys[1]], dual_mesh[vkeys[2]],
                      dt)

    # Interpolate data
    if dkey is not None:
        data = dual_mesh[zkeys[0]] + dual_mesh[dkey] * dt
    else:
        data = dual_mesh[zkeys[0]]
    zval = performInterp(coords, ncoords, data, mthd)
    dispTec = dual_mesh[zkeys[1]] - zval

    # Move backward nodes according to displacement
    ncoords = mvNodes(dual_mesh.uxgrid.node_x, dual_mesh.uxgrid.node_y,
                      dual_mesh.uxgrid.node_z, dual_mesh[vkeys[0]],
                      dual_mesh[vkeys[1]], dual_mesh[vkeys[2]],
                      -dt)

    # Interpolate data
    tval = performInterp(coords, ncoords, dispTec, mthd)
    if dkey is not None:
       dual_mesh['tec'] = ('n_node', tval/dt + dual_mesh[dkey])
    else:
       dual_mesh['tec'] = ('n_node', tval/dt)

    return dual_mesh

def generateVTKmesh(points, cells):
    """
    A global VTK mesh is generated to compute the distance between mesh vertices and coastlines position.

    The distance to the coastline for every marine vertices is used to define a maximum shelf slope during deposition. The coastline contours are efficiently obtained from VTK contouring function. This function is performed on a VTK mesh which is built in this function.
    """

    vtkMesh = vtk.vtkUnstructuredGrid()

    # Define mesh vertices
    vtk_points = vtk.vtkPoints()
    vtk_array = numpy_support.numpy_to_vtk(points, deep=True)
    vtk_points.SetData(vtk_array)
    vtkMesh.SetPoints(vtk_points)

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

    vtkMesh.SetCells(
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

    return vtkMesh


def globalCoastsTree(coastXYZ, points, seaID, k_neighbors=1):
    """
    This function takes all local coastline points and computes locally the distance of all marine points to the coastline.

    :arg coastXYZ: local coastline coordinates
    :arg k_neighbors: number of nodes to use when querying the kd-tree
    """

    coastDist = np.zeros(len(points))

    # Get coastlines points globally
    tree = spatial.cKDTree(coastXYZ, leafsize=10)
    coastDist[seaID], _ = tree.query(
        points[seaID, :], k=k_neighbors
    )

    del tree
    gc.collect()

    return coastDist


def distanceCoasts(vtkMesh, points, data, sl, k_neighbors=1):
    """
    This function computes for every marine vertices the distance to the closest coastline.

    .. important::

        The calculation takes advantage of the `vtkContourFilter` function from VTK library
        which is performed on the **global** VTK mesh. Once the coastlines have been extracted,
        the distances are obtained by querying a kd-tree (initialised with the coastal nodes) for
        marine vertices contained within each partition.

    :arg data: local elevation numpy array
    :arg k_neighbors: number of nodes to use when querying the kd-tree
    """

    pointData = vtkMesh.GetPointData()
    array = numpy_support.numpy_to_vtk(data, deep=1)
    array.SetName("elev")
    pointData.AddArray(array)
    vtkMesh.SetFieldData(pointData)

    cf = vtk.vtkContourFilter()
    cf.SetInputData(vtkMesh)
    cf.SetValue(0, sl)
    cf.SetInputArrayToProcess(0, 0, 0, 0, "elev")
    cf.GenerateTrianglesOff()
    cf.Update()
    if cf.GetOutput().GetPoints() is not None:
        coastXYZ = numpy_support.vtk_to_numpy(cf.GetOutput().GetPoints().GetData())
    else:
        coastXYZ = np.zeros((0, 3), dtype=np.float64)

    # Get coastlines points globally
    seaID = np.where(data <= sl)[0]

    coastDist = globalCoastsTree(coastXYZ, points, seaID, k_neighbors)

    del array, pointData, cf, coastXYZ
    gc.collect()

    return coastDist


def getGridCoast(ncgrid, dual_mesh, dcoast, input_path):

    ds_locs = xr.Dataset()

    ds_locs["lon"] = xr.DataArray(
        data=dual_mesh.uxgrid.node_lon.values, dims=("locations")
    )
    ds_locs["lat"] = xr.DataArray(
        data=dual_mesh.uxgrid.node_lat.values, dims=("locations")
    )
    ds_locs["coast"] = xr.DataArray(
        data=dcoast, dims=("locations")
    )

    if not os.path.exists(input_path+"/weights_distcoast.nc"):
        regridder_loc = xe.Regridder(
            ds_locs, ncgrid, "nearest_s2d", locstream_in=True
        )
        regridder_loc.to_netcdf(input_path+'weights_distcoast.nc')
    else:
        regridder_loc = xe.Regridder(ds_locs, ncgrid, "nearest_s2d", 
                                     filename=input_path+"weights_distcoast.nc", reuse_weights=True)

    return regridder_loc(ds_locs)


def getGridCoast2D(ncgrid, dual_mesh, dcoast):

    xi, yi = np.meshgrid(ncgrid.x.values, ncgrid.y.values)
    zi = griddata((dual_mesh.uxgrid.node_x.values,
                   dual_mesh.uxgrid.node_y.values), dcoast, (xi, yi),
                   method='cubic')
    zi[zi < 0] = 0.
    ncgrid['coast'] = (['y', 'x'], zi)

    return ncgrid
