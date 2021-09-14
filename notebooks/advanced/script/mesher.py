import os
import meshio
import meshplex
import jigsawpy
import numpy as np
import pandas as pd
from scipy import ndimage
from netCDF4 import Dataset
from scipy import interpolate
from time import process_time
from scipy.spatial import cKDTree
from gospl._fortran import definegtin


def buildInputs(
    fdata, idata, hfn0=100.0, hfn1=30.0, hfn2=15, visvtk=True, reverse=False
):

    jigsawMeshScotese(fdata[0], idata[0], key="z", res=0.1)

    coords, cells, elev = unstMesh(fdata[3], idata[0], hfn0, hfn1, hfn2)

    gosplElev(coords, cells, elev, idata[1], visvtk)

    lonlat = gosplRain(coords, cells, fdata[1], idata[2], visvtk)

    gosplDisp(coords, lonlat, cells, fdata[2], idata[3], visvtk, reverse)

    return


def jigsawMeshScotese(infile, outfile, key="z", res=0.1):

    t0 = process_time()
    data = Dataset(infile, "r", format="NETCDF4")
    img = np.fliplr(data[key][:, :].T)
    # width = img.shape[0]
    # height = img.shape[1]

    # etopoLon = np.linspace(-180.0, 180, width)
    # etopoLat = np.linspace(-90.0, 90, height)

    res = 0.1
    Lon = np.round(np.arange(-180.0, 180 + res, res), 2)
    Lat = np.round(np.arange(-90.0, 90 + res, res), 2)

    print("Read scotese map (%0.02f seconds)" % (process_time() - t0))

    # t0 = process_time()
    # mLon, mLat = np.meshgrid(Lon, Lat, sparse=False, indexing="ij")
    # print("Building Jigsaw interpolation grid (%0.02f seconds)" % (process_time() - t0))

    value = np.round(img.flatten(), 3)

    t0 = process_time()
    f = open(outfile, "w+")
    f.write("MSHID=2;EUCLIDEAN-GRID\n")
    f.write("NDIMS=2\n")

    f.write("COORD=1;%d\n" % (len(Lon)))
    for k in range(len(Lon)):
        f.write(str(Lon[k]) + "\n")

    f.write("COORD=2;%d\n" % (len(Lat)))
    for k in range(len(Lat)):
        f.write(str(Lat[k]) + "\n")

    f.write("VALUE=%d;1\n" % (len(value)))
    for k in range(len(value)):
        f.write(str(value[k]) + "\n")

    f.close()
    print("Writing Jigsaw input .msh file (%0.02f seconds)" % (process_time() - t0))

    return


def xyz2lonlat(coords, radius=6378137.0):
    """
    Convert x,y,z representation of cartesian points of the
    spherical triangulation to lat / lon (radians).
    """

    gLonLat = np.zeros((len(coords), 2))

    gLonLat[:, 1] = np.arcsin(coords[:, 2] / radius)
    gLonLat[:, 0] = np.arctan2(coords[:, 1], coords[:, 0])
    gLonLat[:, 1] = np.mod(np.degrees(gLonLat[:, 1]) + 90, 180.0)
    gLonLat[:, 0] = np.mod(np.degrees(gLonLat[:, 0]) + 180.0, 360.0)

    return gLonLat


def unstMesh(esph, infile, hfn0=100.0, hfn1=30.0, hfn2=15, tmp="tmp/"):

    opts = jigsawpy.jigsaw_jig_t()
    topo = jigsawpy.jigsaw_msh_t()
    geom = jigsawpy.jigsaw_msh_t()
    mesh = jigsawpy.jigsaw_msh_t()
    hmat = jigsawpy.jigsaw_msh_t()

    opts.geom_file = esph
    opts.jcfg_file = os.path.join(tmp, "topo.jig")
    opts.mesh_file = os.path.join(tmp, "mesh.msh")
    opts.hfun_file = os.path.join(tmp, "spac.msh")

    # define JIGSAW geometry
    geom.mshID = "ellipsoid-mesh"
    geom.radii = np.full(3, 6.371e003, dtype=geom.REALS_t)

    jigsawpy.savemsh(opts.geom_file, geom)

    # define spacing pattern
    jigsawpy.loadmsh(infile, topo)

    hmat.mshID = "ellipsoid-grid"
    hmat.radii = geom.radii

    hmat.xgrid = topo.xgrid * np.pi / 180.0
    hmat.ygrid = topo.ygrid * np.pi / 180.0

    # Define grid resolution for each elevation domain
    hmat.value = np.full(topo.value.shape, hfn0, dtype=hmat.REALS_t)

    hmat.value[topo.value > -1000.0] = hfn1
    hmat.value[topo.value > 0.0] = hfn2

    # Set HFUN gradient-limiter
    hmat.slope = np.full(topo.value.shape, +0.025, dtype=hmat.REALS_t)

    jigsawpy.savemsh(opts.hfun_file, hmat)
    jigsawpy.cmd.marche(opts, hmat)

    # Make mesh using JIGSAW
    opts.hfun_scal = "absolute"
    opts.hfun_hmax = float("inf")  # null HFUN limits
    opts.hfun_hmin = float(+0.00)

    opts.mesh_dims = +2  # 2-dim. simplexes

    opts.optm_qlim = +9.5e-01  # tighter opt. tol
    opts.optm_iter = +32
    opts.optm_qtol = +1.0e-05

    jigsawpy.cmd.tetris(opts, 4, mesh)

    scr2 = jigsawpy.triscr2(  # "quality" metric
        mesh.point["coord"], mesh.tria3["index"]
    )

    apos = jigsawpy.R3toS2(geom.radii, mesh.point["coord"][:])

    apos = apos * 180.0 / np.pi

    zfun = interpolate.RectBivariateSpline(topo.ygrid, topo.xgrid, topo.value)

    mesh.value = zfun(apos[:, 1], apos[:, 0], grid=False)

    coords = (mesh.vert3["coord"] / 6.371e003) * 6378137.0

    del scr2

    return coords, mesh.tria3["index"], mesh.value


def gosplElev(coords, cells, elev, gmesh, visvtk=False):

    Gmesh = meshplex.MeshTri(coords, cells)
    s = Gmesh.idx_hierarchy.shape
    a = np.sort(Gmesh.idx_hierarchy.reshape(s[0], -1).T)

    if meshplex.__version__ >= "0.16.0":
        Gmesh.edges = {"points": np.unique(a, axis=0)}
        ngbNbs, ngbID = definegtin(
            len(coords), Gmesh.cells("points"), Gmesh.edges["points"]
        )
    elif meshplex.__version__ >= "0.14.0":
        Gmesh.edges = {"points": np.unique(a, axis=0)}
        ngbNbs, ngbID = definegtin(
            len(coords), Gmesh.cells["points"], Gmesh.edges["points"]
        )
    else:
        Gmesh.edges = {"nodes": np.unique(a, axis=0)}
        ngbNbs, ngbID = definegtin(
            len(coords), Gmesh.cells["nodes"], Gmesh.edges["nodes"]
        )

    np.savez_compressed(gmesh, v=coords, c=cells, n=ngbID[:, :8].astype(int), z=elev)

    if visvtk:
        paleovtk = gmesh + ".vtk"
        vis_mesh = meshio.Mesh(coords, {"triangle": cells}, point_data={"z": elev})
        meshio.write(paleovtk, vis_mesh)
        print("Writing VTK file {}".format(paleovtk))

    return


def gosplRain(coords, cells, paleorain, rainmesh, visvtk=False, filter=2):

    lonlat = xyz2lonlat(coords, radius=6378137.0)

    data = Dataset(paleorain, "r", format="NETCDF4")
    paleorain = data["z"][:, :].T

    # Map mesh coordinates on ETOPO1 dataset
    ilons = paleorain.shape[0] * lonlat[:, 0] / float(paleorain.shape[0])
    ilats = paleorain.shape[1] * lonlat[:, 1] / float(paleorain.shape[1])

    icoords = np.stack((ilons, ilats))
    paleorain = ndimage.gaussian_filter(paleorain, sigma=filter)

    rlons = icoords[0, :] * paleorain.shape[0] / 360.0
    rlats = icoords[1, :] * paleorain.shape[1] / 180.0

    rcoords = np.zeros(icoords.shape)
    rcoords[0, :] = rlons
    rcoords[1, :] = rlats

    # Interpolate the paleogrid on global mesh
    meshd = ndimage.map_coordinates(paleorain, rcoords, order=2, mode="nearest").astype(
        np.float64,
    )

    # Conversion from mm/day to m/yr
    meshd *= 365.2422 / 1000.0

    # Save the mesh as compressed numpy file for global simulation
    np.savez_compressed(rainmesh, r=meshd)

    if visvtk:
        paleovtk = rainmesh + ".vtk"
        vis_mesh = meshio.Mesh(coords, {"triangle": cells}, point_data={"r": meshd})
        meshio.write(paleovtk, vis_mesh)
        print("Writing VTK file {}".format(paleovtk))

    return icoords


def gosplDisp(coords, lonlat, cells, paleodisp, dispmesh, visvtk=False, reverse=False):

    with open(paleodisp) as f:
        for i, l in enumerate(f):
            pass

    maxlines = i

    # Open gPlates 1 degree 3D displacement maps (xy files)
    data = pd.read_csv(
        paleodisp,
        sep=r"\s+",
        engine="c",
        header=None,
        skiprows=[0, 1, 2, 3, 4, 5, maxlines],
        error_bad_lines=True,
        na_filter=False,
        dtype=np.float64,
        low_memory=False,
    )

    # Build the kdtree
    lon = data.values[:, 0] + 180.0
    lat = data.values[:, 1] + 90.0
    # Create kdtree...
    tree2 = cKDTree(list(zip(lon, lat)))
    d2, inds2 = tree2.query(list(zip(lonlat[0, :], lonlat[1, :])), k=1)

    # Conversion from cm/yr to m/yr
    if reverse:
        tmpx = -data.values[:, 2] / 100.0
        tmpy = -data.values[:, 3] / 100.0
        tmpz = -data.values[:, 4] / 100.0
    else:
        tmpx = data.values[:, 2] / 100.0
        tmpy = data.values[:, 3] / 100.0
        tmpz = data.values[:, 4] / 100.0

    # speed = np.sqrt(np.square(tmpx) + np.square(tmpy) + np.square(tmpz))

    # Interpolate the paleo displacement on global mesh
    dX = tmpx.flatten()[inds2].reshape(lonlat[0, :].shape)
    dY = tmpy.flatten()[inds2].reshape(lonlat[0, :].shape)
    dZ = tmpz.flatten()[inds2].reshape(lonlat[0, :].shape)

    disps = np.stack((dX, dY, dZ)).T

    # Save the mesh as compressed numpy file for global simulation
    np.savez_compressed(dispmesh, xyz=disps)

    if visvtk:
        paleovtk = dispmesh + ".vtk"
        vis_mesh = meshio.Mesh(
            coords,
            {"triangle": cells},
            point_data={"x": disps[:, 0], "y": disps[:, 1], "z": disps[:, 2]},
        )
        meshio.write(paleovtk, vis_mesh)
        print("Writing VTK file {}".format(paleovtk))

    return
