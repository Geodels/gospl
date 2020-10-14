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

import argparse


def buildRegElevationMesh(infile, outfile, key="z"):

    t0 = process_time()
    data = Dataset(infile, "r", format="NETCDF4")
    img = np.fliplr(data[key][:, :].T)
    width = img.shape[0]
    height = img.shape[1]

    Lon = np.linspace(-180.0, 180, width)
    Lat = np.linspace(-90.0, 90, height)
    print("Read scotese map (%0.02f seconds)" % (process_time() - t0))

    value = np.round(img.flatten(), 3)
    t0 = process_time()
    f = open(outfile, "w+")
    f.write("mshid=3;ellipsoid-grid\n")
    f.write("mdims=2\n")

    f.write("coord=1;%d\n" % (len(Lon)))
    for k in range(len(Lon)):
        f.write(str(Lon[k]) + "\n")

    f.write("coord=2;%d\n" % (len(Lat)))
    for k in range(len(Lat)):
        f.write(str(Lat[k]) + "\n")

    f.write("value=%d;1\n" % (len(value)))
    for k in range(len(value)):
        f.write(str(value[k]) + "\n")

    f.close()
    print(
        "Writing Jigsaw topo input .msh file (%0.02f seconds)" % (process_time() - t0)
    )

    return


def getInitialMesh(topofile, meshfile, spacefile, outfile, dst_path, hfn):

    t0 = process_time()
    opts = jigsawpy.jigsaw_jig_t()
    topo = jigsawpy.jigsaw_msh_t()
    geom = jigsawpy.jigsaw_msh_t()
    mesh = jigsawpy.jigsaw_msh_t()
    hmat = jigsawpy.jigsaw_msh_t()

    jigsawpy.loadmsh(topofile, topo)
    print("Load topography grid (%0.02f seconds)" % (process_time() - t0))

    t0 = process_time()
    opts.geom_file = os.path.join(dst_path, "topology.msh")
    opts.jcfg_file = os.path.join(dst_path, "config.jig")
    opts.mesh_file = meshfile
    opts.hfun_file = spacefile

    geom.mshID = "ellipsoid-mesh"
    geom.radii = np.full(3, 6.371e003, dtype=geom.REALS_t)
    jigsawpy.savemsh(opts.geom_file, geom)

    hmat.mshID = "ellipsoid-grid"
    hmat.radii = geom.radii
    hmat.xgrid = topo.xgrid * np.pi / 180.0
    hmat.ygrid = topo.ygrid * np.pi / 180.0

    # Set HFUN gradient-limiter
    hmat.value = np.full(topo.value.shape, hfn[0], dtype=hmat.REALS_t)
    hmat.value[topo.value > -1000] = hfn[1]
    hmat.value[topo.value > 0] = hfn[2]

    hmat.slope = np.full(topo.value.shape, +0.050, dtype=hmat.REALS_t)
    jigsawpy.savemsh(opts.hfun_file, hmat)
    jigsawpy.cmd.marche(opts, hmat)
    print("Build space function (%0.02f seconds)" % (process_time() - t0))

    t0 = process_time()
    opts.hfun_scal = "absolute"
    opts.hfun_hmax = float("inf")  # null HFUN limits
    opts.hfun_hmin = float(+0.00)

    opts.mesh_dims = +2  # 2-dim. simplexes

    opts.optm_qlim = +9.5e-01  # tighter opt. tol
    opts.optm_iter = +32
    opts.optm_qtol = +1.0e-05

    jigsawpy.cmd.tetris(opts, 3, mesh)
    print("Perform triangulation (%0.02f seconds)" % (process_time() - t0))

    t0 = process_time()
    apos = jigsawpy.R3toS2(geom.radii, mesh.point["coord"][:])

    apos = apos * 180.0 / np.pi

    zfun = interpolate.RectBivariateSpline(topo.ygrid, topo.xgrid, topo.value)

    mesh.value = zfun(apos[:, 1], apos[:, 0], grid=False)

    jigsawpy.savevtk(outfile, mesh)
    jigsawpy.savemsh(opts.mesh_file, mesh)
    print("Get unstructured mesh (%0.02f seconds)" % (process_time() - t0))

    return


def gosplElev(coords, cells, elev, gmesh, visvtk=False):

    Gmesh = meshplex.mesh_tri.MeshTri(coords, cells)
    s = Gmesh.idx_hierarchy.shape
    a = np.sort(Gmesh.idx_hierarchy.reshape(s[0], -1).T)
    Gmesh.edges = {"nodes": np.unique(a, axis=0)}
    ngbNbs, ngbID = definegtin(len(coords), Gmesh.cells["nodes"], Gmesh.edges["nodes"])

    np.savez_compressed(gmesh, v=coords, c=cells, n=ngbID[:, :8].astype(int), z=elev)

    if visvtk:
        paleovtk = gmesh + ".vtk"
        vis_mesh = meshio.Mesh(coords, {"triangle": cells}, point_data={"z": elev})
        meshio.write(paleovtk, vis_mesh)
        print("Writing VTK file {}".format(paleovtk))

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


def gosplRain(coords, cells, paleorain, rainmesh, visvtk=False, filter=2):

    lonlat = xyz2lonlat(coords, radius=6378137.0)

    data = Dataset(paleorain, "r", format="NETCDF4")
    paleorain = data["z"][:, :].T

    # Map mesh coordinates
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
        np.float
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

    return


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
        dtype=np.float,
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


def list_str(values):
    return values.split(",")


def main():
    start_time = process_time()

    # Parsing command line arguments
    parser = argparse.ArgumentParser(
        description="Generation of high resolution mesh for goSPL based on JIGSAW.",
        add_help=True,
    )

    parser.add_argument("-t", "--time", help="Time interval to build", required=True)
    parser.add_argument("-d", "--data", help="Forcing input folder", required=True)
    parser.add_argument("-s", "--space", type=list_str, required=True)

    args = parser.parse_args()

    time = args.time
    src_path = args.data
    hfn = np.zeros(3)
    hfn[0] = float(args.space[0])
    hfn[1] = float(args.space[1])
    hfn[2] = float(args.space[2])

    input_path = "input" + str(time)
    if not os.path.exists(input_path):
        os.makedirs(input_path)

    # EXTRACT PALEO-ELEVATION AS A REGULAR GRID
    infile = os.path.join(src_path, "paleomap/" + str(time) + "Ma.nc")
    outfile = os.path.join(input_path, str(time) + "Ma.msh")
    buildRegElevationMesh(infile, outfile, key="z")

    # DEFINE JIGSAW MESHING DATASET
    topofile = os.path.join(input_path, str(time) + "Ma.msh")
    meshfile = os.path.join(input_path, "mesh" + str(time) + ".msh")
    spacefile = os.path.join(input_path, "spac" + str(time) + ".msh")
    outfile = os.path.join(input_path, "mesh" + str(time) + ".vtk")
    getInitialMesh(topofile, meshfile, spacefile, outfile, input_path, hfn)

    # READ UNSTRUCTURED MESH
    umesh = meshio.read(outfile)
    coords = umesh.points
    coords = (coords / 6.371e003) * 6378137.0
    cells = umesh.cells_dict["triangle"]
    elev = umesh.point_data["value"]

    # DEFINE goSPL INPUT FILES

    # 1- elevation mesh
    npzelev = os.path.join(input_path, str(time) + "Ma")
    gosplElev(coords, cells, elev, npzelev, visvtk=False)

    # 2- precipitation mesh
    npzrain = os.path.join(input_path, "rain" + str(time) + "Ma")
    paleoRain = os.path.join(src_path, "precipitation/" + str(time) + "Ma.nc")
    gosplRain(coords, cells, paleoRain, npzrain, visvtk=False, filter=2)

    print("Creation goSPL files took (%0.02f seconds)" % (process_time() - start_time))


if __name__ == "__main__":

    main()
