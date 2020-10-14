import gc
import glob
import h5py
import meshio
import numpy as np
from scipy import spatial
import ruamel.yaml as yaml


class readOutput:
    def __init__(self, path=None, filename=None, step=None, back=False, uplift=True):

        # Check input file exists
        self.path = path
        if path is not None:
            filename = self.path + filename

        try:
            with open(filename) as finput:
                pass
        except IOError:
            print("Unable to open file: ", filename)
            raise IOError("The input file is not found...")

        # Open YAML file
        with open(filename, "r") as finput:
            self.input = yaml.load(finput, Loader=yaml.Loader)

        self.radius = 6378137.0
        self.back = back
        self.lookuplift = uplift
        self._inputParser()

        self.nbCPUs = len(glob.glob1(self.outputDir + "/h5/", "topology.p*"))

        self._readElevationData(step)

        return

    def _inputParser(self):

        try:
            domainDict = self.input["domain"]
        except KeyError:
            print("Key 'domain' is required and is missing in the input file!")
            raise KeyError("Key domain is required in the input file!")

        try:
            self.npdata = domainDict["npdata"]
        except KeyError:
            print(
                "Key 'npdata' is required and is missing in the 'domain' declaration!"
            )
            raise KeyError("Simulation npdata needs to be declared.")

        try:
            timeDict = self.input["time"]
        except KeyError:
            print("Key 'time' is required and is missing in the input file!")
            raise KeyError("Key time is required in the input file!")

        try:
            self.tStart = timeDict["start"]
        except KeyError:
            print("Key 'start' is required and is missing in the 'time' declaration!")
            raise KeyError("Simulation start time needs to be declared.")

        try:
            self.tEnd = timeDict["end"]
        except KeyError:
            print("Key 'end' is required and is missing in the 'time' declaration!")
            raise KeyError("Simulation end time needs to be declared.")

        try:
            outDict = self.input["output"]
            try:
                self.outputDir = outDict["dir"]
            except KeyError:
                self.outputDir = "output"
        except KeyError:
            self.outputDir = "output"

        if self.path is not None:
            self.outputDir = self.path + self.outputDir

        return

    def _xyz2lonlat(self):

        r = np.sqrt(
            self.vertices[:, 0] ** 2
            + self.vertices[:, 1] ** 2
            + self.vertices[:, 2] ** 2
        )
        # h = r - self.radius

        xs = np.array(self.vertices[:, 0])
        ys = np.array(self.vertices[:, 1])
        zs = np.array(self.vertices[:, 2] / r)

        lons = np.arctan2(ys, xs)
        lats = np.arcsin(zs)

        # Convert spherical mesh longitudes and latitudes to degrees
        self.lonlat = np.empty((len(self.vertices[:, 0]), 2))
        self.lonlat[:, 0] = np.mod(np.degrees(lons) + 180.0, 360.0)
        self.lonlat[:, 1] = np.mod(np.degrees(lats) + 90, 180.0)

        self.tree = spatial.cKDTree(self.lonlat, leafsize=10)

        return

    def _getCoordinates(self, step):

        nerodep = None
        nsedLoad = None
        nflowAcc = None
        nhdisp = None
        nuplift = None

        for k in range(self.nbCPUs):
            df = h5py.File("%s/h5/topology.p%s.h5" % (self.outputDir, k), "r")
            coords = np.array((df["/coords"]))

            df2 = h5py.File("%s/h5/gospl.%s.p%s.h5" % (self.outputDir, step, k), "r")
            elev = np.array((df2["/elev"]))
            rain = np.array((df2["/rain"]))
            if not self.back:
                erodep = np.array((df2["/erodep"]))
                sedLoad = np.array((df2["/sedLoad"]))
                flowAcc = np.array((df2["/flowAcc"]))
                if self.lookuplift:
                    uplift = np.array((df2["/uplift"]))
                    hdisp = np.array((df2["/hdisp"]))

            if k == 0:
                x, y, z = np.hsplit(coords, 3)
                nelev = elev
                nrain = rain

                if not self.back:
                    nerodep = erodep
                    nsedLoad = sedLoad
                    nflowAcc = flowAcc
                    if self.lookuplift:
                        nhdisp = hdisp
                        nuplift = uplift
            else:
                x = np.append(x, coords[:, 0])
                y = np.append(y, coords[:, 1])
                z = np.append(z, coords[:, 2])
                nelev = np.append(nelev, elev)
                nrain = np.append(nrain, rain)

                if not self.back:
                    nerodep = np.append(nerodep, erodep)
                    nsedLoad = np.append(nsedLoad, sedLoad)
                    nflowAcc = np.append(nflowAcc, flowAcc)
                    if self.lookuplift:
                        nhdisp = np.append(nhdisp, hdisp)
                        nuplift = np.append(nuplift, uplift)

            df.close()

        self.nbPts = len(x)
        ncoords = np.zeros((self.nbPts, 3))
        ncoords[:, 0] = x.ravel()
        ncoords[:, 1] = y.ravel()
        ncoords[:, 2] = z.ravel()
        if not self.back:
            if self.lookuplift:
                del coords, elev, erodep, hdisp, uplift, sedLoad, flowAcc, rain
            else:
                del coords, elev
        else:
            del coords, elev
        gc.collect()

        self._getCoordinates2(
            ncoords, nelev, nrain, nerodep, nsedLoad, nflowAcc, nhdisp, nuplift
        )

    def _getCoordinates2(
        self, ncoords, nelev, nrain, nerodep, nsedLoad, nflowAcc, nhdisp, nuplift
    ):

        # Load mesh structure
        mesh_struct = np.load(str(self.npdata) + ".npz")
        self.vertices = mesh_struct["v"]
        self.cells = mesh_struct["c"]
        self._xyz2lonlat()

        tree = spatial.cKDTree(ncoords, leafsize=10)
        distances, indices = tree.query(self.vertices, k=3)

        # Inverse weighting distance...
        weights = 1.0 / distances ** 2
        onIDs = np.where(distances[:, 0] == 0)[0]
        if nelev[indices].ndim == 2:
            self.elev = np.sum(weights * nelev[indices][:, :], axis=1) / np.sum(
                weights, axis=1
            )
            self.rain = np.sum(weights * nrain[indices][:, :], axis=1) / np.sum(
                weights, axis=1
            )

            if not self.back:
                self.erodep = np.sum(weights * nerodep[indices][:, :], axis=1) / np.sum(
                    weights, axis=1
                )
                self.sedLoad = np.sum(
                    weights * nsedLoad[indices][:, :], axis=1
                ) / np.sum(weights, axis=1)
                self.flowAcc = np.sum(
                    weights * nflowAcc[indices][:, :], axis=1
                ) / np.sum(weights, axis=1)
                if self.lookuplift:
                    self.hdisp = np.sum(
                        weights * nhdisp[indices][:, :], axis=1
                    ) / np.sum(weights, axis=1)
                    self.uplift = np.sum(
                        weights * nuplift[indices][:, :], axis=1
                    ) / np.sum(weights, axis=1)

        else:
            self.elev = np.sum(weights * nelev[indices][:, :, 0], axis=1) / np.sum(
                weights, axis=1
            )
            self.rain = np.sum(weights * nrain[indices][:, :, 0], axis=1) / np.sum(
                weights, axis=1
            )

            if not self.back:
                self.erodep = np.sum(
                    weights * nerodep[indices][:, :, 0], axis=1
                ) / np.sum(weights, axis=1)
                self.sedLoad = np.sum(
                    weights * nsedLoad[indices][:, :, 0], axis=1
                ) / np.sum(weights, axis=1)
                self.flowAcc = np.sum(
                    weights * nflowAcc[indices][:, :, 0], axis=1
                ) / np.sum(weights, axis=1)
                if self.lookuplift:
                    self.hdisp = np.sum(
                        weights * nhdisp[indices][:, :, 0], axis=1
                    ) / np.sum(weights, axis=1)
                    self.uplift = np.sum(
                        weights * nuplift[indices][:, :, 0], axis=1
                    ) / np.sum(weights, axis=1)

        if len(onIDs) > 0:
            self.elev[onIDs] = nelev[indices[onIDs, 0]]
            self.rain[onIDs] = nrain[indices[onIDs, 0]]

            if not self.back:
                self.erodep[onIDs] = nerodep[indices[onIDs, 0]]
                self.sedLoad[onIDs] = nsedLoad[indices[onIDs, 0]]
                self.flowAcc[onIDs] = nflowAcc[indices[onIDs, 0]]
                if self.lookuplift:
                    self.hdisp[onIDs] = nhdisp[indices[onIDs, 0]]
                    self.uplift[onIDs] = nuplift[indices[onIDs, 0]]

        if not self.back:
            if self.lookuplift:
                del weights, nelev, nrain, distances, indices, ncoords
                del nerodep, nhdisp, nuplift, nsedLoad, nflowAcc
            else:
                del (
                    weights,
                    nelev,
                    nrain,
                    distances,
                    indices,
                    ncoords,
                    nsedLoad,
                    nflowAcc,
                )
        else:
            del weights, nelev, nrain, distances, indices, ncoords

        gc.collect()
        return

    def exportVTK(self, vtkfile):

        vis_mesh = meshio.Mesh(
            self.vertices,
            {"triangle": self.cells},
            point_data={
                "elev": self.elev,
                "erodep": self.erodep,
                "rain": self.rain,
                "FA": np.log(self.flowAcc),
                "SL": self.sedLoad,
            },
        )
        meshio.write(vtkfile, vis_mesh)
        print("Writing VTK file {}".format(vtkfile))

        return

    def _readElevationData(self, step):

        self._getCoordinates(step)

        return

    def buildLonLatMesh(self, res=0.1, nghb=3):

        self.nx = int(360.0 / res) + 1
        self.ny = int(180.0 / res) + 1
        self.lon = np.linspace(0.0, 360.0, self.nx)
        self.lat = np.linspace(0, 180, self.ny)

        self.lon, self.lat = np.meshgrid(self.lon, self.lat)
        xyi = np.dstack([self.lon.flatten(), self.lat.flatten()])[0]
        # self.zi = np.empty((self.ny, self.nx))
        # self.thi = np.empty((self.ny, self.nx))

        distances, indices = self.tree.query(xyi, k=nghb)
        weights = 1.0 / distances ** 2
        denum = 1.0 / np.sum(weights, axis=1)
        onIDs = np.where(distances[:, 0] == 0)[0]

        zi = np.sum(weights * self.elev[indices], axis=1) * denum
        raini = np.sum(weights * self.rain[indices], axis=1) * denum
        if not self.back:
            erodepi = np.sum(weights * self.erodep[indices], axis=1) * denum

            if self.lookuplift:
                sedLoadi = np.sum(weights * self.sedLoad[indices], axis=1) * denum
                uplifti = np.sum(weights * self.uplift[indices], axis=1) * denum
                flowAcci = np.sum(weights * self.flowAcc[indices], axis=1) * denum
                hdispi = np.sum(weights * self.hdisp[indices], axis=1) * denum

        if len(onIDs) > 0:
            zi[onIDs] = self.elev[indices[onIDs, 0]]
            raini[onIDs] = self.rain[indices[onIDs, 0]]
            if not self.back:
                erodepi[onIDs] = self.erodep[indices[onIDs, 0]]

                if self.lookuplift:
                    sedLoadi[onIDs] = self.sedLoad[indices[onIDs, 0]]
                    uplifti[onIDs] = self.uplift[indices[onIDs, 0]]
                    flowAcci[onIDs] = self.flowAcc[indices[onIDs, 0]]
                    hdispi[onIDs] = self.hdisp[indices[onIDs, 0]]

        self.rain = np.reshape(raini, (self.ny, self.nx))
        self.z = np.reshape(zi, (self.ny, self.nx))

        if not self.back:
            self.th = np.reshape(erodepi, (self.ny, self.nx))

            if self.lookuplift:
                self.sl = np.reshape(sedLoadi, (self.ny, self.nx))
                self.vdisp = np.reshape(uplifti, (self.ny, self.nx))
                self.hdisp = np.reshape(hdispi, (self.ny, self.nx))
                self.fa = np.reshape(flowAcci, (self.ny, self.nx))

        if not self.back:

            if self.lookuplift:
                del (
                    weights,
                    denum,
                    onIDs,
                    zi,
                    raini,
                    erodepi,
                    xyi,
                    sedLoadi,
                    uplifti,
                    flowAcci,
                    hdispi,
                )
            else:
                del (weights, denum, onIDs, raini, zi)
        else:
            del (weights, denum, onIDs, raini, zi)
        gc.collect()

        return
