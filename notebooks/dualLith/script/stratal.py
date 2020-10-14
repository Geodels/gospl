import gc
import sys
import glob
import h5py
import numpy as np
from scipy import spatial
import ruamel.yaml as yaml
from pyevtk.hl import gridToVTK
from scipy.ndimage import gaussian_filter


class stratal:
    def __init__(self, path=None, filename=None, layer=None):

        self.path = path
        if path is not None:
            filename = self.path + filename

        # Check input file exists
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
        self._inputParser()

        self.nbCPUs = len(glob.glob1(self.outputDir + "/h5/", "topology.p*"))

        if layer is not None:
            self.layNb = layer
        else:
            self.layNb = int((self.tEnd - self.tStart) / self.strat)

        print("Maximum layer number:", self.layNb)

        self.nbfile = len(glob.glob1(self.outputDir + "/h5/", "stratal.*.p0.h5"))

        return

    def _inputParser(self):

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
            self.strat = timeDict["strat"]
        except KeyError:
            print(
                "Key 'strat' is required to build the stratigraphy in the input file!"
            )
            raise KeyError("Simulation stratal time needs to be declared.")

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

        r = np.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
        # h = r - self.radius

        xs = np.array(self.x)
        ys = np.array(self.y)
        zs = np.array(self.z / r)

        lons = np.arctan2(ys, xs)
        lats = np.arcsin(zs)

        # Convert spherical mesh longitudes and latitudes to degrees
        self.lonlat = np.empty((len(self.x), 2))
        if lons.ndim == 1:
            self.lonlat[:, 0] = np.mod(np.degrees(lons) + 180.0, 360.0)
            self.lonlat[:, 1] = np.mod(np.degrees(lats) + 90, 180.0)
        else:
            self.lonlat[:, 0] = np.mod(np.degrees(lons[:, 0]) + 180.0, 360.0)
            self.lonlat[:, 1] = np.mod(np.degrees(lats[:, 0]) + 90, 180.0)

        self.tree = spatial.cKDTree(self.lonlat, leafsize=10)

        return

    def _lonlat2xyz(self, lon, lat, elev):
        """
        Convert lon / lat (radians) for the spherical triangulation into x,y,z
        on the unit sphere
        """

        xs = np.cos(lat) * np.cos(lon) * self.radius
        ys = np.cos(lat) * np.sin(lon) * self.radius
        zs = np.sin(lat) * self.radius + elev

        return xs, ys, zs

    def _getCoordinates(self):

        for k in range(self.nbCPUs):
            df = h5py.File("%s/h5/topology.p%s.h5" % (self.outputDir, k), "r")
            coords = np.array((df["/coords"]))
            if k == 0:
                self.x, self.y, self.z = np.hsplit(coords, 3)
            else:
                self.x = np.append(self.x, coords[:, 0])
                self.y = np.append(self.y, coords[:, 1])
                self.z = np.append(self.z, coords[:, 2])
            df.close()

        del coords
        self.nbPts = len(self.x)
        self._xyz2lonlat()

        gc.collect()

        return

    def _getData(self, nbfile):

        for k in range(self.nbCPUs):
            sf = h5py.File("%s/h5/stratal.%s.p%s.h5" % (self.outputDir, nbfile, k), "r")
            if k == 0:
                elev = np.array(sf["/stratZ"])
                th = np.array(sf["/stratH"])
                fine = np.array(sf["/stratF"])
                phiS = np.array(sf["/phiS"])
                phiF = np.array(sf["/phiF"])
                # self.sea[layNb] = np.array((sf["/sea"]))
            else:
                elev = np.append(elev, sf["/stratZ"], axis=0)
                th = np.append(th, sf["/stratH"], axis=0)
                fine = np.append(fine, sf["/stratF"], axis=0)
                phiS = np.append(phiS, sf["/phiS"], axis=0)
                phiF = np.append(phiF, sf["/phiF"], axis=0)
            sf.close()
        self.curLay = th.shape[1]
        print("Current number of sedimentary:", self.curLay)

        self.elev = elev
        self.th = th
        self.fine = fine
        self.phiS = phiS
        self.phiF = phiF

        return

    def readStratalData(self, fileNb=None):

        if fileNb is None:
            fileNb = self.nbfile - 1

        self._getCoordinates()
        # self.sea = np.empty(self.layNb-1)
        # self.elev = np.empty((self.nbPts, layNb-1))
        # self.erodep = np.empty((self.nbPts, layNb-1))

        # Find associated file:
        self._getData(fileNb)

        return

    def _test_progress(self, job_title, progress):

        length = 20
        block = int(round(length * progress))
        msg = "\r{0}: [{1}] {2}%".format(
            job_title, "#" * block + "-" * (length - block), round(progress * 100, 2)
        )
        if progress >= 1:
            msg += " DONE\r\n"
        sys.stdout.write(msg)
        sys.stdout.flush()

    def buildLonLatMesh(self, res=0.1, nghb=3):

        self.nx = int(360.0 / res)
        self.ny = int(180.0 / res)
        self.lon = np.linspace(0.0, 360.0, self.nx)
        self.lat = np.linspace(0, 180, self.ny)

        self.lon, self.lat = np.meshgrid(self.lon, self.lat)
        xyi = np.dstack([self.lon.flatten(), self.lat.flatten()])[0]
        self.zi = np.empty((self.curLay, self.ny, self.nx))
        self.thi = np.empty((self.curLay, self.ny, self.nx))
        self.finei = np.empty((self.curLay, self.ny, self.nx))
        self.phiSi = np.empty((self.curLay, self.ny, self.nx))
        self.phiFi = np.empty((self.curLay, self.ny, self.nx))

        distances, indices = self.tree.query(xyi, k=nghb)
        weights = 1.0 / distances ** 2
        denum = 1.0 / np.sum(weights, axis=1)
        onIDs = np.where(distances[:, 0] == 0)[0]

        print("Start building regular stratigraphic arrays")

        for k in range(self.curLay):

            zz = self.elev[:, k]
            th = self.th[:, k]
            fine = self.fine[:, k]
            phiS = self.phiS[:, k]
            phiF = self.phiF[:, k]
            self._test_progress("Percentage of arrays built ", (k + 1) / self.curLay)
            zi = np.sum(weights * zz[indices], axis=1) * denum
            thi = np.sum(weights * th[indices], axis=1) * denum
            finei = np.sum(weights * fine[indices], axis=1) * denum
            phiSi = np.sum(weights * phiS[indices], axis=1) * denum
            phiFi = np.sum(weights * phiF[indices], axis=1) * denum

            if len(onIDs) > 0:
                zi[onIDs] = zz[indices[onIDs, 0]]
                thi[onIDs] = th[indices[onIDs, 0]]
                finei[onIDs] = fine[indices[onIDs, 0]]
                phiSi[onIDs] = phiS[indices[onIDs, 0]]
                phiFi[onIDs] = phiF[indices[onIDs, 0]]

            self.zi[k, :, :] = np.reshape(zi, (self.ny, self.nx))
            self.thi[k, :, :] = np.reshape(thi, (self.ny, self.nx))
            self.finei[k, :, :] = np.reshape(finei, (self.ny, self.nx))
            self.phiSi[k, :, :] = np.reshape(phiSi, (self.ny, self.nx))
            self.phiFi[k, :, :] = np.reshape(phiFi, (self.ny, self.nx))

        del weights, denum, onIDs, zz, zi, th, xyi, thi, finei, phiSi, phiFi
        gc.collect()

        return

    def writeMesh(self, vtkfile="mesh", lons=None, lats=None, sigma=0.0):
        """
        Create a vtk unstructured grid based on current time step stratal parameters.
        """

        lon = np.linspace(0.0, 360.0, self.nx)
        lat = np.linspace(0.0, 180.0, self.ny)

        if lons is None:
            lons = [lon[0], lon[1]]
        else:
            lons[0] = np.where(self.lon[0, :] >= lons[0] + 180)[0].min()
            lons[1] = np.where(self.lon[0, :] >= lons[1] + 180)[0].min()
        if lats is None:
            lats = [lat[0], lat[1]]
        else:
            lats[0] = np.where(self.lat[:, 0] >= lats[0] + 90)[0].min()
            lats[1] = np.where(self.lat[:, 0] >= lats[1] + 90)[0].min()

        x = np.zeros((lons[1] - lons[0], lats[1] - lats[0], self.curLay))
        y = np.zeros((lons[1] - lons[0], lats[1] - lats[0], self.curLay))
        z = np.zeros((lons[1] - lons[0], lats[1] - lats[0], self.curLay))
        e = np.zeros((lons[1] - lons[0], lats[1] - lats[0], self.curLay))
        h = np.zeros((lons[1] - lons[0], lats[1] - lats[0], self.curLay))
        t = np.zeros((lons[1] - lons[0], lats[1] - lats[0], self.curLay))
        f = np.zeros((lons[1] - lons[0], lats[1] - lats[0], self.curLay))
        ps = np.zeros((lons[1] - lons[0], lats[1] - lats[0], self.curLay))
        pf = np.zeros((lons[1] - lons[0], lats[1] - lats[0], self.curLay))

        zz = self.zi[-1, lats[0] : lats[1], lons[0] : lons[1]]
        zz = gaussian_filter(zz, sigma)

        res = 360.0 / self.nx
        hscale = 11000.0 * res
        xmax = -1.0e12
        ymax = -1.0e12

        for k in range(self.curLay - 1, -1, -1):
            th = gaussian_filter(self.thi[k, :, :], sigma)
            th[th < 0] = 0.0
            if k < self.curLay - 1:
                thu = gaussian_filter(self.thi[k + 1, :, :], sigma)
                thu[thu < 0] = 0.0
            for j in range(lats[1] - lats[0]):
                for i in range(lons[1] - lons[0]):
                    x[i, j, k] = (lon[i + lons[0]] - lon[lons[0]]) * hscale
                    y[i, j, k] = (lat[j + lats[0]] - lat[lats[0]]) * hscale
                    xmax = max(xmax, x[i, j, k])
                    ymax = max(ymax, y[i, j, k])
                    if k == self.curLay - 1:
                        z[i, j, k] = zz[j, i]
                    else:
                        z[i, j, k] = z[i, j, k + 1] - thu[j + lats[0], i + lons[0]]
                    e[i, j, k] = self.zi[k, j + lats[0], i + lons[0]]
                    h[i, j, k] = th[j + lats[0], i + lons[0]]
                    f[i, j, k] = self.finei[k, j + lats[0], i + lons[0]]
                    ps[i, j, k] = self.phiSi[k, j + lats[0], i + lons[0]]
                    pf[i, j, k] = self.phiFi[k, j + lats[0], i + lons[0]]
                    t[i, j, k] = k

        gridToVTK(
            vtkfile,
            x,
            y,
            z,
            pointData={
                "dep elev": e,
                "th": h,
                "layID": t,
                "percfine": f,
                "phiC": ps,
                "phiF": pf,
            },
        )

        return [xmax, ymax]
