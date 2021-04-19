import numpy as np
import pandas as pd
from netCDF4 import Dataset
from time import process_time
from scipy.spatial import cKDTree


class moveElev:
    def __init__(
        self, elevfile=None, dispfile=None, elev=None, erodep=None, scotese=True
    ):

        self.elevfile = elevfile
        self.dispfile = dispfile

        self.coords = None
        self.disps = None
        self.elev = None
        self.erodep = None
        self.lonlat = None
        self.new_elev = None
        self.new_erodep = None

        if elevfile is not None:
            if scotese:
                self.getElevScotese()
            else:
                self.getElevNPZ()
        else:
            self.elev = elev
            self.erodep = erodep
            self.getElev()

        self.getDisplacements()

        return

    def _lonlat2xyz(self, lon, lat, radius=6378137.0):

        rlon = np.radians(lon)
        rlat = np.radians(lat)

        self.coords = np.zeros((len(lon), 3))
        self.coords[:, 0] = np.cos(rlat) * np.cos(rlon) * radius
        self.coords[:, 1] = np.cos(rlat) * np.sin(rlon) * radius
        self.coords[:, 2] = np.sin(rlat) * radius

        return

    def _xyz2lonlat(self, coords):

        r = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2 + coords[:, 2] ** 2)

        xs = np.array(coords[:, 0])
        ys = np.array(coords[:, 1])
        zs = np.array(coords[:, 2] / r)

        lons = np.arctan2(ys, xs)
        lats = np.arcsin(zs)

        # Convert spherical mesh longitudes and latitudes to degrees
        lonlat = np.empty((len(coords[:, 0]), 2))
        lonlat[:, 0] = np.degrees(lons)
        lonlat[:, 1] = np.degrees(lats)

        return lonlat

    def getElevScotese(self):

        t0 = process_time()

        data = Dataset(self.elevfile, "r", format="NETCDF4")
        self.elev = np.flipud(data["z"][:, :])
        width = self.elev.shape[0]
        height = self.elev.shape[1]

        lon = np.linspace(-180.0, 180, height)
        lat = np.linspace(-90.0, 90, width)

        mlon, mlat = np.meshgrid(lon, lat)
        self.lonlat = np.dstack([mlon.flatten(), mlat.flatten()])[0]

        self._lonlat2xyz(mlon.flatten(), mlat.flatten())

        print("Read scotese map (%0.02f seconds)" % (process_time() - t0))

        return

    def getElev(self):

        t0 = process_time()

        lon = np.linspace(-180.0, 180, 3601)
        lat = np.linspace(-90.0, 90, 1801)

        mlon, mlat = np.meshgrid(lon, lat)
        self.lonlat = np.dstack([mlon.flatten(), mlat.flatten()])[0]

        self._lonlat2xyz(mlon.flatten(), mlat.flatten())

        print("Read npz elevation map (%0.02f seconds)" % (process_time() - t0))

        return

    def getElevNPZ(self):

        t0 = process_time()

        data = np.load(self.elevfile)
        self.elev = data["z"].reshape((1801, 3601))
        self.erodep = data["erodep"].reshape((1801, 3601))

        lon = np.linspace(-180.0, 180, 3601)
        lat = np.linspace(-90.0, 90, 1801)

        mlon, mlat = np.meshgrid(lon, lat)
        self.lonlat = np.dstack([mlon.flatten(), mlat.flatten()])[0]

        self._lonlat2xyz(mlon.flatten(), mlat.flatten())

        print("Read npz elevation map (%0.02f seconds)" % (process_time() - t0))

        return

    def getDisplacements(self):

        t0 = process_time()

        with open(self.dispfile) as f:
            for i, l in enumerate(f):
                pass

        maxlines = i

        # Open gPlates 1 degree 3D displacement maps (xy files)
        data = pd.read_csv(
            self.dispfile,
            sep=r"\s+",
            engine="c",
            header=None,
            skiprows=[0, 1, 2, 3, 4, 5, maxlines],
            error_bad_lines=True,
            na_filter=False,
            dtype=np.float64,
            low_memory=False,
        )

        dlon = data.values[:, 0]
        dlat = data.values[:, 1]

        # Conversion from cm/yr to m/yr
        tmpx = data.values[:, 2] / 100.0
        tmpy = data.values[:, 3] / 100.0
        tmpz = data.values[:, 4] / 100.0

        print("Read displacement map (%0.02f seconds)" % (process_time() - t0))

        t0 = process_time()

        # Create kdtree...
        tree = cKDTree(list(zip(dlon, dlat)))
        distances, indices = tree.query(self.lonlat, k=3)

        # Inverse weighting distance...
        weights = np.divide(
            1.0, distances ** 2, out=np.zeros_like(distances), where=distances != 0
        )
        sumweights = np.sum(weights, axis=1)

        onIDs = np.where(sumweights == 0)[0]
        sumweights[sumweights == 0] = 1.0e-4

        dispx = np.sum(weights * tmpx[indices][:, :], axis=1) / sumweights
        dispy = np.sum(weights * tmpy[indices][:, :], axis=1) / sumweights
        dispz = np.sum(weights * tmpz[indices][:, :], axis=1) / sumweights

        if len(onIDs) > 0:
            dispx[onIDs] = tmpx[indices[onIDs, 0]]
            dispy[onIDs] = tmpy[indices[onIDs, 0]]
            dispz[onIDs] = tmpz[indices[onIDs, 0]]

        self.disps = np.stack((dispx, dispy, dispz)).T

        print("Interpolate displacement map (%0.02f seconds)" % (process_time() - t0))

        return

    def applyDisplacements(self, dt=1.0e6):

        t0 = process_time()
        mvcoords = self.coords + self.disps * dt
        mvlonlat = self._xyz2lonlat(mvcoords)
        felev = self.elev.flatten()
        if self.erodep is not None:
            ferodep = self.erodep.flatten()

        tree = cKDTree(mvlonlat)
        distances, indices = tree.query(self.lonlat, k=3)

        # Inverse weighting distance...
        weights = np.divide(
            1.0, distances ** 2, out=np.zeros_like(distances), where=distances != 0
        )
        sumweights = np.sum(weights, axis=1)

        onIDs = np.where(sumweights == 0)[0]
        sumweights[sumweights == 0] = 1.0

        self.new_elev = np.sum(weights * felev[indices][:, :], axis=1) / sumweights
        if self.erodep is not None:
            self.new_erodep = (
                np.sum(weights * ferodep[indices][:, :], axis=1) / sumweights
            )

        if len(onIDs) > 0:
            self.new_elev[onIDs] = felev[indices[onIDs, 0]]
            if self.erodep is not None:
                self.new_erodep[onIDs] = ferodep[indices[onIDs, 0]]

        self.new_elev = self.new_elev.reshape(self.elev.shape)
        if self.erodep is not None:
            self.new_erodep = self.new_erodep.reshape(self.erodep.shape)

        print("Interpolate new elevation map (%0.02f seconds)" % (process_time() - t0))

        return
