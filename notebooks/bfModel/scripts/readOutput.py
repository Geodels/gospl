import gc
import glob
import h5py
import meshio
import numpy as np
from scipy import spatial
import ruamel.yaml as yaml


class readOutput:
    def __init__(self, filename=None, step=None):

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
            erodep = np.array((df2["/erodep"]))
            sedLoad = np.array((df2["/sedLoad"]))
            flowAcc = np.array((df2["/flowAcc"]))
            uplift = np.array((df2["/uplift"]))
            hdisp = np.array((df2["/hdisp"]))

            if k == 0:
                x, y, z = np.hsplit(coords, 3)
                nelev = elev
                nrain = rain

                nerodep = erodep
                nsedLoad = sedLoad
                nflowAcc = flowAcc
                nhdisp = hdisp
                nuplift = uplift
            else:
                x = np.append(x, coords[:, 0])
                y = np.append(y, coords[:, 1])
                z = np.append(z, coords[:, 2])
                nelev = np.append(nelev, elev)
                nrain = np.append(nrain, rain)
                nerodep = np.append(nerodep, erodep)
                nsedLoad = np.append(nsedLoad, sedLoad)
                nflowAcc = np.append(nflowAcc, flowAcc)
                nhdisp = np.append(nhdisp, hdisp)
                nuplift = np.append(nuplift, uplift)

            df.close()

        self.nbPts = len(x)
        ncoords = np.zeros((self.nbPts, 3))
        ncoords[:, 0] = x.ravel()
        ncoords[:, 1] = y.ravel()
        ncoords[:, 2] = z.ravel()
        del coords, elev, erodep, hdisp, uplift, sedLoad, flowAcc, rain
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

            self.erodep = np.sum(weights * nerodep[indices][:, :], axis=1) / np.sum(
                weights, axis=1
            )
            self.sedLoad = np.sum(weights * nsedLoad[indices][:, :], axis=1) / np.sum(
                weights, axis=1
            )
            self.flowAcc = np.sum(weights * nflowAcc[indices][:, :], axis=1) / np.sum(
                weights, axis=1
            )
            self.hdisp = np.sum(weights * nhdisp[indices][:, :], axis=1) / np.sum(
                weights, axis=1
            )
            self.uplift = np.sum(weights * nuplift[indices][:, :], axis=1) / np.sum(
                weights, axis=1
            )

        else:
            self.elev = np.sum(weights * nelev[indices][:, :, 0], axis=1) / np.sum(
                weights, axis=1
            )
            self.rain = np.sum(weights * nrain[indices][:, :, 0], axis=1) / np.sum(
                weights, axis=1
            )

            self.erodep = np.sum(weights * nerodep[indices][:, :, 0], axis=1) / np.sum(
                weights, axis=1
            )
            self.sedLoad = np.sum(
                weights * nsedLoad[indices][:, :, 0], axis=1
            ) / np.sum(weights, axis=1)
            self.flowAcc = np.sum(
                weights * nflowAcc[indices][:, :, 0], axis=1
            ) / np.sum(weights, axis=1)
            self.hdisp = np.sum(weights * nhdisp[indices][:, :, 0], axis=1) / np.sum(
                weights, axis=1
            )
            self.uplift = np.sum(weights * nuplift[indices][:, :, 0], axis=1) / np.sum(
                weights, axis=1
            )

        if len(onIDs) > 0:
            self.elev[onIDs] = nelev[indices[onIDs, 0]]
            self.rain[onIDs] = nrain[indices[onIDs, 0]]

            self.erodep[onIDs] = nerodep[indices[onIDs, 0]]
            self.sedLoad[onIDs] = nsedLoad[indices[onIDs, 0]]
            self.flowAcc[onIDs] = nflowAcc[indices[onIDs, 0]]
            self.hdisp[onIDs] = nhdisp[indices[onIDs, 0]]
            self.uplift[onIDs] = nuplift[indices[onIDs, 0]]

        del weights, nelev, nrain, distances, indices, ncoords
        del nerodep, nhdisp, nuplift, nsedLoad, nflowAcc

        gc.collect()
        return

    def _getCoordinatesOld(self, step):

        for k in range(self.nbCPUs):
            df = h5py.File("%s/h5/topology.p%s.h5" % (self.outputDir, k), "r")
            coords = np.array((df["/coords"]))

            df2 = h5py.File("%s/h5/gospl.%s.p%s.h5" % (self.outputDir, step, k), "r")
            elev = np.array((df2["/elev"]))

            if k == 0:
                x, y, z = np.hsplit(coords, 3)
                nelev = elev
            else:
                x = np.append(x, coords[:, 0])
                y = np.append(y, coords[:, 1])
                z = np.append(z, coords[:, 2])
                nelev = np.append(nelev, elev)
            df.close()

        self.nbPts = len(x)
        ncoords = np.zeros((self.nbPts, 3))
        ncoords[:, 0] = x.ravel()
        ncoords[:, 1] = y.ravel()
        ncoords[:, 2] = z.ravel()
        del coords, elev
        gc.collect()

        # Load mesh structure
        mesh_struct = np.load(str(self.npdata) + ".npz")
        self.vertices = mesh_struct["v"]
        self.cells = mesh_struct["c"]

        tree = spatial.cKDTree(ncoords, leafsize=10)
        distances, indices = tree.query(self.vertices, k=3)

        # Inverse weighting distance...
        weights = 1.0 / distances ** 2
        onIDs = np.where(distances[:, 0] == 0)[0]
        if nelev[indices].ndim == 2:
            self.elev = np.sum(weights * nelev[indices][:, :], axis=1) / np.sum(
                weights, axis=1
            )
        else:
            self.elev = np.sum(weights * nelev[indices][:, :, 0], axis=1) / np.sum(
                weights, axis=1
            )

        if len(onIDs) > 0:
            self.elev[onIDs] = nelev[indices[onIDs, 0]]
        del weights, nelev, distances, indices, ncoords
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
                "vtec": self.uplift,
            },
        )
        meshio.write(vtkfile, vis_mesh)
        print("Writing VTK file {}".format(vtkfile))

        return

    def _readElevationData(self, step):

        self._getCoordinates(step)

        return
