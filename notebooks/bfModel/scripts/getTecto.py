import glob
import h5py
import numpy as np
import ruamel.yaml as yaml
from scipy import spatial


class getTecto:
    def __init__(self, filename=None):

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

        self.nbfile = len(glob.glob1(self.outputDir + "/h5/", "gospl.*.p0.h5"))

        mesh_struct = np.load(self.meshFile)
        self.vertices = mesh_struct["v"]

        return

    def _inputParser(self):

        try:
            domainDict = self.input["domain"]
        except KeyError:
            print(
                "Key 'domain' is required and is missing in the input file!", flush=True
            )
            raise KeyError("Key domain is required in the input file!")

        try:
            meshFile = domainDict["npdata"]
        except KeyError:
            print(
                "Key 'npdata' is required and is missing in the 'domain' declaration!",
                flush=True,
            )
            raise KeyError("Compressed numpy dataset definition is not defined!")

        self.meshFile = meshFile + ".npz"

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
            self.tout = timeDict["tout"]
        except KeyError:
            print(
                "Key 'tout' is required to build the vertical displacement in the input file!"
            )
            raise KeyError("Simulation output time needs to be declared.")

        try:
            outDict = self.input["output"]
            try:
                self.outputDir = outDict["dir"]
            except KeyError:
                self.outputDir = "output"
        except KeyError:
            self.outputDir = "output"

        return

    def _getCoordinates(self):

        for k in range(self.nbCPUs):
            df = h5py.File("%s/h5/topology.p%s.h5" % (self.outputDir, k), "r")
            coords = np.array((df["/coords"]))
            if k == 0:
                x, y, z = np.hsplit(coords, 3)
            else:
                x = np.append(x, coords[:, 0])
                y = np.append(y, coords[:, 1])
                z = np.append(z, coords[:, 2])
            df.close()

        del coords
        self.nbPts = len(x)

        self.coords = np.zeros((self.nbPts, 3))
        self.coords[:, 0] = x.ravel()
        self.coords[:, 1] = y.ravel()
        self.coords[:, 2] = z.ravel()

        self.tree = spatial.cKDTree(self.coords, leafsize=10)

        return

    def _getH5Data(self, nbfile):
        for k in range(self.nbCPUs):
            sf = h5py.File("%s/h5/gospl.%s.p%s.h5" % (self.outputDir, nbfile, k), "r")
            if k == 0:
                uplift = np.array(sf["/uplift"])
            else:
                uplift = np.append(uplift, sf["/uplift"])
            sf.close()
        print("Read uplift/subsidence for step:", nbfile)

        return uplift

    def readData(self, out=None, time=None):

        self._getCoordinates()
        distances, id = self.tree.query(self.vertices, k=1)

        for p in range(0, self.nbfile - 1):
            outfile = out + str(time[p]) + "Ma"
            uplift = self._getH5Data(p + 1)

            nuplift = uplift[id]

            # Save the mesh as compressed numpy file for global simulation
            np.savez_compressed(outfile, z=nuplift)

        return
