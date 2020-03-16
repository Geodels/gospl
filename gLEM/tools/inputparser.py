import os
import errno
import numpy as np
import pandas as pd
import ruamel.yaml as yaml
from operator import itemgetter
from scipy.interpolate import interp1d

from mpi4py import MPI
import sys
import petsc4py

from petsc4py import PETSc

petsc4py.init(sys.argv)
MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIsize = PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class ReadYaml(object):
    """
    Parsing YAML input file

    Args:
        filename: input filename (.yml YAML file)
    """

    def __init__(self, filename):

        if self.showlog:
            self.log = PETSc.Log()
            self.log.begin()

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

        if MPIrank == 0 and "name" in self.input.keys() and self.verbose:
            print("The following model will be run:     {}".format(self.input["name"]))

        # Read simulation parameters
        self._readDomain()
        self._readTime()
        self._readSPL()
        self._readHillslope()
        self._readSealevel()
        self._readTectonic()
        self._readRain()
        self._readPaleo()
        self._readFlex()
        self._readOut()
        self._readForcePaleo()
        self._paleoFit()

        self.gravity = 9.81
        self.tNow = self.tStart
        self.saveTime = self.tNow
        if self.strat > 0:
            self.saveStrat = self.tNow + self.strat
        else:
            self.saveStrat = self.tEnd + self.tout

        return

    def _readDomain(self):
        """
        Read domain definition, boundary conditions and flow direction parameters.
        """

        try:
            domainDict = self.input["domain"]
        except KeyError:
            print("Key 'domain' is required and is missing in the input file!")
            raise KeyError("Key domain is required in the input file!")

        try:
            self.radius = domainDict["radius"]
        except KeyError:
            self.radius = 6378137.0

        try:
            self.flowDir = domainDict["flowdir"]
        except KeyError:
            self.flowDir = 6

        try:
            self.reflevel = domainDict["refinement"]
        except KeyError:
            self.reflevel = 0

        try:
            meshFile = domainDict["npdata"]
        except KeyError:
            print(
                "Key 'npdata' is required and is missing in the 'domain' declaration!"
            )
            raise KeyError("Compressed numpy dataset definition is not defined!")

        if self.reflevel > 0:
            self.meshFile = meshFile + str(self.reflevel) + ".npz"
        else:
            self.meshFile = meshFile + ".npz"

        try:
            with open(self.meshFile) as meshfile:
                meshfile.close()
                pass
        except IOError:
            print("Unable to open numpy dataset: ", self.meshFile)
            raise IOError("The numpy dataset is not found...")

        try:
            self.fast = domainDict["fast"]
        except KeyError:
            self.fast = False

        try:
            self.backward = domainDict["backward"]
        except KeyError:
            self.backward = False

        return

    def _readTime(self):
        """
        Read simulation time declaration.
        """

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

        if self.tEnd <= self.tStart:
            raise ValueError("Simulation end/start times do not make any sense!")

        try:
            self.dt = timeDict["dt"]
        except KeyError:
            print("Key 'dt' is required and is missing in the 'time' declaration!")
            raise KeyError("Simulation discretisation time step needs to be declared.")

        try:
            self.tout = timeDict["tout"]
        except KeyError:
            self.tout = self.tEnd - self.tStart
            print(
                "Output time interval 'tout' has been set to {} years".format(self.tout)
            )

        self._addTime(timeDict)

        return

    def _addTime(self, timeDict):
        """
        Read additional time parameters.
        """

        try:
            self.rStep = timeDict["rstep"]
        except KeyError:
            self.rStep = 0

        if self.tout < self.dt:
            self.tout = self.dt
            print(
                "Output time interval was changed to {} years to match the time step dt".format(
                    self.dt
                )
            )

        try:
            self.tecStep = timeDict["tec"]
        except KeyError:
            self.tecStep = self.tout

        try:
            self.strat = timeDict["strat"]
        except KeyError:
            self.strat = 0

        return

    def _readSPL(self):
        """
        Read surface processes bedrock parameters.
        """

        try:
            splDict = self.input["spl"]
            try:
                self.K = splDict["K"]
            except KeyError:
                print(
                    "When using the Surface Process Model definition of coefficient Kb is required."
                )
                raise ValueError("Surface Process Model: Kb coefficient not found.")
            try:
                self.frac_fine = splDict["Ff"]
            except KeyError:
                self.frac_fine = 0.0
            try:
                # `wgth` is the percentage of upstream sediment flux that will be deposited on each cell...
                self.wgth = splDict["wgth"]
                if self.wgth >= 1.0:
                    self.wgth = 0.999
            except KeyError:
                self.wgth = 0.0

        except KeyError:
            self.K = 1.0e-12
            self.wgth = 0.0
            self.frac_fine = 0.0

        return

    def _readHillslope(self):
        """
        Read hillslope parameters.
        # """

        try:
            hillDict = self.input["diffusion"]
            try:
                self.Cd = hillDict["hillslopeK"]
            except KeyError:
                print(
                    "When declaring diffusion processes, the coefficient hillslopeK is required."
                )
                raise ValueError("Hillslope: Cd coefficient not found.")
            try:
                self.sedimentK = hillDict["sedimentK"]
            except KeyError:
                self.sedimentK = 10.0
        except KeyError:
            self.Cd = 0.0
            self.sedimentK = 10.0

        return

    def _readSealevel(self):
        """
        Define sealevel evolution.
        """

        seafile = None
        sealevel = 0.0
        self.seafunction = None
        try:
            seaDict = self.input["sea"]
            try:
                sealevel = seaDict["position"]
                try:
                    seafile = seaDict["curve"]
                except KeyError:
                    seafile = None
            except KeyError:
                try:
                    seafile = seaDict["curve"]
                except KeyError:
                    seafile = None
        except KeyError:
            sealevel = 0.0

        if seafile is not None:
            try:
                with open(seafile) as fsea:
                    fsea.close()
                    try:
                        seadata = pd.read_csv(
                            seafile,
                            sep=r",",
                            engine="c",
                            header=None,
                            na_filter=False,
                            dtype=np.float,
                            low_memory=False,
                        )
                        pass
                    except ValueError:
                        try:
                            seadata = pd.read_csv(
                                seafile,
                                sep=r"\s+",
                                engine="c",
                                header=None,
                                na_filter=False,
                                dtype=np.float,
                                low_memory=False,
                            )
                            pass
                        except ValueError:
                            print(
                                "The sea-level file is not well formed: it should be comma or tab separated"
                            )
                            raise ValueError("Wrong formating of sea-level file.")
            except IOError:
                print("Unable to open file: ", seafile)
                raise IOError("The sealevel file is not found...")

            seadata[1] += sealevel
            if seadata[0].min() > self.tStart:
                tmpS = []
                tmpS.insert(0, {0: self.tStart, 1: seadata[1].iloc[0]})
                seadata = pd.concat([pd.DataFrame(tmpS), seadata], ignore_index=True)
            if seadata[0].max() < self.tEnd:
                tmpE = []
                tmpE.insert(0, {0: self.tEnd, 1: seadata[1].iloc[-1]})
                seadata = pd.concat([seadata, pd.DataFrame(tmpE)], ignore_index=True)
            self.seafunction = interp1d(
                seadata[0], seadata[1] + sealevel, kind="linear"
            )
        else:
            year = np.linspace(self.tStart, self.tEnd + self.dt, num=11, endpoint=True)
            seaval = np.full(len(year), sealevel)
            self.seafunction = interp1d(year, seaval, kind="linear")

        return

    def _storeTectonic(self, k, tecStart, zMap, tMap, tStep, tEnd, tecdata):

        if tMap is not None:
            if self.meshFile != tMap:
                try:
                    with open(tMap) as tecfile:
                        tecfile.close()
                        pass
                except IOError:
                    print("Unable to open tectonic file: ", tMap)
                    raise IOError(
                        "The tectonic file {} is not found for climatic event {}.".format(
                            tMap, k
                        )
                    )
        else:
            tMap = "empty"

        if zMap is not None:
            if self.meshFile != zMap:
                try:
                    with open(zMap) as tecfile:
                        tecfile.close()
                        pass
                except IOError:
                    print("Unable to open tectonic file: ", zMap)
                    raise IOError(
                        "The tectonic file {} is not found for climatic event {}.".format(
                            zMap, k
                        )
                    )
        else:
            zMap = "empty"

        if tMap == "empty" and zMap == "empty":
            print("For each tectonic event a tectonic grid (mapH or mapV) is required.")
            raise ValueError("Tectonic event {} has no tectonic map (map).".format(k))

        tmpTec = []
        tmpTec.insert(0, {"start": tecStart, "tMap": tMap, "zMap": zMap})

        if k == 0:
            tecdata = pd.DataFrame(tmpTec, columns=["start", "tMap", "zMap"])
        else:
            tecdata = pd.concat(
                [tecdata, pd.DataFrame(tmpTec, columns=["start", "tMap", "zMap"])],
                ignore_index=True,
            )

        if tStep is not None:
            if tEnd is not None:
                tectime = tecStart + tStep
                while tectime < tEnd:
                    tmpTec = []
                    tmpTec.insert(0, {"start": tectime, "tMap": tMap, "zMap": zMap})
                    tecdata = pd.concat(
                        [
                            tecdata,
                            pd.DataFrame(tmpTec, columns=["start", "tMap", "zMap"]),
                        ],
                        ignore_index=True,
                    )
                    tectime = tectime + tStep

        return tecdata

    def _defineTectonic(self, k, tecSort, tecdata):

        tecStart = None
        tEnd = None
        tStep = None
        tMap = None
        zMap = None

        try:
            tecStart = tecSort[k]["start"]
        except Exception:
            print("For each tectonic event a start time is required.")
            raise ValueError("Tectonic event {} has no parameter start".format(k))
        try:
            tMap = tecSort[k]["mapH"] + ".npz"
        except Exception:
            pass
        try:
            zMap = tecSort[k]["mapV"] + ".npz"
        except Exception:
            pass
        try:
            tStep = self.tecStep
        except Exception:
            pass
        try:
            tEnd = tecSort[k]["end"]
        except Exception:
            pass

        tecdata = self._storeTectonic(k, tecStart, zMap, tMap, tStep, tEnd, tecdata)

        return tecdata

    def _readTectonic(self):
        """
        Parse tectonic forcing conditions.
        """

        tecdata = None
        try:
            tecDict = self.input["tectonic"]
            tecSort = sorted(tecDict, key=itemgetter("start"))
            for k in range(len(tecSort)):
                tecdata = self._defineTectonic(k, tecSort, tecdata)

            if tecdata["start"][0] > self.tStart:
                tmpTec = []
                tmpTec.insert(
                    0, {"start": self.tStart, "tMap": "empty", "zMap": "empty"}
                )
                tecdata = pd.concat(
                    [pd.DataFrame(tmpTec, columns=["start", "tMap", "zMap"]), tecdata],
                    ignore_index=True,
                )
            self.tecdata = tecdata

        except KeyError:
            self.tecdata = None
            pass

        return

    def _defineRain(self, k, rStart, rMap, rUniform, raindata):

        if rMap is None and rUniform is None:
            print(
                "For each climate event a rainfall value (uniform) or a rainfall grid (map) is required."
            )
            raise ValueError(
                "Climate event {} has no rainfall value (uniform) or a rainfall map (map).".format(
                    k
                )
            )

        tmpRain = []
        if rMap is None:
            tmpRain.insert(
                0, {"start": rStart, "rUni": rUniform, "rMap": None, "rKey": None},
            )
        else:
            tmpRain.insert(
                0,
                {
                    "start": rStart,
                    "rUni": None,
                    "rMap": rMap[0] + ".npz",
                    "rKey": rMap[1],
                },
            )

        if k == 0:
            raindata = pd.DataFrame(tmpRain, columns=["start", "rUni", "rMap", "rKey"])
        else:
            raindata = pd.concat(
                [
                    raindata,
                    pd.DataFrame(tmpRain, columns=["start", "rUni", "rMap", "rKey"]),
                ],
                ignore_index=True,
            )

        return raindata

    def _readRain(self):
        """
        Parse rain forcing conditions.
        """

        raindata = None
        try:
            rainDict = self.input["climate"]
            rainSort = sorted(rainDict, key=itemgetter("start"))
            for k in range(len(rainSort)):
                rStart = None
                rUniform = None
                rMap = None
                try:
                    rStart = rainSort[k]["start"]
                except Exception:
                    print("For each climate event a start time is required.")
                    raise ValueError(
                        "Climate event {} has no parameter start".format(k)
                    )
                try:
                    rUniform = rainSort[k]["uniform"]
                except Exception:
                    pass
                try:
                    rMap = rainSort[k]["map"]
                except Exception:
                    pass

                if rMap is not None:
                    if self.meshFile != rMap[0] + ".npz":
                        try:
                            with open(rMap[0] + ".npz") as rainfile:
                                rainfile.close()
                                pass
                        except IOError:
                            print("Unable to open rain file: ", rMap[0] + ".npz")
                            raise IOError(
                                "The rain file {} is not found for climatic event {}.".format(
                                    rMap[0] + ".npz", k
                                )
                            )
                        mdata = np.load(rMap[0] + ".npz")
                        rainSet = mdata.files
                    else:
                        mdata = np.load(self.meshFile)
                        rainSet = mdata.files
                    try:
                        rainKey = mdata[rMap[1]]
                        if rainKey is not None:
                            pass
                    except KeyError:
                        print(
                            "Field name {} is missing from rain file {}".format(
                                rMap[1], rMap[0] + ".npz"
                            )
                        )
                        print("The following fields are available: ", rainSet)
                        print("Check your rain file fields definition...")
                        raise KeyError(
                            "Field name for rainfall is not defined correctly or does not exist!"
                        )

                    raindata = self._defineRain(k, rStart, rMap, rUniform, raindata)

            if raindata["start"][0] > self.tStart:
                tmpRain = []
                tmpRain.insert(
                    0, {"start": self.tStart, "rUni": 0.0, "rMap": None, "rKey": None}
                )
                raindata = pd.concat(
                    [
                        pd.DataFrame(
                            tmpRain, columns=["start", "rUni", "rMap", "rKey"]
                        ),
                        raindata,
                    ],
                    ignore_index=True,
                )
            self.raindata = raindata

        except KeyError:
            self.raindata = None
            pass

        return

    def _readPaleo(self):
        """
        Parse paleomap conditions.
        """
        try:
            paleoDict = self.input["paleomap"]
            paleoSort = sorted(paleoDict, key=itemgetter("time"))
            for k in range(len(paleoSort)):
                pTime = None
                pMap = None
                try:
                    pTime = paleoSort[k]["time"]
                except Exception:
                    print("For each paleomap a given time is required.")
                    raise ValueError("Paleomap {} has no parameter time".format(k))
                try:
                    pMap = paleoSort[k]["npdata"]
                except Exception:
                    pass

                if pMap is not None:

                    try:
                        with open(pMap + ".npz") as meshfile:
                            meshfile.close()
                            pass
                    except IOError:
                        print("Unable to open numpy dataset: ", pMap + ".npz")
                        raise IOError("The numpy dataset is not found...")

                tmpPaleo = []
                tmpPaleo.insert(0, {"time": pTime, "pMap": pMap + ".npz"})

                if k == 0:
                    paleodata = pd.DataFrame(tmpPaleo, columns=["time", "pMap"])
                else:
                    paleodata = pd.concat(
                        [paleodata, pd.DataFrame(tmpPaleo, columns=["time", "pMap"])],
                        ignore_index=True,
                    )

            self.paleodata = paleodata
            self.paleoNb = len(paleodata)

        except KeyError:
            self.paleodata = None
            self.paleoNb = 0
            pass

        return

    def _readForcePaleo(self):
        """
        Parse paleomap forcing.
        """

        try:
            fpaleoDict = self.input["forcepaleo"]

            try:
                self.forceDir = fpaleoDict["dir"]
                if not os.path.exists(self.forceDir):
                    print("Forcing paleo directory does not exist!")
                    raise ValueError("Forcing paleo directory does not exist!")

                if self.tout > self.tecStep:
                    self.tout = self.tecStep
                    print(
                        "Output time interval and tectonic forcing time step \
                          have been adjusted to match each others."
                    )
                elif self.tout < self.tecStep:
                    self.tecStep = self.tout
                    print(
                        "Output time interval and tectonic forcing time step \
                          have been adjusted to match each others."
                    )

                out_nb = int((self.tEnd - self.tStart) / self.tout) + 1
                stepf = np.arange(1, out_nb, dtype=int)
                self.stepb = np.flip(np.arange(0, out_nb - 1, dtype=int))
                self.alpha = stepf.astype(float) / (out_nb - 1)
                self.forceStep = 0
            except Exception:
                print("A directory is required to force the model with paleodata.")
                raise ValueError("forcepaleo key requires a directory")

        except KeyError:
            self.forceDir = None
            self.forceStep = -1
            pass

        return

    def _defineFlex(self, flexDict):
        """
        Define isostatic flexure conditions.
        """

        flexSort = sorted(flexDict, key=itemgetter("time"))
        for k in range(len(flexSort)):
            fTime = None
            fMap = None
            fKey = None
            try:
                fTime = flexSort[k]["time"]
            except Exception:
                print("For each elastic thickness map a given time is required.")
                raise ValueError("elastic thickness {} has no parameter time".format(k))
            try:
                fMap = flexSort[k]["npdata"]
            except Exception:
                pass

            if fMap is not None:
                try:
                    with open(fMap + ".npz") as meshfile:
                        meshfile.close()
                        pass
                except IOError:
                    print("Unable to open numpy dataset: ", fMap + ".npz")
                    raise IOError("The numpy dataset is not found...")

            try:
                fKey = flexSort[k]["key"]
            except Exception:
                pass

            tmpFlex = []
            tmpFlex.insert(0, {"time": fTime, "fMap": fMap + ".npz", "fKey": fKey})

            if k == 0:
                flexdata = pd.DataFrame(tmpFlex, columns=["time", "fMap", "fKey"])
            else:
                flexdata = pd.concat(
                    [
                        flexdata,
                        pd.DataFrame(tmpFlex, columns=["time", "fMap", "fKey"]),
                    ],
                    ignore_index=True,
                )
        self.flexdata = flexdata

        return

    def _readFlex(self):
        """
        Parse isostatic flexure conditions.
        """
        try:
            flexDict = self.input["flexure"]
            self.flexure = True

            try:
                self.dmantle = flexDict["dmantle"]
            except KeyError:
                self.dmantle = 3350.0

            try:
                self.dfill = flexDict["dfill"]
            except KeyError:
                self.dfill = 2750.0

            try:
                self.young = flexDict["young"]
            except KeyError:
                self.young = 1.0e11

            try:
                self.fpts = flexDict["npts"]
            except KeyError:
                print(
                    "The number of points under the influence of a given point load is required."
                )
                raise ValueError("The nbPts parameter needs to be be defined")

            self._defineFlex(flexDict)

        except KeyError:
            self.flexure = False
            pass

        return

    def _readOut(self):
        """
        Parse output directory.
        """

        try:
            outDict = self.input["output"]
            try:
                self.outputDir = outDict["dir"]
            except KeyError:
                self.outputDir = "output"
            try:
                self.makedir = outDict["makedir"]
            except KeyError:
                self.makedir = True
        except KeyError:
            self.outputDir = "output"
            self.makedir = True

        if self.rStep > 0:
            self.makedir = False

        return

    def _defineFit(self, paleoDict):
        """
        Define paleotopography files parameters.
        """

        try:
            self.paleoTopo = paleoDict["paleotopo"]
            self.paleoTopo = self.paleoTopo + ".npz"
        except KeyError:
            print(
                "Key 'paleotopo' is required and is missing in the 'paleofit' declaration!"
            )
            raise KeyError("NPZ paleotopography file is not defined!")
        try:
            with open(self.paleoTopo) as paleonpz:
                paleonpz.close()
                pass
        except IOError:
            print("Unable to open NPZ paleo file: ", self.paleoTopo)
            raise IOError("The NPZ paleotopography file is not found...")

        try:
            self.lonlat = paleoDict["lonlat"]
            self.lonlat = self.lonlat + ".npz"
        except KeyError:
            print(
                "Key 'lonlat' is required and is missing in the 'paleofit' declaration!"
            )
            raise KeyError("Longitude/latitude coordinates file is not defined!")
        try:
            with open(self.paleoTopo) as paleonpz:
                paleonpz.close()
                pass
        except IOError:
            print("Unable to open lonlat file: ", self.lonlat)
            raise IOError("The lonlat paleotopography file is not found...")

        try:
            self.paleoDir = paleoDict["vdispdir"]
        except KeyError:
            print(
                "Key 'vdispdir' is required and is missing in the 'paleofit' declaration!"
            )
            raise KeyError(
                "Paleo-tectonic directory to store vertical displacement is not defined!"
            )

        try:
            os.makedirs(self.paleoDir, exist_ok=True)
        except FileExistsError:
            # Directory already exists
            pass

        return

    def _paleoFit(self):
        """
        Parse paleotopography convergence parameters.
        """

        try:
            paleoDict = self.input["paleofit"]
            try:
                self.paleostep = paleoDict["step"]
            except KeyError:
                self.paleostep = 10
            try:
                self.cvglimit = paleoDict["cvglimit"]
            except KeyError:
                self.cvglimit = 0.1
            try:
                self.erange = paleoDict["erange"]
            except KeyError:
                self.erange = 100.0

            self._defineFit(paleoDict)

        except KeyError:
            self.paleostep = 0
            self.cvglimit = 0.0
            self.erange = 100.0
            self.paleoTopo = None
            self.paleoDir = None

        return
