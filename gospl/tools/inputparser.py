import os
import sys
import petsc4py
import numpy as np
import pandas as pd
import ruamel.yaml as yaml

from operator import itemgetter
from scipy.interpolate import interp1d

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class ReadYaml(object):
    """
    Class for reading simulation input file and initialising model parameters.

    Definition of input parameters is provided in the `User Documentation <https://gospl.readthedocs.io/en/latest/inputfile.html>`_
    """

    def __init__(self, filename):
        """
        Parsing YAML file.

        :arg filename: input filename (.yml YAML file)
        """

        if self.showlog:
            self.log = petsc4py.PETSc.Log()
            self.log.begin()

        # Check input file exists
        try:
            with open(filename) as finput:
                pass
        except IOError:
            print("Unable to open file: ", filename, flush=True)
            raise IOError("The input file is not found...")

        # Open YAML file
        with open(filename, "r") as finput:
            self.input = yaml.load(finput, Loader=yaml.Loader)

        if MPIrank == 0 and "name" in self.input.keys() and self.verbose:
            print(
                "The following model will be run:     {}".format(self.input["name"]),
                flush=True,
            )

        # Read simulation parameters
        self._readDomain()
        self._readTime()
        self._readSPL()
        self._readHillslope()
        self._readSealevel()
        self._readTectonic()
        self._readRain()
        self._readCompaction()

        if self.raindata is not None:
            if self.rStep > 0:
                self.raindata = self.raindata[
                    self.raindata["start"] >= self.tStart + self.rStep * self.tout
                ]
                self.raindata.reset_index(drop=True, inplace=True)
                self.rainNb = len(self.raindata)

        self._readBackwardPaleo()
        self._readOut()
        self._readForcePaleo()

        self.radius = 6378137.0
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
            print(
                "Key 'domain' is required and is missing in the input file!", flush=True
            )
            raise KeyError("Key domain is required in the input file!")

        # try:
        #     self.radius = domainDict["radius"]
        # except KeyError:
        #     self.radius = 6378137.0

        try:
            self.flowDir = domainDict["flowdir"]
        except KeyError:
            self.flowDir = 6

        try:
            self.flowExp = domainDict["flowexp"]
        except KeyError:
            self.flowExp = 0.42

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
            with open(self.meshFile) as meshfile:
                meshfile.close()

        except IOError:
            print("Unable to open numpy dataset: {}".format(self.meshFile), flush=True)
            raise IOError("The numpy dataset is not found...")

        try:
            self.fast = domainDict["fast"]
        except KeyError:
            self.fast = False

        try:
            self.backward = domainDict["backward"]
        except KeyError:
            self.backward = False

        try:
            self.interp = domainDict["interp"]
        except KeyError:
            self.interp = 1

        self._extraDomain()

        return

    def _extraDomain(self):
        """
        Read domain additional information.
        """

        domainDict = self.input["domain"]

        try:
            self.overlap = domainDict["overlap"]
        except KeyError:
            self.overlap = 1

        try:
            dataFile = domainDict["nperodep"]
            self.dataFile = dataFile + ".npz"
            with open(self.dataFile) as dataFile:
                dataFile.close()
        except KeyError:
            self.dataFile = None

        try:
            strataFile = domainDict["npstrata"]
            self.strataFile = strataFile + ".npz"
            with open(self.strataFile) as strataFile:
                strataFile.close()
        except KeyError:
            self.strataFile = None

        return

    def _readTime(self):
        """
        Read simulation time declaration.
        """

        try:
            timeDict = self.input["time"]
        except KeyError:
            print(
                "Key 'time' is required and is missing in the input file!", flush=True
            )
            raise KeyError("Key time is required in the input file!")

        try:
            self.tStart = timeDict["start"]
        except KeyError:
            print(
                "Key 'start' is required and is missing in the 'time' declaration!",
                flush=True,
            )
            raise KeyError("Simulation start time needs to be declared.")

        try:
            self.tEnd = timeDict["end"]
        except KeyError:
            print(
                "Key 'end' is required and is missing in the 'time' declaration!",
                flush=True,
            )
            raise KeyError("Simulation end time needs to be declared.")

        if self.tEnd <= self.tStart:
            raise ValueError("Simulation end/start times do not make any sense!")

        try:
            self.dt = timeDict["dt"]
        except KeyError:
            print(
                "Key 'dt' is required and is missing in the 'time' declaration!",
                flush=True,
            )
            raise KeyError("Simulation discretisation time step needs to be declared.")

        try:
            self.tout = timeDict["tout"]
        except KeyError:
            self.tout = self.tEnd - self.tStart
            print(
                "Output time interval 'tout' has been set to {} years".format(
                    self.tout
                ),
                flush=True,
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
                ),
                flush=True,
            )

        try:
            self.tecStep = timeDict["tec"]
        except KeyError:
            self.tecStep = self.tout

        try:
            self.strat = timeDict["strat"]
        except KeyError:
            self.strat = 0

        if self.tout < self.tecStep:
            self.tecStep = self.tout
            print(
                "Output time interval and tectonic forcing time step \
                 have been adjusted to match each others.",
                flush=True,
            )

        if self.tout < self.strat:
            self.strat = self.tout
            print(
                "Output time interval and stratal time step \
                 have been adjusted to match each others.",
                flush=True,
            )

        if self.tecStep > 0:
            if self.tout % self.tecStep != 0:
                print(
                    "When declaring tectonic time interval, the value should be divisible by the output time interval.",
                    flush=True,
                )
                raise ValueError("Tectonic time interval definition is wrong!")

        if self.strat > 0:
            if self.tout % self.strat != 0:
                print(
                    "When declaring stratal time interval, the value should be divisible by the output time interval.",
                    flush=True,
                )
                raise ValueError("Stratal time interval definition is wrong!")
            self.stratNb = int((self.tEnd - self.tStart) / self.strat) + 1
        else:
            self.stratNb = 0

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
                    "When using the Surface Process Model definition of coefficient Kb is required.",
                    flush=True,
                )
                raise ValueError("Surface Process Model: Kb coefficient not found.")
            try:
                self.waterfill = splDict["wfill"]
            except KeyError:
                self.waterfill = 100.0
            try:
                self.sedfill = splDict["sfill"]
            except KeyError:
                self.sedfill = 50.0
            try:
                self.coeffd = splDict["d"]
            except KeyError:
                self.coeffd = 0.0
            # try:
            #     # `wght` is the percentage of upstream sediment flux
            #     # that will be deposited on each cell...
            #     self.wght = splDict["wght"]
            #     if self.wght >= 1.0:
            #         self.wght = 0.999
            # except KeyError:
            #     self.wght = 0.0

        except KeyError:
            self.K = 1.0e-12
            self.waterfill = 100.0
            self.sedfill = 50.0
            # self.wght = 0.0
            self.coeffd = 0.0

        return

    def _readHillslope(self):
        """
        Read hillslope parameters.
        """

        try:
            hillDict = self.input["diffusion"]

            try:
                self.Cda = hillDict["hillslopeKa"]
            except KeyError:
                print(
                    "When declaring diffusion processes, the coefficient hillslopeKa is required.",
                    flush=True,
                )
                raise ValueError(
                    "Hillslope: Cd coefficient for aerial environment not found."
                )
            try:
                self.Cdm = hillDict["hillslopeKm"]
            except KeyError:
                print(
                    "When declaring diffusion processes, the coefficient hillslopeKm is required.",
                    flush=True,
                )
                raise ValueError(
                    "Hillslope: Cd coefficient for marine environment not found."
                )
            try:
                self.sedimentK = hillDict["sedK"]
            except KeyError:
                self.sedimentK = 10000.0
            try:
                self.sedimentKf = hillDict["sedKf"]
            except KeyError:
                self.sedimentKf = 20000.0
            try:
                self.sedimentKw = hillDict["sedKw"]
            except KeyError:
                self.sedimentKw = 30000.0
            try:
                self.diffstep = hillDict["dstep"]
            except KeyError:
                self.diffstep = 3
        except KeyError:
            self.Cda = 0.0
            self.Cdm = 0.0
            self.sedimentK = 10000.0
            self.sedimentKf = 20000.0
            self.sedimentKw = 30000.0
            self.diffstep = 3

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
                            dtype=np.float64,
                            low_memory=False,
                        )

                    except ValueError:
                        try:
                            seadata = pd.read_csv(
                                seafile,
                                sep=r"\s+",
                                engine="c",
                                header=None,
                                na_filter=False,
                                dtype=np.float64,
                                low_memory=False,
                            )

                        except ValueError:
                            print(
                                "The sea-level file is not well formed: it should be comma or tab separated",
                                flush=True,
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
            self.seafunction = interp1d(seadata[0], seadata[1], kind="linear")
        else:
            year = np.linspace(self.tStart, self.tEnd + self.dt, num=11, endpoint=True)
            seaval = np.full(len(year), sealevel)
            self.seafunction = interp1d(year, seaval, kind="linear")

        return

    def _storeTectonic(self, k, tecStart, zMap, tMap, tStep, tEnd, tecdata):
        """
        Record tectonic conditions.

        :arg k: tectonic event number
        :arg tecStart: tectonic event start time
        :arg zMap: horizontal tectonic displacement file
        :arg tMap: vertical tectonic displacement file
        :arg tStep: tectonic time step
        :arg tEnd: tectonic event end time
        :arg tecdata: pandas dataframe storing each tectonic event

        :return: appended tecdata
        """

        if tMap is not None:
            if self.meshFile != tMap:
                try:
                    with open(tMap) as tecfile:
                        tecfile.close()

                except IOError:
                    print("Unable to open tectonic file: {}".format(tMap), flush=True)
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

                except IOError:
                    print("Unable to open tectonic file: {}".format(zMap), flush=True)
                    raise IOError(
                        "The tectonic file {} is not found for climatic event {}.".format(
                            zMap, k
                        )
                    )
        else:
            zMap = "empty"

        if tMap == "empty" and zMap == "empty":
            print(
                "For each tectonic event a tectonic grid (mapH or mapV) is required.",
                flush=True,
            )
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
        """
        Define tectonic conditions.

        :arg k: tectonic event number
        :arg tecSort: sorted tectonic event
        :arg tecdata: pandas dataframe storing each tectonic event
        :return: appended tecdata
        """

        tecStart = None
        tEnd = None
        tStep = None
        tMap = None
        zMap = None

        try:
            tecStart = tecSort[k]["start"]
        except Exception:
            print("For each tectonic event a start time is required.", flush=True)
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
            self.tecdata = tecdata[tecdata["start"] >= self.tStart]
            if self.rStep > 0:
                self.tecdata = tecdata[
                    tecdata["start"] >= self.tStart + self.rStep * self.tout
                ]
            self.tecdata.reset_index(drop=True, inplace=True)

        except KeyError:
            self.tecdata = None

        return

    def _defineRain(self, k, rStart, rMap, rUniform, raindata):
        """
        Define precipitation conditions.

        :arg k: precipitation event number
        :arg rStart: precipitation event start time
        :arg rMap: precipitation map file event
        :arg rUniform: precipitation uniform value event
        :arg raindata: pandas dataframe storing each precipitation event
        :return: appended raindata
        """

        if rMap is None and rUniform is None:
            print(
                "For each climate event a rainfall value (uniform) or a rainfall \
                grid (map) is required.",
                flush=True,
            )
            raise ValueError(
                "Climate event {} has no rainfall value (uniform) or a rainfall \
                map (map).".format(
                    k
                )
            )

        tmpRain = []
        if rMap is None:
            tmpRain.insert(
                0,
                {"start": rStart, "rUni": rUniform, "rMap": None, "rKey": None},
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
                    print(
                        "For each climate event a start time is required.", flush=True
                    )
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

                        except IOError:
                            print(
                                "Unable to open rain file: {}.npz".format(rMap[0]),
                                flush=True,
                            )
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
                            "Field name {} is missing from rain file {}.npz".format(
                                rMap[1], rMap[0]
                            ),
                            flush=True,
                        )
                        print(
                            "The following fields are available: {}".format(rainSet),
                            flush=True,
                        )
                        print("Check your rain file fields definition...", flush=True)
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
            self.raindata = raindata[raindata["start"] >= self.tStart]
            self.raindata.reset_index(drop=True, inplace=True)
            self.rainNb = len(self.raindata)

        except KeyError:
            self.raindata = None

        return

    def _readBackwardPaleo(self):
        """
        Force model with backward paleomaps.
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
                    print("For each paleomap a given time is required.", flush=True)
                    raise ValueError("Paleomap {} has no parameter time".format(k))
                try:
                    pMap = paleoSort[k]["npdata"]
                except Exception:
                    pass

                if pMap is not None:
                    try:
                        with open(pMap + ".npz") as meshfile:
                            meshfile.close()

                    except IOError:
                        print(
                            "Unable to open numpy dataset: {}.npz".format(pMap),
                            flush=True,
                        )
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

            self.paleodata = paleodata[paleodata["time"] >= self.tStart]
            if self.rStep > 0:
                self.paleodata = paleodata[
                    paleodata["start"] >= self.tStart + self.rStep * self.tout
                ]

            self.paleodata.reset_index(drop=True, inplace=True)
            self.paleoNb = len(self.paleodata)

        except KeyError:
            self.paleodata = None
            self.paleoNb = 0

        return

    def _readForcePaleo(self):
        """
        Get series of paleomaps to force the model through time.
        """

        try:
            fpaleoDict = self.input["forcepaleo"]

            try:
                self.forceDir = fpaleoDict["dir"]
                if not os.path.exists(self.forceDir):
                    print("Forcing paleo directory does not exist!", flush=True)
                    raise ValueError("Forcing paleo directory does not exist!")
                try:
                    forceNb = fpaleoDict["steps"]
                except KeyError:
                    print(
                        "New Paleomap loading frequency 'steps' is required",
                        flush=True,
                    )
                outNb = np.sum(forceNb)
                self.stepb = np.flip(np.arange(0, outNb, dtype=int))
                self.alpha = np.zeros(len(self.stepb))
                p = 0
                for k in range(len(forceNb)):
                    out_nb = forceNb[k] + 1
                    stepf = np.arange(1, out_nb, dtype=int)
                    self.alpha[p : p + len(stepf)] = stepf.astype(float) / (out_nb - 1)
                    p += len(stepf)
                self.forceStep = 0

                if self.rStep > 0:
                    steptime = np.arange(0, outNb, dtype=int) * self.tout + self.tStart
                    id = np.where(steptime == self.tStart + self.rStep * self.tout)[0]
                    if len(id) == 0:
                        raise steptime(
                            "Something went wrong with the restart time, it needs to be related to the forcing step."
                        )
                    self.forceStep = id[0]

            except Exception:
                print(
                    "A directory is required to force the model with paleodata.",
                    flush=True,
                )
                raise ValueError("forcepaleo key requires a directory")

        except KeyError:
            self.forceDir = None
            self.forceStep = -1

        return

    def _readCompaction(self):
        """
        Read compaction parameters.
        """

        try:
            compDict = self.input["compaction"]
            try:
                self.phi0s = compDict["phis"]
            except KeyError:
                self.phi0s = 0.49

            try:
                self.phi0f = compDict["phif"]
            except KeyError:
                self.phi0f = 0.63

            try:
                self.phi0w = compDict["phiw"]
            except KeyError:
                self.phi0w = 0.63

            try:
                self.phi0c = compDict["phic"]
            except KeyError:
                self.phi0c = 0.65

        except KeyError:
            self.phi0s = 0.49
            self.phi0f = 0.63
            self.phi0w = 0.63
            self.phi0c = 0.65

        self._extraCompaction()

        return

    def _extraCompaction(self):
        """
        Read compaction additional parameters.
        """

        try:
            compDict = self.input["compaction"]
            try:
                self.z0s = compDict["z0s"]
            except KeyError:
                self.z0s = 3700.0

            try:
                self.z0f = compDict["z0f"]
            except KeyError:
                self.z0f = 1960.0

            try:
                self.z0w = compDict["z0w"]
            except KeyError:
                self.z0w = 1960.0

            try:
                self.z0c = compDict["z0c"]
            except KeyError:
                self.z0c = 2570.0

        except KeyError:
            self.z0s = 3700.0
            self.z0f = 1960.0
            self.z0w = 1960.0
            self.z0c = 2570.0

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
