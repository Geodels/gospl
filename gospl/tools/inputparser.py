import os
import sys
import petsc4py
import numpy as np
import pandas as pd

from operator import itemgetter

if "READTHEDOCS" not in os.environ:
    from ruamel.yaml import YAML
    from scipy.interpolate import interp1d

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
        self.finput = filename
        try:
            with open(filename) as finput:
                pass
        except IOError:
            print("Unable to open file: ", filename, flush=True)
            raise IOError("The input file is not found...")

        # Open YAML file
        with open(filename, "r") as finput:
            yaml = YAML(typ='rt')
            self.input = yaml.load(finput)

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
        self._readSoilInfo()
        self._readSealevel()
        self._readTectonics()
        self._readErofactor()
        self._readRain()
        self._readCompaction()
        self._readIce()
        self._readOrography()
        self._readFlex()
        self._readTeMap()
        self._readOut()

        self.tNow = self.tStart
        self.saveTime = self.tNow
        if self.strat > 0:
            self.saveStrat = self.tNow + self.strat
        else:
            self.saveStrat = self.tEnd + self.tout

        # In case of restarting simulation
        self._restartUpdate()

        return

    def _get_param(self, *keys, default=None):
        """
        Safe traversal into ``self.input``.

        Returns ``default`` (None by default) if any key in the chain is
        missing, or if an intermediate node is not a dict. Never raises
        KeyError.

        Examples::

            self._get_param("name")                  # self.input.get("name", default)
            self._get_param("domain", "flowdir")     # self.input["domain"].get("flowdir", default)
            self._get_param("spl", "K", default=1e-12)

        Use this for direct access through ``self.input``. For inner
        keys on a pre-extracted section dict (e.g. ``splDict =
        self.input["spl"]`` bound on an earlier line), prefer the
        Python builtin ``splDict.get(key, default)`` — it makes the
        data flow more visible and avoids a redundant top-level lookup.

        Out of scope: this helper navigates ``self.input`` only. For
        NPZ-archive accesses elsewhere in this file (``_isKeyinFile``,
        the rainKey/sedKey blocks inside ``_readRain``/``_readErofactor``/
        ``_readTeMap``), keep the existing ``try/except KeyError`` since
        those exceptions carry user-facing diagnostics about missing
        fields in the data file.

        See AGENTS.md > The ``_extra*`` methods are mandatory
        continuations: the call chain (``_readDomain → _extraDomain →
        _extraDomain2``, ``_readHillslope → _extraHillslope``, etc.) is
        load-bearing. The internal ``try/except KeyError`` blocks
        inside those methods are NOT load-bearing and use this helper.
        """
        node = self.input
        for k in keys:
            if not isinstance(node, dict) or k not in node:
                return default
            node = node[k]
        return node

    def _restartUpdate(self):
        """
        Update some forcing parameters in case of a restart.
        """

        if self.rStep > 0:
            rNow = self.tStart + self.rStep * self.tout

            if self.raindata is not None:
                for k in range(len(self.raindata)):
                    if self.raindata["start"][k] < rNow and k < len(self.raindata) - 1:
                        self.raindata.loc[k, ["start"]] = rNow - self.tout
                    elif (
                        self.raindata["start"][k] < rNow and k == len(self.raindata) - 1
                    ):
                        self.raindata.loc[k, ["start"]] = rNow
                self.raindata = self.raindata[self.raindata["start"] >= rNow]
                self.raindata.reset_index(drop=True, inplace=True)
                self.rainNb = len(self.raindata)

            if self.tedata is not None:
                for k in range(len(self.tedata)):
                    if self.tedata["start"][k] < rNow and k < len(self.tedata) - 1:
                        self.tedata.loc[k, ["start"]] = rNow - self.tout
                    elif (
                        self.tedata["start"][k] < rNow and k == len(self.tedata) - 1
                    ):
                        self.tedata.loc[k, ["start"]] = rNow
                self.tedata = self.tedata[self.tedata["start"] >= rNow]
                self.tedata.reset_index(drop=True, inplace=True)
                self.teNb = len(self.tedata)

            if self.sedfacdata is not None:
                for k in range(len(self.sedfacdata)):
                    if self.sedfacdata["start"][k] < rNow and k < len(self.sedfacdata) - 1:
                        self.sedfacdata.loc[k, ["start"]] = rNow - self.tout
                    elif (
                        self.sedfacdata["start"][k] < rNow and k == len(self.sedfacdata) - 1
                    ):
                        self.sedfacdata.loc[k, ["start"]] = rNow
                self.sedfacdata = self.sedfacdata[self.sedfacdata["start"] >= rNow]
                self.sedfacdata.reset_index(drop=True, inplace=True)
                self.sedfactNb = len(self.sedfacdata)

        return

    def _readDomain(self):
        """
        Read domain definition, boundary conditions and flow direction parameters.
        """

        # TODO-REFACTOR: complex except, needs manual review
        try:
            domainDict = self.input["domain"]
        except KeyError:
            print(
                "Key 'domain' is required and is missing in the input file!", flush=True
            )
            raise KeyError("Key domain is required in the input file!")

        self.flowDir = domainDict.get("flowdir", 8)

        # try:
        #     self.flowExp = domainDict["flowexp"]
        # except KeyError:
        self.flowExp = 1.1

        self.boundCond = domainDict.get("bc", '1111')

        # TODO-REFACTOR: complex except, needs manual review
        try:
            meshInfo = domainDict["npdata"]
        except KeyError:
            print(
                "Key 'npdata' is required and is missing in the 'domain' declaration!",
                flush=True,
            )
            raise KeyError("Compressed numpy dataset definition is not defined!")

        self.meshFile = meshInfo[0] + ".npz"
        self.infoCoords = meshInfo[1]
        self.infoCells = meshInfo[2]
        self.infoElev = meshInfo[3]

        try:
            with open(self.meshFile) as meshfile:
                meshfile.close()

        except IOError:
            print("Unable to open numpy dataset: {}".format(self.meshFile), flush=True)
            raise IOError("The numpy dataset is not found...")

        self.fast = domainDict.get("fast", False)

        self._extraDomain()

        return

    def _extraDomain(self):
        """
        Read domain additional information.
        """

        domainDict = self.input["domain"]

        self.seaDepo = domainDict.get("seadepo", True)
        self.overlap = domainDict.get("overlap", 1)

        nperodep = domainDict.get("nperodep")
        if nperodep is None:
            self.dataFile = None
        else:
            self.dataFile = nperodep + ".npz"
            with open(self.dataFile) as fh:
                pass

        self.nodep = domainDict.get("nodep", False)

        npstrata = domainDict.get("npstrata")
        if npstrata is None:
            self.strataFile = None
        else:
            self.strataFile = npstrata + ".npz"
            with open(self.strataFile) as fh:
                pass

        self._extraDomain2()

        return

    def _extraDomain2(self):
        """
        Read domain additional information.
        """

        domainDict = self.input["domain"]

        advscheme = domainDict.get("advect")
        if advscheme is None:
            self.advscheme = 1
        elif advscheme == 'iioe1':
            self.advscheme = 2
        elif advscheme == 'iioe2':
            self.advscheme = 3
        elif advscheme == 'upwind':
            self.advscheme = 1
        elif advscheme == 'interp':
            self.advscheme = 0
        else:
            raise ValueError(
                "Unknown advect scheme '%s'; expected one of "
                "iioe1 / iioe2 / upwind / interp." % advscheme
            )

        self.radius = domainDict.get("radius", 6378137.0)
        self.gravity = domainDict.get("gravity", 9.81)

        return

    def _readTime(self):
        """
        Read simulation time declaration.
        """

        # TODO-REFACTOR: complex except, needs manual review
        try:
            timeDict = self.input["time"]
        except KeyError:
            print(
                "Key 'time' is required and is missing in the input file!", flush=True
            )
            raise KeyError("Key time is required in the input file!")

        # TODO-REFACTOR: complex except, needs manual review
        try:
            self.tStart = timeDict["start"]
        except KeyError:
            print(
                "Key 'start' is required and is missing in the 'time' declaration!",
                flush=True,
            )
            raise KeyError("Simulation start time needs to be declared.")

        # TODO-REFACTOR: complex except, needs manual review
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

        # TODO-REFACTOR: complex except, needs manual review
        try:
            self.dt = timeDict["dt"]
        except KeyError:
            print(
                "Key 'dt' is required and is missing in the 'time' declaration!",
                flush=True,
            )
            raise KeyError("Simulation discretisation time step needs to be declared.")

        # TODO-REFACTOR: complex except, needs manual review (default + print side effect)
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

        self.rStep = timeDict.get("rstep", 0)

        if self.tStart + self.tout > self.tEnd:
            self.tout = self.tEnd - self.tStart
            print(
                "Output time interval was changed to {} years to match the end time".format(
                    self.tout
                ),
                flush=True,
            )

        if self.tout < self.dt:
            self.tout = self.dt
            print(
                "Output time interval was changed to {} years to match the time step dt".format(
                    self.dt
                ),
                flush=True,
            )

        self.strat = timeDict.get("strat", 0)

        if self.tout < self.strat:
            self.strat = self.tout
            print(
                "Output time interval and stratal time step \
                 have been adjusted to match each others.",
                flush=True,
            )

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
        Read surface processes erosion and deposition laws parameters.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section multi-default on missing "spl")
        try:
            splDict = self.input["spl"]
            # TODO-REFACTOR: complex except, needs manual review (K is required inside spl)
            try:
                self.K = splDict["K"]
            except KeyError:
                print(
                    "When using the Surface Process Model definition of coefficient K is required.",
                    flush=True,
                )
                raise ValueError("Surface Process Model: K coefficient not found.")
            self.coeffd = splDict.get("d", 0.0)
            self.fDepa = splDict.get("G", 0.0)
            self.spl_m = splDict.get("m", 0.5)
            self.spl_n = splDict.get("n", 1.0)
        except KeyError:
            self.K = 1.0e-12
            self.coeffd = 0.0
            self.fDepa = 0.0
            self.spl_m = 0.5
            self.spl_n = 1.0

        return

    def _readHillslope(self):
        """
        Read hillslope parameters.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section multi-default on missing "diffusion")
        try:
            hillDict = self.input["diffusion"]

            # TODO-REFACTOR: complex except, needs manual review (hillslopeKa is required inside diffusion)
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
            # TODO-REFACTOR: complex except, needs manual review (hillslopeKm is required inside diffusion)
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
            self.K_nl = hillDict.get("hillslopenl", 1.0)
            self.K_sc = hillDict.get("hillslopeSc", 0.0)
            self.K_nb = int(hillDict.get("hillslopeNb", 0))
            self.oFill = hillDict.get("oFill", -6000.0)

        except KeyError:
            self.Cda = 0.0
            self.Cdm = 0.0
            self.K_nl = 1.0
            self.K_sc = 0.0
            self.K_nb = 0
            self.oFill = -6000.0

        self._extraHillslope()

        return

    def _extraHillslope(self):
        """
        Read extra hillslope parameters.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section multi-default on missing "diffusion")
        try:
            hillDict = self.input["diffusion"]
            self.nlK = hillDict.get("nonlinKm", 10.0)
            self.clinSlp = hillDict.get("clinSlp", 1.0e-6)
            self.tsStep = hillDict.get("tsSteps", 2000)
            self.Gmar = hillDict.get("Gmar", 0.0)
            self.offshore = hillDict.get("offshore", 100.e5)
            # Lake / depression non-linear diffusion thresholds: only run the
            # expensive marine-style nonlinear diffusion in pits that are both
            # large in volume AND deep, otherwise fall back to bottom-up fill.
            self.nl_pit_volume = hillDict.get("nlPitVolume", 1.0e9)   # 1 km^3
            self.nl_pit_depth = hillDict.get("nlPitDepth", 100.0)     # 100 m
            self.nl_pit_K = hillDict.get("nlPitK", self.nlK)
            # Fraction of each pit's deposit concentrated at the inlets
            # (delta seed); the remainder is distributed as a bathymetric
            # bottom-up baseline. 0.0 = pure bowl fill, 1.0 = original
            # inlet-only spike.
            self.nl_pit_inlet_bias = hillDict.get("pitInletBias", 0.10)
        except KeyError:
            self.nlK = 10.0
            self.clinSlp = 1.0e-6
            self.Gmar = 0.
            self.tsStep = 2000
            self.offshore = 100.e5
            self.nl_pit_volume = 1.0e9
            self.nl_pit_depth = 100.0
            self.nl_pit_K = self.nlK
            self.nl_pit_inlet_bias = 0.50

        self.clinSlp = max(1.0e-6, self.clinSlp)
        self.nl_pit_inlet_bias = min(1.0, max(0.0, self.nl_pit_inlet_bias))

        return

    def _readSoilInfo(self):
        """
        Read soil information parameters.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section multi-default on missing "soil" + cptSoil flag)
        try:
            soilDict = self.input["soil"]
            self.cptSoil = True
            # TODO-REFACTOR: complex except, needs manual review (soilK is required inside soil)
            try:
                self.Ksoil = soilDict["soilK"]
            except KeyError:
                print(
                    "When declaring soil production, the erodibility soilK is required.",
                    flush=True,
                )
                raise ValueError(
                    "Soil: Erodibility coefficient for soil not found."
                )
            # soil production maximum rate (50 m/Myr)
            self.P0 = soilDict.get("maxProd", 50.e-6)
            # soil production decay depth
            self.Hs = soilDict.get("depthProd", 0.0)
            # roughness length_scale
            self.h_star = float(soilDict.get("roughnessL", 1.0))
            # soil transport decay depth for diffusion
            self.H0 = soilDict.get("decayDepth", 0.7)
            # soil / bedrock transition limit ratio factor of production
            self.Sperc = soilDict.get("bedrockConv", 0.0001)
            # initial soil thickness
            self.cstSoilH = soilDict.get("uniform", 1.0)
            # TODO-REFACTOR: complex except, needs manual review (except sets BOTH local soilfile and self.soilFile)
            try:
                soilfile = soilDict["soilMap"]
            except KeyError:
                soilfile = None
                self.soilFile = None
            # TODO-REFACTOR: complex except, needs manual review (except sets BOTH local tempfile and self.tempFile)
            try:
                tempfile = soilDict["tempMap"]
            except KeyError:
                tempfile = None
                self.tempFile = None

            if soilfile is not None:
                self.soilFile = soilfile[0] + ".npz"
                self.soilData = soilfile[1]
                try:
                    with open(self.soilFile) as sinfo:
                        sinfo.close()
                except IOError:
                    print("Unable to open numpy dataset: {}".format(self.soilFile), flush=True)
                    raise IOError("The numpy dataset is not found...")
            
            if tempfile is not None:
                self.tempFile = tempfile[0] + ".npz"
                self.tempData = tempfile[1]
                try:
                    with open(self.tempFile) as tinfo:
                        tinfo.close()
                except IOError:
                    print("Unable to open numpy dataset: {}".format(self.tempFile), flush=True)
                    raise IOError("The numpy dataset is not found...")
            # Activation energy (J/mol)
            self.energyAct = soilDict.get("activation", 40.e3)
            # Reference temperature in Celsius
            self.tempRef = soilDict.get("tempRef", 15.0)

        except KeyError:
            self.cptSoil = False
            self.Ksoil = 0.0
            self.P0 = 0.0
            self.Hs = 0.0
            self.h_star = 1.0
            self.H0 = 1.0
            self.Sperc = 0.0
            self.cstSoilH = 0.0
            self.soilFile = None
            self.soilData = None
            self.energyAct = 0.0
            self.tempRef = 15.0
            self.tempFile = None
            self.tempData = None

        return

    def _readSealevel(self):
        """
        Define sealevel evolution.
        """

        seafile = None
        sealevel = 0.0
        self.seafunction = None
        # TODO-REFACTOR: complex except, needs manual review (nested position/curve fallback logic)
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

    def _isKeyinFile(self, dmap):
        '''
        Check if a numpy compressed file exists and that the corresponding keys are present in it.
        '''

        if dmap is not None:
            try:
                with open(dmap[0] + ".npz") as file:
                    file.close()
            except IOError:
                print(
                    "Unable to open file: {}.npz".format(dmap[0]),
                    flush=True,
                )
                raise IOError(
                    "The following file {} is not found.".format(
                        dmap[0] + ".npz"
                    )
                )
            mdata = np.load(dmap[0] + ".npz")

            # TODO-REFACTOR: complex except, needs manual review (out of scope: accesses np.load() archive, not self.input; carries a user-facing diagnostic about missing NPZ field)
            try:
                vkey = mdata[dmap[1]]
                if vkey is not None:
                    pass
            except KeyError:
                print(
                    "Field name {} is missing from file {}.npz".format(
                        dmap[1], dmap[0]
                    ),
                    flush=True,
                )
            del mdata

        return

    def _storeTectonics(self, k, tecStart, hMap, tMap, zMap, tecEnd, tecdata):
        """
        Record tectonic conditions.

        :arg k: tectonic event number
        :arg tecStart: tectonic event start time
        :arg hMap: horizontal tectonic information
        :arg tMap: vertical tectonic displacement information
        :arg zMap: elevation information
        :arg tEnd: tectonic event end time
        :arg tecdata: pandas dataframe storing each tectonic event

        :return: appended tecdata
        """

        if tMap is None:
            tMap = "empty"
        if hMap is None:
            hMap = "empty"
        if zMap is None:
            zMap = "empty"

        tmpTec = []
        tmpTec.insert(0, {"start": tecStart, "end": tecEnd, "tMap": tMap, "zMap": zMap, "hMap": hMap})

        if k == 0:
            tecdata = pd.DataFrame(tmpTec, columns=["start", "end", "tMap", "zMap", "hMap"])
        else:
            tecdata = pd.concat(
                [tecdata, pd.DataFrame(tmpTec, columns=["start", "end", "tMap", "zMap", "hMap"])],
                ignore_index=True,
            )

        return tecdata

    def _defineTectonics(self, k, tecSort, tecdata):
        """
        Define tectonics conditions.

        :arg k: tectonic event number
        :arg tecSort: sorted tectonic event
        :arg tecdata: pandas dataframe storing each tectonic event
        :return: appended tecdata
        """

        tecStart = None
        tecEnd = None
        zMap = None
        tMap = None
        hMap = None

        # TODO-REFACTOR: complex except, needs manual review (tecStart is required)
        try:
            tecStart = tecSort[k]["start"]
        except KeyError:
            print("For each tectonic event a start time is required.", flush=True)
            raise ValueError("Tectonic event {} has no parameter start".format(k))

        # TODO-REFACTOR: complex except, needs manual review (tecEnd is required)
        try:
            tecEnd = tecSort[k]["end"]
        except KeyError:
            print("For each tectonic event an end time is required.", flush=True)
            raise ValueError("Tectonic event {} has no parameter end".format(k))

        tMap = tecSort[k].get("upsub")
        self._isKeyinFile(tMap)

        hMap = tecSort[k].get("hdisp")
        self._isKeyinFile(hMap)

        zMap = tecSort[k].get("zfit")
        self._isKeyinFile(zMap)

        tecdata = self._storeTectonics(k, tecStart, hMap, tMap, zMap, tecEnd, tecdata)

        return tecdata

    def _readTectonics(self):
        """
        Parse tectonics forcing conditions.
        """

        tecdata = None
        # TODO-REFACTOR: complex except, needs manual review (outer-section: sets self.tecdata = None on missing "tectonics")
        try:
            tecDict = self.input["tectonics"]

            tecSort = sorted(tecDict, key=itemgetter("start"))
            for k in range(len(tecSort)):
                tecdata = self._defineTectonics(k, tecSort, tecdata)

            if tecdata["start"][0] > self.tStart:
                tmpTec = []
                tmpTec.insert(
                    0, {"start": self.tStart, "end": tecdata["start"][0], "tMap": "empty", "zMap": "empty", "hMap": "empty"}
                )
                tecdata = pd.concat(
                    [pd.DataFrame(tmpTec, columns=["start", "end", "tMap", "zMap", "hMap"]), tecdata],
                    ignore_index=True,
                )
            self.tecdata = tecdata[tecdata["start"] >= self.tStart]

            if self.rStep > 0:
                rNow = self.tStart + self.rStep * self.tout
                for k in range(len(self.tecdata)):
                    if self.tecdata["start"][k] < rNow and k < len(self.tecdata) - 1:
                        self.tecdata.loc[k, ["start"]] = rNow - self.tout
                    elif self.tecdata["start"][k] < rNow and k == len(self.tecdata) - 1:
                        self.tecdata.loc[k, ["start"]] = rNow
                self.tecdata = self.tecdata[self.tecdata["start"] >= rNow]
                self.tecdata.reset_index(drop=True, inplace=True)

        except KeyError:
            self.tecdata = None

        return

    def _defineErofactor(self, k, sStart, sMap, sUniform, sedfacdata):
        """
        Define sediment surface erodibility factor conditions.

        :arg k: erodibility factor map number
        :arg sStart: erodibility factor map start time
        :arg sMap: erodibility factor map file event
        :arg sUniform: erodibility factor uniform value event
        :arg sedfacdata: pandas dataframe storing each erodibility factor map
        :return: appended sedfacdata
        """

        if sMap is None and sUniform is None:
            print(
                "For each erodibility factor map a factor value (uniform) or a factor \
                grid (map) is required.",
                flush=True,
            )
            raise ValueError(
                "Sediment erodibility factor {} has no value (uniform) or a \
                map (map).".format(
                    k
                )
            )

        tmpErof = []
        if sMap is None:
            tmpErof.insert(
                0,
                {"start": sStart, "sUni": sUniform, "sMap": None, "sKey": None},
            )
        else:
            tmpErof.insert(
                0,
                {
                    "start": sStart,
                    "sUni": None,
                    "sMap": sMap[0] + ".npz",
                    "sKey": sMap[1],
                },
            )

        if k == 0:
            sedfacdata = pd.DataFrame(tmpErof, columns=["start", "sUni", "sMap", "sKey"])
        else:
            sedfacdata = pd.concat(
                [
                    sedfacdata,
                    pd.DataFrame(tmpErof, columns=["start", "sUni", "sMap", "sKey"]),
                ],
                ignore_index=True,
            )

        return sedfacdata

    def _readErofactor(self):
        """
        Parse erodibility factor based on surface geology.
        """

        sedfacdata = None
        # TODO-REFACTOR: complex except, needs manual review (outer-section: sets self.sedfactNb=0, self.sedfacdata=None on missing "sedfactor")
        try:
            sedDict = self.input["sedfactor"]
            sedSort = sorted(sedDict, key=itemgetter("start"))
            for k in range(len(sedSort)):
                sStart = None
                sUniform = None
                sMap = None
                # TODO-REFACTOR: complex except, needs manual review (sStart is required)
                try:
                    sStart = sedSort[k]["start"]
                except KeyError:
                    print(
                        "For each sediment factor a start time is required.", flush=True
                    )
                    raise ValueError(
                        "Sediment factor map {} has no parameter start".format(k)
                    )
                sUniform = sedSort[k].get("uniform")
                sMap = sedSort[k].get("map")

                if sMap is not None:
                    try:
                        with open(sMap[0] + ".npz") as sedfacfile:
                            sedfacfile.close()

                    except IOError:
                        print(
                            "Unable to open sediment factor file: {}.npz".format(sMap[0]),
                            flush=True,
                        )
                        raise IOError(
                            "The sediment factor file {} is not found for event {}.".format(
                                sMap[0] + ".npz", k
                            )
                        )
                    mdata = np.load(sMap[0] + ".npz")
                    sedfacSet = mdata.files

                    # TODO-REFACTOR: complex except, needs manual review (out of scope: accesses np.load() archive, not self.input; carries a user-facing diagnostic about missing NPZ field)
                    try:
                        sedKey = mdata[sMap[1]]
                        if sedKey is not None:
                            pass
                    except KeyError:
                        print(
                            "Field name {} is missing from sediment factor file {}.npz".format(
                                sMap[1], sMap[0]
                            ),
                            flush=True,
                        )
                        print(
                            "The following fields are available: {}".format(sedfacSet),
                            flush=True,
                        )
                        print("Check your sediment factor file fields definition...", flush=True)
                        raise KeyError(
                            "Field name for sediment factor is not defined correctly or does not exist!"
                        )

                sedfacdata = self._defineErofactor(k, sStart, sMap, sUniform, sedfacdata)

            if sedfacdata["start"][0] > self.tStart:
                tmpSedF = []
                tmpSedF.insert(
                    0, {"start": self.tStart, "sUni": 1.0, "sMap": None, "sKey": None}
                )
                sedfacdata = pd.concat(
                    [
                        pd.DataFrame(
                            tmpSedF, columns=["start", "sUni", "sMap", "sKey"]
                        ),
                        sedfacdata,
                    ],
                    ignore_index=True,
                )
            self.sedfacdata = sedfacdata.copy()
            self.sedfacdata.reset_index(drop=True, inplace=True)
            self.sedfactNb = len(self.sedfacdata)

        except KeyError:
            self.sedfactNb = 0
            self.sedfacdata = None

        return

    def _getTe(self, k, tStart, tMap, tUniform, tedata):
        """
        Define elastic map.

        :arg k: elastic map event number
        :arg teStart: elastic map event start time
        :arg teMap: elastic map file event
        :arg teUniform: elastic map uniform thickness value event
        :arg tedata: pandas dataframe storing each elastic map event
        :return: appended tedata
        """

        if tMap is None and tUniform is None:
            print(
                "For each elastic map a thickness value (uniform) or a elastic \
                grid (map) is required.",
                flush=True,
            )
            raise ValueError(
                "Elastic event {} has no thickness value (uniform) or a elastic \
                map (map).".format(
                    k
                )
            )

        tmpTe = []
        if tMap is None:
            tmpTe.insert(
                0,
                {"start": tStart, "tUni": tUniform, "tMap": None, "tKey": None},
            )
        else:
            tmpTe.insert(
                0,
                {
                    "start": tStart,
                    "tUni": 0.,
                    "tMap": tMap[0] + ".npz",
                    "tKey": tMap[1],
                },
            )

        if k == 0:
            tedata = pd.DataFrame(tmpTe, columns=["start", "tUni", "tMap", "tKey"])
        else:
            tedata = pd.concat(
                [
                    tedata,
                    pd.DataFrame(tmpTe, columns=["start", "tUni", "tMap", "tKey"]),
                ],
                ignore_index=True,
            )

        return tedata

    def _readTeMap(self):
        """
        Parse elastic map forcing conditions.
        """

        tedata = None
        # TODO-REFACTOR: complex except, needs manual review (outer-section: sets self.tedata = None on missing "temap")
        try:
            teDict = self.input["temap"]
            teSort = sorted(teDict, key=itemgetter("start"))
            for k in range(len(teSort)):
                rStart = None
                rUniform = None
                rMap = None
                # TODO-REFACTOR: complex except, needs manual review (rStart is required)
                try:
                    rStart = teSort[k]["start"]
                except KeyError:
                    print(
                        "For each elastic map event a start time is required.", flush=True
                    )
                    raise ValueError(
                        "Elastic event {} has no parameter start".format(k)
                    )
                rUniform = teSort[k].get("uniform")
                rMap = teSort[k].get("map")

                if rMap is not None:
                    if self.meshFile != rMap[0] + ".npz":
                        try:
                            with open(rMap[0] + ".npz") as rainfile:
                                rainfile.close()

                        except IOError:
                            print(
                                "Unable to open elastic file: {}.npz".format(rMap[0]),
                                flush=True,
                            )
                            raise IOError(
                                "The elastic file {} is not found for elastic event {}.".format(
                                    rMap[0] + ".npz", k
                                )
                            )
                        mdata = np.load(rMap[0] + ".npz")
                        teSet = mdata.files
                    else:
                        mdata = np.load(self.meshFile)
                        teSet = mdata.files
                    # TODO-REFACTOR: complex except, needs manual review (out of scope: accesses np.load() archive, not self.input; carries a user-facing diagnostic about missing NPZ field)
                    try:
                        rainKey = mdata[rMap[1]]
                        if rainKey is not None:
                            pass
                    except KeyError:
                        print(
                            "Field name {} is missing from elastic file {}.npz".format(
                                rMap[1], rMap[0]
                            ),
                            flush=True,
                        )
                        print(
                            "The following fields are available: {}".format(teSet),
                            flush=True,
                        )
                        print("Check your elastic file fields definition...", flush=True)
                        raise KeyError(
                            "Field name for elastic is not defined correctly or does not exist!"
                        )

                tedata = self._getTe(k, rStart, rMap, rUniform, tedata)

            if tedata["start"][0] > self.tStart:
                tmpT = []
                tmpT.insert(
                    0, {"start": self.tStart, "tUni": 0.0, "tMap": None, "tKey": None}
                )
                tedata = pd.concat(
                    [
                        pd.DataFrame(
                            tmpT, columns=["start", "tUni", "tMap", "tKey"]
                        ),
                        tedata,
                    ],
                    ignore_index=True,
                )
            self.tedata = tedata.copy()
            self.tedata.reset_index(drop=True, inplace=True)
            self.teNb = len(self.tedata)

        except KeyError:
            self.tedata = None

        return

    def _defineRain(self, k, rStart, rMap, rUniform, rZscale,raindata):
        """
        Define precipitation conditions.

        :arg k: precipitation event number
        :arg rStart: precipitation event start time
        :arg rMap: precipitation map file event
        :arg rUniform: precipitation uniform value event
        :arg rZscale: precipitation scaled with elevation value event
        :arg raindata: pandas dataframe storing each precipitation event
        :return: appended raindata
        """

        if rMap is None and rUniform is None and rZscale is None:
            print(
                "For each climate event a rainfall value (uniform), a rainfall "
                "grid (map), or an elevation-scaling pair (zscale: [A, B]) is required.",
                flush=True,
            )
            raise ValueError(
                "Climate event {} has no rainfall value (uniform), rainfall map "
                "(map), or elevation-scaling pair (zscale).".format(k)
            )

        tmpRain = []
        if rMap is None:
            if rUniform is not None:
                tmpRain.insert(
                    0,
                    {"start": rStart, "rUni": rUniform, "rzA": None, "rzB": None, "rMap": None, "rKey": None},
                )
            else:
                tmpRain.insert(
                    0,
                    {"start": rStart, "rUni": None, "rzA": rZscale[0], "rzB": rZscale[1],
                     "rMap": None, "rKey": None},
                )
        else:
            tmpRain.insert(
                0,
                {
                    "start": rStart,
                    "rUni": None,
                    "rzA": None,
                    "rzB": None,
                    "rMap": rMap[0] + ".npz",
                    "rKey": rMap[1],
                },
            )

        if k == 0:
            raindata = pd.DataFrame(tmpRain, columns=["start", "rUni", "rzA", "rzB", "rMap", "rKey"])
        else:
            raindata = pd.concat(
                [
                    raindata,
                    pd.DataFrame(tmpRain, columns=["start", "rUni", "rzA", "rzB", "rMap", "rKey"]),
                ],
                ignore_index=True,
            )

        return raindata

    def _readRain(self):
        """
        Parse rain forcing conditions.
        """

        raindata = None
        # TODO-REFACTOR: complex except, needs manual review (outer-section: sets self.raindata = None on missing "climate")
        try:
            rainDict = self.input["climate"]
            rainSort = sorted(rainDict, key=itemgetter("start"))
            for k in range(len(rainSort)):
                rStart = None
                rUniform = None
                rZscale = None
                rMap = None
                # TODO-REFACTOR: complex except, needs manual review (rStart is required)
                try:
                    rStart = rainSort[k]["start"]
                except KeyError:
                    print(
                        "For each climate event a start time is required.", flush=True
                    )
                    raise ValueError(
                        "Climate event {} has no parameter start".format(k)
                    )
                rUniform = rainSort[k].get("uniform")
                rZscale = rainSort[k].get("zscale")
                rMap = rainSort[k].get("map")

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
                    # TODO-REFACTOR: complex except, needs manual review (out of scope: accesses np.load() archive, not self.input; carries a user-facing diagnostic about missing NPZ field)
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

                raindata = self._defineRain(k, rStart, rMap, rUniform, rZscale, raindata)

            if raindata["start"][0] > self.tStart:
                tmpRain = []
                tmpRain.insert(
                    0, {"start": self.tStart, "rUni": 0.0, "rzA":None, "rzB":None, "rMap": None, "rKey": None}
                )
                raindata = pd.concat(
                    [
                        pd.DataFrame(
                            tmpRain, columns=["start", "rUni", "rzA", "rzB", "rMap", "rKey"]
                        ),
                        raindata,
                    ],
                    ignore_index=True,
                )
            self.raindata = raindata.copy()  # [raindata["start"] >= self.tStart]
            self.raindata.reset_index(drop=True, inplace=True)
            self.rainNb = len(self.raindata)

        except KeyError:
            self.raindata = None

        return

    def _readCompaction(self):
        """
        Read compaction parameters.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section multi-default on missing "compaction")
        try:
            compDict = self.input["compaction"]
            self.phi0s = compDict.get("phis", 0.49)
            self.z0s = compDict.get("z0s", 3700.0)

        except KeyError:
            self.phi0s = 0.49
            self.z0s = 3700.0

        return

    def _readFlex(self):
        """
        Parse flexural isostasy variables.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section: sets flexOn=False, flex_method='FD' on missing "flexure")
        try:
            flexDict = self.input["flexure"]
            self.flexOn = True
            self.flex_method = flexDict.get("method", 'FD')

            if self.flex_method != 'global' and self.flex_method != 'FD' and self.flex_method != 'FFT':
                print(
                    "Method {} is not in the list of possible methods".format(
                        self.flex_method), flush=True,
                )
                raise ValueError(
                    "Method name for flexure is not recognised choices are FD, FFT or global."
                )
            self.reg_dx = flexDict.get("regdx", 1000.0)
            # raise ValueError("Flexure definition: regular grid spacing is required.")
            self.rgrd_interp = flexDict.get("ninterp", 4)
            self.flex_rhoa = flexDict.get("rhoa", 3300.0)
            self.flex_eet = flexDict.get("thick", 10000.0)
            self.flex_rhos = flexDict.get("rhoc", 2300.0)
            self.young = flexDict.get("young", 65e9)

        except KeyError:
            self.flexOn = False
            self.flex_method = 'FD'

        self._extraFlex()

        return

    def _extraFlex(self):
        """
        Read flexure additional information.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section: sets flexOn=False on missing "flexure")
        try:
            flexDict = self.input["flexure"]
            self.nu = flexDict.get("nu", 0.25)
            # Resolution at which the SH expansion is performed. Input is
            # regridded to this; output is interpolated back to the input grid.
            # Default 0.25 deg -> lmax = 359, well below the elastic cutoff.
            self.flex_res_deg = flexDict.get("res_deg", 0.25)
            self.flex_bcN = flexDict.get("bcN", "0Slope0Shear")
            self.flex_bcS = flexDict.get("bcS", "0Slope0Shear")
            self.flex_bcE = flexDict.get("bcE", "0Slope0Shear")
            self.flex_bcW = flexDict.get("bcW", "0Slope0Shear")
        except KeyError:
            self.flexOn = False

        return

    def _readOrography(self):
        """
        Parse orographic precipitation variables.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section: sets oroOn=False on missing "orography")
        try:
            oroDict = self.input["orography"]
            self.oroOn = True

            # TODO-REFACTOR: complex except, needs manual review (regdx is required when "orography" is present)
            try:
                self.reg_dx = oroDict["regdx"]
            except KeyError:
                raise ValueError("Orographic definition: regular grid spacing is required.")
            # TODO-REFACTOR: complex except, needs manual review (try body has bounds-check + raise ValueError)
            try:
                self.wind_latitude = oroDict["latitude"]
                if self.wind_latitude > 90 or self.wind_latitude < -90:
                    print(
                        "Latitude for orographic rain needs to be between -90 and 90.",
                        flush=True,
                    )
                    raise ValueError("Latitude value not appropriately set.")
            except KeyError:
                self.wind_latitude = 0.0
            self.wind_speed = oroDict.get("wind_speed", 10.)
            self.wind_dir = oroDict.get("wind_dir", 0.0)
            self.oro_nm = oroDict.get("nm", 0.01)
            self.oro_hw = oroDict.get("hw", 3400.0)

            self.rgrd_interp = 4
            self._extraOrography(oroDict)

        except KeyError:
            self.oroOn = False

        return

    def _extraOrography(self, oroDict):
        """
        Read domain additional information.
        """

        lapse_rate = oroDict.get("env_lapse_rate", -4.0)
        lapse_rate_m = oroDict.get("moist_lapse_rate", -7.0)
        ref_density = oroDict.get("ref_density", 7.4e-3)
        if lapse_rate == 0:
            raise ValueError(
                "Orographic precipitation: env_lapse_rate must be non-zero."
            )
        self.oro_cw = ref_density * lapse_rate_m / lapse_rate
        self.oro_conv_time = oroDict.get("conv_time", 1000.0)
        self.oro_fall_time = oroDict.get("fall_time", 1000.0)
        self.oro_precip_base = oroDict.get("precip_base", 7.0)
        self.oro_precip_min = oroDict.get("precip_min", 0.01)
        self.rainfall_frequency = oroDict.get("rainfall_frequency", 1)

        return

    def _readIce(self):
        """
        Parse ice flow variables.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section: sets iceOn=False + 4 local defaults on missing "ice")
        try:
            iceDict = self.input["ice"]
            self.iceOn = True
            self.gaussIce = iceDict.get("diff", 10.0)
            elaH = iceDict.get("hela", 2000.0)
            iceH = iceDict.get("hice", 2400.0)
            iceT = iceDict.get("hterm", 1800.0)  # glacier terminus
            icefile = iceDict.get("evol")
            self.Kice = iceDict.get("Ki", 0.0)
            self.iceDir = iceDict.get("icedir", 1)
            self.meltfac = iceDict.get("melt", 10.)  # Melting factor adjustment
            self.icewf = iceDict.get("fwidth", 1.5)  # width factor (a in Eq. (8))
            # thickness-to-width ratio (delta in Eq. (9))
            self.icewe = iceDict.get("eheight", 0.25)
        except KeyError:
            self.iceOn = False
            icefile = None
            elaH = None
            iceH = None
            iceT = None

        self._extraIce(icefile, elaH, iceH, iceT)

        return

    def _extraIce(self, icefile, elaH, iceH, iceT):

        if icefile is not None:
            try:
                with open(icefile) as fice:
                    fice.close()
                    try:
                        icedata = pd.read_csv(
                            icefile,
                            sep=r",",
                            engine="c",
                            header=None,
                            na_filter=False,
                            dtype=np.float64,
                            low_memory=False,
                        )

                    except ValueError:
                        try:
                            icedata = pd.read_csv(
                                icefile,
                                sep=r"\s+",
                                engine="c",
                                header=None,
                                na_filter=False,
                                dtype=np.float64,
                                low_memory=False,
                            )

                        except ValueError:
                            print(
                                "The ice evolution file is not well formed: it should be comma or tab separated",
                                flush=True,
                            )
                            raise ValueError("Wrong formating of ice evolution file.")
            except IOError:
                print("Unable to open file: ", icefile)
                raise IOError("The ice evolution file is not found...")

            if icedata[0].min() > self.tStart:
                tmpS = []
                tmpS.insert(0, {0: self.tStart, 1: icedata[1].iloc[0], 2: icedata[2].iloc[0], 3: icedata[3].iloc[0]})
                icedata = pd.concat([pd.DataFrame(tmpS), icedata], ignore_index=True)
            if icedata[0].max() < self.tEnd:
                tmpE = []
                tmpE.insert(0, {0: self.tEnd, 1: icedata[1].iloc[-1], 2: icedata[2].iloc[-1], 3: icedata[3].iloc[-1]})
                icedata = pd.concat([icedata, pd.DataFrame(tmpE)], ignore_index=True)
            self.iceT = interp1d(icedata[0], icedata[1], kind="linear")
            self.elaH = interp1d(icedata[0], icedata[2], kind="linear")
            self.iceH = interp1d(icedata[0], icedata[3], kind="linear")
        elif self.iceOn:
            year = np.linspace(self.tStart, self.tEnd + self.dt, num=11, endpoint=True)
            iceval = np.full(len(year), iceT)
            self.iceT = interp1d(year, iceval, kind="linear")
            iceval = np.full(len(year), elaH)
            self.elaH = interp1d(year, iceval, kind="linear")
            iceval = np.full(len(year), iceH)
            self.iceH = interp1d(year, iceval, kind="linear")

        return

    def _readOut(self):
        """
        Parse output directory.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section multi-default on missing "output")
        try:
            outDict = self.input["output"]
            self.outputDir = outDict.get("dir", "output")
            self.makedir = outDict.get("makedir", False)
        except KeyError:
            self.outputDir = "output"
            self.makedir = False

        if self.rStep > 0:
            self.makedir = False

        return
