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

# Sentinel for an unprescribed glacier terminus altitude (``ice.hterm``): the
# effective terminus floor is then the sea-level position (see
# iceplex._iceFlowMFD). Any value safely below the lowest plausible sea
# level works.
TERMINUS_UNSET = -1.0e10


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

        # Boundary conditions, one char per edge in order [North, East, South,
        # West] — N at y=ymax, E at x=xmax, S at y=ymin, W at x=xmin:
        #   'o' = open    — deep base-level outlet (edge forced below sea level);
        #   'f' = fixed   — fixed base-level outlet at the natural edge elevation
        #                   (water/sediment leave the domain — NOT a no-flux wall);
        #   'w' = wall    — true no-flux wall: flow and sediment are contained
        #                   (the edge nodes behave as interior, sediment deposits
        #                   against them, mass is conserved);
        #   'c' = cyclic  — periodic (requires a periodic/cylinder input mesh).
        # Legacy digits map '0' -> 'o' and '1' -> 'f' (preserving old behaviour).
        # Cyclic edges must come as an OPPOSITE PAIR (North+South or East+West)
        # and at most ONE pair may be cyclic — i.e. up to two periodic edges,
        # never all four (a full torus has no boundary and would be mis-detected
        # as a global model).
        bc = str(domainDict.get("bc", 'oooo'))
        bc = bc.replace('0', 'o').replace('1', 'f')   # legacy digits
        self.boundCond = bc
        if len(self.boundCond) != 4 or any(c not in 'ofwc' for c in self.boundCond):
            raise ValueError(
                "domain.bc must be 4 characters from {'o' open, 'f' fixed, "
                "'w' wall, 'c' cyclic} (legacy '0'/'1' accepted) in [North, "
                "East, South, West] order; got '%s'." % str(domainDict.get("bc"))
            )
        bcN, bcE, bcS, bcW = self.boundCond
        if (bcN == 'c') != (bcS == 'c'):
            raise ValueError(
                "Cyclic (periodic) boundaries must be paired: set BOTH North and "
                "South to 'c' (got bc='%s')." % self.boundCond
            )
        if (bcE == 'c') != (bcW == 'c'):
            raise ValueError(
                "Cyclic (periodic) boundaries must be paired: set BOTH East and "
                "West to 'c' (got bc='%s')." % self.boundCond
            )
        if bcN == 'c' and bcE == 'c':
            raise ValueError(
                "At most one pair of edges may be cyclic (North/South OR "
                "East/West, not both); got bc='%s'." % self.boundCond
            )

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
            # Non-linear SPL (spl_n != 1) SNES controls. The transport-limited
            # (G>0) solver is otherwise a bare ngmres accelerator that stalls on
            # the stiff residual, so the primary defaults to 'qn' (L-BFGS) with
            # the same complementary-fallback robustness net as the soil solver.
            self.snes_maxit = int(splDict.get("maxIter", 500))
            self.snes_rtol = float(splDict.get("rtol", 1.0e-6))
            self.snes_atol = float(splDict.get("atol", 1.0e-6))
            self.nlspl_solver = str(splDict.get("solver", "qn"))
            self.nlspl_pc = str(splDict.get("pcType", "hypre"))
        except KeyError:
            self.K = 1.0e-12
            self.coeffd = 0.0
            self.fDepa = 0.0
            self.spl_m = 0.5
            self.spl_n = 1.0
            self.snes_maxit = 500
            self.snes_rtol = 1.0e-6
            self.snes_atol = 1.0e-6
            self.nlspl_solver = "qn"
            self.nlspl_pc = "hypre"

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
            # Marine/lake non-linear-diffusion solver. 'ts' (default) = adaptive
            # non-linear PETSc TS (rosw). 'picard' = lagged-diffusivity backward-
            # Euler with linear solves (faster, robust at the C_d kink), an
            # opt-in approximation; picardSub sub-steps over dt, picardIts Picard
            # iterations per sub-step. See sed/hillslope._diffuseImplicitPicard.
            self.marineSolver = str(hillDict.get("marineSolver", "ts")).lower()
            self.picardSub = int(hillDict.get("picardSub", 10))
            self.picardIts = int(hillDict.get("picardIts", 2))
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
            self.marineSolver = "ts"
            self.picardSub = 10
            self.picardIts = 2

        self.clinSlp = max(1.0e-6, self.clinSlp)
        self.nl_pit_inlet_bias = min(1.0, max(0.0, self.nl_pit_inlet_bias))
        if self.marineSolver not in ("ts", "picard"):
            self.marineSolver = "ts"
        self.picardSub = max(1, self.picardSub)
        self.picardIts = max(1, self.picardIts)

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

            # Soil-SPL nonlinear-solver (SNES) controls. The soil-aware SPL
            # residual is stiffer than the bedrock case (soil-production
            # coupling), so the iteration budget defaults to the same value as
            # nlSPL (500, up from the previous 100) and the tolerances and
            # preconditioner are exposed for tuning at scale.
            self.soil_maxit = int(soilDict.get("maxIter", 500))
            self.soil_rtol = float(soilDict.get("rtol", 1.0e-6))
            self.soil_atol = float(soilDict.get("atol", 1.0e-6))
            # Preconditioner for the soil SNES Krylov solve: 'hypre'
            # (BoomerAMG, default), 'gamg', 'bjacobi', 'asm', ...
            self.soil_pc = str(soilDict.get("pcType", "hypre"))
            # Primary nonlinear solver: 'qn' (limited-memory quasi-Newton /
            # L-BFGS, default — ~2.4x faster than ngmres at the same tolerance
            # and solution on a global soil model) or 'ngmres' (accelerator +
            # multigrid PC). Whichever is chosen, the other is the fallback.
            self.soil_solver = str(soilDict.get("solver", "qn"))

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
            self.soil_maxit = 500
            self.soil_rtol = 1.0e-6
            self.soil_atol = 1.0e-6
            self.soil_pc = "hypre"
            self.soil_solver = "qn"

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

    def _defineEvap(self, k, eStart, eMap, eUniform, evapdata):
        """
        Define evaporation conditions for one climate event.

        Mirrors `_defineRain` but for evaporation. Evaporation is opt-in
        per-row: if both `eMap` and `eUniform` are None for this event,
        `evapdata` is returned unchanged (no row appended).

        :arg k: climate event number
        :arg eStart: event start time
        :arg eMap: evaporation map (tuple of path+key) or None
        :arg eUniform: evaporation uniform scalar (m/yr) or None
        :arg evapdata: pandas dataframe storing each evaporation event (or None on first call)
        :return: appended evapdata (or unchanged if this row has no evap)
        """

        if eMap is None and eUniform is None:
            return evapdata

        tmpEvap = []
        if eMap is None:
            tmpEvap.insert(
                0,
                {"start": eStart, "eUni": eUniform, "eMap": None, "eKey": None},
            )
        else:
            tmpEvap.insert(
                0,
                {
                    "start": eStart,
                    "eUni": None,
                    "eMap": eMap[0] + ".npz",
                    "eKey": eMap[1],
                },
            )

        if evapdata is None:
            evapdata = pd.DataFrame(tmpEvap, columns=["start", "eUni", "eMap", "eKey"])
        else:
            evapdata = pd.concat(
                [
                    evapdata,
                    pd.DataFrame(tmpEvap, columns=["start", "eUni", "eMap", "eKey"]),
                ],
                ignore_index=True,
            )

        return evapdata

    def _readRain(self):
        """
        Parse rain and evaporation forcing conditions.

        Both share the same `[climate]` YAML block. Each climate event row
        may declare rainfall (`uniform`/`map`/`zscale`) and, optionally,
        evaporation (`evap_uniform`/`evap_map`). Evaporation is opt-in: if
        no row has an evap field, `self.evapdata` ends up as None and the
        downstream solver bypasses both evap hooks.
        """

        raindata = None
        evapdata = None
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
                eUniform = rainSort[k].get("evap_uniform")
                eMap = rainSort[k].get("evap_map")

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

                if eMap is not None:
                    if self.meshFile != eMap[0] + ".npz":
                        try:
                            with open(eMap[0] + ".npz") as evapfile:
                                evapfile.close()
                        except IOError:
                            print(
                                "Unable to open evaporation file: {}.npz".format(eMap[0]),
                                flush=True,
                            )
                            raise IOError(
                                "The evaporation file {} is not found for climatic event {}.".format(
                                    eMap[0] + ".npz", k
                                )
                            )
                        edata = np.load(eMap[0] + ".npz")
                        evapSet = edata.files
                    else:
                        edata = np.load(self.meshFile)
                        evapSet = edata.files
                    if eMap[1] not in evapSet:
                        print(
                            "Field name {} is missing from evaporation file {}.npz".format(
                                eMap[1], eMap[0]
                            ),
                            flush=True,
                        )
                        print(
                            "The following fields are available: {}".format(evapSet),
                            flush=True,
                        )
                        print("Check your evaporation file fields definition...", flush=True)
                        raise KeyError(
                            "Field name for evaporation is not defined correctly or does not exist!"
                        )

                raindata = self._defineRain(k, rStart, rMap, rUniform, rZscale, raindata)
                evapdata = self._defineEvap(k, rStart, eMap, eUniform, evapdata)

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

            if evapdata is not None:
                # Prepend a zero-evap row at tStart if the first declared event
                # is later, so _updateEvap has a defined value from the start.
                if evapdata["start"][0] > self.tStart:
                    tmpEvap = [
                        {"start": self.tStart, "eUni": 0.0, "eMap": None, "eKey": None}
                    ]
                    evapdata = pd.concat(
                        [
                            pd.DataFrame(
                                tmpEvap, columns=["start", "eUni", "eMap", "eKey"]
                            ),
                            evapdata,
                        ],
                        ignore_index=True,
                    )
                self.evapdata = evapdata.copy()
                self.evapdata.reset_index(drop=True, inplace=True)
            else:
                self.evapdata = None

        except KeyError:
            self.raindata = None
            self.evapdata = None

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

        self._extraStrata()

        return

    def _extraStrata(self):
        """
        Parse the optional dual-lithology (coarse/fine) stratigraphy block.

        Continuation of ``_readCompaction`` (the coarse porosity curve
        defaults to the single-fraction ``compaction`` values, so a model
        with ``strata: dual: False`` — or no ``strata`` block at all — runs
        the existing single-fraction path bitwise-unchanged).

        Sets ``self.stratLith`` (master opt-in flag) plus the per-lithology
        parameters consumed by later phases (porosity-depth curves, bedrock
        composition, marine fine transport efficiency, per-fraction lake
        inlet bias, and subaerial/subaqueous diffusivities). All of these
        are inert while ``self.stratLith`` is False.

        Dual lithology is only meaningful when stratigraphy recording is on;
        if requested with stratigraphy disabled (``self.stratNb == 0``) the
        flag is forced back to False with a rank-0 warning.

        See ``docs/DESIGN_DUAL_LITHOLOGY.md`` and AGENTS.md > The ``_extra*``
        methods are mandatory continuations.
        """

        # Coarse lithology porosity-depth defaults to the single-fraction
        # compaction curve so the dual-off path is unchanged. Fine lithology
        # defaults to a higher surface porosity / shallower decay; both fine
        # values are inert unless dual lithology is enabled.
        self.phi0c = self.phi0s
        self.z0c = self.z0s
        self.phi0f = 0.63
        self.z0f = 1960.0
        # Fine erodibility relative to coarse: K_fine = K_coarse * fine_k_factor.
        # 1.0 (default) = no lithology contrast in the SPL erodibility blend.
        self.fine_k_factor = 1.0
        self.bedrock_coarse_frac = 0.5
        self.fine_efficiency = 0.5
        self.pit_inlet_bias_coarse = 0.50
        self.pit_inlet_bias_fine = 0.0
        # Fine diffusivity relative to coarse: a multiplier on the existing
        # hillslope coefficients (Cda/Cdm/nlK), NOT an absolute coefficient.
        # Cd_eff = Cd_base * (fc + ff * fine_diff_factor). 1.0 = no contrast.
        self.fine_diff_factor = 1.0

        # TODO-REFACTOR: complex except, needs manual review (outer-section: sets stratLith=False on missing "strata")
        try:
            strataDict = self.input["strata"]
            self.stratLith = bool(strataDict.get("dual", False))

            coarseDict = strataDict.get("coarse", {})
            self.phi0c = coarseDict.get("phi0", self.phi0s)
            self.z0c = coarseDict.get("z0", self.z0s)

            fineDict = strataDict.get("fine", {})
            self.phi0f = fineDict.get("phi0", self.phi0f)
            self.z0f = fineDict.get("z0", self.z0f)
            self.fine_k_factor = fineDict.get("k_factor", self.fine_k_factor)

            self.bedrock_coarse_frac = strataDict.get(
                "bedrock_coarse_frac", self.bedrock_coarse_frac
            )
            self.fine_efficiency = strataDict.get("fine_efficiency", self.fine_efficiency)

            biasDict = strataDict.get("pitInletBias", {})
            if isinstance(biasDict, dict):
                self.pit_inlet_bias_coarse = biasDict.get(
                    "coarse", self.pit_inlet_bias_coarse
                )
                self.pit_inlet_bias_fine = biasDict.get("fine", self.pit_inlet_bias_fine)

            self.fine_diff_factor = strataDict.get(
                "fine_diff_factor", self.fine_diff_factor
            )

        except KeyError:
            self.stratLith = False

        # Dual lithology only makes sense when stratigraphy recording is on.
        if self.stratLith and self.stratNb == 0:
            if MPIrank == 0:
                print(
                    "Dual-lithology (strata: dual) requires stratigraphy to be "
                    "enabled (set a positive `strat` interval in the time block); "
                    "falling back to single-fraction sediment.",
                    flush=True,
                )
            self.stratLith = False

        # Clamp inlet-bias fractions to [0, 1] (mirrors nl_pit_inlet_bias).
        self.pit_inlet_bias_coarse = min(1.0, max(0.0, self.pit_inlet_bias_coarse))
        self.pit_inlet_bias_fine = min(1.0, max(0.0, self.pit_inlet_bias_fine))

        self._extraProvenance()

        return

    def _extraProvenance(self):
        """
        Parse the optional in-model **sediment-provenance tracer** block. N
        source-rock classes are carried (as a passive label — no physics
        feedback) through erosion, transport, deposition and the stratigraphic
        record, giving conservation-exact, recycling-aware provenance per layer.
        See ``docs/DESIGN_PROVENANCE.md`` §6.

        Sets ``self.provOn`` (master opt-in), ``self.provNb`` (class count), the
        per-vertex source-class source (``uniform`` scalar or ``source``
        ``[file, key]`` map, resolved post-mesh in ``readStratLayers``), and an
        optional ``cu_weight``. All inert while ``provOn`` is False; like dual
        lithology it requires stratigraphy (``stratNb > 0``).
        """
        self.provOn = False
        self.provNb = 0
        self._provSourceUniform = None
        self._provSourceMap = None
        self.prov_cu_weight = None
        try:
            provDict = self.input["provenance"]
        except KeyError:
            return

        self.provNb = int(provDict.get("classes", 0))
        self._provSourceUniform = provDict.get("uniform")
        self._provSourceMap = provDict.get("source")          # [file, key] or None
        cw = provDict.get("cu_weight")
        if cw is not None:
            self.prov_cu_weight = np.asarray(cw, dtype=np.float64)

        self.provOn = self.provNb > 0
        if self.provOn and self.stratNb == 0:
            if MPIrank == 0:
                print(
                    "Provenance tracers (provenance:) require stratigraphy to be "
                    "enabled (set a positive `strat` interval in the time block); "
                    "disabling provenance.",
                    flush=True,
                )
            self.provOn = False
        # Default to a single source class everywhere if none specified.
        if self.provOn and self._provSourceUniform is None and self._provSourceMap is None:
            self._provSourceUniform = 0
        return

    def _readFlex(self):
        """
        Parse flexural isostasy variables.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section: sets flexOn=False, flex_method='fem' on missing "flexure")
        try:
            flexDict = self.input["flexure"]
            self.flexOn = True
            # 'fem'  -> parallel FV biharmonic on the flat DMPlex (default)
            # 'global' -> spherical-harmonic solve (global models)
            self.flex_method = flexDict.get("method", 'fem')

            if self.flex_method not in ('global', 'fem'):
                print(
                    "Method {} is not in the list of possible methods".format(
                        self.flex_method), flush=True,
                )
                raise ValueError(
                    "Method name for flexure is not recognised: choices are 'fem' (flat) or 'global'."
                )
            # Number of source points for the global ('global') DH-grid
            # interpolation (KDTree in _buildDHGrid); unused by 'fem'.
            self.rgrd_interp = flexDict.get("ninterp", 4)
            self.flex_rhoa = flexDict.get("rhoa", 3300.0)
            self.flex_eet = flexDict.get("thick", 10000.0)
            self.flex_rhos = flexDict.get("rhoc", 2300.0)
            self.young = flexDict.get("young", 65e9)

        except KeyError:
            self.flexOn = False
            self.flex_method = 'fem'

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
            # Global varying-Te flexure: iterative (Anderson/Picard) solve
            # controls + how often flexure is applied. Defaults reproduce the
            # previous hard-coded behaviour (every step, max_iter=50, tol=5e-4,
            # relax=1.0). See tools/addprocess._cmptFlexGlobal and the flexure
            # phase in model.runProcesses.
            self.flex_max_iter = int(flexDict.get("maxIter", 50))
            self.flex_tol = float(flexDict.get("tol", 5.0e-4))
            self.flex_relax = float(flexDict.get("relax", 1.0))
            # Apply flexure every `interval` goSPL steps (the load accumulates in
            # between). 1 = every step (default, unchanged behaviour).
            self.flex_interval = int(flexDict.get("interval", 1))
        except KeyError:
            self.flexOn = False

        if getattr(self, "flexOn", False):
            self.flex_max_iter = max(1, self.flex_max_iter)
            self.flex_tol = max(1.0e-12, self.flex_tol)
            self.flex_relax = min(1.0, max(1.0e-3, self.flex_relax))
            self.flex_interval = max(1, self.flex_interval)

        return

    def _readOrography(self):
        """
        Parse orographic precipitation variables.
        """

        try:
            oroDict = self.input["orography"]
            self.oroOn = True

            # Orographic precipitation is solved directly on the unstructured
            # mesh (advection-relaxation, no regular grid / FFT). Only the
            # uniform wind and the moisture parameters are needed.
            self.wind_speed = oroDict.get("wind_speed", 10.)
            self.wind_dir = oroDict.get("wind_dir", 0.0)

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

        # When an `ice` section is present, goSPL runs the diagnostic glacial
        # model: the ELA accumulation is routed downhill into an ice discharge
        # from which a thickness and a basal sliding velocity are derived (no
        # ice-dynamics solve), driving glacial abrasion, till transport and ice
        # loading. See docs/DESIGN_ICE_SHEET.md.
        #
        # The glacier-geometry altitudes (terminus / ELA / ice-cap) can each be
        # a uniform scalar OR a per-vertex map ``[file, key]`` — the latter for
        # global models where the ELA varies strongly with latitude (tropical
        # vs polar). They may also be a TIME SERIES via the optional `glaciers`
        # list (mirroring the precipitation `climate` block): each entry has a
        # `start` time and uniform-or-map hela/hice/hterm, stepped over time by
        # _updateIce. `self._iceTimeSeries` carries that series (None = legacy
        # uniform/`evol` scalars handled by _extraIce).
        self._iceTimeSeries = None
        self._iceSeriesIdx = -1
        self.termMesh = None
        self.elaMesh = None
        self.iceMesh = None
        useSeries = False
        try:
            iceDict = self.input["ice"]
            self.iceOn = True
            icefile = iceDict.get("evol")
            glaciers = iceDict.get("glaciers")
            elaH = iceDict.get("hela", 2000.0)
            iceH = iceDict.get("hice", 2400.0)
            # Terminus floor = max(hterm, sea level) (applied in iceplex): ice is
            # never kept below the (possibly time-varying) sea surface. The
            # sentinel means "not prescribed" -> the sea-level position.
            iceT = iceDict.get("hterm", TERMINUS_UNSET)

            # Diagnostic glacial driver parameters. The ELA accumulation is routed
            # downhill (`icedir` flow directions) into an ice discharge, from which
            # a Bahr thickness (`eheight`*`fwidth`*Q^0.3) and a basal sliding
            # velocity (Glen sliding law, `slide`/`glen`) are derived; these drive
            # the glacial abrasion / till / loading machinery with no ice-dynamics
            # solve (robust and physical at any resolution).
            self.iceDir = int(iceDict.get("icedir", 1))      # MFD flow directions
            self.ice_meltfac = iceDict.get("melt", 10.0)     # ablation amplifier
            self.icewf = iceDict.get("fwidth", 1.5)          # Bahr width factor
            self.icewe = iceDict.get("eheight", 0.25)        # Bahr thickness factor
            self.ice_slide = iceDict.get("slide", 1.0e-3)    # basal sliding coeff
            self.ice_glen = iceDict.get("glen", 3.0)         # Glen sliding exponent
            # Glacial-meltwater model for the river coupling. True (default):
            # discharge-conserving — the accumulation (water that fell as ice) is
            # routed down-glacier and released as meltwater where the ice melts
            # out, so Σ meltwater == Σ accumulation (the steady-state, long-
            # timescale assumption; closes the glacial water budget). False: the
            # local precipitation-scaled ablation rate (loses water downstream).
            self.ice_melt_conserve = bool(iceDict.get("melt_conserve", True))
            # Surface-mass-balance accumulation controls (applied to the positive
            # ELA ramp only; ablation is unchanged). `accum_factor` is a
            # precipitation->ice conversion fraction (full precipitation is rarely
            # all snow/ice); `accum_max` caps the accumulation rate (m ice/yr) at a
            # realistic ceiling. Defaults (1.0, no cap) keep the raw precip rate.
            self.ice_accum_factor = float(iceDict.get("accum_factor", 1.0))
            am = iceDict.get("accum_max", None)
            self.ice_accum_max = None if am is None else float(am)

            abrDict = iceDict.get("abrasion", {})
            self.ice_Kg = abrDict.get("Kg", 0.0)             # abrasion coeff (0 = off)
            self.ice_abr_l = abrDict.get("l", 1.0)           # sliding-velocity exponent
            # Lateral (valley-wall) glacial erosion — widens glaciated valleys
            # toward a U-profile. Off by default (Kl = 0). `lat_l` defaults to the
            # vertical-abrasion exponent `l`.
            self.ice_Kl = abrDict.get("Kl", 0.0)
            self.ice_lat_l = abrDict.get("lat_l", self.ice_abr_l)

            tillDict = iceDict.get("till", {})
            # Till on by default: when glacial abrasion is enabled (Kg > 0) the
            # eroded rock is carried as till and deposited as moraine — the
            # complete glacial sediment cycle. Set False to send abrasion
            # straight to the fluvial system instead. (No cost when Kg == 0.)
            self.ice_till_on = bool(tillDict.get("on", True))
            # Catchment-aware till routing on by default: the till is routed down
            # the ice-surface flow network and melts out toward each terminus, so
            # deposition stays connected to the upstream erosion. Set False for
            # the older melt-weighted spreading, which pools the GLOBAL abraded
            # volume across all melt cells (it decouples erosion and deposition
            # across separate ice masses — misleading on multi-glacier domains).
            self.ice_till_route = bool(tillDict.get("route", True))

            # Use the per-vertex / time-series path when a `glaciers` series is
            # given or any top-level altitude is a map.
            topIsMap = any(
                isinstance(v, (list, tuple)) for v in (elaH, iceH, iceT)
            )
            useSeries = glaciers is not None or topIsMap
            if useSeries and icefile is not None:
                if MPIrank == 0:
                    print(
                        "Ice: an `evol` time series and per-vertex / `glaciers` "
                        "geometry were both given; using `evol` (the maps are "
                        "ignored).",
                        flush=True,
                    )
                useSeries = False
            if useSeries:
                self._iceTimeSeries = self._buildIceSeries(
                    glaciers, elaH, iceH, iceT
                )
        except KeyError:
            self.iceOn = False
            icefile = None
            elaH = None
            iceH = None
            iceT = None
            self.iceDir = 1
            self.ice_meltfac = 10.0
            self.icewf = 1.5
            self.icewe = 0.25
            self.ice_slide = 1.0e-3
            self.ice_glen = 3.0
            self.ice_melt_conserve = True
            self.ice_accum_factor = 1.0
            self.ice_accum_max = None
            self.ice_Kg = 0.0
            self.ice_abr_l = 1.0
            self.ice_Kl = 0.0
            self.ice_lat_l = 1.0
            self.ice_till_on = True
            self.ice_till_route = True

        # Legacy uniform / `evol` path builds the scalar time functions. In
        # series mode the geometry comes from _iceTimeSeries via _updateIce, so
        # _extraIce is skipped.
        if not useSeries:
            self._extraIce(icefile, elaH, iceH, iceT)

        return

    def _iceGeomField(self, val):
        """
        Split a glacier-geometry input (``hela``/``hice``/``hterm``) into a
        scalar fallback and an optional ``[file, key]`` map spec. A list/tuple
        value is a per-vertex npz map; anything else is a uniform scalar.

        :return: (scalar_or_None, map_spec_or_None)
        """
        if isinstance(val, (list, tuple)):
            return None, list(val)
        return val, None

    def _buildIceSeries(self, glaciers, elaTop, iceTop, termTop):
        """
        Build the glacier-geometry time series consumed by ``_updateIce``.

        ``glaciers`` (when given) is a list of ``{start, hela, hice, hterm}``
        events (each altitude uniform-or-map); otherwise the top-level
        ``hela``/``hice``/``hterm`` define a single interval starting at
        ``tStart``. Each interval stores, per field, a ``(scalar, map_spec)``
        pair (exactly one is non-None). Map files/keys are validated here; the
        arrays are loaded lazily on interval change in ``_updateIce``.

        :return: list of intervals sorted by start time.
        """
        if glaciers is not None:
            events = sorted(glaciers, key=itemgetter("start"))
        else:
            events = [
                {"start": self.tStart, "hela": elaTop, "hice": iceTop, "hterm": termTop}
            ]
        series = []
        for ev in events:
            interval = {"start": ev["start"]}
            for fld, default in (("hela", 2000.0), ("hice", 2400.0), ("hterm", TERMINUS_UNSET)):
                sc, spec = self._iceGeomField(ev.get(fld, default))
                if spec is not None:
                    self._checkMap(spec, "ice %s" % fld)
                interval[fld] = (sc, spec)
            series.append(interval)
        return series

    def _extraIce(self, icefile, elaH, iceH, iceT):
        """
        Legacy uniform glacier geometry: build the scalar time functions
        ``self.iceT`` / ``self.elaH`` / ``self.iceH`` from an ``evol`` CSV or
        from constant ``hterm`` / ``hela`` / ``hice``. The per-vertex /
        time-series map path is handled separately by ``_buildIceSeries`` +
        ``_updateIce``.
        """

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

    def _checkMap(self, spec, name):
        """
        Validate that a ``[file, key]`` npz map exists and contains the field,
        without loading the (potentially large) array. Used at parse time for
        the glacier-geometry maps; the array is loaded lazily by ``_loadIceMap``
        on interval change.
        """
        fname, key = spec[0], spec[1]
        try:
            with open(fname + ".npz") as f:
                f.close()
        except IOError:
            print("Unable to open %s map file: %s.npz" % (name, fname), flush=True)
            raise IOError("The %s map file is not found." % name)
        mdata = np.load(fname + ".npz")
        if key not in mdata:
            print(
                "Field '%s' is missing from %s map file %s.npz" % (key, name, fname),
                flush=True,
            )
            raise ValueError("Missing field in %s map file." % name)
        del mdata
        return

    def _loadIceMap(self, spec, name):
        """
        Load a per-vertex glacier-geometry map (full-mesh array) from a
        validated ``[file, key]`` npz spec — the same convention as the
        precipitation/tectonic maps.
        """
        fname, key = spec[0], spec[1]
        mdata = np.load(fname + ".npz")
        arr = mdata[key].astype(np.float64)
        del mdata
        return arr

    def _readOut(self):
        """
        Parse output directory.
        """

        # TODO-REFACTOR: complex except, needs manual review (outer-section multi-default on missing "output")
        try:
            outDict = self.input["output"]
            self.outputDir = outDict.get("dir", "output")
            self.makedir = outDict.get("makedir", False)
            # Opt-in wall-clock phase profiler (see gospl/tools/profiler.py).
            self.profileFlag = outDict.get("profile", False)
            # HDF5 output compression. `gzip` (default) keeps the historical
            # behaviour; an int sets the gzip level (0-9); `none`/False writes
            # uncompressed (faster I/O, larger files) — see outmesh._h5opts.
            self.outCompress = outDict.get("compression", "gzip")
        except KeyError:
            self.outputDir = "output"
            self.makedir = False
            self.profileFlag = False
            self.outCompress = "gzip"

        if self.rStep > 0:
            self.makedir = False

        return
