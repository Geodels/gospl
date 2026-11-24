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
        self._readGravity()
        self._readIce()
        self._readGroundwater()
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

    @staticmethod
    def _concatRows(df, rows, columns, prepend=False):
        """
        Add event rows to a forcing DataFrame (``tecdata``, ``raindata``,
        ``evapdata``, ``sedfacdata``, ``tedata``).

        Rebuilds the frame from records instead of calling ``pd.concat``:
        these frames have columns that are all-None for some events (e.g.
        ``rzA``/``rzB`` for a rain map), and pandas 2.x emits a
        FutureWarning, once per MPI rank, when concatenating those. Same
        values, same column names, fresh 0..n-1 index.

        :arg df: existing DataFrame
        :arg rows: list of dicts (one per new row)
        :arg columns: column names (part of the forcing API)
        :arg prepend: put ``rows`` before ``df`` instead of after
        :return: new DataFrame
        """
        old = df.to_dict("records")
        return pd.DataFrame(rows + old if prepend else old + rows, columns=columns)

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
        # Cyclic boundaries are SERIAL-ONLY for now. The periodic cylinder seam,
        # once split by the DMPlex partitioner, makes the finite-volume operators
        # ill-conditioned: the pit-deposition diffusion does not converge (the run
        # hangs) and the horizontal-advection solve diverges (elevation collapses
        # to ~0). Flow and pit-filling themselves are fine, but the sediment and
        # advection solves are not — so a parallel cyclic run produces a hang or
        # garbage rather than a wrong-but-finite answer. Fail fast with a clear
        # message instead. (A parallel-safe cylinder-seam discretisation is a
        # tracked follow-up.)
        if 'c' in self.boundCond and petsc4py.PETSc.COMM_WORLD.Get_size() > 1:
            raise ValueError(
                "Cyclic (periodic) boundary conditions (bc='%s') are currently "
                "supported only in SERIAL (one MPI rank). The periodic mesh seam "
                "does not yet partition correctly, so a parallel run would hang or "
                "produce invalid results. Run this model on a single rank, or use "
                "open/fixed/wall boundaries for parallel runs." % self.boundCond
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

        # Optional remap interval for `advect: interp`, in mesh spacings: the
        # displacement of each tectonic event is applied in remaps of about
        # this many cell widths instead of once at the event's end (None:
        # once per event, the historical behaviour). Remaps of half a cell or
        # more carry the drainage with the motion; sub-cell remaps (and the
        # per-step Eulerian schemes) do not (2026-11 lateral-motion test).
        self.interpRemap = domainDict.get("interp_remap")
        if self.interpRemap is not None:
            self.interpRemap = float(self.interpRemap)
            if self.interpRemap <= 0.0:
                raise ValueError("domain: interp_remap must be > 0 (mesh spacings).")
            if self.advscheme != 0:
                raise ValueError(
                    "domain: interp_remap only applies to advect: 'interp' "
                    "(got advect: '%s')." % advscheme
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
            # Optional channel threshold (m3/yr): land cells whose discharge
            # reaches it are rivers, held fixed in the hillslope diffusion; the
            # material hillslopes deliver to them joins the routed river load
            # (hillslope._channelDiffusion). None: diffusion everywhere, as
            # before. With rain P, channelQ / P is a drainage area.
            self.hillChannelQ = hillDict.get("channelQ")
        except KeyError:
            self.hillChannelQ = None
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
        if self.hillChannelQ is not None:
            self.hillChannelQ = float(self.hillChannelQ)
            if self.hillChannelQ <= 0.0:
                raise ValueError("diffusion: channelQ must be > 0 (m3/yr).")

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
            # soil accounting mode (DESIGN_SOIL_REGOLITH.md Option 2.5 step 2):
            #   "lumped"   (default) -> Lsoil is a soft cover that also absorbs
            #                           deposited sediment (legacy behaviour;
            #                           byte-identical to pre-step-2 runs).
            #   "regolith"           -> Lsoil is the WEATHERING-produced regolith
            #                           only; deposited sediment lives in the
            #                           stratigraphy (its erodibility via stratK),
            #                           not in Lsoil.
            self.regolithSoil = str(soilDict.get("mode", "lumped")).lower() == "regolith"
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
            self.regolithSoil = False
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
            tecdata = self._concatRows(tecdata, tmpTec, ["start", "end", "tMap", "zMap", "hMap"])

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
                tecdata = self._concatRows(tecdata, tmpTec, ["start", "end", "tMap", "zMap", "hMap"], prepend=True)
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
            sedfacdata = self._concatRows(sedfacdata, tmpErof, ["start", "sUni", "sMap", "sKey"])

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
                sedfacdata = self._concatRows(sedfacdata, tmpSedF, ["start", "sUni", "sMap", "sKey"], prepend=True)
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
            tedata = self._concatRows(tedata, tmpTe, ["start", "tUni", "tMap", "tKey"])

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
                tedata = self._concatRows(tedata, tmpT, ["start", "tUni", "tMap", "tKey"], prepend=True)
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
            raindata = self._concatRows(raindata, tmpRain, ["start", "rUni", "rzA", "rzB", "rMap", "rKey"])

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
            evapdata = self._concatRows(evapdata, tmpEvap, ["start", "eUni", "eMap", "eKey"])

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
                # Opt-in "losing stream": evaporation debited from the
                # accumulated discharge (not just local runoff). A global
                # numeric mode — enabled if ANY climate event requests it.
                if rainSort[k].get("evap_stream"):
                    self.evapStream = True

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
                raindata = self._concatRows(raindata, tmpRain, ["start", "rUni", "rzA", "rzB", "rMap", "rKey"], prepend=True)
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
                    evapdata = self._concatRows(evapdata, tmpEvap, ["start", "eUni", "eMap", "eKey"], prepend=True)
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
        # Insert a dedicated infinite-bedrock sentinel layer BENEATH the
        # file-provided initial layers (only consulted when an `npstrata` file
        # is given; the no-file path always builds a sentinel). Default False
        # keeps existing file runs unchanged, where the deepest file layer is
        # itself the erosion floor.
        self.bedrock_sentinel = False
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
            self.bedrock_sentinel = bool(
                strataDict.get("bedrock_sentinel", self.bedrock_sentinel)
            )

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

        self._readLithologies()
        self._readWaves()
        self._extraProvenance()

        return

    def _readLithologies(self):
        """
        Build the lithology table (docs/DESIGN_CARBONATES.md, section 2.7).

        Coarse clastic sediment is the implicit remainder of every layer. The
        other lithologies each carry their own pile, porosity curve,
        erodibility and diffusivity factors and grain-size behaviour
        (``size``: ``fine`` or ``coarse``, which decides how pits, the marine
        spreading and the gravity currents treat them). Under dual lithology
        the fine sediment is lithology 0, built from the ``fine:`` block and
        the dual keys, so a dual run is unchanged. An optional
        ``strata: lithologies:`` list appends more entries::

            lithologies:
              - {name: carb_mud, size: fine, phi0: 0.7, z0: 1500.,
                 k_factor: 1., diff_factor: 1., bedrock_frac: 0.}

        Sets ``litNb`` (0 without dual lithology), ``litName``, ``litFine``
        (bool, fine-like), ``litPhi0``, ``litZ0``, ``litK``, ``litD`` and
        ``litBed`` (share of the infinite-bedrock sentinel).
        """
        self.litNb = 0
        self.litName, self.litFine = [], []
        self.litPhi0, self.litZ0, self.litK, self.litD, self.litBed = [], [], [], [], []
        self.carbOn = False
        self.pelOn = False
        self.sstOn = False
        if not self.stratLith:
            if "pelagic" in self.input:
                raise ValueError(
                    "pelagic: needs dual lithology (strata: dual: True) and "
                    "stratigraphy (a positive `strat` interval): the rain is "
                    "recorded as the carb_mud and fine lithologies."
                )
            if "carbonate" in self.input:
                raise ValueError(
                    "carbonate: needs dual lithology (strata: dual: True) and "
                    "stratigraphy (a positive `strat` interval): the carbonate "
                    "sediment is recorded as its own lithologies."
                )
            return

        def add(name, fine, phi0, z0, k, d, bed):
            self.litName.append(name)
            self.litFine.append(bool(fine))
            self.litPhi0.append(phi0)
            self.litZ0.append(z0)
            self.litK.append(k)
            self.litD.append(d)
            self.litBed.append(bed)

        add("fine", True, self.phi0f, self.z0f, self.fine_k_factor,
            self.fine_diff_factor, 1.0 - self.bedrock_coarse_frac)

        extra = self.input["strata"].get("lithologies", None) or []
        if not isinstance(extra, (list, tuple)):
            raise ValueError("strata: lithologies must be a list of mappings.")
        for k, ent in enumerate(extra):
            if not isinstance(ent, dict) or "name" not in ent:
                raise ValueError(
                    "strata: lithologies entry %d needs at least a 'name'." % k
                )
            name = str(ent["name"])
            if name in self.litName or name == "coarse":
                raise ValueError("strata: lithology name '%s' is already used." % name)
            size = str(ent.get("size", "fine")).lower()
            if size not in ("fine", "coarse"):
                raise ValueError(
                    "strata: lithology '%s' size must be 'fine' or 'coarse'." % name
                )
            add(name, size == "fine", float(ent.get("phi0", self.phi0f)),
                float(ent.get("z0", self.z0f)), float(ent.get("k_factor", 1.0)),
                float(ent.get("diff_factor", 1.0)), float(ent.get("bedrock_frac", 0.0)))
        self._readSST()
        self._readCarbonate(add)
        self._readPelagic(add)
        self._readCarbSupply()
        if sum(self.litBed) > 1.0 + 1.0e-12:
            raise ValueError(
                "strata: the bedrock shares of the lithologies (fine: "
                "1 - bedrock_coarse_frac, others: bedrock_frac) sum to more than 1."
            )
        self.litNb = len(self.litName)

        return

    # Default carbonate factories (docs/DESIGN_CARBONATES.md 2.1): the three
    # Bosscher & Schlager factories of CarboKitten (Hidding et al. 2026,
    # Table 1). The framework and sand fractions are illustrative defaults.
    CARB_FACTORIES = (
        {"name": "euphotic", "growth": 500.0, "extinction": 0.8, "saturation": 60.0,
         "framework": 0.5, "sand": 0.7},
        {"name": "oligophotic", "growth": 400.0, "extinction": 0.1, "saturation": 60.0,
         "framework": 0.1, "sand": 0.5},
        {"name": "aphotic", "growth": 100.0, "extinction": 0.005, "saturation": 60.0,
         "framework": 0.0, "sand": 0.3},
    )
    # Default carbonate lithologies: surface porosity and e-folding depth of
    # the porosity-depth curve (grainstone, lime mud, boundstone; orders of
    # magnitude from published carbonate compaction trends, to calibrate).
    CARB_LITHOLOGIES = {
        "carb_sand": {"size": "coarse", "phi0": 0.45, "z0": 2500.0, "k_factor": 1.0,
                      "diff_factor": 1.0},
        "carb_mud": {"size": "fine", "phi0": 0.65, "z0": 1500.0, "k_factor": 1.0,
                     "diff_factor": 1.0},
        "reef": {"size": "coarse", "phi0": 0.2, "z0": 5000.0, "k_factor": 0.1,
                 "diff_factor": 0.1},
    }
    CARB_CONTROLS = ("depth", "slope", "sedrate", "exposure", "temperature", "wave",
                     "solute")

    def _readWaves(self):
        """
        Parse the optional ``waves:`` block (wave energy and shadows,
        docs/DESIGN_CARBONATES.md 2.10, phase C7)::

            waves:
                climate:                         # directions the waves come FROM
                  - {from: 90., weight: 0.7}     # degrees, 0 = N, 90 = E
                  - {from: 180., weight: 0.3}    # or from: zonal (global belts),
                                                 # or from_map: [file, key] (degrees)
                fetch: 3.e5          # growth to full energy over open water (m);
                                     # required on a global mesh
                decay: 1.e5          # attenuation length in deep water (m)
                shoal_depth: 20.     # depth below which the waves feel the bottom (m)
                shoal_length: 2.e3   # attenuation length at the shoreline (m)

        The weights are normalised to sum to 1. Sets ``waveOn``,
        ``waveClimate`` (list of ``(from, weight)``), ``waveDecay``,
        ``waveShoalDepth``, ``waveShoalLength`` and ``waveFetch`` (m or None).
        A climate entry's direction is a bearing (degrees), ``"zonal"`` (trade
        winds, westerlies, polar easterlies by latitude; global meshes) or a
        ``(file, key)`` map of bearings per vertex.
        """
        self.waveOn = False
        wd = self.input.get("waves", None)
        if wd is None:
            return
        if not isinstance(wd, dict):
            raise ValueError("waves: must be a mapping with a `climate` list.")
        clim = wd.get("climate", None)
        if not isinstance(clim, (list, tuple)) or len(clim) == 0:
            raise ValueError("waves: `climate` must be a non-empty list of "
                             "{from: degrees, weight: w}.")
        out = []
        for k, c in enumerate(clim):
            if not isinstance(c, dict) or ("from" not in c and "from_map" not in c):
                raise ValueError("waves: climate entry %d needs `from` (degrees, or "
                                 "`zonal`) or `from_map` ([file, key])." % k)
            w = float(c.get("weight", 1.0))
            if w < 0.0:
                raise ValueError("waves: climate entry %d weight must be >= 0." % k)
            if "from_map" in c:
                fm = c["from_map"]
                if not isinstance(fm, (list, tuple)) or len(fm) != 2:
                    raise ValueError("waves: climate entry %d from_map must be "
                                     "[file, key]." % k)
                frm = (str(fm[0]), str(fm[1]))
            elif str(c["from"]).lower() == "zonal":
                frm = "zonal"
            else:
                frm = float(c["from"]) % 360.0
            out.append((frm, w))
        tot = sum(w for _, w in out)
        if tot <= 0.0:
            raise ValueError("waves: the climate weights sum to zero.")
        self.waveClimate = [(f, w / tot) for f, w in out]
        self.waveDecay = float(wd.get("decay", 1.0e5))
        self.waveShoalDepth = float(wd.get("shoal_depth", 20.0))
        self.waveShoalLength = float(wd.get("shoal_length", 2.0e3))
        fe = wd.get("fetch", None)
        self.waveFetch = None if fe is None else float(fe)
        if self.waveFetch is not None and self.waveFetch <= 0.0:
            raise ValueError("waves: fetch must be positive (m).")
        for key, v in (("decay", self.waveDecay), ("shoal_depth", self.waveShoalDepth),
                       ("shoal_length", self.waveShoalLength)):
            if v <= 0.0:
                raise ValueError("waves: %s must be positive (m)." % key)
        self.waveOn = True

    def _readCarbSupply(self):
        """
        Parse the optional ``carbonate_supply:`` block (supply-limited
        carbonate burial, docs/DESIGN_CARBONATES.md 2.10, phase C6)::

            carbonate_supply:
                species: carbonate   # geochemistry species feeding the burial
                yield: 1.0           # m3 of carbonate per unit of solute
                pelagic: True        # the pelagic rain shares the supply

        The marine reservoir of the species (groundwater geochemistry with
        ``river_load`` and ``marine_coupling``) caps the carbonate buried each
        step. ``yield`` follows the convention of the duricrust
        ``solid_volume``: volume of rock per unit of solute, in the
        concentration units of the geochemistry. Sets ``carbSupply`` (dict or
        None).
        """
        self.carbSupply = None
        cs = self.input.get("carbonate_supply", None)
        if cs is None:
            return
        if not isinstance(cs, dict) or "species" not in cs:
            raise ValueError("carbonate_supply: needs a `species` (geochemistry "
                             "tracer name) and a `yield`.")
        y = float(cs.get("yield", 0.0))
        if y <= 0.0:
            raise ValueError("carbonate_supply: yield must be positive (m3 of "
                             "carbonate per unit of solute).")
        if not (self.carbOn or self.pelOn):
            raise ValueError("carbonate_supply: needs a `carbonate` or a `pelagic` "
                             "block.")
        self.carbSupply = {"species": str(cs["species"]), "yield": y,
                           "pelagic": bool(cs.get("pelagic", True))}

    def _readSST(self):
        """
        Parse the optional ``sst:`` block (sea-surface temperature, degrees C,
        docs/DESIGN_CARBONATES.md 2.10, phase C5), used by the ``temperature``
        response of the carbonate factories and the pelagic rain::

            sst:
              events:                       # maps (global) or uniform values
                - {start: 0., map: ['sst0', 't']}
                - {start: 1.e6, uniform: 24.}
              curve: 'sst_curve.csv'        # time (yr), SST (deg C)

        The events are a step function in time (each holds until the next
        starts); the curve is interpolated linearly (constant beyond its
        ends). With both, the curve is an anomaly added to the event values;
        with the curve alone it is the (uniform) SST. Sets ``sstOn``,
        ``sstEvents`` (list of ``{start, val}``, ``val`` a float or
        ``(file, key)``) and ``sstCurve`` (``(t, v)`` or None).
        """
        sd = self.input.get("sst", None)
        self.sstEvents, self.sstCurve = [], None
        if sd is None:
            return
        if not isinstance(sd, dict):
            raise ValueError("sst: must be a mapping with `events` and/or `curve`.")
        for k, e in enumerate(sorted(sd.get("events", None) or [],
                                     key=lambda x: float(x.get("start", 0.0)))):
            if not isinstance(e, dict) or "start" not in e:
                raise ValueError("sst: event %d needs a start time." % k)
            if "map" in e:
                mp = e["map"]
                if not isinstance(mp, (list, tuple)) or len(mp) != 2:
                    raise ValueError("sst: event %d map must be [file, key]." % k)
                val = (str(mp[0]), str(mp[1]))
            elif "uniform" in e:
                val = float(e["uniform"])
            else:
                raise ValueError("sst: event %d needs `uniform` or `map`." % k)
            self.sstEvents.append({"start": float(e["start"]), "val": val})
        if sd.get("curve", None) is not None:
            try:
                data = np.loadtxt(str(sd["curve"]), delimiter=None if not
                                  str(sd["curve"]).endswith(".csv") else ",", ndmin=2)
            except Exception as exc:
                raise ValueError("sst: cannot read the curve %s (%s)." % (sd["curve"], exc))
            if data.shape[1] < 2 or data.shape[0] < 1 or np.any(np.diff(data[:, 0]) <= 0.0):
                raise ValueError(
                    "sst: the curve needs two columns (time in yr, increasing; SST "
                    "in deg C).")
            self.sstCurve = (data[:, 0].copy(), data[:, 1].copy())
        if not self.sstEvents and self.sstCurve is None:
            raise ValueError("sst: give `events`, a `curve`, or both.")
        self.sstOn = True

    def _readCarbonate(self, add):
        """
        Parse the optional ``carbonate:`` block (benthic carbonate production,
        docs/DESIGN_CARBONATES.md, phase C1) and add its lithologies to the
        table through ``add`` (called from ``_readLithologies``).

        Each factory ``f`` produces, below sea level, at the Bosscher &
        Schlager (1992) rate ``growth * tanh(I0 exp(-extinction w) /
        saturation)`` (``w`` the water depth), times an optional habitat map
        and response tables. Its ``framework`` fraction is built in place
        (lithology ``reef``), the rest is loose sediment split into
        ``carb_sand`` (``sand`` fraction, coarse-like) and ``carb_mud``
        (fine-like), handed to the marine routing. Rates are read in m/Myr
        and kept in m/yr.

        Sets ``carbOn``, ``carbI0``, ``carbMinDepth``, ``carbExpLength`` and
        ``carbFact`` (list
        of dicts with ``g`` (m/yr), ``k``, ``Ik``, ``phi`` (framework),
        ``sand``, ``habitat`` ([file, key] or None) and ``resp`` (list of
        (control, x, f))).
        """
        cdict = self.input.get("carbonate", None)
        if cdict is None:
            return
        if cdict is True or cdict == {}:
            cdict = {}
        if not isinstance(cdict, dict):
            raise ValueError("carbonate: must be a mapping.")
        self.carbOn = True
        self.carbI0 = float(cdict.get("insolation", 400.0))
        self.carbMinDepth = float(cdict.get("min_depth", 0.0))
        # Length (m) of the open-water index `exposure` (elevation minus its
        # smoothing over that length; positive on margins and rims, negative
        # in lagoons). Needed only when a response uses it.
        el = cdict.get("exposure_length", None)
        self.carbExpLength = None if el is None else float(el)
        if self.carbExpLength is not None and self.carbExpLength <= 0.0:
            raise ValueError("carbonate: exposure_length must be positive (m).")
        if self.carbI0 <= 0.0:
            raise ValueError("carbonate: insolation must be positive.")
        facts = cdict.get("factories", None)
        if facts is None:
            facts = [dict(f) for f in self.CARB_FACTORIES]
        if not isinstance(facts, (list, tuple)) or len(facts) == 0:
            raise ValueError("carbonate: factories must be a non-empty list.")
        self.carbFact = []
        names = set()
        for k, f in enumerate(facts):
            if not isinstance(f, dict):
                raise ValueError("carbonate: factory %d must be a mapping." % k)
            name = str(f.get("name", "factory%d" % k))
            if name in names:
                raise ValueError("carbonate: factory name '%s' is used twice." % name)
            names.add(name)
            g = float(f.get("growth", 0.0))
            kx = float(f.get("extinction", 0.1))
            ik = float(f.get("saturation", 60.0))
            phi = float(f.get("framework", 0.0))
            sand = float(f.get("sand", 0.5))
            if g < 0.0 or kx < 0.0 or ik <= 0.0:
                raise ValueError(
                    "carbonate: factory '%s' needs growth >= 0, extinction >= 0 "
                    "and saturation > 0." % name)
            if not (0.0 <= phi <= 1.0 and 0.0 <= sand <= 1.0):
                raise ValueError(
                    "carbonate: factory '%s' framework and sand must be between "
                    "0 and 1." % name)
            resp = []
            solute = None
            for r in f.get("response", None) or []:
                ctl = str(r.get("control", ""))
                if ctl not in self.CARB_CONTROLS:
                    raise ValueError(
                        "carbonate: factory '%s' response control '%s' is not one "
                        "of %s." % (name, ctl, ", ".join(self.CARB_CONTROLS)))
                x = np.asarray(r.get("x", []), dtype=float)
                fv = np.asarray(r.get("f", []), dtype=float)
                if x.size < 2 or x.size != fv.size or np.any(np.diff(x) <= 0.0):
                    raise ValueError(
                        "carbonate: factory '%s' response '%s' needs x and f of "
                        "the same length (>= 2), x increasing." % (name, ctl))
                if np.any(fv < 0.0):
                    raise ValueError(
                        "carbonate: factory '%s' response '%s' factors must be "
                        ">= 0." % (name, ctl))
                if ctl == "solute":
                    # coastal plume of a geochemistry species (phase C6)
                    if solute is not None:
                        raise ValueError(
                            "carbonate: factory '%s' has more than one solute "
                            "response." % name)
                    if "species" not in r or float(r.get("length", 0.0)) <= 0.0:
                        raise ValueError(
                            "carbonate: factory '%s' solute response needs a "
                            "`species` and a positive `length` (m)." % name)
                    solute = (str(r["species"]), float(r["length"]))
                resp.append((ctl, x, fv))
            hab = f.get("habitat", None)
            if hab is not None and (not isinstance(hab, (list, tuple)) or len(hab) != 2):
                raise ValueError(
                    "carbonate: factory '%s' habitat must be [file, key]." % name)
            self.carbFact.append({"name": name, "g": g * 1.0e-6, "k": kx, "Ik": ik,
                                  "phi": phi, "sand": sand, "habitat": hab,
                                  "resp": resp, "solute": solute})

        if self.carbExpLength is None and any(
                ctl == "exposure" for f in self.carbFact for ctl, _, _ in f["resp"]):
            raise ValueError(
                "carbonate: a response uses `exposure`; set carbonate: "
                "exposure_length (m), the length over which the open-water "
                "index compares the seafloor with its surroundings.")
        if "waves" not in self.input and any(
                ctl == "wave" for f in self.carbFact for ctl, _, _ in f["resp"]):
            raise ValueError(
                "carbonate: a response uses `wave`; add a `waves:` block (the "
                "wave climate).")
        if not self.sstOn and any(
                ctl == "temperature" for f in self.carbFact for ctl, _, _ in f["resp"]):
            raise ValueError(
                "carbonate: a response uses `temperature`; add an `sst:` block "
                "(sea-surface temperature maps or curve).")
        over = cdict.get("lithologies", None) or {}
        if not isinstance(over, dict):
            raise ValueError("carbonate: lithologies must be a mapping by name.")
        needed = ["carb_sand", "carb_mud"]
        if any(f["phi"] > 0.0 for f in self.carbFact):
            needed.append("reef")
        for key in over:
            if key not in self.CARB_LITHOLOGIES:
                raise ValueError(
                    "carbonate: lithologies has no '%s' (carb_sand, carb_mud, reef)."
                    % key)
        self.carbLit = {}
        for name in needed:
            if name in self.litName:
                raise ValueError(
                    "carbonate: lithology name '%s' is reserved for the carbonate "
                    "sediment; rename the strata: lithologies entry." % name)
            prop = dict(self.CARB_LITHOLOGIES[name])
            prop.update(over.get(name, {}) or {})
            self.carbLit[name] = len(self.litName)
            add(name, prop["size"] == "fine", float(prop["phi0"]), float(prop["z0"]),
                float(prop["k_factor"]), float(prop["diff_factor"]), 0.0)

    def _readPelagic(self, add):
        """
        Parse the optional ``pelagic:`` block (settled pelagic carbonate and
        hemipelagic clay, docs/DESIGN_CARBONATES.md 2.8, phase C2): a list of
        events sorted by ``start`` (yr), each holding until the next one::

            pelagic:
              - start: 0.
                carbonate: 5.        # pelagic carbonate rain (m/Myr) or [file, key]
                lysocline: 3500.     # depth (m) where dissolution starts
                ccd: 4500.           # carbonate compensation depth (m)
                clay: 1.             # hemipelagic clay (m/Myr) or [file, key]
                clay_decay: 5.e4     # optional: clay * exp(-d/L) from the coast (m)
                resuspend: False     # optional: rain into the marine diffusion

        The carbonate rain settles as ``carb_mud`` with the fraction
        preserved going linearly from 1 at the lysocline to 0 at the CCD; the
        clay settles as ``fine``. Adds ``carb_mud`` to the lithology table if
        the carbonate block did not. Sets ``pelOn`` and ``pelEvents`` (list of
        dicts: ``start``, ``carb``/``clay`` (m/yr scalar or [file, key] with
        the m/Myr scale), ``lys``, ``ccd``, ``decay`` (m or None) and
        ``resuspend``).
        """
        ev = self.input.get("pelagic", None)
        if ev is None:
            return
        if not isinstance(ev, (list, tuple)) or len(ev) == 0:
            raise ValueError("pelagic: must be a list of events with a start time.")
        for k, e in enumerate(ev):
            if not isinstance(e, dict) or "start" not in e:
                raise ValueError("pelagic: event %d needs a start time." % k)
        self.pelEvents = []
        anyCarb = False

        def rate(val, key, k):
            if val is None:
                return 0.0
            if isinstance(val, (list, tuple)):
                if len(val) != 2:
                    raise ValueError(
                        "pelagic: event %d %s map must be [file, key]." % (k, key))
                return (str(val[0]), str(val[1]))
            v = float(val)
            if v < 0.0:
                raise ValueError("pelagic: event %d %s must be >= 0." % (k, key))
            return v * 1.0e-6

        for k, e in enumerate(sorted(ev, key=lambda x: float(x["start"]))):
            carb = rate(e.get("carbonate"), "carbonate", k)
            clay = rate(e.get("clay"), "clay", k)
            lys = float(e.get("lysocline", 1.0e12))
            ccd = float(e.get("ccd", lys))
            hasCarb = isinstance(carb, tuple) or carb > 0.0
            if hasCarb and "lysocline" not in e:
                raise ValueError(
                    "pelagic: event %d has carbonate rain but no lysocline." % k)
            if ccd < lys:
                raise ValueError(
                    "pelagic: event %d needs ccd >= lysocline (depths in m)." % k)
            anyCarb = anyCarb or hasCarb
            decay = e.get("clay_decay", None)
            if decay is not None:
                decay = float(decay)
                if decay <= 0.0:
                    raise ValueError(
                        "pelagic: event %d clay_decay must be a positive length (m)." % k)
            temp = e.get("temperature", None)
            if temp is not None:
                if not self.sstOn:
                    raise ValueError(
                        "pelagic: event %d has a temperature table; add an `sst:` "
                        "block." % k)
                tx = np.asarray(temp.get("x", []), dtype=float)
                tf = np.asarray(temp.get("f", []), dtype=float)
                if tx.size < 2 or tx.size != tf.size or np.any(np.diff(tx) <= 0.0) \
                        or np.any(tf < 0.0):
                    raise ValueError(
                        "pelagic: event %d temperature needs x (increasing) and f "
                        "(>= 0) of the same length (>= 2)." % k)
                temp = (tx, tf)
            self.pelEvents.append({"start": float(e["start"]), "carb": carb,
                                   "clay": clay, "lys": lys, "ccd": ccd,
                                   "decay": decay, "temp": temp,
                                   "resuspend": bool(e.get("resuspend", False))})
        self.pelOn = True
        if anyCarb and "carb_mud" not in self.litName:
            prop = self.CARB_LITHOLOGIES["carb_mud"]
            if not hasattr(self, "carbLit"):
                self.carbLit = {}
            self.carbLit["carb_mud"] = len(self.litName)
            add("carb_mud", True, float(prop["phi0"]), float(prop["z0"]),
                float(prop["k_factor"]), float(prop["diff_factor"]), 0.0)

    def _extraProvenance(self):
        """
        Parse the optional in-model **sediment-provenance tracer** block. N
        source-rock classes are carried (as a passive label — no physics
        feedback) through erosion, transport, deposition and the stratigraphic
        record, giving conservation-exact, recycling-aware provenance per layer.
        See ``docs/DESIGN_PROVENANCE.md`` §6.

        Sets ``self.provOn`` (master opt-in), ``self.provNb`` (class count) and
        the per-vertex source-class source (``uniform`` scalar or ``source``
        ``[file, key]`` map, resolved post-mesh in ``readStratLayers``). All
        inert while ``provOn`` is False; like dual lithology it requires
        stratigraphy (``stratNb > 0``).

        Note: the copper-fertility weighting (``cu_weight``) is **not** a model
        input — it is a post-processing parameter supplied to the standalone
        ``gospl.analyse.provenance`` tool (its ``--cu-weights`` argument), since
        the fertility is a what-if weighting applied to the recorded per-class
        deposits and is best varied without re-running the model.
        """
        self.provOn = False
        self.provNb = 0
        self._provSourceUniform = None
        self._provSourceMap = None
        try:
            provDict = self.input["provenance"]
        except KeyError:
            return

        self.provNb = int(provDict.get("classes", 0))
        self._provSourceUniform = provDict.get("uniform")
        self._provSourceMap = provDict.get("source")          # [file, key] or None

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

    def _readGravity(self):
        """
        Parse the optional marine **gravity-current** block (``gravity:``), see
        ``docs/DESIGN_GRAVITY_CURRENTS.md``.

        Sets ``self.gravOn`` (master opt-in, False when the block is absent so
        every gravity path is gated out), the run-out lengths ``gravRunout``
        (coarse or single fraction, m) and ``gravRunoutFine`` (fine fraction
        under dual lithology, m), the bypass conditions ``gravBypassDepth`` (m
        below sea level, default 50) and ``gravBypassSlope`` (off by default),
        the routing network (``gravNetwork`` 'raw' or 'smoothed',
        ``gravReceivers``, ``gravFlowExp``), the slope erosion
        (``gravErosionLength``, off when absent, ``gravErosionSlope``,
        ``gravMaxErosion`` m per step, ``gravErosionExp``/``gravErosionFlux`` non-linear flux
        dependence), and the trigger: a fixed fraction
        ``gravFraction`` of the river load, or (when ``fraction`` is absent) the
        flood-fraction parameters ``gravCcrit`` (kg/m3), ``gravRating`` (b > 1),
        ``gravSigma`` (log-std of daily discharge) and ``gravRhoS`` (kg/m3).
        """
        self.gravOn = False
        gravDict = self.input.get("gravity")
        if gravDict is None:
            return
        if not isinstance(gravDict, dict):
            raise ValueError("gravity: expects a block of keys (runout, fraction or flood, ...)")

        self.gravRunout = float(gravDict.get("runout", 5.0e4))
        self.gravRunoutFine = float(gravDict.get("runout_fine", 4.0 * self.gravRunout))
        if self.gravRunout <= 0.0 or self.gravRunoutFine <= 0.0:
            raise ValueError("gravity: runout and runout_fine must be positive (m)")

        # Bypass: no deposition shallower than bypass_depth (m below sea level,
        # e.g. storm wave base; 0 = off) or on slopes steeper than bypass_slope
        # (dimensionless; None = off, slopes depend on the mesh resolution).
        self.gravBypassDepth = float(gravDict.get("bypass_depth", 50.0))
        bs = gravDict.get("bypass_slope")
        self.gravBypassSlope = None if bs is None else float(bs)
        if self.gravBypassDepth < 0.0:
            raise ValueError("gravity: bypass_depth must be >= 0 (m below sea level)")
        if self.gravBypassSlope is not None and self.gravBypassSlope <= 0.0:
            raise ValueError("gravity: bypass_slope must be positive")

        # Gravity network (phase G3b): 'raw' = own routing on the eps-filled,
        # UNSMOOTHED bathymetry (concentrates on slopes, spreads on flat floors),
        # 'smoothed' = the clinoform cascade's network (dMat1).
        self.gravNetwork = str(gravDict.get("network", "raw")).lower()
        if self.gravNetwork not in ("raw", "smoothed"):
            raise ValueError("gravity: network must be 'raw' or 'smoothed'")
        self.gravReceivers = int(gravDict.get("receivers", 6))
        self.gravFlowExp = float(gravDict.get("flow_exp", 2.0))
        if not 1 <= self.gravReceivers <= 12:
            raise ValueError("gravity: receivers must lie between 1 and 12")
        if self.gravFlowExp < 0.0:
            raise ValueError("gravity: flow_exp must be >= 0")
        # Fine fraction on its own network (dual lithology): coarse sediment
        # follows concentrated pathways (channel belts), fine spreads (lobes).
        self.gravReceiversFine = int(gravDict.get("receivers_fine", self.gravReceivers))
        self.gravFlowExpFine = float(gravDict.get("flow_exp_fine", self.gravFlowExp))
        if not 1 <= self.gravReceiversFine <= 12 or self.gravFlowExpFine < 0.0:
            raise ValueError("gravity: receivers_fine in 1..12 and flow_exp_fine >= 0")
        # Smoothing passes of the routing surface and of the slope (tapers and
        # erosion), so pathways and deposits follow the seafloor at the scale
        # of a few cells rather than cell-to-cell noise.
        self.gravNetworkSmooth = int(gravDict.get("network_smooth", 0))
        self.gravSlopeSmooth = int(gravDict.get("slope_smooth", 0))
        # Lateral spreading length (m) of each step's gravity deposit: gives
        # concentrated belts a width (number of passes set from the cell size).
        # Absent: two mean cell sizes (set at run time); 0 turns it off.
        ds = gravDict.get("deposit_spread")
        self.gravDepositSpread = None if ds is None else float(ds)
        if self.gravDepositSpread is not None and self.gravDepositSpread < 0.0:
            raise ValueError("gravity: deposit_spread must be >= 0 (m)")
        # Share of the diverted fraction applied to the coarse load (dual
        # lithology): gravity currents carry mostly the suspended load.
        self.gravCoarseFactor = float(gravDict.get("coarse_factor", 0.3))
        if not 0.0 <= self.gravCoarseFactor <= 1.0:
            raise ValueError("gravity: coarse_factor must lie between 0 and 1")
        if self.gravNetworkSmooth < 0 or self.gravSlopeSmooth < 0:
            raise ValueError("gravity: network_smooth and slope_smooth must be >= 0")

        # Slope erosion (phase G3b): off unless erosion_length (m) is given.
        el = gravDict.get("erosion_length")
        self.gravErosionLength = None if el is None else float(el)
        self.gravErosionSlope = float(gravDict.get("erosion_slope", 0.02))
        self.gravMaxErosion = float(gravDict.get("max_erosion", 10.0))
        if self.gravErosionLength is not None and self.gravErosionLength <= 0.0:
            raise ValueError("gravity: erosion_length must be positive (m)")
        if self.gravErosionSlope <= 0.0 or self.gravMaxErosion <= 0.0:
            raise ValueError("gravity: erosion_slope and max_erosion must be positive")
        # Smoothed slope (passes of neighbour averaging) and non-linear flux
        # dependence eps ~ (Q/Q_e)^(a-1) (a = erosion_exp, Q_e = erosion_flux m3/yr).
        self.gravErosionExp = float(gravDict.get("erosion_exp", 1.0))
        self.gravErosionFlux = float(gravDict.get("erosion_flux", 1.0e4))
        if self.gravErosionExp < 1.0 or self.gravErosionFlux <= 0.0:
            raise ValueError("gravity: erosion_exp must be >= 1 and erosion_flux positive (m3/yr)")

        # Slope-failure trigger (phase G3): off unless failure_slope is given.
        fs = gravDict.get("failure_slope")
        self.gravFailureSlope = None if fs is None else float(fs)
        self.gravMaxFailure = float(gravDict.get("max_failure", 10.0))
        if self.gravFailureSlope is not None and self.gravFailureSlope <= 0.0:
            raise ValueError("gravity: failure_slope must be positive")
        if self.gravMaxFailure <= 0.0:
            raise ValueError("gravity: max_failure must be positive (m per step)")

        frac = gravDict.get("fraction")
        self.gravFraction = None if frac is None else float(frac)
        if self.gravFraction is not None and not 0.0 <= self.gravFraction <= 1.0:
            raise ValueError("gravity: fraction must lie between 0 and 1")

        flood = gravDict.get("flood") or {}
        self.gravCcrit = float(flood.get("c_crit", 40.0))
        self.gravRating = float(flood.get("rating_exp", 1.8))
        self.gravSigma = float(flood.get("sigma", 1.2))
        self.gravRhoS = float(flood.get("rho_s", 2650.0))
        if self.gravFraction is None:
            if self.gravRating <= 1.0:
                raise ValueError("gravity: flood rating_exp must be > 1")
            if self.gravSigma <= 0.0 or self.gravCcrit <= 0.0 or self.gravRhoS <= 0.0:
                raise ValueError("gravity: flood c_crit, sigma and rho_s must be positive")

        self.gravOn = True
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
            self.ice_meltfac = iceDict.get("melt", 0.0)      # ablation in net balance (0 = accumulation-only; 1 = true net)
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
            self.ice_meltfac = 0.0
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

    def _readGroundwater(self):
        """
        Parse the optional ``groundwater:`` block (opt-in water table + generic
        duricrust). See ``docs/DESIGN_WATERTABLE_DURICRUST.md``. When the block is
        absent everything is inert (``gwOn = False``) and goSPL is byte-identical
        to a run without it. Phase 0 = parse + flags + defaults only (no solve).
        """

        try:
            gwDict = self.input["groundwater"]
            self.gwOn = True
            # Hydrology (implicit Dupuit-Boussinesq head solve).
            self.gwKsat = float(gwDict.get("Ksat", 1.0))            # K_h (m/yr)
            self.gwSpecificYield = float(gwDict.get("specific_yield", 0.1))  # S
            # z_bed depth below surface (m): a scalar, the string 'from_soil'
            # (tie to lHbed; requires soilSPL), OR a per-vertex map ``[file, key]``
            # (loaded in _GWMesh, which has locIDs — physically the regolith /
            # weathering-front depth varies in space). Kept raw here; resolved
            # later.
            base = gwDict.get("aquifer_base", 50.0)
            if isinstance(base, (list, tuple)):
                self._gwAquiferBaseMap = tuple(base)
                self.gwAquiferBase = None
            else:
                self._gwAquiferBaseMap = None
                self.gwAquiferBase = base
            self.gwBedrockDepth = float(gwDict.get("bedrock_depth", 0.0))
            self.gwMinSatThick = float(gwDict.get("min_sat_thickness", 1.0))
            # Infiltration fraction f_infil: a scalar, OR a per-vertex map
            # ``[file, key]`` (loaded in _GWMesh, which has locIDs). Physically it
            # varies with lithology / regolith / slope; see DESIGN §6.
            infil = gwDict.get("infiltration", 0.3)
            if isinstance(infil, (list, tuple)):
                self._gwInfilMap = tuple(infil)
                self.gwInfiltration = None
            else:
                self._gwInfilMap = None
                self.gwInfiltration = float(infil)
            self.gwConserveBaseflow = bool(gwDict.get("conserve_baseflow", True))
            # Opt-in lake ↔ aquifer VOLUME coupling (default off — the standard
            # fixed-head treatment). When on, the across-bed groundwater flux
            # debits/credits each lake's fill budget (DESIGN §15).
            self.gwLakeExchange = bool(gwDict.get("lake_exchange", False))
            # Recharge refinements (all opt-in, defaults = current behaviour):
            #   subglacial_recharge — fraction of glacial meltwater (iceMeltRiverL)
            #     that infiltrates under ice (0 = the ice gate zeroes recharge);
            #   fine_infil_factor — multiplier on f_infil for the fine end-member
            #     (dual lithology; 1 = no lithology dependence);
            #   infil_slope_ref — reference slope for f/(1+slope/ref) reduction on
            #     steep terrain (0 = off, no slope dependence).
            self.gwSubglacial = float(gwDict.get("subglacial_recharge", 0.0))
            self.gwFineInfilFactor = float(gwDict.get("fine_infil_factor", 1.0))
            self.gwInfilSlopeRef = float(gwDict.get("infil_slope_ref", 0.0))
            self.gwPicardIts = int(gwDict.get("picard_its", 3))
            self.gwSeepagePasses = int(gwDict.get("seepage_passes", 4))

            # Generic duricrust (capillary-fringe induration that armors K).
            duri = gwDict.get("duricrust", {}) or {}
            self.duriOn = bool(duri) if isinstance(duri, dict) else False
            self.duriFormRate = float(duri.get("form_rate", 1.0e-4))       # k_form (m/yr)
            self.duriMaxThick = float(duri.get("max_thickness", 5.0))      # duriH_max (m)
            self.duriFringeDepth = float(duri.get("fringe_depth", 3.0))    # d0 (m)
            self.duriFringeWidth = float(duri.get("fringe_width", 2.0))    # w (m)
            # Absolute-accumulation gate: restrict crust to groundwater DISCHARGE
            # zones (lateral flow converges, delivering dissolved load that emerges
            # and precipitates) rather than forming wherever the water table is
            # shallow (the in-situ / relative-accumulation default). Off ->
            # backwards-compatible (weight 1 everywhere). See DESIGN §3b.
            self.duriDischargeGate = bool(duri.get("discharge_gate", False))
            self.duriSupplyExp = float(duri.get("supply_exp", 1.0))        # p on (rain-evap)
            self.duriWeatherEa = float(duri.get("weather_Ea", 0.0))        # Arrhenius (0=off)
            self.duriArmorMax = float(duri.get("armor_max", 0.9))          # max K reduction 0..1
            self.duriArmorDiffusion = bool(duri.get("armor_diffusion", False))
            self.duriBreakRate = float(duri.get("break_rate", 1.0))        # k_break
            self.duriDecayRate = float(duri.get("decay_rate", 1.0e-6))     # k_decay (1/yr)
            # Weathering-supply coupling (proxy | rate | prodsoil) — Level-A knobs
            # parsed here, used from Phase 3 (see DESIGN §3a).
            weath = duri.get("weathering", {}) or {}
            self.duriWeatherMode = str(weath.get("mode", "proxy")).lower()
            self.duriWeatherCeq = float(weath.get("C_eq", 1.0))
            self.duriWeatherDw = float(weath.get("Dw", 1.0))
            self.duriWeatherL = float(weath.get("path_length", 20.0))
            self.duriWeatherability = weath.get("weatherability", 1.0)

            # Level-B conservative geochemistry (opt-in; DESIGN_WATERTABLE_GEOCHEM.md).
            # The solute is an `n_species` array (default 1): single-tracer ships
            # first, multi-tracer (calcrete/silcrete/ferricrete typing) is just a
            # longer `species:` list. Coupled multi-species aqueous equilibrium is
            # OUT of scope. Per-species params kept as lists here; gwplex builds the
            # numpy arrays + state. State alloc only in G0 (nothing solved yet).
            geo = gwDict.get("geochem", {}) or {}
            self.gwGeochemOn = bool(geo)
            species = geo.get("species", []) or []
            if self.gwGeochemOn and not species:
                species = [{}]                          # default single tracer
            self.gwNspecies = max(1, len(species))
            self.gwGeoName = [
                str(s.get("name", "solute%d" % i)) for i, s in enumerate(species)
            ] or ["solute0"]
            # weatherability: scalar OR a per-vertex `[file, key]` map (a), so
            # lithology can control which species a region yields; kept RAW here
            # (gwplex resolves it). `weatherability_by_class` + `weatherability_from`
            # give the table form (c): a per-(class, species) value gathered by a
            # per-vertex lithology label (currently `source_class`, i.e. provenance).
            self.gwGeoWeather = [s.get("weatherability", 1.0) for s in species] or [1.0]
            self.gwGeoWeatherByClass = [
                s.get("weatherability_by_class", None) for s in species
            ] or [None]
            self.gwWeatherFrom = geo.get("weatherability_from", None)
            # (b) a standalone per-vertex integer lithology map `[file, key]`
            # (independent of provenance) — the label for `weatherability_by_class`.
            self._gwLithoMap = geo.get("lithology", None)
            self.gwGeoCsat = [float(s.get("c_sat", 1.0)) for s in species] or [1.0]
            self.gwGeoPrecip = [float(s.get("precip_rate", 1.0)) for s in species] or [1.0]
            self.gwGeoVsolid = [float(s.get("solid_volume", 1.0)) for s in species] or [1.0]
            # Dissolvable source-rock reservoir per unit cell area (the weatherable
            # rock mass, default 1e6). Dissolution debits it, so a small pool
            # exhausts and the tracer goes inert; set it large where weathering must
            # stay rate-limited for the whole run (e.g. a thick saprolite source).
            self.gwGeoSourcePool = [float(s.get("source_pool", 1.0e6)) for s in species] or [1.0e6]
            self.gwGeochemConserve = bool(geo.get("conserve", True))
            # ext 2: route the groundwater-exported (seepage/baseflow) solute
            # DOWN the surface drainage network to the shoreline (river dissolved
            # load), reusing the flow-accumulation matrix. Opt-in. Per-species.
            self.gwRiverLoad = bool(geo.get("river_load", False))
            # ext 2 refinement: per-species first-order IN-TRANSIT loss along the
            # river (0 = conservative). And marine coupling: accumulate the
            # delivered coastal flux into a per-species marine reservoir.
            self.gwGeoRiverDecay = [
                float(s.get("river_decay", 0.0)) for s in species
            ] or [0.0]
            self.gwMarineCoupling = bool(geo.get("marine_coupling", False))

        except KeyError:
            self.gwOn = False
            self.duriOn = False
            self.gwKsat = 1.0
            self.gwSpecificYield = 0.1
            self.gwAquiferBase = 50.0
            self._gwAquiferBaseMap = None
            self.gwBedrockDepth = 0.0
            self.gwMinSatThick = 1.0
            self.gwInfiltration = 0.3
            self._gwInfilMap = None
            self.gwConserveBaseflow = True
            self.gwLakeExchange = False
            self.gwSubglacial = 0.0
            self.gwFineInfilFactor = 1.0
            self.gwInfilSlopeRef = 0.0
            self.gwPicardIts = 3
            self.gwSeepagePasses = 4
            self.duriFormRate = 0.0
            self.duriMaxThick = 5.0
            self.duriFringeDepth = 3.0
            self.duriFringeWidth = 2.0
            self.duriDischargeGate = False
            self.duriSupplyExp = 1.0
            self.duriWeatherEa = 0.0
            self.duriArmorMax = 0.0
            self.duriArmorDiffusion = False
            self.duriBreakRate = 0.0
            self.duriDecayRate = 0.0
            self.duriWeatherMode = "proxy"
            self.duriWeatherCeq = 1.0
            self.duriWeatherDw = 1.0
            self.duriWeatherL = 20.0
            self.duriWeatherability = 1.0
            self.gwGeochemOn = False
            self.gwNspecies = 1
            self.gwGeoName = ["solute0"]
            self.gwGeoWeather = [1.0]
            self.gwGeoWeatherByClass = [None]
            self.gwWeatherFrom = None
            self._gwLithoMap = None
            self.gwGeoCsat = [1.0]
            self.gwGeoPrecip = [1.0]
            self.gwGeoVsolid = [1.0]
            self.gwGeochemConserve = True
            self.gwRiverLoad = False
            self.gwGeoRiverDecay = [0.0]
            self.gwMarineCoupling = False

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
