import os
import gc
import sys
import petsc4py
import numpy as np

from mpi4py import MPI
from time import process_time

from gospl.tools.constants import BEDROCK_SENTINEL

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import strataonesed

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()


class STRAMesh(object):
    """
    This class encapsulates all the functions related to underlying stratigraphic information. 
    
    .. note::
        Sediment compaction in stratigraphic layers geometry and properties change are also considered.

    """

    _LIT_TABLE = ("litNb", "litName", "litFine", "litPhi0", "litZ0", "litK",
                  "litD", "litBed")

    def __getattr__(self, name):
        # Only reached when the attribute is missing: an object built without
        # the parser (unit tests create a bare STRAMesh and set the dual
        # attributes by hand) gets the one-lithology table of dual lithology.
        if name in STRAMesh._LIT_TABLE or name in ("_litEroded", "_litDeposited",
                                                    "litShare"):
            self._defaultLithTable()
            return self.__dict__[name]
        raise AttributeError(name)

    def _defaultLithTable(self):
        """Lithology table from the dual-lithology attributes (lithology 0 =
        fine), or an empty one without dual lithology."""
        d = self.__dict__
        if "litNb" not in d:
            if d.get("stratLith", False):
                d["litNb"] = 1
                d["litName"] = ["fine"]
                d["litFine"] = [True]
                d["litPhi0"] = [d.get("phi0f", 0.63)]
                d["litZ0"] = [d.get("z0f", 1960.0)]
                d["litK"] = [d.get("fine_k_factor", 1.0)]
                d["litD"] = [d.get("fine_diff_factor", 1.0)]
                d["litBed"] = [1.0 - d.get("bedrock_coarse_frac", 0.5)]
            else:
                d["litNb"] = 0
                for k in STRAMesh._LIT_TABLE[1:]:
                    d[k] = []
        d.setdefault("_litEroded", np.zeros(d["litNb"]))
        d.setdefault("_litDeposited", np.zeros(d["litNb"]))
        d.setdefault("litShare", None)

    # Lithology 0 (the dual-lithology fine fraction) under its historical
    # names. Reads and in-place writes go to stratHL[0] / phiL[0].
    @property
    def stratHf(self):
        hl = self.__dict__.get("stratHL")
        return None if hl is None else hl[0]

    @stratHf.setter
    def stratHf(self, value):
        if value is None:
            return
        if self.__dict__.get("stratHL") is None:
            self.stratHL = [value]
            if self.__dict__.get("phiL") is None:
                self.phiL = [None]
        else:
            self.stratHL[0] = value

    @property
    def phiF(self):
        pl = self.__dict__.get("phiL")
        return None if pl is None else pl[0]

    @phiF.setter
    def phiF(self, value):
        if value is None:
            return
        if self.__dict__.get("phiL") is None:
            self.phiL = [value]
            if self.__dict__.get("stratHL") is None:
                self.stratHL = [None]
        else:
            self.phiL[0] = value

    def _lsum(self, arrays):
        """Sum of a list of arrays; the single array itself when there is one
        (so a one-lithology run does exactly the dual arithmetic)."""
        if len(arrays) == 1:
            return arrays[0]
        out = arrays[0] + arrays[1]
        for a in arrays[2:]:
            out = out + a
        return out

    def _fineGroup(self, arrays):
        """Sum over the fine-like lithologies of a per-lithology list."""
        return self._lsum([a for a, f in zip(arrays, self.litFine) if f])

    def _litSoleFine(self, l):
        """True when lithology ``l`` is the only fine-like one: its share of
        the fine group is exactly 1."""
        return self.litFine[l] and sum(self.litFine) == 1

    def _groupShares(self):
        """``litShare`` made consistent: the fine-like shares of a node sum to
        exactly 1 (equal split where they carry nothing), the coarse-like ones
        to at most 1 (coarse clastic takes the rest). ``None`` for dual
        lithology. Rank-local."""
        if self.litNb <= 1:
            return None
        sh = np.clip(self.litShare, 0.0, None)
        fine = np.array(self.litFine)
        fsum = sh[fine].sum(axis=0)
        nf = int(fine.sum())
        if nf > 0:
            sh[fine] = np.where(fsum > 0.0, sh[fine] / np.where(fsum > 0.0, fsum, 1.0),
                                1.0 / nf)
        if (~fine).any():
            csum = sh[~fine].sum(axis=0)
            sh[~fine] /= np.maximum(csum, 1.0)
        return sh

    def _litDepoFrac(self, l, ff, sh=None):
        """Fraction of a deposit that is lithology ``l``, given the deposit's
        fine-group fraction ``ff``: ``ff * share`` (fine-like) or
        ``(1 - ff) * share`` (coarse-like), ``share`` being the lithology's
        part of its group (``sh`` from ``_groupShares``). Exactly ``ff`` for
        the sole fine lithology (dual lithology)."""
        if self._litSoleFine(l):
            return ff
        if sh is None:
            sh = self._groupShares()
        g = ff if self.litFine[l] else (1.0 - ff)
        return g * sh[l]

    def __init__(self):
        """
        The initialisation of `STRAMesh` class related to stratigraphic informations.
        """

        self.stratH = None
        self.stratZ = None
        self.phiS = None
        # Lithology piles (DESIGN_CARBONATES.md 2.7). Allocated only when
        # self.stratLith is True. stratHL[l] = bulk thickness of lithology l
        # per layer, phiL[l] its porosity; coarse clastic is the remainder
        # (stratH - sum_l stratHL[l]) with porosity phiS. Lithology 0 is the
        # dual-lithology fine fraction, also reachable as stratHf / phiF.
        self.stratHL = None
        self.phiL = None
        # Per-lithology share of a deposit inside its size group (fine-like or
        # coarse-like), (litNb, lpoints); None for a lithology alone in the
        # fine group, whose share is exactly 1 (the dual case).
        self.litShare = None
        self._litEroded = None
        self._litDeposited = None
        # Dual-lithology fine mass-balance diagnostics (m^3, owned-node running
        # totals over the run; reduced across ranks by the conservation test).
        # _fineEroded accumulates fine solid removed from the pile (erodeStrat);
        # _fineDeposited accumulates fine deposited back (deposeStrat, both the
        # continental and marine calls). On a closed sphere they must match to
        # within the floor/transit budget — the fine-specific analogue of the
        # total-cumED mass-conservation check.
        self._fineEroded = 0.0
        self._fineDeposited = 0.0
        # Per-layer erodibility multiplier (lpoints, stratNb). Effective K
        # at the surface = self.K * stratK[<top non-empty layer>]. Fresh
        # deposits and the bedrock sentinel default to 1.0 (no scaling).
        # An optional `stratK` key in the npstrata file lets the user
        # impose a non-uniform multiplier on the initial layers.
        self.stratK = None
        # In-model provenance tracers (opt-in `provenance:`; see
        # DESIGN_PROVENANCE.md §6). stratP[node, layer, class] = thickness of
        # each source-rock class in each layer (Σ over classes == stratH);
        # source_class = per-vertex bedrock class. Allocated only when provOn.
        self.stratP = None
        self.source_class = None
        # Per-layer diagenetic induration degree ∈ [0,1] (the archived duricrust
        # armoring record; DESIGN_WATERTABLE_DURICRUST.md §9). An INTENSIVE
        # property (like phiS): advects with the pile, compaction-neutral, and —
        # unlike stratK — 0 is a legitimate value (uncemented), so it is NEVER
        # forward-filled. Allocated only when the water table is on (its
        # producer) and stratigraphy is recorded (stratNb>0).
        self.stratDuri = None

        # Per-layer duricrust chemistry archive (Level-B geochemistry). Compact
        # categorical records — the DOMINANT solute species (`stratCrustType`)
        # and the DOMINANT source-rock region (`stratCrustSource`) of the crust
        # written into each layer — stored as float64 integer codes (-1 = no
        # crust). Frozen on burial, nearest-neighbour advected, compaction-
        # neutral. Allocated only with the geochemistry on (type) and in-model
        # provenance on (source). Consistent with `stratDuri` (the degree).
        self.stratCrustType = None
        self.stratCrustSource = None

        return

    def readStratLayers(self):
        """
        When stratigraphic layers are turned on, this function reads any initial stratigraphic layers provided in the input file (key: `npstrata`).

        The following variables will be read from the file:

        - thickness of each stratigrapic layer `strataH` accounting for both erosion & deposition events.
        - elevation at time of deposition, considered to be to the current elevation for the top stratigraphic layer `strataZ`.
        - porosity of coarse sediment `phiS` in each stratigraphic layer computed at center of each layer.

        With dual lithology enabled (`strata: dual: True`), two optional keys
        let each initial layer carry its own coarse/fine composition:

        - `strataHf`: fine-fraction bulk thickness per layer (coarse =
          `strataH - strataHf`); absent -> all-coarse.
        - `phiF`: fine porosity per layer; absent -> defaults to `phi0f`.
        """

        if self.strataFile is not None:
            fileData = np.load(self.strataFile)
            # Fail fast on a malformed file (missing required field / wrong
            # shape) and warn about an under-specified dual-lithology setup.
            nlay = self._checkStrataFile(fileData)
            has_hf = "strataHf" in fileData.files
            has_phiF = "phiF" in fileData.files
            has_K = "stratK" in fileData.files
            stratVal = fileData["strataH"]
            # Optionally reserve layer 0 as a dedicated infinite-bedrock
            # sentinel BELOW the file-provided layers (strata.bedrock_sentinel).
            # The file's `nlay` layers then occupy indices [off : off+nlay] and
            # layer 0 is the frozen 1e6 m reservoir. off=0 is the legacy
            # behaviour: the deepest file layer is itself the erosion floor.
            off = 1 if getattr(self, "bedrock_sentinel", False) else 0
            self.initLay = nlay + off
            self.stratNb += self.initLay
            lo, hi = off, off + nlay

            # Create stratigraphic arrays
            self.stratH = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratH[:, lo:hi] = stratVal[self.locIDs, 0:nlay]

            stratVal = fileData["strataZ"]
            self.stratZ = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratZ[:, lo:hi] = stratVal[self.locIDs, 0:nlay]

            stratVal = fileData["phiS"]
            self.phiS = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.phiS[:, lo:hi] = stratVal[self.locIDs, 0:nlay]

            # Per-layer K multiplier. Loaded from the optional `stratK`
            # key in the npstrata file (same (mpoints, nlay) shape as
            # the other layer fields); defaults to 1.0 (use self.K as-is).
            self.stratK = np.ones((self.lpoints, self.stratNb), dtype=np.float64)
            if "stratK" in fileData.files:
                stratVal = fileData["stratK"]
                self.stratK[:, lo:hi] = stratVal[self.locIDs, 0:nlay]

            # Dual-lithology fine-fraction layer fields. Optional in the
            # npstrata file: `strataHf` (per-layer fine bulk thickness, so each
            # initial layer carries its own coarse/fine composition; coarse =
            # strataH - strataHf) and `phiF` (per-layer fine porosity). Absent
            # `strataHf` -> all-coarse, so a single-fraction file still loads.
            if self.stratLith:
                self._allocLithologies()
                for l in range(self.litNb):
                    hkey, pkey = self._litFileKeys(l)
                    if hkey in fileData.files:
                        stratVal = fileData[hkey]
                        self.stratHL[l][:, lo:hi] = stratVal[self.locIDs, 0:nlay]
                        # Keep the partition physical: 0 <= pile <= total thickness.
                        np.clip(
                            self.stratHL[l][:, lo:hi],
                            0.0,
                            self.stratH[:, lo:hi],
                            out=self.stratHL[l][:, lo:hi],
                        )
                    if pkey in fileData.files:
                        stratVal = fileData[pkey]
                        self.phiL[l][:, lo:hi] = stratVal[self.locIDs, 0:nlay]
                    else:
                        # No porosity supplied: default to the lithology's
                        # surface porosity so initial layers holding it compact
                        # sensibly (irrelevant where its pile is empty).
                        self.phiL[l][:, lo:hi] = self.litPhi0[l]
                if self.litNb > 1:
                    # The piles together cannot exceed the layer thickness.
                    self._capLithologies(slice(lo, hi))
            elif "strataHf" in fileData.files and MPIrank == 0:
                print(
                    "Warning: the npstrata file provides 'strataHf' (per-layer "
                    "fine composition) but dual lithology is disabled. Set "
                    "`strata: dual: True` to use it; the field is ignored.",
                    flush=True,
                )

            if off:
                # Dedicated infinite-bedrock sentinel beneath the file layers
                # (mirrors the no-file path): a 1e6 m reservoir at layer 0,
                # split into coarse/fine by bedrock_coarse_frac under dual
                # lithology. bedrockLay > 0 freezes it in compaction
                # (_depthPorosity), and erodeStrat's transient inflate/deflate
                # of layer 0 makes it an un-erodable floor with this defined
                # composition. stratZ/stratK keep their 0.0/1.0 init (frozen,
                # so the values are inert). Provenance is seeded for every
                # layer (incl. this one) by _initProvenance below.
                self.stratH[:, 0] = BEDROCK_SENTINEL
                self.phiS[:, 0] = self.phi0s
                if self.stratLith:
                    for l in range(self.litNb):
                        self.stratHL[l][:, 0] = BEDROCK_SENTINEL * self.litBed[l]
                        self.phiL[l][:, 0] = self.litPhi0[l]
                self.bedrockLay = 1
            else:
                # All layers in the file are real sediment; no bedrock sentinel.
                self.bedrockLay = 0

            self._logStratInit(
                nfile_lay=nlay, has_hf=has_hf, has_phiF=has_phiF, has_K=has_K
            )

            if self.memclear:
                del fileData, stratVal
                gc.collect()
        else:
            self.stratH = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.phiS = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratZ = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratK = np.ones((self.lpoints, self.stratNb), dtype=np.float64)
            # Treat layer 0 as an effectively infinite bedrock reservoir when
            # no initial stratigraphy file is provided. The 1e6 m sentinel
            # cancels out in erodeStrat (cumThick / eroVal share the offset)
            # so it does not contaminate eroded volumes.
            self.stratH[:, 0] = BEDROCK_SENTINEL
            self.phiS[:, 0] = self.phi0s
            self.bedrockLay = 1          # layer 0 is the infinite-bedrock sentinel

            # Dual-lithology fine-fraction fields. The infinite-bedrock
            # sentinel layer is split into coarse/fine by bedrock_coarse_frac
            # so material eroded from bedrock inherits a defined composition.
            if self.stratLith:
                self._allocLithologies()
                for l in range(self.litNb):
                    self.stratHL[l][:, 0] = BEDROCK_SENTINEL * self.litBed[l]
                    self.phiL[l][:, 0] = self.litPhi0[l]

            self._logStratInit()

        if getattr(self, "provOn", False):
            self._initProvenance()

        # Diagenetic induration archive (§9) — allocated uncemented (0) only when
        # the groundwater/duricrust feature is on. Initial layers start with no
        # crust; the record is grown by _recordInduration and restored on restart.
        if getattr(self, "gwOn", False):
            self.stratDuri = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)

        # Per-layer crust chemistry archive (Level-B geochemistry). Integer
        # codes stored as float64, -1 = no crust in that layer. `stratCrustType`
        # (dominant solute species) is grown whenever geochemistry is on;
        # `stratCrustSource` (dominant source region) additionally needs in-model
        # provenance. Written by _recordInduration, restored on restart.
        if getattr(self, "gwGeochemOn", False):
            self.stratCrustType = np.full(
                (self.lpoints, self.stratNb), -1.0, dtype=np.float64
            )
            if getattr(self, "provOn", False):
                self.stratCrustSource = np.full(
                    (self.lpoints, self.stratNb), -1.0, dtype=np.float64
                )

        return

    def _allocLithologies(self):
        """Allocate the lithology piles and their bookkeeping (zeroed)."""
        shape = (self.lpoints, self.stratNb)
        self.stratHL = [np.zeros(shape, dtype=np.float64) for _ in range(self.litNb)]
        self.phiL = [np.zeros(shape, dtype=np.float64) for _ in range(self.litNb)]
        self._litEroded = np.zeros(self.litNb)
        self._litDeposited = np.zeros(self.litNb)
        if self.litNb > 1:
            self.litShare = np.zeros((self.litNb, self.lpoints), dtype=np.float64)

    def _litFileKeys(self, l):
        """npstrata / output dataset names of lithology ``l``'s pile and
        porosity: the historical ``strataHf`` / ``phiF`` for lithology 0,
        ``strataH_<name>`` / ``phi_<name>`` for the others."""
        if l == 0:
            return "strataHf", "phiF"
        return "strataH_" + self.litName[l], "phi_" + self.litName[l]

    def _litOutKeys(self, l):
        """Output (and restart) dataset names of lithology ``l``: ``stratHf`` /
        ``phiF`` for lithology 0, ``stratH_<name>`` / ``phi_<name>`` for the
        others."""
        if l == 0:
            return "stratHf", "phiF"
        return "stratH_" + self.litName[l], "phi_" + self.litName[l]

    def _capLithologies(self, sl):
        """Scale the lithology piles down where, together, they exceed the
        layer thickness (columns ``sl``). Only used with several lithologies;
        dual lithology keeps its own ``np.minimum``."""
        tot = self._lsum([h[:, sl] for h in self.stratHL])
        H = self.stratH[:, sl]
        over = tot > H
        if not over.any():
            return
        scale = np.ones_like(tot)
        scale[over] = np.divide(H[over], tot[over], out=np.zeros_like(H[over]),
                                where=tot[over] > 0)
        for h in self.stratHL:
            h[:, sl] *= scale

    def _checkStrataFile(self, fileData):
        """
        Validate an ``npstrata`` initial-stratigraphy file before it is loaded.

        Complains (raises ``ValueError``) when a required field is missing or
        an array shape is inconsistent, so a malformed file fails fast with a
        clear message instead of a cryptic ``KeyError`` / broadcast error deep
        in the loader. Under dual lithology it also warns (rank 0) when the
        per-layer coarse/fine composition (``strataHf``) is absent, since the
        initial pile then silently defaults to all-coarse.

        :arg fileData: the opened ``np.load`` archive.
        :return: the number of initial layers ``nlay`` (``strataH`` columns).
        """
        avail = set(fileData.files)

        # 1) Fields required for ANY npstrata file.
        required = ("strataH", "strataZ", "phiS")
        missing = [k for k in required if k not in avail]
        if missing:
            raise ValueError(
                "npstrata file '%s' is missing required field(s): %s. An "
                "initial stratigraphy file must provide 'strataH' (per-layer "
                "thickness), 'strataZ' (deposition elevation) and 'phiS' "
                "(coarse porosity), each shaped (mesh_points, n_layers)."
                % (self.strataFile, ", ".join(missing))
            )

        # 2) Shape consistency: strataH is (mesh_points, n_layers); every other
        # layer field must match it (arrays are indexed by global vertex id).
        ref = fileData["strataH"].shape
        if len(ref) != 2 or ref[0] != self.mpoints:
            raise ValueError(
                "npstrata file '%s': 'strataH' has shape %s but must be "
                "(mesh_points=%d, n_layers); each row is one mesh vertex."
                % (self.strataFile, ref, self.mpoints)
            )
        litKeys = []
        for l in range(getattr(self, "litNb", 0)):
            litKeys += list(self._litFileKeys(l))
        for k in ("strataZ", "phiS", "stratK", "strataHf", "phiF", *litKeys):
            if k in avail and fileData[k].shape != ref:
                raise ValueError(
                    "npstrata file '%s': field '%s' has shape %s but must "
                    "match 'strataH' %s."
                    % (self.strataFile, k, fileData[k].shape, ref)
                )

        # 3) Dual-lithology composition: warn if the fine split is unspecified.
        if self.stratLith and "strataHf" not in avail and MPIrank == 0:
            print(
                "Warning: dual lithology is enabled but the npstrata file '%s' "
                "provides no 'strataHf' (per-layer fine bulk thickness); the "
                "initial layers will be treated as ALL-COARSE. Add a 'strataHf' "
                "array (and optionally 'phiF'), shaped like 'strataH', to set a "
                "per-vertex initial coarse/fine distribution." % self.strataFile,
                flush=True,
            )

        return ref[1]

    def _logStratInit(self, nfile_lay=0, has_hf=False, has_phiF=False, has_K=False):
        """
        Rank-0 verbose summary of the stratigraphy setup (run with ``-v``).

        Reports the single- vs dual-lithology mode and, under dual lithology,
        the per-fraction compaction curves and contrast knobs, the source of
        the initial layers (file vs bedrock-only), whether a per-vertex
        coarse/fine distribution was supplied, and the bedrock-floor model.
        """
        if MPIrank != 0 or not getattr(self, "verbose", False):
            return

        if self.stratLith:
            print(
                "Dual-lithology stratigraphy ON — coarse phi0/z0 = %.3g/%.0f m, "
                "fine phi0/z0 = %.3g/%.0f m, fine K x%.3g, fine diffusivity x%.3g"
                % (
                    self.phi0c, self.z0c, self.phi0f, self.z0f,
                    getattr(self, "fine_k_factor", 1.0),
                    getattr(self, "fine_diff_factor", 1.0),
                ),
                flush=True,
            )
        else:
            print("Single-lithology stratigraphy ON", flush=True)

        if self.strataFile is not None:
            if not self.stratLith:
                comp = "single-fraction"
            elif has_hf:
                comp = "per-vertex coarse/fine from 'strataHf'" + (
                    "" if has_phiF else " (phiF defaulted to phi0f)"
                )
            else:
                comp = "all-coarse (no 'strataHf' in file)"
            print(
                "   initial layers from '%s': %d layer(s); composition: %s%s"
                % (
                    self.strataFile, nfile_lay, comp,
                    "; per-layer K from 'stratK'" if has_K else "",
                ),
                flush=True,
            )
            if self.bedrockLay > 0:
                print(
                    "   + dedicated infinite-bedrock sentinel beneath them"
                    + (
                        " (bedrock_coarse_frac = %.3g)" % self.bedrock_coarse_frac
                        if self.stratLith else ""
                    ),
                    flush=True,
                )
            else:
                print(
                    "   deepest file layer is the infinite-bedrock floor "
                    "(set strata.bedrock_sentinel: True for a dedicated one)",
                    flush=True,
                )
        else:
            print(
                "   infinite-bedrock sentinel only (no npstrata file)"
                + (
                    "; bedrock_coarse_frac = %.3g" % self.bedrock_coarse_frac
                    if self.stratLith else ""
                ),
                flush=True,
            )
        return

    def logDualState(self):
        """
        Rank-0 per-step verbose line summarising the dual-lithology sediment
        state (run with ``-v``): the area-weighted mean surface fine fraction
        and the fine share of the whole recorded pile (excluding the bedrock
        sentinel). Called from the time loop so dual lithology stays visible
        throughout the run, not only in the one-off init summary
        (:meth:`_logStratInit`). No-op unless dual lithology is on and verbose.

        The reductions are COLLECTIVE, so the method runs on every rank (gated
        only by the uniform ``stratLith``/``verbose``); the line prints on rank
        0.
        """
        if not self.stratLith or not getattr(self, "verbose", False):
            return

        owned = self.inIDs == 1
        a = self.larea
        ff = 1.0 - self._surfaceComposition()             # surface fine fraction
        comm = MPI.COMM_WORLD
        surfNum = comm.allreduce(float(np.sum((ff * a)[owned])), op=MPI.SUM)
        area = comm.allreduce(float(np.sum(a[owned])), op=MPI.SUM)
        # Pile fine share over the recorded (non-sentinel) layers.
        b = self.bedrockLay
        H = self.stratH[:, b : self.stratStep + 1].sum(axis=1)
        Hf = self._fineGroup(
            [h[:, b : self.stratStep + 1] for h in self.stratHL]
        ).sum(axis=1)
        vH = comm.allreduce(float(np.sum((H * a)[owned])), op=MPI.SUM)
        vHf = comm.allreduce(float(np.sum((Hf * a)[owned])), op=MPI.SUM)
        if MPIrank == 0:
            surf = surfNum / area if area > 0.0 else 0.0
            pile = vHf / vH if vH > 0.0 else 0.0
            print(
                "  [dual] surface fine fraction %.3f | recorded-pile fine "
                "share %.3f" % (surf, pile),
                flush=True,
            )
        return

    def _logProvInit(self):
        """
        Rank-0 ``-v`` setup summary for the in-model provenance tracers: the
        number of source classes and where the per-vertex bedrock class comes
        from (a uniform constant or a ``[file, key]`` map). Descriptive only
        (no reductions); the evolving per-class shares are reported each step
        by :meth:`logProvState`.
        """
        if MPIrank != 0 or not getattr(self, "verbose", False):
            return
        if getattr(self, "_provSourceMap", None) is not None:
            src = "per-vertex map '%s':'%s'" % (
                self._provSourceMap[0], self._provSourceMap[1]
            )
        else:
            src = "uniform class %d (single source)" % int(self._provSourceUniform)
        print(
            "Provenance tracers ON — %d source class(es); bedrock source: %s"
            % (self.provNb, src),
            flush=True,
        )
        return

    def logProvState(self):
        """
        Rank-0 per-step verbose line: the source-class composition of the
        recorded stratigraphic pile (area-weighted volume share per class,
        excluding the bedrock sentinel) — i.e. where the deposited sediment has
        come from so far. No-op unless provenance is on and verbose.

        The reductions are COLLECTIVE (run on every rank, gated only by the
        uniform ``provOn``/``verbose``); the line prints on rank 0.
        """
        if not getattr(self, "provOn", False) or not getattr(self, "verbose", False):
            return

        owned = self.inIDs == 1
        a = self.larea
        b = self.bedrockLay
        # Per-class and total recorded-pile volume (area-weighted, owned nodes).
        P = self.stratP[:, b : self.stratStep + 1, :].sum(axis=1)   # (lpoints, prov)
        H = self.stratH[:, b : self.stratStep + 1].sum(axis=1)      # (lpoints,)
        comm = MPI.COMM_WORLD
        volP = np.ascontiguousarray(
            (P * a[:, None])[owned].sum(axis=0), dtype=np.float64
        )
        comm.Allreduce(MPI.IN_PLACE, volP, op=MPI.SUM)
        vH = comm.allreduce(float(np.sum((H * a)[owned])), op=MPI.SUM)
        if MPIrank == 0:
            if vH > 0.0:
                shares = " ".join(
                    "c%d %.3f" % (c, volP[c] / vH) for c in range(self.provNb)
                )
            else:
                shares = "(empty pile)"
            print("  [prov] recorded-pile source shares: %s" % shares, flush=True)
        return

    def _initProvenance(self):
        """
        Allocate and seed the provenance state (opt-in `provenance:`; Phase 0).

        ``source_class`` (per-vertex bedrock class) is read from the ``uniform``
        scalar or the ``source`` ``[file, key]`` map; ``stratP[node, layer,
        class]`` is the per-class thickness in each layer, seeded so every
        initial layer (and the bedrock sentinel) carries the node's bedrock
        source class (Σ over classes == stratH). A passive tracer — no physics
        feedback — so a provenance-on run is identical to provenance-off until
        the erosion/transport/deposition hooks of later phases are added.
        """
        if self._provSourceMap is not None:
            fname, key = self._provSourceMap[0], self._provSourceMap[1]
            data = np.load(fname + ".npz")
            self.source_class = data[key][self.locIDs].astype(np.int64)
            del data
        else:
            self.source_class = np.full(
                self.lpoints, int(self._provSourceUniform), dtype=np.int64
            )
        np.clip(self.source_class, 0, self.provNb - 1, out=self.source_class)

        # Per-layer per-class thickness; each initial layer's full thickness is
        # assigned to the node's bedrock source class.
        self.stratP = np.zeros(
            (self.lpoints, self.stratNb, self.provNb), dtype=np.float64
        )
        self.stratP[np.arange(self.lpoints), :, self.source_class] = self.stratH

        self._logProvInit()
        return

    def _fillZeroPorosity(self, phiS):
        """
        Where ``phiS == 0`` (layer has been emptied or never deposited),
        inherit the value from the nearest underlying layer with non-zero
        porosity. Vectorised forward-fill from low → high index along axis 1.
        Leading zeros at the column base stay zero (no valid layer below).
        """
        mask = phiS > 0
        idx = np.where(mask, np.arange(phiS.shape[1])[None, :], 0)
        np.maximum.accumulate(idx, axis=1, out=idx)
        return np.take_along_axis(phiS, idx, axis=1)

    def _surfaceK(self):
        """
        Return the per-node erodibility multiplier of the topmost
        non-empty stratigraphic layer. Used by the SPL flavours to scale
        the scalar ``self.K`` according to the bedrock currently exposed
        at the surface.

        Returns 1.0 everywhere when stratigraphy is disabled
        (``stratNb == 0``) or when the column is fully empty at a node,
        so the SPL evaluation falls back to the YAML-default K.
        """
        if self.stratNb == 0 or self.stratK is None:
            return np.ones(self.lpoints, dtype=np.float64)

        sK = self.stratK[:, : self.stratStep + 1]
        mask = sK > 0
        # For each row, locate the highest column index where the layer
        # has a non-zero multiplier. Reverse-then-argmax picks the first
        # non-zero starting from the top of the column.
        rev = mask[:, ::-1]
        any_valid = rev.any(axis=1)
        top_from_right = np.argmax(rev, axis=1)
        top_idx = sK.shape[1] - 1 - top_from_right

        out = np.ones(self.lpoints, dtype=np.float64)
        rows = np.arange(sK.shape[0])
        out[any_valid] = sK[rows[any_valid], top_idx[any_valid]]
        return out

    def _surfaceComposition(self):
        """
        Return the per-node coarse fraction ``fc`` (by bulk thickness) of
        the topmost non-empty stratigraphic layer; the fine fraction is
        ``1 - fc``.

        Returns all-ones (all coarse) when dual lithology is disabled or the
        fine pile is unallocated, so single-fraction callers are unchanged.
        Shared hook for the SPL erodibility blend (`_surfaceLithoK`) and the
        hillslope diffusivity blend (see DESIGN_DUAL_LITHOLOGY.md).
        """
        if not self.stratLith or self.stratHf is None or self.stratNb == 0:
            return np.ones(self.lpoints, dtype=np.float64)

        H = self.stratH[:, : self.stratStep + 1]
        Hf = self._fineGroup([h[:, : self.stratStep + 1] for h in self.stratHL])
        mask = H > 0
        # Topmost non-empty layer per node (reverse-then-argmax, as _surfaceK).
        rev = mask[:, ::-1]
        any_valid = rev.any(axis=1)
        top_from_right = np.argmax(rev, axis=1)
        top_idx = H.shape[1] - 1 - top_from_right

        fc = np.ones(self.lpoints, dtype=np.float64)
        rows = np.arange(H.shape[0])
        Htop = H[rows, top_idx]
        Hftop = Hf[rows, top_idx]
        valid = any_valid & (Htop > 0)
        fc[valid] = 1.0 - (Hftop[valid] / Htop[valid])
        np.clip(fc, 0.0, 1.0, out=fc)
        return fc

    def _surfaceLithFractions(self):
        """
        Per-lithology bulk fraction of the topmost non-empty layer, a list of
        ``litNb`` arrays (zeros where a node has no layer). The coarse clastic
        fraction is ``1 - sum``. Used by the erodibility and diffusivity
        blends when there are several lithologies.
        """
        top = self.stratStep + 1
        H = self.stratH[:, :top]
        mask = H > 0
        rev = mask[:, ::-1]
        any_valid = rev.any(axis=1)
        top_idx = H.shape[1] - 1 - np.argmax(rev, axis=1)
        rows = np.arange(H.shape[0])
        Htop = H[rows, top_idx]
        valid = any_valid & (Htop > 0)
        out = []
        for h in self.stratHL:
            f = np.zeros(self.lpoints, dtype=np.float64)
            f[valid] = h[rows, top_idx][valid] / Htop[valid]
            np.clip(f, 0.0, 1.0, out=f)
            out.append(f)
        return out

    def _surfaceLithoBlend(self, factors):
        """``f_clastic + sum_l f_l * factors[l]`` on the surface layer, the
        clastic fraction being the remainder (several lithologies)."""
        fl = self._surfaceLithFractions()
        fsum = self._lsum(fl)
        out = np.clip(1.0 - fsum, 0.0, 1.0)
        for f, k in zip(fl, factors):
            out = out + f * k
        return out

    def _surfaceArmoringK(self):
        """
        Duricrust erodibility-armoring multiplier (≤ 1) — the **single hook**
        (DESIGN_WATERTABLE_DURICRUST.md §5) through which the generic duricrust
        modulates ALL three eroders (SPL / nlSPL / soilSPL) with **no branching**
        in them: it composes multiplicatively inside ``_surfaceLithoK``.

        Returns the scalar ``1.0`` (a byte-identical no-op) when the duricrust is
        off. When on, ``1 − armor_max · duriF``: a fully indurated cell
        (``duriF = 1``) cuts erodibility by ``armor_max`` (e.g. 0.9 ⇒ 10× more
        resistant, producing relief inversion); ``duriF = 0`` leaves ``K``
        unchanged. Rate-only — no elevation change, so flow/routing are untouched
        this step (like dual-lithology deposition being composition-only).
        """
        if not getattr(self, "duriOn", False):
            return 1.0
        return 1.0 - self.duriArmorMax * self.duriF

    def _surfaceLithoK(self):
        """
        Per-node erodibility multiplier from the exposed surface
        composition — ``(fc + (1 - fc) * fine_k_factor)`` — times the duricrust
        armoring factor (``_surfaceArmoringK``).

        Equals 1.0 everywhere when dual lithology **and** the duricrust are off
        (or ``fine_k_factor == 1`` with no crust), so it composes multiplicatively
        with ``_surfaceK`` in the SPL flavours without altering single-fraction
        behaviour.
        """
        if self.stratLith and self.litNb > 1:
            litho = self._surfaceLithoBlend(self.litK)
        elif self.stratLith:
            fc = self._surfaceComposition()
            litho = fc + (1.0 - fc) * self.fine_k_factor
        else:
            litho = np.ones(self.lpoints, dtype=np.float64)
        return litho * self._surfaceArmoringK()

    def _surfaceLithoD(self):
        """
        Per-node diffusivity multiplier from the exposed surface composition:
        ``fc + (1 - fc) * fine_diff_factor``.

        Equals 1.0 everywhere when dual lithology is off (or
        ``fine_diff_factor == 1``, i.e. no contrast), so it composes
        multiplicatively with the base hillslope coefficients (``Cda``/``Cdm``)
        without altering single-fraction behaviour. Fines diffuse faster when
        ``fine_diff_factor > 1`` (see DESIGN_DUAL_LITHOLOGY.md Section 7).
        A cemented duricrust also resists creep when ``armor_diffusion`` is on
        (off by default), scaling ``Cd`` by the same ``1 − armor_max · duriF``.
        """
        if self.stratLith and self.litNb > 1:
            litho = self._surfaceLithoBlend(self.litD)
        elif self.stratLith:
            fc = self._surfaceComposition()
            litho = fc + (1.0 - fc) * self.fine_diff_factor
        else:
            litho = np.ones(self.lpoints, dtype=np.float64)
        if getattr(self, "duriOn", False) and self.duriArmorDiffusion:
            litho = litho * (1.0 - self.duriArmorMax * self.duriF)
        return litho

    def deposeStrat(self):
        """
        Add deposition on top of an existing stratigraphic layer. The following variables will be recorded:

        - thickness of each stratigrapic layer `stratH` accounting for both erosion & deposition events.
        - porosity of sediment `phiS` in each stratigraphic layer computed at center of each layer.

        In dual-lithology mode the deposit is split into a coarse and a fine
        fraction using the **per-node deposit fine fraction** ``self.depoFineFrac``.
        It starts each step as ``self.fineFrac`` (the composition of the routed
        sediment arriving at each node, from ``sedplex._getSedFlux``) and is
        refined inside continental depressions by ``_pitFineFraction`` (3b:
        fine biased to the depocenter, coarse to the inlet/margins).
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        depo = self.tmpL.getArray().copy()
        depo[depo < 1.0e-4] = 0.0
        self.stratH[:, self.stratStep] += depo
        ids = np.where(depo > 0)[0]
        if self.stratLith:
            # Each lithology's part of the deposit: the fine-group fraction
            # depoFineFrac (or its complement for a coarse-like lithology)
            # times its share of the group (exactly depoFineFrac for the
            # dual-lithology fine fraction).
            own = self.inIDs == 1
            sh = self._groupShares()
            for l in range(self.litNb):
                litDepo = depo * self._litDepoFrac(l, self.depoFineFrac, sh)
                self.stratHL[l][:, self.stratStep] += litDepo
                vol = float(np.sum((litDepo * self.larea)[own]))
                self._litDeposited[l] += vol
                if self.litFine[l]:
                    # Mass-balance diagnostic of the fine group (m^3).
                    self._fineDeposited += vol
                self.phiL[l][ids, self.stratStep] = self.litPhi0[l]
            # Fresh deposit carries each lithology's surface porosity.
            self.phiS[ids, self.stratStep] = self.phi0c
        else:
            self.phiS[ids, self.stratStep] = self.phi0s
        if getattr(self, "provOn", False):
            # Lay the deposit into the layer's provenance composition (the
            # arriving routed composition; a passive tracer, so no sorting).
            # Keeps Σ over classes == stratH for the deposited thickness.
            frac = self.depoProvFrac.copy()
            fsum = frac.sum(axis=1)
            # Depositing nodes with no arriving composition (Σ_c ≈ 0) — e.g.
            # off-channel hillslope-creep deposits, which have no river
            # through-flux so `provFrac` is zero there. Without a fallback the
            # layer gets thickness but zero provenance (Σ_c stratP < stratH), so
            # the post-processed cell shows no source class (`dominant == -1`).
            # Attribute such locally-derived deposits to the in-situ bedrock
            # `source_class` (the best label available without threading
            # provenance through hillslope diffusion). Then renormalise every
            # depositing node so the layer is exactly partitioned (Σ_c == depo);
            # routed/pit/marine fractions already sum to 1, so this is a no-op
            # for them and only repairs the hillslope holes.
            holes = (depo > 0.0) & (fsum < 1.0e-6)
            if holes.any():
                frac[holes, :] = 0.0
                frac[holes, self.source_class[holes]] = 1.0
                fsum = frac.sum(axis=1)
            good = fsum > 0.0
            frac[good] /= fsum[good, None]
            # The composition actually laid down, for _netRoutedSource.
            self._lastDepoProvFrac = frac
            provDepo = depo[:, None] * frac
            self.stratP[:, self.stratStep, :] += provDepo
            self._provDeposited += np.sum(
                (provDepo * self.larea[:, None])[self.inIDs == 1], axis=0
            )
        # Freshly deposited sediment erodibility multiplier (top layer).
        # - lumped mode: 1.0 (no multiplier) — the deposit becomes soil and gets
        #   its soft erodibility (Ksoil) via updateSoilThickness.
        # - regolith mode: the deposit stays in the stratigraphy, so it must
        #   carry its own SOFT erodibility here. Set stratK = Ksoil/K so the SPL
        #   bedrock term Kbr*stratK = Ksoil — fresh unconsolidated sediment erodes
        #   like soil, reusing the already-defined soilK (regolith mode implies
        #   cptSoil, so Ksoil is defined). See DESIGN_SOIL_REGOLITH.md §5.
        if getattr(self, "regolithSoil", False) and self.K > 0.0:
            self.stratK[ids, self.stratStep] = self.Ksoil / self.K
        else:
            self.stratK[ids, self.stratStep] = 1.0

        # Cleaning arrays
        if self.memclear:
            del depo, ids
            gc.collect()

        return

    def _netRoutedSource(self):
        """
        Turn the eroded-sediment source produced by ``erodeStrat`` into the
        NET source that must be routed downstream, by subtracting what the
        same eroder just deposited in place. Call right after the eroder's
        ``erodeStrat`` + ``deposeStrat`` pair (SPL, nlSPL, soilSPL).

        Why: ``erodeStrat`` turns only the INCISING part of the step's signed
        thickness change (``self.tmp`` = ``Eb*dt``) into ``thCoarse`` /
        ``thFine`` / ``provEro``, while ``deposeStrat`` lays the DEPOSITING
        part (transport-limited SPL, ``spl: G > 0``) into the pile.
        ``sedplex._getSedFlux`` then routed the gross erosion, so the in-place
        deposit was sent downstream as well and laid down a second time, in a
        pit or the sea. Without stratigraphy ``_getSedFlux`` routes ``-Eb`` (the
        net rate) and was always right; with it, ``stratigraphic_record``'s
        ``input-stratiG`` (``G = 3``) put 38% more sediment into the sea than
        it eroded. The subtraction uses the full ``max(Eb*dt, 0)`` (what changed
        the elevation, including deposits below ``deposeStrat``'s 0.1 mm
        stratigraphic floor), split by the composition ``deposeStrat`` used:
        ``depoFineFrac`` per fraction, ``_lastDepoProvFrac`` per class (the
        bedrock ``source_class`` where none was recorded). The sources become
        negative at depositing nodes, exactly as ``-Eb`` is.

        Rank-local (no collective). Scratch Vecs: reads ``self.tmp`` through
        ``self.tmpL``.
        """
        self.dm.globalToLocal(self.tmp, self.tmpL)
        dep = np.maximum(self.tmpL.getArray(), 0.0) / self.dt
        if not np.any(dep > 0.0):
            return
        if self.stratLith:
            ff = np.clip(self.depoFineFrac, 0.0, 1.0)
            self.thFine = self.thFine - dep * ff
            self.thCoarse = self.thCoarse - dep * (1.0 - ff)
            if self.litNb > 1:
                sh = self._groupShares()
                for l in range(self.litNb):
                    self.thL[l] = self.thL[l] - dep * self._litDepoFrac(l, ff, sh)
        else:
            self.thCoarse = self.thCoarse - dep
        if getattr(self, "provOn", False):
            pf = getattr(self, "_lastDepoProvFrac", None)
            if pf is None or pf.shape != (self.lpoints, self.provNb):
                pf = np.zeros((self.lpoints, self.provNb))
            pf = pf.copy()
            none = (dep > 0.0) & (pf.sum(axis=1) < 1.0e-6)
            if none.any():
                pf[none, :] = 0.0
                pf[none, self.source_class[none]] = 1.0
            self.provEro = self.provEro - dep[:, None] * pf

        return

    def erodeStrat(self):
        """
        This function removes eroded sediment thicknesses from the stratigraphic pile. The function takes into account the porosity values of considered lithologies in each eroded stratigraphic layers.

        It follows the assumptions:

        - Eroded thicknesses from stream power law and hillslope diffusion are considered to encompass both the solid and void phase.
        - Only the solid phase will be moved dowstream by surface processes.
        - The corresponding deposit thicknesses for those freshly eroded sediments correspond to uncompacted thicknesses based on the porosity at surface given from the input file.
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        ero = self.tmpL.getArray().copy()
        ero[ero > 0] = 0.0

        # Nodes experiencing erosion
        nids = np.where(ero < 0)[0]
        if len(nids) == 0:
            self.thCoarse = np.zeros(self.lpoints)
            if self.stratLith:
                self.thFine = np.zeros(self.lpoints)
                self.thL = [np.zeros(self.lpoints) for _ in range(self.litNb)]
            if getattr(self, "provOn", False):
                self.provEro = np.zeros((self.lpoints, self.provNb))
            return

        # Cumulative thickness for each node
        if self.stratLith:
            # Inflate the fine pile in proportion to layer 0's CURRENT fine
            # fraction so the sentinel never alters that layer's composition
            # (for the infinite-bedrock case this fraction is 1-bedrock_coarse_frac;
            # for a real basal layer it is its own fine ratio, i.e. 0 if all coarse).
            # Same for every lithology pile.
            H0 = self.stratH[nids, 0]
            sentL = [
                BEDROCK_SENTINEL * np.divide(
                    h[nids, 0], H0, out=np.zeros_like(H0), where=H0 > 0
                )
                for h in self.stratHL
            ]
        self.stratH[nids, 0] += BEDROCK_SENTINEL
        if self.stratLith:
            for h, sl in zip(self.stratHL, sentL):
                h[nids, 0] += sl
        provOn = getattr(self, "provOn", False)
        if provOn:
            # The infinite-bedrock sentinel is pure bedrock -> the node's source
            # class (keeps Σ over classes == stratH for layer 0).
            self.stratP[nids, 0, self.source_class[nids]] += BEDROCK_SENTINEL
        cumThick = np.cumsum(self.stratH[nids, self.stratStep :: -1], axis=1)[:, ::-1]
        boolMask = cumThick < -ero[nids].reshape((len(nids), 1))
        mask = boolMask.astype(int)
        if provOn:
            # Per-class bulk eroded from the FULLY consumed layers (captured
            # before they are cleared below).
            provMaskedBulk = (
                self.stratP[nids, : self.stratStep + 1, :] * mask[:, :, None]
            ).sum(axis=1)

        thickS = self.stratH[nids, 0 : self.stratStep + 1]
        if self.stratLith:
            # Each layer's bulk thickness splits into coarse and fine; each
            # fraction yields its own solid phase (1 - its own porosity).
            # Coarse clastic is the remainder of the lithology piles.
            litTh = [h[nids, 0 : self.stratStep + 1] for h in self.stratHL]
            coarseTh = thickS - self._lsum(litTh)
            thCoarse = coarseTh * (1.0 - self.phiS[nids, 0 : self.stratStep + 1])
            thLit = [
                np.sum(
                    (th * (1.0 - p[nids, 0 : self.stratStep + 1])) * mask, axis=1
                )
                for th, p in zip(litTh, self.phiL)
            ]
            thCoarse = np.sum((thCoarse * mask), axis=1)
        else:
            thCoarse = thickS * (1.0 - self.phiS[nids, 0 : self.stratStep + 1])
            thCoarse = np.sum((thCoarse * mask), axis=1)

        # Clear all stratigraphy points which are eroded
        cumThick[boolMask] = 0.0
        tmp = self.stratH[nids, : self.stratStep + 1]
        tmp[boolMask] = 0
        self.stratH[nids, : self.stratStep + 1] = tmp
        if self.stratLith:
            for h in self.stratHL:
                tmpf = h[nids, : self.stratStep + 1]
                tmpf[boolMask] = 0
                h[nids, : self.stratStep + 1] = tmpf
        if provOn:
            tmpp = self.stratP[nids, : self.stratStep + 1, :]
            tmpp[boolMask] = 0.0
            self.stratP[nids, : self.stratStep + 1, :] = tmpp

        # Erode remaining stratal layers
        # Get non-zero top layer number. `minlength=len(nids)` is required: a
        # node whose entire column was consumed (every layer in `boolMask`, i.e.
        # the erosion demand exceeded the total pile *including* the
        # BEDROCK_SENTINEL infinite-bedrock floor) contributes no nonzero row to
        # `np.nonzero(cumThick)[0]`. Without minlength, `bincount` silently
        # truncates such an all-zero row when it is the LAST node, so `eroLayNb`
        # comes back one short and the fancy-index below raises a broadcast
        # IndexError (a middle all-zero row already gets a 0 count, so only
        # trailing ones crash). The sentinel makes a fully-consumed column
        # impossible for physical erosion, but a non-converged upstream solve
        # (e.g. glacial abrasion feeding `dz_ero`) can produce
        # a pathological demand that reaches here; padding keeps `erodeStrat`
        # crash-safe and treats those rows like any other fully-eroded node.
        eroLayNb = np.bincount(np.nonzero(cumThick)[0], minlength=len(nids)) - 1
        eroVal = cumThick[np.arange(len(nids)), eroLayNb] + ero[nids]

        self.thCoarse = np.zeros(self.lpoints)
        # From the partially eroded top layer extract the solid phase removed.
        H_old = self.stratH[nids, eroLayNb]
        tmp = H_old - eroVal
        tmp[tmp < 1.0e-8] = 0.0  # TODO-REFACTOR: value matches DISCHARGE_FLOOR but distinct role (thickness numerical-noise floor); do not replace
        if self.stratLith:
            Hl_old = [h[nids, eroLayNb] for h in self.stratHL]
            # Well-mixed layer: each lithology's bulk is removed in proportion
            # to the layer composition; the remainder keeps the same ratio.
            frac = np.divide(tmp, H_old, out=np.zeros_like(tmp), where=H_old > 0)
            coarse_rm = (H_old - self._lsum(Hl_old)) * frac
            thCoarse = thCoarse + coarse_rm * (1.0 - self.phiS[nids, eroLayNb])
            keep = np.divide(eroVal, H_old, out=np.zeros_like(eroVal), where=H_old > 0)
            # Uncompacted thickness deposited downstream, per lithology, using
            # each lithology's surface porosity.
            self.thCoarse[nids] = thCoarse / (1.0 - self.phi0c)
            self.thCoarse[self.thCoarse < 0.0] = 0.0
            self.thL = []
            for l in range(self.litNb):
                thl = thLit[l] + (Hl_old[l] * frac) * (1.0 - self.phiL[l][nids, eroLayNb])
                th = np.zeros(self.lpoints)
                th[nids] = thl / (1.0 - self.litPhi0[l])
                th[th < 0.0] = 0.0
                self.thL.append(th)
                # Remaining pile in the partial layer scales with the remaining
                # total.
                self.stratHL[l][nids, eroLayNb] = Hl_old[l] * keep
        else:
            # Define the uncompacted sand thickness that will be deposited dowstream
            thCoarse += tmp * (1.0 - self.phiS[nids, eroLayNb])
            self.thCoarse[nids] = thCoarse / (1.0 - self.phi0s)
            self.thCoarse[self.thCoarse < 0.0] = 0.0

        if provOn:
            # Per-class bulk removed from the partial top layer (proportional to
            # the bulk removed), plus the fully-eroded layers captured earlier.
            fracP = np.divide(
                H_old - eroVal, H_old, out=np.zeros_like(H_old), where=H_old > 0
            )
            P_part = self.stratP[nids, eroLayNb, :].copy()
            provBulkEro = provMaskedBulk + P_part * fracP[:, None]   # (n, classes)
            keepP = np.divide(eroVal, H_old, out=np.zeros_like(eroVal), where=H_old > 0)
            self.stratP[nids, eroLayNb, :] = P_part * keepP[:, None]

        # Update thickness of top stratigraphic layer
        self.stratH[nids, eroLayNb] = eroVal
        self.stratH[nids, 0] -= BEDROCK_SENTINEL
        if self.stratLith:
            for h, sl in zip(self.stratHL, sentL):
                h[nids, 0] -= sl
        if provOn:
            self.stratP[nids, 0, self.source_class[nids]] -= BEDROCK_SENTINEL
        neg = self.stratH < 0
        self.stratH[neg] = 0.0
        self.phiS[neg] = 0.0
        self.stratK[neg] = 0.0
        # Induration archive: emptied layers lose their crust record. NOT
        # forward-filled (unlike stratK/phiS) — stratDuri=0 is a valid value
        # (uncemented), and the exhumation re-arm (_recordInduration) reads the
        # top NON-EMPTY layer directly, so no fill is needed.
        if getattr(self, "stratDuri", None) is not None:
            self.stratDuri[neg] = 0.0
        # Emptied layers also lose their crust chemistry code (-1 = no crust).
        if getattr(self, "stratCrustType", None) is not None:
            self.stratCrustType[neg] = -1.0
        if getattr(self, "stratCrustSource", None) is not None:
            self.stratCrustSource[neg] = -1.0
        self.phiS[:, : self.stratStep + 1] = self._fillZeroPorosity(
            self.phiS[:, : self.stratStep + 1]
        )
        # Same forward-fill for the K multiplier so an emptied layer
        # inherits the bedrock K of the layer below it once exposed.
        self.stratK[:, : self.stratStep + 1] = self._fillZeroPorosity(
            self.stratK[:, : self.stratStep + 1]
        )
        if self.stratLith:
            top = self.stratStep + 1
            for h, p in zip(self.stratHL, self.phiL):
                h[neg] = 0.0
                p[neg] = 0.0
                np.clip(h, 0.0, None, out=h)
            if self.litNb == 1:
                # Fine bulk can never exceed the layer total thickness.
                np.minimum(
                    self.stratHf[:, :top], self.stratH[:, :top],
                    out=self.stratHf[:, :top],
                )
            else:
                self._capLithologies(slice(0, top))
            for p in self.phiL:
                p[:, :top] = self._fillZeroPorosity(p[:, :top])
        self.thCoarse /= self.dt
        if self.stratLith:
            own = self.inIDs == 1
            for l in range(self.litNb):
                self.thL[l] /= self.dt
                # Mass-balance diagnostic: solid removed this step (owned
                # nodes, m^3), as the uncompacted deposit-equivalent rate, the
                # same basis transported and deposited downstream.
                vol = float(np.sum((self.thL[l] * self.dt * self.larea)[own]))
                self._litEroded[l] += vol
                if self.litFine[l]:
                    self._fineEroded += vol
            # Routed fractions: the fine group, and the coarse group (coarse
            # clastic plus the coarse-like lithologies).
            self.thFine = self._fineGroup(self.thL).copy()
            for l in range(self.litNb):
                if not self.litFine[l]:
                    self.thCoarse = self.thCoarse + self.thL[l]

        if provOn:
            top = self.stratStep + 1
            np.clip(self.stratP[:, :top, :], 0.0, None, out=self.stratP[:, :top, :])
            # Re-impose Σ over classes == stratH per layer (absorbs the sentinel
            # round-off and the neg-thickness clamp); empty layers stay zero.
            psum = self.stratP[:, :top, :].sum(axis=2)
            scale = np.divide(
                self.stratH[:, :top], psum, out=np.ones_like(psum), where=psum > 0
            )
            self.stratP[:, :top, :] *= scale[:, :, None]

            # Per-class eroded sediment as the uncompacted deposit-equivalent
            # RATE: the total (self.thCoarse[+thFine]) split by the eroded bulk's
            # provenance composition. Σ over classes == the total, so the class
            # sub-fluxes routed in transport sum to the total flux.
            thTot = self.thCoarse.copy()
            if self.stratLith:
                thTot = thTot + self.thFine
            totBulk = provBulkEro.sum(axis=1)
            provComp = np.divide(
                provBulkEro, totBulk[:, None],
                out=np.zeros_like(provBulkEro), where=totBulk[:, None] > 0,
            )
            self.provEro = np.zeros((self.lpoints, self.provNb))
            self.provEro[nids] = provComp * thTot[nids][:, None]
            self._provEroded += np.sum(
                (self.provEro * self.dt * self.larea[:, None])[self.inIDs == 1], axis=0
            )

        return

    def elevStrat(self):
        """
        This function updates the current stratigraphic layer elevation.
        """

        self.stratZ[:, self.stratStep] = self.hLocal.getArray()

        return
    
    def _depthPorosity(self, depth):
        """
        This function uses the depth-porosity relationships to compute the porosities for each lithology and then the solid phase to get each layer thickness changes due to compaction.

        .. note::

            We assume that porosity cannot increase after unloading.

        :arg depth: depth below basement for each sedimentary layer

        :return: newH updated sedimentary layer thicknesses after compaction
        """

        if self.stratLith:
            return self._depthPorosityDual(depth)

        # Depth-porosity functions
        phiS = self.phi0s * np.exp(depth / self.z0s)
        phiS = np.minimum(phiS, self.phiS[:, : self.stratStep + 1])

        # Compute the solid phase in each layers
        tmpS = self.stratH[:, : self.stratStep + 1].copy()
        tmpS *= 1.0 - self.phiS[:, : self.stratStep + 1]
        solidPhase = tmpS

        # Get new layer thickness after porosity change
        tot = 1.0 - phiS[:, : self.stratStep + 1]

        ids = np.where(tot > 0.0)
        newH = np.zeros(tot.shape)
        newH[ids] = solidPhase[ids] / tot[ids]
        newH[newH <= 0] = 0.0
        phiS[newH <= 0] = 0.0

        # Freeze the bedrock sentinel so it neither compacts nor drives surface drop
        if self.bedrockLay > 0:
            b = self.bedrockLay
            newH[:, :b] = self.stratH[:, :b]      # thickness unchanged (1e6 sentinel)
            phiS[:, :b] = self.phiS[:, :b]        # porosity unchanged

        # Update porosities in each sedimentary layer
        phiS = self._fillZeroPorosity(phiS)
        self.phiS[:, : self.stratStep + 1] = phiS

        if self.memclear:
            del phiS, solidPhase
            del ids, tmpS, tot
            gc.collect()

        return newH

    def _depthPorosityDual(self, depth):
        """
        Dual-lithology depth-porosity / compaction. Each fraction (coarse,
        fine) compacts on its own porosity-depth curve; the layer's new
        thickness is the sum of the two recompacted bulk thicknesses, and the
        fine pile ``stratHf`` is updated consistently.

        This is the headline physics of dual lithology: fine and coarse have
        different compaction, so the same solid load yields a different
        preserved thickness depending on the layer's grain mix.

        :arg depth: depth below basement for each sedimentary layer
        :return: newH updated layer thicknesses after compaction
        """

        top = self.stratStep + 1
        phiS_cur = self.phiS[:, :top]
        phiL_cur = [p[:, :top] for p in self.phiL]
        H = self.stratH[:, :top]
        HL = [h[:, :top] for h in self.stratHL]
        Hc = H - self._lsum(HL)

        # Per-lithology depth-porosity, capped at the current value (porosity
        # cannot increase on unloading — same assumption as single-fraction).
        phiS_new = np.minimum(self.phi0c * np.exp(depth / self.z0c), phiS_cur)
        phiL_new = [
            np.minimum(self.litPhi0[l] * np.exp(depth / self.litZ0[l]), phiL_cur[l])
            for l in range(self.litNb)
        ]

        # Solid phase preserved per lithology.
        coarseSolid = Hc * (1.0 - phiS_cur)
        litSolid = [HL[l] * (1.0 - phiL_cur[l]) for l in range(self.litNb)]

        # Recompacted bulk thickness per lithology.
        totC = 1.0 - phiS_new
        Hc_new = np.zeros_like(Hc)
        idc = totC > 0.0
        Hc_new[idc] = coarseSolid[idc] / totC[idc]
        Hc_new[Hc_new <= 0] = 0.0
        HL_new = []
        for l in range(self.litNb):
            totL = 1.0 - phiL_new[l]
            hn = np.zeros_like(HL[l])
            idl = totL > 0.0
            hn[idl] = litSolid[l][idl] / totL[idl]
            hn[hn <= 0] = 0.0
            HL_new.append(hn)
        newH = Hc_new + HL_new[0]
        for hn in HL_new[1:]:
            newH = newH + hn

        # Emptied layers carry no porosity.
        empty = newH <= 0
        phiS_new[empty] = 0.0
        for pn in phiL_new:
            pn[empty] = 0.0

        # Freeze the bedrock sentinel layers (thickness + porosities unchanged).
        if self.bedrockLay > 0:
            b = self.bedrockLay
            newH[:, :b] = H[:, :b]
            phiS_new[:, :b] = phiS_cur[:, :b]
            for l in range(self.litNb):
                HL_new[l][:, :b] = HL[l][:, :b]
                phiL_new[l][:, :b] = phiL_cur[l][:, :b]

        phiS_new = self._fillZeroPorosity(phiS_new)
        self.phiS[:, :top] = phiS_new
        for l in range(self.litNb):
            self.phiL[l][:, :top] = self._fillZeroPorosity(phiL_new[l])
            self.stratHL[l][:, :top] = HL_new[l]

        if self.memclear:
            del phiS_new, phiL_new, coarseSolid, litSolid
            del Hc, HL, Hc_new, HL_new, totC
            gc.collect()

        return newH

    def getCompaction(self):
        """
        This function computes the change in sedimentary layers porosities and thicknesses due to compaction.

        .. note::

            We assume a simple depth-porosiy relationship.
        """

        t0 = process_time()
        topZ = np.vstack(self.hLocal.getArray())
        totH = np.sum(self.stratH[:, : self.stratStep + 1], axis=1)

        # Height of the sediment column above the center of each layer is given by
        cumZ = -np.cumsum(self.stratH[:, self.stratStep :: -1], axis=1) + topZ
        elev = np.append(topZ, cumZ[:, :-1], axis=1)
        zlay = np.fliplr(elev - np.fliplr(self.stratH[:, : self.stratStep + 1] / 2.0))

        # Compute lithologies porosities for each depositional layers
        # Get depth below basement
        depth = zlay - topZ

        # Now using depth-porosity relationships we compute the porosities
        newH = self._depthPorosity(depth)

        # Get the total thickness changes induced by compaction and
        # update the elevation accordingly
        dz = totH - np.sum(newH, axis=1)
        dz[dz <= 0] = 0.0
        self.hLocal.setArray(topZ.flatten() - dz.flatten())
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        # Update each layer thicknesses
        self.stratH[:, : self.stratStep + 1] = newH

        # Provenance: compaction only expels pore water — it is composition-
        # neutral (every source class's solid grains are preserved in the same
        # proportion), but it shrinks the layer thickness. Rescale stratP by the
        # same per-layer ratio so Σ over classes == the compacted stratH stays
        # exact; without this stratP keeps the pre-compaction (larger) thickness
        # and the post-processed fraction stratP[c]/stratH exceeds 1. Same renorm
        # idiom as erodeStrat and the advection step. (stratHf is already
        # recompacted consistently inside _depthPorosityDual.)
        if getattr(self, "provOn", False):
            top = self.stratStep + 1
            psum = self.stratP[:, :top, :].sum(axis=2)
            scale = np.divide(
                self.stratH[:, :top], psum, out=np.ones_like(psum), where=psum > 0
            )
            self.stratP[:, :top, :] *= scale[:, :, None]

        if self.memclear:
            del dz, newH, totH, topZ
            del depth, zlay, cumZ, elev
            gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Lithology Porosity Values (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return

    def _stratFields(self):
        """
        The per-layer stratigraphic fields present in this run, as a dict of
        local views ``name -> array[:, :stratStep(, provNb)]``. The single list
        of what a lateral move of the column must carry (used by
        `stratalRecord` and by plate advection, `tectonics._advectStrati`).
        """
        top = self.stratStep
        f = {"stratH": self.stratH[:, :top], "stratZ": self.stratZ[:, :top],
             "phiS": self.phiS[:, :top]}
        if getattr(self, "stratK", None) is not None:
            f["stratK"] = self.stratK[:, :top]
        if getattr(self, "stratLith", False):
            for l in range(self.litNb):
                f["stratHL%d" % l] = self.stratHL[l][:, :top]
                f["phiL%d" % l] = self.phiL[l][:, :top]
        for name in ("stratDuri", "stratCrustType", "stratCrustSource"):
            if getattr(self, name, None) is not None:
                f[name] = getattr(self, name)[:, :top]
        if getattr(self, "provOn", False):
            f["stratP"] = self.stratP[:, :top, :]
        return f

    def stratalRecord(self, indices, weights, onIDs, source=None):
        """
        Once the interpolation has been performed, the following function updates the stratigraphic records based on the advected mesh.

        The function relies on the fortran subroutine strataonesed. With
        lithology piles each pile (`stratHL[l]`, `phiL[l]`) is advected with
        its own strataonesed call — a pile is a bulk thickness, so the
        plain weighted interpolation applies directly (the unused `stratathreesed`
        kernel instead treats its extra fields as 0–1 fractions and renormalises
        them, which does not match the thickness representation here).

        Intensive per-layer properties (`phiS`, `phiL`, `stratDuri`, `stratK`)
        are interpolated with the same weights; the categorical crust codes take
        the record of the primary (highest-weight) source node.

        Scratch used: ``self.tmp`` (global) and ``self.tmpL`` (local), for the
        halo sync of every advected column; both hold the last synced column
        on exit.

        :arg indices: (lpoints, 3) rows of ``source`` used for each local node
        :arg weights: (lpoints, 3) interpolation weights
        :arg onIDs: local nodes that coincide with a source node (copied)
        :arg source: dict ``name -> rows`` as returned by `_stratFields`; the
            rows ``indices`` points into. Default: this rank's own layers
            (indices are local node ids). Plate advection passes the rows
            fetched from their owners (`tectonics._advectStrati`).
        """

        if source is None:
            source = self._stratFields()
        nsrc = source["stratH"].shape[0]
        top = self.stratStep
        has = lambda name: name in source  # noqa: E731

        def interp(field, ref_h=None, ref_z=None):
            # strataonesed interpolates three columns per call; the H/Z slots
            # are reused when only one intensive field is wanted.
            h = source["stratH"] if ref_h is None else ref_h
            z = source["stratZ"] if ref_z is None else ref_z
            a, b, c = strataonesed(nsrc, top, indices, weights, h, z, field)
            return a, b, c

        def pin(new, src_field):
            if len(onIDs) > 0:
                new[onIDs, :] = src_field[indices[onIDs, 0], :]
            return new

        # Get stratal dataset after displacements
        loc_stratH = source["stratH"]
        loc_stratZ = source["stratZ"]
        loc_phiS = source["phiS"]
        nstratH, nstratZ, nphiS = interp(loc_phiS)
        pin(nstratZ, loc_stratZ)
        pin(nstratH, loc_stratH)
        pin(nphiS, loc_phiS)

        if has("stratK"):
            _, _, nK = interp(source["stratK"])
            pin(nK, source["stratK"])

        nlit = sum(1 for k in source if k.startswith("stratHL"))
        nstratHL, nphiL = [], []
        for l in range(nlit):
            # Advect each lithology pile with the same interpolation (pile as
            # thickness, its porosity as the intensive field); the
            # re-interpolated elevation is identical to nstratZ above and
            # discarded.
            hk, pk = "stratHL%d" % l, "phiL%d" % l
            nh, _, npl = interp(source[pk], ref_h=source[hk])
            nstratHL.append(pin(nh, source[hk]))
            nphiL.append(pin(npl, source[pk]))

        if has("stratDuri"):
            # Induration archive: an INTENSIVE per-layer property (like phiS).
            _, _, nDuri = interp(source["stratDuri"])
            pin(nDuri, source["stratDuri"])

        if has("stratP"):
            # Each provenance class's per-layer thickness, with the same
            # interpolation (a bulk thickness, like stratHf). Linear weights
            # preserve Σ over classes == stratH (re-imposed below after sync).
            nstratP = np.zeros((indices.shape[0], top, self.provNb))
            for c in range(self.provNb):
                loc_P = np.ascontiguousarray(source["stratP"][:, :, c])
                nP, _, _ = interp(loc_phiS, ref_h=loc_P)
                nstratP[:, :, c] = pin(nP, loc_P)

        # Categorical crust chemistry codes advect by NEAREST neighbour — an
        # integer species/source label cannot be linearly interpolated, so each
        # column takes the record of its primary (highest-weight) source node.
        if has("stratCrustType"):
            nCrustType = source["stratCrustType"][indices[:, 0], :]
        if has("stratCrustSource"):
            nCrustSource = source["stratCrustSource"][indices[:, 0], :]

        # Make the records identical on every rank sharing a node: the OWNED
        # rows go into the global scratch Vec and are scattered back to the
        # local (ghost-including) arrays. (Until 2026-11 the full local column
        # was set into the global Vec, which is only the right size at np=1.)
        own = self.glIDs

        def sync(new_col, dst, k, rint=False):
            self.tmp.setArray(new_col[own])
            self.dm.globalToLocal(self.tmp, self.tmpL)
            v = self.tmpL.getArray()
            dst[:, k] = np.rint(v) if rint else v.copy()

        for k in range(top):
            sync(nstratZ[:, k], self.stratZ, k)
            sync(nstratH[:, k], self.stratH, k)
            sync(nphiS[:, k], self.phiS, k)
            if has("stratK"):
                sync(nK[:, k], self.stratK, k)
            for l in range(nlit):
                sync(nstratHL[l][:, k], self.stratHL[l], k)
                sync(nphiL[l][:, k], self.phiL[l], k)
            if has("stratDuri"):
                sync(nDuri[:, k], self.stratDuri, k)
            if has("stratCrustType"):
                sync(nCrustType[:, k], self.stratCrustType, k, rint=True)
            if has("stratCrustSource"):
                sync(nCrustSource[:, k], self.stratCrustSource, k, rint=True)
            if has("stratP"):
                for c in range(self.provNb):
                    self.tmp.setArray(nstratP[own, k, c])
                    self.dm.globalToLocal(self.tmp, self.tmpL)
                    self.stratP[:, k, c] = self.tmpL.getArray().copy()

        if has("stratP"):
            # Re-impose Σ over classes == stratH after the interpolation drift.
            np.clip(self.stratP[:, :top, :], 0.0, None, out=self.stratP[:, :top, :])
            psum = self.stratP[:, :top, :].sum(axis=2)
            scale = np.divide(
                self.stratH[:, :top], psum, out=np.ones_like(psum), where=psum > 0
            )
            self.stratP[:, :top, :] *= scale[:, :, None]

        return
