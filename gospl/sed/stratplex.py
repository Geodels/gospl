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

    def __init__(self):
        """
        The initialisation of `STRAMesh` class related to stratigraphic informations.
        """

        self.stratH = None
        self.stratZ = None
        self.phiS = None
        # Dual-lithology (coarse/fine) layer fields. Allocated only when
        # self.stratLith is True (see _extraStrata / DESIGN_DUAL_LITHOLOGY.md).
        # stratHf = fine-fraction layer thickness (coarse = stratH - stratHf);
        # phiF = fine porosity per layer (phiS then holds the coarse porosity).
        self.stratHf = None
        self.phiF = None
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
                self.stratHf = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
                if "strataHf" in fileData.files:
                    stratVal = fileData["strataHf"]
                    self.stratHf[:, lo:hi] = stratVal[self.locIDs, 0:nlay]
                    # Keep the partition physical: 0 <= fine <= total thickness.
                    np.clip(
                        self.stratHf[:, lo:hi],
                        0.0,
                        self.stratH[:, lo:hi],
                        out=self.stratHf[:, lo:hi],
                    )
                self.phiF = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
                if "phiF" in fileData.files:
                    stratVal = fileData["phiF"]
                    self.phiF[:, lo:hi] = stratVal[self.locIDs, 0:nlay]
                else:
                    # No fine porosity supplied: default to the fine surface
                    # porosity so fine-bearing initial layers compact sensibly
                    # (irrelevant where strataHf == 0).
                    self.phiF[:, lo:hi] = self.phi0f
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
                    self.stratHf[:, 0] = BEDROCK_SENTINEL * (
                        1.0 - self.bedrock_coarse_frac
                    )
                    self.phiF[:, 0] = self.phi0f
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
                self.stratHf = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
                self.phiF = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
                self.stratHf[:, 0] = BEDROCK_SENTINEL * (1.0 - self.bedrock_coarse_frac)
                self.phiF[:, 0] = self.phi0f

            self._logStratInit()

        if getattr(self, "provOn", False):
            self._initProvenance()

        return

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
        for k in ("strataZ", "phiS", "stratK", "strataHf", "phiF"):
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
        Hf = self.stratHf[:, b : self.stratStep + 1].sum(axis=1)
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
        Hf = self.stratHf[:, : self.stratStep + 1]
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

    def _surfaceLithoK(self):
        """
        Per-node erodibility multiplier from the exposed surface
        composition: ``fc + (1 - fc) * fine_k_factor``.

        Equals 1.0 everywhere when dual lithology is off (or
        ``fine_k_factor == 1``, i.e. no lithology contrast), so it composes
        multiplicatively with ``_surfaceK`` in the SPL flavours without
        altering single-fraction behaviour.
        """
        if not self.stratLith:
            return np.ones(self.lpoints, dtype=np.float64)
        fc = self._surfaceComposition()
        return fc + (1.0 - fc) * self.fine_k_factor

    def _surfaceLithoD(self):
        """
        Per-node diffusivity multiplier from the exposed surface composition:
        ``fc + (1 - fc) * fine_diff_factor``.

        Equals 1.0 everywhere when dual lithology is off (or
        ``fine_diff_factor == 1``, i.e. no contrast), so it composes
        multiplicatively with the base hillslope coefficients (``Cda``/``Cdm``)
        without altering single-fraction behaviour. Fines diffuse faster when
        ``fine_diff_factor > 1`` (see DESIGN_DUAL_LITHOLOGY.md Section 7).
        """
        if not self.stratLith:
            return np.ones(self.lpoints, dtype=np.float64)
        fc = self._surfaceComposition()
        return fc + (1.0 - fc) * self.fine_diff_factor

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
            fineDepo = depo * self.depoFineFrac
            self.stratHf[:, self.stratStep] += fineDepo
            # Mass-balance diagnostic (owned nodes only, m^3).
            self._fineDeposited += float(
                np.sum((fineDepo * self.larea)[self.inIDs == 1])
            )
            # Fresh deposit carries each lithology's surface porosity.
            self.phiS[ids, self.stratStep] = self.phi0c
            self.phiF[ids, self.stratStep] = self.phi0f
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
            provDepo = depo[:, None] * frac
            self.stratP[:, self.stratStep, :] += provDepo
            self._provDeposited += np.sum(
                (provDepo * self.larea[:, None])[self.inIDs == 1], axis=0
            )
        # Freshly deposited sediment carries the default erodibility (no
        # multiplier). If you want re-deposited sediment to keep its
        # source-layer K, this is the line to revisit.
        self.stratK[ids, self.stratStep] = 1.0

        # Cleaning arrays
        if self.memclear:
            del depo, ids
            gc.collect()

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
            if getattr(self, "provOn", False):
                self.provEro = np.zeros((self.lpoints, self.provNb))
            return

        # Cumulative thickness for each node
        if self.stratLith:
            # Inflate the fine pile in proportion to layer 0's CURRENT fine
            # fraction so the sentinel never alters that layer's composition
            # (for the infinite-bedrock case this fraction is 1-bedrock_coarse_frac;
            # for a real basal layer it is its own fine ratio, i.e. 0 if all coarse).
            H0 = self.stratH[nids, 0]
            sentFine = BEDROCK_SENTINEL * np.divide(
                self.stratHf[nids, 0], H0, out=np.zeros_like(H0), where=H0 > 0
            )
        self.stratH[nids, 0] += BEDROCK_SENTINEL
        if self.stratLith:
            self.stratHf[nids, 0] += sentFine
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
            fineTh = self.stratHf[nids, 0 : self.stratStep + 1]
            coarseTh = thickS - fineTh
            thCoarse = coarseTh * (1.0 - self.phiS[nids, 0 : self.stratStep + 1])
            thFine = fineTh * (1.0 - self.phiF[nids, 0 : self.stratStep + 1])
            thCoarse = np.sum((thCoarse * mask), axis=1)
            thFine = np.sum((thFine * mask), axis=1)
        else:
            thCoarse = thickS * (1.0 - self.phiS[nids, 0 : self.stratStep + 1])
            thCoarse = np.sum((thCoarse * mask), axis=1)

        # Clear all stratigraphy points which are eroded
        cumThick[boolMask] = 0.0
        tmp = self.stratH[nids, : self.stratStep + 1]
        tmp[boolMask] = 0
        self.stratH[nids, : self.stratStep + 1] = tmp
        if self.stratLith:
            tmpf = self.stratHf[nids, : self.stratStep + 1]
            tmpf[boolMask] = 0
            self.stratHf[nids, : self.stratStep + 1] = tmpf
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
            self.thFine = np.zeros(self.lpoints)
            Hf_old = self.stratHf[nids, eroLayNb]
            # Well-mixed layer: coarse/fine bulk are removed in proportion to
            # the layer composition; the remainder keeps the same ratio.
            frac = np.divide(tmp, H_old, out=np.zeros_like(tmp), where=H_old > 0)
            coarse_rm = (H_old - Hf_old) * frac
            fine_rm = Hf_old * frac
            thCoarse = thCoarse + coarse_rm * (1.0 - self.phiS[nids, eroLayNb])
            thFine = thFine + fine_rm * (1.0 - self.phiF[nids, eroLayNb])
            # Uncompacted thickness deposited downstream, per fraction, using
            # each lithology's surface porosity.
            self.thCoarse[nids] = thCoarse / (1.0 - self.phi0c)
            self.thFine[nids] = thFine / (1.0 - self.phi0f)
            self.thCoarse[self.thCoarse < 0.0] = 0.0
            self.thFine[self.thFine < 0.0] = 0.0
            # Remaining fine in the partial layer scales with remaining total.
            self.stratHf[nids, eroLayNb] = Hf_old * np.divide(
                eroVal, H_old, out=np.zeros_like(eroVal), where=H_old > 0
            )
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
            self.stratHf[nids, 0] -= sentFine
        if provOn:
            self.stratP[nids, 0, self.source_class[nids]] -= BEDROCK_SENTINEL
        neg = self.stratH < 0
        self.stratH[neg] = 0.0
        self.phiS[neg] = 0.0
        self.stratK[neg] = 0.0
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
            self.stratHf[neg] = 0.0
            self.phiF[neg] = 0.0
            np.clip(self.stratHf, 0.0, None, out=self.stratHf)
            # Fine bulk can never exceed the layer total thickness.
            np.minimum(
                self.stratHf[:, :top], self.stratH[:, :top],
                out=self.stratHf[:, :top],
            )
            self.phiF[:, :top] = self._fillZeroPorosity(self.phiF[:, :top])
        self.thCoarse /= self.dt
        if self.stratLith:
            self.thFine /= self.dt
            # Mass-balance diagnostic: fine solid removed this step (owned
            # nodes, m^3). thFine is the uncompacted deposit-equivalent rate,
            # the same basis transported and deposited downstream.
            self._fineEroded += float(
                np.sum((self.thFine * self.dt * self.larea)[self.inIDs == 1])
            )

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
        phiF_cur = self.phiF[:, :top]
        H = self.stratH[:, :top]
        Hf = self.stratHf[:, :top]
        Hc = H - Hf

        # Per-fraction depth-porosity, capped at the current value (porosity
        # cannot increase on unloading — same assumption as single-fraction).
        phiS_new = np.minimum(self.phi0c * np.exp(depth / self.z0c), phiS_cur)
        phiF_new = np.minimum(self.phi0f * np.exp(depth / self.z0f), phiF_cur)

        # Solid phase preserved per fraction.
        coarseSolid = Hc * (1.0 - phiS_cur)
        fineSolid = Hf * (1.0 - phiF_cur)

        # Recompacted bulk thickness per fraction.
        totC = 1.0 - phiS_new
        totF = 1.0 - phiF_new
        Hc_new = np.zeros_like(Hc)
        Hf_new = np.zeros_like(Hf)
        idc = totC > 0.0
        idf = totF > 0.0
        Hc_new[idc] = coarseSolid[idc] / totC[idc]
        Hf_new[idf] = fineSolid[idf] / totF[idf]
        Hc_new[Hc_new <= 0] = 0.0
        Hf_new[Hf_new <= 0] = 0.0
        newH = Hc_new + Hf_new

        # Emptied layers carry no porosity.
        empty = newH <= 0
        phiS_new[empty] = 0.0
        phiF_new[empty] = 0.0

        # Freeze the bedrock sentinel layers (thickness + porosities unchanged).
        if self.bedrockLay > 0:
            b = self.bedrockLay
            newH[:, :b] = H[:, :b]
            Hf_new[:, :b] = Hf[:, :b]
            phiS_new[:, :b] = phiS_cur[:, :b]
            phiF_new[:, :b] = phiF_cur[:, :b]

        phiS_new = self._fillZeroPorosity(phiS_new)
        phiF_new = self._fillZeroPorosity(phiF_new)
        self.phiS[:, :top] = phiS_new
        self.phiF[:, :top] = phiF_new
        self.stratHf[:, :top] = Hf_new

        if self.memclear:
            del phiS_new, phiF_new, coarseSolid, fineSolid
            del Hc, Hf, Hc_new, Hf_new, totC, totF
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

    def stratalRecord(self, indices, weights, onIDs):
        """
        Once the interpolation has been performed, the following function updates the stratigraphic records based on the advected mesh.

        The function relies on the fortran subroutine strataonesed. In
        dual-lithology mode the fine pile (`stratHf`, `phiF`) is advected with
        a second strataonesed call — `stratHf` is a bulk thickness, so the
        plain weighted interpolation applies directly (the unused `stratathreesed`
        kernel instead treats its extra fields as 0–1 fractions and renormalises
        them, which does not match the thickness representation here).

        :arg indices: indices of the closest nodes used for interpolation
        :arg weights: weights based on the distances to closest nodes
        :arg onIDs: index of nodes remaining at the same position.

        """

        # Get local stratal dataset after displacements
        loc_stratH = self.stratH[:, : self.stratStep]
        loc_stratZ = self.stratZ[:, : self.stratStep]
        loc_phiS = self.phiS[:, : self.stratStep]
        nstratH, nstratZ, nphiS = strataonesed(
            self.lpoints,
            self.stratStep,
            indices,
            weights,
            loc_stratH,
            loc_stratZ,
            loc_phiS,
        )

        if len(onIDs) > 0:
            nstratZ[onIDs, :] = loc_stratZ[indices[onIDs, 0], :]
            nstratH[onIDs, :] = loc_stratH[indices[onIDs, 0], :]
            nphiS[onIDs, :] = loc_phiS[indices[onIDs, 0], :]

        if self.stratLith:
            # Advect the fine pile with the same interpolation (stratHf as
            # thickness, phiF as porosity); the re-interpolated elevation is
            # identical to nstratZ above and discarded.
            loc_stratHf = self.stratHf[:, : self.stratStep]
            loc_phiF = self.phiF[:, : self.stratStep]
            nstratHf, _, nphiF = strataonesed(
                self.lpoints,
                self.stratStep,
                indices,
                weights,
                loc_stratHf,
                loc_stratZ,
                loc_phiF,
            )
            if len(onIDs) > 0:
                nstratHf[onIDs, :] = loc_stratHf[indices[onIDs, 0], :]
                nphiF[onIDs, :] = loc_phiF[indices[onIDs, 0], :]

        provOn = getattr(self, "provOn", False)
        if provOn:
            # Advect each provenance class's per-layer thickness with the same
            # interpolation (a bulk thickness, like stratHf). Linear weights
            # preserve Σ over classes == stratH (re-imposed below after sync).
            nstratP = np.zeros((self.lpoints, self.stratStep, self.provNb))
            for c in range(self.provNb):
                loc_P = np.ascontiguousarray(self.stratP[:, : self.stratStep, c])
                nP, _, _ = strataonesed(
                    self.lpoints, self.stratStep, indices, weights,
                    loc_P, loc_stratZ, loc_phiS,
                )
                if len(onIDs) > 0:
                    nP[onIDs, :] = loc_P[indices[onIDs, 0], :]
                nstratP[:, :, c] = nP

        # Updates stratigraphic records after mesh advection on the edges of each partition
        # to ensure that all stratigraphic information on the adjacent nodes of the neighbouring
        # partition are equals on all processors sharing a common number of nodes.
        for k in range(self.stratStep):
            self.tmp.setArray(nstratZ[:, k])
            self.dm.globalToLocal(self.tmp, self.tmpL)
            self.stratZ[:, k] = self.tmpL.getArray().copy()

            self.tmp.setArray(nstratH[:, k])
            self.dm.globalToLocal(self.tmp, self.tmpL)
            self.stratH[:, k] = self.tmpL.getArray().copy()
            self.tmp.setArray(nphiS[:, k])
            self.dm.globalToLocal(self.tmp, self.tmpL)
            self.phiS[:, k] = self.tmpL.getArray().copy()

            if self.stratLith:
                self.tmp.setArray(nstratHf[:, k])
                self.dm.globalToLocal(self.tmp, self.tmpL)
                self.stratHf[:, k] = self.tmpL.getArray().copy()
                self.tmp.setArray(nphiF[:, k])
                self.dm.globalToLocal(self.tmp, self.tmpL)
                self.phiF[:, k] = self.tmpL.getArray().copy()

            if provOn:
                for c in range(self.provNb):
                    self.tmp.setArray(nstratP[:, k, c])
                    self.dm.globalToLocal(self.tmp, self.tmpL)
                    self.stratP[:, k, c] = self.tmpL.getArray().copy()

        if provOn:
            # Re-impose Σ over classes == stratH after the interpolation drift.
            top = self.stratStep
            np.clip(self.stratP[:, :top, :], 0.0, None, out=self.stratP[:, :top, :])
            psum = self.stratP[:, :top, :].sum(axis=2)
            scale = np.divide(
                self.stratH[:, :top], psum, out=np.ones_like(psum), where=psum > 0
            )
            self.stratP[:, :top, :] *= scale[:, :, None]

        return
