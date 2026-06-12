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
        # Per-layer erodibility multiplier (lpoints, stratNb). Effective K
        # at the surface = self.K * stratK[<top non-empty layer>]. Fresh
        # deposits and the bedrock sentinel default to 1.0 (no scaling).
        # An optional `stratK` key in the npstrata file lets the user
        # impose a non-uniform multiplier on the initial layers.
        self.stratK = None

        return

    def readStratLayers(self):
        """
        When stratigraphic layers are turned on, this function reads any initial stratigraphic layers provided in the input file (key: `npstrata`).

        The following variables will be read from the file:

        - thickness of each stratigrapic layer `strataH` accounting for both erosion & deposition events.
        - elevation at time of deposition, considered to be to the current elevation for the top stratigraphic layer `strataZ`.
        - porosity of coarse sediment `phiS` in each stratigraphic layer computed at center of each layer.
        """

        if self.strataFile is not None:
            fileData = np.load(self.strataFile)
            stratVal = fileData["strataH"]
            self.initLay = stratVal.shape[1]
            self.stratNb += self.initLay

            # Create stratigraphic arrays
            self.stratH = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratH[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            stratVal = fileData["strataZ"]
            self.stratZ = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.stratZ[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            stratVal = fileData["phiS"]
            self.phiS = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
            self.phiS[:, 0 : self.initLay] = stratVal[self.locIDs, 0 : self.initLay]

            # Per-layer K multiplier. Loaded from the optional `stratK`
            # key in the npstrata file (same (mpoints, initLay) shape as
            # the other layer fields); defaults to 1.0 (use self.K as-is).
            self.stratK = np.ones((self.lpoints, self.stratNb), dtype=np.float64)
            if "stratK" in fileData.files:
                stratVal = fileData["stratK"]
                self.stratK[:, 0 : self.initLay] = stratVal[
                    self.locIDs, 0 : self.initLay
                ]

            # Dual-lithology fine-fraction layer fields. Optional in the
            # npstrata file (`strataHf`, `phiF`); absent -> all-coarse
            # (fine thickness 0), so a single-fraction file still loads.
            if self.stratLith:
                self.stratHf = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
                if "strataHf" in fileData.files:
                    stratVal = fileData["strataHf"]
                    self.stratHf[:, 0 : self.initLay] = stratVal[
                        self.locIDs, 0 : self.initLay
                    ]
                self.phiF = np.zeros((self.lpoints, self.stratNb), dtype=np.float64)
                if "phiF" in fileData.files:
                    stratVal = fileData["phiF"]
                    self.phiF[:, 0 : self.initLay] = stratVal[
                        self.locIDs, 0 : self.initLay
                    ]

            # All layers in the file are real sediment; no bedrock sentinel.
            self.bedrockLay = 0

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
            self.stratHf[:, self.stratStep] += depo * self.depoFineFrac
            # Fresh deposit carries each lithology's surface porosity.
            self.phiS[ids, self.stratStep] = self.phi0c
            self.phiF[ids, self.stratStep] = self.phi0f
        else:
            self.phiS[ids, self.stratStep] = self.phi0s
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
        cumThick = np.cumsum(self.stratH[nids, self.stratStep :: -1], axis=1)[:, ::-1]
        boolMask = cumThick < -ero[nids].reshape((len(nids), 1))
        mask = boolMask.astype(int)

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

        # Erode remaining stratal layers
        # Get non-zero top layer number
        eroLayNb = np.bincount(np.nonzero(cumThick)[0]) - 1
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

        # Update thickness of top stratigraphic layer
        self.stratH[nids, eroLayNb] = eroVal
        self.stratH[nids, 0] -= BEDROCK_SENTINEL
        if self.stratLith:
            self.stratHf[nids, 0] -= sentFine
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

        return
