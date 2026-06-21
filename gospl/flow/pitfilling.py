import os
import gc
import sys
import petsc4py
import numpy as np
import pandas as pd
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

from gospl.tools.constants import GRAPH_OUTLIER_CAP, MISSING_DATA_SENTINEL

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import edge_tile
    from gospl._fortran import fill_tile
    from gospl._fortran import fill_edges
    from gospl._fortran import fill_dir
    from gospl._fortran import nghb_dir
    from gospl._fortran import getpitvol
    from gospl._fortran import pits_cons
    from gospl._fortran import fill_depressions
    from gospl._fortran import graph_nodes
    from gospl._fortran import combine_edges
    from gospl._fortran import label_pits
    from gospl._fortran import spill_pts
    from gospl._fortran import sort_ids

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()


class PITFill(object):
    """
    Depression filling is an important preconditioning step to many landscape evolution models.

    This class implements a linearly-scaling parallel priority-flood depression-filling algorithm based on `Barnes (2016) <https://arxiv.org/pdf/1606.06204.pdf>`_ algorithm.

    .. note::

        Unlike previous algorithms, `Barnes (2016) <https://arxiv.org/pdf/1606.06204.pdf>`_ approach guarantees a fixed number of memory access and communication events per processors.

    The approach proposed in goSPL handles irregular meshes, allowing for complex distributed meshes to be used as long as a clear definition of inter-mesh connectivities is available. It also creates directions over flat regions allowing for downstream flows in cases where the entire volume of a depression is filled.

    For inter-mesh connections and message passing, the approach relies on PETSc DMPlex functions.

    **Parallel (MPI) structure — the three Barnes (2016) phases:**

    1. *Per-tile flood (parallel).* Each rank priority-floods its own partition
       (``fill_tile``) and emits a small spillover graph linking its local
       depression labels by their minimum spill elevation (``graph_nodes`` /
       ``combine_edges``).
    2. *Global graph solve (serial, rank 0).* The per-rank graphs are gathered,
       labels describing the same basin across a partition border are unified
       (:meth:`_transferIDs` → :meth:`_unifyLabels`), and a single priority-flood
       on the combined graph (:meth:`_fillFromEdges` → ``fill_edges``) computes
       each depression's final water level. The graph scales with the partition
       *perimeter* (``O(sqrt(N))`` per tile), so — as in Barnes (2016) — this
       master solve is intentionally serial and cheap even on large meshes.
    3. *Apply (parallel).* Each rank raises its cells to the per-label water
       levels (``fill_depressions``).

    .. note::

        Cross-rank label unification (:meth:`_unifyLabels`) is a single
        deterministic **union-find** pass that collapses each connected component
        of the equivalence graph to its minimum label — so the depression
        identifiers (and the filled surface) are **independent of the MPI
        partition**. The serial graph solve uses a compact (densely re-indexed)
        label space and a binary-heap priority queue, so its cost tracks the
        number of distinct spillover basins rather than the rank-dependent raw
        (globally-offset) label range.

    The main functions return the following parameters:

    - the elevation of the filled surface,
    - the information for each depression (e.g., a unique global ID, its spillover local points and related processor),
    - the description of each depression (total volume and maximum filled depth).
    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `PITFill` class consists in the declaration of PETSc vectors, matrices and each partition internals edge vertices.
        """

        # Petsc vectors. Domain spill-over outlets are the DRAINING borders
        # (`outletIDs` = open/fixed edges); true-wall edges are excluded so a
        # depression touching a wall fills and is contained rather than spilling
        # out of the domain. (`idLBounds` are the parallel partition-local
        # boundaries — unrelated to the physical domain boundary — and always
        # participate in the spill graph.)
        edges = -np.ones((self.lpoints, 2), dtype=int)
        edges[self.idLBounds, 0] = self.idLBounds
        edges[self.outletIDs, 0] = self.outletIDs
        edges[self.idLBounds, 1] = 0
        edges[self.outletIDs, 1] = 1
        self.borders = edges

        self.outEdges = np.zeros(self.lpoints, dtype=int)
        self.outEdges[self.ghostIDs] = 1

        return

    def _buildPitDataframe(self, label1, label2):
        """
        Definition of a Pandas Dataframe used to find unique pit ID across processors.

        :arg label1: depression ID in a given processor
        :arg label2: same depression ID in a neighbouring mesh

        :return: df (sorted Dataframe of pit ID between processors)
        """

        data = {
            "p1": label1,
            "p2": label2,
        }
        df = pd.DataFrame(data, columns=["p1", "p2"])
        df = df.drop_duplicates().sort_values(["p2", "p1"], ascending=(False, False))

        return df

    def _sortingPits(self, df):
        """
        Sorts depressions number before combining them to ensure no depression index is changed in an unsorted way.

        :arg df: Pandas Dataframe containing depression numbers which have to be combined.

        :return: df sorted pandas dataframe containing depression numbers.

        .. note::
            Superseded by :meth:`_unifyLabels` (union-find). Kept for reference
            / the standalone sort_ids kernel; no longer on the hot path.
        """

        df["p2"] = sort_ids(df["p1"].values.astype(int), df["p2"].values.astype(int))
        df = df.drop_duplicates().sort_values(["p2", "p1"], ascending=(False, False))

        return df

    def _unifyLabels(self, df, label):
        """
        Collapse cross-processor depression-label equivalence pairs into one
        canonical id per connected component using union-find (disjoint-set),
        then apply the mapping to ``label``.

        ``df`` holds the global equivalence pairs ``(p1, p2)`` with ``p1 < p2``
        (two labels that are the same depression across a tile border). Each
        connected component is collapsed to its MINIMUM label by keeping the
        smaller label as the set root — matching the previous iterative
        ``sort_ids`` fixpoint, but in a single deterministic O(N α(N)) pass.
        Determinism (root = min, independent of pair order) is what makes the
        result identical regardless of the domain partition.

        :arg df: Dataframe of equivalence pairs, columns ``p1`` < ``p2``.
        :arg label: local depression-label array (globally-offset ids).

        :return: ``label`` with every label remapped to its component minimum.
        """

        pairs = df[["p1", "p2"]].values.astype(np.int64)

        # Union-find over the (sparse, globally-offset) labels that appear in a
        # pair. Most mesh labels never cross a border, so the structure is small.
        parent = {}

        def find(x):
            root = x
            while parent[root] != root:
                root = parent[root]
            while parent[x] != root:  # path compression
                parent[x], x = root, parent[x]
            return root

        for a, b in pairs:
            a, b = int(a), int(b)
            if a not in parent:
                parent[a] = a
            if b not in parent:
                parent[b] = b
            ra, rb = find(a), find(b)
            if ra != rb:
                # smaller label becomes the root → component minimum is canonical
                if ra < rb:
                    parent[rb] = ra
                else:
                    parent[ra] = rb

        if not parent:
            return label

        # Vectorised remap: only labels appearing in a pair change; all others
        # (the vast majority, plus the -1 border sentinel) pass through.
        keys = np.fromiter(sorted(parent), dtype=np.int64, count=len(parent))
        vals = np.fromiter((find(int(k)) for k in keys), dtype=np.int64,
                           count=len(keys))
        idx = np.searchsorted(keys, label)
        idx_clip = np.clip(idx, 0, len(keys) - 1)
        hit = keys[idx_clip] == label
        out = label.copy()
        out[hit] = vals[idx_clip[hit]]

        return out

    def _offsetGlobal(self, lgth):
        """
        Computes the offset between processors to ensure unique number for considered indices.

        :arg lgth: local length of the data to distribute

        :return: cumulative sum and sum of the labels to add to each processor
        """

        label_offset = np.zeros(MPIsize + 1, dtype=int)
        label_offset[MPIrank + 1] = lgth
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, label_offset, op=MPI.MAX)

        return np.cumsum(label_offset), np.sum(label_offset)

    def _fillFromEdges(self, mgraph):
        """
        Combine local meshes by joining their edges based on local spillover graphs.

        :arg mgraph: Numpy Array containing local mesh edges information

        :arg ggraph: Numpy Array containing filled elevation values based on other processors values
        """

        # Get bidirectional edges connections (pandas drop_duplicates is
        # hash-based and faster than any numpy sort-based dedup here — verified
        # by micro-benchmark; this dedup is <~30 ms even at 200k edges and is
        # NOT the pitgraph hot spot).
        cgraph = pd.DataFrame(
            mgraph, columns=["source", "target", "elev", "spill", "rank"]
        )
        cgraph = cgraph.sort_values("elev")
        cgraph = cgraph.drop_duplicates(["source", "target"], keep="first")
        cgraph = cgraph.values

        # --- Pit-label compaction --------------------------------------------
        # source/target are globally-offset pit labels (disjoint per rank), so
        # the label space is SPARSE and its maximum grows ~linearly with the
        # rank count even though the number of distinct spillover pits is far
        # smaller. fill_edges sizes its O(nb) / O(nb*maxnghbs) work arrays from
        # that maximum (nb = max(label)+2), so the serial solve grows with
        # ranks. Relabel to a dense [0, nlab) space, solve there, then scatter
        # the results back onto the global label index the downstream code
        # expects. Label 0 (the global outlet, == node 1, the priority-flood
        # seed) is the smallest label so it always maps to compact 0; union it
        # in explicitly in case no edge references it directly. This is a pure
        # relabeling: the fill result is identical to solving in global space.
        nb_global = int(cgraph[:, 1].max()) + 2
        src = cgraph[:, 0].astype(np.int64)
        tgt = cgraph[:, 1].astype(np.int64)
        uniq = np.unique(np.concatenate((src, tgt)))
        if uniq[0] != 0:
            uniq = np.concatenate((np.array([0], dtype=uniq.dtype), uniq))
        nlab = len(uniq)
        ccgraph = cgraph.copy()
        ccgraph[:, 0] = np.searchsorted(uniq, src)
        ccgraph[:, 1] = np.searchsorted(uniq, tgt)
        c12 = ccgraph[:, :2].astype(int).ravel()
        cmax = np.max(np.bincount(c12)) + 1

        # Filling the bidirectional graph (serial priority-flood; timed
        # separately so we can see whether the Fortran graph solve or the
        # surrounding MPI Reduce/bcast dominates the serial pitgraph cost).
        # nlab+1 always covers every compact label 0..nlab-1 (the +1 keeps the
        # one-node slack the global formula had).
        self.profiler.start("flow_pit_edges")
        elev, rank, nodes, spillID = fill_edges(nlab + 1, ccgraph, cmax)
        self.profiler.stop("flow_pit_edges")

        # Scatter compact results back to the global label index. Labels with
        # no spillover edge are absent from the compact graph; in global space
        # fill_edges would have left them at its nelev sentinel (1e8), which the
        # caller rewrites to MISSING_DATA_SENTINEL — reproduce that default so
        # behaviour is identical. spillID is itself a (compact) label, so map it
        # back to the global label too (nodes/rank are mesh ids, not labels).
        ggraph = -np.ones((nb_global, 5))
        ggraph[:, 0] = 1.0e8
        ggraph[uniq, 0] = elev[:nlab]
        ggraph[uniq, 1] = nodes[:nlab]
        ggraph[uniq, 2] = rank[:nlab]
        gspill = spillID[:nlab].astype(np.int64)
        valid = (gspill >= 0) & (gspill < nlab)
        gspill[valid] = uniq[gspill[valid]]
        ggraph[uniq, 3] = gspill
        if self.memclear:
            del elev, nodes, rank, spillID
            del src, tgt, uniq, c12, gspill, ccgraph

        return ggraph

    def _transferIDs(self, pitIDs):
        """
        This function transfers local depression IDs along local borders and combines them with a unique identifier.

        :arg pitIDs: local depression index.

        :return: number of depressions.
        """

        # Define globally unique watershed index
        t0 = process_time()
        fillIDs = pitIDs >= 0
        offset, _ = self._offsetGlobal(np.amax(pitIDs))
        pitIDs[fillIDs] += offset[MPIrank]
        if MPIrank == 0 and self.verbose:
            print(
                "Offset Global pit ids (%0.02f seconds)" % (process_time() - t0)
            )
        t0 = process_time()

        # Transfer depression IDs along local borders
        self.tmpL.setArray(pitIDs)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        label = self.tmpL.getArray().copy().astype(int)

        ids = label < pitIDs
        df = self._buildPitDataframe(label[ids], pitIDs[ids])
        ids = label > pitIDs
        df2 = self._buildPitDataframe(pitIDs[ids], label[ids])
        df = pd.concat([df, df2], ignore_index=True)
        df = df.drop_duplicates().sort_values(["p2", "p1"], ascending=(False, False))
        df = df[(df["p1"] >= 0) & (df["p2"] >= 0)]
        if MPIrank == 0 and self.verbose:
            print(
                "Build pit dataframe (%0.02f seconds)" % (process_time() - t0)
            )
        t0 = process_time()

        # Gather every rank's equivalence pairs onto all ranks. Allgatherv of
        # the real pairs (Barnes-style) replaces a padded Allreduce(MAX) over a
        # (total, 2) array — same content, no -1 padding to move. All ranks then
        # run the identical union-find.
        offset, total = self._offsetGlobal(len(df))
        sendbuf = (
            df.values.astype(np.int64).ravel()
            if len(df) > 0
            else np.empty(0, dtype=np.int64)
        )
        counts = (np.diff(offset) * 2).astype(int)
        recvbuf = np.empty(int(total) * 2, dtype=np.int64)
        MPI.COMM_WORLD.Allgatherv(sendbuf, [recvbuf, counts])
        combIds = recvbuf.reshape(-1, 2)
        df = self._buildPitDataframe(combIds[:, 0], combIds[:, 1])
        df = df[(df["p1"] >= 0) & (df["p2"] >= 0)]
        if MPIrank == 0 and self.verbose:
            print(
                "Combine pit dataframe (%0.02f seconds)" % (process_time() - t0)
            )
        t0 = process_time()

        # Collapse the cross-processor label-equivalence pairs (p1 < p2) into a
        # single canonical id per connected component with union-find. This is
        # Barnes-style label unification: the canonical id is the MINIMUM label
        # in each component, identical to the previous iterative sort_ids
        # fixpoint but in a single deterministic pass — no convergence loop (the
        # old code iterated up to 1000 times) and no partition-dependent
        # ordering, which removes the instability that forced the serial path.
        if len(df) > 0:
            label = self._unifyLabels(df, label)
        if MPIrank == 0 and self.verbose:
            print(
                "Sorting pits (%0.02f seconds)" % (process_time() - t0)
            )
        t0 = process_time()

        # Transfer depression IDs along local borders
        self.tmpL.setArray(label)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)

        # At this point all pits have a unique IDs across processors
        self.pitIDs = self.tmpL.getArray().astype(int)

        # Lets make consecutive indices
        pitnbs = np.zeros(1, dtype=int)
        pitnbs[0] = np.max(self.pitIDs) + 1
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, pitnbs, op=MPI.MAX)
        fillIDs = self.pitIDs >= 0
        valpit = -np.ones(pitnbs[0], dtype=int)
        unique, idx_groups = npi.group_by(
            self.pitIDs[fillIDs], np.arange(len(self.pitIDs[fillIDs]))
        )
        valpit[unique] = 1
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, valpit, op=MPI.MAX)
        pitNb = np.where(valpit > 0)[0]
        if MPIrank == 0 and self.verbose:
            print(
                "Define consecutive pit ids (%0.02f seconds)" % (process_time() - t0)
            )
        t0 = process_time()

        self.pitIDs = pits_cons(self.pitIDs, pitNb)

        return pitNb

    def _dirFlats(self):
        """
        This function finds routes to spillover points on filled depressions to ensure downstream distribution if they are overfilled.
        """

        # Find local receivers to direct flow to spillover nodes
        ids = self.pitIDs > -1
        if ids.any():
            pdir = fill_dir(self.lspillIDs, self.pitIDs, self.lFill)
        else:
            pdir = -np.ones(self.lpoints, dtype=np.int32)

        # Iteratively propagate flat-region directions from the spill points
        # across partition boundaries. Each pass advances the spill->interior
        # BFS front one cell (plus one halo layer across partitions). The OLD
        # fixed 100-pass cap stopped early on a large flat region (big endorheic
        # basin / flat shelf on a high-res global mesh) or one cut unfavourably
        # across ranks, leaving its interior with inconsistent cross-partition
        # `ptdir` that fill_rcvs (run independently per rank) wired into
        # mutually-pointing receivers -> a cycle -> a singular (I - W^T)
        # sub-block that diverged the flow/sediment KSP (P=48/96/192, not 144 --
        # partition-dependent). Instead of guessing a cap, run until either
        # every pit node is directed (remaining == 0) or a pass directs no new
        # node anywhere (remaining stops shrinking -> the rest is unreachable
        # from any spill seed; spinning further cannot help). Same one collective
        # per pass as before (a SUM count instead of a MIN), with a high
        # absolute backstop. Any leftover undirected node becomes a self-sink in
        # fill_rcvs, which is safe (non-singular).
        remaining = np.array([-1.0])
        prev_remaining = -1.0
        stall = 0
        stp = 0
        while True:
            self.tmp.setArray(pdir[self.glIDs])
            self.dm.globalToLocal(self.tmp, self.tmpL, 3)
            pdir = self.tmpL.getArray().copy().astype(int)
            pdir = nghb_dir(self.pitIDs, self.lFill, pdir)

            # Global count of pit nodes still without a direction.
            if ids.any():
                remaining[0] = float(np.count_nonzero(pdir[ids] < 0))
            else:
                remaining[0] = 0.0
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, remaining, op=MPI.SUM)
            stp += 1

            if remaining[0] == 0.0:
                break  # fully directed

            # Stagnation: two consecutive passes with no net progress means the
            # remainder is unreachable from any spill seed (one stall pass can be
            # a transient waiting on a halo hop, so require two).
            if prev_remaining >= 0.0 and remaining[0] >= prev_remaining:
                stall += 1
            else:
                stall = 0
            prev_remaining = remaining[0]

            if stall >= 2 or stp > 10000:
                # A handful of flat nodes can't be directed to a spill and fall
                # back to sinks. The pinpoint diagnostic (P=240 sweep) showed
                # these are ISOLATED filled pockets: single cells / tiny clusters
                # whose neighbours are non-pit (higher) terrain or a different,
                # higher pit -- they share a watershed's pitID+spill but have no
                # pit-connected path to it and no downhill exit at the fill level.
                # A sink is therefore the physically-correct local outcome (they
                # pond; mass stays in the pit). This is a depression-merge
                # labelling artefact, NOT a routing failure -- count it (owned
                # only) so the magnitude stays visible, but don't call it a bug.
                owned = self.inIDs == 1
                unreached = owned & (self.pitIDs > -1) & (pdir < 0)
                ncnt = np.array([float(np.count_nonzero(unreached))])
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, ncnt, op=MPI.SUM)
                if MPIrank == 0 and self.verbose:
                    print(
                        "[pitfill] %d isolated flat node(s) had no pit-path to a "
                        "spill -> left as sinks (pond locally; mass conserved)"
                        % int(ncnt[0]),
                        flush=True,
                    )
                # Verbose-only pinpoint diagnostic: dump the first unreachable
                # owned node's local neighbourhood (same-pit? directed? owned vs
                # ghost?) -- used to characterise the isolated-pocket pattern.
                if self.verbose and not getattr(self, "_pitDiagDone", False):
                    uidx = np.where(unreached)[0]
                    if uidx.size > 0:
                        try:
                            n0 = int(uidx[0])
                            bits = []
                            for c in self.FVmesh_ngbID[n0]:
                                c = int(c)
                                if c < 0:
                                    continue
                                bits.append(
                                    "g%d/pit%d/pdir%d/%s"
                                    % (
                                        self.gid[c],
                                        self.pitIDs[c],
                                        pdir[c],
                                        "own" if self.inIDs[c] == 1 else "ghost",
                                    )
                                )
                            print(
                                "[pitfill-diag] rank%d node g%d pit%d pdir%d "
                                "spill_g%d | nbrs: %s"
                                % (
                                    MPIrank,
                                    self.gid[n0],
                                    self.pitIDs[n0],
                                    pdir[n0],
                                    self.pitInfo[self.pitIDs[n0], 0],
                                    " ".join(bits),
                                ),
                                flush=True,
                            )
                            self._pitDiagDone = True
                        except Exception:
                            pass
                break

        # Final halo refresh so fill_rcvs sees `ptdir` consistent across
        # partition boundaries. The last nghb_dir advanced owned nodes only,
        # leaving ghost directions one pass stale; without this refresh two ranks
        # can pick mutually-pointing receivers on a shared boundary -> a 2-cycle
        # in the assembled matrix.
        self.tmp.setArray(pdir[self.glIDs])
        self.dm.globalToLocal(self.tmp, self.tmpL, 3)
        pdir = self.tmpL.getArray().copy().astype(int)

        # Define receiver nodes on each depression — PARTITION-INVARIANT.
        #
        # Each flat node routes to a neighbour strictly closer to the spill
        # (smaller `pdir`, not higher: lFill[c] <= lFill[i]). `pdir` is itself
        # partition-invariant, but in a flat region MANY neighbours tie at the
        # minimum distance, and the fortran `fill_rcvs` broke that tie by
        # "first in FVnID neighbour order" — which is partition-DEPENDENT (local
        # numbering/halo), so the same node picked a different (equally-valid)
        # receiver at different rank counts (~0.5% of nodes on a 370k earth
        # mesh). That made the assembled (I - W^T) operator partition-dependent
        # → the near-singular pinches / KSP failures that varied with the cut.
        #
        # We resolve the tie on `locIDs` (the input-mesh id), which IS
        # partition-invariant (unlike `self.gid`, PETSc's per-partition
        # numbering). Among the strictly-closer, not-higher neighbours we pick
        # the smallest `pdir`, ties broken by smallest `locIDs` — identical
        # across any decomposition. Replaces the fortran `fill_rcvs` tie-break;
        # non-tied cases are unchanged.
        ngb = self.FVmesh_ngbID                       # (n, K) local ids, -1 pad
        n = self.lpoints
        valid = ngb >= 0
        cidx = np.where(valid, ngb, 0)
        cpd = pdir[cidx]
        ch = self.lFill[cidx]
        clid = self.locIDs[cidx].astype(np.int64)
        cand = (
            valid
            & (cpd > -1)
            & (cpd < pdir[:, None])
            & (ch <= self.lFill[:, None])
            & (self.pitIDs[:, None] > -1)
        )
        BIG = np.iinfo(np.int64).max
        nid = np.int64(self.locIDs.max()) + 1         # lexicographic (pdir, locID)
        comp = np.where(cand, cpd.astype(np.int64) * nid + clid, BIG)
        kbest = comp.argmin(axis=1)
        rows = np.arange(n)
        has = comp[rows, kbest] < BIG
        self.flatDirs = np.full(n, -1, dtype=np.int32)
        ispit = self.pitIDs > -1
        # pit node with a candidate -> chosen receiver; pit node w/o -> self (sink)
        self.flatDirs[ispit] = np.where(
            has[ispit], ngb[rows, kbest][ispit], rows[ispit]
        ).astype(np.int32)
        self.flatDirs[self.lspillIDs] = -1

        return

    def _getPitParams(self, hl, nbpits):
        """
        Define depression global parameters:

        - volume of each depression
        - maximum filled depth

        :arg hl: Numpy Array of unfilled surface elevation
        :arg nbpits: number of depression in the global mesh
        """

        # Get pit parameters (volume and maximum filled elevation)
        ids = self.inIDs == 1
        grp = npi.group_by(self.pitIDs[ids])
        uids = grp.unique
        _, vol = grp.sum((self.lFill[ids] - hl[ids]) * self.larea[ids])
        _, hh = grp.max(self.lFill[ids])
        _, dh = grp.max(self.lFill[ids] - hl[ids])
        totv = np.zeros(nbpits, dtype=np.float64)
        hmax = np.full(nbpits, MISSING_DATA_SENTINEL, dtype=np.float64)
        diffh = np.zeros(nbpits, dtype=np.float64)
        ids = uids > -1
        totv[uids[ids]] = vol[ids]
        hmax[uids[ids]] = hh[ids]
        diffh[uids[ids]] = dh[ids]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, totv, op=MPI.SUM)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, hmax, op=MPI.MAX)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, diffh, op=MPI.MAX)

        self.pitParams = np.empty((nbpits, 3), dtype=np.float64)
        self.pitParams[:, 0] = totv
        self.pitParams[:, 1] = hmax
        self.pitParams[:, 2] = diffh

        return

    def _performFilling(self, hl, level, sed):
        """
        This functions implements the linearly-scaling parallel priority-flood depression-filling algorithm from `Barnes (2016) <https://arxiv.org/pdf/1606.06204.pdf>`_ but adapted to unstructured meshes.

        :arg hl: local elevation.
        :arg level: minimal elevation above which the algorithm is performed.
        :arg sed: boolean specifying if the pits are filled with water or sediments.
        """

        t0 = process_time()

        # Get local meshes edges and communication nodes
        ledges = edge_tile(0.0, self.borders, hl)
        out = np.where(ledges >= 0)[0]
        localEdges = np.empty((len(out), 2), dtype=int)
        localEdges[:, 0] = np.arange(self.lpoints)[out].astype(int)
        localEdges[:, 1] = ledges[out]
        out = np.where(ledges == -2)[0]
        inIDs = self.inIDs.copy()
        inIDs[out] = 0
        outEdges = self.outEdges.copy()
        outEdges[out] = 0
        out = np.where((ledges == 0) | (ledges == 2))[0]
        gBounds = np.zeros(self.lpoints, dtype=int)
        gBounds[out] = 1

        # Local pit filling
        lFill, label, gnb = fill_tile(localEdges, hl, inIDs)

        # Graph associates label pairs with the minimum spillover elevation
        lgth = 0
        if gnb > 0:
            graph = graph_nodes(gnb)
            lgth = np.amax(graph[:, :2])

        # Define globally unique watershed index
        offset, _ = self._offsetGlobal(lgth)
        label += offset[MPIrank]
        if lgth > 0:
            graph[:, 0] += offset[MPIrank]
            ids = np.where(graph[:, 1] > 0)[0]
            graph[ids, 1] += offset[MPIrank]

        # Transfer watershed values along local borders
        self.tmpL.setArray(label)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        label = self.tmpL.getArray().copy().astype(int)

        # Transfer filled values along the local borders
        self.tmpL.setArray(lFill)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        lFill = self.tmpL.getArray()

        # Combine tiles edges
        cgraph, graphnb = combine_edges(lFill, label, localEdges[:, 0], outEdges)

        lgrph = 0
        if graphnb > 0 and lgth > 0:
            cgraph = np.concatenate((graph, cgraph[:graphnb]))
            lgrph = len(cgraph)
        elif graphnb > 0 and lgth == 0:
            cgraph = cgraph[:graphnb]
            lgrph = len(cgraph)
        elif graphnb == 0 and lgth > 0:
            cgraph = graph
            lgrph = len(cgraph)

        # Add processor number to the graph
        offset, total = self._offsetGlobal(lgrph)

        if total == 0:
            # Fully closed domain (e.g. every edge a wall): there is no
            # spill-over graph at all, so every depression is endorheic and the
            # local fill from fill_tile is already final. A ggraph left at the
            # MISSING_DATA_SENTINEL default makes fill_depressions keep that
            # local fill (no global merge, no outflow) — and skips the serial
            # rank-0 solve / Reduce / bcast entirely.
            nlab = int(label.max()) + 2 if label.size else 2
            graph = -np.ones((nlab, 5), dtype=float)
        else:
            graph = -np.ones((total, 5), dtype=float)
            if lgrph > 0:
                graph[offset[MPIrank] : offset[MPIrank] + lgrph, :4] = cgraph
                graph[offset[MPIrank] : offset[MPIrank] + lgrph, 4] = MPIrank

            # Build global spillover graph on master. SERIAL on rank 0 (the
            # depression-merge graph) + a Reduce-to-root and a bcast — the
            # suspected non-scaling part of `flow`. Timed as "flow_pitgraph"
            # (no-op when profiling off) to confirm whether it (vs the parallel
            # solves) drives the high-rank flow plateau.
            self.profiler.start("flow_pitgraph")
            if MPIrank == 0:
                mgraph = -np.ones((total, 5), dtype=float)
            else:
                mgraph = None
            self.profiler.start("flow_pit_reduce")
            MPI.COMM_WORLD.Reduce(graph, mgraph, op=MPI.MAX, root=0)
            self.profiler.stop("flow_pit_reduce")
            # _fillFromEdges (serial rank-0 graph solve) self-times
            # "flow_pit_edges".
            if MPIrank == 0:
                ggraph = self._fillFromEdges(mgraph)
            else:
                ggraph = None

            # Send filled graph dataset to each processors
            self.profiler.start("flow_pit_bcast")
            graph = MPI.COMM_WORLD.bcast(ggraph, root=0)
            self.profiler.stop("flow_pit_bcast")
            self.profiler.stop("flow_pitgraph")

        # Drain pit on local boundaries and towards mesh edges
        keep = graph[:, 2].astype(int) == MPIrank
        proc = -np.ones(len(graph))
        proc[keep] = graph[keep, 1]
        keep = proc > -1
        proc[keep] = gBounds[proc[keep].astype(int)]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, proc, op=MPI.MAX)
        ids = np.where(proc == 1)[0]
        ids2 = np.where(graph[ids, 0] == graph[graph[ids, 3].astype(int), 0])[0]
        graph[ids[ids2], 3] = 0.0
        # Sentinel filter: physical Earth elevations are bounded by ~9 km;
        # entries outside [MISSING_DATA_SENTINEL, GRAPH_OUTLIER_CAP] m come
        # from the global-graph initialiser (MISSING_DATA_SENTINEL) and any
        # extreme outliers, both flagged as MISSING_DATA_SENTINEL.
        graph[graph[:, 0] < MISSING_DATA_SENTINEL, 0] = MISSING_DATA_SENTINEL
        graph[graph[:, 0] > GRAPH_OUTLIER_CAP, 0] = MISSING_DATA_SENTINEL

        # Define global solution by combining depressions/flat together
        lFill = fill_depressions(0.0, hl, lFill, label.astype(int), graph[:, 0])

        # Define filling in land and enclosed seas only
        if not sed:
            id = lFill < self.sealevel - level
            lFill[id] = hl[id]
        self.tmpL.setArray(lFill)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        self.lFill = self.tmpL.getArray().copy()

        if MPIrank == 0 and self.verbose:
            print("Remove depressions (%0.02f seconds)" % (process_time() - t0))

        return

    def _pitInformation(self, hl, level, sed=False):
        """
        This function extracts depression information available to all processors. It stores the following:

        - the information for each depression (*e.g.*, unique global ID, its spillover local points and related processor),
        - the description of each depression (total volume and maximum filled depth).

        :arg hl: local elevation.
        :arg sed: boolean specifying if the pits are filled with water or sediments.
        """

        t0 = process_time()

        # Combine pits locally to get a unique local ID per depression
        if sed:
            pitIDs = label_pits(level, self.lFill)
        else:
            pitIDs = label_pits(self.sealevel, self.lFill)
        if MPIrank == 0 and self.verbose:
            print(
                "Define pit labels (%0.02f seconds)" % (process_time() - t0)
            )

        t0 = process_time()
        pitIDs[self.outletIDs] = -1   # draining borders aren't pits; walls can be
        pitNb = self._transferIDs(pitIDs)
        if MPIrank == 0 and self.verbose:
            print(
                "Define transfer IDs (%0.02f seconds)" % (process_time() - t0)
            )
        t0 = process_time()
        pitnbs = len(pitNb) + 1
        spillIDs, lspill, rank = spill_pts(
            MPIrank, pitnbs, self.lFill, self.pitIDs, self.borders[:, 1]
        )
        if MPIrank == 0 and self.verbose:
            print(
                "Define spill points (%0.02f seconds)" % (process_time() - t0)
            )
        t0 = process_time()
        self.lspillIDs = np.where(lspill == 1)[0]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, rank, op=MPI.MAX)
        spillIDs[rank != MPIrank] = -1
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, spillIDs, op=MPI.MAX)
        if MPIrank == 0 and self.verbose:
            print(
                "Define spill points (%0.02f seconds)" % (process_time() - t0)
            )
        t0 = process_time()

        # Get depression information
        self.pitInfo = np.zeros((pitnbs, 2), dtype=int)
        self.pitInfo[:, 0] = spillIDs
        self.pitInfo[:, 1] = rank

        # Transfer depression IDs along local borders
        self.tmpL.setArray(self.pitIDs)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        self.pitIDs = self.tmpL.getArray().copy().astype(int)
        self.pitIDs[self.outletIDs] = -1   # draining borders aren't pits; walls can be

        # Get directions on flats
        self._dirFlats()
        id = self.lFill <= self.sealevel

        h = hl.copy()
        self.lFill[id] = hl[id]
        self.flatDirs[id] = -1
        # Only compute the water volume for incoming water fluxes above sea level
        if sed:
            self.lFill[id] = hl[id]
            self.flatDirs[id] = -1
        else:
            h[h < self.sealevel] = self.sealevel
        if MPIrank == 0 and self.verbose:
            print(
                "Define flat directions (%0.02f seconds)" % (process_time() - t0)
            )
        t0 = process_time()

        # Get pit parameters
        self._getPitParams(h, pitnbs)

        if MPIrank == 0 and self.verbose:
            print(
                "Define depressions parameters (%0.02f seconds)" % (process_time() - t0)
            )

        return

    def fillElevation(self, sed=False):
        """
        This function is the main entry point to perform pit filling.

        It relies on the following private functions:

        - _performFilling
        - _pitInformation

        :arg sed: boolean specifying if the pits are filled with water or sediments.
        """

        tfill = process_time()

        hl = self.hLocal.getArray().copy()
        minh = self.hGlobal.min()[1]
        if not self.flatModel:
            minh += 1.0e-3  # TODO-REFACTOR: value matches DEPOSIT_FLOOR but distinct role (minh epsilon nudge for fillElevation); do not replace
        level = max(minh, self.sealevel + self.oFill)

        self._performFilling(hl - level, level, sed)
        self.lFill += level
        self._pitInformation(hl, level, sed)

        # Compute discrete depth-volume curves per pit. Used by water filling
        # (sed=False) for routing residual flux, and by the bottom-up
        # sediment fill in _updateSinks (sed=True) to find the lake-surface
        # level that matches the deposited volume.
        ids = self.pitParams[:, 0] > 0.0
        dh = np.zeros((len(self.pitInfo), 6), dtype=np.float64)
        dh[ids, 0] = self.pitParams[ids, 1] - self.pitParams[ids, 2]
        dh[ids, 1:] = np.expand_dims(self.pitParams[ids, 2] / 5.0, axis=1)
        self.filled_lvl = np.cumsum(dh, axis=1)[:, 1:]

        self.filled_vol = np.zeros((len(self.pitInfo), 5), dtype=np.float64)
        hl_clip = hl.copy()
        if not sed:
            hl_clip[hl_clip < self.sealevel] = self.sealevel
        self.filled_vol[:, :-1] = getpitvol(
            self.filled_lvl[:, :-1], hl_clip, self.pitIDs, self.inIDs
        )
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, self.filled_vol, op=MPI.SUM)
        self.filled_vol[:, -1] = self.pitParams[:, 0]

        if MPIrank == 0 and self.verbose:
            print(
                "Handling depressions over the surface (%0.02f seconds)"
                % (process_time() - tfill)
            )

        if self.memclear:
            del hl, minh
            gc.collect()

        return

    def fillIceElevation(self, hl):
        """
        This function is the main entry point to perform pit filling for glacier.

        It relies on the following private functions:

        - _performFilling
        - _pitInformation

        :arg sed: boolean specifying if the pits are filled with water or sediments.
        """

        tfill = process_time()
        minh = self.hGlobal.min()[1]
        if not self.flatModel:
            minh += 1.0e-3  # TODO-REFACTOR: value matches DEPOSIT_FLOOR but distinct role (minh epsilon nudge for fillIceElevation); do not replace
        level = max(minh, self.sealevel + self.oFill)

        self._performFilling(hl - level, level, False)
        self.lFill += level
        self._pitInformation(hl, level, False)

        # Define specific filling levels for unfilled water depressions
        # if not False:
        ids = self.pitParams[:, 0] > 0.0
        dh = np.zeros((len(self.pitInfo), 6), dtype=np.float64)
        dh[ids, 0] = self.pitParams[ids, 1] - self.pitParams[ids, 2]
        dh[ids, 1:] = np.expand_dims(self.pitParams[ids, 2] / 5.0, axis=1)
        self.filled_lvl = np.cumsum(dh, axis=1)[:, 1:]

        self.filled_vol = np.zeros((len(self.pitInfo), 5), dtype=np.float64)
        hl[hl < self.sealevel] = self.sealevel
        self.filled_vol[:, :-1] = getpitvol(
            self.filled_lvl[:, :-1], hl, self.pitIDs, self.inIDs
        )
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, self.filled_vol, op=MPI.SUM)
        self.filled_vol[:, -1] = self.pitParams[:, 0]

        if MPIrank == 0 and self.verbose:
            print(
                "Handling ice depressions over the smooth surface (%0.02f seconds)"
                % (process_time() - tfill)
            )

        if self.memclear:
            del minh
            gc.collect()

        return