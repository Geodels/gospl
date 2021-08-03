import os
import gc
import sys
import warnings
import petsc4py
import numpy as np
import pandas as pd
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import fillpit
    from gospl._fortran import edge_tile
    from gospl._fortran import fill_tile
    from gospl._fortran import fill_edges
    from gospl._fortran import fill_depressions
    from gospl._fortran import graph_nodes
    from gospl._fortran import combine_edges
    from gospl._fortran import label_pits
    from gospl._fortran import pit_nodes
    from gospl._fortran import spill_pts

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()


class PITFill(object):
    """
    Depression filling is an important preconditioning step to many landscape evolution models.

    This class implements a linearly-scaling parallel priority-flood depression-filling algorithm based on `Barnes (2016) <https://arxiv.org/pdf/1606.06204.pdf>`_ algorithm.

    .. note::

        Unlike previous algorithms, `Barnes (2016) <https://arxiv.org/pdf/1606.06204.pdf>`_ approach guarantees a fixed number of memory access and communication events per processors. As mentionned in his paper based on comparison testing, it runs generally faster while using fewer resources than previous methods.

    The approach proposed here is more general than the one in the initial paper. First, it handles both regular and irregular meshes, allowing for complex distributed meshes to be  used as long as a clear definition of inter-mesh connectivities is available. Secondly, to prevent iteration over unnecessary vertices (such as marine regions), it is possible to define a minimal elevation (i.e. sea-level position) above which the algorithm is performed.

    For inter-mesh connections and message passing, the approach relies on PETSc DMPlex functions.

    The main functions return the following parameters:

    - the elevation of the filled surface,
    - the information for each depression (e.g., a unique global ID, its spillover local points and related processor),
    - the description of each depression (total volume and maximum filled depth).

    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `PITFill` class consists in the declaration of PETSc vectors, matrices and each partition internals edge vertices.
        """

        # Petsc vectors
        self.fZg = self.dm.createGlobalVector()  # fill elevation global
        self.fZl = self.dm.createLocalVector()  # fill elevation local
        self.lbg = self.dm.createGlobalVector()  # watershed label global
        self.lbl = self.dm.createLocalVector()  # watershed label local

        edges = -np.ones((self.lpoints, 2), dtype=int)
        edges[self.idLBounds, 0] = self.idLBounds
        edges[self.idBorders, 0] = self.idBorders
        edges[self.idLBounds, 1] = 0
        edges[self.idBorders, 1] = 1
        self.borders = edges
        out = np.where(edges[:, 1] > -1)[0]
        self.localEdges = edges[out, :]

        self.outEdges = np.zeros(self.lpoints, dtype=int)
        self.outEdges[self.shadow_lOuts] = 1

        self.rankIDs = None

        return

    def _buildPitDataframe(self, label1, label2):
        """
        Definition of a Pandas data frame used to find a unique pit ID between processors.

        :arg label1: depression ID in a given processors
        :arg label2: same depression ID in a neighbouring mesh

        :return: df (sorted dataframe of pit ID between processors)
        """

        data = {
            "p1": label1,
            "p2": label2,
        }
        df = pd.DataFrame(data, columns=["p1", "p2"])
        df = df.drop_duplicates().sort_values(["p2", "p1"], ascending=(False, False))

        return df

    def _offsetGlobal(self, lgth):
        """
        Compute the offset between processors to ensure a unique number for considered indices.

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

        :arg mgraph: numpy array containing local mesh edges information

        :arg ggraph: numpy array containing filled elevation values based on other processors values
        """

        # Get bidirectional edges connections
        cgraph = pd.DataFrame(
            mgraph, columns=["source", "target", "elev", "spill", "rank"]
        )
        cgraph = cgraph.sort_values("elev")
        cgraph = cgraph.drop_duplicates(["source", "target"], keep="first")
        c12 = np.concatenate((cgraph["source"].values, cgraph["target"].values))
        cmax = np.max(np.bincount(c12.astype(int))) + 1

        # Filling the bidirectional graph
        cgraph = cgraph.values
        elev, rank, nodes, spillID = fill_edges(
            int(max(cgraph[:, 1]) + 2), cgraph, cmax
        )
        ggraph = -np.ones((len(elev), 5))
        ggraph[:, 0] = elev
        ggraph[:, 1] = nodes
        ggraph[:, 2] = rank
        ggraph[:, 3] = spillID
        if self.memclear:
            del elev, nodes, rank
            del spillID, c12

        return ggraph

    def _getPitParams(self, lFill, hl, nbpits):
        """
        Define depression global parameters:

        - volume of each depression
        - maximum filled depth

        :arg lFill: numpy array of filled elevation
        :arg hl: numpy array of unfilled surface elevation
        :arg nbpits: number of depression in the global mesh
        """

        # Get pit parameters (volume and maximum filled elevation)
        ids = np.where(self.inIDs == 1)[0]
        grp = npi.group_by(self.pitIDs[ids])
        uids = grp.unique
        _, vol = grp.sum((lFill[ids] - hl[ids]) * self.larea[ids])
        _, height = grp.max(lFill[ids] - hl[ids])
        totv = np.zeros(nbpits, dtype=np.float64)
        toth = np.zeros(nbpits, dtype=np.float64)
        for k in range(len(uids)):
            if uids[k] > -1:
                totv[uids[k]] = vol[k]
                toth[uids[k]] = height[k]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, totv, op=MPI.SUM)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, toth, op=MPI.MAX)

        self.pitParams = np.empty((len(totv[1:]), 2), dtype=np.float64)
        self.pitParams[:, 0] = totv[1:]
        self.pitParams[:, 1] = toth[1:]

        return

    def serialFilling(self, level):
        """
        This function performs the depression-filling in serial using the *priority-flood* algorithm proposed in `Barnes et al. (2014) <https://doi.org/10.1016/j.cageo.2013.04.024>`_.

        :arg level: minimal elevation above which the algorithm is performed.
        """

        if self.rankIDs is None:
            self.rankIDs = np.empty(self.mpoints, dtype=int)
            self.rankIDs[self.locIDs] = MPIrank
            self.rankIDs[self.outIDs] = -1
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, self.rankIDs, op=MPI.MAX)

        hl = self.hLocal.getArray().copy()
        gZ = np.empty(self.mpoints, dtype=np.float64)
        gZ[self.locIDs] = hl
        gZ[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)
        if self.flatModel:
            mingZ = gZ.min()

        # Perform pit filling on process rank 0
        if MPIrank == 0:
            if self.flatModel:
                if mingZ > level:
                    hFill, pits = self._fill2D(mingZ, gZ, 1.0e6)
                else:
                    hFill, pits = fillpit(level, gZ, 1.0e6)
            else:
                hFill, pits = fillpit(level, gZ, 1.0e6)

        else:
            hFill = np.empty(self.mpoints, dtype=np.float64)
            pits = np.empty((self.mpoints, 2), dtype=np.float64)
        hFill = MPI.COMM_WORLD.bcast(hFill, root=0)
        pits = MPI.COMM_WORLD.bcast(pits, root=0)

        self.pitIDs = pits[self.locIDs, 0]

        # Find spillover points
        upits, ids = np.unique(pits[:, 0], return_index=True)
        self.pitInfo = np.zeros((len(upits), 3), dtype=int)
        gpts = pits[ids, 1]
        pitrank = self.rankIDs[gpts]
        ids = np.in1d(self.locIDs, gpts).nonzero()
        loc = np.arange(self.lpoints, dtype=int)[ids]
        ids = np.where(self.inIDs[loc] == 1)[0]
        loc = loc[ids]
        gpts = -np.ones(len(upits), dtype=int)
        gpts[self.pitIDs[loc]] = loc
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gpts, op=MPI.MAX)

        # Define pit information
        lFill = hFill[self.locIDs]
        self.fZl.setArray(lFill)
        self.dm.localToGlobal(self.fZl, self.fZg)
        self.dm.globalToLocal(self.fZg, self.fZl)

        self.pitInfo[:, 0] = upits
        self.pitInfo[:, 1] = gpts
        self.pitInfo[:, 2] = pitrank
        id = np.where(upits > -1)[0]
        self.pitInfo = self.pitInfo[id, :]

        # Get pit parameters (volume and maximum filled elevation)
        self._getPitParams(lFill, hl, len(upits))

        if self.memclear:
            del hl, gZ, ids, gpts, loc, upits, pits, lFill
            gc.collect()

        return

    def parallelFilling(self, level):
        """
        This functions implements the linearly-scaling parallel priority-flood depression-filling algorithm from `Barnes (2016) <https://arxiv.org/pdf/1606.06204.pdf>`_ but adapted to unstructured meshes.

        :arg level: minimal elevation above which the algorithm is performed.
        """

        t0 = process_time()
        hl = self.hLocal.getArray().copy()

        # Get local meshes edges and communication nodes
        ledges = edge_tile(level, self.borders, hl)
        out = np.where(ledges >= 0)[0]
        localEdges = np.empty((len(out), 2), dtype=int)
        localEdges[:, 0] = np.arange(self.lpoints)[out].astype(int)
        localEdges[:, 1] = ledges[out]
        out = np.where(ledges == -2)[0]
        inIDs = self.inIDs.copy()
        inIDs[out] = 0
        outEdges = self.outEdges.copy()
        outEdges[out] = 0
        out = np.where((ledges == 0) & (ledges == 2))[0]
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
        self.lbl.setArray(label.astype(int))
        self.dm.localToGlobal(self.lbl, self.lbg)
        self.dm.globalToLocal(self.lbg, self.lbl)
        label = self.lbl.getArray()

        # Transfer filled values along the local borders
        self.fZl.setArray(lFill)
        self.dm.localToGlobal(self.fZl, self.fZg)
        self.dm.globalToLocal(self.fZg, self.fZl)
        lFill = self.fZl.getArray()

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
        offset, sum = self._offsetGlobal(lgrph)
        graph = -np.ones((sum, 5), dtype=float)
        if lgrph > 0:
            graph[offset[MPIrank] : offset[MPIrank] + lgrph, :4] = cgraph
            graph[offset[MPIrank] : offset[MPIrank] + lgrph, 4] = MPIrank

        # Build global spillover graph on master
        if MPIrank == 0:
            mgraph = -np.ones((sum, 5), dtype=float)
        else:
            mgraph = None
        MPI.COMM_WORLD.Reduce(graph, mgraph, op=MPI.MAX, root=0)
        if MPIrank == 0:
            ggraph = self._fillFromEdges(mgraph)
        else:
            ggraph = None

        # Send filled graph dataset to each processors
        graph = MPI.COMM_WORLD.bcast(ggraph, root=0)

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
        graph[graph[:, 0] < -1.0e8, 0] = -1.0e8
        graph[graph[:, 0] > 1.0e7, 0] = -1.0e8

        # Define global solution by combining depressions/flat together
        lFill = fill_depressions(level, hl, lFill, label.astype(int), graph[:, 0])

        self.fZl.setArray(lFill)
        self.dm.localToGlobal(self.fZl, self.fZg)
        self.dm.globalToLocal(self.fZg, self.fZl)

        if self.memclear:
            del hl, label, gnb, graph, lFill
            del label_offset, offset, proc, keep
            del cgraph, outs, mgraph, ggraph
            gc.collect()

        if MPIrank == 0 and self.verbose:
            print("Compute pit filling (%0.02f seconds)" % (process_time() - t0))

        return

    def pitInformation(self):
        """
        This functions extracts depression informations available to all processors. It stores the following things:

        - the information for each depression (e.g., a unique global ID, its spillover local points and related processor),
        - the description of each depression (total volume and maximum filled depth).

        :arg level: minimal elevation above which the algorithm is performed.
        """

        t0 = process_time()

        # Get filled information
        hl = self.hLocal.getArray().copy()
        lFill = self.fZl.getArray().copy()
        fillIDs = np.where(lFill > hl)[0]

        # Combine pits locally to get a unique ID per depression
        t1 = process_time()
        pitIDs, pitnbs = label_pits(hl, lFill)
        pitArray = pit_nodes(pitnbs)
        df = self._buildPitDataframe(pitArray[:, 0], pitArray[:, 1])
        for k in range(len(df)):
            id2 = df["p2"].iloc[k]
            id1 = df["p1"].iloc[k]
            pitIDs[pitIDs == id2] = id1

        # Define globally unique watershed index
        offset, _ = self._offsetGlobal(np.amax(pitIDs))
        pitIDs[fillIDs] += offset[MPIrank]

        # Transfer depression IDs along local borders
        self.lbl.setArray(pitIDs)
        self.dm.localToGlobal(self.lbl, self.lbg)
        self.dm.globalToLocal(self.lbg, self.lbl)
        label = self.lbl.getArray()

        ids = np.where(label < pitIDs)[0]
        df = self._buildPitDataframe(label[ids], pitIDs[ids])
        ids = np.where(label > pitIDs)[0]
        df2 = self._buildPitDataframe(pitIDs[ids], label[ids])
        df = df.append(df2, ignore_index=True)
        df = df.drop_duplicates().sort_values(["p2", "p1"], ascending=(False, False))
        df = df[(df["p1"] >= 0) & (df["p2"] >= 0)]

        # Send depression IDs globally
        offset, _ = self._offsetGlobal(len(df))
        combIds = -np.ones((np.amax(offset), 2), dtype=int)
        if len(df) > 0:
            combIds[offset[MPIrank] : offset[MPIrank + 1], :] = df.values
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, combIds, op=MPI.MAX)
        df = self._buildPitDataframe(combIds[:, 0], combIds[:, 1])
        df = df[(df["p1"] >= 0) & (df["p2"] >= 0)]
        for k in range(len(df)):
            id1 = df["p1"].iloc[k]
            if k == 0:
                id2 = df["p2"].iloc[0]
            else:
                if df["p2"].iloc[k] == df["p2"].iloc[k - 1]:
                    id2 = df["p1"].iloc[k - 1]
                else:
                    id2 = df["p2"].iloc[k]
            label[label == id2] = id1

        # Transfer depression IDs along local borders
        self.lbl.setArray(label)
        self.dm.localToGlobal(self.lbl, self.lbg)
        self.dm.globalToLocal(self.lbg, self.lbl)

        # At this point all pits have a unique IDs across processors
        self.pitIDs = self.lbl.getArray().astype(int)
        if MPIrank == 0 and self.verbose:
            print(
                "Define unique depression IDs (%0.02f seconds)" % (process_time() - t1)
            )

        # Get spill over points
        t1 = process_time()
        pitnbs = np.zeros(1, dtype=int)
        pitnbs[0] = np.max(self.pitIDs)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, pitnbs, op=MPI.MAX)
        rank = -np.ones(pitnbs[0], dtype=int)
        spillIDs = -np.ones(pitnbs[0], dtype=np.int32)
        if len(fillIDs) > 0:
            spillIDs = spill_pts(
                pitnbs[0],
                fillIDs,
                self.pitIDs[fillIDs],
                lFill,
                self.borders[:, 1],
                self.inIDs,
            )
            id = np.where(spillIDs > -1)[0]
            rank[id] = MPIrank
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, rank, op=MPI.MAX)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, spillIDs, op=MPI.MAX)

        # Get depression informations:
        pitids = np.arange(pitnbs[0]) + 1
        self.pitInfo = np.zeros((pitnbs[0], 3), dtype=int)
        self.pitInfo[:, 0] = pitids
        self.pitInfo[:, 1] = spillIDs
        self.pitInfo[:, 2] = rank
        id = np.where(rank > -1)[0]
        self.pitInfo = self.pitInfo[id, :]

        # Get pit parameters
        self._getPitParams(lFill, hl, pitnbs[0] + 1)
        self.pitParams = self.pitParams[id, :]

        if MPIrank == 0 and self.verbose:
            print("Spill points (%0.02f seconds)" % (process_time() - t1))

        if self.memclear:
            del hl, label, df, df2, data, lFill
            del label_offset, offset, pitIDs
            del fillIDs, combIds, pitArray
            gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Define depression parameters (%0.02f seconds)" % (process_time() - t0)
            )

        return
