import os
import gc
import sys
import petsc4py
import numpy as np
import pandas as pd
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import edge_tile
    from gospl._fortran import fill_tile
    from gospl._fortran import fill_edges
    from gospl._fortran import fill_dir
    from gospl._fortran import nghb_dir
    from gospl._fortran import fill_rcvs
    from gospl._fortran import getpitvol
    from gospl._fortran import pits_cons
    from gospl._fortran import fill_depressions
    from gospl._fortran import graph_nodes
    from gospl._fortran import combine_edges
    from gospl._fortran import label_pits
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

    The approach proposed here is more general than the one in the initial paper. First, it handles both regular and irregular meshes, allowing for complex distributed meshes to be used as long as a clear definition of inter-mesh connectivities is available. Secondly, to prevent iteration over unnecessary vertices (such as marine regions), it is possible to define a minimal elevation (i.e. sea-level position) above which the algorithm is performed. Finally, it creates directions over flat regions allowing for downstream flows in cases where the entire volume of a depression is filled.

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
        edges = -np.ones((self.lpoints, 2), dtype=int)
        edges[self.idLBounds, 0] = self.idLBounds
        edges[self.idBorders, 0] = self.idBorders
        edges[self.idLBounds, 1] = 0
        edges[self.idBorders, 1] = 1
        self.borders = edges

        self.outEdges = np.zeros(self.lpoints, dtype=int)
        self.outEdges[self.ghostIDs] = 1

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

    def _sortingPits(self, df):
        """
        Sorts depressions number before combining them to ensure no depression index is
        changed in an unsorted way.

        :arg df: pandas dataframe containing depression numbers which have to be combined.

        :return: df sorted pandas dataframe containing depression numbers.
        """

        p1 = []
        p2 = []
        for k in range(len(df)):
            id1 = df["p1"].iloc[k]
            if k == 0:
                id2 = df["p2"].iloc[0]
            else:
                if df["p2"].iloc[k] == df["p2"].iloc[k - 1]:
                    id2 = df["p1"].iloc[k - 1]
                else:
                    id2 = df["p2"].iloc[k]
            p1.append(id1)
            p2.append(id2)
        data = {
            "p1": p1,
            "p2": p2,
        }
        df = pd.DataFrame(data, columns=["p1", "p2"])
        df = df.drop_duplicates().sort_values(["p2", "p1"], ascending=(False, False))

        return df

    def _offsetGlobal(self, lgth):
        """
        Computes the offset between processors to ensure a unique number for considered indices.

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

    def _transferIDs(self, pitIDs):
        """
        This function transfers local depression IDs along local borders and combines them with a unique identifier.

        :arg pitIDs: local depression index.

        :return: number of depressions.
        """

        # Define globally unique watershed index
        fillIDs = pitIDs >= 0
        offset, _ = self._offsetGlobal(np.amax(pitIDs))
        pitIDs[fillIDs] += offset[MPIrank]

        # Transfer depression IDs along local borders
        self.tmpL.setArray(pitIDs)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        label = self.tmpL.getArray().copy().astype(int)

        ids = label < pitIDs
        df = self._buildPitDataframe(label[ids], pitIDs[ids])
        ids = label > pitIDs
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

        # Sorting label transfer between processors
        sorting = True
        while sorting:
            df2 = self._sortingPits(df)
            sorting = not df.equals(df2)
            df = df2.copy()
        for k in range(len(df)):
            label[label == df["p2"].iloc[k]] = df["p1"].iloc[k]

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

        # Transfer filled values along the local borders
        keep = -np.ones(1)
        stp = 0
        while keep[0] < 0:
            self.tmp.setArray(pdir[self.glIDs])
            self.dm.globalToLocal(self.tmp, self.tmpL, 3)
            pdir = self.tmpL.getArray().copy().astype(int)
            pdir = nghb_dir(self.pitIDs, self.lFill, pdir)

            # Check if there are pit nodes without directions
            if ids.any():
                keep[0] = pdir[ids].min()
            else:
                keep[0] = 1.0
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, keep, op=MPI.MIN)
            stp += 1
            if stp > 100:
                break

        # Define receiver nodes on each depression
        self.flatDirs = fill_rcvs(self.pitIDs, self.lFill, pdir)
        self.flatDirs[self.lspillIDs] = -1

        return

    def _getPitParams(self, hl, nbpits):
        """
        Define depression global parameters:

        - volume of each depression
        - maximum filled depth

        :arg hl: numpy array of unfilled surface elevation
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
        hmax = -np.ones(nbpits, dtype=np.float64) * 1.0e8
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
        This function extracts depression informations available to all processors. It stores the following things:

        - the information for each depression (e.g., a unique global ID, its spillover local points and related processor),
        - the description of each depression (total volume and maximum filled depth).

        :arg hl: local elevation.
        :arg sed: boolean specifying if the pits are filled with water or sediments.
        """

        t0 = process_time()

        # Combine pits locally to get an unique local ID per depression
        if sed:
            pitIDs = label_pits(level, self.lFill)
        else:
            pitIDs = label_pits(self.sealevel, self.lFill)
        pitIDs[self.idBorders] = -1
        pitNb = self._transferIDs(pitIDs)
        pitnbs = len(pitNb) + 1
        spillIDs, lspill, rank = spill_pts(
            MPIrank, pitnbs, self.lFill, self.pitIDs, self.borders[:, 1]
        )
        self.lspillIDs = np.where(lspill == 1)[0]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, rank, op=MPI.MAX)
        spillIDs[rank != MPIrank] = -1
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, spillIDs, op=MPI.MAX)

        # Get depression informations:
        self.pitInfo = np.zeros((pitnbs, 2), dtype=int)
        self.pitInfo[:, 0] = spillIDs
        self.pitInfo[:, 1] = rank

        # Transfer depression IDs along local borders
        self.tmpL.setArray(self.pitIDs)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        self.pitIDs = self.tmpL.getArray().copy().astype(int)
        self.pitIDs[self.idBorders] = -1

        # Get directions on flats
        self._dirFlats()
        id = self.lFill <= self.sealevel

        h = hl.copy()
        if sed:
            self.oceanFill = self.lFill.copy()
            self.lFill[id] = hl[id]
            self.flatOcean = self.flatDirs.copy()
            self.flatDirs[id] = -1
        else:
            self.lFill[id] = hl[id]
            self.flatDirs[id] = -1
            # Only compute the water volume for incoming water fluxes above sea level
            h[h < self.sealevel] = self.sealevel

        # Get pit parameters
        self._getPitParams(h, pitnbs)

        if MPIrank == 0 and self.verbose:
            print(
                "Define depressions parameters (%0.02f seconds)" % (process_time() - t0)
            )

        return

    def fillElevation(self, sed=False):
        """
        This functions is the main entry point to perform pit filling.

        It relies on the following private functions:

        - _performFilling
        - _pitInformation

        :arg sed: boolean specifying if the pits are filled with water or sediments.
        """

        tfill = process_time()

        hl = self.hLocal.getArray().copy()
        minh = self.hGlobal.min()[1]
        if not self.flatModel:
            minh += 1.0e-3
        level = max(minh, self.sealevel - 6000.0)

        self._performFilling(hl - level, level, sed)
        self.lFill += level
        self._pitInformation(hl, level, sed)

        # Define specific filling levels for unfilled water depressions
        if not sed:
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
                "Handling depressions over the surface (%0.02f seconds)"
                % (process_time() - tfill)
            )

        if self.memclear:
            del hl, minh
            gc.collect()

        return
