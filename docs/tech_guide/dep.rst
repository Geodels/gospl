.. _dep:

==================================
Inland depressions & deposition
==================================

Priority-flood
---------------------------------

In most landscape evolution models, internally-draining regions (*e.g.*, depressions and pits) are usually filled before the calculation of flow discharge and erosion-deposition rates. This ensures that all flows conveniently reach the coast or the boundary of the simulated domain. In models intended to simulate purely erosional features, such depressions are usually treated as transient features and often ignored.

.. note::

  However, goSPL is designed to not only address erosion problems but also to simulate source-to-sink transfer and sedimentary basins formation and evolution in potentially complex tectonic settings. In such cases, depressions may be formed at different periods during runtime and may be filled or remain internally drained (*e.g.*, endorheic basins) depending on the volume of sediment transported by upstream catchments.

Depression filling approaches have received some attention in recent years with the development of new and more efficient algorithms such as the work from `Barnes et al. (2016) <https://arxiv.org/pdf/1606.06204.pdf>`_. These methods based on **priority-flood**  offer a time complexity of the order of :math:`\mathrm{O(Nlog(N))}`.

Priority-flood algorithms consist in finding the minimum elevation a cell needs to be raised to (*e.g.*, spill elevation of a cell) to prevent downstream ascending path to occur. They rely on priority queue data structure used to efficiently find the lowest spill elevation in a grid.

.. figure:: ../images/flood.png
  :scale: 60 %
  :align: center

  The Priority-Flood begins by adding all of the edge cells to the priority queue. Queued cells are represented by a black circle. Each edge cell is the mouth of its own watershed, represented with different colours here. The queue's lowest cell c is dequeued and its neighbours added to the queue; the neighbours inherit c's watershed label. Depressions are filled in. When two different watersheds meet, the maximum elevation of the two meeting cells is noted: here there are five distinct elevation levels and the two watersheds meet at an elevation of 5. If this noted elevation is the lowest of any meeting of the two watersheds, it is retained as the watersheds' spillover elevation (adapted from Barnes et al. (2014)).


In goSPL, the priority-flood algorithm proposed in `Barnes et al. (2016) <https://arxiv.org/pdf/1606.06204.pdf>`_ is implemented. It provides a solution to remove automatically flat surfaces.

The approach proposed in goSPL is more general than the one in the initial paper. First, it handles both regular and irregular meshes, allowing for complex distributed meshes to be used as long as a clear definition of inter-mesh connectivities is available. Secondly, to prevent iteration over unnecessary vertices (such as marine regions), it is possible to define a minimal elevation (i.e. sea-level position) above which the algorithm is performed. Finally, it creates directions over flat regions allowing for downstream flows in cases where the entire volume of a depression is filled.

.. important::

  The priority-flood algorithm returns the **flooded elevation** and for each depression its **volume**, its **spill over node**  and a **basin unique identifier** assigned to each of the nodes belonging to the depression.

Depression filling
---------------------------------

Once the volumes of these depressions have been obtained, their subsequent filling is dependent of the sediment fluxes calculation defined in the :ref:`previous section <ero>` (Fig. 4a).

.. figure:: ../images/filldep.png
  :align: center

  Illustration of the two cases that may arise depending on the volume of sediment entering an internally drained depression (panel **a**). The red line shows the limit of the depression at the minimal spillover elevation. **b)** The volume of sediment (:math:`\mathrm{V_s^{in}}`) is lower than the depression volume :math:`\mathrm{V_{pit}}`. In this case all sediments are deposited and no additional calculation is required. **c)** If :math:`\mathrm{V_s^{in}\ge V_{pit}}`, the depression is filled up to depression filling elevation (priority-flood + :math:`\mathrm{\epsilon}`), the flow calculation needs to be recalculated and the excess sediment flux (:math:`\mathrm{Q_s^{ex}}`) is transported to downstream nodes.

In cases where the incoming sediment volume is lower than the depression volume (Fig. 4b), all sediments are deposited and the elevation at node :math:`i` in the depression is increased by a thickness :math:`\mathrm{\delta_i}` such that:

.. math::

    \mathrm{\delta_{i}} = \mathrm{\Upsilon (\eta^{f}_{i}-\eta_{i})}


where :math:`\mathrm{\eta^{f}_{i}}` is the filling elevation of node :math:`\mathrm{i}` obtained with the priority-flood + :math:`\mathrm{\epsilon}` algorithm and  the ratio :math:`\mathrm{\Upsilon}` is set to :math:`\mathrm{V_s^{in}/V_{pit}}`.

If the cumulative sediment volume  transported by the rivers draining in a specific depression is above the volume of the depression (:math:`\mathrm{V_s^{in} \ge V_{pit}}` - Fig. 4c) the elevation of each node :math:`\mathrm{i}` is increased to its filling elevation (:math:`\mathrm{\eta^{f}_{i}}`) and the excess sediment volume is allocated to the spillover node (Fig. 4c).

The updated elevation field is then used to compute the flow accumulation following the approach presented in section :ref:`1 <flow>` and :ref:`2 <ero>`. The sediment fluxes are initially set to zero except on the spillover nodes and the excess sediments are transported downstream.

During a specific time step, the processed described above is iteratively repeated until all sediments are deposited in inlands depressions or have entered the marine environment.
