.. _flow:

=========================
River Discharge
=========================



SFD to MFD
-----------

.. image:: ../images/flowpath.png
   :align: center


The flow discharge at node i (:math:`\mathrm{q_i}`) is determined as follows:


.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-body">

.. math::

  \mathrm{q_i} = \mathrm{b_i} + \mathrm{\sum_{d=1}^{N_d} q_d}

.. raw:: html

                </div>
                </div>
            </div>
        </div>
    </div>


where :math:`\mathrm{b_i}` is the local volume of water :math:`\mathrm{\Omega_i P_i}` where :math:`\mathrm{\Omega_i}` is the voronoi area and :math:`\mathrm{P_i}` the local precipitation value available for runoff during a given time step. :math:`\mathrm{N_d}` is the number of donors with a donor defined as a node that drains into i (as an example the donor of vertex 5 in the SFD sketch in the above figure is 1). To find the donors of each node, the method consists in finding their receivers first. Then, the receivers of each donor is saved into a receiver matrix, noting that the nodes, which are local minima, are their own receivers.

Finally the transpose of the matrix is used to get the donor matrix. When the previous equation is applied to all nodes and considering the MFD case illustrated above, the following relations are obtained:



.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-body">

.. math::

  \mathrm{q_1} &= \mathrm{b_1} \\
  \mathrm{q_2} &= \mathrm{b_2 + q_1 w_{1,2}}  \\
  \mathrm{q_3} &= \mathrm{b_3 + q_2 w_{2,3} + q_4 w_{4,3} } \\
  \mathrm{q_4} &= \mathrm{b_4 +  q_1 w_{1,4} + q_2 w_{2,4}}  \\
  \mathrm{q_5} &= \mathrm{b_5 + q_1 w_{1,5} + q_4 w_{4,5}} \\
  \mathrm{q_6} &= \mathrm{b_6 + q_4 w_{4,6} + q_5 w_{5,6} + q_7 w_{7,6}}  \\
  \mathrm{q_7} &= \mathrm{b_7 + q_{10} w_{10,7}} \\
  \mathrm{q_8} &= \mathrm{b_8 + q_3 w_{3,8} + q_4 w_{4,8} + q_6 w_{6,8} + q_7 w_{7,8} + q_{10} w_{10,8}}\\
  \mathrm{q_9} &= \mathrm{b_9 + q_3 w_{3,9} + q_8 w_{8,9} + q_{10} w_{10,9}}

.. raw:: html

                </div>
                </div>
            </div>
        </div>
    </div>


Algorithm
---------------
