.. _tech_guide:

================
Technical Guide
================

The Technical Guide covers the implicit, iterative approaches used to solve the multiple flow direction water routing and the erosion deposition processes main algorithms implemented in :mod:`gospl`.

:mod:`gospl` is mostly written in ``Python`` with some functions in ``Fortran`` and takes advantage of ``PETSc`` solvers over parallel computing architectures using ``MPI``.

Further information on any specific method can be obtained in the
:ref:`api_ref`.


.. raw:: html

    <div class="shadow gs-callout gs-callout-remember">
        <h4>Short Description</h4>

The code is primarily a **parallel global scale landscape evolution model**, built to simulate **topography and basins** dynamics. It has the ability to track **two types of clastic sediment size**. The following processes are considered:

- **river incision** using stream power law,
- inland **river deposition** in depressions,
- **marine deposition** at river mouth,
- **hillslope processes** in both marine and inland areas, and
- **sediment compaction** as stratigraphic layers geometry and properties change.

It can be forced with spatially and temporally varying tectonics (horizontal and vertical displacements) and climatic forces (temporal and spatial precipitation changes and sea-level fluctuations).

.. raw:: html

    </div>




.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    River flow implementation
                </div>
                <div class="card-body">
                    <p class="card-text">Based on a parallel implicit
                      drainage area (IDA) method.
                      Want to gain insights on the implemented approach?</p>

.. container:: custom-button

    :ref:`Learn more <flow>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Erosion rate and sediment flux
                </div>
                <div class="card-body">
                    <p class="card-text">Based on the stream power law
                      (SPL), river erosion depends on local slope,
                      discharge and erodibility coefficient.</p>

.. container:: custom-button

    :ref:`Learn more <ero>`

.. raw:: html

                </div>
                </div>
            </div>
        </div>
    </div>



.. toctree::
    :maxdepth: 3
    :hidden:

    flow
    ero
    dep
    mar
    hill
