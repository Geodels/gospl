:notoc:

.. gospl documentation master file, created by

.. module:: gospl

********************
gospl documentation
********************

**Date**: |today| **Version**: |version|


**Useful links**:
`Binary Installers <https://pypi.org/project/gospl>`__ |
`Source Repository <https://github.com/Geodels/gospl>`__ |
`Issues & Ideas <https://github.com/Geodels/gospl/issues>`__ |
`Q&A Support <https://stackoverflow.com/questions/tagged/gospl>`__ |
`Mailing List <https://groups.google.com/forum/#!forum/gospl>`__


:mod:`gospl` is an open source, GPL-licensed library providing a scalable parallelised Python-based numerical model to simulate landscapes and basins reconstruction at global scale.


.. raw:: html

    <div class="shadow gs-callout gs-callout-need">
        <h4>Statement of need</h4>

Since the '90s, many software have been designed to estimate long-term catchment dynamic, drainage evolution as well as sedimentary basins formation  in response to various mechanisms such as tectonic or climatic forcing. These models rely on a set of mathematical and physical expressions that simulate sediment erosion, transport and deposition and can reproduce the first order complexity of Earth's surface geomorphological evolution.

Yet, we were still missing a tool to evaluate **global scale** evolution of Earth surface and its interaction with the atmosphere, the hydrosphere, the tectonic and mantle dynamics. :mod:`gospl` is the first model designed to address this gap. It can be used to better characterise many aspects of the Earth system ranging from the role of atmospheric circulation on physical denudation, from the influence of erosion and deposition of sediments on mantle convection, from the location and abundance of natural resources to the evolution of life.

.. raw:: html

    </div>


.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="_static/index_getting_started.svg" class="card-img-top" alt="getting started with gospl action icon" height="52">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Getting started</h5>
                    <p class="card-text">New to <code>gospl</code>? Check out the getting started guides. They
                    contain an introduction to the code installation.</p>

.. container:: custom-button

    :ref:`To the getting started guides<getting_started>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="_static/index_user_guide.svg" class="card-img-top" alt="gospl physics action icon" height="52">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Technical guide</h5>
                    <p class="card-text">The technical guide provides in-depth information on the
                    underlying physics of <code>gospl</code>.</p>

.. container:: custom-button

    :ref:`To the technical guide<tech_guide>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="_static/index_contribute.svg" class="card-img-top" alt="running gospl action icon" height="52">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Hands-on guide</h5>
                    <p class="card-text">Learning how to use <code>gospl</code> by running some pre- and post processing examples available as
                    <a href="https://jupyter.org">Jupyter notebooks</a>.</p>

.. container:: custom-button

    :ref:`To the examples<example>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="_static/index_api.svg" class="card-img-top" alt="api of gospl action icon" height="52">
                <div class="card-body flex-fill">
                    <h5 class="card-title">API reference</h5>
                    <p class="card-text">This guide contains a detailed description of
                    <code>gospl</code> API. It describes how methods work and functions have
                    been declared. </p>

.. container:: custom-button

    :ref:`To the reference guide<api_ref>`

.. raw:: html

                </div>
                </div>
            </div>
        </div>
    </div>


.. toctree::
    :maxdepth: 3
    :hidden:
    :titlesonly:

    getting_started/index
    tech_guide/index
    examples.rst
    api_ref/index
