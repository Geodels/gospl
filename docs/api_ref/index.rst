.. _api_ref:

==============
API reference
==============

This page gives an overview of all :mod:`gospl` objects, functions and methods.

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
  :target: https://github.com/psf/black
  :alt: Code style: black

.. image:: https://img.shields.io/lgtm/alerts/g/Geodels/gospl.svg?logo=lgtm&logoWidth=18
  :target: https://lgtm.com/projects/g/Geodels/gospl/alerts/
  :alt: Total alerts

.. image:: https://img.shields.io/lgtm/grade/python/g/Geodels/gospl.svg?logo=lgtm&logoWidth=18
  :target: https://lgtm.com/projects/g/Geodels/gospl/context:python
  :alt: Language grade: Python

Main classes
--------------

.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Class Model
                </div>
                <div class="card-body">
                  <p class="card-text"> Instantiates model <code>object</code> and performs <code>surface processes</code> evolution. </p>

.. container:: custom-button

    :ref:`See functions and source code <model_ref>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Class Mesh
                </div>
                <div class="card-body">
                  <p class="card-text"> Defines spherical mesh characteristics and builds <code>PETSc DMPlex</code>. </p>


.. container:: custom-button

    :ref:`See functions and source code <mesh_ref>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Class Flow
                </div>
                <div class="card-body">
                  <p class="card-text"> <code>Flow accumulation</code> computation for global unstructured mesh. </p>

.. container:: custom-button

    :ref:`See functions and source code <flow_ref>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Class Sediment
                </div>
                <div class="card-body">
                  <p class="card-text"> Functions related to <code>sediment transport</code> and <code>deposition</code>. </p>


.. container:: custom-button

    :ref:`See functions and source code <sed_ref>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Class Stratigraphy
                </div>
                <div class="card-body">
                  <p class="card-text"> Functions related to <code>stratigraphic</code> architecture and <code>compaction</code>. </p>


.. container:: custom-button

    :ref:`See functions and source code <stra_ref>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-12 d-flex install-block">
                <div class="card install-card shadow w-100">
                <div class="card-header">
                    Class Input & Class Output
                </div>
                <div class="card-body">
                  <p class="card-text"> Input/output methods declaration. </p>


.. container:: custom-button

    See :ref:`input <in_ref>` & :ref:`output<out_ref>` functions

.. raw:: html

                </div>
                </div>
            </div>
        </div>
    </div>

.. toctree::
    :maxdepth: 3
    :hidden:

    model_ref
    mesh_ref
    flow_ref
    sed_ref
    stra_ref
    in_ref
    out_ref
