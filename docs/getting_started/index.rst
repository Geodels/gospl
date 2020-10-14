.. _getting_started:

===============
Getting started
===============


Python version support
----------------------

Officially Python 3.7.1 and above, 3.8, and 3.9.

Installation instructions
-------------------------

The easiest way to start using :mod:`gospl` is to install it as part of the
`Docker <https://hub.docker.com/r/geodels/gospl>`__ distribution, a cross platform
distribution for data analysis and scientific computing. This is the recommended
installation method for most users.

Instructions for installing from PyPI, Anaconda, or source are also provided.


.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-lg-6 col-md-6 col-sm-12 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/docker.svg" class="card-img-top" alt="Docker logo" height="72">
                <div class="card-body flex-fill">
                    <p class="card-text">Lightweight virtualisation with <a href="https://hub.docker.com/r/geodels/gospl">Docker container</a> provides a simple approach for running <code>gospl</code> model.</p>

.. container:: custom-button

    :ref:`Learn more <installDocker>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-12 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/Conda_logo.svg" class="card-img-top" alt="Conda logo" height="42">
                <div class="card-body flex-fill">
                    <p class="card-text">For <code>Anaconda</code> users, an environment with specific collection of conda packages and <code>gospl</code> dependencies is available.</p>

.. container:: custom-button

    :ref:`Learn more <installConda>`

.. raw:: html

                    </div>
                    </div>
                </div>
                <div class="col-lg-6 col-md-6 col-sm-12 col-xs-12 d-flex">
                    <div class="card text-center intro-card shadow">
                    <img src="../_static/PyPI_logo.svg" class="card-img-top" alt="PyPI logo" height="72">
                    <div class="card-body flex-fill">
                        <p class="card-text">For Linux user we provide a step-by-step guide on how to install <code>gospl</code> and its dependencies using <code>PyPI</code>.</p>

.. container:: custom-button

    :ref:`Learn more <installPypi>`

.. raw:: html

                    </div>
                    </div>
                </div>
                <div class="col-lg-6 col-md-6 col-sm-12 col-xs-12 d-flex">
                    <div class="card text-center intro-card shadow">
                    <img src="../_static/Pythonlogo.svg" class="card-img-top" alt="Python logo" height="62">
                    <div class="card-body flex-fill">
                        <p class="card-text"><code>gospl</code> can also be installed directly from source but this requires installation of all other dependencies.</p>

.. container:: custom-button

    :ref:`Learn more <installSrc>`

.. raw:: html

                    </div>
                    </div>
                </div>
        </div>
    </div>


.. toctree::
    :maxdepth: 2
    :hidden:

    installDocker
    installConda
    installPypi
    installSrc
