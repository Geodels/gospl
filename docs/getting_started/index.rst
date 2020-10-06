.. _getting_started:

===============
Getting started
===============

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
                <img src="../_static/docker.svg" class="card-img-top" alt="Docker logo" height="42">
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
                <img src="../_static/Conda_logo.svg" class="card-img-top" alt="SQL logo" height="32">
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
                    <img src="../_static/PyPI_logo.svg" class="card-img-top" alt="STATA logo" height="62">
                    <div class="card-body flex-fill">
                        <p class="card-text">The <code>data set</code> included in the
                            <a href="https://en.wikipedia.org/wiki/Stata">STATA</a> statistical software suite corresponds
                            to the pandas <code>dataframe</code>. Many of the operations known from STATA have an equivalent
                            in pandas.</p>

.. container:: custom-button

    :ref:`Learn more <_installPypi>`

.. raw:: html

                    </div>
                    </div>
                </div>
                <div class="col-lg-6 col-md-6 col-sm-12 col-xs-12 d-flex">
                    <div class="card text-center intro-card shadow">
                    <img src="../_static/Python_logo.svg" class="card-img-top" alt="SAS logo" height="52">
                    <div class="card-body flex-fill">
                        <p class="card-text">The  <a href="https://en.wikipedia.org/wiki/SAS_(software)">SAS</a> statistical software suite
                            also provides the <code>data set</code> corresponding to the pandas <code>dataframe</code>.
                            Also SAS vectorized operations, filtering, string processing operations, and more have similar
                            functions in pandas.</p>

.. container:: custom-button

    :ref:`Learn more <_installSrc>`

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
