.. _examples:

========
Examples
========

Ready-to-run models from the :examples:`goSPL-examples <.>`
repository, from small regional
set-ups to global runs. Each page gives the example's description (what it
shows, the set-up block by block, how to run it) and its results notebook
(``view_Results.ipynb``) as it was last run, with its figures.

To run an example yourself, clone the repository, open the example's folder
and follow its description: build the inputs (``scripts/make_inputs.py`` or
the setup notebook), run the model with ``gospl -i input.yml`` (or the
folder's ``run_all.sh``), then open ``view_Results.ipynb``. The examples use
the ``gospl-smoke`` environment of the repository (``environment.yml``), also
available as a Docker image (see :ref:`getting_started`).

.. include:: generated/gallery.rst
