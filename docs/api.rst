##################
API Documentation
##################

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
  :target: https://github.com/psf/black
  :alt: Code style: black

.. image:: https://img.shields.io/lgtm/alerts/g/Geodels/gospl.svg?logo=lgtm&logoWidth=18
  :target: https://lgtm.com/projects/g/Geodels/gospl/alerts/
  :alt: Total alerts

.. image:: https://img.shields.io/lgtm/grade/python/g/Geodels/gospl.svg?logo=lgtm&logoWidth=18
  :target: https://lgtm.com/projects/g/Geodels/gospl/context:python
  :alt: Language grade: Python


.. warning::
  It is worth mentioning that a set of **fortran** functions are part of *gospl* code but are not described in this API...


Model class
------------

.. automodule:: model
    :members:

Flow class
------------

.. automodule:: flow
    :members:

Flow path & erosion
^^^^^^^^^^^^^^^^^^^^^

.. automodule:: flow.flowplex
    :members:

Sediment class
---------------

.. automodule:: sed
    :members:

Continental & marine sedimentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: sed.sedplex
    :members:


Meshing class
--------------

.. automodule:: mesher
    :members:

External forcing over unstructured mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mesher.unstructuredmesh
    :members:


I/O classes
------------

.. automodule:: tools
    :members:

Input parser
^^^^^^^^^^^^^^^^^

.. automodule:: tools.inputparser
    :members:

Output mesh
^^^^^^^^^^^^^^^^^

.. automodule:: tools.outmesh
    :members:
