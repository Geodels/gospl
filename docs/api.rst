##################
API Documentation
##################

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
  :target: https://github.com/psf/black

  
**gLEM** is a global-scale surface evolution model built to simulate landscape development, sediment transport over Earth surface.


.. warning::
  It is worth mentioning that a set of **fortran** functions are part of **gLEM** code but are not described in this API...


Model Class
------------

.. automodule:: model
    :members:

Flow Model
------------

.. automodule:: flow
    :members:

Surface Processes
^^^^^^^^^^^^^^^^^^

.. automodule:: flow.surfprocplex
    :members:


Mesher Model
------------

.. automodule:: mesher
    :members:

Unstructured Mesh
^^^^^^^^^^^^^^^^^^^^

.. automodule:: mesher.unstructuredmesh
    :members:


Tools
---------

.. automodule:: tools
    :members:

Input Parser
^^^^^^^^^^^^^^^^^

.. automodule:: tools.inputparser
    :members:

Output Mesh
^^^^^^^^^^^^^^^^^

.. automodule:: tools.outmesh
    :members:
