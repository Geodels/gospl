.. _hill:

====================================
Hillslope processes and stratigraphy
====================================


Hillslope: soil creep
-----------------------

Hillslope processes are known to strongly influence catchment morphology and drainage density and several formulations of hillslope transport laws have been proposed. Most of these formulations are based on a mass conservation equation and assume that a layer of soil available for transport is always present (*i.e.* precluding  case of bare exposed bedrock) and that dissolution and mass transport in solution can be neglected.

Under such assumptions and via the Exner's law, the mass conservation equation widely applied in landscape modelling is of the form:

.. math::

  \mathrm{\frac{\partial \eta}{\partial t}} = -\mathrm{\nabla \cdot {q_{ds}}}

where :math:`\mathrm{q_{ds}}` is the volumetric soil flux of transportable sediment per unit width of the land surface. In its simplest form, :math:`\mathrm{q_{ds}}` obeys the Culling model and hypothesises a proportional relationship to local hillslope gradient (*i.e.* :math:`\mathrm{q_{ds}=-D\nabla \eta}`, also referred to as the **creep diffusion equation**):

.. math::

  \mathrm{\frac{\partial \eta}{\partial t}} = \mathrm{D \Delta \eta}


in which :math:`\mathrm{D}` is the diffusion coefficient that encapsulates a variety of processes operating on the superficial soil layer. As an example, :math:`\mathrm{D}` may vary as a function of substrate, lithology, soil depth, climate and biological activity.

In :mod:`gospl`, hillslope processes rely on this approximation even though field evidence suggest that the creep approximation is only rarely appropriate.

For a discrete element, considering a node :math:`\mathrm{i}` the implicit finite volume representation of  the above equation is:

.. math::

  \mathrm{\frac{\partial \eta_i}{\partial t}} = \mathrm{\frac{\eta_i^{t+\Delta t}-\eta_i^t}{\Delta t} = D \sum_{j=1}^N \frac{  \chi_{i,j}(\eta_j^{t+\Delta t} - \eta_i^{t+\Delta t}) }{\Omega_i \lambda_{i,j}} }


:math:`\mathrm{N}` is the number of neighbours surrounding node :math:`\mathrm{i}`, :math:`\mathrm{\Omega_i}` is the voronoi area,  :math:`\mathrm{\lambda_{i,j}}` is the length of the edge connecting the considered nodes and :math:`\mathrm{\chi_{i,j}}` is the length of voronoi face shared by nodes :math:`\mathrm{i}` and :math:`\mathrm{j}`.

Applied to the entire domain, the equation above can be rewritten as a matrix system:

.. math::

  \mathrm{\mathbf Q \boldsymbol\eta^{t+\Delta t}} = \mathrm{\boldsymbol\eta^{t}}

where :math:`\mathrm{\mathbf Q}` is sparse. The matrix terms  only depend on the diffusion coefficient :math:`\mathrm{D}`, the grid parameters and voronoi variables (:math:`\mathrm{\chi_{i,j}}`,  :math:`\mathrm{\lambda_{i,j}}`, :math:`\mathrm{\Omega_i}`).

In :mod:`gospl`, these parameters remain fixed  during a model run and therefore :math:`\mathrm{\mathbf Q}` needs to be created only once at initialisation. At each iteration, hillslope induced changes in elevation :math:`\mathrm{\boldsymbol \eta}` are then obtained in a similar way as for the solution of the other systems using `PETSc <https://www.mcs.anl.gov/petsc/>`_ *Richardson solver* and *block Jacobi* preconditioning.


Stratigraphy evolution
-------------------------


Dual-lithology cases
^^^^^^^^^^^^^^^^^^^^^



Porosity and compaction
^^^^^^^^^^^^^^^^^^^^^^^^
