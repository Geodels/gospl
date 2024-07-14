.. _hill:

==============================================
Hillslope and marine deposition
==============================================

Hillslope: soil creep
-----------------------

Hillslope processes are known to strongly influence catchment morphology and drainage density and several formulations of hillslope transport laws have been proposed. Most of these formulations are based on a mass conservation equation and assume that a layer of soil available for transport is always present (*i.e.* precluding  case of bare exposed bedrock) and that dissolution and mass transport in solution can be neglected.

Under such assumptions, the mass conservation equation widely applied in landscape modelling is of the form:

.. math::

  \mathrm{\frac{\partial \eta}{\partial t}} = -\mathrm{\nabla \cdot {q_{ds}}}

where :math:`\mathrm{q_{ds}}` is the volumetric soil flux of transportable sediment per unit width of the land surface. In its simplest form, :math:`\mathrm{q_{ds}}` obeys the Culling model and hypothesises a proportional relationship to local hillslope gradient (*i.e.* :math:`\mathrm{q_{ds}=-D\nabla \eta}`, also referred to as the **creep diffusion equation**):

.. math::

  \mathrm{\frac{\partial \eta}{\partial t}} = \mathrm{D \nabla^2 \eta}


in which :math:`\mathrm{D}` is the diffusion coefficient that encapsulates a variety of processes operating on the superficial soil layer. As an example, :math:`\mathrm{D}` may vary as a function of substrate, lithology, soil depth, climate and biological activity.

In goSPL, hillslope processes rely on this approximation even though field evidence suggest that the creep approximation is only rarely appropriate.

For a discrete element, considering a node :math:`\mathrm{i}` the implicit finite volume representation of the above equation is:

.. math::

  \mathrm{\frac{\partial \eta_i}{\partial t}} = \mathrm{\frac{\eta_i^{t+\Delta t}-\eta_i^t}{\Delta t} = D \sum_{j=1}^N \frac{  \chi_{i,j}(\eta_j^{t+\Delta t} - \eta_i^{t+\Delta t}) }{\Omega_i \lambda_{i,j}} }


:math:`\mathrm{N}` is the number of neighbours surrounding node :math:`\mathrm{i}`, :math:`\mathrm{\Omega_i}` is the voronoi area,  :math:`\mathrm{\lambda_{i,j}}` is the length of the edge connecting the considered nodes and :math:`\mathrm{\chi_{i,j}}` is the length of voronoi face shared by nodes :math:`\mathrm{i}` and :math:`\mathrm{j}`.

Applied to the entire domain, the equation above can be rewritten as a matrix system:

.. math::

  \mathrm{\mathbf Q \boldsymbol\eta^{t+\Delta t}} = \mathrm{\boldsymbol\eta^{t}}

where :math:`\mathrm{\mathbf Q}` is sparse. The matrix terms  only depend on the diffusion coefficient :math:`\mathrm{D}`, the grid parameters and voronoi variables (:math:`\mathrm{\chi_{i,j}}`,  :math:`\mathrm{\lambda_{i,j}}`, :math:`\mathrm{\Omega_i}`).

In goSPL, these parameters remain fixed  during a model run and therefore :math:`\mathrm{\mathbf Q}` needs to be created only once at initialisation. At each iteration, hillslope induced changes in elevation :math:`\mathrm{\boldsymbol \eta}` are then obtained in a similar way as for the solution of the other systems using `PETSc <https://www.mcs.anl.gov/petsc/>`_ *Richardson solver* and *block Jacobi* preconditioning.

Marine deposition
--------------------

In the marine realm, sediment transport is modelled through nonlinear diffusion equation. The rate of elevation change in the marine environment is governed by:


.. math::

  \mathrm{\frac{\partial \eta}{\partial t}} = \mathrm{\nabla \cdot \left( K_m(\eta) \nabla \eta \right)} + Q_{sr}
  

where :math:`\mathrm{K_m}` is the marine sediment transport coefficient (m2/yr), and :math:`\mathrm{Q_{sr}}` is the sediment flux coming at the river mouth. As the model progresses over time so does the shoreline position due to both offshore sedimentation and prescribed eustatic sea-level variations.

During a single time step, marine sediment transport is performed until all available sediment transported by rivers to the ocean have been diffused and the accumulations on these specific nodes remains below water depth or below a prescribed slope (i.e., clinoform slope) computed based on local water depth and distance to the nearest coastline.

The implicit finite volume approach is implemented using `PETSc <https://www.mcs.anl.gov/petsc/>`_ SNES functionality. The nonlinear system at each time step (``arkimex``) is solved iteratively with time stepping and the SNES solution is based on a Nonlinear Generalized Minimum Residual method (``NGMRES``) and the linear solver uses a ``richardson`` KSP with a ``bjacobi`` preconditioner.  (See PETSC documentation for more details about the solver and preconditoner options and settings). 
