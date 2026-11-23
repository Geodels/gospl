.. _hill:

==============================================
Hillslope and marine deposition
==============================================

Hillslope: soil creep
-----------------------

Linear creep law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


Non-linear creep laws
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A more complex version of the creep law involves a non-linear relationship between soil flux and topographic gradient. Several nonlinear creep-transport laws have been suggested in the literature. 

.. math::

  \frac{\partial \eta}{\partial t}= -\nabla \cdot \left( D(\eta) \nabla \eta \right)


and in goSPL, the user might choose between three laws: the truncated Taylor series below, the full critical-slope law it approximates (``hillslopeNb: 0``), and a power law. 


The first one follows is derived by Ganti et al. (2012) for the Andrews-Bucknam (see `Barnhart et al. (2019) <https://gmd.copernicus.org/articles/12/1267/2019/gmd-12-1267-2019.pdf>`_). Where the non-linear formulation captures accelerated creep and shallow landsliding as the gradient approaches an effective angle of repose for loose granular material. It uses a truncated Taylor series formulation for soil flux and takes the following form:

.. math::

  \mathrm{q_h} = \mathrm{-D S \left[  1 + \left(\frac{S}{S_c}\right)^2 + \left(\frac{S}{S_c}\right)^4 + ... + + \left(\frac{S}{S_c}\right)^{2(N-1)}  \right] }

where :math:`\mathrm{S = −\nabla \eta}` is topographic gradient (positive downhill), :math:`\mathrm{D}` is the transport efficiency factor, and :math:`\mathrm{S_c}` is a critical gradient. The user specifies the number of terms :math:`\mathrm{N}` to be used in the approximation. Setting :math:`\mathrm{N = 0}` uses the full critical-slope law instead, :math:`\mathrm{q_h = -D S / (1 - (S/S_c)^2)}` (Roering et al., 1999), with :math:`\mathrm{S/S_c}` capped at 0.8 so the flux stays finite. The nonlinear flux rule results in convex-up topography for shallow slopes and transitions to linear hillslopes for steeper slopes (see `Barnhart et al. (2019) <https://gmd.copernicus.org/articles/12/1267/2019/gmd-12-1267-2019.pdf>`_).


The second approach uses a non-critical hillslope model following the work of `Wang et al. (2024) <https://www.sciencedirect.com/science/article/pii/S0169555X24001053>`_. This hillslope model is compatible with the existing linear and non-linear hillslope models and is capable of emulating topographic evolution under arbitrary hillslope gradients. It takes the following form:

.. math::

  \mathrm{q_h} = \mathrm{-D S^n}

:math:`\mathrm{n}` is a positive constant of the hillslope gradient exponent. When :math:`\mathrm{n=1}`, the equation simplifies to a linear hillslope model; when :math:`\mathrm{n>1}`, it becomes a non-linear hillslope model. 

**Coupling with rivers and resolution.** Fluvial incision and sediment routing act first in each step, then the diffusion acts on every land cell, channel cells included: what diffusion brings into a valley is deposited in the channel cells and removed by the next fluvial step (as eroded surface in the detachment-limited mode, as part of the routed load with ``G > 0``). A river is a line of single cells, so the diffusive flux delivered into it grows as the cells shrink and the balance between the two processes depends on the mesh resolution (Hergarten and Pietrek, 2023): with ``hillslopeKa: 0.2`` the steady channel steepness of an uplifted block rose from 30.5 to 66.4 (instead of 25) between 200 and 50 m spacing, while without diffusion it stayed at 25. Setting ``channelQ`` (m³/yr) removes this dependence: cells whose discharge reaches it are rivers and are held at their elevation during the diffusion solve (identity rows, like the outlets), so hillslopes grade to them and their slope is set by the fluvial law alone. The volume the hillslopes deliver to each channel cell, :math:`\mathrm{d = h_{old} - (A h)}` with the unmasked operator row :math:`\mathrm{A}` (or :math:`\mathrm{\Delta t}` times the flux divergence for the non-linear laws), is added to the routed sediment source of the same step: with ``channelQ`` and deposition on, the hillslope step runs between the fluvial incision and the sediment routing, so the rivers carry the hillslope supply in the step it is produced. With stratigraphy the routed source is the deposit-equivalent thickness of the eroded material: the supply is converted with the ratio of deposit-equivalent to in-situ thickness of what its upslope neighbours eroded, and split by their fine fraction and source classes (flux-weighted by their creep into the channel cell); a channel cell that would lose material keeps that loss, so the volume budget is exact up to the solver tolerance. With the soil model, whose creep is integrated by an adaptive TS, the channel rows are held fixed in its residual and Jacobian; the volume the non-channel cells lose is exactly what reached the channels, and it is shared among them in proportion to their inflow at the end of the step. With ``channelQ: 1.e6`` the channel steepness of the same test is 25.1, 25.2 and 25.4. Without ``channelQ``, choose the land diffusivity for the mesh resolution (see the user guide, *Changing the mesh resolution*).

**Discretisation.** On each Voronoi face the flux is the linear finite-volume flux :math:`\mathrm{D\,\Delta\eta/L}` multiplied by the law's factor :math:`\mathrm{q_h/(D S)}`, evaluated on the face average of the gradient magnitude :math:`\mathrm{|\nabla\eta|}` of the two cells (a Green-Gauss gradient on the Voronoi cell). Using the slope along the face's own edge instead would under-estimate :math:`\mathrm{S}` on every edge not aligned with the gradient and leave a steady ridge a few percent too high at any resolution. The implicit step is solved with a PETSc SNES (quasi-Newton L-BFGS, options prefix ``hillnl_``). All three laws are validated against their exact steady ridge profiles in ``benchmarks/test_hillslope_nonlinear.py``.

.. note::
  
  Since they share a common mathematical form, this last formulation can vary smoothly between linear and non-linear behaviors. Depending on the specific transport mechanisms involved, n is suggested to vary approximately between 1. and 3. (`Wang et al. (2024) <https://www.sciencedirect.com/science/article/pii/S0169555X24001053>`_).


Marine deposition
--------------------

In the marine realm, sediment transport is modelled through nonlinear diffusion equation. The rate of elevation change in the marine environment is governed by:


.. math::

  \mathrm{\frac{\partial \eta}{\partial t}} = \mathrm{\nabla \cdot \left( K_m(\eta) \nabla \eta \right)} + Q_{sr}

where :math:`\mathrm{K_m}` is the marine sediment transport coefficient (m2/yr), and :math:`\mathrm{Q_{sr}}` is the sediment flux coming at the river mouth. As the model progresses over time so does the shoreline position due to both offshore sedimentation and prescribed eustatic sea-level variations.

During a single time step, marine sediment transport is performed until all available sediment transported by rivers to the ocean have been diffused and the accumulations on these specific nodes remains below water depth or below a prescribed slope (i.e., clinoform slope) computed based on local water depth and distance to the nearest coastline.

The transport coefficient depends on the thickness :math:`\mathrm{h_d}` of the fresh deposit, :math:`\mathrm{K_m(h_d) = K_{min} + K_m (1 - e^{-0.05 h_d})}` (with :math:`\mathrm{K_m}` = ``nonlinKm``, :math:`\mathrm{K_{min}} = 10^{-4}` m²/yr, and :math:`\mathrm{h_d}` set to zero below 0.1 m), and the equation is solved for the absolute surface (bed plus deposit) with an implicit finite-volume scheme. Across each cell face the coefficient is taken from the **higher** of the two cells, the one sediment moves out of. A cell carrying no fresh deposit therefore keeps only :math:`\mathrm{K_{min}}` and cannot shed its pre-existing bed into a neighbouring deposit. (An average of the two cells, used until November 2026, let bare seafloor above a deposit feed material into it; the deposit then carried about twice its delivered volume on steep bathymetry and a final rescale halved every deposit.) The ``picard`` integrator also limits supply: over a sub-step a cell cannot export more than the deposit it holds (its outgoing fluxes are scaled down when needed), so a thin deposit on a steep slope cannot cut into its bed on fine meshes. Every face flux is computed once, so the scheme conserves volume; the remaining small negative values are clipped, a final volume rescale corrects the difference (it should be close to 1), and the result does not depend on the number of processors.

Two time integrators are available (``marineSolver``). The default ``ts`` is an adaptive linearly-implicit Rosenbrock-W scheme (PETSc ``TS rosw``) whose stage systems are solved by GMRES with a GASM preconditioner; its error is controlled on the deposit thickness, not the absolute elevation. The ``picard`` option takes ``picardSub`` backward-Euler sub-steps with the coefficient lagged over ``picardIts`` iterations, each a linear solve. If the adaptive ``ts`` integrator fails (too many rejected steps), goSPL switches to ``picard`` for the rest of the run.

Before the diffusion, river sediment is routed across the sea floor along flow directions computed on a smoothed bathymetry whose depressions are filled by a priority flood. The fill starts at ``sea level + oFill``; on a global model it starts no lower than the deepest point of the smoothed surface plus 0.1 m, so the whole ocean floor drains to a few terminal sinks instead of hundreds of closed pockets.

Gravity currents (optional)
---------------------------

With a ``gravity`` block, part of the river load reaching the sea travels as a gravity current (hyperpycnal flow or turbidity current) instead of being spread by the clinoform scheme. A fraction :math:`\mathrm{f_h}` of the load arriving at each marine node is diverted and routed down the same seafloor flow network as the clinoform cascade (receiver matrix :math:`\mathrm{W}`). A fraction :math:`\mathrm{\Lambda_i}` of the flux passing a node deposits there and the rest moves on, so the flux :math:`\mathrm{q}` solves the linear system

.. math::

  \mathrm{\left(I - W^T \mathrm{diag}(1 - \Lambda)\right) q = s}, \qquad \mathrm{deposit_i = \Lambda_i q_i}, \qquad \mathrm{\Lambda_i = 1 - e^{-\ell_i / L}}

where :math:`\mathrm{s}` is the diverted supply, :math:`\mathrm{\ell_i}` the cell size and :math:`\mathrm{L}` the run-out length (``runout``; ``runout_fine`` for the fine fraction under dual lithology). On a uniform slope the deposit decays as :math:`\mathrm{e^{-x/L}}` whatever the mesh spacing. Terminal sinks keep what reaches them and open boundaries export it, so the supplied volume equals the deposited plus the exported volume to the solver tolerance. The system is solved once per fraction and step with the routing solver of the flow accumulation. The gravity deposit is added to the marine deposit after the diffusion, and under dual lithology the fine fraction of each node is the volume-weighted blend of the clinoform deposit and the coarse and fine gravity deposits, so the fine volume is conserved. With provenance tracers, each source class is routed by the same operator (its supply is the arriving composition times the gravity supply), so the class deposits sum to the gravity deposit and the deposit composition is blended in the same way.

The diverted fraction is either fixed (``fraction``) or computed from floods. goSPL's discharge :math:`\mathrm{\bar Q}` and load :math:`\mathrm{\bar Q_s}` are long-term means, and the mean concentration :math:`\mathrm{\bar C = \rho_s \bar Q_s/\bar Q}` is almost always far below the hyperpycnal threshold :math:`\mathrm{C_c}` (about 40 kg/m³, Mulder and Syvitski, 1995). With log-normal daily discharges :math:`\mathrm{X = Q/\bar Q}` (log-standard deviation :math:`\mathrm{\sigma}`, mean 1, :math:`\mathrm{\mu = -\sigma^2/2}`) and a rating curve :math:`\mathrm{Q_s \propto Q^b}`, the fraction of the load carried above the threshold is

.. math::

  \mathrm{f_h = \Phi\left(\frac{\mu + b\sigma^2 - \ln X_0}{\sigma}\right)}, \qquad \mathrm{X_0 = \left(\frac{C_c\, E[X^b]}{\bar C}\right)^{1/(b-1)}}, \qquad \mathrm{E[X^b] = e^{b\mu + b^2\sigma^2/2}}

where :math:`\mathrm{\Phi}` is the standard normal distribution function. The design, limitations and planned extensions (slope failure) are in ``docs/DESIGN_GRAVITY_CURRENTS.md``.
