.. _tech_gravity:

==============================
Marine gravity currents
==============================

With a ``gravity`` block, part of the river load reaching the sea travels as a gravity current (hyperpycnal flow or turbidity current) instead of being spread by the clinoform scheme. The currents are not resolved as flows (a depth-averaged turbidity-current model needs time steps of seconds); they are represented by a routing of sediment down the seafloor with deposition over a run-out length, which is linear, solved once per step, and exactly conservative. The user keys are listed in the user guide; this section gives the equations.

**Position in the time step.** The gravity supply is taken from the river load before the clinoform cascade; the cascade spreads the rest; the currents are routed on the seafloor as it was before the marine deposition of the step; the gravity deposit is added after the marine diffusion (it is not diffused in the step it forms), and elevation, cumulative erosion-deposition and stratigraphy then receive the combined marine deposit.

Supply
------

At each marine node :math:`\mathrm{i}` reached by rivers, the load :math:`\mathrm{Q_{s,i}}` (m³ over the step; :math:`\mathrm{Q_{s,i}^c}` coarse and :math:`\mathrm{Q_{s,i}^f}` fine under dual lithology) is split into a gravity supply

.. math::

  \mathrm{s_i^f = f_{h,i}\, Q_{s,i}^f}, \qquad \mathrm{s_i^c = \gamma\, f_{h,i}\, Q_{s,i}^c}

and a remainder that goes to the clinoform cascade. :math:`\mathrm{\gamma}` is ``coarse_factor`` (0.3 by default; gravity currents carry mostly the suspended load). The fraction :math:`\mathrm{f_h}` is either fixed (``fraction``) or carried by floods. goSPL's discharge :math:`\mathrm{\bar Q}` and load :math:`\mathrm{\bar Q_s}` are long-term means, and the mean concentration :math:`\mathrm{\bar C = \rho_s \bar Q_s/\bar Q}` is almost always far below the hyperpycnal threshold :math:`\mathrm{C_c}` (about 40 kg/m³, Mulder and Syvitski, 1995). With log-normal daily discharges :math:`\mathrm{X = Q/\bar Q}` (log-standard deviation :math:`\mathrm{\sigma}`, mean 1, so :math:`\mathrm{\mu = -\sigma^2/2}`) and a rating curve :math:`\mathrm{Q_s \propto Q^b}`, the fraction of the load carried by days above the threshold is

.. math::

  \mathrm{f_h = \Phi\left(\frac{\mu + b\sigma^2 - \ln X_0}{\sigma}\right)}, \qquad \mathrm{X_0 = \left(\frac{C_c\, E[X^b]}{\bar C}\right)^{1/(b-1)}}, \qquad \mathrm{E[X^b] = e^{b\mu + b^2\sigma^2/2}}

where :math:`\mathrm{\Phi}` is the standard normal distribution function and :math:`\mathrm{X_0}` the discharge ratio at which the daily concentration :math:`\mathrm{\bar C X^{b-1}/E[X^b]}` reaches :math:`\mathrm{C_c}`.

Routing network
---------------

The currents follow a receiver matrix :math:`\mathrm{W}` (:math:`\mathrm{W_{ij}}` is the share of node :math:`\mathrm{i}`'s outflow sent to node :math:`\mathrm{j}`, each row summing to 1, or 0 at a sink). With ``network: smoothed`` it is the matrix of the clinoform cascade. With ``network: raw`` (default) it is built every step on the seafloor :math:`\mathrm{z}` as follows.

1. Optional smoothing (``network_smooth`` = :math:`\mathrm{n_s}` passes), each pass an area-weighted average over the node and its neighbours :math:`\mathrm{N(i)}`:

   .. math::

     \mathrm{z_i \leftarrow \frac{a_i z_i + \sum_{j \in N(i)} a_j z_j}{a_i + \sum_{j \in N(i)} a_j}}

   with :math:`\mathrm{a_i}` the cell area.
2. Depression filling of the (smoothed) surface by a priority flood from the open boundaries (or the open ocean on a global mesh) and the deepest cells, with an infinitesimal increment so that every filled cell has a lower neighbour: :math:`\mathrm{\tilde z}`.
3. Receivers: among the neighbours strictly lower than :math:`\mathrm{i}`, the :math:`\mathrm{n_r}` (``receivers``) with the steepest slopes :math:`\mathrm{S_{ij} = (\tilde z_i - \tilde z_j)/d_{ij}}` (:math:`\mathrm{d_{ij}}` the distance between nodes; exact ties broken by the input-mesh node index), with weights

   .. math::

     \mathrm{W_{ij} = \frac{S_{ij}^{\,p}}{\sum_{k \in R(i)} S_{ik}^{\,p}}}

   where :math:`\mathrm{p}` is ``flow_exp`` and :math:`\mathrm{R(i)}` the selected receivers. A large :math:`\mathrm{p}` and few receivers concentrate the flow on the steepest path; on gentle floors the slopes are alike and the flow spreads. A node with no lower neighbour is a terminal sink. Open boundary nodes have no receivers and export.

Under dual lithology the fine fraction can have its own network (``receivers_fine``, ``flow_exp_fine``): with a concentrated coarse network and a spreading fine one, coarse sediment builds belts along the main pathways and fine sediment broad lobes beyond.

Deposition over a run-out length
--------------------------------

A fraction :math:`\mathrm{\Lambda_i}` of the flux :math:`\mathrm{q_i}` through node :math:`\mathrm{i}` deposits there and the rest moves on:

.. math::

  \mathrm{q_i = s_i + \sum_j W_{ji}\,(1 - \Lambda_j)\, q_j} \quad\Longleftrightarrow\quad \mathrm{\left(I - W^T \mathrm{diag}(1 - \Lambda)\right) q = s}, \qquad \mathrm{D_i = \Lambda_i q_i}

with

.. math::

  \mathrm{\Lambda_i = \left(1 - e^{-\ell_i / L}\right)\, T_d(d_i)\, T_S(S_i)}, \qquad \mathrm{T_d(d) = \min\left(1, \max\left(0, \frac{d - d_b}{d_b}\right)\right)}, \qquad \mathrm{T_S(S) = \frac{1}{1 + (S/S_b)^4}}

:math:`\mathrm{\ell_i = \sqrt{a_i}}` is the cell size, :math:`\mathrm{L}` the run-out length (``runout``, or ``runout_fine`` for the fine fraction), :math:`\mathrm{d_i}` the water depth, :math:`\mathrm{d_b}` ``bypass_depth``, :math:`\mathrm{S_i}` the seafloor gradient (Green-Gauss on the Voronoi cell, optionally smoothed ``slope_smooth`` times with the same average as above) and :math:`\mathrm{S_b}` ``bypass_slope`` (:math:`\mathrm{T_S = 1}` when it is not set). :math:`\mathrm{\Lambda = 1}` at terminal sinks (they keep what reaches them) and 0 at open boundaries (export). Along a path of cells the flux decays as :math:`\mathrm{e^{-x/L}}` whatever the mesh spacing, so :math:`\mathrm{L}` is resolution independent (``benchmarks/test_gravity.py`` checks it to round-off on square lattices at two spacings, with one and six receivers). The tapers are smooth on purpose: hard thresholds on a gradient that is noisy from cell to cell made neighbouring cells alternate between deposition and bypass.

**Accommodation.** A gravity deposit cannot build a cell above :math:`\mathrm{\eta - d_b}` (:math:`\mathrm{\eta}` the sea level). With :math:`\mathrm{r_i = \max\left(0, (\eta - d_b - z_i)\, a_i - D^{cl}_i\right)}` the room left after the clinoform deposit :math:`\mathrm{D^{cl}_i}`, every node with :math:`\mathrm{\Lambda_i q_i > r_i}` gets :math:`\mathrm{\Lambda_i = r_i / q_i}` (it deposits exactly its room and passes the rest on) and the system is solved again, until no node overfills (:math:`\mathrm{\Lambda}` only decreases, so a few passes suffice). The fine fraction gets the room the coarse one leaves.

Slope erosion
-------------

With ``erosion_length`` :math:`\mathrm{L_e}`, the currents pick up sediment where the seafloor is steep, in proportion to the flux:

.. math::

  \mathrm{E_i = \epsilon_i q_i}, \qquad \mathrm{\epsilon^0_i = \min\left(1,\; \frac{\ell_i}{L_e} \frac{(S_i/S_e)^4}{1 + (S_i/S_e)^4}\, T_d(d_i)\right)}, \qquad \mathrm{\epsilon_i = \min\left(1,\; \epsilon^0_i \left(\frac{q_i^{(0)}}{Q_e\, \Delta t}\right)^{a-1}\right)}

with :math:`\mathrm{S_e}` ``erosion_slope``, :math:`\mathrm{a}` ``erosion_exp`` (1 by default) and :math:`\mathrm{Q_e}` ``erosion_flux`` (m³/yr). The total flux then solves :math:`\mathrm{\left(I - W^T \mathrm{diag}(1 - \Lambda + \epsilon)\right) q = s}`. A node either erodes (:math:`\mathrm{\epsilon_i > \Lambda_i}`, and then :math:`\mathrm{\Lambda_i = 0}`) or deposits (:math:`\mathrm{\epsilon_i = 0}`). Where :math:`\mathrm{E_i}` would exceed ``max_erosion`` × :math:`\mathrm{a_i}`, :math:`\mathrm{\epsilon_i}` is lowered to fit and the solve repeated. For :math:`\mathrm{a > 1}`, :math:`\mathrm{q^{(0)}}` is the flux of the linear solve and the factor is evaluated once per step: iterating on :math:`\mathrm{q}` within a step diverges, because the pick-up multiplies the flux by :math:`\mathrm{\prod (1 + \epsilon)}` along a path. For the same reason :math:`\mathrm{\epsilon_i \le 1}` and :math:`\mathrm{L_e}` must span many cells. The eroded volume is removed from the stratigraphy (which gives its fine fraction and provenance) and added, with that composition, to the supply of each fraction at the eroding node, so the routing below treats it as an explicit source.

Slope erosion is an expert option: it needs submarine channels resolved by the mesh (several cells across channels a few hundred metres to a few kilometres wide), and goSPL warns when it is used on cells coarser than 200 m.

Slope failure
-------------

With ``failure_slope`` :math:`\mathrm{S_f}`, sediment on marine cells steeper than :math:`\mathrm{S_f}` fails and joins the gravity currents. The failed thickness at node :math:`\mathrm{i}` brings its gradient back to about :math:`\mathrm{S_f}`:

.. math::

  \mathrm{\Delta z_i = T_d(d_i)\, \min\left(\frac{(S_i - S_f)\, \ell_i}{2},\; z^{max}_f,\; D_i\right)}

with :math:`\mathrm{S_i}` the gradient (smoothed ``slope_smooth`` times), :math:`\mathrm{z^{max}_f}` ``max_failure`` (m per step) and :math:`\mathrm{D_i}` the sediment deposited during the run (the positive cumulative erosion-deposition, so bedrock never fails); lowering a cell by :math:`\mathrm{(S - S_f)\ell/2}` takes about :math:`\mathrm{S - S_f}` off its gradient. The depth taper keeps failures below the bypass depth. The failed volume follows the same path as the slope erosion: it is removed from the stratigraphy, which gives its composition, added to the supply at that node, and the node receives no gravity deposit in that step. It is reported apart (``failed``) and counted with the eroded volume in the budget below.

Lateral spreading of the deposit
--------------------------------

A concentrated network delivers the coarse load along paths one cell wide, which would build one-cell ridges. With ``deposit_spread`` :math:`\mathrm{L_s}` (m; two mean cell sizes by default, 0 to turn it off), each step's gravity deposit :math:`\mathrm{v_i}` (m³) is spread by a conservative diffusion of its thickness :math:`\mathrm{t_i = v_i / a_i}`:

.. math::

  \mathrm{v_i \leftarrow v_i + \sum_{j \in N(i)} w_{ij} (t_j - t_i)}, \qquad \mathrm{w_{ij} = w_{ji} = \frac{1}{2} \frac{\min(a_i, a_j)}{\max(n_i, n_j)}}

repeated :math:`\mathrm{\lceil 2 (L_s / \overline{\ell})^2 \rceil}` times (:math:`\mathrm{n_i}` the number of neighbours, :math:`\mathrm{\overline{\ell}}` the mean cell size), only between cells deeper than :math:`\mathrm{d_b}` that are not eroding or failing in that step. The exchange is antisymmetric, so volume is conserved exactly; a uniform thickness is left unchanged on any mesh; and :math:`\mathrm{\sum_j w_{ij} \le a_i/2}` keeps every pass positive. Each pass widens a deposit by about half a cell size squared in variance, so the same :math:`\mathrm{L_s}` gives the same width at any resolution. (An operator that shares each donor's volume in proportion to its neighbours' areas also conserves volume but does not keep a uniform layer uniform: on an irregular mesh it grew chains of metre-scale bumps at cells larger or smaller than their neighbours.)

Composition and budget
----------------------

The coarse and fine supplies (and, with provenance tracers, one column per source class, whose supply is the arriving composition times the gravity supply, plus the eroded classes) are routed with the same operator, so the class deposits sum to the gravity deposit. The fine fraction recorded in the stratigraphy at node :math:`\mathrm{i}` is the volume-weighted blend of the clinoform deposit :math:`\mathrm{D^{cl}_i}` (fine fraction :math:`\mathrm{F^{cl}_i}`) and the coarse and fine gravity deposits :math:`\mathrm{D^c_i}`, :math:`\mathrm{D^f_i}`:

.. math::

  \mathrm{F_i = \frac{F^{cl}_i D^{cl}_i + D^f_i}{D^{cl}_i + D^c_i + D^f_i}}

and the provenance composition is blended the same way, so the fine and class volumes are conserved. The stratigraphy receives the deposits only (the slope erosion has already been removed from it), while elevation and erosion-deposition rates receive the net change. Over each step

.. math::

  \mathrm{\sum_i s_i + \sum_i E_i = \sum_i (D^c_i + D^f_i) + X}

with :math:`\mathrm{X}` the volume exported through open boundaries; the closure is reported with verbose output and in the ``--summary`` record, and holds to the solver tolerance.

Numerics and cost
-----------------

Each system is triangular in flow order (receivers are strictly lower), so it is well posed; it is solved with the routing solver of the flow accumulation (``fgmres`` with block Jacobi, and an exact-factor rescue). One solve per fraction, plus one per accommodation or erosion-cap pass and one per provenance class. The matrices are assembled on owned rows only and every derived field (gradient, :math:`\mathrm{\Lambda}`, :math:`\mathrm{\epsilon}`, smoothed values, areas used by the spreading) is synchronised from its owner to the ghost nodes, so the result does not depend on the number of processors. The depression filling is done on one processor (the global surface is gathered on it, filled, and each processor receives back only its own nodes), as for the clinoform network; ``network: raw`` adds a second such fill each step, while ``smoothed`` reuses the clinoform one. On a 52,805-node regional mesh with 4 processors the module adds about 30 % to a step (routing 0.20 s, network 0.04 s, spreading 0.04 s). The design, the tests and the history of the choices are in ``docs/DESIGN_GRAVITY_CURRENTS.md``.

The user keys and the set-up by resolution are in :ref:`the user guide
<gravity>`.
