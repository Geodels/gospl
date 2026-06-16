.. _ice:

==============================================
Ice sheets, glacial erosion and meltwater
==============================================

goSPL represents glaciers with a **Shallow-Ice-Approximation (SIA) ice-sheet
model**: the ice thickness is evolved as a non-linear diffusion of the ice
surface, and the resulting flow drives glacial abrasion, till transport and
flexural loading. The model is enabled by adding an ``ice`` section to the
input file.

Surface mass balance
--------------------

The accumulation/ablation source :math:`\dot{m}` (m ice/yr) at each cell is a
fraction of the local precipitation set by the surface elevation relative to two
time-varying altitudes — the equilibrium-line altitude (``hela``) and the
ice-cap altitude (``hice``):

.. math::

   \mathrm{\dot{m}_i} = \mathrm{P_i} \cdot \min\!\left(1, \frac{\eta_i - h_{ela}}{h_{ice} - h_{ela}}\right)

Above the ice-cap altitude the cell accumulates the full local precipitation as
ice; between the ELA and the ice cap the fraction ramps linearly from 0 to 1;
below the ELA the ramp goes negative and the source becomes ablation (melt).

A degenerate-configuration guard zeroes the ice and returns immediately when
``hice <= hela`` or when the maximum surface elevation lies below the ELA.

.. note::

   For **global** models a single ELA is unphysical — it ranges from ~5000–6000 m
   in the tropics to near sea level at the poles. ``hela``, ``hice`` and
   ``hterm`` may therefore each be a **per-vertex map** rather than a scalar, and
   may vary **in time** through a ``glaciers`` time series (mirroring the
   precipitation ``climate`` block). The mass-balance ramp is then evaluated
   per node with its local ELA / ice-cap altitude, so the same run can grow
   tropical summit glaciers and polar ice sheets simultaneously. See the
   :ref:`user guide <surfproc>` for the input syntax.

Ice dynamics (SIA)
------------------

The ice thickness :math:`H` evolves as a non-linear diffusion of the ice surface
:math:`s = z_\mathrm{bed} + H`:

.. math::

   \frac{\partial H}{\partial t} = \dot{m} - \nabla \cdot \mathbf{q},
   \qquad
   \mathbf{q} = -\left[\frac{2A(\rho_i g)^n}{n+2}\,H^{\,n+2} + a_s\,H^{\,n}\right]
                |\nabla s|^{\,n-1}\,\nabla s

The flux :math:`\mathbf{q}` combines Glen's-law internal **deformation** (rate
factor :math:`A` = ``Aglen``, exponent :math:`n` = ``glen``, ice density
:math:`\rho_i = 910\ \mathrm{kg\,m^{-3}}`) and basal **sliding** (coefficient
:math:`a_s` = ``slide``). Its divergence is evaluated on the unstructured
finite-volume mesh by the ``ice_flux`` Fortran kernel, using the same
mass-conservative edge fluxes as the hillslope diffusion operator. Because the
flux is computed on the *total* surface (bed + ice) it naturally fills and
overflows closed basins.

Explicit thickness solve
------------------------

Each goSPL time step the thickness is advanced by an **explicit, mass-conserving,
positivity-preserving** integration of

.. math::

   \frac{\partial H}{\partial t} = \dot{m} - \nabla \cdot \mathbf{q}(H),

split into one or more substeps :math:`H \leftarrow H + \Delta t_\mathrm{sub}
(\dot{m} - \nabla \cdot \mathbf{q})`. The flux divergence comes from the
``ice_flux_limiter`` / ``ice_flux_rscaled`` kernels, which cap each cell's
outflux to the ice it actually holds,

.. math::

   R(k) = \min\!\Big(1,\; \frac{H_k A_k}{\Delta t_\mathrm{sub}\,\text{outflux}_k}\Big),

and scales every edge flux by the source cell's :math:`R`. Because the scaled
flux is applied equal-and-opposite across each shared face, ice volume is
conserved to machine precision; because a cell never exports more ice than it
has, the free boundary :math:`H \ge 0` is preserved for **any** substep with no
mass-injecting clamp.

.. note::

   The ice margin is a *free-boundary (obstacle) problem*, which makes an
   implicit ``F(H)=0`` thickness solve ill-posed at the margin (it was tried,
   and diverged on real runs). The flux limiter handles the boundary natively.
   Positivity being unconditional, the substep size is an **accuracy** choice
   set by ``sia.cfl`` (a fraction of the per-cell time to lose its ice), capped
   by ``sia.max_substeps``. The stable substep at goSPL's km /
   :math:`10^2`–:math:`10^4` yr resolution is :math:`10^3`–:math:`10^4` yr, so a
   step typically needs only a handful of substeps (often one).

Ice is removed below the terminus floor :math:`\max(h_\mathrm{term}, \text{sea
level})` — so it never persists below the (possibly time-varying) sea surface,
and an ``hterm`` below sea level is raised to it. When ``hterm`` is omitted the
floor is simply the sea-level position, leaving the dynamics and ablation to set
the terminus. The solver is validated against the analytical SIA dome: with zero
surface mass balance the flux conserves ice volume to the numerical floor, even
for a thick dome at a :math:`5\times10^4` yr step.

From the converged thickness, goSPL derives the **basal sliding speed**
:math:`u_b` (``ice_velocity`` kernel), which is written to the output as
``iceUb`` and drives glacial abrasion.

By default the ice grows in from zero. A **pre-existing ice thickness** can be
seeded at the first step with ``hinit`` (a uniform scalar or a per-vertex map);
the SIA solve then evolves it. The flexural reference is taken after seeding, so
a pre-existing ice load does not shock the plate, and a restart restores the
evolved thickness instead of re-seeding.

Meltwater re-injection into the river network
---------------------------------------------

Where the surface mass balance is negative and ice is present, the ablation rate
is captured as liquid meltwater (m\ :sup:`3`/yr):

.. math::

   \mathrm{m_i} = \max(-\dot{m}_i, 0)\;\Omega_i\;[H_i > 0]

The river accumulation step in ``flowAccumulation`` removes the precipitation
that was captured as ice (above the ELA) from the river source and then adds
this meltwater back, so cells downstream of glacier termini see the
corresponding meltwater discharge instead of losing it.

Glacial abrasion
----------------

When ``abrasion.Kg > 0``, sliding ice abrades the bed at the rate

.. math::

   E_g = K_g\,|u_b|^{\,l}

(``Kg`` and ``l`` in the ``ice`` block). The abrasion is masked to subaerial,
ice-covered cells (no marine abrasion; :math:`u_b` is already zero where there
is no ice). When **till handling is off**, :math:`E_g` is added directly to the
erosion–deposition rate as an incision and flows into the fluvial sediment
system through the standard erosion path in all SPL flavours (``SPL``,
``nlSPL``, ``soilSPL``).

Glacial till and moraine deposition
-----------------------------------

When ``till.on`` is set, the abraded material is not released straight into the
rivers but carried as **till** by the ice and deposited where the ice melts out
— the ablation zone / terminal moraine. The bed is lowered under fast-sliding
ice and raised in the ablation zone, weighted by the local meltwater rate. Two
modes are used depending on whether stratigraphy is active:

- **Bulk bed** (no stratigraphy): the abraded volume is redistributed as a
  bed-to-bed transport, so the net bed-volume change is exactly zero (rock is
  moved, not created or destroyed).
- **Stratigraphic**: the till is removed from the stratigraphic layers it was
  abraded from and re-deposited as a fresh moraine layer in the ablation zone.
  Under **dual lithology** the moraine is split into coarse and fine fractions
  carrying the abraded (ice-mixed) fine fraction, so the per-fraction solid
  budget stays balanced. The bed bulks up by the porosity contrast between the
  uncompacted till and the compacted source rock.

Independently of the bulk/stratigraphic split, the **deposition distribution**
has two options. By default the till is spread across the ablation zone weighted
by the local meltwater rate — adequate when a mesh cell aggregates a whole
glacier (continental/global resolution). With ``till.route: True`` the till is
instead **routed down the ice-surface flow network**: it is transported along
steepest descent of the surface :math:`s = z_\mathrm{bed} + H` and melts out
progressively toward each catchment's terminus, building moraine at the actual
ice margins. This is a transport-with-loss solved with the same MPI-correct
flow-matrix / KSP machinery as the river accumulation,

.. math::
   L_i = A_i + \sum_{u \to i} w_{ui}\,(1-f_u)\,L_u, \qquad D_i = f_i\,L_i,

where :math:`A_i` is the local abraded volume, :math:`L_i` the till load passing
through cell :math:`i`, and :math:`f_i = \min(1, \dot a_i \Delta t / H_i)` the
melt-out fraction (forced to 1 at the ice margin so no till leaks onto bare
ground). Routing is appropriate for high-resolution (sub-km) regional runs where
individual glacier catchments and termini are resolved; at coarse resolution it
reduces to near-local deposition. Both options conserve mass and are protected by
dedicated tests.

Ice loading and flexure
-----------------------

When ``flexure`` is enabled, the SIA ice thickness is used as the ice load in
the flexural-isostasy computation: the change in ice thickness between steps is
converted to an equivalent load (scaled by :math:`\rho_i / \rho_c`) and applied
through the existing flexure path, so a growing ice sheet drives isostatic
subsidence and deglaciation drives rebound.

Output
------

Four ice fields are written to the HDF5 / XDMF output whenever the ice model is
on:

- ``iceH`` — ice thickness (m); restored on restart;
- ``iceUb`` — basal sliding speed (m/yr);
- ``iceMelt`` — ablation meltwater (m\ :sup:`3`/yr) re-injected into the rivers
  (the glacial contribution to downstream discharge);
- ``iceAbr`` — glacial abrasion rate :math:`E_g = K_g|u_b|^{l}` (m/yr); zero
  where abrasion is off (``Kg = 0``).
