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

.. note::

   The diffusivity scales as :math:`H^{\,n+2}|\nabla s|^{\,n-1}`, so an explicit
   integration would require a CFL time step orders of magnitude smaller than
   goSPL's km / :math:`10^2`–:math:`10^4` yr resolution. goSPL therefore solves
   the thickness **implicitly**.

Implicit thickness solve
-------------------------

Each goSPL time step the thickness is advanced with a single implicit
(semi-implicit) solve of the backward-Euler residual

.. math::

   F(H) = H - H_\mathrm{old} - \Delta t\,\big(\dot{m} - \nabla \cdot \mathbf{q}(H)\big) = 0,

using a cached, Jacobian-free PETSc ``SNES`` (``ngmres`` with a CG / HYPRE
BoomerAMG inner solve) — the same non-linear-diffusion machinery as the
non-linear hillslope solver. The scheme is unconditionally stable, so it takes
the **full goSPL time step in one solve**. The free boundary :math:`H \ge 0` is
enforced by clamping after convergence, and ice is removed below the glacier
terminus elevation (``hterm``). The solver is validated against the analytical
SIA dome: with zero surface mass balance the flux conserves ice volume to the
numerical floor.

From the converged thickness, goSPL derives the **basal sliding speed**
:math:`u_b` (``ice_velocity`` kernel), which is written to the output as
``iceUb`` and drives glacial abrasion.

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

.. note::

   Till is deposited across the ablation zone weighted by the meltwater rate
   rather than routed per ice-catchment to each individual terminus. Both modes
   are protected by a dedicated mass-conservation test.

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
