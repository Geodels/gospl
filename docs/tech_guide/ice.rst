.. _ice:

==============================================
Ice sheets, glacial erosion and meltwater
==============================================

goSPL represents glaciers with a cheap, robust **diagnostic glacial-erosion
model**. Each step the ELA surface mass balance is routed downhill into an ice
discharge, from which an ice thickness, a basal sliding velocity, glacial
abrasion, till/moraine deposition, meltwater and an ice load are derived — one
linear solve, no time integration. The model is enabled by adding an ``ice``
section to the input file.

.. note::

   **Why a diagnostic, not full ice dynamics.** A true Shallow-Ice-Approximation
   (SIA) thickness solve was implemented and then removed. The ice margin
   :math:`H \ge 0` is a *free-boundary (obstacle) problem*, on which an implicit
   ``F(H)=0`` thickness solve diverged on real continental runs; and an
   ice-dynamics solve is stiff and **over-thickens km-scale continental ice at
   coarse resolution** over goSPL's long (:math:`10^2`–:math:`10^4` yr)
   timesteps. The diagnostic below is physical and robust at any resolution and
   over those long steps, which is why it is the only model goSPL ships.

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
below the ELA the ramp goes negative and the source becomes ablation (melt). The
accumulation part is scaled by ``accum_factor`` (a precipitation→ice conversion
fraction) and optionally capped at ``accum_max`` (m ice/yr); ablation is
amplified by ``melt``.

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

Ice routing and discharge
-------------------------

The **net mass balance** — accumulation above the ELA minus the ``melt``-scaled
ablation below it — is **routed downhill** on a *drainage-conditioned* bed using
a multiple-flow-direction (MFD) algorithm — the same flow-matrix / KSP machinery
that builds the river flow accumulation, with the number of directions set by
``icedir``. A single linear solve turns the per-cell net balance into an **ice
discharge** :math:`Q` (m\ :sup:`3`/yr) at each cell: the volume of ice passing
through it per year. Where downstream ablation has consumed all the upstream
accumulation the flux falls to zero — that is the glacier snout, emerging from
the mass balance rather than a hard cut-off. The ``melt`` factor sets how much
ablation is subtracted (``0`` = accumulation-only, ``1`` = true net balance,
``> 1`` = shorter/thinner tongues). There is no time integration of an
ice-thickness PDE.

Because the discharge is a single :math:`(\mathbf{I}-\mathbf{W}^{\!\top})`
solve, the routing surface must let **every cell above the glacier terminus
strictly drain** — no closed depressions and no flats — or the operator is
singular and the solve collapses to zero discharge (no ice). This conditioning
is built **in parallel**, anchored at the terminus (the physical ice outlet,
where ice melts out): a parallel priority-flood removes depressions and a
parallel flat-router gives every filled-flat cell a spill-ward direction, while
cells at or below the terminus and the domain edges act as absorbing outlets.
The result is partition-invariant (identical at any processor count) and
avoids gathering the whole mesh to one rank, so it scales to large meshes.

Ice thickness (Bahr scaling)
----------------------------

The **ice thickness** :math:`H` follows from a Bahr-type width–area (here
discharge) scaling of the discharge:

.. math::

   H = \mathrm{eheight}\,\cdot\,\mathrm{fwidth}\,\cdot\,Q^{\,0.3}

with the Bahr thickness factor ``eheight`` and width factor ``fwidth``. Thicker
ice forms where the routed discharge is larger (the trunk of the ice stream),
tapering to zero where there is no ice.

Basal sliding velocity (Glen sliding law)
-----------------------------------------

From that thickness and the bed-surface slope, goSPL derives the **basal sliding
velocity** :math:`u_b` from Glen's sliding law (the ``ice_velocity`` Fortran
kernel):

.. math::

   u_b \;\propto\; H^{\,n-1}\,|\nabla s|^{\,n-1}\,\nabla s,
   \qquad s = z_\mathrm{bed} + H

with sliding coefficient ``slide`` and Glen exponent :math:`n` = ``glen`` (usually
3). The velocity is physically bounded. It is written to the output as ``iceUb``
and drives glacial abrasion.

Glacial abrasion (vertical and lateral)
---------------------------------------

When ``abrasion.Kg > 0``, sliding ice abrades the bed (deepening valleys) at the
rate

.. math::

   E_g = K_g\,|u_b|^{\,l}

(``Kg`` and ``l`` in the ``abrasion`` block). The abrasion is masked to
subaerial, ice-covered cells (no marine abrasion; :math:`u_b` is already zero
where there is no ice).

An optional **lateral (valley-wall) erosion** term widens glaciated valleys
toward a U-profile. With ``abrasion.Kl > 0``, each wall cell flanking fast ice is
eroded at :math:`K_l\,|u_{b,\mathrm{neighbour}}|^{\,\mathrm{lat\_l}}`, tapered by
how much of the wall is in contact with the neighbouring ice column. The eroded
wall rock joins the same conserved till → moraine budget as the vertical
abrasion. ``lat_l`` (default = ``l``) is its velocity exponent. Both ``Kg`` and
``Kl`` default to ``0`` (abrasion off).

When **till handling is off**, the abraded material is added directly to the
erosion–deposition rate as an incision and flows into the fluvial sediment system
through the standard erosion path in all SPL flavours (``SPL``, ``nlSPL``,
``soilSPL``).

Glacial till and moraine deposition
-----------------------------------

When ``till.on`` is set (the default), the abraded rock is not released straight
into the rivers but carried as **till** by the ice and deposited as **moraine**
where the ice melts out — the ablation zone / terminus. Rock volume is conserved
(rock is moved, not created or destroyed). Two modes are used depending on
whether stratigraphy is active:

- **Bulk bed** (no stratigraphy): the abraded volume is redistributed as a
  bed-to-bed transport, so the net bed-volume change is exactly zero.
- **Stratigraphic**: the till is removed from the stratigraphic layers it was
  abraded from and re-deposited as a fresh moraine layer in the ablation zone.
  Under **dual lithology** the moraine is split into coarse and fine fractions
  carrying the abraded (ice-mixed) fine fraction, so the per-fraction solid
  budget stays balanced. The bed bulks up by the porosity contrast between the
  uncompacted till and the compacted source rock.

Independently of the bulk/stratigraphic split, the **deposition distribution**
has two options. With ``till.route: True`` (the default) the till is **routed
down the ice-surface flow network**: it is transported along steepest descent of
the surface :math:`s = z_\mathrm{bed} + H` and melts out progressively toward
each catchment's terminus, building moraine at the actual ice margins. This is a
transport-with-loss solved with the same MPI-correct flow-matrix / KSP machinery
as the river accumulation,

.. math::
   L_i = A_i + \sum_{u \to i} w_{ui}\,(1-f_u)\,L_u, \qquad D_i = f_i\,L_i,

where :math:`A_i` is the local abraded volume, :math:`L_i` the till load passing
through cell :math:`i`, and :math:`f_i = \min(1, \dot a_i \Delta t / H_i)` the
melt-out fraction (forced to 1 at the ice margin so no till leaks onto bare
ground). With ``till.route: False`` the **global** abraded volume is instead
spread across the whole ablation zone weighted by the local meltwater rate —
cheaper (no extra solve) but it decouples erosion and deposition across separate
ice masses. Both options conserve mass and are protected by dedicated tests.

Meltwater re-injection into the river network
---------------------------------------------

By default (``melt_conserve: True``) the glacial meltwater delivered to the
rivers is **discharge-conserving**. The precipitation that fell as ice above the
ELA is routed down-glacier and released as liquid meltwater where the ice melts
out, so the **total meltwater equals the total accumulation** — closing the
glacial water budget so downstream basins do not under-predict discharge. The
river accumulation step in ``flowAccumulation`` removes the precipitation that
was captured as ice from the river source and adds this meltwater back, so cells
downstream of glacier termini see the corresponding meltwater discharge instead
of losing it.

With ``melt_conserve: False`` goSPL reverts to the local precipitation-scaled
ablation rate,

.. math::

   \mathrm{m_i} = \max(-\dot{m}_i, 0)\;\Omega_i\;[H_i > 0]

which is cheaper but generally returns less water than fell as ice.

Ice loading and flexure
-----------------------

When ``flexure`` is enabled, the diagnostic ice thickness is used as the ice load
in the flexural-isostasy computation: the change in ice thickness between steps is
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
