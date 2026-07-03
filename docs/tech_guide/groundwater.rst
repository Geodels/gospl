.. _groundwater:

==============================================
Water table and generic duricrust
==============================================

goSPL can evolve a **water table** and a generic **duricrust** — near-surface
hydrology coupled to a chemical armoring of the erodibility. Each step a shallow
groundwater head is solved on the mesh from the infiltrated rainfall; where the
water-table depth lies in a capillary fringe an indurated crust precipitates and
resists erosion, producing relief inversion and — with stratigraphy on — stacked,
buried and re-exhumed duricrusts. The feature is enabled by a ``groundwater``
section in the input file and is a strict no-op when absent. The full design
(algorithm, invariants, caveats) is in ``docs/DESIGN_WATERTABLE_DURICRUST.md``.

The whole update runs once per step, **after flow accumulation and before
erosion**, so each process reads the previous process's state (soil → groundwater
→ duricrust → erosion). Every step of that loop is slow (:math:`10^4`–:math:`10^6`
yr) and the coupling is explicit/sequential within a step, so there is no stiff
intra-step feedback — the same stability argument as goSPL's other explicit
couplings.

Water-table head — Dupuit–Boussinesq
------------------------------------

The saturated groundwater flow under a free surface is described by the
unconfined (Dupuit–Boussinesq) approximation. The water-table head :math:`h`
(elevation of the saturated surface) obeys

.. math::

   S\,\frac{\partial h}{\partial t} = \nabla\!\cdot\!\big(T(h)\,\nabla h\big) + R,
   \qquad T(h) = K_h\,\max(h - z_{bed},\, b_{min}),

with :math:`S` the specific yield (drainable porosity), :math:`T` the
transmissivity, :math:`K_h` the hydraulic conductivity, :math:`z_{bed}` the
impermeable base, :math:`b_{min}` a minimum saturated thickness (so
:math:`T>0`), and :math:`R` the net recharge. Recharge is the fraction of
effective rainfall that infiltrates,

.. math::

   R = f_{infil}\,\max(0,\, \mathrm{rain} - \mathrm{evap}),

held at zero under standing water (sea, ponded lakes) and under ice, where the
head is pinned to the surface or precipitation does not infiltrate.

The infiltration fraction :math:`f_{infil}` is a scalar or a per-vertex map, and
three optional (default-off) modulations refine it: **subglacial recharge** lets
a fraction of the glacial meltwater (the ice model's ``iceMeltRiverL``) infiltrate
under ice — the one recharge path allowed there; **lithology** scales it by the
exposed coarse fraction so coarse infiltrates more than fine (dual lithology); and
**slope** reduces it on steeper terrain as :math:`f/(1+\mathrm{slope}/\mathrm{slope}_{ref})`.

Implicit discretisation (and why)
---------------------------------

The groundwater response time :math:`\tau_{gw}` is short relative to goSPL's
century-to-millennium steps :math:`\Delta t`, so an explicit update would be
severely time-step limited. goSPL instead takes **one implicit (backward-Euler)
solve per step** on the mesh's finite-volume Laplacian :math:`L(T)` (the same
operator the marine diffusion uses):

.. math::

   \Big(I + \tfrac{\Delta t}{S}\,L(T)\Big)\,h = h_{old} + \tfrac{\Delta t}{S}\,R.

As :math:`\Delta t/S` grows this tends to the steady elliptic limit
:math:`L(T)\,h = R`, so the solve is quasi-steady at long steps — exactly the
regime goSPL runs in. The operator is a stiff 2-D elliptic system: it is solved
with a Krylov method preconditioned by **algebraic multigrid** (fgmres + hypre
BoomerAMG). A simple block-Jacobi/ILU preconditioner does *not* converge on it
and silently drifts the water table to the surface, so multigrid is required —
the same lesson as the flexure biharmonic.

.. note::

   The implicit solve was validated against the exact steady 1-D Dupuit parabola
   :math:`s(x)^2 = s_d^2 + (R/K)(2Wx - x^2)` on a drained strip (the
   ``benchmarks/test_dupuit.py`` fixture): the steady head matches the analytic
   profile to a root-mean-square error of 0.00 %.

Two non-linearities are handled with cheap inner iterations:

- **Unconfined transmissivity** :math:`T(h)`: a few Picard sweeps
  (``picard_its``) lag :math:`T` on the previous iterate.
- **Seepage free boundary** :math:`h \le z`: the head is pinned to the surface
  (Dirichlet) at the base seepage set — sea, ponded lakes and open outlets —
  then, after each Picard block, any cell whose head over-topped the surface is
  added to that set and the solve repeated, up to ``seepage_passes`` passes. A
  dry-aquifer floor :math:`h \ge z_{bed}` is imposed by clipping.

The aquifer base :math:`z_{bed}` is prescribed (``aquifer_base`` scalar or
per-vertex map), or — with ``aquifer_base: from_soil`` — tied to the bedrock
elevation :math:`z_{bed} = l_{Hbed} - d_{bedrock}` (the *permeable regolith over
impermeable bedrock* model of laterite terrains), which requires soil tracking.
In a **depositional basin** the porous sediment fill is itself the aquifer, so the
base is taken as the deeper of :math:`l_{Hbed} - d_{bedrock}` and the bottom of the
stratigraphic pile :math:`z - \sum H` — uplands keep the thin-regolith base, basins
deepen to the fill bottom.

Baseflow closure
----------------

With ``conserve_baseflow`` on, the recharge that does not go into aquifer storage
returns to the surface-water network as **baseflow** at the seepage nodes,

.. math::

   Q_{seep} = \sum_{owned}\!\big(R\,A\big) - \sum_{owned} S\,\frac{h - h_{old}}{\Delta t}\,A ,

(written as the ``baseflow`` output). It is also **re-injected into the surface-flow
source**: the runoff becomes :math:`b = \mathrm{rain}\,A - R\,A + Q_{seep}`, so the
infiltrated recharge leaves surface runoff and re-emerges at the seepage nodes —
rivers are physically baseflow-fed and total discharge stays
:math:`\approx \mathrm{rain} - \mathrm{evap}` (net-neutral globally; in the steady
limit :math:`\sum \mathrm{baseflow} \approx \sum \mathrm{recharge}`). It is applied
once per step (the single runoff-source reset) so the two per-step flow-accumulation
solves do not double-count it.

Lake ↔ aquifer volume coupling
------------------------------

Lakes, rivers and the sea are fixed-head boundaries (:math:`h = z`) for the head
solve. With the opt-in ``lake_exchange`` the **lake volume budget** is additionally
debited/credited by the across-bed groundwater flux. After the head solve the signed
per-node exchange :math:`\nabla\!\cdot\!(T\nabla h)\,A = -(L h)\,A` is formed
(positive = the aquifer discharges *into* the lake, negative = the lake *leaks* into
the aquifer), aggregated per lake, and fed into the lake's inflow — gains added,
leakage debited and clamped at the available water (mirroring the lake-surface
evaporation). Off by default (the conventional fixed-head treatment).

Duricrust formation and armoring
--------------------------------

A generic capillary-fringe crust grows from the water-table depth
:math:`wt = z - h`. Its thickness :math:`duriH` evolves as

.. math::

   \frac{\mathrm{d}\,duriH}{\mathrm{d}t}
     = k_{form}\,\Phi\,\Psi\,\Big(1 - \frac{duriH}{duriH_{max}}\Big),
   \qquad
   \Phi = \exp\!\left[-\Big(\frac{wt - d_0}{w}\Big)^2\right],

self-limiting to :math:`duriH_{max}`. The **fringe favourability** :math:`\Phi`
is a Gaussian band on the water-table depth centred on the mean fringe depth
:math:`d_0` (half-width :math:`w`): the crust forms only where the table sits at
the fringe. The **weathering supply** :math:`\Psi` has three modes: a default
climate/temperature *proxy* :math:`\max(0, \mathrm{rain}-\mathrm{evap})^{p}`, an
explicit *rate* (a Maher–Chamberlain kinetic–thermodynamic law driven by the
recharge), or *prodsoil* (reuse of the soil production rate). When soil is
tracked, formation is additionally **regolith-limited** — chemical crust growth
cannot outpace physical regolith production. A breakdown term strips the crust
under surface incision and relaxes it away from the fringe.

The induration degree :math:`duriF = duriH/duriH_{max} \in [0,1]` is the single
control on erodibility. It enters as a multiplicative **armoring factor** at the
shared surface-erodibility hook,

.. math::

   K_{eff} = K \cdot \mathrm{surfLithoK} \cdot \big(1 - armor_{max}\,duriF\big),

so a fully indurated cell (:math:`duriF = 1`) is :math:`armor_{max}` less
erodible (e.g. 0.9 ⇒ 10× more resistant). This modulates **all three eroders**
(SPL, non-linear SPL and soil-aware SPL) with no branching in them. Armoring is
**rate-only**: it changes no elevation the step it forms, so flow routing and
mass balance are untouched — relief inversion emerges through the normal erosion
pathway because crusted cells simply erode more slowly.

Stratigraphic induration record
--------------------------------

When stratigraphy is recorded (a positive ``strat`` interval), the induration is
archived per layer as ``stratDuri`` — an intensive property in :math:`[0,1]`,
distinct from the depositional erodibility ``stratK`` (they multiply at the
exposed surface). Each step the live crust is written down across its **depth
range** — every layer whose top lies within the crust thickness :math:`duriH`
below the surface is raised toward the live induration, so a thick crust indurates
the several thin layers it physically spans (**formation**); a fresh deposit buries
it with an uncemented (0) layer (**preservation**); and when erosion exhumes a
previously buried indurated layer its preserved value re-arms the surface, which
therefore resists incision over the crust's full thickness (**exhumation**). This
is the multi-cycle, stacked-duricrust / relief-inversion behaviour of cratonic
(e.g. Australian laterite) landscapes. ``stratDuri`` advects with the pile like the porosity, is
compaction-neutral, and is written to / restored from the stratal output; it is
exposed per layer as the ``induration`` field of ``gospl-strata-volume``.

Conservative geochemistry (Level B, opt-in)
-------------------------------------------

By default the duricrust supply is a *local, non-conservative* rate (§ above).
Adding a ``groundwater: geochem:`` block turns on a **conservative solute cycle**:
one or more lumped tracers are **dissolved**, **transported** with the groundwater
flux, **precipitated** at the fringe (feeding the crust) and **exported** to the
rivers, with a closed mass budget. See ``docs/DESIGN_WATERTABLE_GEOCHEM.md``.

Each tracer ``c`` (an ``n_species`` array — a single tracer by default, several
for calcrete / silcrete / ferricrete typing) obeys a **steady advection–reaction**
balance, solved once per step (the transport equilibrates within a century step,
just like the head):

.. math::

   \nabla\!\cdot\!(q\,c) = D - P, \qquad q = -T\nabla h,

- **Dissolution** :math:`D` — a chemical-weathering source on subaerial land
  (the climate/temperature supply scaled per tracer by ``weatherability``),
  drawing from a conserved per-node source pool.
- **Transport** — first-order upwind advection by the groundwater flux ``q``,
  built from the head operator's face conductances (geometry-correct on flat and
  global meshes); a diagonal **seepage sink** where the aquifer discharges makes
  the operator a well-posed M-matrix.
- **Precipitation** :math:`P = k_p\,\Phi\,c` — a linear sink at the capillary
  fringe ``Φ`` that **feeds the crust** ``duriH`` (replacing the local proxy
  supply when geochem is on). The saturation threshold ``c_sat`` is a future
  refinement.
- **Export** — the seepage sink removes the solute discharging to the surface;
  because the upwind operator is exactly conservative, **each step
  ``dissolved = precipitated + exported``** (machine precision). The dissolved
  export is the ``soluteflux`` output and a per-tracer flux to the ocean — the
  weathering-derived alkalinity/solute delivery relevant to long-term carbonate
  and climate budgets.

With several tracers, the dominant crust former per node is written as the
``crust_type`` field (different tracers win in different settings). With in-model
:ref:`provenance <surfproc>` on, the solute is additionally transported **per
source-rock class** (by linearity of the transport operator), so the downstream
crust is **attributed to the upgradient region where its solute dissolved** — the
dominant source is the ``crust_source`` field (solute-source provenance). The
crust's dominant species and dominant source region are additionally archived
**per stratigraphic layer** (``stratCrustType`` / ``stratCrustSource``, the
categorical companions to the ``stratDuri`` induration degree) — frozen on
burial and exposed by ``gospl-strata-volume`` — so a cross-section shows *what*
each buried crust is and *where* its chemistry came from. Outputs:
``solute`` (concentration), ``soluteflux`` (export), ``crust_type``,
``crust_source``. Coupled
multi-species aqueous equilibrium (speciation, pH) is out of scope by design —
the lumped multi-tracer model is the scale-appropriate choice at km / My.

Compatibility
-------------

The water table and duricrust compose with the other stratigraphic features:
``stratDuri`` is one more independent, intensive per-layer field alongside the
dual-lithology fine pile (``stratHf``/``phiF``) and the provenance partition
(``stratP``), and the armoring is a multiplicative factor at the same erodibility
hook the others already use. A run with groundwater, duricrust, dual lithology
and provenance all enabled conserves every invariant.
