.. _carb:

==============================
Carbonate platforms and reefs
==============================

This page gives the equations behind the carbonate options: the lithology
list they rely on, benthic production (``carbonate:``), the open-water index,
and the pelagic rain (``pelagic:``). The user keys and the set-up by
resolution are in :ref:`the user guide <carbonates>`. Design notes and the
history of the choices: ``docs/DESIGN_CARBONATES.md``.

.. contents:: On this page
   :local:
   :depth: 1


Lithologies
-----------

With dual lithology a layer holds coarse and fine sediment. The lithology
list generalises it: a layer's bulk thickness :math:`H` is shared between
lithologies :math:`l = 1 \dots n` with piles :math:`H_l` (the dual-lithology
fine fraction is the first one), and coarse clastic sediment takes the rest,
:math:`H_c = H - \sum_l H_l`. Each lithology has

- a porosity-depth curve :math:`\phi_l(z) = \phi_{0,l}\, e^{-z/z_{0,l}}`, so
  each pile compacts on its own: the solid volume
  :math:`H_l (1 - \phi_l)` is kept and the pile recomputed with the porosity
  at its new depth (porosity never increases on unloading);
- factors on the river erodibility and the hillslope diffusivity, applied in
  proportion to its share of the surface layer,
  :math:`K = K_0 \left(f_c + \sum_l f_l k_l\right)`;
- a **size**, fine-like or coarse-like, which decides how it is sorted.

**Sorting by group.** The sorting rules of dual lithology (a filled pit keeps
coarse sediment first and lets fine overspill; the marine fine fraction
increases with water depth; coarse and fine gravity currents have their own
run-out) act on two groups: the fine group (the fine-like lithologies) and
the coarse group (coarse clastic and the coarse-like lithologies). Each
deposit has a fine-group fraction :math:`F`, and lithology :math:`l` takes the
share :math:`s_l` of its group:

.. math::

    \Delta H_l = \Delta H \times \begin{cases} F\, s_l & l \ \text{fine-like}
    \\ (1 - F)\, s_l & l \ \text{coarse-like} \end{cases},
    \qquad \sum_{l\ \text{fine}} s_l = 1, \quad \sum_{l\ \text{coarse}} s_l \le 1 .

The shares are the composition of what arrived, inside its group:

- *rivers*: each lithology's eroded volume is routed with the same implicit
  operator as the total, :math:`s_l = q_l / \sum_{\text{group}} q`;
- *pits*: a pit retains the same fraction of each lithology as of its group,
  so its deposit has the mix of its retained volumes;
- *marine*: each lithology is carried through the marine cascade and
  diffusion with the total (the same machinery as the provenance classes), so
  carbonate from a platform stays on its margin instead of being averaged with
  a delta's mud across the ocean. The depth-biased fine fraction then fixes
  each cell's fine and coarse volumes, and the routed mix is balanced by
  iterative proportional fitting so that each cell keeps its group volume and
  each lithology deposits its share of the group input (what leaves through
  open edges carries the same mix);
- *gravity currents*: each lithology rides its group's route;
- *hillslope creep and channel supply*: the donors' eroded lithologies,
  weighted by their creep flux into the cell and their group volume.

Every lithology is conserved. A fraction split into two identical lithologies
reproduces the dual-lithology run to round-off, which is how this is tested.


Benthic production
------------------

Below sea level, factory :math:`f` produces (Bosscher and Schlager, 1992)

.. math::

    P_f(w) = g_f \tanh\left(\frac{I_0\, e^{-k_f w}}{I_{k,f}}\right)
             M_f \prod_j R_{fj}

where :math:`w` is the water depth, :math:`g_f` the maximum growth rate,
:math:`I_0` the surface light, :math:`k_f` the extinction coefficient,
:math:`I_{k,f}` the saturating light, :math:`M_f` an optional habitat map and
:math:`R_{fj}` piecewise-linear response factors on the depth, the seafloor
slope, the clastic deposition rate of the previous step, the open-water
index and the sea-surface temperature.

**Integration over a step.** Each cell grows with the depth it creates,

.. math::

    \frac{dh}{dt} = \sum_f P_f(\eta - h), \qquad h \le \eta - w_{min},

with :math:`\eta` the sea level. The step is split into second-order (Heun)
sub-steps of about 0.25 m of growth at the largest possible rate (the same
number on every processor, from the inputs only), so the result does not
depend on the model time step: at 1 and 10 kyr steps it matches the exact
solution of this equation to :math:`10^{-4}`.

**Keep-up and give-up.** Under a relative sea-level rise :math:`\sigma` the
depth evolves as :math:`dw/dt = \sigma - P(w)`. Production grows as the water
shallows, so the depth where :math:`P(w^*) = \sigma`,

.. math::

    w^* = \frac{1}{k} \ln\left(\frac{I_0}{I_k\, \mathrm{atanh}(\sigma/g)}\right),

is unstable: a platform shallower than :math:`w^*` reaches sea level and keeps
up with it (its growth then limited to :math:`\sigma` by the accommodation); a
deeper one drowns, and its depth eventually increases at the rate
:math:`\sigma`.

**Framework and loose sediment.** A fraction :math:`\phi_f` of each factory's
production is built in place (lithology ``reef``); the rest is split into
carbonate sand (coarse-like, fraction :math:`\sigma_f`) and mud (fine-like)
and added to the marine cascade of the same step, after the gravity supply is
taken from the river load. Framework thinner than the stratigraphic floor
(0.1 mm) joins the loose sand, so elevation, stratigraphy and budget stay
consistent.


Sea-surface temperature
-----------------------

The ``sst`` block gives :math:`T(\mathbf{x}, t)`: event maps (or uniform
values) as a step function in time, plus a curve interpolated linearly in time
(an anomaly when both are given). The ``temperature`` response is a
piecewise-linear factor :math:`R(T)` on a factory's production, and on the
pelagic carbonate rain. It is evaluated once per step, node by node.


Open-water index
----------------

Reefs grow best at the platform margin, facing open water, and less in
lagoons and behind barriers (the distance-to-open-water penalty of Pastier et
al., 2019, the wave-energy control of pyBadlands). goSPL measures this with
the elevation minus its smoothing over a length :math:`L`
(``exposure_length``):

.. math::

    E = h - s, \qquad (I - L^2 \nabla^2)\, s = h .

:math:`E` is positive on margins, rims and highs and negative in lagoons. It
is solved for :math:`E` directly, :math:`(I - L^2 \nabla^2) E = -L^2 \nabla^2 h`,
so the solver tolerance applies to the relief rather than to the elevation.
For a parabolic ridge :math:`h = -b x^2` the index is exactly :math:`2 b L^2`.
Because :math:`L` is a length, the index does not change with the mesh
resolution (a smoothing set in numbers of cells would).


Wave energy and shadows
-----------------------

For each direction of the wave climate a relative energy :math:`E` travels
along the local wave direction :math:`\hat{\mathbf{d}}`,

.. math::

    \left(\hat{\mathbf{d}} \cdot \nabla + \frac{1}{L_f} + \kappa\right) E
    = \frac{1}{L_f}, \qquad
    \kappa = \frac{1}{L_d} + \frac{1}{L_s}\max\left(0, 1 - \frac{w}{w_s}\right),

with :math:`E = 0` on land (the shadow), :math:`L_f` the fetch over which the
waves grow to full energy (no growth term without it), :math:`L_d` the
attenuation length in deep water, :math:`w_s` the depth below which the waves
feel the bottom and :math:`L_s` the attenuation length at the shoreline (the
shoreward dissipation of Pastier et al., 2019). On a flat mesh the domain edges
facing the waves impose :math:`E = 1`; a global mesh has none and relies on the
fetch. Far from land the energy tends to
:math:`E_\infty = (1/L_f) / (1/L_f + 1/L_d)`; behind an island it grows again
from zero over the fetch.

The direction is a bearing turned into a tangent vector with the local east
and north (on the sphere) or :math:`x, y` (flat): one bearing, a map, or the
zonal belts. The scheme is the first-order upwind finite-volume operator of the
orographic rain and of the tectonic advection, with a velocity per node; the
coefficients of each direction are computed once, the matrix every step (the
coastline moves). The climate is :math:`E = \sum_k w_k E_k`. Along a mesh axis
the scheme is the recurrence
:math:`E_i = (E_{i-1} + \Delta x / L_f) / (1 + (1/L_f + \kappa)\Delta x)`;
the benchmark recovers it to :math:`10^{-14}` from the inflow edge and behind
an island (away from the half cells of the walls and outflow edge), and
:math:`(1 + \kappa\Delta x)^{-i}` without fetch.


*Refraction* (``refraction``). In water shallower than :math:`w_s` the local
direction turns towards the normal :math:`\hat{\mathbf{n}}` of the depth
contours (pointing to shallower water), by Snell's law with the shallow-water
celerity :math:`c = \sqrt{g w}`:

.. math::

    \sin\theta = \sqrt{\frac{w}{w_s}}\, \sin\theta_0,

with :math:`\theta_0` the angle of the offshore direction to
:math:`\hat{\mathbf{n}}`; waves travelling away from the shore are not turned.
:math:`\hat{\mathbf{n}}` comes from a least-squares plane fit of the
(optionally smoothed) depth in the local tangent plane, synchronised from the
owners to the ghost nodes; the upwind coefficients are then recomputed every
step. On a planar beach with straight contours the direction follows Snell's
angle to round-off (benchmark). Caustics are smoothed by the upwind scheme and
diffraction is not represented.


Coastal erosion
---------------

With the wave energy :math:`E` and the erodibility :math:`K` of the exposed
layer (stratigraphic multiplier times the lithology and duricrust factor), a
marine cell at depth :math:`w` erodes over a step

.. math::

    \Delta z = r\, E\, e^{-w/z_0}\, K\, \Delta t ,

and a cliff cell (land no higher than :math:`h_c` above sea level, next to the
sea) erodes :math:`r\, E_{sea}\, K\, \Delta t`, with :math:`E_{sea}` the largest
energy of its sea neighbours, at most down to sea level; both are capped per
step (``max_erosion``). This is the depth-decaying wave erosion of Pastier et
al. (2019) on a 2-D mesh, with the shoreward dissipation carried by the wave
energy itself. The eroded material is removed with the stratigraphic erosion
and added to the marine cascade of the step with its composition. With
``rubble`` the eroded reef volume :math:`V_r` joins the rubble lithology at the
same solid volume, :math:`V_r (1 - \phi_{0,reef}) / (1 - \phi_{0,rubble})`. Benchmarks:
the abrasion of a seafloor at three depths in uniform wave energy is exact to
:math:`10^{-15}` m; a low cliff retreats by :math:`r E \Delta t` and stops at
sea level; nothing erodes in the lee of an island.


Pelagic rain
------------

A carbonate rain :math:`R_c` (export production, uniform or a map) and a clay
rain :math:`R_k` settle on the marine cells each step. The fraction of the
carbonate preserved at water depth :math:`w` is

.. math::

    f(w) = \begin{cases} 1 & w \le z_{lys} \\
                         \dfrac{z_{ccd} - w}{z_{ccd} - z_{lys}} & z_{lys} < w < z_{ccd} \\
                         0 & w \ge z_{ccd} \end{cases}

with :math:`z_{lys}` the lysocline and :math:`z_{ccd}` the carbonate
compensation depth, both given per event. The settled thickness over a step
is :math:`(R_c f(w) + R_k)\, \Delta t`, capped at the water depth; the dissolved
volume :math:`R_c (1 - f)\, \Delta t` is kept in the budget. A deposit thinner
than the stratigraphic floor (0.1 mm) is carried over to the next steps.

**Clay decreasing away from the coast** (``clay_decay``: :math:`L`). The clay
rate is the coastal rate times :math:`c`, the solution of

.. math::

    (I - L^2 \nabla^2)\, c = 0 \ \text{in the sea}, \qquad c = 1 \ \text{on land},

which is :math:`e^{-x/L}` at a distance :math:`x` from a straight coast and
follows the coastline around headlands and islands. On a lattice of spacing
:math:`\Delta x` the discrete solution is exactly :math:`\lambda^i` with
:math:`\cosh(\Delta x / L') = 1 + \Delta x^2 / (2 L^2)`. The operator is rebuilt
each step, as the coastline moves.

**Resuspension** (``resuspend``). The settled thickness is added to the input
of the marine diffusion of the same step instead of being laid in place; the
carbonate mud and the clay keep their lithology through the diffusion (their
own columns of the lithology transport) and are counted in the inputs of the
marine composition, so both are conserved with the rest of the marine
deposit.


Link to the geochemistry
------------------------

*Coastal plume* (``solute`` response). With :math:`q_k` the river dissolved
load of species :math:`k` delivered at the coastal exits per unit area, the
plume is :math:`(I - L^2 \nabla^2)\, c_k = q_k` with land as identity rows
(the operator of the clay decay); the response factor is read at :math:`c_k`.

*Supply-limited burial* (``carbonate_supply``). With :math:`R_k` the marine
reservoir of the species (accumulated from the rivers by the geochemistry) and
:math:`y` the yield, the potential burial of a step :math:`B` (benthic
production plus settled pelagic carbonate) is scaled by

.. math::

    s = \min\left(1, \frac{y R_k}{B}\right),

and the reservoir loses what is buried, :math:`R_k \leftarrow R_k - s B / y`.
Over many steps the burial is bounded by the supply from weathering; when the
supply exceeds the potential production, the reservoir accumulates.


Numerics and parallel consistency
---------------------------------

Production and the pelagic rain are computed node by node, so they do not
depend on the number of processors. The open-water index and the clay decay
are linear solves on the mesh (``fgmres`` with block Jacobi), assembled on
owned rows; the index operator is built once, the clay operator every step.
The per-lithology routing adds one solve per lithology to each river and
cascade routing solve. With one lithology (dual lithology) none of this runs
and the results are unchanged bit for bit.


Benchmarks
----------

``benchmarks/test_carbonate.py``:

.. list-table::
   :header-rows: 1
   :widths: 26 44 30

   * - Test
     - Exact solution
     - Result
   * - Keep-up / give-up
     - half a metre above :math:`w^*` the platform reaches sea level and keeps
       up; half a metre below it drowns
     - both sides as predicted
   * - Drowning
     - depth increases at :math:`\sigma` once production has stopped
     - rate within 2 %
   * - Time-step independence
     - ODE :math:`dw/dt = -P(w)` (tight-tolerance reference)
     - :math:`10^{-4}` at 1 and 10 kyr
   * - Open-water index
     - parabolic ridge, :math:`E = 2 b L^2`
     - :math:`2 \times 10^{-6}` at two spacings and lengths
   * - Clay decay
     - 1-D strip, :math:`\lambda^i` and :math:`e^{-x/L}`
     - :math:`10^{-10}`; second order in :math:`\Delta x/L`
   * - Coastal erosion
     - abrasion at three depths, cliff retreat, island lee
     - :math:`10^{-15}` m; nothing in the lee
   * - Wave refraction
     - planar beach, Snell's angle at every depth
     - :math:`10^{-16}`
   * - Wave energy
     - strip :math:`(1 + \kappa\Delta x)^{-i}`; island shadow; two directions
     - :math:`10^{-14}`; nothing in the lee (``benchmarks/test_waves.py``)

``tests/test_lithologies.py`` and ``tests/test_carbonate.py`` check the
lithology split against dual lithology, the conservation of every lithology,
the carbonate and pelagic budgets, and the agreement of 1 and 3 processors.
