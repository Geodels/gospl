.. _ice:

==============================================
Ice accumulation and meltwater
==============================================

Mass-balance proxy
------------------

goSPL implements a first-order, mass-balance proxy for glacial ice rather than full ice-sheet dynamics (no shear-thinning, no SIA flow). The accumulation source at each cell is a fraction of the local precipitation determined by elevation relative to two time-varying altitudes — the equilibrium-line altitude (``hela``) and the ice-cap altitude (``hice``):

.. math::

   \mathrm{r_i} = \begin{cases}
     \mathrm{P_i} \cdot 1 & \mathrm{\eta_i \ge h_{ice}} \\
     \mathrm{P_i} \cdot \frac{\eta_i - h_{ela}}{h_{ice} - h_{ela}} & \mathrm{h_{ela} < \eta_i < h_{ice}} \\
     \mathrm{P_i} \cdot \frac{\eta_i - h_{ela}}{h_{ice} - h_{ela}} \cdot m_f & \mathrm{\eta_i \le h_{ela}}
   \end{cases}

Above the ice cap altitude the cell accumulates the full local precipitation as ice; between the ELA and the ice cap the fraction ramps linearly from 0 to 1; below the ELA the source becomes negative (ablation) and is amplified by ``melt`` (default ``10``) to make it act as a strong sink in the implicit solver.

A degenerate-configuration guard returns immediately when ``hice <= hela`` or when the maximum surface elevation lies below the ELA.

Ice flow accumulation
---------------------

Ice flux is routed via a multiple-flow-direction (MFD) accumulation on the **epsilon-filled** digital elevation model, with a configurable number of receivers per cell (``icedir``, 1 by default, up to 8). The MFD matrix is built the same way as the river flow matrix but on the filled surface so that ice does not pond.

The implicit drainage area is then solved via PETSc to obtain the ice flow field :math:`\mathrm{F^{ice}}`:

.. math::

   \mathrm{W^{ice}\; F^{ice}} = \mathrm{r}\; \mathrm{\Omega}

Two clamps are applied to the solution: negative values (cells where ablation exceeded upstream supply) are set to zero, and cells below the glacier terminus elevation (``hterm``) are also zeroed.

The clamped field is finally passed through a single linear-diffusion smoothing pass (``_hillSlope(smooth=1)``) to suppress the spiky pattern that MFD on a filled DEM tends to produce.

.. warning::

   The smoothing step is **not mass-conservative** — the published ``iceFA`` field differs slightly from the strict accumulated source. It is intended for visualisation and to give the erosion driver a regular field to integrate against; do not quote it as a flux balance.

Meltwater re-injection into the river network
---------------------------------------------

Without coupling, the clamp at line :math:`\mathrm{F^{ice}_i \ge 0}` silently discards the meltwater generated when ice from upstream reaches sub-ELA cells, and rivers downstream of glaciers therefore under-predict discharge. goSPL closes this loop by capturing the local ablation potential at every cell where ice is actually present:

.. math::

   \mathrm{m_i} = \begin{cases}
     \mathrm{P_i \cdot \frac{h_{ela} - \eta_i}{h_{ice} - h_{ela}}} & \mathrm{if}\; \eta_i < h_{ela}\; \mathrm{and}\; F^{ice}_i > 0 \\
     0 & \mathrm{otherwise}
   \end{cases}

This meltwater field is computed *before* the smoothing pass so the ice-presence gate uses the true (un-smeared) accumulation. The ``meltfac`` amplifier is deliberately **not** applied to :math:`\mathrm{m_i}` because it is a numerical sink-strengthening trick inside the implicit solver, not a physical melt multiplier.

The river accumulation step in ``flowAccumulation`` then adds :math:`\mathrm{m_i}` to the river source term (after subtracting the above-ELA precipitation already routed to ice), so cells downstream of glacier termini see the corresponding meltwater discharge.

.. note::

   The estimate is *local* — it is bounded by the local ablation potential but not strictly capped by the available upstream ice supply. It is a strict improvement over the previous zero-melt behaviour but not a strict mass closure. In practice, errors are small except where ablation potential drastically exceeds local ice flux.

Ice thickness and flexure coupling
----------------------------------

When the ice option is enabled, an empirical Bahr-style width-area scaling is applied to the smoothed ice-flow field to produce a per-cell ice thickness:

.. math::

   \mathrm{H^{ice}_i} = \mathrm{e \cdot a \cdot (F^{ice}_i)^{0.3}}

where :math:`\mathrm{e}` and :math:`\mathrm{a}` are the YAML parameters ``eheight`` and ``fwidth``. The result is smoothed again and written to the output as ``iceH``.

When ``flexure`` is also enabled, this thickness is used as the ice load contribution to the flexural isostatic computation.

Ice-driven erosion
------------------

Ice flux contributes to bedrock erosion in parallel with river-driven erosion in all SPL flavours (``SPL``, ``nlSPL``, ``soilSPL``). The additional coefficient takes the same stream-power form:

.. math::

   \mathrm{K^{ice}_{b,i}} = \mathrm{K_{ice} \cdot \Delta t \cdot (F^{ice}_i)^m}

controlled by ``Ki`` in the YAML ``ice`` block (default ``0``, so ice-driven erosion is off unless explicitly enabled).

Output and partition-edge artefacts
-----------------------------------

``iceFA`` and ``iceH`` are written to the HDF5 output whenever ``iceOn`` is true. The on-disk ``iceFA`` field is **floored at 1.0** to suppress sub-unit partition-edge noise produced by the linear-diffusion smoothing pass; the in-memory ``iceFAL`` (used for meltwater capture and erosion) is unchanged. ``iceH`` is now populated whenever ice is on, not only when flexure is enabled.
