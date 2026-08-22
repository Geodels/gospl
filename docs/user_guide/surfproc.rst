.. _surfproc:

=================================
Surface processes  parameters
=================================

Stream Power Law parameters
---------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
                
        **Declaration example**:

        .. code:: yaml

            spl:
                K: 3.e-8
                d: 0.42
                m: 0.4
                n: 1.0
                G: 1.

        This part of the input file define the parameters for the fluvial surface processes based on the *Stream Power Law* (SPL) and is composed of:

        a. ``K`` representing the erodibility coefficient which is scale-dependent and its value depend on lithology and mean precipitation rate, channel width, flood frequency, channel hydraulics. It is used in the SPL law: :math:`E = K (\bar{P}A)^m S^n`

        The following parameters are **optional**:

        b. Studies have shown that the physical strength of bedrock which varies with the degree of chemical weathering, increases systematically with local rainfall rate. Following `Murphy et al. (2016) <https://doi.org/10.1038/nature17449>`_, the stream power equation could be adapted to explicitly incorporate the effect of local mean annual precipitation rate, P, on erodibility: :math:`E = (K_i P^d) (\bar{P}A)^m S^n`. ``d`` (:math:`d` in the equation) is a positive exponent that has been estimated from field-based relationships to 0.42. Its default value is set to 0.0
        c. ``m`` is the flow accumulation coefficient from the SPL law: :math:`E = K (\bar{P}A)^m S^n` and takes the default value of 0.5.
        d. ``n`` is the slope coefficient from the SPL law: :math:`E = K (\bar{P}A)^m S^n` and takes the default value of 1.0.
        e. ``G`` dimensionless deposition coefficient for continental domain when accounting for sedimentation rate in the SPL following the model of `Yuan et al, 2019 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JF004867>`_. The default value is 0.0 (purely detachment-limited model).

        When ``n`` is **not** equal to 1.0 the SPL is solved with a non-linear PETSc ``SNES`` (the linear ``n == 1`` case uses a direct Krylov solve and needs none of the controls below). The transport-limited (``G > 0``) solver is the dominant cost of such a run, and its behaviour can be tuned (all optional):

        f. ``maxIter`` is the maximum number of non-linear iterations (default ``500``),
        g. ``rtol`` / ``atol`` are the relative / absolute convergence tolerances (default ``1.e-6``),
        h. ``pcType`` is the preconditioner for the ``ngmres`` Krylov solve (default ``'hypre'`` BoomerAMG),
        i. ``solver`` selects the primary non-linear solver: ``'qn'`` (default, limited-memory quasi-Newton / L-BFGS) or ``'ngmres'`` (accelerator + multigrid preconditioner).

        .. tip::

            ``solver: 'qn'`` is the default because on a global model it converged the transport-limited solve in well under 50 iterations versus ~200 for ``ngmres`` (each iteration is also cheaper), cutting the erosion phase by an order of magnitude at the **same** tolerance and solution. Whichever primary solver is chosen, goSPL automatically retries a stalled timestep with the complementary solver before continuing.


Hillslope and marine deposition parameters
-------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
                
        **Declaration example**:

        .. code:: yaml

            diffusion:
                hillslopeKa: 0.02
                hillslopeKm: 0.2
                nonlinKm: 500.0
                Gmar: 1.0
                clinSlp: 5.e-5

        Hillslope processes in goSPL is defined using a classical *diffusion law* in which sediment deposition and erosion depend on slopes (*simple creep*). The marine deposition of freshly deposited sediments by rivers is obtained using a non-linear diffusion and the following parameters can be tuned based on your model resolution:

        a. ``hillslopeKa`` is the diffusion coefficient for the aerial domain,
        b. ``hillslopeKm`` is the diffusion coefficient for the marine domain,
        c. ``nonlinKm`` is the transport coefficient of freshly deposited sediments entering the ocean from rivers (non-linear diffusion),
        d. ``Gmar`` is a dimensionless deposition coefficient for marine domain,
        e. ``clinSlp`` is the maximum slope of clinoforms (needs to be positive), this slope is then used to estimate the top of the marine deposition based on distance to shore.

        The following parameters tune the **marine** non-linear diffusion solver:

        f. ``tsSteps`` is the maximum number of internal time-steps the adaptive controller may take per goSPL step (default: 2000). Increase if the verbose log shows many rejected steps,
        g. ``offshore`` is the distance offshore (m) beyond which the clinoform-distance cap is no longer applied (default: 1.0e7),
        h. ``oFill`` is a depth (m, negative) **relative to sea level** below which the priority-flood algorithm is not applied — used to skip deep ocean cells (default: -6000.0). The cut-off elevation is ``sea level + oFill``, so with a sea level of 0 the value reads directly as an elevation, but on a model with a deep datum (e.g. ``sea: position: -2200``) an ``oFill`` of ``-1500`` means ``-3700`` m, not ``-1500`` m. Keep it well **below** sea level: the depression fill and the marine flow-direction surface are only smoothed at or above this cut-off, so a value above sea level leaves every bathymetric pocket on the shelf in place and sediment gets trapped at the coast.

        *Optional: marine/lake diffusion solver*

        .. code:: yaml

            diffusion:
                # ... above keys ...
                marineSolver: picard   # 'ts' (default) | 'picard'
                picardSub: 10
                picardIts: 2

        By default the marine and lake non-linear diffusion is integrated with an adaptive non-linear PETSc time-stepper (``marineSolver: ts``). On large, stiff marine inputs this can dominate the run (the adaptive controller takes many sub-steps and stalls at the diffusivity threshold). The opt-in ``marineSolver: picard`` uses a **lagged-diffusivity** backward-Euler scheme instead: it takes ``picardSub`` sub-steps over the goSPL step (default 10), freezing the diffusivity within each and solving the resulting *linear* system with ``picardIts`` Picard updates (default 2). Each solve is linear (no non-linear stalls or step rejections), which is markedly faster on production meshes. It is an **approximation** — on small problems it matches the default solver almost exactly, but on large runs the deposit geometry differs slightly; increase ``picardSub`` to converge toward the time-stepper solution. The default remains ``ts``.

        .. note::

           The ``picard`` solver has been **validated on a production earth model**: its marine deposit distribution was confirmed to match the default ``ts`` time-stepper closely (a few percent on the deposit geometry) while running roughly 5–8× faster on the marine diffusion. It is therefore a safe opt-in for large runs where marine diffusion dominates the wall-clock time; the exact adaptive ``ts`` solver remains the default.

        The following parameters control how partial-fill deposition in **inland depressions** is distributed (see :ref:`dep`):

        .. code:: yaml

            diffusion:
                # ... above keys ...
                nlPitVolume: 1.0e9
                nlPitDepth: 100.0
                nlPitK: 10.0
                pitInletBias: 0.10

        i. ``nlPitVolume`` is the volume threshold (m³) above which a partially-filled depression uses the bathymetric-pile + inlet-bias geometry instead of bottom-up fill (default: 1.0e9),
        j. ``nlPitDepth`` is the depth threshold (m) — both ``nlPitVolume`` AND ``nlPitDepth`` must be exceeded to trigger the diffusion path (default: 100.0),
        k. ``nlPitK`` is the non-linear diffusion coefficient (m²/yr) used inside selected pits (default: same as ``nonlinKm``). Higher values let the delta wedge prograde further per step,
        l. ``pitInletBias`` is the fraction (0–1) of each pit's deposit concentrated at the inlets to seed delta progradation; the remainder is distributed as a bathymetric bottom-up baseline. ``0.0`` = pure bowl fill, ``1.0`` = original inlet-only spike. Default: ``0.10``. See :ref:`dep` for the underlying algorithm.
                
        *Optional additions for non-linear diffusion model*

        A more complex version of the creep law involves a non-linear relationship between soil flux and topographic gradient. 
        
        .. note::
    
            Several non-linear creep-transport laws have been suggested in the literature and 2 non-linear formulations are available in goSPL. 

        **Either** by adding the following to the above parameters:

        .. code:: yaml

                hillslopenl: 2.5
                
        e. ``hillslopenl`` is the slope exponent in the non-critical hillslope model defined in the work of `Wang et al. (2024) <https://www.sciencedirect.com/science/article/pii/S0169555X24001053>`_. Here the coefficient of diffusion is set to the values of ``hillslopeKa`` and ``hillslopeKm``.

        **Or** by defining the following two parameters:

        .. code:: yaml

                hillslopeSc: 0.8
                hillslopeNb: 4

        f. ``hillslopeSc`` is the critical slope,
        g. ``hillslopeNb`` is the number of terms used in the truncated Taylor series formulation.

        In this last model, the non-linear creep formulation is described in `Barnhart et al. (2019) <https://gmd.copernicus.org/articles/12/1267/2019/gmd-12-1267-2019.pdf>`_ (section 3.4.3 **EQ. 14**).


Ice sheets and glacial erosion
-----------------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::

        Adding an ``ice`` section turns on goSPL's glacial model, driving glacial
        abrasion, till transport and ice loading. It is a cheap, robust
        **diagnostic glacial-erosion model**: the ELA net mass balance is routed
        downhill into an ice discharge, from which a Bahr ice thickness and a
        bounded Glen-sliding velocity are derived (one linear solve, no
        ice-dynamics time integration), then abrasion, till/moraine deposition,
        meltwater and an ice load follow. It is fast and physical at any
        resolution and over goSPL's long timesteps, suited to the **morphology of
        glacial erosion**. The full algorithm is in the technical guide
        (:ref:`ice`).

        **Declaration example**:

        .. code:: yaml

            ice:
                # melt_conserve: True       # discharge-conserving river meltwater (default)
                # Either constant glacial parameters
                hterm: 1700.0
                hela: 1850.0
                hice: 2100.0
                # Or using a file to characterise glacial evolution
                # evol: 'data/ice_evol.csv'
                # Diagnostic controls (all optional):
                icedir: 1                 # MFD flow directions for ice routing
                eheight: 0.25             # Bahr thickness factor
                fwidth: 1.5               # Bahr width factor
                melt: 1.0                 # ablation in the net mass balance (0 = accumulation-only; 1 = true net)
                slide: 1.0e-3             # basal sliding coefficient (Glen sliding law)
                glen: 3.0                 # Glen sliding exponent n
                # accum_factor: 1.0       # precip->ice accumulation fraction (optional)
                # accum_max: 2.0          # cap accumulation rate, m ice/yr (optional)
                abrasion:
                    Kg: 1.0e-4
                    l: 1.0
                    # Kl: 0.0          # lateral (valley-wall) erosion coeff -> U-shaping
                    # lat_l: 1.0        # lateral velocity exponent (defaults to l)
                till:
                    on: True          # default True (set False -> abrasion goes straight to rivers)
                    route: True       # default True (catchment-routed; False = global melt-spread)

        Each goSPL step the model routes the ELA **net mass balance**
        (accumulation minus the ``melt``-scaled ablation) downhill on a
        **terminus-anchored, drainage-conditioned bed** — filled in parallel so
        that every cell above the glacier terminus strictly drains (no closed
        depression, no flat), which is required for the single ice-discharge
        solve to be well-posed — into an ice discharge :math:`Q`, derives a Bahr
        ice thickness :math:`H = \mathrm{eheight}\cdot\mathrm{fwidth}\cdot Q^{0.3}`
        and a bounded basal sliding velocity from Glen's sliding law, and then
        drives the abrasion / till / loading machinery below — all in a single
        linear solve, with no ice-thickness PDE time integration. The fill is
        partition-invariant, so results are identical at any processor count.

        ``melt_conserve`` (default ``True``) sets how glacial meltwater is
        delivered to the rivers. ``True`` is **discharge-conserving**: the
        precipitation that fell as ice above the ELA is routed down-glacier and
        released as liquid meltwater where the ice melts out, so the total
        meltwater equals the total accumulation — the right assumption over
        goSPL's long (steady-state) timesteps, and it closes the glacial water
        budget so downstream basins don't under-predict discharge. ``False``
        reverts to the local precipitation-scaled ablation rate (cheaper, but it
        generally returns less water than fell as ice).

        The equilibrium-line / ice-cap geometry controls where ice accumulates
        and melts:

        a. ``hterm`` is the glacier terminus elevation (m) — no ice is kept below it. The effective floor is ``max(hterm, sea level)``: ice never persists below the (possibly time-varying) sea surface, and a ``hterm`` below sea level is raised to it. When omitted, the terminus defaults to the **sea-level position**,
        b. ``hela`` is the equilibrium-line altitude (m) — ablation below, accumulation above,
        c. ``hice`` is the ice-cap altitude (m) — full precipitation is captured as ice above it.

        The diagnostic flow controls are all optional, with the defaults shown
        above:

        d. ``icedir`` is the number of MFD flow directions used to route the ice accumulation into the discharge :math:`Q`,
        e. ``eheight`` and ``fwidth`` are the Bahr thickness and width scaling factors in :math:`H = \mathrm{eheight}\cdot\mathrm{fwidth}\cdot Q^{0.3}`. Only their **product** sets the thickness, so tune ``eheight * fwidth`` (the split between them has no separate effect),
        f. ``melt`` controls how much the below-ELA **ablation** is subtracted from the accumulation when forming the discharge — i.e. how far/thick glacier tongues reach. ``0`` (default) ignores ablation (accumulation-only, the legacy behaviour); ``1`` is the true net balance; ``> 1`` amplifies ablation so tongues are shorter and thinner. It changes only the ice extent — the raw ablation still drives the till melt-out and meltwater,
        g. ``slide`` is the basal-sliding coefficient (Glen sliding law) and ``glen`` the Glen sliding exponent :math:`n` (usually 3), setting the bounded basal velocity :math:`u_b \propto H^{n-1}|\nabla s|^{n-1}\nabla s` that drives abrasion,
        h. ``accum_factor`` (default ``1.0``) and ``accum_max`` (default unset) control the **accumulation** part of the surface mass balance only (ablation is untouched). Full precipitation is rarely all snow/ice, so ``accum_factor`` is a precipitation→ice conversion fraction and ``accum_max`` caps the accumulation rate (m ice/yr) at a realistic ceiling (real ice sheets accumulate ~0.1–2 m/yr). **Set these for high-precipitation runs** — converting several m/yr of rainfall directly to ice produces unphysically thick ice.

        The ``abrasion`` sub-block enables velocity-based glacial erosion
        :math:`E_g = K_g\,|u_b|^{l}` (off by default, ``Kg: 0``):

        i. ``Kg`` is the (vertical) abrasion coefficient (default ``0.0`` — set it to enable glacial erosion),
        j. ``l`` is the basal-sliding-velocity exponent (default ``1.0``),
        k. ``Kl`` is the **lateral** (valley-wall) erosion coefficient (default ``0.0`` = off). Vertical abrasion only deepens the trough; ``Kl > 0`` adds erosion of the *walls* flanking fast ice — each wall cell (little ice of its own) is abraded at ``Kl·u_b,neighbour^{lat\_l}``, tapered by how much of the wall is in contact with the neighbouring ice column — which **widens glaciated valleys toward a U-profile**. The eroded wall rock joins the same conserved till → moraine budget. ``lat_l`` (default = ``l``) is its velocity exponent.

        The ``till`` sub-block controls glacial sediment (on by default, but
        inert until abrasion is enabled with ``Kg > 0``):

        l. ``on`` (default ``True``) — abraded rock is carried as **till** and
           deposited as a moraine where the ice melts out (the ablation zone),
           conserving the abraded volume. With stratigraphy on, the till is
           layered into the stratigraphic record and split into the coarse/fine
           lithology fractions when dual lithology is enabled. Set ``False`` to
           instead send abrasion straight into the fluvial sediment system.
        m. ``route`` (default ``True``) — controls how the till is distributed.
           ``True`` **routes the till down the ice-surface flow network** and
           melts it out toward each catchment's terminus, so deposition stays
           connected to the upstream erosion (correct on multi-glacier / global
           domains). ``False`` instead spreads the **global** abraded volume
           across the whole ablation zone weighted by the meltwater rate — cheaper
           (no extra solve) but it decouples erosion and deposition across
           separate ice masses, so prefer it only on a single-ice-mass domain.
           Both conserve mass.

        The glacier geometry can instead be read from a file:

        n. ``evol`` is the glacier characteristics over time (`csv` file). When
           used, ``hterm``, ``hela`` and ``hice`` are not required because they
           are defined in this file.

        When flexural isostasy is enabled, the diagnostic ice thickness is
        automatically used as the ice load contribution to the isostatic
        computation — no extra parameters are needed.

        .. tip::

            **Tuning the glacial model.** A few relationships are worth knowing
            before turning the knobs:

            - **Ice thickness is one knob, not two.** :math:`H` depends only on
              the *product* ``eheight * fwidth`` (in
              :math:`H = \mathrm{eheight}\cdot\mathrm{fwidth}\cdot Q^{0.3}`), so
              raise/lower that product until the ``iceH`` output matches the
              thickness you expect (~10\ :sup:`2`–10\ :sup:`3` m); splitting it
              between the two factors has no separate effect.
            - **Accumulation.** ``accum_factor`` is a *fraction* of precipitation
              that becomes ice — keep it **≤ 1** (it is not an amplifier; values
              > 1 are unphysical). Use ``accum_max`` (m ice/yr, ~0.1–2 for real
              ice) as the realistic ceiling, and raise the *precipitation forcing*
              or the Bahr product — not ``accum_factor`` — if you need thicker
              ice.
            - **Glacier extent / melt.** ``melt`` sets how much ablation eats the
              ice discharge: ``0`` ignores it (accumulation-only), ``1`` is the
              true net balance, ``> 1`` gives shorter/thinner tongues. Use it to
              control how far ice reaches below the ELA; it does not change the
              accumulation-zone thickness.
            - **Erosion: deepen vs widen.** ``Kg`` (vertical abrasion)
              **deepens** the trough; ``Kl`` (lateral wall erosion) **widens** it
              toward a U-profile. Keep ``Kl`` the *same order of magnitude* as
              ``Kg`` (e.g. both ~``1e-3``): an ``Kl`` orders of magnitude larger
              than ``Kg`` makes the walls retreat far faster than the floor
              incises (runaway widening). Increasing ``Kg`` alone deepens — it
              does not widen.

        .. important::

            The glacial evolution file is defined as a 4 columns **csv** file containing in the first column the time in years (it doesn't need to be regularly temporally spaced) and in the second the glacier characteristics for the given time. When goSPL interprets this file, it will interpolate linearly between the defined times to find the values of ``hterm``, ``hela`` and ``hice`` for every time step.

    .. grid-item-card::

        **Spatially-varying ELA (global models)**

        In a **global** run a single ``hela`` cannot be right everywhere — the
        equilibrium-line altitude is ~5000–6000 m in the tropics but near sea
        level at the poles. Each of ``hela``, ``hice`` and ``hterm`` can
        therefore be given as a **per-vertex map** (``[file, key]``) instead of a
        scalar — the same convention as the precipitation maps:

        .. code:: yaml

            ice:
                hela:  ['input/ela', 'ela']    # per-vertex ELA (e.g. latitude-varying)
                hice:  ['input/ela', 'hice']   # ice-cap altitude must track hela
                hterm: 0.                       # a global scalar floor is usually fine
                slide: 1.0e-3
                glen: 3.0

        Only ``hela`` really needs a map; ``hice`` should follow it (it is the
        top of the accumulation band, so ``hice > hela`` everywhere), while
        ``hterm`` is a backstop floor that can stay a global scalar.

        The geometry can also vary **in time** through a ``glaciers`` time series
        (mirroring the precipitation ``climate`` block) — each entry has a
        ``start`` time and uniform-or-map ``hela``/``hice``/``hterm``, stepped as
        the simulation advances. This gives a latitude-**and**-time varying ELA,
        e.g. for glacial–interglacial cycles:

        .. code:: yaml

            ice:
                slide: 1.0e-3
                glen: 3.0
                glaciers:
                  - start: -120000.
                    hela:  ['input/ela_lgm', 'ela']
                    hice:  ['input/ela_lgm', 'hice']
                    hterm: 0.
                  - start: -20000.
                    hela:  ['input/ela_holocene', 'ela']
                    hice:  ['input/ela_holocene', 'hice']
                    hterm: 500.

        .. note::

            The map files follow the standard goSPL ``.npz`` convention: the
            ``key`` field is a per-vertex array over the mesh. A uniform scalar
            and the ``evol`` CSV remain available and unchanged; ``evol`` takes
            precedence over maps if both are given.

        **Deriving the ELA from paleo-climate temperature.** The ``gospl-ela``
        command (see :ref:`Running goSPL <running>`) turns a per-vertex
        temperature map into the ``hela``/``hice`` ``.npz`` maps above by
        lapse-rate inversion (``ELA = z + (T - T_ELA)/Gamma``). It derives the
        ELA *position* only; the ablation magnitude stays precipitation-scaled
        (it is not a degree-day melt model).

        From a **terminal** (run once per climate snapshot; ``--start`` prints a
        ready-to-paste ``glaciers`` entry):

        .. code:: bash

            gospl-ela \
                --temperature climate/t2m_21ka.npz --t-key t2m \
                --reference surface --elevation input/mesh.npz --z-key z \
                --lapse 0.0065 --t-ela -2.0 --band 400 \
                --out input/ela_21ka.npz --start -21000

        ``--help`` lists every option; loop over snapshots with a shell ``for``.
        (``python -m gospl.tools.ela_from_temperature ...`` works without the
        console-script install.)

        From a **Jupyter notebook**, either shell out with the ``!`` magic:

        .. code:: bash

            !gospl-ela \
                --temperature climate/t2m_21ka.npz --t-key t2m \
                --reference surface --elevation input/mesh.npz --z-key z \
                --lapse 0.0065 --t-ela -2.0 --band 400 --out input/ela_21ka.npz

        or import the conversion and build the ``glaciers`` list directly:

        .. code:: python

            import numpy as np
            from gospl.tools.ela_from_temperature import derive_ela

            z = np.load("input/mesh.npz")["z"]
            glaciers = []
            for t, f in [(-21000, "climate/t2m_21ka.npz"),
                         (-18000, "climate/t2m_18ka.npz")]:
                T = np.load(f)["t2m"]
                hela, hice = derive_ela(T, lapse=0.0065, t_ela=-2.0, band=400.0,
                                        reference="surface", elevation=z)
                out = f"input/ela_{int(-t/1000)}ka.npz"
                np.savez(out, hela=hela, hice=hice)
                glaciers.append({"start": float(t),
                                 "hela": [out[:-4], "hela"],
                                 "hice": [out[:-4], "hice"]})


Soil production, erosion, transport and deposition
-----------------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            soil:
                mode: 'lumped'
                soilK: 4.e-6
                maxProd: 50.e-6
                depthProd: 0.5
                roughnessL: 1.0
                decayDepth: 0.7
                bedrockConv: 0.0001
                uniform: 0.5
                soilMap: ['test_mesh8/hsoil', 'soil']
                tempMap: ['test_mesh8/temperature', 'temp']
                activation: 40.e3
                tempRef: 15.0

        a. ``soilK`` is the erodibility coefficient for soil,
        b. ``maxProd`` is the soil production maximum rate (m/yr),
        c. ``depthProd`` is the soil production decay depth (m),
        d. ``roughnessL`` is the roughness length scale,
        e. ``decayDepth`` is the soil transport decay depth for non-linear diffusion where the coefficient of diffusion is set to the values of ``hillslopeKa`` and ``hillslopeKm``,
        f. ``bedrockConv`` is the soil to bedrock conversion fraction, bedrock begins where soil production is a very small fraction of the maximum soil production (optional, default ``0.0001``). It also sets the **maximum soil thickness** to ``-ln(bedrockConv) * depthProd`` (the depth at which production has decayed to that fraction; e.g. ~4.6 m for the default with ``depthProd = 0.5``). Setting ``bedrockConv: 0`` means *no* bedrock-conversion depth and therefore **no maximum-soil cap** (soil thickness is unbounded).
        g. ``mode`` selects how ``soil`` (the near-surface layer) is accounted (optional, default ``'lumped'``):

           - ``'lumped'`` — the soil layer is a soft surface cover that absorbs **both** weathering-produced regolith **and** deposited sediment (fluvial, lake/pit, marine). This is the historical behaviour (unchanged).
           - ``'regolith'`` — the soil layer is the **weathering-produced regolith only**; deposited sediment is instead kept in the stratigraphy, where freshly deposited layers are given a soft erodibility (``stratK = soilK/K``, i.e. fresh sediment erodes like soil). Underwater (marine or ponded lake) there is no soil, and under ice the regolith is frozen (preserved).

        .. important::

            ``mode: 'regolith'`` needs **stratigraphic recording turned on** so the deposited
            sediment has somewhere to live — i.e. a stratal time step ``strat`` in the
            ``time`` block (which sets the number of stratigraphic layers). No ``strata:``
            block is required (that is only for *initial* layers). Without a ``strat`` time
            step goSPL prints a warning and freshly deposited sediment erodes at bedrock
            erodibility.

        Then the user can specify the initial soil thickness if any by setting **either**:

        h. ``uniform`` a uniform soil thickness on the entire surface (m),

        **or**:

        i. ``soilMap`` a soil thickness map given as ``[file, key]``.

        .. important::

            When defining a soil thickness grid, one needs to use the **npz** format and needs to specify the key corresponding to the soil thickness value in the file. In the ``soilMap`` example above the file is ``test_mesh8/hsoil.npz`` and this key is ``'soil'``. The soil grid needs to define values for all vertices in the mesh in metres.

        Soil production can optionally be made **temperature-dependent** via an
        Arrhenius scaling of ``maxProd`` (warmer ⇒ faster weathering). This is
        activated by supplying an annual-mean surface-temperature map:

        j. ``tempMap`` a temperature map given as ``[file, key]`` — an **npz** file whose ``key`` holds the annual-mean surface temperature (in **degrees Celsius**) at every mesh vertex. When present, the production rate becomes ``maxProd * exp( Ea/Rg * (1/T_ref - 1/T) )`` (temperatures internally converted to Kelvin; ``Rg = 8.314`` J/mol/K). When omitted, production uses the constant ``maxProd`` everywhere.
        k. ``activation`` the Arrhenius activation energy ``Ea`` (J/mol, optional, default ``40.e3``); only used when ``tempMap`` is set.
        l. ``tempRef`` the reference temperature (degrees Celsius, optional, default ``15.0``) at which the production rate equals ``maxProd``; only used when ``tempMap`` is set.

        The soil-aware non-linear SPL is solved with a PETSc ``SNES``. Its
        behaviour can be tuned (all optional) with:

        m. ``maxIter`` is the maximum number of non-linear iterations (default ``500``),
        n. ``rtol`` / ``atol`` are the relative / absolute convergence tolerances (default ``1.e-6``),
        o. ``pcType`` is the preconditioner for the ``ngmres`` Krylov solve (default ``'hypre'`` BoomerAMG; ``'gamg'``, ``'bjacobi'`` or ``'asm'`` can help on heavily-decomposed / ocean-dominated partitions),
        p. ``solver`` selects the primary non-linear solver: ``'qn'`` (default, limited-memory quasi-Newton / L-BFGS) or ``'ngmres'`` (accelerator + multigrid preconditioner).

        .. tip::

            The soil solve is usually the dominant cost of a soil-enabled run. The default ``solver: 'qn'`` was chosen because on a global model it cut the soil-solve wall time roughly in half (or more) versus ``'ngmres'`` at the **same** tolerance and solution — L-BFGS reaches a comparable iteration count but each iteration is far cheaper than an ``ngmres``/multigrid sweep. If ``'qn'`` struggles on a particular configuration, set ``solver: 'ngmres'``. Prefer switching ``solver`` over relaxing ``rtol``: a loose ``rtol`` can leave the elevation field under-resolved and **destabilise the downstream sediment-routing solver**, which is both slower and less accurate.

        .. tip::

            **Validated fast preset for a stiff regional soil run.** Where the default ``'qn'`` stalls, ::

                soil:
                    solver: 'ngmres'
                    pcType: 'hypre'
                    rtol: 1.e-5
                    atol: 1.e-5

            converged cleanly in ~30 non-linear iterations and cut the soil solve by roughly 3–5× versus the un-tuned ``ngmres`` at ``rtol = 1.e-6``, with no artefacts in the elevation or deposition fields. Two findings worth keeping in mind when tuning from here: do **not** loosen ``rtol``/``atol`` to ``1.e-4`` — at that tolerance the global residual norm stops while a single node remains under-resolved, producing an isolated elevation **spike** ("needle"); and at ``1.e-5`` the stronger ``pcType: 'hypre'`` BoomerAMG preconditioner *beat* the cheaper ``'bjacobi'`` on wall-clock (fewer Krylov iterations per ``ngmres`` step more than offset its higher per-application cost). The cheaper ``'bjacobi'`` only wins at the over-loose ``1.e-4`` that should not be used. So: tighten the tolerance to ``1.e-5`` and keep the stronger preconditioner.

        .. note::

            If the chosen primary solver stalls on a stiff soil-production residual, goSPL automatically retries that timestep with the *complementary* solver (quasi-Newton ⇄ ``ngmres`` multigrid accelerator) before continuing. If **both** solvers still diverge — which happens when a node whose duricrust/armour has just eroded off suddenly captures a large drainage, so the demanded single-step incision becomes very large and the slope term overshoots — the fluvial solve is retried with **adaptive sub-stepping**: the step is split into ``N = 4 → 8 → 16`` increments of :math:`\Delta t/N` (each cuts the demanded per-solve incision by ``N`` until the residual is smooth enough to converge), and the increments sum to the full step's erosion. Erosion is therefore *retained*, not skipped; only a step still diverging after ``N = 16`` is reverted to its prior elevation (no fluvial erosion that step) as a last resort. This makes the solve robust to drainage-capture / de-armouring transients that would otherwise spike the topography.


Sediment surface erodibility factor
-------------------------------------


.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            sedfactor:
                - start: 200000.
                  uniform: 3
                - start: 400000.
                  map: ['facEro','fsed']

        One could choose to impose variable erodibility factors through space and time to reflect different surficial rock composition. For example, those maps could be set to represent different rock erodibility index as proposed in `Mossdorf et al. (2018) <https://www.sciencedirect.com/science/article/abs/pii/S0143622817306859>`_. The factor are then used in front of the erodibility coefficient (``K`` in the SPL).

        .. important::

            When defining your variable erodibility factors grid, you needs to use the **npz** format and your factors would be specified by a key corresponding to the factor values for each vertice of the mesh. In the above example this key is ``'fsed'``. 


Compaction & porosity variables definition
------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::  
        
        **Declaration example**:

        .. code:: yaml

            compaction:
                phis: 0.49
                z0s: 3700.0

        We assume  a depth-porosity relationship for the sediment compaction based on the following parameters:

        a. porosity at the surface ``phis``, default value is set to 0.49,       
        b. e-folding depth ``z0s`` (in metres), default value is set to 3700.       

        .. note::

            See the technical `documentation <https://gospl.readthedocs.io/en/latest/tech_guide/strat.html>`_ for more information.


Dual-lithology (coarse/fine) variables definition
--------------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::

        **Declaration example**:

        .. code:: yaml

            strata:
                dual: True
                coarse: {phi0: 0.49, z0: 3700.}
                fine:   {phi0: 0.63, z0: 1960., k_factor: 1.5}
                bedrock_coarse_frac: 0.6
                fine_diff_factor: 2.0
                pitInletBias: {coarse: 0.5, fine: 0.0}

        Optional. When stratigraphy is enabled (a positive ``strat`` interval in
        the ``time`` block), set ``dual: True`` to track coarse (sand) and fine
        (silt/clay) sediment separately. **Omitting the block — or ``dual:
        False`` — reproduces the single-fraction behaviour exactly.** The keys:

        a. ``coarse`` / ``fine`` — per-lithology depth-porosity curves
           (``phi0`` surface porosity, ``z0`` e-folding depth in metres). The
           coarse curve defaults to the ``compaction`` ``phis``/``z0s`` values.
        b. ``fine.k_factor`` — fine erodibility relative to coarse
           (``K_fine = K_coarse * k_factor``); ``1.0`` (default) = no contrast.
        c. ``fine_diff_factor`` — fine hillslope/soil diffusivity relative to
           coarse; a **multiplier** on the existing ``hillslopeKa`` /
           ``hillslopeKm`` / ``nonlinKm`` coefficients (it does not replace
           them). ``> 1`` makes fines diffuse (and so travel) farther.
        d. ``bedrock_coarse_frac`` — coarse fraction of the underlying bedrock
           (default 0.5), i.e. the composition of material eroded from bedrock.
        e. ``bedrock_sentinel`` — ``True`` inserts a dedicated infinite-bedrock
           reservoir *beneath* the initial layers supplied by an ``npstrata``
           file, with the ``bedrock_coarse_frac`` composition (default
           ``False``). Without it the deepest file layer is itself the
           un-erodable floor; with it those layers are finite and erosion that
           cuts through them exposes the bedrock reservoir below. Only relevant
           when an ``npstrata`` file is given — the no-file path always builds a
           sentinel.
        f. ``pitInletBias`` — per-fraction lake/depression deposition bias
           (coarse builds inlet deltas, fine settles in the depocenter). It is
           the **contrast** ``coarse − fine`` that sets how strongly the
           pit-deposit composition segregates with bathymetric depth: equal
           values give a uniform composition (no segregation), ``coarse > fine``
           (e.g. ``{coarse: 0.5, fine: 0.0}``) concentrates fine in the deep
           depocenter, and ``coarse − fine == 1`` is the maximum (fully
           depth-proportional) split. Conserves each pit's fine volume at any
           strength.

        Fines preferentially reach the **distal** parts of the system: the
        depocenter of lakes/depressions and the deeper/distal marine shelf,
        while coarse stays proximal.

        .. note::

            See the technical `documentation <https://gospl.readthedocs.io/en/latest/tech_guide/strat.html>`_ (Dual lithology section) for the physics and current limitations.


Sediment provenance tracers
-----------------------------------------------------

.. grid:: 1
    :padding: 3

    .. grid-item-card::

        **Declaration example**:

        .. code:: yaml

            provenance:
                classes: 3
                source: ['input/source', 'rock']    # per-vertex source class
                # uniform: 0                          #   ... or one class everywhere

        When stratigraphy is enabled (a positive ``strat`` interval in the
        ``time`` block), adding a ``provenance`` block carries **N source-rock
        classes** through erosion, transport, deposition and the stratigraphic
        record, so the composition of each layer (and hence per-pixel /
        per-basin provenance) is tracked conservatively. Provenance is a passive
        label — it has no effect on erodibility, diffusivity or deposition, so a
        model without the block is unchanged.

        The recorded composition is **conservation-exact for any number of
        classes** (the per-class thicknesses always sum to the layer thickness),
        and both depositional sinks carry the exact delivered source mix: marine
        deposits take the basin-delivered composition and intracontinental
        pit/lake deposits take each lake's cascade-retained mix (overspill
        between chained lakes is mixed exactly).

        a. ``classes`` — number of source classes,
        b. ``source`` — a per-vertex integer class map ``[file, key]`` (values
           in ``[0, classes)``), **or**
        c. ``uniform`` — a single source class everywhere.

        The per-layer composition is written to the stratal output (``stratP``).
        The optional **copper-fertility weighting** is *not* a model input: it
        is supplied at post-processing to the standalone
        :mod:`gospl.analyse.provenance` tool (its ``--cu-weights`` argument),
        which reads the recorded ``stratP`` and blends a Cu-sourced fraction —
        a what-if weighting best varied without re-running the model.

        .. note::

            See the technical `documentation <https://gospl.readthedocs.io/en/latest/tech_guide/provenance.html>`_ (provenance section) for the algorithm, the standalone post-processing tool, and the copper-prospectivity scope.


Groundwater & duricrust
-----------------------

Adding an optional ``groundwater`` section turns on a **water table** and a
generic **duricrust** — near-surface hydrology plus chemical armoring of the
erodibility. It is fully **opt-in**: with no ``groundwater`` block the model is
byte-identical to a run without it.

The physical picture: rainfall that does not run off **infiltrates** and feeds a
water table, solved each step as an implicit Dupuit–Boussinesq head on the mesh
(draining to the sea, lakes and rivers). Where the water-table depth sits in a
shallow **capillary fringe**, an indurated crust (calcrete / silcrete / ferricrete,
treated generically) precipitates and **hardens the surface**, so crusted cells
erode more slowly — producing relief inversion and, when stratigraphy is on,
stacked duricrusts that are buried and later re-exposed. The duricrust runs
**soil-independent** by default and **couples to** :ref:`soil <surfproc>`
production when soil is tracked.

.. grid:: 1
    :padding: 3

    .. grid-item-card::

        **Declaration example**:

        .. code:: yaml

            groundwater:
                Ksat: 3.65e4          # hydraulic conductivity K_h (m/yr)
                specific_yield: 0.1   # S (drainable porosity)
                aquifer_base: 50.0    # z_bed depth below surface (m)
                bedrock_depth: 0.0    # permeable rock below lHbed (from_soil only)
                min_sat_thickness: 1.0
                infiltration: 0.3     # fraction of (rain − evap) recharging
                conserve_baseflow: True
                picard_its: 3
                seepage_passes: 4
                duricrust:
                    form_rate: 1.0e-4     # k_form (m/yr at Φ=Ψ=1)
                    max_thickness: 5.0    # duriH_max (m)
                    fringe_depth: 3.0     # d0 — fringe centre below surface (m)
                    fringe_width: 2.0     # w — Gaussian half-width (m)
                    discharge_gate: False # crust only in gw discharge zones
                    supply_exp: 1.0       # p on (rain − evap) in the proxy
                    weather_Ea: 0.0       # Arrhenius activation energy (0 ⇒ off)
                    armor_max: 0.9        # max fractional K reduction (0..1)
                    armor_diffusion: False  # also armor hillslope Cd
                    break_rate: 1.0       # k_break per unit incision
                    decay_rate: 1.0e-6    # k_decay disequilibrium (1/yr)
                    weathering:
                        mode: proxy       # proxy | rate | prodsoil
                        C_eq: 1.0
                        Dw: 1.0
                        path_length: 20.0
                        weather_Ea: 0.0
                        weatherability: 1.0

The **water-table (hydrology) keys** are:

a. ``Ksat`` — saturated hydraulic conductivity ``K_h`` (m/yr); a scalar, a per-vertex map ``[file, key]``, or per-lithology.
b. ``specific_yield`` — drainable porosity ``S`` (the storage coefficient linking recharge to head change), default ``0.1``.
c. ``aquifer_base`` — depth of the impermeable base ``z_bed`` below the surface (m): a **scalar**, a per-vertex **map** ``[file, key]``, or the string ``from_soil`` (tie the base to the bedrock elevation ``z_bed = lHbed − bedrock_depth`` — requires soil tracking; in a depositional basin the base deepens to the bottom of the porous sediment fill). Default ``50.0``.
d. ``bedrock_depth`` — permeable weathered/fractured-rock thickness below ``lHbed`` (m), used **only** with ``aquifer_base: from_soil`` (default ``0``).
e. ``min_sat_thickness`` — floor ``b_min`` on the saturated thickness so the transmissivity stays positive near the base (m, default ``1.0``).
f. ``infiltration`` — fraction ``f_infil`` of ``max(0, rain − evap)`` that recharges the aquifer; a scalar or a per-vertex map ``[file, key]`` (default ``0.3``).
g. ``conserve_baseflow`` — return the seepage discharge to the river network so total river discharge stays ``≈ rain − evap`` (default ``True``). The infiltrated recharge leaves surface runoff and is **re-injected as baseflow** at the seepage nodes (rivers become baseflow-fed); also writes the ``baseflow`` output.
h. ``lake_exchange`` — opt-in lake ↔ aquifer **volume** coupling (default ``False``). When on, the signed across-bed groundwater flux debits/credits each lake's fill budget — a lake ringed by a higher water table gains groundwater, one ringed by a lower table leaks. Off ⇒ lakes are fixed-head only (unchanged).
i. ``subglacial_recharge`` — fraction of the glacial meltwater (``iceMeltRiverL``) that infiltrates the aquifer where the ice melts out (default ``0`` — under ice the rain path is gated off; this is the one recharge path allowed there).
j. ``fine_infil_factor`` — multiplier on ``f_infil`` for the fine end-member (dual lithology); ``< 1`` makes clay/fine surfaces infiltrate less than coarse/sand (default ``1`` — no lithology dependence).
k. ``infil_slope_ref`` — reference slope for a ``f/(1 + slope/infil_slope_ref)`` reduction of infiltration on steep terrain (default ``0`` — off; slope is the steepest-descent gradient).
l. ``picard_its`` / ``seepage_passes`` — inner iteration counts for the unconfined non-linearity ``T(h)`` and the seepage free-boundary discovery (defaults ``3`` / ``4``).

.. important::

    ``aquifer_base: from_soil`` needs :ref:`soil production <surfproc>` tracked (it
    reads the bedrock elevation ``lHbed``). Without soil it falls back to the
    surface as the base with a warning.

The nested **duricrust keys** (omit the ``duricrust:`` block for a water table with no crust) are:

a. ``form_rate`` — crust formation rate ``k_form`` (m/yr at full favourability and supply).
b. ``max_thickness`` — maximum crust thickness ``duriH_max`` (m); the induration degree is ``duriF = duriH/max_thickness``.
c. ``fringe_depth`` / ``fringe_width`` — centre ``d0`` and Gaussian half-width ``w`` (m) of the capillary-fringe favourability band ``Φ`` on the water-table depth.
d. ``discharge_gate`` — restrict crust formation to groundwater **discharge zones** (default ``False``). When ``True`` the favourability ``Φ`` is multiplied by ``G = (−∇·q)⁺ / ((−∇·q)⁺ + R)``, the fraction of the cell's upward discharge imported by **lateral convergence** of the groundwater flux (vs local recharge ``R``): ``G→1`` at a convergent valley floor / seepage face, ``G→0`` at a divergent recharge rise. This is the **absolute-accumulation** (lateral, valley/footslope ferricrete) style, as opposed to the default **relative-accumulation** style that indurates wherever the water table is shallow (the plateau/*bowal* cuirasses that blanket flat uplands). Off ⇒ weight ``1`` everywhere (unchanged).
e. ``supply_exp`` — exponent ``p`` on ``(rain − evap)`` in the default climate proxy supply.
f. ``weather_Ea`` — Arrhenius activation energy (J/mol) for an optional temperature scaling of the supply (``0`` ⇒ off; reuses the soil ``tempMap`` when present).
g. ``armor_max`` — maximum fractional erodibility reduction (0–1); a fully indurated cell (``duriF = 1``) has its ``K`` multiplied by ``1 − armor_max`` (e.g. ``0.9`` ⇒ 10× more resistant).
h. ``armor_diffusion`` — also armor the hillslope diffusivity ``Cd`` by the same factor (default ``False``).
i. ``break_rate`` / ``decay_rate`` — breakdown per unit surface incision ``k_break`` and the slow disequilibrium decay ``k_decay`` (1/yr) away from the fringe.

The optional ``weathering:`` sub-block selects the **solute supply** ``Ψ`` feeding
formation: ``mode: proxy`` (default, climate/temperature stand-in), ``rate`` (an
explicit Maher–Chamberlain chemical-weathering rate driven by the recharge, keys
``C_eq``/``Dw``/``path_length``/``weatherability``), or ``prodsoil`` (reuse the
soil production rate). ``rate`` and ``prodsoil`` fall back to the proxy when soil
is off.

A further optional ``geochem:`` sub-block turns on **conservative solute
geochemistry** (Level B). One or more lumped **tracers** (a "species" = a solute
class such as carbonate, silica or iron) are, each step: **dissolved** on
subaerial land, **transported** down the groundwater flux, **precipitated** where
the water table sits in the capillary fringe (feeding the crust ``duriH``) and
the remainder **exported** to the rivers — a closed mass budget
(``dissolved = precipitated + exported``) with a dissolved flux to the ocean.
When ``geochem:`` is on the transported precipitation **replaces** the Level-A
proxy/rate supply as the crust source — so the duricrust ``form_rate`` is **not
used**; the per-species ``precip_rate`` (below) governs crust growth instead.
Coupled multi-species aqueous equilibrium (speciation, pH, activity) is
deliberately **out of scope** — the lumped multi-tracer model is the
scale-appropriate choice at km / My.

A single tracer ships by default; several tracers give calcrete / silcrete /
ferricrete typing (the dominant one per node is the ``crust_type`` output). Full
declaration with sensible values:

.. code:: yaml

    groundwater:
        # ... hydrology + duricrust keys ...
        geochem:
            conserve: True                # keep/report the closed solute budget
            weatherability_from: source_class   # (optional) lithology label for the tables below
            lithology: [litho, rock_class]      # (optional) standalone per-vertex lithology map
            river_load: True              # (optional) route the export down the rivers
            marine_coupling: True         # (optional) accumulate the coastal delivery
            species:
              - name: carbonate           # tracer label (names the per-species outputs)
                weatherability: 1.0        # relative dissolution rate (scalar, map, or by-class)
                precip_rate: 1.0           # fringe-precipitation RATE (1/yr); scale down for large dt
                solid_volume: 1.0          # crust volume per unit precipitated solute
                c_sat: 1.0                 # saturation threshold (reserved; see note)
                river_decay: 0.0           # in-transit river loss (0 = conservative)
              - name: silica
                weatherability: 0.5
                precip_rate: 0.8
                river_decay: 0.3

**Block-level** ``geochem:`` keys:

.. list-table::
    :header-rows: 1
    :widths: 22 58 20

    * - Key
      - Meaning
      - Default
    * - ``species``
      - list of tracers (one mapping per species; keys below). Empty / absent ⇒ one default tracer.
      - one tracer
    * - ``conserve``
      - keep the per-step mass budget (``dissolved = precipitated + exported``) and report it.
      - ``True``
    * - ``weatherability_from``
      - lithology label for the ``weatherability_by_class`` tables: ``source_class`` (static bedrock) or ``surface_class`` (dynamic top stratigraphic layer). Needs ``provenance:``.
      - *(none)*
    * - ``lithology``
      - ``[file, key]`` per-vertex **integer** lithology map (an alternative label for the by-class tables, independent of provenance).
      - *(none)*
    * - ``river_load``
      - route the exported solute **down the drainage network** to the coast (per species) — the river dissolved load ``riverSolute``.
      - ``False``
    * - ``marine_coupling``
      - accumulate the delivered coastal flux into a per-species marine reservoir + write ``marineSoluteInput``. Needs ``river_load``.
      - ``False``

**Per-species** keys (inside each ``species:`` entry). ``weatherability``,
``precip_rate`` and ``solid_volume`` are **relative, dimensionless** tuning knobs
— what matters is their ratio between species, not the absolute value:

.. list-table::
    :header-rows: 1
    :widths: 20 60 20

    * - Key
      - Meaning
      - Default
    * - ``name``
      - tracer label; names the per-species outputs (e.g. ``solute_carbonate``).
      - ``solute<i>``
    * - ``weatherability``
      - relative **dissolution** rate (scales the Level-A weathering supply for this tracer). A scalar, a per-vertex map ``[file, key]``, or a per-lithology ``weatherability_by_class`` list. ``0`` ⇒ this species does not weather here.
      - ``1.0``
    * - ``weatherability_by_class``
      - list of per-lithology weatherabilities gathered by the ``lithology`` map or the provenance class (see the three forms below).
      - *(none)*
    * - ``precip_rate``
      - **fringe-precipitation rate** ``k_p`` (per year) — how fast the transported solute precipitates into crust in the capillary fringe. Precipitation is **self-limiting** (``·(1−duriH/max_thickness)``): a crust at ``max_thickness`` rejects further solute, which is exported instead. ``k_p`` is a **rate**, so scale it with the timestep — at a My ``dt`` a value near 1 saturates the crust to ``max_thickness`` in one step; use a small value (e.g. ``1e-6``) for gradual growth.
      - ``1.0``
    * - ``solid_volume``
      - crust **volume produced per unit precipitated solute** (a molar-volume-like factor turning precipitated mass into ``duriH``).
      - ``1.0``
    * - ``source_pool``
      - dissolvable **source-rock reservoir** per unit cell area (the weatherable rock mass). Dissolution debits it, so a pool too small for the run's weathering rate **exhausts** and the tracer goes inert (a one-time ``[gw] geochem: source pool exhausting`` warning is printed). Set it **large** (a thick saprolite source) to keep weathering *rate-limited* over long, fast-weathering runs, or small to model a finite, exhaustible weathering front (e.g. a wet-phase iron pulse).
      - ``1.0e6``
    * - ``c_sat``
      - saturation threshold (**reserved** — precipitation is currently a linear fringe sink; a hard ``c > c_sat`` gate is a future nonlinear refinement).
      - ``1.0``
    * - ``river_decay``
      - first-order **in-transit loss** along the river (in-channel precipitation / uptake), per drainage step. ``0`` = conservative; a small value (``0.1``–``0.5``) for a reactive species.
      - ``0.0``

**Spatial (lithology-driven) weatherability.** The ``weatherability`` may vary in
space so *lithology* controls which species each region yields (mafic rock →
Fe / silica, a carbonate platform → carbonate). Three forms, on top of the scalar
default:

.. code:: yaml

        geochem:
            lithology: [litho, rock_class]         # (b) per-vertex integer lithology map
            # weatherability_from: source_class    # (c) OR reuse the provenance regions
            species:
              - {name: carbonate, weatherability_by_class: [1.0, 0.0]}   # per rock class
              - {name: silica,    weatherability: [litho, sil_wab]}      # (a) per-vertex map

* **(a)** ``weatherability: [file, key]`` — a **per-vertex map** for that species
  (loaded like ``infiltration``).
* **(b)** ``weatherability_by_class: [...]`` with a standalone ``lithology:
  [file, key]`` integer map — a **per-(class, species) table** gathered by that
  per-vertex lithology label (independent of provenance).
* **(c)** the same ``weatherability_by_class`` with ``weatherability_from:``
  a provenance class — ``source_class`` (**static** bedrock) or ``surface_class``
  (**dynamic** — the dominant provenance of the top stratigraphic layer,
  re-derived each step, so the weatherability tracks the rock actually exposed as
  erosion exhumes deeper layers or deposition buries the surface). Needs
  ``provenance:``; no extra input.

The label for (b)/(c) is the ``lithology:`` map when given, otherwise the
provenance class.

**River dissolved load** (``river_load: True``). The exported (seepage / baseflow)
solute is routed **down the drainage network** to the shoreline — the dominant
natural pathway of weathering products to the sea — **per species** on the flow
matrix. ``riverSolute`` (m³/yr) grows downstream, is delivered at the coast, and
is trapped in closed continental basins (evaporite behaviour). Two refinements:

* per-species ``river_decay`` — a first-order **in-transit loss** (in-channel
  precipitation / uptake; ``0`` = conservative);
* ``marine_coupling: True`` — accumulate the delivered coastal flux into a
  per-species marine reservoir + write the per-node ``marineSoluteInput``
  (where weathering solute enters the sea).

Extract the per-basin dissolved / species flux at river mouths with
``gospl-catchment`` (see :ref:`running`).

.. note::

    **Outputs.** Water table: ``recharge``, ``wtable``, ``wtdepth``, and
    ``baseflow`` (with ``conserve_baseflow``). Duricrust: ``duricrust``,
    ``induration``, ``Karmor``; with stratigraphy, the per-layer archive
    ``stratDuri`` (degree) plus ``stratCrustType`` / ``stratCrustSource`` (which
    crust and from where), all shown by ``gospl-strata-volume``. Geochemistry:
    ``solute`` (concentration), ``soluteflux`` (dissolved export), ``crust_type``
    (dominant former, several tracers), ``crust_source`` (with ``provenance:``),
    ``riverSolute`` (with ``river_load``) and ``marineSoluteInput`` (with
    ``marine_coupling``). With several tracers every aggregate is **also written
    per species** (``solute_<name>``, ``soluteflux_<name>``, ``crust_<name>``,
    ``riverSolute_<name>``). See the technical
    `groundwater documentation <https://gospl.readthedocs.io/en/latest/tech_guide/groundwater.html>`_
    for the formulation.
