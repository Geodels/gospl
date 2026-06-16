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
        h. ``oFill`` is the minimum elevation (m) below which the priority-flood algorithm is not applied — used to skip deep ocean cells (default: -6000.0).

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
        abrasion, till transport and ice loading. Two flow models are available
        via ``flow_model`` — the full **Shallow-Ice-Approximation (SIA)**
        (``sia``, default): an explicit, mass-conserving non-linear diffusion of
        the ice thickness; or a cheap, stable **diagnostic** (``mfd``): the ELA
        accumulation is routed downhill into an ice discharge from which a Bahr
        thickness and a balance velocity are derived — no dynamics solve. Use
        ``mfd`` when the **morphology of glacial erosion** matters more than the
        ice dynamics themselves (it is fast and robust at any resolution, where
        the SIA becomes stiff for km-thick continental ice). The full algorithm
        is described in the technical guide (:ref:`ice`).

        **Declaration example**:

        .. code:: yaml

            ice:
                # flow_model: sia          # 'sia' (default) | 'mfd' (diagnostic)
                # Either constant glacial parameters
                hterm: 1700.0
                hela: 1850.0
                hice: 2100.0
                # Or using a file to characterise glacial evolution
                # evol: 'data/ice_evol.csv'
                # hinit: ['input/ice0', 'H']   # optional pre-existing ice
                # Diagnostic ('mfd') controls (ignored under flow_model: sia):
                # icedir: 1                 # MFD flow directions for ice routing
                # eheight: 0.25             # Bahr thickness factor
                # fwidth: 1.5               # Bahr width factor
                # melt: 10.                 # ablation amplifier
                sia:
                    Aglen: 1.0e-16
                    slide: 1.0e-3
                    glen: 3.0
                    # cfl: 0.5            # explicit substep accuracy (optional)
                    # max_substeps: 500   # substep cap per goSPL step (optional)
                    # accum_factor: 1.0   # precip->ice accumulation fraction (optional)
                    # accum_max: 2.0      # cap accumulation rate, m ice/yr (optional)
                abrasion:
                    Kg: 1.0e-4
                    l: 1.0
                till:
                    on: True
                    route: False

        ``flow_model`` selects how ice is computed: ``sia`` (default) solves the
        Shallow-Ice-Approximation thickness; ``mfd`` is the diagnostic proxy
        (route the accumulation into an ice discharge, then a Bahr thickness and
        a balance velocity — no dynamics solve). Both then drive the *same*
        abrasion / till / loading machinery below. The diagnostic ``mfd`` adds
        four optional controls: ``icedir`` (number of MFD flow directions for the
        ice routing), ``eheight`` and ``fwidth`` (the Bahr thickness/width
        scaling factors :math:`H = \mathrm{eheight}\cdot\mathrm{fwidth}\cdot Q^{0.3}`),
        and ``melt`` (an ablation amplifier). They are ignored under ``sia``.

        The equilibrium-line / ice-cap geometry controls where ice accumulates
        and melts:

        a. ``hterm`` is the glacier terminus elevation (m) — no ice is kept below it. The effective floor is ``max(hterm, sea level)``: ice never persists below the (possibly time-varying) sea surface, and a ``hterm`` below sea level is raised to it. When omitted, the terminus defaults to the **sea-level position**, so the SIA dynamics and ablation set the actual terminus,
        b. ``hela`` is the equilibrium-line altitude (m) — ablation below, accumulation above,
        c. ``hice`` is the ice-cap altitude (m) — full precipitation is captured as ice above it.

        The ``sia`` sub-block sets the flow physics (all optional, with the
        defaults shown above):

        d. ``Aglen`` is the Glen's-law rate factor (ice softness) controlling internal deformation,
        e. ``slide`` is the basal-sliding coefficient,
        f. ``glen`` is the Glen's-law exponent :math:`n` (usually 3),
        g. ``cfl`` (default ``0.5``) sets the explicit-substep size as a fraction of the time for a cell to shed (or gain) its ice — an accuracy/cost knob, not a stability limit (the flux limiter keeps the solve stable and positive at any size); ``max_substeps`` (default ``500``) caps the substeps per goSPL step. Thick (km-scale) ice is genuinely stiff: it needs more substeps for accuracy, so lower ``cfl`` (e.g. ``0.1``) for thick ice, accepting a higher cost.
        h. ``accum_factor`` (default ``1.0``) and ``accum_max`` (default unset) control the **accumulation** part of the surface mass balance only (ablation is untouched). Full precipitation is rarely all snow/ice, so ``accum_factor`` is a precipitation→ice conversion fraction and ``accum_max`` caps the accumulation rate (m ice/yr) at a realistic ceiling (real ice sheets accumulate ~0.1–2 m/yr). **Set these for high-precipitation runs** — converting several m/yr of rainfall directly to ice produces unphysically thick, expensive-to-solve ice.

        The ``abrasion`` sub-block enables velocity-based glacial erosion
        :math:`E_g = K_g\,|u_b|^{l}` (off by default, ``Kg: 0``):

        i. ``Kg`` is the abrasion coefficient (default ``0.0`` — set it to enable glacial erosion),
        j. ``l`` is the basal-sliding-velocity exponent (default ``1.0``).

        The ``till`` sub-block controls glacial sediment (default off):

        k. ``on`` — when ``True``, abraded rock is carried as **till** and
           deposited as a moraine where the ice melts out (the ablation zone),
           conserving the abraded volume. With stratigraphy on, the till is
           layered into the stratigraphic record and split into the coarse/fine
           lithology fractions when dual lithology is enabled.
        l. ``route`` (default ``False``) — controls how the till is distributed.
           ``False`` spreads it across the whole ablation zone weighted by the
           meltwater rate (appropriate when a cell aggregates a glacier, i.e.
           continental/global resolution). ``True`` instead **routes the till
           down the ice-surface flow network** and melts it out toward each
           catchment's terminus, building moraine at the actual ice margins — for
           high-resolution (sub-km) regional runs where individual glacier
           catchments and termini are resolved. Both conserve mass.

        The glacier geometry can instead be read from a file:

        m. ``evol`` is the glacier characteristics over time (`csv` file). When
           used, ``hterm``, ``hela`` and ``hice`` are not required because they
           are defined in this file.
        n. ``hinit`` (optional) is a **pre-existing ice thickness** (m) — a
           uniform scalar or a per-vertex ``[file, key]`` map — used to seed the
           ice at the start of the run; the SIA solve then evolves it. Without
           it the ice grows in from zero. The equilibrium-line geometry
           (``hela``/``hice``/``hterm``) is still required, as it drives the
           evolution. On a restart the evolved ice thickness is restored
           instead, so ``hinit`` only applies to a fresh start.

        When flexural isostasy is enabled, the SIA ice thickness is automatically
        used as the ice load contribution to the isostatic computation — no extra
        parameters are needed.

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
                sia: {Aglen: 1.0e-16, slide: 1.0e-3, glen: 3.0}

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
                sia: {Aglen: 1.0e-16, slide: 1.0e-3, glen: 3.0}
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

        **Deriving the ELA from paleo-climate temperature.** The helper
        ``scripts/ela_from_temperature.py`` turns a per-vertex temperature map
        into the ``hela``/``hice`` ``.npz`` maps above by lapse-rate inversion
        (``ELA = z + (T - T_ELA)/Gamma``). It derives the ELA *position* only;
        the ablation magnitude stays precipitation-scaled (it is not a
        degree-day melt model).

        From a **terminal** (run once per climate snapshot; ``--start`` prints a
        ready-to-paste ``glaciers`` entry):

        .. code:: bash

            python scripts/ela_from_temperature.py \
                --temperature climate/t2m_21ka.npz --t-key t2m \
                --reference surface --elevation input/mesh.npz --z-key z \
                --lapse 0.0065 --t-ela -2.0 --band 400 \
                --out input/ela_21ka.npz --start -21000

        ``--help`` lists every option; loop over snapshots with a shell ``for``.

        From a **Jupyter notebook**, either shell out with the ``!`` magic:

        .. code:: python

            !python scripts/ela_from_temperature.py \
                --temperature climate/t2m_21ka.npz --t-key t2m \
                --reference surface --elevation input/mesh.npz --z-key z \
                --lapse 0.0065 --t-ela -2.0 --band 400 --out input/ela_21ka.npz

        or import the conversion and build the ``glaciers`` list directly:

        .. code:: python

            import sys, numpy as np
            sys.path.append("scripts")
            from ela_from_temperature import derive_ela

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
                soilK: 4.e-6
                maxProd: 50.e-6
                depthProd: 0.5
                roughnessL: 1.0
                decayDepth: 0.7
                bedrockConv: 0.0001
                uniform: 0.5
                map: ['test_mesh8/hsoil', 'soil']

        a. ``soilK`` is the erodibility coefficient for soil,
        b. ``maxProd`` is the soil production maximum rate (m/yr),
        c. ``depthProd`` is the soil production decay depth (m),
        d. ``roughnessL`` is the roughness length scale,
        e. ``decayDepth`` is the soil transport decay depth for non-linear diffusion where the coefficient of diffusion is set to the values of ``hillslopeKa`` and ``hillslopeKm``,
        f. ``bedrockConv`` is the soil to bedrock conversion fraction, bedrock begins where soil production is a very small fraction of the maximum soil production (optional). 

        Then the user can specify the initial soil thickness if any by setting **either**:

        g. ``uniform`` a uniform soil thickness on the entire surface (m),

        **or**:

        h. ``map`` a soil thickness map. 

        .. important::

            When defining a soil thickness grid, one needs to use the **npz** format and needs to specify the key corresponding to the soil thickness value in the file. In the above example this key is ``'soil'``. The soil grid needs to define values for all vertices in the mesh in metres.

        The soil-aware non-linear SPL is solved with a PETSc ``SNES``. Its
        behaviour can be tuned (all optional) with:

        i. ``maxIter`` is the maximum number of non-linear iterations (default ``500``),
        j. ``rtol`` / ``atol`` are the relative / absolute convergence tolerances (default ``1.e-6``),
        k. ``pcType`` is the preconditioner for the ``ngmres`` Krylov solve (default ``'hypre'`` BoomerAMG; ``'gamg'``, ``'bjacobi'`` or ``'asm'`` can help on heavily-decomposed / ocean-dominated partitions),
        l. ``solver`` selects the primary non-linear solver: ``'qn'`` (default, limited-memory quasi-Newton / L-BFGS) or ``'ngmres'`` (accelerator + multigrid preconditioner).

        .. tip::

            The soil solve is usually the dominant cost of a soil-enabled run. The default ``solver: 'qn'`` was chosen because on a global model it cut the soil-solve wall time roughly in half (or more) versus ``'ngmres'`` at the **same** tolerance and solution — L-BFGS reaches a comparable iteration count but each iteration is far cheaper than an ``ngmres``/multigrid sweep. If ``'qn'`` struggles on a particular configuration, set ``solver: 'ngmres'``. Prefer switching ``solver`` over relaxing ``rtol``: a loose ``rtol`` can leave the elevation field under-resolved and **destabilise the downstream sediment-routing solver**, which is both slower and less accurate.

        .. note::

            If the chosen primary solver stalls on a stiff soil-production residual, goSPL automatically retries that timestep with the *complementary* solver (quasi-Newton ⇄ ``ngmres`` multigrid accelerator) at a relaxed tolerance before continuing.


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
        e. ``pitInletBias`` — per-fraction lake-deposition bias (coarse builds
           inlet deltas, fine settles in the depocenter).

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
                cu_weight: [1.0, 0.0, 0.3]            # optional copper fertility

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
        c. ``uniform`` — a single source class everywhere,
        d. ``cu_weight`` (optional) — copper fertility per class, for a
           Cu-sourced fraction diagnostic.

        The per-layer composition is written to the stratal output (``stratP``).

        .. note::

            See the technical `documentation <https://gospl.readthedocs.io/en/latest/tech_guide/provenance.html>`_ (provenance section) for the algorithm, the standalone post-processing tool, and the copper-prospectivity scope.
