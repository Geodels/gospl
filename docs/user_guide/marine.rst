.. _marine_systems:

==============================
Marine sedimentary systems
==============================

Beyond the marine deposition of river sediment (the clinoform cascade and
the marine diffusion, see :ref:`the surface processes parameters <surfproc>`),
two opt-in modules represent marine sedimentary systems:

.. grid:: 1 1 2 2
    :gutter: 2

    .. grid-item-card:: Marine gravity currents
        :link: gravity
        :link-type: ref

        Part of the river load travels down the slope as hyperpycnal flows or
        turbidity currents and builds coarse belts at the foot of the slope
        and fine lobes on the basin floor; slope failure and erosion.

    .. grid-item-card:: Carbonate platforms and reefs
        :link: carbonates
        :link-type: ref

        Carbonate factories produce reef framework and loose carbonate sand
        and mud on the seafloor; pelagic carbonate and hemipelagic clay
        settle in the ocean basins.

Both work on regional and global meshes, with or without dual lithology
(the carbonates need it), and together: carbonate on the slope can fail and
be carried by the gravity currents.

**Order within a time step.** River erosion and continental deposition;
carbonate production (framework built in place, loose carbonate kept for the
marine step); the marine step (gravity supply taken from the river load,
clinoform cascade and marine diffusion of the river and loose carbonate
sediment, gravity currents); pelagic rain (unless resuspended into the marine
diffusion); hillslope diffusion.

.. toctree::
    :maxdepth: 2
    :hidden:

    gravity
    carbonates
