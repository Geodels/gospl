.. _carb:

==============================
Carbonate production
==============================

Benthic carbonate production (``carbonate:`` block) follows Bosscher and
Schlager (1992), with one or more carbonate factories as in CarboKitten
(Hidding et al., 2026). Below sea level, factory :math:`f` produces at

.. math::

    P_f(w) = g_f \tanh\left(\frac{I_0\, e^{-k_f w}}{I_{k,f}}\right)
             M_f \prod_j R_{fj}

where :math:`w` is the water depth, :math:`g_f` the maximum growth rate,
:math:`I_0` the surface light, :math:`k_f` the extinction coefficient,
:math:`I_{k,f}` the saturating light, :math:`M_f` an optional habitat map and
:math:`R_{fj}` piecewise-linear response factors on the depth, the seafloor
slope and the clastic deposition rate of the previous step.

**Integration over a step.** Each cell grows with the depth it creates,

.. math::

    \frac{dh}{dt} = \sum_f P_f(\eta - h), \qquad h \le \eta - w_{min},

with :math:`\eta` the sea level. The step is split into second-order (Heun)
sub-steps of about 0.25 m of growth at the largest possible rate, so the result
does not depend on the model time step: at 1 and 10 kyr steps it matches the
exact solution of this ODE to :math:`10^{-4}` (``benchmarks/test_carbonate.py``).

**Keep-up and give-up.** Under a relative sea-level rise :math:`\sigma`, the
depth evolves as :math:`dw/dt = \sigma - P(w)`. Because production grows as the
water shallows, the depth :math:`w^*` where :math:`P(w^*) = \sigma`,

.. math::

    w^* = \frac{1}{k} \ln\left(\frac{I_0}{I_k\, \mathrm{atanh}(\sigma/g)}\right),

is unstable: a platform shallower than :math:`w^*` reaches sea level and keeps
up with it (growth limited to :math:`\sigma` by the accommodation), a deeper
one drowns. The benchmark checks both sides of :math:`w^*` and the drowning rate.

**Framework and loose sediment.** A fraction :math:`\phi_f` of each factory's
production is built in place (lithology ``reef``); the rest is split into
carbonate sand (coarse-like) and mud (fine-like) and added to the marine
routing of the same step, so carbonate and river sediment share the marine
cascade, diffusion and stratigraphy. Each carbonate lithology keeps its own
compaction curve and erodibility.
