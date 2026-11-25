.. list-table::
   :header-rows: 1
   :widths: 28 13 9 9 9 32

   * - Option
     - added to
     - 1 rank
     - 4 ranks
     - own phases
     - notes
   * - base
     -
     - 100 (2.47 s)
     - 100 (0.71 s)
     -
     - mostly the marine deposition and the flow accumulation
   * - **Erosion and hillslope laws**
     -
     -
     -
     -
     -
   * - non-linear SPL
     - base
     - 53
     - 58
     - 0
     - changes the relief: the marine step does less work here
   * - transport-limited SPL
     - base
     - 112
     - 124
     - 25
     - 
   * - non-linear hillslope
     - base
     - 98
     - 111
     - 0
     - 
   * - soil
     - base
     - 200
     - 173
     - 78
     - soil-aware SPL and soil diffusion (SNES, TS)
   * - **Climate, ice and the lithosphere**
     -
     -
     -
     -
     -
   * - flexure
     - base
     - 93
     - 100
     - 1
     - operator factorised once
   * - orography
     - base
     - 49
     - 53
     - 2
     - lowers the rain here: less flow and marine work
   * - ice
     - base
     - 100
     - 105
     - 8
     - 
   * - **Stratigraphy**
     -
     -
     -
     -
     -
   * - stratigraphy
     - base
     - 102
     - 103
     - 2
     - memory grows with the layers (below)
   * - dual lithology
     - stratigraphy
     - 111
     - 115
     - 0
     - a second (mud) routing in the cascades
   * - provenance (3 classes)
     - dual
     - 150
     - 159
     - 0
     - one routed column per class
   * - **Groundwater and geochemistry**
     -
     -
     -
     -
     -
   * - groundwater
     - base
     - 114
     - 129
     - 29
     - AMG head solve
   * - + duricrust
     - groundwater
     - 112
     - 130
     - 28
     - 
   * - + geochemistry
     - duricrust
     - 112
     - 170
     - 95
     - one direct (MUMPS) solve per species
   * - **Marine processes**
     -
     -
     -
     -
     -
   * - gravity currents
     - dual
     - 127
     - 127
     - 0
     - its own network and solve, inside the marine step
   * - pelagic rain
     - dual
     - 151
     - 158
     - 9
     - clay decay solve, and carb_mud added to the marine routing
   * - carbonates
     - dual
     - 195
     - 182
     - 2
     - three more lithologies in the marine routing
   * - + exposure
     - carbonates
     - 202
     - 192
     - 4
     - one screened-Poisson solve
   * - + turbidity
     - carbonates
     - 191
     - 188
     - 11
     - 
   * - + waves
     - carbonates
     - 186
     - 187
     - 8
     - one upwind solve per wave direction
   * - + coastal erosion
     - carbonates
     - 189
     - 188
     - 2
     - 
   * - + karst
     - carbonates
     - 185
     - 182
     - 1
     - 
