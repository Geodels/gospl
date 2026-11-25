.. list-table::
   :header-rows: 1
   :widths: 28 13 9 9 9 32

   * - Option
     - added to
     - 1 rank
     - 4 ranks
     - added time
     - notes
   * - base
     -
     - 100 (1.67 s)
     - 100 (0.53 s)
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
     - 78
     - 82
     - 0
     - changes the relief: the marine step does less work here
   * - transport-limited SPL
     - base
     - 137
     - 137
     - 33
     - 
   * - non-linear hillslope
     - base
     - 109
     - 96
     - 0
     - 
   * - soil
     - base
     - 127
     - 120
     - 25
     - soil-aware SPL and soil diffusion (SNES, TS)
   * - **Climate, ice and the lithosphere**
     -
     -
     -
     -
     -
   * - flexure
     - base
     - 110
     - 100
     - 1
     - operator factorised once
   * - orography
     - base
     - 84
     - 87
     - 2
     - lowers the rain here: less flow and marine work
   * - ice
     - base
     - 107
     - 104
     - 11
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
     - 110
     - 5
     - a second (mud) routing in the cascades
   * - provenance (3 classes)
     - dual lithology
     - 143
     - 144
     - 23
     - one routed column per class
   * - **Groundwater and geochemistry**
     -
     -
     -
     -
     -
   * - groundwater
     - base
     - 124
     - 140
     - 37
     - AMG head solve
   * - + duricrust
     - groundwater
     - 124
     - 135
     - 0
     - 
   * - + geochemistry
     - duricrust
     - 140
     - 154
     - 26
     - one direct (MUMPS) solve per species
   * - **Marine processes**
     -
     -
     -
     -
     -
   * - gravity currents
     - dual lithology
     - 137
     - 131
     - 21
     - its own network and solve, inside the marine step
   * - pelagic rain
     - dual lithology
     - 159
     - 153
     - 40
     - clay decay solve, and carb_mud added to the marine routing
   * - carbonates
     - dual lithology
     - 187
     - 173
     - 56
     - three more lithologies in the marine routing
   * - + exposure
     - carbonates
     - 177
     - 175
     - 3
     - one screened-Poisson solve
   * - + turbidity
     - carbonates
     - 185
     - 186
     - 9
     - 
   * - + waves
     - carbonates
     - 198
     - 177
     - 7
     - one upwind solve per wave direction
   * - + coastal erosion
     - carbonates
     - 203
     - 188
     - 0
     - 
   * - + karst
     - carbonates
     - 199
     - 187
     - 1
     - 
