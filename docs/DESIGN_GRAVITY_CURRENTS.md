# Design: marine gravity currents (opt-in `gravity:` block)

Status: phases G1 and G2 done, bypass and accommodation cap added after the first
regional test (2026-11-23); G3 to do.

## 1. Purpose and scope

goSPL distributes the river load at sea with the clinoform cascade (`seaplex._distOcean`,
capacity set by water depth and distance from the shore) followed by the marine
diffusion. Neither depends on discharge, and the dual-lithology marine composition puts
fine sediment in the deepest water (`_marineFineFraction`). So the model cannot carry
sediment past the foresets to the base of the slope, which is what hyperpycnal flows,
turbidity currents and slope failures do. That is the origin of coarse bottomsets and
basin-floor fans.

This feature adds a **second, routed pathway** for the part of the river load that
travels as a gravity current. It is opt-in and, when the `gravity:` block is absent,
every path is gated out and results are byte-identical.

Out of scope: a depth-averaged turbidity-current solver (Parker et al. type). It is
explicit and limited to time steps of seconds to minutes within single events, which
does not fit goSPL's kyr steps or global meshes. It can be used offline to calibrate
the run-out lengths below.

## 2. Governing equations

### 2.1 Routing with deposition

The gravity flux `q` (m3 over the step) travels down the marine flow network already
built for the clinoform cascade (`dMat1 = W^T`, weights from the smoothed, eps-filled
bathymetry). At each node a fraction `Lambda_i` of the flux passing through it deposits
and the rest moves on:

    q_i = s_i + sum_j W_ji (1 - Lambda_j) q_j        (I - W^T diag(1 - Lambda)) q = s
    deposit_i = Lambda_i q_i

`s` is the gravity supply (section 2.3). The operator is an M-matrix (column sums of
`W^T diag(1-Lambda)` are at most 1), so the solve is well posed; it is solved with the
shared routing KSP (`flowplex._solve_KSP(..., route=True)`), like the river-solute
routing with in-transit loss. No new cached solver.

Boundary nodes:
- terminal sinks (`is_sink_local`, zero column in `W^T`): `Lambda = 1`, the flux
  deposits there (as the cascade's force-deposit);
- open outlets of a flat model (`outletIDs`): nothing deposits (`Lambda = 0`) and the
  flux leaves the domain (export).

Budget (machine precision, checked every step): supplied = deposited + exported.

### 2.2 Deposition law: run-out length

    Lambda_i = 1 - exp(-l_i / L)

with `l_i = sqrt(Omega_i)` the cell size and `L` a run-out length (m). On a uniform
1-D slope the deposit then decays as `exp(-x/L)` whatever the mesh spacing (as long as
`l << L`): the same `L` gives the same deposit at regional and global resolution. This
is the default and the only law of phase G1. Slope-dependent deposition (bypass on steep
slopes) is phase G3, because slopes depend on the mesh (global meshes smear the shelf
break).

With dual lithology the coarse and fine fractions each have their own run-out length
(`runout`, `runout_fine`; coarse shorter), so a single flow builds a coarse proximal
lobe and a finer distal one. Each fraction is one solve with the same `W^T`.

### 2.2b Bypass and accommodation (added after the first regional test)

With the run-out law alone, flows started at the river mouths deposited a fixed share per
cell from the coast, so most of the load dropped on the inner shelf, starved the
clinoform front (the shelf break retreated ~15 km on the `dual_lithology` margin) and
built mounds up to sea level; and the deposit could only decrease away from the source,
so no lobe could form at the base of the slope. Two changes, both setting `Lambda = 0`
(so the solve stays linear and the budget exact):

- **depth bypass** (`bypass_depth` = `d_b`, m below sea level, default 50; 0 = off):
  `Lambda` is multiplied by `clip((d - d_b)/d_b, 0, 1)` (water depth `d`): no
  deposition on the shelf above the wave base, full deposition from `2 d_b`.
  Resolution-independent.
- **slope bypass** (`bypass_slope` = `S_b`, off by default): `Lambda` is multiplied by
  `1/(1 + (S/S_b)^4)` (gradient magnitude `cellslope`, synced owner -> ghost; one half
  at `S_b`), so flows cross the foresets and deposit at their toe. Mesh-dependent: for
  regional runs.

The bypasses were first hard switches (`Lambda = 0` beyond the threshold). On a rough
foreset the gradient is noisy from cell to cell, so neighbouring cells alternated
between deposit and bypass and the deposits flattened their own cells: a speckled,
lumpy deposit (seafloor roughness, RMS departure from a 2 km smoothed surface over the
slope and basin, 6.4 m against 4.1 m without gravity currents; gravity deposit 9.6 m).
The smooth tapers give 3.6 and 3.7 m with the same coarse toe deposit.

And an **accommodation cap**: a gravity deposit cannot raise a cell above
`sealevel - bypass_depth` (room = that height minus the cascade deposit, times the cell
area; the fine fraction gets what the coarse one leaves). It is applied by iterating the
total solve, setting `Lambda = room / q` wherever a node overfills, until none does
(the provenance classes are then solved with the final `Lambda`). Terminal sinks are
exempt (mass must stay in the domain).

Regional test (`dual_lithology`, 250 kyr, run-out 5/20 km): with `fraction: 0.3` and
both bypasses, a coarse body (coarse fraction ~0.6) forms at the foreset toe and a 5 to
40 m deposit covers the basin floor; with the flood trigger (`sigma: 2.0`) the wedge is
unchanged and a coarse bottomset (~0.58) forms at the toe. Depth bypass alone gives only
a thin fine apron. With the default flood parameters (`sigma: 1.2`) the diverted share is
0.005 to 0.2 % on this margin.

### 2.3 Trigger: which part of the river load becomes a gravity current

The supply is a fraction `f_h` of the river load arriving at each marine node
(`sedFlux`, and `fineFlux` with dual lithology): `s = f_h * sedFlux`. Two options:

1. **Fixed fraction** (`fraction: f`): `f_h = f` everywhere. For regional calibration.
2. **Flood fraction** (`flood:` sub-block, default when `fraction` is absent). goSPL's
   fluxes are long-term means, and the mean concentration `C = rho_s Qs/Q`
   (`rho_s E/P` for a catchment eroding at `E` under rain `P`) is far below the
   hyperpycnal threshold of about 40 kg/m3 (Mulder and Syvitski, 1995) almost
   everywhere. Gravity flows come from floods. So `f_h` is the fraction of the long-term
   load carried in flows above a threshold concentration `C_c`, for daily discharges
   `X = Q/Qbar` log-normal with log-standard deviation `sigma` (mean 1,
   `mu = -sigma^2/2`) and a rating curve `Qs ~ Q^b`:

        C(X) = Cbar X^(b-1) / E[X^b],   E[X^b] = exp(b mu + b^2 sigma^2 / 2)
        X0   = (C_c E[X^b] / Cbar)^(1/(b-1))
        f_h  = E[X^b 1{X > X0}] / E[X^b] = Phi((mu + b sigma^2 - ln X0) / sigma)

   with `Cbar = rho_s * Qs / Q` at the node (`Qs = vSed`, `Q = FA`). Analytic, one
   value per node, no solve. Parameters `c_crit` (kg/m3), `rating_exp` (b > 1),
   `sigma`, `rho_s`. `f_h` is set to 0 where `Q` or `Qs` is zero.

Limitation (document in the user guide): on coarse global meshes, small steep coastal
catchments, the most hyperpycnal-prone rivers, are merged into larger outlets with lower
variability, so gravity flows from them are underestimated. A catchment-area dependence
of `sigma` is a possible phase-G3 refinement.

## 3. Position in the step and coupling

In `seaplex.seaChange`, only when the block is present:

1. split the arriving load: `sedFlux -> (1 - f_h) sedFlux` for the clinoform cascade and
   `f_h sedFlux` for the gravity pathway (same for `fineFlux` with dual lithology);
2. the cascade, `_depMarineSystem` (global, `Gmar > 0`) and `_diffuseOcean` run as
   before on the reduced load;
3. the gravity routing runs on the same `dMat1` (before it is destroyed); its deposit
   thickness is added to the marine deposit `self.tmp` after the diffusion (gravity
   deposits are not diffused in the step they form);
4. elevation, `cumED`, `EbLocal`, soil and `deposeStrat` then see the combined deposit.

Dual lithology: `_marineFineFraction` computes the fine fraction of the cascade deposit
only (optional `mdep` argument, unchanged call when off), and the node fine fraction is
the volume-weighted blend of that deposit and the gravity coarse/fine deposits, so the
fine volume is conserved exactly.

Provenance (phase G2, done): each class supply is the arriving composition (`provFrac`)
times the gravity supply, shared between the coarse and fine parts (a passive label has
no grain size). The classes are solved with the same operator as the total, so their
deposits sum to it to the KSP tolerance, and `_gravProvFraction` blends the cascade
composition with the per-class gravity deposits, volume-weighted (as for the fine
fraction).

## 4. Inputs (YAML)

    gravity:
      runout: 5.0e4        # m, run-out length (coarse, or single-fraction)
      runout_fine: 2.0e5   # m, fine fraction (dual lithology; default 4 x runout)
      fraction: 0.2        # optional fixed fraction; if absent the flood trigger is used
      flood:
        c_crit: 40.        # kg/m3 hyperpycnal threshold
        rating_exp: 1.8    # b in Qs ~ Q^b (> 1)
        sigma: 1.2         # log-std of daily discharge
        rho_s: 2650.       # kg/m3

## 5. Outputs

- `gravDep`: cumulative gravity-current deposit thickness (m), per vertex; restored on
  restart (zero when absent from older outputs).
- Verbose line per step: supplied, deposited and exported volumes (km3) and the
  closure error.
- `--summary` record: a `gravity` entry (supplied, deposited, exported, m3) per step;
  `gospl-inspect --summary` reports a `gravity_budget` issue when it does not close.

## 6. Tests (phase G1)

- parse + allocation, and off-path byte-identity of `seaChange` outputs;
- budget: supplied = deposited + exported to round-off, on the marine fixture;
- analytical: 1-D uniform channel of cells, `Lambda` constant -> deposit
  `s (1-Lambda)^k Lambda` at the k-th cell (exact discrete solution), and the run-out
  length independent of the mesh spacing;
- flood fraction: `f_h` against a Monte-Carlo estimate of the same log-normal model;
- dual lithology: fine volume conserved; the coarse deposit nearer the source than the
  fine one;
- np=1 vs np=3 invariance (`mpi`).

## 6b. Why the deposits are still sheets with dots (motivation for G3b)

After the bypass, cap and tapers, the `fraction: 0.3` regional run still shows isolated
dots in the gravity deposit (48 cells more than 50 % above their neighbours, up to
80 m on one cell in one step at the foreset toe) and no channel. An instrumented run
(every gravity step from 190 kyr) rules out terminal sinks, the cap (100-200 m of room
left), local supply and `Lambda` (0.02-0.1, like the neighbours). The cause is the flux:
the marine network funnels the gravity flux into single-cell paths where the foreset is
gullied in front of the river lobes, with 5-10 times the neighbours' flux, and the
deposit is a fixed share of the flux. Two things are missing:

1. **Concentrated flow should bypass or erode, not deposit.** A confined, faster current
   carries its load further and cuts its bed. Without erosion there is no positive
   feedback (a deeper path captures more flow and is cut deeper), so no channel forms.
2. **The network spreads the flow by design.** `_matOcean` computes directions on a
   SMOOTHED bathymetry with up to 12 receivers, made to spread river sediment at sea.
   A channel needs flow concentrated on the actual seafloor; a lobe needs it to spread
   where the channel ends.

## 6c. Phase G3b: gravity network, slope erosion, channels and lobes

### Gravity network

The gravity currents get their own routing matrix `W_g`, built each step like
`_matOcean` but on the eps-filled ACTUAL bathymetry (no smoothing), with `n_g`
receivers (`gravity: receivers`, default 6) and a flow-partition exponent `p_g`
(`gravity: flow_exp`, default 2; weights proportional to `S^p_g`, see `mfdrcvrs`). On
the slope the steepest path dominates (concentration into channels); on the gentle basin
floor the gradients are similar in every direction and the flow spreads (lobes). The
fill reuses the serial `epsfill` path of `_matOcean` (one more gather on rank 0 per step;
acceptable regionally, to be measured globally, where `network: smoothed` keeps the
current matrix).

### Slope erosion, linear in the flux

A current picks up sediment where it is steep and strong. With `l_i` the cell size, an
erosion length `L_e` (`gravity: erosion_length`, m; the shorter, the stronger) and a
slope dependence `f(S) = (S/S_e)^4 / (1 + (S/S_e)^4)` (`gravity: erosion_slope` = `S_e`;
the mirror of the deposition taper: 0 on flat ground, 1/2 at `S_e`, ~1 above), the
volume picked up at node `i` is

    e_i = eps_i q_i,   eps_i = (l_i / L_e) f(S_i)

proportional to the flux through the node, so the flux grows downstream on steep ground
and is self-reinforcing where it concentrates. The balance becomes

    q_i = s_i + sum_j W_ji (1 - Lambda_j + eps_j) q_j
    net deposit_i = (Lambda_i - eps_i) q_i          (negative = erosion)

still linear and exactly conservative (supplied + eroded = deposited + exported). The
operator is no longer diagonally dominant where `eps > Lambda`, but the eps-filled
network is a DAG (triangular in flow order), so the system stays uniquely solvable; the
shared routing KSP and its exact-factor rescue handle it.

Limits (both iterated like the cap: lower `eps_i` where exceeded and re-solve):
- **erodible thickness**: at most `max_erosion` (m per step, default 10) per cell, which
  also guards against the runaway seen with the explicit glacial abrasion (erosion on
  slopes steepens them);
- **no erosion above the bypass depth** (the shelf is the river and wave domain).

### Composition (dual lithology, provenance)

Erosion removes material from the stratigraphy through `erodeStrat` (thickness in
`self.tmp`), which gives the eroded coarse/fine and class split; the river step's
`thCoarse`/`thFine`/`provEro` are saved and restored around the call (the
`hillslopeBeforeRouting` pattern). Because the fractions have their own run-out lengths,
the coupled total is solved first (Picard, 2-3 passes: route each fraction with the
eroded volume of the previous pass added as a source in the eroded composition, update
`e = eps q_total`), and the final pass routes every fraction and class with the final
erosion as an explicit source. Conservation is exact whatever the number of passes,
because the source added to the routing is exactly what `erodeStrat` removed.

### As implemented (2026-11-23)

Three points differ from, or add to, the proposal above.

- **mfdrcvrs with fewer than 12 receivers.** The kernel pads the unused slots with -1
  and keeps the first `n_g` weights without renormalising them: a node with 7 downslope
  neighbours kept 6 weights of 1/7 and 1/7 of its flux vanished (2 to 90 % of the
  supply lost on the test mesh). `_gravNetwork` points the padding at the node itself
  with zero weight and renormalises every row.
- **Receivers in numpy, not `mfdrcvrs` (after the first regional run).** `mfdrcvrs`
  gives every selected receiver the SAME weight `1/k` (the slope-weighted line is
  commented out), so `flow_exp` had no effect and the flux spread evenly over up to 6
  neighbours; and its branch for nodes with more downslope neighbours than receivers
  left ~120 nodes with lower neighbours flagged as sinks. Sinks keep everything (no
  cap), so they grew into isolated mounds up to 420 m in one step. `_gravReceivers`
  now keeps the `receivers` steepest strictly-downhill neighbours with weights
  `S^flow_exp`, normalised, with the same input-id tie break.
- **Pits below the fill cut-off.** `epsfill` starts at `sealevel + oFill`, so a closed
  deep below it keeps its raw elevation, and `mfdrcvrs` hands its flux to all its
  (higher) neighbours: an uphill cycle and a singular solve. Only receivers strictly
  lower than the node are kept; a node with none is a terminal sink (a closed deep keeps
  what reaches it).
- **Erode or deposit, never both, and no Picard loop.** At each node the larger of
  `eps` and `Lambda` (coarse) wins: an eroding node has `Lambda = 0`, a depositing node
  `eps = 0`. The erosion is computed once, on the total gravity supply
  (`_gravErosion`: `eps q`, capped at `max_erosion` by iteration), removed from the
  stratigraphy with `erodeStrat` (river arrays saved and restored), and its coarse,
  fine and class parts are added to the routed supply as an explicit source at the
  eroding nodes, which do not deposit (`nodep`). The budget is exact (supplied + eroded
  = deposited + exported) because the source is exactly what `erodeStrat` removed; the
  Picard passes of the proposal were not needed. The stratigraphy then receives the
  deposits only (`deposeStrat`), while elevation, `cumED`, `EDrate` and `gravDep` get
  the net change; giving `deposeStrat` the net change made the eroded fine count twice
  (5 % fine excess in the test).

### Channel belts without resolved channels (2026-11-23, regional tests)

Slope erosion (linear, on a smoothed slope, or non-linear in the flux) only roughened
the toe of the slope on the 615 m regional mesh with 5 kyr steps: channels are not
resolved there. What worked is to represent channelised transport by its deposits:

- **Separate networks per fraction** (`receivers`/`flow_exp` for coarse,
  `receivers_fine`/`flow_exp_fine` for fine): coarse belts along the main pathways
  (2 receivers, exponent 4), fine lobes beyond (6, exponent 1).
- **Smoothed routing surface and slope** (`network_smooth`, `slope_smooth`, passes of
  area-weighted neighbour averaging, owner -> ghost sync per pass): removed the
  hummocky toe.
- **Lateral spreading of each step's deposit** (`deposit_spread`, a LENGTH in m; passes
  `ceil(2 (L/dx)^2)`): a conservative diffusion of thickness with symmetric weights
  `w_ij = 0.5 min(a_i, a_j)/max(n_i, n_j)`. A first donor-share operator (`v_j` split in
  proportion to the recipients' areas) conserved volume but not a uniform thickness,
  and grew chains of 2-5 m dots at odd-sized cells along the mesh seams (correlation
  0.55 with the cell-area anomaly). Measured: dots 767 -> 65, cell-scale roughness
  1.8 -> 0.5 m, spreading cost 37 ms/step.
- **`coarse_factor`**: with a fixed 30 % fraction, 43 km3 of coarse per 100 kyr went
  into the gravity currents and ~100 % of it deposited at the foot of the slope.
  Turbidity currents carry mostly the suspended load; `coarse_factor: 0.3` brings the
  lobes to a mixed composition.
- Tried and dropped: pit filling by the currents (its iteration never converged, 21
  solves per route, 64 % of the run time, no effect on the dots) and receiver weights
  times the Voronoi face length (no effect).

Cost (regional, 52,805 nodes, 4 ranks): 0.96 s/step without gravity, 1.25 s/step with
(network 41 ms, routing 196 ms, spreading 37 ms).

### Expected behaviour and tests

- inclined plane with a small gully: the flux concentrates, the gully deepens over steps
  (channel), and a lobe forms where it reaches the flat floor;
- budget: supplied + eroded = deposited + exported, every step;
- dual and provenance conservation (`_fineEroded`/`_fineDeposited`, `stratP`);
- erosion only where `S > ~S_e` and deeper than the bypass depth; `max_erosion` respected;
- partition invariance (`ab_partition --steps 1`);
- the regional margin: an incised slope channel and a terminal lobe, no dots.

Effort: about 6-9 working days (network 1, erosion and limits 2, composition 2-3,
tests and the regional example 2-3). Main risks: stability of slope erosion over kyr
steps (hence `max_erosion` and a stability note, as for the glacial abrasion), and the
extra serial fill on global meshes.

## 7. Phases

- **G1 (this branch):** split + routed solve + run-out length + fixed/flood trigger,
  single fraction and dual lithology, `gravDep` output, tests, docs.
- **G2 (done):** provenance classes in lockstep; budget in the `--summary` JSON.
- **G3b (coded, see 6c):** gravity network on the actual bathymetry, linear slope
  erosion with limits, composition of the eroded material. Next: the regional
  channel/lobe test.
- **G3:** slope-failure trigger (remobilise sediment above a critical slope), optional
  slope-dependent deposition, catchment-area-dependent `sigma`; example in
  goSPL-examples (delta and slope, regional; one global run).
