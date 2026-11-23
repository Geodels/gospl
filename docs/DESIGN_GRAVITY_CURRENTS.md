# Design: marine gravity currents (opt-in `gravity:` block)

Status: phases G1 and G2 done on branch `feat/gravity-currents` (2026-11-23); G3 to do.

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

## 7. Phases

- **G1 (this branch):** split + routed solve + run-out length + fixed/flood trigger,
  single fraction and dual lithology, `gravDep` output, tests, docs.
- **G2 (done):** provenance classes in lockstep; budget in the `--summary` JSON.
- **G3:** slope-failure trigger (remobilise sediment above a critical slope), optional
  slope-dependent deposition, catchment-area-dependent `sigma`; example in
  goSPL-examples (delta and slope, regional; one global run).
