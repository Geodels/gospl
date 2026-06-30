import os 
import time 
from scripts import analysis as anlys
from gospl.model import Model as sim

def run_model(inputfile):

    # Reading input file
    model = sim(inputfile, False, False)

    # Running forward model
    model.runProcesses()

    # Cleaning model
    model.destroy()

Lx = 1e3   # 1 km
Ly = 1e3   # 1 km
nx = 21
ny = 21
noise_amp = 0.1

# ══════════════════════════════════════════════════════════════════════════════
# Model parameters — edit these to match your goSPL run
# ══════════════════════════════════════════════════════════════════════════════
U     = 5e-4   # m/yr  uplift rate
KAPPA = 0.1    # m²/yr hillslope diffusivity
L     = 1000.0 # m     domain length (x direction)

# ══════════════════════════════════════════════════════════════════════════════
# Test tolerances
# ══════════════════════════════════════════════════════════════════════════════
TOL_TAU_PCT      = 10.0   # % — fitted τ vs theoretical τ
TOL_PEAK_PCT     = 5.0    # % — asymptotic z_ss vs analytical z_max
                           #     (relaxed to account for node offset from x=L/2)
TOL_PROFILE_R2   = 0.990  # R² of spatial profile vs parabola
TOL_CONVERGENCE  = 99.0   # % — model must reach this % of its own asymptote
TOL_MONO         = True   # z_max must increase monotonically
TOL_EARLY_PCT    = 8.0    # % — early-time growth rate (z_ss/τ from the fit).
                           #     Ratio of two already-checked quantities, so it
                           #     compounds the peak (±5%) and τ (±10%) tolerances;
                           #     with the true no-flux 'w' wall the peak sits ~3.5%
                           #     low and τ ~2.5% high, summing to ~6% on the rate.

anlys.configure(U,KAPPA,L)
anlys.tolerances(TOL_TAU_PCT, TOL_PEAK_PCT, TOL_PROFILE_R2, TOL_CONVERGENCE, TOL_MONO, TOL_EARLY_PCT)  

# If initial conditions are needed, make them with the following command and then edit the input path in the next line to match your goSPL run
# input_path = 'based'
# anlys.makeInputs(Lx, Ly, nx, ny,noise_amp, input_path)

outfolder = 'sims'
yml_files = ['hillslope_base.yml', 'hillslope_dx25.yml', 'hillslope_dx10.yml', 
            'hillslope_dt_5000.yml', 'hillslope_dt_500.yml', 
            'hillslope_kappa2.yml', 'hillslope_U2.yml']

os.makedirs('sims_outputs', exist_ok=True)

time_start = time.time()

for file in yml_files:
    run_model(f'{outfolder}/{file}')

time_end = time.time()
print(f"Total execution time: {time_end - time_start} seconds")

# ── Tolerances (set once) ─────────────────────────────────────────────────────
anlys.tolerances(TOL_TAU_PCT, TOL_PEAK_PCT, TOL_PROFILE_R2,
                 TOL_CONVERGENCE, TOL_MONO, TOL_EARLY_PCT)

# ── Step 1: extract each run (configure before each group) ───────────────────
anlys.configure(U=5e-4, KAPPA=0.1, L=1000.0)
r_base   = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_base",    plotdir="results")
r_dx25   = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_dx25",    plotdir="results")
r_dx10   = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_dx10",    plotdir="results")
r_dt5000 = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_dt_5000", plotdir="results")
r_dt500  = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_dt_500",  plotdir="results")

anlys.configure(U=1e-3, KAPPA=0.1, L=1000.0)
r_2U     = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_U2",      plotdir="results")

anlys.configure(U=5e-4, KAPPA=0.2, L=1000.0)
r_kappa2 = anlys.extractOutputs(80, tol=20, model_path="sims_outputs/hillslope_kappa2",  plotdir="results")

# ── Step 2: combined summary ──────────────────────────────────────────────────
anlys.configure(U=5e-4, KAPPA=0.1, L=1000.0)   # restore base for any _set reads
results = anlys.runFullBenchmark(
    results_base   = r_base,
    results_2U     = r_2U,
    results_kappa2 = r_kappa2,
    results_dx25   = r_dx25,
    results_dx10   = r_dx10,
    results_dt5000 = r_dt5000,
    results_dt500  = r_dt500,
    plotdir        = "results",
)

os.system("rm -rf sims_outputs")