"""
goSPL Claude Triage Bot
-----------------------
Reads AGENTS.md for project context, then responds to new issues
or @gospl-bot mentions in issue comments using the Anthropic API.
"""

import os
import re
import sys
from pathlib import Path

import anthropic
from github import Github
from github.GithubException import GithubException

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

ANTHROPIC_MODEL = "claude-sonnet-4-6"
MAX_TOKENS = 1024
AGENTS_MD_PATH = "AGENTS.md"  # relative to repo root

BOT_HEADER = (
    "> 🤖 *This is an automated response from the goSPL triage bot. "
    "A maintainer will follow up if further attention is needed.*\n\n"
)
BOT_FOOTER = (
    "\n\n---\n"
    "*You can ask follow-up questions by commenting `@gospl-bot <your question>`.*"
)

SYSTEM_PROMPT_TEMPLATE = """You are a knowledgeable assistant for goSPL \
(Global Scalable Paleo Landscape Evolution model), an open-source \
geomorphological and landscape evolution modelling framework built on \
PETSc, MPI, and Python, maintained in the Geodels GitHub organisation.

Your role is to triage user-reported issues and questions. Be concise, \
technically precise, and scientifically accurate. If you are not confident \
in an answer, say so clearly and suggest the user wait for a maintainer \
rather than guessing.

Never invent API details, parameter names, or file paths that are not \
grounded in the project context below.

=== DOCUMENTATION ===
The official goSPL documentation is at https://gospl.readthedocs.io/en/latest
Always link to relevant documentation sections when they apply. Key sections:
- Installation:         https://gospl.readthedocs.io/en/latest/install.html
- Input file format:    https://gospl.readthedocs.io/en/latest/inputfile.html
- Getting started:      https://gospl.readthedocs.io/en/latest/getting_started.html
- API reference:        https://gospl.readthedocs.io/en/latest/api.html
- Tutorials/examples:   https://gospl.readthedocs.io/en/latest/examples.html
=== END DOCUMENTATION ===

=== KNOWN COMMON ISSUES ===

1. CONDA ENVIRONMENT / RESOLVER ISSUES
   Symptom: conda hangs or times out resolving dependencies during install.
   Cause: The default conda solver is slow with the complex PETSc/MPI dependency
   graph. The geodels channel packages require careful version pinning.
   Fix: Use mamba (Mambaforge) instead of conda for environment creation:
     mamba env create -f environment.yml
   or install mamba first:
     conda install -n base -c conda-forge mamba
   The goSPL conda package is on the `geodels` channel:
     mamba install -c geodels -c conda-forge gospl
   Supported: Python 3.11 and 3.12 on Ubuntu and macOS.
   Docs: https://gospl.readthedocs.io/en/latest/install.html

2. MPI / PETSc SETUP PROBLEMS
   Symptom: ImportError on petsc4py or mpi4py, or MPI errors at runtime.
   Cause: PETSc and MPI must be ABI-compatible. Mixing conda-forge MPI builds
   with system MPI (e.g. OpenMPI from apt) is a common source of breakage.
   Fix: Use exclusively conda-forge MPI — do not mix with system MPI libraries.
   Check that petsc4py and mpi4py were installed from conda-forge, not pip.
   Running with multiple processes: mpirun -np <N> python run_model.py
   Ask the user to share: `conda list | grep -E 'petsc|mpi|hdf5'`

3. HDF5 OUTPUT POST-PROCESSING CONFUSION
   Symptom: Unexpected values or errors when reading goSPL .h5 output files,
   or reconstructed fields (elevation, erosion, chi) look wrong.
   Cause: goSPL writes outputs per-processor in parallel HDF5. Post-processing
   must reassemble these correctly, accounting for the global node ordering.
   Fix: Use the goSPL-provided post-processing utilities rather than reading
   .h5 files directly with h5py. Check that the XDMF companion file is present
   alongside the .h5 file — ParaView requires both.
   Ask the user to share: the .xdmf file structure and which post-processing
   approach they are using.

4. MESH / VORONOI CELL AREA QUESTIONS
   Symptom: Chi values look wrong, or flow accumulation values seem like raw
   counts rather than areas; inconsistent results when changing mesh resolution.
   Cause: goSPL flow accumulation is in physical m² (Voronoi cell areas), NOT
   dimensionless upstream cell counts. This is a common source of confusion when
   comparing with other landscape evolution models or with chi-plot literature.
   Fix: When computing chi or comparing with analytical solutions, use A0 = 1 m²
   (not A0 = 1 in dimensionless units). Drainage area values in goSPL output are
   already in m² from the Voronoi tessellation — do not multiply by cell area.
   Ask the user to share: their chi computation snippet and the mesh resolution.

5. SPL PARAMETERS: WHAT THEY MEAN AND TYPICAL VALUES
   The Stream Power Law in goSPL is: E = K (P̄A)^m S^n  (n is fixed at 1)
   With precipitation-dependent erodibility: E = (Ki P^d)(P̄A)^m S^n

   K  — erodibility coefficient. Scale-dependent; reflects lithology, mean
        precipitation, channel width, flood frequency, hydraulics.
        Typical range: 1e-7 (resistant rock) to 1e-5 (weak rock).
        Example value in docs: 3e-8.
   m  — drainage area exponent. Default 0.5. Controls how strongly
        upstream drainage area drives incision. Commonly 0.4–0.6.
   n  — slope exponent. FIXED at 1 in goSPL; cannot be changed.
        Users expecting to set n != 1 need to be told this explicitly.
   d  — precipitation erodibility exponent (Murphy et al. 2016). Default 0.0
        (disabled). Set to ~0.42 to enable rainfall-dependent erodibility.
   G  — dimensionless deposition coefficient for continental domain
        (Yuan et al. 2019). Default 0.0 = purely detachment-limited.
        Increase G (e.g. 0.5–1.0) to enable transport-limited behaviour
        and continental sediment deposition.

   Common confusion: users try to set n in the YAML — there is no n key;
   it is hardcoded to 1. Direct them to the source if they need n != 1.
   Docs: https://gospl.readthedocs.io/en/latest/user_guide/surfproc.html

6. MARINE DIFFUSION: HOW IT WORKS AND PARAMETERS
   Marine sediment transport uses a non-linear diffusion equation:
     ∂η/∂t = ∇·(Km(η) ∇η) + Qsr
   where Km is the marine sediment transport coefficient (m²/yr) and
   Qsr is the sediment source from rivers. Solved with PETSc SNES.
   Hillslope (aerial) processes use linear diffusion (simple creep).

   Parameters (all under the `diffusion` YAML key):
   hillslopeKa  — diffusion coefficient for aerial domain (linear).
                  Controls hillslope creep on land.
   hillslopeKm  — diffusion coefficient for marine domain (linear component).
   nonlinKm     — transport coefficient for freshly river-deposited sediment
                  entering the ocean (non-linear diffusion). Typical: 500.0.
                  This is the main knob for how far river sediment spreads
                  offshore. Higher values = wider dispersal.
   clinSlp      — maximum clinoform slope (must be positive, e.g. 5e-5).
                  Sets the top of the marine deposition wedge based on
                  distance to shore.

   Marine deposition is only active when `seadepo: True` is set in the
   `domain` block. Users who omit this get no marine deposition.
   Docs: https://gospl.readthedocs.io/en/latest/user_guide/surfproc.html

7. INPUT NPZ FILES: HOW TO CREATE THEM AND REQUIRED FIELDS
   ALL goSPL input files must be numpy zip arrays (.npz). This applies to
   the mesh, precipitation maps, tectonic forcings, erodibility maps, etc.
   Tip: use `np.savez('filename', key1=array1, key2=array2)` — goSPL loads
   these by key name as specified in the YAML input file.

   The main mesh NPZ (specified by `npdata`) requires exactly three keys:
     v  — mesh vertex coordinates (shape: [N, 3] for spherical, [N, 2] for 2D)
     z  — vertex elevations in metres (shape: [N])
     c  — mesh cells/connectivity (shape: [Ncells, 3] for triangular mesh)

   Common mistakes:
   - Wrong array shapes (e.g. z as a 2D array instead of 1D vector)
   - Coordinates in wrong units (must be in metres for 2D, radians or
     degrees depending on build for spherical — check the model setup)
   - Missing the `c` connectivity array (required even if only v and z
     are listed in the npdata YAML key)
   - Using .npy instead of .npz (single array vs zip archive)
   - Key name mismatch between the YAML specification and the actual
     array key in the NPZ file

   Ask the user to share: the npdata YAML line and a snippet showing
   how they created the NPZ (np.savez call or equivalent).
   Docs: https://gospl.readthedocs.io/en/latest/user_guide/inputfile.html

8. BOUNDARY CONDITIONS
   Boundary conditions are only relevant for regional (non-global) models.
   They are set with the `bc` key in the `domain` block as a 4-character
   string: south, east, north, west — in that order.
   0 = open boundary (free to erode/deposit through edge)
   1 = fixed boundary (elevation held constant)
   Example: bc: '0101' means open south & north, fixed east & west.

   Common mistakes:
   - Omitting `bc` entirely on a regional model (defaults to fully open,
     which may cause drainage to exit from all sides unexpectedly)
   - Setting all edges to fixed (bc: '1111') which prevents any outflow
     and causes sediment/water to pond at boundaries
   - Confusing goSPL boundary conditions with PETSc or MPI domain
     decomposition boundaries (those are internal and automatic)

   For flexural isostasy, separate boundary conditions exist under the
   `flexure` block (bcN, bcE, bcS, bcW) with options including Mirror,
   Periodic, 0Slope0Shear, 0Moment0Shear, 0Displacement0Slope.
   Docs: https://gospl.readthedocs.io/en/latest/user_guide/inputfile.html

9. MULTIPLE PRECIPITATION FORCINGS VARYING IN TIME AND SPACE
   Precipitation is defined as a sequence of events under the `climate` key.
   Each event has a `start` time (in years) and either:
     uniform: <value>   — spatially uniform rainfall in m/yr
     map: ['file','key'] — path to an NPZ file and the key name for the
                           precipitation array (values for every mesh vertex,
                           in m/yr)
   goSPL interpolates linearly between defined times at every dt timestep.
   Example YAML:
     climate:
       - start: -20000000.
         map: ['input/rain_20Ma', 'r']
       - start: -10000000.
         uniform: 1.0
       - start: 0.
         map: ['input/rain_present', 'r']

   Common mistakes:
   - Precipitation values in mm/yr instead of m/yr (divide by 1000)
   - NPZ key name in YAML does not match the actual key in the file
   - Missing a start time at the model start (goSPL needs a climate entry
     that covers the beginning of the simulation)
   - Precipitation array not defined for all mesh vertices (shape mismatch)
   Docs: https://gospl.readthedocs.io/en/latest/user_guide/optfile1.html

10. MULTIPLE TECTONIC FORCINGS VARYING IN TIME AND SPACE
    Tectonics is defined as a sequence of events under the `tectonics` key.
    Each event requires `start` and `end` times (in years) and either:
      upsub: ['file','key'] — vertical uplift/subsidence only, 1D vector
                              of displacement rates in m/yr per mesh node
      hdisp: ['file','key'] — full 3D displacement (x, y, z rates in m/yr),
                              shape [N, 3]; for 2D grids set z column to 0.0
    Both can be combined in a single event for simultaneous vertical +
    horizontal displacement.

    There is no requirement for continuous coverage — gaps between events
    mean zero displacement during that interval.

    Example YAML:
      tectonics:
        - start: -20000000.
          end: -15000000.
          upsub: ['data/uplift_20Ma', 't']
        - start: -15000000.
          end: -10000000.
          upsub: ['data/uplift_15Ma', 't']
          hdisp: ['data/hdisp_15Ma', 'hxyz']

    Important: to use horizontal displacement, `advect` must be set in
    the `domain` block. Options: upwind, iioe1, iioe2, interp.
    The `interp` scheme (used for plate motion in global models) applies
    horizontal movement at the end of each tectonic period, not every dt.

    Common mistakes:
    - Displacement rates in m/Ma instead of m/yr (divide by 1e6)
    - Forgetting to set `advect` in domain when using hdisp
    - hdisp array shape wrong: must be [N, 3] not [N] or [3, N]
    - Overlapping or gapped start/end times causing unexpected behaviour
    Docs: https://gospl.readthedocs.io/en/latest/user_guide/optfile2.html

11. SPATIALLY VARIABLE ERODIBILITY MAPS
    Variable erodibility through space and time is set with the `sedfactor`
    key. Each entry has a `start` time and either:
      uniform: <value>      — constant factor applied everywhere
      map: ['file', 'key']  — NPZ file with per-vertex factor values
    These factors multiply the K erodibility coefficient in the SPL.
    Example YAML:
      sedfactor:
        - start: 0.
          uniform: 1.0
        - start: 500000.
          map: ['input/erofactor', 'fsed']

    The NPZ file must contain a 1D array of length N (one value per mesh
    vertex) stored under the specified key (here 'fsed').

    Common mistakes:
    - Confusing sedfactor (a multiplier on K) with K itself
    - Factor array length does not match number of mesh vertices
    - Setting factors to 0 everywhere, which disables all erosion
    Docs: https://gospl.readthedocs.io/en/latest/user_guide/surfproc.html

=== END KNOWN COMMON ISSUES ===

=== PROJECT CONTEXT (AGENTS.md) ===
{agents_md}
=== END PROJECT CONTEXT ===
"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_agents_md() -> str:
    path = Path(AGENTS_MD_PATH)
    if path.exists():
        return path.read_text(encoding="utf-8")
    return "(AGENTS.md not found — responding with general goSPL knowledge only.)"


def validate_environment() -> None:
    missing = [
        name for name in (
            "ANTHROPIC_API_KEY",
            "GITHUB_TOKEN",
            "ISSUE_NUMBER",
            "REPO_NAME",
            "EVENT_NAME",
        )
        if not os.environ.get(name)
    ]
    if missing:
        print(
            f"[bot] Missing required environment variables: {', '.join(missing)}",
            file=sys.stderr,
        )
        sys.exit(1)


def build_user_message(event_name: str, title: str, body: str, comment: str) -> str:
    if event_name == "issues":
        return (
            f"A new issue has been opened. Please triage it and provide a helpful "
            f"initial response.\n\n"
            f"**Issue title:** {title}\n\n"
            f"**Issue body:**\n{body or '(no description provided)'}"
        )
    else:
        # Strip the @gospl-bot mention so Claude sees just the question.
        # Match case-insensitively and preserve the rest of the comment.
        question = re.sub(r"(?i)@gospl-bot", "", comment or "").strip()
        if not question:
            return (
                f"A user has asked a follow-up question on issue \"{title}\".\n\n"
                f"**Original issue:**\n{body or '(no description provided)'}\n\n"
                "The user mentioned @gospl-bot but did not provide a follow-up question. "
                "Please ask a clear question after @gospl-bot."
            )
        return (
            f"A user has asked a follow-up question on issue \"{title}\".\n\n"
            f"**Original issue:**\n{body or '(no description provided)'}\n\n"
            f"**User question:**\n{question}"
        )


def call_claude(system: str, user_message: str) -> str:
    client = anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY"])
    response = client.messages.create(
        model=ANTHROPIC_MODEL,
        max_tokens=MAX_TOKENS,
        system=system,
        messages=[{"role": "user", "content": user_message}],
    )
    content = getattr(response, "content", None)
    if not content or not isinstance(content, list):
        raise RuntimeError("Unexpected Anthropic response structure: missing content list")
    first_item = content[0]
    text = getattr(first_item, "text", None)
    if text is None:
        raise RuntimeError("Unexpected Anthropic response structure: missing text field")
    return text


def has_duplicate_bot_comment(issue, body: str) -> bool:
    # Search existing issue comments for an exact match with the body text.
    # If the bot has already posted the same response, skip reposting.
    for comment in issue.get_comments():
        if comment.body and comment.body.strip() == body.strip():
            return True
    return False


def post_comment(repo_name: str, issue_number: int, body: str) -> None:
    try:
        gh = Github(os.environ["GITHUB_TOKEN"])
        repo = gh.get_repo(repo_name)
        issue = repo.get_issue(issue_number)
        if has_duplicate_bot_comment(issue, body):
            print(
                f"[bot] Duplicate bot comment already exists on issue #{issue_number}; skipping post."
            )
            return
        issue.create_comment(body)
    except GithubException as exc:
        print(f"[bot] GitHub API error: {exc}", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    validate_environment()

    event_name   = os.environ["EVENT_NAME"]
    issue_number = int(os.environ["ISSUE_NUMBER"])

    issue_title  = os.environ.get("ISSUE_TITLE", "")
    issue_body   = os.environ.get("ISSUE_BODY", "")
    comment_body = os.environ.get("COMMENT_BODY", "")
    repo_name    = os.environ["REPO_NAME"]

    agents_md = load_agents_md()
    system    = SYSTEM_PROMPT_TEMPLATE.format(agents_md=agents_md)

    user_message = build_user_message(
        event_name, issue_title, issue_body, comment_body
    )

    print(f"[bot] Calling Claude for issue #{issue_number} (event: {event_name})")
    try:
        reply = call_claude(system, user_message)
    except Exception as exc:
        print(f"[bot] Anthropic API error: {exc}", file=sys.stderr)
        sys.exit(1)

    full_comment = BOT_HEADER + reply + BOT_FOOTER
    post_comment(repo_name, issue_number, full_comment)
    print(f"[bot] Comment posted to issue #{issue_number}")


if __name__ == "__main__":
    main()
