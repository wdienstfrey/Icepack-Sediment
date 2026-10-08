# Icepack-Sediment

Can prescribed terminus advance over a shallowing bed slow and thicken Eqalorutsit
Kangilliit Sermiat? This minimal UW Icepack/Firedrake SSA experiment implements
the comparison in [Dachauer et al. (2026), Figure 13](https://doi.org/10.1017/jog.2026.10172).
Sediment evolution is the next step, after the geometric response is established.

**Status:** implementation and synthetic solver checks; no calibrated historical
EKaS reproduction or agreement claim. See `docs/results.pdf` for validation
and limitations. Exact regional inputs still need preparation; see [data/README.md](data/README.md).
The requested `OOI-Image` template could not be located; its AGENTS.md and
plotting style must be supplied before the template requirement can be completed.

The paper initializes from 2021 geometry and velocity, inverts A and slipperiness,
relaxes for two years and re-inverts. It then truncates to the 1953 front, relaxes
for 25 years, and compares a fixed-front control with advance to 2021. Here both
runs independently execute an identical 25-year spinup on the truncated domain;
validation checks the start and post-spinup fields. Time before 1953 is a model
relaxation clock, not a reconstruction of 1928–1953 climate.

`src/model.py` solves `IceStream` stress balance and ice continuity. Units are
meters, years and megapascals. Glen n=3, prescribed A and the prepared basal
friction field are fixed. The mixed Coulomb–Weertman law approximates Ua A.2
through a quadrature dissipation potential. Effective pressure is ice overburden
minus ocean-connected hydrostatic pressure at the ice base, bounded at zero;
floating ice therefore has zero basal drag. Icepack supplies the grounded/floating
surface and ocean traction. The fixed MAR 2021 annual SMB is identical in both
runs. Upstream velocity/thickness are prescribed from inputs; side walls have
Icepack's no-normal-flow penalty and configurable sliding drag. No basal melt,
calving law or explicit sediments are included.

Front offsets live in [data/terminus.csv](data/terminus.csv). They are approximate
integrals of the Figure 13 advance rates, and are editable. Monotone PCHIP makes
front motion differentiable at observations; set `smooth_front` false for the
paper's piecewise-linear position. The initial front is x=`length_1953_m`.
A fresh, identically ordered rectangular triangular mesh stretches in x each
step. Nodal thickness is multiplied by the old/new length ratio to conserve
volume on transfer, then transported with ice velocity minus mesh velocity.
This first-order ALE split requires convergence checks. It creates no arbitrary
new ice slab, and records the implied ice discharge relative to the moving front.
It assumes a straight constant-width corridor; it does not reproduce the
paper's curved planform fronts or full catchment.

Install Firedrake/PETSc using its [official installation instructions](https://www.firedrakeproject.org/install.html),
then install UW [Icepack](https://icepack.github.io/install/) in that environment.
Do not install the unrelated CICE sea-ice package named Icepack. API reference
checked against icepack commit `c9a29780cd0f7d068d206cb5a170fa367a7655b0`.
`environment.yml` provides a separate preprocessing/plotting environment and
**does not install Firedrake**. Add numpy/scipy/matplotlib/rasterio to the
Firedrake environment if needed. Run from the repository root, on one MPI rank:

```bash
python scripts/prepare_inputs.py data/raw/manifest.json
python scripts/run_fixed_front.py
python scripts/run_advancing_front.py
python scripts/validate_results.py
python scripts/make_figure13.py
jupyter lab notebooks/figure13_reproduction.ipynb
```

Scientific outputs are ignored by git. CSV histories include gate thickness,
speed and flux at x=1953-front minus 1 km, domain volume, minimum thickness and
mass residual. NPZ snapshots on the shared 1953 domain allow spatial differences.
The figure script draws the Figure 13 time-series comparison plus thickness,
speed and spatial-difference panels. It refuses missing results. A separate
terminus panel can be generated without simulations; no modeled values are
filled in when results do not exist.

Fast checks and a short, explicitly synthetic solver run:

```bash
python -m unittest discover -s tests -v
python tests/check_firedrake_geometry.py  # in the Firedrake environment
python scripts/prepare_smoke_inputs.py
python scripts/run_fixed_front.py --inputs data/processed/synthetic.npz --config data/processed/smoke_config.json --output results/smoke_fixed
python scripts/run_advancing_front.py --inputs data/processed/synthetic.npz --config data/processed/smoke_config.json --output results/smoke_advancing
python scripts/validate_results.py --fixed results/smoke_fixed --advancing results/smoke_advancing --allow-short
```

For scientific validation, halve dt and mesh spacing and compare gate curves,
mass budgets and maps. Reject nonfinite/nonpositive thickness rather than clipping.
The complete-run validator checks slowdown and thickening relative to control;
a failed expectation is a scientific result to investigate, not something to
force through parameter constraints. Short smoke runs do not establish historical
stability or the expected 68-year response.

Known differences: no A/C inversion or re-inversion; truncated trunk inflow rather
than ice-divide boundaries; uniform A by default; stock/user-corrected bed rather
than the author's kriged product; approximate scalar front offsets; corridor
stretching rather than Ua remeshing; regularized quadrature sliding potential;
smooth rather than piecewise-linear front motion by default. These differences
can alter both the sign and magnitude of the response. No synthetic run is
presented as EKaS or used to infer agreement with the paper.

Planned next step: add explicit sediment thickness and bank geometry to the bed
and ocean-contact boundary, with mass conservation and independently constrained
sediment parameters, after establishing the geometric control experiment.

Rebuild the report with `cd docs` and `pdflatex results.tex` (twice for links).
