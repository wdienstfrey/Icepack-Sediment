# Icepack-Sediment

Prescribed terminus advance changes glacier geometry and frontal stress.
This minimal UW Icepack/Firedrake SSA experiment tests the resulting slowdown
and thickening of Eqalorutsit Kangilliit Sermiat and implements
the comparison in [Dachauer et al. (2026), Figure 13](https://doi.org/10.1017/jog.2026.10172).
Sediment evolution is the next step, after the geometric response is established.

**Status:** implementation and synthetic solver checks; no calibrated historical
EKaS reproduction or agreement claim. See [docs/results.pdf](docs/results.pdf) for validation
and limitations. Regional inputs and the corridor domain still need preparation; see [data/README.md](data/README.md).
[OMP-OOI](https://github.com/wdienstfrey/OMP-OOI) supplies the byte-identical
`AGENTS.md` and `brad-lipovsky-academic-style-guide.md`. It has no standalone
plotting style file; `src/plotting.py` reuses the exact rcParams block from its
`nlayer-sensitivity/energy-likelihood/scripts/analyze.py`, with 220 dpi exports.
The repository is public; all changes remain on the feature branch in PR #1.

The paper initializes from 2021 geometry and velocity, inverts A and slipperiness,
relaxes for two years and re-inverts. It then truncates to the 1953 front, relaxes
for 25 years, and compares a fixed-front control with advance to 2021. Here both
runs independently execute an identical 25-year spinup on the truncated domain;
validation checks the start and post-spinup fields. Time before 1953 is a model
relaxation clock, not a reconstruction of 1928–1953 climate.

`src/model.py` solves `IceStream` stress balance and ice continuity. Units are
meters, years and megapascals. Glen n=3, prescribed A and the prepared basal
friction field are fixed. The mixed Coulomb–Weertman law approximates Ua A.2
through a 16-point quadrature dissipation potential (m=3). Effective pressure is ice overburden
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

Additional checks use cached regional geometry without claiming a reconstruction:

```bash
python scripts/check_regional_geometry.py /path/to/BedMachineGreenland-v6.nc /path/to/ekas_flowline_shiver_2023.csv
```

This probes a candidate straight corridor against the ice mask. The modern
flowline and approximate historical offset do not define a surveyed 1953 front.
The tested candidate includes land with zero ice thickness, so it is unsuitable
for the present solver. A curved ice-conforming mesh or a verified narrower
trunk domain is needed before interpreting regional simulations. See the report
for extended synthetic runs, timestep failures and data limitations. Partial
histories and failure details are saved when a transient aborts.

To reproduce the extended synthetic stability checks, first generate the synthetic
inputs, then run each case with `--inputs data/processed/synthetic.npz` and
`--config configs/synthetic_coarse.json` or `configs/synthetic_refined.json`.
Use distinct `--output` directories. Both configurations retain the 25-year
spinup and 1953–2021 chronology; only the timestep differs (1 versus 0.25 year).
These experiments assess the numerical implementation and remain separate from
the historical EKaS comparison.

The default spatial quadrature degree is 2, matching Icepack for the CG1 spaces,
and can be changed through the model hook.
An experimental analytic-potential comparison at degrees 12 and 24 agreed
to about 5e-13 in short-run gate speed and flux. It was not adopted in production
because a runtime improvement was not established. This is not a mesh or
timestep convergence study. Traction
checks cover 16 pressure/speed combinations, including zero effective pressure.
The one-year timestep test fails despite closing mass budgets: continuous CG1
transport and the split stress/continuity update do not guarantee positive
thickness. No minimum-thickness constraint is used to hide that failure.

The retained full-duration synthetic pair passed: advancing-minus-control gate
thickness +67.83 m, speed −515.54 m/year and flux −0.18842 Gt/year in 2021.
Both remained positive over all 93 model years; this does not validate regional EKaS magnitudes.
