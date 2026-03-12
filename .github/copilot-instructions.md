## Quick orientation for AI assistants

This repo contains a small, self-contained analysis pipeline that implements
Proxy-based Treatment Debiasing (PTD) and an adversarial debiasing baseline
(Gordon et al. 2026) applied to a mines → agriculture regression dataset.

Keep the instructions concise and actionable: below are the minimal facts an
assistant needs to make safe, correct edits or extend the analysis.

### Where to start (entry points)
- `run_all.R` — top-level runner that sources the main scripts in `code/`.
- `code/01_prepare-data-ddml-minesriversyields.R` — builds `data/ptd_labelled.RDS`
  and `data/ptd_unlabelled.RDS` from raw inputs. This MUST run before the report.
- `code/02_analysis-ddml-minesriversyields.R` — reads prepared RDS files,
  runs `ptd_debias()` / `adversarial_debias()` and writes `output/plots/ptd_overlay.pdf`.
- `code/helpers_ddml.R` — canonical implementations of PTD, adversarial loss,
  and bootstrap helpers. Prefer changes here when modifying estimator behaviour.
- `code/sim_ddml.R` — synthetic validation / unit-like simulations (no external data).

### Required data files (scripts will stop early with a clear message)
- `data/df_reg.RDS` (main panel)
- `data/df_reg_basins_2023.gpkg` (basins geometry used for spatial joins)
- `data/additional/lsms_plot.dta` and `data/additional/63_evi-mean_c_broad_6km2.csv`
  (used to estimate the θ slope and build `y_pred_theta`).

If any of these are missing, `code/01_prepare...R` explicitly stops with an
informative error; don't try to silently mock them unless writing dedicated
tests or fast CI stubs.

### Naming & data conventions to follow
- Proxy prediction column: `y_pred_theta` (constructed as θ̂ × EVI). Keep this name.
- Labelled (survey) dataset saved to `data/ptd_labelled.RDS` and contains
  `log_yield_usd`, `y_pred_theta`, `downstream`, `mine_basin`, `year`.
- Unlabelled (basin panel) saved to `data/ptd_unlabelled.RDS` and contains
  `HYBAS_ID`, `mine_basin`, `year`, `y_pred_theta`, `downstream`, plus covs.
- Key covariate names used across scripts: `max_EVI_16_cci_c_broad`,
  `elevation`, `slope`, `soilgrid_grouped`, `tmax_tc`, `pre_chirps`,
  `accessibility_to_cities_2015`, `pop_2015`. Note some are rescaled in code
  (e.g. `elevation = elevation / 100`). Mirror this scaling when adding tests.

### Important API / functions (in `code/helpers_ddml.R`)
- ptd_debias(y_true_name, y_pred_name, x_name, labelled_data, unlabelled_data,
  W_names = NULL, use_feols = FALSE, cluster = NULL, fe = NULL, fe_lab = NULL)
  — returns `beta_raw`, `gamma_bias`, `gamma_se`, `beta_corrected` and models.
  Use `fe_lab` when the labelled-data FE would be collinear (this file already
  follows that convention).
- bootstrap_ptd_se(...) — resamples labelled & unlabelled jointly; default
  scripts call with B=100 (faster) but production runs may increase B.
- adversarial_debias(...) — performs optim(BFGS) on labelled data then
  predicts debiased outputs for the unlabelled panel. Check `fit$convergence`.
- bootstrap_adversarial_se(...) — bootstraps both datasets to capture
  optimiser + regression uncertainty (B=100 in the report scripts).

When implementing new estimators, prefer adding small helper functions in
`helpers_ddml.R` and call them from `02_analysis...R` to preserve the single
report entry point.

### Fixed effects, clustering & inference
- `fixest::feols()` is used for FE regressions and clustered SEs. Caller code
  builds formulas like `y ~ downstream + covs | year` and passes `vcov = ~mine_basin`.
- In the labelled-data bias regression the code purposely allows a reduced FE
  spec (`fe_lab`) to avoid perfect collinearity — keep this pattern.

### Developer workflows & common commands
Run the whole pipeline from the project root (or the `ptd_analysis/` folder if
you've moved the scripts) so relative paths resolve:

```zsh
cd /path/to/repo
Rscript run_all.R
```

Run pieces interactively in R/RStudio (useful while editing):

```r
setwd("/path/to/repo")
source("code/01_prepare-data-ddml-minesriversyields.R")   # create RDS files
source("code/02_analysis-ddml-minesriversyields.R")     # run estimators + plot
source("code/sim_ddml.R")                               # quick synthetic checks
```

Notes: bootstrapping (B ≥ 100) can take minutes; CI code in `02_analysis...R`
uses `B = 100` for speed. If you change B, update any runtime docs or CI budgets.

### Examples of patterns to preserve in edits
- Always preserve column names used by helpers (`y_pred_theta`, `log_yield_usd`, `downstream`, `mine_basin`, `HYBAS_ID`).
- When adding covariates, add them to the `covs` vector in both prepare and
  analysis scripts, and ensure they exist in both labelled & unlabelled datasets
  before including them in `W_names`.
- If modifying FE or clustering logic, mirror changes between `ptd_debias()` and
  `adversarial_debias()` so the inference remains comparable.

### Outputs to check after edits
- `data/ptd_labelled.RDS`, `data/ptd_unlabelled.RDS` (size/columns)
- `output/plots/ptd_overlay.pdf` (visual sanity check)
- console messages the scripts print (they intentionally surface missing-file
  problems and optimisation convergence warnings)

If anything here is unclear or you want more project-specific rules (naming
conventions, acceptable library versions, or CI commands), tell me which area
to expand and I'll update this file.
