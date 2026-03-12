# PTD Bias-Correction Analysis

Self-contained folder for the Proxy-based Treatment Debiasing (PTD) pipeline
(Gordon et al. 2026), applied to the mines → agriculture RDD project.

## Structure

```
ptd_analysis/
├── run_all.R                      # entry point — runs all three scripts in order
├── code/
│   ├── helpers_remote_sensing.R   # ptd_debias(), adversarial_debias(), bootstrap_*()
│   ├── remote_sensing_script.R    # synthetic-data validation (no external data needed)
│   ├── 49_ptd_prepare_data.R      # builds ptd_labelled.RDS + ptd_unlabelled.RDS
│   └── 49_ptd_report.R            # PTD + adversarial estimates, comparison plot
├── data/                          # symlinks to main project data (not duplicated)
│   ├── df_reg.RDS
│   ├── df_reg_basins_2023.gpkg
│   └── additional/
│       ├── lsms_plot.dta
│       └── 63_evi-mean_c_broad_6km2.csv
└── output/
    ├── plots/                     # ptd_overlay.pdf written here
    └── tables/
```

## How to run

**Working directory must be `ptd_analysis/`** for all relative paths to resolve.

### From RStudio
Open `run_all.R` and source it — it sets the working directory automatically.

### From the terminal
```zsh
cd /Users/thomasbombarde/Documents/R/minesriverstyields/minesriversyields/ptd_analysis
Rscript run_all.R
```

### Individual scripts
```r
setwd(".../ptd_analysis")
source("code/49_ptd_prepare_data.R")   # only needs to run once; saves RDS files
source("code/49_ptd_report.R")         # reads the RDS files, runs estimation, saves plot
source("code/remote_sensing_script.R") # fully synthetic, no data needed
```

## Required R packages
```r
install.packages(c("tidyverse", "fixest", "sf", "haven", "broom", "psych", "purrr", "ggplot2"))
```
(`torch` is loaded in `remote_sensing_script.R` but not actively used — safe to skip if not installed.)

## Key outputs

| File | Description |
|---|---|
| `data/ptd_labelled.RDS` | 567 LSMS plots with `log_yield_usd`, `y_pred_theta`, `downstream` |
| `data/ptd_unlabelled.RDS` | ~93k basin-year rows with `y_pred_theta`, covariates |
| `output/plots/ptd_overlay.pdf` | Coefficient comparison: fe-OLS vs adversarial debiasing |

## Methods summary

- **θ model**: `log(yield_USD) ~ max_EVI_16_cci_c_broad | country×year + main_crop` on LSMS data → slope θ̂ used to construct `y_pred_theta = θ̂ × EVI` for all basin-years
- **PTD**: `β_corrected = β_raw − γ̂`, where γ̂ = regression of prediction errors on `downstream` in labelled data; SEs via bootstrap (B=100)
- **Adversarial**: minimises `MSE + α·Cov(ν, downstream)²` over prediction weights on labelled data, then regresses debiased predictions on `downstream | year + mine_basin`; SEs via bootstrap (B=100)
