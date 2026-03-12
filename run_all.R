#!/usr/bin/env Rscript
# ── Entry point: run all PTD analysis from the project root ──────────────────
# Set working directory to ptd_analysis/ before running, or run via:
#   Rscript ptd_analysis/run_all.R   (from minesriversyields/)
#
# Order:
#   1. remote_sensing_script.R  -- synthetic-data validation (no external data needed)
#   2. 49_ptd_prepare_data.R    -- builds data/ptd_labelled.RDS + ptd_unlabelled.RDS
#   3. 49_ptd_report.R          -- runs PTD + adversarial, saves output/plots/ptd_overlay.pdf

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # works in RStudio
# If running via Rscript, uncomment the next line instead:
# setwd("/Users/thomasbombarde/Documents/R/double-debiased-ml")

message("\n====== Step 1: Synthetic validation ======")
source("code/sim_ddml.R")

message("\n====== Step 2: Prepare PTD data ======")
source("code/01_prepare-data-ddml-minesriversyields.R")

message("\n====== Step 3: PTD + Adversarial report ======")
source("code/02_analysis-ddml-minesriversyields.R")

message("\nAll done. Outputs in output/plots/ and output/tables/")
