#!/usr/bin/env Rscript
# Prepare labelled and unlabelled datasets for PTD helpers
# Saves:
#  - data/ptd_labelled.RDS  (labelled LSMS plots with y_pred_theta, log_yield_usd, downstream)
#  - data/ptd_unlabelled.RDS (basin-level unlabelled rows with y_pred_theta, downstream, covs, year)
#
# y_pred_theta = theta_hat * max_EVI_16_cci_c_broad, where theta_hat is the
# EVI-to-log-yield slope from:
#   log(yield_value_USD) ~ max_EVI_16_cci_c_broad | country^year + main_crop
# estimated on the LSMS labelled data (same spec as 42_survey.R lines 224-226).
# Only the slope is used for out-of-sample prediction; country×year and crop
# fixed effects are LSMS-specific and cannot be transported to basin observations.

library("tidyverse")
library("sf")
library("fixest")
library("haven")

message("Preparing labelled and unlabelled PTD datasets...")

if (!file.exists("data/df_reg.RDS")) stop("data/df_reg.RDS missing")
df_reg_raw <- readRDS("data/df_reg.RDS")

# Apply same rescaling as 31_regression_main.R
df_reg <- df_reg_raw |>
  dplyr::mutate(
    elevation                    = elevation / 100,
    pre_chirps                   = pre_chirps / 100,
    accessibility_to_cities_2015 = accessibility_to_cities_2015 / 100,
    pop_2015                     = pop_2015 / 1e5
  )

# Covariate names (must match 31_regression_main.R)
covs <- c("elevation", "slope", "soilgrid_grouped",
          "tmax_tc", "pre_chirps",
          "accessibility_to_cities_2015", "pop_2015")

# ── Read LSMS plot data ───────────────────────────────────────────────────────
pd <- tryCatch(haven::read_dta("data/additional/lsms_plot.dta"), error = function(e) {
  message("Could not read lsms_plot.dta: ", conditionMessage(e)); NULL
})

# Derive year from harvest_interview_month (Stata months since 1960)
if (!is.null(pd)) {
  pd <- pd %>% dplyr::mutate(year = 1960L + as.integer(harvest_interview_month %/% 12))
}

evi_csv <- tryCatch(readr::read_csv("data/additional/63_evi-mean_c_broad_6km2.csv",
                                    show_col_types = FALSE), error = function(e) {
  message("Could not read 63_evi-mean_c_broad_6km2.csv: ", conditionMessage(e)); NULL
})
basins_sf <- tryCatch(sf::st_read("data/df_reg_basins_2023.gpkg", quiet = TRUE), error = function(e) {
  message("Could not read df_reg_basins_2023.gpkg: ", conditionMessage(e)); NULL
})

if (is.null(pd) || is.null(evi_csv) || is.null(basins_sf)) {
  stop("Missing required LSMS/EVI/basins files; cannot prepare labelled/unlabelled datasets.")
}

# ── Build EVI per geocoords_id × year (max 16-day composite, same as 42_survey.R) ──
evi_lsms_cluster_data <- evi_csv %>%
  filter(!is.na(geocoords_id)) %>%
  transmute(geocoords_id,
            year                   = lubridate::year(image_date),
            max_EVI_16_cci_c_broad = mean_EVI) %>%
  arrange(geocoords_id, year) %>%
  group_by(geocoords_id, year) %>%
  slice_max(max_EVI_16_cci_c_broad, n = 1, na_rm = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

lsms_evi_discs <- dplyr::left_join(pd, evi_lsms_cluster_data, by = c("geocoords_id", "year"))

# ── Fit theta model: same spec as 42_survey.R lines 224-226 ──────────────────
# log(yield_value_USD) ~ max_EVI_16_cci_c_broad | country^year + main_crop
# This gives theta_hat = the EVI-to-log-yield elasticity
lsms_model_data <- lsms_evi_discs %>%
  dplyr::filter(
    !is.na(max_EVI_16_cci_c_broad),
    !is.na(yield_value_USD),
    is.finite(yield_value_USD),
    yield_value_USD > 0,
    !is.na(country), !is.na(year), !is.na(main_crop)
  ) %>%
  dplyr::mutate(
    yield_value_USD = {
      x <- yield_value_USD
      qs <- quantile(x, probs = c(0.01, 0.99), na.rm = TRUE)
      x[!is.na(x) & x < qs[1]] <- qs[1]
      x[!is.na(x) & x > qs[2]] <- qs[2]
      x
    },
    log_yield_usd   = log(yield_value_USD)
  )

theta_mod <- fixest::feols(
  log_yield_usd ~ max_EVI_16_cci_c_broad | country^year + main_crop,
  data = lsms_model_data
)

theta_est <- coef(theta_mod)[["max_EVI_16_cci_c_broad"]]
theta_se  <- se(theta_mod)[["max_EVI_16_cci_c_broad"]]
message(sprintf("Fitted theta = %.4f (SE %.4f) from %d LSMS observations",
                theta_est, theta_se, nobs(theta_mod)))

# ── Generate y_pred_theta = theta_hat * EVI ──────────────────────────────────
# Slope-only prediction: country×year and crop FEs from the LSMS model are
# unit-specific and not available for basin observations, so they are dropped.
# The PTD proxy is simply theta_hat * EVI; the dropped FEs contribute only an
# additive constant which cancels out in the debiasing regression.
df_reg <- df_reg %>%
  dplyr::mutate(y_pred_theta = theta_est * max_EVI_16_cci_c_broad)

lsms_model_data <- lsms_model_data %>%
  dplyr::mutate(y_pred_theta = theta_est * max_EVI_16_cci_c_broad)

# ── Spatial join: LSMS plots → basins (RDD window dist_order ∈ {-1, +1}) ─────
coords_pd <- lsms_evi_discs %>%
  dplyr::select(geocoords_id, lat_modified, lon_modified) %>%
  distinct(geocoords_id, .keep_all = TRUE) %>%
  filter(!is.na(lat_modified), !is.na(lon_modified))

coords_sf    <- sf::st_as_sf(coords_pd, coords = c("lon_modified", "lat_modified"), crs = 4326)
coords_sf    <- sf::st_transform(coords_sf, crs = sf::st_crs(basins_sf))
coords_basin <- sf::st_join(
    coords_sf,
    basins_sf %>% dplyr::filter(dist_order %in% c(-1, 1)),
    join = sf::st_within, left = FALSE
  ) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(geocoords_id, HYBAS_ID, mine_basin, dist_order, downstream)

pd_basin <- lsms_model_data %>%
  dplyr::inner_join(coords_basin, by = "geocoords_id")

message(sprintf("PTD labelled: %d plots in RDD window (downstream=1: %d, upstream=0: %d)",
                nrow(pd_basin),
                sum(pd_basin$downstream == 1, na.rm = TRUE),
                sum(pd_basin$downstream == 0, na.rm = TRUE)))

# ── Labelled dataset ──────────────────────────────────────────────────────────
# Join basin-level covariates from df_reg onto LSMS plots via HYBAS_ID.
# Basin covariates are time-invariant so we take one row per HYBAS_ID.
basin_covs <- df_reg %>%
  dplyr::distinct(HYBAS_ID, .keep_all = TRUE) %>%
  dplyr::select(HYBAS_ID, dplyr::any_of(covs))

# Identify which covs were successfully extracted (guards against name mismatches)
covs_in_basin <- intersect(covs, names(basin_covs))
if (length(covs_in_basin) < length(covs)) {
  message(sprintf("Note: %d covs not found in df_reg and will be absent from ptd_labelled: %s",
                  length(covs) - length(covs_in_basin),
                  paste(setdiff(covs, covs_in_basin), collapse = ", ")))
}

ptd_labelled <- pd_basin %>%
  dplyr::filter(!is.na(downstream), !is.na(log_yield_usd), !is.na(y_pred_theta)) %>%
  dplyr::left_join(basin_covs, by = "HYBAS_ID") %>%
  dplyr::select(geocoords_id, HYBAS_ID, mine_basin, year,
                max_EVI_16_cci_c_broad, y_pred_theta,
                log_yield_usd, downstream,
                dplyr::any_of(covs))

# ── Unlabelled dataset: full basin panel from df_reg ─────────────────────────
ptd_unlabelled <- df_reg %>%
  dplyr::select(HYBAS_ID, mine_basin, year, max_EVI_16_cci_c_broad, downstream,
                dist_order, y_pred_theta, dplyr::all_of(covs)) %>%
  dplyr::filter(!is.na(max_EVI_16_cci_c_broad), !is.na(downstream))

# ── Save ──────────────────────────────────────────────────────────────────────
if (!dir.exists("data")) dir.create("data", recursive = TRUE)
saveRDS(ptd_labelled,  file = "data/ptd_labelled.RDS")
saveRDS(ptd_unlabelled, file = "data/ptd_unlabelled.RDS")

message(sprintf("Saved labelled:   %d rows -> data/ptd_labelled.RDS",  nrow(ptd_labelled)))
message(sprintf("Saved unlabelled: %d rows -> data/ptd_unlabelled.RDS", nrow(ptd_unlabelled)))

invisible(list(labelled = ptd_labelled, unlabelled = ptd_unlabelled,
               theta_est = theta_est, theta_se = theta_se, theta_mod = theta_mod))
