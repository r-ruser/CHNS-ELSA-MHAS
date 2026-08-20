options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(dplyr); library(arrow)})

set.seed(20260817)
project <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
src <- file.path(project, "data", "causal_landmark")
out <- file.path(project, "causal_dml", "data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

imp <- readRDS(file.path(src, "landmark_MICE5_no_outcomes_internal_restricted.rds"))
stopifnot(length(imp) == 5L)

# Outcome access occurs only here, after the baseline/MICE objects already exist.
outcomes <- readRDS(file.path(project, "data", "causal", "causal_long_full_internal_restricted.rds")) |>
  filter(wave == 2009) |>
  select(Idind, HD_z, KDM_BAA) |>
  distinct(Idind, .keep_all = TRUE)

id_map <- imp[[1]] |>
  transmute(Idind, dml_id = sprintf("DML%06d", row_number()))
write.csv(id_map, file.path(out, "internal_id_crosswalk_restricted.csv"), row.names = FALSE)

feature_names <- c(
  "age2004", "sex", "province", "fuel2004", "previous_fuel_state",
  "number_prior_fuel_waves", "proportion_prior_solid_only", "proportion_prior_mixed",
  "ever_clean_before_2006", "ever_mixed_before_2006", "number_fuel_transitions",
  "number_clean_to_mixed_reversals", "number_mixed_to_solid_reversals",
  "years_since_first_clean_observation",
  "duration_weighted_solid_equivalent_history_to_2004", "last_two_wave_pattern",
  "urban_recent", "education_recent", "income_recent", "income_mean", "income_change",
  "assets_recent", "assets_mean", "assets_change", "smoking_recent", "alcohol_recent"
)

manifest <- list()
for (j in seq_along(imp)) {
  z <- imp[[j]] |>
    left_join(id_map, by = "Idind") |>
    left_join(outcomes, by = "Idind") |>
    mutate(
      R = as.integer(!is.na(biomarker2009) & biomarker2009 == 1 &
                     (!is.na(HD_z) | !is.na(KDM_BAA))),
      A = as.integer(A),
      commid2004 = as.character(commid2004),
      imputation = j
    ) |>
    select(dml_id, imputation, commid2004, A, R, all_of(feature_names), HD_z, KDM_BAA)
  stopifnot(nrow(z) == 5366L, !anyDuplicated(z$dml_id), all(z$A %in% 0:1), all(z$R %in% 0:1))
  write_parquet(z, file.path(out, sprintf("dml_imputation_%d.parquet", j)), compression = "zstd")
  manifest[[j]] <- data.frame(
    imputation = j, n = nrow(z), clean = sum(z$A == 1), any_solid = sum(z$A == 0),
    outcome_observed = sum(z$R == 1), communities = n_distinct(z$commid2004)
  )
}

write.csv(bind_rows(manifest), file.path(out, "dml_input_manifest.csv"), row.names = FALSE)
writeLines(feature_names, file.path(out, "frozen_feature_list.txt"), useBytes = TRUE)
writeLines(c(
  "DML_INPUT_PREPARATION=PASS",
  "Eligibility=2004 any-solid; treatment=2006 clean-only versus continued any-solid; outcome=2009.",
  "Five MICE baseline datasets; outcomes joined only after baseline feature list was frozen.",
  "Public analytical exports use dml_id; the Idind crosswalk is internal restricted."
), file.path(out, "20_prepare_dml_inputs.log"), useBytes = TRUE)

