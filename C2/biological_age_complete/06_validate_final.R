options(stringsAsFactors = FALSE, cli.unicode = FALSE)
suppressPackageStartupMessages({library(png); library(tiff)})

project <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(project) != "C2") stop("Run from C2")
formal <- file.path(project, "results", "formal")
fig <- file.path(project, "figures", "Figure2_fuel_trajectory_biological_age")
logcon <- file(file.path(project, "logs", "06_validate_final.log"), open = "wt", encoding = "UTF-8")
sink(logcon, split = TRUE)
sink(logcon, type = "message")
on.exit({sink(type = "message"); sink(); close(logcon)}, add = TRUE)

formal_files <- c(
  "formal_sample_flow.csv", "table1.csv", "main_secondary_results.csv",
  "sensitivity_results.csv", "extended_sensitivity_results.csv",
  "sex_interaction_tests.csv", "selection_balance.csv",
  "selection_weight_diagnostics.csv", "fuel_reclassification_matrix_counts.csv",
  "fuel_reclassification_individual_deidentified.csv", "sap_deviations.csv",
  "covariate_balance.csv", "propensity_weight_diagnostics.csv",
  "positivity_overlap.csv", "kdm_crossfit_diagnostics.csv",
  "figure2_source_data.csv", "model_n_integrity.csv",
  "cluster_bootstrap_replicates.csv", "cluster_bootstrap_summary.csv",
  "cluster_bootstrap_diagnostics.csv", "cluster_bootstrap_warnings.csv",
  "cluster_bootstrap_failures.csv", "cluster_bootstrap_contract.csv"
)
required <- c(
  file.path(project, "data", c(
    "analysis_dataset.rds", "analysis_dataset_deidentified.csv",
    "analysis_dataset_deidentified.rds"
  )),
  file.path(formal, formal_files),
  paste0(fig, c(".svg", ".pdf", ".tiff", ".png")),
  file.path(project, c("study_report.md", "protocol.md", "analysis_plan.md", "go_no_go.md"))
)
stopifnot(all(file.exists(required)), all(file.info(required)$size > 0))

d <- readRDS(file.path(project, "data", "analysis_dataset.rds"))
res <- read.csv(file.path(formal, "main_secondary_results.csv"))
integ <- read.csv(file.path(formal, "model_n_integrity.csv"))
src <- read.csv(file.path(formal, "figure2_source_data.csv"))
stopifnot(nrow(d) == 6544, all(integ$pass), all(integ$value == 4514))
stopifnot(nrow(src) == 10, all(src$n == 4514), all(src$lower <= src$estimate & src$estimate <= src$upper))

hd <- subset(res, outcome == "HD_z" & analysis == "Fully adjusted")
ky <- subset(res, outcome == "KDM_BAA" & analysis == "Fully adjusted")
stopifnot(nrow(hd) == 1, nrow(ky) == 1)
stopifnot(abs(hd$estimate - 0.0495033833) < 1e-8, abs(ky$estimate - 0.6770059474) < 1e-8)

# Exactly 1,000 successful community-cluster bootstrap replicates, with fixed
# constructed outcomes and directly reproducible percentile intervals.
boot <- read.csv(file.path(formal, "cluster_bootstrap_replicates.csv"))
bsum <- read.csv(file.path(formal, "cluster_bootstrap_summary.csv"))
bdiag <- read.csv(file.path(formal, "cluster_bootstrap_diagnostics.csv"))
bwarn <- read.csv(file.path(formal, "cluster_bootstrap_warnings.csv"))
bfail <- read.csv(file.path(formal, "cluster_bootstrap_failures.csv"))
bcontract <- read.csv(file.path(formal, "cluster_bootstrap_contract.csv"))
stopifnot(
  nrow(boot) == 1000L,
  identical(boot$replicate, seq_len(1000L)),
  all(is.finite(boot$HD_z_estimate)),
  all(is.finite(boot$KDM_BAA_estimate)),
  all(boot$n_cluster_draws == 193L),
  all(boot$n_new_cluster_ids == 193L),
  nrow(bsum) == 2L,
  all(bsum$successful_replicates == 1000L),
  all(bsum$attempts == 1000L),
  all(bsum$outcome_construction == "Fixed from original analysis; not re-estimated within bootstrap"),
  all(bdiag$pass),
  nrow(bwarn) == 0L,
  nrow(bfail) == 0L,
  any(bcontract$item == "Primary uncertainty" & grepl("sandwich", bcontract$specification, ignore.case = TRUE)),
  any(bcontract$item == "Outcome construction" & grepl("not re-estimated", bcontract$specification, fixed = TRUE))
)
q_hd <- unname(quantile(boot$HD_z_estimate, c(.025, .975), type = 7))
q_kdm <- unname(quantile(boot$KDM_BAA_estimate, c(.025, .975), type = 7))
s_hd <- subset(bsum, outcome == "HD_z")
s_kdm <- subset(bsum, outcome == "KDM_BAA")
stopifnot(
  max(abs(q_hd - c(s_hd$lower, s_hd$upper))) < 1e-12,
  max(abs(q_kdm - c(s_kdm$lower, s_kdm$upper))) < 1e-12,
  nrow(subset(res, analysis == "Fully adjusted; community cluster bootstrap")) == 2L,
  sum(src$analysis == "Fully adjusted; community cluster bootstrap") == 2L
)

# Cross-fitting and all extended sensitivity outputs must be internally valid.
kd <- read.csv(file.path(formal, "kdm_crossfit_diagnostics.csv"))
stopifnot(nrow(kd) == 10, all(kd$n_train >= 200), all(is.finite(kd$age_rmse)), all(kd$age_correlation > 0.7))
sens <- read.csv(file.path(formal, "sensitivity_results.csv"))
ext <- read.csv(file.path(formal, "extended_sensitivity_results.csv"))
inter <- read.csv(file.path(formal, "sex_interaction_tests.csv"))
selbal <- read.csv(file.path(formal, "selection_balance.csv"))
selwd <- read.csv(file.path(formal, "selection_weight_diagnostics.csv"))
recmat <- read.csv(file.path(formal, "fuel_reclassification_matrix_counts.csv"))
recind <- read.csv(file.path(formal, "fuel_reclassification_individual_deidentified.csv"))
sapdev <- read.csv(file.path(formal, "sap_deviations.csv"))
check_effect_rows <- function(z) {
  stopifnot(nrow(z) > 0, all(is.finite(z$estimate)), all(is.finite(z$lower)),
            all(is.finite(z$upper)), all(z$lower <= z$estimate & z$estimate <= z$upper),
            all(z$n > 0))
}
check_effect_rows(sens)
check_effect_rows(ext)
expected_ext <- c(
  ">=3 prior effective fuel waves", "First-last state definition",
  "Sex-stratified exploratory: Male", "Sex-stratified exploratory: Female",
  "Exploratory inverse selection weighted"
)
stopifnot(
  all(expected_ext %in% ext$analysis),
  nrow(inter) == 2L,
  all(is.finite(inter$interaction_estimate)),
  all(inter$p_interaction >= 0 & inter$p_interaction <= 1),
  nrow(selbal) > 0,
  all(is.finite(selbal$smd_selected_vs_excluded)),
  nrow(selwd) == 1L,
  selwd$ESS_trimmed > 0,
  sum(recmat$n) == 6544L,
  nrow(recind) == 6544L,
  nrow(sapdev) >= 1L
)

# De-identification checks: public analysis/reclassification/bootstrap tables must
# not contain the original person, household, community, stratum, or internal keys.
forbidden <- c("Idind", "hhid", "commid", "stratum", "Idind_key")
deid_names <- names(read.csv(file.path(project, "data", "analysis_dataset_deidentified.csv"), nrows = 1))
rec_names <- names(recind)
boot_names <- names(boot)
stopifnot(
  !any(forbidden %in% deid_names),
  !any(forbidden %in% rec_names),
  !any(forbidden %in% boot_names),
  "record_id" %in% deid_names,
  "record_id" %in% rec_names,
  !anyDuplicated(recind$record_id)
)

# R-only Nature export checks.
pdim <- dim(png::readPNG(paste0(fig, ".png")))
tdim <- dim(tiff::readTIFF(paste0(fig, ".tiff")))
stopifnot(all(pdim[1:2] == c(2362, 4322)), all(tdim[1:2] == c(2362, 4322)))
svg <- paste(readLines(paste0(fig, ".svg"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
stopifnot(
  grepl("<text", svg, fixed = TRUE),
  grepl("KDM biological age deviation", svg, fixed = TRUE),
  grepl("Fully adjusted (bootstrap)", svg, fixed = TRUE)
)

script <- paste(readLines(file.path(project, "scripts", "03_formal_analysis.R"), warn = FALSE), collapse = "\n")
stopifnot(grepl("wave <= 2009", script, fixed = TRUE), !grepl("wave == 2011", script, fixed = TRUE), !grepl("wave == 2015", script, fixed = TRUE))

# Text quality scan limited to user-facing Markdown/CSV/log/TXT files.
textfiles <- list.files(project, pattern = "\\.(md|csv|log|txt)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
badpat <- c("\uFFFD", "????", "閸?", "閺?", "閳?", "閴?")
bad <- character()
for (f in textfiles) {
  z <- paste(suppressWarnings(readLines(f, warn = FALSE, encoding = "UTF-8", skipNul = TRUE)), collapse = "\n")
  z <- suppressWarnings(iconv(z, from = "UTF-8", to = "UTF-8", sub = ""))
  if (any(vapply(badpat, function(p) grepl(p, z, fixed = TRUE), logical(1)))) bad <- c(bad, f)
}
stopifnot(length(bad) == 0)

manifest <- data.frame(
  file = sub(paste0("^", project, "/"), "", normalizePath(required, winslash = "/", mustWork = TRUE)),
  bytes = file.info(required)$size,
  md5 = unname(tools::md5sum(required))
)
write.csv(manifest, file.path(formal, "final_manifest_md5.csv"), row.names = FALSE)
cat("PASS: final C2 validation with community cluster bootstrap\n")
cat("analysis n=6544; primary model n=4514; KDM folds=10; bootstrap=1000/1000; warnings=0; failures=0\n")
cat("bootstrap percentile CIs reproducible; extended sensitivities PASS; de-identification PASS; Figure2=4322x2362\n")
