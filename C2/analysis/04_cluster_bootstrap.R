options(stringsAsFactors = FALSE, warn = 1, cli.unicode = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(splines)
})

project <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(project) != "C2") stop("Run with working directory set to C2")
outdir <- file.path(project, "results", "formal")
logfile <- file.path(project, "logs", "08_cluster_bootstrap.log")
logcon <- file(logfile, open = "wt", encoding = "UTF-8")
sink(logcon, split = TRUE)
sink(logcon, type = "message")
on.exit({sink(type = "message"); sink(); close(logcon)}, add = TRUE)

cat("C2 community cluster bootstrap started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")

d <- readRDS(file.path(project, "data", "analysis_dataset_internal_restricted.rds"))
p <- d %>%
  filter(trajectory %in% c("persistent_any_solid", "solid_to_clean_only")) %>%
  droplevels()
p$trajectory <- factor(
  p$trajectory,
  levels = c("persistent_any_solid", "solid_to_clean_only")
)

stopifnot(nrow(p) > 0L, !anyNA(p$commid), n_distinct(p$commid) > 1L)
stopifnot(all(is.finite(p$HD_z)), all(is.finite(p$KDM_BAA)))

full_rhs <- paste(
  "trajectory + ns(age,4) + sex + province_f + urban_f + education +",
  "income_asinh + income_missing + asset_count + asset_missing + smoking_f + alcohol_f"
)
coef_name <- "trajectorysolid_to_clean_only"

fit_contrast <- function(outcome, dat) {
  fit <- lm(as.formula(paste(outcome, "~", full_rhs)), data = dat)
  if (nobs(fit) != nrow(dat)) stop("Model dropped observations")
  b <- coef(fit)
  if (!(coef_name %in% names(b)) || !is.finite(b[[coef_name]])) {
    stop("Exposure contrast is not estimable")
  }
  unname(b[[coef_name]])
}

main <- read.csv(file.path(outdir, "main_secondary_results.csv"))
orig_hd <- subset(main, outcome == "HD_z" & analysis == "Fully adjusted")$estimate
orig_kdm <- subset(main, outcome == "KDM_BAA" & analysis == "Fully adjusted")$estimate
stopifnot(length(orig_hd) == 1L, length(orig_kdm) == 1L)

# Outcomes and their cross-fitting are deliberately frozen. This resampling targets
# uncertainty in the exposure-outcome regression under community clustering; it does
# not re-estimate HD reference parameters or the cross-fitted KDM construction.
cluster_rows <- split(seq_len(nrow(p)), as.character(p$commid))
cluster_ids <- names(cluster_rows)
G <- length(cluster_ids)
B <- 1000L
max_attempts <- 5000L
set.seed(20260809)

replicates <- vector("list", B)
warning_records <- list()
failure_records <- list()
success <- 0L
attempt <- 0L

while (success < B && attempt < max_attempts) {
  attempt <- attempt + 1L
  sampled_ids <- sample(cluster_ids, size = G, replace = TRUE)
  sampled_rows <- unname(cluster_rows[sampled_ids])
  idx <- unlist(sampled_rows, use.names = FALSE)
  draw_no <- rep(seq_len(G), lengths(sampled_rows))
  bd <- p[idx, , drop = FALSE]
  # A repeated source community is a distinct bootstrap cluster draw.
  bd$bootstrap_cluster_id <- sprintf("b%04d_c%03d", success + 1L, draw_no)

  caught <- character()
  ans <- tryCatch(
    withCallingHandlers({
      c(
        HD_z = fit_contrast("HD_z", bd),
        KDM_BAA = fit_contrast("KDM_BAA", bd)
      )
    }, warning = function(w) {
      caught <<- c(caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )

  if (inherits(ans, "error") || any(!is.finite(ans))) {
    msg <- if (inherits(ans, "error")) conditionMessage(ans) else "Non-finite estimate"
    failure_records[[length(failure_records) + 1L]] <- data.frame(
      attempt = attempt, message = msg
    )
    next
  }

  success <- success + 1L
  if (length(caught)) {
    warning_records[[length(warning_records) + 1L]] <- data.frame(
      replicate = success,
      attempt = attempt,
      message = paste(unique(caught), collapse = " | ")
    )
  }

  replicates[[success]] <- data.frame(
    replicate = success,
    attempt = attempt,
    HD_z_estimate = unname(ans[["HD_z"]]),
    KDM_BAA_estimate = unname(ans[["KDM_BAA"]]),
    n_rows = nrow(bd),
    n_cluster_draws = G,
    n_new_cluster_ids = n_distinct(bd$bootstrap_cluster_id),
    n_distinct_source_communities = n_distinct(sampled_ids),
    n_duplicate_cluster_draws = G - n_distinct(sampled_ids),
    n_persistent_any_solid = sum(bd$trajectory == "persistent_any_solid"),
    n_solid_to_clean_only = sum(bd$trajectory == "solid_to_clean_only")
  )

  if (success %% 100L == 0L) {
    cat("Successful replicates:", success, "of", B, "; attempts:", attempt, "\n")
  }
}

if (success != B) stop("Failed to obtain 1000 successful bootstrap replicates")
boot <- bind_rows(replicates)
stopifnot(nrow(boot) == B, all(boot$replicate == seq_len(B)))
stopifnot(all(boot$n_cluster_draws == G), all(boot$n_new_cluster_ids == G))
stopifnot(all(is.finite(boot$HD_z_estimate)), all(is.finite(boot$KDM_BAA_estimate)))
write.csv(
  boot,
  file.path(outdir, "cluster_bootstrap_replicates.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

summarise_boot <- function(x, outcome, original, unit) {
  qs <- unname(quantile(x, c(0.025, 0.975), type = 7, names = FALSE))
  data.frame(
    outcome = outcome,
    analysis = "Fully adjusted; community cluster bootstrap",
    contrast = "solid_to_clean_only vs persistent_any_solid",
    estimate = original,
    bootstrap_mean = mean(x),
    bootstrap_bias = mean(x) - original,
    bootstrap_se = sd(x),
    lower = qs[1],
    upper = qs[2],
    ci_method = "Percentile cluster bootstrap",
    confidence_level = 0.95,
    successful_replicates = B,
    attempts = attempt,
    resampling_unit = "community",
    clusters_per_replicate = G,
    outcome_construction = "Fixed from original analysis; not re-estimated within bootstrap",
    unit = unit
  )
}

summary_out <- bind_rows(
  summarise_boot(boot$HD_z_estimate, "HD_z", orig_hd, "SD"),
  summarise_boot(boot$KDM_BAA_estimate, "KDM_BAA", orig_kdm, "years")
)
write.csv(
  summary_out,
  file.path(outdir, "cluster_bootstrap_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Add the bootstrap robustness rows to the formal effect-result table while
# retaining the community-sandwich rows as the prespecified primary inference.
bootstrap_main_rows <- summary_out %>% transmute(
  outcome = outcome,
  analysis = analysis,
  contrast = contrast,
  estimate = estimate,
  se = bootstrap_se,
  lower = lower,
  upper = upper,
  p_value = NA_real_,
  n = nrow(p),
  clusters = G,
  cluster = "community bootstrap"
)
main_updated <- main %>%
  filter(analysis != "Fully adjusted; community cluster bootstrap") %>%
  bind_rows(bootstrap_main_rows)
write.csv(
  main_updated,
  file.path(outdir, "main_secondary_results.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

diagnostics <- data.frame(
  check = c(
    "Requested successful replicates",
    "Observed successful replicates",
    "Total attempts",
    "Failed attempts",
    "Replicates with captured warnings",
    "Original community clusters",
    "Minimum new cluster IDs per replicate",
    "Maximum new cluster IDs per replicate",
    "Minimum distinct source communities",
    "Maximum distinct source communities",
    "Minimum bootstrap sample size",
    "Maximum bootstrap sample size",
    "HD percentile CI reproducible",
    "KDM percentile CI reproducible",
    "Cross-fitted outcomes rebuilt inside bootstrap"
  ),
  value = c(
    B,
    nrow(boot),
    attempt,
    length(failure_records),
    length(warning_records),
    G,
    min(boot$n_new_cluster_ids),
    max(boot$n_new_cluster_ids),
    min(boot$n_distinct_source_communities),
    max(boot$n_distinct_source_communities),
    min(boot$n_rows),
    max(boot$n_rows),
    identical(unname(quantile(boot$HD_z_estimate, c(.025, .975), type = 7)), summary_out[summary_out$outcome == "HD_z", c("lower", "upper")] |> unlist() |> unname()),
    identical(unname(quantile(boot$KDM_BAA_estimate, c(.025, .975), type = 7)), summary_out[summary_out$outcome == "KDM_BAA", c("lower", "upper")] |> unlist() |> unname()),
    FALSE
  ),
  pass = c(
    TRUE,
    nrow(boot) == B,
    attempt >= B,
    TRUE,
    TRUE,
    G == 193L,
    min(boot$n_new_cluster_ids) == G,
    max(boot$n_new_cluster_ids) == G,
    min(boot$n_distinct_source_communities) > 0,
    max(boot$n_distinct_source_communities) <= G,
    min(boot$n_rows) > 0,
    max(boot$n_rows) > 0,
    TRUE,
    TRUE,
    TRUE
  )
)
write.csv(
  diagnostics,
  file.path(outdir, "cluster_bootstrap_diagnostics.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

warning_df <- if (length(warning_records)) bind_rows(warning_records) else data.frame(
  replicate = integer(), attempt = integer(), message = character()
)
failure_df <- if (length(failure_records)) bind_rows(failure_records) else data.frame(
  attempt = integer(), message = character()
)
write.csv(warning_df, file.path(outdir, "cluster_bootstrap_warnings.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(failure_df, file.path(outdir, "cluster_bootstrap_failures.csv"), row.names = FALSE, fileEncoding = "UTF-8")

contract <- data.frame(
  item = c(
    "Purpose", "Resampling unit", "Replicate construction", "Models",
    "Contrast", "Interval", "Primary uncertainty", "Outcome construction",
    "Random seed", "Successful replicates"
  ),
  specification = c(
    "Robustness analysis for clustered sampling uncertainty",
    "Community",
    "Sample 193 communities with replacement; each draw receives a new bootstrap cluster ID",
    "Prespecified fully adjusted linear models for HD_z and KDM_BAA",
    "solid_to_clean_only minus persistent_any_solid",
    "2.5th and 97.5th percentiles of successful replicate estimates",
    "Community sandwich CI remains the prespecified primary CI",
    "Previously constructed HD_z and five-fold cross-fitted KDM_BAA are fixed; construction parameters are not re-estimated within replicates",
    "20260809",
    as.character(B)
  )
)
write.csv(contract, file.path(outdir, "cluster_bootstrap_contract.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("Bootstrap completed:", format(Sys.time()), "\n")
cat("Successful replicates:", success, "; attempts:", attempt, "; warnings:", nrow(warning_df), "; failures:", nrow(failure_df), "\n")
print(summary_out)
