options(stringsAsFactors = FALSE, scipen = 999)

root_dir <- normalizePath(getwd(), mustWork = TRUE)
project_dir <- file.path(root_dir, "SDC1_AAD_HARRELL_BOOTCAL_20260728")
run_label <- "FINAL_B300_M50_NESTED_LOCKED"
source_dir <- file.path(
  project_dir, "04_results", "source_data", "discovery_nested_selection", run_label
)
stat_dir <- file.path(
  project_dir, "04_results", "statistical_results", "discovery_nested_selection", run_label
)
qc_dir <- file.path(
  project_dir, "07_quality_control", "discovery_nested_selection", run_label
)

draws <- read.csv(
  file.path(source_dir, "nested_mi_selection_bootstrap_draws.csv"),
  check.names = FALSE
)
reported <- read.csv(
  file.path(stat_dir, "nested_mi_sdc1_selection_stability.csv"),
  check.names = FALSE
)
tasks <- readRDS(file.path(source_dir, "bootstrap_sample_indices.rds"))
raw <- readRDS(file.path(
  root_dir,
  "SDC1_AAD_REANALYSIS_V6_1_REVIEWER_REFINED_20260715",
  "02_analysis_data", "proteomics_raw_deidentified.rds"
))

B <- nrow(draws)
retained <- sum(draws$sdc1_retained)
rank1_all <- sum(draws$sdc1_rank == 1L, na.rm = TRUE)
rank1_retained <- sum(draws$sdc1_rank == 1L, na.rm = TRUE)
retain_ci <- binom.test(retained, B)$conf.int
rank1_all_ci <- binom.test(rank1_all, B)$conf.int
rank1_retained_ci <- binom.test(rank1_retained, retained)$conf.int

checks <- data.frame(
  check = c(
    "draw_count_300",
    "all_imputation_counts_50",
    "retained_count_matches",
    "rank1_all_count_matches",
    "rank1_retained_count_matches",
    "retained_ci_matches",
    "rank1_all_ci_matches",
    "rank1_retained_ci_matches",
    "every_resample_preserves_stratum_sizes",
    "only_one_retained_non_rank1_draw",
    "retained_non_rank1_is_rank2"
  ),
  pass = c(
    B == 300L,
    all(draws$imputations_repeated == 50L),
    retained == reported$retained_count,
    rank1_all == reported$rank1_count_all,
    rank1_retained == reported$rank1_count_retained,
    max(abs(retain_ci - c(
      reported$retained_ci_lower, reported$retained_ci_upper
    ))) < 1e-12,
    max(abs(rank1_all_ci - c(
      reported$rank1_all_ci_lower, reported$rank1_all_ci_upper
    ))) < 1e-12,
    max(abs(rank1_retained_ci - c(
      reported$rank1_retained_ci_lower, reported$rank1_retained_ci_upper
    ))) < 1e-12,
    all(vapply(tasks, function(task) {
      identical(
        sort(as.integer(table(raw$sample_metadata$group[task$sampled]))),
        sort(as.integer(table(raw$sample_metadata$group)))
      )
    }, logical(1))),
    sum(draws$sdc1_retained & draws$sdc1_rank != 1L, na.rm = TRUE) == 1L,
    unique(draws$sdc1_rank[
      draws$sdc1_retained & draws$sdc1_rank != 1L
    ]) == 2L
  ),
  stringsAsFactors = FALSE
)

write.csv(
  checks,
  file.path(qc_dir, "independent_nested_selection_audit.csv"),
  row.names = FALSE
)
if (!all(checks$pass)) stop("Independent nested-selection audit failed")

cat("Independent nested-selection audit PASS\n")
cat("Retained:", retained, "/", B, "\n")
cat("Rank 1 among all:", rank1_all, "/", B, "\n")
cat("Rank 1 conditional on retention:", rank1_retained, "/", retained, "\n")

