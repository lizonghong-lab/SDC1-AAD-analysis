options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
source_dir <- file.path(root, "04_results", "source_data")
result_dir <- file.path(root, "04_results", "statistical_results")
qa_dir <- file.path(root, "07_quality_control")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

patients <- read.csv(
  file.path(source_dir, "fig5_standalone_biomarker_roc_source_data.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stored_auc <- read.csv(
  file.path(result_dir, "fig5_standalone_biomarker_roc_auc.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stored_delta <- read.csv(
  file.path(result_dir, "fig5_standalone_biomarker_roc_paired_difference.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
draws <- read.csv(
  file.path(source_dir, "fig5_standalone_biomarker_roc_bootstrap_draws.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pairwise_auc <- function(y, marker) {
  case_values <- marker[y == 1]
  control_values <- marker[y == 0]
  comparisons <- outer(case_values, control_values, FUN = "-")
  wins <- sum(comparisons > 0)
  ties <- sum(comparisons == 0)
  losses <- sum(comparisons < 0)
  pairs <- length(case_values) * length(control_values)
  score <- wins + 0.5 * ties
  list(
    auc = score / pairs,
    wins = wins,
    ties = ties,
    losses = losses,
    score = score,
    pairs = pairs,
    n_case = length(case_values),
    n_control = length(control_values)
  )
}

populations <- list(
  "All symptomatic comparators" = c("AAD", "AMI", "Angina", "PE"),
  "AAD versus PE" = c("AAD", "PE")
)
markers <- c(D_dimer = "d_dimer_mg_l", SDC1 = "sdc1_source_value")

audit_rows <- list()
for (population_name in names(populations)) {
  subset_data <- patients[patients$group %in% populations[[population_name]], , drop = FALSE]
  y <- as.integer(subset_data$group == "AAD")
  for (biomarker_name in names(markers)) {
    marker_column <- markers[[biomarker_name]]
    result <- pairwise_auc(y, subset_data[[marker_column]])
    stored <- stored_auc[
      stored_auc$population == population_name & stored_auc$biomarker == biomarker_name,
      ,
      drop = FALSE
    ]
    marker_draws <- draws[
      draws$population == population_name & draws$biomarker == biomarker_name,
      "auc"
    ]
    interval <- as.numeric(quantile(marker_draws, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
    audit_rows[[length(audit_rows) + 1]] <- data.frame(
      population = population_name,
      biomarker = biomarker_name,
      n_case = result$n_case,
      n_control = result$n_control,
      pair_count = result$pairs,
      correctly_ordered_pairs = result$wins,
      tied_pairs = result$ties,
      incorrectly_ordered_pairs = result$losses,
      concordance_score = result$score,
      independently_recomputed_auc = result$auc,
      stored_auc = stored$estimate,
      auc_absolute_difference = abs(result$auc - stored$estimate),
      recomputed_lower_95 = interval[1],
      stored_lower_95 = stored$lower_95,
      recomputed_upper_95 = interval[2],
      stored_upper_95 = stored$upper_95,
      stringsAsFactors = FALSE
    )
  }
}
audit <- do.call(rbind, audit_rows)

delta_rows <- list()
for (population_name in names(populations)) {
  d_auc <- audit$independently_recomputed_auc[
    audit$population == population_name & audit$biomarker == "D_dimer"
  ]
  s_auc <- audit$independently_recomputed_auc[
    audit$population == population_name & audit$biomarker == "SDC1"
  ]
  wide_draws <- reshape(
    draws[draws$population == population_name, c("bootstrap_replicate", "biomarker", "auc")],
    idvar = "bootstrap_replicate",
    timevar = "biomarker",
    direction = "wide"
  )
  delta_draws <- wide_draws$auc.SDC1 - wide_draws$auc.D_dimer
  interval <- as.numeric(quantile(delta_draws, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
  stored <- stored_delta[stored_delta$population == population_name, , drop = FALSE]
  pair_count <- unique(audit$pair_count[audit$population == population_name])
  delta_rows[[length(delta_rows) + 1]] <- data.frame(
    population = population_name,
    independently_recomputed_delta_auc = s_auc - d_auc,
    stored_delta_auc = stored$difference,
    pairwise_concordance_score_difference = (s_auc - d_auc) * pair_count,
    recomputed_lower_95 = interval[1],
    stored_lower_95 = stored$lower_95,
    recomputed_upper_95 = interval[2],
    stored_upper_95 = stored$upper_95,
    stringsAsFactors = FALSE
  )
}
delta_audit <- do.call(rbind, delta_rows)

group_counts <- table(factor(patients$group, levels = c("AAD", "AMI", "Angina", "PE")))
missing_sdc1 <- sum(!is.finite(patients$sdc1_source_value))
missing_d_dimer <- sum(!is.finite(patients$d_dimer_mg_l))
numeric_tolerance <- 1e-12
pass <-
  nrow(patients) == 329 &&
  identical(as.integer(group_counts), c(99L, 102L, 85L, 43L)) &&
  missing_sdc1 == 0 &&
  missing_d_dimer == 0 &&
  all(audit$auc_absolute_difference < numeric_tolerance) &&
  all(abs(audit$recomputed_lower_95 - audit$stored_lower_95) < numeric_tolerance) &&
  all(abs(audit$recomputed_upper_95 - audit$stored_upper_95) < numeric_tolerance) &&
  all(abs(delta_audit$independently_recomputed_delta_auc - delta_audit$stored_delta_auc) < numeric_tolerance) &&
  all(abs(delta_audit$recomputed_lower_95 - delta_audit$stored_lower_95) < numeric_tolerance) &&
  all(abs(delta_audit$recomputed_upper_95 - delta_audit$stored_upper_95) < numeric_tolerance)

write.csv(audit, file.path(qa_dir, "FIGURE5_AUC_INDEPENDENT_AUDIT.csv"), row.names = FALSE, na = "")
write.csv(delta_audit, file.path(qa_dir, "FIGURE5_DELTA_AUC_INDEPENDENT_AUDIT.csv"), row.names = FALSE, na = "")

fmt <- function(x) sprintf("%.6f", x)
report <- c(
  "# Figure 5 AUROC independent audit",
  "",
  paste0("- Overall status: **", if (pass) "PASS" else "FAIL", "**"),
  "- Method: direct case-control pair enumeration from the 329-row Figure 5 patient-level Source Data; the stored ROC result table was not used to calculate AUROC.",
  paste0("- Group counts: AAD ", group_counts[["AAD"]], ", AMI ", group_counts[["AMI"]], ", Angina ", group_counts[["Angina"]], ", PE ", group_counts[["PE"]], "."),
  paste0("- Missing values: SDC1 ", missing_sdc1, ", D-dimer ", missing_d_dimer, "."),
  "",
  "| Population | Biomarker | Pair count | Concordance score | AUROC | 95% bootstrap CI |",
  "|---|---|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(audit))) {
  report <- c(
    report,
    paste0(
      "| ", audit$population[i], " | ", audit$biomarker[i], " | ",
      audit$pair_count[i], " | ", audit$concordance_score[i], " | ",
      fmt(audit$independently_recomputed_auc[i]), " | ",
      fmt(audit$recomputed_lower_95[i]), "-", fmt(audit$recomputed_upper_95[i]), " |"
    )
  )
}
report <- c(
  report,
  "",
  "| Population | Delta AUROC (SDC1-D-dimer) | Pairwise concordance-score difference | 95% paired bootstrap CI |",
  "|---|---:|---:|---:|"
)
for (i in seq_len(nrow(delta_audit))) {
  report <- c(
    report,
    paste0(
      "| ", delta_audit$population[i], " | ",
      fmt(delta_audit$independently_recomputed_delta_auc[i]), " | ",
      sprintf("%.1f", delta_audit$pairwise_concordance_score_difference[i]), " | ",
      fmt(delta_audit$recomputed_lower_95[i]), "-", fmt(delta_audit$recomputed_upper_95[i]), " |"
    )
  )
}
report <- c(
  report,
  "",
  "The all-comparator AUROCs are not mathematically identical. They differ by 0.000220, equivalent to five additional concordant AAD/non-AAD pair scores for SDC1 across 22,770 possible pairs; both round to 0.902 at three decimals. In AAD versus PE, SDC1 gains 708 concordant pair scores across 4,257 pairs, producing a Delta AUROC of 0.166314."
)
writeLines(report, file.path(qa_dir, "FIGURE5_AUC_INDEPENDENT_AUDIT.md"), useBytes = TRUE)

cat(paste(report, collapse = "\n"), "\n")
if (!pass) quit(status = 1)
