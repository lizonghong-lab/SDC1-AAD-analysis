options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
data_file <- file.path(root, "02_analysis_data", "paired_sdc1_subset.csv")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
qa_dir <- file.path(root, "07_quality_control")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(data_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
data$group <- factor(data$group, levels = c("Healthy control", "AAD", "AMI", "Angina", "PE"))
for (column in c("age_years", "female", "hypertension", "total_sdc1", "ev_sdc1", "shed_sdc1")) {
  data[[column]] <- as.numeric(data[[column]])
}

saved_summary <- read.csv(file.path(result_dir, "fig3_paired_compartment_group_summary.csv"), check.names = FALSE)
saved_effects <- read.csv(file.path(result_dir, "fig3_paired_compartment_adjusted_effects.csv"), check.names = FALSE)
saved_correlations <- read.csv(file.path(result_dir, "fig3_paired_compartment_correlations.csv"), check.names = FALSE)
saved_draws <- read.csv(file.path(source_dir, "fig3_paired_compartment_correlation_bootstrap_draws.csv"), check.names = FALSE)

assays <- c("Total plasma" = "total_sdc1", "EV-associated" = "ev_sdc1", "Shed" = "shed_sdc1")
expected_counts <- c(`Healthy control` = 19L, AAD = 21L, AMI = 17L, Angina = 20L, PE = 19L)

recomputed_medians <- do.call(rbind, lapply(names(assays), function(assay) {
  do.call(rbind, lapply(levels(data$group), function(group_name) {
    values <- data[[assays[[assay]]]][data$group == group_name]
    data.frame(assay = assay, group = group_name, n = length(values), median = median(values))
  }))
}))
summary_merged <- merge(
  recomputed_medians,
  saved_summary[, c("assay", "group", "n", "median")],
  by = c("assay", "group"), suffixes = c("_audit", "_saved"), sort = FALSE
)

hc3_covariance <- function(fit) {
  x <- model.matrix(fit)
  e <- residuals(fit)
  h <- hatvalues(fit)
  bread <- solve(crossprod(x))
  bread %*% crossprod(x, x * (e / (1 - h))^2) %*% bread
}

comparators <- c("Healthy control", "AMI", "Angina", "PE")
audit_effects <- list()
for (assay in names(assays)) {
  model_data <- data.frame(
    outcome = log2(data[[assays[[assay]]]]),
    group = data$group,
    age_years = data$age_years,
    female = data$female,
    hypertension = data$hypertension
  )
  fit <- lm(outcome ~ group + age_years + female + hypertension, data = model_data)
  covariance <- hc3_covariance(fit)
  beta <- coef(fit)
  for (comparator in comparators) {
    contrast <- setNames(rep(0, length(beta)), names(beta))
    contrast["groupAAD"] <- 1
    if (comparator != "Healthy control") contrast[paste0("group", comparator)] <- -1
    estimate <- sum(contrast * beta)
    se <- sqrt(drop(t(contrast) %*% covariance %*% contrast))
    audit_effects[[length(audit_effects) + 1L]] <- data.frame(
      assay = assay,
      comparator = comparator,
      estimate_audit = estimate,
      lower_audit = estimate - qnorm(0.975) * se,
      upper_audit = estimate + qnorm(0.975) * se
    )
  }
}
audit_effects <- do.call(rbind, audit_effects)
effect_merged <- merge(
  audit_effects,
  saved_effects[, c("assay", "comparator", "adjusted_log2_fold_change", "lower_95", "upper_95")],
  by = c("assay", "comparator"), sort = FALSE
)

pair_specifications <- list(
  "Total plasma vs EV-associated" = c("total_sdc1", "ev_sdc1"),
  "Total plasma vs shed" = c("total_sdc1", "shed_sdc1"),
  "EV-associated vs shed" = c("ev_sdc1", "shed_sdc1")
)
audit_correlations <- do.call(rbind, lapply(names(pair_specifications), function(pair_name) {
  variables <- pair_specifications[[pair_name]]
  data.frame(
    pair = pair_name,
    rho_audit = suppressWarnings(cor(log2(data[[variables[1]]]), log2(data[[variables[2]]]), method = "spearman"))
  )
}))
correlation_merged <- merge(
  audit_correlations,
  saved_correlations[, c("pair", "spearman_rho", "lower_95", "upper_95")],
  by = "pair", sort = FALSE
)

draw_intervals <- do.call(rbind, lapply(split(saved_draws$spearman_rho, saved_draws$pair), function(values) {
  interval <- unname(quantile(values, c(0.025, 0.975), na.rm = TRUE))
  data.frame(lower_draws = interval[1], upper_draws = interval[2])
}))
draw_intervals$pair <- rownames(draw_intervals)
rownames(draw_intervals) <- NULL
correlation_interval_merged <- merge(correlation_merged, draw_intervals, by = "pair", sort = FALSE)

tolerance <- 1e-10
checks <- c(
  records_96 = nrow(data) == 96L,
  unique_patients_96 = length(unique(data$cluster_id)) == 96L,
  group_counts_match = identical(as.integer(table(data$group)), unname(expected_counts)),
  three_assays_complete = all(complete.cases(data[, unname(assays)])),
  summary_rows_match = nrow(summary_merged) == 15L,
  group_ns_match = all(summary_merged$n_audit == summary_merged$n_saved),
  group_medians_match = max(abs(summary_merged$median_audit - summary_merged$median_saved)) < tolerance,
  adjusted_effect_rows_match = nrow(effect_merged) == 12L,
  adjusted_estimates_match = max(abs(effect_merged$estimate_audit - effect_merged$adjusted_log2_fold_change)) < tolerance,
  adjusted_lower_limits_match = max(abs(effect_merged$lower_audit - effect_merged$lower_95)) < tolerance,
  adjusted_upper_limits_match = max(abs(effect_merged$upper_audit - effect_merged$upper_95)) < tolerance,
  correlation_rows_match = nrow(correlation_merged) == 3L,
  spearman_estimates_match = max(abs(correlation_merged$rho_audit - correlation_merged$spearman_rho)) < tolerance,
  bootstrap_draws_6000 = nrow(saved_draws) == 6000L,
  bootstrap_lower_limits_match = max(abs(correlation_interval_merged$lower_draws - correlation_interval_merged$lower_95)) < tolerance,
  bootstrap_upper_limits_match = max(abs(correlation_interval_merged$upper_draws - correlation_interval_merged$upper_95)) < tolerance
)

qa <- data.frame(check = names(checks), pass = as.logical(checks))
write.csv(qa, file.path(qa_dir, "FINAL_FIGURE3_PAIRED_INDEPENDENT_AUDIT.csv"), row.names = FALSE)
write.csv(effect_merged, file.path(qa_dir, "FINAL_FIGURE3_PAIRED_EFFECT_AUDIT.csv"), row.names = FALSE)
write.csv(correlation_interval_merged, file.path(qa_dir, "FINAL_FIGURE3_PAIRED_CORRELATION_AUDIT.csv"), row.names = FALSE)

cat("Figure 3 independent numerical audit:", if (all(checks)) "PASS" else "FAIL", "\n")
cat("Checks passed:", sum(checks), "of", length(checks), "\n")
if (!all(checks)) quit(status = 1)
