options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "02_analysis_data")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
figure_dir <- file.path(root, "04_results", "figures", "main")
panel_dir <- file.path(root, "04_results", "figures", "panels", "figure_7_clinical")
qa_dir <- file.path(root, "07_quality_control")

read_final <- function(file) {
  data <- read.csv(file.path(data_dir, file), check.names = FALSE, fileEncoding = "UTF-8-BOM")
  data$group <- factor(data$group, levels = c("Healthy control", "AAD", "AMI", "Angina", "PE"))
  data
}
discovery <- read_final("discovery_clinical_deidentified.csv")
clinical <- read_final("clinical_evaluation_all.csv")
paired <- read_final("paired_sdc1_subset.csv")

saved_correlations <- read.csv(file.path(result_dir, "fig7_clinical_correlation_summary.csv"), check.names = FALSE)
saved_draws <- read.csv(file.path(source_dir, "fig7_clinical_correlation_bootstrap_draws.csv"), check.names = FALSE)
saved_onset <- read.csv(file.path(result_dir, "fig7_aad_onset_associations.csv"), check.names = FALSE)
saved_onset_draws <- read.csv(file.path(source_dir, "fig7_aad_onset_bootstrap_draws.csv"), check.names = FALSE)
saved_outcomes <- read.csv(file.path(result_dir, "fig7_clinical_outcome_associations.csv"), check.names = FALSE)

marker_map <- c(
  "D-dimer" = "d_dimer_mg_l", "CRP" = "crp_mg_l", "cTnT" = "ctnt_source_units",
  "BNP" = "bnp_source_units", "CK-MB" = "ckmb_source_units",
  "Creatinine" = "creatinine_umol_l", "WBC" = "wbc_10e9_l", "Platelets" = "platelets_10e9_l"
)

get_analysis_data <- function(row_id, family) {
  if (row_id == "Discovery EV") return(list(data = discovery[discovery$group != "Healthy control", , drop = FALSE], sdc1 = "ev_sdc1"))
  if (row_id == "Clinical total") return(list(data = clinical[clinical$group != "Healthy control", , drop = FALSE], sdc1 = "total_sdc1"))
  if (row_id == "Paired total") return(list(data = paired[paired$group != "Healthy control", , drop = FALSE], sdc1 = "total_sdc1"))
  if (row_id == "Paired EV") return(list(data = paired[paired$group != "Healthy control", , drop = FALSE], sdc1 = "ev_sdc1"))
  if (row_id == "Paired shed") return(list(data = paired[paired$group != "Healthy control", , drop = FALSE], sdc1 = "shed_sdc1"))
  if (row_id == "Discovery AAD") return(list(data = discovery[discovery$group == "AAD", , drop = FALSE], sdc1 = "ev_sdc1"))
  if (row_id == "Clinical AAD") return(list(data = clinical[clinical$group == "AAD", , drop = FALSE], sdc1 = "total_sdc1"))
  stop("Unknown row_id: ", row_id)
}

audit_correlation <- function(data, x, y, adjusters) {
  working <- data[complete.cases(data[, c(x, y, adjusters), drop = FALSE]), c(x, y, adjusters), drop = FALSE]
  rank_x <- rank(working[[x]], ties.method = "average")
  rank_y <- rank(working[[y]], ties.method = "average")
  raw <- cor(rank_x, rank_y)
  working$rank_x <- rank_x
  working$rank_y <- rank_y
  rhs <- paste(adjusters, collapse = " + ")
  residual_x <- residuals(lm(as.formula(paste("rank_x ~", rhs)), data = working))
  residual_y <- residuals(lm(as.formula(paste("rank_y ~", rhs)), data = working))
  c(n = nrow(working), raw = raw, adjusted = cor(residual_x, residual_y))
}

audit_rows <- vector("list", nrow(saved_correlations))
for (index in seq_len(nrow(saved_correlations))) {
  row <- saved_correlations[index, ]
  specification <- get_analysis_data(row$row_id, row$family)
  adjusters <- if (row$family == "AAD only") c("age_years", "female", "hypertension") else c("group", "age_years", "female", "hypertension")
  estimate <- audit_correlation(specification$data, specification$sdc1, marker_map[[row$marker]], adjusters)
  audit_rows[[index]] <- data.frame(
    family = row$family,
    row_id = row$row_id,
    marker = row$marker,
    n_audit = estimate[["n"]],
    raw_audit = estimate[["raw"]],
    adjusted_audit = estimate[["adjusted"]]
  )
}
audit_correlations <- do.call(rbind, audit_rows)
correlation_merged <- merge(
  saved_correlations,
  audit_correlations,
  by = c("family", "row_id", "marker"),
  sort = FALSE
)

draw_keys <- unique(saved_draws[, c("family", "row_id", "marker")])
draw_interval_rows <- vector("list", nrow(draw_keys))
for (index in seq_len(nrow(draw_keys))) {
  key <- draw_keys[index, ]
  selected <- saved_draws$family == key$family & saved_draws$row_id == key$row_id & saved_draws$marker == key$marker
  raw_interval <- unname(quantile(saved_draws$raw_rho[selected], c(0.025, 0.975), na.rm = TRUE))
  adjusted_interval <- unname(quantile(saved_draws$adjusted_rho[selected], c(0.025, 0.975), na.rm = TRUE))
  draw_interval_rows[[index]] <- data.frame(
    family = key$family,
    row_id = key$row_id,
    marker = key$marker,
    raw_lower_draws = raw_interval[1],
    raw_upper_draws = raw_interval[2],
    adjusted_lower_draws = adjusted_interval[1],
    adjusted_upper_draws = adjusted_interval[2]
  )
}
draw_intervals <- do.call(rbind, draw_interval_rows)
correlation_merged <- merge(correlation_merged, draw_intervals, by = c("family", "row_id", "marker"), sort = FALSE)

onset_audit_rows <- list()
for (index in seq_len(nrow(saved_onset))) {
  cohort <- saved_onset$cohort[index]
  if (cohort == "Discovery cohort") {
    data <- discovery[discovery$group == "AAD" & discovery$onset_hours > 0, , drop = FALSE]
    sdc1 <- "ev_sdc1"
  } else {
    data <- clinical[clinical$group == "AAD" & clinical$onset_hours > 0, , drop = FALSE]
    sdc1 <- "total_sdc1"
  }
  estimate <- audit_correlation(data, sdc1, "onset_hours", c("age_years", "female", "hypertension"))
  selected_draws <- saved_onset_draws$cohort == cohort
  interval <- unname(quantile(saved_onset_draws$partial_spearman_rho[selected_draws], c(0.025, 0.975), na.rm = TRUE))
  onset_audit_rows[[index]] <- data.frame(
    cohort = cohort,
    n_audit = estimate[["n"]],
    rho_audit = estimate[["adjusted"]],
    lower_draws = interval[1],
    upper_draws = interval[2]
  )
}
onset_audit <- do.call(rbind, onset_audit_rows)
onset_merged <- merge(saved_onset, onset_audit, by = "cohort", sort = FALSE)

hc3_covariance <- function(fit) {
  x <- model.matrix(fit)
  e <- residuals(fit)
  h <- hatvalues(fit)
  bread <- solve(crossprod(x))
  bread %*% crossprod(x, x * (e / (1 - h))^2) %*% bread
}

outcome_audit <- function(data, sdc1, outcome, adjusters) {
  working <- data[complete.cases(data[, c(sdc1, outcome, adjusters), drop = FALSE]), c(sdc1, outcome, adjusters), drop = FALSE]
  working$outcome_binary <- as.numeric(working[[outcome]])
  events <- sum(working$outcome_binary == 1)
  if (events < 2L) return(c(n = nrow(working), events = events, estimate = NA, lower = NA, upper = NA))
  working$z <- as.numeric(scale(log2(working[[sdc1]])))
  fit <- lm(as.formula(paste("z ~ outcome_binary +", paste(adjusters, collapse = " + "))), data = working)
  covariance <- hc3_covariance(fit)
  estimate <- coef(fit)[["outcome_binary"]]
  se <- sqrt(covariance["outcome_binary", "outcome_binary"])
  c(n = nrow(working), events = events, estimate = estimate, lower = estimate - qnorm(0.975) * se, upper = estimate + qnorm(0.975) * se)
}

outcome_map <- list(
  "Discovery cohort|Symptomatic|Cardiogenic shock" = list(data = discovery[discovery$group != "Healthy control", , drop = FALSE], sdc1 = "ev_sdc1", outcome = "cardiogenic_shock", adjusters = c("group", "age_years", "female", "hypertension")),
  "Clinical-evaluation cohort|Symptomatic|Cardiogenic shock" = list(data = clinical[clinical$group != "Healthy control", , drop = FALSE], sdc1 = "total_sdc1", outcome = "cardiogenic_shock", adjusters = c("group", "age_years", "female", "hypertension")),
  "Discovery cohort|AAD|Cardiac death" = list(data = discovery[discovery$group == "AAD", , drop = FALSE], sdc1 = "ev_sdc1", outcome = "cardiac_death", adjusters = c("age_years", "female", "hypertension")),
  "Discovery cohort|AAD|Malperfusion" = list(data = discovery[discovery$group == "AAD", , drop = FALSE], sdc1 = "ev_sdc1", outcome = "malperfusion", adjusters = c("age_years", "female", "hypertension")),
  "Discovery cohort|AAD|Heart failure" = list(data = discovery[discovery$group == "AAD", , drop = FALSE], sdc1 = "ev_sdc1", outcome = "heart_failure", adjusters = c("age_years", "female", "hypertension")),
  "Clinical-evaluation cohort|AAD|Cardiac death" = list(data = clinical[clinical$group == "AAD", , drop = FALSE], sdc1 = "total_sdc1", outcome = "cardiac_death", adjusters = c("age_years", "female", "hypertension")),
  "Clinical-evaluation cohort|AAD|Malperfusion" = list(data = clinical[clinical$group == "AAD", , drop = FALSE], sdc1 = "total_sdc1", outcome = "malperfusion", adjusters = c("age_years", "female", "hypertension")),
  "Clinical-evaluation cohort|AAD|Heart failure" = list(data = clinical[clinical$group == "AAD", , drop = FALSE], sdc1 = "total_sdc1", outcome = "heart_failure", adjusters = c("age_years", "female", "hypertension"))
)

outcome_audit_rows <- vector("list", nrow(saved_outcomes))
for (index in seq_len(nrow(saved_outcomes))) {
  row <- saved_outcomes[index, ]
  key <- paste(row$cohort, row$population, row$outcome, sep = "|")
  specification <- outcome_map[[key]]
  estimate <- outcome_audit(specification$data, specification$sdc1, specification$outcome, specification$adjusters)
  outcome_audit_rows[[index]] <- data.frame(
    cohort = row$cohort,
    population = row$population,
    outcome = row$outcome,
    n_audit = estimate[["n"]],
    events_audit = estimate[["events"]],
    estimate_audit = estimate[["estimate"]],
    lower_audit = estimate[["lower"]],
    upper_audit = estimate[["upper"]]
  )
}
outcome_audit_table <- do.call(rbind, outcome_audit_rows)
outcome_merged <- merge(saved_outcomes, outcome_audit_table, by = c("cohort", "population", "outcome"), sort = FALSE)

tolerance <- 1e-10
finite_outcomes <- is.finite(outcome_merged$adjusted_sd_difference)
checks <- c(
  cohort_rows_match = nrow(discovery) == 307L && nrow(clinical) == 415L && nrow(paired) == 96L,
  correlation_rows_56 = nrow(correlation_merged) == 56L,
  correlation_n_match = all(correlation_merged$n == correlation_merged$n_audit),
  raw_correlations_match = max(abs(correlation_merged$raw_rho - correlation_merged$raw_audit)) < tolerance,
  adjusted_correlations_match = max(abs(correlation_merged$adjusted_rho - correlation_merged$adjusted_audit)) < tolerance,
  correlation_draws_56000 = nrow(saved_draws) == 56000L,
  raw_interval_lower_match = max(abs(correlation_merged$raw_lower_95 - correlation_merged$raw_lower_draws)) < tolerance,
  raw_interval_upper_match = max(abs(correlation_merged$raw_upper_95 - correlation_merged$raw_upper_draws)) < tolerance,
  adjusted_interval_lower_match = max(abs(correlation_merged$adjusted_lower_95 - correlation_merged$adjusted_lower_draws)) < tolerance,
  adjusted_interval_upper_match = max(abs(correlation_merged$adjusted_upper_95 - correlation_merged$adjusted_upper_draws)) < tolerance,
  onset_rows_2 = nrow(onset_merged) == 2L,
  onset_n_match = all(onset_merged$n == onset_merged$n_audit),
  onset_correlations_match = max(abs(onset_merged$partial_spearman_rho - onset_merged$rho_audit)) < tolerance,
  onset_interval_lower_match = max(abs(onset_merged$lower_95 - onset_merged$lower_draws)) < tolerance,
  onset_interval_upper_match = max(abs(onset_merged$upper_95 - onset_merged$upper_draws)) < tolerance,
  outcome_rows_8 = nrow(outcome_merged) == 8L,
  outcome_n_match = all(outcome_merged$n == outcome_merged$n_audit),
  outcome_event_counts_match = all(outcome_merged$events == outcome_merged$events_audit),
  outcome_estimates_match = max(abs(outcome_merged$adjusted_sd_difference[finite_outcomes] - outcome_merged$estimate_audit[finite_outcomes])) < tolerance,
  outcome_lower_limits_match = max(abs(outcome_merged$lower_95[finite_outcomes] - outcome_merged$lower_audit[finite_outcomes])) < tolerance,
  outcome_upper_limits_match = max(abs(outcome_merged$upper_95[finite_outcomes] - outcome_merged$upper_audit[finite_outcomes])) < tolerance,
  main_exports_present = all(file.exists(file.path(figure_dir, paste0("Figure_7_clinical_associations_across_cohorts", c(".pdf", ".svg", ".tiff", ".png"))))),
  panel_exports_present = all(vapply(LETTERS[1:6], function(panel) all(file.exists(file.path(panel_dir, paste0("Figure_7", panel, c(".pdf", ".svg", ".tiff", ".png"))))), logical(1)))
)

qa <- data.frame(check = names(checks), pass = as.logical(checks))
write.csv(qa, file.path(qa_dir, "FINAL_FIGURE7_CLINICAL_ASSOCIATIONS_INDEPENDENT_AUDIT.csv"), row.names = FALSE)
write.csv(correlation_merged, file.path(qa_dir, "FINAL_FIGURE7_CORRELATION_AUDIT.csv"), row.names = FALSE)
write.csv(onset_merged, file.path(qa_dir, "FINAL_FIGURE7_ONSET_AUDIT.csv"), row.names = FALSE)
write.csv(outcome_merged, file.path(qa_dir, "FINAL_FIGURE7_OUTCOME_AUDIT.csv"), row.names = FALSE)

cat("Figure 7 independent numerical audit:", if (all(checks)) "PASS" else "FAIL", "\n")
cat("Checks passed:", sum(checks), "of", length(checks), "\n")
if (!all(checks)) quit(status = 1)
