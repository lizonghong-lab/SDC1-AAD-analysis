options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "02_analysis_data")
qa_dir <- file.path(root, "07_quality_control")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

read_final <- function(file) {
  data <- read.csv(file.path(data_dir, file), check.names = FALSE, fileEncoding = "UTF-8-BOM")
  for (name in names(data)) {
    if (is.character(data[[name]]) && all(is.na(suppressWarnings(as.numeric(data[[name]]))) == is.na(data[[name]]))) {
      data[[name]] <- suppressWarnings(as.numeric(data[[name]]))
    }
  }
  data
}

discovery <- read_final("discovery_clinical_deidentified.csv")
clinical <- read_final("clinical_evaluation_all.csv")

clinical_markers <- c(
  "D-dimer" = "d_dimer_mg_l",
  "CRP" = "crp_mg_l",
  "cTnT" = "ctnt_source_units",
  "BNP" = "bnp_source_units",
  "CK-MB" = "ckmb_source_units",
  "Creatinine" = "creatinine_umol_l",
  "WBC" = "wbc_10e9_l",
  "Platelets" = "platelets_10e9_l"
)

partial_spearman <- function(data, x, y, adjusters) {
  selected <- complete.cases(data[, c(x, y, adjusters), drop = FALSE])
  working <- data[selected, c(x, y, adjusters), drop = FALSE]
  if (nrow(working) < 15L || length(unique(working[[x]])) < 4L || length(unique(working[[y]])) < 4L) {
    return(c(n = nrow(working), rho = NA_real_, p = NA_real_))
  }
  working$rank_x <- rank(working[[x]], ties.method = "average")
  working$rank_y <- rank(working[[y]], ties.method = "average")
  rhs <- paste(adjusters, collapse = " + ")
  residual_x <- residuals(lm(as.formula(paste("rank_x ~", rhs)), data = working))
  residual_y <- residuals(lm(as.formula(paste("rank_y ~", rhs)), data = working))
  test <- cor.test(residual_x, residual_y, method = "pearson")
  c(n = nrow(working), rho = unname(test$estimate), p = test$p.value)
}

correlation_rows <- list()
cohort_specs <- list(
  discovery = list(data = discovery, sdc1 = "ev_sdc1", sdc1_label = "EV-associated SDC1"),
  clinical_evaluation = list(data = clinical, sdc1 = "total_sdc1", sdc1_label = "Total-plasma SDC1")
)

for (cohort_name in names(cohort_specs)) {
  specification <- cohort_specs[[cohort_name]]
  data <- specification$data
  data$group <- factor(data$group)
  for (population in c("All participants", "AAD only")) {
    analysis_data <- if (population == "AAD only") data[data$group == "AAD", , drop = FALSE] else data
    adjusters <- if (population == "AAD only") c("age_years", "female", "hypertension") else c("group", "age_years", "female", "hypertension")
    for (marker_label in names(clinical_markers)) {
      marker <- clinical_markers[[marker_label]]
      estimate <- partial_spearman(analysis_data, specification$sdc1, marker, adjusters)
      correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(
        cohort = cohort_name,
        sdc1_measure = specification$sdc1_label,
        population = population,
        marker = marker_label,
        marker_variable = marker,
        n = estimate[["n"]],
        partial_spearman_rho = estimate[["rho"]],
        p_value = estimate[["p"]],
        stringsAsFactors = FALSE
      )
    }
  }
}
correlations <- do.call(rbind, correlation_rows)
correlations$bh_q_within_population <- ave(correlations$p_value, correlations$cohort, correlations$population, FUN = function(x) p.adjust(x, method = "BH"))

outcomes <- c(
  "Cardiogenic shock" = "cardiogenic_shock",
  "Heart failure" = "heart_failure",
  "Cardiac death" = "cardiac_death",
  "Malperfusion" = "malperfusion"
)

event_rows <- list()
for (cohort_name in names(cohort_specs)) {
  specification <- cohort_specs[[cohort_name]]
  aad <- specification$data[specification$data$group == "AAD", , drop = FALSE]
  for (outcome_label in names(outcomes)) {
    outcome <- outcomes[[outcome_label]]
    available <- !is.na(aad[[outcome]]) & as.character(aad[[outcome]]) != ""
    values <- suppressWarnings(as.numeric(aad[[outcome]][available]))
    event_rows[[length(event_rows) + 1L]] <- data.frame(
      cohort = cohort_name,
      sdc1_measure = specification$sdc1_label,
      outcome = outcome_label,
      n_available = sum(available),
      events = sum(values == 1, na.rm = TRUE),
      non_events = sum(values == 0, na.rm = TRUE),
      other_values = sum(!values %in% c(0, 1), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}
events <- do.call(rbind, event_rows)

symptomatic_event_rows <- list()
for (cohort_name in names(cohort_specs)) {
  specification <- cohort_specs[[cohort_name]]
  symptomatic <- specification$data[specification$data$group != "Healthy control", , drop = FALSE]
  for (outcome_label in names(outcomes)) {
    outcome <- outcomes[[outcome_label]]
    available <- !is.na(symptomatic[[outcome]]) & as.character(symptomatic[[outcome]]) != ""
    values <- suppressWarnings(as.numeric(symptomatic[[outcome]][available]))
    symptomatic_event_rows[[length(symptomatic_event_rows) + 1L]] <- data.frame(
      cohort = cohort_name,
      sdc1_measure = specification$sdc1_label,
      population = "Symptomatic participants",
      outcome = outcome_label,
      n_available = sum(available),
      events = sum(values == 1, na.rm = TRUE),
      non_events = sum(values == 0, na.rm = TRUE),
      other_values = sum(!values %in% c(0, 1), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}
symptomatic_events <- do.call(rbind, symptomatic_event_rows)

hc3_covariance <- function(fit) {
  design <- model.matrix(fit)
  residual <- residuals(fit)
  leverage <- hatvalues(fit)
  bread <- solve(crossprod(design))
  bread %*% crossprod(design, design * (residual / (1 - leverage))^2) %*% bread
}

adjusted_difference <- function(data, sdc1, outcome, adjusters) {
  selected <- complete.cases(data[, c(sdc1, outcome, adjusters), drop = FALSE])
  working <- data[selected, c(sdc1, outcome, adjusters), drop = FALSE]
  working$outcome_binary <- as.numeric(working[[outcome]])
  events <- sum(working$outcome_binary == 1)
  non_events <- sum(working$outcome_binary == 0)
  if (events < 2L || non_events < 2L) {
    return(c(n = nrow(working), events = events, estimate = NA, lower = NA, upper = NA, p = NA))
  }
  working$log2_sdc1_z <- as.numeric(scale(log2(working[[sdc1]])))
  formula <- as.formula(paste("log2_sdc1_z ~ outcome_binary +", paste(adjusters, collapse = " + ")))
  fit <- lm(formula, data = working)
  covariance <- hc3_covariance(fit)
  estimate <- coef(fit)[["outcome_binary"]]
  standard_error <- sqrt(covariance["outcome_binary", "outcome_binary"])
  c(
    n = nrow(working),
    events = events,
    estimate = estimate,
    lower = estimate - qnorm(0.975) * standard_error,
    upper = estimate + qnorm(0.975) * standard_error,
    p = 2 * pnorm(abs(estimate / standard_error), lower.tail = FALSE)
  )
}

outcome_association_rows <- list()
for (cohort_name in names(cohort_specs)) {
  specification <- cohort_specs[[cohort_name]]
  full_data <- specification$data
  symptomatic <- full_data[full_data$group != "Healthy control", , drop = FALSE]
  aad <- full_data[full_data$group == "AAD", , drop = FALSE]

  shock_result <- adjusted_difference(
    symptomatic,
    specification$sdc1,
    "cardiogenic_shock",
    c("group", "age_years", "female", "hypertension")
  )
  outcome_association_rows[[length(outcome_association_rows) + 1L]] <- data.frame(
    cohort = cohort_name,
    sdc1_measure = specification$sdc1_label,
    population = "Symptomatic participants",
    outcome = "Cardiogenic shock",
    n = shock_result[["n"]],
    events = shock_result[["events"]],
    adjusted_sdc1_sd_difference = shock_result[["estimate"]],
    lower_95 = shock_result[["lower"]],
    upper_95 = shock_result[["upper"]],
    p_value = shock_result[["p"]],
    stringsAsFactors = FALSE
  )

  for (outcome_label in c("Cardiac death", "Malperfusion", "Heart failure")) {
    outcome_variable <- outcomes[[outcome_label]]
    result <- adjusted_difference(
      aad,
      specification$sdc1,
      outcome_variable,
      c("age_years", "female", "hypertension")
    )
    outcome_association_rows[[length(outcome_association_rows) + 1L]] <- data.frame(
      cohort = cohort_name,
      sdc1_measure = specification$sdc1_label,
      population = "AAD only",
      outcome = outcome_label,
      n = result[["n"]],
      events = result[["events"]],
      adjusted_sdc1_sd_difference = result[["estimate"]],
      lower_95 = result[["lower"]],
      upper_95 = result[["upper"]],
      p_value = result[["p"]],
      stringsAsFactors = FALSE
    )
  }
}
outcome_associations <- do.call(rbind, outcome_association_rows)
outcome_associations$eligible_primary_figure <- outcome_associations$events >= 10L

time_rows <- list()
for (cohort_name in names(cohort_specs)) {
  specification <- cohort_specs[[cohort_name]]
  aad <- specification$data[specification$data$group == "AAD", , drop = FALSE]
  available <- complete.cases(aad[, c(specification$sdc1, "onset_hours", "age_years", "female", "hypertension")])
  working <- aad[available, , drop = FALSE]
  estimate <- partial_spearman(working, specification$sdc1, "onset_hours", c("age_years", "female", "hypertension"))
  time_rows[[cohort_name]] <- data.frame(
    cohort = cohort_name,
    sdc1_measure = specification$sdc1_label,
    n = estimate[["n"]],
    partial_spearman_rho = estimate[["rho"]],
    p_value = estimate[["p"]],
    onset_min = min(working$onset_hours, na.rm = TRUE),
    onset_median = median(working$onset_hours, na.rm = TRUE),
    onset_max = max(working$onset_hours, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
time_summary <- do.call(rbind, time_rows)

write.csv(correlations, file.path(qa_dir, "FIGURE7_EXPLORATORY_PARTIAL_CORRELATIONS.csv"), row.names = FALSE)
write.csv(events, file.path(qa_dir, "FIGURE7_AAD_EVENT_COUNTS.csv"), row.names = FALSE)
write.csv(symptomatic_events, file.path(qa_dir, "FIGURE7_SYMPTOMATIC_EVENT_COUNTS.csv"), row.names = FALSE)
write.csv(outcome_associations, file.path(qa_dir, "FIGURE7_EXPLORATORY_OUTCOME_ASSOCIATIONS.csv"), row.names = FALSE)
write.csv(time_summary, file.path(qa_dir, "FIGURE7_AAD_ONSET_ASSOCIATION_PROFILE.csv"), row.names = FALSE)

cat("\nPartial correlations:\n")
print(correlations, row.names = FALSE)
cat("\nAAD event counts:\n")
print(events, row.names = FALSE)
cat("\nSymptomatic-cohort event counts:\n")
print(symptomatic_events, row.names = FALSE)
cat("\nExploratory outcome associations:\n")
print(outcome_associations, row.names = FALSE)
cat("\nAAD onset-time profile:\n")
print(time_summary, row.names = FALSE)
