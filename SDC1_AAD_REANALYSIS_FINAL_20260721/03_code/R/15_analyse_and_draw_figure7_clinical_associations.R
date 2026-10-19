options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
})

data_dir <- file.path(root, "02_analysis_data")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
figure_dir <- file.path(root, "04_results", "figures", "main")
panel_dir <- file.path(root, "04_results", "figures", "panels", "figure_7_clinical")
qa_dir <- file.path(root, "07_quality_control")
log_dir <- file.path(root, "04_results", "logs")
for (directory in c(result_dir, source_dir, figure_dir, panel_dir, qa_dir, log_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

read_final <- function(file) {
  read.csv(file.path(data_dir, file), check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

discovery <- read_final("discovery_clinical_deidentified.csv")
clinical <- read_final("clinical_evaluation_all.csv")
paired <- read_final("paired_sdc1_subset.csv")

group_levels <- c("Healthy control", "AAD", "AMI", "Angina", "PE")
for (object_name in c("discovery", "clinical", "paired")) {
  object <- get(object_name)
  object$group <- factor(object$group, levels = group_levels)
  assign(object_name, object)
}

stopifnot(
  nrow(discovery) == 307L,
  nrow(clinical) == 415L,
  nrow(paired) == 96L,
  length(unique(discovery$cluster_id)) == 307L,
  length(unique(clinical$cluster_id)) == 415L,
  length(unique(paired$cluster_id)) == 96L
)

numeric_variables <- c(
  "age_years", "female", "hypertension", "ev_sdc1", "total_sdc1", "shed_sdc1",
  "d_dimer_mg_l", "crp_mg_l", "ctnt_source_units", "bnp_source_units",
  "ckmb_source_units", "creatinine_umol_l", "wbc_10e9_l", "platelets_10e9_l",
  "onset_hours", "cardiogenic_shock", "heart_failure", "cardiac_death", "malperfusion"
)
for (object_name in c("discovery", "clinical", "paired")) {
  object <- get(object_name)
  for (variable in intersect(numeric_variables, names(object))) {
    object[[variable]] <- suppressWarnings(as.numeric(object[[variable]]))
  }
  assign(object_name, object)
}

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

symptomatic_specs <- list(
  list(
    row_id = "Discovery EV",
    row_label = "Discovery\nEV-associated",
    cohort = "Discovery cohort",
    assay = "EV-associated SDC1",
    data = discovery[discovery$group != "Healthy control", , drop = FALSE],
    sdc1 = "ev_sdc1"
  ),
  list(
    row_id = "Clinical total",
    row_label = "Clinical evaluation\nTotal plasma",
    cohort = "Clinical-evaluation cohort",
    assay = "Total-plasma SDC1",
    data = clinical[clinical$group != "Healthy control", , drop = FALSE],
    sdc1 = "total_sdc1"
  ),
  list(
    row_id = "Paired total",
    row_label = "Nested paired subset\nTotal plasma",
    cohort = "Nested paired subset",
    assay = "Total-plasma SDC1",
    data = paired[paired$group != "Healthy control", , drop = FALSE],
    sdc1 = "total_sdc1"
  ),
  list(
    row_id = "Paired EV",
    row_label = "Nested paired subset\nEV-associated",
    cohort = "Nested paired subset",
    assay = "EV-associated SDC1",
    data = paired[paired$group != "Healthy control", , drop = FALSE],
    sdc1 = "ev_sdc1"
  ),
  list(
    row_id = "Paired shed",
    row_label = "Nested paired subset\nShed",
    cohort = "Nested paired subset",
    assay = "Shed SDC1",
    data = paired[paired$group != "Healthy control", , drop = FALSE],
    sdc1 = "shed_sdc1"
  )
)

aad_specs <- list(
  list(
    row_id = "Discovery AAD",
    row_label = "Discovery AAD\nEV-associated",
    cohort = "Discovery cohort",
    assay = "EV-associated SDC1",
    data = discovery[discovery$group == "AAD", , drop = FALSE],
    sdc1 = "ev_sdc1"
  ),
  list(
    row_id = "Clinical AAD",
    row_label = "Clinical-evaluation AAD\nTotal plasma",
    cohort = "Clinical-evaluation cohort",
    assay = "Total-plasma SDC1",
    data = clinical[clinical$group == "AAD", , drop = FALSE],
    sdc1 = "total_sdc1"
  )
)

rank_association <- function(data, x, y, adjusters) {
  required <- c(x, y, adjusters)
  selected <- complete.cases(data[, required, drop = FALSE])
  working <- data[selected, required, drop = FALSE]
  working <- working[is.finite(working[[x]]) & is.finite(working[[y]]), , drop = FALSE]
  n <- nrow(working)
  if (n < 15L || length(unique(working[[x]])) < 4L || length(unique(working[[y]])) < 4L) {
    return(list(n = n, raw = NA_real_, adjusted = NA_real_, p_raw = NA_real_, p_adjusted = NA_real_, working = working))
  }

  rank_x <- rank(working[[x]], ties.method = "average")
  rank_y <- rank(working[[y]], ties.method = "average")
  raw <- cor(rank_x, rank_y)
  raw_df <- n - 2L
  raw_t <- raw * sqrt(raw_df / max(1e-12, 1 - raw^2))
  p_raw <- 2 * pt(abs(raw_t), df = raw_df, lower.tail = FALSE)

  design <- model.matrix(reformulate(adjusters), data = working)
  decomposition <- qr(design)
  residual_x <- qr.resid(decomposition, rank_x)
  residual_y <- qr.resid(decomposition, rank_y)
  adjusted <- cor(residual_x, residual_y)
  adjusted_df <- n - decomposition$rank - 1L
  adjusted_t <- adjusted * sqrt(adjusted_df / max(1e-12, 1 - adjusted^2))
  p_adjusted <- 2 * pt(abs(adjusted_t), df = adjusted_df, lower.tail = FALSE)

  list(
    n = n,
    raw = raw,
    adjusted = adjusted,
    p_raw = p_raw,
    p_adjusted = p_adjusted,
    working = working
  )
}

resample_indices <- function(data, stratify) {
  if (is.null(stratify)) return(sample.int(nrow(data), nrow(data), replace = TRUE))
  strata <- split(seq_len(nrow(data)), droplevels(data[[stratify]]))
  unlist(lapply(strata, function(indices) sample(indices, length(indices), replace = TRUE)), use.names = FALSE)
}

bootstrap_rank_association <- function(working, x, y, adjusters, stratify, repetitions, seed) {
  set.seed(seed)
  draws <- matrix(NA_real_, nrow = repetitions, ncol = 2L, dimnames = list(NULL, c("raw", "adjusted")))
  for (replicate_index in seq_len(repetitions)) {
    sampled <- working[resample_indices(working, stratify), , drop = FALSE]
    estimate <- rank_association(sampled, x, y, adjusters)
    draws[replicate_index, ] <- c(estimate$raw, estimate$adjusted)
  }
  as.data.frame(draws)
}

n_bootstrap_correlations <- 1000L
correlation_rows <- list()
correlation_draw_rows <- list()

run_correlation_family <- function(specifications, adjusters, stratify, family_name, seed_offset) {
  family_summary <- list()
  family_draws <- list()
  counter <- 0L
  for (specification in specifications) {
    for (marker_label in names(clinical_markers)) {
      counter <- counter + 1L
      marker <- clinical_markers[[marker_label]]
      estimate <- rank_association(specification$data, specification$sdc1, marker, adjusters)
      draws <- bootstrap_rank_association(
        estimate$working,
        specification$sdc1,
        marker,
        adjusters,
        stratify,
        n_bootstrap_correlations,
        seed = 20260721L + seed_offset + counter
      )
      raw_interval <- unname(quantile(draws$raw, c(0.025, 0.975), na.rm = TRUE))
      adjusted_interval <- unname(quantile(draws$adjusted, c(0.025, 0.975), na.rm = TRUE))
      family_summary[[counter]] <- data.frame(
        family = family_name,
        row_id = specification$row_id,
        row_label = specification$row_label,
        cohort = specification$cohort,
        assay = specification$assay,
        marker = marker_label,
        marker_variable = marker,
        n = estimate$n,
        raw_rho = estimate$raw,
        raw_lower_95 = raw_interval[1],
        raw_upper_95 = raw_interval[2],
        raw_p_value = estimate$p_raw,
        adjusted_rho = estimate$adjusted,
        adjusted_lower_95 = adjusted_interval[1],
        adjusted_upper_95 = adjusted_interval[2],
        adjusted_p_value = estimate$p_adjusted,
        stringsAsFactors = FALSE
      )
      family_draws[[counter]] <- data.frame(
        family = family_name,
        row_id = specification$row_id,
        marker = marker_label,
        bootstrap_replicate = seq_len(n_bootstrap_correlations),
        raw_rho = draws$raw,
        adjusted_rho = draws$adjusted,
        stringsAsFactors = FALSE
      )
    }
  }
  list(summary = do.call(rbind, family_summary), draws = do.call(rbind, family_draws))
}

symptomatic_correlations <- run_correlation_family(
  symptomatic_specs,
  adjusters = c("group", "age_years", "female", "hypertension"),
  stratify = "group",
  family_name = "Symptomatic participants",
  seed_offset = 1000L
)

aad_correlations <- run_correlation_family(
  aad_specs,
  adjusters = c("age_years", "female", "hypertension"),
  stratify = NULL,
  family_name = "AAD only",
  seed_offset = 100000L
)

correlation_summary <- rbind(symptomatic_correlations$summary, aad_correlations$summary)
correlation_draws <- rbind(symptomatic_correlations$draws, aad_correlations$draws)
correlation_summary$raw_bh_q <- ave(
  correlation_summary$raw_p_value,
  correlation_summary$family,
  correlation_summary$row_id,
  FUN = function(values) p.adjust(values, method = "BH")
)
correlation_summary$adjusted_bh_q <- ave(
  correlation_summary$adjusted_p_value,
  correlation_summary$family,
  correlation_summary$row_id,
  FUN = function(values) p.adjust(values, method = "BH")
)

build_marker_patient_source <- function(specifications, family_name, adjusters) {
  rows <- list()
  counter <- 0L
  for (specification in specifications) {
    for (marker_label in names(clinical_markers)) {
      counter <- counter + 1L
      marker <- clinical_markers[[marker_label]]
      data <- specification$data
      required <- c(specification$sdc1, marker, adjusters)
      rows[[counter]] <- data.frame(
        family = family_name,
        row_id = specification$row_id,
        cohort = specification$cohort,
        assay = specification$assay,
        row_uid = data$row_uid,
        cluster_id = data$cluster_id,
        group = as.character(data$group),
        age_years = data$age_years,
        female = data$female,
        hypertension = data$hypertension,
        sdc1_value = data[[specification$sdc1]],
        marker = marker_label,
        marker_value = data[[marker]],
        analysis_complete = complete.cases(data[, required, drop = FALSE]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

clinical_marker_patient_data <- rbind(
  build_marker_patient_source(
    symptomatic_specs,
    "Symptomatic participants",
    c("group", "age_years", "female", "hypertension")
  ),
  build_marker_patient_source(
    aad_specs,
    "AAD only",
    c("age_years", "female", "hypertension")
  )
)

onset_specs <- list(
  list(cohort = "Discovery cohort", assay = "EV-associated SDC1", data = discovery[discovery$group == "AAD", , drop = FALSE], sdc1 = "ev_sdc1"),
  list(cohort = "Clinical-evaluation cohort", assay = "Total-plasma SDC1", data = clinical[clinical$group == "AAD", , drop = FALSE], sdc1 = "total_sdc1")
)
n_bootstrap_onset <- 2000L
onset_summary_rows <- list()
onset_draw_rows <- list()
onset_patient_rows <- list()
for (index in seq_along(onset_specs)) {
  specification <- onset_specs[[index]]
  estimate <- rank_association(
    specification$data,
    specification$sdc1,
    "onset_hours",
    c("age_years", "female", "hypertension")
  )
  working <- estimate$working[estimate$working$onset_hours > 0 & estimate$working[[specification$sdc1]] > 0, , drop = FALSE]
  estimate <- rank_association(working, specification$sdc1, "onset_hours", c("age_years", "female", "hypertension"))
  draws <- bootstrap_rank_association(
    estimate$working,
    specification$sdc1,
    "onset_hours",
    c("age_years", "female", "hypertension"),
    NULL,
    n_bootstrap_onset,
    20261721L + index
  )
  interval <- unname(quantile(draws$adjusted, c(0.025, 0.975), na.rm = TRUE))
  onset_summary_rows[[index]] <- data.frame(
    cohort = specification$cohort,
    assay = specification$assay,
    n = estimate$n,
    partial_spearman_rho = estimate$adjusted,
    lower_95 = interval[1],
    upper_95 = interval[2],
    p_value = estimate$p_adjusted,
    stringsAsFactors = FALSE
  )
  onset_draw_rows[[index]] <- data.frame(
    cohort = specification$cohort,
    bootstrap_replicate = seq_len(n_bootstrap_onset),
    partial_spearman_rho = draws$adjusted,
    stringsAsFactors = FALSE
  )
  onset_patient_rows[[index]] <- data.frame(
    cohort = specification$cohort,
    assay = specification$assay,
    onset_hours = estimate$working$onset_hours,
    sdc1_pg_ml = estimate$working[[specification$sdc1]],
    log10_onset_hours = log10(estimate$working$onset_hours),
    log2_sdc1 = log2(estimate$working[[specification$sdc1]]),
    stringsAsFactors = FALSE
  )
}
onset_summary <- do.call(rbind, onset_summary_rows)
onset_summary$bh_q_value <- p.adjust(onset_summary$p_value, method = "BH")
onset_draws <- do.call(rbind, onset_draw_rows)
onset_patient_data <- do.call(rbind, onset_patient_rows)

hc3_covariance <- function(fit) {
  design <- model.matrix(fit)
  residual <- residuals(fit)
  leverage <- hatvalues(fit)
  bread <- solve(crossprod(design))
  bread %*% crossprod(design, design * (residual / (1 - leverage))^2) %*% bread
}

adjusted_outcome_difference <- function(data, sdc1, outcome, adjusters) {
  required <- c(sdc1, outcome, adjusters)
  working <- data[complete.cases(data[, required, drop = FALSE]), required, drop = FALSE]
  working$outcome_binary <- as.numeric(working[[outcome]])
  events <- sum(working$outcome_binary == 1)
  non_events <- sum(working$outcome_binary == 0)
  if (events < 2L || non_events < 2L) {
    return(c(n = nrow(working), events = events, estimate = NA, lower = NA, upper = NA, p = NA))
  }
  working$standardised_log2_sdc1 <- as.numeric(scale(log2(working[[sdc1]])))
  fit <- lm(
    as.formula(paste("standardised_log2_sdc1 ~ outcome_binary +", paste(adjusters, collapse = " + "))),
    data = working
  )
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

outcome_specifications <- list(
  list(cohort = "Discovery cohort", assay = "EV-associated SDC1", population = "Symptomatic", data = discovery[discovery$group != "Healthy control", , drop = FALSE], sdc1 = "ev_sdc1", outcome = "cardiogenic_shock", label = "Cardiogenic shock", adjusters = c("group", "age_years", "female", "hypertension")),
  list(cohort = "Clinical-evaluation cohort", assay = "Total-plasma SDC1", population = "Symptomatic", data = clinical[clinical$group != "Healthy control", , drop = FALSE], sdc1 = "total_sdc1", outcome = "cardiogenic_shock", label = "Cardiogenic shock", adjusters = c("group", "age_years", "female", "hypertension")),
  list(cohort = "Discovery cohort", assay = "EV-associated SDC1", population = "AAD", data = discovery[discovery$group == "AAD", , drop = FALSE], sdc1 = "ev_sdc1", outcome = "cardiac_death", label = "Cardiac death", adjusters = c("age_years", "female", "hypertension")),
  list(cohort = "Discovery cohort", assay = "EV-associated SDC1", population = "AAD", data = discovery[discovery$group == "AAD", , drop = FALSE], sdc1 = "ev_sdc1", outcome = "malperfusion", label = "Malperfusion", adjusters = c("age_years", "female", "hypertension")),
  list(cohort = "Discovery cohort", assay = "EV-associated SDC1", population = "AAD", data = discovery[discovery$group == "AAD", , drop = FALSE], sdc1 = "ev_sdc1", outcome = "heart_failure", label = "Heart failure", adjusters = c("age_years", "female", "hypertension")),
  list(cohort = "Clinical-evaluation cohort", assay = "Total-plasma SDC1", population = "AAD", data = clinical[clinical$group == "AAD", , drop = FALSE], sdc1 = "total_sdc1", outcome = "cardiac_death", label = "Cardiac death", adjusters = c("age_years", "female", "hypertension")),
  list(cohort = "Clinical-evaluation cohort", assay = "Total-plasma SDC1", population = "AAD", data = clinical[clinical$group == "AAD", , drop = FALSE], sdc1 = "total_sdc1", outcome = "malperfusion", label = "Malperfusion", adjusters = c("age_years", "female", "hypertension")),
  list(cohort = "Clinical-evaluation cohort", assay = "Total-plasma SDC1", population = "AAD", data = clinical[clinical$group == "AAD", , drop = FALSE], sdc1 = "total_sdc1", outcome = "heart_failure", label = "Heart failure", adjusters = c("age_years", "female", "hypertension"))
)

outcome_rows <- list()
for (index in seq_along(outcome_specifications)) {
  specification <- outcome_specifications[[index]]
  estimate <- adjusted_outcome_difference(
    specification$data,
    specification$sdc1,
    specification$outcome,
    specification$adjusters
  )
  outcome_rows[[index]] <- data.frame(
    cohort = specification$cohort,
    assay = specification$assay,
    population = specification$population,
    outcome = specification$label,
    n = estimate[["n"]],
    events = estimate[["events"]],
    adjusted_sd_difference = estimate[["estimate"]],
    lower_95 = estimate[["lower"]],
    upper_95 = estimate[["upper"]],
    p_value = estimate[["p"]],
    eligible_main_figure = estimate[["events"]] >= 10,
    stringsAsFactors = FALSE
  )
}
outcome_summary <- do.call(rbind, outcome_rows)
eligible_outcomes <- outcome_summary[outcome_summary$eligible_main_figure, , drop = FALSE]
eligible_outcomes$bh_q_value <- p.adjust(eligible_outcomes$p_value, method = "BH")
eligible_outcomes$display_label <- paste0(
  ifelse(eligible_outcomes$cohort == "Discovery cohort", "Discovery", "Clinical evaluation"),
  ": ", eligible_outcomes$outcome, "\n", eligible_outcomes$events, "/", eligible_outcomes$n, " events"
)

outcome_patient_rows <- list()
for (index in seq_along(outcome_specifications)) {
  specification <- outcome_specifications[[index]]
  data <- specification$data
  required <- c(specification$sdc1, specification$outcome, specification$adjusters)
  outcome_patient_rows[[index]] <- data.frame(
    cohort = specification$cohort,
    assay = specification$assay,
    population = specification$population,
    outcome = specification$label,
    row_uid = data$row_uid,
    cluster_id = data$cluster_id,
    group = as.character(data$group),
    age_years = data$age_years,
    female = data$female,
    hypertension = data$hypertension,
    sdc1_value = data[[specification$sdc1]],
    outcome_value = data[[specification$outcome]],
    analysis_complete = complete.cases(data[, required, drop = FALSE]),
    main_figure_eligible = outcome_summary$eligible_main_figure[index],
    stringsAsFactors = FALSE
  )
}
outcome_patient_data <- do.call(rbind, outcome_patient_rows)

write.csv(correlation_summary, file.path(result_dir, "fig7_clinical_correlation_summary.csv"), row.names = FALSE)
write.csv(correlation_draws, file.path(source_dir, "fig7_clinical_correlation_bootstrap_draws.csv"), row.names = FALSE)
write.csv(clinical_marker_patient_data, file.path(source_dir, "fig7_clinical_marker_patient_data.csv"), row.names = FALSE)
write.csv(onset_summary, file.path(result_dir, "fig7_aad_onset_associations.csv"), row.names = FALSE)
write.csv(onset_draws, file.path(source_dir, "fig7_aad_onset_bootstrap_draws.csv"), row.names = FALSE)
write.csv(onset_patient_data, file.path(source_dir, "fig7_aad_onset_patient_data.csv"), row.names = FALSE)
write.csv(outcome_summary, file.path(result_dir, "fig7_clinical_outcome_associations.csv"), row.names = FALSE)
write.csv(eligible_outcomes, file.path(source_dir, "fig7_clinical_outcome_main_panel.csv"), row.names = FALSE)
write.csv(outcome_patient_data, file.path(source_dir, "fig7_clinical_outcome_patient_data.csv"), row.names = FALSE)

theme_pub <- function() {
  theme_classic(base_family = "Arial", base_size = 6.5) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#202124"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#202124"),
      axis.text = element_text(size = 5.9, colour = "#30343B"),
      axis.title = element_text(size = 6.3, colour = "#202124"),
      plot.title = element_text(size = 7.3, face = "bold", colour = "#202124", margin = margin(b = 1.5)),
      plot.subtitle = element_text(size = 5.9, colour = "#70777D", margin = margin(b = 3)),
      strip.text = element_text(size = 6.2, face = "bold"),
      legend.title = element_text(size = 5.9),
      legend.text = element_text(size = 5.7),
      plot.margin = margin(4, 5, 4, 5)
    )
}

marker_levels <- names(clinical_markers)
symptomatic_row_levels <- vapply(symptomatic_specs, `[[`, character(1), "row_label")
aad_row_levels <- vapply(aad_specs, `[[`, character(1), "row_label")

make_heatmap_data <- function(summary, family, estimate_name, q_name, row_levels) {
  data <- summary[summary$family == family, , drop = FALSE]
  data$rho <- data[[estimate_name]]
  data$q_value <- data[[q_name]]
  data$marker <- factor(data$marker, levels = marker_levels)
  data$row_label <- factor(data$row_label, levels = rev(row_levels))
  data$label <- sprintf("%.2f%s", data$rho, ifelse(data$q_value < 0.05, "*", ""))
  data
}

raw_heatmap <- make_heatmap_data(
  correlation_summary,
  "Symptomatic participants",
  "raw_rho",
  "raw_bh_q",
  symptomatic_row_levels
)
adjusted_heatmap <- make_heatmap_data(
  correlation_summary,
  "Symptomatic participants",
  "adjusted_rho",
  "adjusted_bh_q",
  symptomatic_row_levels
)
aad_heatmap <- make_heatmap_data(
  correlation_summary,
  "AAD only",
  "adjusted_rho",
  "adjusted_bh_q",
  aad_row_levels
)

heatmap_scale <- scale_fill_gradient2(
  low = "#2A9D8F",
  mid = "#F7F7F5",
  high = "#BF3F4C",
  midpoint = 0,
  limits = c(-0.60, 0.60),
  oob = scales::squish,
  name = "rho"
)

make_heatmap <- function(data, title, subtitle, show_legend) {
  ggplot(data, aes(marker, row_label, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.55) +
    geom_text(aes(label = label), family = "Arial", size = 2.05, colour = "#202124") +
    heatmap_scale +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal(base_family = "Arial", base_size = 6.5) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 32, hjust = 1, size = 5.6, colour = "#30343B"),
      axis.text.y = element_text(size = 5.8, colour = "#30343B"),
      plot.title = element_text(size = 7.3, face = "bold", colour = "#202124", margin = margin(b = 1.5)),
      plot.subtitle = element_text(size = 5.9, colour = "#70777D", margin = margin(b = 3)),
      legend.position = if (show_legend) "bottom" else "none",
      legend.key.width = grid::unit(18, "mm"),
      legend.key.height = grid::unit(2.2, "mm"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

p_a <- make_heatmap(
  raw_heatmap,
  "Unadjusted clinical-marker correlations",
  "Symptomatic participants; * BH q<0.05 within each assay-cohort row",
  FALSE
)
p_b <- make_heatmap(
  adjusted_heatmap,
  "Adjusted clinical-marker correlations",
  "Partial Spearman rho; adjusted for group and fixed covariates",
  TRUE
)
p_c <- make_heatmap(
  aad_heatmap,
  "AAD-only clinical-marker correlations",
  "Partial Spearman rho adjusted for age, sex and hypertension",
  FALSE
)

onset_palette <- c("Discovery cohort" = "#2A9D8F", "Clinical-evaluation cohort" = "#BF3F4C")
make_onset_plot <- function(cohort_name, title, y_title) {
  plot_data <- onset_patient_data[onset_patient_data$cohort == cohort_name, , drop = FALSE]
  estimate <- onset_summary[onset_summary$cohort == cohort_name, , drop = FALSE]
  annotation <- sprintf(
    "partial rho = %.2f\n95%% CI %.2f to %.2f\nBH q = %.3f",
    estimate$partial_spearman_rho,
    estimate$lower_95,
    estimate$upper_95,
    estimate$bh_q_value
  )
  ggplot(plot_data, aes(log10_onset_hours, log2_sdc1)) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 0.85, colour = onset_palette[[cohort_name]], fill = scales::alpha(onset_palette[[cohort_name]], 0.16), linewidth = 0.65) +
    geom_point(colour = onset_palette[[cohort_name]], alpha = 0.52, size = 1.05) +
    annotate("label", x = -Inf, y = Inf, label = annotation, hjust = -0.05, vjust = 1.08, family = "Arial", size = 2.05, linewidth = 0, fill = scales::alpha("white", 0.85)) +
    scale_x_continuous(
      breaks = log10(c(0.5, 1, 3, 6, 12, 24, 72, 168, 336)),
      labels = c("0.5", "1", "3", "6", "12", "24", "72", "168", "336")
    ) +
    labs(
      x = "Time from symptom onset (hours; log scale)",
      y = y_title,
      title = title,
      subtitle = paste0("AAD only; n=", estimate$n, "; continuous time without binning")
    ) +
    theme_pub()
}

p_d <- make_onset_plot("Discovery cohort", "Discovery AAD time course", "EV-associated SDC1 (log2 pg/mL)")
p_e <- make_onset_plot("Clinical-evaluation cohort", "Clinical-evaluation AAD time course", "Total-plasma SDC1 (log2 pg/mL)")

eligible_outcomes$cohort_short <- ifelse(eligible_outcomes$cohort == "Discovery cohort", "Discovery", "Clinical evaluation")
eligible_outcomes$display_label <- factor(eligible_outcomes$display_label, levels = rev(eligible_outcomes$display_label))
outcome_palette <- c("Discovery" = "#2A9D8F", "Clinical evaluation" = "#BF3F4C")
p_f <- ggplot(eligible_outcomes, aes(adjusted_sd_difference, display_label, colour = cohort_short)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = "#8E969C") +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), orientation = "y", width = 0, linewidth = 0.6) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = outcome_palette, name = NULL) +
  labs(
    x = "Adjusted standardised log2 SDC1 difference (event minus no event)",
    y = NULL,
    title = "No clear association with adverse clinical phenotypes",
    subtitle = "Only combinations with at least 10 events are displayed; 95% HC3 robust intervals"
  ) +
  theme_pub() +
  theme(legend.position = "bottom")

figure7 <- (p_a | p_b) / p_c / (p_d | p_e) / p_f +
  plot_layout(heights = c(1.05, 0.62, 1.02, 0.72)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 8, family = "Arial", margin = margin(r = 3, b = 2)))
  )

save_figure <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = 300)
  print(plot)
  dev.off()
}

figure_stem <- file.path(figure_dir, "Figure_7_clinical_associations_across_cohorts")
save_figure(figure7, figure_stem, width_mm = 183, height_mm = 230, dpi = 600)

panel_list <- list(A = p_a, B = p_b, C = p_c, D = p_d, E = p_e, F = p_f)
panel_dimensions <- list(A = c(90, 92), B = c(90, 92), C = c(183, 60), D = c(90, 82), E = c(90, 82), F = c(183, 68))
for (panel_name in names(panel_list)) {
  dimensions <- panel_dimensions[[panel_name]]
  save_figure(
    panel_list[[panel_name]],
    file.path(panel_dir, paste0("Figure_7", panel_name)),
    width_mm = dimensions[1],
    height_mm = dimensions[2],
    dpi = 600
  )
}

qa_checks <- c(
  discovery_records_307 = nrow(discovery) == 307L,
  clinical_records_415 = nrow(clinical) == 415L,
  paired_records_96 = nrow(paired) == 96L,
  symptomatic_correlation_rows_40 = nrow(symptomatic_correlations$summary) == 40L,
  aad_correlation_rows_16 = nrow(aad_correlations$summary) == 16L,
  correlation_bootstrap_draws_56000 = nrow(correlation_draws) == 56000L,
  all_correlation_estimates_finite = all(is.finite(correlation_summary$raw_rho)) && all(is.finite(correlation_summary$adjusted_rho)),
  onset_rows_2 = nrow(onset_summary) == 2L,
  onset_bootstrap_draws_4000 = nrow(onset_draws) == 4000L,
  eligible_outcome_rows_4 = nrow(eligible_outcomes) == 4L,
  ineligible_outcomes_retained = sum(!outcome_summary$eligible_main_figure) == 4L,
  marker_patient_source_present = file.exists(file.path(source_dir, "fig7_clinical_marker_patient_data.csv")),
  outcome_patient_source_present = file.exists(file.path(source_dir, "fig7_clinical_outcome_patient_data.csv")),
  main_exports_present = all(file.exists(paste0(figure_stem, c(".pdf", ".svg", ".tiff", ".png")))),
  standalone_exports_present = all(vapply(LETTERS[1:6], function(panel_name) {
    all(file.exists(file.path(panel_dir, paste0("Figure_7", panel_name, c(".pdf", ".svg", ".tiff", ".png")))))
  }, logical(1)))
)
qa_table <- data.frame(check = names(qa_checks), pass = as.logical(qa_checks))
write.csv(qa_table, file.path(qa_dir, "FINAL_FIGURE7_CLINICAL_ASSOCIATIONS_QA.csv"), row.names = FALSE)

summary_lines <- c(
  "Figure 7 clinical associations across cohorts",
  paste0("Overall computational QA: ", if (all(qa_checks)) "PASS" else "FAIL"),
  "Discovery: n=307, EV-associated SDC1",
  "Clinical evaluation: n=415, total-plasma SDC1",
  "Nested paired subset: n=96; not counted as independent",
  "Eight prespecified clinical markers; no data-driven marker selection",
  "No pooling across cohorts because assay types differ and patient overlap is unresolved",
  "Sparse outcomes with fewer than 10 events are retained in results but excluded from the main panel",
  paste0("Main figure: ", figure_stem)
)
writeLines(summary_lines, file.path(log_dir, "fig7_clinical_associations_summary.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "fig7_clinical_associations_sessionInfo.txt"), useBytes = TRUE)
cat(paste(summary_lines, collapse = "\n"), "\n")
if (!all(qa_checks)) quit(status = 1)
