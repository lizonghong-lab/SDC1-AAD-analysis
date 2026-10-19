options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
})

data_file <- file.path(root, "02_analysis_data", "paired_sdc1_subset.csv")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
figure_dir <- file.path(root, "04_results", "figures", "main")
panel_dir <- file.path(root, "04_results", "figures", "panels", "figure_3_paired")
qa_dir <- file.path(root, "07_quality_control")
log_dir <- file.path(root, "04_results", "logs")
for (directory in c(result_dir, source_dir, figure_dir, panel_dir, qa_dir, log_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

set.seed(20260721)

data <- read.csv(data_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
required <- c(
  "row_uid", "cluster_id", "group", "age_years", "female", "hypertension",
  "total_sdc1", "ev_sdc1", "shed_sdc1"
)
missing_columns <- setdiff(required, names(data))
if (length(missing_columns)) stop("Missing required columns: ", paste(missing_columns, collapse = ", "))

group_order <- c("Healthy control", "AAD", "AMI", "Angina", "PE")
expected_counts <- c("Healthy control" = 19L, "AAD" = 21L, "AMI" = 17L, "Angina" = 20L, "PE" = 19L)
data$group <- factor(data$group, levels = group_order)

numeric_columns <- c("age_years", "female", "hypertension", "total_sdc1", "ev_sdc1", "shed_sdc1")
data[numeric_columns] <- lapply(data[numeric_columns], as.numeric)

if (nrow(data) != 96L) stop("Expected 96 paired records; found ", nrow(data))
if (length(unique(data$cluster_id)) != 96L) stop("Expected 96 unique patient keys")
if (!identical(as.integer(table(data$group)), unname(expected_counts))) stop("Unexpected paired-subset group counts")
if (any(!complete.cases(data[, required]))) stop("Required paired analysis fields are incomplete")
if (any(data[, c("total_sdc1", "ev_sdc1", "shed_sdc1")] <= 0)) stop("SDC1 measurements must be positive")

data$log2_total_sdc1 <- log2(data$total_sdc1)
data$log2_ev_sdc1 <- log2(data$ev_sdc1)
data$log2_shed_sdc1 <- log2(data$shed_sdc1)

group_palette <- c(
  "Healthy control" = "#4B8F8A",
  "AAD" = "#BF3F4C",
  "AMI" = "#3F78A8",
  "Angina" = "#D78B2B",
  "PE" = "#7B5AA6"
)
assay_palette <- c(
  "Total plasma" = "#BF3F4C",
  "EV-associated" = "#2A9D8F",
  "Shed" = "#D78B2B"
)

assay_map <- data.frame(
  assay = c("Total plasma", "EV-associated", "Shed"),
  raw_variable = c("total_sdc1", "ev_sdc1", "shed_sdc1"),
  log2_variable = c("log2_total_sdc1", "log2_ev_sdc1", "log2_shed_sdc1"),
  stringsAsFactors = FALSE
)

long_rows <- lapply(seq_len(nrow(assay_map)), function(index) {
  data.frame(
    row_uid = data$row_uid,
    cluster_id = data$cluster_id,
    group = data$group,
    assay = assay_map$assay[index],
    sdc1_pg_ml = data[[assay_map$raw_variable[index]]],
    log2_sdc1 = data[[assay_map$log2_variable[index]]],
    stringsAsFactors = FALSE
  )
})
long_data <- do.call(rbind, long_rows)
long_data$assay <- factor(long_data$assay, levels = assay_map$assay)
long_data$group <- factor(long_data$group, levels = group_order)

summary_rows <- list()
for (assay_name in assay_map$assay) {
  for (group_name in group_order) {
    values <- long_data$sdc1_pg_ml[
      long_data$assay == assay_name & long_data$group == group_name
    ]
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      assay = assay_name,
      group = group_name,
      n = length(values),
      median = median(values),
      q1 = unname(quantile(values, 0.25)),
      q3 = unname(quantile(values, 0.75)),
      minimum = min(values),
      maximum = max(values),
      stringsAsFactors = FALSE
    )
  }
}
group_summary <- do.call(rbind, summary_rows)

hc3_vcov <- function(fit) {
  design <- model.matrix(fit)
  residual <- residuals(fit)
  leverage <- hatvalues(fit)
  inverse_information <- solve(crossprod(design))
  scaled_residual_sq <- (residual / (1 - leverage))^2
  meat <- crossprod(design, design * scaled_residual_sq)
  inverse_information %*% meat %*% inverse_information
}

comparators <- c("Healthy control", "AMI", "Angina", "PE")
effect_rows <- list()
for (index in seq_len(nrow(assay_map))) {
  outcome_variable <- assay_map$log2_variable[index]
  fit_data <- data.frame(
    log2_sdc1 = data[[outcome_variable]],
    group = data$group,
    age_years = data$age_years,
    female = data$female,
    hypertension = data$hypertension
  )
  fit <- lm(log2_sdc1 ~ group + age_years + female + hypertension, data = fit_data)
  robust_vcov <- hc3_vcov(fit)
  coefficients <- coef(fit)

  for (comparator in comparators) {
    contrast <- setNames(rep(0, length(coefficients)), names(coefficients))
    contrast["groupAAD"] <- 1
    if (comparator != "Healthy control") contrast[paste0("group", comparator)] <- -1
    estimate <- sum(contrast * coefficients)
    standard_error <- sqrt(drop(t(contrast) %*% robust_vcov %*% contrast))
    z_value <- estimate / standard_error
    p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
    effect_rows[[length(effect_rows) + 1L]] <- data.frame(
      assay = assay_map$assay[index],
      comparison = paste0("AAD vs ", ifelse(comparator == "Healthy control", "healthy controls", comparator)),
      comparator = comparator,
      adjusted_log2_fold_change = estimate,
      standard_error_hc3 = standard_error,
      lower_95 = estimate - qnorm(0.975) * standard_error,
      upper_95 = estimate + qnorm(0.975) * standard_error,
      geometric_mean_ratio = 2^estimate,
      p_value = p_value,
      model_n = nrow(fit_data),
      stringsAsFactors = FALSE
    )
  }
}
adjusted_effects <- do.call(rbind, effect_rows)
adjusted_effects$bh_q_value <- p.adjust(adjusted_effects$p_value, method = "BH")
comparison_levels <- c("AAD vs PE", "AAD vs Angina", "AAD vs AMI", "AAD vs healthy controls")
adjusted_effects$comparison <- factor(adjusted_effects$comparison, levels = comparison_levels)
adjusted_effects$assay <- factor(adjusted_effects$assay, levels = assay_map$assay)

correlation_pairs <- list(
  list(pair = "Total plasma vs EV-associated", x = "log2_total_sdc1", y = "log2_ev_sdc1"),
  list(pair = "Total plasma vs shed", x = "log2_total_sdc1", y = "log2_shed_sdc1"),
  list(pair = "EV-associated vs shed", x = "log2_ev_sdc1", y = "log2_shed_sdc1")
)
n_bootstrap <- 2000L
correlation_rows <- list()
correlation_draw_rows <- vector("list", length(correlation_pairs))
scatter_rows <- list()

for (pair_index in seq_along(correlation_pairs)) {
  specification <- correlation_pairs[[pair_index]]
  x <- data[[specification$x]]
  y <- data[[specification$y]]
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(20263721L + pair_index)
  draws <- replicate(n_bootstrap, {
    sampled <- sample.int(nrow(data), nrow(data), replace = TRUE)
    suppressWarnings(cor(x[sampled], y[sampled], method = "spearman"))
  })
  interval <- unname(quantile(draws, c(0.025, 0.975), na.rm = TRUE))
  correlation_rows[[pair_index]] <- data.frame(
    pair = specification$pair,
    n = nrow(data),
    spearman_rho = rho,
    lower_95 = interval[1],
    upper_95 = interval[2],
    stringsAsFactors = FALSE
  )
  correlation_draw_rows[[pair_index]] <- data.frame(
    pair = specification$pair,
    bootstrap_replicate = seq_len(n_bootstrap),
    spearman_rho = draws,
    stringsAsFactors = FALSE
  )
  scatter_rows[[pair_index]] <- data.frame(
    row_uid = data$row_uid,
    cluster_id = data$cluster_id,
    group = data$group,
    pair = specification$pair,
    x = x,
    y = y,
    stringsAsFactors = FALSE
  )
}
correlation_summary <- do.call(rbind, correlation_rows)
correlation_draws <- do.call(rbind, correlation_draw_rows)
scatter_data <- do.call(rbind, scatter_rows)
pair_levels <- vapply(correlation_pairs, `[[`, character(1), "pair")
correlation_summary$pair <- factor(correlation_summary$pair, levels = pair_levels)
correlation_draws$pair <- factor(correlation_draws$pair, levels = pair_levels)
scatter_data$pair <- factor(scatter_data$pair, levels = pair_levels)
scatter_data$group <- factor(scatter_data$group, levels = group_order)
correlation_summary$label <- sprintf(
  "Spearman rho = %.2f\n95%% CI %.2f to %.2f",
  correlation_summary$spearman_rho,
  correlation_summary$lower_95,
  correlation_summary$upper_95
)

write.csv(
  data[, c("row_uid", "cluster_id", "group", "age_years", "female", "hypertension", "total_sdc1", "ev_sdc1", "shed_sdc1")],
  file.path(source_dir, "fig3_paired_compartment_patient_data.csv"),
  row.names = FALSE
)
write.csv(long_data, file.path(source_dir, "fig3_paired_compartment_long_data.csv"), row.names = FALSE)
write.csv(group_summary, file.path(result_dir, "fig3_paired_compartment_group_summary.csv"), row.names = FALSE)
write.csv(adjusted_effects, file.path(result_dir, "fig3_paired_compartment_adjusted_effects.csv"), row.names = FALSE)
write.csv(correlation_summary, file.path(result_dir, "fig3_paired_compartment_correlations.csv"), row.names = FALSE)
write.csv(correlation_draws, file.path(source_dir, "fig3_paired_compartment_correlation_bootstrap_draws.csv"), row.names = FALSE)
write.csv(scatter_data, file.path(source_dir, "fig3_paired_compartment_scatter_source_data.csv"), row.names = FALSE)

theme_pub <- function() {
  theme_classic(base_size = 6.7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#202124"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#202124"),
      axis.title = element_text(size = 6.5),
      axis.text = element_text(size = 5.9, colour = "#202124"),
      legend.title = element_blank(),
      legend.text = element_text(size = 5.6),
      strip.text = element_text(size = 6.2, face = "bold"),
      plot.title = element_text(size = 7.2, face = "bold"),
      plot.subtitle = element_text(size = 5.8, colour = "#687078"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

p_a <- ggplot() +
  annotate("rect", xmin = 0.3, xmax = 2.45, ymin = 1.05, ymax = 2.95, fill = "#F1F3F4", colour = "#60676D", linewidth = 0.45) +
  annotate("text", x = 1.375, y = 2.30, label = "Clinical-evaluation cohort", family = "Arial", fontface = "bold", size = 2.6) +
  annotate("text", x = 1.375, y = 1.72, label = "n = 415", family = "Arial", size = 2.55) +
  annotate("segment", x = 2.55, xend = 3.10, y = 2.0, yend = 2.0, arrow = arrow(length = grid::unit(2.2, "mm")), linewidth = 0.5, colour = "#60676D") +
  annotate("text", x = 2.82, y = 2.45, label = "nested", family = "Arial", size = 2.2, colour = "#687078") +
  annotate("rect", xmin = 3.20, xmax = 5.85, ymin = 0.70, ymax = 3.30, fill = "#E8F2F1", colour = "#2A9D8F", linewidth = 0.65) +
  annotate("text", x = 4.525, y = 2.72, label = "Paired-compartment subset", family = "Arial", fontface = "bold", size = 2.7) +
  annotate("text", x = 4.525, y = 2.17, label = "same 96 participants", family = "Arial", fontface = "bold", size = 2.55) +
  annotate("text", x = 4.525, y = 1.56, label = "HC 19 | AAD 21 | AMI 17", family = "Arial", size = 2.25) +
  annotate("text", x = 4.525, y = 1.12, label = "Angina 20 | PE 19", family = "Arial", size = 2.25) +
  annotate("segment", x = 5.85, xend = 6.35, y = 2.0, yend = 2.0, linewidth = 0.45, colour = "#60676D") +
  annotate("segment", x = 6.35, xend = 6.35, y = 0.68, yend = 3.32, linewidth = 0.45, colour = "#60676D") +
  annotate("segment", x = 6.35, xend = 6.72, y = 0.95, yend = 0.95, arrow = arrow(length = grid::unit(1.8, "mm")), linewidth = 0.45, colour = "#60676D") +
  annotate("segment", x = 6.35, xend = 6.72, y = 2.00, yend = 2.00, arrow = arrow(length = grid::unit(1.8, "mm")), linewidth = 0.45, colour = "#60676D") +
  annotate("segment", x = 6.35, xend = 6.72, y = 3.05, yend = 3.05, arrow = arrow(length = grid::unit(1.8, "mm")), linewidth = 0.45, colour = "#60676D") +
  annotate("rect", xmin = 6.78, xmax = 9.75, ymin = 2.63, ymax = 3.47, fill = "#FBEAEC", colour = assay_palette[["Total plasma"]], linewidth = 0.5) +
  annotate("rect", xmin = 6.78, xmax = 9.75, ymin = 1.58, ymax = 2.42, fill = "#E4F3F1", colour = assay_palette[["EV-associated"]], linewidth = 0.5) +
  annotate("rect", xmin = 6.78, xmax = 9.75, ymin = 0.53, ymax = 1.37, fill = "#FAEEDC", colour = assay_palette[["Shed"]], linewidth = 0.5) +
  annotate("text", x = 8.265, y = 3.05, label = "Total-plasma SDC1  n = 96", family = "Arial", fontface = "bold", size = 2.45) +
  annotate("text", x = 8.265, y = 2.00, label = "EV-associated SDC1  n = 96", family = "Arial", fontface = "bold", size = 2.45) +
  annotate("text", x = 8.265, y = 0.95, label = "Shed SDC1  n = 96", family = "Arial", fontface = "bold", size = 2.45) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.3, 3.7), clip = "off") +
  labs(
    title = "Same-participant assay-compartment design",
    subtitle = "The paired subset is nested within the clinical-evaluation cohort; it is not an independent validation cohort"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 7.4, face = "bold", margin = margin(b = 2)),
    plot.subtitle = element_text(size = 5.9, colour = "#687078"),
    plot.margin = margin(4, 6, 2, 6)
  )

group_axis_labels <- c(
  "Healthy control" = "HC\nn=19",
  "AAD" = "AAD\nn=21",
  "AMI" = "AMI\nn=17",
  "Angina" = "Angina\nn=20",
  "PE" = "PE\nn=19"
)

iqr_summary <- function(values) {
  data.frame(
    y = median(values),
    ymin = unname(quantile(values, 0.25)),
    ymax = unname(quantile(values, 0.75))
  )
}

make_distribution_plot <- function(assay_name, title, y_title) {
  plot_data <- long_data[long_data$assay == assay_name, , drop = FALSE]
  ggplot(plot_data, aes(group, sdc1_pg_ml, colour = group)) +
    geom_point(position = position_jitter(width = 0.13, height = 0, seed = 20260721), size = 1.15, alpha = 0.62) +
    stat_summary(fun.data = iqr_summary, geom = "errorbar", width = 0.18, linewidth = 0.55, colour = "#202124") +
    stat_summary(fun = median, geom = "point", shape = 21, size = 2.0, stroke = 0.55, fill = "white", colour = "#202124") +
    scale_colour_manual(values = group_palette, guide = "none") +
    scale_x_discrete(labels = group_axis_labels) +
    scale_y_log10(labels = scales::label_number(big.mark = ",")) +
    labs(x = NULL, y = y_title, title = title, subtitle = "Median and interquartile range") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

p_b <- make_distribution_plot("Total plasma", "Total-plasma SDC1", "SDC1 (pg/mL; log scale)")
p_c <- make_distribution_plot("EV-associated", "EV-associated SDC1", "SDC1 (pg/mL; log scale)")
p_d <- make_distribution_plot("Shed", "Shed SDC1", "SDC1 (pg/mL; log scale)")

p_e <- ggplot(adjusted_effects, aes(adjusted_log2_fold_change, comparison, colour = assay)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#8F969C", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    orientation = "y",
    position = position_dodge(width = 0.54),
    width = 0,
    linewidth = 0.55
  ) +
  geom_point(position = position_dodge(width = 0.54), size = 1.9) +
  scale_colour_manual(values = assay_palette) +
  labs(
    x = "Adjusted log2 fold change (95% HC3 CI)",
    y = NULL,
    title = "AAD-versus-comparator effects",
    subtitle = "Adjusted for age, sex and hypertension"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(colour = guide_legend(nrow = 1))

p_f <- ggplot(scatter_data, aes(x, y, colour = group)) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.45, colour = "#8F969C") +
  geom_point(size = 0.95, alpha = 0.60) +
  geom_label(
    data = correlation_summary,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.08,
    size = 2.05,
    linewidth = 0,
    fill = scales::alpha("white", 0.82)
  ) +
  facet_wrap(
    ~pair,
    nrow = 1,
    scales = "free",
    labeller = as_labeller(c(
      "Total plasma vs EV-associated" = "Total plasma vs\nEV-associated",
      "Total plasma vs shed" = "Total plasma vs\nshed",
      "EV-associated vs shed" = "EV-associated vs\nshed"
    ))
  ) +
  scale_colour_manual(values = group_palette) +
  labs(
    x = "First SDC1 measure (log2 pg/mL)",
    y = "Second SDC1 measure (log2 pg/mL)",
    title = "Same-participant cross-compartment associations",
    subtitle = "Spearman correlations with 2,000 patient-bootstrap intervals"
  ) +
  theme_pub() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    strip.background = element_rect(fill = "#F1F3F4", colour = NA),
    panel.spacing.x = grid::unit(1.6, "mm")
  ) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 1.5)))

middle_row <- p_b | p_c | p_d
figure3 <- p_a / middle_row / p_e / p_f +
  plot_layout(heights = c(0.62, 1.02, 0.78, 1.18)) +
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

figure_stem <- file.path(figure_dir, "Figure_3_paired_compartment_SDC1")
save_figure(figure3, figure_stem, width_mm = 183, height_mm = 230, dpi = 600)

panel_list <- list(A = p_a, B = p_b, C = p_c, D = p_d, E = p_e, F = p_f)
panel_dimensions <- list(
  A = c(183, 52),
  B = c(61, 75),
  C = c(61, 75),
  D = c(61, 75),
  E = c(183, 68),
  F = c(183, 88)
)
for (panel_name in names(panel_list)) {
  dimensions <- panel_dimensions[[panel_name]]
  save_figure(
    panel_list[[panel_name]],
    file.path(panel_dir, paste0("Figure_3", panel_name)),
    width_mm = dimensions[1],
    height_mm = dimensions[2],
    dpi = 600
  )
}

qa_checks <- c(
  paired_records_96 = nrow(data) == 96L,
  unique_patient_keys_96 = length(unique(data$cluster_id)) == 96L,
  final_group_counts_match = identical(as.integer(table(data$group)), unname(expected_counts)),
  all_three_sdc1_layers_complete = all(complete.cases(data[, c("total_sdc1", "ev_sdc1", "shed_sdc1")])),
  all_sdc1_values_positive = all(data[, c("total_sdc1", "ev_sdc1", "shed_sdc1")] > 0),
  adjusted_effect_rows_12 = nrow(adjusted_effects) == 12L,
  adjusted_effects_finite = all(is.finite(unlist(adjusted_effects[, c("adjusted_log2_fold_change", "lower_95", "upper_95", "p_value", "bh_q_value")]))),
  correlation_pairs_3 = nrow(correlation_summary) == 3L,
  correlation_bootstrap_draws_6000 = nrow(correlation_draws) == 6000L,
  correlation_intervals_finite = all(is.finite(unlist(correlation_summary[, c("spearman_rho", "lower_95", "upper_95")]))),
  main_exports_present = all(file.exists(paste0(figure_stem, c(".pdf", ".svg", ".tiff", ".png")))),
  standalone_exports_present = all(vapply(
    LETTERS[1:6],
    function(panel_name) all(file.exists(file.path(panel_dir, paste0("Figure_3", panel_name, c(".pdf", ".svg", ".tiff", ".png"))))),
    logical(1)
  ))
)
qa_table <- data.frame(check = names(qa_checks), pass = as.logical(qa_checks))
write.csv(qa_table, file.path(qa_dir, "FINAL_FIGURE3_PAIRED_QA.csv"), row.names = FALSE)

summary_lines <- c(
  "Figure 3 paired-compartment SDC1 analysis",
  paste0("Overall computational QA: ", if (all(qa_checks)) "PASS" else "FAIL"),
  "Population: 96 participants nested within the 415-record clinical-evaluation cohort",
  "Final counts: HC 19; AAD 21; AMI 17; Angina 20; PE 19",
  "All 96 participants have total-plasma, EV-associated and shed SDC1 measurements",
  "Adjusted effects: age, sex and hypertension; HC3 confidence intervals",
  "Correlations: Spearman rho with 2,000 patient-level bootstrap intervals",
  "No ROC, diagnostic cutoff, prediction model or independent-validation claim",
  paste0("Main figure: ", figure_stem)
)
writeLines(summary_lines, file.path(log_dir, "fig3_paired_analysis_summary.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "fig3_paired_sessionInfo.txt"), useBytes = TRUE)

cat(paste(summary_lines, collapse = "\n"), "\n")
if (!all(qa_checks)) quit(status = 1)
