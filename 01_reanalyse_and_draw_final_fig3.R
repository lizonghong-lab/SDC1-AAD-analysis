options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(Sys.getenv("SDC1_FINAL_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
legacy_rlib <- Sys.getenv(
  "SDC1_RLIB",
  unset = ""
)
if (dir.exists(legacy_rlib)) .libPaths(c(legacy_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
})

set.seed(20260721)

data_path <- file.path(root, "02_analysis_data", "discovery_clinical_deidentified.csv")
result_dir <- file.path(root, "04_results", "statistical_results")
source_dir <- file.path(root, "04_results", "source_data")
figure_dir <- file.path(root, "04_results", "figures", "main")
panel_dir <- file.path(root, "04_results", "figures", "panels", "figure_3")
log_dir <- file.path(root, "04_results", "logs")
for (directory in c(result_dir, source_dir, figure_dir, panel_dir, log_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")

hc3_vcov <- function(fit) {
  x <- model.matrix(fit)
  residual <- residuals(fit)
  leverage <- pmin(hatvalues(fit), 1 - 1e-8)
  omega <- (residual / (1 - leverage))^2
  bread <- solve(crossprod(x))
  meat <- crossprod(x, sweep(x, 1, omega, "*"))
  bread %*% meat %*% bread
}

contrast_result <- function(fit, vcov_hc3, contrast) {
  beta <- coef(fit)
  estimate <- as.numeric(crossprod(contrast, beta))
  variance <- as.numeric(t(contrast) %*% vcov_hc3 %*% contrast)
  standard_error <- sqrt(max(variance, 0))
  statistic <- estimate / standard_error
  p_value <- 2 * pnorm(abs(statistic), lower.tail = FALSE)
  data.frame(
    estimate_log2 = estimate,
    robust_se = standard_error,
    lower_log2_95 = estimate - qnorm(0.975) * standard_error,
    upper_log2_95 = estimate + qnorm(0.975) * standard_error,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}

adjusted_predictions <- function(fit, vcov_hc3, newdata) {
  terms_no_response <- delete.response(terms(fit))
  x_new <- model.matrix(terms_no_response, newdata, contrasts.arg = fit$contrasts)
  estimate <- as.numeric(x_new %*% coef(fit))
  standard_error <- sqrt(pmax(diag(x_new %*% vcov_hc3 %*% t(x_new)), 0))
  data.frame(
    estimate_log2 = estimate,
    robust_se = standard_error,
    lower_log2_95 = estimate - qnorm(0.975) * standard_error,
    upper_log2_95 = estimate + qnorm(0.975) * standard_error,
    stringsAsFactors = FALSE
  )
}

prepare_covariates <- function(data) {
  data$z_age <- as.numeric(scale(data$age_years))
  data$z_log1p_wbc <- as.numeric(scale(log1p(data$wbc_10e9_l)))
  data$z_log1p_crp <- as.numeric(scale(log1p(data$crp_mg_l)))
  data
}

format_probability <- function(value) {
  if (!is.finite(value)) return("NA")
  if (value < 0.0001) return("<0.0001")
  if (value < 0.001) return(sprintf("%.4f", value))
  sprintf("%.3f", value)
}

discovery <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
numeric_columns <- c("ev_sdc1", "d_dimer_mg_l", "ctnt_source_units", "age_years", "female", "hypertension", "wbc_10e9_l", "crp_mg_l")
discovery[numeric_columns] <- lapply(discovery[numeric_columns], as.numeric)
discovery$group <- factor(discovery$group, levels = c("Healthy control", "AAD", "AMI"))

biomarker_specs <- list(
  ev_sdc1 = list(title = "AF4-EV SDC1", y = "AF4-EV SDC1 (pg/mL)", display = "original"),
  d_dimer_mg_l = list(title = "D-dimer", y = "log2 D-dimer (mg/L)", display = "log2"),
  ctnt_source_units = list(title = "hs-cTnT", y = "log2 hs-cTnT (source units)", display = "log2")
)

group_fits <- list()
group_means <- list()
group_contrasts <- list()
group_plot_data <- list()
completeness <- list()

for (outcome in names(biomarker_specs)) {
  required <- c(outcome, "group", "age_years", "female", "hypertension", "wbc_10e9_l", "crp_mg_l")
  keep <- complete.cases(discovery[, required]) & is.finite(discovery[[outcome]]) & discovery[[outcome]] > 0
  data <- droplevels(discovery[keep, , drop = FALSE])
  data <- prepare_covariates(data)
  data$log2_outcome <- log2(data[[outcome]])
  fit <- lm(log2_outcome ~ group + z_age + female + hypertension + z_log1p_wbc + z_log1p_crp, data = data)
  vcov_hc3 <- hc3_vcov(fit)

  common <- data.frame(
    group = factor(levels(discovery$group), levels = levels(discovery$group)),
    z_age = 0,
    female = mean(data$female),
    hypertension = mean(data$hypertension),
    z_log1p_wbc = 0,
    z_log1p_crp = 0
  )
  prediction <- cbind(common["group"], adjusted_predictions(fit, vcov_hc3, common))
  prediction$outcome <- outcome
  prediction$n_complete <- nrow(data)
  prediction$n_group <- as.integer(table(data$group)[as.character(prediction$group)])
  prediction$estimate_display <- if (biomarker_specs[[outcome]]$display == "original") 2^prediction$estimate_log2 else prediction$estimate_log2
  prediction$lower_display_95 <- if (biomarker_specs[[outcome]]$display == "original") 2^prediction$lower_log2_95 else prediction$lower_log2_95
  prediction$upper_display_95 <- if (biomarker_specs[[outcome]]$display == "original") 2^prediction$upper_log2_95 else prediction$upper_log2_95

  x_new <- model.matrix(delete.response(terms(fit)), common, contrasts.arg = fit$contrasts)
  rownames(x_new) <- as.character(common$group)
  planned <- list(
    "AAD versus Healthy control" = x_new["AAD", ] - x_new["Healthy control", ],
    "AAD versus AMI" = x_new["AAD", ] - x_new["AMI", ]
  )
  contrasts <- do.call(rbind, lapply(names(planned), function(label) {
    result <- contrast_result(fit, vcov_hc3, planned[[label]])
    result$outcome <- outcome
    result$contrast <- label
    result$n_complete <- nrow(data)
    result$geometric_mean_ratio <- 2^result$estimate_log2
    result$ratio_lower_95 <- 2^result$lower_log2_95
    result$ratio_upper_95 <- 2^result$upper_log2_95
    result
  }))

  data$plot_value <- if (biomarker_specs[[outcome]]$display == "original") data[[outcome]] else data$log2_outcome
  data$outcome <- outcome
  group_fits[[outcome]] <- fit
  group_means[[outcome]] <- prediction
  group_contrasts[[outcome]] <- contrasts
  group_plot_data[[outcome]] <- data[, c("row_uid", "source_id", "group", "hypertension", "plot_value", "outcome")]
  completeness[[outcome]] <- data.frame(
    analysis = outcome,
    total_records = nrow(discovery),
    complete_records = nrow(data),
    healthy_control = sum(data$group == "Healthy control"),
    aad = sum(data$group == "AAD"),
    ami = sum(data$group == "AMI"),
    adjustment_set = "Age, sex, hypertension, log1p(WBC), log1p(CRP)",
    stringsAsFactors = FALSE
  )
}

group_means <- do.call(rbind, group_means)
group_contrasts <- do.call(rbind, group_contrasts)
group_contrasts$bh_q_across_six <- p.adjust(group_contrasts$p_value, method = "BH")
group_plot_data <- do.call(rbind, group_plot_data)
completeness <- do.call(rbind, completeness)

htn_means <- list()
htn_effects <- list()
htn_plot_data <- list()

for (group_name in levels(discovery$group)) {
  required <- c("ev_sdc1", "hypertension", "age_years", "female", "wbc_10e9_l", "crp_mg_l")
  subset <- discovery[discovery$group == group_name, , drop = FALSE]
  keep <- complete.cases(subset[, required]) & subset$ev_sdc1 > 0
  data <- subset[keep, , drop = FALSE]
  data <- prepare_covariates(data)
  data$log2_ev_sdc1 <- log2(data$ev_sdc1)
  fit <- lm(log2_ev_sdc1 ~ hypertension + z_age + female + z_log1p_wbc + z_log1p_crp, data = data)
  vcov_hc3 <- hc3_vcov(fit)

  common <- data.frame(
    hypertension = c(0, 1),
    z_age = 0,
    female = mean(data$female),
    z_log1p_wbc = 0,
    z_log1p_crp = 0
  )
  prediction <- adjusted_predictions(fit, vcov_hc3, common)
  prediction$hypertension <- common$hypertension
  prediction$hypertension_label <- factor(ifelse(common$hypertension == 1, "Hypertension", "No hypertension"), levels = c("No hypertension", "Hypertension"))
  prediction$group <- group_name
  prediction$n_group <- c(sum(data$hypertension == 0), sum(data$hypertension == 1))
  prediction$estimate_display <- 2^prediction$estimate_log2
  prediction$lower_display_95 <- 2^prediction$lower_log2_95
  prediction$upper_display_95 <- 2^prediction$upper_log2_95

  coefficient <- rep(0, length(coef(fit)))
  names(coefficient) <- names(coef(fit))
  coefficient["hypertension"] <- 1
  effect <- contrast_result(fit, vcov_hc3, coefficient)
  effect$group <- group_name
  effect$n_complete <- nrow(data)
  effect$n_no_hypertension <- sum(data$hypertension == 0)
  effect$n_hypertension <- sum(data$hypertension == 1)
  effect$geometric_mean_ratio <- 2^effect$estimate_log2
  effect$ratio_lower_95 <- 2^effect$lower_log2_95
  effect$ratio_upper_95 <- 2^effect$upper_log2_95
  effect$adjustment_set <- "Age, sex, log1p(WBC), log1p(CRP)"

  data$hypertension_label <- factor(ifelse(data$hypertension == 1, "Hypertension", "No hypertension"), levels = c("No hypertension", "Hypertension"))
  data$group <- group_name
  htn_means[[group_name]] <- prediction
  htn_effects[[group_name]] <- effect
  htn_plot_data[[group_name]] <- data[, c("row_uid", "source_id", "group", "hypertension", "hypertension_label", "ev_sdc1")]
}

htn_means <- do.call(rbind, htn_means)
htn_effects <- do.call(rbind, htn_effects)
htn_effects$bh_q_across_three <- p.adjust(htn_effects$p_value, method = "BH")
htn_plot_data <- do.call(rbind, htn_plot_data)

write_csv(group_means, file.path(result_dir, "fig3_adjusted_group_means.csv"))
write_csv(group_contrasts, file.path(result_dir, "fig3_adjusted_group_contrasts.csv"))
write_csv(htn_means, file.path(result_dir, "fig3_adjusted_hypertension_means.csv"))
write_csv(htn_effects, file.path(result_dir, "fig3_adjusted_hypertension_effects.csv"))
write_csv(completeness, file.path(result_dir, "fig3_complete_case_counts.csv"))
write_csv(group_plot_data, file.path(source_dir, "fig3_panels_abc_source_data.csv"))
write_csv(htn_plot_data, file.path(source_dir, "fig3_panels_def_source_data.csv"))

palette <- c("Healthy control" = "#4B8F8A", "AAD" = "#BF3F4C", "AMI" = "#3F78A8")
htn_palette <- c("No hypertension" = "#B8BDC2", "Hypertension" = "#5C6470")

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6.2, colour = "black"),
      plot.title = element_text(size = 7, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = 6.1, colour = "#555555", margin = margin(b = 3)),
      plot.margin = margin(6, 5, 5, 5),
      legend.position = "none"
    )
}

add_bracket <- function(plot, x_start, x_end, y, height, label) {
  plot +
    annotate("segment", x = x_start, xend = x_end, y = y, yend = y, linewidth = 0.35) +
    annotate("segment", x = x_start, xend = x_start, y = y, yend = y - height, linewidth = 0.35) +
    annotate("segment", x = x_end, xend = x_end, y = y, yend = y - height, linewidth = 0.35) +
    annotate("text", x = (x_start + x_end) / 2, y = y + 0.02 * abs(y), label = label, size = 2.15, vjust = 0)
}

make_group_panel <- function(outcome) {
  raw <- group_plot_data[group_plot_data$outcome == outcome, , drop = FALSE]
  raw$group <- factor(raw$group, levels = levels(discovery$group))
  mean_data <- group_means[group_means$outcome == outcome, , drop = FALSE]
  mean_data$group <- factor(mean_data$group, levels = levels(discovery$group))
  comparison <- group_contrasts[group_contrasts$outcome == outcome, , drop = FALSE]
  n_lookup <- setNames(mean_data$n_group, as.character(mean_data$group))
  labels <- setNames(paste0(c("HC", "AAD", "AMI"), "\n(n=", n_lookup[levels(discovery$group)], ")"), levels(discovery$group))
  range_value <- range(c(raw$plot_value, mean_data$lower_display_95, mean_data$upper_display_95), na.rm = TRUE)
  span <- diff(range_value)
  if (!is.finite(span) || span == 0) span <- max(abs(range_value), 1)
  y_one <- range_value[2] + 0.10 * span
  y_two <- range_value[2] + 0.26 * span
  bracket_height <- 0.04 * span
  q_hc <- comparison$bh_q_across_six[comparison$contrast == "AAD versus Healthy control"]
  q_ami <- comparison$bh_q_across_six[comparison$contrast == "AAD versus AMI"]

  p <- ggplot(raw, aes(group, plot_value, colour = group)) +
    geom_jitter(width = 0.11, size = 0.72, alpha = 0.48) +
    geom_errorbar(
      data = mean_data,
      aes(x = group, ymin = lower_display_95, ymax = upper_display_95),
      inherit.aes = FALSE,
      width = 0.14, linewidth = 0.62, colour = "#202428"
    ) +
    geom_point(
      data = mean_data,
      aes(x = group, y = estimate_display),
      inherit.aes = FALSE, shape = 21, size = 2.3, stroke = 0.65,
      fill = "white", colour = "#202428"
    ) +
    scale_colour_manual(values = palette) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.34))) +
    labs(x = NULL, y = biomarker_specs[[outcome]]$y, title = biomarker_specs[[outcome]]$title, subtitle = "Adjusted for age, sex, hypertension, WBC and CRP") +
    coord_cartesian(clip = "off") +
    theme_pub()
  p <- add_bracket(p, 1, 2, y_one, bracket_height, paste0("BH q", format_probability(q_hc)))
  add_bracket(p, 2, 3, y_two, bracket_height, paste0("BH q", format_probability(q_ami)))
}

make_htn_panel <- function(group_name) {
  raw <- htn_plot_data[htn_plot_data$group == group_name, , drop = FALSE]
  mean_data <- htn_means[htn_means$group == group_name, , drop = FALSE]
  effect <- htn_effects[htn_effects$group == group_name, , drop = FALSE]
  labels <- setNames(paste0(c("No HTN", "HTN"), "\n(n=", mean_data$n_group, ")"), levels(raw$hypertension_label))
  range_value <- range(c(raw$ev_sdc1, mean_data$lower_display_95, mean_data$upper_display_95), na.rm = TRUE)
  span <- diff(range_value)
  if (!is.finite(span) || span == 0) span <- max(abs(range_value), 1)
  y <- range_value[2] + 0.13 * span
  bracket_height <- 0.05 * span
  p <- ggplot(raw, aes(hypertension_label, ev_sdc1, colour = hypertension_label)) +
    geom_jitter(width = 0.10, size = 0.75, alpha = 0.50) +
    geom_errorbar(
      data = mean_data,
      aes(x = hypertension_label, ymin = lower_display_95, ymax = upper_display_95),
      inherit.aes = FALSE, width = 0.14, linewidth = 0.62, colour = "#202428"
    ) +
    geom_point(
      data = mean_data,
      aes(x = hypertension_label, y = estimate_display),
      inherit.aes = FALSE, shape = 21, size = 2.3, stroke = 0.65,
      fill = "white", colour = "#202428"
    ) +
    scale_colour_manual(values = htn_palette) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.26))) +
    labs(x = NULL, y = "AF4-EV SDC1 (pg/mL)", title = group_name, subtitle = "Adjusted for age, sex, WBC and CRP") +
    coord_cartesian(clip = "off") +
    theme_pub()
  add_bracket(p, 1, 2, y, bracket_height, paste0("BH q", format_probability(effect$bh_q_across_three)))
}

p_a <- make_group_panel("ev_sdc1")
p_b <- make_group_panel("d_dimer_mg_l")
p_c <- make_group_panel("ctnt_source_units")
p_d <- make_htn_panel("AAD")
p_e <- make_htn_panel("AMI")
p_f <- make_htn_panel("Healthy control")

figure <- (p_a | p_b | p_c) / (p_d | p_e | p_f) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 8, family = "Arial"))) +
  plot_layout(heights = c(1, 1))

save_figure <- function(plot, stem, width_mm = 183, height_mm = 132, dpi = 600) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = cairo_pdf, family = "Arial")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, device = svglite::svglite)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300, device = ragg::agg_png)
  ggsave(paste0(stem, ".tiff"), plot, width = width, height = height, dpi = dpi, compression = "lzw", device = ragg::agg_tiff)
}

stem <- file.path(figure_dir, "Figure_3_discovery_covariate_adjusted")
save_figure(figure, stem)

panel_list <- list(A = p_a, B = p_b, C = p_c, D = p_d, E = p_e, F = p_f)
for (panel_name in names(panel_list)) {
  save_figure(
    panel_list[[panel_name]],
    file.path(panel_dir, paste0("Figure_3", panel_name)),
    width_mm = 61,
    height_mm = 66,
    dpi = 600
  )
}

log_lines <- c(
  paste("Input discovery records:", nrow(discovery)),
  paste("A-C complete-case records:", paste(completeness$analysis, completeness$complete_records, sep = "=", collapse = "; ")),
  paste("AF4-EV SDC1 AAD-vs-HC GMR:", round(group_contrasts$geometric_mean_ratio[group_contrasts$outcome == "ev_sdc1" & group_contrasts$contrast == "AAD versus Healthy control"], 4)),
  paste("AF4-EV SDC1 AAD-vs-AMI GMR:", round(group_contrasts$geometric_mean_ratio[group_contrasts$outcome == "ev_sdc1" & group_contrasts$contrast == "AAD versus AMI"], 4)),
  paste("D-F hypertension GMRs:", paste(htn_effects$group, round(htn_effects$geometric_mean_ratio, 4), sep = "=", collapse = "; ")),
  "Input authority: FINAL_20260721 lock and its derived discovery analysis table.",
  "Available-case rule: each model includes every record complete for that model's prespecified variables.",
  "Panel G removed: no ROC curve or AUROC is produced by this script.",
  "hs-cTnT unit remains source units pending author verification."
)
writeLines(log_lines, file.path(log_dir, "fig3_adjusted_summary.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "fig3_sessionInfo.txt"), useBytes = TRUE)
cat(paste(log_lines, collapse = "\n"), "\n")
