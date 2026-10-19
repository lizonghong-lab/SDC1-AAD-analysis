options(stringsAsFactors = FALSE, scipen = 999)

project_root <- normalizePath(
  Sys.getenv("SDC1_AAD_ONLY_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)
workspace_root <- normalizePath(file.path(project_root, ".."), winslash = "/", mustWork = TRUE)
input_file <- file.path(
  workspace_root,
  "SDC1_AAD_REANALYSIS_FINAL_20260721",
  "02_analysis_data",
  "paired_sdc1_subset.csv"
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
})

data_dir <- file.path(project_root, "01_data")
result_dir <- file.path(project_root, "02_results")
code_dir <- file.path(project_root, "03_code", "R")
figure_dir <- file.path(project_root, "04_figures")
qa_dir <- file.path(project_root, "05_qa")
for (directory in c(data_dir, result_dir, code_dir, figure_dir, qa_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

set.seed(20260727)

paired <- read.csv(input_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
required <- c("row_uid", "cluster_id", "group", "total_sdc1", "ev_sdc1", "shed_sdc1")
missing_columns <- setdiff(required, names(paired))
if (length(missing_columns)) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

aad <- paired[paired$group == "AAD", required, drop = FALSE]
for (column in c("total_sdc1", "ev_sdc1", "shed_sdc1")) {
  aad[[column]] <- as.numeric(aad[[column]])
}

if (nrow(aad) != 21L) stop("Expected 21 AAD participants; found ", nrow(aad))
if (length(unique(aad$cluster_id)) != 21L) stop("Expected 21 unique AAD patient keys")
if (any(!complete.cases(aad))) stop("AAD paired-compartment fields are incomplete")
if (any(aad[, c("total_sdc1", "ev_sdc1", "shed_sdc1")] <= 0)) {
  stop("All SDC1 measurements must be positive")
}

aad$participant_index <- sprintf("AAD%02d", seq_len(nrow(aad)))
aad$log2_total_sdc1 <- log2(aad$total_sdc1)
aad$log2_ev_sdc1 <- log2(aad$ev_sdc1)
aad$log2_shed_sdc1 <- log2(aad$shed_sdc1)

pair_specifications <- list(
  list(
    pair = "Total plasma vs EV-associated",
    x = "log2_total_sdc1",
    y = "log2_ev_sdc1",
    x_label = "Total-plasma SDC1 (log2 pg/mL)",
    y_label = "EV-associated SDC1 (log2 pg/mL)"
  ),
  list(
    pair = "Total plasma vs shed",
    x = "log2_total_sdc1",
    y = "log2_shed_sdc1",
    x_label = "Total-plasma SDC1 (log2 pg/mL)",
    y_label = "Shed SDC1 (log2 pg/mL)"
  ),
  list(
    pair = "EV-associated vs shed",
    x = "log2_ev_sdc1",
    y = "log2_shed_sdc1",
    x_label = "EV-associated SDC1 (log2 pg/mL)",
    y_label = "Shed SDC1 (log2 pg/mL)"
  )
)

n_bootstrap <- 2000L
summary_rows <- vector("list", length(pair_specifications))
draw_rows <- vector("list", length(pair_specifications))
scatter_rows <- vector("list", length(pair_specifications))

for (index in seq_along(pair_specifications)) {
  specification <- pair_specifications[[index]]
  x <- aad[[specification$x]]
  y <- aad[[specification$y]]
  rho <- suppressWarnings(cor(x, y, method = "spearman"))

  set.seed(20264727L + index)
  draws <- replicate(n_bootstrap, {
    sampled <- sample.int(nrow(aad), nrow(aad), replace = TRUE)
    suppressWarnings(cor(x[sampled], y[sampled], method = "spearman"))
  })
  interval <- unname(quantile(draws, c(0.025, 0.975), na.rm = TRUE))

  summary_rows[[index]] <- data.frame(
    population = "AAD",
    pair = specification$pair,
    n = nrow(aad),
    spearman_rho = rho,
    lower_95 = interval[1],
    upper_95 = interval[2],
    interval_method = "2,000 patient-level percentile bootstrap resamples",
    stringsAsFactors = FALSE
  )
  draw_rows[[index]] <- data.frame(
    population = "AAD",
    pair = specification$pair,
    bootstrap_replicate = seq_len(n_bootstrap),
    spearman_rho = draws,
    stringsAsFactors = FALSE
  )
  scatter_rows[[index]] <- data.frame(
    participant_index = aad$participant_index,
    pair = specification$pair,
    x = x,
    y = y,
    stringsAsFactors = FALSE
  )
}

correlation_summary <- do.call(rbind, summary_rows)
bootstrap_draws <- do.call(rbind, draw_rows)
scatter_data <- do.call(rbind, scatter_rows)

write.csv(
  aad[, c(
    "participant_index", "total_sdc1", "ev_sdc1", "shed_sdc1",
    "log2_total_sdc1", "log2_ev_sdc1", "log2_shed_sdc1"
  )],
  file.path(data_dir, "AAD_only_paired_SDC1_patient_data.csv"),
  row.names = FALSE
)
write.csv(
  aad[, c("participant_index", "row_uid", "cluster_id")],
  file.path(qa_dir, "AAD_only_internal_patient_key_audit.csv"),
  row.names = FALSE
)
write.csv(
  scatter_data,
  file.path(data_dir, "AAD_only_paired_SDC1_scatter_source_data.csv"),
  row.names = FALSE
)
write.csv(
  correlation_summary,
  file.path(result_dir, "AAD_only_paired_SDC1_correlation_summary.csv"),
  row.names = FALSE
)
write.csv(
  bootstrap_draws,
  file.path(result_dir, "AAD_only_paired_SDC1_bootstrap_draws.csv"),
  row.names = FALSE
)

theme_pub <- function() {
  theme_classic(base_size = 7.2, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.42, colour = "#202124"),
      axis.ticks = element_line(linewidth = 0.38, colour = "#202124"),
      axis.title = element_text(size = 7.0, colour = "#202124"),
      axis.text = element_text(size = 6.3, colour = "#202124"),
      plot.title = element_text(size = 7.8, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 6.2, colour = "#687078", hjust = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(5, 6, 4, 6)
    )
}

aad_colour <- "#BF3F4C"
aad_outline <- "#9E3340"
trend_colour <- "#8F969C"

make_panel <- function(index) {
  specification <- pair_specifications[[index]]
  panel_data <- scatter_data[scatter_data$pair == specification$pair, , drop = FALSE]
  panel_summary <- correlation_summary[
    correlation_summary$pair == specification$pair,
    ,
    drop = FALSE
  ]
  x_range <- range(panel_data$x)
  y_range <- range(panel_data$y)
  x_text <- x_range[1] + 0.20 * diff(x_range)
  y_text <- y_range[2] - 0.02 * diff(y_range)
  label <- sprintf(
    "Spearman rho = %.2f\n95%% CI %.2f to %.2f",
    panel_summary$spearman_rho,
    panel_summary$lower_95,
    panel_summary$upper_95
  )

  ggplot(panel_data, aes(x, y)) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.55,
      colour = trend_colour
    ) +
    geom_point(
      shape = 21,
      size = 1.55,
      stroke = 0.28,
      fill = aad_colour,
      colour = aad_outline,
      alpha = 0.68
    ) +
    annotate(
      "text",
      x = x_text,
      y = y_text,
      label = label,
      hjust = 0,
      vjust = 1,
      family = "Arial",
      size = 2.35,
      lineheight = 1.18,
      colour = "#111111"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    labs(
      x = specification$x_label,
      y = specification$y_label,
      title = specification$pair,
      subtitle = "AAD only, n = 21"
    ) +
    theme_pub()
}

p_total_ev <- make_panel(1)
p_total_shed <- make_panel(2)
p_ev_shed <- make_panel(3)

two_panel <- p_total_ev | p_total_shed
two_panel <- two_panel +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 8.5,
        family = "Arial",
        colour = "#111111"
      )
    )
  )

three_panel <- p_total_ev | p_total_shed | p_ev_shed
three_panel <- three_panel +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 8.5,
        family = "Arial",
        colour = "#111111"
      )
    )
  )

save_figure <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"),
    width = width_in,
    height = height_in,
    family = "Arial"
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"),
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi,
    compression = "lzw"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    paste0(stem, ".png"),
    width = width_in,
    height = height_in,
    units = "in",
    res = 300
  )
  print(plot)
  dev.off()
}

two_panel_stem <- file.path(figure_dir, "AAD_only_paired_SDC1_correlations_two_panel")
three_panel_stem <- file.path(figure_dir, "AAD_only_paired_SDC1_correlations_three_panel")
save_figure(two_panel, two_panel_stem, width_mm = 183, height_mm = 76)
save_figure(three_panel, three_panel_stem, width_mm = 183, height_mm = 69)

panel_names <- c("total_vs_EV_associated", "total_vs_shed", "EV_associated_vs_shed")
panel_plots <- list(p_total_ev, p_total_shed, p_ev_shed)
for (index in seq_along(panel_plots)) {
  save_figure(
    panel_plots[[index]],
    file.path(figure_dir, paste0("AAD_only_", panel_names[index])),
    width_mm = 89,
    height_mm = 76
  )
}

recomputed_rho <- vapply(pair_specifications, function(specification) {
  suppressWarnings(cor(aad[[specification$x]], aad[[specification$y]], method = "spearman"))
}, numeric(1))
saved_rho <- correlation_summary$spearman_rho
expected_exports <- unlist(lapply(
  c(two_panel_stem, three_panel_stem, file.path(figure_dir, paste0("AAD_only_", panel_names))),
  function(stem) paste0(stem, c(".svg", ".pdf", ".tiff", ".png"))
))

qa_checks <- c(
  aad_records_21 = nrow(aad) == 21L,
  unique_aad_patient_keys_21 = length(unique(aad$cluster_id)) == 21L,
  three_assays_complete = all(complete.cases(aad[, c("total_sdc1", "ev_sdc1", "shed_sdc1")])),
  all_values_positive = all(aad[, c("total_sdc1", "ev_sdc1", "shed_sdc1")] > 0),
  correlation_pairs_3 = nrow(correlation_summary) == 3L,
  bootstrap_draws_6000 = nrow(bootstrap_draws) == 6000L,
  bootstrap_intervals_finite = all(is.finite(unlist(correlation_summary[, c("lower_95", "upper_95")]))),
  spearman_recalculation_matches = max(abs(recomputed_rho - saved_rho)) < 1e-12,
  all_exports_present = all(file.exists(expected_exports))
)
qa_table <- data.frame(check = names(qa_checks), pass = as.logical(qa_checks))
write.csv(qa_table, file.path(qa_dir, "AAD_only_paired_SDC1_QA.csv"), row.names = FALSE)

summary_lines <- c(
  "AAD-only paired-compartment SDC1 correlations",
  paste0("Overall QA: ", if (all(qa_checks)) "PASS" else "FAIL"),
  "Population: 21 AAD participants from the locked paired-compartment subset",
  "Measurements: total-plasma, EV-associated and shed SDC1",
  "Statistics: Spearman rho with 2,000 patient-level percentile-bootstrap resamples",
  "Grey lines: ordinary least-squares visual guides only",
  "No causal, kinetic, prognostic or compartment-equivalence claim"
)
writeLines(summary_lines, file.path(qa_dir, "AAD_only_paired_SDC1_summary.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(qa_dir, "sessionInfo.txt"), useBytes = TRUE)

print(correlation_summary)
cat(paste(summary_lines, collapse = "\n"), "\n")
if (!all(qa_checks)) quit(status = 1)
