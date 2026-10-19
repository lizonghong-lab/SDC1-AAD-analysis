options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

root <- Sys.getenv("SDC1_PROJECT_ROOT", unset = normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]), "..", ".."), mustWork = FALSE))
if (!dir.exists(root)) root <- normalizePath(file.path(getwd()), mustWork = TRUE)
source_dir <- file.path(root, "04_results", "source_data")
out_dir <- file.path(root, "04_results", "figures", "response")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

palette <- c(
  "Single marker" = "#6F8F88",
  "116 dimensions" = "#D69032",
  "1,837 proteins" = "#B6404C",
  "Complete data" = "#2F6FA3",
  "Observed detection" = "#B6404C"
)

theme_response <- theme_classic(base_size = 8.2, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(size = 9.2, face = "bold", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 7.5, colour = "#555555", margin = margin(b = 6)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7.2),
    plot.tag = element_text(size = 10, face = "bold")
  )

curve <- read.csv(file.path(source_dir, "proteomics_power_false_negative_curve.csv"), check.names = FALSE)
curve <- subset(curve, contrast == "AAD versus AMI")
curve$scenario <- factor(
  curve$multiplicity_scenario,
  levels = c("Prespecified single candidate", "Effective-dimension family-wise control", "Protein-level family-wise control"),
  labels = c("Single marker", "116 dimensions", "1,837 proteins")
)

p1 <- ggplot(curve, aes(standardised_mean_difference, false_negative_probability, colour = scenario)) +
  geom_hline(yintercept = 0.20, linetype = 2, colour = "#888888", linewidth = 0.4) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = palette[c("Single marker", "116 dimensions", "1,837 proteins")]) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 1.15)) +
  labs(
    x = "True standardised mean difference",
    y = "False-negative probability",
    title = "Detectable-effect limits depend on multiplicity",
    subtitle = "AAD versus AMI; dashed line marks 20% false-negative probability"
  ) + theme_response

bh <- read.csv(file.path(source_dir, "proteomics_bh_simulation_power_v4.csv"), check.names = FALSE)
bh <- subset(bh, standardised_mean_difference <= 1.4 & grepl("^AAD versus AMI;", scenario))
bh$scenario_label <- ifelse(
  grepl("observed SDC1 detection counts", bh$scenario),
  "Observed detection",
  "Complete data"
)
bh$fn_lower <- 1 - bh$upper_95
bh$fn_upper <- 1 - bh$lower_95

p2 <- ggplot(bh, aes(standardised_mean_difference, simulated_false_negative_probability, colour = scenario_label, fill = scenario_label)) +
  geom_hline(yintercept = 0.20, linetype = 2, colour = "#888888", linewidth = 0.4) +
  geom_ribbon(aes(ymin = fn_lower, ymax = fn_upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = palette[c("Complete data", "Observed detection")]) +
  scale_fill_manual(values = palette[c("Complete data", "Observed detection")]) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "True standardised mean difference",
    y = "Simulated false-negative probability",
    title = "BH-adjusted discovery can miss moderate effects",
    subtitle = "1,837 tests; ribbons are binomial 95% intervals"
  ) + theme_response

figure <- (p1 | p2) + plot_annotation(tag_levels = "a")

save_response <- function(plot, stem, width_mm = 183, height_mm = 88, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::svg(file.path(out_dir, paste0(stem, ".svg")), width = width_in, height = height_in, family = "Arial")
  print(plot)
  dev.off()
  grDevices::cairo_pdf(file.path(out_dir, paste0(stem, ".pdf")), width = width_in, height = height_in, family = "Arial")
  print(plot)
  dev.off()
  grDevices::tiff(file.path(out_dir, paste0(stem, ".tiff")), width = width_in, height = height_in, units = "in", res = dpi, compression = "lzw", type = "cairo")
  print(plot)
  dev.off()
  grDevices::png(file.path(out_dir, paste0(stem, ".png")), width = width_in, height = height_in, units = "in", res = 300, type = "cairo")
  print(plot)
  dev.off()
}

save_response(figure, "Response_Figure_1_false_negative_limits")
cat("Generated response-only false-negative figure in", out_dir, "\n")
