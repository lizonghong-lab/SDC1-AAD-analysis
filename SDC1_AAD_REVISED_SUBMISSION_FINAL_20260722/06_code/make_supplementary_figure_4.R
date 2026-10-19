options(stringsAsFactors = FALSE)

input <- file.path(getwd(), "SDC1_AAD_REANALYSIS_FINAL_20260721", "02_analysis_data", "paired_sdc1_subset.csv")
output_stem <- file.path(getwd(), "SDC1_AAD_REVISED_SUBMISSION_FINAL_20260722", "04_figures", "supplementary", "Supplementary_Figure_5_EV_vs_shed_SDC1")

d <- read.csv(input, check.names = FALSE)
d <- d[is.finite(d$ev_sdc1) & d$ev_sdc1 > 0 & is.finite(d$shed_sdc1) & d$shed_sdc1 > 0, ]
d$x <- log2(d$ev_sdc1)
d$y <- log2(d$shed_sdc1)

groups <- c("Healthy control", "AAD", "AMI", "Angina", "PE")
raw_groups <- c("Healthy control", "AAD", "AMI", "Angina", "PE")
palette <- c("#5B9C95", "#C44755", "#4A83B3", "#D88B2B", "#7F5AA6")
names(palette) <- raw_groups

draw_figure <- function() {
  par(mar = c(5.2, 5.4, 2.4, 1.2), family = "sans", las = 1)
  plot(d$x, d$y,
       type = "n",
       xlab = "EV-associated SDC1 (log2 pg/mL)",
       ylab = "Shed SDC1 (log2 pg/mL)",
       cex.lab = 1.25,
       cex.axis = 1.05,
       bty = "l")
  for (g in raw_groups) {
    idx <- d$group == g
    points(d$x[idx], d$y[idx], pch = 21, bg = palette[g], col = palette[g], cex = 1.15)
  }
  fit <- lm(y ~ x, data = d)
  abline(fit, col = "#7A7F87", lwd = 2.2)
  legend("topright", legend = groups, pt.bg = palette, col = palette, pch = 21,
         bty = "n", cex = 0.95, title = "Diagnostic group")
  text(par("usr")[1] + 0.03 * diff(par("usr")[1:2]),
       par("usr")[4] - 0.07 * diff(par("usr")[3:4]),
       labels = "Spearman rho = 0.20\n95% CI -0.02 to 0.42",
       adj = c(0, 1), cex = 1.05)
}

png(paste0(output_stem, ".png"), width = 2400, height = 1800, res = 300, bg = "white")
draw_figure()
dev.off()

cairo_pdf(paste0(output_stem, ".pdf"), width = 8, height = 6, family = "Arial")
draw_figure()
dev.off()

svg(paste0(output_stem, ".svg"), width = 8, height = 6, family = "Arial")
draw_figure()
dev.off()

tiff(paste0(output_stem, ".tiff"), width = 8, height = 6, units = "in", res = 600,
     compression = "lzw", bg = "white")
draw_figure()
dev.off()
