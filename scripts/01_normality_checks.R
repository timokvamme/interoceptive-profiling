# =============================================================================
# Normality Checks for LPA Assumptions
# =============================================================================
#
# Output: analysis_output/01_output.txt
# Plots:  1_s_8_normality_histograms.png, 1_s_9_normality_qq_plots.png
#
# =============================================================================

library(tidyverse)

setwd("C:/code/projects/mi/analyses/soma")
output_dir <- "C:/code/projects/intero_mod/analysis_output"
plots_dir <- "C:/code/projects/intero_mod/plots/supplementary_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Redirect output to file
sink(file.path(output_dir, "01_output.txt"), split = TRUE)

# Load data
df <- read.csv("results/dfc_vviq_q_k.csv")
cat("N =", nrow(df), "\n\n")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("NORMALITY ASSESSMENT FOR LPA ASSUMPTIONS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# =============================================================================
# DESCRIPTIVE STATISTICS
# =============================================================================

cat("--- IAS (Interoceptive Accuracy) ---\n")
ias <- na.omit(df$ias)
cat("N:", length(ias), "\n")
cat("Mean:", round(mean(ias), 2), "| SD:", round(sd(ias), 2), "\n")
cat("Min:", round(min(ias), 2), "| Max:", round(max(ias), 2), "\n")

# Skewness and Kurtosis (manual calculation)
ias_skew <- (sum((ias - mean(ias))^3) / length(ias)) / (sd(ias)^3)
ias_kurt <- (sum((ias - mean(ias))^4) / length(ias)) / (sd(ias)^4) - 3

cat("Skewness:", round(ias_skew, 3))
if (abs(ias_skew) < 0.5) { cat(" (excellent: |skew| < 0.5)\n")
} else if (abs(ias_skew) < 1) { cat(" (acceptable: |skew| < 1)\n")
} else { cat(" (problematic: |skew| > 1)\n") }

cat("Excess Kurtosis:", round(ias_kurt, 3))
if (abs(ias_kurt) < 1) { cat(" (excellent: |kurt| < 1)\n")
} else if (abs(ias_kurt) < 2) { cat(" (acceptable: |kurt| < 2)\n")
} else { cat(" (problematic: |kurt| > 2)\n") }

# Shapiro-Wilk test
sw_ias <- shapiro.test(ias)
cat("Shapiro-Wilk W:", round(sw_ias$statistic, 4), "| p =", format.pval(sw_ias$p.value, digits = 3), "\n")

# Kolmogorov-Smirnov test
ks_ias <- ks.test(ias, "pnorm", mean(ias), sd(ias))
cat("Kolmogorov-Smirnov D:", round(ks_ias$statistic, 4), "| p =", format.pval(ks_ias$p.value, digits = 3), "\n\n")

# -----------------------------------------------------------------------------

cat("--- IATS (Interoceptive Attention) ---\n")
iats <- na.omit(df$iats)
cat("N:", length(iats), "\n")
cat("Mean:", round(mean(iats), 2), "| SD:", round(sd(iats), 2), "\n")
cat("Min:", round(min(iats), 2), "| Max:", round(max(iats), 2), "\n")

# Skewness and Kurtosis
iats_skew <- (sum((iats - mean(iats))^3) / length(iats)) / (sd(iats)^3)
iats_kurt <- (sum((iats - mean(iats))^4) / length(iats)) / (sd(iats)^4) - 3

cat("Skewness:", round(iats_skew, 3))
if (abs(iats_skew) < 0.5) { cat(" (excellent: |skew| < 0.5)\n")
} else if (abs(iats_skew) < 1) { cat(" (acceptable: |skew| < 1)\n")
} else { cat(" (problematic: |skew| > 1)\n") }

cat("Excess Kurtosis:", round(iats_kurt, 3))
if (abs(iats_kurt) < 1) { cat(" (excellent: |kurt| < 1)\n")
} else if (abs(iats_kurt) < 2) { cat(" (acceptable: |kurt| < 2)\n")
} else { cat(" (problematic: |kurt| > 2)\n") }

# Shapiro-Wilk test
sw_iats <- shapiro.test(iats)
cat("Shapiro-Wilk W:", round(sw_iats$statistic, 4), "| p =", format.pval(sw_iats$p.value, digits = 3), "\n")

# Kolmogorov-Smirnov test
ks_iats <- ks.test(iats, "pnorm", mean(iats), sd(iats))
cat("Kolmogorov-Smirnov D:", round(ks_iats$statistic, 4), "| p =", format.pval(ks_iats$p.value, digits = 3), "\n\n")

# =============================================================================
# INTERPRETATION
# =============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("INTERPRETATION FOR LPA\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Note: Statistical tests (Shapiro-Wilk, KS) are overly sensitive with\n")
cat("large samples (n=833). With large N, even trivial deviations from\n")
cat("normality will be statistically significant.\n\n")

cat("For LPA, focus on practical metrics:\n")
cat("  1. Skewness: |value| < 1 is acceptable\n")
cat("  2. Kurtosis: |value| < 2 is acceptable\n")
cat("  3. Visual inspection (Q-Q plots)\n\n")

cat("SUMMARY:\n")
cat("  IAS: skew =", round(ias_skew, 3), ", kurtosis =", round(ias_kurt, 3), "\n")
cat("  IATS: skew =", round(iats_skew, 3), ", kurtosis =", round(iats_kurt, 3), "\n\n")

if (abs(ias_skew) < 1 && abs(ias_kurt) < 2 && abs(iats_skew) < 1 && abs(iats_kurt) < 2) {
  cat("CONCLUSION: Both variables show acceptable normality for LPA.\n")
  cat("LPA assumptions are reasonably met.\n")
} else { cat("CONCLUSION: Some deviation from normality detected.\n")
  cat("Consider robust LPA variants or transformation.\n") }

# =============================================================================
# Q-Q PLOTS
# =============================================================================

png(file.path(plots_dir, "1_s_9_normality_qq_plots.png"), width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))

qqnorm(ias, main = "Q-Q Plot: IAS", pch = 20, col = rgb(0, 0, 0, 0.3))
qqline(ias, col = "red", lwd = 2)

qqnorm(iats, main = "Q-Q Plot: IATS", pch = 20, col = rgb(0, 0, 0, 0.3))
qqline(iats, col = "red", lwd = 2)

dev.off()

cat("\nQ-Q plots saved: plots/supplementary_plots/1_s_9_normality_qq_plots.png\n")

# =============================================================================
# HISTOGRAMS WITH NORMAL OVERLAY AND NORMALITY LABELS
# =============================================================================

png(file.path(plots_dir, "1_s_8_normality_histograms.png"), width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# Determine normality assessment
ias_normal <- abs(ias_skew) < 1 && abs(ias_kurt) < 2
iats_normal <- abs(iats_skew) < 1 && abs(iats_kurt) < 2

# IAS histogram
hist(ias, breaks = 30, probability = TRUE, main = "IAS Distribution",
     xlab = "IAS Score", col = "lightblue", border = "white")
curve(dnorm(x, mean = mean(ias), sd = sd(ias)), add = TRUE, col = "red", lwd = 2)

# Add normality label for IAS
ias_label <- paste0(
  "Skewness: ", round(ias_skew, 2), "\n",
  "Kurtosis: ", round(ias_kurt, 2), "\n",
  "SW W = ", round(sw_ias$statistic, 3), "\n",
  ifelse(ias_normal, "NORMAL", "NON-NORMAL")
)
ias_color <- ifelse(ias_normal, "darkgreen", "red")
legend("topright", legend = ias_label, bty = "n", cex = 0.9,
       text.col = ias_color, text.font = 2)

# IATS histogram
hist(iats, breaks = 30, probability = TRUE, main = "IATS Distribution",
     xlab = "IATS Score", col = "lightgreen", border = "white")
curve(dnorm(x, mean = mean(iats), sd = sd(iats)), add = TRUE, col = "red", lwd = 2)

# Add normality label for IATS
iats_label <- paste0(
  "Skewness: ", round(iats_skew, 2), "\n",
  "Kurtosis: ", round(iats_kurt, 2), "\n",
  "SW W = ", round(sw_iats$statistic, 3), "\n",
  ifelse(iats_normal, "NORMAL", "NON-NORMAL")
)
iats_color <- ifelse(iats_normal, "darkgreen", "red")
legend("topright", legend = iats_label, bty = "n", cex = 0.9,
       text.col = iats_color, text.font = 2)

dev.off()

cat("Histograms saved: plots/supplementary_plots/1_s_8_normality_histograms.png\n")
cat("\nDone!\n")

# Close output file
sink()
