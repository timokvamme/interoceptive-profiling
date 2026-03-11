# =============================================================================
# Basic Correlation: IAS vs IATS
# =============================================================================
#
# Output: analysis_output/02_output.txt
# Plots:  2_s_1_item_mh_correlations.png
#
# =============================================================================

setwd("C:/code/projects/interoceptive-profiling")

output_dir <- "C:/code/projects/interoceptive-profiling/analysis_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Redirect output to file
sink(file.path(output_dir, "02_output.txt"), split = TRUE)

# Load data
dfc <- read.csv('dfc_interoception_profiling.csv')

library(tidyverse)

# Directories
supplementary_data_dir <- "C:/code/projects/interoceptive-profiling/supplementary_data"
supplementary_plots_dir <- "C:/code/projects/interoceptive-profiling/plots/supplementary_plots"
dir.create(supplementary_plots_dir, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================================\n")
cat("INTEROCEPTION AND MENTAL HEALTH: CORRELATION ANALYSIS\n")
cat("=============================================================================\n\n")

cat("N =", nrow(dfc), "\n\n")

# =============================================================================
# PART 1: IAS vs IATS CORRELATION
# =============================================================================

cat("=============================================================================\n")
cat("PART 1: IAS vs IATS SCALE CORRELATION\n")
cat("=============================================================================\n\n")

cor_test <- cor.test(dfc$ias, dfc$iats, use = "complete.obs")

cat("IAS (Interoceptive Accuracy Scale): 21 items, sum score\n")
cat("  Mean =", round(mean(dfc$ias, na.rm = TRUE), 2),
    "| SD =", round(sd(dfc$ias, na.rm = TRUE), 2), "\n")
cat("  Range:", round(min(dfc$ias, na.rm = TRUE), 1), "-",
    round(max(dfc$ias, na.rm = TRUE), 1), "\n\n")

cat("IATS (Interoceptive Attention Scale): 21 items, sum score\n")
cat("  Mean =", round(mean(dfc$iats, na.rm = TRUE), 2),
    "| SD =", round(sd(dfc$iats, na.rm = TRUE), 2), "\n")
cat("  Range:", round(min(dfc$iats, na.rm = TRUE), 1), "-",
    round(max(dfc$iats, na.rm = TRUE), 1), "\n\n")

cat("IAS-IATS Correlation:\n")
cat("  r =", round(cor_test$estimate, 3), "\n")
cat("  p =", format.pval(cor_test$p.value, digits = 3), "\n")
cat("  95% CI:", round(cor_test$conf.int[1], 3), "to", round(cor_test$conf.int[2], 3), "\n")
cat("  n =", sum(complete.cases(dfc$ias, dfc$iats)), "\n\n")

cat("Interpretation: IAS and IATS are weakly negatively correlated,\n")
cat("suggesting they measure distinct constructs.\n\n")

# =============================================================================
# PART 2: MENTAL HEALTH COMPOSITE CONSTRUCTION
# =============================================================================

cat("=============================================================================\n")
cat("PART 2: MENTAL HEALTH COMPOSITE CONSTRUCTION\n")
cat("=============================================================================\n\n")

outcome_vars <- c("tas", "phq9", "gad7", "stai", "sss8", "pcs")
outcome_labels <- c("TAS (Alexithymia)", "PHQ-9 (Depression)", "GAD-7 (Anxiety)",
                    "STAI (Trait Anxiety)", "SSS-8 (Somatic Symptoms)", "PCS (Pain Catastrophizing)")

cat("The MH Composite is constructed as:\n")
cat("  MH_Composite = rowMeans(z-scored outcomes)\n\n")
cat("Components (all z-scored before averaging):\n")
for (i in seq_along(outcome_vars)) {
  v <- outcome_vars[i]
  cat(sprintf("  %d. %s: Mean = %.2f, SD = %.2f\n",
              i, outcome_labels[i], mean(dfc[[v]], na.rm = TRUE), sd(dfc[[v]], na.rm = TRUE)))
}

# Create composite
dfc$mh_composite <- rowMeans(scale(dfc[, outcome_vars]), na.rm = TRUE)

cat("\nMH Composite:\n")
cat("  Mean =", round(mean(dfc$mh_composite, na.rm = TRUE), 3), "(by construction ~0)\n")
cat("  SD =", round(sd(dfc$mh_composite, na.rm = TRUE), 3), "\n")
cat("  Range:", round(min(dfc$mh_composite, na.rm = TRUE), 2), "to",
    round(max(dfc$mh_composite, na.rm = TRUE), 2), "\n\n")

cat("Higher MH Composite = worse mental health\n\n")

# =============================================================================
# PART 3: SCALE TOTALS vs MH COMPOSITE
# =============================================================================

cat("=============================================================================\n")
cat("PART 3: SCALE TOTALS vs MH COMPOSITE\n")
cat("=============================================================================\n\n")

ias_cor <- cor.test(dfc$ias, dfc$mh_composite, use = "complete.obs")
iats_cor <- cor.test(dfc$iats, dfc$mh_composite, use = "complete.obs")

cat("IAS Total vs MH Composite:\n")
cat("  r =", round(ias_cor$estimate, 3), "\n")
cat("  p =", format.pval(ias_cor$p.value, digits = 3), "\n")
cat("  95% CI:", round(ias_cor$conf.int[1], 3), "to", round(ias_cor$conf.int[2], 3), "\n")
cat("  Interpretation: NEGATIVE correlation - higher accuracy = BETTER mental health\n\n")

cat("IATS Total vs MH Composite:\n")
cat("  r =", round(iats_cor$estimate, 3), "\n")
cat("  p =", format.pval(iats_cor$p.value, digits = 3), "\n")
cat("  95% CI:", round(iats_cor$conf.int[1], 3), "to", round(iats_cor$conf.int[2], 3), "\n")
cat("  Interpretation: POSITIVE correlation - higher attention = WORSE mental health\n\n")

# =============================================================================
# PART 4: ITEM-LEVEL CORRELATIONS WITH MH COMPOSITE
# =============================================================================

cat("=============================================================================\n")
cat("PART 4: ITEM-LEVEL CORRELATIONS WITH MH COMPOSITE\n")
cat("=============================================================================\n\n")

# Calculate correlations for all IAS items
ias_items <- paste0("ias_", 1:21)
iats_items <- paste0("iats_", 1:21)

cat("IAS Items (negative r = protective):\n")
cat(sprintf("%-8s %8s %8s %s\n", "Item", "r", "p", "Direction"))
cat(paste(rep("-", 60), collapse = ""), "\n")

ias_cors <- data.frame(item = ias_items, r = NA, p = NA)
for (i in 1:21) {
  item <- ias_items[i]
  ct <- cor.test(dfc[[item]], dfc$mh_composite, use = "complete.obs")
  ias_cors$r[i] <- ct$estimate
  ias_cors$p[i] <- ct$p.value
  direction <- ifelse(ct$estimate < 0, "Protective", "Risk")
  sig <- ifelse(ct$p.value < 0.001, "***", ifelse(ct$p.value < 0.01, "**", ifelse(ct$p.value < 0.05, "*", "")))
  cat(sprintf("%-8s %8.3f %8s %s\n", item, ct$estimate, paste0(format.pval(ct$p.value, digits = 2), sig), direction))
}

cat("\nIAS Summary:\n")
cat("  Items with protective effect (r < 0):", sum(ias_cors$r < 0), "/ 21\n")
cat("  Significant protective (p < .05):", sum(ias_cors$r < 0 & ias_cors$p < 0.05), "/ 21\n")
cat("  Mean r =", round(mean(ias_cors$r), 3), "\n\n")

cat("IATS Items (positive r = risk):\n")
cat(sprintf("%-8s %8s %8s %s\n", "Item", "r", "p", "Direction"))
cat(paste(rep("-", 60), collapse = ""), "\n")

iats_cors <- data.frame(item = iats_items, r = NA, p = NA)
for (i in 1:21) {
  item <- iats_items[i]
  ct <- cor.test(dfc[[item]], dfc$mh_composite, use = "complete.obs")
  iats_cors$r[i] <- ct$estimate
  iats_cors$p[i] <- ct$p.value
  direction <- ifelse(ct$estimate > 0, "Risk", "Protective")
  sig <- ifelse(ct$p.value < 0.001, "***", ifelse(ct$p.value < 0.01, "**", ifelse(ct$p.value < 0.05, "*", "")))
  cat(sprintf("%-8s %8.3f %8s %s\n", item, ct$estimate, paste0(format.pval(ct$p.value, digits = 2), sig), direction))
}

cat("\nIATS Summary:\n")
cat("  Items with risk effect (r > 0):", sum(iats_cors$r > 0), "/ 21\n")
cat("  Significant risk (p < .05):", sum(iats_cors$r > 0 & iats_cors$p < 0.05), "/ 21\n")
cat("  Mean r =", round(mean(iats_cors$r), 3), "\n\n")

# =============================================================================
# PART 5: SCALE CORRELATIONS WITH INDIVIDUAL OUTCOMES
# =============================================================================

cat("=============================================================================\n")
cat("PART 5: SCALE CORRELATIONS WITH INDIVIDUAL MH OUTCOMES\n")
cat("=============================================================================\n\n")

cat("IAS correlations with MH outcomes:\n")
cat(sprintf("%-25s %8s %8s\n", "Outcome", "r", "p"))
cat(paste(rep("-", 45), collapse = ""), "\n")
for (i in seq_along(outcome_vars)) {
  ct <- cor.test(dfc$ias, dfc[[outcome_vars[i]]], use = "complete.obs")
  cat(sprintf("%-25s %8.3f %8s\n", outcome_labels[i], ct$estimate, format.pval(ct$p.value, digits = 3)))
}

cat("\nIATS correlations with MH outcomes:\n")
cat(sprintf("%-25s %8s %8s\n", "Outcome", "r", "p"))
cat(paste(rep("-", 45), collapse = ""), "\n")
for (i in seq_along(outcome_vars)) {
  ct <- cor.test(dfc$iats, dfc[[outcome_vars[i]]], use = "complete.obs")
  cat(sprintf("%-25s %8.3f %8s\n", outcome_labels[i], ct$estimate, format.pval(ct$p.value, digits = 3)))
}

cat("\n")

# =============================================================================
# SUPPLEMENTARY FIGURE
# =============================================================================

cat("=============================================================================\n")
cat("GENERATING SUPPLEMENTARY FIGURE\n")
cat("=============================================================================\n\n")

# Load item data with correlations
item_data <- read.csv(file.path(supplementary_data_dir, "interoception_item_data.csv"))

# Extract item number for ordering
item_data <- item_data %>%
  mutate(
    item_num = as.numeric(gsub(".*_", "", item)),
    display_label = paste0(item_num, ". ", questionnaire_text)
  )

ias_total_cor <- cor(dfc$ias, dfc$mh_composite, use = "complete.obs")
iats_total_cor <- cor(dfc$iats, dfc$mh_composite, use = "complete.obs")

ias_total_cor <- cor(dfc$ias, dfc$mh_composite, use = "complete.obs")
iats_total_cor <- cor(dfc$iats, dfc$mh_composite, use = "complete.obs")

# Create scale total rows
scale_totals <- data.frame(
  item = c("ias_total", "iats_total"),
  item_label = c("IAS (Total)", "IATS (Total)"),
  scale = c("IAS", "IATS"),
  questionnaire_text = c("IAS (Total)", "IATS (Total)"),
  correlation_with_mh_composite = c(ias_total_cor, iats_total_cor),
  item_num = c(22, 22),  # Place at end of each section
  display_label = c("IAS (Total)", "IATS (Total)")
)

# Combine items and totals
plot_data <- item_data %>%
  select(item, item_label, scale, questionnaire_text,
         correlation_with_mh_composite, item_num, display_label) %>%
  bind_rows(scale_totals) %>%
  # Order: IAS items 1-21 then total, IATS items 1-21 then total
  arrange(desc(scale), item_num) %>%
  mutate(
    y_order = row_number(),
    # Color based on correlation direction
    fill_color = ifelse(correlation_with_mh_composite < 0, "protective", "risk"),
    cor_label = sprintf("%.2f", correlation_with_mh_composite)
  )

# Create the visualization
p <- ggplot(plot_data, aes(x = 1, y = reorder(display_label, -y_order))) +
  geom_tile(aes(fill = correlation_with_mh_composite), width = 0.9, height = 0.9) +
  geom_text(aes(label = cor_label), color = "black", size = 3, fontface = "bold") +
  scale_fill_gradient2(
    low = "#4DAF4A",      # Green for negative (protective)
    mid = "white",
    high = "#E41A1C",     # Red for positive (risk)
    midpoint = 0,
    limits = c(-0.5, 0.5),
    name = "Correlation\nwith MH"
  ) +
  facet_grid(scale ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title = "Interoception Items: Correlations with Mental Health Composite",
    subtitle = "Green = protective (better MH) | Red = risk (worse MH)",
    x = "MH Composite",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7, hjust = 1),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
    strip.placement = "outside",
    panel.grid = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
    legend.position = "right"
  )

# Save the plot
output_file <- file.path(supplementary_plots_dir, "2_s_1_item_mh_correlations.png")
ggsave(output_file, p, width = 12, height = 16, dpi = 300, bg = "white")
cat("\nSupplementary figure saved:", output_file, "\n")

cat("\nDone.\n")

# Close output file
sink()
