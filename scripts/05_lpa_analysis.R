# =============================================================================
# LATENT PROFILE ANALYSIS (LPA)
# =============================================================================
#
# Output: analysis_output/05_output.txt
# Plots:  1_s_4_lpa_profiles.png
#
# PART 1: Model comparison (1-6 profiles) and k-means validation
# PART 2: Profile visualization and mental health outcomes
#
# =============================================================================

library(tidyverse)
library(tidyLPA)
library(mclust)
library(gridExtra)

setwd("C:/code/projects/interoceptive-profiling")
output_dir <- "C:/code/projects/interoceptive-profiling"
analysis_output_dir <- "C:/code/projects/interoceptive-profiling/analysis_output"
supplementary_data_dir <- "C:/code/projects/interoceptive-profiling/supplementary_data"
plots_dir <- "C:/code/projects/interoceptive-profiling/plots"
dir.create(analysis_output_dir, showWarnings = FALSE, recursive = TRUE)

# Redirect output to file
sink(file.path(analysis_output_dir, "05_output.txt"), split = TRUE)

# =============================================================================
# PART 1: LPA MODEL FITTING AND COMPARISON
# =============================================================================

df <- read.csv("dfc_vviq_q_k.csv")
cat("Loaded", nrow(df), "participants\n\n")

lpa_data <- df %>%
  select(ias, iats) %>%
  drop_na()

cat("Data for LPA:", nrow(lpa_data), "complete cases\n\n")

# Run LPA with different profile numbers
cat("Running LPA models for 1-6 profiles...\n\n")

lpa_results <- lpa_data %>%
  estimate_profiles(1:6,
                    variances = c("equal", "varying"),
                    covariances = c("zero", "varying"))

fit_stats <- get_fit(lpa_results)
write.csv(fit_stats, file.path(supplementary_data_dir, "lpa_fit_statistics.csv"), row.names = FALSE)

# Model comparison
model1_fits <- fit_stats %>%
  filter(Model == 1) %>%
  select(Classes, AIC, BIC, Entropy, BLRT_p)

cat("MODEL COMPARISON (BIC - lower is better):\n")
print(model1_fits)

best_bic <- model1_fits %>% filter(BIC == min(BIC, na.rm = TRUE))
cat("\nBest model by BIC:", best_bic$Classes, "profiles\n")

# Fit 3-profile solution for k-means comparison
lpa_3 <- lpa_data %>%
  estimate_profiles(3, variances = "equal", covariances = "zero")

lpa_3_results <- get_data(lpa_3)

# Compare with k-means
kmeans_data <- read.csv(file.path(supplementary_data_dir, "data_with_k_means_clusters.csv"))
lpa_3_results$kmeans_cluster <- kmeans_data$cluster

cross_tab <- table(LPA = lpa_3_results$Class, Kmeans = lpa_3_results$kmeans_cluster)
ari <- adjustedRandIndex(lpa_3_results$Class, lpa_3_results$kmeans_cluster)

cat("\nCross-tabulation LPA vs K-means:\n")
print(cross_tab)
cat("\nAdjusted Rand Index:", round(ari, 3), "\n")

# Save comparison
lpa_output <- lpa_3_results %>%
  select(ias, iats, Class, starts_with("CPROB")) %>%
  rename(lpa_profile = Class)
lpa_output$kmeans_cluster <- kmeans_data$cluster
write.csv(lpa_output, file.path(supplementary_data_dir, "lpa_vs_kmeans_comparison.csv"), row.names = FALSE)

# =============================================================================
# PART 2: PROFILE VISUALIZATION AND OUTCOMES
# =============================================================================

dfc <- df
dfc$lpa_profile <- lpa_3_results$Class

# Profile labels based on IAS/IATS patterns
profile_stats <- dfc %>%
  group_by(lpa_profile) %>%
  summarise(n = n(), IAS_mean = mean(ias), IATS_mean = mean(iats), .groups = "drop")

# Determine profile labels based on actual IAS/IATS means
# Profile with lowest IATS = Disengaged, highest IATS = Attuned, middle = Average
profile_order <- profile_stats %>%
  arrange(IATS_mean) %>%
  mutate(label = c("Disengaged", "Average", "Attuned"))

short_labels <- setNames(profile_order$label, as.character(profile_order$lpa_profile))
color_map_short <- c("Disengaged" = "#4DAF4A", "Average" = "#FF7F00", "Attuned" = "#E41A1C")

dfc$profile_short <- factor(short_labels[as.character(dfc$lpa_profile)],
                             levels = c("Disengaged", "Average", "Attuned"))

# Mental health outcomes
outcomes <- dfc %>%
  group_by(lpa_profile) %>%
  summarise(
    n = n(),
    TAS_mean = round(mean(tas, na.rm = TRUE), 1),
    PHQ9_mean = round(mean(phq9, na.rm = TRUE), 1),
    GAD7_mean = round(mean(gad7, na.rm = TRUE), 1),
    STAI_mean = round(mean(stai, na.rm = TRUE), 1),
    SSS8_mean = round(mean(sss8, na.rm = TRUE), 1),
    PCS_mean = round(mean(pcs, na.rm = TRUE), 1),
    .groups = "drop"
  )
outcomes$label <- short_labels[as.character(outcomes$lpa_profile)]

cat("\nLPA Profile Outcomes:\n")
print(outcomes)

# Scatter plot
p_scatter <- ggplot(dfc, aes(x = ias, y = iats, color = profile_short)) +
  geom_point(alpha = 0.5, size = 2) +
  stat_ellipse(level = 0.90, linewidth = 1.2) +
  scale_color_manual(values = color_map_short, na.value = "grey50") +
  labs(
    title = "A. LPA Interoceptive Profiles",
    x = "Interoceptive Accuracy (IAS)",
    y = "Interoceptive Attention (IATS)",
    color = "Profile"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# Bar plot function
create_bar_plot <- function(var, var_label, panel_letter, data, color_map) {
  summary_data <- data %>%
    group_by(profile_short) %>%
    summarize(
      mean = mean(get(var), na.rm = TRUE),
      se = sd(get(var), na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  ggplot(summary_data, aes(x = profile_short, y = mean, fill = profile_short)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_fill_manual(values = color_map) +
    labs(title = paste0(panel_letter, ". ", var_label), x = NULL, y = "Mean (+/- SE)") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
}

# Create bar plots
p_tas <- create_bar_plot("tas", "Alexithymia (TAS)", "B", dfc, color_map_short)
p_sss8 <- create_bar_plot("sss8", "Somatic (SSS-8)", "C", dfc, color_map_short)
p_pcs <- create_bar_plot("pcs", "Pain Catast. (PCS)", "D", dfc, color_map_short)
p_phq9 <- create_bar_plot("phq9", "Depression (PHQ-9)", "E", dfc, color_map_short)
p_gad7 <- create_bar_plot("gad7", "Anxiety (GAD-7)", "F", dfc, color_map_short)
p_stai <- create_bar_plot("stai", "Trait Anxiety (STAI)", "G", dfc, color_map_short)

# Combine
bar_grid <- arrangeGrob(p_tas, p_sss8, p_pcs, p_phq9, p_gad7, p_stai, ncol = 3)
combined_figure <- arrangeGrob(p_scatter, bar_grid, ncol = 1, heights = c(1.2, 1.5))

# Save figure (renamed to match numbering pattern)
ggsave(file.path(plots_dir, "supplementary_plots/1_s_4_lpa_profiles.png"), combined_figure,
       width = 11, height = 12, dpi = 300)

write.csv(outcomes, file.path(supplementary_data_dir, "lpa_profile_outcomes.csv"), row.names = FALSE)

# Summary
cat("\n================================================================================\n")
cat("LATENT PROFILE ANALYSIS - SUMMARY\n")
cat("================================================================================\n\n")
cat("Best model by BIC:", best_bic$Classes, "profiles\n")
cat("Adjusted Rand Index (vs k-means):", round(ari, 3), "\n\n")
cat("Profile sizes:\n")
print(table(dfc$lpa_profile))
cat("\n\nProfile outcomes:\n")
print(outcomes)

cat("\nLPA analysis complete. Output saved to:", output_dir, "\n")

# Close output file
sink()
