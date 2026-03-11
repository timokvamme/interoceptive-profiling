# =============================================================================
# COMPREHENSIVE MODERATED MEDIATION FIGURE (Supplementary Figure 19)
# =============================================================================
#
# Panels:
#   A. Conditional indirect effects at -1SD, Mean, +1SD of IATS
#   B. Profile-specific indirect effects from multi-group SEM
#   C. Forest plot of multi-group path coefficients (NEW - replaces text panel)
#
# Output: plots/supplementary_plots/3_s_7_comprehensive_moderated_mediation.png
#
# =============================================================================

setwd("C:/code/projects/interoceptive-profiling")

library(dplyr)
library(lavaan)
library(ggplot2)
library(gridExtra)
library(grid)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dfc <- read.csv("C:/code/projects/interoceptive-profiling/dfc_interoception_profiling.csv")

# Save full-scale totals for clustering (consistent with 04_cluster_analysis.R)
ias_full <- dfc$ias
iats_full <- dfc$iats

# Compute composite outcomes
dfc$mh <- rowMeans(cbind(scale(dfc$phq9), scale(dfc$gad7), scale(dfc$stai)), na.rm = TRUE)
dfc$somatic <- rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE)

# Standardize (full-scale IAS/IATS for SEM, consistent with 12_output.txt)
dfc$ias_z <- scale(dfc$ias)[,1]
dfc$iats_z <- scale(dfc$iats)[,1]
dfc$tas_z <- scale(dfc$tas)[,1]
dfc$somatic_z <- scale(dfc$somatic)[,1]
dfc$mh_z <- scale(dfc$mh)[,1]
dfc$ias_x_iats <- dfc$ias_z * dfc$iats_z

# Assign profiles via clustering on FULL-SCALE totals
# (matching 04_cluster_analysis.R: seed 123, nstart 50)
set.seed(123)
km <- kmeans(scale(cbind(ias_full, iats_full)), centers = 3, nstart = 50)
dfc$cluster <- km$cluster

# Label clusters using full-scale z-scores
dfc$ias_full_z <- scale(ias_full)[,1]
dfc$iats_full_z <- scale(iats_full)[,1]
cluster_means <- aggregate(cbind(ias_full_z, iats_full_z) ~ cluster, data = dfc, mean)
cluster_means$profile <- NA
cluster_means$profile[which.max(cluster_means$ias_full_z + cluster_means$iats_full_z)] <- "Hypervigilant"
cluster_means$profile[which.min(cluster_means$ias_full_z)] <- "Uncertain"
cluster_means$profile[is.na(cluster_means$profile)] <- "Efficient"
dfc$profile <- cluster_means$profile[match(dfc$cluster, cluster_means$cluster)]
dfc$ias_full_z <- NULL
dfc$iats_full_z <- NULL

cat("Sample sizes by profile:\n")
print(table(dfc$profile))

# Profile-typical IATS values
profile_values <- aggregate(iats_z ~ profile, data = dfc, mean)

# =============================================================================
# 2. PANEL A: CONDITIONAL INDIRECT EFFECTS AT -1SD, MEAN, +1SD IATS
# =============================================================================

cat("\nFitting SEM for Panel A (conditional indirect effects)...\n")

model_interaction <- '
  tas_z ~ a1*ias_z + a2*iats_z + a3*ias_x_iats
  somatic_z ~ b1*ias_z + b2*iats_z + b3*ias_x_iats
  mh_z ~ c1*tas_z + c2*somatic_z + d1*ias_z + d2*iats_z + d3*ias_x_iats

  # Conditional indirect at -1SD IATS
  cond_tas_low := (a1 + a3*(-1)) * c1
  cond_som_low := (b1 + b3*(-1)) * c2

  # Conditional indirect at mean IATS (=0)
  cond_tas_mean := a1 * c1
  cond_som_mean := b1 * c2

  # Conditional indirect at +1SD IATS
  cond_tas_high := (a1 + a3*(1)) * c1
  cond_som_high := (b1 + b3*(1)) * c2
'

fit_int <- sem(model_interaction, data = dfc)
params_int <- parameterEstimates(fit_int, ci = TRUE)

# Extract conditional indirect effects
cond_labels <- c("cond_tas_low", "cond_som_low",
                 "cond_tas_mean", "cond_som_mean",
                 "cond_tas_high", "cond_som_high")
cond_effects <- params_int[params_int$label %in% cond_labels,
                           c("label", "est", "ci.lower", "ci.upper", "pvalue")]

# Build Panel A data
panel_a_data <- data.frame(
  IATS_Level = factor(rep(c("-1 SD\n(Low)", "Mean", "+1 SD\n(High)"), each = 2),
                      levels = c("-1 SD\n(Low)", "Mean", "+1 SD\n(High)")),
  Pathway = rep(c("Somatic Pathway\n(IAS\u2192Som\u2192MH)", "TAS Pathway\n(IAS\u2192TAS\u2192MH)"), 3),
  Est = c(
    cond_effects$est[cond_effects$label == "cond_som_low"],
    cond_effects$est[cond_effects$label == "cond_tas_low"],
    cond_effects$est[cond_effects$label == "cond_som_mean"],
    cond_effects$est[cond_effects$label == "cond_tas_mean"],
    cond_effects$est[cond_effects$label == "cond_som_high"],
    cond_effects$est[cond_effects$label == "cond_tas_high"]
  ),
  CI_low = c(
    cond_effects$ci.lower[cond_effects$label == "cond_som_low"],
    cond_effects$ci.lower[cond_effects$label == "cond_tas_low"],
    cond_effects$ci.lower[cond_effects$label == "cond_som_mean"],
    cond_effects$ci.lower[cond_effects$label == "cond_tas_mean"],
    cond_effects$ci.lower[cond_effects$label == "cond_som_high"],
    cond_effects$ci.lower[cond_effects$label == "cond_tas_high"]
  ),
  CI_high = c(
    cond_effects$ci.upper[cond_effects$label == "cond_som_low"],
    cond_effects$ci.upper[cond_effects$label == "cond_tas_low"],
    cond_effects$ci.upper[cond_effects$label == "cond_som_mean"],
    cond_effects$ci.upper[cond_effects$label == "cond_tas_mean"],
    cond_effects$ci.upper[cond_effects$label == "cond_som_high"],
    cond_effects$ci.upper[cond_effects$label == "cond_tas_high"]
  ),
  p = c(
    cond_effects$pvalue[cond_effects$label == "cond_som_low"],
    cond_effects$pvalue[cond_effects$label == "cond_tas_low"],
    cond_effects$pvalue[cond_effects$label == "cond_som_mean"],
    cond_effects$pvalue[cond_effects$label == "cond_tas_mean"],
    cond_effects$pvalue[cond_effects$label == "cond_som_high"],
    cond_effects$pvalue[cond_effects$label == "cond_tas_high"]
  )
)
panel_a_data$Significant <- ifelse(panel_a_data$p < 0.05, "Yes", "No")

p_a <- ggplot(panel_a_data, aes(x = IATS_Level, y = Est, fill = Pathway, alpha = Significant)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(width = 0.7),
                width = 0.2, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Somatic Pathway\n(IAS\u2192Som\u2192MH)" = "#7C3AED",
                                "TAS Pathway\n(IAS\u2192TAS\u2192MH)" = "#F59E0B")) +
  scale_alpha_manual(values = c("Yes" = 1.0, "No" = 0.4), name = "Significant") +
  labs(
    title = "A. Conditional Indirect Effects of IAS on Mental Health",
    subtitle = "How accuracy's protective effect changes with attention level",
    x = "IATS Level (Interoceptive Attention)",
    y = "Indirect Effect on MH"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "gray40"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2))

# =============================================================================
# 3. PANEL B: PROFILE-SPECIFIC INDIRECT EFFECTS (from multi-group SEM)
# =============================================================================

cat("Fitting multi-group SEM for Panels B and C...\n")

# Multi-group model (no interaction term - same as 12c)
model_mg <- '
  tas_z ~ ias_z + iats_z
  somatic_z ~ ias_z + iats_z
  mh_z ~ tas_z + somatic_z + ias_z + iats_z
'

fit_configural <- sem(model_mg, data = dfc, group = "profile")
params_mg <- parameterEstimates(fit_configural, ci = TRUE, standardized = TRUE)

# Compute indirect effects by profile
groups <- c("Efficient", "Hypervigilant", "Uncertain")
indirect_data <- data.frame()

for (g in 1:3) {
  gp <- params_mg[params_mg$group == g & params_mg$op == "~", ]

  a1 <- gp$est[gp$lhs == "tas_z" & gp$rhs == "ias_z"]
  b1 <- gp$est[gp$lhs == "somatic_z" & gp$rhs == "ias_z"]
  c1 <- gp$est[gp$lhs == "mh_z" & gp$rhs == "tas_z"]
  c2 <- gp$est[gp$lhs == "mh_z" & gp$rhs == "somatic_z"]

  # Standard errors for products (delta method approximation)
  a1_se <- gp$se[gp$lhs == "tas_z" & gp$rhs == "ias_z"]
  b1_se <- gp$se[gp$lhs == "somatic_z" & gp$rhs == "ias_z"]
  c1_se <- gp$se[gp$lhs == "mh_z" & gp$rhs == "tas_z"]
  c2_se <- gp$se[gp$lhs == "mh_z" & gp$rhs == "somatic_z"]

  ind_tas <- a1 * c1
  ind_som <- b1 * c2

  # Approximate SE for product using delta method: SE(ab) ~ sqrt(a^2*SE_b^2 + b^2*SE_a^2)
  se_ind_tas <- sqrt(a1^2 * c1_se^2 + c1^2 * a1_se^2)
  se_ind_som <- sqrt(b1^2 * c2_se^2 + c2^2 * b1_se^2)

  indirect_data <- rbind(indirect_data, data.frame(
    Profile = groups[g],
    Pathway = c("IAS\u2192Somatic\u2192MH", "IAS\u2192TAS\u2192MH"),
    Est = c(ind_som, ind_tas),
    CI_low = c(ind_som - 1.96 * se_ind_som, ind_tas - 1.96 * se_ind_tas),
    CI_high = c(ind_som + 1.96 * se_ind_som, ind_tas + 1.96 * se_ind_tas)
  ))
}

indirect_data$Profile <- factor(indirect_data$Profile,
                                levels = c("Efficient", "Hypervigilant", "Uncertain"))

p_b <- ggplot(indirect_data, aes(x = Profile, y = Est, fill = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(width = 0.7),
                width = 0.2, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("IAS\u2192Somatic\u2192MH" = "#7C3AED",
                                "IAS\u2192TAS\u2192MH" = "#F59E0B")) +
  labs(
    title = "B. Profile-Specific Indirect Effects",
    subtitle = "Multi-group SEM results (no interaction term)",
    x = "Interoceptive Profile",
    y = "Indirect Effect on MH\n(Product of Path Coefficients)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "gray40"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# 4. PANEL C: FOREST PLOT OF MULTI-GROUP PATH COEFFICIENTS (NEW)
# =============================================================================

cat("Creating forest plot for Panel C...\n")

# Get standardized estimates with CIs
std_solution <- standardizedSolution(fit_configural)

# Extract the 4 key mediation paths for each group
key_paths <- data.frame()
path_specs <- list(
  list(lhs = "tas_z", rhs = "ias_z", label = "IAS \u2192 TAS"),
  list(lhs = "somatic_z", rhs = "ias_z", label = "IAS \u2192 Somatic"),
  list(lhs = "mh_z", rhs = "tas_z", label = "TAS \u2192 MH"),
  list(lhs = "mh_z", rhs = "somatic_z", label = "Somatic \u2192 MH")
)

# Also get the Wald test p-values for annotation
wald_p <- c(0.8740, 0.1768, 0.3314, 0.1456)

for (i in seq_along(path_specs)) {
  spec <- path_specs[[i]]
  for (g in 1:3) {
    row <- std_solution[std_solution$lhs == spec$lhs &
                        std_solution$rhs == spec$rhs &
                        std_solution$group == g &
                        std_solution$op == "~", ]
    if (nrow(row) > 0) {
      key_paths <- rbind(key_paths, data.frame(
        Path = spec$label,
        Profile = groups[g],
        Est = row$est.std[1],
        CI_low = row$ci.lower[1],
        CI_high = row$ci.upper[1],
        p = row$pvalue[1],
        Wald_p = wald_p[i]
      ))
    }
  }
}

# Order paths and profiles
key_paths$Path <- factor(key_paths$Path,
                         levels = rev(c("IAS \u2192 TAS", "IAS \u2192 Somatic",
                                        "TAS \u2192 MH", "Somatic \u2192 MH")))
key_paths$Profile <- factor(key_paths$Profile,
                            levels = c("Efficient", "Hypervigilant", "Uncertain"))

# Add Wald annotation labels
wald_labels <- data.frame(
  Path = factor(c("IAS \u2192 TAS", "IAS \u2192 Somatic",
                   "TAS \u2192 MH", "Somatic \u2192 MH"),
                levels = rev(c("IAS \u2192 TAS", "IAS \u2192 Somatic",
                               "TAS \u2192 MH", "Somatic \u2192 MH"))),
  label = c("Wald p = .874", "Wald p = .177", "Wald p = .331", "Wald p = .146")
)

profile_colors <- c("Efficient" = "#22C55E", "Hypervigilant" = "#F59E0B", "Uncertain" = "#EF4444")

p_c <- ggplot(key_paths, aes(x = Est, y = Path, color = Profile, shape = Profile)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high),
                  position = position_dodge(width = 0.6),
                  size = 0.6, linewidth = 0.7) +
  # Highlight the IAS->Somatic row
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 2.5, ymax = 3.5,
           fill = "#7C3AED", alpha = 0.08) +
  # Re-add points on top of highlight
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high),
                  position = position_dodge(width = 0.6),
                  size = 0.6, linewidth = 0.7) +
  # Add Wald test annotations on the right
  geom_text(data = wald_labels,
            aes(x = 0.85, y = Path, label = label),
            inherit.aes = FALSE,
            size = 2.8, color = "gray30", hjust = 0) +
  scale_color_manual(values = profile_colors) +
  scale_shape_manual(values = c("Efficient" = 16, "Hypervigilant" = 17, "Uncertain" = 15)) +
  scale_x_continuous(limits = c(-0.5, 1.1), breaks = seq(-0.4, 0.8, 0.2)) +
  labs(
    title = "C. Multi-Group Path Coefficients Across Profiles",
    subtitle = expression(paste("Omnibus invariance test: ", chi^2, "(16) = 16.40, ",
                                italic("p"), " = .425")),
    x = "Standardized Coefficient (\u03b2)",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "gray40"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 9, face = "bold")
  )

# =============================================================================
# 5. COMBINE AND SAVE
# =============================================================================

cat("Combining panels...\n")

# Top row: A and B side by side
top_row <- arrangeGrob(p_a, p_b, ncol = 2, widths = c(1, 1))

# Bottom row: C (full width)
combined <- arrangeGrob(
  top_row,
  p_c,
  nrow = 2,
  heights = c(1, 0.8),
  top = NULL
)

ggsave("plots/supplementary_plots/3_s_7_comprehensive_moderated_mediation.png",
       combined, width = 12, height = 10, dpi = 300, bg = "white")

cat("\nFigure saved: plots/supplementary_plots/3_s_7_comprehensive_moderated_mediation.png\n")
cat("DONE\n")
