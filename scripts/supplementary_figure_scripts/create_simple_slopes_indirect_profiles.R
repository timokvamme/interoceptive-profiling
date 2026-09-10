# =============================================================================
# Simple Slopes: Indirect Effects on Mental Health by Interoceptive Profile
# Shows indirect effects through Somatic and TAS for each cluster profile
# Indirect effect = (IAS -> Mediator) x (Mediator -> MH)
# =============================================================================

# Paths are relative to the repository root; run scripts from there.

library(tidyverse)
library(gridExtra)
library(grid)

# =============================================================================
# 1. LOAD DATA, CREATE PROFILES, DEFINE B-PATHS
# =============================================================================

dfc <- read.csv("dfc_interoception_profiling.csv")

dfc$ias_z <- scale(dfc$ias)[,1]
dfc$iats_z <- scale(dfc$iats)[,1]
dfc$tas_z <- scale(dfc$tas)[,1]
dfc$somatic_z <- scale(rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE))[,1]

# Create cluster profiles
set.seed(42)
cluster_data <- cbind(dfc$ias_z, dfc$iats_z)
km <- kmeans(cluster_data, centers = 3, nstart = 25)
dfc$cluster <- km$cluster

# Label clusters
cluster_means <- data.frame(
  cluster = 1:3,
  IAS_mean = km$centers[, 1],
  IATS_mean = km$centers[, 2]
)
cluster_means$profile <- NA
cluster_means$profile[which.min(cluster_means$IATS_mean)] <- "Efficient"
cluster_means$profile[which.max(cluster_means$IATS_mean)] <- "Hypervigilant"
cluster_means$profile[is.na(cluster_means$profile)] <- "Uncertain"

dfc$profile <- cluster_means$profile[match(dfc$cluster, cluster_means$cluster)]

# B-paths from SEM (Mediator -> Mental Health)
b_somatic_mh <- 0.61  # Somatic -> MH
b_tas_mh <- 0.31      # TAS -> MH

cat("B-paths from SEM:\n")
cat("  Somatic -> MH: \u03B2 =", b_somatic_mh, "\n")
cat("  TAS -> MH: \u03B2 =", b_tas_mh, "\n")

# =============================================================================
# 2. CALCULATE INDIRECT EFFECTS BY PROFILE
# =============================================================================

# Profile order and colors (consistent with figure_1_profiles.png)
profile_order <- c("Hypervigilant", "Uncertain", "Efficient")
profile_colors <- c("Hypervigilant" = "#e67e22", "Uncertain" = "#e74c3c", "Efficient" = "#2ecc71")

# Function to format significance
format_sig <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return(" (n.s.)")
}

ias_range <- seq(-2, 2, by = 0.1)

# Calculate indirect effects through Somatic by profile
indirect_somatic <- data.frame()
results_somatic <- data.frame()

for (prof in profile_order) {
  subset_data <- dfc[dfc$profile == prof, ]

  # Fit a-path model for this profile
  mod <- lm(somatic_z ~ ias_z, data = subset_data)
  a_path <- coef(mod)["ias_z"]
  p_val <- summary(mod)$coefficients["ias_z", "Pr(>|t|)"]

  # Indirect effect = a * b
  indirect_effect <- a_path * b_somatic_mh

  # Predicted indirect contribution across IAS range
  pred <- indirect_effect * ias_range

  indirect_somatic <- rbind(indirect_somatic, data.frame(
    IAS = ias_range,
    Indirect_Effect = pred,
    Profile = prof,
    stringsAsFactors = FALSE
  ))

  results_somatic <- rbind(results_somatic, data.frame(
    Profile = prof,
    a_path = round(a_path, 3),
    indirect = round(indirect_effect, 3),
    p = p_val,
    sig = format_sig(p_val),
    n = nrow(subset_data),
    label = paste0("\u03B2 = ", round(indirect_effect, 3), format_sig(p_val)),
    stringsAsFactors = FALSE
  ))
}

indirect_somatic$Profile <- factor(indirect_somatic$Profile, levels = profile_order)

# Calculate indirect effects through TAS by profile
indirect_tas <- data.frame()
results_tas <- data.frame()

for (prof in profile_order) {
  subset_data <- dfc[dfc$profile == prof, ]

  # Fit a-path model for this profile
  mod <- lm(tas_z ~ ias_z, data = subset_data)
  a_path <- coef(mod)["ias_z"]
  p_val <- summary(mod)$coefficients["ias_z", "Pr(>|t|)"]

  # Indirect effect = a * b
  indirect_effect <- a_path * b_tas_mh

  pred <- indirect_effect * ias_range

  indirect_tas <- rbind(indirect_tas, data.frame(
    IAS = ias_range,
    Indirect_Effect = pred,
    Profile = prof,
    stringsAsFactors = FALSE
  ))

  results_tas <- rbind(results_tas, data.frame(
    Profile = prof,
    a_path = round(a_path, 3),
    indirect = round(indirect_effect, 3),
    p = p_val,
    sig = format_sig(p_val),
    n = nrow(subset_data),
    label = paste0("\u03B2 = ", round(indirect_effect, 3), format_sig(p_val)),
    stringsAsFactors = FALSE
  ))
}

indirect_tas$Profile <- factor(indirect_tas$Profile, levels = profile_order)

cat("\n=== INDIRECT THROUGH SOMATIC: IAS -> Somatic -> MH by Profile ===\n")
print(results_somatic)

cat("\n=== INDIRECT THROUGH TAS: IAS -> TAS -> MH by Profile ===\n")
print(results_tas)

# =============================================================================
# 3. CREATE PLOTS
# =============================================================================

label_y_positions_som <- c(0.38, 0.30, 0.22)  # Top right for somatic
label_y_positions_tas <- c(0.38, 0.30, 0.22)  # Top right for TAS

# Panel A: Indirect through Somatic (MODERATED by profile)
p_somatic <- ggplot(indirect_somatic, aes(x = IAS, y = Indirect_Effect, color = Profile)) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  scale_color_manual(values = profile_colors, name = "Interoceptive\nProfile") +
  scale_x_continuous(breaks = -2:2, limits = c(-2.2, 2.2)) +
  scale_y_continuous(limits = c(-0.45, 0.45)) +
  # Add labels in top right
  annotate("text", x = 1.9, y = label_y_positions_som[1],
           label = results_somatic$label[results_somatic$Profile == "Hypervigilant"],
           color = profile_colors["Hypervigilant"], size = 3.2, fontface = "bold", hjust = 1) +
  annotate("text", x = 1.9, y = label_y_positions_som[2],
           label = results_somatic$label[results_somatic$Profile == "Uncertain"],
           color = profile_colors["Uncertain"], size = 3.2, fontface = "bold", hjust = 1) +
  annotate("text", x = 1.9, y = label_y_positions_som[3],
           label = results_somatic$label[results_somatic$Profile == "Efficient"],
           color = profile_colors["Efficient"], size = 3.2, fontface = "bold", hjust = 1) +
  labs(
    title = "A. Indirect via Somatic: Moderation by Interoceptive Profile",
    subtitle = paste0("Indirect = (IAS \u2192 Somatic) \u00d7 (Somatic \u2192 MH), b-path = ", b_somatic_mh),
    x = "Interoceptive Accuracy (IAS)",
    y = "Indirect Effect on\nMental Health Symptoms"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# Panel B: Indirect through TAS (NO MODERATION by profile)
p_tas <- ggplot(indirect_tas, aes(x = IAS, y = Indirect_Effect, color = Profile)) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  scale_color_manual(values = profile_colors, name = "Interoceptive\nProfile") +
  scale_x_continuous(breaks = -2:2, limits = c(-2.2, 2.2)) +
  scale_y_continuous(limits = c(-0.45, 0.45)) +
  # Add labels in top right
  annotate("text", x = 1.9, y = label_y_positions_tas[1],
           label = results_tas$label[results_tas$Profile == "Hypervigilant"],
           color = profile_colors["Hypervigilant"], size = 3.2, fontface = "bold", hjust = 1) +
  annotate("text", x = 1.9, y = label_y_positions_tas[2],
           label = results_tas$label[results_tas$Profile == "Uncertain"],
           color = profile_colors["Uncertain"], size = 3.2, fontface = "bold", hjust = 1) +
  annotate("text", x = 1.9, y = label_y_positions_tas[3],
           label = results_tas$label[results_tas$Profile == "Efficient"],
           color = profile_colors["Efficient"], size = 3.2, fontface = "bold", hjust = 1) +
  labs(
    title = "B. Indirect via TAS: No Moderation by Interoceptive Profile",
    subtitle = paste0("Indirect = (IAS \u2192 TAS) \u00d7 (TAS \u2192 MH), b-path = ", b_tas_mh),
    x = "Interoceptive Accuracy (IAS)",
    y = "Indirect Effect on\nMental Health Symptoms"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# =============================================================================
# 4. COMBINE AND SAVE
# =============================================================================

combined <- grid.arrange(
  p_somatic, p_tas,
  ncol = 2,
  top = textGrob("Indirect Effects on Mental Health: Profile Moderates Only the Somatic Pathway",
                 gp = gpar(fontsize = 14, fontface = "bold")),
  bottom = textGrob(sprintf("Note: Indirect effect = a-path \u00d7 b-path. Positive values indicate higher MH symptoms. ***p < .001, **p < .01, *p < .05. N = %d.", nrow(dfc)),
                    gp = gpar(fontsize = 9, col = "gray40"))
)

ggsave("plots/supplementary_plots/3_s_5_simple_slopes_indirect_profiles.png",
       combined, width = 13, height = 6, dpi = 300, bg = "white")

cat("\nFigure saved: figure_simple_slopes_indirect_profiles.png\n")

# =============================================================================
# 5. PRINT SUMMARY
# =============================================================================

cat("\n=== SUMMARY: INDIRECT EFFECTS BY PROFILE ===\n")
cat("\nVIA SOMATIC (IAS -> Somatic -> MH):\n")
for (i in 1:nrow(results_somatic)) {
  cat(sprintf("  %s (n=%d): a = %.3f, indirect %s (p = %.4f)\n",
              results_somatic$Profile[i], results_somatic$n[i],
              results_somatic$a_path[i], results_somatic$label[i], results_somatic$p[i]))
}

cat("\nVIA TAS (IAS -> TAS -> MH):\n")
for (i in 1:nrow(results_tas)) {
  cat(sprintf("  %s (n=%d): a = %.3f, indirect %s (p = %.4f)\n",
              results_tas$Profile[i], results_tas$n[i],
              results_tas$a_path[i], results_tas$label[i], results_tas$p[i]))
}

som_range <- max(results_somatic$indirect) - min(results_somatic$indirect)
tas_range <- max(results_tas$indirect) - min(results_tas$indirect)
cat(sprintf("\nIndirect effect range - Somatic: %.3f | TAS: %.3f\n", som_range, tas_range))

cat("\n=== INTERPRETATION ===\n")
cat("Negative indirect = IAS reduces MH symptoms through this pathway\n")
cat("Positive indirect = IAS increases MH symptoms through this pathway\n")
