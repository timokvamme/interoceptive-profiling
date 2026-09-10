# =============================================================================
# Simple Slopes: Indirect Effects on Mental Health (Continuous IATS Moderator)
# Shows indirect effects through Somatic and TAS at Low/Mean/High IATS
# Indirect effect = (IAS -> Mediator) x (Mediator -> MH)
# =============================================================================

# Paths are relative to the repository root; run scripts from there.

library(tidyverse)
library(gridExtra)
library(grid)

# =============================================================================
# 1. LOAD DATA AND DEFINE B-PATHS FROM SEM
# =============================================================================

dfc <- read.csv("dfc_interoception_profiling.csv")

dfc$ias_z <- scale(dfc$ias)[,1]
dfc$iats_z <- scale(dfc$iats)[,1]
dfc$tas_z <- scale(dfc$tas)[,1]
dfc$somatic_z <- scale(rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE))[,1]

# B-paths from SEM (Mediator -> Mental Health)
# From sem_diagram_results.txt:
b_somatic_mh <- 0.61  # Somatic -> MH
b_tas_mh <- 0.31      # TAS -> MH

cat("B-paths from SEM:\n")
cat("  Somatic -> MH: \u03B2 =", b_somatic_mh, "\n")
cat("  TAS -> MH: \u03B2 =", b_tas_mh, "\n")

# Create interaction term
dfc$ias_x_iats <- dfc$ias_z * dfc$iats_z

# Fit a-path models (IAS -> Mediators)
mod_somatic <- lm(somatic_z ~ ias_z + iats_z + ias_x_iats, data = dfc)
mod_tas <- lm(tas_z ~ ias_z + iats_z + ias_x_iats, data = dfc)

coef_som <- coef(mod_somatic)
coef_tas <- coef(mod_tas)

# =============================================================================
# 2. CALCULATE INDIRECT EFFECTS AT DIFFERENT IATS LEVELS
# =============================================================================

ias_range <- seq(-2, 2, by = 0.1)

# IATS levels and styling
iats_levels <- c(-1, 0, 1)
iats_labels <- c("Low IATS (-1 SD)", "Mean IATS", "High IATS (+1 SD)")
iats_colors <- c("Low IATS (-1 SD)" = "#3498db", "Mean IATS" = "#7f8c8d", "High IATS (+1 SD)" = "#e74c3c")

# Function to format significance
format_sig <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return(" (n.s.)")
}

# Calculate indirect effects through Somatic
indirect_somatic <- data.frame()
results_somatic <- data.frame()

for (i in 1:length(iats_levels)) {
  iats_val <- iats_levels[i]

  # A-path: simple slope of IAS -> Somatic at this IATS level
  a_path <- coef_som["ias_z"] + coef_som["ias_x_iats"] * iats_val

  # Indirect effect = a * b
  indirect_effect <- a_path * b_somatic_mh

  # Test a-path significance (indirect significance depends on a-path)
  vcov_mat <- vcov(mod_somatic)
  se_a <- sqrt(vcov_mat["ias_z", "ias_z"] +
               iats_val^2 * vcov_mat["ias_x_iats", "ias_x_iats"] +
               2 * iats_val * vcov_mat["ias_z", "ias_x_iats"])
  t_val <- a_path / se_a
  p_val <- 2 * pt(-abs(t_val), df = nrow(dfc) - 4)

  # For visualization: show effect on MH as function of IAS
  # At fixed IATS: indirect contribution = indirect_effect * IAS
  pred <- indirect_effect * ias_range

  indirect_somatic <- rbind(indirect_somatic, data.frame(
    IAS = ias_range,
    Indirect_Effect = pred,
    IATS_level = iats_labels[i],
    stringsAsFactors = FALSE
  ))

  results_somatic <- rbind(results_somatic, data.frame(
    IATS_level = iats_labels[i],
    a_path = round(a_path, 3),
    indirect = round(indirect_effect, 3),
    p = p_val,
    sig = format_sig(p_val),
    label = paste0("\u03B2 = ", round(indirect_effect, 3), format_sig(p_val)),
    stringsAsFactors = FALSE
  ))
}

indirect_somatic$IATS_level <- factor(indirect_somatic$IATS_level, levels = iats_labels)

# Calculate indirect effects through TAS
indirect_tas <- data.frame()
results_tas <- data.frame()

for (i in 1:length(iats_levels)) {
  iats_val <- iats_levels[i]

  # A-path: simple slope of IAS -> TAS at this IATS level
  a_path <- coef_tas["ias_z"] + coef_tas["ias_x_iats"] * iats_val

  # Indirect effect = a * b
  indirect_effect <- a_path * b_tas_mh

  # Test a-path significance
  vcov_mat <- vcov(mod_tas)
  se_a <- sqrt(vcov_mat["ias_z", "ias_z"] +
               iats_val^2 * vcov_mat["ias_x_iats", "ias_x_iats"] +
               2 * iats_val * vcov_mat["ias_z", "ias_x_iats"])
  t_val <- a_path / se_a
  p_val <- 2 * pt(-abs(t_val), df = nrow(dfc) - 4)

  pred <- indirect_effect * ias_range

  indirect_tas <- rbind(indirect_tas, data.frame(
    IAS = ias_range,
    Indirect_Effect = pred,
    IATS_level = iats_labels[i],
    stringsAsFactors = FALSE
  ))

  results_tas <- rbind(results_tas, data.frame(
    IATS_level = iats_labels[i],
    a_path = round(a_path, 3),
    indirect = round(indirect_effect, 3),
    p = p_val,
    sig = format_sig(p_val),
    label = paste0("\u03B2 = ", round(indirect_effect, 3), format_sig(p_val)),
    stringsAsFactors = FALSE
  ))
}

indirect_tas$IATS_level <- factor(indirect_tas$IATS_level, levels = iats_labels)

cat("\n=== INDIRECT THROUGH SOMATIC: IAS -> Somatic -> MH ===\n")
print(results_somatic)

cat("\n=== INDIRECT THROUGH TAS: IAS -> TAS -> MH ===\n")
print(results_tas)

# =============================================================================
# 3. CREATE PLOTS
# =============================================================================

label_y_positions <- c(0.55, 0.45, 0.35)

# Panel A: Indirect through Somatic (MODERATED)
p_somatic <- ggplot(indirect_somatic, aes(x = IAS, y = Indirect_Effect, color = IATS_level)) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  scale_color_manual(values = iats_colors, name = "Interoceptive\nAttention (IATS)") +
  scale_x_continuous(breaks = -2:2, limits = c(-2.2, 2.2)) +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  # Add labels
  annotate("text", x = 1.9, y = label_y_positions[1],
           label = results_somatic$label[results_somatic$IATS_level == "High IATS (+1 SD)"],
           color = iats_colors["High IATS (+1 SD)"], size = 3.2, fontface = "bold", hjust = 1) +
  annotate("text", x = 1.9, y = label_y_positions[2],
           label = results_somatic$label[results_somatic$IATS_level == "Mean IATS"],
           color = iats_colors["Mean IATS"], size = 3.2, fontface = "bold", hjust = 1) +
  annotate("text", x = 1.9, y = label_y_positions[3],
           label = results_somatic$label[results_somatic$IATS_level == "Low IATS (-1 SD)"],
           color = iats_colors["Low IATS (-1 SD)"], size = 3.2, fontface = "bold", hjust = 1) +
  labs(
    title = "A. Indirect via Somatic: Moderated by IATS",
    subtitle = paste0("Indirect = (IAS \u2192 Somatic) \u00d7 (Somatic \u2192 MH), b-path = ", b_somatic_mh),
    x = "Interoceptive Accuracy (IAS)",
    y = "Indirect Effect on\nMental Health Symptoms"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = c(0.18, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# Panel B: Indirect through TAS (UNMODERATED)
p_tas <- ggplot(indirect_tas, aes(x = IAS, y = Indirect_Effect, color = IATS_level)) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  scale_color_manual(values = iats_colors, name = "Interoceptive\nAttention (IATS)") +
  scale_x_continuous(breaks = -2:2, limits = c(-2.2, 2.2)) +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  # Add labels
  annotate("text", x = 1.9, y = label_y_positions[1],
           label = results_tas$label[results_tas$IATS_level == "High IATS (+1 SD)"],
           color = iats_colors["High IATS (+1 SD)"], size = 3.2, fontface = "bold", hjust = 1) +
  annotate("text", x = 1.9, y = label_y_positions[2],
           label = results_tas$label[results_tas$IATS_level == "Mean IATS"],
           color = iats_colors["Mean IATS"], size = 3.2, fontface = "bold", hjust = 1) +
  annotate("text", x = 1.9, y = label_y_positions[3],
           label = results_tas$label[results_tas$IATS_level == "Low IATS (-1 SD)"],
           color = iats_colors["Low IATS (-1 SD)"], size = 3.2, fontface = "bold", hjust = 1) +
  labs(
    title = "B. Indirect via TAS: No Moderation by IATS",
    subtitle = paste0("Indirect = (IAS \u2192 TAS) \u00d7 (TAS \u2192 MH), b-path = ", b_tas_mh),
    x = "Interoceptive Accuracy (IAS)",
    y = "Indirect Effect on\nMental Health Symptoms"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = c(0.18, 0.85),
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
  top = textGrob("Indirect Effects on Mental Health: IATS Moderates Only the Somatic Pathway",
                 gp = gpar(fontsize = 14, fontface = "bold")),
  bottom = textGrob(sprintf("Note: Indirect effect = a-path \u00d7 b-path. Positive values indicate higher MH symptoms. ***p < .001, **p < .01, *p < .05. N = %d.", nrow(dfc)),
                    gp = gpar(fontsize = 9, col = "gray40"))
)

ggsave("plots/supplementary_plots/3_s_4_simple_slopes_indirect_continuous.png",
       combined, width = 13, height = 6, dpi = 300, bg = "white")

cat("\nFigure saved: figure_simple_slopes_indirect_continuous.png\n")

# =============================================================================
# 5. PRINT SUMMARY
# =============================================================================

cat("\n=== SUMMARY: INDIRECT EFFECTS ===\n")
cat("\nVIA SOMATIC (IAS -> Somatic -> MH):\n")
for (i in 1:nrow(results_somatic)) {
  cat(sprintf("  %s: a = %.3f, indirect %s (p = %.4f)\n",
              results_somatic$IATS_level[i], results_somatic$a_path[i],
              results_somatic$label[i], results_somatic$p[i]))
}

cat("\nVIA TAS (IAS -> TAS -> MH):\n")
for (i in 1:nrow(results_tas)) {
  cat(sprintf("  %s: a = %.3f, indirect %s (p = %.4f)\n",
              results_tas$IATS_level[i], results_tas$a_path[i],
              results_tas$label[i], results_tas$p[i]))
}

som_range <- max(results_somatic$indirect) - min(results_somatic$indirect)
tas_range <- max(results_tas$indirect) - min(results_tas$indirect)
cat(sprintf("\nIndirect effect range - Somatic: %.3f | TAS: %.3f\n", som_range, tas_range))
