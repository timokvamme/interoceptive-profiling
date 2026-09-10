# =============================================================================
# Simple Slopes: Direct Effects on Mediators (Continuous IATS Moderator)
# Shows IAS effect on Somatic and TAS at Low/Mean/High IATS
# =============================================================================

# Paths are relative to the repository root; run scripts from there.

library(tidyverse)
library(gridExtra)
library(grid)

# =============================================================================
# 1. LOAD DATA AND FIT MODELS
# =============================================================================

dfc <- read.csv("dfc_interoception_profiling.csv")

dfc$ias_z <- scale(dfc$ias)[,1]
dfc$iats_z <- scale(dfc$iats)[,1]
dfc$tas_z <- scale(dfc$tas)[,1]
dfc$somatic_z <- scale(rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE))[,1]

# Create interaction term
dfc$ias_x_iats <- dfc$ias_z * dfc$iats_z

# Fit regression models with interaction
mod_somatic <- lm(somatic_z ~ ias_z + iats_z + ias_x_iats, data = dfc)
mod_tas <- lm(tas_z ~ ias_z + iats_z + ias_x_iats, data = dfc)

# Extract coefficients
coef_som <- coef(mod_somatic)
coef_tas <- coef(mod_tas)

# Get p-values for interaction
p_som_int <- summary(mod_somatic)$coefficients["ias_x_iats", "Pr(>|t|)"]
p_tas_int <- summary(mod_tas)$coefficients["ias_x_iats", "Pr(>|t|)"]

cat("Somatic interaction: b =", round(coef_som["ias_x_iats"], 3), ", p =", round(p_som_int, 4), "\n")
cat("TAS interaction: b =", round(coef_tas["ias_x_iats"], 3), ", p =", round(p_tas_int, 4), "\n")

# =============================================================================
# 2. CREATE SIMPLE SLOPES DATA
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

# Calculate simple slopes and test significance for Somatic
slopes_somatic <- data.frame()
results_somatic <- data.frame()

for (i in 1:length(iats_levels)) {
  iats_val <- iats_levels[i]

  # Simple slope: b1 + b3*IATS
  slope <- coef_som["ias_z"] + coef_som["ias_x_iats"] * iats_val
  intercept <- coef_som["(Intercept)"] + coef_som["iats_z"] * iats_val

  # Test simple slope significance
  # SE of simple slope = sqrt(Var(b1) + IATS^2*Var(b3) + 2*IATS*Cov(b1,b3))
  vcov_mat <- vcov(mod_somatic)
  se_slope <- sqrt(vcov_mat["ias_z", "ias_z"] +
                   iats_val^2 * vcov_mat["ias_x_iats", "ias_x_iats"] +
                   2 * iats_val * vcov_mat["ias_z", "ias_x_iats"])
  t_val <- slope / se_slope
  p_val <- 2 * pt(-abs(t_val), df = nrow(dfc) - 4)

  pred <- intercept + slope * ias_range

  slopes_somatic <- rbind(slopes_somatic, data.frame(
    IAS = ias_range,
    Outcome = pred,
    IATS_level = iats_labels[i],
    stringsAsFactors = FALSE
  ))

  results_somatic <- rbind(results_somatic, data.frame(
    IATS_level = iats_labels[i],
    Slope = round(slope, 3),
    p = p_val,
    sig = format_sig(p_val),
    label = paste0("\u03B2 = ", round(slope, 3), format_sig(p_val)),
    stringsAsFactors = FALSE
  ))
}

slopes_somatic$IATS_level <- factor(slopes_somatic$IATS_level, levels = iats_labels)

# Calculate simple slopes for TAS
slopes_tas <- data.frame()
results_tas <- data.frame()

for (i in 1:length(iats_levels)) {
  iats_val <- iats_levels[i]

  slope <- coef_tas["ias_z"] + coef_tas["ias_x_iats"] * iats_val
  intercept <- coef_tas["(Intercept)"] + coef_tas["iats_z"] * iats_val

  # Test simple slope significance
  vcov_mat <- vcov(mod_tas)
  se_slope <- sqrt(vcov_mat["ias_z", "ias_z"] +
                   iats_val^2 * vcov_mat["ias_x_iats", "ias_x_iats"] +
                   2 * iats_val * vcov_mat["ias_z", "ias_x_iats"])
  t_val <- slope / se_slope
  p_val <- 2 * pt(-abs(t_val), df = nrow(dfc) - 4)

  pred <- intercept + slope * ias_range

  slopes_tas <- rbind(slopes_tas, data.frame(
    IAS = ias_range,
    Outcome = pred,
    IATS_level = iats_labels[i],
    stringsAsFactors = FALSE
  ))

  results_tas <- rbind(results_tas, data.frame(
    IATS_level = iats_labels[i],
    Slope = round(slope, 3),
    p = p_val,
    sig = format_sig(p_val),
    label = paste0("\u03B2 = ", round(slope, 3), format_sig(p_val)),
    stringsAsFactors = FALSE
  ))
}

slopes_tas$IATS_level <- factor(slopes_tas$IATS_level, levels = iats_labels)

cat("\n=== SOMATIC: IAS -> Somatic at IATS levels ===\n")
print(results_somatic)

cat("\n=== TAS: IAS -> TAS at IATS levels ===\n")
print(results_tas)

# =============================================================================
# 3. CREATE PLOTS
# =============================================================================

# Label positions
label_y_positions <- c(1.05, 0.90, 0.75)

# Panel A: Somatic (MODERATED)
p_somatic <- ggplot(slopes_somatic, aes(x = IAS, y = Outcome, color = IATS_level)) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  scale_color_manual(values = iats_colors, name = "Interoceptive\nAttention (IATS)") +
  scale_x_continuous(breaks = -2:2, limits = c(-2.2, 2.2)) +
  scale_y_continuous(limits = c(-1.0, 1.2)) +
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
    title = "A. Somatic Pathway: Moderation by Interoceptive Attention",
    subtitle = paste0("IAS \u00d7 IATS interaction: \u03B2 = ", round(coef_som["ias_x_iats"], 3),
                     ", p = ", ifelse(p_som_int < 0.001, "< .001", round(p_som_int, 3))),
    x = "Interoceptive Accuracy (IAS)",
    y = "Somatic Symptoms\n(SSS-8 + PCS)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = c(0.18, 0.20),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# Panel B: TAS (UNMODERATED)
p_tas <- ggplot(slopes_tas, aes(x = IAS, y = Outcome, color = IATS_level)) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  scale_color_manual(values = iats_colors, name = "Interoceptive\nAttention (IATS)") +
  scale_x_continuous(breaks = -2:2, limits = c(-2.2, 2.2)) +
  scale_y_continuous(limits = c(-1.2, 1.2)) +
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
    title = "B. TAS Pathway: No Moderation by Interoceptive Attention",
    subtitle = paste0("IAS \u00d7 IATS interaction: \u03B2 = ", round(coef_tas["ias_x_iats"], 3),
                     ", p = ", round(p_tas_int, 2), " (n.s.)"),
    x = "Interoceptive Accuracy (IAS)",
    y = "Alexithymia\n(TAS-20)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = c(0.18, 0.20),
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
  top = textGrob("Direct Effects of IAS on Mediators: IATS Moderates Only the Somatic Pathway",
                 gp = gpar(fontsize = 14, fontface = "bold")),
  bottom = textGrob(sprintf("Note: Simple slopes showing IAS effect on each mediator at different IATS levels. ***p < .001, **p < .01, *p < .05. N = %d.", nrow(dfc)),
                    gp = gpar(fontsize = 9, col = "gray40"))
)

ggsave("plots/supplementary_plots/3_s_3_simple_slopes_continuous.png",
       combined, width = 13, height = 6, dpi = 300, bg = "white")

cat("\nFigure saved: figure_simple_slopes_continuous.png\n")

# =============================================================================
# 5. PRINT SUMMARY
# =============================================================================

cat("\n=== SUMMARY ===\n")
cat("\nSOMATIC (IAS -> Somatic):\n")
for (i in 1:nrow(results_somatic)) {
  cat(sprintf("  %s: %s (p = %.4f)\n", results_somatic$IATS_level[i],
              results_somatic$label[i], results_somatic$p[i]))
}

cat("\nTAS (IAS -> TAS):\n")
for (i in 1:nrow(results_tas)) {
  cat(sprintf("  %s: %s (p = %.4f)\n", results_tas$IATS_level[i],
              results_tas$label[i], results_tas$p[i]))
}

som_range <- max(results_somatic$Slope) - min(results_somatic$Slope)
tas_range <- max(results_tas$Slope) - min(results_tas$Slope)
cat(sprintf("\nSlope range - Somatic: %.3f | TAS: %.3f\n", som_range, tas_range))
