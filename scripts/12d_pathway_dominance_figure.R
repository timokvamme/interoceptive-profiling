# =============================================================================
# Panel C: Pathway Dominance with Error Bars
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# Data from bootstrap analysis (12b_pathway_dominance_test.R)
# These are the conditional indirect effect percentages with 95% CIs

pathway_data <- data.frame(
  Profile = factor(c("Hypervigilant", "Uncertain", "Efficient"),
                   levels = c("Hypervigilant", "Uncertain", "Efficient")),
  TAS_pct = c(57.2, 68.7, 78.6),
  TAS_ci_low = c(39.4, 53.1, 59.7),
  TAS_ci_high = c(84.6, 89.2, 98.2),
  Som_pct = c(42.8, 31.3, 21.4),
  Som_ci_low = c(15.4, 10.8, 1.8),
  Som_ci_high = c(60.6, 46.9, 40.3)
)

# Reshape for plotting
plot_data <- pathway_data %>%
  pivot_longer(
    cols = c(TAS_pct, Som_pct),
    names_to = "Pathway",
    values_to = "Percentage"
  ) %>%
  mutate(
    CI_low = ifelse(Pathway == "TAS_pct", TAS_ci_low, Som_ci_low),
    CI_high = ifelse(Pathway == "TAS_pct", TAS_ci_high, Som_ci_high),
    Pathway = factor(ifelse(Pathway == "TAS_pct", "TAS (Alexithymia)", "Somatic"),
                     levels = c("TAS (Alexithymia)", "Somatic"))
  ) %>%
  select(Profile, Pathway, Percentage, CI_low, CI_high)

# Colors matching the SEM diagram
colors <- c("TAS (Alexithymia)" = "#E57373", "Somatic" = "#F59E0B")

# =============================================================================
# VERSION 1: Grouped Bar Chart with Error Bars (Side by Side)
# =============================================================================

p1 <- ggplot(plot_data, aes(x = Profile, y = Percentage, fill = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(width = 0.8),
                width = 0.25, linewidth = 0.6) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title = "C. Pathway Dominance by Profile",
    subtitle = "IAS → MH indirect effect: TAS vs Somatic route (with 95% CIs)",
    x = NULL,
    y = "Percentage of Indirect Effect",
    fill = "Pathway"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  # Add percentage labels
  geom_text(aes(label = sprintf("%.0f%%", Percentage)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.5, fontface = "bold")

# Output disabled to avoid creating non-approved plots.

# =============================================================================
# VERSION 2: Stacked Bar with CI annotations (like original but with uncertainty)
# =============================================================================

# For stacked bars, we can add CI as text annotations
p2 <- ggplot(pathway_data, aes(x = Profile)) +
  # Somatic on bottom
  geom_bar(aes(y = Som_pct, fill = "Somatic"), stat = "identity",
           width = 0.7, color = "black", linewidth = 0.3) +
  # TAS on top
  geom_bar(aes(y = TAS_pct, fill = "TAS (Alexithymia)"), stat = "identity",
           width = 0.7, color = "black", linewidth = 0.3,
           position = position_stack()) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("Somatic" = "#F59E0B", "TAS (Alexithymia)" = "#E57373"),
                    name = "Pathway") +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title = "C. Pathway Dominance by Profile",
    subtitle = "IAS → MH indirect effect (bootstrap 95% CIs shown)",
    x = NULL,
    y = "Percentage of Indirect Effect"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    legend.position = "top"
  ) +
  # Add percentage labels with CIs
  geom_text(aes(y = Som_pct/2,
                label = sprintf("%.0f%%\n[%.0f-%.0f]", Som_pct, Som_ci_low, Som_ci_high)),
            color = "white", size = 3.2, fontface = "bold", lineheight = 0.9) +
  geom_text(aes(y = Som_pct + TAS_pct/2,
                label = sprintf("%.0f%%\n[%.0f-%.0f]", TAS_pct, TAS_ci_low, TAS_ci_high)),
            color = "white", size = 3.2, fontface = "bold", lineheight = 0.9)

# Output disabled to avoid creating non-approved plots.

# =============================================================================
# VERSION 3: Point + Range plot (Cleveland dot plot style)
# =============================================================================

p3 <- ggplot(plot_data, aes(x = Percentage, y = Profile, color = Pathway)) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "gray50") +
  geom_linerange(aes(xmin = CI_low, xmax = CI_high),
                 position = position_dodge(width = 0.5),
                 linewidth = 2, alpha = 0.4) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_text(aes(label = sprintf("%.0f%%", Percentage)),
            position = position_dodge(width = 0.5),
            vjust = -1.2, size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = colors) +
  scale_x_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title = "C. Pathway Dominance by Profile",
    subtitle = "IAS → MH indirect effect: point estimates with 95% bootstrap CIs",
    x = "Percentage of Indirect Effect",
    y = NULL,
    color = "Pathway"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# Output disabled to avoid creating non-approved plots.
