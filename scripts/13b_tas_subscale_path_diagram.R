# =============================================================================
# TAS SUBSCALE PATH DIAGRAMS: Three-Panel SEM Diagram (DIF | DDF | EOT)
# =============================================================================
#
# Same layout as the dual-pathway SEM diagram (Figure 3B / Supp Fig S1)
# but with each TAS subscale replacing TAS total. The visual contrast in
# arrow thickness immediately shows DIF is the dominant mediator.
#
# Output: plots/supplementary_plots/3_s_9_tas_subscale_path_diagrams.png
#
# =============================================================================

setwd("C:/code/projects/intero_mod")

library(dplyr)
library(lavaan)
library(ggplot2)
library(ggforce)
library(gridExtra)
library(grid)

# =============================================================================
# 1. LOAD DATA AND COMPUTE SUBSCALES
# =============================================================================

dfc <- read.csv("C:/code/projects/mi/analyses/soma/results/dfc_interoception_profiling.csv")

# TAS-20 subscales (items already reverse-scored)
dfc$dif <- rowSums(dfc[, paste0("tas_", c(1, 3, 6, 7, 9, 13, 14))], na.rm = TRUE)
dfc$ddf <- rowSums(dfc[, paste0("tas_", c(2, 4, 11, 12, 17))], na.rm = TRUE)
dfc$eot <- rowSums(dfc[, paste0("tas_", c(5, 8, 10, 15, 16, 18, 19, 20))], na.rm = TRUE)

# Full-scale IAS/IATS
dfc$ias_full <- rowMeans(dfc[, paste0("ias_", 1:21)], na.rm = TRUE)
dfc$iats_full <- rowMeans(dfc[, paste0("iats_", 1:21)], na.rm = TRUE)

# Composites
dfc$mh_composite <- rowMeans(cbind(scale(dfc$phq9), scale(dfc$gad7),
                                    scale(dfc$stai)), na.rm = TRUE)
dfc$somatic_combined <- rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE)

cat("Data loaded. N =", nrow(dfc), "\n")

# =============================================================================
# 2. FIT SEPARATE DUAL-PATHWAY SEMS
# =============================================================================

fit_subscale_sem <- function(data, subscale_var) {
  data$ias_z <- scale(data$ias_full)[, 1]
  data$iats_z <- scale(data$iats_full)[, 1]
  data$sub_z <- scale(data[[subscale_var]])[, 1]
  data$somatic_z <- scale(data$somatic_combined)[, 1]
  data$mh_z <- scale(data$mh_composite)[, 1]
  data$ias_x_iats <- data$ias_z * data$iats_z

  model <- '
    sub_z ~ a1*ias_z + a2*iats_z + a3*ias_x_iats
    somatic_z ~ b1*iats_z + b2*ias_z + b3*ias_x_iats
    mh_z ~ c1*sub_z + c2*somatic_z + d1*ias_z + d2*iats_z + d3*ias_x_iats
    ind_ias_sub := a1 * c1
    ind_ias_som := b2 * c2
  '

  fit <- sem(model, data = data)
  params <- parameterEstimates(fit, standardized = TRUE)
  r2 <- inspect(fit, "rsquare")

  get_val <- function(lbl) {
    row <- params[params$label == lbl, ]
    if (nrow(row) > 0) return(list(est = round(row$std.all[1], 2),
                                    p = row$pvalue[1]))
    return(list(est = NA, p = NA))
  }

  list(
    ias_sub = get_val("a1"), iats_sub = get_val("a2"), ixi_sub = get_val("a3"),
    ias_som = get_val("b2"), iats_som = get_val("b1"), ixi_som = get_val("b3"),
    sub_mh = get_val("c1"), som_mh = get_val("c2"),
    ias_mh = get_val("d1"), iats_mh = get_val("d2"), ixi_mh = get_val("d3"),
    ind_sub = get_val("ind_ias_sub"), ind_som = get_val("ind_ias_som"),
    r2_sub = round(r2["sub_z"] * 100), r2_som = round(r2["somatic_z"] * 100),
    r2_mh = round(r2["mh_z"] * 100)
  )
}

cat("Fitting SEMs...\n")
res_dif <- fit_subscale_sem(dfc, "dif")
res_ddf <- fit_subscale_sem(dfc, "ddf")
res_eot <- fit_subscale_sem(dfc, "eot")
cat("Done.\n")

# =============================================================================
# 3. DRAW PATH DIAGRAM FUNCTION
# =============================================================================

format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("< .001")
  if (p < 0.01) return("< .01")
  if (p < 0.05) return("< .05")
  return(sprintf("%.3f", p))
}

draw_sem_panel <- function(res, sub_name, sub_color, pct_alex, panel_label) {

  # Arrow color/thickness based on coefficient
  arrow_style <- function(beta, p, is_direct = FALSE) {
    if (is_direct) return(list(color = "gray50", lw = 0.6, lt = "dotted"))
    if (is.na(p) || p >= 0.05) return(list(color = "gray60", lw = 0.5, lt = "dashed"))
    clr <- if (beta > 0) "#DC2626" else "#10B981"
    lw <- max(0.6, min(3.0, abs(beta) * 5))
    list(color = clr, lw = lw, lt = "solid")
  }

  # Node coordinates
  ias_x <- 1.5; ias_y <- 5.8
  iats_x <- 1.5; iats_y <- 1.8
  sub_x <- 4.8; sub_y <- 7.0
  som_x <- 4.8; som_y <- 1.0
  mh_x <- 8.0; mh_y <- 4.0
  ind_x <- 10.7; ind_w <- 0.6

  p <- ggplot() + theme_void() +
    xlim(0, 12) + ylim(-0.5, 9.5) +
    theme(plot.margin = margin(2, 2, 2, 2, "pt"))

  # --- Title ---
  p <- p +
    annotate("text", x = 6, y = 9.3,
             label = sprintf("%s. %s (%s of alexithymia pathway)",
                             panel_label, sub_name, pct_alex),
             size = 4.2, fontface = "bold") +
    annotate("text", x = 6, y = 9.0,
             label = "Red = worsens MH | Green = improves MH | N = 833",
             size = 2.8, color = "gray40")

  # --- Nodes ---
  # IAS
  p <- p +
    geom_ellipse(aes(x0 = ias_x, y0 = ias_y, a = 0.9, b = 0.9, angle = 0),
                 fill = "#0891B2", color = "black", linewidth = 0.8) +
    annotate("text", x = ias_x, y = ias_y + 0.2, label = "IAS",
             color = "white", size = 5.5, fontface = "bold") +
    annotate("text", x = ias_x, y = ias_y - 0.12, label = "Interoceptive",
             color = "white", size = 3.8) +
    annotate("text", x = ias_x, y = ias_y - 0.38, label = "Accuracy",
             color = "white", size = 3.8)

  # IATS
  p <- p +
    geom_ellipse(aes(x0 = iats_x, y0 = iats_y, a = 0.9, b = 0.9, angle = 0),
                 fill = "#7C3AED", color = "black", linewidth = 0.8) +
    annotate("text", x = iats_x, y = iats_y + 0.2, label = "IATS",
             color = "white", size = 5.5, fontface = "bold") +
    annotate("text", x = iats_x, y = iats_y - 0.12, label = "Interoceptive",
             color = "white", size = 3.8) +
    annotate("text", x = iats_x, y = iats_y - 0.38, label = "Attention",
             color = "white", size = 3.8)

  # Subscale node
  p <- p +
    geom_ellipse(aes(x0 = sub_x, y0 = sub_y, a = 0.85, b = 0.85, angle = 0),
                 fill = sub_color, color = "black", linewidth = 0.8) +
    annotate("text", x = sub_x, y = sub_y + 0.15, label = sub_name,
             size = 5.5, fontface = "bold") +
    annotate("text", x = sub_x, y = sub_y - 0.2, label = "(Alexithymia)",
             size = 3.5) +
    annotate("text", x = sub_x, y = 8.1,
             label = sprintf("R\u00b2 = %d%%", res$r2_sub),
             size = 3.8, color = "gray30")

  # Somatic node
  p <- p +
    geom_ellipse(aes(x0 = som_x, y0 = som_y, a = 0.85, b = 0.85, angle = 0),
                 fill = "#F59E0B", color = "black", linewidth = 0.8) +
    annotate("text", x = som_x, y = som_y + 0.15, label = "Somatic",
             size = 5.5, fontface = "bold") +
    annotate("text", x = som_x, y = som_y - 0.25, label = "(SSS + PCS)",
             size = 3.5) +
    annotate("text", x = som_x, y = -0.1,
             label = sprintf("R\u00b2 = %d%%", res$r2_som),
             size = 3.8, color = "gray30")

  # Mental Health node
  p <- p +
    geom_ellipse(aes(x0 = mh_x, y0 = mh_y, a = 1.1, b = 1.1, angle = 0),
                 fill = "#22C55E", color = "black", linewidth = 0.8) +
    annotate("text", x = mh_x, y = mh_y + 0.15, label = "Mental",
             size = 5.5, fontface = "bold", color = "white") +
    annotate("text", x = mh_x, y = mh_y - 0.2, label = "Health",
             size = 5.5, fontface = "bold", color = "white") +
    annotate("text", x = mh_x, y = 2.5,
             label = sprintf("R\u00b2 = %d%%", res$r2_mh),
             size = 3.8, color = "gray30")

  # MH indicators
  p <- p +
    annotate("rect", xmin = 10.3, xmax = 11.4, ymin = 5.8, ymax = 6.3,
             fill = "#22C55E", color = "black") +
    annotate("text", x = 10.85, y = 6.05, label = "PHQ-9",
             size = 4, color = "white", fontface = "bold") +
    annotate("rect", xmin = 10.3, xmax = 11.4, ymin = 3.7, ymax = 4.2,
             fill = "#22C55E", color = "black") +
    annotate("text", x = 10.85, y = 3.95, label = "GAD-7",
             size = 4, color = "white", fontface = "bold") +
    annotate("rect", xmin = 10.3, xmax = 11.4, ymin = 1.6, ymax = 2.1,
             fill = "#22C55E", color = "black") +
    annotate("text", x = 10.85, y = 1.85, label = "STAI",
             size = 4, color = "white", fontface = "bold") +
    # MH → indicator arrows
    annotate("segment", x = 8.9, y = 4.8, xend = 10.3, yend = 5.9,
             arrow = arrow(length = unit(0.12, "cm"), ends = "both"), linewidth = 0.7) +
    annotate("segment", x = 9.1, y = 4.0, xend = 10.3, yend = 3.95,
             arrow = arrow(length = unit(0.12, "cm"), ends = "both"), linewidth = 0.7) +
    annotate("segment", x = 8.9, y = 3.2, xend = 10.3, yend = 2.0,
             arrow = arrow(length = unit(0.12, "cm"), ends = "both"), linewidth = 0.7)

  # --- Helper: draw labeled arrow ---
  add_arrow <- function(p, x0, y0, x1, y1, beta, pval, lbl_x, lbl_y,
                         is_direct = FALSE, show_ns = FALSE) {
    s <- arrow_style(beta, pval, is_direct)
    # Don't show n.s. label prefix for direct effects
    lbl <- if (!is.na(pval) && pval >= 0.05 && !is_direct) {
      sprintf("%.2f\n(n.s.)", beta)
    } else {
      sprintf("%.2f", beta)
    }

    p <- p +
      annotate("segment", x = x0, y = y0, xend = x1, yend = y1,
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               linewidth = s$lw, color = s$color, linetype = s$lt)

    # Coefficient label
    box_color <- if (is_direct) NA else s$color
    p <- p +
      annotate("rect",
               xmin = lbl_x - 0.38, xmax = lbl_x + 0.38,
               ymin = lbl_y - 0.2, ymax = lbl_y + 0.2,
               fill = "white", color = box_color, linewidth = 0.3) +
      annotate("text", x = lbl_x, y = lbl_y, label = lbl,
               size = 3.5, color = s$color, fontface = "bold")

    p
  }

  # --- Pathway arrows ---

  # IAS → Subscale
  p <- add_arrow(p, 2.4, 6.3, 4.0, 6.8,
                 res$ias_sub$est, res$ias_sub$p, 3.0, 6.8)

  # IATS → Subscale
  p <- add_arrow(p, 2.3, 2.5, 4.1, 6.2,
                 res$iats_sub$est, res$iats_sub$p, 3.5, 4.6)

  # IAS → Somatic
  p <- add_arrow(p, 2.3, 5.1, 4.1, 1.7,
                 res$ias_som$est, res$ias_som$p, 2.8, 3.2)

  # IATS → Somatic
  p <- add_arrow(p, 2.4, 1.3, 4.0, 1.1,
                 res$iats_som$est, res$iats_som$p, 3.1, 0.9)

  # Subscale → MH
  p <- add_arrow(p, 5.6, 6.4, 7.1, 4.9,
                 res$sub_mh$est, res$sub_mh$p, 6.5, 5.8)

  # Somatic → MH
  p <- add_arrow(p, 5.6, 1.6, 7.1, 3.2,
                 res$som_mh$est, res$som_mh$p, 6.5, 2.1)

  # Direct: IAS → MH (dotted)
  p <- add_arrow(p, 2.4, 5.6, 7.0, 4.5,
                 res$ias_mh$est, res$ias_mh$p, 4.8, 5.2, is_direct = TRUE)

  # Direct: IATS → MH (dotted)
  p <- add_arrow(p, 2.4, 2.1, 7.0, 3.5,
                 res$iats_mh$est, res$iats_mh$p, 4.8, 2.6, is_direct = TRUE)

  # --- Indirect effect annotation ---
  ind_sub_val <- res$ind_sub$est
  ind_som_val <- res$ind_som$est

  p <- p +
    annotate("text", x = 6, y = -0.45,
             label = sprintf("Indirect via %s: %.3f  |  Indirect via Somatic: %.3f",
                             sub_name, ind_sub_val, ind_som_val),
             size = 2.8, color = "gray30")

  p
}

# =============================================================================
# 4. CREATE THREE-PANEL FIGURE
# =============================================================================

cat("Drawing panels...\n")

p_dif <- draw_sem_panel(res_dif, "DIF", "#E64B35", "67.2%", "A")
p_ddf <- draw_sem_panel(res_ddf, "DDF", "#4DBBD5", "26.9%", "B")
p_eot <- draw_sem_panel(res_eot, "EOT", "#00A087", "5.9%", "C")

fig <- arrangeGrob(
  p_dif, p_ddf, p_eot,
  ncol = 3,
  top = textGrob(
    "TAS-20 Subscale Decomposition: Dual-Pathway SEM",
    gp = gpar(fontface = "bold", fontsize = 16)
  ),
  bottom = textGrob(
    "DIF = Difficulty Identifying Feelings  |  DDF = Difficulty Describing Feelings  |  EOT = Externally Oriented Thinking",
    gp = gpar(fontsize = 10, col = "gray30")
  )
)

fig_path <- file.path("plots/supplementary_plots",
                       "3_s_9_tas_subscale_path_diagrams.png")
ggsave(fig_path, fig, width = 24, height = 8, dpi = 300, bg = "white")

cat(sprintf("Figure saved: %s\n", fig_path))
cat("Done.\n")
