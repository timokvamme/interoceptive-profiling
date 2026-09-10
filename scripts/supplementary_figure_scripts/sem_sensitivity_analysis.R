# =============================================================================
# SEM Sensitivity Analysis: Full vs Minimal IAS/IATS Scales
# Compares SEM results using all 42 items vs top 21 minimal items
# =============================================================================

# Paths are relative to the repository root; run scripts from there.

library(tidyverse)
library(lavaan)
library(ggforce)
library(gridExtra)
library(grid)

# =============================================================================
# 1. DEFINE ITEM SETS
# =============================================================================

# Full scales (all 21 items each)
ias_full <- paste0("ias_", 1:21)
iats_full <- paste0("iats_", 1:21)

# Minimal scales (top 21 items from consensus ranking)
# 10 IAS + 11 IATS items
# Same minimal battery as 12_sem_pathway_analysis.R: the top 21 items of the
# consensus ranking written by 09_item_selection_analysis.R.
consensus_items <- read.csv("supplementary_data/item_consensus_ranking.csv")$item[1:21]
ias_minimal <- consensus_items[grepl("^ias_", consensus_items)]
iats_minimal <- consensus_items[grepl("^iats_", consensus_items)]

# =============================================================================
# 2. LOAD DATA AND COMPUTE SCALE SCORES
# =============================================================================

dfc <- read.csv("dfc_interoception_profiling.csv")

# Function to compute scale score from items
compute_scale <- function(data, items) {
  rowMeans(data[, items], na.rm = TRUE)
}

# Compute both versions
dfc$ias_full <- compute_scale(dfc, ias_full)
dfc$iats_full <- compute_scale(dfc, iats_full)
dfc$ias_minimal <- compute_scale(dfc, ias_minimal)
dfc$iats_minimal <- compute_scale(dfc, iats_minimal)

# Other variables (same for both)
dfc$tas <- compute_scale(dfc, paste0("tas_", 1:20))
dfc$somatic <- rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE)
dfc$psych_mh <- rowMeans(cbind(scale(dfc$phq9), scale(dfc$gad7), scale(dfc$stai)), na.rm = TRUE)

# =============================================================================
# 3. FIT SEM MODELS
# =============================================================================

fit_sem_model <- function(data, ias_var, iats_var, label_suffix = "") {
  # Standardize variables
  data$ias_z <- scale(data[[ias_var]])[,1]
  data$iats_z <- scale(data[[iats_var]])[,1]
  data$tas_z <- scale(data$tas)[,1]
  data$somatic_z <- scale(data$somatic)[,1]
  data$psych_mh_z <- scale(data$psych_mh)[,1]

  model <- '
    tas_z ~ a1*ias_z + a2*iats_z
    somatic_z ~ b1*iats_z + b2*ias_z
    psych_mh_z ~ c1*tas_z + c2*somatic_z + d1*ias_z + d2*iats_z

    ind_ias_tas := a1 * c1
    ind_iats_som := b1 * c2
    ind_ias_som := b2 * c2
    ind_iats_tas := a2 * c1
  '

  fit <- sem(model, data = data)
  params <- parameterEstimates(fit, standardized = TRUE)
  r2 <- inspect(fit, "rsquare")

  # Extract coefficients
  get_std <- function(label) round(params$std.all[params$label == label], 2)
  get_p <- function(label) params$pvalue[params$label == label]
  get_ind <- function(label) round(params$est[params$label == label], 3)
  get_ind_p <- function(label) params$pvalue[params$label == label]

  list(
    fit = fit,
    params = params,
    r2 = r2,
    coefficients = list(
      a1 = get_std("a1"), a2 = get_std("a2"),
      b1 = get_std("b1"), b2 = get_std("b2"),
      c1 = get_std("c1"), c2 = get_std("c2"),
      d1 = get_std("d1"), d2 = get_std("d2")
    ),
    pvalues = list(
      a1 = get_p("a1"), a2 = get_p("a2"),
      b1 = get_p("b1"), b2 = get_p("b2"),
      c1 = get_p("c1"), c2 = get_p("c2"),
      d1 = get_p("d1"), d2 = get_p("d2")
    ),
    indirect = list(
      ias_tas = get_ind("ind_ias_tas"),
      iats_som = get_ind("ind_iats_som"),
      ias_som = get_ind("ind_ias_som"),
      iats_tas = get_ind("ind_iats_tas")
    ),
    indirect_p = list(
      ias_tas = get_ind_p("ind_ias_tas"),
      iats_som = get_ind_p("ind_iats_som"),
      ias_som = get_ind_p("ind_ias_som"),
      iats_tas = get_ind_p("ind_iats_tas")
    ),
    data = data
  )
}

# Fit both models
cat("Fitting Full Scale Model (42 items)...\n")
fit_full <- fit_sem_model(dfc, "ias_full", "iats_full")

cat("Fitting Minimal Scale Model (21 items)...\n")
fit_minimal <- fit_sem_model(dfc, "ias_minimal", "iats_minimal")

# =============================================================================
# 4. COMPARISON TABLE
# =============================================================================

comparison <- data.frame(
  Path = c("IAS → TAS", "IATS → TAS", "IAS → Somatic", "IATS → Somatic",
           "TAS → MH", "Somatic → MH", "IAS → MH (direct)", "IATS → MH (direct)"),
  Full_beta = c(fit_full$coefficients$a1, fit_full$coefficients$a2,
                fit_full$coefficients$b2, fit_full$coefficients$b1,
                fit_full$coefficients$c1, fit_full$coefficients$c2,
                fit_full$coefficients$d1, fit_full$coefficients$d2),
  Minimal_beta = c(fit_minimal$coefficients$a1, fit_minimal$coefficients$a2,
                   fit_minimal$coefficients$b2, fit_minimal$coefficients$b1,
                   fit_minimal$coefficients$c1, fit_minimal$coefficients$c2,
                   fit_minimal$coefficients$d1, fit_minimal$coefficients$d2)
)
comparison$Difference <- comparison$Minimal_beta - comparison$Full_beta
comparison$Pct_Change <- round((comparison$Minimal_beta - comparison$Full_beta) / abs(comparison$Full_beta) * 100, 1)

cat("\n=============================================================================\n")
cat("PATH COEFFICIENT COMPARISON: Full (42 items) vs Minimal (21 items)\n")
cat("=============================================================================\n\n")
print(comparison, row.names = FALSE)

# R² comparison
cat("\n\nR² Comparison:\n")
cat(sprintf("  TAS:     Full = %.1f%%, Minimal = %.1f%%\n",
            fit_full$r2["tas_z"]*100, fit_minimal$r2["tas_z"]*100))
cat(sprintf("  Somatic: Full = %.1f%%, Minimal = %.1f%%\n",
            fit_full$r2["somatic_z"]*100, fit_minimal$r2["somatic_z"]*100))
cat(sprintf("  MH:      Full = %.1f%%, Minimal = %.1f%%\n",
            fit_full$r2["psych_mh_z"]*100, fit_minimal$r2["psych_mh_z"]*100))

# =============================================================================
# 5. DYNAMIC SEM DIAGRAM FUNCTION
# =============================================================================

create_sem_diagram <- function(sem_result, title_suffix = "", subtitle_extra = "") {

  coef <- sem_result$coefficients
  pval <- sem_result$pvalues
  r2 <- sem_result$r2
  data <- sem_result$data

  # Calculate MH indicator loadings
  mh_composite <- rowMeans(cbind(scale(data$phq9), scale(data$gad7), scale(data$stai)), na.rm = TRUE)
  loading_phq <- cor(data$phq9, mh_composite)
  loading_gad <- cor(data$gad7, mh_composite)
  loading_stai <- cor(data$stai, mh_composite)

  # Node definitions
  nodes <- data.frame(
    name = c("IAS", "IATS", "TAS", "Somatic", "MH", "PHQ9", "GAD7", "STAI"),
    label = c("IAS\nInteroceptive\nAccuracy", "IATS\nInteroceptive\nAttention", "TAS\n(Alexithymia)",
              "Somatic\n(SSS + PCS)", "Mental\nHealth", "PHQ-9", "GAD-7", "STAI"),
    x = c(1, 1, 3, 3, 5, 6.5, 6.5, 6.5),
    y = c(4, 1.5, 4.5, 1, 2.75, 4.2, 2.75, 1.3),
    color = c("#00BCD4", "#9b59b6", "#FF9800", "#e67e22", "#27ae60", "#2ecc71", "#2ecc71", "#2ecc71"),
    radius = c(0.55, 0.55, 0.55, 0.55, 0.65, 0.4, 0.4, 0.4),
    type = c("predictor", "predictor", "mediator", "mediator", "latent", "indicator", "indicator", "indicator"),
    stringsAsFactors = FALSE
  )

  # Path definitions
  paths <- data.frame(
    from = c("IAS", "IATS", "IAS", "IATS", "TAS", "Somatic", "IAS", "IATS"),
    to = c("TAS", "Somatic", "Somatic", "TAS", "MH", "MH", "MH", "MH"),
    label = c("a1", "b1", "b2", "a2", "c1", "c2", "d1", "d2"),
    type = c("mediated", "mediated", "mediated", "mediated", "mediated", "mediated", "direct", "direct"),
    stringsAsFactors = FALSE
  )

  paths$beta <- c(coef$a1, coef$b1, coef$b2, coef$a2, coef$c1, coef$c2, coef$d1, coef$d2)
  paths$p <- c(pval$a1, pval$b1, pval$b2, pval$a2, pval$c1, pval$c2, pval$d1, pval$d2)
  paths$sig <- paths$p < 0.05

  # Path coordinate calculator
  calc_path_coords <- function(from_name, to_name, nodes, curve_offset = 0, is_direct = FALSE) {
    from <- nodes[nodes$name == from_name, ]
    to <- nodes[nodes$name == to_name, ]

    dx <- to$x - from$x
    dy <- to$y - from$y
    dist <- sqrt(dx^2 + dy^2)

    from_pad <- from$radius + 0.12
    to_pad <- if (is_direct) to$radius + 0.25 else to$radius + 0.18

    x1 <- from$x + (dx / dist) * from_pad
    y1 <- from$y + (dy / dist) * from_pad
    x2 <- to$x - (dx / dist) * to_pad
    y2 <- to$y - (dy / dist) * to_pad

    return(list(x1 = x1, y1 = y1, x2 = x2, y2 = y2,
                midx = (x1+x2)/2 + curve_offset, midy = (y1+y2)/2))
  }

  # Build path data
  path_data <- data.frame()
  for (i in 1:nrow(paths)) {
    offset <- 0
    is_direct <- paths$type[i] == "direct"
    if (is_direct) {
      if (paths$from[i] == "IAS") offset <- 0.15
      if (paths$from[i] == "IATS") offset <- -0.15
    }

    coords <- calc_path_coords(paths$from[i], paths$to[i], nodes, offset, is_direct)
    path_data <- rbind(path_data, data.frame(
      x = coords$x1, y = coords$y1,
      xend = coords$x2, yend = coords$y2,
      midx = coords$midx, midy = coords$midy,
      beta = paths$beta[i],
      p = paths$p[i],
      sig = paths$sig[i],
      type = paths$type[i],
      from = paths$from[i],
      to = paths$to[i]
    ))
  }

  # Label adjustments
  path_data$label_x <- path_data$midx
  path_data$label_y <- path_data$midy

  path_data$label_y[path_data$from == "IAS" & path_data$to == "Somatic"] <-
    path_data$label_y[path_data$from == "IAS" & path_data$to == "Somatic"] - 0.25
  path_data$label_y[path_data$from == "IATS" & path_data$to == "TAS"] <-
    path_data$label_y[path_data$from == "IATS" & path_data$to == "TAS"] + 0.25
  path_data$label_y[path_data$from == "IAS" & path_data$to == "MH"] <-
    path_data$label_y[path_data$from == "IAS" & path_data$to == "MH"] + 0.20
  path_data$label_y[path_data$from == "IATS" & path_data$to == "MH"] <-
    path_data$label_y[path_data$from == "IATS" & path_data$to == "MH"] - 0.20

  # Create plot
  p_sem <- ggplot() +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40")
    ) +
    coord_fixed(ratio = 1, xlim = c(-0.2, 7.5), ylim = c(0, 5.5))

  # Draw mediated paths first
  for (i in 1:nrow(path_data)) {
    row <- path_data[i, ]
    if (row$type == "direct") next

    arrow_color <- ifelse(row$beta > 0, "#e74c3c", "#27ae60")
    line_width <- abs(row$beta) * 4 + 1

    p_sem <- p_sem +
      geom_segment(
        data = row,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = arrow_color, linewidth = line_width, linetype = "solid",
        arrow = arrow(length = unit(0.15, "inches"), type = "closed")
      )

    label_text <- sprintf("%.2f", row$beta)
    if (!row$sig) label_text <- paste0(label_text, " (n.s.)")

    p_sem <- p_sem +
      geom_label(
        data = row, aes(x = label_x, y = label_y),
        label = label_text,
        fill = "white", color = arrow_color, size = 4, fontface = "bold",
        label.padding = unit(0.25, "lines")
      )
  }

  # Draw direct paths on top
  for (i in 1:nrow(path_data)) {
    row <- path_data[i, ]
    if (row$type != "direct") next

    arrow_color <- ifelse(row$beta > 0, "#e74c3c", "#27ae60")
    if (row$p > 0.1) arrow_color <- "gray60"

    p_sem <- p_sem +
      geom_segment(
        data = row,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = arrow_color, linewidth = 1.5, linetype = "dotted",
        arrow = arrow(length = unit(0.15, "inches"), type = "closed")
      )

    label_text <- sprintf("%.2f", row$beta)
    if (!row$sig) label_text <- paste0(label_text, " (n.s.)")

    p_sem <- p_sem +
      geom_label(
        data = row, aes(x = label_x, y = label_y),
        label = label_text,
        fill = "white", color = arrow_color, size = 4, fontface = "bold",
        label.padding = unit(0.25, "lines")
      )
  }

  # Draw main nodes
  main_nodes <- nodes[nodes$type != "indicator", ]
  for (i in 1:nrow(main_nodes)) {
    node <- main_nodes[i, ]

    p_sem <- p_sem +
      geom_circle(
        data = data.frame(x0 = node$x, y0 = node$y, r = node$radius),
        aes(x0 = x0, y0 = y0, r = r),
        fill = node$color, color = "gray30", linewidth = 1.5
      ) +
      geom_text(
        data = data.frame(x = node$x, y = node$y, label = node$label),
        aes(x = x, y = y, label = label),
        color = "white", size = 4, fontface = "bold", lineheight = 0.85
      )
  }

  # Draw indicator nodes
  indicator_nodes <- nodes[nodes$type == "indicator", ]
  for (i in 1:nrow(indicator_nodes)) {
    node <- indicator_nodes[i, ]

    p_sem <- p_sem +
      geom_rect(
        data = data.frame(xmin = node$x - 0.35, xmax = node$x + 0.35,
                          ymin = node$y - 0.3, ymax = node$y + 0.3),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = node$color, color = "gray30", linewidth = 1
      ) +
      geom_text(
        data = data.frame(x = node$x, y = node$y, label = node$label),
        aes(x = x, y = y, label = label),
        color = "white", size = 3.5, fontface = "bold"
      )
  }

  # Draw bidirectional arrows to indicators
  mh_node <- nodes[nodes$name == "MH", ]
  loadings <- c(PHQ9 = loading_phq, GAD7 = loading_gad, STAI = loading_stai)

  for (ind_name in c("PHQ9", "GAD7", "STAI")) {
    ind_node <- nodes[nodes$name == ind_name, ]
    loading_val <- loadings[ind_name]

    dx <- ind_node$x - mh_node$x
    dy <- ind_node$y - mh_node$y
    dist <- sqrt(dx^2 + dy^2)

    mh_edge_x <- mh_node$x + (dx / dist) * (mh_node$radius + 0.1)
    mh_edge_y <- mh_node$y + (dy / dist) * (mh_node$radius + 0.1)
    ind_edge_x <- ind_node$x - 0.4
    ind_edge_y <- ind_node$y

    p_sem <- p_sem +
      geom_segment(
        data = data.frame(x = mh_edge_x, y = mh_edge_y, xend = ind_edge_x, yend = ind_edge_y),
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "black", linewidth = 1.5,
        arrow = arrow(length = unit(0.10, "inches"), type = "closed", ends = "both")
      ) +
      geom_label(
        data = data.frame(x = (mh_edge_x + ind_edge_x) / 2, y = (mh_edge_y + ind_edge_y) / 2),
        aes(x = x, y = y),
        label = sprintf("%.2f", loading_val),
        fill = "white", color = "black", size = 3.5, fontface = "bold",
        label.padding = unit(0.15, "lines")
      )
  }

  # Add R² labels
  p_sem <- p_sem +
    annotate("text", x = 3, y = 5.3,
             label = sprintf("R² = %.0f%%", r2["tas_z"] * 100),
             size = 3.5, color = "gray30") +
    annotate("text", x = 3, y = 0.2,
             label = sprintf("R² = %.0f%%", r2["somatic_z"] * 100),
             size = 3.5, color = "gray30") +
    annotate("text", x = 5, y = 1.9,
             label = sprintf("R² = %.0f%%", r2["psych_mh_z"] * 100),
             size = 3.5, color = "gray30")

  subtitle <- "Red = worsens MH | Green = improves MH | Dotted = direct"
  if (subtitle_extra != "") subtitle <- paste0(subtitle, "\n", subtitle_extra)

  p_sem <- p_sem +
    labs(
      title = paste0("SEM: Dual-Pathway Model", title_suffix),
      subtitle = subtitle
    )

  return(p_sem)
}

# =============================================================================
# 6. CREATE COMPARISON FIGURE
# =============================================================================

cat("\nCreating comparison figure...\n")

p_full <- create_sem_diagram(fit_full, " (Full: 42 items)",
                              sprintf("IAS: 21 items | IATS: 21 items | N = %d", nrow(dfc)))

p_minimal <- create_sem_diagram(fit_minimal, " (Minimal: 21 items)",
                                 sprintf("IAS: 10 items | IATS: 11 items | N = %d", nrow(dfc)))

# Combine side by side
combined <- grid.arrange(p_full, p_minimal, ncol = 2)

ggsave("plots/supplementary_plots/3_s_1_sem_sensitivity_comparison.png",
       combined, width = 22, height = 10, dpi = 200, bg = "white")

cat("Figure saved: figure_sem_sensitivity_comparison.png\n")

# =============================================================================
# 7. SAVE COMPARISON RESULTS
# =============================================================================

output_file <- "analysis_output/sem_sensitivity_sem_sensitivity_results.txt"
sink(output_file)

cat("=============================================================================\n")
cat("SEM SENSITIVITY ANALYSIS: Full vs Minimal Interoceptive Scales\n")
cat("=============================================================================\n")
cat(sprintf("N = %d\n", nrow(dfc)))
cat(sprintf("Date: %s\n\n", Sys.time()))

cat("SCALE DEFINITIONS:\n")
cat("  Full:    IAS (21 items) + IATS (21 items) = 42 items total\n")
cat("  Minimal: IAS (10 items) + IATS (11 items) = 21 items total\n")
cat("           Selected via consensus ranking (LASSO + RF + Network centrality)\n\n")

cat("=============================================================================\n")
cat("PATH COEFFICIENT COMPARISON\n")
cat("=============================================================================\n\n")

cat(sprintf("%-20s %10s %10s %10s %10s\n", "Path", "Full", "Minimal", "Diff", "% Change"))
cat(sprintf("%s\n", paste(rep("-", 65), collapse = "")))
for (i in 1:nrow(comparison)) {
  cat(sprintf("%-20s %10.2f %10.2f %10.2f %9.1f%%\n",
              comparison$Path[i], comparison$Full_beta[i], comparison$Minimal_beta[i],
              comparison$Difference[i], comparison$Pct_Change[i]))
}

cat("\n=============================================================================\n")
cat("VARIANCE EXPLAINED (R²) COMPARISON\n")
cat("=============================================================================\n\n")

cat(sprintf("%-15s %10s %10s %10s\n", "Variable", "Full", "Minimal", "Difference"))
cat(sprintf("%s\n", paste(rep("-", 50), collapse = "")))
cat(sprintf("%-15s %9.1f%% %9.1f%% %9.1f%%\n", "TAS",
            fit_full$r2["tas_z"]*100, fit_minimal$r2["tas_z"]*100,
            (fit_minimal$r2["tas_z"] - fit_full$r2["tas_z"])*100))
cat(sprintf("%-15s %9.1f%% %9.1f%% %9.1f%%\n", "Somatic",
            fit_full$r2["somatic_z"]*100, fit_minimal$r2["somatic_z"]*100,
            (fit_minimal$r2["somatic_z"] - fit_full$r2["somatic_z"])*100))
cat(sprintf("%-15s %9.1f%% %9.1f%% %9.1f%%\n", "Mental Health",
            fit_full$r2["psych_mh_z"]*100, fit_minimal$r2["psych_mh_z"]*100,
            (fit_minimal$r2["psych_mh_z"] - fit_full$r2["psych_mh_z"])*100))

cat("\n=============================================================================\n")
cat("INDIRECT EFFECTS COMPARISON\n")
cat("=============================================================================\n\n")

indirect_comparison <- data.frame(
  Pathway = c("IAS → TAS → MH", "IATS → Somatic → MH", "IAS → Somatic → MH", "IATS → TAS → MH"),
  Full = c(fit_full$indirect$ias_tas, fit_full$indirect$iats_som,
           fit_full$indirect$ias_som, fit_full$indirect$iats_tas),
  Minimal = c(fit_minimal$indirect$ias_tas, fit_minimal$indirect$iats_som,
              fit_minimal$indirect$ias_som, fit_minimal$indirect$iats_tas)
)

cat(sprintf("%-25s %10s %10s\n", "Pathway", "Full", "Minimal"))
cat(sprintf("%s\n", paste(rep("-", 50), collapse = "")))
for (i in 1:nrow(indirect_comparison)) {
  cat(sprintf("%-25s %10.3f %10.3f\n",
              indirect_comparison$Pathway[i],
              indirect_comparison$Full[i],
              indirect_comparison$Minimal[i]))
}

cat("\n=============================================================================\n")
cat("INTERPRETATION\n")
cat("=============================================================================\n\n")

# Calculate correlation between full and minimal scales
cor_ias <- cor(dfc$ias_full, dfc$ias_minimal)
cor_iats <- cor(dfc$iats_full, dfc$iats_minimal)

cat(sprintf("Scale correlations (Full vs Minimal):\n"))
cat(sprintf("  IAS:  r = %.3f\n", cor_ias))
cat(sprintf("  IATS: r = %.3f\n\n", cor_iats))

avg_change <- mean(abs(comparison$Pct_Change))
cat(sprintf("Average absolute change in path coefficients: %.1f%%\n\n", avg_change))

if (avg_change < 10) {
  cat("CONCLUSION: The minimal 21-item scale produces highly similar results\n")
  cat("to the full 42-item scale. The SEM pathway structure is robust to\n")
  cat("item reduction, supporting the use of the minimal scale in future studies.\n")
} else if (avg_change < 25) {
  cat("CONCLUSION: The minimal scale produces moderately different results.\n")
  cat("Some pathways show meaningful changes. Consider which version is most\n")
  cat("appropriate for your research question.\n")
} else {
  cat("CONCLUSION: Substantial differences between full and minimal scales.\n")
  cat("The item selection may have altered the construct coverage.\n")
  cat("Careful consideration needed before using minimal scale.\n")
}

sink()

cat(sprintf("Results saved to: %s\n", output_file))

# =============================================================================
# 8. SAVE INDIVIDUAL MODEL RESULTS
# =============================================================================

# Function to write detailed results for a single model
write_model_results <- function(sem_result, output_path, model_name, n_ias, n_iats, ias_items, iats_items) {
  coef <- sem_result$coefficients
  pval <- sem_result$pvalues
  r2 <- sem_result$r2
  ind <- sem_result$indirect
  ind_p <- sem_result$indirect_p

  format_p <- function(p) {
    if (p < 0.001) return("< .001")
    if (p < 0.01) return("< .01")
    if (p < 0.05) return("< .05")
    return(sprintf("%.3f", p))
  }

  sink(output_path)

  cat("=============================================================================\n")
  cat(sprintf("DUAL-PATHWAY SEM: %s\n", model_name))
  cat("=============================================================================\n")
  cat(sprintf("N = %d\n", nrow(dfc)))
  cat(sprintf("Date: %s\n\n", Sys.time()))

  cat("SCALE COMPOSITION:\n")
  cat(sprintf("  IAS:  %d items\n", n_ias))
  cat(sprintf("  IATS: %d items\n", n_iats))
  cat(sprintf("  Total: %d items\n\n", n_ias + n_iats))

  cat("IAS Items:\n")
  cat(sprintf("  %s\n\n", paste(ias_items, collapse = ", ")))
  cat("IATS Items:\n")
  cat(sprintf("  %s\n\n", paste(iats_items, collapse = ", ")))

  cat("=============================================================================\n")
  cat("VARIANCE EXPLAINED (R²)\n")
  cat("=============================================================================\n")
  cat(sprintf("  TAS (Alexithymia): %.1f%%\n", r2["tas_z"] * 100))
  cat(sprintf("  Somatic Symptoms:  %.1f%%\n", r2["somatic_z"] * 100))
  cat(sprintf("  Mental Health:     %.1f%%\n\n", r2["psych_mh_z"] * 100))

  cat("=============================================================================\n")
  cat("DIRECT EFFECTS (Standardized β)\n")
  cat("=============================================================================\n\n")

  cat("TO MEDIATORS:\n")
  cat(sprintf("  IAS → TAS:        β = %.2f, %s\n", coef$a1, format_p(pval$a1)))
  cat(sprintf("  IATS → TAS:       β = %.2f, %s\n", coef$a2, format_p(pval$a2)))
  cat(sprintf("  IATS → Somatic:   β = %.2f, %s\n", coef$b1, format_p(pval$b1)))
  cat(sprintf("  IAS → Somatic:    β = %.2f, %s\n\n", coef$b2, format_p(pval$b2)))

  cat("TO MENTAL HEALTH:\n")
  cat(sprintf("  TAS → MH:          β = %.2f, %s\n", coef$c1, format_p(pval$c1)))
  cat(sprintf("  Somatic → MH:      β = %.2f, %s\n", coef$c2, format_p(pval$c2)))
  sig_d1 <- if (pval$d1 < 0.05) "" else " ** NOT SIGNIFICANT **"
  sig_d2 <- if (pval$d2 < 0.05) "" else " ** NOT SIGNIFICANT **"
  cat(sprintf("  IAS → MH (direct): β = %.2f, %s%s\n", coef$d1, format_p(pval$d1), sig_d1))
  cat(sprintf("  IATS → MH (direct): β = %.2f, %s%s\n\n", coef$d2, format_p(pval$d2), sig_d2))

  cat("=============================================================================\n")
  cat("INDIRECT EFFECTS (Mediation)\n")
  cat("=============================================================================\n\n")

  # Calculate total effect for percentages
  ind_effects <- c(abs(ind$ias_tas), abs(ind$iats_som), abs(ind$ias_som), abs(ind$iats_tas))
  direct_effects <- c(abs(coef$d1), abs(coef$d2))
  total_effect <- sum(ind_effects) + sum(direct_effects)

  cat(sprintf("%-25s %10s %8s %10s\n", "Pathway", "Effect", "% Total", "p-value"))
  cat(sprintf("%s\n", paste(rep("-", 55), collapse = "")))
  cat(sprintf("%-25s %10.3f %7.1f%% %10s\n", "IAS → TAS → MH",
              ind$ias_tas, abs(ind$ias_tas)/total_effect*100, format_p(ind_p$ias_tas)))
  cat(sprintf("%-25s %10.3f %7.1f%% %10s\n", "IATS → TAS → MH",
              ind$iats_tas, abs(ind$iats_tas)/total_effect*100, format_p(ind_p$iats_tas)))
  cat(sprintf("%-25s %10.3f %7.1f%% %10s\n", "IAS → Somatic → MH",
              ind$ias_som, abs(ind$ias_som)/total_effect*100, format_p(ind_p$ias_som)))
  cat(sprintf("%-25s %10.3f %7.1f%% %10s\n", "IATS → Somatic → MH",
              ind$iats_som, abs(ind$iats_som)/total_effect*100, format_p(ind_p$iats_som)))

  cat(sprintf("\n%s\n", paste(rep("-", 55), collapse = "")))
  cat(sprintf("Total indirect effect: %.1f%% of total effect on MH\n",
              sum(ind_effects) / total_effect * 100))
  cat(sprintf("Total direct effect:   %.1f%% of total effect on MH\n\n",
              sum(direct_effects) / total_effect * 100))

  cat("=============================================================================\n")
  cat("KEY FINDING: PATHWAY DOMINANCE\n")
  cat("=============================================================================\n\n")

  # IAS pathway specificity
  ias_tas_pct <- abs(ind$ias_tas) / (abs(ind$ias_tas) + abs(ind$ias_som)) * 100
  iats_som_pct <- abs(ind$iats_som) / (abs(ind$iats_som) + abs(ind$iats_tas)) * 100

  cat(sprintf("IAS operates primarily through TAS (Alexithymia):\n"))
  cat(sprintf("  - IAS → TAS → MH:     %.1f%% of IAS's indirect effect\n", ias_tas_pct))
  cat(sprintf("  - IAS → Somatic → MH: %.1f%% of IAS's indirect effect\n\n", 100 - ias_tas_pct))

  cat(sprintf("IATS operates primarily through Somatic:\n"))
  cat(sprintf("  - IATS → Somatic → MH: %.1f%% of IATS's indirect effect\n", iats_som_pct))
  cat(sprintf("  - IATS → TAS → MH:     %.1f%% of IATS's indirect effect\n", 100 - iats_som_pct))

  sink()
}

# Write results for full model
write_model_results(
  fit_full,
  "analysis_output/sem_sensitivity_sem_results_full_42_items.txt",
  "Full Scale (42 Items)",
  n_ias = 21, n_iats = 21,
  ias_items = ias_full,
  iats_items = iats_full
)
cat("Full model results saved: sem_results_full_42_items.txt\n")

# Write results for minimal model
write_model_results(
  fit_minimal,
  "analysis_output/sem_sensitivity_sem_results_minimal_21_items.txt",
  "Minimal Scale (21 Items)",
  n_ias = 10, n_iats = 11,
  ias_items = ias_minimal,
  iats_items = iats_minimal
)
cat("Minimal model results saved: sem_results_minimal_21_items.txt\n")

# =============================================================================
# 9. SAVE COMPARISON DATA
# =============================================================================

write.csv(comparison,
          "analysis_output/sem_sensitivity_sem_sensitivity_coefficients.csv",
          row.names = FALSE)

cat("Coefficient comparison saved: sem_sensitivity_coefficients.csv\n")
