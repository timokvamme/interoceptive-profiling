# =============================================================================
# Publication Figure Generation
# =============================================================================
#
# CONSOLIDATED SCRIPT: Creates publication-quality combined figures from
# the item selection and network analysis results.
#
# READS FROM (run these scripts first):
#   - 1_item_selection_analysis.R outputs
#   - 2_network_analysis.R outputs
#
# OUTPUTS:
#   - figure_combined.png (main publication figure)
#   - fig2_A_item_importance.png
#   - fig2_B_saturation_curve.png
#   - figure_minimal_sets.png
#   - fig2_C_centrality_vs_prediction.png
#   - figure_network.png (internal network visualization)
#   - top_10_items_full_text.csv
#
# =============================================================================

setwd("C:/code/projects/interoceptive-profiling")

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(igraph)
library(ggrepel)

supplementary_data_dir <- "C:/code/projects/interoceptive-profiling/supplementary_data"
plots_dir <- "C:/code/projects/interoceptive-profiling/plots"
sub_plots_dir <- file.path(plots_dir, "sub_plots")
supplementary_plots_dir <- file.path(plots_dir, "supplementary_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sub_plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supplementary_plots_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD RESULTS FROM PREVIOUS SCRIPTS
# =============================================================================

consensus_df <- read.csv(file.path(supplementary_data_dir, "item_consensus_ranking.csv"))
saturation_results <- read.csv(file.path(supplementary_data_dir, "saturation_results.csv"))
config_results <- read.csv(file.path(supplementary_data_dir, "minimal_set_comparison.csv"))
centrality_df <- read.csv(file.path(supplementary_data_dir, "network_centrality.csv"))
dfc <- read.csv("C:/code/projects/interoceptive-profiling/dfc_interoception_profiling.csv")

# =============================================================================
# 2. FULL QUESTIONNAIRE TEXT
# =============================================================================

ias_full_text <- c(
  "ias_1" = "I can always accurately perceive when my heart is beating fast",
  "ias_2" = "I can always accurately perceive when I am hungry",
  "ias_3" = "I can always accurately perceive when I am breathing fast",
  "ias_4" = "I can always accurately perceive when I am thirsty",
  "ias_5" = "I can always accurately perceive when I need to urinate",
  "ias_6" = "I can always accurately perceive when I need to defecate",
  "ias_7" = "I can always accurately perceive when I encounter different tastes",
  "ias_8" = "I can always accurately perceive when I am going to vomit",
  "ias_9" = "I can always accurately perceive when I am going to sneeze",
  "ias_10" = "I can always accurately perceive when I am going to cough",
  "ias_11" = "I can always accurately perceive when I am hot/cold",
  "ias_12" = "I can always accurately perceive when I am sexually aroused",
  "ias_13" = "I can always accurately perceive when I am going to pass wind",
  "ias_14" = "I can always accurately perceive when I am going to burp",
  "ias_15" = "I can always accurately perceive when my muscles are tired/sore",
  "ias_16" = "I can always accurately perceive when I am going to get a bruise",
  "ias_17" = "I can always accurately perceive when I am in pain",
  "ias_18" = "I can always accurately perceive when my blood sugar is low",
  "ias_19" = "I can always accurately perceive when someone is touching me affectionately",
  "ias_20" = "I can always accurately perceive when something is going to be ticklish",
  "ias_21" = "I can always accurately perceive when something is going to be itchy"
)

iats_full_text <- c(
  "iats_1" = "Most of the time my attention is focused on whether my heart is beating fast",
  "iats_2" = "Most of the time my attention is focused on whether I am hungry",
  "iats_3" = "Most of the time my attention is focused on whether I am breathing fast",
  "iats_4" = "Most of the time my attention is focused on whether I am thirsty",
  "iats_5" = "Most of the time my attention is focused on whether I need to urinate",
  "iats_6" = "Most of the time my attention is focused on whether I need to defecate",
  "iats_7" = "Most of the time my attention is focused on different tastes when eating",
  "iats_8" = "Most of the time my attention is focused on whether I am nauseated",
  "iats_9" = "Most of the time my attention is focused on whether I need to sneeze",
  "iats_10" = "Most of the time my attention is focused on whether I need to cough",
  "iats_11" = "Most of the time my attention is focused on my body temperature",
  "iats_12" = "Most of the time my attention is focused on whether I am sexually aroused",
  "iats_13" = "Most of the time my attention is focused on whether I need to pass wind",
  "iats_14" = "Most of the time my attention is focused on whether I need to burp",
  "iats_15" = "Most of the time my attention is focused on whether my muscles are tired/sore",
  "iats_16" = "Most of the time my attention is focused on pain after injury",
  "iats_17" = "Most of the time my attention is focused on pain NOT caused by injury",
  "iats_18" = "Most of the time my attention is focused on whether my blood sugar is low",
  "iats_19" = "Most of the time my attention is focused on whether touch is pleasant",
  "iats_20" = "Most of the time my attention is focused on whether touch feels ticklish",
  "iats_21" = "Most of the time my attention is focused on whether my body feels itchy"
)

all_full_text <- c(ias_full_text, iats_full_text)

scale_names <- c(
  "IAS" = "Interoceptive Accuracy Scale (IAS)",
  "IATS" = "Interoceptive Attention Scale (IATS)"
)

# =============================================================================
# 3. FIGURE A: ITEM IMPORTANCE WITH FULL TEXT
# =============================================================================

cat("Creating Panel A: Item importance...\n")

plot_df <- consensus_df %>%
  mutate(
    importance = (43 - avg_rank) / 42,
    full_text = all_full_text[item],
    scale_full = ifelse(scale == "IAS", scale_names["IAS"], scale_names["IATS"])
  ) %>%
  arrange(desc(importance))

# Take top 21 items
plot_df_top <- head(plot_df, 21)

# Create labels with item codes + full text
plot_df_top$item_code <- toupper(gsub("_", " ", plot_df_top$item))
plot_df_top$wrapped_text <- sapply(seq_len(nrow(plot_df_top)), function(i) {
  code <- plot_df_top$item_code[i]
  text <- plot_df_top$full_text[i]
  paste0(code, ": ", text)
})

scale_colors <- c(
  "Interoceptive Accuracy Scale (IAS)" = "#2E86AB",
  "Interoceptive Attention Scale (IATS)" = "#A23B72"
)

cutoff_importance <- 0.65

n_ias_above <- sum(plot_df_top$scale == "IAS" & plot_df_top$importance > cutoff_importance)
n_iats_above <- sum(plot_df_top$scale == "IATS" & plot_df_top$importance > cutoff_importance)
n_total_above <- n_ias_above + n_iats_above

p_importance <- ggplot(plot_df_top, aes(x = reorder(wrapped_text, importance),
                                         y = importance, fill = scale_full)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_hline(yintercept = cutoff_importance, linetype = "dashed",
             color = "#E74C3C", linewidth = 0.8) +
  annotate("text", x = 2, y = cutoff_importance + 0.012,
           label = paste0("10-item knee: ", n_total_above, " items (",
                          n_ias_above, " IAS + ", n_iats_above, " IATS)"),
           color = "#E74C3C", size = 5.5, hjust = 0) +
  scale_fill_manual(values = scale_colors, name = "Scale") +
  coord_flip(ylim = c(0.43, 0.95)) +
  labs(
    title = "A. Item Importance for Mental Health Prediction",
    subtitle = "Top 21 items shown (maximum R² solution) | Red line = 10-item knee cutoff",
    x = NULL,
    y = "Relative Importance Score"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 21, hjust = 0.5),
    plot.subtitle = element_text(size = 14, color = "gray40", hjust = 0.5),
    plot.title.position = "plot",
    axis.text.y = element_text(size = 13, lineheight = 0.9),
    axis.text.x = element_text(size = 13),
    axis.title = element_text(size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# 4. FIGURE B: SATURATION CURVE
# =============================================================================

cat("Creating Panel B: Saturation curve...\n")

full_scale_r2 <- config_results$cv_r2[config_results$config == "Full scales (totals)"]

p_saturation <- ggplot(saturation_results, aes(x = n_items, y = r2)) +
  geom_line(color = "#2E86AB", linewidth = 1.2) +
  geom_point(color = "#2E86AB", size = 3) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "#E74C3C",
             linewidth = 0.6, alpha = 0.7) +
  geom_vline(xintercept = 21, linetype = "dashed", color = "#9C27B0",
             linewidth = 0.6, alpha = 0.7) +
  geom_hline(yintercept = full_scale_r2, linetype = "dashed",
             color = "#27AE60", linewidth = 0.8) +
  annotate("text", x = 38, y = full_scale_r2 + 0.025,
           label = "Traditional:\nsum scores", color = "#27AE60", size = 5.5,
           lineheight = 0.9) +
  annotate("label", x = 10, y = 0.27, label = "10 items\n(knee)",
           color = "#E74C3C", fill = "white", size = 5.5, hjust = 0.5,
           lineheight = 0.9, label.size = 0.3, fontface = "bold") +
  annotate("label", x = 21, y = 0.27, label = "21 items\n(maximum)",
           color = "#9C27B0", fill = "white", size = 5.5, hjust = 0.5,
           lineheight = 0.9, label.size = 0.3, fontface = "bold") +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 10, 15, 21, 30, 42)) +
  labs(
    title = "B. Saturation Curve",
    subtitle = "Cross-validated R² for mental health composite",
    x = "Number of Items",
    y = "Cross-Validated R²"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 21),
    plot.subtitle = element_text(size = 14, color = "gray40"),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# 5. FIGURE C: CENTRALITY vs PREDICTION SCATTER
# =============================================================================

cat("Creating Panel C: Centrality vs Prediction...\n")

scatter_df <- consensus_df %>%
  left_join(centrality_df, by = "item") %>%
  mutate(
    prediction_importance = (43 - avg_rank) / 42,
    scale_full = ifelse(scale == "IAS", "Interoceptive Accuracy Scale (IAS)",
                        "Interoceptive Attention Scale (IATS)"),
    short_label = toupper(gsub("_", "", item))
  ) %>%
  group_by(scale) %>%
  mutate(
    network_centrality_norm = (strength - min(strength, na.rm=T)) /
      (max(strength, na.rm=T) - min(strength, na.rm=T))
  ) %>%
  ungroup() %>%
  mutate(
    quadrant = case_when(
      prediction_importance > cutoff_importance & network_centrality_norm < 0.5 ~ "Peripheral but Predictive",
      prediction_importance > cutoff_importance & network_centrality_norm >= 0.5 ~ "Central AND Predictive",
      prediction_importance <= cutoff_importance & network_centrality_norm >= 0.5 ~ "Central but Less Predictive",
      TRUE ~ "Peripheral and Less Predictive"
    )
  )

key_predictors <- c("ias_19", "iats_17", "iats_21", "iats_3")
peripheral_predictive_iats <- scatter_df %>%
  filter(scale == "IATS", network_centrality_norm < 0.5, prediction_importance > 0.75) %>%
  pull(item)
central_ias_items <- scatter_df %>%
  filter(scale == "IAS", network_centrality_norm > 0.7, prediction_importance < cutoff_importance,
         prediction_importance > 0.4, item != "ias_17") %>%
  pull(item)

key_items <- c(key_predictors, peripheral_predictive_iats, central_ias_items)
scatter_df$is_key <- scatter_df$item %in% key_items

centrality_cutoff <- 0.5
moderate_importance <- 0.45

p_scatter <- ggplot(scatter_df, aes(x = network_centrality_norm, y = prediction_importance)) +
  annotate("rect", xmin = -0.1, xmax = centrality_cutoff, ymin = cutoff_importance, ymax = 1.05,
           fill = "#FFCDD2", alpha = 0.6) +
  annotate("rect", xmin = centrality_cutoff, xmax = 1.1, ymin = cutoff_importance, ymax = 1.05,
           fill = "#C8E6C9", alpha = 0.6) +
  annotate("rect", xmin = -0.1, xmax = centrality_cutoff, ymin = moderate_importance, ymax = cutoff_importance,
           fill = "#FFEBEE", alpha = 0.4) +
  annotate("rect", xmin = centrality_cutoff, xmax = 1.1, ymin = moderate_importance, ymax = cutoff_importance,
           fill = "#E8F5E9", alpha = 0.4) +
  geom_point(aes(color = scale_full, size = prediction_importance), alpha = 0.7) +
  geom_label_repel(
    data = filter(scatter_df, is_key),
    aes(label = short_label),
    size = 5, box.padding = 0.5, point.padding = 0.3,
    segment.color = "gray50", min.segment.length = 0, max.overlaps = 20
  ) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "gray50") +
  geom_hline(yintercept = cutoff_importance, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = scale_colors, name = "Scale") +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  labs(
    title = "C. Network Centrality vs. Predictive Importance",
    subtitle = "Centrality normalized within each scale",
    x = "Network Centrality (within-scale normalized)",
    y = "Predictive Importance"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 21),
    plot.subtitle = element_text(size = 14, color = "gray40"),
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(0.0, 1.05))

# =============================================================================
# 6. FIGURE D: NETWORK VISUALIZATION
# =============================================================================

cat("Creating Panel D: Network...\n")

ias_items <- paste0("ias_", 1:21)
iats_items <- paste0("iats_", 1:21)
all_items <- c(ias_items, iats_items)
X <- as.matrix(dfc[, all_items])

cor_matrix <- cor(X, use = "pairwise.complete.obs")
adj_matrix <- abs(cor_matrix)
threshold <- 0.35
adj_matrix[adj_matrix < threshold] <- 0
diag(adj_matrix) <- 0

g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

set.seed(123)
communities <- cluster_louvain(g, weights = E(g)$weight)
membership <- membership(communities)

importance_map <- setNames(scatter_df$prediction_importance, scatter_df$item)

set.seed(42)
fr_layout <- layout_with_fr(g, niter = 500)

layout_df <- data.frame(
  x = fr_layout[, 1],
  y = fr_layout[, 2],
  item = V(g)$name,
  community = as.factor(membership),
  importance = importance_map[V(g)$name],
  scale = ifelse(grepl("^ias", V(g)$name), "IAS", "IATS")
)

layout_df$label <- sapply(layout_df$item, function(x) {
  if (grepl("^ias", x)) paste0("IAS", gsub("ias_", "", x))
  else paste0("IATS", gsub("iats_", "", x))
})

top_21_items <- head(plot_df$item, 21)
top_10_items <- head(plot_df$item, 10)
layout_df$is_top_21 <- layout_df$item %in% top_21_items
layout_df$is_top_10 <- layout_df$item %in% top_10_items

edge_list <- as_edgelist(g)
edge_df <- data.frame(
  x = layout_df$x[match(edge_list[,1], layout_df$item)],
  y = layout_df$y[match(edge_list[,1], layout_df$item)],
  xend = layout_df$x[match(edge_list[,2], layout_df$item)],
  yend = layout_df$y[match(edge_list[,2], layout_df$item)],
  weight = E(g)$weight
)

community_colors <- c(
  "1" = "#7CB342", "2" = "#F4A460", "3" = "#6A5ACD",
  "4" = "#FFD700", "5" = "#FF6B6B", "6" = "#00CED1"
)

community_names <- c(
  "1" = "Core IAS (11 items)",
  "2" = "Autonomic reflexes (IAS)",
  "3" = "Sensory sensitivity",
  "4" = "IAS19 alone (isolated)",
  "5" = "Arousal IATS (9 items)",
  "6" = "General IATS (11 items)"
)

ias19_pos <- layout_df[layout_df$item == "ias_19", ]
y_range <- diff(range(layout_df$y))

n_ias_10 <- sum(head(plot_df, 10)$scale == "IAS")
n_iats_10 <- sum(head(plot_df, 10)$scale == "IATS")

p_network <- ggplot() +
  geom_segment(data = edge_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "gray70", alpha = 0.4, linewidth = edge_df$weight * 1.5) +
  geom_point(data = filter(layout_df, is_top_10),
             aes(x = x, y = y, size = importance * 1.25),
             shape = 21, fill = NA, color = "#E74C3C", stroke = 2.5) +
  geom_point(data = layout_df,
             aes(x = x, y = y, fill = community, size = importance),
             shape = 21, color = "white", stroke = 0.5) +
  geom_point(data = filter(layout_df, is_top_21),
             aes(x = x, y = y, size = importance),
             shape = 21, fill = NA, color = "black", stroke = 1.5) +
  geom_text(data = layout_df,
            aes(x = x, y = y, label = label),
            size = 5.5, fontface = "bold") +
  annotate("label",
           x = ias19_pos$x + y_range * 0.15,
           y = ias19_pos$y - y_range * 0.08,
           label = "IAS19: Peripheral\nbut Top Predictor",
           size = 5, fill = "#FFE4E1", color = "#C0392B",
           label.size = 0.3, lineheight = 0.9) +
  annotate("segment",
           x = ias19_pos$x + y_range * 0.08,
           y = ias19_pos$y - y_range * 0.04,
           xend = ias19_pos$x + y_range * 0.02,
           yend = ias19_pos$y - y_range * 0.01,
           arrow = arrow(length = unit(0.15, "cm")), color = "#C0392B") +
  scale_fill_manual(values = community_colors, labels = community_names, name = "Community") +
  scale_size_continuous(range = c(8, 18), guide = "none") +
  coord_cartesian(
    xlim = c(min(layout_df$x) - y_range * 0.1, max(layout_df$x) + y_range * 0.1),
    ylim = c(min(layout_df$y) - y_range * 0.15, max(layout_df$y) + y_range * 0.15),
    clip = "off"
  ) +
  labs(
    title = "D. Item Network Structure",
    subtitle = paste0("Node size = predictive importance | Red border = top 10 items (", n_ias_10, " IAS + ", n_iats_10, " IATS) | Black border = top 21 items\n",
                      "Colors = network communities (detected via Louvain algorithm)")
  ) +
  theme_void(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 21, hjust = 0.5),
    plot.subtitle = element_text(size = 14, color = "gray40", hjust = 0.5, lineheight = 1.2),
    legend.position = "bottom",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_legend(nrow = 2))

# =============================================================================
# 7. SAVE INDIVIDUAL FIGURES
# =============================================================================

ggsave(file.path(plots_dir, "sub_plots/fig2_A_item_importance.png"), p_importance,
       width = 10, height = 8, dpi = 300)

ggsave(file.path(plots_dir, "sub_plots/fig2_B_saturation_curve.png"), p_saturation,
       width = 7, height = 5, dpi = 300)

ggsave(file.path(plots_dir, "sub_plots/fig2_C_centrality_vs_prediction.png"), p_scatter,
       width = 8, height = 7, dpi = 300)

# Save community network to supplementary plots
ggsave(file.path(plots_dir, "supplementary_plots/2_s_2_network_communities.png"), p_network,
       width = 10, height = 9, dpi = 300)

cat("All individual figures saved.\n")

# =============================================================================
# 8. COMBINED FIGURE (using MDS network with mental health)
# =============================================================================

cat("Creating combined publication figure...\n")

# Load the MDS figure with mental health as Panel D
library(png)
mds_img <- readPNG(file.path(plots_dir, "sub_plots/fig2_D_mds_network.png"))
p_mds <- rasterGrob(mds_img, interpolate = TRUE)

# Layout:
# Top row: A (Item importance) - full width
# Middle row: B (Saturation) + C (Centrality vs Prediction)
# Bottom row: D (MDS Network with MH) - full width

middle_row <- arrangeGrob(p_saturation, p_scatter, ncol = 2)

final_figure <- arrangeGrob(
  p_importance,
  middle_row,
  p_mds,
  ncol = 1,
  heights = c(1.2, 0.9, 1.4)
)

ggsave(file.path(plots_dir, "figure_2_minimal_subtyping_network.png"), final_figure,
       width = 16, height = 18, dpi = 300)

cat("Combined figure saved: figure_2_minimal_subtyping_network.png\n")

# =============================================================================
# 9. SAVE TOP 10 ITEMS WITH FULL TEXT
# =============================================================================

top_10_balanced <- rbind(
  head(filter(plot_df, scale == "IAS"), 5),
  head(filter(plot_df, scale == "IATS"), 5)
) %>%
  arrange(desc(importance)) %>%
  select(item, scale, full_text, importance)

write.csv(top_10_balanced, file.path(supplementary_data_dir, "top_10_items_full_text.csv"),
          row.names = FALSE)

cat("\nTop 10 items (5 IAS + 5 IATS) saved with full questionnaire text.\n")

# =============================================================================
# 10. PRINT SUMMARY
# =============================================================================

cat("\n========================================\n")
cat("PUBLICATION FIGURES GENERATED\n")
cat("========================================\n\n")

cat("Individual panels (in sub_plots/):\n")
cat("  - fig2_A_item_importance.png (Panel A)\n")
cat("  - fig2_B_saturation_curve.png (Panel B)\n")
cat("  - fig2_C_centrality_vs_prediction.png (Panel C)\n")
cat("  - fig2_D_mds_network.png (Panel D)\n")

cat("\nCombined figure:\n")
cat("  - figure_2_minimal_subtyping_network.png (4-panel publication figure)\n")

cat("\nSupplementary data:\n")
cat("  - top_10_items_full_text.csv\n")

cat("\nOutputs: supplementary_data/, plots/\n")
cat("\nVisualization complete.\n")
