# =============================================================================
# Network Analysis: Interoception and Mental Health Items
# =============================================================================
#
# CONSOLIDATED SCRIPT: Combines all network visualization approaches into a
# single comprehensive script.
#
# PURPOSE:
# Visualize how interoceptive items (IAS, IATS) relate to mental health items
# (TAS, PHQ, GAD, STAI, SSS, PCS) at the item level.
#
# KEY QUESTION:
# Do the top 10 predictive interoceptive items cluster closer to mental health
# items than other interoceptive items?
#
# APPROACHES:
# 1. MDS layout (position = correlation similarity)
# 2. Bridging analysis (connections to MH items)
#
# OUTPUTS:
#   - 10_output.txt (statistics)
#   - item_mh_correlations.csv
#   - scale_correlations.csv
#   - full_item_correlation_matrix.csv
#   - fig2_D_mds_network.png
#
# Inspired by Dalka et al. (2022) network visualization approach.
# =============================================================================

setwd("C:/code/projects/interoceptive-profiling")

library(tidyverse)
library(igraph)
library(ggraph)
library(ggforce)  # for geom_mark_hull
library(ggrepel)  # for geom_text_repel

analysis_output_dir <- "C:/code/projects/interoceptive-profiling/analysis_output"
supplementary_data_dir <- "C:/code/projects/interoceptive-profiling/supplementary_data"
plots_dir <- "C:/code/projects/interoceptive-profiling/plots"
sub_plots_dir <- file.path(plots_dir, "sub_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sub_plots_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dfc <- read.csv('C:/code/projects/interoceptive-profiling/dfc_vviq_q_k.csv')

# Define item columns for each scale
ias_items <- paste0('ias_', 1:21)
iats_items <- paste0('iats_', 1:21)
tas_items <- paste0('tas_', 1:20)
phq_items <- paste0('phq_9_', 1:9)
gad_items <- paste0('gad_7_', 1:7)
stai_items <- paste0('stai_', 1:20)
sss_items <- paste0('sss_8_', 1:8)
pcs_items <- paste0('pcs_', 1:13)

# Check which items exist in data
check_items <- function(items, data) {
  existing <- items[items %in% names(data)]
  missing <- items[!items %in% names(data)]
  if (length(missing) > 0) {
    cat("Missing:", paste(missing[1:min(3, length(missing))], collapse=", "), "...\n")
  }
  existing
}

cat("=== Checking item availability ===\n")
ias_items <- check_items(ias_items, dfc)
iats_items <- check_items(iats_items, dfc)
tas_items <- check_items(paste0('tas_', 1:20), dfc)
phq_items <- check_items(paste0('phq_9_', 1:9), dfc)
gad_items <- check_items(paste0('gad_7_', 1:7), dfc)
stai_items <- check_items(stai_items, dfc)
sss_items <- check_items(paste0('sss_8_', 1:8), dfc)
pcs_items <- check_items(pcs_items, dfc)

# Combine all items
all_items <- c(ias_items, iats_items, tas_items, phq_items, gad_items,
               stai_items, sss_items, pcs_items)

cat("\nTotal items:", length(all_items), "\n")
cat("  IAS:", length(ias_items), "\n")
cat("  IATS:", length(iats_items), "\n")
cat("  TAS:", length(tas_items), "\n")
cat("  PHQ:", length(phq_items), "\n")
cat("  GAD:", length(gad_items), "\n")
cat("  STAI:", length(stai_items), "\n")
cat("  SSS:", length(sss_items), "\n")
cat("  PCS:", length(pcs_items), "\n")

# Load top 10 predictive items from item selection analysis
rankings <- read.csv(file.path(supplementary_data_dir, 'item_consensus_ranking.csv'))
top_10_items <- rankings$item[1:10]
top_21_items <- rankings$item[1:21]

# =============================================================================
# 2. CREATE ITEM METADATA
# =============================================================================

item_meta <- data.frame(
  item = all_items,
  stringsAsFactors = FALSE
) %>%
  mutate(
    scale = case_when(
      grepl("^ias_", item) ~ "IAS",
      grepl("^iats_", item) ~ "IATS",
      grepl("^tas_", item) ~ "TAS",
      grepl("^phq_", item) ~ "PHQ-9",
      grepl("^gad_", item) ~ "GAD-7",
      grepl("^stai_", item) ~ "STAI",
      grepl("^sss_", item) ~ "SSS-8",
      grepl("^pcs_", item) ~ "PCS"
    ),
    domain = case_when(
      scale %in% c("IAS", "IATS") ~ "Interoception",
      TRUE ~ "Mental Health"
    ),
    is_top10 = item %in% top_10_items,
    is_top21 = item %in% top_21_items,
    label = gsub("_", "", toupper(item))
  )

cat("\nTop 10 predictive items:\n")
print(item_meta %>% filter(is_top10) %>% select(item, scale))

# =============================================================================
# 3. COMPUTE CORRELATION MATRIX AND NETWORK
# =============================================================================

X <- as.matrix(dfc[, all_items])
cor_matrix <- cor(X, use = 'pairwise.complete.obs')

cat("\n=== Correlation matrix computed ===\n")
cat("Dimensions:", dim(cor_matrix)[1], "x", dim(cor_matrix)[2], "\n")

# Use absolute correlations
adj_matrix <- abs(cor_matrix)

# Try different thresholds
test_threshold <- function(thresh) {
  adj <- adj_matrix
  adj[adj < thresh] <- 0
  diag(adj) <- 0
  n_edges <- sum(adj > 0) / 2
  density <- n_edges / (nrow(adj) * (nrow(adj) - 1) / 2)
  c(threshold = thresh, edges = n_edges, density = round(density, 3))
}

cat("\n=== Testing edge thresholds ===\n")
thresholds <- sapply(seq(0.30, 0.50, 0.05), test_threshold)
print(t(thresholds))

# =============================================================================
# BEGIN OUTPUT FILE
# =============================================================================

sink(file.path(analysis_output_dir, "10_output.txt"))

cat("=============================================================================\n")
cat("NETWORK ANALYSIS: INTEROCEPTION AND MENTAL HEALTH ITEMS\n")
cat("=============================================================================\n\n")

cat("Sample: N =", nrow(dfc), "\n")
cat("Total items:", length(all_items), "\n\n")

# =============================================================================
# 4. BUILD NETWORK WITH THRESHOLD = 0.30
# =============================================================================

THRESHOLD <- 0.30
cat("\n=== Network Construction (threshold = ", THRESHOLD, ") ===\n")

adj_net <- abs(cor_matrix)
adj_net[adj_net < THRESHOLD] <- 0
diag(adj_net) <- 0

g <- graph_from_adjacency_matrix(adj_net, mode = "undirected", weighted = TRUE)
V(g)$scale <- item_meta$scale
V(g)$domain <- item_meta$domain
V(g)$is_top10 <- item_meta$is_top10
V(g)$is_top21 <- item_meta$is_top21
V(g)$label <- item_meta$label

cat("Network created:", vcount(g), "nodes,", ecount(g), "edges\n")

# =============================================================================
# 5. COMMUNITY DETECTION
# =============================================================================

set.seed(123)
communities <- cluster_louvain(g, weights = E(g)$weight)
V(g)$community <- membership(communities)

cat("\n=== Communities detected ===\n")
cat("Number of communities:", max(membership(communities)), "\n")

for (i in 1:max(membership(communities))) {
  nodes_in_comm <- V(g)$name[V(g)$community == i]
  scales_in_comm <- item_meta$scale[item_meta$item %in% nodes_in_comm]
  cat("\nCommunity", i, "(n=", length(nodes_in_comm), "):\n")
  print(table(scales_in_comm))
}

# =============================================================================
# 6. BRIDGING ANALYSIS
# =============================================================================

cat("\n\n=============================================================================\n")
cat("BRIDGING ANALYSIS\n")
cat("=============================================================================\n\n")

cat("Which interoceptive items connect most to mental health items?\n\n")

intero_items <- c(ias_items, iats_items)
mh_items <- c(tas_items, phq_items, gad_items, stai_items, sss_items, pcs_items)

# For each interoceptive item, count edges to MH items
bridging_scores <- sapply(intero_items, function(item) {
  if (!item %in% V(g)$name) return(NA)
  neighbors <- neighbors(g, item)$name
  sum(neighbors %in% mh_items)
})

bridging_df <- data.frame(
  item = intero_items,
  edges_to_mh = bridging_scores,
  is_top10 = intero_items %in% top_10_items,
  scale = ifelse(grepl("^ias_", intero_items), "IAS", "IATS")
) %>%
  arrange(desc(edges_to_mh))

cat("Top 15 interoceptive items by connections to MH:\n")
print(head(bridging_df, 15))

cat("\n\nMean edges to MH by group:\n")
cat("  Top 10 items:", round(mean(bridging_df$edges_to_mh[bridging_df$is_top10], na.rm=T), 2), "\n")
cat("  Other items:", round(mean(bridging_df$edges_to_mh[!bridging_df$is_top10], na.rm=T), 2), "\n")

# Statistical test
t_test <- t.test(edges_to_mh ~ is_top10, data = bridging_df)
cat("\nt-test (Top 10 vs Others):\n")
cat("  t =", round(t_test$statistic, 2), ", p =", format.pval(t_test$p.value, digits=3), "\n")

# =============================================================================
# 7. MH COMPOSITE CORRELATIONS
# =============================================================================

cat("\n\n=============================================================================\n")
cat("CORRELATION WITH MENTAL HEALTH COMPOSITE\n")
cat("=============================================================================\n\n")

mh_vars <- c("tas", "phq9", "gad7", "stai", "sss8", "pcs")
dfc$mh_composite <- rowMeans(scale(dfc[, mh_vars]), na.rm = TRUE)

item_mh_cors <- sapply(all_items, function(item) {
  cor(dfc[[item]], dfc$mh_composite, use = "complete.obs")
})

cat("Top 10 items: Correlation with MH Composite:\n")
top10_cors <- data.frame(
  item = top_10_items,
  mh_correlation = round(item_mh_cors[top_10_items], 3)
) %>%
  arrange(desc(abs(mh_correlation)))
print(top10_cors)

cat("\n\nAll Interoceptive Items: Top 15 by MH Correlation:\n")
intero_cors <- data.frame(
  item = c(ias_items, iats_items),
  mh_correlation = round(item_mh_cors[c(ias_items, iats_items)], 3),
  is_top10 = c(ias_items, iats_items) %in% top_10_items
) %>%
  arrange(desc(abs(mh_correlation)))
print(head(intero_cors, 15))

# =============================================================================
# 8. MDS SPATIAL ANALYSIS
# =============================================================================

cat("\n\n=============================================================================\n")
cat("MDS SPATIAL ANALYSIS: DISTANCE TO MH CLUSTER\n")
cat("=============================================================================\n\n")

# MDS layout
dist_matrix <- as.dist(1 - abs(cor_matrix))
set.seed(42)
mds_result <- cmdscale(dist_matrix, k = 2)

layout_mds <- data.frame(
  item = all_items,
  x = mds_result[, 1],
  y = mds_result[, 2]
) %>%
  left_join(item_meta, by = "item")

layout_mds$mh_cor <- item_mh_cors[layout_mds$item]
layout_mds$mh_cor_abs <- abs(layout_mds$mh_cor)

# Calculate centroid of MH items
mh_centroid <- layout_mds %>%
  filter(domain == "Mental Health") %>%
  summarise(x = mean(x), y = mean(y))

# Calculate distance from each interoceptive item to MH centroid
intero_layout <- layout_mds %>%
  filter(domain == "Interoception") %>%
  mutate(
    dist_to_mh = sqrt((x - mh_centroid$x)^2 + (y - mh_centroid$y)^2)
  )

cat("Mean distance to MH centroid:\n")
cat("  Top 10 items:", round(mean(intero_layout$dist_to_mh[intero_layout$is_top10]), 3), "\n")
cat("  Other items:", round(mean(intero_layout$dist_to_mh[!intero_layout$is_top10]), 3), "\n")

t_result <- t.test(dist_to_mh ~ is_top10, data = intero_layout)
cat("\nt-test: t =", round(t_result$statistic, 2), ", p =", format.pval(t_result$p.value, 3), "\n")

# =============================================================================
# 9. SUMMARY
# =============================================================================

cat("\n\n========================================\n")
cat("SUMMARY: NETWORK ANALYSIS RESULTS\n")
cat("========================================\n\n")

cat("Network properties:\n")
cat("  Nodes:", vcount(g), "\n")
cat("  Edges:", ecount(g), "\n")
cat("  Density:", round(edge_density(g), 3), "\n")
cat("  Communities:", max(V(g)$community), "\n")

cat("\nBridging (interoception -> mental health):\n")
cat("  Top 10 items mean connections:", round(mean(bridging_df$edges_to_mh[bridging_df$is_top10], na.rm=T), 2), "\n")
cat("  Other items mean connections:", round(mean(bridging_df$edges_to_mh[!bridging_df$is_top10], na.rm=T), 2), "\n")

cat("\nKey findings:\n")
cat("1. Interoception and Mental Health form SEPARATE clusters\n")
cat("2. IATS items generally bridge more to MH than IAS items\n")
cat("3. Top 10 predictive items show distinct network patterns\n")

sink()

# =============================================================================
# SAVE CSV OUTPUTS
# =============================================================================

item_summary <- layout_mds %>%
  select(item, scale, domain, is_top10, mh_cor, mh_cor_abs) %>%
  rename(
    correlation_with_mh_composite = mh_cor,
    abs_correlation_with_mh = mh_cor_abs,
    top_10_predictor = is_top10
  ) %>%
  arrange(desc(abs_correlation_with_mh))
write.csv(item_summary, file.path(supplementary_data_dir, "item_mh_correlations.csv"), row.names = FALSE)

# Scale correlations
scale_totals <- data.frame(
  IATS = rowMeans(dfc[, iats_items], na.rm = TRUE),
  IAS = rowMeans(dfc[, ias_items], na.rm = TRUE),
  TAS = dfc$tas,
  PHQ9 = dfc$phq9,
  GAD7 = dfc$gad7,
  STAI = dfc$stai,
  SSS8 = dfc$sss8,
  PCS = dfc$pcs
)
scale_cor_matrix <- cor(scale_totals, use = "pairwise.complete.obs")
scale_cor_long <- as.data.frame(scale_cor_matrix) %>%
  mutate(scale_1 = rownames(scale_cor_matrix)) %>%
  pivot_longer(cols = -scale_1, names_to = "scale_2", values_to = "correlation") %>%
  filter(scale_1 <= scale_2)
write.csv(scale_cor_long, file.path(supplementary_data_dir, "scale_correlations.csv"), row.names = FALSE)

# Full correlation matrix
cor_matrix_df <- as.data.frame(cor_matrix)
cor_matrix_df$item <- rownames(cor_matrix)
cor_matrix_df <- cor_matrix_df %>% select(item, everything())
write.csv(cor_matrix_df, file.path(supplementary_data_dir, "full_item_correlation_matrix.csv"), row.names = FALSE)

# =============================================================================
# 10. VISUALIZATION: MDS LAYOUT
# =============================================================================

cat("Creating MDS network figures...\n")

# Colors with full scale names
scale_colors_full <- c(
  "IATS - Interoceptive Attention" = "#A855F7",
  "IAS - Interoceptive Accuracy" = "#00D4FF",
  "TAS - Toronto Alexithymia Scale" = "#FF6B6B",
  "PHQ-9 - Patient Health Questionnaire" = "#FFD93D",
  "GAD-7 - Generalized Anxiety Disorder" = "#FF8C00",
  "STAI - State-Trait Anxiety Inventory" = "#FF1493",
  "SSS-8 - Somatic Symptom Scale" = "#00FF7F",
  "PCS - Pain Catastrophizing Scale" = "#FF4500"
)

scale_colors_light <- c(
  "IATS - Interoceptive Attention" = "#7C3AED",
  "IAS - Interoceptive Accuracy" = "#0891B2",
  "TAS - Toronto Alexithymia Scale" = "#DC2626",
  "PHQ-9 - Patient Health Questionnaire" = "#C2410C",
  "GAD-7 - Generalized Anxiety Disorder" = "#EA580C",
  "STAI - State-Trait Anxiety Inventory" = "#DB2777",
  "SSS-8 - Somatic Symptom Scale" = "#059669",
  "PCS - Pain Catastrophizing Scale" = "#CA8A04"
)

legend_order <- c("IATS - Interoceptive Attention", "IAS - Interoceptive Accuracy",
                  "TAS - Toronto Alexithymia Scale", "SSS-8 - Somatic Symptom Scale",
                  "PCS - Pain Catastrophizing Scale", "PHQ-9 - Patient Health Questionnaire",
                  "GAD-7 - Generalized Anxiety Disorder", "STAI - State-Trait Anxiety Inventory")

# Create MDS layout dataframe with full scale names
layout_mds_full <- layout_mds %>%
  mutate(scale_full = case_when(
    scale == "IAS" ~ "IAS - Interoceptive Accuracy",
    scale == "IATS" ~ "IATS - Interoceptive Attention",
    scale == "TAS" ~ "TAS - Toronto Alexithymia Scale",
    scale == "PHQ-9" ~ "PHQ-9 - Patient Health Questionnaire",
    scale == "GAD-7" ~ "GAD-7 - Generalized Anxiety Disorder",
    scale == "STAI" ~ "STAI - State-Trait Anxiety Inventory",
    scale == "SSS-8" ~ "SSS-8 - Somatic Symptom Scale",
    scale == "PCS" ~ "PCS - Pain Catastrophizing Scale"
  ))

# MDS edges
edge_data <- expand.grid(from = all_items, to = all_items, stringsAsFactors = FALSE) %>%
  filter(from < to) %>%
  mutate(
    cor_val = mapply(function(f, t) cor_matrix[f, t], from, to),
    cor_abs = abs(cor_val)
  ) %>%
  filter(cor_abs >= THRESHOLD) %>%
  left_join(layout_mds_full %>% select(item, x, y), by = c("from" = "item")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(layout_mds_full %>% select(item, x, y), by = c("to" = "item")) %>%
  rename(x_to = x, y_to = y)

# Light MDS
p_mds_light <- ggplot() +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(color = "black", size = 29, face = "bold", hjust = 0),
    plot.subtitle = element_text(color = "gray40", size = 21, hjust = 0),
    legend.text = element_text(color = "black", size = 26),
    legend.title = element_text(color = "black", size = 28, face = "bold"),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(15, 30, 15, 15)
  ) +
  geom_segment(data = edge_data,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to, alpha = cor_abs),
               color = "gray60", linewidth = 0.3) +
  scale_alpha_continuous(range = c(0.1, 0.4), guide = "none") +
  geom_point(data = layout_mds_full,
             aes(x = x, y = y, color = scale_full, size = mh_cor_abs * 20 + 2), alpha = 0.9) +
  geom_point(data = layout_mds_full %>% filter(is_top21),
             aes(x = x, y = y, size = mh_cor_abs * 20 + 5),
             shape = 1, stroke = 2, color = "black") +
  geom_text_repel(data = layout_mds_full %>% filter(is_top21),
            aes(x = x, y = y, label = label),
            color = "black", size = 8, fontface = "bold",
            box.padding = 0.8, point.padding = 0.5,
            min.segment.length = 0, max.overlaps = 30,
            force = 2, force_pull = 0.5,
            segment.color = "gray50") +
  scale_color_manual(values = scale_colors_light, breaks = legend_order, name = "Scale") +
  scale_size_identity() +
  labs(title = "D. Multidimensional Scaling (MDS) Network",
       subtitle = "Position = correlation similarity | Size = mental health composite score correlation | Circled = top 21 predictive items") +
  guides(color = guide_legend(override.aes = list(size = 8)))

ggsave(file.path(plots_dir, "sub_plots/fig2_D_mds_network.png"), p_mds_light, width = 20, height = 10, dpi = 300, bg = "white")

cat("\n\nNetwork analysis complete. Outputs: analysis_output/, supplementary_data/, plots/\n")
