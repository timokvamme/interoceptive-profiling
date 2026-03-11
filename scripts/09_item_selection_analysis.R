# =============================================================================
# Minimal Interoceptive Subtyping: Item Selection Analysis
# =============================================================================
#
# CONSOLIDATED SCRIPT: Combines item importance, saturation, and cluster
# preservation analyses into a single comprehensive script.
#
# RESEARCH QUESTION:
# Can we reduce the 42-item IAS+IATS battery (21 items each) to a minimal set
# while retaining predictive power for mental health outcomes?
#
# APPROACHES:
# 1. LASSO regression - penalized regression for automatic item selection
# 2. Random Forest - non-linear feature importance
# 3. Network centrality - items with highest connectivity (NALS-inspired)
# 4. Incremental R² - saturation curve showing diminishing returns
# 5. Cluster preservation - do reduced item sets preserve 3-cluster structure?
# 6. Cluster flow - who gets reclassified with reduced items?
#
# READS FROM:
#   - Main project data: C:/code/projects/interoceptive-profiling/data_with_k_means_clusters.csv
#   - Original data: C:/code/projects/interoceptive-profiling/dfc_interoception_profiling.csv
#
# OUTPUTS:
#   - 09_output.txt (all statistics)
#   - item_consensus_ranking.csv
#   - saturation_results.csv
#   - minimal_set_comparison.csv
#   - lasso_coefficients.csv
#   - network_centrality.csv
#   - cluster_preservation_summary.csv
#   - bridging_scores.csv
#   - 2_s_4_transition_heatmap.png
#   - 2_s_3_cluster_preservation.png
#   - figure_reclassified.png
#
# =============================================================================

# Set working directory to project root
setwd("C:/code/projects/interoceptive-profiling")

library(glmnet)      # LASSO/Elastic Net
library(caret)       # Cross-validation framework
library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)      # Network analysis
library(cluster)     # Silhouette
library(mclust)      # Adjusted Rand Index
library(patchwork)

# Output directories
analysis_output_dir <- "C:/code/projects/interoceptive-profiling/analysis_output"
supplementary_data_dir <- "C:/code/projects/interoceptive-profiling/supplementary_data"
plots_dir <- "C:/code/projects/interoceptive-profiling/plots"
supplementary_plots_dir <- file.path(plots_dir, "supplementary_plots")
dir.create(supplementary_plots_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 0. LOAD DATA
# =============================================================================

# Load original data with all items
dfc <- read.csv("C:/code/projects/interoceptive-profiling/dfc_interoception_profiling.csv")
cat("Sample size:", nrow(dfc), "\n")

# Define item names
ias_items <- paste0("ias_", 1:21)
iats_items <- paste0("iats_", 1:21)
all_items <- c(ias_items, iats_items)

# Mental health outcomes
outcome_vars <- c("tas", "phq9", "gad7", "stai", "sss8", "pcs")
outcome_labels <- c("Alexithymia", "Depression", "Anxiety", "Trait Anxiety",
                    "Somatic Symptoms", "Pain Catastrophizing")

# Create composite mental health score (z-scored average)
dfc$mh_composite <- rowMeans(scale(dfc[, outcome_vars]), na.rm = TRUE)

# Extract item data
X <- as.matrix(dfc[, all_items])
X <- scale(X)  # Standardize

# =============================================================================
# BEGIN OUTPUT FILE
# =============================================================================

sink(file.path(analysis_output_dir, "09_output.txt"))

cat("=============================================================================\n")
cat("MINIMAL INTEROCEPTIVE SUBTYPING: ITEM SELECTION ANALYSIS\n")
cat("=============================================================================\n\n")

cat("Sample: N =", nrow(dfc), "\n")
cat("IAS items: 21\n")
cat("IATS items: 21\n")
cat("Total items: 42\n\n")

# =============================================================================
# 1. LASSO REGRESSION FOR ITEM SELECTION
# =============================================================================

cat("=============================================================================\n")
cat("APPROACH 1: LASSO REGRESSION\n")
cat("=============================================================================\n\n")

cat("Using LASSO (L1 penalty) to identify items that survive regularization.\n")
cat("Items with non-zero coefficients are most predictive.\n\n")

# LASSO for composite mental health
set.seed(123)
cv_lasso <- cv.glmnet(X, dfc$mh_composite, alpha = 1, nfolds = 10)

# Get coefficients at lambda.1se (more parsimonious)
lasso_coef <- coef(cv_lasso, s = "lambda.1se")
lasso_df <- data.frame(
  item = rownames(lasso_coef)[-1],
  coef = as.numeric(lasso_coef)[-1]
)
lasso_df$abs_coef <- abs(lasso_df$coef)
lasso_df <- lasso_df[order(-lasso_df$abs_coef), ]
lasso_df$selected <- lasso_df$coef != 0

cat("LASSO Results (Mental Health Composite):\n")
cat("Lambda (1se):", cv_lasso$lambda.1se, "\n")
cat("Items with non-zero coefficients:", sum(lasso_df$selected), "\n\n")

cat("Top 10 items by |coefficient|:\n")
top_lasso <- head(lasso_df[lasso_df$selected, ], 10)
for (i in 1:nrow(top_lasso)) {
  cat(sprintf("  %d. %s: %.4f\n", i, top_lasso$item[i], top_lasso$coef[i]))
}

# Run LASSO for each outcome
lasso_results <- list()
for (i in seq_along(outcome_vars)) {
  y <- dfc[[outcome_vars[i]]]
  cv_fit <- cv.glmnet(X, y, alpha = 1, nfolds = 10)
  coefs <- coef(cv_fit, s = "lambda.1se")
  lasso_results[[outcome_vars[i]]] <- data.frame(
    item = rownames(coefs)[-1],
    coef = as.numeric(coefs)[-1]
  )
}

# Count how often each item is selected across outcomes
item_selection_counts <- data.frame(item = all_items)
for (outcome in outcome_vars) {
  item_selection_counts[[outcome]] <- lasso_results[[outcome]]$coef != 0
}
item_selection_counts$total_selections <- rowSums(item_selection_counts[, -1])
item_selection_counts <- item_selection_counts[order(-item_selection_counts$total_selections), ]

cat("\n\nItems selected across multiple outcomes:\n")
for (i in 1:min(15, nrow(item_selection_counts))) {
  if (item_selection_counts$total_selections[i] > 0) {
    cat(sprintf("  %s: selected in %d/6 outcomes\n",
                item_selection_counts$item[i],
                item_selection_counts$total_selections[i]))
  }
}

# =============================================================================
# 2. RANDOM FOREST FEATURE IMPORTANCE
# =============================================================================

cat("\n\n=============================================================================\n")
cat("APPROACH 2: RANDOM FOREST FEATURE IMPORTANCE\n")
cat("=============================================================================\n\n")

if (requireNamespace("randomForest", quietly = TRUE)) {
  library(randomForest)

  set.seed(123)
  rf_model <- randomForest(x = X, y = dfc$mh_composite,
                           ntree = 500, importance = TRUE)

  rf_importance <- data.frame(
    item = rownames(importance(rf_model)),
    importance = importance(rf_model)[, 1]  # %IncMSE
  )
  rf_importance <- rf_importance[order(-rf_importance$importance), ]

  cat("Random Forest Variable Importance (% Inc MSE):\n\n")
  cat("Top 15 items:\n")
  for (i in 1:15) {
    cat(sprintf("  %d. %s: %.2f\n", i, rf_importance$item[i], rf_importance$importance[i]))
  }

} else {
  cat("randomForest package not available. Skipping RF analysis.\n")
  rf_importance <- NULL
}

# =============================================================================
# 3. NETWORK CENTRALITY (NALS-INSPIRED)
# =============================================================================

cat("\n\n=============================================================================\n")
cat("APPROACH 3: NETWORK CENTRALITY (NALS-Inspired)\n")
cat("=============================================================================\n\n")

cat("Building correlation network of item responses.\n")
cat("Items with high centrality are 'hub' items that connect to many others.\n\n")

# Compute correlation matrix
cor_matrix <- cor(X, use = "pairwise.complete.obs")

# Create adjacency matrix (threshold correlations)
adj_matrix <- abs(cor_matrix)
adj_matrix[adj_matrix < 0.3] <- 0
diag(adj_matrix) <- 0

# Create graph
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

# Calculate weighted degree (strength)
strength <- strength(g)

# Calculate betweenness centrality
betweenness_cent <- betweenness(g, weights = 1/E(g)$weight)

# Combine centrality measures
centrality_df <- data.frame(
  item = names(strength),
  strength = strength,
  betweenness = betweenness_cent
)
centrality_df$combined_centrality <- scale(centrality_df$strength) + scale(centrality_df$betweenness)
centrality_df <- centrality_df[order(-centrality_df$combined_centrality), ]

cat("Network Centrality (Weighted Degree + Betweenness):\n\n")
cat("Top 15 most central items:\n")
for (i in 1:15) {
  cat(sprintf("  %d. %s: strength=%.2f, betweenness=%.1f\n",
              i, centrality_df$item[i], centrality_df$strength[i], centrality_df$betweenness[i]))
}

# =============================================================================
# 4. COMBINE RANKINGS INTO CONSENSUS
# =============================================================================

cat("\n\n=============================================================================\n")
cat("CONSENSUS RANKING\n")
cat("=============================================================================\n\n")

# Create combined ranking
consensus_df <- data.frame(item = all_items)

# Add LASSO rank
lasso_rank <- match(consensus_df$item, lasso_df$item)
consensus_df$lasso_rank <- lasso_rank

# Add RF rank (if available)
if (!is.null(rf_importance)) {
  rf_rank <- match(consensus_df$item, rf_importance$item)
  consensus_df$rf_rank <- rf_rank
} else {
  consensus_df$rf_rank <- NA
}

# Add centrality rank
centrality_rank <- match(consensus_df$item, centrality_df$item)
consensus_df$centrality_rank <- centrality_rank

# Add selection count
sel_idx <- match(consensus_df$item, item_selection_counts$item)
consensus_df$selection_count <- item_selection_counts$total_selections[sel_idx]

# Calculate average rank (lower is better)
consensus_df$avg_rank <- rowMeans(consensus_df[, c("lasso_rank", "rf_rank", "centrality_rank")],
                                   na.rm = TRUE)

# Add selection count bonus (items selected in more outcomes get lower rank)
consensus_df$final_score <- consensus_df$avg_rank - (consensus_df$selection_count * 3)
consensus_df <- consensus_df[order(consensus_df$final_score), ]

# Add scale indicator
consensus_df$scale <- ifelse(grepl("^ias_", consensus_df$item), "IAS", "IATS")

cat("Combined consensus ranking (averaging LASSO, RF, and Network methods):\n\n")
cat("Rank  Item        Scale  LASSO  RF   Network  #Outcomes\n")
cat("----  ----------  -----  -----  ---  -------  ---------\n")
for (i in 1:20) {
  cat(sprintf(" %2d   %-10s  %s    %2d     %2d    %2d        %d\n",
              i, consensus_df$item[i], consensus_df$scale[i],
              consensus_df$lasso_rank[i], consensus_df$rf_rank[i],
              consensus_df$centrality_rank[i], consensus_df$selection_count[i]))
}

# =============================================================================
# 5. INCREMENTAL R² ANALYSIS (SATURATION CURVE)
# =============================================================================

cat("\n\n=============================================================================\n")
cat("SATURATION ANALYSIS: HOW MANY ITEMS DO WE NEED?\n")
cat("=============================================================================\n\n")

# Sort items by consensus
top_items <- consensus_df$item

# Function to calculate cross-validated R²
cv_r2 <- function(X_subset, y, nfolds = 5) {
  set.seed(123)
  folds <- createFolds(y, k = nfolds)
  r2_vals <- numeric(nfolds)

  for (i in seq_along(folds)) {
    train_idx <- unlist(folds[-i])
    test_idx <- folds[[i]]

    X_train <- X_subset[train_idx, , drop = FALSE]
    X_test <- X_subset[test_idx, , drop = FALSE]
    y_train <- y[train_idx]
    y_test <- y[test_idx]

    model <- lm(y_train ~ X_train)

    # Predict on test set
    new_data <- cbind(1, X_test)
    pred <- new_data %*% coef(model)

    # Calculate R²
    ss_res <- sum((y_test - pred)^2)
    ss_tot <- sum((y_test - mean(y_test))^2)
    r2_vals[i] <- 1 - ss_res/ss_tot
  }

  return(mean(r2_vals, na.rm = TRUE))
}

# Calculate incremental R² for composite MH
cat("Incremental R² for Mental Health Composite:\n\n")

n_items_to_test <- c(1, 2, 3, 4, 5, 6, 8, 10, 15, 21, 30, 42)
saturation_results <- data.frame()

for (n in n_items_to_test) {
  if (n > length(top_items)) next

  items_subset <- top_items[1:n]
  X_subset <- X[, items_subset, drop = FALSE]

  r2 <- cv_r2(X_subset, dfc$mh_composite)

  saturation_results <- rbind(saturation_results, data.frame(
    n_items = n,
    r2 = r2
  ))

  cat(sprintf("  %2d items: CV R² = %.3f\n", n, r2))
}

# Calculate R² for full scale totals (baseline)
full_scale_r2 <- cv_r2(as.matrix(dfc[, c("ias", "iats")]), dfc$mh_composite)
cat(sprintf("\n  Full IAS+IATS totals: CV R² = %.3f\n", full_scale_r2))

# Find saturation point (90% of max R²)
max_r2 <- max(saturation_results$r2)
saturation_threshold <- 0.90 * max_r2
saturation_point <- min(saturation_results$n_items[saturation_results$r2 >= saturation_threshold])

cat(sprintf("\nSATURATION POINT: %d items achieve 90%% of maximum predictive power\n", saturation_point))

# =============================================================================
# 6. RECOMMENDED MINIMAL ITEM SETS
# =============================================================================

cat("\n\n=============================================================================\n")
cat("RECOMMENDED MINIMAL ITEM SETS\n")
cat("=============================================================================\n\n")

# Get top IAS and IATS items separately
top_ias <- consensus_df$item[consensus_df$scale == "IAS"][1:5]
top_iats <- consensus_df$item[consensus_df$scale == "IATS"][1:5]

cat("TOP IAS ITEMS (by consensus ranking):\n")
for (i in 1:5) {
  cat(sprintf("  %d. %s\n", i, top_ias[i]))
}

cat("\nTOP IATS ITEMS (by consensus ranking):\n")
for (i in 1:5) {
  cat(sprintf("  %d. %s\n", i, top_iats[i]))
}

# Test different minimal configurations
cat("\n\nMINIMAL SET PERFORMANCE COMPARISON:\n\n")

configs <- list(
  "1 IAS + 1 IATS" = c(top_ias[1], top_iats[1]),
  "2 IAS + 2 IATS" = c(top_ias[1:2], top_iats[1:2]),
  "3 IAS + 3 IATS" = c(top_ias[1:3], top_iats[1:3]),
  "4 IAS + 4 IATS" = c(top_ias[1:4], top_iats[1:4]),
  "5 IAS + 5 IATS" = c(top_ias[1:5], top_iats[1:5]),
  "Full scales (totals)" = c("ias", "iats")
)

config_results <- data.frame()

for (config_name in names(configs)) {
  items <- configs[[config_name]]

  if (config_name == "Full scales (totals)") {
    X_config <- as.matrix(dfc[, items])
  } else {
    X_config <- X[, items, drop = FALSE]
  }

  r2 <- cv_r2(X_config, dfc$mh_composite)

  config_results <- rbind(config_results, data.frame(
    config = config_name,
    n_items = length(items),
    cv_r2 = r2
  ))

  cat(sprintf("  %-20s (%d items): CV R² = %.3f\n", config_name, length(items), r2))
}

# Percent of full-scale performance
full_r2 <- config_results$cv_r2[config_results$config == "Full scales (totals)"]
config_results$pct_of_full <- round(config_results$cv_r2 / full_r2 * 100, 1)

cat("\n\nPERCENT OF FULL SCALE PERFORMANCE:\n")
for (i in 1:nrow(config_results)) {
  cat(sprintf("  %-20s: %.1f%% of full scale performance\n",
              config_results$config[i], config_results$pct_of_full[i]))
}

# =============================================================================
# 7. CLUSTER PRESERVATION ANALYSIS
# =============================================================================

cat("\n\n=============================================================================\n")
cat("CLUSTER PRESERVATION ANALYSIS\n")
cat("=============================================================================\n\n")

cat("Question: Do the 10-item and 21-item reduced sets preserve the 3-cluster\n")
cat("structure found with the full 42-item scales?\n\n")

# Define the three item sets
# Top 10 items (by final_score - lower is better)
top_10 <- consensus_df$item[1:10]
top_10_ias <- top_10[grepl("ias", top_10)]
top_10_iats <- top_10[grepl("iats", top_10)]

# Top 21 items
top_21 <- consensus_df$item[1:21]
top_21_ias <- top_21[grepl("ias", top_21)]
top_21_iats <- top_21[grepl("iats", top_21)]

cat("Item Sets:\n")
cat("  Top 10 items: IAS (", length(top_10_ias), "):", paste(top_10_ias, collapse=", "), "\n")
cat("                IATS (", length(top_10_iats), "):", paste(top_10_iats, collapse=", "), "\n")
cat("  Top 21 items: IAS (", length(top_21_ias), "):", paste(top_21_ias, collapse=", "), "\n")
cat("                IATS (", length(top_21_iats), "):", paste(top_21_iats, collapse=", "), "\n\n")

# Create subscores for each item set
dfc$ias_full <- rowMeans(dfc[, ias_items], na.rm = TRUE)
dfc$iats_full <- rowMeans(dfc[, iats_items], na.rm = TRUE)
dfc$ias_10 <- rowMeans(dfc[, top_10_ias], na.rm = TRUE)
dfc$iats_10 <- rowMeans(dfc[, top_10_iats], na.rm = TRUE)
dfc$ias_21 <- rowMeans(dfc[, top_21_ias], na.rm = TRUE)
dfc$iats_21 <- rowMeans(dfc[, top_21_iats], na.rm = TRUE)

# Clustering function
run_clustering <- function(ias_col, iats_col, data, label) {
  ias_z <- scale(data[[ias_col]])
  iats_z <- scale(data[[iats_col]])
  cluster_data <- data.frame(ias = ias_z, iats = iats_z)
  cluster_data <- cluster_data[complete.cases(cluster_data), ]

  set.seed(123)
  km <- kmeans(cluster_data, centers = 3, nstart = 25)

  sil <- silhouette(km$cluster, dist(cluster_data))
  avg_sil <- mean(sil[, 3])

  profiles <- data.frame(
    cluster = 1:3,
    ias_z = km$centers[, 1],
    iats_z = km$centers[, 2],
    n = as.vector(table(km$cluster))
  )

  # Label profiles based on IAS/IATS pattern
  profiles <- profiles %>%
    mutate(
      label = case_when(
        ias_z > 0.3 & iats_z < -0.5 ~ "Efficient",
        iats_z > 0.5 ~ "Hypervigilant",
        ias_z < -0.5 ~ "Uncertain",
        TRUE ~ "Mixed"
      )
    )

  list(
    km = km,
    profiles = profiles,
    silhouette = avg_sil,
    cluster_data = cluster_data,
    cluster_assignments = km$cluster
  )
}

# Run clustering for each item set
result_full <- run_clustering("ias_full", "iats_full", dfc, "Full (42 items)")
result_21 <- run_clustering("ias_21", "iats_21", dfc, "21 items")
result_10 <- run_clustering("ias_10", "iats_10", dfc, "10 items")

cat("--- Full 42 Items ---\n")
cat("Silhouette:", round(result_full$silhouette, 3), "\n")
print(result_full$profiles %>% arrange(desc(ias_z)))

cat("\n--- Top 21 Items ---\n")
cat("Silhouette:", round(result_21$silhouette, 3), "\n")
print(result_21$profiles %>% arrange(desc(ias_z)))

cat("\n--- Top 10 Items ---\n")
cat("Silhouette:", round(result_10$silhouette, 3), "\n")
print(result_10$profiles %>% arrange(desc(ias_z)))

# Compare cluster assignments (Adjusted Rand Index)
comparison <- data.frame(
  full = result_full$cluster_assignments,
  items_21 = result_21$cluster_assignments,
  items_10 = result_10$cluster_assignments
)

ari_21_vs_full <- adjustedRandIndex(comparison$full, comparison$items_21)
ari_10_vs_full <- adjustedRandIndex(comparison$full, comparison$items_10)
ari_10_vs_21 <- adjustedRandIndex(comparison$items_21, comparison$items_10)

cat("\n\n=== Cluster Agreement Analysis ===\n")
cat("\nAdjusted Rand Index (1.0 = perfect agreement, 0 = random):\n")
cat("  21 items vs Full:", round(ari_21_vs_full, 3), "\n")
cat("  10 items vs Full:", round(ari_10_vs_full, 3), "\n")
cat("  10 items vs 21 items:", round(ari_10_vs_21, 3), "\n")

# Summary table
summary_table <- data.frame(
  Item_Set = c("Full (42)", "Top 21", "Top 10"),
  N_IAS = c(21, length(top_21_ias), length(top_10_ias)),
  N_IATS = c(21, length(top_21_iats), length(top_10_iats)),
  Silhouette = c(result_full$silhouette, result_21$silhouette, result_10$silhouette),
  ARI_vs_Full = c(1.0, ari_21_vs_full, ari_10_vs_full)
)

cat("\n\n=== SUMMARY: Cluster Preservation ===\n")
print(summary_table)

cat("\n\nInterpretation:\n")
cat("- ARI > 0.80: Excellent agreement (clusters preserved)\n")
cat("- ARI 0.65-0.80: Good agreement (similar structure)\n")
cat("- ARI 0.50-0.65: Moderate agreement (some differences)\n")
cat("- ARI < 0.50: Poor agreement (different structure)\n")

# =============================================================================
# 8. CLUSTER FLOW ANALYSIS (RECLASSIFICATION)
# =============================================================================

cat("\n\n=============================================================================\n")
cat("CLUSTER FLOW ANALYSIS: WHO GETS RECLASSIFIED?\n")
cat("=============================================================================\n\n")

# Function to run clustering and return profile labels
run_and_label_clusters <- function(ias_col, iats_col, data) {
  ias_z <- scale(data[[ias_col]])
  iats_z <- scale(data[[iats_col]])
  cluster_data <- data.frame(ias = ias_z, iats = iats_z)

  set.seed(123)
  km <- kmeans(cluster_data, centers = 3, nstart = 25)

  centers <- data.frame(
    cluster = 1:3,
    ias_z = km$centers[, 1],
    iats_z = km$centers[, 2]
  )

  centers <- centers %>%
    mutate(
      profile = case_when(
        ias_z > 0.3 & iats_z < -0.3 ~ "Efficient",
        iats_z > 0.5 ~ "Hypervigilant",
        ias_z < -0.5 ~ "Uncertain",
        TRUE ~ "Mixed"
      )
    )

  cluster_to_profile <- setNames(centers$profile, centers$cluster)
  cluster_to_profile[as.character(km$cluster)]
}

# Get profile assignments for each item set
dfc_complete <- dfc[complete.cases(dfc[, c("ias_full", "iats_full")]), ]

dfc_complete$profile_full <- run_and_label_clusters("ias_full", "iats_full", dfc_complete)
dfc_complete$profile_21 <- run_and_label_clusters("ias_21", "iats_21", dfc_complete)
dfc_complete$profile_10 <- run_and_label_clusters("ias_10", "iats_10", dfc_complete)

# Transition matrices
cat("=== Transition Matrices ===\n")

cat("\nFull 42 -> Top 21:\n")
t1 <- table(dfc_complete$profile_full, dfc_complete$profile_21)
print(t1)
cat("% staying same profile:", round(100 * sum(diag(t1)) / sum(t1), 1), "%\n")

cat("\nTop 21 -> Top 10:\n")
t2 <- table(dfc_complete$profile_21, dfc_complete$profile_10)
print(t2)
cat("% staying same profile:", round(100 * sum(diag(t2)) / sum(t2), 1), "%\n")

cat("\nFull 42 -> Top 10:\n")
t3 <- table(dfc_complete$profile_full, dfc_complete$profile_10)
print(t3)
cat("% staying same profile:", round(100 * sum(diag(t3)) / sum(t3), 1), "%\n")

# Analyze who gets reclassified
dfc_complete$ias_z_full <- scale(dfc_complete$ias_full)[,1]
dfc_complete$iats_z_full <- scale(dfc_complete$iats_full)[,1]
dfc_complete$changed_21 <- dfc_complete$profile_full != dfc_complete$profile_21
dfc_complete$changed_10 <- dfc_complete$profile_full != dfc_complete$profile_10

cat("\n\n=== Who Gets Reclassified? ===\n")
cat("\nMean |IAS z-score| for those who changed vs. stayed (Full->21):\n")
cat("  Changed:", round(mean(abs(dfc_complete$ias_z_full[dfc_complete$changed_21])), 2), "\n")
cat("  Stayed:", round(mean(abs(dfc_complete$ias_z_full[!dfc_complete$changed_21])), 2), "\n")

cat("\nMean |IATS z-score| for those who changed vs. stayed (Full->21):\n")
cat("  Changed:", round(mean(abs(dfc_complete$iats_z_full[dfc_complete$changed_21])), 2), "\n")
cat("  Stayed:", round(mean(abs(dfc_complete$iats_z_full[!dfc_complete$changed_21])), 2), "\n")

cat("\n-> People who get reclassified tend to be closer to cluster boundaries\n")
cat("   (smaller absolute z-scores = less extreme on IAS/IATS)\n")

sink()

# =============================================================================
# SAVE CSV OUTPUTS
# =============================================================================

write.csv(consensus_df, file.path(supplementary_data_dir, "item_consensus_ranking.csv"), row.names = FALSE)
write.csv(saturation_results, file.path(supplementary_data_dir, "saturation_results.csv"), row.names = FALSE)
write.csv(config_results, file.path(supplementary_data_dir, "minimal_set_comparison.csv"), row.names = FALSE)
write.csv(lasso_df, file.path(supplementary_data_dir, "lasso_coefficients.csv"), row.names = FALSE)
write.csv(centrality_df, file.path(supplementary_data_dir, "network_centrality.csv"), row.names = FALSE)
write.csv(summary_table, file.path(supplementary_data_dir, "cluster_preservation_summary.csv"), row.names = FALSE)

# =============================================================================
# CREATE FIGURES
# =============================================================================

# Figure: Transition heatmaps
trans_full_21 <- as.data.frame(t1) %>%
  rename(From = Var1, To = Var2, Count = Freq) %>%
  group_by(From) %>%
  mutate(Pct = Count / sum(Count) * 100) %>%
  ungroup()

trans_full_10 <- as.data.frame(t3) %>%
  rename(From = Var1, To = Var2, Count = Freq) %>%
  group_by(From) %>%
  mutate(Pct = Count / sum(Count) * 100) %>%
  ungroup()

p_heat1 <- ggplot(trans_full_21, aes(x = To, y = From, fill = Pct)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = paste0(round(Pct), "%\n(", Count, ")")),
            size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "#FFF3E0", high = "#E65100",
                       midpoint = 50, limits = c(0, 100),
                       name = "% of\noriginal") +
  labs(title = "Full 42 -> Top 21 Items", x = "Profile with Top 21", y = "Profile with Full 42") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10)) +
  coord_fixed()

p_heat2 <- ggplot(trans_full_10, aes(x = To, y = From, fill = Pct)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = paste0(round(Pct), "%\n(", Count, ")")),
            size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "#FFF3E0", high = "#E65100",
                       midpoint = 50, limits = c(0, 100),
                       name = "% of\noriginal") +
  labs(title = "Full 42 -> Top 10 Items", x = "Profile with Top 10", y = "Profile with Full 42") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10)) +
  coord_fixed()

p_heatmaps <- p_heat1 + p_heat2 +
  plot_annotation(
    title = "Profile Reclassification When Using Reduced Item Sets",
    subtitle = "Diagonal shows % remaining in same profile; off-diagonal shows reclassification",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40"))
  )

ggsave(file.path(supplementary_plots_dir, "2_s_4_transition_heatmap.png"),
       p_heatmaps, width = 12, height = 5.5, dpi = 300)

# Figure: Cluster preservation
profile_colors <- c("Efficient" = "#4DAF4A", "Hypervigilant" = "#FF7F00",
                    "Uncertain" = "#E41A1C", "Mixed" = "#3498DB")

plot_clusters <- function(result, ias_col, iats_col, data, title) {
  plot_data <- data.frame(
    ias = scale(data[[ias_col]])[,1],
    iats = scale(data[[iats_col]])[,1],
    cluster = result$cluster_assignments
  )
  plot_data <- plot_data[complete.cases(plot_data), ]

  cluster_to_label <- setNames(result$profiles$label, result$profiles$cluster)
  plot_data$profile <- factor(cluster_to_label[as.character(plot_data$cluster)],
                               levels = c("Efficient", "Hypervigilant", "Uncertain", "Mixed"))

  center_df <- data.frame(
    x = result$profiles$ias_z,
    y = result$profiles$iats_z,
    profile = factor(result$profiles$label,
                     levels = c("Efficient", "Hypervigilant", "Uncertain", "Mixed")),
    label = result$profiles$label
  )

  ggplot(plot_data, aes(x = ias, y = iats, color = profile)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_point(data = center_df, aes(x = x, y = y, color = profile), size = 5, shape = 18) +
    geom_label(data = center_df, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, vjust = -1, size = 2.5, fontface = "bold",
               fill = "white", alpha = 0.8, label.size = 0.3) +
    scale_color_manual(values = profile_colors, drop = FALSE) +
    labs(title = title, subtitle = paste0("Silhouette = ", round(result$silhouette, 3)),
         x = "IAS (z-score)", y = "IATS (z-score)") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "gray40")) +
    coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))
}

p1 <- plot_clusters(result_full, "ias_full", "iats_full", dfc, "Full 42 Items")
p2 <- plot_clusters(result_21, "ias_21", "iats_21", dfc, "Top 21 Items")
p3 <- plot_clusters(result_10, "ias_10", "iats_10", dfc, "Top 10 Items")

combined <- p1 + p2 + p3 +
  plot_annotation(
    title = "Cluster Preservation: Do Reduced Item Sets Recover the Same Profiles?",
    subtitle = "Comparing 3-cluster solutions across full and reduced item sets",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray40"))
  )

ggsave(file.path(supplementary_plots_dir, "2_s_3_cluster_preservation.png"),
       combined, width = 14, height = 5, dpi = 300)

# Figure: Who gets reclassified scatter
# (Plot output disabled to keep the plot set limited to approved figures.)
p_scatter <- ggplot(dfc_complete, aes(x = ias_z_full, y = iats_z_full)) +
  geom_point(aes(color = changed_10), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("FALSE" = "#95A5A6", "TRUE" = "#E74C3C"),
                     labels = c("Same profile", "Reclassified"),
                     name = "Full 42 -> Top 10") +
  labs(title = "Who Gets Reclassified with Reduced Items?",
       subtitle = "Red points = assigned different profile with 10 items vs. 42 items",
       x = "IAS (z-score, full scale)", y = "IATS (z-score, full scale)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        legend.position = "bottom")

cat("\n\nItem selection analysis complete.\n")
cat("Outputs: analysis_output/, supplementary_data/, plots/\n")
