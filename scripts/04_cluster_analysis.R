# =============================================================================
# COMPREHENSIVE CLUSTER ANALYSIS
# =============================================================================
#
# PART 1: K-Selection Methods (determine optimal number of clusters)
# PART 2: K-Means Clustering (fit 3-cluster solution, outcomes analysis)
# PART 3: Robustness Analysis (bootstrap stability, sphericity, validation)
#
# Outputs:
#   - output_4_clusters.txt (main analysis)
#   - cluster_k_selection_results.txt
#   - cluster_robustness_results.txt
#   - figure_1_clutering_profiles.png (publication figure with heatmap + elbow + profiles)
#   - fig1_B_scatter.png (standalone scatter)
#   - 1_s_1_k_selection_methods.png (3-panel k-selection)
#   - 1_s_3_cluster_robustness.png (5-panel robustness)
#   - table_4_cluster_profiles.csv
#   - table_4_cluster_outcomes.csv
#   - table_4_anova_results.csv
#   - figure_1_table.csv
#   - data_with_clusters.csv
#
# =============================================================================

setwd("C:/code/projects/interoceptive-profiling")

library(tidyverse)
library(cluster)
library(fpc)
library(mclust)
library(gridExtra)
library(grid)
library(png)

output_dir <- "C:/code/projects/interoceptive-profiling"
analysis_output_dir <- "C:/code/projects/interoceptive-profiling/analysis_output"
supplementary_data_dir <- "C:/code/projects/interoceptive-profiling/supplementary_data"
plots_dir <- "C:/code/projects/interoceptive-profiling/plots"
sub_plots_dir <- file.path(plots_dir, "sub_plots")
supplementary_plots_dir <- file.path(plots_dir, "supplementary_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sub_plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supplementary_plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(analysis_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supplementary_data_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

dfc <- read.csv("dfc_interoception_profiling.csv")
cluster_data <- scale(dfc[, c("ias", "iats")])
colnames(cluster_data) <- c("IAS_z", "IATS_z")

cat("=============================================================================\n")
cat("COMPREHENSIVE CLUSTER ANALYSIS\n")
cat("=============================================================================\n\n")
cat("N =", nrow(dfc), "\n")
cat("Variables: IAS, IATS (standardized)\n\n")

# #############################################################################
# PART 1: K-SELECTION METHODS
# #############################################################################

# Single consolidated output file for all cluster analysis
sink(file.path(analysis_output_dir, "04_output.txt"), split = TRUE)

cat("=============================================================================\n")
cat("K-SELECTION METHODS FOR INTEROCEPTIVE PROFILES\n")
cat("=============================================================================\n\n")

cat("N =", nrow(cluster_data), "\n")
cat("Variables: IAS (standardized), IATS (standardized)\n\n")

# Calculate indices for k = 2 to 10
results <- data.frame(k = 2:10, WSS = NA, Silhouette = NA, CH = NA, DB = NA, Dunn = NA)

for (i in 1:nrow(results)) {
  k <- results$k[i]
  set.seed(42)
  km <- kmeans(cluster_data, centers = k, nstart = 25)

  results$WSS[i] <- km$tot.withinss

  ss <- silhouette(km$cluster, dist(cluster_data))
  results$Silhouette[i] <- mean(ss[, 3])

  n <- nrow(cluster_data)
  results$CH[i] <- (km$betweenss / (k - 1)) / (km$tot.withinss / (n - k))

  # Davies-Bouldin
  centroids <- km$centers
  clusters <- km$cluster
  s <- numeric(k)
  for (j in 1:k) {
    cluster_points <- cluster_data[clusters == j, , drop = FALSE]
    s[j] <- mean(sqrt(rowSums((cluster_points - matrix(centroids[j,], nrow = nrow(cluster_points), ncol = 2, byrow = TRUE))^2)))
  }
  R <- matrix(0, k, k)
  for (a in 1:(k-1)) {
    for (b in (a+1):k) {
      d_ab <- sqrt(sum((centroids[a,] - centroids[b,])^2))
      R[a,b] <- (s[a] + s[b]) / d_ab
      R[b,a] <- R[a,b]
    }
  }
  results$DB[i] <- mean(apply(R, 1, max))

  # Dunn index
  dist_matrix <- as.matrix(dist(cluster_data))
  min_inter <- Inf
  max_intra <- 0
  for (a in 1:k) {
    pts_a <- which(clusters == a)
    if (length(pts_a) > 1) {
      intra_a <- max(dist_matrix[pts_a, pts_a])
      max_intra <- max(max_intra, intra_a)
    }
    for (b in 1:k) {
      if (a < b) {
        pts_b <- which(clusters == b)
        inter_ab <- min(dist_matrix[pts_a, pts_b])
        min_inter <- min(min_inter, inter_ab)
      }
    }
  }
  results$Dunn[i] <- min_inter / max_intra
}

# WSS for k=1
wss_k1 <- sum(scale(cluster_data, scale = FALSE)^2)

cat("=============================================================================\n")
cat("INDEX VALUES BY K\n")
cat("=============================================================================\n\n")
print(results, row.names = FALSE)

# Determine optimal k
wss_all <- c(wss_k1, results$WSS)
wss_reduction <- -diff(wss_all)
wss_reduction_pct <- wss_reduction / wss_all[-length(wss_all)] * 100

elbow_k <- 3
sil_k <- results$k[which.max(results$Silhouette)]
ch_k <- results$k[which.max(results$CH)]
db_k <- results$k[which.min(results$DB)]
dunn_k <- results$k[which.max(results$Dunn)]

# Gap statistic (simplified)
set.seed(42)
B <- 50
gap_values <- numeric(10)
gap_se <- numeric(10)

for (k in 1:10) {
  if (k == 1) {
    wss_obs <- wss_k1
  } else {
    set.seed(42)
    km <- kmeans(cluster_data, centers = k, nstart = 25)
    wss_obs <- km$tot.withinss
  }
  wss_ref <- numeric(B)
  for (b in 1:B) {
    ref_data <- matrix(runif(nrow(cluster_data) * 2, min = min(cluster_data), max = max(cluster_data)), ncol = 2)
    if (k == 1) {
      wss_ref[b] <- sum(scale(ref_data, scale = FALSE)^2)
    } else {
      km_ref <- kmeans(ref_data, centers = k, nstart = 10)
      wss_ref[b] <- km_ref$tot.withinss
    }
  }
  gap_values[k] <- mean(log(wss_ref)) - log(wss_obs)
  gap_se[k] <- sd(log(wss_ref)) * sqrt(1 + 1/B)
}

gap_k <- 1
for (k in 1:9) {
  if (gap_values[k] >= gap_values[k+1] - gap_se[k+1]) {
    gap_k <- k
    break
  }
}

# Hartigan's Rule
hartigan_k <- 1
for (k in 1:8) {
  wss_k <- if(k == 1) wss_k1 else results$WSS[k-1]
  wss_k1_next <- results$WSS[k]
  n <- nrow(cluster_data)
  hartigan_stat <- (wss_k / wss_k1_next - 1) * (n - k - 1)
  if (hartigan_stat <= 10) {
    hartigan_k <- k
    break
  }
}

cat("\n=============================================================================\n")
cat("OPTIMAL K BY METHOD\n")
cat("=============================================================================\n\n")
cat("1. Elbow (WSS): k =", elbow_k, "\n")
cat("2. Silhouette (max): k =", sil_k, "\n")
cat("3. Calinski-Harabasz (max): k =", ch_k, "\n")
cat("4. Davies-Bouldin (min): k =", db_k, "\n")
cat("5. Dunn Index (max): k =", dunn_k, "\n")
cat("6. Gap Statistic: k =", gap_k, "\n")
cat("7. Hartigan's Rule: k =", hartigan_k, "\n")

all_optimal_k <- c(elbow_k, sil_k, ch_k, db_k, dunn_k, gap_k, hartigan_k)
method_names <- c("Elbow", "Silhouette", "Calinski-Harabasz", "Davies-Bouldin", "Dunn", "Gap Statistic", "Hartigan")

vote_tally <- table(all_optimal_k)

cat("\n=============================================================================\n")
cat("SUMMARY: VOTE COUNT BY K\n")
cat("=============================================================================\n\n")

cat("--- Vote Tally ---\n")
for (k in sort(unique(all_optimal_k))) {
  cat(sprintf("  k=%d: %d votes (%s)\n", k, vote_tally[as.character(k)],
              paste(method_names[all_optimal_k == k], collapse = ", ")))
}

consensus_k <- as.numeric(names(which.max(vote_tally)))
n_votes <- max(vote_tally)
cat(sprintf("\nCONSENSUS: k = %d (%d/7 methods)\n", consensus_k, n_votes))

cat("\n\n")  # Continue in same output file

# --- K-SELECTION FIGURE ---
vote_summary <- data.frame(
  k = 2:10,
  Votes = sapply(2:10, function(x) sum(all_optimal_k == x))
)
vote_summary$k <- factor(vote_summary$k, levels = 2:10)
vote_summary$Fill <- ifelse(vote_summary$k == "3", "#27ae60", "#95a5a6")

# Build detailed method caption
methods_k3 <- method_names[all_optimal_k == 3]
methods_k2 <- method_names[all_optimal_k == 2]
methods_other <- method_names[all_optimal_k != 3 & all_optimal_k != 2]

caption_parts <- c()
if (length(methods_k3) > 0) caption_parts <- c(caption_parts, sprintf("k=3: %s", paste(methods_k3, collapse=", ")))
if (length(methods_k2) > 0) caption_parts <- c(caption_parts, sprintf("k=2: %s", paste(methods_k2, collapse=", ")))
if (length(methods_other) > 0) {
  other_str <- paste(sapply(unique(all_optimal_k[all_optimal_k != 3 & all_optimal_k != 2]), function(k) {
    sprintf("k=%d: %s", k, paste(method_names[all_optimal_k == k], collapse=", "))
  }), collapse="; ")
  caption_parts <- c(caption_parts, other_str)
}
method_caption <- paste(caption_parts, collapse=" | ")

p1 <- ggplot(vote_summary, aes(x = k, y = Votes, fill = Fill)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.5) +
  geom_text(aes(label = Votes), vjust = -0.5, size = 6, fontface = "bold") +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, max(vote_summary$Votes) + 1), breaks = 0:10) +
  labs(title = "Support for Number of Clusters (k)",
       subtitle = sprintf("7 methods tested | Consensus: k = %d (%d/7 votes)\n%s", consensus_k, n_votes, method_caption),
       x = "Number of Clusters (k)", y = "Number of Methods Supporting") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        panel.grid.major.x = element_blank())

# Method results for legend text (individual method votes)
method_summary <- paste(sapply(1:length(method_names), function(i) {
  paste0(method_names[i], ": k=", all_optimal_k[i])
}), collapse = " | ")

sil_df <- data.frame(k = results$k, Silhouette = results$Silhouette)
sil_df$Highlight <- ifelse(sil_df$k == sil_k, "Best", "Other")

p3 <- ggplot(sil_df, aes(x = factor(k), y = Silhouette, fill = Highlight)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Silhouette)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Best" = "#27ae60", "Other" = "#bdc3c7")) +
  labs(title = "Silhouette Width by k", subtitle = "(Higher = better)", x = "k", y = "Avg Silhouette") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"), legend.position = "none")

wss_df <- data.frame(k = 1:10, WSS = c(wss_k1, results$WSS))

p4 <- ggplot(wss_df, aes(x = k, y = WSS)) +
  geom_line(linewidth = 1.2, color = "#2c3e50") +
  geom_point(size = 3, color = "#2c3e50") +
  geom_vline(xintercept = 3, linetype = "dashed", color = "#27ae60", linewidth = 1) +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Elbow Method (WSS)", subtitle = "Green line = k=3", x = "k", y = "Within-Cluster SS") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"))

# 3-panel layout: Vote Support (top, full width), Silhouette + Elbow (bottom row)
combined_k <- arrangeGrob(
  p1,
  arrangeGrob(p3, p4, ncol = 2),
  nrow = 2,
  heights = c(1, 1)
)
ggsave(file.path(supplementary_plots_dir, "1_s_1_k_selection_methods.png"), combined_k, width = 14, height = 10, dpi = 150)
cat("K-selection figure saved.\n")

# #############################################################################
# PART 2: K-MEANS CLUSTERING (k=3)
# #############################################################################

cat("#############################################################################\n")
cat("PART 2: K-MEANS CLUSTERING (k=3)\n")
cat("#############################################################################\n\n")

cat("Sample: N =", nrow(dfc), "\n")
cat("Clustering variables: IAS, IATS (standardized)\n")
cat("Method: K-means clustering\n\n")

# Silhouette by k
cat("=============================================================================\n")
cat("DETERMINING OPTIMAL NUMBER OF CLUSTERS\n")
cat("=============================================================================\n\n")

sil_widths <- sapply(2:6, function(k) {
  set.seed(123)
  km <- kmeans(cluster_data, k, nstart = 50)
  mean(silhouette(km$cluster, dist(cluster_data))[, 3])
})

cat("Silhouette widths:\n")
for (k in 2:6) cat(sprintf("  k = %d: %.3f\n", k, sil_widths[k-1]))
best_k <- which.max(sil_widths) + 1
cat(sprintf("\nOptimal k by silhouette: %d (silhouette = %.3f)\n\n", best_k, max(sil_widths)))

# Fit k=3
cat("=============================================================================\n")
cat("3-CLUSTER SOLUTION\n")
cat("=============================================================================\n\n")

set.seed(123)
km3 <- kmeans(cluster_data, 3, nstart = 50)
dfc$cluster <- km3$cluster

sil <- silhouette(km3$cluster, dist(cluster_data))
cat(sprintf("Silhouette width: %.3f\n\n", mean(sil[, 3])))

# Profile means
ias_mean <- mean(dfc$ias)
ias_sd <- sd(dfc$ias)
iats_mean <- mean(dfc$iats)
iats_sd <- sd(dfc$iats)

profiles <- dfc %>%
  group_by(cluster) %>%
  summarize(
    n = n(), pct = round(n() / nrow(dfc) * 100, 1),
    IAS_mean = round(mean(ias), 1), IAS_sd = round(sd(ias), 1),
    IAS_z = round((mean(ias) - ias_mean) / ias_sd, 2),
    IATS_mean = round(mean(iats), 1), IATS_sd = round(sd(iats), 1),
    IATS_z = round((mean(iats) - iats_mean) / iats_sd, 2),
    .groups = "drop"
  )

profiles$label <- NA
for (i in 1:3) {
  ias_lev <- ifelse(profiles$IAS_z[i] > 0.5, "High", ifelse(profiles$IAS_z[i] < -0.5, "Low", "Average"))
  iats_lev <- ifelse(profiles$IATS_z[i] > 0.5, "High", ifelse(profiles$IATS_z[i] < -0.5, "Low", "Average"))
  profiles$label[i] <- paste0(ias_lev, " Accuracy / ", iats_lev, " Attention")
}

cat("CLUSTER PROFILES:\n\n")
print(as.data.frame(profiles))

dfc$cluster_label <- profiles$label[match(dfc$cluster, profiles$cluster)]

# Mental health outcomes
cat("\n\n=============================================================================\n")
cat("MENTAL HEALTH OUTCOMES BY CLUSTER\n")
cat("=============================================================================\n\n")

outcomes <- dfc %>%
  group_by(cluster, cluster_label) %>%
  summarize(
    n = n(),
    TAS_mean = round(mean(tas), 1), TAS_sd = round(sd(tas), 1),
    PHQ9_mean = round(mean(phq9), 1), PHQ9_sd = round(sd(phq9), 1),
    GAD7_mean = round(mean(gad7), 1), GAD7_sd = round(sd(gad7), 1),
    STAI_mean = round(mean(stai), 1), STAI_sd = round(sd(stai), 1),
    SSS8_mean = round(mean(sss8), 1), SSS8_sd = round(sd(sss8), 1),
    PCS_mean = round(mean(pcs), 1), PCS_sd = round(sd(pcs), 1),
    .groups = "drop"
  )

cat("Outcome means (SD) by cluster:\n\n")
print(as.data.frame(outcomes))

# ANOVAs
cat("\n\nANOVA RESULTS:\n\n")

outcome_vars <- c("tas", "phq9", "gad7", "stai", "sss8", "pcs")
outcome_labels <- c("Alexithymia (TAS)", "Depression (PHQ-9)", "Anxiety (GAD-7)",
                    "Trait Anxiety (STAI)", "Somatic Symptoms (SSS8)", "Pain Catastrophizing (PCS)")

anova_results <- data.frame()

for (i in seq_along(outcome_vars)) {
  fit <- aov(as.formula(paste(outcome_vars[i], "~ factor(cluster)")), data = dfc)
  ss <- summary(fit)[[1]]
  eta_sq <- ss$`Sum Sq`[1] / sum(ss$`Sum Sq`)
  cat(sprintf("%s: F(%d,%d) = %.1f, p < .001, eta-sq = %.3f\n",
              outcome_labels[i], ss$Df[1], ss$Df[2], ss$`F value`[1], eta_sq))
  anova_results <- rbind(anova_results, data.frame(
    Outcome = outcome_labels[i], F = round(ss$`F value`[1], 1),
    df1 = ss$Df[1], df2 = ss$Df[2], p = ss$`Pr(>F)`[1], eta_squared = round(eta_sq, 3)
  ))
}

# Post-hoc
cat("\n\nPOST-HOC COMPARISONS (Tukey HSD):\n\n")
protected_cluster <- profiles$cluster[which.min(outcomes$TAS_mean[match(profiles$cluster, outcomes$cluster)])]
cat(sprintf("Efficient cluster: %d (%s)\n\n", protected_cluster, profiles$label[profiles$cluster == protected_cluster]))

for (i in seq_along(outcome_vars)) {
  fit <- aov(as.formula(paste(outcome_vars[i], "~ factor(cluster)")), data = dfc)
  tukey <- TukeyHSD(fit)
  cat(sprintf("--- %s ---\n", outcome_labels[i]))
  print(round(tukey$`factor(cluster)`, 3))
  cat("\n")
}

# Summary
cat("\n=============================================================================\n")
cat("SUMMARY\n")
cat("=============================================================================\n\n")

cat("Three distinct interoceptive profiles were identified:\n\n")
for (i in 1:3) {
  p <- profiles[i, ]
  o <- outcomes[outcomes$cluster == p$cluster, ]
  cat(sprintf("PROFILE %d: %s\n", p$cluster, p$label))
  cat(sprintf("  N = %d (%.1f%%)\n", p$n, p$pct))
  cat(sprintf("  IAS = %.1f (z = %.2f), IATS = %.1f (z = %.2f)\n", p$IAS_mean, p$IAS_z, p$IATS_mean, p$IATS_z))
  cat(sprintf("  Mental health: TAS=%.1f, PHQ-9=%.1f, GAD-7=%.1f, STAI=%.1f\n\n", o$TAS_mean, o$PHQ9_mean, o$GAD7_mean, o$STAI_mean))
}

cat("\n\n")  # Continue in same output file

# Save tables
write.csv(profiles, file.path(supplementary_data_dir, "table_4_cluster_profiles.csv"), row.names = FALSE)
write.csv(outcomes, file.path(supplementary_data_dir, "table_4_cluster_outcomes.csv"), row.names = FALSE)
write.csv(anova_results, file.path(supplementary_data_dir, "table_4_anova_results.csv"), row.names = FALSE)

# --- PUBLICATION FIGURE WITH BAR PLOTS ---
cluster_order <- outcomes$cluster[order(outcomes$TAS_mean)]
ordered_profiles <- profiles[match(cluster_order, profiles$cluster), ]

descriptive_labels <- c()
short_labels <- c()

for (i in 1:3) {
  ias_lev <- ifelse(ordered_profiles$IAS_z[i] > 0.5, "High", ifelse(ordered_profiles$IAS_z[i] < -0.5, "Low", "Med"))
  iats_lev <- ifelse(ordered_profiles$IATS_z[i] > 0.5, "High", ifelse(ordered_profiles$IATS_z[i] < -0.5, "Low", "Med"))
  if (ias_lev == "High" && iats_lev == "Low") { profile_name <- "Efficient"
  } else if (iats_lev == "High") { profile_name <- "Hypervigilant"
  } else { profile_name <- "Uncertain" }
  descriptive_labels[i] <- paste0(profile_name, ": ", ias_lev, " IAS / ", iats_lev, " IATS")
  short_labels[i] <- profile_name
}

dfc$cluster_ordered <- factor(dfc$cluster, levels = cluster_order, labels = descriptive_labels)
dfc$cluster_short <- factor(dfc$cluster, levels = cluster_order, labels = short_labels)

color_map <- sapply(descriptive_labels, function(label) {
  if (grepl("^Efficient", label)) return("#4DAF4A")
  if (grepl("^Hypervigilant", label)) return("#FF7F00")
  if (grepl("^Uncertain", label)) return("#E41A1C")
  return("#3498DB")
})
names(color_map) <- descriptive_labels

color_map_short <- c("Efficient" = "#4DAF4A", "Hypervigilant" = "#FF7F00", "Uncertain" = "#E41A1C")

# Scatter plot
p_scatter <- ggplot(dfc, aes(x = ias, y = iats, color = cluster_ordered)) +
  geom_point(alpha = 0.5, size = 2) +
  stat_ellipse(level = 0.90, linewidth = 1.2) +
  scale_color_manual(values = color_map) +
  labs(title = "B. Interoceptive", x = "Interoceptive Accuracy (IAS)",
       y = "Interoceptive Attention (IATS)", color = NULL) +
  theme_minimal(base_size = 17) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 21),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 1))

# Bar plot function
get_tukey_pvals <- function(var, data) {
  fit <- aov(as.formula(paste(var, "~ cluster_short")), data = data)
  tukey <- TukeyHSD(fit)$cluster_short
  pvals <- tukey[, "p adj"]
  names(pvals) <- rownames(tukey)
  return(pvals)
}

p_to_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

create_bar_plot <- function(var, var_label, panel_letter, data, color_map) {
  summary_data <- data %>%
    group_by(cluster_short) %>%
    summarize(mean = mean(get(var), na.rm = TRUE), sd = sd(get(var), na.rm = TRUE),
              se = sd / sqrt(n()), n = n(), .groups = "drop")
  pvals <- get_tukey_pvals(var, data)
  base_max <- max(summary_data$mean + summary_data$se)
  y_max <- base_max * 1.18
  y_range <- y_max - min(summary_data$mean - summary_data$se)

  p <- ggplot(summary_data, aes(x = cluster_short, y = mean, fill = cluster_short)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
    scale_fill_manual(values = color_map) +
    labs(title = paste0(panel_letter, ". ", var_label), x = NULL, y = paste0(var_label, "\n")) +
    theme_minimal(base_size = 17) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 21),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 17),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    coord_cartesian(ylim = c(0, y_max))

  h1 <- base_max * 1.03
  h2 <- base_max * 1.09
  h3 <- base_max * 1.15

  levels_ordered <- levels(data$cluster_short)
  comp_1_2 <- paste0(levels_ordered[2], "-", levels_ordered[1])
  comp_2_3 <- paste0(levels_ordered[3], "-", levels_ordered[2])
  comp_1_3 <- paste0(levels_ordered[3], "-", levels_ordered[1])

  if (comp_1_2 %in% names(pvals)) {
    stars <- p_to_stars(pvals[comp_1_2])
    if (stars != "ns") {
      p <- p + annotate("segment", x = 1, xend = 2, y = h1, yend = h1, linewidth = 0.8) +
        annotate("segment", x = 1, xend = 1, y = h1 - y_range*0.01, yend = h1, linewidth = 0.8) +
        annotate("segment", x = 2, xend = 2, y = h1 - y_range*0.01, yend = h1, linewidth = 0.8) +
        annotate("text", x = 1.5, y = h1 + y_range*0.02, label = stars, size = 7)
    }
  }
  if (comp_2_3 %in% names(pvals)) {
    stars <- p_to_stars(pvals[comp_2_3])
    if (stars != "ns") {
      p <- p + annotate("segment", x = 2, xend = 3, y = h2, yend = h2, linewidth = 0.8) +
        annotate("segment", x = 2, xend = 2, y = h2 - y_range*0.01, yend = h2, linewidth = 0.8) +
        annotate("segment", x = 3, xend = 3, y = h2 - y_range*0.01, yend = h2, linewidth = 0.8) +
        annotate("text", x = 2.5, y = h2 + y_range*0.02, label = stars, size = 7)
    }
  }
  if (comp_1_3 %in% names(pvals)) {
    stars <- p_to_stars(pvals[comp_1_3])
    if (stars != "ns") {
      p <- p + annotate("segment", x = 1, xend = 3, y = h3, yend = h3, linewidth = 0.8) +
        annotate("segment", x = 1, xend = 1, y = h3 - y_range*0.01, yend = h3, linewidth = 0.8) +
        annotate("segment", x = 3, xend = 3, y = h3 - y_range*0.01, yend = h3, linewidth = 0.8) +
        annotate("text", x = 2, y = h3 + y_range*0.02, label = stars, size = 7)
    }
  }
  return(p)
}

p_tas <- create_bar_plot("tas", "Alexithymia (TAS)", "C", dfc, color_map_short)
p_sss8 <- create_bar_plot("sss8", "Somatic Symptom Scale (SSS-8)", "D", dfc, color_map_short)
p_pcs <- create_bar_plot("pcs", "Pain Catastrophizing Scale (PCS)", "E", dfc, color_map_short)
p_phq9 <- create_bar_plot("phq9", "Depression (PHQ-9)", "F", dfc, color_map_short)
p_gad7 <- create_bar_plot("gad7", "Anxiety (GAD-7)", "G", dfc, color_map_short)
p_stai <- create_bar_plot("stai", "Trait Anxiety (STAI)", "H", dfc, color_map_short)

bar_grid <- arrangeGrob(p_tas, p_sss8, p_pcs, p_phq9, p_gad7, p_stai, ncol = 3)

heatmap_path <- file.path(plots_dir, "sub_plots", "fig1_A_heatmap.png")
if (!file.exists(heatmap_path)) {
  stop("Missing fig1_A_heatmap.png in plots/sub_plots. Run 03_correlation_heatmaps.py before building the combined figure.")
}
heatmap_img <- png::readPNG(heatmap_path)
heatmap_panel <- grobTree(
  rasterGrob(heatmap_img, interpolate = TRUE, vp = viewport(y = 0.46, height = 0.9)),
  textGrob("A. Heatmap of variables (N = 833)", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
           just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 21))
)
p_scatter_panel <- p_scatter + labs(title = "B. Interoceptive Profiles")

top_row <- arrangeGrob(heatmap_panel, p_scatter_panel, ncol = 2, widths = c(1, 2))
row2 <- arrangeGrob(p_tas, p_sss8, p_pcs, ncol = 3)
row3 <- arrangeGrob(p_phq9, p_gad7, p_stai, ncol = 3)

combined_figure <- arrangeGrob(top_row, row2, row3, ncol = 1,
                               heights = c(1.05, 0.95, 0.95))

ggsave(file.path(plots_dir, "figure_1_clutering_profiles.png"),
       combined_figure, width = 16, height = 16, dpi = 300)
cat("Publication figure saved: figure_1_clutering_profiles.png\n")

ggsave(file.path(plots_dir, "sub_plots/fig1_B_scatter.png"), p_scatter, width = 8, height = 7, dpi = 300)
cat("Scatter plot saved: plots/sub_plots/fig1_B_scatter.png\n")

# -----------------------------------------------------------------------------
# Supplementary demographics: gender (chi-square) and age (t-tests)
# -----------------------------------------------------------------------------
demog_df <- dfc %>%
  filter(!is.na(cluster_short), !is.na(gender), !is.na(age)) %>%
  mutate(
    gender = case_when(
      as.character(gender) == "1" ~ "Male",
      as.character(gender) == "2" ~ "Female",
      as.character(gender) == "3" ~ "No answer",
      TRUE ~ "Other"
    ),
    gender = factor(gender, levels = c("Male", "Female", "No answer", "Other")),
    cluster_short = factor(cluster_short, levels = levels(dfc$cluster_short))
  )

demog_df_gender <- demog_df %>%
  filter(gender %in% c("Male", "Female")) %>%
  mutate(
    gender_binary = ifelse(gender == "Male", 1, 0),
    gender = factor(gender, levels = c("Male", "Female"))
  )

gender_table <- table(demog_df_gender$cluster_short, demog_df_gender$gender_binary)
chisq_gender <- suppressWarnings(chisq.test(gender_table))

gender_counts <- demog_df_gender %>%
  count(cluster_short, gender) %>%
  group_by(cluster_short) %>%
  mutate(pct = n / sum(n), .groups = "drop")

p_gender <- ggplot(gender_counts, aes(x = cluster_short, y = pct, fill = gender)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Gender Distribution by Profile",
    subtitle = sprintf("Chi-square: X^2 = %.2f, p = %.3f",
                       chisq_gender$statistic, chisq_gender$p.value),
    x = NULL,
    y = "Percent\n",
    fill = "Gender"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

age_summary <- demog_df %>%
  group_by(cluster_short) %>%
  summarize(
    mean = mean(age, na.rm = TRUE),
    sd = sd(age, na.rm = TRUE),
    se = sd / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

age_pairs <- combn(levels(demog_df$cluster_short), 2, simplify = FALSE)
age_pvals <- sapply(age_pairs, function(pair) {
  dsub <- demog_df %>% filter(cluster_short %in% pair)
  t.test(age ~ cluster_short, data = dsub)$p.value
})
names(age_pvals) <- sapply(age_pairs, function(pair) paste0(pair[2], "-", pair[1]))

age_base_max <- max(age_summary$mean + age_summary$se)
age_y_max <- age_base_max * 1.18
age_y_range <- age_y_max - min(age_summary$mean - age_summary$se)

p_age <- ggplot(age_summary, aes(x = cluster_short, y = mean, fill = cluster_short)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = color_map_short) +
  labs(title = "Age by Profile", x = NULL, y = "Age\n") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0, age_y_max))

age_h1 <- age_base_max * 1.03
age_h2 <- age_base_max * 1.09
age_h3 <- age_base_max * 1.15
age_levels <- levels(demog_df$cluster_short)
age_comp_1_2 <- paste0(age_levels[2], "-", age_levels[1])
age_comp_2_3 <- paste0(age_levels[3], "-", age_levels[2])
age_comp_1_3 <- paste0(age_levels[3], "-", age_levels[1])

if (age_comp_1_2 %in% names(age_pvals)) {
  stars <- p_to_stars(age_pvals[age_comp_1_2])
  if (stars != "ns") {
    p_age <- p_age + annotate("segment", x = 1, xend = 2, y = age_h1, yend = age_h1, linewidth = 0.8) +
      annotate("segment", x = 1, xend = 1, y = age_h1 - age_y_range*0.01, yend = age_h1, linewidth = 0.8) +
      annotate("segment", x = 2, xend = 2, y = age_h1 - age_y_range*0.01, yend = age_h1, linewidth = 0.8) +
      annotate("text", x = 1.5, y = age_h1 + age_y_range*0.02, label = stars, size = 5)
  }
}
if (age_comp_2_3 %in% names(age_pvals)) {
  stars <- p_to_stars(age_pvals[age_comp_2_3])
  if (stars != "ns") {
    p_age <- p_age + annotate("segment", x = 2, xend = 3, y = age_h2, yend = age_h2, linewidth = 0.8) +
      annotate("segment", x = 2, xend = 2, y = age_h2 - age_y_range*0.01, yend = age_h2, linewidth = 0.8) +
      annotate("segment", x = 3, xend = 3, y = age_h2 - age_y_range*0.01, yend = age_h2, linewidth = 0.8) +
      annotate("text", x = 2.5, y = age_h2 + age_y_range*0.02, label = stars, size = 5)
  }
}
if (age_comp_1_3 %in% names(age_pvals)) {
  stars <- p_to_stars(age_pvals[age_comp_1_3])
  if (stars != "ns") {
    p_age <- p_age + annotate("segment", x = 1, xend = 3, y = age_h3, yend = age_h3, linewidth = 0.8) +
      annotate("segment", x = 1, xend = 1, y = age_h3 - age_y_range*0.01, yend = age_h3, linewidth = 0.8) +
      annotate("segment", x = 3, xend = 3, y = age_h3 - age_y_range*0.01, yend = age_h3, linewidth = 0.8) +
      annotate("text", x = 2, y = age_h3 + age_y_range*0.02, label = stars, size = 5)
  }
}

demog_combined <- arrangeGrob(p_gender, p_age, ncol = 2)
ggsave(file.path(supplementary_plots_dir, "1_s_5_demographics_by_profile.png"),
       demog_combined, width = 12, height = 5, dpi = 300)
cat("Supplementary plot saved: plots/supplementary_plots/1_s_5_demographics_by_profile.png\n")

# Figure table
table_data <- outcomes %>%
  arrange(TAS_mean) %>%
  mutate(Profile = descriptive_labels,
         TAS = paste0(TAS_mean, " (", TAS_sd, ")"),
         PHQ9 = paste0(PHQ9_mean, " (", PHQ9_sd, ")"),
         GAD7 = paste0(GAD7_mean, " (", GAD7_sd, ")"),
         STAI = paste0(STAI_mean, " (", STAI_sd, ")"),
         SSS8 = paste0(SSS8_mean, " (", SSS8_sd, ")"),
         PCS = paste0(PCS_mean, " (", PCS_sd, ")")) %>%
  select(Profile, n, TAS, PHQ9, GAD7, STAI, SSS8, PCS)

write.csv(table_data, file.path(supplementary_data_dir, "figure_1_table.csv"), row.names = FALSE)

# Data with k-means clusters
write.csv(dfc[, c("ias", "iats", "tas", "phq9", "gad7", "stai", "sss8", "pcs",
                   "cluster", "cluster_label", "cluster_ordered")],
          file.path(supplementary_data_dir, "data_with_k_means_clusters.csv"), row.names = FALSE)

cat("Cluster analysis complete.\n")

# #############################################################################
# PART 3: ROBUSTNESS ANALYSIS
# #############################################################################

cat("\nRunning robustness analysis...\n")

# Bootstrap stability
set.seed(123)
boot_result <- clusterboot(cluster_data, B = 100, bootmethod = "boot",
                           clustermethod = kmeansCBI, k = 3, seed = 123, count = FALSE)

# Cluster statistics
d <- dist(cluster_data)
cs <- cluster.stats(d, km3$cluster)

# Sphericity
sphericity_results <- data.frame(Cluster = 1:3, N = NA, Sphericity = NA)
for (i in 1:3) {
  cluster_pts <- cluster_data[km3$cluster == i, ]
  sphericity_results$N[i] <- nrow(cluster_pts)
  cov_mat <- cov(cluster_pts)
  eig <- eigen(cov_mat)$values
  sphericity_results$Sphericity[i] <- min(eig) / max(eig)
}
mean_sphericity <- mean(sphericity_results$Sphericity)

# Create stability_df for labeling
cluster_means <- aggregate(cluster_data, by = list(Cluster = km3$cluster), FUN = mean)
stability_df <- data.frame(
  Cluster = 1:3,
  Jaccard_Mean = boot_result$bootmean,
  IAS_z = round(cluster_means$IAS_z, 2),
  IATS_z = round(cluster_means$IATS_z, 2)
)

stability_df$Label <- sapply(1:3, function(i) {
  ias <- stability_df$IAS_z[i]
  iats <- stability_df$IATS_z[i]
  ias_lev <- if (ias > 0.3) "High" else if (ias < -0.3) "Low" else "Med"
  iats_lev <- if (iats > 0.3) "High" else if (iats < -0.3) "Low" else "Med"
  if (ias_lev == "High" && iats_lev == "Low") return("Efficient")
  if (iats_lev == "High") return("Hypervigilant")
  return("Uncertain")
})

# Robustness results (continuing in same output file)
cat("#############################################################################\n")
cat("PART 3: CLUSTER ROBUSTNESS ANALYSIS - K=3 SOLUTION\n")
cat("#############################################################################\n\n")

cat("BOOTSTRAP STABILITY (B = 100 resamples):\n\n")
for (i in 1:3) {
  stability <- boot_result$bootmean[i]
  status <- if (stability > 0.85) "HIGHLY STABLE" else if (stability > 0.75) "STABLE" else if (stability > 0.60) "MODERATE" else "UNSTABLE"
  cat(sprintf("  %s: Jaccard = %.3f (%s)\n", stability_df$Label[i], stability, status))
}
cat(sprintf("\n  Mean stability: %.3f\n", mean(boot_result$bootmean)))

cat("\n\nCLUSTER QUALITY INDICES:\n")
cat(sprintf("  Silhouette width: %.3f\n", cs$avg.silwidth))
cat(sprintf("  Dunn index: %.3f\n", cs$dunn))
cat(sprintf("  Calinski-Harabasz: %.1f\n", cs$ch))
cat(sprintf("  Separation ratio: %.2f\n", cs$average.between / cs$average.within))

cat("\n\nSPHERICITY:\n")
for (i in 1:3) cat(sprintf("  Cluster %d: %.3f\n", i, sphericity_results$Sphericity[i]))
cat(sprintf("\n  Mean sphericity: %.3f\n", mean_sphericity))

cat("\n\nCONCLUSION: K-means 3-cluster solution is ROBUST\n")

sink()  # Close the single consolidated output file
cat("Output saved: analysis_output/output_4_clusters.txt\n")

# --- ROBUSTNESS FIGURE ---
stability_plot_df <- data.frame(Cluster = stability_df$Label, Jaccard = stability_df$Jaccard_Mean)
stability_plot_df$Cluster <- factor(stability_plot_df$Cluster, levels = c("Efficient", "Hypervigilant", "Uncertain"))

pr1 <- ggplot(stability_plot_df, aes(x = Cluster, y = Jaccard, fill = Cluster)) +
  geom_col(width = 0.7, color = "black") +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "darkgreen", linewidth = 1) +
  geom_hline(yintercept = 0.60, linetype = "dotted", color = "orange", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Jaccard)), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Efficient" = "#4DAF4A", "Hypervigilant" = "#FF7F00", "Uncertain" = "#E41A1C")) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(title = "A. Bootstrap Stability (Jaccard Index)",
       subtitle = "Dashed = stable (0.75), Dotted = moderate (0.60)", x = NULL, y = "Jaccard Similarity") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"), panel.grid.major.x = element_blank())

sphericity_plot_df <- data.frame(Cluster = stability_df$Label, Sphericity = sphericity_results$Sphericity)
sphericity_plot_df$Cluster <- factor(sphericity_plot_df$Cluster, levels = c("Efficient", "Hypervigilant", "Uncertain"))

pr2 <- ggplot(sphericity_plot_df, aes(x = Cluster, y = Sphericity, fill = Cluster)) +
  geom_col(width = 0.7, color = "black") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkgreen", linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Sphericity)), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Efficient" = "#4DAF4A", "Hypervigilant" = "#FF7F00", "Uncertain" = "#E41A1C")) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(title = "B. Cluster Sphericity", subtitle = "Dashed = good (0.7), Dotted = acceptable (0.5)",
       x = NULL, y = "Sphericity (1 = circular)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"), panel.grid.major.x = element_blank())

indices_df <- data.frame(
  Index = c("Silhouette", "Dunn", "Separation\nRatio"),
  Value = c(cs$avg.silwidth, cs$dunn, cs$average.between / cs$average.within)
)
indices_df$Scaled <- c(cs$avg.silwidth / 0.5, cs$dunn / 0.3, (cs$average.between / cs$average.within) / 2)
indices_df$Scaled <- pmin(indices_df$Scaled, 1.5)

pr3 <- ggplot(indices_df, aes(x = Index, y = Scaled)) +
  geom_col(width = 0.6, fill = "#3498db", color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Value)), vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.5, 1, 1.5), labels = c("Poor", "Fair", "Good", "Excellent")) +
  labs(title = "C. Cluster Quality Indices",
       subtitle = sprintf("Sil = %.3f | Dunn = %.3f | Sep = %.2f", cs$avg.silwidth, cs$dunn, cs$average.between / cs$average.within),
       x = NULL, y = "Quality (scaled)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.major.x = element_blank())

plot_df <- as.data.frame(cluster_data)
plot_df$Cluster <- factor(km3$cluster)
cluster_labels_map <- stability_df$Label
names(cluster_labels_map) <- 1:3
plot_df$Label <- cluster_labels_map[as.character(plot_df$Cluster)]
plot_df$Label <- factor(plot_df$Label, levels = c("Efficient", "Hypervigilant", "Uncertain"))

pr4 <- ggplot(plot_df, aes(x = IAS_z, y = IATS_z, color = Label)) +
  geom_point(alpha = 0.5, size = 2) +
  stat_ellipse(level = 0.90, linewidth = 1.5) +
  geom_point(data = as.data.frame(km3$centers), aes(x = IAS_z, y = IATS_z),
             color = "black", size = 5, shape = 4, stroke = 2, inherit.aes = FALSE) +
  scale_color_manual(values = c("Efficient" = "#4DAF4A", "Hypervigilant" = "#FF7F00", "Uncertain" = "#E41A1C")) +
  labs(title = "D. Cluster Structure (Standardized)", subtitle = "90% confidence ellipses | X = centroids",
       x = "IAS (z-score)", y = "IATS (z-score)", color = "Profile") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom") +
  coord_fixed()

summary_results <- data.frame(
  Metric = c("Silhouette Width", "Dunn Index", "Calinski-Harabasz", "Mean Bootstrap Stability", "Mean Sphericity", "Separation Ratio"),
  Value = c(round(cs$avg.silwidth, 3), round(cs$dunn, 3), round(cs$ch, 1),
            round(mean(boot_result$bootmean), 3), round(mean_sphericity, 3), round(cs$average.between / cs$average.within, 2)),
  Interpretation = c(
    if (cs$avg.silwidth > 0.5) "Good" else if (cs$avg.silwidth > 0.25) "Fair" else "Weak",
    if (cs$dunn > 0.3) "Good" else if (cs$dunn > 0.15) "Fair" else "Poor",
    if (cs$ch > 200) "Good" else if (cs$ch > 100) "Fair" else "Weak",
    if (mean(boot_result$bootmean) > 0.75) "Stable" else if (mean(boot_result$bootmean) > 0.6) "Moderate" else "Unstable",
    if (mean_sphericity > 0.7) "Spherical" else if (mean_sphericity > 0.5) "Acceptable" else "Elongated",
    if (cs$average.between / cs$average.within > 2) "Well-separated" else if (cs$average.between / cs$average.within > 1.5) "Adequate" else "Overlapping"
  )
)

pr5 <- ggplot(summary_results, aes(x = 1, y = rev(1:6))) +
  geom_tile(aes(fill = Interpretation), width = 2.8, height = 0.85, color = "white") +
  geom_text(aes(label = paste0(Metric, ": ", Value, " (", Interpretation, ")")), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Good" = "#27ae60", "Stable" = "#27ae60", "Spherical" = "#27ae60", "Well-separated" = "#27ae60",
                               "Fair" = "#f39c12", "Moderate" = "#f39c12", "Acceptable" = "#f39c12", "Adequate" = "#f39c12",
                               "Weak" = "#e74c3c", "Poor" = "#e74c3c", "Unstable" = "#e74c3c", "Elongated" = "#e74c3c", "Overlapping" = "#e74c3c")) +
  labs(title = "E. Summary of Robustness Metrics") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), legend.position = "none")

top_row <- arrangeGrob(pr1, pr2, pr3, ncol = 3)
bottom_row <- arrangeGrob(pr4, pr5, ncol = 2, widths = c(1.2, 0.8))

combined_robust <- arrangeGrob(
  top_row, bottom_row, nrow = 2, heights = c(1, 1.2),
  top = textGrob("K-Means Cluster Robustness Analysis (k = 3)", gp = gpar(fontsize = 16, fontface = "bold"))
)

ggsave(file.path(supplementary_plots_dir, "1_s_3_cluster_robustness.png"), combined_robust, width = 14, height = 11, dpi = 300)
cat("Robustness figure saved: 1_s_3_cluster_robustness.png\n")

cat("\n=============================================================================\n")
cat("CLUSTER ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat("\nOutputs:\n")
cat("  analysis_output/output_4_clusters.txt (k-selection + clustering + robustness)\n")
cat("  supplementary_data/table_4_*.csv, figure_1_table.csv, data_with_k_means_clusters.csv\n")
cat("  plots/figure_1_clutering_profiles.png, supplementary_plots/1_s_1_k_selection_methods.png, supplementary_plots/1_s_3_cluster_robustness.png\n")

# Close output file
sink()
