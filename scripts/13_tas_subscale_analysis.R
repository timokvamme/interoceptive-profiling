# =============================================================================
# TAS-20 SUBSCALE DECOMPOSITION: Which Alexithymia Facet Drives the Effect?
# =============================================================================
#
# Motivated by reviewer suggestion (Luminet): "Would it be possible to run
# your analyses for each alexithymia facet separately, because the pattern
# of association usually varies as a function of the alexithymia facets?"
#
# TAS-20 SUBSCALES (Bagby, Parker, & Taylor, 1994):
#   DIF = Difficulty Identifying Feelings (items 1, 3, 6, 7, 9, 13, 14)
#   DDF = Difficulty Describing Feelings (items 2, 4, 11, 12, 17)
#   EOT = Externally Oriented Thinking  (items 5, 8, 10, 15, 16, 18, 19, 20)
#
# NOTE: Items 4, 5, 10, 18, 19 are already reverse-scored in this dataset
#       (verified: raw item sum == TAS total exactly)
#
# RESEARCH QUESTIONS:
#   1. Which subscale(s) mediate IAS -> Mental Health? (the TAS pathway)
#   2. Does the null IAS x IATS -> TAS interaction hold for all subscales?
#   3. Which subscale carries the dominant indirect effect?
#
# OUTPUTS:
#   - analysis_output/13_output.txt
#   - supplementary_data/tas_subscale_results.xlsx
#   - plots/supplementary_plots/3_s_8_tas_subscale_decomposition.png
#
# =============================================================================

setwd("C:/code/projects/intero_mod")

library(dplyr)
library(tidyr)
library(lavaan)
library(ggplot2)
library(gridExtra)
library(grid)
library(openxlsx)

# Output directories
analysis_output_dir <- "C:/code/projects/intero_mod/analysis_output"
supplementary_data_dir <- "C:/code/projects/intero_mod/supplementary_data"
plots_dir <- "C:/code/projects/intero_mod/plots"
sub_plots_dir <- file.path(plots_dir, "supplementary_plots")
dir.create(sub_plots_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD DATA AND COMPUTE SUBSCALES
# =============================================================================

dfc <- read.csv("C:/code/projects/mi/analyses/soma/results/dfc_interoception_profiling.csv")

cat("Loading data...\n")
cat("Sample: N =", nrow(dfc), "\n\n")

# TAS-20 subscale definitions (items already reverse-scored in dataset)
dif_items <- c(1, 3, 6, 7, 9, 13, 14)        # 7 items
ddf_items <- c(2, 4, 11, 12, 17)              # 5 items
eot_items <- c(5, 8, 10, 15, 16, 18, 19, 20)  # 8 items

dfc$dif <- rowSums(dfc[, paste0("tas_", dif_items)], na.rm = TRUE)
dfc$ddf <- rowSums(dfc[, paste0("tas_", ddf_items)], na.rm = TRUE)
dfc$eot <- rowSums(dfc[, paste0("tas_", eot_items)], na.rm = TRUE)

# Validate: subscale sum must equal TAS total
subscale_sum <- dfc$dif + dfc$ddf + dfc$eot
stopifnot(all(abs(dfc$tas - subscale_sum) < 0.01))

# Full-scale IAS/IATS (consistent with script 12)
dfc$ias_full <- rowMeans(dfc[, paste0("ias_", 1:21)], na.rm = TRUE)
dfc$iats_full <- rowMeans(dfc[, paste0("iats_", 1:21)], na.rm = TRUE)

# Outcome composites (consistent with script 12)
dfc$mh_composite <- rowMeans(cbind(scale(dfc$phq9), scale(dfc$gad7),
                                    scale(dfc$stai)), na.rm = TRUE)
dfc$somatic_combined <- rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE)

cat("Subscales computed. DIF + DDF + EOT = TAS total (verified).\n\n")

# =============================================================================
# BEGIN OUTPUT FILE
# =============================================================================

sink(file.path(analysis_output_dir, "13_output.txt"))

cat("=============================================================================\n")
cat("TAS-20 SUBSCALE DECOMPOSITION ANALYSIS\n")
cat("=============================================================================\n\n")
cat("Sample: N =", nrow(dfc), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
cat("Research question: Which TAS-20 subscale (DIF, DDF, EOT) drives the\n")
cat("alexithymia pathway in the dual-pathway model of interoception?\n\n")

# Initialize xlsx workbook
wb <- createWorkbook()

# =============================================================================
# 2. DESCRIPTIVES AND INTERNAL CONSISTENCY
# =============================================================================

cat("#############################################################################\n")
cat("PART 1: SUBSCALE DESCRIPTIVES AND INTERNAL CONSISTENCY\n")
cat("#############################################################################\n\n")

# Cronbach's alpha function
cronbach_alpha <- function(items_df) {
  items_df <- items_df[complete.cases(items_df), ]
  k <- ncol(items_df)
  item_vars <- apply(items_df, 2, var)
  total_var <- var(rowSums(items_df))
  alpha <- (k / (k - 1)) * (1 - sum(item_vars) / total_var)
  return(alpha)
}

subscale_names <- c("DIF", "DDF", "EOT", "TAS Total")
subscale_vars <- c("dif", "ddf", "eot", "tas")
n_items <- c(7, 5, 8, 20)
item_lists <- list(dif_items, ddf_items, eot_items, 1:20)

desc_table <- data.frame(
  Subscale = subscale_names,
  Items = n_items,
  Mean = NA, SD = NA, Min = NA, Max = NA, Alpha = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(subscale_vars)) {
  x <- dfc[[subscale_vars[i]]]
  desc_table$Mean[i] <- round(mean(x, na.rm = TRUE), 2)
  desc_table$SD[i] <- round(sd(x, na.rm = TRUE), 2)
  desc_table$Min[i] <- min(x, na.rm = TRUE)
  desc_table$Max[i] <- max(x, na.rm = TRUE)
  desc_table$Alpha[i] <- round(cronbach_alpha(dfc[, paste0("tas_", item_lists[[i]])]), 3)
}

cat("--- Subscale Descriptives ---\n\n")
cat(sprintf("%-12s %5s %8s %8s %5s %5s %8s\n",
            "Subscale", "Items", "Mean", "SD", "Min", "Max", "Alpha"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in 1:nrow(desc_table)) {
  cat(sprintf("%-12s %5d %8.2f %8.2f %5d %5d %8.3f\n",
              desc_table$Subscale[i], desc_table$Items[i],
              desc_table$Mean[i], desc_table$SD[i],
              desc_table$Min[i], desc_table$Max[i], desc_table$Alpha[i]))
}

# Inter-subscale correlations
cat("\n--- Inter-Subscale Correlations ---\n\n")
sub_cors <- cor(dfc[, c("dif", "ddf", "eot", "tas")], use = "complete.obs")
rownames(sub_cors) <- colnames(sub_cors) <- c("DIF", "DDF", "EOT", "TAS")
cat(sprintf("%-5s %8s %8s %8s %8s\n", "", "DIF", "DDF", "EOT", "TAS"))
for (i in 1:4) {
  cat(sprintf("%-5s %8.3f %8.3f %8.3f %8.3f\n",
              rownames(sub_cors)[i],
              sub_cors[i, 1], sub_cors[i, 2], sub_cors[i, 3], sub_cors[i, 4]))
}

addWorksheet(wb, "Descriptives")
writeData(wb, "Descriptives", desc_table)

# =============================================================================
# 3. CORRELATIONS WITH INTEROCEPTIVE AND CLINICAL MEASURES
# =============================================================================

cat("\n\n#############################################################################\n")
cat("PART 2: CORRELATIONS - SUBSCALES x INTEROCEPTIVE AND CLINICAL MEASURES\n")
cat("#############################################################################\n\n")

cor_vars <- c("ias_full", "iats_full", "sss8", "pcs", "phq9", "gad7", "stai",
              "mh_composite", "somatic_combined")
cor_labels <- c("IAS", "IATS", "SSS-8", "PCS", "PHQ-9", "GAD-7", "STAI",
                "MH Composite", "Somatic Composite")

cor_matrix <- matrix(NA, nrow = length(cor_vars), ncol = 4)
p_matrix <- matrix(NA, nrow = length(cor_vars), ncol = 4)
rownames(cor_matrix) <- rownames(p_matrix) <- cor_labels

for (j in 1:4) {
  subscale <- c("dif", "ddf", "eot", "tas")[j]
  for (i in seq_along(cor_vars)) {
    ct <- cor.test(dfc[[subscale]], dfc[[cor_vars[i]]], use = "complete.obs")
    cor_matrix[i, j] <- ct$estimate
    p_matrix[i, j] <- ct$p.value
  }
}
colnames(cor_matrix) <- colnames(p_matrix) <- c("DIF", "DDF", "EOT", "TAS")

cat("--- Pearson Correlations ---\n\n")
cat(sprintf("%-20s %8s %8s %8s %8s\n", "Variable", "DIF", "DDF", "EOT", "TAS"))
cat(paste(rep("-", 55), collapse = ""), "\n")
for (i in 1:nrow(cor_matrix)) {
  sig_markers <- sapply(1:4, function(j) {
    p <- p_matrix[i, j]
    if (p < 0.001) "***" else if (p < 0.01) "** " else if (p < 0.05) "*  " else "   "
  })
  cat(sprintf("%-20s %5.3f%s %5.3f%s %5.3f%s %5.3f%s\n",
              cor_labels[i],
              cor_matrix[i, 1], sig_markers[1],
              cor_matrix[i, 2], sig_markers[2],
              cor_matrix[i, 3], sig_markers[3],
              cor_matrix[i, 4], sig_markers[4]))
}
cat("\nNote: *** p < .001, ** p < .01, * p < .05\n")

# Identify strongest subscale correlate for each variable
cat("\n--- Strongest Subscale Correlate per Variable ---\n\n")
for (i in 1:nrow(cor_matrix)) {
  abs_cors <- abs(cor_matrix[i, 1:3])
  best <- which.max(abs_cors)
  cat(sprintf("  %-20s → %s (r = %.3f)\n",
              cor_labels[i], c("DIF", "DDF", "EOT")[best], cor_matrix[i, best]))
}

cor_df <- as.data.frame(cor_matrix)
cor_df$Variable <- cor_labels
addWorksheet(wb, "Correlations")
writeData(wb, "Correlations", cor_df)

# =============================================================================
# 4. MODERATION: IAS x IATS INTERACTION ON EACH SUBSCALE
# =============================================================================

cat("\n\n#############################################################################\n")
cat("PART 3: MODERATION - IAS x IATS INTERACTION ON EACH SUBSCALE\n")
cat("#############################################################################\n\n")

cat("Testing whether the null IAS x IATS -> TAS interaction holds for all\n")
cat("subscales, or whether one subscale shows a hidden interaction effect.\n\n")

# Standardize predictors
dfc$ias_z <- scale(dfc$ias_full)[, 1]
dfc$iats_z <- scale(dfc$iats_full)[, 1]
dfc$ias_x_iats <- dfc$ias_z * dfc$iats_z

mod_outcomes <- c("dif", "ddf", "eot", "tas")
mod_labels <- c("DIF", "DDF", "EOT", "TAS Total")

mod_results <- data.frame()

for (i in seq_along(mod_outcomes)) {
  dv <- mod_outcomes[i]

  fit_add <- lm(as.formula(paste(dv, "~ ias_z + iats_z")), data = dfc)
  fit_int <- lm(as.formula(paste(dv, "~ ias_z * iats_z")), data = dfc)

  int_coef <- coef(fit_int)["ias_z:iats_z"]
  int_se <- summary(fit_int)$coefficients["ias_z:iats_z", 2]
  int_p <- summary(fit_int)$coefficients["ias_z:iats_z", 4]
  delta_r2 <- summary(fit_int)$r.squared - summary(fit_add)$r.squared

  # Main effects from interaction model
  ias_b <- coef(fit_int)["ias_z"]
  ias_p <- summary(fit_int)$coefficients["ias_z", 4]
  iats_b <- coef(fit_int)["iats_z"]
  iats_p <- summary(fit_int)$coefficients["iats_z", 4]

  mod_results <- rbind(mod_results, data.frame(
    subscale = mod_labels[i],
    ias_beta = round(ias_b, 3),
    ias_p = round(ias_p, 4),
    iats_beta = round(iats_b, 3),
    iats_p = round(iats_p, 4),
    interaction_beta = round(int_coef, 3),
    interaction_se = round(int_se, 3),
    interaction_p = round(int_p, 4),
    delta_r2 = round(delta_r2, 4),
    significant = int_p < 0.05,
    stringsAsFactors = FALSE
  ))
}

cat(sprintf("%-12s %8s %8s %8s %8s %8s %8s  %s\n",
            "Subscale", "IAS b", "IAS p", "IATS b", "IATS p", "IxI b", "IxI p", "Sig"))
cat(paste(rep("-", 78), collapse = ""), "\n")
for (i in 1:nrow(mod_results)) {
  cat(sprintf("%-12s %8.3f %8.4f %8.3f %8.4f %8.3f %8.4f  %s\n",
              mod_results$subscale[i],
              mod_results$ias_beta[i], mod_results$ias_p[i],
              mod_results$iats_beta[i], mod_results$iats_p[i],
              mod_results$interaction_beta[i], mod_results$interaction_p[i],
              ifelse(mod_results$significant[i], "***", "")))
}

cat("\nDelta R-squared for interaction term:\n")
for (i in 1:nrow(mod_results)) {
  cat(sprintf("  %-12s: Delta R2 = %.4f\n", mod_results$subscale[i], mod_results$delta_r2[i]))
}

# Interpretation
any_sig <- any(mod_results$significant[1:3])
if (any_sig) {
  sig_subs <- mod_results$subscale[mod_results$significant[1:3]]
  cat(sprintf("\nKEY FINDING: IAS x IATS interaction IS significant for %s,\n",
              paste(sig_subs, collapse = " and ")))
  cat("suggesting the null effect on TAS total masks a subscale-specific interaction.\n")
} else {
  cat("\nKEY FINDING: IAS x IATS interaction is non-significant for ALL subscales,\n")
  cat("confirming that the null interaction on TAS total is not masking\n")
  cat("subscale-specific effects. The alexithymia pathway operates additively.\n")
}

addWorksheet(wb, "Moderation_Subscales")
writeData(wb, "Moderation_Subscales", mod_results)

# =============================================================================
# 5. SEPARATE DUAL-PATHWAY SEMS (ONE PER SUBSCALE)
# =============================================================================

cat("\n\n#############################################################################\n")
cat("PART 4: SEPARATE DUAL-PATHWAY SEMS (ONE PER SUBSCALE)\n")
cat("#############################################################################\n\n")

cat("Each model replaces TAS total with one subscale in the dual-pathway SEM.\n")
cat("This shows the total (non-unique) mediation through each subscale.\n\n")

run_subscale_sem <- function(data, subscale_var, subscale_label) {
  data$ias_z <- scale(data$ias_full)[, 1]
  data$iats_z <- scale(data$iats_full)[, 1]
  data$sub_z <- scale(data[[subscale_var]])[, 1]
  data$somatic_z <- scale(data$somatic_combined)[, 1]
  data$mh_z <- scale(data$mh_composite)[, 1]
  data$ias_x_iats <- data$ias_z * data$iats_z

  model <- '
    # Paths to mediators
    sub_z ~ a1*ias_z + a2*iats_z + a3*ias_x_iats
    somatic_z ~ b1*iats_z + b2*ias_z + b3*ias_x_iats

    # Mediators to outcome
    mh_z ~ c1*sub_z + c2*somatic_z + d1*ias_z + d2*iats_z + d3*ias_x_iats

    # Indirect effects
    ind_ias_sub := a1 * c1
    ind_iats_sub := a2 * c1
    ind_ias_som := b2 * c2
    ind_iats_som := b1 * c2
  '

  fit <- sem(model, data = data, se = "bootstrap", bootstrap = 1000,
             iseed = 42)
  params <- parameterEstimates(fit, ci = TRUE, standardized = TRUE)

  get_est <- function(lbl) {
    row <- params[params$label == lbl, ]
    if (nrow(row) > 0) {
      return(list(
        est = row$std.all[1], p = row$pvalue[1],
        ci_lo = row$ci.lower[1], ci_hi = row$ci.upper[1]
      ))
    }
    return(list(est = NA, p = NA, ci_lo = NA, ci_hi = NA))
  }

  list(
    label = subscale_label,
    fit = fit,
    params = params,
    # Path coefficients
    ias_to_sub = get_est("a1"),
    iats_to_sub = get_est("a2"),
    ixi_to_sub = get_est("a3"),
    ias_to_som = get_est("b2"),
    iats_to_som = get_est("b1"),
    ixi_to_som = get_est("b3"),
    sub_to_mh = get_est("c1"),
    som_to_mh = get_est("c2"),
    # Indirect effects
    ind_ias_sub = get_est("ind_ias_sub"),
    ind_iats_sub = get_est("ind_iats_sub"),
    ind_ias_som = get_est("ind_ias_som"),
    ind_iats_som = get_est("ind_iats_som")
  )
}

cat("Running bootstrap SEMs (1000 replicates each)...\n")
cat("This may take a few minutes.\n\n")

sem_dif <- run_subscale_sem(dfc, "dif", "DIF")
cat("  DIF model complete.\n")
sem_ddf <- run_subscale_sem(dfc, "ddf", "DDF")
cat("  DDF model complete.\n")
sem_eot <- run_subscale_sem(dfc, "eot", "EOT")
cat("  EOT model complete.\n")
sem_tas_total <- run_subscale_sem(dfc, "tas", "TAS Total")
cat("  TAS Total model complete.\n\n")

# Print results for each
sep_models <- list(sem_dif, sem_ddf, sem_eot)

for (m in sep_models) {
  cat("---", m$label, "as Mediator (replacing TAS total) ---\n\n")

  cat("  Paths to", m$label, "(alexithymia subscale):\n")
  cat(sprintf("    IAS → %-3s:       beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
              m$label, m$ias_to_sub$est, m$ias_to_sub$p,
              m$ias_to_sub$ci_lo, m$ias_to_sub$ci_hi))
  cat(sprintf("    IATS → %-3s:      beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
              m$label, m$iats_to_sub$est, m$iats_to_sub$p,
              m$iats_to_sub$ci_lo, m$iats_to_sub$ci_hi))
  cat(sprintf("    IASxIATS → %-3s:  beta = %6.3f (p = %.4f) [%.3f, %.3f] %s\n",
              m$label, m$ixi_to_sub$est, m$ixi_to_sub$p,
              m$ixi_to_sub$ci_lo, m$ixi_to_sub$ci_hi,
              ifelse(m$ixi_to_sub$p < 0.05, "***", "")))

  cat(sprintf("\n  %-3s → MH:          beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
              m$label, m$sub_to_mh$est, m$sub_to_mh$p,
              m$sub_to_mh$ci_lo, m$sub_to_mh$ci_hi))
  cat(sprintf("  Somatic → MH:      beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
              m$som_to_mh$est, m$som_to_mh$p,
              m$som_to_mh$ci_lo, m$som_to_mh$ci_hi))

  cat(sprintf("\n  Indirect: IAS → %-3s → MH = %.4f (p = %.4f) [%.4f, %.4f]\n",
              m$label, m$ind_ias_sub$est, m$ind_ias_sub$p,
              m$ind_ias_sub$ci_lo, m$ind_ias_sub$ci_hi))
  cat(sprintf("  Indirect: IAS → Som → MH  = %.4f (p = %.4f) [%.4f, %.4f]\n\n",
              m$ind_ias_som$est, m$ind_ias_som$p,
              m$ind_ias_som$ci_lo, m$ind_ias_som$ci_hi))
}

# =============================================================================
# 6. PARALLEL MEDIATION SEM (ALL SUBSCALES SIMULTANEOUSLY)
# =============================================================================

cat("\n#############################################################################\n")
cat("PART 5: PARALLEL MEDIATION SEM (ALL SUBSCALES + SOMATIC SIMULTANEOUSLY)\n")
cat("#############################################################################\n\n")

cat("This model includes DIF, DDF, EOT, and Somatic as simultaneous mediators.\n")
cat("It reveals the UNIQUE contribution of each subscale after controlling\n")
cat("for the others - the critical test of which subscale drives the effect.\n\n")

# Standardize all variables
dfc$ias_z <- scale(dfc$ias_full)[, 1]
dfc$iats_z <- scale(dfc$iats_full)[, 1]
dfc$dif_z <- scale(dfc$dif)[, 1]
dfc$ddf_z <- scale(dfc$ddf)[, 1]
dfc$eot_z <- scale(dfc$eot)[, 1]
dfc$somatic_z <- scale(dfc$somatic_combined)[, 1]
dfc$mh_z <- scale(dfc$mh_composite)[, 1]
dfc$ias_x_iats <- dfc$ias_z * dfc$iats_z

parallel_model <- '
  # Paths from predictors to DIF
  dif_z ~ a1_dif*ias_z + a2_dif*iats_z + a3_dif*ias_x_iats

  # Paths from predictors to DDF
  ddf_z ~ a1_ddf*ias_z + a2_ddf*iats_z + a3_ddf*ias_x_iats

  # Paths from predictors to EOT
  eot_z ~ a1_eot*ias_z + a2_eot*iats_z + a3_eot*ias_x_iats

  # Paths from predictors to Somatic
  somatic_z ~ b1*iats_z + b2*ias_z + b3*ias_x_iats

  # All mediators to outcome
  mh_z ~ c_dif*dif_z + c_ddf*ddf_z + c_eot*eot_z + c_som*somatic_z +
          d1*ias_z + d2*iats_z + d3*ias_x_iats

  # Allow subscale residuals to correlate
  dif_z ~~ ddf_z
  dif_z ~~ eot_z
  ddf_z ~~ eot_z

  # Indirect effects through each subscale (IAS pathway)
  ind_ias_dif := a1_dif * c_dif
  ind_ias_ddf := a1_ddf * c_ddf
  ind_ias_eot := a1_eot * c_eot
  ind_ias_som := b2 * c_som

  # Indirect effects through each subscale (IATS pathway)
  ind_iats_dif := a2_dif * c_dif
  ind_iats_ddf := a2_ddf * c_ddf
  ind_iats_eot := a2_eot * c_eot
  ind_iats_som := b1 * c_som

  # Total alexithymia indirect (sum of subscale indirects)
  total_alex_ias := ind_ias_dif + ind_ias_ddf + ind_ias_eot
  total_alex_iats := ind_iats_dif + ind_iats_ddf + ind_iats_eot

  # Contrasts: pairwise differences between subscale indirect effects
  diff_dif_ddf := ind_ias_dif - ind_ias_ddf
  diff_dif_eot := ind_ias_dif - ind_ias_eot
  diff_ddf_eot := ind_ias_ddf - ind_ias_eot
'

cat("Running parallel mediation SEM with bootstrap (1000 replicates)...\n")
fit_parallel <- sem(parallel_model, data = dfc, se = "bootstrap",
                    bootstrap = 1000, iseed = 42)
cat("  Parallel model complete.\n\n")

params_par <- parameterEstimates(fit_parallel, ci = TRUE, standardized = TRUE)

get_par <- function(lbl) {
  row <- params_par[params_par$label == lbl, ]
  if (nrow(row) > 0) {
    return(list(est = row$std.all[1], p = row$pvalue[1],
                ci_lo = row$ci.lower[1], ci_hi = row$ci.upper[1]))
  }
  return(list(est = NA, p = NA, ci_lo = NA, ci_hi = NA))
}

# --- Print paths to each subscale ---

cat("=== A-Paths: Predictors → Alexithymia Subscales ===\n\n")
cat(sprintf("%-20s %8s %8s %8s %8s %8s\n",
            "Path", "beta", "p", "CI_lo", "CI_hi", "Sig"))
cat(paste(rep("-", 65), collapse = ""), "\n")

a_paths <- list(
  c("a1_dif", "IAS → DIF"), c("a2_dif", "IATS → DIF"), c("a3_dif", "IxI → DIF"),
  c("a1_ddf", "IAS → DDF"), c("a2_ddf", "IATS → DDF"), c("a3_ddf", "IxI → DDF"),
  c("a1_eot", "IAS → EOT"), c("a2_eot", "IATS → EOT"), c("a3_eot", "IxI → EOT")
)

for (ap in a_paths) {
  v <- get_par(ap[1])
  cat(sprintf("%-20s %8.3f %8.4f %8.3f %8.3f %8s\n",
              ap[2], v$est, v$p, v$ci_lo, v$ci_hi,
              ifelse(v$p < 0.05, "***", "")))
  if (grepl("IxI", ap[2])) cat("\n")
}

# --- Print b-paths: mediators → MH ---

cat("\n=== B-Paths: Mediators → Mental Health (unique effects) ===\n\n")
cat(sprintf("%-20s %8s %8s %8s %8s %8s\n",
            "Path", "beta", "p", "CI_lo", "CI_hi", "Sig"))
cat(paste(rep("-", 65), collapse = ""), "\n")

b_labels <- list(
  c("c_dif", "DIF → MH"), c("c_ddf", "DDF → MH"),
  c("c_eot", "EOT → MH"), c("c_som", "Somatic → MH")
)

for (bl in b_labels) {
  v <- get_par(bl[1])
  cat(sprintf("%-20s %8.3f %8.4f %8.3f %8.3f %8s\n",
              bl[2], v$est, v$p, v$ci_lo, v$ci_hi,
              ifelse(v$p < 0.05, "***", "")))
}

# --- Print indirect effects ---

cat("\n\n=== Indirect Effects Through Each Mediator (IAS → Mediator → MH) ===\n\n")
cat(sprintf("%-25s %10s %8s %10s %10s %8s\n",
            "Pathway", "Effect", "p", "CI_lo", "CI_hi", "Sig"))
cat(paste(rep("-", 75), collapse = ""), "\n")

ind_labels <- list(
  c("ind_ias_dif", "IAS → DIF → MH"),
  c("ind_ias_ddf", "IAS → DDF → MH"),
  c("ind_ias_eot", "IAS → EOT → MH"),
  c("ind_ias_som", "IAS → Somatic → MH"),
  c("total_alex_ias", "IAS → All Alex → MH")
)

for (il in ind_labels) {
  v <- get_par(il[1])
  cat(sprintf("%-25s %10.4f %8.4f %10.4f %10.4f %8s\n",
              il[2], v$est, v$p, v$ci_lo, v$ci_hi,
              ifelse(v$p < 0.05, "***", "")))
}

# --- Print contrasts ---

cat("\n\n=== Pairwise Contrasts: IAS Indirect Effect Differences ===\n\n")
cat(sprintf("%-25s %10s %8s %10s %10s %8s\n",
            "Contrast", "Diff", "p", "CI_lo", "CI_hi", "Sig"))
cat(paste(rep("-", 75), collapse = ""), "\n")

contrasts <- list(
  c("diff_dif_ddf", "DIF - DDF"),
  c("diff_dif_eot", "DIF - EOT"),
  c("diff_ddf_eot", "DDF - EOT")
)

for (ct in contrasts) {
  v <- get_par(ct[1])
  cat(sprintf("%-25s %10.4f %8.4f %10.4f %10.4f %8s\n",
              ct[2], v$est, v$p, v$ci_lo, v$ci_hi,
              ifelse(v$p < 0.05, "***", "")))
}

# --- Proportion of alexithymia pathway carried by each subscale ---

cat("\n\n=== Proportion of Alexithymia Pathway by Subscale ===\n\n")

ind_dif <- get_par("ind_ias_dif")$est
ind_ddf <- get_par("ind_ias_ddf")$est
ind_eot <- get_par("ind_ias_eot")$est
total_alex <- get_par("total_alex_ias")$est

# Use absolute values for proportions
abs_total <- abs(ind_dif) + abs(ind_ddf) + abs(ind_eot)
if (abs_total > 0) {
  pct_dif <- abs(ind_dif) / abs_total * 100
  pct_ddf <- abs(ind_ddf) / abs_total * 100
  pct_eot <- abs(ind_eot) / abs_total * 100
} else {
  pct_dif <- pct_ddf <- pct_eot <- NA
}

cat(sprintf("  DIF: %.4f (%.1f%% of alexithymia pathway)\n", ind_dif, pct_dif))
cat(sprintf("  DDF: %.4f (%.1f%% of alexithymia pathway)\n", ind_ddf, pct_ddf))
cat(sprintf("  EOT: %.4f (%.1f%% of alexithymia pathway)\n", ind_eot, pct_eot))
cat(sprintf("  Total alexithymia indirect: %.4f\n", total_alex))

# Add to xlsx
par_results <- data.frame(
  parameter = params_par$label,
  est_std = round(params_par$std.all, 4),
  p = round(params_par$pvalue, 4),
  ci_lower = round(params_par$ci.lower, 4),
  ci_upper = round(params_par$ci.upper, 4)
)
par_results <- par_results[par_results$parameter != "", ]

addWorksheet(wb, "Parallel_SEM")
writeData(wb, "Parallel_SEM", par_results)

# =============================================================================
# 7. SUMMARY AND INTERPRETATION
# =============================================================================

cat("\n\n#############################################################################\n")
cat("SUMMARY: KEY FINDINGS\n")
cat("#############################################################################\n\n")

# Determine which subscale has the largest indirect effect
ind_effects <- c(DIF = ind_dif, DDF = ind_ddf, EOT = ind_eot)
dominant <- names(which.max(abs(ind_effects)))
dominant_pct <- max(abs(pct_dif), abs(pct_ddf), abs(pct_eot), na.rm = TRUE)

# Check which subscales have significant indirect effects
sig_dif <- get_par("ind_ias_dif")$p < 0.05
sig_ddf <- get_par("ind_ias_ddf")$p < 0.05
sig_eot <- get_par("ind_ias_eot")$p < 0.05
sig_subs_indirect <- c("DIF", "DDF", "EOT")[c(sig_dif, sig_ddf, sig_eot)]

cat("1. SUBSCALE STRUCTURE:\n")
cat(sprintf("   - DIF-DDF correlation: r = %.3f (shared affective alexithymia)\n",
            sub_cors["DIF", "DDF"]))
cat(sprintf("   - DIF-EOT correlation: r = %.3f\n", sub_cors["DIF", "EOT"]))
cat(sprintf("   - DDF-EOT correlation: r = %.3f\n", sub_cors["DDF", "EOT"]))

cat("\n2. MODERATION (IAS x IATS → subscale):\n")
for (i in 1:3) {
  cat(sprintf("   - IAS x IATS → %s: b = %.3f (p = %.4f) %s\n",
              mod_results$subscale[i], mod_results$interaction_beta[i],
              mod_results$interaction_p[i],
              ifelse(mod_results$significant[i], "*** SIGNIFICANT", "n.s.")))
}

cat("\n3. UNIQUE MEDIATION (parallel model):\n")
if (length(sig_subs_indirect) > 0) {
  cat(sprintf("   - Significant mediators: %s\n", paste(sig_subs_indirect, collapse = ", ")))
} else {
  cat("   - No individual subscale shows a significant unique indirect effect\n")
}
cat(sprintf("   - Dominant subscale: %s (%.1f%% of alexithymia pathway)\n",
            dominant, dominant_pct))

cat("\n4. SUBSCALE-SPECIFIC PATH COEFFICIENTS (unique, from parallel model):\n")
for (sub in c("dif", "ddf", "eot")) {
  a1 <- get_par(paste0("a1_", sub))
  c1 <- get_par(paste0("c_", sub))
  sub_upper <- toupper(sub)
  cat(sprintf("   - IAS → %s: beta = %.3f (p = %.4f); %s → MH: beta = %.3f (p = %.4f)\n",
              sub_upper, a1$est, a1$p, sub_upper, c1$est, c1$p))
}

cat("\n5. CLINICAL IMPLICATION:\n")
cat(sprintf("   %s is the primary alexithymia facet linking interoceptive\n", dominant))
cat("   accuracy to mental health outcomes in the dual-pathway model.\n")

sink()

# =============================================================================
# 8. SUPPLEMENTARY FIGURE: TAS SUBSCALE DECOMPOSITION (with TAS Total)
# =============================================================================

cat("Creating supplementary figure...\n")

# --- Collect separate-model coefficients for all 4 mediators (DIF, DDF, EOT, TAS) ---

all_models <- list(
  list(m = sem_dif, name = "DIF"),
  list(m = sem_ddf, name = "DDF"),
  list(m = sem_eot, name = "EOT"),
  list(m = sem_tas_total, name = "TAS Total")
)

# --- Panel A: Path coefficients forest plot ---

coef_data <- data.frame(
  path = character(), mediator = character(), beta = numeric(),
  ci_lo = numeric(), ci_hi = numeric(), p = numeric(),
  stringsAsFactors = FALSE
)

for (ml in all_models) {
  nm <- ml$name; m <- ml$m
  # IAS → mediator
  coef_data <- rbind(coef_data, data.frame(
    path = paste0("IAS -> ", nm), mediator = nm,
    beta = m$ias_to_sub$est, ci_lo = m$ias_to_sub$ci_lo,
    ci_hi = m$ias_to_sub$ci_hi, p = m$ias_to_sub$p))
  # IATS → mediator
  coef_data <- rbind(coef_data, data.frame(
    path = paste0("IATS -> ", nm), mediator = nm,
    beta = m$iats_to_sub$est, ci_lo = m$iats_to_sub$ci_lo,
    ci_hi = m$iats_to_sub$ci_hi, p = m$iats_to_sub$p))
  # IxI → mediator
  coef_data <- rbind(coef_data, data.frame(
    path = paste0("IASxIATS -> ", nm), mediator = nm,
    beta = m$ixi_to_sub$est, ci_lo = m$ixi_to_sub$ci_lo,
    ci_hi = m$ixi_to_sub$ci_hi, p = m$ixi_to_sub$p))
  # mediator → MH
  coef_data <- rbind(coef_data, data.frame(
    path = paste0(nm, " -> MH"), mediator = nm,
    beta = m$sub_to_mh$est, ci_lo = m$sub_to_mh$ci_lo,
    ci_hi = m$sub_to_mh$ci_hi, p = m$sub_to_mh$p))
}

coef_data$sig <- ifelse(coef_data$p < 0.05, "p < .05", "n.s.")

# Order paths: group by path type, TAS Total first as reference
path_levels <- rev(c(
  "IAS -> TAS Total", "IAS -> DIF", "IAS -> DDF", "IAS -> EOT",
  "IATS -> TAS Total", "IATS -> DIF", "IATS -> DDF", "IATS -> EOT",
  "IASxIATS -> TAS Total", "IASxIATS -> DIF", "IASxIATS -> DDF", "IASxIATS -> EOT",
  "TAS Total -> MH", "DIF -> MH", "DDF -> MH", "EOT -> MH"
))
coef_data$path <- factor(coef_data$path, levels = path_levels)

coef_data$mediator <- factor(coef_data$mediator,
  levels = c("TAS Total", "DIF", "DDF", "EOT"))

med_colors <- c("TAS Total" = "#333333", "DIF" = "#E64B35",
                "DDF" = "#4DBBD5", "EOT" = "#00A087")

p_a <- ggplot(coef_data, aes(x = beta, y = path, color = mediator, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), size = 0.5, linewidth = 0.7) +
  scale_color_manual(values = med_colors, name = "Mediator") +
  scale_shape_manual(values = c("p < .05" = 16, "n.s." = 1), name = "") +
  labs(x = expression("Standardized " * beta * " [95% Bootstrap CI]"),
       y = NULL,
       title = "A. Path Coefficients: TAS Total vs. Subscales (separate models)") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# --- Panel B: Indirect effects comparison ---

ind_data <- data.frame(
  pathway = character(), mediator = character(),
  effect = numeric(), ci_lo = numeric(), ci_hi = numeric(), p = numeric(),
  stringsAsFactors = FALSE
)

for (ml in all_models) {
  nm <- ml$name; m <- ml$m
  ind_data <- rbind(ind_data, data.frame(
    pathway = paste0("IAS -> ", nm, " -> MH"), mediator = nm,
    effect = m$ind_ias_sub$est, ci_lo = m$ind_ias_sub$ci_lo,
    ci_hi = m$ind_ias_sub$ci_hi, p = m$ind_ias_sub$p))
}
# Add somatic pathway (same across models, take from TAS Total model)
ind_data <- rbind(ind_data, data.frame(
  pathway = "IAS -> Somatic -> MH", mediator = "Somatic",
  effect = sem_tas_total$ind_ias_som$est,
  ci_lo = sem_tas_total$ind_ias_som$ci_lo,
  ci_hi = sem_tas_total$ind_ias_som$ci_hi,
  p = sem_tas_total$ind_ias_som$p))

ind_data$sig <- ifelse(ind_data$p < 0.05, "p < .05", "n.s.")

ind_data$pathway <- factor(ind_data$pathway,
  levels = rev(c("IAS -> TAS Total -> MH", "IAS -> DIF -> MH",
                 "IAS -> DDF -> MH", "IAS -> EOT -> MH",
                 "IAS -> Somatic -> MH")))
ind_data$mediator <- factor(ind_data$mediator,
  levels = c("TAS Total", "DIF", "DDF", "EOT", "Somatic"))

ind_colors <- c(med_colors, "Somatic" = "#7E6148")

p_b <- ggplot(ind_data, aes(x = effect, y = pathway, color = mediator, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), size = 0.7, linewidth = 0.9) +
  scale_color_manual(values = ind_colors, name = "Mediator") +
  scale_shape_manual(values = c("p < .05" = 16, "n.s." = 1), name = "") +
  labs(x = expression("Indirect Effect [95% Bootstrap CI]"),
       y = NULL,
       title = "B. Indirect Effects: IAS -> Mediator -> MH") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# --- Panel C: Subscale proportion bar chart (with TAS Total reference) ---

prop_data <- data.frame(
  subscale = c("TAS Total", "DIF", "DDF", "EOT"),
  pct = c(100, pct_dif, pct_ddf, pct_eot),
  stringsAsFactors = FALSE
)
prop_data$subscale <- factor(prop_data$subscale,
  levels = c("TAS Total", "DIF", "DDF", "EOT"))

prop_colors <- c("TAS Total" = "#333333", "DIF" = "#E64B35",
                 "DDF" = "#4DBBD5", "EOT" = "#00A087")

p_c <- ggplot(prop_data, aes(x = subscale, y = pct, fill = subscale)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = prop_colors) +
  labs(x = NULL, y = "% of Alexithymia Pathway",
       title = "C. Subscale Decomposition of Alexithymia Mediation") +
  ylim(0, 115) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Combine panels
fig <- arrangeGrob(
  p_a, p_b, p_c,
  layout_matrix = rbind(c(1, 1), c(2, 3)),
  heights = c(1.2, 1)
)

fig_path <- file.path(sub_plots_dir, "3_s_8_tas_subscale_decomposition.png")
ggsave(fig_path, fig, width = 14, height = 10, dpi = 300, bg = "white")

cat(sprintf("Figure saved: %s\n", fig_path))

# Save xlsx
saveWorkbook(wb, file.path(supplementary_data_dir, "tas_subscale_results.xlsx"),
             overwrite = TRUE)

cat("Output saved: analysis_output/13_output.txt\n")
cat("Excel saved: supplementary_data/tas_subscale_results.xlsx\n")
cat("Done.\n")
