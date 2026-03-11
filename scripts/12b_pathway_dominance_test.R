# =============================================================================
# PATHWAY DOMINANCE ANALYSIS: Statistical Testing of TAS vs Somatic Pathways
# =============================================================================
#
# PURPOSE: Properly compute and statistically test pathway dominance by profile
#
# This analysis:
#   1. Decomposes ALL variance from IAS, IATS, and their interaction
#   2. Computes percentages for each pathway (TAS vs Somatic vs Direct)
#   3. Computes conditional indirect effects at profile-typical values
#   4. Bootstraps everything for proper confidence intervals
#   5. Formally tests whether pathways differ (overall and by profile)
#
# KEY OUTPUTS:
#   - Pathway percentages with 95% CIs
#   - INDEX of moderated mediation (formal test)
#   - Profile-specific conditional indirect effects
#   - Statistical tests for pathway dominance
#
# =============================================================================

setwd("C:/code/projects/intero_mod")

library(dplyr)
library(tidyr)
library(lavaan)
library(ggplot2)
library(gridExtra)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dfc <- read.csv("C:/code/projects/mi/analyses/soma/results/dfc_vviq_q_k.csv")

cat("=============================================================================\n")
cat("PATHWAY DOMINANCE ANALYSIS: STATISTICAL TESTING\n")
cat("=============================================================================\n\n")
cat("Sample: N =", nrow(dfc), "\n\n")

# Save full-scale totals for clustering (consistent with 04_cluster_analysis.R)
ias_full <- dfc$ias
iats_full <- dfc$iats

# Compute composite outcomes
dfc$mh <- rowMeans(cbind(scale(dfc$phq9), scale(dfc$gad7), scale(dfc$stai)), na.rm = TRUE)
dfc$somatic <- rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE)

# Standardize all variables (full-scale IAS/IATS for SEM, consistent with 12_output.txt)
dfc$ias_z <- scale(dfc$ias)[,1]
dfc$iats_z <- scale(dfc$iats)[,1]
dfc$tas_z <- scale(dfc$tas)[,1]
dfc$somatic_z <- scale(dfc$somatic)[,1]
dfc$mh_z <- scale(dfc$mh)[,1]

# Create interaction term
dfc$ias_x_iats <- dfc$ias_z * dfc$iats_z

# =============================================================================
# 2. GET PROFILE-TYPICAL VALUES FROM CLUSTER ANALYSIS
# =============================================================================

# Assign profiles via clustering on FULL-SCALE totals
# (matching 04_cluster_analysis.R: seed 123, nstart 50)
set.seed(123)
km <- kmeans(scale(cbind(ias_full, iats_full)), centers = 3, nstart = 50)
dfc$cluster <- km$cluster

# Identify which cluster is which based on full-scale IAS/IATS patterns
dfc$ias_full_z <- scale(ias_full)[,1]
dfc$iats_full_z <- scale(iats_full)[,1]
cluster_means <- aggregate(cbind(ias_full_z, iats_full_z) ~ cluster, data = dfc, mean)
cat("Cluster Centers:\n")
print(cluster_means)

# Assign labels based on pattern
# Hypervigilant: High IAS, High IATS
# Uncertain: Low IAS, variable IATS
# Efficient: Medium-High IAS, Low IATS
cluster_means$profile <- NA
cluster_means$profile[which.max(cluster_means$ias_full_z + cluster_means$iats_full_z)] <- "Hypervigilant"
cluster_means$profile[which.min(cluster_means$ias_full_z)] <- "Uncertain"
cluster_means$profile[is.na(cluster_means$profile)] <- "Efficient"

# Map to data
dfc$profile <- cluster_means$profile[match(dfc$cluster, cluster_means$cluster)]
dfc$ias_full_z <- NULL
dfc$iats_full_z <- NULL

# Profile-typical z-scores for conditional effects (using minimal-scale z-scores for SEM)
profile_values <- aggregate(cbind(ias_z, iats_z) ~ profile, data = dfc, mean)
names(profile_values) <- c("profile", "IAS_typical", "IATS_typical")
cat("\nProfile-Typical Values (z-scored):\n")
print(profile_values)

# =============================================================================
# 3. COMPREHENSIVE SEM MODEL WITH ALL INDIRECT EFFECTS
# =============================================================================

cat("\n\n=============================================================================\n")
cat("PART 1: FULL VARIANCE DECOMPOSITION SEM\n")
cat("=============================================================================\n\n")

# Get profile-typical values for conditional effects
hyper_iats <- profile_values$IATS_typical[profile_values$profile == "Hypervigilant"]
uncert_iats <- profile_values$IATS_typical[profile_values$profile == "Uncertain"]
effic_iats <- profile_values$IATS_typical[profile_values$profile == "Efficient"]

hyper_ias <- profile_values$IAS_typical[profile_values$profile == "Hypervigilant"]
uncert_ias <- profile_values$IAS_typical[profile_values$profile == "Uncertain"]
effic_ias <- profile_values$IAS_typical[profile_values$profile == "Efficient"]

# Build model string with profile-specific values inserted
model_full <- sprintf('
  # =================================================================
  # STRUCTURAL PATHS
  # =================================================================

  # Paths to TAS (Alexithymia)
  tas_z ~ a1*ias_z + a2*iats_z + a3*ias_x_iats

  # Paths to Somatic
  somatic_z ~ b1*ias_z + b2*iats_z + b3*ias_x_iats

  # Paths to Mental Health
  mh_z ~ c1*tas_z + c2*somatic_z + d1*ias_z + d2*iats_z + d3*ias_x_iats

  # =================================================================
  # SIMPLE INDIRECT EFFECTS (at mean levels, interaction = 0)
  # =================================================================

  # Through TAS pathway
  ind_ias_tas := a1 * c1
  ind_iats_tas := a2 * c1
  ind_int_tas := a3 * c1

  # Through Somatic pathway
  ind_ias_som := b1 * c2
  ind_iats_som := b2 * c2
  ind_int_som := b3 * c2

  # =================================================================
  # TOTAL EFFECTS BY PATHWAY (signed sums)
  # =================================================================

  total_tas_signed := ind_ias_tas + ind_iats_tas + ind_int_tas
  total_som_signed := ind_ias_som + ind_iats_som + ind_int_som
  total_direct := d1 + d2 + d3

  # =================================================================
  # ABSOLUTE VALUE SUMS FOR PERCENTAGE CALCULATION
  # Note: We compute abs() manually in post-processing
  # Here we define components for later calculation
  # =================================================================

  # =================================================================
  # PATHWAY DIFFERENCE TESTS (signed)
  # =================================================================

  # Test: Is Somatic pathway larger than TAS pathway?
  diff_som_minus_tas := total_som_signed - total_tas_signed

  # =================================================================
  # INDEX OF MODERATED MEDIATION
  # Tests whether the indirect effect differs based on moderator (IATS)
  # =================================================================

  # For Somatic pathway: Does IAS->Somatic->MH depend on IATS?
  index_mod_med_som := b3 * c2

  # For TAS pathway: Does IAS->TAS->MH depend on IATS?
  index_mod_med_tas := a3 * c1

  # =================================================================
  # CONDITIONAL INDIRECT EFFECTS BY PROFILE
  # IAS effect on MH, conditional on IATS level
  # =================================================================

  # --- HYPERVIGILANT (IATS_z = %s) ---
  # IAS -> TAS -> MH at Hypervigilant IATS level
  cond_ias_tas_hyper := (a1 + a3*%s) * c1
  # IAS -> Somatic -> MH at Hypervigilant IATS level
  cond_ias_som_hyper := (b1 + b3*%s) * c2

  # --- UNCERTAIN (IATS_z = %s) ---
  cond_ias_tas_uncert := (a1 + a3*%s) * c1
  cond_ias_som_uncert := (b1 + b3*%s) * c2

  # --- EFFICIENT (IATS_z = %s) ---
  cond_ias_tas_effic := (a1 + a3*%s) * c1
  cond_ias_som_effic := (b1 + b3*%s) * c2

  # =================================================================
  # CONDITIONAL PATHWAY COMPARISONS BY PROFILE
  # =================================================================

  # Difference (Somatic - TAS) for each profile
  diff_hyper := cond_ias_som_hyper - cond_ias_tas_hyper
  diff_uncert := cond_ias_som_uncert - cond_ias_tas_uncert
  diff_effic := cond_ias_som_effic - cond_ias_tas_effic

', round(hyper_iats, 4), round(hyper_iats, 4), round(hyper_iats, 4),
   round(uncert_iats, 4), round(uncert_iats, 4), round(uncert_iats, 4),
   round(effic_iats, 4), round(effic_iats, 4), round(effic_iats, 4))

# =============================================================================
# 4. FIT MODEL WITH BOOTSTRAP
# =============================================================================

cat("Fitting SEM with bootstrap (5000 iterations)...\n")
cat("This may take a few minutes.\n\n")

set.seed(42)
fit_boot <- sem(model_full, data = dfc, se = "bootstrap", bootstrap = 5000)

# Extract parameter estimates with bootstrap CIs
params <- parameterEstimates(fit_boot, ci = TRUE, standardized = TRUE)

# =============================================================================
# 5. EXTRACT AND DISPLAY RESULTS
# =============================================================================

cat("\n=============================================================================\n")
cat("RESULTS: PATH COEFFICIENTS\n")
cat("=============================================================================\n\n")

# Direct paths
direct_paths <- params[params$op == "~", c("lhs", "rhs", "est", "se", "ci.lower", "ci.upper", "pvalue")]
names(direct_paths) <- c("DV", "IV", "Est", "SE", "CI_low", "CI_high", "p")
cat("Direct Path Coefficients:\n")
print(direct_paths, row.names = FALSE)

cat("\n=============================================================================\n")
cat("RESULTS: INDIRECT EFFECTS\n")
cat("=============================================================================\n\n")

# Indirect effects
indirect_labels <- c("ind_ias_tas", "ind_iats_tas", "ind_int_tas",
                     "ind_ias_som", "ind_iats_som", "ind_int_som")
indirect <- params[params$label %in% indirect_labels,
                   c("label", "est", "se", "ci.lower", "ci.upper", "pvalue")]

cat("Individual Indirect Effects:\n")
cat(sprintf("%-15s %8s %8s %12s %8s\n", "Effect", "Est", "SE", "95% CI", "p"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in 1:nrow(indirect)) {
  cat(sprintf("%-15s %8.4f %8.4f [%6.4f, %6.4f] %8.4f %s\n",
              indirect$label[i], indirect$est[i], indirect$se[i],
              indirect$ci.lower[i], indirect$ci.upper[i], indirect$pvalue[i],
              ifelse(indirect$pvalue[i] < 0.05, "*", "")))
}

cat("\n=============================================================================\n")
cat("RESULTS: PATHWAY TOTALS AND COMPARISON\n")
cat("=============================================================================\n\n")

# Pathway totals
pathway_labels <- c("total_tas_signed", "total_som_signed", "total_direct", "diff_som_minus_tas")
pathway_results <- params[params$label %in% pathway_labels,
                          c("label", "est", "se", "ci.lower", "ci.upper", "pvalue")]

cat("Pathway Totals (Signed):\n")
print(pathway_results, row.names = FALSE)

cat("\n*** KEY TEST: Somatic vs TAS Pathway Difference ***\n")
diff_row <- params[params$label == "diff_som_minus_tas", ]
cat(sprintf("Difference (Somatic - TAS): %.4f, 95%% CI [%.4f, %.4f], p = %.4f\n",
            diff_row$est, diff_row$ci.lower, diff_row$ci.upper, diff_row$pvalue))
if (diff_row$pvalue < 0.05) {
  cat("CONCLUSION: Somatic pathway is SIGNIFICANTLY stronger than TAS pathway\n")
} else {
  cat("CONCLUSION: No significant difference between pathways\n")
}

cat("\n=============================================================================\n")
cat("RESULTS: INDEX OF MODERATED MEDIATION\n")
cat("=============================================================================\n\n")

index_labels <- c("index_mod_med_som", "index_mod_med_tas")
index_results <- params[params$label %in% index_labels,
                        c("label", "est", "se", "ci.lower", "ci.upper", "pvalue")]

cat("INDEX of Moderated Mediation:\n")
cat("(Tests whether indirect effect depends on IATS level)\n\n")
print(index_results, row.names = FALSE)

cat("\n*** Interpretation ***\n")
som_index <- params[params$label == "index_mod_med_som", ]
tas_index <- params[params$label == "index_mod_med_tas", ]

if (!is.na(som_index$ci.lower) && !is.na(som_index$ci.upper)) {
  if (som_index$ci.lower > 0 || som_index$ci.upper < 0) {
    cat("Somatic pathway: MODERATED (CI excludes zero)\n")
    cat("  -> IAS->Somatic->MH effect VARIES by IATS level (profile)\n")
  } else {
    cat("Somatic pathway: Not significantly moderated (CI includes zero)\n")
  }
}

if (!is.na(tas_index$ci.lower) && !is.na(tas_index$ci.upper)) {
  if (tas_index$ci.lower > 0 || tas_index$ci.upper < 0) {
    cat("TAS pathway: MODERATED (CI excludes zero)\n")
  } else {
    cat("TAS pathway: Not significantly moderated (CI includes zero)\n")
    cat("  -> IAS->TAS->MH effect is CONSISTENT across profiles\n")
  }
}

cat("\n=============================================================================\n")
cat("RESULTS: CONDITIONAL INDIRECT EFFECTS BY PROFILE\n")
cat("=============================================================================\n\n")

cond_labels <- c("cond_ias_tas_hyper", "cond_ias_som_hyper",
                 "cond_ias_tas_uncert", "cond_ias_som_uncert",
                 "cond_ias_tas_effic", "cond_ias_som_effic")
cond_results <- params[params$label %in% cond_labels,
                       c("label", "est", "se", "ci.lower", "ci.upper", "pvalue")]

cat("Conditional Indirect Effect of IAS on MH, by Profile:\n\n")

profiles <- c("Hypervigilant", "Uncertain", "Efficient")
for (prof in profiles) {
  prof_lower <- tolower(substr(prof, 1, 5))
  if (prof == "Uncertain") prof_lower <- "uncert"

  tas_row <- cond_results[cond_results$label == paste0("cond_ias_tas_", prof_lower), ]
  som_row <- cond_results[cond_results$label == paste0("cond_ias_som_", prof_lower), ]

  cat(sprintf("--- %s ---\n", prof))
  cat(sprintf("  Through TAS:     %.4f [%.4f, %.4f]\n",
              tas_row$est, tas_row$ci.lower, tas_row$ci.upper))
  cat(sprintf("  Through Somatic: %.4f [%.4f, %.4f]\n",
              som_row$est, som_row$ci.lower, som_row$ci.upper))

  # Compute percentage (using absolute values)
  abs_tas <- abs(tas_row$est)
  abs_som <- abs(som_row$est)
  total <- abs_tas + abs_som
  pct_tas <- round(abs_tas / total * 100, 1)
  pct_som <- round(abs_som / total * 100, 1)
  cat(sprintf("  Percentage: TAS = %.1f%%, Somatic = %.1f%%\n\n", pct_tas, pct_som))
}

cat("\n=============================================================================\n")
cat("RESULTS: PROFILE PATHWAY DIFFERENCE TESTS (with Bonferroni Correction)\n")
cat("=============================================================================\n\n")

diff_labels <- c("diff_hyper", "diff_uncert", "diff_effic")
diff_results <- params[params$label %in% diff_labels,
                       c("label", "est", "se", "ci.lower", "ci.upper", "pvalue")]

# Bonferroni correction for 3 tests
n_tests <- 3
alpha_corrected <- 0.05 / n_tests

cat("Test: Is TAS pathway stronger than Somatic pathway for each profile?\n")
cat("(Difference = Somatic - TAS; positive diff + sig = TAS has larger absolute effect)\n")
cat(sprintf("Bonferroni-corrected α = 0.05 / %d = %.4f\n\n", n_tests, alpha_corrected))

cat(sprintf("%-15s %10s %22s %10s %12s %s\n",
            "Profile", "Diff", "95% CI", "p-value", "p < 0.0167?", "Conclusion"))
cat(paste(rep("-", 85), collapse = ""), "\n")

for (i in 1:nrow(diff_results)) {
  prof <- gsub("diff_", "", diff_results$label[i])
  prof <- tools::toTitleCase(prof)

  p_val <- diff_results$pvalue[i]
  sig_bonf <- p_val < alpha_corrected

  # Determine conclusion
  # Both indirect effects are negative (protective). Diff = Som - TAS.
  # If diff > 0, Som is less negative than TAS, meaning TAS has larger absolute effect.
  # If diff < 0, TAS is less negative than Som, meaning Som has larger absolute effect.
  if (sig_bonf) {
    if (diff_results$est[i] > 0) {
      conclusion <- "TAS > Somatic ***"
    } else {
      conclusion <- "Somatic > TAS ***"
    }
  } else {
    conclusion <- "Balanced (n.s.)"
  }

  cat(sprintf("%-15s %10.4f [%7.4f, %7.4f] %10.4f %12s %s\n",
              prof, diff_results$est[i], diff_results$ci.lower[i],
              diff_results$ci.upper[i], p_val,
              ifelse(sig_bonf, "Yes", "No"),
              conclusion))
}

cat("\n*** = survives Bonferroni correction (p < 0.0167)\n")

# Store for later use in summary
diff_results$bonferroni_sig <- diff_results$pvalue < alpha_corrected

# =============================================================================
# 6. COMPUTE PERCENTAGES WITH BOOTSTRAP CIs
# =============================================================================

cat("\n\n=============================================================================\n")
cat("PART 2: PATHWAY PERCENTAGES WITH BOOTSTRAP CIs\n")
cat("=============================================================================\n\n")

# Extract bootstrap samples for percentage calculation
boot_est <- lavInspect(fit_boot, "boot")

# Function to compute percentages from parameter vector
compute_percentages <- function(pars, par_names) {
  # Get indices for each parameter
  get_par <- function(name) {
    idx <- which(par_names == name)
    if (length(idx) == 0) return(NA)
    pars[idx]
  }

  # Get path coefficients
  a1 <- get_par("a1"); a2 <- get_par("a2"); a3 <- get_par("a3")
  b1 <- get_par("b1"); b2 <- get_par("b2"); b3 <- get_par("b3")
  c1 <- get_par("c1"); c2 <- get_par("c2")
  d1 <- get_par("d1"); d2 <- get_par("d2"); d3 <- get_par("d3")

  # Indirect effects
  ind_ias_tas <- a1 * c1
  ind_iats_tas <- a2 * c1
  ind_int_tas <- a3 * c1
  ind_ias_som <- b1 * c2
  ind_iats_som <- b2 * c2
  ind_int_som <- b3 * c2

  # Absolute sums for percentages
  abs_tas <- abs(ind_ias_tas) + abs(ind_iats_tas) + abs(ind_int_tas)
  abs_som <- abs(ind_ias_som) + abs(ind_iats_som) + abs(ind_int_som)
  abs_direct <- abs(d1) + abs(d2) + abs(d3)

  total <- abs_tas + abs_som + abs_direct

  pct_tas <- abs_tas / total * 100
  pct_som <- abs_som / total * 100
  pct_direct <- abs_direct / total * 100

  c(pct_tas = pct_tas, pct_som = pct_som, pct_direct = pct_direct)
}

# Get parameter names
par_names <- names(coef(fit_boot))

# Apply to all bootstrap samples
boot_pcts <- t(apply(boot_est, 1, compute_percentages, par_names = par_names))

# Compute CIs
pct_results <- data.frame(
  Pathway = c("TAS (Alexithymia)", "Somatic", "Direct"),
  Estimate = colMeans(boot_pcts, na.rm = TRUE),
  CI_lower = apply(boot_pcts, 2, quantile, probs = 0.025, na.rm = TRUE),
  CI_upper = apply(boot_pcts, 2, quantile, probs = 0.975, na.rm = TRUE)
)

cat("Pathway Contribution Percentages (with 95% Bootstrap CIs):\n\n")
cat(sprintf("%-20s %10s %15s\n", "Pathway", "Estimate", "95% CI"))
cat(paste(rep("-", 50), collapse = ""), "\n")
for (i in 1:nrow(pct_results)) {
  cat(sprintf("%-20s %9.1f%% [%5.1f%%, %5.1f%%]\n",
              pct_results$Pathway[i], pct_results$Estimate[i],
              pct_results$CI_lower[i], pct_results$CI_upper[i]))
}

# =============================================================================
# 7. PROFILE-SPECIFIC PERCENTAGES
# =============================================================================

cat("\n\n=============================================================================\n")
cat("PART 3: PROFILE-SPECIFIC PATHWAY PERCENTAGES\n")
cat("=============================================================================\n\n")

# Function to compute profile-specific percentages
compute_profile_pcts <- function(pars, par_names, iats_val) {
  get_par <- function(name) {
    idx <- which(par_names == name)
    if (length(idx) == 0) return(NA)
    pars[idx]
  }

  a1 <- get_par("a1"); a3 <- get_par("a3")
  b1 <- get_par("b1"); b3 <- get_par("b3")
  c1 <- get_par("c1"); c2 <- get_par("c2")

  # Conditional indirect effects (for IAS effect at given IATS level)
  cond_tas <- (a1 + a3 * iats_val) * c1
  cond_som <- (b1 + b3 * iats_val) * c2

  # Percentages
  total <- abs(cond_tas) + abs(cond_som)
  pct_tas <- abs(cond_tas) / total * 100
  pct_som <- abs(cond_som) / total * 100

  c(pct_tas = pct_tas, pct_som = pct_som)
}

# Bootstrap for each profile
profile_pcts <- list()
for (prof in profiles) {
  if (prof == "Hypervigilant") iats_val <- hyper_iats
  else if (prof == "Uncertain") iats_val <- uncert_iats
  else iats_val <- effic_iats

  boot_prof <- t(apply(boot_est, 1, compute_profile_pcts,
                       par_names = par_names, iats_val = iats_val))

  # Handle case where result is a matrix or vector
  if (is.null(dim(boot_prof))) {
    # Single row result - convert to matrix
    boot_prof <- matrix(boot_prof, nrow = 1)
    colnames(boot_prof) <- c("pct_tas", "pct_som")
  }

  profile_pcts[[prof]] <- data.frame(
    Profile = prof,
    TAS_pct = mean(boot_prof[, 1], na.rm = TRUE),
    TAS_ci_low = quantile(boot_prof[, 1], 0.025, na.rm = TRUE),
    TAS_ci_high = quantile(boot_prof[, 1], 0.975, na.rm = TRUE),
    Som_pct = mean(boot_prof[, 2], na.rm = TRUE),
    Som_ci_low = quantile(boot_prof[, 2], 0.025, na.rm = TRUE),
    Som_ci_high = quantile(boot_prof[, 2], 0.975, na.rm = TRUE)
  )
}

profile_pct_df <- do.call(rbind, profile_pcts)
rownames(profile_pct_df) <- NULL

cat("Profile-Specific Pathway Percentages:\n")
cat("(For IAS effect on MH, conditional on profile-typical IATS)\n\n")

cat(sprintf("%-15s %20s %20s\n", "Profile", "TAS Pathway", "Somatic Pathway"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in 1:nrow(profile_pct_df)) {
  cat(sprintf("%-15s %5.1f%% [%4.1f, %4.1f] %5.1f%% [%4.1f, %4.1f]\n",
              profile_pct_df$Profile[i],
              profile_pct_df$TAS_pct[i], profile_pct_df$TAS_ci_low[i], profile_pct_df$TAS_ci_high[i],
              profile_pct_df$Som_pct[i], profile_pct_df$Som_ci_low[i], profile_pct_df$Som_ci_high[i]))
}

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

# Save to file
sink("analysis_output/12b_output.txt")
cat("=============================================================================\n")
cat("PATHWAY DOMINANCE ANALYSIS: STATISTICAL TESTING\n")
cat("=============================================================================\n\n")
cat("Sample: N =", nrow(dfc), "\n")
cat("Bootstrap iterations: 5000\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("Profile-Typical IATS Values (z-scored):\n")
cat(sprintf("  Hypervigilant: %.3f\n", hyper_iats))
cat(sprintf("  Uncertain:     %.3f\n", uncert_iats))
cat(sprintf("  Efficient:     %.3f\n\n", effic_iats))

cat("=============================================================================\n")
cat("KEY RESULTS SUMMARY\n")
cat("=============================================================================\n\n")

cat("1. OVERALL PATHWAY PERCENTAGES (with 95% CIs):\n")
for (i in 1:nrow(pct_results)) {
  cat(sprintf("   %-20s: %5.1f%% [%5.1f%%, %5.1f%%]\n",
              pct_results$Pathway[i], pct_results$Estimate[i],
              pct_results$CI_lower[i], pct_results$CI_upper[i]))
}

cat("\n2. INDEX OF MODERATED MEDIATION:\n")
cat(sprintf("   Somatic pathway: %.4f, 95%% CI [%.4f, %.4f]\n",
            som_index$est, som_index$ci.lower, som_index$ci.upper))
cat(sprintf("   TAS pathway:     %.4f, 95%% CI [%.4f, %.4f]\n",
            tas_index$est, tas_index$ci.lower, tas_index$ci.upper))

cat("\n3. PROFILE-SPECIFIC PATHWAY PERCENTAGES (with 95% CIs):\n")
for (i in 1:nrow(profile_pct_df)) {
  cat(sprintf("   %s: TAS = %.1f%% [%.1f, %.1f], Somatic = %.1f%% [%.1f, %.1f]\n",
              profile_pct_df$Profile[i],
              profile_pct_df$TAS_pct[i], profile_pct_df$TAS_ci_low[i], profile_pct_df$TAS_ci_high[i],
              profile_pct_df$Som_pct[i], profile_pct_df$Som_ci_low[i], profile_pct_df$Som_ci_high[i]))
}

cat("\n4. WITHIN-PROFILE PATHWAY DOMINANCE TESTS (Bonferroni-corrected):\n")
cat("   Bonferroni α = 0.05 / 3 = 0.0167\n\n")
cat(sprintf("   %-15s %10s %12s %s\n", "Profile", "p-value", "Significant?", "Conclusion"))
cat("   ", paste(rep("-", 55), collapse = ""), "\n")

# Get the diff_results with p-values
for (i in 1:nrow(diff_results)) {
  prof <- gsub("diff_", "", diff_results$label[i])
  prof <- tools::toTitleCase(prof)
  p_val <- diff_results$pvalue[i]
  sig_bonf <- p_val < 0.0167

  # Both effects are negative (protective). Diff = Som - TAS.
  # Positive diff + sig means TAS has larger absolute effect
  if (sig_bonf) {
    if (diff_results$est[i] > 0) {
      conclusion <- "TAS > Somatic ***"
    } else {
      conclusion <- "Somatic > TAS ***"
    }
  } else {
    conclusion <- "Balanced (n.s.)"
  }

  cat(sprintf("   %-15s %10.4f %12s %s\n",
              prof, p_val,
              ifelse(sig_bonf, "Yes", "No"),
              conclusion))
}

cat("\n5. OVERALL PATHWAY DIFFERENCE (TAS vs Somatic):\n")
cat("   Somatic > TAS overall: ")
if (diff_row$pvalue < 0.05) {
  cat(sprintf("YES (p = %.4f)\n", diff_row$pvalue))
} else {
  cat(sprintf("NO (p = %.4f)\n", diff_row$pvalue))
}

cat("\n*** = survives Bonferroni correction (p < 0.0167)\n")

cat("\n6. INDIVIDUAL INDIRECT EFFECTS:\n")
for (i in 1:nrow(indirect)) {
  cat(sprintf("   %s: %.4f [%.4f, %.4f], p = %.4f\n",
              indirect$label[i], indirect$est[i],
              indirect$ci.lower[i], indirect$ci.upper[i], indirect$pvalue[i]))
}

cat("\n7. PREDICTOR-SPECIFIC PATHWAY SPLITS:\n")
iats_som_val <- abs(indirect$est[indirect$label == "ind_iats_som"])
iats_tas_val <- abs(indirect$est[indirect$label == "ind_iats_tas"])
cat(sprintf("   IATS through Somatic: %.1f%%\n", iats_som_val / (iats_som_val + iats_tas_val) * 100))
cat(sprintf("   IATS through TAS: %.1f%%\n", iats_tas_val / (iats_som_val + iats_tas_val) * 100))
ias_tas_val <- abs(indirect$est[indirect$label == "ind_ias_tas"])
ias_som_val <- abs(indirect$est[indirect$label == "ind_ias_som"])
cat(sprintf("   IAS through TAS: %.1f%%\n", ias_tas_val / (ias_tas_val + ias_som_val) * 100))
cat(sprintf("   IAS through Somatic: %.1f%%\n", ias_som_val / (ias_tas_val + ias_som_val) * 100))

cat(sprintf("\n8. TOTAL MEDIATION PERCENTAGE:\n"))
cat(sprintf("   TAS pathway: %.1f%%\n", pct_results$Estimate[1]))
cat(sprintf("   Somatic pathway: %.1f%%\n", pct_results$Estimate[2]))
cat(sprintf("   Total mediated: %.1f%%\n", pct_results$Estimate[1] + pct_results$Estimate[2]))
cat(sprintf("   Direct: %.1f%%\n", pct_results$Estimate[3]))

sink()

cat("\n\nResults saved to: analysis_output/12b_output.txt\n")

# =============================================================================
# 9. VISUALIZATION
# =============================================================================

cat("\nCreating visualization...\n")

# Prepare data for plotting
plot_data <- data.frame(
  Profile = factor(profile_pct_df$Profile,
                   levels = c("Hypervigilant", "Uncertain", "Efficient")),
  TAS = profile_pct_df$TAS_pct,
  TAS_low = profile_pct_df$TAS_ci_low,
  TAS_high = profile_pct_df$TAS_ci_high,
  Somatic = profile_pct_df$Som_pct,
  Som_low = profile_pct_df$Som_ci_low,
  Som_high = profile_pct_df$Som_ci_high
)

# Long format for stacked bar
plot_long <- pivot_longer(plot_data,
                          cols = c(TAS, Somatic),
                          names_to = "Pathway",
                          values_to = "Percentage")
plot_long$Pathway <- factor(plot_long$Pathway, levels = c("TAS", "Somatic"))

# Create plot with proper statistical backing
p_dominance <- ggplot(plot_long, aes(x = Profile, y = Percentage, fill = Pathway)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("TAS" = "#E57373", "Somatic" = "#F59E0B")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  labs(
    title = "Pathway Dominance by Profile (with Bootstrap 95% CIs)",
    subtitle = "Conditional indirect effect of IAS on MH through each pathway",
    x = NULL,
    y = "Percentage of Total Indirect Effect",
    caption = sprintf("N = %d, Bootstrap = 5000", nrow(dfc))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  ) +
  # Add percentage labels with CIs
  geom_text(aes(label = sprintf("%.0f%%", Percentage)),
            position = position_stack(vjust = 0.5),
            color = "white", size = 5, fontface = "bold")

ggsave("plots/supplementary_plots/3_s_2_pathway_dominance_statistical.png", p_dominance,
       width = 8, height = 6, dpi = 300, bg = "white")

cat("Figure saved: plots/supplementary_plots/3_s_2_pathway_dominance_statistical.png\n")

cat("\n=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
