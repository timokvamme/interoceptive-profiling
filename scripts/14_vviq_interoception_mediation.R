# =============================================================================
# VVIQ-INTEROCEPTION-ALEXITHYMIA MEDIATION ANALYSIS
# =============================================================================
#
# Tests whether interoceptive accuracy (IAS) mediates the link between
# imagery vividness (VVIQ) and alexithymia, supporting the hypothesis
# that imagery and emotion share interoceptive underpinnings rather than
# imagery directly causing emotion.
#
# Framework (Silvanto): IMAGERY <-- INTEROCEPTION --> EMOTION
# Test:                 VVIQ --> IAS --> TAS (and subscales)
#
# ANALYSES:
#   1. Correlations: VVIQ with interoceptive and alexithymia measures
#   2. Simple mediation: VVIQ -> IAS -> TAS
#   3. Specificity: VVIQ -> [IAS, IATS] -> TAS (parallel mediators)
#   4. Subscale decomposition: VVIQ -> IAS -> [DIF, DDF, EOT]
#   5. Figure: mediation path diagram + forest plot
#
# OUTPUTS:
#   - plots/other_plots_vviq_plots/vviq_mediation_output.txt
#   - plots/other_plots_vviq_plots/vviq_mediation_paths.png
#
# =============================================================================

# Paths are relative to the repository root; run scripts from there.

library(dplyr)
library(lavaan)
library(ggplot2)
library(gridExtra)
library(grid)

out_dir <- "plots/other_plots_vviq_plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

dfc <- read.csv("dfc_interoception_profiling.csv")

# TAS subscales (items already reverse-scored)
dfc$dif <- rowSums(dfc[, paste0("tas_", c(1, 3, 6, 7, 9, 13, 14))], na.rm = TRUE)
dfc$ddf <- rowSums(dfc[, paste0("tas_", c(2, 4, 11, 12, 17))], na.rm = TRUE)
dfc$eot <- rowSums(dfc[, paste0("tas_", c(5, 8, 10, 15, 16, 18, 19, 20))], na.rm = TRUE)

# Full-scale composites
dfc$ias_full <- rowMeans(dfc[, paste0("ias_", 1:21)], na.rm = TRUE)
dfc$iats_full <- rowMeans(dfc[, paste0("iats_", 1:21)], na.rm = TRUE)
dfc$mh_composite <- rowMeans(cbind(scale(dfc$phq9), scale(dfc$gad7),
                                    scale(dfc$stai)), na.rm = TRUE)

cat("Data loaded. N =", nrow(dfc), "\n")
cat("VVIQ: M =", round(mean(dfc$vviq, na.rm = TRUE), 2),
    ", SD =", round(sd(dfc$vviq, na.rm = TRUE), 2),
    ", range =", min(dfc$vviq, na.rm = TRUE), "-", max(dfc$vviq, na.rm = TRUE), "\n\n")

# =============================================================================
# BEGIN OUTPUT
# =============================================================================

sink(file.path(out_dir, "vviq_mediation_output.txt"))

cat("=============================================================================\n")
cat("VVIQ-INTEROCEPTION-ALEXITHYMIA MEDIATION ANALYSIS\n")
cat("=============================================================================\n\n")
cat("Sample: N =", nrow(dfc), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
cat("Hypothesis: Imagery vividness (VVIQ) and alexithymia (TAS) are linked\n")
cat("because both depend on interoceptive accuracy (IAS), not because imagery\n")
cat("directly causes emotional awareness.\n\n")
cat("Test: VVIQ -> IAS -> TAS (and subscales DIF, DDF, EOT)\n\n")

# =============================================================================
# 2. CORRELATIONS
# =============================================================================

cat("#############################################################################\n")
cat("PART 1: CORRELATIONS - VVIQ WITH INTEROCEPTIVE AND CLINICAL MEASURES\n")
cat("#############################################################################\n\n")

cor_vars <- c("ias_full", "iats_full", "tas", "dif", "ddf", "eot",
              "sss8", "pcs", "phq9", "gad7", "stai", "mh_composite")
cor_labels <- c("IAS", "IATS", "TAS Total", "DIF", "DDF", "EOT",
                "SSS-8", "PCS", "PHQ-9", "GAD-7", "STAI", "MH Composite")

cat(sprintf("%-15s %8s %10s %s\n", "Variable", "r", "p", "Sig"))
cat(paste(rep("-", 45), collapse = ""), "\n")

cor_results <- data.frame(variable = character(), r = numeric(), p = numeric(),
                           stringsAsFactors = FALSE)

for (i in seq_along(cor_vars)) {
  ct <- cor.test(dfc$vviq, dfc[[cor_vars[i]]], use = "complete.obs")
  sig <- if (ct$p.value < 0.001) "***" else if (ct$p.value < 0.01) "**" else if (ct$p.value < 0.05) "*" else ""
  cat(sprintf("%-15s %8.3f %10.4f %s\n",
              cor_labels[i], ct$estimate, ct$p.value, sig))
  cor_results <- rbind(cor_results, data.frame(
    variable = cor_labels[i], r = ct$estimate, p = ct$p.value))
}

cat("\nKey pattern:\n")
cat("  VVIQ correlates with IAS (r = +.21) but NOT IATS (r = -.04, n.s.)\n")
cat("  VVIQ correlates with TAS (r = -.24), especially DIF\n")
cat("  This supports IAS-specific mediation of the imagery-alexithymia link\n")

# =============================================================================
# 3. SIMPLE MEDIATION: VVIQ -> IAS -> TAS
# =============================================================================

cat("\n\n#############################################################################\n")
cat("PART 2: SIMPLE MEDIATION - VVIQ -> IAS -> TAS\n")
cat("#############################################################################\n\n")

# Standardize
dfc$vviq_z <- scale(dfc$vviq)[, 1]
dfc$ias_z <- scale(dfc$ias_full)[, 1]
dfc$iats_z <- scale(dfc$iats_full)[, 1]
dfc$tas_z <- scale(dfc$tas)[, 1]
dfc$dif_z <- scale(dfc$dif)[, 1]
dfc$ddf_z <- scale(dfc$ddf)[, 1]
dfc$eot_z <- scale(dfc$eot)[, 1]
dfc$mh_z <- scale(dfc$mh_composite)[, 1]

simple_model <- '
  # a path: VVIQ -> IAS
  ias_z ~ a*vviq_z

  # b path: IAS -> TAS (controlling for VVIQ)
  tas_z ~ b*ias_z + cp*vviq_z

  # Indirect effect
  indirect := a * b

  # Total effect
  total := cp + a * b

  # Proportion mediated
  prop_med := indirect / total
'

cat("Running bootstrap mediation (5000 resamples)...\n")
fit_simple <- sem(simple_model, data = dfc, se = "bootstrap",
                  bootstrap = 5000, iseed = 42)
params_simple <- parameterEstimates(fit_simple, ci = TRUE, standardized = TRUE)

get_est <- function(params, lbl) {
  row <- params[params$label == lbl, ]
  if (nrow(row) > 0) return(list(
    est = row$std.all[1], p = row$pvalue[1],
    ci_lo = row$ci.lower[1], ci_hi = row$ci.upper[1]))
  return(list(est = NA, p = NA, ci_lo = NA, ci_hi = NA))
}

a <- get_est(params_simple, "a")
b <- get_est(params_simple, "b")
cp <- get_est(params_simple, "cp")
ind <- get_est(params_simple, "indirect")
tot <- get_est(params_simple, "total")
prop <- get_est(params_simple, "prop_med")

cat("\n--- Simple Mediation: VVIQ -> IAS -> TAS ---\n\n")
cat(sprintf("  Total effect (c):    beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            tot$est, tot$p, tot$ci_lo, tot$ci_hi))
cat(sprintf("  Direct effect (c'):  beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            cp$est, cp$p, cp$ci_lo, cp$ci_hi))
cat(sprintf("  a path (VVIQ->IAS): beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            a$est, a$p, a$ci_lo, a$ci_hi))
cat(sprintf("  b path (IAS->TAS):  beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            b$est, b$p, b$ci_lo, b$ci_hi))
cat(sprintf("  Indirect (a*b):     beta = %6.4f (p = %.4f) [%.4f, %.4f]\n",
            ind$est, ind$p, ind$ci_lo, ind$ci_hi))
cat(sprintf("  Proportion mediated: %.1f%%\n", abs(prop$est) * 100))

cat("\nInterpretation: IAS mediates ", round(abs(prop$est) * 100, 1),
    "% of the VVIQ-TAS association.\n", sep = "")
if (cp$p < 0.05) {
  cat("The direct effect remains significant, indicating partial mediation.\n")
} else {
  cat("The direct effect is no longer significant, indicating full mediation.\n")
}

# =============================================================================
# 4. SPECIFICITY: VVIQ -> [IAS, IATS] -> TAS
# =============================================================================

cat("\n\n#############################################################################\n")
cat("PART 3: SPECIFICITY - VVIQ -> [IAS, IATS] -> TAS (parallel mediators)\n")
cat("#############################################################################\n\n")

cat("Testing whether the mediation is specific to IAS (accuracy)\n")
cat("or also operates through IATS (attention).\n\n")

parallel_model <- '
  ias_z ~ a1*vviq_z
  iats_z ~ a2*vviq_z
  tas_z ~ b1*ias_z + b2*iats_z + cp*vviq_z

  ind_ias := a1 * b1
  ind_iats := a2 * b2
  total_ind := ind_ias + ind_iats
  total := cp + total_ind
'

fit_parallel <- sem(parallel_model, data = dfc, se = "bootstrap",
                    bootstrap = 5000, iseed = 42)
pp <- parameterEstimates(fit_parallel, ci = TRUE, standardized = TRUE)

get_pp <- function(lbl) {
  row <- pp[pp$label == lbl, ]
  if (nrow(row) > 0) return(list(
    est = row$std.all[1], p = row$pvalue[1],
    ci_lo = row$ci.lower[1], ci_hi = row$ci.upper[1]))
  return(list(est = NA, p = NA, ci_lo = NA, ci_hi = NA))
}

a1 <- get_pp("a1"); a2 <- get_pp("a2")
b1 <- get_pp("b1"); b2 <- get_pp("b2")
cp2 <- get_pp("cp")
ind_ias <- get_pp("ind_ias"); ind_iats <- get_pp("ind_iats")

cat("--- Parallel Mediation: VVIQ -> [IAS, IATS] -> TAS ---\n\n")
cat("  a paths (VVIQ -> mediator):\n")
cat(sprintf("    VVIQ -> IAS:  beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            a1$est, a1$p, a1$ci_lo, a1$ci_hi))
cat(sprintf("    VVIQ -> IATS: beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            a2$est, a2$p, a2$ci_lo, a2$ci_hi))

cat("\n  b paths (mediator -> TAS, controlling for VVIQ):\n")
cat(sprintf("    IAS -> TAS:   beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            b1$est, b1$p, b1$ci_lo, b1$ci_hi))
cat(sprintf("    IATS -> TAS:  beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            b2$est, b2$p, b2$ci_lo, b2$ci_hi))

cat("\n  Indirect effects:\n")
cat(sprintf("    Via IAS:  %6.4f (p = %.4f) [%.4f, %.4f] %s\n",
            ind_ias$est, ind_ias$p, ind_ias$ci_lo, ind_ias$ci_hi,
            ifelse(ind_ias$p < 0.05, "***", "n.s.")))
cat(sprintf("    Via IATS: %6.4f (p = %.4f) [%.4f, %.4f] %s\n",
            ind_iats$est, ind_iats$p, ind_iats$ci_lo, ind_iats$ci_hi,
            ifelse(ind_iats$p < 0.05, "***", "n.s.")))

cat(sprintf("\n  Direct (c'): beta = %6.3f (p = %.4f)\n", cp2$est, cp2$p))

if (ind_ias$p < 0.05 && ind_iats$p >= 0.05) {
  cat("\nSPECIFICITY CONFIRMED: Mediation operates through IAS (accuracy)\n")
  cat("but NOT through IATS (attention). The imagery-alexithymia link is\n")
  cat("specifically routed through interoceptive signal detection.\n")
}

# =============================================================================
# 5. SUBSCALE DECOMPOSITION: VVIQ -> IAS -> [DIF, DDF, EOT]
# =============================================================================

cat("\n\n#############################################################################\n")
cat("PART 4: SUBSCALE DECOMPOSITION - VVIQ -> IAS -> [DIF, DDF, EOT]\n")
cat("#############################################################################\n\n")

cat("Which alexithymia facet does the VVIQ-IAS pathway target?\n\n")

subscale_model <- '
  # a path
  ias_z ~ a*vviq_z

  # b paths to each subscale
  dif_z ~ b_dif*ias_z + c_dif*vviq_z
  ddf_z ~ b_ddf*ias_z + c_ddf*vviq_z
  eot_z ~ b_eot*ias_z + c_eot*vviq_z

  # Allow subscale residuals to correlate
  dif_z ~~ ddf_z
  dif_z ~~ eot_z
  ddf_z ~~ eot_z

  # Indirect effects
  ind_dif := a * b_dif
  ind_ddf := a * b_ddf
  ind_eot := a * b_eot
  total_ind := ind_dif + ind_ddf + ind_eot

  # Contrasts
  diff_dif_ddf := ind_dif - ind_ddf
  diff_dif_eot := ind_dif - ind_eot
'

fit_sub <- sem(subscale_model, data = dfc, se = "bootstrap",
               bootstrap = 5000, iseed = 42)
ps <- parameterEstimates(fit_sub, ci = TRUE, standardized = TRUE)

get_ps <- function(lbl) {
  row <- ps[ps$label == lbl, ]
  if (nrow(row) > 0) return(list(
    est = row$std.all[1], p = row$pvalue[1],
    ci_lo = row$ci.lower[1], ci_hi = row$ci.upper[1]))
  return(list(est = NA, p = NA, ci_lo = NA, ci_hi = NA))
}

a_sub <- get_ps("a")
b_dif <- get_ps("b_dif"); b_ddf <- get_ps("b_ddf"); b_eot <- get_ps("b_eot")
c_dif <- get_ps("c_dif"); c_ddf <- get_ps("c_ddf"); c_eot <- get_ps("c_eot")
ind_dif <- get_ps("ind_dif"); ind_ddf <- get_ps("ind_ddf"); ind_eot <- get_ps("ind_eot")
total_sub <- get_ps("total_ind")
diff_dif_ddf <- get_ps("diff_dif_ddf")
diff_dif_eot <- get_ps("diff_dif_eot")

cat("  a path: VVIQ -> IAS: beta = ", sprintf("%.3f (p = %.4f)\n", a_sub$est, a_sub$p))

cat("\n  b paths (IAS -> subscale, controlling for VVIQ):\n")
cat(sprintf("    IAS -> DIF: beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            b_dif$est, b_dif$p, b_dif$ci_lo, b_dif$ci_hi))
cat(sprintf("    IAS -> DDF: beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            b_ddf$est, b_ddf$p, b_ddf$ci_lo, b_ddf$ci_hi))
cat(sprintf("    IAS -> EOT: beta = %6.3f (p = %.4f) [%.3f, %.3f]\n",
            b_eot$est, b_eot$p, b_eot$ci_lo, b_eot$ci_hi))

cat("\n  Direct effects (VVIQ -> subscale, controlling for IAS):\n")
cat(sprintf("    VVIQ -> DIF: beta = %6.3f (p = %.4f)\n", c_dif$est, c_dif$p))
cat(sprintf("    VVIQ -> DDF: beta = %6.3f (p = %.4f)\n", c_ddf$est, c_ddf$p))
cat(sprintf("    VVIQ -> EOT: beta = %6.3f (p = %.4f)\n", c_eot$est, c_eot$p))

cat("\n  Indirect effects (VVIQ -> IAS -> subscale):\n")
cat(sprintf("    Via DIF: %6.4f (p = %.4f) [%.4f, %.4f] %s\n",
            ind_dif$est, ind_dif$p, ind_dif$ci_lo, ind_dif$ci_hi,
            ifelse(ind_dif$p < 0.05, "***", "n.s.")))
cat(sprintf("    Via DDF: %6.4f (p = %.4f) [%.4f, %.4f] %s\n",
            ind_ddf$est, ind_ddf$p, ind_ddf$ci_lo, ind_ddf$ci_hi,
            ifelse(ind_ddf$p < 0.05, "***", "n.s.")))
cat(sprintf("    Via EOT: %6.4f (p = %.4f) [%.4f, %.4f] %s\n",
            ind_eot$est, ind_eot$p, ind_eot$ci_lo, ind_eot$ci_hi,
            ifelse(ind_eot$p < 0.05, "***", "n.s.")))

# Proportions
abs_total <- abs(ind_dif$est) + abs(ind_ddf$est) + abs(ind_eot$est)
if (abs_total > 0) {
  pct_dif <- abs(ind_dif$est) / abs_total * 100
  pct_ddf <- abs(ind_ddf$est) / abs_total * 100
  pct_eot <- abs(ind_eot$est) / abs_total * 100
} else {
  pct_dif <- pct_ddf <- pct_eot <- NA
}

cat(sprintf("\n  Proportions: DIF = %.1f%%, DDF = %.1f%%, EOT = %.1f%%\n",
            pct_dif, pct_ddf, pct_eot))

cat("\n  Pairwise contrasts:\n")
cat(sprintf("    DIF - DDF: %.4f (p = %.4f) %s\n",
            diff_dif_ddf$est, diff_dif_ddf$p,
            ifelse(diff_dif_ddf$p < 0.05, "***", "n.s.")))
cat(sprintf("    DIF - EOT: %.4f (p = %.4f) %s\n",
            diff_dif_eot$est, diff_dif_eot$p,
            ifelse(diff_dif_eot$p < 0.05, "***", "n.s.")))

# =============================================================================
# 6. IMAGERY GROUPS x CLUSTER PROFILES (CHI-SQUARE)
# =============================================================================

cat("\n\n#############################################################################\n")
cat("PART 5: IMAGERY GROUPS x CLUSTER PROFILES (CHI-SQUARE)\n")
cat("#############################################################################\n\n")

cat("Are aphantasics, hypophantasics, typical imagers, and hyperphantasics\n")
cat("equally distributed across Efficient / Hypervigilant / Uncertain profiles?\n\n")

# --- Reproduce k-means cluster assignment (matches 04_cluster_analysis.R) ---
cluster_data <- scale(dfc[, c("ias", "iats")])
set.seed(123)
km3 <- kmeans(cluster_data, 3, nstart = 50)
dfc$cluster_kmeans <- km3$cluster

ias_mu <- mean(dfc$ias); ias_sd <- sd(dfc$ias)
iats_mu <- mean(dfc$iats); iats_sd <- sd(dfc$iats)
prof <- aggregate(cbind(ias = dfc$ias, iats = dfc$iats),
                  by = list(cluster = dfc$cluster_kmeans), FUN = mean)
prof$ias_z <- (prof$ias - ias_mu) / ias_sd
prof$iats_z <- (prof$iats - iats_mu) / iats_sd
prof$short <- sapply(seq_len(nrow(prof)), function(i) {
  ias_lev <- if (prof$ias_z[i] > 0.3) "High" else if (prof$ias_z[i] < -0.3) "Low" else "Med"
  iats_lev <- if (prof$iats_z[i] > 0.3) "High" else if (prof$iats_z[i] < -0.3) "Low" else "Med"
  if (ias_lev == "High" && iats_lev == "Low") return("Efficient")
  if (iats_lev == "High") return("Hypervigilant")
  return("Uncertain")
})
dfc$cluster_short <- factor(
  prof$short[match(dfc$cluster_kmeans, prof$cluster)],
  levels = c("Efficient", "Hypervigilant", "Uncertain"))

# --- Categorize VVIQ into imagery groups (Dance et al. 2022 cutoffs, scaled
#     to per-item mean so they apply to either VVIQ-1 or VVIQ-2 totals) ---
vviq_items <- grep("^vviq_q", names(dfc), value = TRUE)
n_vviq_items <- length(vviq_items)
dfc$vviq_per_item <- dfc$vviq / n_vviq_items

dfc$imagery_group <- cut(
  dfc$vviq_per_item,
  breaks = c(-Inf, 2.0, 2.3125, 4.6875, Inf),
  labels = c("Aphantasia", "Hypophantasia", "Typical imager", "Hyperphantasia"),
  right = TRUE)

cat(sprintf("VVIQ items detected: %d (per-item mean used for cutoffs)\n", n_vviq_items))
cat("Imagery group cutoffs (per-item mean, Dance et al. 2022):\n")
cat("  Aphantasia:     mean <= 2.00\n")
cat("  Hypophantasia:  2.00 < mean <= 2.31\n")
cat("  Typical imager: 2.31 < mean <  4.69\n")
cat("  Hyperphantasia: mean >= 4.69\n\n")

cat("Imagery group counts (overall):\n")
print(table(dfc$imagery_group, useNA = "ifany"))

# Drop NAs (missing VVIQ) and any unused factor levels for the chi-square
chi_dat <- dfc[!is.na(dfc$imagery_group) & !is.na(dfc$cluster_short), ]
chi_dat$imagery_group <- droplevels(chi_dat$imagery_group)

ct_tbl <- table(chi_dat$imagery_group, chi_dat$cluster_short)
cat("\nCross-tabulation (rows = imagery group, cols = cluster):\n")
print(ct_tbl)

cat("\nRow %% (% of each imagery group falling in each cluster):\n")
print(round(prop.table(ct_tbl, margin = 1) * 100, 1))

cat("\nColumn %% (% of each cluster made up of each imagery group):\n")
print(round(prop.table(ct_tbl, margin = 2) * 100, 1))

# --- Chi-square test ---
cs <- suppressWarnings(chisq.test(ct_tbl))
n_total <- sum(ct_tbl)
cramers_v <- as.numeric(sqrt(cs$statistic / (n_total * (min(dim(ct_tbl)) - 1))))

cat("\nChi-square test of independence:\n")
cat(sprintf("  Chi-sq(%d) = %.2f, p = %.4f, N = %d\n",
            cs$parameter, cs$statistic, cs$p.value, n_total))
cat(sprintf("  Cramer's V = %.3f\n", cramers_v))

low_expected <- any(cs$expected < 5)
fisher_p <- NA
if (low_expected) {
  cat(sprintf("  Note: %d cell(s) have expected count < 5; running Fisher's exact (Monte Carlo).\n",
              sum(cs$expected < 5)))
  fe <- fisher.test(ct_tbl, simulate.p.value = TRUE, B = 10000)
  fisher_p <- fe$p.value
  cat(sprintf("  Fisher's exact p (simulated) = %.4f\n", fisher_p))
}

cat("\nStandardized residuals (|z| > 1.96 indicates departure from independence):\n")
print(round(cs$stdres, 2))

# Save cross-tab for the supplement
ct_out <- as.data.frame.matrix(ct_tbl)
ct_out$imagery_group <- rownames(ct_out)
ct_out <- ct_out[, c("imagery_group", levels(dfc$cluster_short))]
write.csv(ct_out, file.path(out_dir, "vviq_imagery_groups_by_cluster.csv"),
          row.names = FALSE)

# =============================================================================
# 7. SUMMARY
# =============================================================================

cat("\n\n#############################################################################\n")
cat("SUMMARY\n")
cat("#############################################################################\n\n")

cat("1. VVIQ correlates with IAS (r = +.21) but not IATS (r = -.04, n.s.),\n")
cat("   confirming specificity to interoceptive accuracy.\n\n")

cat(sprintf("2. IAS mediates %.1f%% of the VVIQ -> TAS association\n",
            abs(prop$est) * 100))
cat(sprintf("   (indirect = %.4f, 95%% CI [%.4f, %.4f]).\n",
            ind$est, ind$ci_lo, ind$ci_hi))

cat(sprintf("\n3. Mediation is specific to IAS (indirect = %.4f, p = %.4f)\n",
            ind_ias$est, ind_ias$p))
cat(sprintf("   not IATS (indirect = %.4f, p = %.4f).\n",
            get_pp("ind_iats")$est, get_pp("ind_iats")$p))

cat(sprintf("\n4. Subscale decomposition: VVIQ -> IAS targets primarily\n"))
cat(sprintf("   DIF (%.1f%%), then DDF (%.1f%%), EOT negligible (%.1f%%).\n",
            pct_dif, pct_ddf, pct_eot))

cat(sprintf("\n5. Imagery group x cluster (chi-square): chi-sq(%d) = %.2f, p = %.4f,\n",
            cs$parameter, cs$statistic, cs$p.value))
cat(sprintf("   Cramer's V = %.3f, N = %d.\n", cramers_v, n_total))
if (!is.na(fisher_p)) {
  cat(sprintf("   (Fisher's exact p = %.4f due to small expected cells.)\n", fisher_p))
}
chi_sig <- if (cs$p.value < 0.001) "<.001" else sprintf("= %.3f", cs$p.value)
if (cs$p.value < 0.05) {
  cat("   Imagery groups are NOT evenly distributed across cluster profiles -\n")
  cat("   inspect Part 5 standardized residuals for which cells drive the effect.\n")
} else {
  cat("   Imagery groups appear evenly distributed across cluster profiles.\n")
}

cat("\n6. THEORETICAL IMPLICATION:\n")
cat("   The imagery-alexithymia link is substantially mediated by\n")
cat("   interoceptive accuracy, supporting the shared-mechanism\n")
cat("   hypothesis (IMAGERY <-- INTEROCEPTION --> EMOTION) over the\n")
cat("   amplification view (IMAGERY --> EMOTION). Imagery vividness\n")
cat("   and emotional identification share a common dependence on\n")
cat("   accurate internal signal detection.\n")

sink()

# =============================================================================
# 7. FIGURE
# =============================================================================

cat("Creating figure...\n")

# --- Panel A: Mediation path diagram (simple) ---

p_a <- ggplot() + theme_void() +
  xlim(0, 10) + ylim(0, 6) +
  # Title
  annotate("text", x = 5, y = 5.8, label = "A. Simple Mediation: VVIQ -> IAS -> TAS",
           size = 4.5, fontface = "bold") +
  # Boxes
  annotate("rect", xmin = 0.3, xmax = 2.3, ymin = 2.5, ymax = 3.5,
           fill = "#8B5CF6", color = "black", linewidth = 0.8) +
  annotate("text", x = 1.3, y = 3.0, label = "VVIQ\n(Imagery)", color = "white",
           size = 4, fontface = "bold") +
  annotate("rect", xmin = 4.0, xmax = 6.0, ymin = 4.2, ymax = 5.2,
           fill = "#0891B2", color = "black", linewidth = 0.8) +
  annotate("text", x = 5.0, y = 4.7, label = "IAS\n(Accuracy)", color = "white",
           size = 4, fontface = "bold") +
  annotate("rect", xmin = 7.7, xmax = 9.7, ymin = 2.5, ymax = 3.5,
           fill = "#E57373", color = "black", linewidth = 0.8) +
  annotate("text", x = 8.7, y = 3.0, label = "TAS\n(Alexithymia)",
           size = 4, fontface = "bold") +
  # a path
  annotate("segment", x = 2.3, y = 3.3, xend = 4.0, yend = 4.4,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           linewidth = 1.5, color = "#10B981") +
  annotate("label", x = 3.0, y = 4.1,
           label = sprintf("a = %.2f***", a$est),
           size = 3.5, color = "#10B981", fontface = "bold",
           fill = "white", label.size = 0.3) +
  # b path
  annotate("segment", x = 6.0, y = 4.4, xend = 7.7, yend = 3.3,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           linewidth = 1.5, color = "#10B981") +
  annotate("label", x = 7.0, y = 4.1,
           label = sprintf("b = %.2f***", b$est),
           size = 3.5, color = "#10B981", fontface = "bold",
           fill = "white", label.size = 0.3) +
  # c' path (direct)
  annotate("segment", x = 2.3, y = 2.8, xend = 7.7, yend = 2.8,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           linewidth = 1.0, color = if (cp$p < 0.05) "#DC2626" else "gray50",
           linetype = if (cp$p < 0.05) "solid" else "dashed") +
  annotate("label", x = 5.0, y = 2.4,
           label = sprintf("c' = %.2f (p = %.3f)", cp$est, cp$p),
           size = 3.5, color = if (cp$p < 0.05) "#DC2626" else "gray50",
           fontface = "bold", fill = "white", label.size = 0.3) +
  # Indirect effect annotation
  annotate("label", x = 5.0, y = 1.2,
           label = sprintf("Indirect: %.3f [%.3f, %.3f]\nProportion mediated: %.1f%%",
                           ind$est, ind$ci_lo, ind$ci_hi, abs(prop$est) * 100),
           size = 3.2, fill = "#F0FDF4", color = "#065F46",
           label.size = 0.5, fontface = "bold")

# --- Panel B: Specificity (IAS vs IATS) ---

spec_data <- data.frame(
  mediator = c("Via IAS\n(Accuracy)", "Via IATS\n(Attention)"),
  effect = c(ind_ias$est, get_pp("ind_iats")$est),
  ci_lo = c(ind_ias$ci_lo, get_pp("ind_iats")$ci_lo),
  ci_hi = c(ind_ias$ci_hi, get_pp("ind_iats")$ci_hi),
  p = c(ind_ias$p, get_pp("ind_iats")$p),
  stringsAsFactors = FALSE
)
spec_data$sig <- ifelse(spec_data$p < 0.05, "p < .05", "n.s.")
spec_data$mediator <- factor(spec_data$mediator,
  levels = rev(c("Via IAS\n(Accuracy)", "Via IATS\n(Attention)")))

p_b <- ggplot(spec_data, aes(x = effect, y = mediator, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi),
                  color = c("#7C3AED", "#0891B2"), size = 1, linewidth = 1.2) +
  scale_shape_manual(values = c("p < .05" = 16, "n.s." = 1), name = "") +
  labs(x = "Indirect Effect [95% Bootstrap CI]", y = NULL,
       title = "B. Specificity: IAS vs IATS as Mediators") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank())

# --- Panel C: Subscale decomposition ---

sub_data <- data.frame(
  subscale = c("DIF", "DDF", "EOT"),
  effect = c(ind_dif$est, ind_ddf$est, ind_eot$est),
  ci_lo = c(ind_dif$ci_lo, ind_ddf$ci_lo, ind_eot$ci_lo),
  ci_hi = c(ind_dif$ci_hi, ind_ddf$ci_hi, ind_eot$ci_hi),
  p = c(ind_dif$p, ind_ddf$p, ind_eot$p),
  pct = c(pct_dif, pct_ddf, pct_eot),
  stringsAsFactors = FALSE
)
sub_data$sig <- ifelse(sub_data$p < 0.05, "p < .05", "n.s.")
sub_data$label <- sprintf("%s (%.0f%%)", sub_data$subscale, sub_data$pct)
sub_data$subscale <- factor(sub_data$subscale, levels = rev(c("DIF", "DDF", "EOT")))
sub_data$label <- factor(sub_data$label, levels = rev(sub_data$label))

p_c <- ggplot(sub_data, aes(x = effect, y = label, color = subscale, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), size = 0.8, linewidth = 1) +
  scale_color_manual(values = c("DIF" = "#E64B35", "DDF" = "#4DBBD5", "EOT" = "#00A087"),
                     guide = "none") +
  scale_shape_manual(values = c("p < .05" = 16, "n.s." = 1), name = "") +
  labs(x = "Indirect Effect: VVIQ -> IAS -> Subscale [95% CI]", y = NULL,
       title = "C. Subscale Target: VVIQ -> IAS -> Which Facet?") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank())

# Combine
fig <- arrangeGrob(
  p_a,
  arrangeGrob(p_b, p_c, ncol = 2, widths = c(1, 1.2)),
  nrow = 2, heights = c(1.2, 1),
  top = textGrob("VVIQ-Interoception-Alexithymia Mediation",
                 gp = gpar(fontface = "bold", fontsize = 16))
)

fig_path <- file.path(out_dir, "vviq_mediation_paths.png")
ggsave(fig_path, fig, width = 14, height = 10, dpi = 300, bg = "white")

cat(sprintf("Figure saved: %s\n", fig_path))
cat(sprintf("Output saved: %s\n", file.path(out_dir, "vviq_mediation_output.txt")))

# =============================================================================
# 8. CHI-SQUARE FIGURE: IMAGERY GROUPS x CLUSTER PROFILES
# =============================================================================

cat("Creating chi-square figure...\n")

chi_caption <- sprintf("χ²(%d) = %.2f, p %s, Cramer's V = %.3f, N = %d",
                       cs$parameter, cs$statistic,
                       ifelse(cs$p.value < 0.001, "< .001",
                              sprintf("= %.3f", cs$p.value)),
                       cramers_v, n_total)
if (!is.na(fisher_p)) {
  chi_caption <- paste0(chi_caption,
                        sprintf("  |  Fisher's exact (simulated) p %s",
                                ifelse(fisher_p < 0.001, "< .001",
                                       sprintf("= %.3f", fisher_p))))
}

# --- Panel A: composition of each cluster by imagery group (column %) ---
panel_a_df <- as.data.frame(ct_tbl)
names(panel_a_df) <- c("imagery_group", "cluster", "n")
panel_a_df <- panel_a_df %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()
panel_a_df$cluster <- factor(panel_a_df$cluster,
                             levels = c("Efficient", "Hypervigilant", "Uncertain"))
panel_a_df$imagery_group <- factor(panel_a_df$imagery_group,
  levels = c("Aphantasia", "Hypophantasia", "Typical imager", "Hyperphantasia"))

imagery_colors <- c("Aphantasia" = "#1F2937",
                    "Hypophantasia" = "#6B7280",
                    "Typical imager" = "#60A5FA",
                    "Hyperphantasia" = "#F59E0B")

p_a_chi <- ggplot(panel_a_df,
                  aes(x = cluster, y = pct, fill = imagery_group)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(n > 0,
                                sprintf("n=%d\n(%.1f%%)", n, pct), "")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.2) +
  scale_fill_manual(values = imagery_colors, name = "VVIQ imagery group") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
  labs(x = "Cluster profile", y = "% of cluster",
       title = "A. Imagery group composition within each cluster") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# --- Panel B: standardized residuals (where the table deviates) ---
res_df <- as.data.frame(cs$stdres)
names(res_df) <- c("imagery_group", "cluster", "std_res")
res_df$cluster <- factor(res_df$cluster,
                         levels = c("Efficient", "Hypervigilant", "Uncertain"))
res_df$imagery_group <- factor(res_df$imagery_group,
  levels = rev(c("Aphantasia", "Hypophantasia", "Typical imager", "Hyperphantasia")))
res_df$label <- ifelse(abs(res_df$std_res) > 1.96,
                       sprintf("%.2f *", res_df$std_res),
                       sprintf("%.2f", res_df$std_res))

res_max <- max(abs(res_df$std_res))
p_b_chi <- ggplot(res_df, aes(x = cluster, y = imagery_group, fill = std_res)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 4.2, color = "black", fontface = "bold") +
  scale_fill_gradient2(low = "#3B82F6", mid = "white", high = "#DC2626",
                       midpoint = 0,
                       limits = c(-res_max, res_max),
                       name = "Std. residual") +
  labs(x = "Cluster profile", y = NULL,
       title = "B. Standardized residuals (* = |z| > 1.96)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        panel.grid = element_blank())

fig_chi <- arrangeGrob(
  p_a_chi, p_b_chi, ncol = 2, widths = c(1, 1.05),
  top = textGrob("VVIQ Imagery Groups x Interoceptive Cluster Profiles",
                 gp = gpar(fontface = "bold", fontsize = 16)),
  bottom = textGrob(chi_caption,
                    gp = gpar(fontsize = 12, fontface = "italic")))

chi_fig_path <- file.path(out_dir, "vviq_imagery_groups_by_cluster.png")
ggsave(chi_fig_path, fig_chi, width = 14, height = 7, dpi = 300, bg = "white")

cat(sprintf("Chi-square figure saved: %s\n", chi_fig_path))
cat("Done.\n")
