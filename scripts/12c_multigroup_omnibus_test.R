# =============================================================================
# MULTI-GROUP SEM: Omnibus Test for Pathway Differences Across Profiles
# =============================================================================
#
# Output: analysis_output/12c_output.txt
#
# This is the "ANOVA equivalent" for SEM - tests whether path coefficients
# differ significantly across profile groups using chi-square difference test.
#
# =============================================================================

setwd("C:/code/projects/interoceptive-profiling")

library(dplyr)
library(lavaan)

analysis_output_dir <- "C:/code/projects/interoceptive-profiling/analysis_output"
dir.create(analysis_output_dir, showWarnings = FALSE, recursive = TRUE)

# Redirect output to file
sink(file.path(analysis_output_dir, "12c_output.txt"), split = TRUE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dfc <- read.csv("C:/code/projects/interoceptive-profiling/dfc_interoception_profiling.csv")

cat("=============================================================================\n")
cat("MULTI-GROUP SEM: OMNIBUS TEST FOR PATHWAY DIFFERENCES\n")
cat("=============================================================================\n\n")

# Save full-scale totals for clustering (consistent with 04_cluster_analysis.R)
ias_full <- dfc$ias
iats_full <- dfc$iats

# Compute composite outcomes
dfc$mh <- rowMeans(cbind(scale(dfc$phq9), scale(dfc$gad7), scale(dfc$stai)), na.rm = TRUE)
dfc$somatic <- rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE)

# Standardize (full-scale IAS/IATS for SEM, consistent with 12_output.txt)
dfc$ias_z <- scale(dfc$ias)[,1]
dfc$iats_z <- scale(dfc$iats)[,1]
dfc$tas_z <- scale(dfc$tas)[,1]
dfc$somatic_z <- scale(dfc$somatic)[,1]
dfc$mh_z <- scale(dfc$mh)[,1]

# Assign profiles via clustering on FULL-SCALE totals
# (matching 04_cluster_analysis.R: seed 123, nstart 50)
set.seed(123)
km <- kmeans(scale(cbind(ias_full, iats_full)), centers = 3, nstart = 50)
dfc$cluster <- km$cluster

# Label clusters using full-scale z-scores
dfc$ias_full_z <- scale(ias_full)[,1]
dfc$iats_full_z <- scale(iats_full)[,1]
cluster_means <- aggregate(cbind(ias_full_z, iats_full_z) ~ cluster, data = dfc, mean)
cluster_means$profile <- NA
cluster_means$profile[which.max(cluster_means$ias_full_z + cluster_means$iats_full_z)] <- "Hypervigilant"
cluster_means$profile[which.min(cluster_means$ias_full_z)] <- "Uncertain"
cluster_means$profile[is.na(cluster_means$profile)] <- "Efficient"
dfc$profile <- cluster_means$profile[match(dfc$cluster, cluster_means$cluster)]
dfc$ias_full_z <- NULL
dfc$iats_full_z <- NULL

cat("Sample sizes by profile:\n")
print(table(dfc$profile))
cat("\n")

# =============================================================================
# 2. DEFINE BASE MODEL (No labels - for configural model)
# =============================================================================

# Model without labels (for configural - paths free across groups)
model_free <- '
  # Paths to mediators
  tas_z ~ ias_z + iats_z
  somatic_z ~ ias_z + iats_z

  # Paths to outcome
  mh_z ~ tas_z + somatic_z + ias_z + iats_z
'

# Model with group-specific labels (for Wald tests)
model_labeled <- '
  # Paths to mediators (group-specific labels)
  tas_z ~ c(a1_E, a1_H, a1_U)*ias_z + c(a2_E, a2_H, a2_U)*iats_z
  somatic_z ~ c(b1_E, b1_H, b1_U)*ias_z + c(b2_E, b2_H, b2_U)*iats_z

  # Paths to outcome
  mh_z ~ c(c1_E, c1_H, c1_U)*tas_z + c(c2_E, c2_H, c2_U)*somatic_z +
         c(d1_E, d1_H, d1_U)*ias_z + c(d2_E, d2_H, d2_U)*iats_z
'

# =============================================================================
# 3. FIT MULTI-GROUP MODELS
# =============================================================================

cat("=============================================================================\n")
cat("FITTING MULTI-GROUP SEM MODELS\n")
cat("=============================================================================\n\n")

# Model 1: Configural invariance (all paths free across groups)
cat("Fitting Model 1: Configural (all paths FREE across profiles)...\n")
fit_configural <- sem(model_free, data = dfc, group = "profile")

# Model 2: Constrain regression paths equal across groups
cat("Fitting Model 2: Regression paths CONSTRAINED equal...\n")
fit_equal_paths <- sem(model_free, data = dfc, group = "profile",
                       group.equal = c("regressions"))

# Model 3: Fit labeled model for Wald tests
cat("Fitting Model 3: Labeled model for Wald tests...\n")
fit_labeled <- sem(model_labeled, data = dfc, group = "profile")

cat("\n")

# =============================================================================
# 4. CHI-SQUARE DIFFERENCE TESTS (Omnibus Tests)
# =============================================================================

cat("=============================================================================\n")
cat("OMNIBUS TESTS: Chi-Square Difference (Like ANOVA F-test)\n")
cat("=============================================================================\n\n")

# Compare configural vs constrained paths
comparison1 <- anova(fit_configural, fit_equal_paths)

cat("TEST 1: Do path coefficients differ across profiles?\n")
cat("        (Configural vs. Equal Regressions)\n\n")
print(comparison1)

# Extract key values
chi_diff <- comparison1$`Chisq diff`[2]
df_diff <- comparison1$`Df diff`[2]
p_value <- comparison1$`Pr(>Chisq)`[2]

cat("\n*** OMNIBUS TEST RESULT ***\n")
cat(sprintf("Chi-square difference: %.3f\n", chi_diff))
cat(sprintf("df difference: %d\n", df_diff))
cat(sprintf("p-value: %.4f\n", p_value))

if (!is.na(p_value) && p_value < 0.05) {
  cat("\nCONCLUSION: Path coefficients SIGNIFICANTLY DIFFER across profiles (p < .05)\n")
  cat("            At least one profile has different pathway strengths.\n")
} else if (!is.na(p_value)) {
  cat("\nCONCLUSION: Path coefficients do NOT significantly differ across profiles\n")
  cat("            Pathway dominance is CONSISTENT across all profiles.\n")
} else {
  cat("\nCONCLUSION: Could not compute p-value (models may be equivalent)\n")
}

# =============================================================================
# 5. PROFILE-SPECIFIC ESTIMATES
# =============================================================================

cat("\n\n=============================================================================\n")
cat("PROFILE-SPECIFIC PATH ESTIMATES (from Configural Model)\n")
cat("=============================================================================\n\n")

# Get parameter estimates by group
params <- parameterEstimates(fit_configural, standardized = TRUE)

# Extract regression paths for each group (using op == "~" instead of labels)
groups <- c("Efficient", "Hypervigilant", "Uncertain")

for (g in 1:3) {
  group_name <- groups[g]
  cat(sprintf("--- %s (n = %d) ---\n", group_name, sum(dfc$profile == group_name)))

  group_params <- params[params$group == g & params$op == "~",
                         c("lhs", "rhs", "est", "std.all", "pvalue")]

  cat(sprintf("%-12s → %-12s %8s %8s %10s\n", "From", "To", "Est", "Std", "p"))
  for (i in 1:nrow(group_params)) {
    cat(sprintf("%-12s → %-12s %8.3f %8.3f %10.4f\n",
                group_params$rhs[i], group_params$lhs[i],
                group_params$est[i], group_params$std.all[i],
                group_params$pvalue[i]))
  }
  cat("\n")
}

# =============================================================================
# 6. INDIRECT EFFECTS BY PROFILE (Computed from path estimates)
# =============================================================================

cat("=============================================================================\n")
cat("INDIRECT EFFECTS BY PROFILE (computed from path products)\n")
cat("=============================================================================\n\n")

# Compute indirect effects manually from the configural model estimates
cat(sprintf("%-15s %15s %15s\n", "Profile", "IAS→TAS→MH", "IAS→Som→MH"))
cat(paste(rep("-", 50), collapse = ""), "\n")

for (g in 1:3) {
  group_name <- groups[g]
  gp <- params[params$group == g & params$op == "~", ]

  # Get path estimates
  a1 <- gp$std.all[gp$lhs == "tas_z" & gp$rhs == "ias_z"]
  b1 <- gp$std.all[gp$lhs == "somatic_z" & gp$rhs == "ias_z"]
  c1 <- gp$std.all[gp$lhs == "mh_z" & gp$rhs == "tas_z"]
  c2 <- gp$std.all[gp$lhs == "mh_z" & gp$rhs == "somatic_z"]

  ind_tas <- a1 * c1
  ind_som <- b1 * c2

  cat(sprintf("%-15s %15.4f %15.4f\n", group_name, ind_tas, ind_som))
}
cat("\n")

# =============================================================================
# 7. WALD TEST FOR SPECIFIC PATH EQUALITY
# =============================================================================

cat("\n=============================================================================\n")
cat("WALD TESTS: Specific Path Comparisons Across Profiles\n")
cat("=============================================================================\n\n")

# fit_labeled already fitted above with group-specific labels

# Wald test: Is IAS->Somatic (b1) equal across profiles?
cat("Test: Is IAS → Somatic path equal across profiles?\n")
wald_b1 <- lavTestWald(fit_labeled, constraints = "b1_E == b1_H; b1_H == b1_U")
cat(sprintf("  Wald chi-sq = %.3f, df = %d, p = %.4f\n\n",
            wald_b1$stat, wald_b1$df, wald_b1$p.value))

# Wald test: Is IAS->TAS (a1) equal across profiles?
cat("Test: Is IAS → TAS path equal across profiles?\n")
wald_a1 <- lavTestWald(fit_labeled, constraints = "a1_E == a1_H; a1_H == a1_U")
cat(sprintf("  Wald chi-sq = %.3f, df = %d, p = %.4f\n\n",
            wald_a1$stat, wald_a1$df, wald_a1$p.value))

# Wald test: Is TAS->MH (c1) equal across profiles?
cat("Test: Is TAS → MH path equal across profiles?\n")
wald_c1 <- lavTestWald(fit_labeled, constraints = "c1_E == c1_H; c1_H == c1_U")
cat(sprintf("  Wald chi-sq = %.3f, df = %d, p = %.4f\n\n",
            wald_c1$stat, wald_c1$df, wald_c1$p.value))

# Wald test: Is Somatic->MH (c2) equal across profiles?
cat("Test: Is Somatic → MH path equal across profiles?\n")
wald_c2 <- lavTestWald(fit_labeled, constraints = "c2_E == c2_H; c2_H == c2_U")
cat(sprintf("  Wald chi-sq = %.3f, df = %d, p = %.4f\n\n",
            wald_c2$stat, wald_c2$df, wald_c2$p.value))

# =============================================================================
# 8. SUMMARY
# =============================================================================

cat("\n=============================================================================\n")
cat("SUMMARY: OMNIBUS TEST RESULTS\n")
cat("=============================================================================\n\n")

cat("1. OVERALL CHI-SQUARE DIFFERENCE TEST:\n")
if (!is.na(p_value)) {
  cat(sprintf("   χ²(%d) = %.3f, p = %.4f\n", df_diff, chi_diff, p_value))
  if (p_value < 0.05) {
    cat("   → Paths DIFFER significantly across profiles\n\n")
  } else {
    cat("   → Paths do NOT differ significantly across profiles\n\n")
  }
} else {
  cat("   Could not compute (see comparison above)\n\n")
}

cat("2. SPECIFIC PATH TESTS (Wald):\n")
cat(sprintf("   IAS → TAS:     p = %.4f %s\n", wald_a1$p.value,
            ifelse(wald_a1$p.value < 0.05, "***", "")))
cat(sprintf("   IAS → Somatic: p = %.4f %s\n", wald_b1$p.value,
            ifelse(wald_b1$p.value < 0.05, "***", "")))
cat(sprintf("   TAS → MH:      p = %.4f %s\n", wald_c1$p.value,
            ifelse(wald_c1$p.value < 0.05, "***", "")))
cat(sprintf("   Somatic → MH:  p = %.4f %s\n", wald_c2$p.value,
            ifelse(wald_c2$p.value < 0.05, "***", "")))

cat("\n3. INTERPRETATION:\n")
all_ns <- wald_b1$p.value >= 0.05 && wald_a1$p.value >= 0.05 &&
          wald_c1$p.value >= 0.05 && wald_c2$p.value >= 0.05
if (all_ns) {
  cat("   The pathway model is INVARIANT across profiles.\n")
  cat("   Pathway dominance (TAS vs Somatic) does not significantly\n")
  cat("   depend on which profile a person belongs to.\n")
} else {
  cat("   Some paths differ across profiles. See specific tests above.\n")
}

cat("\n=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n")

# Close output file
sink()
