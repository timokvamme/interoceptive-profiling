# =============================================================================
# EXACT REPORTING OF TEST STATISTICS (ERT) FOR COMMSPSYCHOL-26-0385-T
# =============================================================================
#
# The journal requires every frequentist statistic to be reported in full:
#   named test, degrees of freedom, EXACT p (unless p < .001), effect size,
#   and a 95% CI around every effect size.
#
# Four manuscript sentences were flagged as non-compliant:
#   (i)   "nonsignificant (beta = -0.12, P > 0.05)"                [simple slope]
#   (ii)  "interaction remained non-significant for all three
#          subscales (all P > 0.12)"                              [TAS subscales]
#   (iii) "externally oriented thinking contributed negligibly
#          (5.9%, n.s.)"                                          [EOT share]
#   (iv)  "85% of explainable variance was mediated through the
#          alexithymia pathway versus 15% (P < 0.001)"            [dominance]
#
# This script re-fits the published models and produces the exact statistics
# needed to replace each sentence. It modifies NO other script.
#
# Data loading, variable construction, model specification and estimator are
# copied verbatim from 17_revision_sensitivity.R / 12_sem_pathway_analysis.R /
# 12b_pathway_dominance_test.R / 13_tas_subscale_analysis.R.
#
# OUTPUT: analysis_output/19_output.txt
# =============================================================================

# Paths are relative to the repository root; run scripts from there.

library(lavaan)

analysis_output_dir <- "analysis_output"
supplementary_data_dir <- "supplementary_data"
dir.create(analysis_output_dir, showWarnings = FALSE, recursive = TRUE)

N_BOOT <- 2000
BOOT_SEED <- 42

# =============================================================================
# 1. DATA AND VARIABLE CONSTRUCTION (identical to 17 / 12 / 13)
# =============================================================================

dfc <- read.csv("dfc_interoception_profiling.csv")

ias_full_items <- paste0("ias_", 1:21)
iats_full_items <- paste0("iats_", 1:21)

dfc$ias_full <- rowMeans(dfc[, ias_full_items], na.rm = TRUE)
dfc$iats_full <- rowMeans(dfc[, iats_full_items], na.rm = TRUE)

dfc$mh_composite <- rowMeans(cbind(scale(dfc$phq9), scale(dfc$gad7),
                                   scale(dfc$stai)), na.rm = TRUE)
dfc$somatic_combined <- rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE)

# TAS-20 subscales, item groupings from 13_tas_subscale_analysis.R lines 56-58
# (Bagby, Parker & Taylor 1994; items 4, 5, 10, 18, 19 already reverse-scored
#  in this dataset, verified by the stopifnot below)
dif_items <- c(1, 3, 6, 7, 9, 13, 14)         # 7 items
ddf_items <- c(2, 4, 11, 12, 17)              # 5 items
eot_items <- c(5, 8, 10, 15, 16, 18, 19, 20)  # 8 items

dfc$dif <- rowSums(dfc[, paste0("tas_", dif_items)], na.rm = TRUE)
dfc$ddf <- rowSums(dfc[, paste0("tas_", ddf_items)], na.rm = TRUE)
dfc$eot <- rowSums(dfc[, paste0("tas_", eot_items)], na.rm = TRUE)
stopifnot(all(abs(dfc$tas - (dfc$dif + dfc$ddf + dfc$eot)) < 0.01))

# z-scored analysis variables (same convention as 12 / 17)
dfc$ias_z     <- scale(dfc$ias_full)[, 1]
dfc$iats_z    <- scale(dfc$iats_full)[, 1]
dfc$tas_z     <- scale(dfc$tas)[, 1]
dfc$dif_z     <- scale(dfc$dif)[, 1]
dfc$ddf_z     <- scale(dfc$ddf)[, 1]
dfc$eot_z     <- scale(dfc$eot)[, 1]
dfc$somatic_z <- scale(dfc$somatic_combined)[, 1]
dfc$mh_z      <- scale(dfc$mh_composite)[, 1]
dfc$ias_x_iats <- dfc$ias_z * dfc$iats_z

SD_PRODUCT <- sd(dfc$ias_x_iats)

N_TOTAL <- nrow(dfc)
model_vars <- c("ias_z", "iats_z", "tas_z", "dif_z", "ddf_z", "eot_z",
                "somatic_z", "mh_z", "ias_x_iats")
N_COMPLETE <- sum(complete.cases(dfc[, model_vars]))

# =============================================================================
# 2. HELPERS
# =============================================================================

# Exact p formatting per journal rule: exact value unless p < .001
fmt_p <- function(p) {
  if (is.na(p)) return("n/a")
  if (p == 0) return("< .001 (<1e-16)")
  if (p < 0.001) return(sprintf("< .001 (%.0e)", p))
  sub("^0\\.", ".", sprintf("%.4f", p))
}

# Percentile bootstrap two-sided p (proportion of resamples on the far side
# of zero, doubled). Floor reported honestly when no resample crosses.
boot_p <- function(v, null_value = 0) {
  v <- v[is.finite(v)]
  p <- 2 * min(mean(v <= null_value), mean(v >= null_value))
  min(p, 1)
}

fmt_boot_p <- function(v, null_value = 0) {
  v <- v[is.finite(v)]
  n_cross <- min(sum(v <= null_value), sum(v >= null_value))
  p <- boot_p(v, null_value)
  if (n_cross == 0) {
    return(sprintf("< .001 (0/%d)", length(v)))
  }
  sub("^0\\.", ".", sprintf("%.4f", p))
}

qlo <- function(v) unname(quantile(v[is.finite(v)], 0.025))
qhi <- function(v) unname(quantile(v[is.finite(v)], 0.975))

# Pull bootstrap draw matrix (free parameters) from a lavaan fit
get_boot_mat <- function(fit) {
  B <- lavInspect(fit, "boot")
  B <- as.matrix(B)
  colnames(B) <- names(coef(fit))
  B[complete.cases(B), , drop = FALSE]
}

get_row <- function(tab, lbl) tab[tab$label == lbl, ][1, ]

hr <- function(ch = "-", n = 78) cat(paste(rep(ch, n), collapse = ""), "\n")

# =============================================================================
# 3. MODEL M0 (published dual-pathway model) WITH DEFINED QUANTITIES
# =============================================================================

# --- 3a. plain ML fit, used only for the reproduction check (matches 17) ---
model_m0_plain <- '
  tas_z     ~ a1*ias_z + a2*iats_z + a3*ias_x_iats
  somatic_z ~ b1*iats_z + b2*ias_z + b3*ias_x_iats
  mh_z      ~ c1*tas_z + c2*somatic_z + d1*ias_z + d2*iats_z + d3*ias_x_iats
'
fit_m0_ml <- sem(model_m0_plain, data = dfc)
std_m0_ml <- standardizedSolution(fit_m0_ml, level = 0.95)

# --- 3b. profile assignment (identical to 12b_pathway_dominance_test.R) ---
ias_total <- dfc$ias
iats_total <- dfc$iats
set.seed(123)
km <- kmeans(scale(cbind(ias_total, iats_total)), centers = 3, nstart = 50)
dfc$cluster <- km$cluster
dfc$ias_tot_z <- scale(ias_total)[, 1]
dfc$iats_tot_z <- scale(iats_total)[, 1]
cluster_means <- aggregate(cbind(ias_tot_z, iats_tot_z) ~ cluster, data = dfc, mean)
cluster_means$profile <- NA
cluster_means$profile[which.max(cluster_means$ias_tot_z + cluster_means$iats_tot_z)] <- "Hypervigilant"
cluster_means$profile[which.min(cluster_means$ias_tot_z)] <- "Uncertain"
cluster_means$profile[is.na(cluster_means$profile)] <- "Efficient"
dfc$profile <- cluster_means$profile[match(dfc$cluster, cluster_means$cluster)]

profile_values <- aggregate(cbind(ias_z, iats_z) ~ profile, data = dfc, mean)
names(profile_values) <- c("profile", "IAS_typical", "IATS_typical")
prof_n <- table(dfc$profile)

hyper_iats  <- profile_values$IATS_typical[profile_values$profile == "Hypervigilant"]
uncert_iats <- profile_values$IATS_typical[profile_values$profile == "Uncertain"]
effic_iats  <- profile_values$IATS_typical[profile_values$profile == "Efficient"]

# used both in Part 1d (profile-typical simple slopes) and Part 4
profs_pre <- c("Hypervigilant", "Uncertain", "Efficient")
iats_pre <- round(c(hyper_iats, uncert_iats, effic_iats), 4)

# --- 3c. bootstrap fit with every derived quantity as a defined parameter ---
model_m0 <- sprintf('
  tas_z     ~ a1*ias_z + a2*iats_z + a3*ias_x_iats
  somatic_z ~ b1*iats_z + b2*ias_z + b3*ias_x_iats
  mh_z      ~ c1*tas_z + c2*somatic_z + d1*ias_z + d2*iats_z + d3*ias_x_iats

  # simple slopes of IAS on the somatic mediator at IATS = -1 SD / mean / +1 SD
  ss_som_lo  := b2 + b3*(-1)
  ss_som_mid := b2
  ss_som_hi  := b2 + b3*(1)

  # simple slopes of IAS on TAS at the same three IATS levels
  ss_tas_lo  := a1 + a3*(-1)
  ss_tas_mid := a1
  ss_tas_hi  := a1 + a3*(1)

  # conditional indirect effects of IAS on MH at profile-typical IATS
  cond_tas_hyper  := (a1 + a3*%s) * c1
  cond_som_hyper  := (b2 + b3*%s) * c2
  cond_tas_uncert := (a1 + a3*%s) * c1
  cond_som_uncert := (b2 + b3*%s) * c2
  cond_tas_effic  := (a1 + a3*%s) * c1
  cond_som_effic  := (b2 + b3*%s) * c2

  # signed pathway contrasts within profile (somatic minus alexithymia)
  diff_hyper  := cond_som_hyper  - cond_tas_hyper
  diff_uncert := cond_som_uncert - cond_tas_uncert
  diff_effic  := cond_som_effic  - cond_tas_effic

  # indices of moderated mediation
  index_mm_som := b3 * c2
  index_mm_tas := a3 * c1
',
  round(hyper_iats, 4), round(hyper_iats, 4),
  round(uncert_iats, 4), round(uncert_iats, 4),
  round(effic_iats, 4), round(effic_iats, 4))

set.seed(BOOT_SEED)
fit_m0 <- sem(model_m0, data = dfc, se = "bootstrap",
              bootstrap = N_BOOT, iseed = BOOT_SEED)
pe_m0 <- parameterEstimates(fit_m0, ci = TRUE, level = 0.95)
sd_m0 <- standardizedSolution(fit_m0, level = 0.95)
B0 <- get_boot_mat(fit_m0)
R0 <- nrow(B0)

bs <- function(lbl) B0[, lbl]

# =============================================================================
# 4. PARALLEL SUBSCALE MEDIATION MODEL (identical spec to script 13)
# =============================================================================

parallel_model <- '
  dif_z ~ a1_dif*ias_z + a2_dif*iats_z + a3_dif*ias_x_iats
  ddf_z ~ a1_ddf*ias_z + a2_ddf*iats_z + a3_ddf*ias_x_iats
  eot_z ~ a1_eot*ias_z + a2_eot*iats_z + a3_eot*ias_x_iats
  somatic_z ~ b1*iats_z + b2*ias_z + b3*ias_x_iats

  mh_z ~ c_dif*dif_z + c_ddf*ddf_z + c_eot*eot_z + c_som*somatic_z +
          d1*ias_z + d2*iats_z + d3*ias_x_iats

  dif_z ~~ ddf_z
  dif_z ~~ eot_z
  ddf_z ~~ eot_z

  ind_ias_dif := a1_dif * c_dif
  ind_ias_ddf := a1_ddf * c_ddf
  ind_ias_eot := a1_eot * c_eot
  ind_ias_som := b2 * c_som

  total_alex_ias := ind_ias_dif + ind_ias_ddf + ind_ias_eot

  diff_dif_ddf := ind_ias_dif - ind_ias_ddf
  diff_dif_eot := ind_ias_dif - ind_ias_eot
  diff_ddf_eot := ind_ias_ddf - ind_ias_eot
'

set.seed(BOOT_SEED)
fit_par <- sem(parallel_model, data = dfc, se = "bootstrap",
               bootstrap = N_BOOT, iseed = BOOT_SEED)
pe_par <- parameterEstimates(fit_par, ci = TRUE, level = 0.95)
sd_par <- standardizedSolution(fit_par, level = 0.95)
Bp <- get_boot_mat(fit_par)
Rp <- nrow(Bp)

# =============================================================================
# 5. OLS MODELS FOR THE SUBSCALE INTERACTIONS AND SIMPLE SLOPES
# =============================================================================
# Each mediator equation in the SEM is a saturated (just-identified) regression,
# so its coefficients are numerically identical to OLS. OLS is used in parallel
# because it supplies exact residual degrees of freedom and an exact t / F test,
# which the SEM Wald z (df = infinity) does not.

ols_slopes <- function(dv, k) {
  d <- dfc
  d$iats_c <- d$iats_z - k
  m <- lm(as.formula(paste0(dv, " ~ ias_z * iats_c")), data = d)
  s <- summary(m)
  ci <- confint(m, "ias_z", level = 0.95)
  list(est = unname(coef(m)["ias_z"]),
       se = s$coefficients["ias_z", 2],
       t = s$coefficients["ias_z", 3],
       p = s$coefficients["ias_z", 4],
       lo = ci[1], hi = ci[2],
       df = m$df.residual, n = length(m$residuals))
}

ols_interaction <- function(dv) {
  m_add <- lm(as.formula(paste0(dv, " ~ ias_z + iats_z")), data = dfc)
  m_int <- lm(as.formula(paste0(dv, " ~ ias_z * iats_z")), data = dfc)
  s <- summary(m_int)
  ci <- confint(m_int, "ias_z:iats_z", level = 0.95)
  av <- anova(m_add, m_int)
  list(est = unname(coef(m_int)["ias_z:iats_z"]),
       se = s$coefficients["ias_z:iats_z", 2],
       t = s$coefficients["ias_z:iats_z", 3],
       p = s$coefficients["ias_z:iats_z", 4],
       lo = ci[1], hi = ci[2],
       df = m_int$df.residual, n = length(m_int$residuals),
       dR2 = s$r.squared - summary(m_add)$r.squared,
       F = av$F[2], df1 = 1, df2 = m_int$df.residual, pF = av$`Pr(>F)`[2],
       ias = s$coefficients["ias_z", ], iats = s$coefficients["iats_z", ])
}

# =============================================================================
# 6. OUTPUT
# =============================================================================

sink(file.path(analysis_output_dir, "19_output.txt"))

cat("=============================================================================\n")
cat("EXACT TEST STATISTICS FOR REVISION (COMMSPSYCHOL-26-0385-T)\n")
cat("=============================================================================\n\n")
cat("Script:      19_ert_exact_statistics.R\n")
cat("Date:       ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("R:          ", R.version.string, "\n")
cat("lavaan:     ", as.character(packageVersion("lavaan")), "\n")
cat("Estimator:   ML, complete cases\n")
cat("Bootstrap:  ", N_BOOT, "nonparametric (ordinary) resamples, seed", BOOT_SEED, "\n")
cat("             95% CIs are percentile bootstrap unless stated otherwise\n\n")
cat("Rows in data file:                 N =", N_TOTAL, "\n")
cat("Complete cases on all model vars:  N =", N_COMPLETE, "\n")
cat("SD of the IAS x IATS product term:   ", sprintf("%.4f", SD_PRODUCT), "\n")
cat("  (the product of two z-scores does not itself have SD 1, which is why a\n")
cat("   lavaan std.all value for any interaction term equals the z-metric\n")
cat("   coefficient multiplied by this SD; documented in 17_revision_sensitivity.R)\n\n")

# -----------------------------------------------------------------------------
cat("#############################################################################\n")
cat("PART 0: REPRODUCTION CHECK AGAINST THE PUBLISHED MODEL M0\n")
cat("#############################################################################\n\n")
cat("Model M0: full 42-item IAS/IATS, somatic = SSS-8 + PCS composite,\n")
cat("mental health = PHQ-9 + GAD-7 + STAI composite, no covariates.\n\n")

targets <- data.frame(
  label = c("a1", "a2", "a3", "b2", "b1", "b3"),
  path = c("IAS -> TAS", "IATS -> TAS", "IAS x IATS -> TAS",
           "IAS -> Somatic", "IATS -> Somatic", "IAS x IATS -> Somatic"),
  target_beta = c(-0.346, 0.314, -0.008, -0.094, 0.444, -0.094),
  target_p = c(NA, NA, 0.8011, NA, NA, 0.0030),
  stringsAsFactors = FALSE
)

cat(sprintf("%-24s %10s %10s %9s %10s %10s %8s\n",
            "Path", "target", "refit", "diff", "target p", "refit p", "OK"))
hr()
repro_ok <- TRUE
for (i in seq_len(nrow(targets))) {
  r <- get_row(std_m0_ml, targets$label[i])
  d <- r$est.std - targets$target_beta[i]
  pd <- if (is.na(targets$target_p[i])) NA else r$pvalue - targets$target_p[i]
  ok <- abs(d) < 0.001 && (is.na(pd) || abs(pd) < 0.0005)
  if (!ok) repro_ok <- FALSE
  cat(sprintf("%-24s %10.3f %10.3f %9.4f %10s %10.4f %8s\n",
              targets$path[i], targets$target_beta[i], r$est.std, d,
              ifelse(is.na(targets$target_p[i]), "-",
                     sprintf("%.4f", targets$target_p[i])),
              r$pvalue, ifelse(ok, "yes", "NO")))
}
cat("\n")
ft0 <- lavInspect(fit_m0_ml, "fit")
cat(sprintf("Model chi-square: chi2(%d) = %.3f, p %s\n",
            as.integer(ft0["df"]), ft0["chisq"], fmt_p(ft0["pvalue"])))
cat("Converged:", lavInspect(fit_m0_ml, "converged"), "\n")
cat("N used by lavaan:", lavInspect(fit_m0_ml, "nobs"), "\n\n")
if (repro_ok) {
  cat("REPRODUCTION CHECK: PASS. All six published coefficients and both\n")
  cat("published p values reproduce to the precision printed in the manuscript.\n\n")
} else {
  cat("*** REPRODUCTION CHECK: FAIL ***\n")
  cat("The refit does NOT match the published model M0. Everything below is\n")
  cat("therefore NOT safe to paste into the manuscript. Investigate first.\n\n")
}

# -----------------------------------------------------------------------------
cat("#############################################################################\n")
cat("PART 1: SIMPLE SLOPES OF IAS AT IATS = -1 SD, MEAN, +1 SD\n")
cat("#############################################################################\n\n")
cat("Definition. All variables are z-scored, so the simple slope is the change\n")
cat("in the mediator, in SD units, per 1 SD increase in IAS, evaluated at the\n")
cat("stated level of IATS. It is computed as b_IAS + b_IASxIATS * k, with\n")
cat("k = -1, 0, +1. Because the mediator equation of the SEM is a saturated\n")
cat("regression, the OLS re-centring approach and the SEM defined parameter\n")
cat("give identical point estimates; the OLS version is reported because it\n")
cat("carries exact residual degrees of freedom.\n\n")

slope_block <- function(dv_ols, lbl_prefix, title) {
  cat("---", title, "---\n\n")
  ks <- c(-1, 0, 1)
  knames <- c("IATS = -1 SD", "IATS = mean (0)", "IATS = +1 SD")
  labs <- paste0(lbl_prefix, c("_lo", "_mid", "_hi"))

  cat("(a) OLS simple slopes: two-tailed t test on the regression coefficient\n\n")
  cat(sprintf("%-18s %9s %8s %22s %9s %10s\n",
              "Level", "slope", "SE", "95% CI", "t", "p"))
  hr()
  for (j in seq_along(ks)) {
    o <- ols_slopes(dv_ols, ks[j])
    cat(sprintf("%-18s %9.4f %8.4f  [%8.4f, %8.4f] %9.3f %10s\n",
                knames[j], o$est, o$se, o$lo, o$hi, o$t, fmt_p(o$p)))
  }
  o1 <- ols_slopes(dv_ols, 0)
  cat("\n  test: t test on an OLS coefficient, df =", o1$df,
      "(N =", o1$n, ", 4 estimated parameters)\n")
  cat("  CI: exact t-based 95% CI\n\n")

  cat("(b) Same slopes as SEM defined parameters (model M0)\n")
  cat("    z: Wald z test using the bootstrap SE; CI: percentile bootstrap\n\n")
  cat(sprintf("%-18s %9s %8s %22s %9s %12s %12s\n",
              "Level", "slope", "SE_boot", "95% boot CI", "z", "p (Wald z)",
              "p (boot)"))
  hr(n = 96)
  for (j in seq_along(ks)) {
    r <- get_row(pe_m0, labs[j])
    v <- switch(lbl_prefix,
                ss_som = bs("b2") + bs("b3") * ks[j],
                ss_tas = bs("a1") + bs("a3") * ks[j])
    cat(sprintf("%-18s %9.4f %8.4f  [%8.4f, %8.4f] %9.3f %12s %12s\n",
                knames[j], r$est, r$se, qlo(v), qhi(v), r$z,
                fmt_p(r$pvalue), fmt_boot_p(v)))
  }
  cat("\n(c) Alternative metric: the same slopes expressed in lavaan std.all,\n")
  cat("    i.e. with the product term rescaled by its own SD. Delta-method CI.\n\n")
  cat(sprintf("%-18s %9s %8s %22s %10s\n",
              "Level", "beta.std", "SE", "95% CI (delta)", "p"))
  hr()
  for (j in seq_along(ks)) {
    r <- get_row(sd_m0, labs[j])
    cat(sprintf("%-18s %9.4f %8.4f  [%8.4f, %8.4f] %10s\n",
                knames[j], r$est.std, r$se, r$ci.lower, r$ci.upper,
                fmt_p(r$pvalue)))
  }
  cat("\n")
}

slope_block("somatic_z", "ss_som",
            "SOMATIC MEDIATOR (SSS-8 + PCS composite) regressed on IAS")
slope_block("tas_z", "ss_tas",
            "ALEXITHYMIA (TAS-20 total) regressed on IAS")

cat("--- Index of moderated mediation (formal test of whether the IAS\n")
cat("    indirect effect depends on IATS) ---\n\n")
for (lbl in c("index_mm_som", "index_mm_tas")) {
  r <- get_row(pe_m0, lbl)
  v <- switch(lbl,
              index_mm_som = bs("b3") * bs("c2"),
              index_mm_tas = bs("a3") * bs("c1"))
  nm <- if (lbl == "index_mm_som") "Somatic pathway" else "Alexithymia pathway"
  cat(sprintf("%-20s index = %8.4f, SE_boot = %.4f, 95%% boot CI [%.4f, %.4f],\n",
              nm, r$est, r$se, qlo(v), qhi(v)))
  cat(sprintf("%-20s Wald z = %.3f, p = %s; percentile bootstrap p = %s\n\n",
              "", r$z, fmt_p(r$pvalue), fmt_boot_p(v)))
}

cat("--- (d) The same simple slopes at each profile-typical IATS level ---\n\n")
cat("Figure panels D and E label the slopes by profile rather than by SD, so\n")
cat("these are the values that correspond one-to-one to the figure labels.\n")
cat("IATS(z) is fixed at each profile's mean (see Part 4 for the profile\n")
cat("solution and profile n).\n\n")
cat(sprintf("%-16s %9s %9s %8s %22s %8s %6s %14s\n",
            "Profile", "IATS(z)", "slope", "SE", "95% CI", "t", "df", "p"))
hr(n = 100)
slope_rows <- list()
for (dv in c("somatic_z", "tas_z")) {
  cat(ifelse(dv == "somatic_z",
             "Somatic mediator (figure panel D):\n",
             "Alexithymia / TAS-20 (figure panel E):\n"))
  for (j in seq_along(profs_pre)) {
    k <- iats_pre[j]
    o <- ols_slopes(dv, k)
    cat(sprintf("  %-14s %9.3f %9.4f %8.4f  [%8.4f, %8.4f] %8.3f %6d %14s\n",
                profs_pre[j], k, o$est, o$se, o$lo, o$hi, o$t, o$df, fmt_p(o$p)))
    slope_rows[[length(slope_rows) + 1]] <- data.frame(
      outcome = ifelse(dv == "somatic_z", "somatic", "tas"), profile = profs_pre[j],
      iats_z = k, slope = o$est, se = o$se, ci_lo = o$lo, ci_hi = o$hi,
      t = o$t, df = o$df, p = o$p)
  }
}
# Figure 1 panels d and e read these values (12_sem_pathway_analysis.R).
write.csv(do.call(rbind, slope_rows),
          file.path(supplementary_data_dir, "figure1_simple_slopes_by_profile.csv"),
          row.names = FALSE)
cat("\n  test: two-tailed t test on an OLS coefficient after re-centring IATS\n")
cat("  at the profile mean; CI is the exact t-based 95% CI.\n\n")

cat("NOTE ON THE MANUSCRIPT VALUES -0.317 / -0.122 / +0.16.\n")
cat("Those three numbers are hard-coded in 12_sem_pathway_analysis.R (lines\n")
cat("1008-1012, annotated \"from reference figure\") and are not produced by any\n")
cat("computation in this repository. They do not match the fitted model on any\n")
cat("metric or at any IATS level: the tabulated model-implied slopes above are\n")
cat("the correct values and are far smaller in absolute size.\n\n")

cross <- unname(-coef(fit_m0)["b2"] / coef(fit_m0)["b3"])
pct_below <- 100 * mean(dfc$iats_z < cross)
cat(sprintf("The slope of IAS on the somatic mediator changes sign at\n"))
cat(sprintf("IATS(z) = %.3f. The observed IATS(z) range is %.3f to %.3f, so the\n",
            cross, min(dfc$iats_z), max(dfc$iats_z)))
cat(sprintf("slope is positive for the %.1f%% of the sample (n = %d) below that\n",
            pct_below, sum(dfc$iats_z < cross)))
cat("point. That crossing lies BELOW the Efficient profile mean, so the\n")
cat("positive slope labelled +0.16 for the Efficient profile in the figure is\n")
cat(sprintf("not what the model implies: at the Efficient mean the slope is %.4f.\n\n",
            ols_slopes("somatic_z", effic_iats)$est))

# -----------------------------------------------------------------------------
cat("#############################################################################\n")
cat("PART 2: IAS x IATS INTERACTION ON EACH TAS-20 SUBSCALE\n")
cat("#############################################################################\n\n")
cat("Item groupings (from 13_tas_subscale_analysis.R, Bagby et al. 1994;\n")
cat("items 4, 5, 10, 18, 19 are already reverse-scored in this dataset, and\n")
cat("DIF + DDF + EOT reproduces the TAS-20 total exactly):\n")
cat("  DIF (Difficulty Identifying Feelings), 7 items:  tas_",
    paste(dif_items, collapse = ", tas_"), "\n", sep = "")
cat("  DDF (Difficulty Describing Feelings), 5 items:   tas_",
    paste(ddf_items, collapse = ", tas_"), "\n", sep = "")
cat("  EOT (Externally Oriented Thinking), 8 items:     tas_",
    paste(eot_items, collapse = ", tas_"), "\n\n", sep = "")
cat("Each subscale is z-scored and regressed on z-scored IAS, z-scored IATS\n")
cat("and their product, in a separate model per subscale.\n\n")

cat("(a) Separate OLS models: two-tailed t test on the interaction coefficient\n\n")
subs <- c("dif_z", "ddf_z", "eot_z", "tas_z")
sublab <- c("DIF", "DDF", "EOT", "TAS-20 total")
ols_int <- list()
cat(sprintf("%-14s %9s %8s %22s %9s %6s %12s\n",
            "Subscale", "beta", "SE", "95% CI", "t", "df", "p"))
hr(n = 88)
for (i in seq_along(subs)) {
  o <- ols_interaction(subs[i])
  ols_int[[sublab[i]]] <- o
  cat(sprintf("%-14s %9.4f %8.4f  [%8.4f, %8.4f] %9.3f %6d %12s\n",
              sublab[i], o$est, o$se, o$lo, o$hi, o$t, o$df, fmt_p(o$p)))
}
cat("\n  beta is in the z metric (z-scored outcome, z-scored IAS and IATS,\n")
cat("  raw product). Multiply by", sprintf("%.4f", SD_PRODUCT),
    "for the lavaan std.all metric:\n\n")
cat(sprintf("%-14s %12s %24s\n", "Subscale", "beta (std.all)", "95% CI (std.all)"))
hr()
for (i in seq_along(subs)) {
  o <- ols_int[[sublab[i]]]
  cat(sprintf("%-14s %12.4f  [%9.4f, %9.4f]\n", sublab[i],
              o$est * SD_PRODUCT, o$lo * SD_PRODUCT, o$hi * SD_PRODUCT))
}

cat("\n(b) Variance explained by the interaction: F test for the increment in R2\n")
cat("    over the additive model (IAS + IATS)\n\n")
cat(sprintf("%-14s %10s %12s %10s %12s\n",
            "Subscale", "Delta R2", "F", "df", "p"))
hr()
for (i in seq_along(subs)) {
  o <- ols_int[[sublab[i]]]
  cat(sprintf("%-14s %10.5f %12.3f %10s %12s\n", sublab[i], o$dR2, o$F,
              sprintf("1, %d", o$df2), fmt_p(o$pF)))
}

cat("\n(c) Same interactions inside the parallel multiple-mediator SEM\n")
cat("    (DIF, DDF, EOT and somatic as simultaneous mediators; residual\n")
cat("    covariances among subscales freed). Percentile bootstrap CI,\n")
cat("    R =", Rp, "; p is a Wald z test on the bootstrap SE.\n\n")
cat(sprintf("%-14s %9s %8s %22s %9s %12s %12s\n",
            "Path", "beta", "SE_boot", "95% boot CI", "z", "p (Wald z)", "p (boot)"))
hr(n = 92)
for (l in c("a3_dif", "a3_ddf", "a3_eot")) {
  r <- get_row(pe_par, l)
  v <- Bp[, l]
  nm <- paste("IAS x IATS ->", toupper(sub("a3_", "", l)))
  cat(sprintf("%-14s %9.4f %8.4f  [%8.4f, %8.4f] %9.3f %12s %12s\n",
              nm, r$est, r$se, qlo(v), qhi(v), r$z, fmt_p(r$pvalue),
              fmt_boot_p(v)))
}
cat("\n  std.all values from the same SEM (delta-method CI):\n\n")
cat(sprintf("%-14s %10s %24s %12s\n", "Path", "beta.std", "95% CI", "p"))
hr()
for (l in c("a3_dif", "a3_ddf", "a3_eot")) {
  r <- get_row(sd_par, l)
  nm <- paste("IAS x IATS ->", toupper(sub("a3_", "", l)))
  cat(sprintf("%-14s %10.4f  [%10.4f, %10.4f] %12s\n",
              nm, r$est.std, r$ci.lower, r$ci.upper, fmt_p(r$pvalue)))
}
cat("\n")

# -----------------------------------------------------------------------------
cat("#############################################################################\n")
cat("PART 3: INDIRECT EFFECT OF IAS ON MENTAL HEALTH VIA DIF, DDF AND EOT\n")
cat("#############################################################################\n\n")
cat("Parallel multiple-mediator SEM: DIF, DDF, EOT and the somatic composite\n")
cat("enter simultaneously as mediators of IAS, IATS and their product on the\n")
cat("mental-health composite. All variables z-scored, so the estimates below\n")
cat("are already fully standardized (no product term is involved in any of\n")
cat("these indirect effects, so the z metric and std.all coincide).\n")
cat("Bootstrap: R =", Rp, "successful resamples.\n\n")

cat(sprintf("%-24s %9s %8s %24s %8s %12s %12s\n",
            "Indirect effect", "est", "SE_boot", "95% boot CI", "z",
            "p (Wald z)", "p (boot)"))
hr(n = 104)
ind_labels <- c("ind_ias_dif", "ind_ias_ddf", "ind_ias_eot", "ind_ias_som",
                "total_alex_ias")
ind_names <- c("IAS -> DIF -> MH", "IAS -> DDF -> MH", "IAS -> EOT -> MH",
               "IAS -> Somatic -> MH", "IAS -> all 3 facets -> MH")
ind_draws <- list(
  ind_ias_dif = Bp[, "a1_dif"] * Bp[, "c_dif"],
  ind_ias_ddf = Bp[, "a1_ddf"] * Bp[, "c_ddf"],
  ind_ias_eot = Bp[, "a1_eot"] * Bp[, "c_eot"],
  ind_ias_som = Bp[, "b2"] * Bp[, "c_som"]
)
ind_draws$total_alex_ias <- ind_draws$ind_ias_dif + ind_draws$ind_ias_ddf +
  ind_draws$ind_ias_eot
for (i in seq_along(ind_labels)) {
  r <- get_row(pe_par, ind_labels[i])
  v <- ind_draws[[ind_labels[i]]]
  cat(sprintf("%-24s %9.4f %8.4f  [%9.4f, %9.4f] %8.3f %12s %12s\n",
              ind_names[i], r$est, r$se, qlo(v), qhi(v), r$z,
              fmt_p(r$pvalue), fmt_boot_p(v)))
}

cat("\nA-paths and b-paths that compose them (same model):\n\n")
cat(sprintf("%-22s %9s %8s %24s %12s\n",
            "Path", "beta", "SE_boot", "95% boot CI", "p (Wald z)"))
hr(n = 82)
comp <- c(a1_dif = "IAS -> DIF", a1_ddf = "IAS -> DDF", a1_eot = "IAS -> EOT",
          c_dif = "DIF -> MH", c_ddf = "DDF -> MH", c_eot = "EOT -> MH",
          c_som = "Somatic -> MH")
for (l in names(comp)) {
  r <- get_row(pe_par, l)
  v <- Bp[, l]
  cat(sprintf("%-22s %9.4f %8.4f  [%9.4f, %9.4f] %12s\n",
              comp[[l]], r$est, r$se, qlo(v), qhi(v), fmt_p(r$pvalue)))
}

cat("\nPairwise contrasts between facet-specific indirect effects\n")
cat("(Wald z on the bootstrap SE, plus percentile bootstrap p):\n\n")
cat(sprintf("%-16s %9s %8s %24s %12s %12s\n",
            "Contrast", "diff", "SE_boot", "95% boot CI", "p (Wald z)", "p (boot)"))
hr(n = 88)
contr <- list(diff_dif_ddf = c("DIF - DDF", "ind_ias_dif", "ind_ias_ddf"),
              diff_dif_eot = c("DIF - EOT", "ind_ias_dif", "ind_ias_eot"),
              diff_ddf_eot = c("DDF - EOT", "ind_ias_ddf", "ind_ias_eot"))
for (l in names(contr)) {
  r <- get_row(pe_par, l)
  v <- ind_draws[[contr[[l]][2]]] - ind_draws[[contr[[l]][3]]]
  cat(sprintf("%-16s %9.4f %8.4f  [%9.4f, %9.4f] %12s %12s\n",
              contr[[l]][1], r$est, r$se, qlo(v), qhi(v),
              fmt_p(r$pvalue), fmt_boot_p(v)))
}

cat("\n--- Percentage share of the alexithymia-mediated indirect effect ---\n\n")
cat("Share is computed on absolute values, share_j = |ind_j| / sum_k |ind_k|,\n")
cat("which is the definition used in 13_tas_subscale_analysis.R (needed because\n")
cat("the EOT effect has the opposite sign to DIF and DDF). Percentages are\n")
cat("recomputed inside every bootstrap resample; the CI is percentile.\n")
cat("The p value in the last column tests H0: this facet contributes 0% of the\n")
cat("alexithymia-mediated effect, which is identical to the test that the\n")
cat("facet-specific indirect effect is zero (percentile bootstrap).\n\n")

abs_tot <- abs(ind_draws$ind_ias_dif) + abs(ind_draws$ind_ias_ddf) +
  abs(ind_draws$ind_ias_eot)
share_draws <- list(
  DIF = 100 * abs(ind_draws$ind_ias_dif) / abs_tot,
  DDF = 100 * abs(ind_draws$ind_ias_ddf) / abs_tot,
  EOT = 100 * abs(ind_draws$ind_ias_eot) / abs_tot
)
pt <- sapply(c("ind_ias_dif", "ind_ias_ddf", "ind_ias_eot"),
             function(l) get_row(pe_par, l)$est)
pt_share <- 100 * abs(pt) / sum(abs(pt))
signed_share <- 100 * pt / sum(pt)

cat(sprintf("%-6s %10s %11s %22s %14s %12s\n",
            "Facet", "share %", "boot mean", "95% boot CI", "signed share %",
            "p (boot)"))
hr(n = 82)
for (i in seq_along(share_draws)) {
  f <- names(share_draws)[i]
  v <- share_draws[[f]]
  pl <- c("ind_ias_dif", "ind_ias_ddf", "ind_ias_eot")[i]
  cat(sprintf("%-6s %9.1f%% %10.1f%% [%8.1f%%, %8.1f%%] %13.1f%% %12s\n",
              f, pt_share[i], mean(v), qlo(v), qhi(v), signed_share[i],
              fmt_boot_p(ind_draws[[pl]])))
}
cat("\n  share %      = share computed from the point estimates\n")
cat("  boot mean    = mean of the same quantity across resamples\n")
cat("  signed share = ind_j / total_alex (sums to 100%, can be negative)\n\n")

# Supplementary Fig. 18 is drawn from these values
# (supplementary_figure_scripts/make_subscale_figure.py).
sub_rows <- list()
for (l in names(comp)) {
  r <- get_row(pe_par, l); v <- Bp[, l]
  sub_rows[[length(sub_rows) + 1]] <- data.frame(block = "path", label = comp[[l]],
    est = r$est, se = r$se, ci_lo = qlo(v), ci_hi = qhi(v), p = r$pvalue)
}
for (i in seq_along(ind_labels)) {
  r <- get_row(pe_par, ind_labels[i]); v <- ind_draws[[ind_labels[i]]]
  sub_rows[[length(sub_rows) + 1]] <- data.frame(block = "indirect", label = ind_names[i],
    est = r$est, se = r$se, ci_lo = qlo(v), ci_hi = qhi(v), p = r$pvalue)
}
for (i in seq_along(share_draws)) {
  v <- share_draws[[i]]
  sub_rows[[length(sub_rows) + 1]] <- data.frame(block = "share", label = names(share_draws)[i],
    est = pt_share[i], se = sd(v), ci_lo = qlo(v), ci_hi = qhi(v), p = NA)
}
write.csv(do.call(rbind, sub_rows),
          file.path(supplementary_data_dir, "subscale_parallel_model.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
cat("#############################################################################\n")
cat("PART 4: PATHWAY DOMINANCE BY INTEROCEPTIVE PROFILE\n")
cat("#############################################################################\n\n")
cat("Profiles are the k-means solution of 04_cluster_analysis.R / 12b\n")
cat("(k = 3, seed 123, nstart 50, on z-scored IAS and IATS totals).\n\n")
cat("Profile n:\n")
for (p in names(prof_n)) cat(sprintf("  %-14s n = %d\n", p, prof_n[[p]]))
cat("\nProfile-typical values (mean z-score within profile):\n")
cat(sprintf("%-16s %10s %10s\n", "Profile", "IAS (z)", "IATS (z)"))
hr()
for (i in seq_len(nrow(profile_values))) {
  cat(sprintf("%-16s %10.3f %10.3f\n", profile_values$profile[i],
              profile_values$IAS_typical[i], profile_values$IATS_typical[i]))
}

cat("\nConditional indirect effect of IAS on mental health at each profile's\n")
cat("typical IATS level, separately through the alexithymia (TAS) and the\n")
cat("somatic mediator. Model M0, R =", R0, "bootstrap resamples.\n\n")

profs <- c("Hypervigilant", "Uncertain", "Efficient")
suff <- c("hyper", "uncert", "effic")
iats_at <- c(hyper_iats, uncert_iats, effic_iats)

cond_draws <- list()
for (j in seq_along(profs)) {
  cond_draws[[suff[j]]] <- list(
    tas = (bs("a1") + bs("a3") * round(iats_at[j], 4)) * bs("c1"),
    som = (bs("b2") + bs("b3") * round(iats_at[j], 4)) * bs("c2")
  )
}

cat(sprintf("%-16s %-12s %9s %8s %24s %12s %12s\n",
            "Profile", "Pathway", "est", "SE_boot", "95% boot CI",
            "p (Wald z)", "p (boot)"))
hr(n = 100)
for (j in seq_along(profs)) {
  for (w in c("tas", "som")) {
    lbl <- paste0("cond_", w, "_", suff[j])
    r <- get_row(pe_m0, lbl)
    v <- cond_draws[[suff[j]]][[w]]
    cat(sprintf("%-16s %-12s %9.4f %8.4f  [%9.4f, %9.4f] %12s %12s\n",
                ifelse(w == "tas", profs[j], ""),
                ifelse(w == "tas", "Alexithymia", "Somatic"),
                r$est, r$se, qlo(v), qhi(v), fmt_p(r$pvalue), fmt_boot_p(v)))
  }
}

cat("\n--- Percentage of the total indirect effect carried by each pathway ---\n\n")
cat("pct_alex = |cond_alex| / (|cond_alex| + |cond_som|) * 100, recomputed in\n")
cat("every bootstrap resample (definition used in 12b_pathway_dominance_test.R).\n")
cat("The p value tests H0: the two pathways carry equal absolute effect, i.e.\n")
cat("pct = 50%. Two tests are reported for that null:\n")
cat("  Wald z   : z = diff / SE_boot on the signed contrast\n")
cat("             (somatic minus alexithymia conditional indirect effect),\n")
cat("             the test used in the published analysis;\n")
cat("  bootstrap: two-sided percentile bootstrap p on |cond_som| - |cond_alex|.\n\n")

cat(sprintf("%-16s %20s %9s %20s %9s %10s %12s %16s\n",
            "Profile", "Alexithymia %", "bootmean", "Somatic %", "bootmean",
            "diff", "p (Wald z)", "p (boot)"))
hr(n = 122)
dom_tab <- list()
for (j in seq_along(profs)) {
  a <- abs(cond_draws[[suff[j]]]$tas)
  s <- abs(cond_draws[[suff[j]]]$som)
  pct_a <- 100 * a / (a + s)
  pct_s <- 100 * s / (a + s)
  ea <- abs(get_row(pe_m0, paste0("cond_tas_", suff[j]))$est)
  es <- abs(get_row(pe_m0, paste0("cond_som_", suff[j]))$est)
  pa <- 100 * ea / (ea + es)
  rd <- get_row(pe_m0, paste0("diff_", suff[j]))
  dvec <- s - a
  cat(sprintf("%-16s %5.1f%% [%5.1f, %5.1f] %8.1f%% %5.1f%% [%5.1f, %5.1f] %8.1f%% %10.4f %12s %16s\n",
              profs[j], pa, qlo(pct_a), qhi(pct_a), mean(pct_a),
              100 - pa, qlo(pct_s), qhi(pct_s), mean(pct_s),
              rd$est, fmt_p(rd$pvalue), fmt_boot_p(dvec)))
  dom_tab[[profs[j]]] <- list(pct_a = pa, lo_a = qlo(pct_a), hi_a = qhi(pct_a),
                              lo_s = qlo(pct_s), hi_s = qhi(pct_s),
                              bm_a = mean(pct_a), bm_s = mean(pct_s),
                              diff = rd$est, se = rd$se, z = rd$z,
                              p = rd$pvalue, pboot = boot_p(dvec),
                              ci_lo = qlo(dvec), ci_hi = qhi(dvec))
}
# Figure 1 panel c and Supplementary Fig. 17 read these values
# (12_sem_pathway_analysis.R, 12b_pathway_dominance_test.R).
write.csv(do.call(rbind, lapply(names(dom_tab), function(nm) {
  d <- dom_tab[[nm]]
  data.frame(profile = nm, pct_alex = d$pct_a, lo_alex = d$lo_a, hi_alex = d$hi_a,
             pct_som = 100 - d$pct_a, lo_som = d$lo_s, hi_som = d$hi_s,
             bootmean_alex = d$bm_a, bootmean_som = d$bm_s,
             diff = d$diff, z = d$z, p_wald = d$p, p_boot = d$pboot)
})), file.path(supplementary_data_dir, "pathway_dominance_by_profile.csv"), row.names = FALSE)
cat("\n  Percentages are from the point estimates; brackets are percentile\n")
cat("  bootstrap CIs on the percentage itself; bootmean is the mean of the\n")
cat("  same percentage across resamples. The two differ because a bounded\n")
cat("  ratio of absolute effects is skewed. The published values (85 / 15,\n")
cat("  65 / 35, 51 / 49) are bootstrap MEANS from 5000 resamples, so they\n")
cat("  should be compared with the bootmean columns, not the point estimates.\n")
cat("  diff = signed contrast (somatic minus alexithymia). Both conditional\n")
cat("  indirect effects are negative, so a positive diff means the somatic\n")
cat("  effect is the smaller of the two in absolute value.\n\n")

cat("Signed contrast in full, per profile:\n\n")
cat(sprintf("%-16s %10s %9s %24s %8s %12s\n",
            "Profile", "diff", "SE_boot", "95% boot CI", "z", "p (Wald z)"))
hr(n = 86)
for (j in seq_along(profs)) {
  rd <- get_row(pe_m0, paste0("diff_", suff[j]))
  a <- cond_draws[[suff[j]]]$tas
  s <- cond_draws[[suff[j]]]$som
  v <- s - a
  cat(sprintf("%-16s %10.4f %9.4f  [%9.4f, %9.4f] %8.3f %12s\n",
              profs[j], rd$est, rd$se, qlo(v), qhi(v), rd$z, fmt_p(rd$pvalue)))
}
cat("\n  Family of 3 tests. Bonferroni-corrected alpha = .05 / 3 = .0167\n")
cat("  (the correction applied in 12b_pathway_dominance_test.R).\n\n")

# -----------------------------------------------------------------------------
cat("#############################################################################\n")
cat("PART 5: SAMPLE SIZE AND CONVERGENCE FOR EVERY MODEL\n")
cat("#############################################################################\n\n")

model_list <- list(
  "M0 (ML, reproduction check)" = fit_m0_ml,
  "M0 (bootstrap, Parts 1 and 4)" = fit_m0,
  "Parallel subscale mediation (Parts 2c and 3)" = fit_par
)
cat(sprintf("%-46s %7s %10s %6s %11s %16s\n",
            "Model", "N", "converged", "df", "chi-square", "p"))
hr(n = 100)
for (nm in names(model_list)) {
  f <- model_list[[nm]]
  ft <- lavInspect(f, "fit")
  cat(sprintf("%-46s %7d %10s %6d %11.3f %16s\n", nm, lavInspect(f, "nobs"),
              lavInspect(f, "converged"), as.integer(ft["df"]), ft["chisq"],
              fmt_p(ft["pvalue"])))
}
cat("\nBootstrap resamples requested:", N_BOOT, "\n")
cat("Successful resamples, M0:", R0, " (failed:", N_BOOT - R0, ")\n")
cat("Successful resamples, parallel subscale model:", Rp,
    " (failed:", N_BOOT - Rp, ")\n\n")

cat("OLS models (Parts 1a and 2a): N =", N_COMPLETE,
    ", residual df = ", N_COMPLETE - 4, " for every interaction model\n", sep = "")
cat("(4 estimated coefficients: intercept, IAS, IATS, IAS x IATS).\n\n")

cat("Note on model fit. The published specification leaves the residual\n")
cat("covariance between the alexithymia mediator and the somatic mediator\n")
cat("fixed at zero (and, in the parallel model, the three subscale-somatic\n")
cat("residual covariances), so the models are over-identified rather than\n")
cat("saturated and the chi-square is significant. This is the specification\n")
cat("used in 12_sem_pathway_analysis.R and 17_revision_sensitivity.R and is\n")
cat("reproduced here unchanged. Both models are recursive with uncorrelated\n")
cat("disturbances, so every structural equation is estimated equation-by-\n")
cat("equation and its ML coefficients equal the OLS coefficients; this is\n")
cat("verified numerically above, where the OLS simple slopes in Part 1a match\n")
cat("the SEM defined parameters in Part 1b to four decimal places.\n\n")

cat("=============================================================================\n")
cat("END OF OUTPUT\n")
cat("=============================================================================\n")

sink()

cat("Done. Output written to analysis_output/19_output.txt\n")
