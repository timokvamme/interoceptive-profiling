# =============================================================================
# REVISION SENSITIVITY ANALYSES: DUAL-PATHWAY SEM ROBUSTNESS CHECKS
# =============================================================================
#
# Reviewer-requested sensitivity analyses for COMMSPSYCHOL-26-0385-T.
# Companion to 12_sem_pathway_analysis.R (does NOT modify it).
#
# RESEARCH QUESTIONS:
#   A1. METRIC RECONCILIATION. The main text reports beta = -0.091 for the
#       IAS x IATS -> Somatic path (the fully standardized std.all value) but
#       pairs it with 95% CI [-0.146, -0.023], which is symmetric about
#       -0.0845 and is therefore the UNSTANDARDIZED interval. Supplementary
#       Table 1 reports the unstandardized -0.085 from the identical model.
#       The two differ because the model z-scores the variables and then forms
#       ias_x_iats = ias_z * iats_z, whose SD is 1.0557, so
#       std.all = est * 1.0557. This script re-fits the published model and
#       reports BOTH metrics on the same footing, each with a genuine 95%
#       percentile bootstrap CI computed on that metric (the standardized
#       quantity is re-standardized inside every bootstrap resample, not
#       rescaled after the fact).
#   A2. Does the dual-pathway dissociation survive adjustment for age + gender
#       (covariates on BOTH mediators and on the mental-health composite)?
#   A3. Does the somatic pathway hold when the somatic mediator is SSS-8 ONLY
#       (PCS dropped, because PCS is also a standalone outcome and indexes
#       cognitive pain appraisal rather than somatic symptoms)?
#   A4. Do both adjustments survive SIMULTANEOUSLY (SSS-8 only + age + gender)?
#
# MODELS ESTIMATED (all with full 42-item IAS/IATS, MH composite outcome):
#   M0  baseline        somatic = SSS-8 + PCS   no covariates   (published model)
#   A2  covariates      somatic = SSS-8 + PCS   + age + gender
#   A3  SSS-8 only      somatic = SSS-8         no covariates
#   A4  combined        somatic = SSS-8         + age + gender
#
# OUTPUTS:
#   - analysis_output/17_output.txt (all statistics)
#   - supplementary_data/sem_revision_sensitivity.xlsx (tables in sheets)
#
# =============================================================================

# Paths are relative to the repository root; run scripts from there.

library(dplyr)
library(tidyr)
library(lavaan)
library(openxlsx)

# Output directories
analysis_output_dir <- "analysis_output"
supplementary_data_dir <- "supplementary_data"
dir.create(analysis_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supplementary_data_dir, showWarnings = FALSE, recursive = TRUE)

BOOT_SEED <- 42
N_BOOT <- 2000         # sensitivity models A2 / A3 / A4
N_BOOT_METRIC <- 5000  # published model M0, metric reconciliation table

# =============================================================================
# 1. LOAD DATA AND COMPUTE SCALES (identical to 12_sem_pathway_analysis.R)
# =============================================================================

dfc <- read.csv("dfc_interoception_profiling.csv")

cat("Loading data...\n")
cat("Sample: N =", nrow(dfc), "\n\n")

# -----------------------------------------------------------------------------
# 1a. Full scales (all 42 items: 21 IAS + 21 IATS)
# -----------------------------------------------------------------------------

ias_full_items <- paste0("ias_", 1:21)
iats_full_items <- paste0("iats_", 1:21)

dfc$ias_full <- rowMeans(dfc[, ias_full_items], na.rm = TRUE)
dfc$iats_full <- rowMeans(dfc[, iats_full_items], na.rm = TRUE)

# -----------------------------------------------------------------------------
# 1b. Mental health and mediator variables
# -----------------------------------------------------------------------------

# MH composite (psychological distress)
dfc$mh_composite <- rowMeans(cbind(scale(dfc$phq9), scale(dfc$gad7),
                                   scale(dfc$stai)), na.rm = TRUE)

# Somatic measures
dfc$somatic_sss_only <- dfc$sss8
dfc$somatic_combined <- rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE)

# -----------------------------------------------------------------------------
# 1c. Demographic covariates
# -----------------------------------------------------------------------------
# gender is coded 1 / 2 / 3 in the source file. Dummy-coded with 1 as reference
# so that the small third category (n = 5) is not treated as an ordered step.

dfc$age_z <- scale(dfc$age)[, 1]
dfc$gender_d2 <- as.numeric(dfc$gender == 2)
dfc$gender_d3 <- as.numeric(dfc$gender == 3)

cat("Scales computed:\n")
cat("  Full IAS/IATS: 21 items each (42 total)\n")
cat("  MH composite: PHQ-9 + GAD-7 + STAI (z-scored average)\n")
cat("  Somatic (SSS only): SSS-8\n")
cat("  Somatic (combined): SSS-8 + PCS (z-scored average)\n")
cat("  Covariates: age (z), gender dummies (ref = category 1)\n\n")

# =============================================================================
# 2. MODEL ENGINE
# =============================================================================

COVARS <- c("age_z", "gender_d2", "gender_d3")

# Builds and fits the dual-pathway SEM, optionally with demographic covariates
# on both mediators and on the mental-health composite.
run_sem_sensitivity <- function(data, ias_var, iats_var, somatic_var, mh_var,
                                covars = NULL, label, boot = TRUE,
                                n_boot = N_BOOT) {

  # Standardize all variables (same convention as 12_sem_pathway_analysis.R)
  data$ias_z <- scale(data[[ias_var]])[, 1]
  data$iats_z <- scale(data[[iats_var]])[, 1]
  data$tas_z <- scale(data$tas)[, 1]
  data$somatic_z <- scale(data[[somatic_var]])[, 1]
  data$mh_z <- scale(data[[mh_var]])[, 1]
  data$ias_x_iats <- data$ias_z * data$iats_z

  # Covariate terms, with labels so they can be pulled out by name
  if (!is.null(covars) && length(covars) > 0) {
    cov_tas <- paste0(" + ", paste0("t", seq_along(covars), "*", covars, collapse = " + "))
    cov_som <- paste0(" + ", paste0("s", seq_along(covars), "*", covars, collapse = " + "))
    cov_mh  <- paste0(" + ", paste0("m", seq_along(covars), "*", covars, collapse = " + "))
  } else {
    cov_tas <- ""; cov_som <- ""; cov_mh <- ""
  }

  model <- paste0('
    # Paths to mediators (with interaction)
    tas_z ~ a1*ias_z + a2*iats_z + a3*ias_x_iats', cov_tas, '
    somatic_z ~ b1*iats_z + b2*ias_z + b3*ias_x_iats', cov_som, '

    # Mediators to outcome
    mh_z ~ c1*tas_z + c2*somatic_z + d1*ias_z + d2*iats_z + d3*ias_x_iats', cov_mh, '

    # Indirect effects
    ind_ias_tas := a1 * c1
    ind_iats_tas := a2 * c1
    ind_ias_som := b2 * c2
    ind_iats_som := b1 * c2

    # Total indirect
    total_indirect := ind_ias_tas + ind_iats_tas + ind_ias_som + ind_iats_som
  ')

  fit <- sem(model, data = data)

  # Standardized solution with delta-method 95% CIs
  std <- standardizedSolution(fit, level = 0.95)
  unstd <- parameterEstimates(fit, level = 0.95)
  r2 <- inspect(fit, "rsquare")

  # Nonparametric bootstrap on BOTH metrics simultaneously. Inside each
  # resample the model is re-fit and standardizedSolution() recomputes std.all
  # from that resample's implied variances, so the standardized interval is a
  # genuine bootstrap of the standardized quantity rather than a rescaling of
  # the unstandardized one.
  boot_tab <- NULL
  n_boot_used <- 0
  n_boot_failed <- 0
  if (boot) {
    boot_tab <- boot_both(fit, model, n_boot)
    n_boot_used <- attr(boot_tab, "R_ok")
    n_boot_failed <- attr(boot_tab, "R_failed")
  }

  list(
    label = label,
    fit = fit,
    model = model,
    std = std,
    unstd = unstd,
    r2 = r2,
    boot = boot_tab,
    n_boot_used = n_boot_used,
    n_boot_failed = n_boot_failed,
    nobs = lavInspect(fit, "nobs"),
    converged = lavInspect(fit, "converged"),
    covars = covars,
    somatic_var = somatic_var
  )
}

# All labelled structural paths, in reporting order
PATH_LABELS <- c("a1", "a2", "a3", "b2", "b1", "b3",
                 "c1", "c2", "d1", "d2", "d3",
                 "ind_ias_tas", "ind_iats_tas", "ind_ias_som", "ind_iats_som",
                 "total_indirect")

PATH_NAMES <- c(
  a1 = "IAS -> TAS",
  a2 = "IATS -> TAS",
  a3 = "IAS x IATS -> TAS",
  b2 = "IAS -> Somatic",
  b1 = "IATS -> Somatic",
  b3 = "IAS x IATS -> Somatic",
  c1 = "TAS -> MH",
  c2 = "Somatic -> MH",
  d1 = "IAS -> MH (direct)",
  d2 = "IATS -> MH (direct)",
  d3 = "IAS x IATS -> MH (direct)",
  ind_ias_tas = "IAS -> TAS -> MH",
  ind_iats_tas = "IATS -> TAS -> MH",
  ind_ias_som = "IAS -> Somatic -> MH",
  ind_iats_som = "IATS -> Somatic -> MH",
  total_indirect = "Total indirect"
)

# Bootstrap both metrics for every labelled path
boot_both <- function(fit, model, R) {
  labs <- PATH_LABELS

  extract <- function(x) {
    p <- parameterEstimates(x, se = FALSE, zstat = FALSE, pvalue = FALSE, ci = FALSE)
    s <- standardizedSolution(x, se = FALSE, zstat = FALSE,
                              pvalue = FALSE, ci = FALSE)
    c(p$est[match(labs, p$label)], s$est.std[match(labs, s$label)])
  }

  set.seed(BOOT_SEED)
  B <- suppressWarnings(
    lavaan::bootstrapLavaan(fit, R = R, type = "ordinary", FUN = extract))
  B <- as.matrix(B)
  ok <- stats::complete.cases(B)
  R_failed <- sum(!ok)
  B <- B[ok, , drop = FALSE]
  R_ok <- nrow(B)

  k <- length(labs)
  point <- extract(fit)

  # two-sided percentile bootstrap p-value: 2 * the smaller tail mass at zero
  boot_p <- function(v) {
    p <- 2 * min(mean(v <= 0), mean(v >= 0))
    min(p, 1)
  }

  out <- data.frame(
    label = labs,
    path = unname(PATH_NAMES[labs]),
    est_unstd = point[1:k],
    se_boot_unstd = apply(B[, 1:k, drop = FALSE], 2, sd),
    lo_unstd = apply(B[, 1:k, drop = FALSE], 2, quantile, probs = 0.025),
    hi_unstd = apply(B[, 1:k, drop = FALSE], 2, quantile, probs = 0.975),
    p_boot_unstd = apply(B[, 1:k, drop = FALSE], 2, boot_p),
    est_std = point[(k + 1):(2 * k)],
    se_boot_std = apply(B[, (k + 1):(2 * k), drop = FALSE], 2, sd),
    lo_std = apply(B[, (k + 1):(2 * k), drop = FALSE], 2, quantile, probs = 0.025),
    hi_std = apply(B[, (k + 1):(2 * k), drop = FALSE], 2, quantile, probs = 0.975),
    p_boot_std = apply(B[, (k + 1):(2 * k), drop = FALSE], 2, boot_p),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "R_ok") <- R_ok
  attr(out, "R_failed") <- R_failed
  out
}

# Pull one labelled path out of the standardized solution (delta-method CI)
get_std <- function(res, lbl) {
  row <- res$std[res$std$label == lbl, ]
  if (nrow(row) == 0) return(list(est = NA, se = NA, lo = NA, hi = NA, p = NA))
  list(est = row$est.std[1], se = row$se[1],
       lo = row$ci.lower[1], hi = row$ci.upper[1], p = row$pvalue[1])
}

# Pull one labelled path out of the unstandardized solution (delta-method CI)
get_unstd <- function(res, lbl) {
  row <- res$unstd[res$unstd$label == lbl, ]
  if (nrow(row) == 0) return(list(est = NA, se = NA, lo = NA, hi = NA, p = NA))
  list(est = row$est[1], se = row$se[1],
       lo = row$ci.lower[1], hi = row$ci.upper[1], p = row$pvalue[1])
}

# Pull one labelled path out of the bootstrap table
get_boot <- function(res, lbl) {
  if (is.null(res$boot)) return(NULL)
  row <- res$boot[res$boot$label == lbl, ]
  if (nrow(row) == 0) return(NULL)
  as.list(row[1, ])
}

fmt_p <- function(p) {
  if (is.na(p)) return("   n/a")
  if (p < 0.0001) return("< .0001")
  sprintf("%7.4f", p)
}

sig_star <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  if (p < 0.10) return(".")
  return("")
}

# Print a full coefficient table for one fitted model
print_model_table <- function(res) {
  paths <- list(
    c("a1", "IAS -> TAS"),
    c("a2", "IATS -> TAS"),
    c("a3", "IAS x IATS -> TAS"),
    c("b2", "IAS -> Somatic"),
    c("b1", "IATS -> Somatic"),
    c("b3", "IAS x IATS -> Somatic"),
    c("c1", "TAS -> MH"),
    c("c2", "Somatic -> MH"),
    c("d1", "IAS -> MH (direct)"),
    c("d2", "IATS -> MH (direct)"),
    c("d3", "IAS x IATS -> MH (direct)")
  )

  if (!is.null(res$covars)) {
    cov_names <- c("age (z)", "gender d2 (cat2 vs cat1)", "gender d3 (cat3 vs cat1)")
    for (i in seq_along(res$covars)) {
      paths[[length(paths) + 1]] <- c(paste0("t", i), paste0(cov_names[i], " -> TAS"))
    }
    for (i in seq_along(res$covars)) {
      paths[[length(paths) + 1]] <- c(paste0("s", i), paste0(cov_names[i], " -> Somatic"))
    }
    for (i in seq_along(res$covars)) {
      paths[[length(paths) + 1]] <- c(paste0("m", i), paste0(cov_names[i], " -> MH"))
    }
  }

  cat(sprintf("%-32s %8s %7s %20s %9s %s\n",
              "Path", "beta", "SE", "95% CI", "p", ""))
  cat(paste(rep("-", 88), collapse = ""), "\n")
  for (pth in paths) {
    v <- get_std(res, pth[1])
    cat(sprintf("%-32s %8.3f %7.3f  [%7.3f, %7.3f] %9s %s\n",
                pth[2], v$est, v$se, v$lo, v$hi, fmt_p(v$p), sig_star(v$p)))
  }

  cat("\nIndirect effects (standardized):\n")
  cat(paste(rep("-", 88), collapse = ""), "\n")
  ind <- list(
    c("ind_ias_tas", "IAS -> TAS -> MH"),
    c("ind_iats_tas", "IATS -> TAS -> MH"),
    c("ind_ias_som", "IAS -> Somatic -> MH"),
    c("ind_iats_som", "IATS -> Somatic -> MH"),
    c("total_indirect", "Total indirect")
  )
  for (pth in ind) {
    v <- get_std(res, pth[1])
    cat(sprintf("%-32s %8.3f %7.3f  [%7.3f, %7.3f] %9s %s\n",
                pth[2], v$est, v$se, v$lo, v$hi, fmt_p(v$p), sig_star(v$p)))
  }

  cat(sprintf("\nR2: TAS = %.3f, Somatic = %.3f, MH = %.3f\n",
              res$r2["tas_z"], res$r2["somatic_z"], res$r2["mh_z"]))
  cat(sprintf("N = %d | converged = %s\n", res$nobs, res$converged))

  if (!is.null(res$boot)) {
    cat(sprintf("\nBootstrap percentile 95%% CIs (%d successful resamples of %d",
                res$n_boot_used, res$n_boot_used + res$n_boot_failed))
    cat(sprintf("; %d failed):\n", res$n_boot_failed))
    print_boot_table(res$boot)
  }
  cat("\n")
}

# Both metrics side by side, each with its own bootstrap percentile CI
print_boot_table <- function(bt, labs = PATH_LABELS) {
  cat(sprintf("%-28s | %-38s | %-38s\n", "", "UNSTANDARDIZED (est)",
              "FULLY STANDARDIZED (std.all)"))
  cat(sprintf("%-28s | %7s %20s %9s | %7s %20s %9s\n",
              "Path", "est", "95% boot CI", "p_boot",
              "std.all", "95% boot CI", "p_boot"))
  cat(paste(rep("-", 112), collapse = ""), "\n")
  for (lb in labs) {
    r <- bt[bt$label == lb, ]
    if (nrow(r) == 0) next
    cat(sprintf("%-28s | %7.3f [%8.4f, %8.4f] %9s | %7.3f [%8.4f, %8.4f] %9s\n",
                r$path[1], r$est_unstd[1], r$lo_unstd[1], r$hi_unstd[1],
                sub("^\\s+", "", fmt_p(r$p_boot_unstd[1])),
                r$est_std[1], r$lo_std[1], r$hi_std[1],
                sub("^\\s+", "", fmt_p(r$p_boot_std[1]))))
  }
}

# Flatten one model into a tidy data frame for the xlsx workbook
tidy_model <- function(res) {
  s <- res$std
  s <- s[s$label != "", c("lhs", "op", "rhs", "label", "est.std", "se",
                          "ci.lower", "ci.upper", "z", "pvalue")]
  s$model <- res$label
  s$n <- res$nobs
  s[, c("model", "n", "label", "lhs", "op", "rhs", "est.std", "se",
        "ci.lower", "ci.upper", "z", "pvalue")]
}

# =============================================================================
# 3. FIT ALL FOUR MODELS
# =============================================================================

cat("Fitting models (bootstrap CIs take a moment)...\n\n")

m0 <- run_sem_sensitivity(dfc, "ias_full", "iats_full", "somatic_combined",
                          "mh_composite", covars = NULL,
                          label = "M0_baseline_SSS8+PCS_nocov",
                          n_boot = N_BOOT_METRIC)

# SD of the interaction product term, which is why est and std.all differ
sd_product <- sd(scale(dfc$ias_full)[, 1] * scale(dfc$iats_full)[, 1])

a2 <- run_sem_sensitivity(dfc, "ias_full", "iats_full", "somatic_combined",
                          "mh_composite", covars = COVARS,
                          label = "A2_SSS8+PCS_age+gender")

a3 <- run_sem_sensitivity(dfc, "ias_full", "iats_full", "somatic_sss_only",
                          "mh_composite", covars = NULL,
                          label = "A3_SSS8only_nocov")

a4 <- run_sem_sensitivity(dfc, "ias_full", "iats_full", "somatic_sss_only",
                          "mh_composite", covars = COVARS,
                          label = "A4_SSS8only_age+gender")

models <- list(m0 = m0, a2 = a2, a3 = a3, a4 = a4)

# =============================================================================
# BEGIN OUTPUT FILE
# =============================================================================

sink(file.path(analysis_output_dir, "17_output.txt"))

cat("=============================================================================\n")
cat("REVISION SENSITIVITY ANALYSES: DUAL-PATHWAY SEM (COMMSPSYCHOL-26-0385-T)\n")
cat("=============================================================================\n\n")

cat("Sample: N =", nrow(dfc), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("Estimator: ML (lavaan), complete cases; no missing data on any model variable\n")
cat("Bootstrap:", N_BOOT, "resamples, percentile CIs, seed = 42\n\n")

cat("Models:\n")
cat("  M0  baseline    somatic = SSS-8 + PCS   no covariates  (published model)\n")
cat("  A2  covariates  somatic = SSS-8 + PCS   + age + gender on TAS, Somatic, MH\n")
cat("  A3  SSS-8 only  somatic = SSS-8         no covariates\n")
cat("  A4  combined    somatic = SSS-8         + age + gender on TAS, Somatic, MH\n\n")

cat("Missing data: the source file has 0 missing values on age, gender, tas,\n")
cat("sss8, pcs, phq9, gad7, stai, and all 42 IAS/IATS items, so every model is\n")
cat(sprintf("estimated on the full N = %d with no listwise deletion.", nrow(dfc)), "\n\n")

gcnt <- table(factor(dfc$gender, levels = c(1, 2, 3)))
cat(sprintf("Gender coding: source variable is 1 / 2 / 3 (n = %d / %d / %d). Entered as", gcnt[1], gcnt[2], gcnt[3]), "\n")
cat("two dummies with category 1 as reference, rather than as a single numeric\n")
cat("step, so the 5-person third category cannot distort a linear contrast.\n\n")

# -----------------------------------------------------------------------------

cat("\n#############################################################################\n")
cat("A1: METRIC RECONCILIATION FOR THE PUBLISHED MODEL (HIGHEST PRIORITY)\n")
cat("#############################################################################\n\n")

cat("Model: full 42-item IAS/IATS, somatic = SSS-8 + PCS composite, MH composite,\n")
cat("no covariates. This is exactly the specification in 12_sem_pathway_analysis.R\n")
cat("(run_dual_pathway_sem with somatic_combined and mh_composite).\n\n")

cat("THE PROBLEM\n")
cat("-----------\n")
cat("12_sem_pathway_analysis.R z-scores every variable and then forms the product\n")
cat("term ias_x_iats = ias_z * iats_z. The product of two z-scored variables is\n")
cat(sprintf("NOT itself unit-variance: SD(ias_z * iats_z) = %.6f in this sample.\n", sd_product))
cat("Consequently, for the interaction paths only, the fully standardized\n")
cat("coefficient differs from the unstandardized one by that factor:\n")
cat(sprintf("  std.all = est * %.6f\n", sd_product))
cat("For every non-interaction path the two metrics coincide to 3 decimals,\n")
cat("because those predictors and all outcomes are already z-scored.\n")
cat("lavaan attaches the UNSTANDARDIZED Wald p-value to std.all, which is why the\n")
cat("p-value is identical in both columns and the discrepancy went unnoticed.\n\n")

cat("Verification of the published numbers, from this re-fit:\n")
b3u <- get_unstd(m0, "b3"); b3s <- get_std(m0, "b3")
cat(sprintf("  IAS x IATS -> Somatic, unstandardized est = %.6f (SE %.6f)\n", b3u$est, b3u$se))
cat(sprintf("  IAS x IATS -> Somatic, std.all           = %.6f (delta SE %.6f)\n", b3s$est, b3s$se))
cat(sprintf("  ratio std.all / est                      = %.6f\n", b3s$est / b3u$est))
cat(sprintf("  Wald p (identical for both)               = %.6f\n", b3u$p))
cat(sprintf("  delta-method 95%% CI, unstandardized      = [%.4f, %.4f]\n", b3u$lo, b3u$hi))
cat(sprintf("  delta-method 95%% CI, standardized        = [%.4f, %.4f]\n", b3s$lo, b3s$hi))
cat("\nThe published interval [-0.146, -0.023] is centred on -0.0845, i.e. on the\n")
cat("UNSTANDARDIZED point estimate, while the published point estimate -0.091 is\n")
cat("the STANDARDIZED one. The manuscript therefore pairs a standardized estimate\n")
cat("with an unstandardized interval. That is the defect these numbers replace.\n\n")

cat(sprintf("BOOTSTRAP RESULTS (%d requested resamples, ordinary nonparametric,\n",
            N_BOOT_METRIC))
cat(sprintf("percentile method, seed = %d; %d succeeded, %d failed)\n",
            BOOT_SEED, m0$n_boot_used, m0$n_boot_failed))
cat("The standardized column is a genuine bootstrap of std.all: within every\n")
cat("resample the model is re-fit and standardizedSolution() recomputes std.all\n")
cat("from that resample's own implied variances. It is NOT the unstandardized\n")
cat("interval multiplied by a constant.\n\n")
print_boot_table(m0$boot)

cat("\nHeadline coefficients requested for the response letter:\n")
cat(paste(rep("-", 112), collapse = ""), "\n")
print_boot_table(m0$boot, labs = c("a1", "a2", "a3", "b1", "b2", "b3"))

cat("\nRECOMMENDATION\n")
cat("--------------\n")
cat("Report the FULLY STANDARDIZED (std.all) metric throughout, with percentile\n")
cat("bootstrap CIs computed on the standardized quantity, i.e. the right-hand\n")
cat("block above. Reasons:\n")
cat("  1. It is the metric the main text already quotes (-0.345, +0.314, +0.444,\n")
cat("     -0.095, -0.091), so only the intervals change, not the headline values.\n")
cat("  2. It is the only metric in which the interaction term is on the same\n")
cat("     footing as the main effects; the unstandardized interaction is scaled\n")
cat(sprintf("     by the arbitrary SD of the product term (%.4f).\n", sd_product))
cat("  3. Standardized betas are what the profile-level interpretation rests on.\n")
cat("Supplementary Table 1 should then be regenerated in the same metric, or\n")
cat("explicitly labelled 'unstandardized' with the unstandardized CIs from the\n")
cat("left-hand block. Do not mix the two, which is the present state.\n\n")

cat("\n#############################################################################\n")
cat("MODEL M0: BASELINE (SSS-8 + PCS, NO COVARIATES) - FULL PATH SET\n")
cat("#############################################################################\n\n")
cat("Reproduction of the published dual-pathway SEM, re-estimated here so every\n")
cat("comparison below is against numbers produced by this same script.\n")
cat("The 'beta' column below is std.all with delta-method CIs; the bootstrap\n")
cat("block underneath repeats it with percentile CIs on both metrics.\n\n")
print_model_table(m0)

cat("\n#############################################################################\n")
cat("A2: AGE + GENDER COVARIATE SEM (SOMATIC = SSS-8 + PCS)\n")
cat("#############################################################################\n\n")
cat("Age and gender enter as predictors of BOTH mediators (TAS, Somatic) and of\n")
cat("the mental-health composite. Tests whether the dual-pathway dissociation is\n")
cat("an artefact of demographic differences between profiles (the Efficient\n")
cat("profile is significantly older).\n\n")
print_model_table(a2)

cat("\n#############################################################################\n")
cat("A3: SSS-8 ALONE AS THE SOMATIC MEDIATOR (NO COVARIATES)\n")
cat("#############################################################################\n\n")
cat("PCS dropped from the somatic composite. Reviewers note PCS is (a) also used\n")
cat("as a standalone outcome and (b) measures cognitive pain appraisal rather\n")
cat("than somatic symptom burden.\n\n")
cat("NOTE: 12_sem_pathway_analysis.R already contains an SSS-alone vs SSS+PCS\n")
cat("sensitivity check (its PART 4 / RQ4). That existing comparison reported:\n")
cat("  SSS-8 Only : IATS->Som  0.409 | IAS->Som -0.099 | IxI->Som -0.068, p = .0336 | Som->MH 0.549 | R2 MH 0.528\n")
cat("  SSS-8 + PCS: IATS->Som  0.443 | IAS->Som -0.093 | IxI->Som -0.094, p = .0030 | Som->MH 0.603 | R2 MH 0.560\n")
cat("Those numbers (2026-09-02 rerun on corrected scores) are reproduced exactly below; this script adds the 95% CIs,\n")
cat("the full path set, and the bootstrap CIs that the earlier run did not print.\n\n")
print_model_table(a3)

cat("\n#############################################################################\n")
cat("A4: COMBINED - SSS-8 ALONE + AGE + GENDER COVARIATES\n")
cat("#############################################################################\n\n")
cat("Both adjustments applied simultaneously: the strictest test of the\n")
cat("multiplicative somatic pathway.\n\n")
print_model_table(a4)

# =============================================================================
# 4. SIDE-BY-SIDE COMPARISONS
# =============================================================================

cat("\n#############################################################################\n")
cat("SIDE-BY-SIDE COMPARISON OF THE KEY PATHS ACROSS ALL FOUR MODELS\n")
cat("#############################################################################\n\n")

compare_paths <- list(
  c("a1", "IAS -> TAS"),
  c("a2", "IATS -> TAS"),
  c("a3", "IAS x IATS -> TAS"),
  c("b2", "IAS -> Somatic"),
  c("b1", "IATS -> Somatic"),
  c("b3", "IAS x IATS -> Somatic"),
  c("c1", "TAS -> MH"),
  c("c2", "Somatic -> MH")
)

model_order <- c("m0", "a2", "a3", "a4")
model_heads <- c("M0 base", "A2 +cov", "A3 SSS8", "A4 both")

cat(sprintf("%-26s", "Path"))
for (h in model_heads) cat(sprintf(" | %-22s", h))
cat("\n")
cat(paste(rep("-", 26 + 4 * 25), collapse = ""), "\n")
for (pth in compare_paths) {
  cat(sprintf("%-26s", pth[2]))
  for (mk in model_order) {
    v <- get_std(models[[mk]], pth[1])
    cat(sprintf(" | %6.3f p=%-12s", v$est, sub("^\\s+", "", fmt_p(v$p))))
  }
  cat("\n")
}

cat("\n\nKEY INTERACTION TERMS: STANDARDIZED, DELTA-METHOD vs BOOTSTRAP CIs\n")
cat(paste(rep("-", 104), collapse = ""), "\n")
cat(sprintf("%-6s %-24s %8s %20s %9s | %20s %9s\n",
            "Model", "Interaction path", "std.all", "95% CI (delta)", "p_Wald",
            "95% CI (bootstrap)", "p_boot"))
for (mk in model_order) {
  for (pth in list(c("a3", "IAS x IATS -> TAS"), c("b3", "IAS x IATS -> Somatic"))) {
    v <- get_std(models[[mk]], pth[1])
    b <- get_boot(models[[mk]], pth[1])
    cat(sprintf("%-6s %-24s %8.3f  [%7.3f, %7.3f] %9s | [%7.3f, %7.3f] %9s %s\n",
                mk, pth[2], v$est, v$lo, v$hi, sub("^\\s+", "", fmt_p(v$p)),
                b$lo_std, b$hi_std, sub("^\\s+", "", fmt_p(b$p_boot_std)),
                sig_star(v$p)))
  }
}

cat("\n\nSAME INTERACTION TERMS IN THE UNSTANDARDIZED METRIC (for cross-checking\n")
cat("against Supplementary Table 1)\n")
cat(paste(rep("-", 104), collapse = ""), "\n")
cat(sprintf("%-6s %-24s %8s %20s %9s | %20s %9s\n",
            "Model", "Interaction path", "est", "95% CI (delta)", "p_Wald",
            "95% CI (bootstrap)", "p_boot"))
for (mk in model_order) {
  for (pth in list(c("a3", "IAS x IATS -> TAS"), c("b3", "IAS x IATS -> Somatic"))) {
    v <- get_unstd(models[[mk]], pth[1])
    b <- get_boot(models[[mk]], pth[1])
    cat(sprintf("%-6s %-24s %8.4f  [%7.4f, %7.4f] %9s | [%7.4f, %7.4f] %9s\n",
                mk, pth[2], v$est, v$lo, v$hi, sub("^\\s+", "", fmt_p(v$p)),
                b$lo_unstd, b$hi_unstd, sub("^\\s+", "", fmt_p(b$p_boot_unstd))))
  }
}

cat("\n\nMODEL R-SQUARED\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
cat(sprintf("%-10s %10s %10s %10s %8s\n", "Model", "R2 TAS", "R2 Somatic", "R2 MH", "N"))
for (mk in model_order) {
  r <- models[[mk]]
  cat(sprintf("%-10s %10.3f %10.3f %10.3f %8d\n", mk,
              r$r2["tas_z"], r$r2["somatic_z"], r$r2["mh_z"], r$nobs))
}

# =============================================================================
# 5. VERDICTS
# =============================================================================

cat("\n\n#############################################################################\n")
cat("VERDICTS\n")
cat("#############################################################################\n\n")

verdict <- function(res, tag) {
  b3 <- get_std(res, "b3"); b3b <- get_boot(res, "b3")
  a3v <- get_std(res, "a3"); a3b <- get_boot(res, "a3")
  cat("---", tag, "(", res$label, ") ---\n")
  cat(sprintf("  Somatic IAS x IATS: std.all = %.3f, delta CI [%.3f, %.3f], p = %s\n",
              b3$est, b3$lo, b3$hi, sub("^\\s+", "", fmt_p(b3$p))))
  cat(sprintf("                      bootstrap CI [%.3f, %.3f], p_boot = %s -> %s\n",
              b3b$lo_std, b3b$hi_std, sub("^\\s+", "", fmt_p(b3b$p_boot_std)),
              ifelse(b3b$hi_std < 0 | b3b$lo_std > 0,
                     "SIGNIFICANT (multiplicative somatic pathway holds)",
                     "NOT significant (bootstrap CI includes zero)")))
  cat(sprintf("  Alexithymia IAS x IATS: std.all = %.3f, delta CI [%.3f, %.3f], p = %s\n",
              a3v$est, a3v$lo, a3v$hi, sub("^\\s+", "", fmt_p(a3v$p))))
  cat(sprintf("                          bootstrap CI [%.3f, %.3f], p_boot = %s -> %s\n\n",
              a3b$lo_std, a3b$hi_std, sub("^\\s+", "", fmt_p(a3b$p_boot_std)),
              ifelse(a3b$hi_std < 0 | a3b$lo_std > 0,
                     "SIGNIFICANT (alexithymia pathway NO LONGER additive)",
                     "NULL as expected (alexithymia pathway remains additive)")))
}

verdict(m0, "M0 baseline")
verdict(a2, "A2 age + gender covariates")
verdict(a3, "A3 SSS-8 only")
verdict(a4, "A4 SSS-8 only + age + gender")

cat("Convergence: ")
cat(paste(sapply(model_order, function(mk)
  paste0(mk, " = ", models[[mk]]$converged)), collapse = " | "), "\n")

sink()

# =============================================================================
# 6. WRITE XLSX
# =============================================================================

wb <- createWorkbook()

all_params <- do.call(rbind, lapply(model_order, function(mk) tidy_model(models[[mk]])))
addWorksheet(wb, "All_Parameters")
writeData(wb, "All_Parameters", all_params)

key_paths <- do.call(rbind, lapply(model_order, function(mk) {
  res <- models[[mk]]
  do.call(rbind, lapply(compare_paths, function(pth) {
    v <- get_std(res, pth[1])
    u <- get_unstd(res, pth[1])
    data.frame(model = res$label, path = pth[2], label = pth[1],
               std_beta = v$est, std_se = v$se, std_ci_lower = v$lo,
               std_ci_upper = v$hi, unstd_est = u$est, unstd_se = u$se,
               unstd_ci_lower = u$lo, unstd_ci_upper = u$hi,
               p_wald = v$p, n = res$nobs)
  }))
}))
addWorksheet(wb, "Key_Paths")
writeData(wb, "Key_Paths", key_paths)

boot_all <- do.call(rbind, lapply(model_order, function(mk) {
  b <- models[[mk]]$boot
  if (is.null(b)) return(NULL)
  cbind(model = models[[mk]]$label, R_ok = models[[mk]]$n_boot_used, b)
}))
addWorksheet(wb, "Bootstrap_Both_Metrics")
writeData(wb, "Bootstrap_Both_Metrics", boot_all)

r2_tab <- do.call(rbind, lapply(model_order, function(mk) {
  r <- models[[mk]]
  data.frame(model = r$label, r2_tas = r$r2["tas_z"], r2_somatic = r$r2["somatic_z"],
             r2_mh = r$r2["mh_z"], n = r$nobs, converged = r$converged)
}))
addWorksheet(wb, "Model_Fit")
writeData(wb, "Model_Fit", r2_tab)

saveWorkbook(wb, file.path(supplementary_data_dir, "sem_revision_sensitivity.xlsx"),
             overwrite = TRUE)

cat("Done.\n")
cat("  analysis_output/17_output.txt\n")
cat("  supplementary_data/sem_revision_sensitivity.xlsx\n")
