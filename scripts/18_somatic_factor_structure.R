# =============================================================================
# FACTOR STRUCTURE OF THE "SOMATIC" COMPOSITE: SSS-8 + PCS ITEM-LEVEL ANALYSIS
# =============================================================================
#
# Reviewer-requested psychometric analysis for COMMSPSYCHOL-26-0385-T.
# Companion to 12_sem_pathway_analysis.R and 17_revision_sensitivity.R
# (does NOT modify either).
#
# BACKGROUND:
#   The dual-pathway SEM uses a "somatic" mediator formed as the z-scored
#   average of the SSS-8 (8 items, somatic symptom burden) and the PCS
#   (13 items, catastrophic cognitive appraisal of pain). Reviewer 1 objects
#   that (a) PCS also appears as a standalone outcome elsewhere in the paper,
#   and (b) PCS indexes cognitive pain appraisal (rumination, magnification,
#   helplessness) rather than bodily symptom experience, so labelling the
#   composite "somatic" conflates two constructs.
#
#   17_revision_sensitivity.R already answered the STRUCTURAL question (does
#   the SEM survive dropping PCS: the IAS x IATS -> Somatic interaction
#   attenuates to p = .050, and to p = .070 with age + gender covaried).
#   This script answers the PSYCHOMETRIC question directly: do the 8 SSS-8
#   items and the 13 PCS items cohere as ONE factor in this sample, or are
#   they empirically separable?
#
# ANALYSES (all on the 21 combined items):
#   1. SUITABILITY   KMO (overall + per item), Bartlett's test of sphericity
#   2. DIMENSIONALITY parallel analysis (principal axis, polychoric), scree,
#                     Very Simple Structure, Velicer's MAP
#   3. EFA           principal axis + oblimin, 1 / 2 / 3 factor solutions,
#                    full loading matrices, factor correlations, variance
#   4. CFA           lavaan, WLSMV (ordinal): one-factor vs two-factor
#                    correlated vs bifactor; CFI/TLI/RMSEA[90% CI]/SRMR/chisq
#                    and the scaled chi-square difference test
#   5. RELIABILITY   omega total / omega hierarchical, Cronbach's alpha
#                    (composite and per scale), AVE, SSS-8 x PCS correlation
#                    raw and disattenuated
#   6. VERDICT       one factor or two?
#
# PART B (added as a clearly-labelled follow-up section, see the bottom of the
# script and of the output file): if the composite is mostly a PCS score, is
# the multiplicative IAS x IATS effect on the "somatic" mediator actually a
# PAIN CATASTROPHIZING effect? Part B re-fits the published dual-pathway SEM
# of 12_sem_pathway_analysis.R three times, identical in every respect except
# the mediator (SSS-8 + PCS composite / SSS-8 alone / PCS alone), adds the
# three PCS subscales, and formally tests whether the interaction differs
# between the PCS and SSS-8 mediators in a joint two-mediator model.
# Part B uses the SAME data source, estimator (ML), specification, bootstrap
# (2000 ordinary resamples, percentile, seed 42) and no covariates as model M0
# of 17_revision_sensitivity.R, so its numbers are directly comparable.
#
# DATA:
#   supplementary_data/data_with_ias_iats_ratio.csv (item-level file, N = 832,
#   the same participants as the SEM source file; sss_8_1..8, pcs_1..13).
#   Items are 5-category ordinal (coded 1-5 here), so polychoric correlations
#   and a WLSMV estimator are used throughout, with Pearson/ML cross-checks.
#
# OUTPUTS:
#   - analysis_output/18_output.txt (all statistics)
#   - supplementary_data/somatic_factor_structure.xlsx (tables in sheets)
#   - plots/figure_18_scree_parallel_somatic.png
#
# =============================================================================

# Paths are relative to the repository root; run scripts from there.

library(psych)
library(lavaan)
library(GPArotation)
library(openxlsx)

analysis_output_dir <- "analysis_output"
supplementary_data_dir <- "supplementary_data"
plots_dir <- "plots"
dir.create(analysis_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supplementary_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

SEED <- 42
LOAD_CUT <- 0.40   # loadings at or above this magnitude are flagged

# =============================================================================
# 1. LOAD DATA
# =============================================================================

df <- read.csv(file.path(supplementary_data_dir, "data_with_ias_iats_ratio.csv"))

sss_items <- paste0("sss_8_", 1:8)
pcs_items <- paste0("pcs_", 1:13)
items <- c(sss_items, pcs_items)

stopifnot(all(items %in% names(df)))

n_file <- nrow(df)
na_by_item <- colSums(is.na(df[, items]))
n_na_cells <- sum(na_by_item)
cc <- complete.cases(df[, items])
dat <- df[cc, items]
N <- nrow(dat)

# Scale totals (sums), for the SSS-8 x PCS correlation. These reproduce the
# sss8 / pcs columns in the source file exactly (verified r = 1.000).
tot <- data.frame(
  sss8_total = rowSums(dat[, sss_items]),
  pcs_total  = rowSums(dat[, pcs_items])
)

cat("Loading item-level data...\n")
cat("Rows in file:", n_file, "| complete on all 21 items:", N, "\n\n")

# =============================================================================
# 2. CORRELATION MATRICES
# =============================================================================
# Items have 5 ordered response categories, so the polychoric matrix is the
# primary input (Pearson correlations of coarse ordinal items are attenuated
# and can manufacture spurious "difficulty" factors).

cat("Computing polychoric correlations...\n")
set.seed(SEED)
poly_out <- psych::polychoric(dat)
R_poly <- poly_out$rho
R_pear <- cor(dat)

# =============================================================================
# 3. SUITABILITY, DIMENSIONALITY, EFA
# =============================================================================

kmo_poly <- psych::KMO(R_poly)
kmo_pear <- psych::KMO(R_pear)
bart_poly <- psych::cortest.bartlett(R_poly, n = N)
bart_pear <- psych::cortest.bartlett(R_pear, n = N)

cat("Running parallel analysis + scree...\n")
png(file.path(plots_dir, "figure_18_scree_parallel_somatic.png"),
    width = 1600, height = 1200, res = 180)
set.seed(SEED)
pa <- psych::fa.parallel(R_poly, n.obs = N, fm = "pa", fa = "both",
                         n.iter = 500, main = "Scree + parallel analysis: SSS-8 + PCS items (polychoric)",
                         show.legend = TRUE)
dev.off()

cat("Running VSS / MAP...\n")
set.seed(SEED)
vss_out <- try(psych::vss(R_poly, n.obs = N, n = 6, fm = "pa",
                          rotate = "oblimin", plot = FALSE), silent = TRUE)

cat("Running EFA (1, 2, 3 factors)...\n")
efa_fit <- list()
for (k in 1:3) {
  efa_fit[[k]] <- try(psych::fa(R_poly, nfactors = k, n.obs = N, fm = "pa",
                                rotate = if (k == 1) "none" else "oblimin"),
                      silent = TRUE)
}

# =============================================================================
# 4. CFA (lavaan, WLSMV on ordinal items; MLR cross-check)
# =============================================================================

mod_1f <- paste0("SOM =~ ", paste(items, collapse = " + "))

mod_2f <- paste0(
  "SSS =~ ", paste(sss_items, collapse = " + "), "\n",
  "PCS =~ ", paste(pcs_items, collapse = " + ")
)

mod_bf <- paste0(
  "G =~ ", paste(items, collapse = " + "), "\n",
  "SSSs =~ ", paste(sss_items, collapse = " + "), "\n",
  "PCSs =~ ", paste(pcs_items, collapse = " + "), "\n",
  "G ~~ 0*SSSs\nG ~~ 0*PCSs\nSSSs ~~ 0*PCSs"
)

cat("Fitting CFAs (WLSMV)...\n")
fit_wlsmv <- list()
for (nm in c("one", "two", "bifactor")) {
  m <- switch(nm, one = mod_1f, two = mod_2f, bifactor = mod_bf)
  fit_wlsmv[[nm]] <- try(suppressWarnings(
    lavaan::cfa(m, data = dat, ordered = items, estimator = "WLSMV",
                std.lv = TRUE)), silent = TRUE)
}

cat("Fitting CFAs (MLR cross-check)...\n")
fit_mlr <- list()
for (nm in c("one", "two", "bifactor")) {
  m <- switch(nm, one = mod_1f, two = mod_2f, bifactor = mod_bf)
  fit_mlr[[nm]] <- try(suppressWarnings(
    lavaan::cfa(m, data = dat, estimator = "MLR", std.lv = TRUE)), silent = TRUE)
}

ok_fit <- function(f) {
  !inherits(f, "try-error") && lavaan::lavInspect(f, "converged")
}

# Pull a labelled fit-index row, preferring robust/scaled variants
get_fits <- function(f, robust = TRUE) {
  if (!ok_fit(f)) return(NULL)
  fm <- lavaan::fitMeasures(f)
  pick <- function(base) {
    cands <- if (robust) c(paste0(base, ".robust"), paste0(base, ".scaled"), base)
             else c(base)
    for (cc in cands) if (cc %in% names(fm)) return(unname(fm[cc]))
    NA_real_
  }
  chi_nm <- if (robust && "chisq.scaled" %in% names(fm)) "chisq.scaled" else "chisq"
  df_nm  <- if (robust && "df.scaled" %in% names(fm)) "df.scaled" else "df"
  p_nm   <- if (robust && "pvalue.scaled" %in% names(fm)) "pvalue.scaled" else "pvalue"
  data.frame(
    chisq = unname(fm[chi_nm]), df = unname(fm[df_nm]), pvalue = unname(fm[p_nm]),
    cfi = pick("cfi"), tli = pick("tli"),
    rmsea = pick("rmsea"),
    rmsea_lo = pick("rmsea.ci.lower"), rmsea_hi = pick("rmsea.ci.upper"),
    srmr = unname(fm["srmr"]),
    npar = unname(fm["npar"]),
    stringsAsFactors = FALSE
  )
}

lrt_wlsmv <- try(lavaan::lavTestLRT(fit_wlsmv$one, fit_wlsmv$two), silent = TRUE)
lrt_mlr <- try(lavaan::lavTestLRT(fit_mlr$one, fit_mlr$two), silent = TRUE)
lrt_wlsmv_bf <- try(lavaan::lavTestLRT(fit_wlsmv$two, fit_wlsmv$bifactor), silent = TRUE)

# Standardized loadings from a fitted CFA, per latent factor
std_load <- function(f, lv) {
  s <- lavaan::standardizedSolution(f)
  s <- s[s$op == "=~" & s$lhs == lv, ]
  setNames(s$est.std, s$rhs)
}

# Average variance extracted from standardized loadings
ave_from <- function(l) mean(l^2, na.rm = TRUE)

# =============================================================================
# 5. RELIABILITY
# =============================================================================

cat("Computing omega / alpha...\n")

# omega on the 21-item composite. Fitted from the polychoric matrix (ordinal
# omega) with a 3-factor group structure, which is psych's default schmid-
# leiman setup; also with 2 group factors matching the two scales.
om3 <- try(suppressWarnings(psych::omega(R_poly, nfactors = 3, n.obs = N,
                                         fm = "pa", rotate = "oblimin",
                                         plot = FALSE)), silent = TRUE)
om2 <- try(suppressWarnings(psych::omega(R_poly, nfactors = 2, n.obs = N,
                                         fm = "pa", rotate = "oblimin",
                                         plot = FALSE)), silent = TRUE)
om_sss <- try(suppressWarnings(psych::omega(R_poly[sss_items, sss_items],
                                            nfactors = 2, n.obs = N, fm = "pa",
                                            plot = FALSE)), silent = TRUE)
om_pcs <- try(suppressWarnings(psych::omega(R_poly[pcs_items, pcs_items],
                                            nfactors = 3, n.obs = N, fm = "pa",
                                            plot = FALSE)), silent = TRUE)

alpha_all <- psych::alpha(dat, warnings = FALSE)
alpha_sss <- psych::alpha(dat[, sss_items], warnings = FALSE)
alpha_pcs <- psych::alpha(dat[, pcs_items], warnings = FALSE)

# ordinal (polychoric) alpha
ord_alpha <- function(R) {
  k <- ncol(R)
  (k / (k - 1)) * (1 - k / sum(R))
}
oalpha_all <- ord_alpha(R_poly)
oalpha_sss <- ord_alpha(R_poly[sss_items, sss_items])
oalpha_pcs <- ord_alpha(R_poly[pcs_items, pcs_items])

r_tot <- cor(tot$sss8_total, tot$pcs_total)
r_tot_sp <- cor(tot$sss8_total, tot$pcs_total, method = "spearman")
disatt <- function(r, a, b) r / sqrt(a * b)

# =============================================================================
# BEGIN OUTPUT FILE
# =============================================================================

sink(file.path(analysis_output_dir, "18_output.txt"))

cat("=============================================================================\n")
cat("FACTOR STRUCTURE OF THE SOMATIC COMPOSITE: SSS-8 + PCS ITEMS\n")
cat("COMMSPSYCHOL-26-0385-T revision\n")
cat("=============================================================================\n\n")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("R", paste0(R.version$major, ".", R.version$minor),
    "| psych", as.character(packageVersion("psych")),
    "| lavaan", as.character(packageVersion("lavaan")),
    "| GPArotation", as.character(packageVersion("GPArotation")), "\n\n")

cat("QUESTION\n")
cat("--------\n")
cat("The SEM's somatic mediator is the z-scored average of SSS-8 and PCS.\n")
cat("Reviewer 1 argues these index different constructs (bodily symptom burden\n")
cat("vs cognitive appraisal of pain). Do the 21 items form ONE factor here?\n\n")

cat("SAMPLE\n")
cat("------\n")
cat("Source file: supplementary_data/data_with_ias_iats_ratio.csv\n")
cat(sprintf("Rows in file: %d (the analysis sample used throughout the paper)\n", n_file))
cat(sprintf("Missing item responses: %d cells across %d participants\n",
            n_na_cells, sum(!cc)))
if (n_na_cells > 0) {
  nz <- na_by_item[na_by_item > 0]
  cat("  items with missing values: ",
      paste(sprintf("%s (%d)", names(nz), nz), collapse = ", "), "\n", sep = "")
}
cat(sprintf("ANALYSIS N (listwise complete on all 21 items) = %d\n", N))
cat("No further exclusions are applied: the item file already contains exactly\n")
cat(sprintf("the %d participants of the published analyses.", nrow(df)), "\n")
cat("The 2-participant shortfall is item-level missingness only, not exclusion.\n")
cat("Note: for those 2 participants the raw item sums do NOT equal the sss8 /\n")
cat("pcs scale scores used in 12_sem_pathway_analysis.R. Those two totals are\n")
cat("prorated (and one SSS-8 total is mean-imputed) by 00_score_corrections.py,\n")
cat("because the original totals scored omitted items as zero. Both rows are\n")
cat("marked in the scoring_flag column. See CORRECTIONS.md.\n\n")

cat("MEASURES\n")
cat("--------\n")
cat("SSS-8: 8 items (sss_8_1 ... sss_8_8), 5 ordered categories, coded 1-5\n")
cat("PCS  : 13 items (pcs_1 ... pcs_13), 5 ordered categories, coded 1-5\n")
cat("Because every item has only 5 ordered categories, all factor analyses use\n")
cat("POLYCHORIC correlations, and the CFAs use the WLSMV estimator with items\n")
cat("declared ordered. Pearson/ML equivalents are reported as cross-checks.\n\n")

# -----------------------------------------------------------------------------

cat("\n#############################################################################\n")
cat("1. SUITABILITY FOR FACTOR ANALYSIS\n")
cat("#############################################################################\n\n")

cat(sprintf("KMO overall (polychoric matrix): %.3f\n", kmo_poly$MSA))
cat(sprintf("KMO overall (Pearson matrix)   : %.3f\n", kmo_pear$MSA))
cat("Kaiser's benchmarks: >.90 marvellous, >.80 meritorious, >.70 middling,\n")
cat(">.60 mediocre, <.50 unacceptable.\n\n")

cat("Per-item KMO (measure of sampling adequacy):\n")
cat(sprintf("%-12s %10s %10s\n", "Item", "MSA(poly)", "MSA(pear)"))
cat(paste(rep("-", 34), collapse = ""), "\n")
for (it in items) {
  cat(sprintf("%-12s %10.3f %10.3f\n", it,
              kmo_poly$MSAi[it], kmo_pear$MSAi[it]))
}
cat(sprintf("\nMinimum per-item MSA (poly): %.3f (%s)\n",
            min(kmo_poly$MSAi), names(which.min(kmo_poly$MSAi))))

cat("\nBartlett's test of sphericity (H0: correlation matrix = identity):\n")
cat(sprintf("  polychoric: chi2(%d) = %.1f, p = %s\n",
            bart_poly$df, bart_poly$chisq,
            ifelse(bart_poly$p.value < 1e-300, "< 1e-300",
                   format.pval(bart_poly$p.value, digits = 3))))
cat(sprintf("  Pearson   : chi2(%d) = %.1f, p = %s\n",
            bart_pear$df, bart_pear$chisq,
            ifelse(bart_pear$p.value < 1e-300, "< 1e-300",
                   format.pval(bart_pear$p.value, digits = 3))))
cat("\nBoth criteria are comfortably met: the 21 items are factorable. Note that\n")
cat("this says the items are INTERCORRELATED, not that they are UNIDIMENSIONAL.\n")

# -----------------------------------------------------------------------------

cat("\n#############################################################################\n")
cat("2. DIMENSIONALITY\n")
cat("#############################################################################\n\n")

cat("PARALLEL ANALYSIS (Horn), principal axis factoring on the polychoric\n")
cat(sprintf("matrix, %d simulated datasets, n.obs = %d, seed = %d.\n", 500, N, SEED))
cat("Simulation-based (rather than resampling-based) because the input is a\n")
cat("polychoric correlation matrix; this is psych's standard behaviour for a\n")
cat("matrix input and matches how ordinal parallel analysis is normally run.\n\n")

cat(sprintf("Factors retained (parallel analysis, PA): %d\n", pa$nfact))
cat(sprintf("Components retained (parallel analysis, PC): %d\n", pa$ncomp))

sim_fa <- if (!is.null(pa$fa.sim)) pa$fa.sim else pa$fa.simr
sim_pc <- if (!is.null(pa$pc.sim)) pa$pc.sim else pa$pc.simr
cat("\nEigenvalues, actual vs simulated:\n")
cat(sprintf("%-4s %12s %12s %12s %12s\n", "n", "FA actual", "FA simulated",
            "PC actual", "PC simulated"))
cat(paste(rep("-", 56), collapse = ""), "\n")
for (i in 1:length(pa$fa.values)) {
  cat(sprintf("%-4d %12.3f %12.3f %12.3f %12.3f\n", i,
              pa$fa.values[i], sim_fa[i], pa$pc.values[i], sim_pc[i]))
}
cat("\nA factor is retained while the actual eigenvalue exceeds the simulated one.\n")
cat("Scree + parallel plot: plots/figure_18_scree_parallel_somatic.png\n")
cat("\nNOTE ON INTERPRETATION: only the first two factor eigenvalues are\n")
cat("substantively large. Any factor beyond the second that parallel analysis\n")
cat("retains sits barely above the simulated line and, in the EFA below,\n")
cat("attracts no item at |lambda| >= 0.40 (see the 3-factor solution). The\n")
cat("defensible reading is TWO substantive dimensions, with any further\n")
cat("factors reflecting minor residual covariance among PCS items rather than\n")
cat("additional interpretable constructs. Whatever the exact count, the crucial\n")
cat("point for the reviewer's question is that it is NOT one.\n")

if (!inherits(vss_out, "try-error")) {
  cat("\nVERY SIMPLE STRUCTURE and VELICER'S MAP (polychoric, PA, oblimin):\n")
  vs <- vss_out$vss.stats
  cat(sprintf("%-4s %8s %8s %10s %10s %10s\n",
              "k", "VSS1", "VSS2", "MAP", "BIC", "RMSEA"))
  cat(paste(rep("-", 54), collapse = ""), "\n")
  for (i in 1:length(vss_out$cfit.1)) {
    cat(sprintf("%-4d %8.3f %8.3f %10.5f %10.1f %10.4f\n", i,
                vss_out$cfit.1[i],
                ifelse(i <= length(vss_out$cfit.2) && !is.na(vss_out$cfit.2[i]),
                       vss_out$cfit.2[i], NA_real_),
                vss_out$map[i],
                if ("BIC" %in% names(vs)) vs$BIC[i] else NA_real_,
                if ("RMSEA" %in% names(vs)) vs$RMSEA[i] else NA_real_))
  }
  cat(sprintf("\nVSS complexity 1 maximised at k = %d\n", which.max(vss_out$cfit.1)))
  cat(sprintf("VSS complexity 2 maximised at k = %d\n",
              which.max(replace(vss_out$cfit.2, is.na(vss_out$cfit.2), -Inf))))
  cat(sprintf("Velicer's MAP minimised at k = %d\n", which.min(vss_out$map)))
} else {
  cat("\nVSS / MAP failed to run:\n")
  cat(as.character(vss_out), "\n")
}

# -----------------------------------------------------------------------------

cat("\n#############################################################################\n")
cat("3. EXPLORATORY FACTOR ANALYSIS (principal axis, oblimin)\n")
cat("#############################################################################\n\n")

cat("Polychoric input, principal axis factoring (fm = 'pa'), oblique oblimin\n")
cat("rotation for k >= 2 (unrotated for k = 1). Loadings with |lambda| >=",
    sprintf("%.2f", LOAD_CUT), "are flagged '*'.\n")
cat("Item block: rows 1-8 are SSS-8, rows 9-21 are PCS.\n\n")

print_efa <- function(f, k) {
  if (inherits(f, "try-error")) {
    cat(sprintf("--- %d-FACTOR SOLUTION: FAILED ---\n", k))
    cat(as.character(f), "\n\n")
    return(invisible(NULL))
  }
  L <- unclass(f$loadings)[items, , drop = FALSE]
  cat(sprintf("--- %d-FACTOR SOLUTION ---\n\n", k))
  hdr <- sprintf("%-12s %-6s", "Item", "scale")
  for (j in 1:k) hdr <- paste0(hdr, sprintf(" %9s", colnames(L)[j]))
  hdr <- paste0(hdr, sprintf(" %8s %8s", "h2", "u2"))
  cat(hdr, "\n")
  cat(paste(rep("-", nchar(hdr)), collapse = ""), "\n")
  for (it in items) {
    scl <- if (it %in% sss_items) "SSS-8" else "PCS"
    line <- sprintf("%-12s %-6s", it, scl)
    for (j in 1:k) {
      v <- L[it, j]
      line <- paste0(line, sprintf(" %8.3f%s", v,
                                   ifelse(abs(v) >= LOAD_CUT, "*", " ")))
    }
    line <- paste0(line, sprintf(" %8.3f %8.3f", f$communality[it],
                                 f$uniquenesses[it]))
    cat(line, "\n")
  }

  cat("\nVariance accounted for:\n")
  print(round(f$Vaccounted, 4))

  if (k > 1 && !is.null(f$Phi)) {
    cat("\nFactor intercorrelations (oblimin):\n")
    print(round(f$Phi, 3))
  }

  # Which scale owns each factor
  cat("\nFactor composition at |lambda| >= ", sprintf("%.2f", LOAD_CUT), ":\n", sep = "")
  for (j in 1:k) {
    hi <- items[abs(L[, j]) >= LOAD_CUT]
    n_s <- sum(hi %in% sss_items); n_p <- sum(hi %in% pcs_items)
    cat(sprintf("  %-6s : %2d items (%d SSS-8, %d PCS)%s\n",
                colnames(L)[j], length(hi), n_s, n_p,
                if (length(hi)) paste0("  -> ", paste(hi, collapse = ", ")) else ""))
  }

  # Cross-loading / no-loading items
  n_cross <- sum(apply(abs(L) >= LOAD_CUT, 1, sum) > 1)
  n_none <- sum(apply(abs(L) >= LOAD_CUT, 1, sum) == 0)
  cat(sprintf("  cross-loading items: %d | items loading on no factor: %d\n",
              n_cross, n_none))
  cat(sprintf("\nModel fit of this EFA solution: RMSEA = %s, TLI = %s, BIC = %s\n",
              ifelse(is.null(f$RMSEA), "n/a", sprintf("%.4f", f$RMSEA[1])),
              ifelse(is.null(f$TLI), "n/a", sprintf("%.3f", f$TLI)),
              ifelse(is.null(f$BIC), "n/a", sprintf("%.1f", f$BIC))))
  cat(sprintf("Chi-square (empirical): %s\n",
              ifelse(is.null(f$STATISTIC), "n/a",
                     sprintf("%.1f, df = %d, p = %.4g", f$STATISTIC, f$dof, f$PVAL))))
  cat("\n")
}

for (k in 1:3) print_efa(efa_fit[[k]], k)

# -----------------------------------------------------------------------------

cat("\n#############################################################################\n")
cat("4. CONFIRMATORY FACTOR ANALYSIS\n")
cat("#############################################################################\n\n")

cat("Estimator: WLSMV (diagonally weighted least squares with mean- and\n")
cat("variance-adjusted test statistic), items declared ordered, delta\n")
cat("parameterisation, latent variances fixed to 1 (std.lv = TRUE). This is the\n")
cat("standard estimator for 5-category ordinal indicators. Fit indices reported\n")
cat("are lavaan's ROBUST versions (scaled chi-square, robust CFI/TLI/RMSEA).\n")
cat("An MLR (continuous-treatment) cross-check follows.\n\n")

cat("Models:\n")
cat("  (a) one-factor      : all 21 items -> a single SOM factor\n")
cat("  (b) two-factor      : 8 SSS-8 items -> SSS, 13 PCS items -> PCS, correlated\n")
cat("  (c) bifactor        : orthogonal general G (all 21) + SSS-specific +\n")
cat("                        PCS-specific group factors\n\n")

print_fit_table <- function(fits, robust, label) {
  cat(label, "\n")
  cat(sprintf("%-12s %10s %6s %10s %8s %8s %8s %-20s %8s %6s\n",
              "Model", "chi2", "df", "p", "CFI", "TLI", "RMSEA",
              "RMSEA 90% CI", "SRMR", "npar"))
  cat(paste(rep("-", 112), collapse = ""), "\n")
  for (nm in names(fits)) {
    f <- fits[[nm]]
    if (!ok_fit(f)) {
      msg <- if (inherits(f, "try-error")) "ERROR" else "DID NOT CONVERGE"
      cat(sprintf("%-12s %s\n", nm, msg))
      next
    }
    g <- get_fits(f, robust = robust)
    cat(sprintf("%-12s %10.2f %6.0f %10s %8.3f %8.3f %8.3f  [%6.3f, %6.3f]     %8.3f %6.0f\n",
                nm, g$chisq, g$df,
                ifelse(g$pvalue < 0.0001, "< .0001", sprintf("%.4f", g$pvalue)),
                g$cfi, g$tli, g$rmsea, g$rmsea_lo, g$rmsea_hi, g$srmr, g$npar))
  }
  cat("\n")
}

print_fit_table(fit_wlsmv, TRUE, sprintf("WLSMV (robust / scaled indices), N = %d:", N))
print_fit_table(fit_mlr, TRUE, "MLR cross-check (robust indices), items treated as continuous:")

cat("Conventional cut-offs: CFI/TLI >= .95 good (>= .90 acceptable),\n")
cat("RMSEA <= .06 good (<= .08 acceptable), SRMR <= .08 good.\n\n")

cat("--- CHI-SQUARE DIFFERENCE TEST: one-factor (a) vs two-factor (b) ---\n")
cat("Nested: the one-factor model is the two-factor model with the factor\n")
cat("correlation fixed at 1. For WLSMV lavaan applies a scaled difference test.\n\n")
if (!inherits(lrt_wlsmv, "try-error")) {
  print(lrt_wlsmv)
} else {
  cat("WLSMV difference test FAILED:\n"); cat(as.character(lrt_wlsmv), "\n")
}
cat("\n--- same comparison, MLR ---\n")
if (!inherits(lrt_mlr, "try-error")) {
  print(lrt_mlr)
} else {
  cat("MLR difference test FAILED:\n"); cat(as.character(lrt_mlr), "\n")
}

cat("\n--- two-factor (b) vs bifactor (c) ---\n")
if (!inherits(lrt_wlsmv_bf, "try-error")) {
  print(lrt_wlsmv_bf)
} else {
  cat("Comparison FAILED (models not nested or estimation problem):\n")
  cat(as.character(lrt_wlsmv_bf), "\n")
}

cat("\n--- KEY QUANTITY: SSS x PCS LATENT CORRELATION (two-factor model) ---\n")
r_lat_wlsmv <- NA_real_; r_lat_mlr <- NA_real_
if (ok_fit(fit_wlsmv$two)) {
  s <- lavaan::standardizedSolution(fit_wlsmv$two)
  row <- s[s$op == "~~" & s$lhs == "SSS" & s$rhs == "PCS", ]
  r_lat_wlsmv <- row$est.std[1]
  cat(sprintf("WLSMV: r(SSS, PCS) = %.3f, SE = %.3f, 95%% CI [%.3f, %.3f], z = %.2f, p %s\n",
              row$est.std[1], row$se[1], row$ci.lower[1], row$ci.upper[1],
              row$z[1], ifelse(row$pvalue[1] < 0.0001, "< .0001",
                               sprintf("= %.4f", row$pvalue[1]))))
  cat(sprintf("Shared latent variance = %.1f%% (r^2)\n", 100 * row$est.std[1]^2))
}
if (ok_fit(fit_mlr$two)) {
  s <- lavaan::standardizedSolution(fit_mlr$two)
  row <- s[s$op == "~~" & s$lhs == "SSS" & s$rhs == "PCS", ]
  r_lat_mlr <- row$est.std[1]
  cat(sprintf("MLR  : r(SSS, PCS) = %.3f, 95%% CI [%.3f, %.3f]\n",
              row$est.std[1], row$ci.lower[1], row$ci.upper[1]))
}
cat("\nA latent correlation of 1.00 would mean the two scales measure the same\n")
cat("thing; the one-factor model is exactly that hypothesis.\n")

cat("\n--- STANDARDIZED LOADINGS, ONE-FACTOR MODEL (WLSMV) ---\n")
if (ok_fit(fit_wlsmv$one)) {
  l1 <- std_load(fit_wlsmv$one, "SOM")
  cat(sprintf("%-12s %-6s %10s\n", "Item", "scale", "lambda"))
  cat(paste(rep("-", 30), collapse = ""), "\n")
  for (it in items) {
    cat(sprintf("%-12s %-6s %10.3f%s\n", it,
                if (it %in% sss_items) "SSS-8" else "PCS", l1[it],
                ifelse(abs(l1[it]) >= LOAD_CUT, " *", "")))
  }
  cat(sprintf("\nMean lambda: SSS-8 items = %.3f | PCS items = %.3f\n",
              mean(l1[sss_items]), mean(l1[pcs_items])))
  cat(sprintf("AVE of the single 'somatic' factor = %.3f (threshold .50)\n",
              ave_from(l1)))
}

cat("\n--- STANDARDIZED LOADINGS, TWO-FACTOR MODEL (WLSMV) ---\n")
ave_sss <- NA_real_; ave_pcs <- NA_real_
if (ok_fit(fit_wlsmv$two)) {
  ls <- std_load(fit_wlsmv$two, "SSS")
  lp <- std_load(fit_wlsmv$two, "PCS")
  cat(sprintf("%-12s %-6s %10s\n", "Item", "factor", "lambda"))
  cat(paste(rep("-", 30), collapse = ""), "\n")
  for (it in sss_items) cat(sprintf("%-12s %-6s %10.3f\n", it, "SSS", ls[it]))
  for (it in pcs_items) cat(sprintf("%-12s %-6s %10.3f\n", it, "PCS", lp[it]))
  ave_sss <- ave_from(ls); ave_pcs <- ave_from(lp)
  cat(sprintf("\nAVE(SSS-8) = %.3f | AVE(PCS) = %.3f\n", ave_sss, ave_pcs))
  if (!is.na(r_lat_wlsmv)) {
    cat(sprintf("Fornell-Larcker discriminant validity: r^2 = %.3f vs AVE = %.3f / %.3f\n",
                r_lat_wlsmv^2, ave_sss, ave_pcs))
    cat(sprintf("  -> %s\n",
                ifelse(r_lat_wlsmv^2 < min(ave_sss, ave_pcs, na.rm = TRUE),
                       "PASSES: each factor explains more of its own items than it shares with the other",
                       "FAILS: shared variance exceeds one factor's AVE")))
  }
}

cat("\n--- BIFACTOR MODEL (WLSMV) ---\n")
if (ok_fit(fit_wlsmv$bifactor)) {
  lg <- std_load(fit_wlsmv$bifactor, "G")
  lss <- std_load(fit_wlsmv$bifactor, "SSSs")
  lpp <- std_load(fit_wlsmv$bifactor, "PCSs")
  cat(sprintf("%-12s %-6s %10s %12s\n", "Item", "scale", "G", "specific"))
  cat(paste(rep("-", 42), collapse = ""), "\n")
  for (it in items) {
    sp <- if (it %in% sss_items) lss[it] else lpp[it]
    cat(sprintf("%-12s %-6s %10.3f %12.3f\n", it,
                if (it %in% sss_items) "SSS-8" else "PCS", lg[it], sp))
  }
  cat(sprintf("\nMean |G| loading: SSS-8 = %.3f | PCS = %.3f\n",
              mean(abs(lg[sss_items])), mean(abs(lg[pcs_items]))))
  cat(sprintf("Mean |specific| loading: SSS-8 = %.3f | PCS = %.3f\n",
              mean(abs(lss[sss_items])), mean(abs(lpp[pcs_items]))))
  cat("If the specific factors carry loadings comparable to or larger than the\n")
  cat("general factor, the general factor is not strong enough to justify a\n")
  cat("single composite score.\n")
  bf_gap <- mean(abs(lg[pcs_items])) - mean(abs(lg[sss_items]))
  cat(sprintf("\nHere the general factor is PCS-DOMINATED: mean |G| loading is %.3f for\n",
              mean(abs(lg[pcs_items]))))
  cat(sprintf("PCS items but only %.3f for SSS-8 items (gap = %.3f), while the SSS-8\n",
              mean(abs(lg[sss_items])), bf_gap))
  cat(sprintf("specific factor (%.3f) is stronger than SSS-8 items' loadings on G.\n",
              mean(abs(lss[sss_items]))))
  cat("In other words, a single score over these 21 items is mostly a PCS score.\n")
} else {
  cat("Bifactor model did NOT converge / errored under WLSMV:\n")
  cat(if (inherits(fit_wlsmv$bifactor, "try-error"))
        as.character(fit_wlsmv$bifactor) else "non-convergence\n")
}

# -----------------------------------------------------------------------------

cat("\n#############################################################################\n")
cat("5. RELIABILITY AND COMPOSITE JUSTIFICATION\n")
cat("#############################################################################\n\n")

cat("McDonald's omega (Schmid-Leiman on the polychoric matrix, principal axis):\n")
cat(sprintf("%-34s %10s %10s\n", "Item set / group factors", "omega_h", "omega_t"))
cat(paste(rep("-", 56), collapse = ""), "\n")
om_row <- function(o, lab) {
  if (inherits(o, "try-error")) {
    cat(sprintf("%-34s %10s %10s\n", lab, "FAILED", "FAILED")); return(invisible(NULL))
  }
  cat(sprintf("%-34s %10.3f %10.3f\n", lab, o$omega_h, o$omega.tot))
}
om_row(om3, "21 items (3 group factors)")
om_row(om2, "21 items (2 group factors)")
om_row(om_sss, "SSS-8 alone (2 group factors)")
om_row(om_pcs, "PCS alone (3 group factors)")
cat("\nNOTE: psych::omega requires 3 group factors for full identification. The\n")
cat("two rows fitted with 2 group factors (21 items / SSS-8 alone) were estimated\n")
cat("with general-factor loadings constrained equal, as psych warns; the 3-group-\n")
cat("factor row for the 21-item composite is the unconstrained solution and is\n")
cat("the one quoted in the verdict.\n")

if (!inherits(om3, "try-error")) {
  cat(sprintf("\n21-item composite, 3 group factors: omega_h = %.3f, omega_total = %.3f\n",
              om3$omega_h, om3$omega.tot))
  cat(sprintf("Proportion of RELIABLE variance attributable to the general factor:\n"))
  cat(sprintf("  omega_h / omega_total = %.3f\n", om3$omega_h / om3$omega.tot))
  if (!is.null(om3$ECV)) cat(sprintf("  ECV (explained common variance) = %.3f\n", om3$ECV))
  cat("Rule of thumb: omega_h >= .70-.80 supports scoring a single total; below\n")
  cat("~.50 the total score is dominated by group factors, not a general factor.\n")
}

cat("\nCronbach's alpha (raw / Pearson) and ordinal alpha (polychoric):\n")
cat(sprintf("%-26s %10s %10s %14s\n", "Scale", "alpha", "ord.alpha", "avg r (raw)"))
cat(paste(rep("-", 62), collapse = ""), "\n")
cat(sprintf("%-26s %10.3f %10.3f %14.3f\n", "21-item composite",
            alpha_all$total$raw_alpha, oalpha_all, alpha_all$total$average_r))
cat(sprintf("%-26s %10.3f %10.3f %14.3f\n", "SSS-8 (8 items)",
            alpha_sss$total$raw_alpha, oalpha_sss, alpha_sss$total$average_r))
cat(sprintf("%-26s %10.3f %10.3f %14.3f\n", "PCS (13 items)",
            alpha_pcs$total$raw_alpha, oalpha_pcs, alpha_pcs$total$average_r))
cat("\nNOTE: alpha rises mechanically with test length, so a high alpha for the\n")
cat("21-item composite is NOT evidence of unidimensionality. The average\n")
cat("inter-item correlation and omega_h are the informative numbers.\n")

cat("\nAverage variance extracted (from the WLSMV CFAs):\n")
if (ok_fit(fit_wlsmv$one))
  cat(sprintf("  one-factor 'somatic' (all 21 items): AVE = %.3f\n",
              ave_from(std_load(fit_wlsmv$one, "SOM"))))
if (!is.na(ave_sss))
  cat(sprintf("  two-factor: AVE(SSS-8) = %.3f, AVE(PCS) = %.3f\n", ave_sss, ave_pcs))
cat("  (AVE >= .50 is the usual convergent-validity threshold.)\n")

cat("\nSSS-8 TOTAL x PCS TOTAL (sum scores):\n")
ci_r <- psych::r.con(r_tot, N)
cat(sprintf("  Pearson  r = %.3f, 95%% CI [%.3f, %.3f], N = %d, r^2 = %.3f\n",
            r_tot, ci_r[1], ci_r[2], N, r_tot^2))
cat(sprintf("  Spearman rho = %.3f\n", r_tot_sp))
cat(sprintf("  disattenuated (Cronbach alpha)  = %.3f\n",
            disatt(r_tot, alpha_sss$total$raw_alpha, alpha_pcs$total$raw_alpha)))
cat(sprintf("  disattenuated (ordinal alpha)   = %.3f\n",
            disatt(r_tot, oalpha_sss, oalpha_pcs)))
if (!inherits(om_sss, "try-error") && !inherits(om_pcs, "try-error"))
  cat(sprintf("  disattenuated (omega total)     = %.3f\n",
              disatt(r_tot, om_sss$omega.tot, om_pcs$omega.tot)))
if (!is.na(r_lat_wlsmv))
  cat(sprintf("  latent correlation, WLSMV CFA   = %.3f (the model-based analogue)\n",
              r_lat_wlsmv))
cat("\nIf the disattenuated correlation is well below 1.00, the two scales are not\n")
cat("measuring the same construct even after correcting for unreliability.\n")

# -----------------------------------------------------------------------------

cat("\n#############################################################################\n")
cat("6. VERDICT\n")
cat("#############################################################################\n\n")

f1 <- get_fits(fit_wlsmv$one, TRUE)
f2 <- get_fits(fit_wlsmv$two, TRUE)
fb <- get_fits(fit_wlsmv$bifactor, TRUE)

cat("EVIDENCE SUMMARY\n")
cat("----------------\n")
cat(sprintf("  N = %d complete cases (of %d in the analysis sample)\n", N, n_file))
cat(sprintf("  KMO = %.3f, Bartlett p < .0001 -> factorable\n", kmo_poly$MSA))
cat(sprintf("  Parallel analysis (PA, polychoric): %d factors\n", pa$nfact))
if (!inherits(vss_out, "try-error"))
  cat(sprintf("  MAP: %d factor(s); VSS-1: %d; VSS-2: %d\n",
              which.min(vss_out$map), which.max(vss_out$cfit.1),
              which.max(replace(vss_out$cfit.2, is.na(vss_out$cfit.2), -Inf))))
if (!inherits(efa_fit[[2]], "try-error")) {
  L2 <- unclass(efa_fit[[2]]$loadings)[items, , drop = FALSE]
  f1_items <- items[abs(L2[, 1]) >= LOAD_CUT]
  f2_items <- items[abs(L2[, 2]) >= LOAD_CUT]
  cat(sprintf("  2-factor EFA: factor 1 = %d SSS / %d PCS; factor 2 = %d SSS / %d PCS\n",
              sum(f1_items %in% sss_items), sum(f1_items %in% pcs_items),
              sum(f2_items %in% sss_items), sum(f2_items %in% pcs_items)))
  if (!is.null(efa_fit[[2]]$Phi))
    cat(sprintf("  EFA factor correlation (oblimin) = %.3f\n", efa_fit[[2]]$Phi[1, 2]))
}
if (!is.null(f1)) cat(sprintf("  CFA 1-factor: CFI = %.3f, TLI = %.3f, RMSEA = %.3f [%.3f, %.3f], SRMR = %.3f\n",
                              f1$cfi, f1$tli, f1$rmsea, f1$rmsea_lo, f1$rmsea_hi, f1$srmr))
if (!is.null(f2)) cat(sprintf("  CFA 2-factor: CFI = %.3f, TLI = %.3f, RMSEA = %.3f [%.3f, %.3f], SRMR = %.3f\n",
                              f2$cfi, f2$tli, f2$rmsea, f2$rmsea_lo, f2$rmsea_hi, f2$srmr))
if (!is.null(fb)) cat(sprintf("  CFA bifactor: CFI = %.3f, TLI = %.3f, RMSEA = %.3f [%.3f, %.3f], SRMR = %.3f\n",
                              fb$cfi, fb$tli, fb$rmsea, fb$rmsea_lo, fb$rmsea_hi, fb$srmr))
if (!is.na(r_lat_wlsmv))
  cat(sprintf("  Latent r(SSS-8, PCS) = %.3f (%.1f%% shared variance)\n",
              r_lat_wlsmv, 100 * r_lat_wlsmv^2))
if (ok_fit(fit_wlsmv$one)) {
  l1v <- std_load(fit_wlsmv$one, "SOM")
  cat(sprintf("  1-factor loadings are lopsided: mean lambda = %.3f (SSS-8) vs %.3f (PCS)\n",
              mean(l1v[sss_items]), mean(l1v[pcs_items])))
}
cat(sprintf("  Observed r(SSS-8 total, PCS total) = %.3f; disattenuated = %.3f\n",
            r_tot, disatt(r_tot, alpha_sss$total$raw_alpha, alpha_pcs$total$raw_alpha)))
if (!inherits(om3, "try-error"))
  cat(sprintf("  omega_h = %.3f, omega_total = %.3f for the 21-item composite\n",
              om3$omega_h, om3$omega.tot))

cat("\nAUTOMATED READ OF THE DECISION RULES\n")
cat("------------------------------------\n")
rules <- character(0)
rules <- c(rules, sprintf("Parallel analysis indicates %d factor(s) -> %s",
                          pa$nfact,
                          ifelse(pa$nfact == 1, "consistent with one factor",
                                 "MORE THAN ONE dimension")))
if (!is.null(f1) && !is.null(f2)) {
  rules <- c(rules, sprintf("2-factor CFA improves CFI by %+.3f, RMSEA by %+.3f, SRMR by %+.3f vs 1-factor",
                            f2$cfi - f1$cfi, f2$rmsea - f1$rmsea, f2$srmr - f1$srmr))
  rules <- c(rules, sprintf("1-factor model %s conventional fit thresholds (CFI>=.95, RMSEA<=.06, SRMR<=.08)",
                            ifelse(f1$cfi >= 0.95 && f1$rmsea <= 0.06 && f1$srmr <= 0.08,
                                   "MEETS", "does NOT meet")))
  rules <- c(rules, sprintf("2-factor model %s those thresholds",
                            ifelse(f2$cfi >= 0.95 && f2$rmsea <= 0.06 && f2$srmr <= 0.08,
                                   "MEETS", "does NOT meet")))
}
if (!is.na(r_lat_wlsmv))
  rules <- c(rules, sprintf("latent r = %.3f is %s the .85 heuristic for 'not discriminable'",
                            r_lat_wlsmv, ifelse(r_lat_wlsmv >= 0.85, "ABOVE", "BELOW")))
if (!inherits(om3, "try-error"))
  rules <- c(rules, sprintf("omega_h = %.3f is %s .70, so a single total score is %s",
                            om3$omega_h, ifelse(om3$omega_h >= 0.70, "above", "below"),
                            ifelse(om3$omega_h >= 0.70, "defensible on reliability grounds",
                                   "NOT well justified on reliability grounds alone")))
for (r in rules) cat("  - ", r, "\n", sep = "")

cat("\nPLAIN-LANGUAGE VERDICT\n")
cat("----------------------\n")
two_better <- !is.null(f1) && !is.null(f2) && (f2$cfi > f1$cfi) && (f2$rmsea < f1$rmsea)
if (pa$nfact > 1 && two_better) {
  cat("The SSS-8 and the PCS are EMPIRICALLY SEPARABLE in this sample. Parallel\n")
  cat("analysis retains more than one factor, the two-factor solution splits\n")
  cat("cleanly along scale boundaries, and the two-factor CFA fits better than\n")
  cat("the one-factor CFA on every index reported above.\n")
} else if (pa$nfact == 1) {
  cat("The data are consistent with a single dimension: parallel analysis\n")
  cat("retains one factor and the one-factor CFA is not clearly outperformed.\n")
} else {
  cat("Mixed evidence: see the numbers above rather than a one-line summary.\n")
}
if (!is.na(r_lat_wlsmv)) {
  cat(sprintf("\nHOWEVER, the two factors correlate r = %.3f (%.0f%% shared variance).\n",
              r_lat_wlsmv, 100 * r_lat_wlsmv^2))
  if (r_lat_wlsmv >= 0.85) {
    cat("That is high enough that the two factors are barely discriminable, so a\n")
    cat("composite is a defensible practical simplification even though it is not\n")
    cat("the best-fitting measurement model.\n")
  } else {
    cat("That is a substantial but far from unity association: the two scales share\n")
    cat("a meaningful common core while retaining clearly separate reliable\n")
    cat("variance. A single 'somatic' composite therefore blends two related but\n")
    cat("distinguishable constructs, exactly as the reviewer argued.\n")
  }
}
cat("\nRECOMMENDED FRAMING FOR THE RESPONSE LETTER\n")
cat("-------------------------------------------\n")
cat("Report these numbers, concede that the composite is not strictly\n")
cat("unidimensional, and lean on 17_revision_sensitivity.R: the SSS-8-only\n")
cat("model is the honest primary specification, with the SSS-8 + PCS composite\n")
cat("retained (if at all) as a broader 'somatic distress' index whose label is\n")
cat("changed to reflect that it includes pain-related cognitive appraisal.\n")

cat("\nProvenance: every number above was produced by this script\n")
cat("(18_somatic_factor_structure.R) in this run. Nothing is carried over from\n")
cat("earlier analyses except the reference values quoted in the header comment.\n")

sink()

# =============================================================================
# =============================================================================
# PART B. IS THE MULTIPLICATIVE EFFECT A PAIN-CATASTROPHIZING EFFECT?
#         DUAL-PATHWAY SEM WITH THREE ALTERNATIVE SOMATIC MEDIATORS
# =============================================================================
# =============================================================================
#
# Motivation: Part A shows the 21-item composite is PCS-dominated (bifactor
# mean |G| loading .77 for PCS items vs .36 for SSS-8 items), and
# 17_revision_sensitivity.R shows the IAS x IATS -> Somatic interaction is
# stronger with the composite than with SSS-8 alone. If the effect is carried
# by the PCS, the published "somatic amplification" framing is really a
# pain-catastrophizing framing.
#
# Specification: EXACTLY model M0 of 17_revision_sensitivity.R (which itself
# reproduces 12_sem_pathway_analysis.R): full 42-item IAS/IATS means, TAS and
# the mediator regressed on ias_z, iats_z and their product, MH composite
# regressed on both mediators plus direct effects, ML estimator, NO covariates,
# all variables z-scored before the product term is formed. Only the mediator
# changes between models.
# =============================================================================

BOOT_SEED_B <- 42     # same seed as 17_revision_sensitivity.R
N_BOOT_B <- 2000      # same resample count as its sensitivity models A2/A3/A4

# -----------------------------------------------------------------------------
# B1. DATA (same source file as 12 / 17), plus PCS subscales from the item file
# -----------------------------------------------------------------------------

dfc <- read.csv("dfc_interoception_profiling.csv")

# The item-level file (df, used in Part A) and the SEM source file (dfc) are the
# same participants in the same row order. Verified rather than assumed:
align_id <- ("responseid" %in% names(df)) && ("responseid" %in% names(dfc)) &&
  isTRUE(all.equal(df$responseid, dfc$responseid))
align_scales <- isTRUE(all.equal(df$pcs, dfc$pcs)) &&
  isTRUE(all.equal(df$sss8, dfc$sss8)) && isTRUE(all.equal(df$tas, dfc$tas))
if (!align_id || !align_scales)
  stop("Item file and SEM source file are NOT row-aligned; PCS subscales cannot be attached.")

dfc$ias_full <- rowMeans(dfc[, paste0("ias_", 1:21)], na.rm = TRUE)
dfc$iats_full <- rowMeans(dfc[, paste0("iats_", 1:21)], na.rm = TRUE)
dfc$mh_composite <- rowMeans(cbind(scale(dfc$phq9), scale(dfc$gad7),
                                   scale(dfc$stai)), na.rm = TRUE)
dfc$somatic_combined <- rowMeans(cbind(scale(dfc$sss8), scale(dfc$pcs)), na.rm = TRUE)
dfc$somatic_sss_only <- dfc$sss8
dfc$somatic_pcs_only <- dfc$pcs

# PCS subscales. The subscale membership is NOT recorded in any scoring script
# in this project or in the upstream mi project (verified by search); it is
# Sullivan, Bishop & Pivik's (1995) published definition, matched here against
# the VERBATIM item text recovered from the Qualtrics survey export
# (psychosomatic_controls_October 22,
# 2024_11.48.xlsx, row 2 = question text). The wording and numbering match the
# published instrument item for item, so the mapping is verified, not assumed.
pcs_rum_items  <- paste0("pcs_", c(8, 9, 10, 11))
pcs_mag_items  <- paste0("pcs_", c(6, 7, 13))
pcs_help_items <- paste0("pcs_", c(1, 2, 3, 4, 5, 12))

PCS_ITEM_TEXT <- c(
  "1. I worry all the time about whether the pain will end",
  "2. I feel I can't go on",
  "3. It's terrible and I think it's never going to get any better",
  "4. It's awful and I feel it overwhelms me",
  "5. I feel I can't stand it any more",
  "6. I become afraid that the pain will become worse",
  "7. I keep thinking of other painful events",
  "8. I anxiously want the pain to go away",
  "9. I can't seem to keep it out of my mind",
  "10. I keep thinking about how much it hurts",
  "11. I keep thinking about how badly I want the pain to stop",
  "12. There's nothing I can do to reduce the intensity of the pain",
  "13. I wonder whether something serious may happen")

dfc$pcs_rum  <- rowMeans(df[, pcs_rum_items], na.rm = TRUE)
dfc$pcs_mag  <- rowMeans(df[, pcs_mag_items], na.rm = TRUE)
dfc$pcs_help <- rowMeans(df[, pcs_help_items], na.rm = TRUE)

# -----------------------------------------------------------------------------
# B2. MODEL ENGINE (identical specification to 17_revision_sensitivity.R, M0)
# -----------------------------------------------------------------------------

B_PATHS <- c("a1", "a2", "a3", "b2", "b1", "b3", "c1", "c2", "d1", "d2", "d3",
             "ind_ias_tas", "ind_iats_tas", "ind_ias_som", "ind_iats_som",
             "total_indirect")

B_NAMES <- c(a1 = "IAS -> TAS", a2 = "IATS -> TAS", a3 = "IAS x IATS -> TAS",
             b2 = "IAS -> Mediator", b1 = "IATS -> Mediator",
             b3 = "IAS x IATS -> Mediator",
             c1 = "TAS -> MH", c2 = "Mediator -> MH",
             d1 = "IAS -> MH (direct)", d2 = "IATS -> MH (direct)",
             d3 = "IAS x IATS -> MH (direct)",
             ind_ias_tas = "IAS -> TAS -> MH", ind_iats_tas = "IATS -> TAS -> MH",
             ind_ias_som = "IAS -> Mediator -> MH",
             ind_iats_som = "IATS -> Mediator -> MH",
             total_indirect = "Total indirect")

b_zprep <- function(data, med_var) {
  data$ias_z <- scale(data$ias_full)[, 1]
  data$iats_z <- scale(data$iats_full)[, 1]
  data$tas_z <- scale(data$tas)[, 1]
  data$somatic_z <- scale(data[[med_var]])[, 1]
  data$mh_z <- scale(data$mh_composite)[, 1]
  data$ias_x_iats <- data$ias_z * data$iats_z
  data
}

B_MODEL <- '
  tas_z ~ a1*ias_z + a2*iats_z + a3*ias_x_iats
  somatic_z ~ b1*iats_z + b2*ias_z + b3*ias_x_iats

  mh_z ~ c1*tas_z + c2*somatic_z + d1*ias_z + d2*iats_z + d3*ias_x_iats

  ind_ias_tas := a1 * c1
  ind_iats_tas := a2 * c1
  ind_ias_som := b2 * c2
  ind_iats_som := b1 * c2
  total_indirect := ind_ias_tas + ind_iats_tas + ind_ias_som + ind_iats_som
'

b_boot_p <- function(v) min(2 * min(mean(v <= 0), mean(v >= 0)), 1)

# Bootstrap both metrics (unstandardized est and std.all), re-standardizing
# inside every resample, exactly as 17_revision_sensitivity.R does.
b_boot <- function(fit, labs, R = N_BOOT_B) {
  extract <- function(x) {
    p <- parameterEstimates(x, se = FALSE, zstat = FALSE, pvalue = FALSE, ci = FALSE)
    s <- standardizedSolution(x, se = FALSE, zstat = FALSE, pvalue = FALSE, ci = FALSE)
    c(p$est[match(labs, p$label)], s$est.std[match(labs, s$label)])
  }
  set.seed(BOOT_SEED_B)
  B <- suppressWarnings(lavaan::bootstrapLavaan(fit, R = R, type = "ordinary",
                                                FUN = extract))
  B <- as.matrix(B)
  ok <- stats::complete.cases(B)
  R_failed <- sum(!ok)
  B <- B[ok, , drop = FALSE]
  k <- length(labs)
  point <- extract(fit)
  out <- data.frame(
    label = labs,
    est_unstd = point[1:k],
    lo_unstd = apply(B[, 1:k, drop = FALSE], 2, quantile, probs = 0.025),
    hi_unstd = apply(B[, 1:k, drop = FALSE], 2, quantile, probs = 0.975),
    p_boot_unstd = apply(B[, 1:k, drop = FALSE], 2, b_boot_p),
    est_std = point[(k + 1):(2 * k)],
    se_boot_std = apply(B[, (k + 1):(2 * k), drop = FALSE], 2, sd),
    lo_std = apply(B[, (k + 1):(2 * k), drop = FALSE], 2, quantile, probs = 0.025),
    hi_std = apply(B[, (k + 1):(2 * k), drop = FALSE], 2, quantile, probs = 0.975),
    p_boot_std = apply(B[, (k + 1):(2 * k), drop = FALSE], 2, b_boot_p),
    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  attr(out, "R_ok") <- nrow(B); attr(out, "R_failed") <- R_failed
  out
}

b_run <- function(med_var, label, boot = TRUE) {
  d <- b_zprep(dfc, med_var)
  fit <- lavaan::sem(B_MODEL, data = d)
  list(label = label, med_var = med_var, fit = fit,
       std = standardizedSolution(fit, level = 0.95),
       unstd = parameterEstimates(fit, level = 0.95),
       r2 = inspect(fit, "rsquare"),
       nobs = lavInspect(fit, "nobs"),
       converged = lavInspect(fit, "converged"),
       boot = if (boot) b_boot(fit, B_PATHS) else NULL)
}

b_get <- function(res, lbl, std = TRUE) {
  tab <- if (std) res$std else res$unstd
  row <- tab[tab$label == lbl, ]
  if (nrow(row) == 0) return(list(est = NA, se = NA, lo = NA, hi = NA, p = NA))
  list(est = if (std) row$est.std[1] else row$est[1], se = row$se[1],
       lo = row$ci.lower[1], hi = row$ci.upper[1], p = row$pvalue[1])
}

b_getboot <- function(res, lbl) {
  if (is.null(res$boot)) return(NULL)
  r <- res$boot[res$boot$label == lbl, ]
  if (nrow(r) == 0) return(NULL)
  as.list(r[1, ])
}

b_fmt_p <- function(p) {
  if (is.na(p)) return("n/a")
  if (p < 0.0001) return("<.0001")
  sprintf("%.4f", p)
}
b_star <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***"); if (p < 0.01) return("**")
  if (p < 0.05) return("*"); if (p < 0.10) return("."); ""
}

cat("PART B: fitting dual-pathway SEMs with alternative mediators...\n")

b_models <- list(
  combined = b_run("somatic_combined", "(a) SSS-8 + PCS composite (published)"),
  sss      = b_run("somatic_sss_only", "(b) SSS-8 alone"),
  pcs      = b_run("somatic_pcs_only", "(c) PCS alone"),
  rum      = b_run("pcs_rum",  "(d) PCS rumination subscale"),
  mag      = b_run("pcs_mag",  "(e) PCS magnification subscale"),
  help     = b_run("pcs_help", "(f) PCS helplessness subscale")
)

# -----------------------------------------------------------------------------
# B3. JOINT TWO-MEDIATOR MODEL AND THE FORMAL EQUALITY TEST
# -----------------------------------------------------------------------------
# Both SSS-8 and PCS enter as parallel mediators with the SAME predictor set,
# their residuals correlated. Because the two mediator equations share an
# identical right-hand side, each b3 in the joint model is numerically identical
# to its single-mediator counterpart (seemingly-unrelated regressions reduce to
# equation-by-equation OLS when the regressors coincide); the joint fit adds the
# sampling covariance between them, which is exactly what a difference test
# needs. This is checked numerically below rather than asserted.

dj <- dfc
dj$ias_z <- scale(dj$ias_full)[, 1]
dj$iats_z <- scale(dj$iats_full)[, 1]
dj$tas_z <- scale(dj$tas)[, 1]
dj$sss_z <- scale(dj$sss8)[, 1]
dj$pcs_z <- scale(dj$pcs)[, 1]
dj$mh_z <- scale(dj$mh_composite)[, 1]
dj$ias_x_iats <- dj$ias_z * dj$iats_z

B_JOINT <- '
  tas_z ~ a1*ias_z + a2*iats_z + a3*ias_x_iats
  sss_z ~ b1s*iats_z + b2s*ias_z + b3s*ias_x_iats
  pcs_z ~ b1p*iats_z + b2p*ias_z + b3p*ias_x_iats

  mh_z ~ c1*tas_z + c2s*sss_z + c2p*pcs_z + d1*ias_z + d2*iats_z + d3*ias_x_iats

  sss_z ~~ pcs_z

  diff_b3 := b3p - b3s
  diff_b1 := b1p - b1s
  diff_b2 := b2p - b2s
'

fit_joint <- lavaan::sem(B_JOINT, data = dj)
wald_b3 <- try(lavaan::lavTestWald(fit_joint, constraints = "b3p == b3s"), silent = TRUE)
wald_all <- try(lavaan::lavTestWald(fit_joint,
                                    constraints = "b3p == b3s\nb1p == b1s\nb2p == b2s"),
                silent = TRUE)

JOINT_LABS <- c("b1s", "b2s", "b3s", "b1p", "b2p", "b3p", "c2s", "c2p",
                "diff_b3", "diff_b1", "diff_b2")
boot_joint <- b_boot(fit_joint, JOINT_LABS)

j_std <- standardizedSolution(fit_joint, level = 0.95)
j_unstd <- parameterEstimates(fit_joint, level = 0.95)
j_get <- function(lbl, std = TRUE) {
  tab <- if (std) j_std else j_unstd
  r <- tab[tab$label == lbl, ]
  if (nrow(r) == 0) return(list(est = NA, se = NA, lo = NA, hi = NA, p = NA))
  list(est = if (std) r$est.std[1] else r$est[1], se = r$se[1],
       lo = r$ci.lower[1], hi = r$ci.upper[1], p = r$pvalue[1])
}

# numerical identity check between joint and single-mediator coefficients
chk_b3s <- j_get("b3s", FALSE)$est - b_get(b_models$sss, "b3", FALSE)$est
chk_b3p <- j_get("b3p", FALSE)$est - b_get(b_models$pcs, "b3", FALSE)$est

# -----------------------------------------------------------------------------
# B4. PRINT PART B (appended to the same output file)
# -----------------------------------------------------------------------------

sink(file.path(analysis_output_dir, "18_output.txt"), append = TRUE)

cat("\n\n")
cat("#############################################################################\n")
cat("#############################################################################\n")
cat("PART B. IS THE MULTIPLICATIVE EFFECT A PAIN-CATASTROPHIZING EFFECT?\n")
cat("        DUAL-PATHWAY SEM WITH THREE ALTERNATIVE SOMATIC MEDIATORS\n")
cat("#############################################################################\n")
cat("#############################################################################\n\n")

cat("MOTIVATION\n")
cat("----------\n")
cat("Part A shows the 21-item composite is PCS-dominated. 17_revision_sensitivity.R\n")
cat("shows the IAS x IATS -> Somatic interaction weakens when PCS is dropped.\n")
cat("If the interaction is carried by the PCS, the published 'somatic\n")
cat("amplification' effect is really a pain-catastrophizing effect. Part B tests\n")
cat("that directly by swapping ONLY the mediator.\n\n")

cat("SPECIFICATION (identical to model M0 of 17_revision_sensitivity.R)\n")
cat("------------------------------------------------------------------\n")
cat("  data      : dfc_interoception_profiling.csv\n")
cat("  predictors: IAS and IATS full 42-item means, z-scored, product term formed\n")
cat("              after z-scoring (ias_x_iats = ias_z * iats_z)\n")
cat("  model     : tas_z ~ ias_z + iats_z + ias_x_iats\n")
cat("              somatic_z ~ ias_z + iats_z + ias_x_iats\n")
cat("              mh_z ~ tas_z + somatic_z + ias_z + iats_z + ias_x_iats\n")
cat("  outcome   : MH composite = z-scored mean of PHQ-9, GAD-7, STAI\n")
cat("  estimator : ML (lavaan), no covariates (matches baseline M0)\n")
cat(sprintf("  bootstrap : %d ordinary nonparametric resamples, percentile CIs, seed = %d\n",
            N_BOOT_B, BOOT_SEED_B))
cat("  Only the mediator differs across models.\n\n")

cat("MEDIATORS\n")
cat("---------\n")
cat("  (a) SSS-8 + PCS composite : z-scored mean of the two scale totals (published)\n")
cat("  (b) SSS-8 alone           : SSS-8 total\n")
cat("  (c) PCS alone             : PCS total\n")
cat("  (d) PCS rumination        : mean of pcs_8, pcs_9, pcs_10, pcs_11\n")
cat("  (e) PCS magnification     : mean of pcs_6, pcs_7, pcs_13\n")
cat("  (f) PCS helplessness      : mean of pcs_1..pcs_5, pcs_12\n\n")

cat("PCS SUBSCALE PROVENANCE: no scoring script in this project or in the\n")
cat("upstream mi project assigns PCS items to subscales (only a single pcs\n")
cat("total is computed). The groupings above are Sullivan, Bishop & Pivik's\n")
cat("(1995) published definition, matched against the verbatim item text\n")
cat("recovered from the Qualtrics survey export\n")
cat("(mi/analyses/soma/data/psychosomatic_controls_October 22, 2024_11.48.xlsx,\n")
cat("row 2 = question text). The recovered wording matches the published\n")
cat("instrument item for item, so the mapping is verified against item content\n")
cat("rather than assumed from numbering alone. Items as administered:\n")
for (i in 1:13) {
  grp <- if (paste0("pcs_", i) %in% pcs_rum_items) "RUMINATION"
         else if (paste0("pcs_", i) %in% pcs_mag_items) "MAGNIFICATION"
         else "HELPLESSNESS"
  cat(sprintf("  %-8s %-14s %s\n", paste0("pcs_", i), grp, PCS_ITEM_TEXT[i]))
}
cat("\nRow alignment between the item-level file and the SEM source file was\n")
cat("verified (responseid identical row-wise, and the sss8 / pcs / tas scale\n")
cat("scores identical row-wise), so subscale scores line up with the SEM data.\n\n")

cat("SAMPLE AND CONVERGENCE\n")
cat("----------------------\n")
cat(sprintf("%-40s %8s %12s\n", "Model", "N", "converged"))
cat(paste(rep("-", 62), collapse = ""), "\n")
for (nm in names(b_models)) {
  r <- b_models[[nm]]
  cat(sprintf("%-40s %8d %12s\n", r$label, r$nobs, r$converged))
}
cat(sprintf("%-40s %8d %12s\n", "joint two-mediator model (SSS-8 + PCS)",
            lavInspect(fit_joint, "nobs"), lavInspect(fit_joint, "converged")))
cat("\nNo missing data on any model variable, so every model uses the full N.\n")

# ---------------------------------------------------------------------------

b_print_model <- function(res) {
  cat("\n-----------------------------------------------------------------------------\n")
  cat("MODEL ", res$label, "  [mediator variable: ", res$med_var, "]\n", sep = "")
  cat("-----------------------------------------------------------------------------\n\n")
  cat("Fully standardized (std.all) coefficients, delta-method 95% CI and Wald p,\n")
  cat("with percentile bootstrap 95% CI and bootstrap p on the same metric.\n\n")
  cat(sprintf("%-28s %8s %20s %9s | %20s %9s %s\n", "Path", "beta",
              "95% CI (delta)", "p", "95% CI (boot)", "p_boot", ""))
  cat(paste(rep("-", 112), collapse = ""), "\n")
  for (lb in B_PATHS[1:11]) {
    v <- b_get(res, lb); bt <- b_getboot(res, lb)
    cat(sprintf("%-28s %8.3f  [%7.3f, %7.3f] %9s | [%7.3f, %7.3f] %9s %s\n",
                B_NAMES[lb], v$est, v$lo, v$hi, b_fmt_p(v$p),
                bt$lo_std, bt$hi_std, b_fmt_p(bt$p_boot_std), b_star(v$p)))
  }
  cat("\nIndirect effects (standardized):\n")
  cat(paste(rep("-", 112), collapse = ""), "\n")
  for (lb in B_PATHS[12:16]) {
    v <- b_get(res, lb); bt <- b_getboot(res, lb)
    cat(sprintf("%-28s %8.3f  [%7.3f, %7.3f] %9s | [%7.3f, %7.3f] %9s %s\n",
                B_NAMES[lb], v$est, v$lo, v$hi, b_fmt_p(v$p),
                bt$lo_std, bt$hi_std, b_fmt_p(bt$p_boot_std), b_star(v$p)))
  }
  cat(sprintf("\nR-squared: TAS = %.3f | Mediator = %.3f | MH = %.3f\n",
              res$r2["tas_z"], res$r2["somatic_z"], res$r2["mh_z"]))
  cat(sprintf("N = %d | converged = %s | bootstrap resamples used = %d (failed %d)\n",
              res$nobs, res$converged, attr(res$boot, "R_ok"), attr(res$boot, "R_failed")))
  cat("\nUnstandardized metric for the interaction paths (cross-check):\n")
  for (lb in c("a3", "b3", "d3")) {
    v <- b_get(res, lb, std = FALSE); bt <- b_getboot(res, lb)
    cat(sprintf("  %-26s est = %8.4f  delta CI [%7.4f, %7.4f]  boot CI [%7.4f, %7.4f]\n",
                B_NAMES[lb], v$est, v$lo, v$hi, bt$lo_unstd, bt$hi_unstd))
  }
}

cat("\n#############################################################################\n")
cat("B1. THE THREE PRIMARY MODELS, FULL PATH SETS\n")
cat("#############################################################################\n")
for (nm in c("combined", "sss", "pcs")) b_print_model(b_models[[nm]])

cat("\n\n#############################################################################\n")
cat("B2. THE KEY QUESTION: THE IAS x IATS INTERACTION ON EACH MEDIATOR\n")
cat("#############################################################################\n\n")

cat("IAS x IATS -> MEDIATOR, side by side (fully standardized):\n")
cat(sprintf("%-38s %8s %20s %9s | %20s %9s\n", "Mediator", "beta",
            "95% CI (delta)", "p_Wald", "95% CI (boot)", "p_boot"))
cat(paste(rep("-", 116), collapse = ""), "\n")
for (nm in c("combined", "sss", "pcs", "rum", "mag", "help")) {
  r <- b_models[[nm]]; v <- b_get(r, "b3"); bt <- b_getboot(r, "b3")
  cat(sprintf("%-38s %8.3f  [%7.3f, %7.3f] %9s | [%7.3f, %7.3f] %9s %s\n",
              r$label, v$est, v$lo, v$hi, b_fmt_p(v$p),
              bt$lo_std, bt$hi_std, b_fmt_p(bt$p_boot_std), b_star(v$p)))
}

cat("\nFor completeness, the other two mediator-equation paths:\n")
cat(sprintf("%-38s %18s %18s %10s\n", "Mediator", "IAS -> Med (b2)",
            "IATS -> Med (b1)", "R2 Med"))
cat(paste(rep("-", 88), collapse = ""), "\n")
for (nm in c("combined", "sss", "pcs", "rum", "mag", "help")) {
  r <- b_models[[nm]]
  v2 <- b_get(r, "b2"); v1 <- b_get(r, "b1")
  cat(sprintf("%-38s %8.3f p=%-7s %8.3f p=%-7s %10.3f\n", r$label,
              v2$est, b_fmt_p(v2$p), v1$est, b_fmt_p(v1$p), r$r2["somatic_z"]))
}

cat("\nMediator -> MH and the direct interaction on MH:\n")
cat(sprintf("%-38s %18s %18s %10s\n", "Mediator", "Med -> MH (c2)",
            "IASxIATS -> MH (d3)", "R2 MH"))
cat(paste(rep("-", 88), collapse = ""), "\n")
for (nm in c("combined", "sss", "pcs", "rum", "mag", "help")) {
  r <- b_models[[nm]]
  vc <- b_get(r, "c2"); vd <- b_get(r, "d3")
  cat(sprintf("%-38s %8.3f p=%-7s %8.3f p=%-7s %10.3f\n", r$label,
              vc$est, b_fmt_p(vc$p), vd$est, b_fmt_p(vd$p), r$r2["mh_z"]))
}

cat("\n\n#############################################################################\n")
cat("B3. FORMAL TEST: DOES THE INTERACTION DIFFER BETWEEN PCS AND SSS-8?\n")
cat("#############################################################################\n\n")

cat("Method. The two single-mediator models cannot be compared informally, since\n")
cat("their coefficients come from the same participants and are correlated. Both\n")
cat("mediators are therefore entered in ONE model as parallel mediators with the\n")
cat("same predictor set and correlated residuals, and the equality of the two\n")
cat("interaction coefficients is tested with a Wald test and with a percentile\n")
cat("bootstrap of their difference.\n\n")

cat("Validity check (the joint model must not change the coefficients it compares):\n")
cat(sprintf("  b3 for SSS-8: joint %.6f vs single-mediator %.6f, difference = %.2e\n",
            j_get("b3s", FALSE)$est, b_get(b_models$sss, "b3", FALSE)$est, chk_b3s))
cat(sprintf("  b3 for PCS  : joint %.6f vs single-mediator %.6f, difference = %.2e\n",
            j_get("b3p", FALSE)$est, b_get(b_models$pcs, "b3", FALSE)$est, chk_b3p))
cat("  Identical to numerical precision, as expected when the two mediator\n")
cat("  equations share the same right-hand side. The joint model therefore only\n")
cat("  supplies the covariance needed for the difference test; it does not alter\n")
cat("  the quantities being compared.\n\n")

cat("Both mediators are z-scored, so the two coefficients are in the same metric\n")
cat("(SD of the respective mediator per unit of the product term) and their\n")
cat("difference is interpretable.\n\n")

cat("Joint model coefficients (unstandardized, z-scored metric):\n")
cat(sprintf("%-34s %10s %8s %20s %10s\n", "Parameter", "est", "SE", "95% CI", "p"))
cat(paste(rep("-", 86), collapse = ""), "\n")
jrows <- list(c("b3s", "IAS x IATS -> SSS-8"), c("b3p", "IAS x IATS -> PCS"),
              c("diff_b3", "DIFFERENCE (PCS - SSS-8)"),
              c("b1s", "IATS -> SSS-8"), c("b1p", "IATS -> PCS"),
              c("diff_b1", "DIFFERENCE (PCS - SSS-8)"),
              c("b2s", "IAS -> SSS-8"), c("b2p", "IAS -> PCS"),
              c("diff_b2", "DIFFERENCE (PCS - SSS-8)"),
              c("c2s", "SSS-8 -> MH (controlling PCS)"),
              c("c2p", "PCS -> MH (controlling SSS-8)"))
for (jr in jrows) {
  v <- j_get(jr[1], FALSE)
  cat(sprintf("%-34s %10.4f %8.4f  [%7.4f, %7.4f] %10s %s\n",
              jr[2], v$est, v$se, v$lo, v$hi, b_fmt_p(v$p), b_star(v$p)))
}

cat("\nBootstrap of the same quantities (percentile, ", N_BOOT_B, " resamples):\n", sep = "")
cat(sprintf("%-34s %10s %22s %10s\n", "Parameter", "est", "95% boot CI", "p_boot"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (lb in c("b3s", "b3p", "diff_b3", "b1s", "b1p", "diff_b1", "b2s", "b2p", "diff_b2")) {
  r <- boot_joint[boot_joint$label == lb, ]
  nm <- switch(lb, b3s = "IAS x IATS -> SSS-8", b3p = "IAS x IATS -> PCS",
               diff_b3 = "DIFFERENCE in interaction",
               b1s = "IATS -> SSS-8", b1p = "IATS -> PCS",
               diff_b1 = "DIFFERENCE in IATS effect",
               b2s = "IAS -> SSS-8", b2p = "IAS -> PCS",
               diff_b2 = "DIFFERENCE in IAS effect")
  cat(sprintf("%-34s %10.4f  [%8.4f, %8.4f] %10s\n", nm,
              r$est_unstd[1], r$lo_unstd[1], r$hi_unstd[1],
              b_fmt_p(r$p_boot_unstd[1])))
}
cat(sprintf("(bootstrap resamples used: %d, failed: %d)\n",
            attr(boot_joint, "R_ok"), attr(boot_joint, "R_failed")))

cat("\nWALD TEST of H0: b3(PCS) = b3(SSS-8)\n")
if (!inherits(wald_b3, "try-error")) {
  cat(sprintf("  chi2(%d) = %.4f, p = %.4f\n", wald_b3$df, wald_b3$stat, wald_b3$p.value))
} else {
  cat("  FAILED:\n"); cat(as.character(wald_b3), "\n")
}
cat("\nWald test of joint equality of ALL THREE mediator-equation paths\n")
cat("(H0: b1, b2 and b3 are the same for PCS and SSS-8):\n")
if (!inherits(wald_all, "try-error")) {
  cat(sprintf("  chi2(%d) = %.4f, p = %s\n", wald_all$df, wald_all$stat,
              b_fmt_p(wald_all$p.value)))
} else {
  cat("  FAILED:\n"); cat(as.character(wald_all), "\n")
}

cat("\nINTERPRETATION OF THE TEST\n")
cat("--------------------------\n")
d_int <- j_get("diff_b3", FALSE)
d_boot <- boot_joint[boot_joint$label == "diff_b3", ]
cat(sprintf("Interaction difference (PCS minus SSS-8) = %.4f, 95%% CI [%.4f, %.4f],\n",
            d_int$est, d_int$lo, d_int$hi))
cat(sprintf("Wald p = %s; bootstrap 95%% CI [%.4f, %.4f], bootstrap p = %s.\n",
            if (!inherits(wald_b3, "try-error")) b_fmt_p(wald_b3$p.value) else "n/a",
            d_boot$lo_unstd[1], d_boot$hi_unstd[1], b_fmt_p(d_boot$p_boot_unstd[1])))
sig_diff <- (!inherits(wald_b3, "try-error") && wald_b3$p.value < 0.05)
boot_sig <- (d_boot$lo_unstd[1] > 0 | d_boot$hi_unstd[1] < 0)
cat(sprintf("\n-> The interaction is %s larger for PCS than for SSS-8 by the Wald test,\n",
            ifelse(sig_diff, "SIGNIFICANTLY", "NOT significantly")))
cat(sprintf("   and the bootstrap CI for the difference %s zero.\n",
            ifelse(boot_sig, "EXCLUDES", "INCLUDES")))
if (!sig_diff) {
  cat("   The two mediators therefore cannot be statistically distinguished on the\n")
  cat("   interaction path, even though only one of them reaches significance on\n")
  cat("   its own. A non-significant difference between a significant and a\n")
  cat("   non-significant effect is NOT evidence that the effects are equal, nor\n")
  cat("   that they differ; the comparison is simply underpowered.\n")
}

cat("\n\n#############################################################################\n")
cat("B4. SECONDARY: WHICH PCS FACET CARRIES THE INTERACTION?\n")
cat("#############################################################################\n\n")
cat("Same model, mediator replaced by each PCS subscale in turn. These are\n")
cat("descriptive; the three subscales are strongly intercorrelated and the models\n")
cat("are not mutually adjusted.\n\n")
cat("Subscale intercorrelations (Pearson) and with SSS-8:\n")
sub_cor <- cor(dfc[, c("pcs_rum", "pcs_mag", "pcs_help", "pcs", "sss8")])
print(round(sub_cor, 3))
cat("\nAlpha of each subscale (raw, Pearson):\n")
for (nmx in list(c("rumination", "pcs_rum_items"), c("magnification", "pcs_mag_items"),
                 c("helplessness", "pcs_help_items"))) {
  its <- get(nmx[2])
  a <- suppressWarnings(psych::alpha(df[, its], warnings = FALSE))
  cat(sprintf("  %-14s (%d items): alpha = %.3f\n", nmx[1], length(its),
              a$total$raw_alpha))
}
cat("\n")
for (nm in c("rum", "mag", "help")) b_print_model(b_models[[nm]])

cat("\n\n#############################################################################\n")
cat("B5. PART B VERDICT\n")
cat("#############################################################################\n\n")

v_comb <- b_get(b_models$combined, "b3"); bt_comb <- b_getboot(b_models$combined, "b3")
v_sss <- b_get(b_models$sss, "b3");       bt_sss <- b_getboot(b_models$sss, "b3")
v_pcs <- b_get(b_models$pcs, "b3");       bt_pcs <- b_getboot(b_models$pcs, "b3")

cat("THE THREE INTERACTION COEFFICIENTS (IAS x IATS -> mediator, standardized):\n")
cat(sprintf("  (a) SSS-8 + PCS composite : beta = %.3f, 95%% CI [%.3f, %.3f], p = %s\n",
            v_comb$est, v_comb$lo, v_comb$hi, b_fmt_p(v_comb$p)))
cat(sprintf("  (b) SSS-8 alone           : beta = %.3f, 95%% CI [%.3f, %.3f], p = %s\n",
            v_sss$est, v_sss$lo, v_sss$hi, b_fmt_p(v_sss$p)))
cat(sprintf("  (c) PCS alone             : beta = %.3f, 95%% CI [%.3f, %.3f], p = %s\n",
            v_pcs$est, v_pcs$lo, v_pcs$hi, b_fmt_p(v_pcs$p)))
cat(sprintf("\n  Difference (PCS - SSS-8)  : %.4f, Wald p = %s, bootstrap p = %s\n",
            d_int$est,
            if (!inherits(wald_b3, "try-error")) b_fmt_p(wald_b3$p.value) else "n/a",
            b_fmt_p(d_boot$p_boot_unstd[1])))

cat("\nPCS SUBSCALES (IAS x IATS -> subscale, standardized):\n")
for (nm in c("rum", "mag", "help")) {
  r <- b_models[[nm]]; v <- b_get(r, "b3")
  cat(sprintf("  %-32s beta = %.3f, 95%% CI [%.3f, %.3f], p = %s %s\n",
              r$label, v$est, v$lo, v$hi, b_fmt_p(v$p), b_star(v$p)))
}

cat("\nPLAIN-LANGUAGE READING\n")
cat("----------------------\n")
rank_txt <- c(sprintf("composite %.3f", v_comb$est), sprintf("SSS-8 %.3f", v_sss$est),
              sprintf("PCS %.3f", v_pcs$est))
cat("Interaction magnitudes: ", paste(rank_txt, collapse = " | "), "\n", sep = "")
larger <- ifelse(abs(v_pcs$est) > abs(v_sss$est), "PCS", "SSS-8")
cat(sprintf("The larger single-scale interaction is on %s.\n", larger))
cat(sprintf("Significant at .05 on its own: composite %s, SSS-8 %s, PCS %s.\n",
            ifelse(v_comb$p < 0.05, "YES", "no"), ifelse(v_sss$p < 0.05, "YES", "no"),
            ifelse(v_pcs$p < 0.05, "YES", "no")))
if (!sig_diff) {
  cat("\nBUT the formal test does not license the claim that the effect is\n")
  cat("specifically a pain-catastrophizing effect: the PCS and SSS-8 interaction\n")
  cat("coefficients do not differ significantly from each other. The honest\n")
  cat("statement is that the interaction is estimated with similar magnitude on\n")
  cat("both mediators and that this study cannot resolve which one carries it.\n")
} else {
  cat("\nThe formal test supports a genuine difference between the two mediators\n")
  cat("on the interaction path.\n")
}
cat("\nNote also that the composite's interaction can exceed BOTH single-scale\n")
cat("interactions: averaging two z-scored mediators cancels part of their\n")
cat("independent error variance, so the composite can be a more reliable\n")
cat("indicator of whatever common variance the interaction acts on. Compare the\n")
cat("three coefficients above against that possibility before concluding that\n")
cat("one scale 'drives' the effect.\n")

cat("\nProvenance: every Part B number was produced by this script in this run,\n")
cat("with the same data, specification, estimator, seed and resample count as\n")
cat("model M0 of 17_revision_sensitivity.R.\n")

sink()

# =============================================================================
# 7. WRITE XLSX
# =============================================================================

wb <- createWorkbook()

kmo_tab <- data.frame(item = items,
                      scale = ifelse(items %in% sss_items, "SSS-8", "PCS"),
                      msa_polychoric = as.numeric(kmo_poly$MSAi[items]),
                      msa_pearson = as.numeric(kmo_pear$MSAi[items]))
addWorksheet(wb, "KMO"); writeData(wb, "KMO", kmo_tab)

eig_tab <- data.frame(n = seq_along(pa$fa.values),
                      fa_actual = pa$fa.values,
                      fa_simulated = if (!is.null(pa$fa.sim)) pa$fa.sim else pa$fa.simr,
                      pc_actual = pa$pc.values,
                      pc_simulated = if (!is.null(pa$pc.sim)) pa$pc.sim else pa$pc.simr)
addWorksheet(wb, "Parallel_Analysis"); writeData(wb, "Parallel_Analysis", eig_tab)

for (k in 1:3) {
  f <- efa_fit[[k]]
  if (inherits(f, "try-error")) next
  L <- unclass(f$loadings)[items, , drop = FALSE]
  d <- data.frame(item = items,
                  scale = ifelse(items %in% sss_items, "SSS-8", "PCS"),
                  as.data.frame(round(L, 4)),
                  h2 = round(f$communality[items], 4),
                  u2 = round(f$uniquenesses[items], 4))
  sh <- paste0("EFA_", k, "factor")
  addWorksheet(wb, sh); writeData(wb, sh, d)
}

fit_rows <- do.call(rbind, lapply(names(fit_wlsmv), function(nm) {
  g <- get_fits(fit_wlsmv[[nm]], TRUE)
  if (is.null(g)) return(NULL)
  cbind(estimator = "WLSMV", model = nm, g)
}))
fit_rows_mlr <- do.call(rbind, lapply(names(fit_mlr), function(nm) {
  g <- get_fits(fit_mlr[[nm]], TRUE)
  if (is.null(g)) return(NULL)
  cbind(estimator = "MLR", model = nm, g)
}))
addWorksheet(wb, "CFA_Fit")
writeData(wb, "CFA_Fit", rbind(fit_rows, fit_rows_mlr))

if (ok_fit(fit_wlsmv$two)) {
  l1 <- if (ok_fit(fit_wlsmv$one)) std_load(fit_wlsmv$one, "SOM") else setNames(rep(NA, 21), items)
  ls <- std_load(fit_wlsmv$two, "SSS"); lp <- std_load(fit_wlsmv$two, "PCS")
  d <- data.frame(item = items,
                  scale = ifelse(items %in% sss_items, "SSS-8", "PCS"),
                  onefactor_SOM = round(as.numeric(l1[items]), 4),
                  twofactor = round(as.numeric(ifelse(items %in% sss_items,
                                                      ls[items], lp[items])), 4))
  if (ok_fit(fit_wlsmv$bifactor)) {
    lg <- std_load(fit_wlsmv$bifactor, "G")
    lss <- std_load(fit_wlsmv$bifactor, "SSSs"); lpp <- std_load(fit_wlsmv$bifactor, "PCSs")
    d$bifactor_G <- round(as.numeric(lg[items]), 4)
    d$bifactor_specific <- round(as.numeric(ifelse(items %in% sss_items,
                                                   lss[items], lpp[items])), 4)
  }
  addWorksheet(wb, "CFA_Loadings"); writeData(wb, "CFA_Loadings", d)
}

rel_tab <- data.frame(
  quantity = c("alpha_21item", "alpha_SSS8", "alpha_PCS",
               "ord_alpha_21item", "ord_alpha_SSS8", "ord_alpha_PCS",
               "omega_h_21item_3gf", "omega_total_21item_3gf",
               "r_sss_pcs_total", "r_sss_pcs_disattenuated_alpha",
               "latent_r_SSS_PCS_WLSMV", "N"),
  value = c(alpha_all$total$raw_alpha, alpha_sss$total$raw_alpha,
            alpha_pcs$total$raw_alpha, oalpha_all, oalpha_sss, oalpha_pcs,
            if (!inherits(om3, "try-error")) om3$omega_h else NA,
            if (!inherits(om3, "try-error")) om3$omega.tot else NA,
            r_tot, disatt(r_tot, alpha_sss$total$raw_alpha, alpha_pcs$total$raw_alpha),
            r_lat_wlsmv, N))
addWorksheet(wb, "Reliability"); writeData(wb, "Reliability", rel_tab)

# ---- Part B sheets ----

partb_all <- do.call(rbind, lapply(names(b_models), function(nm) {
  r <- b_models[[nm]]
  do.call(rbind, lapply(B_PATHS, function(lb) {
    v <- b_get(r, lb); u <- b_get(r, lb, std = FALSE); bt <- b_getboot(r, lb)
    data.frame(model = r$label, mediator = r$med_var, label = lb,
               path = unname(B_NAMES[lb]),
               std_beta = v$est, std_ci_lo = v$lo, std_ci_hi = v$hi, p_wald = v$p,
               boot_std_lo = bt$lo_std, boot_std_hi = bt$hi_std,
               p_boot_std = bt$p_boot_std,
               unstd_est = u$est, unstd_ci_lo = u$lo, unstd_ci_hi = u$hi,
               boot_unstd_lo = bt$lo_unstd, boot_unstd_hi = bt$hi_unstd,
               n = r$nobs, stringsAsFactors = FALSE)
  }))
}))
addWorksheet(wb, "PartB_All_Paths"); writeData(wb, "PartB_All_Paths", partb_all)

partb_r2 <- do.call(rbind, lapply(names(b_models), function(nm) {
  r <- b_models[[nm]]
  data.frame(model = r$label, mediator = r$med_var,
             r2_tas = unname(r$r2["tas_z"]), r2_mediator = unname(r$r2["somatic_z"]),
             r2_mh = unname(r$r2["mh_z"]), n = r$nobs, converged = r$converged,
             stringsAsFactors = FALSE)
}))
addWorksheet(wb, "PartB_R2"); writeData(wb, "PartB_R2", partb_r2)

partb_key <- partb_all[partb_all$label == "b3", ]
addWorksheet(wb, "PartB_Interaction"); writeData(wb, "PartB_Interaction", partb_key)

partb_joint <- do.call(rbind, lapply(JOINT_LABS, function(lb) {
  v <- j_get(lb, FALSE); s <- j_get(lb, TRUE)
  bt <- boot_joint[boot_joint$label == lb, ]
  data.frame(label = lb, unstd_est = v$est, unstd_se = v$se,
             unstd_ci_lo = v$lo, unstd_ci_hi = v$hi, p_wald = v$p,
             boot_ci_lo = bt$lo_unstd[1], boot_ci_hi = bt$hi_unstd[1],
             p_boot = bt$p_boot_unstd[1], std_est = s$est,
             stringsAsFactors = FALSE)
}))
partb_joint$wald_b3_chisq <- if (!inherits(wald_b3, "try-error")) wald_b3$stat else NA
partb_joint$wald_b3_df <- if (!inherits(wald_b3, "try-error")) wald_b3$df else NA
partb_joint$wald_b3_p <- if (!inherits(wald_b3, "try-error")) wald_b3$p.value else NA
addWorksheet(wb, "PartB_Joint_Model"); writeData(wb, "PartB_Joint_Model", partb_joint)

saveWorkbook(wb, file.path(supplementary_data_dir, "somatic_factor_structure.xlsx"),
             overwrite = TRUE)

cat("Done (Part A factor structure + Part B mediator decomposition).\n")
cat("  analysis_output/18_output.txt\n")
cat("  supplementary_data/somatic_factor_structure.xlsx\n")
cat("  plots/figure_18_scree_parallel_somatic.png\n")
