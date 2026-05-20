"""
================================================================================
ALEXITHYMIA SUBSCALES AND THE SOMATIC PATHWAY
================================================================================

The canonical literature predicts that alexithymia (particularly DIF) drives
somatisation. In this dataset, the profile with the highest alexithymia
(Disconnected / Uncertain, IAS-low / IATS-average) does NOT show the highest
somatic burden — that belongs to the Hypervigilant profile (IAS-average / IATS-high).

This script tests three concrete questions that disentangle the two pathways:

  Q1. Do the three profiles differ in WHICH facets of alexithymia they show
      (DIF, DDF, EOT)?

  Q2. Does any TAS facet predict somatic burden ABOVE AND BEYOND the
      IAS x IATS interaction? Hierarchical regression: step 1 = DIF/DDF/EOT,
      step 2 = + IAS + IATS, step 3 = + IAS x IATS.

  Q3. Decomposing the Uncertain-vs-Hypervigilant somatic mean difference:
      how much is "explained" by each predictor in the joint model?

Outputs land in supplementary_data/ and a four-panel figure is saved to
plots/supplementary_plots/3_s_10_alexithymia_somatic_pathway.png.
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

warnings.filterwarnings("ignore")
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

PLOTS_DIR = "C:/code/projects/interoceptive-profiling/plots/supplementary_plots"
DATA_OUT = "C:/code/projects/interoceptive-profiling/supplementary_data"
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(DATA_OUT, exist_ok=True)

CANONICAL_COLORS = {
    "Hypervigilant": "#FF7F00",
    "Efficient": "#4DAF4A",
    "Uncertain": "#E41A1C",
}
SUBSCALE_COLORS = {"DIF": "#3B7DD8", "DDF": "#8C5BA3", "EOT": "#888888"}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
df = pd.read_csv("C:/code/projects/interoceptive-profiling/dfc_interoception_profiling.csv")
labels_df = pd.read_csv("C:/code/projects/interoceptive-profiling/supplementary_data/data_with_k_means_clusters.csv")
df["profile"] = labels_df["cluster_label"].values

# Map long labels to canonical short labels
label_map = {
    "High Accuracy / Low Attention": "Efficient",
    "Average Accuracy / High Attention": "Hypervigilant",
    "Low Accuracy / Average Attention": "Uncertain",
}
df["profile_short"] = df["profile"].map(label_map)

# TAS-20 subscale items (items 4, 5, 10, 18, 19 already reverse-scored in dataset)
dif_items = [1, 3, 6, 7, 9, 13, 14]
ddf_items = [2, 4, 11, 12, 17]
eot_items = [5, 8, 10, 15, 16, 18, 19, 20]

df["dif"] = df[[f"tas_{i}" for i in dif_items]].sum(axis=1)
df["ddf"] = df[[f"tas_{i}" for i in ddf_items]].sum(axis=1)
df["eot"] = df[[f"tas_{i}" for i in eot_items]].sum(axis=1)
# sanity: subscales sum to TAS total
assert (np.abs(df["tas"] - (df["dif"] + df["ddf"] + df["eot"])).max() < 0.01)

# Somatic composite as used in the main SEM
df["somatic_combined"] = (
    stats.zscore(df["sss8"]) + stats.zscore(df["pcs"])
) / 2

# z-scored predictors
df["ias_z"] = stats.zscore(df["ias"])
df["iats_z"] = stats.zscore(df["iats"])
df["dif_z"] = stats.zscore(df["dif"])
df["ddf_z"] = stats.zscore(df["ddf"])
df["eot_z"] = stats.zscore(df["eot"])
df["ias_x_iats"] = df["ias_z"] * df["iats_z"]
df["somatic_z"] = stats.zscore(df["somatic_combined"])
df["sss8_z"] = stats.zscore(df["sss8"])
df["pcs_z"] = stats.zscore(df["pcs"])

sub = df.dropna(
    subset=["dif", "ddf", "eot", "ias", "iats", "sss8", "pcs", "profile_short"]
).copy()
N = len(sub)
print(f"N = {N}")

# ===========================================================================
# Q1. TAS facet means by profile
# ===========================================================================
print("\n[Q1] TAS facet means by profile")
facet_means = sub.groupby("profile_short")[["dif", "ddf", "eot", "tas"]].agg(["mean", "std", "count"])
print(facet_means.round(2).to_string())

# Per-facet one-way ANOVA across profiles
facet_anova_rows = []
for facet in ["dif", "ddf", "eot", "tas"]:
    groups = [
        sub.loc[sub["profile_short"] == p, facet].values
        for p in ["Efficient", "Hypervigilant", "Uncertain"]
    ]
    f, p = stats.f_oneway(*groups)
    grand = np.concatenate(groups)
    ss_total = ((grand - grand.mean()) ** 2).sum()
    ss_within = sum(((g - g.mean()) ** 2).sum() for g in groups)
    eta2 = (ss_total - ss_within) / ss_total
    facet_anova_rows.append({
        "facet": facet, "F": round(f, 2), "p": p, "eta2_p": round(eta2, 4),
        "mean_Efficient": round(groups[0].mean(), 2),
        "mean_Hypervigilant": round(groups[1].mean(), 2),
        "mean_Uncertain": round(groups[2].mean(), 2),
    })
facet_anova_df = pd.DataFrame(facet_anova_rows)
print("\nFacet ANOVAs (F, p, eta2_p):")
print(facet_anova_df.to_string(index=False))
facet_anova_df.to_csv(os.path.join(DATA_OUT, "alexithymia_facets_by_profile.csv"), index=False)

# Tukey pairwise for each facet (Uncertain vs Hypervigilant is the focal contrast)
print("\nPairwise (Welch t) Uncertain vs Hypervigilant per facet:")
pair_rows = []
for facet in ["dif", "ddf", "eot", "tas"]:
    a = sub.loc[sub["profile_short"] == "Uncertain", facet].values
    b = sub.loc[sub["profile_short"] == "Hypervigilant", facet].values
    t, p = stats.ttest_ind(a, b, equal_var=False)
    d = (a.mean() - b.mean()) / np.sqrt(((a.std(ddof=1) ** 2 + b.std(ddof=1) ** 2) / 2))
    pair_rows.append({
        "facet": facet, "t": round(t, 3), "p": p,
        "mean_Uncertain": round(a.mean(), 2),
        "mean_Hypervigilant": round(b.mean(), 2),
        "diff": round(a.mean() - b.mean(), 2),
        "cohen_d": round(d, 3),
    })
pair_df = pd.DataFrame(pair_rows)
print(pair_df.to_string(index=False))
pair_df.to_csv(os.path.join(DATA_OUT, "alexithymia_facets_uncertain_vs_hypervigilant.csv"), index=False)

# ===========================================================================
# Q2. Hierarchical regression: somatic_combined ~ steps
# ===========================================================================
print("\n[Q2] Hierarchical regression predicting somatic burden")

steps = {
    "step1_TAS_total":      "somatic_z ~ tas",
    "step1b_facets":        "somatic_z ~ dif_z + ddf_z + eot_z",
    "step2_facets_plus_intero": "somatic_z ~ dif_z + ddf_z + eot_z + ias_z + iats_z",
    "step3_full":           "somatic_z ~ dif_z + ddf_z + eot_z + ias_z + iats_z + ias_x_iats",
}

models = {}
step_rows = []
prev_r2 = 0.0
for name, formula in steps.items():
    m = ols(formula, data=sub).fit()
    models[name] = m
    delta = m.rsquared - prev_r2
    step_rows.append({
        "step": name,
        "formula": formula,
        "n_params": int(m.df_model),
        "R2": round(m.rsquared, 4),
        "adj_R2": round(m.rsquared_adj, 4),
        "delta_R2": round(delta, 4),
        "F": round(m.fvalue, 2),
        "p": m.f_pvalue,
    })
    prev_r2 = m.rsquared
step_df = pd.DataFrame(step_rows)
print(step_df.to_string(index=False))
step_df.to_csv(os.path.join(DATA_OUT, "alexithymia_somatic_hierarchical_r2.csv"), index=False)

# Detailed coefficients from the full model
full = models["step3_full"]
coef_df = pd.DataFrame({
    "predictor": full.params.index,
    "beta": full.params.values.round(4),
    "se": full.bse.values.round(4),
    "t": full.tvalues.values.round(3),
    "p": full.pvalues.values,
    "ci_lower": full.conf_int()[0].values.round(4),
    "ci_upper": full.conf_int()[1].values.round(4),
}).query("predictor != 'Intercept'").reset_index(drop=True)
print("\nFull-model standardised coefficients (somatic_z ~ DIF + DDF + EOT + IAS + IATS + IAS x IATS):")
print(coef_df.to_string(index=False))
coef_df.to_csv(os.path.join(DATA_OUT, "alexithymia_somatic_full_model_betas.csv"), index=False)

# Compare: model with TAS_total + IAS x IATS only (no facets) vs full model
m_tasonly_int = ols("somatic_z ~ tas + ias_z + iats_z + ias_x_iats", data=sub).fit()
m_intero_only = ols("somatic_z ~ ias_z + iats_z + ias_x_iats", data=sub).fit()
m_facets_only = ols("somatic_z ~ dif_z + ddf_z + eot_z", data=sub).fit()
m_tas_only = ols("somatic_z ~ tas", data=sub).fit()
m_dif_only = ols("somatic_z ~ dif_z", data=sub).fit()

compare_rows = [
    {"model": "DIF only", "R2": round(m_dif_only.rsquared, 4),
     "predictors": "DIF"},
    {"model": "TAS total only", "R2": round(m_tas_only.rsquared, 4),
     "predictors": "TAS"},
    {"model": "DIF + DDF + EOT", "R2": round(m_facets_only.rsquared, 4),
     "predictors": "facets"},
    {"model": "IAS + IATS + IAS x IATS", "R2": round(m_intero_only.rsquared, 4),
     "predictors": "intero"},
    {"model": "TAS + IAS + IATS + IAS x IATS", "R2": round(m_tasonly_int.rsquared, 4),
     "predictors": "TAS + intero"},
    {"model": "Full (facets + intero + interaction)", "R2": round(full.rsquared, 4),
     "predictors": "facets + intero"},
]
print("\nModel comparison (R^2 for somatic_z):")
print(pd.DataFrame(compare_rows).to_string(index=False))
pd.DataFrame(compare_rows).to_csv(os.path.join(DATA_OUT, "alexithymia_somatic_model_comparison.csv"), index=False)

# Variance shared/unique decomposition for the question:
#   how much somatic variance is unique to TAS vs to IAS+IATS+interaction?
r2_facets_only = m_facets_only.rsquared
r2_intero_only = m_intero_only.rsquared
r2_full = full.rsquared
unique_facets = r2_full - r2_intero_only        # added by facets over intero
unique_intero = r2_full - r2_facets_only        # added by intero over facets
shared = r2_full - unique_facets - unique_intero
print(f"\nVariance decomposition of somatic_z (R^2):")
print(f"  unique to TAS facets       = {unique_facets:.4f}")
print(f"  unique to IAS+IATS+int     = {unique_intero:.4f}")
print(f"  shared                     = {shared:.4f}")
print(f"  total (full model R^2)     = {r2_full:.4f}")

with open(os.path.join(DATA_OUT, "alexithymia_somatic_variance_decomposition.json"), "w") as f:
    json.dump({
        "r2_facets_only": float(r2_facets_only),
        "r2_intero_only": float(r2_intero_only),
        "r2_full": float(r2_full),
        "unique_to_facets": float(unique_facets),
        "unique_to_intero": float(unique_intero),
        "shared": float(shared),
    }, f, indent=2)

# ===========================================================================
# Q3. Decomposing the Uncertain-vs-Hypervigilant somatic gap
# ===========================================================================
print("\n[Q3] Decomposing the Uncertain - Hypervigilant somatic gap")

# Centred predictors so the regression coefficients act on standardised units;
# compute each profile's mean predictor vector and multiply by full-model betas.
predictors = ["dif_z", "ddf_z", "eot_z", "ias_z", "iats_z", "ias_x_iats"]
beta = full.params[predictors]
profile_means_predictors = sub.groupby("profile_short")[predictors].mean()
profile_somatic_means = sub.groupby("profile_short")["somatic_z"].mean()
print("\nProfile means of predictors (z units):")
print(profile_means_predictors.round(3).to_string())
print("\nProfile means of somatic_z:")
print(profile_somatic_means.round(3).to_string())

# Contribution of each predictor to the Hypervigilant - Uncertain gap
gap_predictors = profile_means_predictors.loc["Hypervigilant"] - profile_means_predictors.loc["Uncertain"]
contributions = gap_predictors * beta
print("\nPredicted contributions to Hypervigilant - Uncertain somatic_z gap:")
gap_df = pd.DataFrame({
    "predictor": predictors,
    "mean_Uncertain": [profile_means_predictors.loc["Uncertain", p] for p in predictors],
    "mean_Hypervigilant": [profile_means_predictors.loc["Hypervigilant", p] for p in predictors],
    "mean_diff": gap_predictors.values,
    "beta_full_model": beta.values,
    "contribution_to_gap": contributions.values,
})
gap_df["pct_of_total_gap"] = (
    gap_df["contribution_to_gap"] /
    gap_df["contribution_to_gap"].sum() * 100
).round(2)
gap_df = gap_df.round(4)
print(gap_df.to_string(index=False))
gap_df.to_csv(os.path.join(DATA_OUT, "alexithymia_somatic_gap_decomposition.csv"), index=False)

predicted_gap = contributions.sum()
observed_gap = profile_somatic_means["Hypervigilant"] - profile_somatic_means["Uncertain"]
print(f"\nObserved Hypervigilant - Uncertain somatic_z gap = {observed_gap:.4f}")
print(f"Predicted gap from joint model = {predicted_gap:.4f}")

# ===========================================================================
# FIGURE: four panels
# ===========================================================================
print("\n[Figure] Building 3_s_10_alexithymia_somatic_pathway.png")

fig = plt.figure(figsize=(15.0, 9.0))
gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.26,
                      left=0.08, right=0.97, top=0.94, bottom=0.08)

# Panel A: TAS facet means by profile
axA = fig.add_subplot(gs[0, 0])
profiles_order = ["Efficient", "Hypervigilant", "Uncertain"]
facets_order = ["DIF", "DDF", "EOT"]
xpos = np.arange(len(profiles_order))
w = 0.27
for i, facet in enumerate(facets_order):
    means = [sub.loc[sub["profile_short"] == p, facet.lower()].mean()
             for p in profiles_order]
    sems = [
        sub.loc[sub["profile_short"] == p, facet.lower()].std(ddof=1) /
        np.sqrt((sub["profile_short"] == p).sum())
        for p in profiles_order
    ]
    axA.bar(xpos + (i - 1) * w, means, w, color=SUBSCALE_COLORS[facet],
            yerr=sems, capsize=3, label=facet, alpha=0.95)
axA.set_xticks(xpos)
axA.set_xticklabels(profiles_order)
axA.set_ylabel("Subscale total")
axA.set_title("a  TAS-20 facets by profile", fontsize=11, weight="bold", loc="left")
axA.legend(fontsize=9, loc="upper left", framealpha=0.85)
axA.grid(alpha=0.25, axis="y")
# Annotate Uncertain vs Hypervigilant DIF gap (the canonical-model expectation)
dif_unc = sub.loc[sub["profile_short"] == "Uncertain", "dif"].mean()
dif_hyp = sub.loc[sub["profile_short"] == "Hypervigilant", "dif"].mean()
axA.annotate(f"Δ DIF\n{dif_unc - dif_hyp:+.1f}",
             xy=(2 + (-1) * w, dif_unc), xytext=(0.9, 24),
             fontsize=9, color="#3B7DD8",
             arrowprops=dict(arrowstyle="->", color="#3B7DD8", alpha=0.7))

# Panel B: forest of standardised betas in the full model
axB = fig.add_subplot(gs[0, 1])
pred_pretty = {
    "dif_z": "DIF",
    "ddf_z": "DDF",
    "eot_z": "EOT",
    "ias_z": "IAS",
    "iats_z": "IATS",
    "ias_x_iats": "IAS × IATS",
}
ord_pred = ["dif_z", "ddf_z", "eot_z", "ias_z", "iats_z", "ias_x_iats"]
betas = [full.params[p] for p in ord_pred]
lows = [full.conf_int().loc[p, 0] for p in ord_pred]
highs = [full.conf_int().loc[p, 1] for p in ord_pred]
ps = [full.pvalues[p] for p in ord_pred]
y = np.arange(len(ord_pred))[::-1]
colors_pred = ["#3B7DD8", "#8C5BA3", "#888888", "#2C7A2C", "#D95F02", "#B22222"]
for yi, b, lo, hi, p, c in zip(y, betas, lows, highs, ps, colors_pred):
    axB.errorbar(b, yi, xerr=[[b - lo], [hi - b]], fmt="o", color=c,
                 markersize=8, capsize=4, linewidth=2)
    star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
    axB.text(hi + 0.01, yi, f"  β={b:+.3f} {star}",
             va="center", fontsize=9, color=c)
axB.axvline(0, color="#888", linestyle=":", linewidth=1)
axB.set_yticks(y)
axB.set_yticklabels([pred_pretty[p] for p in ord_pred])
axB.set_xlabel("Standardised β (predicting somatic burden, z)")
axB.set_title("b  Joint predictors of somatic burden", fontsize=11, weight="bold", loc="left")
axB.set_xlim(-0.05, max(highs) + 0.20)
axB.grid(alpha=0.25, axis="x")

# Panel C: variance decomposition / hierarchical R^2
axC = fig.add_subplot(gs[1, 0])
# show R^2 for nested models, plus shared/unique decomposition
bars = [
    ("DIF only", r2_facets_only - (full.rsquared - r2_intero_only) - shared,
     m_dif_only.rsquared),  # use simple model R^2
    ("TAS total only", None, m_tas_only.rsquared),
    ("Facets only (DIF+DDF+EOT)", None, r2_facets_only),
    ("Intero only (IAS+IATS+IAS×IATS)", None, r2_intero_only),
    ("Full (facets + intero + interaction)", None, full.rsquared),
]
labels = [b[0] for b in bars]
r2_vals = [b[2] for b in bars]
y2 = np.arange(len(bars))[::-1]
axC.barh(y2, r2_vals, color=["#3B7DD8", "#5DA0E0", "#7AB8E8", "#D95F02", "#2C3E50"])
for yi, v in zip(y2, r2_vals):
    axC.text(v + 0.005, yi, f"R²={v:.3f}", va="center", fontsize=9)
axC.set_yticks(y2)
axC.set_yticklabels(labels)
axC.set_xlim(0, max(r2_vals) * 1.18)
axC.set_xlabel(r"R$^2$ for somatic burden")
axC.set_title("c  Variance in somatic burden explained", fontsize=11, weight="bold", loc="left")
axC.grid(alpha=0.25, axis="x")
# Add unique/shared decomposition as inset text
axC.text(
    0.02, 0.04,
    f"Variance decomposition (Full R² = {full.rsquared:.3f}):\n"
    f"  unique to TAS facets       = {unique_facets:.3f}\n"
    f"  unique to IAS+IATS+IAS×IATS = {unique_intero:.3f}\n"
    f"  shared                      = {shared:.3f}",
    transform=axC.transAxes,
    fontsize=8.5, family="monospace",
    bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFF8DC", edgecolor="#999", alpha=0.9),
)

# Panel D: profile-gap decomposition (Hypervigilant - Uncertain)
axD = fig.add_subplot(gs[1, 1])
pretty_labels = [pred_pretty[p] for p in ord_pred]
contribs = [gap_df.loc[gap_df["predictor"] == p, "contribution_to_gap"].values[0] for p in ord_pred]
colors_d = colors_pred
xp = np.arange(len(ord_pred))
axD.bar(xp, contribs, color=colors_d, alpha=0.9)
for xi, c in zip(xp, contribs):
    axD.text(xi, c + (0.005 if c >= 0 else -0.015),
             f"{c:+.3f}", ha="center", fontsize=9,
             color="black" if c >= 0 else "darkred")
axD.axhline(0, color="#888", linewidth=0.8)
axD.axhline(observed_gap, color="#2C3E50", linestyle="--", linewidth=1.5,
            label=f"observed gap = {observed_gap:+.3f}")
axD.set_xticks(xp)
axD.set_xticklabels(pretty_labels, rotation=15)
axD.set_ylabel("Contribution to Hypervigilant−Uncertain somatic_z gap")
axD.set_title("d  What explains the Hypervigilant−Uncertain somatic gap?",
              fontsize=11, weight="bold", loc="left")
axD.legend(fontsize=9, loc="upper left", framealpha=0.85)
axD.grid(alpha=0.25, axis="y")

fig_path = os.path.join(PLOTS_DIR, "3_s_10_alexithymia_somatic_pathway.png")
fig.savefig(fig_path, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Figure saved to {fig_path}")

# ===========================================================================
# Machine-readable summary
# ===========================================================================
summary = {
    "n": int(N),
    "facets_by_profile_anova": facet_anova_df.to_dict(orient="records"),
    "uncertain_vs_hypervigilant_pairwise": pair_df.to_dict(orient="records"),
    "hierarchical_R2": step_df.to_dict(orient="records"),
    "full_model_coefs": coef_df.to_dict(orient="records"),
    "model_comparison_R2": compare_rows,
    "variance_decomposition": {
        "r2_facets_only": float(r2_facets_only),
        "r2_intero_only": float(r2_intero_only),
        "r2_full": float(r2_full),
        "unique_to_facets": float(unique_facets),
        "unique_to_intero": float(unique_intero),
        "shared": float(shared),
    },
    "gap_decomposition": gap_df.to_dict(orient="records"),
    "observed_hyp_minus_unc_somatic_z": float(observed_gap),
    "predicted_hyp_minus_unc_somatic_z": float(predicted_gap),
}
with open(os.path.join(DATA_OUT, "alexithymia_somatic_summary.json"), "w") as f:
    json.dump(summary, f, indent=2, default=str)
print(f"\nSummary JSON: {os.path.join(DATA_OUT, 'alexithymia_somatic_summary.json')}")
print("\nDone.")
