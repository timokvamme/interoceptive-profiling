"""
================================================================================
ITEM-LEVEL CLUSTERING SENSITIVITY ANALYSIS
================================================================================

Tests whether the 3-profile k=3 solution is stable when the clustering algorithm
operates on all 42 IAS+IATS items rather than on the two summary scores. Two
alternatives are evaluated:

  (1) k-means on z-standardized 42 items (treats Likert as continuous)
  (2) PAM with Gower distance on 42 items treated as ordinal/mixed-type

For each alternative the script reports k-selection diagnostics, profile
centroids, cross-tabulation against the primary solution, outcome ANOVAs, and
bootstrap stability. Outputs land in supplementary_data/ and a six-panel figure
is saved to plots/supplementary_plots/.

The primary solution clusters on the two summary scores; this script asks
whether moving to the full item space materially changes the structure or the
mental-health architecture downstream.
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Ellipse
import matplotlib.transforms as mtransforms
from sklearn.cluster import KMeans
from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score,
)
from sklearn.preprocessing import StandardScaler
from scipy import stats
import gower
from sklearn_extra.cluster import KMedoids

warnings.filterwarnings("ignore")
RNG = np.random.default_rng(42)

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
GOWER_COLORS = {
    "Low Attention": "#1F77B4",
    "Average": "#7F7F7F",
    "High Attention": "#9467BD",
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
df = pd.read_csv("C:/code/projects/interoceptive-profiling/dfc_interoception_profiling.csv")
labels_df = pd.read_csv("C:/code/projects/interoceptive-profiling/supplementary_data/data_with_k_means_clusters.csv")

ias_cols = [f"ias_{i}" for i in range(1, 22)]
iats_cols = [f"iats_{i}" for i in range(1, 22)]
all_items = ias_cols + iats_cols
outcomes = ["tas", "sss8", "pcs", "phq9", "gad7", "stai"]
outcome_pretty = {
    "tas": "TAS (alexithymia)",
    "sss8": "SSS-8 (somatic)",
    "pcs": "PCS (catastrophizing)",
    "phq9": "PHQ-9 (depression)",
    "gad7": "GAD-7 (anxiety)",
    "stai": "STAI (trait anxiety)",
}

assert len(df) == len(labels_df), f"row mismatch {len(df)} vs {len(labels_df)}"
df["orig_cluster"] = labels_df["cluster"].values
df["orig_label"] = labels_df["cluster_label"].values
df["orig_ordered"] = labels_df["cluster_ordered"].values

sub = df.dropna(subset=all_items).copy()
N = len(sub)
print(f"N = {N}")

X_items = sub[all_items].values.astype(float)
X_z = StandardScaler().fit_transform(X_items)
sub["ias_z"] = stats.zscore(sub["ias"])
sub["iats_z"] = stats.zscore(sub["iats"])


def cluster_centroid_zspace(labels, ias_z, iats_z):
    """Return DataFrame of (ias_z, iats_z) centroid by cluster id."""
    return pd.DataFrame(
        {
            "ias_z": pd.Series(ias_z).groupby(labels).mean(),
            "iats_z": pd.Series(iats_z).groupby(labels).mean(),
            "n": pd.Series(labels).value_counts().sort_index(),
        }
    )


def canonical_label_from_centroid(ias_z, iats_z):
    """Map (IAS_z, IATS_z) centroid → manuscript profile label."""
    if iats_z < -0.4 and ias_z > 0:
        return "Efficient"
    if iats_z > 0.5:
        return "Hypervigilant"
    return "Uncertain"


def gower_label_from_iats(iats_z):
    if iats_z < -0.4:
        return "Low Attention"
    if iats_z > 0.6:
        return "High Attention"
    return "Average"


# ===========================================================================
# 1) K-MEANS ON 42 ITEMS — K-SELECTION
# ===========================================================================
print("\n[1] k-means on z-scored 42 items — k-selection")
kmeans_ksel = []
for k in range(2, 9):
    km = KMeans(n_clusters=k, n_init=50, random_state=42).fit(X_z)
    lab = km.labels_
    sil = silhouette_score(X_z, lab, sample_size=2000, random_state=42)
    ch = calinski_harabasz_score(X_z, lab)
    db = davies_bouldin_score(X_z, lab)
    sizes = pd.Series(lab).value_counts().sort_index().tolist()
    kmeans_ksel.append(
        {
            "k": k,
            "silhouette": round(sil, 4),
            "CH": round(ch, 2),
            "DB": round(db, 4),
            "inertia": round(km.inertia_, 2),
            "min_max_ratio": round(min(sizes) / max(sizes), 3),
            "sizes": sizes,
        }
    )
kmeans_ksel_df = pd.DataFrame(kmeans_ksel)
kmeans_ksel_df.to_csv(
    os.path.join(DATA_OUT, "item_level_kmeans_kselection.csv"), index=False
)
print(kmeans_ksel_df.to_string(index=False))

# Fit final k=3
km3 = KMeans(n_clusters=3, n_init=500, random_state=42).fit(X_z)
sub["c42_km3_raw"] = km3.labels_

# Re-label by canonical profile from centroid
cent_km = cluster_centroid_zspace(sub["c42_km3_raw"], sub["ias_z"], sub["iats_z"])
cent_km["canonical"] = cent_km.apply(
    lambda r: canonical_label_from_centroid(r["ias_z"], r["iats_z"]), axis=1
)
km_map = cent_km["canonical"].to_dict()
sub["c42_km3"] = sub["c42_km3_raw"].map(km_map)
print("\nk-means 42-item k=3 centroids and label assignment:")
print(cent_km.round(3).to_string())

# ===========================================================================
# 2) GOWER + PAM ON 42 ITEMS — K-SELECTION
# ===========================================================================
print("\n[2] Gower-PAM on 42 ordinal items — k-selection")
D = gower.gower_matrix(sub[all_items].astype(float).values)
gower_ksel = []
for k in range(2, 7):
    pam = KMedoids(
        n_clusters=k,
        metric="precomputed",
        method="pam",
        init="k-medoids++",
        random_state=42,
        max_iter=500,
    ).fit(D)
    lab = pam.labels_
    sil = silhouette_score(D, lab, metric="precomputed")
    sizes = pd.Series(lab).value_counts().sort_index().tolist()
    gower_ksel.append(
        {
            "k": k,
            "silhouette": round(sil, 4),
            "sizes": sizes,
            "min_max_ratio": round(min(sizes) / max(sizes), 3),
        }
    )
gower_ksel_df = pd.DataFrame(gower_ksel)
gower_ksel_df.to_csv(
    os.path.join(DATA_OUT, "item_level_gower_pam_kselection.csv"), index=False
)
print(gower_ksel_df.to_string(index=False))

# Fit final k=3
pam3 = KMedoids(
    n_clusters=3,
    metric="precomputed",
    method="pam",
    init="k-medoids++",
    random_state=42,
    max_iter=500,
).fit(D)
sub["c42_pam3_raw"] = pam3.labels_
cent_pam = cluster_centroid_zspace(sub["c42_pam3_raw"], sub["ias_z"], sub["iats_z"])
cent_pam["label"] = cent_pam.apply(lambda r: gower_label_from_iats(r["iats_z"]), axis=1)
pam_map = cent_pam["label"].to_dict()
sub["c42_pam3"] = sub["c42_pam3_raw"].map(pam_map)
print("\nGower-PAM k=3 centroids and label assignment:")
print(cent_pam.round(3).to_string())

# ===========================================================================
# 3) OUTCOME ANOVAS — FOR EACH SOLUTION
# ===========================================================================
def anova_table(labels_series, outcomes, prefix):
    rows = []
    levels = sorted(labels_series.dropna().unique())
    for o in outcomes:
        groups = [
            sub.loc[labels_series == lvl, o].dropna().values for lvl in levels
        ]
        f, p = stats.f_oneway(*groups)
        grand = np.concatenate(groups)
        ss_total = ((grand - grand.mean()) ** 2).sum()
        ss_within = sum(((g - g.mean()) ** 2).sum() for g in groups)
        eta2 = (ss_total - ss_within) / ss_total
        row = {
            "outcome": o,
            "F": round(f, 2),
            "p": p,
            "eta2_p": round(eta2, 4),
        }
        for lvl, g in zip(levels, groups):
            row[f"mean_{lvl}"] = round(g.mean(), 2)
            row[f"sd_{lvl}"] = round(g.std(ddof=1), 2)
            row[f"n_{lvl}"] = len(g)
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA_OUT, f"{prefix}.csv"), index=False)
    return out, levels


anova_km, km_levels = anova_table(sub["c42_km3"], outcomes, "item_level_kmeans_k3_anova")
print("\nk-means 42-item k=3 ANOVAs:")
print(anova_km.to_string(index=False))

anova_pam, pam_levels = anova_table(
    sub["c42_pam3"], outcomes, "item_level_gower_pam_k3_anova"
)
print("\nGower-PAM k=3 ANOVAs:")
print(anova_pam.to_string(index=False))

# Published-solution ANOVAs (recomputed here for an apples-to-apples comparison)
anova_orig, orig_levels = anova_table(
    sub["orig_label"], outcomes, "item_level_published_k3_anova"
)
print("\nPublished 2-variable k=3 ANOVAs (recomputed):")
print(anova_orig.to_string(index=False))

# ===========================================================================
# 4) CROSS-TABULATIONS
# ===========================================================================
ct_km_orig = pd.crosstab(sub["c42_km3"], sub["orig_label"], margins=True)
ct_km_orig.to_csv(os.path.join(DATA_OUT, "item_level_crosstab_kmeans_vs_published.csv"))
print("\nCross-tab k-means-42 (rows) vs published (cols):")
print(ct_km_orig.to_string())

ct_pam_orig = pd.crosstab(sub["c42_pam3"], sub["orig_label"], margins=True)
ct_pam_orig.to_csv(
    os.path.join(DATA_OUT, "item_level_crosstab_gower_pam_vs_published.csv")
)
print("\nCross-tab Gower-PAM (rows) vs published (cols):")
print(ct_pam_orig.to_string())

# Diagonal agreement: rows use canonical labels, columns use the long
# "Acc / Att" labels — map cols → canonical for the diagonal computation.
col_to_canonical = {
    "High Accuracy / Low Attention": "Efficient",
    "Average Accuracy / High Attention": "Hypervigilant",
    "Low Accuracy / Average Attention": "Uncertain",
}
diag = 0
ct_no_margin = ct_km_orig.drop(index="All", columns="All")
for col_lbl, canon in col_to_canonical.items():
    if col_lbl in ct_no_margin.columns and canon in ct_no_margin.index:
        diag += ct_no_margin.loc[canon, col_lbl]
total = ct_km_orig.loc["All", "All"]
agreement = diag / total
print(f"\nk-means-42 vs published canonical agreement: {diag}/{total} = {agreement:.3f}")

# ===========================================================================
# 5) BOOTSTRAP STABILITY (k-means 42, k=3)
# ===========================================================================
print("\n[5] Bootstrap Jaccard stability — B=100")


def bootstrap_jaccard(X_z, base_labels, B=100, k=3, rng=None):
    rng = rng or np.random.default_rng(0)
    N = len(X_z)
    per_cluster = {c: [] for c in np.unique(base_labels)}
    base_lab = np.asarray(base_labels)
    for b in range(B):
        idx = rng.choice(N, size=N, replace=True)
        Xb = X_z[idx]
        kmb = KMeans(n_clusters=k, n_init=30, random_state=b).fit(Xb)
        # cluster→cluster best-match Jaccard against base
        for c in np.unique(base_lab):
            base_mask = base_lab == c
            best = 0.0
            for cb in np.unique(kmb.labels_):
                boot_mask = np.zeros(N, dtype=bool)
                boot_mask[idx[kmb.labels_ == cb]] = True
                inter = (base_mask & boot_mask).sum()
                union = (base_mask | boot_mask).sum()
                if union > 0:
                    j = inter / union
                    if j > best:
                        best = j
            per_cluster[c].append(best)
    return per_cluster


# stability for the 42-item k-means
stab_km = bootstrap_jaccard(X_z, sub["c42_km3"].values, B=100, k=3, rng=RNG)
# stability for the published (2-variable) solution
X_2var = StandardScaler().fit_transform(sub[["ias", "iats"]].values)
stab_2var = bootstrap_jaccard(X_2var, sub["orig_label"].values, B=100, k=3, rng=RNG)

stab_rows = []
for label, vals in stab_km.items():
    stab_rows.append({"algorithm": "kmeans_42", "cluster": label, "mean_jaccard": round(float(np.mean(vals)), 3)})
for label, vals in stab_2var.items():
    stab_rows.append({"algorithm": "kmeans_2var", "cluster": label, "mean_jaccard": round(float(np.mean(vals)), 3)})
stab_df = pd.DataFrame(stab_rows)
stab_df.to_csv(os.path.join(DATA_OUT, "item_level_bootstrap_stability.csv"), index=False)
print(stab_df.to_string(index=False))

# ===========================================================================
# 6) SAVE PER-PARTICIPANT LABELS
# ===========================================================================
keep = ["orig_cluster", "orig_label", "ias", "iats", "c42_km3", "c42_pam3"]
if "responseid" in sub.columns:
    keep = ["responseid"] + keep
sub[keep].to_csv(os.path.join(DATA_OUT, "item_level_cluster_labels.csv"), index=False)

# ===========================================================================
# 7) SUPPLEMENTARY FIGURE — 2×3 PANEL LAYOUT
# ===========================================================================
print("\n[7] Building supplementary figure")

fig = plt.figure(figsize=(15.5, 9.5))
gs = fig.add_gridspec(2, 3, hspace=0.36, wspace=0.30,
                      left=0.06, right=0.98, top=0.94, bottom=0.08)


def scatter_panel(ax, labels, palette, title, show_legend=True):
    ias_v = sub["ias"].values
    iats_v = sub["iats"].values
    for name, color in palette.items():
        m = (labels == name).values
        ax.scatter(ias_v[m], iats_v[m], s=12, c=color, alpha=0.55,
                   edgecolors="white", linewidths=0.3, label=f"{name} (n={m.sum()})")
        # 90% confidence ellipse
        if m.sum() >= 10:
            x = ias_v[m]
            y = iats_v[m]
            cov = np.cov(x, y)
            pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
            rx = np.sqrt(1 + pearson)
            ry = np.sqrt(1 - pearson)
            chi2_90 = 4.605  # 2-df 90% quantile
            scale_x = np.sqrt(cov[0, 0] * chi2_90)
            scale_y = np.sqrt(cov[1, 1] * chi2_90)
            ell = Ellipse((0, 0), width=rx * 2, height=ry * 2,
                          facecolor="none", edgecolor=color, linewidth=1.6, alpha=0.9)
            transf = (mtransforms.Affine2D()
                      .rotate_deg(45)
                      .scale(scale_x, scale_y)
                      .translate(x.mean(), y.mean()))
            ell.set_transform(transf + ax.transData)
            ax.add_patch(ell)
    ax.set_xlabel("IAS total")
    ax.set_ylabel("IATS total")
    ax.set_title(title, fontsize=11, weight="bold", loc="left")
    if show_legend:
        ax.legend(fontsize=7.5, loc="lower right", framealpha=0.85)
    ax.grid(alpha=0.25)


# A: k-means 42-item k=3 scatter
axA = fig.add_subplot(gs[0, 0])
scatter_panel(axA, sub["c42_km3"], CANONICAL_COLORS,
              "a  k-means on 42 items (k=3)")

# B: Gower-PAM 42-item k=3 scatter
axB = fig.add_subplot(gs[0, 1])
scatter_panel(axB, sub["c42_pam3"], GOWER_COLORS,
              "b  Gower-PAM on 42 items (k=3)")

# C: k-selection silhouette
axC = fig.add_subplot(gs[0, 2])
ks_km = kmeans_ksel_df["k"].values
sil_km = kmeans_ksel_df["silhouette"].values
ks_pam = gower_ksel_df["k"].values
sil_pam = gower_ksel_df["silhouette"].values
axC.plot(ks_km, sil_km, "o-", color="#333333", label="k-means (42 items)", linewidth=1.8, markersize=6)
axC.plot(ks_pam, sil_pam, "s--", color="#5B9BD5", label="Gower-PAM (42 items)", linewidth=1.8, markersize=6)
axC.axvline(3, linestyle=":", color="#888", linewidth=1)
axC.set_xlabel("k")
axC.set_ylabel("Silhouette width")
axC.set_title("c  Cluster validation by k", fontsize=11, weight="bold", loc="left")
axC.set_xticks(range(2, 9))
axC.legend(fontsize=8, loc="upper right")
axC.grid(alpha=0.25)
axC.text(3.05, axC.get_ylim()[1] * 0.05, "k=3 (primary)", fontsize=7.5, color="#666")

# D: η²p comparison across 3 solutions
axD = fig.add_subplot(gs[1, 0])
eta_pub = anova_orig.set_index("outcome").loc[outcomes, "eta2_p"].values
eta_km = anova_km.set_index("outcome").loc[outcomes, "eta2_p"].values
eta_pam = anova_pam.set_index("outcome").loc[outcomes, "eta2_p"].values
xp = np.arange(len(outcomes))
w = 0.27
axD.bar(xp - w, eta_pub, w, color="#2C3E50", label="Primary (2 summary scores)")
axD.bar(xp, eta_km, w, color="#27AE60", label="k-means (42 items)")
axD.bar(xp + w, eta_pam, w, color="#5B9BD5", label="Gower-PAM (42 items)")
axD.set_xticks(xp)
axD.set_xticklabels([outcome_pretty[o].split(" ")[0] for o in outcomes], rotation=20, ha="right")
axD.set_ylabel(r"Partial $\eta^2$")
axD.set_title("d  Outcome separation by solution", fontsize=11, weight="bold", loc="left")
axD.legend(fontsize=8, loc="upper right", framealpha=0.85)
axD.grid(alpha=0.25, axis="y")

# E: Cross-tab heatmap (k-means 42 vs published)
axE = fig.add_subplot(gs[1, 1])
ct_plot = pd.crosstab(sub["c42_km3"], sub["orig_label"])
# normalize by row
row_norm = ct_plot.div(ct_plot.sum(axis=1), axis=0)
order_rows = [c for c in ["Efficient", "Hypervigilant", "Uncertain"] if c in row_norm.index]
order_cols = [c for c in row_norm.columns if c]  # keep natural order
# manually order cols by canonical mapping
col_canon = col_to_canonical
order_cols = sorted(row_norm.columns, key=lambda c: ["Efficient","Hypervigilant","Uncertain"].index(col_canon.get(c, c)))
row_norm = row_norm.loc[order_rows, order_cols]
im = axE.imshow(row_norm.values, cmap="Greens", vmin=0, vmax=1, aspect="auto")
for i in range(row_norm.shape[0]):
    for j in range(row_norm.shape[1]):
        n = ct_plot.loc[row_norm.index[i], row_norm.columns[j]]
        pct = row_norm.values[i, j]
        axE.text(j, i, f"{n}\n({pct:.0%})", ha="center", va="center",
                 fontsize=9, color="white" if pct > 0.5 else "black")
axE.set_xticks(range(len(order_cols)))
axE.set_xticklabels([col_canon[c] for c in order_cols], rotation=20, ha="right", fontsize=9)
axE.set_yticks(range(len(order_rows)))
axE.set_yticklabels(order_rows, fontsize=9)
axE.set_xlabel("Primary 3-profile solution")
axE.set_ylabel("42-item k-means")
axE.set_title(
    f"e  Recovery of profiles ({diag}/{total} = {agreement:.0%} agreement)",
    fontsize=11, weight="bold", loc="left",
)
plt.colorbar(im, ax=axE, fraction=0.04, pad=0.02, label="row proportion")

# F: outcome means by Gower-PAM cluster (showing monotonic IATS pattern)
axF = fig.add_subplot(gs[1, 2])
pam_order = ["Low Attention", "Average", "High Attention"]
n_levels = len(pam_order)
xp2 = np.arange(len(outcomes))
w2 = 0.27
for i, lvl in enumerate(pam_order):
    means = [anova_pam.set_index("outcome").loc[o, f"mean_{lvl}"] for o in outcomes]
    sds = [anova_pam.set_index("outcome").loc[o, f"sd_{lvl}"] for o in outcomes]
    ns = [anova_pam.set_index("outcome").loc[o, f"n_{lvl}"] for o in outcomes]
    sems = [s / np.sqrt(n) for s, n in zip(sds, ns)]
    axF.bar(xp2 + (i - 1) * w2, means, w2, color=GOWER_COLORS[lvl],
            yerr=sems, capsize=2, label=f"{lvl} (n={int(ns[0])})", alpha=0.95)
axF.set_xticks(xp2)
axF.set_xticklabels([outcome_pretty[o].split(" ")[0] for o in outcomes],
                    rotation=20, ha="right")
axF.set_ylabel("Mean outcome score")
axF.set_title("f  Gower-PAM: monotonic IATS gradient", fontsize=11, weight="bold", loc="left")
axF.legend(fontsize=8, loc="upper left", framealpha=0.85)
axF.grid(alpha=0.25, axis="y")

fig_path = os.path.join(PLOTS_DIR, "1_s_9_item_level_clustering_sensitivity.png")
fig.savefig(fig_path, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Figure saved to {fig_path}")

# ===========================================================================
# 8) MACHINE-READABLE SUMMARY
# ===========================================================================
summary = {
    "n": int(N),
    "kmeans_42item": {
        "k_selection": kmeans_ksel_df.to_dict(orient="records"),
        "k3_centroids": cent_km.reset_index().to_dict(orient="records"),
        "k3_anova": anova_km.to_dict(orient="records"),
        "k3_crosstab_vs_published": ct_km_orig.to_dict(),
        "k3_diagonal_agreement": float(agreement),
        "k3_bootstrap_jaccard_mean": {
            str(k): float(np.mean(v)) for k, v in stab_km.items()
        },
    },
    "gower_pam_42item": {
        "k_selection": gower_ksel_df.to_dict(orient="records"),
        "k3_centroids": cent_pam.reset_index().to_dict(orient="records"),
        "k3_anova": anova_pam.to_dict(orient="records"),
        "k3_crosstab_vs_published": ct_pam_orig.to_dict(),
    },
    "primary_solution_recomputed": {
        "k3_anova": anova_orig.to_dict(orient="records"),
        "k3_bootstrap_jaccard_mean": {
            str(k): float(np.mean(v)) for k, v in stab_2var.items()
        },
    },
}
with open(os.path.join(DATA_OUT, "item_level_summary.json"), "w") as f:
    json.dump(summary, f, indent=2, default=str)
print(f"\nSummary JSON: {os.path.join(DATA_OUT, 'item_level_summary.json')}")
print("\nDone.")
