"""Build the per-figure source data workbooks the journal requires.

Writes source_data/Source_Data_Figure_1.xlsx, _2.xlsx and _3.xlsx: one sheet per
panel, holding the numbers behind that panel, taken from the files the figure
scripts read or write. Nothing is re-estimated here except three derived
quantities that the figure scripts compute inline and never save:

  Figure 2A  p-values for the correlation heatmap (t-test on r, n = 832), which
             03_correlation_heatmaps.py computes but does not write
  Figure 2C-H per-profile summary statistics of the box plots, and the Tukey HSD
             contrasts behind the significance markers, computed the same way
             04_cluster_analysis.R computes them
  Figure 3D  the classical MDS coordinates, obtained by calling R's cmdscale on
             the deposited correlation matrix, so they equal the plotted layout

Run from the repository root after the pipeline:
  python scripts/make_source_data.py
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
while not (ROOT / "dfc_interoception_profiling.csv").exists():
    if ROOT.parent == ROOT:
        sys.exit("cannot find the repository root")
    ROOT = ROOT.parent
SD = ROOT / "supplementary_data"
OUT = ROOT / "source_data"
OUT.mkdir(exist_ok=True)

PROFILE = {1: "Uncertain", 2: "Efficient", 3: "Hypervigilant"}   # cluster number -> label, from table_4_cluster_profiles.csv
OUTCOMES = [("tas", "C", "Alexithymia (TAS-20)"), ("sss8", "D", "Somatic symptoms (SSS-8)"),
            ("pcs", "E", "Pain catastrophizing (PCS)"), ("phq9", "F", "Depression (PHQ-9)"),
            ("gad7", "G", "Anxiety (GAD-7)"), ("stai", "H", "Trait anxiety (STAI-T)")]


def readme(rows):
    return pd.DataFrame(rows, columns=["sheet", "panel", "content", "taken from", "produced by"])


def write(path, sheets):
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in sheets:
            df.to_excel(xw, sheet_name=name[:31], index=False)
            ws = xw.sheets[name[:31]]
            for col in ws.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in col)
                ws.column_dimensions[col[0].column_letter].width = min(60, max(10, width + 2))
    print(f"wrote {path.relative_to(ROOT)}  ({len(sheets)} sheets)")


# ------------------------------------------------------------------ Figure 1
def figure_1():
    x = pd.ExcelFile(SD / "sem_pathway_results.xlsx")
    coef = x.parse("SEM_Coefficients")
    indirect = x.parse("Indirect_Effects")
    dom = pd.read_csv(SD / "pathway_dominance_by_profile.csv")
    slopes = pd.read_csv(SD / "figure1_simple_slopes_by_profile.csv")
    sheets = [
        ("README", readme([
            ("1A_1B_path_coefficients", "a, b", "Standardized path coefficients and p-values of the dual-pathway SEM, one row per model specification; the figure shows the composite-outcome model", "supplementary_data/sem_pathway_results.xlsx, sheet SEM_Coefficients", "12_sem_pathway_analysis.R"),
            ("1A_1B_indirect_effects", "a, b", "Indirect effects through each mediator", "sem_pathway_results.xlsx, sheet Indirect_Effects", "12_sem_pathway_analysis.R"),
            ("1C_pathway_dominance", "c", "Percentage of the IAS indirect effect carried by each pathway per profile, with bootstrap 95% CIs (5,000 resamples) and tests", "supplementary_data/pathway_dominance_by_profile.csv", "12b_pathway_dominance_test.R"),
            ("1D_somatic_slopes", "d", "Simple slope of IAS on the somatic mediator at each profile's mean IATS, with SE, 95% CI, t, df, p", "supplementary_data/figure1_simple_slopes_by_profile.csv, outcome = somatic", "19_ert_exact_statistics.R"),
            ("1E_alexithymia_slopes", "e", "Simple slope of IAS on alexithymia at each profile's mean IATS", "figure1_simple_slopes_by_profile.csv, outcome = tas", "19_ert_exact_statistics.R"),
        ])),
        ("1A_1B_path_coefficients", coef),
        ("1A_1B_indirect_effects", indirect),
        ("1C_pathway_dominance", dom),
        ("1D_somatic_slopes", slopes[slopes.outcome == "somatic"].reset_index(drop=True)),
        ("1E_alexithymia_slopes", slopes[slopes.outcome == "tas"].reset_index(drop=True)),
    ]
    write(OUT / "Source_Data_Figure_1.xlsx", sheets)


# ------------------------------------------------------------------ Figure 2
def figure_2():
    n = 832
    cm = pd.read_csv(SD / "correlation_matrix_intero_mh.csv", index_col=0)
    r = cm.values.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = r * np.sqrt((n - 2) / (1 - r ** 2))
    p = 2 * stats.t.sf(np.abs(t), df=n - 2)
    np.fill_diagonal(p, np.nan)
    pm = pd.DataFrame(p, index=cm.index, columns=cm.columns)
    cm_out = cm.reset_index().rename(columns={"index": "variable", "Unnamed: 0": "variable"})
    pm_out = pm.reset_index().rename(columns={"index": "variable"})

    d = pd.read_csv(SD / "data_with_k_means_clusters.csv")
    d["profile"] = d["cluster"].map(PROFILE)
    # Both label forms: the manuscript's profile name and the descriptive label
    # the clustering script writes, so the correspondence is explicit.
    scatter = d[["ias", "iats", "cluster", "cluster_label", "profile"]].copy()
    scatter.insert(0, "participant", range(1, len(d) + 1))
    prof = pd.read_csv(SD / "table_4_cluster_profiles.csv")
    prof.insert(1, "profile", prof["cluster"].map(PROFILE))

    indiv = d[["cluster", "cluster_label", "profile"] + [o for o, _, _ in OUTCOMES]].copy()
    indiv.insert(0, "participant", range(1, len(d) + 1))

    summ_rows, tukey_rows = [], []
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    for var, panel, label in OUTCOMES:
        g = d.groupby("profile")[var]
        for prof_name, s in g:
            q1, med, q3 = s.quantile([0.25, 0.5, 0.75])
            iqr = q3 - q1
            lo_w = s[s >= q1 - 1.5 * iqr].min()
            hi_w = s[s <= q3 + 1.5 * iqr].max()
            summ_rows.append([panel, label, prof_name, int(s.count()), round(s.mean(), 3),
                              round(med, 3), round(q1, 3), round(q3, 3), round(lo_w, 3), round(hi_w, 3),
                              int(((s < lo_w) | (s > hi_w)).sum())])
        tk = pairwise_tukeyhsd(d[var], d["profile"], alpha=0.05)
        for row in tk.summary().data[1:]:
            tukey_rows.append([panel, label] + list(row))
    summ = pd.DataFrame(summ_rows, columns=["panel", "outcome", "profile", "n", "mean", "median",
                                            "q1", "q3", "whisker_low", "whisker_high", "points_beyond_whiskers"])
    tukey = pd.DataFrame(tukey_rows, columns=["panel", "outcome", "group1", "group2", "mean_diff",
                                              "p_adj", "ci_lower", "ci_upper", "reject_at_0.05"])
    anova = pd.read_csv(SD / "table_4_anova_results.csv")

    sheets = [
        ("README", readme([
            ("2A_correlation_r", "a", "Pearson correlations among the eight scale totals (lower triangle plotted)", "supplementary_data/correlation_matrix_intero_mh.csv", "03_correlation_heatmaps.py"),
            ("2A_correlation_p", "a", "Two-sided p-values for those correlations, t-test on r with n = 832, which set the cell colours", "computed here from r and n, as 03_correlation_heatmaps.py does", "03_correlation_heatmaps.py"),
            ("2B_scatter_individuals", "b", "IAS and IATS score and k-means profile of every participant; the 90% ellipses are drawn from these points", "supplementary_data/data_with_k_means_clusters.csv", "04_cluster_analysis.R"),
            ("2B_profile_centroids", "b", "Profile sizes and mean scores", "supplementary_data/table_4_cluster_profiles.csv", "04_cluster_analysis.R"),
            ("2C-2H_individual_values", "c to h", "Every participant's six outcome scores and profile, from which each box plot is drawn", "supplementary_data/data_with_k_means_clusters.csv", "04_cluster_analysis.R"),
            ("2C-2H_box_summary", "c to h", "Per profile: n, mean (white diamond), median, quartiles, whisker ends, points beyond whiskers", "computed here from the individual values", "make_source_data.py"),
            ("2C-2H_tukey", "c to h", "Tukey HSD pairwise contrasts behind the significance markers", "computed here as in 04_cluster_analysis.R (TukeyHSD on the profile factor)", "04_cluster_analysis.R"),
            ("2C-2H_anova", "c to h", "One-way ANOVA per outcome", "supplementary_data/table_4_anova_results.csv", "04_cluster_analysis.R"),
        ])),
        ("2A_correlation_r", cm_out),
        ("2A_correlation_p", pm_out),
        ("2B_scatter_individuals", scatter),
        ("2B_profile_centroids", prof),
        ("2C-2H_individual_values", indiv),
        ("2C-2H_box_summary", summ),
        ("2C-2H_tukey", tukey),
        ("2C-2H_anova", anova),
    ]
    write(OUT / "Source_Data_Figure_2.xlsx", sheets)


# ------------------------------------------------------------------ Figure 3
def r_cmdscale(dist_csv):
    """Classical MDS exactly as the figure script does it: R's cmdscale(k = 2)."""
    rscript = os.environ.get("RSCRIPT", "Rscript")
    out = Path(tempfile.mkdtemp()) / "mds.csv"
    code = (f'm <- as.matrix(read.csv("{dist_csv.as_posix()}", row.names = 1, check.names = FALSE));'
            'd <- as.dist(1 - abs(m)); z <- cmdscale(d, k = 2);'
            f'write.csv(data.frame(item = rownames(m), x = z[,1], y = z[,2]), "{out.as_posix()}", row.names = FALSE)')
    r = subprocess.run([rscript, "--vanilla", "-e", code], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("Rscript failed for cmdscale:\n" + r.stderr[-800:])
    return pd.read_csv(out)


def figure_3():
    cons = pd.read_csv(SD / "item_consensus_ranking.csv")
    cons["prediction_importance"] = (43 - cons["avg_rank"]) / 42
    cons["consensus_rank"] = range(1, len(cons) + 1)
    cons["shown_in_panel_a"] = cons["consensus_rank"] <= 21
    cent = pd.read_csv(SD / "network_centrality.csv")
    sc = cons.merge(cent, on="item", how="left")
    sc["network_centrality_norm"] = sc.groupby("scale")["strength"].transform(
        lambda s: (s - s.min()) / (s.max() - s.min()))
    sat = pd.read_csv(SD / "saturation_results.csv")
    cfg = pd.read_csv(SD / "minimal_set_comparison.csv")
    full_r2 = float(cfg.loc[cfg["config"] == "Full scales (totals)", "cv_r2"].iloc[0])
    sat_lines = pd.DataFrame([["horizontal dashed line: cross-validated R2 of the full-scale totals", full_r2],
                              ["vertical dashed line: 10 items", 10],
                              ["vertical dashed line: 21 items", 21]], columns=["reference line", "value"])

    cm = pd.read_csv(SD / "full_item_correlation_matrix.csv")
    first = cm.columns[0]
    cm = cm.set_index(first)
    cm.index.name = "item"
    mh = pd.read_csv(SD / "item_mh_correlations.csv")
    top21 = set(cons["item"].head(21))
    top10 = set(cons["item"].head(10))
    tmp = Path(tempfile.mkdtemp()) / "cor.csv"
    cm.to_csv(tmp)
    mds = r_cmdscale(tmp)
    nodes = mds.merge(mh[["item", "scale", "domain", "correlation_with_mh_composite",
                          "abs_correlation_with_mh"]], on="item", how="left")
    nodes["black_ring_top21"] = nodes["item"].isin(top21)
    nodes["top10_predictor"] = nodes["item"].isin(top10)
    items = list(cm.index)
    R = cm.values.astype(float)
    edges = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if abs(R[i, j]) >= 0.30:
                edges.append([items[i], items[j], round(R[i, j], 6), round(abs(R[i, j]), 6)])
    edges = pd.DataFrame(edges, columns=["item_1", "item_2", "r", "abs_r"])

    sheets = [
        ("README", readme([
            ("3A_item_importance", "a", "Consensus rank across LASSO, random forest and network centrality for all 42 interoceptive items; importance = (43 - mean rank) / 42; the panel shows the top 21; dashed line at 0.65", "supplementary_data/item_consensus_ranking.csv", "09_item_selection_analysis.R, drawn by 11_figure_combined.R"),
            ("3B_saturation_curve", "b", "Cross-validated R2 as items are added by consensus rank", "supplementary_data/saturation_results.csv", "09_item_selection_analysis.R"),
            ("3B_reference_lines", "b", "Values of the dashed reference lines", "supplementary_data/minimal_set_comparison.csv", "09_item_selection_analysis.R"),
            ("3C_centrality_vs_importance", "c", "Per item: prediction importance and strength centrality, with centrality min-max normalised within scale as plotted", "item_consensus_ranking.csv joined with network_centrality.csv", "11_figure_combined.R"),
            ("3D_nodes", "d", "All 119 items: classical MDS coordinates on 1 - |r| (R cmdscale, k = 2), correlation with the mental health composite (node size), and ring flags", "full_item_correlation_matrix.csv, item_mh_correlations.csv", "10_network_analysis.R"),
            ("3D_edges", "d", "Every plotted edge: item pairs with |r| >= 0.30 and their correlation", "full_item_correlation_matrix.csv", "10_network_analysis.R"),
            ("3D_correlation_matrix", "d", "The full 119 x 119 item correlation matrix the panel is built from", "supplementary_data/full_item_correlation_matrix.csv", "10_network_analysis.R"),
        ])),
        ("3A_item_importance", cons[["consensus_rank", "item", "scale", "lasso_rank", "rf_rank", "centrality_rank",
                                     "avg_rank", "prediction_importance", "shown_in_panel_a"]]),
        ("3B_saturation_curve", sat),
        ("3B_reference_lines", sat_lines),
        ("3C_centrality_vs_importance", sc[["item", "scale", "prediction_importance", "strength", "betweenness",
                                            "combined_centrality", "network_centrality_norm"]]),
        ("3D_nodes", nodes),
        ("3D_edges", edges),
        ("3D_correlation_matrix", cm.reset_index()),
    ]
    write(OUT / "Source_Data_Figure_3.xlsx", sheets)


if __name__ == "__main__":
    figure_1()
    figure_2()
    figure_3()
    print("done")
