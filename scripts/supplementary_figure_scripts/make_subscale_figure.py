"""Supplementary figure: TAS-20 subscale decomposition of the alexithymia pathway.

Draws the PARALLEL multiple-mediator model (DIF, DDF, EOT and somatic distress
entered simultaneously), i.e. the same model whose numbers appear in
Supplementary Table 5 and in the Results. The values come from
supplementary_data/subscale_parallel_model.csv, written by
19_ert_exact_statistics.R (2,000 percentile bootstrap resamples).

Run after 19 has run, from anywhere:
  python scripts/supplementary_figure_scripts/make_subscale_figure.py
or through the pipeline runner:
  python scripts/run_analysis.py s18
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Find the repository root by looking for the data file, walking up from this
# file. Counting parent directories breaks the moment the script moves, and it
# has: this folder sits at the repository root in the working repo and under
# scripts/ in the public release.
ROOT = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(ROOT, "dfc_interoception_profiling.csv")):
    parent = os.path.dirname(ROOT)
    if parent == ROOT:
        raise SystemExit(
            "Cannot find dfc_interoception_profiling.csv above "
            + os.path.abspath(__file__)
            + ". Run this from a complete checkout of the repository.")
    ROOT = parent
SRC = os.path.join(ROOT, "supplementary_data", "subscale_parallel_model.csv")
OUT = os.path.join(ROOT, "plots", "supplementary_plots", "3_s_8_tas_subscale_decomposition.png")

d = pd.read_csv(SRC)
# Sample size for the caption, read from the deposited data rather than
# hard-coded, so it cannot drift when the sample changes.
N = len(pd.read_csv(os.path.join(ROOT, "dfc_interoception_profiling.csv")))
paths = d[d.block == "path"].reset_index(drop=True)
ind = d[d.block == "indirect"].reset_index(drop=True)
share = d[d.block == "share"].reset_index(drop=True)

COL = {"DIF": "#E74C3C", "DDF": "#3498DB", "EOT": "#1ABC9C", "Somatic": "#8D6E63", "TAS": "#333333"}

def colour(label):
    for k in ("DIF", "DDF", "EOT", "Somatic"):
        if k in label:
            return COL[k]
    return COL["TAS"]

fig = plt.figure(figsize=(13, 9.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], width_ratios=[1.15, 1])

# ---- A: a-paths and b-paths -------------------------------------------------
ax = fig.add_subplot(gs[0, :])
order = list(paths.index)[::-1]
for y, i in enumerate(order):
    r = paths.loc[i]
    c = colour(r.label)
    sig = r.p < 0.05
    ax.errorbar(r.est, y, xerr=[[r.est - r.ci_lo], [r.ci_hi - r.est]], fmt="o", color=c,
                mfc=c if sig else "white", mec=c, ms=7, capsize=3, lw=1.8)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([paths.loc[i, "label"].replace("->", "→") for i in order])
ax.axvline(0, color="0.4", ls="--", lw=1)
ax.set_xlabel("Standardized coefficient [95% percentile bootstrap CI]")
ax.set_title("A. Path coefficients in the parallel four-mediator model", loc="left", fontweight="bold")
ax.grid(axis="x", alpha=0.3)

# ---- B: indirect effects -----------------------------------------------------
ax = fig.add_subplot(gs[1, 0])
order = list(ind.index)[::-1]
for y, i in enumerate(order):
    r = ind.loc[i]
    c = colour(r.label)
    sig = r.p < 0.05
    ax.errorbar(r.est, y, xerr=[[r.est - r.ci_lo], [r.ci_hi - r.est]], fmt="o", color=c,
                mfc=c if sig else "white", mec=c, ms=7, capsize=3, lw=1.8)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([ind.loc[i, "label"].replace("->", "→") for i in order])
ax.axvline(0, color="0.4", ls="--", lw=1)
ax.set_xlabel("Indirect effect of IAS on mental health [95% percentile bootstrap CI]")
ax.set_title("B. Indirect effects", loc="left", fontweight="bold")
ax.grid(axis="x", alpha=0.3)

# ---- C: shares ---------------------------------------------------------------
ax = fig.add_subplot(gs[1, 1])
xs = range(len(share))
ax.bar(xs, share.est, color=[COL[l] for l in share.label], edgecolor="black", lw=0.5)
ax.errorbar(xs, share.est, yerr=[share.est - share.ci_lo, share.ci_hi - share.est],
            fmt="none", ecolor="black", capsize=4)
for x, v in zip(xs, share.est):
    ax.text(x, v + 4, f"{v:.1f}%", ha="center", fontweight="bold")
ax.set_xticks(list(xs))
ax.set_xticklabels(share.label)
ax.set_ylim(0, 100)
ax.set_ylabel("Share of the alexithymia-mediated indirect effect (%)")
ax.set_title("C. Facet shares of the alexithymia pathway", loc="left", fontweight="bold")
ax.grid(axis="y", alpha=0.3)

fig.text(0.01, 0.005, f"Filled markers: 95% bootstrap CI excludes zero. N = {N}; 2,000 percentile bootstrap resamples.",
         fontsize=9, color="0.3")
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
print("wrote", OUT)
