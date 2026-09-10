"""Supplementary figure: dual-pathway coefficients by definition of the somatic mediator.

Same communicative pattern as the TAS-20 subscale figure: show the alternatives
side by side and let the reader compare, rather than enumerating models in prose.
Values are the standardized solutions in analysis_output/17_output.txt (composite,
SSS-8 alone, age+gender adjusted) and 18_output.txt (PCS alone), regenerated
2026-09-02 after the score corrections in 00_score_corrections.py.
Run from anywhere: python scripts/supplementary_figure_scripts/make_somatic_mediator_figure.py
or through the pipeline runner: python scripts/run_analysis.py s24
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# path: {model: (beta, lo, hi)}
PATHS = [
    ("IAS \u2192 somatic", {
        "SSS-8 + PCS": (-0.094, -0.154, -0.033),
        "SSS-8 alone": (-0.099, -0.160, -0.037),
        "PCS alone":   (-0.061, -0.125,  0.003)}),
    ("IATS \u2192 somatic", {
        "SSS-8 + PCS": (0.444, 0.390, 0.497),
        "SSS-8 alone": (0.409, 0.354, 0.465),
        "PCS alone":   (0.346, 0.287, 0.406)}),
    ("IAS \u00d7 IATS \u2192 somatic", {
        "SSS-8 + PCS": (-0.094, -0.156, -0.032),
        "SSS-8 alone": (-0.068, -0.131, -0.006),
        "PCS alone":   (-0.092, -0.156, -0.027)}),
    ("somatic \u2192 mental health", {
        "SSS-8 + PCS": (0.611, 0.567, 0.655),
        "SSS-8 alone": (0.554, 0.508, 0.600),
        "PCS alone":   (0.388, 0.336, 0.440)}),
]

# covariate adjustment, interaction path only
COV = [
    ("SSS-8 + PCS, unadjusted",        (-0.094, -0.156, -0.032)),
    ("SSS-8 + PCS, age and gender",    (-0.088, -0.149, -0.027)),
    ("SSS-8 alone, unadjusted",        (-0.068, -0.131, -0.006)),
    ("SSS-8 alone, age and gender",    (-0.063, -0.125, -0.001)),
]

MODELS = ["SSS-8 + PCS", "SSS-8 alone", "PCS alone"]
COLS = {"SSS-8 + PCS": "#1b4965", "SSS-8 alone": "#5fa8d3", "PCS alone": "#e07a5f"}

fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.6),
                         gridspec_kw={"width_ratios": [1.5, 1]})

# ---- panel A ----
ax = axes[0]
ypos, ylabels = [], []
y = 0.0
for label, models in PATHS:
    for m in MODELS:
        b, lo, hi = models[m]
        ax.errorbar(b, y, xerr=[[b - lo], [hi - b]], fmt="o", ms=6, capsize=3,
                    color=COLS[m], lw=1.8,
                    label=m if label == PATHS[0][0] else None)
        y -= 1.0
    ypos.append(y + 2.0)
    ylabels.append(label)
    y -= 0.8
ax.axvline(0, color="0.4", lw=1, ls="--")
ax.set_yticks(ypos)
ax.set_yticklabels(ylabels)
ax.set_xlabel("Standardized coefficient (95% CI)")
ax.set_title("A. Path coefficients by definition of the somatic mediator", loc="left",
             fontsize=11, fontweight="bold")
ax.legend(frameon=False, loc="lower left", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

# ---- panel B ----
ax = axes[1]
for i, (label, (b, lo, hi)) in enumerate(COV):
    col = "#1b4965" if "SSS-8 + PCS" in label else "#5fa8d3"
    ax.errorbar(b, -i, xerr=[[b - lo], [hi - b]], fmt="o", ms=6, capsize=3,
                color=col, lw=1.8)
ax.axvline(0, color="0.4", lw=1, ls="--")
ax.set_yticks([-i for i in range(len(COV))])
ax.set_yticklabels([l for l, _ in COV], fontsize=9)
ax.set_xlabel("IAS \u00d7 IATS \u2192 somatic (95% CI)")
ax.set_title("B. Interaction under demographic adjustment", loc="left",
             fontsize=11, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
import os
# Same root search as make_subscale_figure.py. Counting parent directories was
# correct in only one of the two layouts this file lives in.
ROOT = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(ROOT, "dfc_interoception_profiling.csv")):
    parent = os.path.dirname(ROOT)
    if parent == ROOT:
        raise SystemExit(
            "Cannot find dfc_interoception_profiling.csv above "
            + os.path.abspath(__file__)
            + ". Run this from a complete checkout of the repository.")
    ROOT = parent
out = os.path.join(ROOT, "plots", "supplementary_plots",
                   "1_s_24_somatic_mediator_definitions.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
print("wrote", out)
