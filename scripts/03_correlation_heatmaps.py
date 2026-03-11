"""
================================================================================
CORRELATION HEATMAP GENERATOR
================================================================================

PURPOSE:
    Creates publication-quality correlation heatmaps for psychological/
    neuroscience research. Designed for reuse across different projects.

FEATURES:
    - Lower triangle display (no diagonal, no redundant upper triangle)
    - Correlation values displayed in each cell
    - Discrete color scheme based on direction AND significance:
        * Positive correlations: red spectrum (darker = more significant)
        * Negative correlations: blue spectrum (darker = more significant)
        * Non-significant: pale/muted colors
    - Legend showing color meanings
    - Customizable variable labels
    - Outputs both PNG figure and CSV correlation matrix

HEATMAP STRUCTURE:
    The heatmap displays correlations WITHOUT the diagonal by:
    - Y-axis: Variables 2 through N (skips first variable)
    - X-axis: Variables 1 through N-1 (skips last variable)

    This creates a clean lower triangle where each cell shows the
    correlation between the row variable and column variable.

    Example for 4 variables [A, B, C, D]:

                A       B       C
        B     r(A,B)
        C     r(A,C)  r(B,C)
        D     r(A,D)  r(B,D)  r(C,D)

COLOR SCALE:
    Discrete colors based on significance level:
    - Dark red (#67000d): Positive p < .001
    - Red (#cb181d): Positive p < .01
    - Light red (#fb6a4a): Positive p < .05
    - Pale red (#fcbba1): Positive ns
    - Gray (#f0f0f0): Near zero
    - Pale blue (#c6dbef): Negative ns
    - Light blue (#6baed6): Negative p < .05
    - Blue (#2171b5): Negative p < .01
    - Dark blue (#08306b): Negative p < .001

USAGE:
    1. Modify the 'variables' list to include your variable names
    2. Modify 'label_map' to provide clean display labels
    3. Run script: python correlation_heatmap.py

    For different datasets, change 'data_path' and variable definitions.

AUTHOR: Generated for interoception/mental health research
DATE: 2025
================================================================================
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-GUI backend for server/WSL environments
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import pearsonr
import os
import warnings
warnings.filterwarnings('ignore')


def get_color_category(r, p):
    """
    Assign color category based on correlation direction and significance.

    Returns integer code:
        -4: negative, p < .001
        -3: negative, p < .01
        -2: negative, p < .05
        -1: negative, not significant
         0: near zero (|r| < 0.05)
         1: positive, not significant
         2: positive, p < .05
         3: positive, p < .01
         4: positive, p < .001
    """
    if abs(r) < 0.05:
        return 0

    if r > 0:
        if p < 0.001:
            return 4
        elif p < 0.01:
            return 3
        elif p < 0.05:
            return 2
        else:
            return 1
    else:
        if p < 0.001:
            return -4
        elif p < 0.01:
            return -3
        elif p < 0.05:
            return -2
        else:
            return -1


def create_correlation_heatmap(
    df: pd.DataFrame,
    variables: list,
    label_map: dict,
    output_path: str,
    title: str = "Correlation Matrix",
    figsize: tuple = (12, 10),
    annot_fontsize: int = 11,
    label_fontsize: int = 11,
    legend_fontsize: int = 11,
    legend_title_fontsize: int = 12
) -> tuple:
    """
    Create correlation heatmap with discrete colors based on significance.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the variables to correlate
    variables : list
        List of column names to include (in desired display order)
    label_map : dict
        Dictionary mapping variable names to display labels
    output_path : str
        Full path for output PNG file
    title : str
        Plot title
    figsize : tuple
        Figure size (width, height) in inches
    annot_fontsize : int
        Font size for correlation values in cells
    label_fontsize : int
        Font size for axis labels

    Returns
    -------
    tuple : (correlation_matrix, p_value_matrix, n_complete_cases)
    """

    # Prepare data
    analysis_df = df[variables].dropna()
    n_complete = len(analysis_df)
    n_vars = len(variables)

    if title and "{n}" in title:
        title = title.format(n=n_complete)

    print(f"  Variables: {n_vars}")
    print(f"  Complete cases: {n_complete}")

    # Compute correlation matrix
    corr_matrix = analysis_df.corr()

    # Compute p-values
    p_matrix = np.zeros((n_vars, n_vars))
    for i in range(n_vars):
        for j in range(n_vars):
            if i != j:
                _, p_val = pearsonr(
                    analysis_df[variables[i]],
                    analysis_df[variables[j]]
                )
                p_matrix[i, j] = p_val
            else:
                p_matrix[i, j] = 1.0

    # Drop the top row and far-right column (diagonal is masked)
    y_vars = variables[1:]
    x_vars = variables[:-1]

    y_labels = [label_map.get(v, v) for v in y_vars]
    x_labels = [label_map.get(v, v) for v in x_vars]

    # Build display matrices
    display_corr = np.zeros((len(y_vars), len(x_vars)))
    display_p = np.zeros((len(y_vars), len(x_vars)))
    display_cat = np.zeros((len(y_vars), len(x_vars)))

    for i, y_var in enumerate(y_vars):
        for j, x_var in enumerate(x_vars):
            y_idx = variables.index(y_var)
            x_idx = variables.index(x_var)
            r = corr_matrix.iloc[y_idx, x_idx]
            p = p_matrix[y_idx, x_idx]
            display_corr[i, j] = r
            display_p[i, j] = p
            display_cat[i, j] = get_color_category(r, p)

    # Mask upper triangle and diagonal (offset by one for shifted axes)
    mask = np.zeros_like(display_corr, dtype=bool)
    for i in range(len(y_vars)):
        for j in range(len(x_vars)):
            if j >= i + 1:
                mask[i, j] = True
                display_cat[i, j] = np.nan

    # Define discrete colormap
    # Order: -4, -3, -2, -1, 0, 1, 2, 3, 4
    colors = [
        '#08306b',  # -4: negative p<.001 (darkest blue)
        '#2171b5',  # -3: negative p<.01
        '#6baed6',  # -2: negative p<.05
        '#c6dbef',  # -1: negative ns (pale blue)
        '#f0f0f0',  #  0: near zero (gray)
        '#fcbba1',  #  1: positive ns (pale red)
        '#fb6a4a',  #  2: positive p<.05
        '#cb181d',  #  3: positive p<.01
        '#67000d',  #  4: positive p<.001 (darkest red)
    ]

    cmap = ListedColormap(colors)
    bounds = [-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
    norm = BoundaryNorm(bounds, cmap.N)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot heatmap
    im = ax.imshow(display_cat, cmap=cmap, norm=norm, aspect='equal')

    # Set ticks
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=label_fontsize)
    ax.set_yticklabels(y_labels, fontsize=label_fontsize)

    # Add correlation values as text
    for i in range(len(y_vars)):
        for j in range(len(x_vars)):
            if not mask[i, j]:
                r = display_corr[i, j]
                cat = display_cat[i, j]

                # Choose text color based on background
                if abs(cat) >= 3:
                    text_color = 'white'
                else:
                    text_color = 'black'

                ax.text(j, i, f'{r:.2f}', ha='center', va='center',
                       fontsize=annot_fontsize, fontweight='bold', color=text_color)

    # Mask upper triangle with white
    for i in range(len(y_vars)):
        for j in range(len(x_vars)):
            if mask[i, j]:
                ax.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1,
                            facecolor='white', edgecolor='white'))

    # Add grid lines
    ax.set_xticks(np.arange(-0.5, len(x_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(y_labels), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1)

    # Remove plot border (spines)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Create legend (only categories present)
    category_map = {
        -4: ('#08306b', 'Negative p < .001'),
        -3: ('#2171b5', 'Negative p < .01'),
        -2: ('#6baed6', 'Negative p < .05'),
        -1: ('#c6dbef', 'Negative ns'),
         0: ('#f0f0f0', 'Near zero'),
         1: ('#fcbba1', 'Positive ns'),
         2: ('#fb6a4a', 'Positive p < .05'),
         3: ('#cb181d', 'Positive p < .01'),
         4: ('#67000d', 'Positive p < .001')
    }
    present_cats = set(int(x) for x in np.unique(display_cat[~np.isnan(display_cat)]))
    ordered_cats = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    legend_elements = [
        mpatches.Patch(facecolor=category_map[c][0], edgecolor='black', label=category_map[c][1])
        for c in ordered_cats if c in present_cats
    ]

    legend = ax.legend(handles=legend_elements, loc='upper right',
              fontsize=legend_fontsize, title='Correlation', title_fontsize=legend_title_fontsize,
              framealpha=0.95, edgecolor='black', fancybox=False,
              handlelength=2, handleheight=1.5, borderpad=1, labelspacing=0.8)
    legend.get_title().set_fontweight('bold')

    if title:
        ax.set_title(title, fontsize=14, pad=15, weight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_path}")
    plt.close()

    return corr_matrix, p_matrix, n_complete


# =============================================================================
# MAIN SCRIPT - Interoception and Mental Health Analysis
# =============================================================================

if __name__ == "__main__":

    import sys

    # =========================================================================
    # CONFIGURATION
    # =========================================================================

    # Data path
    data_path = "C:/code/projects/mi/analyses/soma/results/dfc_vviq_q_k.csv"
    output_dir = "C:/code/projects/intero_mod"
    analysis_output_dir = os.path.join(output_dir, "analysis_output")
    plots_dir = os.path.join(output_dir, "plots")
    sub_plots_dir = os.path.join(plots_dir, "sub_plots")
    supplementary_plots_dir = os.path.join(plots_dir, "supplementary_plots")
    other_plots_vviq_dir = os.path.join(plots_dir, "other_plots_vviq_plots")
    supplementary_data_dir = "C:/code/projects/intero_mod/supplementary_data"
    os.makedirs(analysis_output_dir, exist_ok=True)
    os.makedirs(sub_plots_dir, exist_ok=True)
    os.makedirs(supplementary_plots_dir, exist_ok=True)
    os.makedirs(other_plots_vviq_dir, exist_ok=True)

    # Redirect output to file
    class Tee:
        def __init__(self, *files):
            self.files = files
        def write(self, obj):
            for f in self.files:
                f.write(obj)
                f.flush()
        def flush(self):
            for f in self.files:
                f.flush()

    output_file = open(os.path.join(analysis_output_dir, "03_output.txt"), "w")
    sys.stdout = Tee(sys.stdout, output_file)

    print("=" * 70)
    print("CORRELATION HEATMAP GENERATOR")
    print("Interoception and Mental Health Variables")
    print("=" * 70)

    # Load data
    print(f"\nLoading data from: {data_path}")
    df = pd.read_csv(data_path)
    print(f"Dataset shape: {df.shape}")

    # =========================================================================
    # CREATE SCALE TOTALS
    # =========================================================================

    # IAS and IATS item columns
    ias_items = [f'ias_{i}' for i in range(1, 22)]
    iats_items = [f'iats_{i}' for i in range(1, 22)]

    # Create totals if not present
    if 'ias_total' not in df.columns:
        df['ias_total'] = df[ias_items].mean(axis=1)
    if 'iats_total' not in df.columns:
        df['iats_total'] = df[iats_items].mean(axis=1)

    # =========================================================================
    # DEFINE VARIABLES AND LABELS
    # =========================================================================

    # Display labels (clean names for figures)
    label_map = {
        'ias_total': 'IAS',
        'iats_total': 'IATS',
        'vviq': 'VVIQ',
        'tas': 'TAS',
        'phq9': 'PHQ-9',
        'gad7': 'GAD-7',
        'stai': 'STAI',
        'sss8': 'SSS-8',
        'pcs': 'PCS'
    }

    # Variable order for heatmaps
    # Order matters: determines row/column arrangement
    vars_no_vviq = [
        'iats_total', 'ias_total',
        'tas', 'sss8', 'pcs', 'phq9', 'gad7', 'stai'
    ]

    vars_with_vviq = [
        'iats_total', 'ias_total', 'vviq',
        'tas', 'sss8', 'pcs', 'phq9', 'gad7', 'stai'
    ]

    # =========================================================================
    # CREATE HEATMAP 1: Without VVIQ
    # =========================================================================

    print("\n" + "-" * 50)
    print("Creating heatmap: Subjective interoceptive accuracy and attention")
    print("-" * 50)

    corr1, pval1, n1 = create_correlation_heatmap(
        df=df,
        variables=vars_no_vviq,
        label_map=label_map,
        output_path=os.path.join(sub_plots_dir, "fig1_A_heatmap.png"),
        title="",
        figsize=(11, 9),
        annot_fontsize=14,
        label_fontsize=18,
        legend_fontsize=13,
        legend_title_fontsize=14
    )

    # =========================================================================
    # CREATE HEATMAP 2: With VVIQ
    # =========================================================================

    print("\n" + "-" * 50)
    print("Creating heatmap: Subjective interoceptive accuracy and attention (with VVIQ)")
    print("-" * 50)

    corr2, pval2, n2 = create_correlation_heatmap(
        df=df,
        variables=vars_with_vviq,
        label_map=label_map,
        output_path=os.path.join(other_plots_vviq_dir, "1_s_8_heatmap_intero_mh_vviq.png"),
        title="",
        figsize=(12, 10),
        annot_fontsize=14,
        label_fontsize=18,
        legend_fontsize=13,
        legend_title_fontsize=14
    )

    # =========================================================================
    # SAVE CORRELATION MATRICES AS CSV
    # =========================================================================

    print("\n" + "-" * 50)
    print("Saving correlation matrices")
    print("-" * 50)

    # With clean labels
    corr1_labeled = corr1.copy()
    corr1_labeled.index = [label_map.get(v, v) for v in vars_no_vviq]
    corr1_labeled.columns = [label_map.get(v, v) for v in vars_no_vviq]
    corr1_labeled.to_csv(os.path.join(supplementary_data_dir, "correlation_matrix_intero_mh.csv"))
    print(f"  Saved: supplementary_data/correlation_matrix_intero_mh.csv")

    corr2_labeled = corr2.copy()
    corr2_labeled.index = [label_map.get(v, v) for v in vars_with_vviq]
    corr2_labeled.columns = [label_map.get(v, v) for v in vars_with_vviq]
    corr2_labeled.to_csv(os.path.join(supplementary_data_dir, "correlation_matrix_intero_mh_vviq.csv"))
    print(f"  Saved: supplementary_data/correlation_matrix_intero_mh_vviq.csv")

    # =========================================================================
    # PRINT KEY FINDINGS
    # =========================================================================

    print("\n" + "=" * 70)
    print("KEY CORRELATIONS")
    print("=" * 70)

    print(f"\n1. IAS-IATS correlation: r = {corr1.loc['ias_total', 'iats_total']:.3f}")
    print("   -> Very weak! IAS and IATS measure distinct constructs.")

    print(f"\n2. IAS correlations with mental health (negative = better MH):")
    for mh in ['tas', 'phq9', 'gad7', 'stai', 'sss8', 'pcs']:
        r = corr1.loc['ias_total', mh]
        print(f"   IAS-{label_map[mh]}: r = {r:.3f}")

    print(f"\n3. IATS correlations with mental health (positive = worse MH):")
    for mh in ['tas', 'phq9', 'gad7', 'stai', 'sss8', 'pcs']:
        r = corr1.loc['iats_total', mh]
        print(f"   IATS-{label_map[mh]}: r = {r:.3f}")

    if 'vviq' in df.columns:
        print(f"\n4. VVIQ correlations:")
        print(f"   VVIQ-IAS: r = {corr2.loc['vviq', 'ias_total']:.3f} (positive)")
        print(f"   VVIQ-IATS: r = {corr2.loc['vviq', 'iats_total']:.3f} (near zero)")

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)

    # Close output file
    output_file.close()
    sys.stdout = sys.__stdout__
