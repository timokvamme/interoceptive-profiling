# Interoceptive Profiling Identifies Separable Somatic and Alexithymic Contributors to Mental Health

Timo L. Kvamme & Juha Silvanto

## Overview

In 833 adults, we identified three interoceptive profiles based on self-reported interoceptive accuracy (IAS) and attention (IATS) using k-means clustering, and examined their associations with mental health outcomes. Structural equation modeling revealed a dual-pathway dissociation: alexithymia follows an additive deficit model, whereas somatic symptoms follow a multiplicative amplification model in which elevated attention weakens accuracy's protective effect.

### Three interoceptive profiles

- **Efficient** (n = 302, 36%): high accuracy, low attention — lowest psychopathology across all domains
- **Hypervigilant** (n = 273, 33%): moderate accuracy, high attention — highest somatic symptoms and pain catastrophizing
- **Uncertain** (n = 258, 31%): low accuracy, moderate attention — highest alexithymia

### Key findings

- Profile membership explained 6–15% of variance in mental health outcomes (all P < 0.001)
- The alexithymia pathway showed purely additive IAS and IATS effects (no interaction, P = 0.80)
- The somatic pathway showed a significant IAS x IATS interaction (P = 0.004), where accuracy was protective only at high attention levels
- A minimal 2-item subset matched the predictive performance of the full 42-item battery

## Repository structure

```
├── dfc_interoception_profiling.csv    # Dataset (N = 833)
├── scripts/                           # Analysis pipeline (R and Python)
│   ├── 01_normality_checks.R
│   ├── 02_correlation_ias_iats.R
│   ├── 03_correlation_heatmaps.py
│   ├── 04_cluster_analysis.R
│   ├── 05_lpa_analysis.R
│   ├── 06_simpsons_paradox_bayes.py
│   ├── 07_at_risk_subcluster.py
│   ├── 08_cluster_visualizations.py
│   ├── 09_item_selection_analysis.R
│   ├── 10_network_analysis.R
│   ├── 11_figure_combined.R
│   ├── 12_sem_pathway_analysis.R
│   ├── 12b_pathway_dominance_test.R
│   ├── 12c_multigroup_omnibus_test.R
│   ├── 12d_pathway_dominance_figure.R
│   ├── 12e_comprehensive_moderation_figure.R
│   └── run_analysis.py                # Runs all scripts in sequence
├── analysis_output/                   # Statistical output logs
├── plots/
│   ├── main/                          # Main manuscript figures (Figs. 1–3)
│   ├── sub_plots/                     # Individual figure panels
│   └── supplementary_plots/           # Supplementary figures
└── supplementary_data/                # Intermediate data tables and results
```

## Running the analyses

Scripts are numbered in execution order. All scripts assume the dataset (`dfc_interoception_profiling.csv`) is located in the project root (`C:/code/projects/interoceptive-profiling/`).

**Requirements:** R (4.4.0+) with `lavaan`, `glmnet`, `randomForest`, `igraph`, `tidyLPA`, `ggplot2`, `gridExtra`; Python (3.10+) with `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `scikit-learn`.

To run the full pipeline:

```bash
python scripts/run_analysis.py
```

Or run individual scripts:

```bash
Rscript scripts/04_cluster_analysis.R
Rscript scripts/12_sem_pathway_analysis.R
```

## License

Please cite the paper if you use this code or data.
