# Interoceptive Profiling Identifies Separable Somatic and Alexithymic Contributors to Mental Health

Timo L. Kvamme & Juha Silvanto

Data and analysis code for the paper. Everything needed to reproduce every
number and figure is in this repository.

## What this is

In 832 adults we identified three interoceptive profiles from self-reported
interoceptive accuracy (IAS) and interoceptive attention (IATS) using k-means
clustering, then examined how those profiles relate to mental health. Structural
equation modelling showed a dual-pathway dissociation: alexithymia follows an
additive deficit model, while somatic symptoms follow a multiplicative
amplification model in which elevated attention weakens the protective effect of
accuracy.

### Three interoceptive profiles

- **Efficient** (n = 301, 36.2%): high accuracy, low attention. Lowest psychopathology across all domains.
- **Hypervigilant** (n = 273, 32.8%): moderate accuracy, high attention. Highest somatic symptoms and pain catastrophizing.
- **Uncertain** (n = 258, 31.0%): low accuracy, moderate attention. Highest alexithymia.

### Key findings

- Profile membership explained 6 to 15% of variance in mental health outcomes (all P < 0.001).
- The alexithymia pathway showed purely additive IAS and IATS effects, with no interaction (P = 0.80).
- The somatic pathway showed a significant IAS x IATS interaction (P = 0.003), where accuracy was most protective at high attention.
- A minimal 2-item subset exceeded the predictive performance of the two 42-item scale totals (within-sample cross-validation, so optimistic).

### New here? Read these three files in this order

1. **This README**, for what the project is and how to run it.
2. **[CODEBOOK.md](CODEBOOK.md)**, for what every column means, how each
   questionnaire is scored, and which items are reverse-keyed. Read it before
   you compute anything from the item columns. Two things in there catch people
   out: the response anchors are 1-indexed rather than 0-indexed, and two rows
   carry totals that are not the sum of their items.
3. **[CORRECTIONS.md](CORRECTIONS.md)**, for the three scoring defects and the
   duplicate participant that were
   found and fixed, and exactly what each one changed.

## Data

`dfc_interoception_profiling.csv` is the analysis dataset. One row per
participant, N = 832, 185 columns. It contains raw item responses, scale totals,
and demographics.

- Direct and indirect identifiers are **anonymized** in this release. No
  anonymized column is used by any analysis.
- Item columns are stored **after** reverse keying. Do not reverse them again.
- Response anchors are **1-indexed**, so several totals sit one point per item
  above their published scale. This changes no result but does change every
  reported mean. [CODEBOOK.md](CODEBOOK.md) has the per-scale table.
- Ten item responses are missing across two participants. Their affected totals
  are prorated or imputed and both rows are marked in the `scoring_flag` column.

Everything above is documented per column in [CODEBOOK.md](CODEBOOK.md).

## How to reproduce

All paths in this repository are **relative**. Clone it, install the
dependencies, and run from the repository root. Nothing needs editing.

**Requirements**

- R 4.4.0 or newer, with `lavaan`, `glmnet`, `randomForest`, `igraph`,
  `tidyLPA`, `ggplot2`, `gridExtra`, `dplyr`, `tidyr`, `openxlsx`.
- Python 3.10 or newer, with `pandas`, `numpy`, `scipy`, `matplotlib`,
  `seaborn`, `scikit-learn`.
- Script 15 additionally needs `scikit-learn-extra` and `gower`.
  `scikit-learn-extra` is compiled against NumPy 1.x and fails to import under
  NumPy 2.x, and it has not been updated. Run that one script in a Python
  3.9-3.11 environment with NumPy < 2. Every other script runs without it, and
  nothing downstream reads its output. Point the runner at that environment with
  `PYTHON_EXE`; [REQUIREMENTS.md](REQUIREMENTS.md) gives the exact versions and
  the two commands that build it.
- [REQUIREMENTS.md](REQUIREMENTS.md) lists every package version the deposited
  outputs were produced with, and records what does and does not reproduce
  exactly between two runs.

**Run the whole pipeline**

```bash
python scripts/run_analysis.py
```

The runner locates the repository root itself, so it works whichever directory
you call it from.

**Run selected steps, or list them**

```bash
python scripts/run_analysis.py --list
python scripts/run_analysis.py 04 12
```

**Run one script directly**, from the repository root:

```bash
Rscript scripts/04_cluster_analysis.R
Rscript scripts/12_sem_pathway_analysis.R
```

Scripts are numbered in execution order. Step `00b` rebuilds the two derived
files under `supplementary_data/` from the master data file, so later steps
never read a stale snapshot. Run it before any individual script that reads
`supplementary_data/`.

`00_score_corrections.py` is **not** part of the pipeline. It was applied once to
produce the deposited dataset and is kept for provenance. It is idempotent, so
running it against the deposited file reports "No changes written."

**Verify the scoring logic**

Both data-preparation scripts carry a self-check that fails loudly if their
logic breaks:

```bash
python scripts/00_score_corrections.py --self-check
python scripts/00b_build_derived_data.py --self-check
```

## Repository structure

```
├── dfc_interoception_profiling.csv    # Dataset (N = 832). Start at CODEBOOK.md.
├── CODEBOOK.md                        # Every column, every scale, every reverse-keyed item
├── CORRECTIONS.md                     # The four data corrections and what each changed
├── README.md                          # This file
├── scripts/
│   ├── 00_score_corrections.py        # One-off scoring corrections (provenance, idempotent)
│   ├── 00b_build_derived_data.py      # Rebuilds derived supplementary_data files (pipeline step 1)
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
│   ├── 13_tas_subscale_analysis.R
│   ├── 13b_tas_subscale_path_diagram.R
│   ├── 14_vviq_interoception_mediation.R
│   ├── 15_item_level_clustering_sensitivity.py
│   ├── 16_alexithymia_somatic_pathway.py
│   ├── 17_revision_sensitivity.R      # Age/gender covariates, SSS-8-only mediator (Supp. Table 7, Supp. Fig. 24)
│   ├── 18_somatic_factor_structure.R  # Factor structure of the somatic composite (Supp. Table 6)
│   ├── 19_ert_exact_statistics.R      # Exact test statistics quoted in the text
│   ├── supplementary_figure_scripts/  # Scripts for Supp. Figs. 8, 11-14, 18 and 24
│   └── run_analysis.py                # Runs all scripts in sequence
├── analysis_output/                   # Statistical output logs, one per script
├── plots/                             # Main manuscript figures (Figs. 1-3)
│   ├── sub_plots/                     # Individual figure panels
│   ├── other_plots_vviq_plots/        # VVIQ analyses, not printed in the paper
│   └── supplementary_plots/           # Supplementary figures, see FIGURE_MAP.md
└── supplementary_data/                # Intermediate data tables and results
```

The file names of the supplementary figures use a per-analysis prefix, not the
numbers printed in the paper. `plots/supplementary_plots/FIGURE_MAP.md` maps
Supplementary Figures S1 to S24 to their file names and to the script that
writes each one.

`plots/figure_18_scree_parallel_somatic.png` is not a manuscript figure. The
paper prints Figures 1 to 3 only. The 18 is the number of the script that writes
it, `18_somatic_factor_structure.R`. It is the scree and parallel analysis plot
behind Supplementary Table 6.

## Scoring corrections

An item-level audit of every scale during the Communications Psychology revision
found three scoring defects and one duplicate participant in the originally
deposited file. All four are
fixed here, and every output and figure in this repository was regenerated from
the corrected data.

1. **STAI-T reverse keying.** Two of the nine reverse-keyed trait items (items 3
   and 14) had never been reverse-scored.
2. **Missing items scored as zero.** Totals were computed with
   `rowSums(na.rm = TRUE)`, so an omitted item counted as zero. One participant
   who skipped the entire SSS-8 received a total of 0, below the scale floor of 8.
3. **GAD-7 response-coding gap.** The four response options were exported from
   Qualtrics as 1, 3, 4, 5 rather than 1, 2, 3, 4, so the stored total was a
   non-linear transform of the true score.

The full detail, including which rows were affected and exactly what each
correction changed in the results, is in **[CORRECTIONS.md](CORRECTIONS.md)**.
The headline conclusions and the three-profile solution are unchanged, but one
sensitivity analysis moved across the .05 threshold. That is documented rather
than glossed over.

## Verification and known limitations

The scoring was verified against the original Qualtrics export, which holds both
the numeric response codes and the choice labels the participants saw. The check
joined that export to the pre-anonymisation working copy on the response
identifier (the released file carries anonymised identifiers; the same match can be
reproduced from start time, end time and duration), matched all 833 rows of the
file as first deposited (832 after the duplicate completion was removed),
and compared every item. [CODEBOOK.md](CODEBOOK.md) section 5 has
the detail. In summary:

- **The column-to-item mapping is verified.** `stai_3` really is STAI-T item 3,
  and so on for all 145 items. This is not an inference from the column naming.
- **The reverse keying is verified at source.** TAS-20 items 4, 5, 10, 18, 19 and
  STAI-T items 1, 3, 6, 7, 10, 13, 14, 16, 19 are reversed, exactly matching the
  published keys. No other scale has a reverse-keyed item.
- **The VVIQ direction is settled and the analysis is correct.** In this file a
  high value means vivid imagery. The raw survey used the opposite Marks (1973)
  order, and preprocessing reversed all 16 items.
- **The `sass_*` columns are the Somatosensory Amplification Scale** (SSAS;
  Barsky, Wyshak & Klerman, 1990). It has no reverse-keyed items. It was
  collected in recruitment wave 1 only and no analysis script reads it.

One limitation remains:

- **One SSS-8 total is imputed, not measured** (row 511, flagged in
  `scoring_flag`). The published results do not depend on whether that row is
  kept or dropped.

## License

See [LICENSE](LICENSE). Please cite the paper if you use this code or data.
