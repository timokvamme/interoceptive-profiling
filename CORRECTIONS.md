# Scoring corrections

An item-level audit of every scale during the Communications Psychology revision
(COMMSPSYCHOL-26-0385-T) found three scoring defects and one duplicate
participant in the originally deposited
data file. All three are fixed here. Every output and figure in this repository
was regenerated from the corrected data.

This file records what was wrong, how it was found, what was done, and exactly
what changed. Nothing is glossed over: one sensitivity analysis moved across the
.05 threshold, and that is documented in full below.

`scripts/00_score_corrections.py` applies all four corrections. It is
idempotent, so running it against the deposited file reports "No changes
written". It carries a self-check:

```bash
python scripts/00_score_corrections.py --self-check
```

---

## Correction 1: STAI-T reverse keying

**The defect.** The trait scale of the State-Trait Anxiety Inventory (Form Y-2)
has nine reverse-keyed items: 1, 3, 6, 7, 10, 13, 14, 16, 19. Preprocessing
reversed seven of them and missed items 3 ("I feel satisfied with myself") and
14 ("I make decisions easily").

**How it was found.** Both items had strongly negative corrected item-total
correlations, which is the signature of an item that is still scored in the
wrong direction.

| Item | Item-total r before | Item-total r after |
|---|---|---|
| `stai_3` | **-0.734** | **+0.740** |
| `stai_14` | **-0.577** | **+0.583** |

For comparison, the other 18 STAI-T items range from +0.54 to +0.78.

**The fix.** `stai_3` and `stai_14` are reversed (`5 - x`) and the STAI total is
recomputed as the sum of the 20 items.

**What changed.** The STAI total mean moved from 46.21 to 46.77. The corrected
total correlates **0.9864** with the previous one. Every STAI-derived number in
the manuscript moved at the second or third decimal.

**Verified at source.** The reverse set was later confirmed against the original
Qualtrics export, participant by participant: exactly items 1, 3, 6, 7, 10, 13,
14, 16, 19 are reversed in the corrected file, matching the published key. See
[CODEBOOK.md](CODEBOOK.md) section 5.

---

## Correction 2: missing item responses scored as zero

**The defect.** Scale totals were computed with `rowSums(na.rm = TRUE)`, so an
omitted item contributed 0 to the total rather than being treated as missing.
Ten item responses are missing, across two participants.

**How it was found.** One SSS-8 total was 0, which is below the scale floor of 8
and therefore impossible.

**The fix.** Totals with isolated missing items are prorated from the answered
items (person mean times the number of items, rounded). The participant who
omitted the entire SSS-8 is retained, because every other questionnaire is
complete, with the SSS-8 total set to the rounded sample mean of the other 831
participants. Every affected row is marked in the new `scoring_flag` column.

| Row | Column | Before | After | Method |
|---|---|---|---|---|
| 512 | `sss8` | **0** | 17 | mean-imputed; all 8 items missing |
| 512 | `pcs` | 25 | 27 | prorated from 12 of 13 items |
| 515 | `sss8` | 15 | 17 | prorated from 7 of 8 items |

**Row 512's SSS-8 total is imputed, not measured.** That is the one place in this
dataset where a value was invented rather than observed. It is flagged so anyone
recomputing totals from items can identify it. The published results do not
depend on whether the row is kept or dropped, as shown below.

**What changed, and the one significance crossing.** This correction, not the
STAI one, moved a result across the .05 threshold. In the SSS-8-only sensitivity
model, the IAS x IATS interaction on the somatic factor was:

| | Before | After |
|---|---|---|
| Standardized estimate | -0.063 | -0.068 |
| p, standardized solution | .0497 | .0331 |
| p, delta method | .0503 | .0336 |
| p, bootstrap | .0500 | .0330 |
| Delta-method 95% CI (standardized) | [-0.126, -0.000] | [-0.131, -0.005] |

The effect was sitting exactly on the threshold before, and its significance
flipped depending on which estimator was used. It is now below .05 by all three.

**The cause was isolated, not assumed.** Refitting `z(SSS-8) ~ IAS * IATS` by
ordinary least squares on both versions of the file reproduces the move exactly
(p = .0512 before, p = .0343 after). Substituting only the corrected `sss8`
column into the old file reproduces the new result exactly. Substituting only the
old `stai` column into the new file leaves the new result unchanged. So the STAI
correction contributes nothing to this path, and the driver is the row with the
all-missing SSS-8 (row 512 before the duplicate was removed, row 511 in the
released file), whose
SSS-8 total of 0 was a floor artifact created by `rowSums(na.rm = TRUE)`.

**Robustness.** Dropping both affected participants from the analysis entirely
gives p = **.0378**. So the corrected conclusion does not depend on the
missing-data method, and the pre-correction p = .051 was the artifact.

The primary model, which uses the combined SSS-8 + PCS somatic mediator, is
unaffected: the interaction is p = .0041 before and p = .0030 after, significant
in both.

---

## Correction 3: GAD-7 response-coding gap

**The defect.** The seven GAD-7 items were exported from Qualtrics with the
recode values 1, 3, 4, 5. The value 2 never occurred in any item, in any of the
833 rows of the file as first deposited. GAD-7 has four response options, so the
second option ("Several days")
carried the code 3 instead of 2. The stored total was therefore a non-linear
transform of the true score.

**How it was found.** A value-frequency check across every scale. A zero-count
middle category, between counts of 2,495 and 1,835, is not a response pattern.
It was later confirmed in the raw Qualtrics export, where the same gap is
present at source.

**The fix.** The four options are remapped 1, 3, 4, 5 to 1, 2, 3, 4 and the
GAD-7 total is recomputed as the sum of the seven items.

**What changed.**

| | Before | After |
|---|---|---|
| Mean | 17.51 | **13.51** |
| SD | 7.94 | 5.65 |
| Range | 7-35 | 7-28 |

The corrected total correlates **0.9880** with the previous one, so the ordering
of participants barely moves. The mental-health composite, which is built from
z-scores, correlates **0.9971** with its pre-correction version across all three
corrections combined (STAI-T alone 0.9864, GAD-7 alone 0.9880).

**Effect on the results.** The whole pipeline was rerun and every output file was
compared line by line against the pre-correction version.

| Measure | Result |
|---|---|
| Output lines changed | 559 of 4,001 (14.0%) |
| **p-values crossing .05, in either direction** | **0** |
| **Confidence intervals changing whether they contain zero** | **0** |
| k-means solution | unchanged (k = 3) |
| LPA solution | unchanged (3 profiles by BIC) |
| Adjusted Rand index | unchanged (0.266) |
| Primary SEM interaction, IAS x IATS on somatic | p = .0030, unchanged |
| SSS-8-only sensitivity interaction | p = .0336, unchanged |

So the GAD-7 correction moves many numbers at the second or third decimal and
corrects the reported GAD-7 descriptives, but changes no inferential conclusion
anywhere in the analysis. Correlations with other scales move by at most 0.018.

The figures in that table were measured on the 833-row file, before correction 4
removed the duplicate completion. They record what the GAD-7 fix alone changed.
On the final 832-participant sample the adjusted Rand index is **0.269**; the
k-means and LPA solutions and both interaction p-values are unchanged.

**A related gap that was already handled.** The raw export codes the IAS response
"Strongly agree" as 13 rather than 5, leaving the same kind of gap. Preprocessing
remapped it correctly, and this was verified item by item. No action was needed.

---

## Related repository fixes

Two problems found in the same audit were not scoring defects but did affect
reproducibility. Both are fixed.

**Stale derived data files.** `supplementary_data/data_with_ias_iats_ratio.csv`
and `supplementary_data/interoception_item_data.csv` were read by the pipeline
but written by no step in it. They were hand-made snapshots, so they went stale
whenever the master data file changed. The first still carried the broken GAD-7
coding, so `18_somatic_factor_structure.R` was running on uncorrected data. The
second holds each item's correlation with the mental-health composite, which
`02_correlation_ias_iats.R` plots, so a figure was drawn from stale values.
Both are now rebuilt from the master by `scripts/00b_build_derived_data.py`,
which runs as the first step of the pipeline.

**Hard-coded absolute paths.** Every script read the data from an absolute path
on the original author's machine, so no one else could run the pipeline. All
paths are now relative to the repository root.

**A false statement in the output.** `18_somatic_factor_structure.R` printed that
"the item sums reproduce the sss8 / pcs scale scores". After correction 2 that is
no longer true for the two flagged rows. The text now says so explicitly.

**Scripts that chose the wrong working directory.** Scripts 06, 08, 15 and 16
resolved their paths from the folder holding the script, which is correct in the
working repository but wrong in this release, where the scripts live in
`scripts/` and the data sits at the repository root. They now locate the
repository root by looking for the data file, so both layouts work.

---

## Script 15 needs a separate environment

`15_item_level_clustering_sensitivity.py` cannot run under NumPy 2. Its outputs in this repository were regenerated on 2026-09-03 from the fully corrected data in a NumPy 1 environment, so they are current. It needs
`scikit-learn-extra`, whose `KMedoids` is compiled against NumPy 1.x and fails to
import under NumPy 2.x with:

```
ImportError: numpy.core.multiarray failed to import
```

`scikit-learn-extra` has not been updated for NumPy 2, so there is no version of
the package that fixes this. To regenerate this one analysis, run it in an
environment with NumPy < 2:

```bash
pip install scikit-learn-extra gower      # in a Python 3.9-3.11 env with numpy<2
python scripts/15_item_level_clustering_sensitivity.py
```

**The files in `analysis_output/15_output.txt` and `supplementary_data/item_level_*.csv`** were produced after all four corrections.
---

## Correction 4: one participant completed the survey twice

**The defect.** One person completed the whole survey twice. Both responses
passed every data-quality screen and both were in the analysed sample, so that
person contributed two of the 833 rows and the effective sample was 832 people.

**How it was found.** An identity audit of the whole sample. Exactly one
self-entered Prolific ID occurs twice, and the same two rows are the only pair
that shares an IP address. No other pair of rows shares an ID, an ID within one
character, or an IP.

| | First completion | Second completion |
|---|---|---|
| Prolific ID | identical | identical |
| IP address, geolocation | identical | identical |
| Age, gender | 36, same | 36, same |
| Start (Excel serial) | 45574.74073 | 45574.75612 |
| Duration | 1302 s | 679 s |

The second session began **27 seconds after the first one ended**. It is a
repeat completion rather than a duplicated export row: only 51% of the 119 item
responses agree between the two rows (r = .56). Response similarity alone would
not have found it; the identity evidence is what settles it.

**The fix.** The second completion is removed and the first is kept, because the
first was given without prior exposure to the items. The analysed sample is
**832 participants, one row per person**. The upstream preprocessing script now
drops later repeat completions by the same Prolific ID, so a rerun cannot
reintroduce the row.

`scripts/00_score_corrections.py` identifies the row two ways, which select the
same row: by Prolific ID where the identifiable source file is used, and by the
exact start time, end time and duration in this released file, whose identifiers
are anonymised.

**What changed.** Every analysis was rerun on the 832-participant sample. No
conclusion changed. Degrees of freedom fall by one throughout, and coefficients
move in the third decimal.


**Which of their numbers are affected, and which are not.** The script clusters
on the 42 IAS and IATS items only (`all_items = ias_cols + iats_cols`), and no
correction touched those items. So:

- The cluster solutions, k-selection tables, silhouettes, cluster sizes,
  cross-tabulations against the published clusters, and the bootstrap stability
  results are **unaffected and remain valid**.
- Only the outcome ANOVA tables are affected, and within them only the **GAD-7**
  and **STAI** rows, because those two totals changed. Given that the full rerun
  of every other script produced zero p-value crossings, the direction and
  significance of those rows are very unlikely to change, but they have not been
  recomputed and should not be quoted from these files without rerunning.

This analysis is a supplementary sensitivity check. No result in the main
manuscript depends on it.
