# Codebook

Every variable in `dfc_interoception_profiling.csv`, how each questionnaire is
scored, and which items are reverse-keyed.

Read this before you compute anything from the item columns. Two things in this
file surprise people: the response anchors are 1-indexed rather than 0-indexed,
and two rows carry totals that are not the sum of their items.

- N = 832 participants, 185 columns. One row per person: the one participant who
  completed the survey twice is represented by their first completion only
  (see [CORRECTIONS.md](CORRECTIONS.md), correction 4).
- One row per participant. Row order is the analysis order used throughout.
- Corrections applied to this file are documented in [CORRECTIONS.md](CORRECTIONS.md).

---

## 1. Response coding: read this first

**Every questionnaire in this file is stored 1-indexed, exactly as Qualtrics
exported it.** Several of these instruments are published with 0-based anchors.
The stored coding is shifted up by one point per item.

| Scale | Stored anchors | Published anchors | Stored total is inflated by |
|---|---|---|---|
| SSS-8 | 1 to 5 | 0 to 4 | 8 points |
| PCS | 1 to 5 | 0 to 4 | 13 points |
| PHQ-9 | 1 to 4 | 0 to 3 | 9 points |
| GAD-7 | 1 to 4 | 0 to 3 | 7 points |
| TAS-20 | 1 to 5 | 1 to 5 | not shifted |
| STAI-T | 1 to 4 | 1 to 4 | not shifted |
| IAS, IATS, VVIQ, SASS | 1 to 5 | 1 to 5 | not shifted |

What this does and does not affect:

- **It does not affect any result in the paper.** A constant shift per item adds
  a constant to the total. It leaves every correlation, z-score, cluster
  solution, factor loading and path coefficient identical. The mental-health
  composite is built from z-scores, so it is shift-invariant by construction.
- **It does affect every reported mean**, by the item count.
- **It makes published clinical cut-offs inapplicable to the stored totals.**
  To apply, for example, the PHQ-9 cut-off of 10, first subtract 9 from the
  stored `phq9`. Applying the cut-off to the stored total instead would flag
  760 of 832 participants rather than 294. No script in this repository applies
  a clinical cut-off, so nothing in the published analyses depends on this.

---

## 2. Scales

Observed values below are measured from this file, not from the manuals.
Alpha is Cronbach's alpha on complete cases.

| Scale | Item columns | k | Total column | Observed total | Mean (SD) | Alpha | Reverse-keyed items |
|---|---|---|---|---|---|---|---|
| TAS-20 | `tas_1` … `tas_20` | 20 | `tas` | 20-87 | 50.80 (12.74) | .880 | **4, 5, 10, 18, 19** |
| IAS | `ias_1` … `ias_21` | 21 | `ias` | 36-105 | 77.72 (10.79) | .865 | none |
| IATS | `iats_1` … `iats_21` | 21 | `iats` | 21-91 | 47.11 (13.97) | .911 | none |
| SSAS (`sass_*`) | `sass_1` … `sass_10` | 10 | `sass` | 12-49 | 29.70 (6.23) | .712 | none |
| SSS-8 | `sss_8_1` … `sss_8_8` | 8 | `sss8` | 8-39 | 16.97 (6.02) | .797 | none |
| PCS | `pcs_1` … `pcs_13` | 13 | `pcs` | 13-65 | 30.38 (12.51) | .955 | none |
| PHQ-9 | `phq_9_1` … `phq_9_9` | 9 | `phq9` | 9-36 | 17.15 (6.38) | .900 | none |
| GAD-7 | `gad_7_1` … `gad_7_7` | 7 | `gad7` | 7-28 | 13.51 (5.65) | .923 | none |
| STAI-T | `stai_1` … `stai_20` | 20 | `stai` | 20-80 | 46.77 (13.13) | .951 | **1, 3, 6, 7, 10, 13, 14, 16, 19** |
| VVIQ | `vviq_q_2_1_1` … `vviq_q2_4_4` | 16 | `vviq` | 16-80 | 48.91 (16.69) | .968 | none |

**Every total column is the plain sum of its item columns**, with two documented
exceptions listed in section 4.

### Reverse keying: what is already applied

**The item columns in this file are stored POST-reversal.** Do not reverse them
again. A reverse-keyed item's stored value is already `(max + 1) - raw`:

- TAS-20 items 4, 5, 10, 18, 19: stored as `6 - raw`.
- STAI-T items 1, 3, 6, 7, 10, 13, 14, 16, 19: stored as `5 - raw`.
- **All 16 VVIQ items**: stored as `6 - raw`, so a high stored value means vivid
  imagery. The raw survey ran the other way. See section 5.

You can confirm this yourself: every item in every scale has a positive
corrected item-total correlation. If a reverse-keyed item were still stored raw,
its correlation would be strongly negative. That is exactly how the STAI defect
described in [CORRECTIONS.md](CORRECTIONS.md) was found. Section 5 goes further
and verifies the keying against the original survey export, item by item.

### The one weak item

`tas_5` has a corrected item-total correlation of **-0.011** against the full
TAS-20. This is **not** a keying error, and it should be left alone:

- Against its own subscale (EOT) the correlation is **+0.123**.
- If the item were flipped, that within-subscale correlation would become
  **-0.123**, which is worse.
- TAS-20 item 5 is the weakest item of the instrument in most published samples.

The whole EOT subscale is weak in this sample (within-subscale item-total
correlations .12 to .46), which is consistent with the literature.

### TAS-20 subscales

| Subscale | Items |
|---|---|
| DIF (difficulty identifying feelings) | 1, 3, 6, 7, 9, 13, 14 |
| DDF (difficulty describing feelings) | 2, 4, 11, 12, 17 |
| EOT (externally oriented thinking) | 5, 8, 10, 15, 16, 18, 19, 20 |

### Item wordings

The 42 IAS and IATS item wordings are in
`supplementary_data/interoception_item_data.csv`, in the `questionnaire_text`
column, keyed by the same column names used here.

The other instruments are **not** reproduced here. TAS-20, STAI-T, PCS, PHQ-9,
GAD-7, SSS-8 and VVIQ are published, separately licensed instruments; several
are copyrighted and cannot be redistributed. Use the citations below to obtain
the official item text and confirm the numbering against the column names.

| Scale | Source |
|---|---|
| TAS-20 | Bagby, Parker & Taylor (1994), *J Psychosom Res* 38(1), 23-32 |
| IAS | Murphy et al. (2019), *Behav Res Methods* 52, 2287-2299 |
| IATS | Gabriele et al. (2022), *Biol Psychol* 168, 108243 |
| SSAS (`sass_*`) | Barsky, Wyshak & Klerman (1990), *J Psychiatr Res* 24(4), 323-334 |
| SSS-8 | Gierk et al. (2014), *JAMA Intern Med* 174(3), 399-407 |
| PCS | Sullivan, Bishop & Pivik (1995), *Psychol Assess* 7(4), 524-532 |
| PHQ-9 | Kroenke, Spitzer & Williams (2001), *J Gen Intern Med* 16(9), 606-613 |
| GAD-7 | Spitzer et al. (2006), *Arch Intern Med* 166(10), 1092-1097 |
| STAI-T | Spielberger et al. (1983), Form Y-2, Mind Garden (copyrighted) |
| VVIQ | Marks (1973), *Br J Psychol* 64(1), 17-24 |

STAI-T Form Y-2 items 1 to 20 correspond to items 21 to 40 of the full Form Y.
Its nine reverse-keyed items are 21, 23, 26, 27, 30, 33, 34, 36, 39 in full-form
numbering.

---

## 3. Other columns

| Column | Meaning |
|---|---|
| `age`, `gender` | Demographics. |
| `evenodd_score` | Attention-check / even-odd consistency score. |
| `time_used`, `duration (in seconds)` | Completion time. |
| `q1_2` … `q1_8`, `q93` | Recruitment and platform identifiers. **Anonymized in this release**; not used by any analysis. |
| `ipaddress`, `responseid`, `recipient*`, `locationlatitude`, `locationlongitude` | **Anonymized in this release**; not used by any analysis. |
| `startdate`, `enddate`, `recordeddate` | Qualtrics timestamps, serial format. |
| `status`, `progress`, `finished`, `distributionchannel`, `userlanguage`, `q_recaptchascore` | Qualtrics survey metadata. |
| `scoring_flag` | **Added by the corrections.** Empty for all but two rows. See section 4. |

### Derived columns in `supplementary_data/data_with_ias_iats_ratio.csv`

That file is this data file plus six derived columns. It is rebuilt from the
master by `00b_build_derived_data.py`, which runs as the first pipeline step.
Do not edit it by hand; it will be overwritten.

| Column | Definition |
|---|---|
| `cohort` | `control` if the SASS was administered (recruitment wave 1, n = 511), `aphan` otherwise (wave 2, aphantasia-enriched, n = 321). |
| `ias_total` | `ias / 21`, the per-item mean. |
| `iats_total` | `iats / 21`, the per-item mean. |
| `ias_iats_ratio` | `ias_total / iats_total`. |
| `vviq_k`, `vviq_q` | Carried over from the original file. They depend only on the VVIQ items, which no correction touched. |

---

## 4. Missing data and the `scoring_flag` column

Ten item responses are missing, across two participants. Their scale totals are
therefore **not** the plain sum of their items. Both rows are marked in
`scoring_flag` so you can identify them when recomputing:

Row numbers below are 1-based data rows, counting the first participant as row 1.
In pandas, `pd.read_csv(...)` gives them the 0-based index one lower, so they are
`df.loc[511]` and `df.loc[514]`. Do not rely on the number alone: select on
`scoring_flag` instead.

| Row (1-based) | pandas index | `scoring_flag` | What was done |
|---|---|---|---|
| 512 | 511 | `pcs_prorated_from_12_of_13_items;sss8_mean_imputed_all_8_items_missing` | The PCS total is prorated from the 12 answered items. All eight SSS-8 items are missing, so the SSS-8 total is the rounded sample mean of the other 831 participants. |
| 515 | 514 | `sss8_prorated_from_7_of_8_items` | The SSS-8 total is prorated from the 7 answered items. |

Row 512's SSS-8 total is **imputed, not measured**. If your analysis depends on
the SSS-8, decide deliberately whether to keep or drop that row. The published
results do not depend on the choice; see [CORRECTIONS.md](CORRECTIONS.md).

Every other row is complete, and for every other row each total equals the sum
of its items exactly.

---

## 5. Verification against the original survey

Everything in this codebook was checked against the original Qualtrics export,
not inferred from the column names. The export exists in two paired versions,
one holding the numeric response codes and one holding the choice labels, so
each stored value can be traced back to the wording the participant saw. The
export contains identifiers and is not part of this public release.

The check joined the raw export to this data file on the Qualtrics response
identifier, matched all 833 rows of the file as first deposited (832 after the
duplicate completion was removed), and compared every item value.

**The column-to-item mapping is verified.** For every scale, each column matched
its raw counterpart participant by participant. `stai_3` is STAI-T item 3,
`tas_5` is TAS-20 item 5, and so on for all 145 items.

**The reverse keying is verified at source**, not inferred from correlations:

| Scale | Items reversed during preprocessing | Matches the published key |
|---|---|---|
| TAS-20 | 4, 5, 10, 18, 19 | yes |
| STAI-T | 1, 3, 6, 7, 10, 13, 14, 16, 19 | yes |
| VVIQ | all 16 | yes, see below |
| all others | none | yes |

**The VVIQ direction is settled: high means vivid, and the analysis is correct.**
The raw survey used the original Marks (1973) order, where 1 is "Perfectly clear
and as vivid as normal vision" and 5 is "No image at all". Preprocessing reversed
all 16 items, so in **this** file the order is inverted:

| Stored value | What the participant chose |
|---|---|
| 5 | Perfectly clear and as vivid as normal vision |
| 4 | Clear and reasonably vivid |
| 3 | Moderately clear and vivid |
| 2 | Vague and dim |
| 1 | No image at all, you only "know" that you are thinking of an object |

So `14_vviq_interoception_mediation.R` is right to treat a high total as vivid
imagery and a per-item mean at or below 2.0 as aphantasia. The resulting 18.4%
aphantasia rate is far above the 1 to 4% population estimate because recruitment
wave 2 was deliberately aphantasia-enriched (`cohort == "aphan"`), not because
the scale is inverted.

**The IAS recode gap was already handled.** In the raw export the IAS response
"Strongly agree" carried the code 13, not 5, leaving a gap like the GAD-7 one.
Preprocessing remapped 13 to 5 correctly, and this file contains contiguous
values 1 to 5. This was checked item by item.

**The SASS is the Somatosensory Amplification Scale (SSAS).** Barsky, Wyshak &
Klerman (1990), *J Psychiatr Res* 24(4), 323-334. Ten items, anchored 1 ("not at
all true") to 5 ("extremely true"). **The SSAS has no reverse-keyed items**, which
matches the verification above. It was collected in recruitment wave 1 only
(n = 511) and dropped from wave 2 to reduce participant burden. No analysis
script reads the `sass` column, so no published result depends on it. The items:

| Column | Item |
|---|---|
| `sass_1` | When someone else coughs, it makes me cough too. |
| `sass_2` | I can't stand smog, smoke or pollutants in the air. |
| `sass_3` | I am often aware of various things happening within my body. |
| `sass_4` | When I bruise myself, it stays noticeable for a long time. |
| `sass_5` | Sudden loud noises really bother me. |
| `sass_6` | I can sometimes hear my pulse or my heartbeat throbbing in my ear. |
| `sass_7` | I hate to be too hot or too cold. |
| `sass_8` | I am quick to sense hunger contractions in my body. |
| `sass_9` | Even something minor, like an insect bite or a splinter, really bothers me. |
| `sass_10` | I have low tolerance for pain. |

## 6. Remaining limitation

**One SSS-8 total is imputed, not measured.** Row 512 omitted all eight SSS-8
items, so its total is the rounded sample mean of the other 831 participants. It
is marked in `scoring_flag`. If your analysis depends on the SSS-8, decide
deliberately whether to keep or drop that row. The published results do not
depend on the choice; see [CORRECTIONS.md](CORRECTIONS.md).
