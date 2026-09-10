"""
================================================================================
00 - SCORE CORRECTIONS (applied once to the deposited data file)
================================================================================

Three scoring defects and one duplicate completion were found during the
Communications Psychology revision (COMMSPSYCHOL-26-0385-T), in an item-level
audit of every scale and an identity audit of every response.

1. STAI-T reverse keying.
   The trait scale of the State-Trait Anxiety Inventory (Form Y-2, 20 items)
   has nine reverse-keyed items: 1, 3, 6, 7, 10, 13, 14, 16, 19.
   Preprocessing reversed seven of them and missed items 3 ("I feel satisfied
   with myself") and 14 ("I make decisions easily"). Both items showed strongly
   negative corrected item-total correlations (-0.73 and -0.58) in the
   deposited file, which is the signature of an unreversed item.
   Fix: stai_3 and stai_14 are reversed (5 - x) and the STAI total is recomputed
   as the sum of the 20 items.

2. Missing item responses scored as zero.
   Scale totals were computed with rowSums(na.rm = TRUE), so an omitted item
   contributed 0 to the total. Two participants had missing items:
     - one omitted a single SSS-8 item,
     - one omitted all eight SSS-8 items (total recorded as 0, below the
       scale minimum of 8) and a single PCS item.
   Fix: totals with isolated missing items are prorated from the answered
   items (person mean x number of items, rounded to the nearest integer).
   The participant who omitted the whole SSS-8 is retained, because every
   other questionnaire is complete and the profile solution does not depend
   on this scale; their SSS-8 total is set to the rounded sample mean of the
   other 832 participants, which is a neutral value in the standardized
   somatic distress composite.
   Every affected row is marked in the `scoring_flag` column, so a reader who
   recomputes totals from the items can tell an imputed cell from a real one.

3. GAD-7 response-coding gap.
   The seven GAD-7 items were exported from Qualtrics with the recode values
   1, 3, 4, 5. The value 2 never occurs in any item, in any of the 833 rows.
   GAD-7 has four response options, so the second option ("Several days")
   carried the code 3 instead of 2. The stored total was therefore a
   non-linear transform of the true score and its mean was inflated by about
   four points.
   Fix: the four options are remapped 1, 3, 4, 5 -> 1, 2, 3, 4 and the GAD-7
   total is recomputed as the sum of the seven items.

4. Duplicate completion by one participant.
   One person completed the whole survey twice. Both responses passed every
   data-quality screen and both were in the analysed sample, so that person
   contributed two rows and the effective sample was 832 people in 833 rows.
   The two responses carry the same self-entered Prolific ID, the same IP
   address, the same geolocation, the same age (36) and the same gender, and
   the second session began 27 seconds after the first one ended. It is a
   repeat completion, not a duplicated export row: only 51% of the item
   responses agree (r = .56 across the 119 items).
   Fix: the second completion is removed and the first is kept, because the
   first was given without prior exposure to the items. The final sample is
   832 participants, one row per person.
   The row is identified by its Prolific ID where the identifiable source file
   is used, and by its exact start time, end time and duration in the released
   file, whose identifiers are anonymised. Both routes select the same row.

RESPONSE CODING NOTE. Every questionnaire in this file is stored 1-indexed, as
Qualtrics exported it. SSS-8 and PCS run 1-5 where the manuals use 0-4, and
PHQ-9 and GAD-7 run 1-4 where the manuals use 0-3. This is a constant shift of
one point per item. It changes no correlation, z-score, cluster or model
result, but it does shift every reported mean by the item count, and published
clinical cut-offs cannot be applied to these totals without subtracting it.
See CODEBOOK.md for the per-scale table.

The script is idempotent: it builds the corrected table and writes the file
only if something actually changed. Run it once, then run the pipeline
(run_analysis.py).

Usage:
  python 00_score_corrections.py                # corrects public_version/ file
  python 00_score_corrections.py <path.csv>     # corrects the given file
  python 00_score_corrections.py --self-check   # runs the built-in test, writes nothing
================================================================================
"""

import sys

import numpy as np
import pandas as pd

import os

def _default_path():
    """The data file sits next to this script (working repo) or one level up
    (public repo, where the script lives in scripts/); older layouts kept it in
    public_version/. Use whichever exists."""
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, "dfc_interoception_profiling.csv"),
                 os.path.join(here, "..", "dfc_interoception_profiling.csv"),
                 os.path.join(here, "public_version", "dfc_interoception_profiling.csv")):
        if os.path.exists(cand):
            return os.path.normpath(cand)
    raise SystemExit("dfc_interoception_profiling.csv not found; pass the path explicitly")

DEFAULT = _default_path()

STAI_ITEMS = [f"stai_{i}" for i in range(1, 21)]
MISSED_REVERSE = ["stai_3", "stai_14"]
STAI_MAX_PLUS_ONE = 5  # 1-4 scale, so a reversal is 5 - x

SSS8_ITEMS = [f"sss_8_{i}" for i in range(1, 9)]
PCS_ITEMS = [f"pcs_{i}" for i in range(1, 14)]

GAD7_ITEMS = [f"gad_7_{i}" for i in range(1, 8)]
GAD7_GAP_MAP = {1: 1, 3: 2, 4: 3, 5: 4}
GAD7_BROKEN_VALUES = {1, 3, 4, 5}
GAD7_FIXED_VALUES = {1, 2, 3, 4}

FLAG_COL = "scoring_flag"

# The second completion of the one participant who took the survey twice.
# Signature = (startdate, enddate, duration in seconds), which survives the
# anonymisation of the released file. Derived from the Prolific ID in the
# identifiable source file; the ID-based rule below selects the same row.
DUPLICATE_SECOND_COMPLETIONS = [(45574.75612, 45574.76398, 679)]
PROLIFIC_ID_COL = "q1_2"
ANONYMISED = {"anonymized", "anonymised", "nan", ""}


def item_total_r(df, item, items, total):
    """Corrected item-total correlation, item against the sum of the others."""
    rest = df[total] - df[item]
    return np.corrcoef(df[item], rest)[0, 1]


def fix_stai(df, log):
    """Reverse the two items preprocessing missed, then recompute the total."""
    r3 = item_total_r(df, "stai_3", STAI_ITEMS, "stai")
    r14 = item_total_r(df, "stai_14", STAI_ITEMS, "stai")
    log(f"STAI item-total r before: stai_3 = {r3:+.2f}, stai_14 = {r14:+.2f}")
    if not (r3 < 0 and r14 < 0):
        log("  STAI already corrected; nothing to do")
        return df
    for it in MISSED_REVERSE:
        df[it] = STAI_MAX_PLUS_ONE - df[it]
    old_mean = df["stai"].mean()
    df["stai"] = df[STAI_ITEMS].sum(axis=1).astype(int)
    r3 = item_total_r(df, "stai_3", STAI_ITEMS, "stai")
    r14 = item_total_r(df, "stai_14", STAI_ITEMS, "stai")
    log(f"  reversed stai_3 and stai_14; STAI total mean {old_mean:.2f} -> "
        f"{df['stai'].mean():.2f}; item-total r now {r3:+.2f}, {r14:+.2f}")
    return df


def fix_gad7(df, log):
    """Close the Qualtrics recode gap, then recompute the total.

    The four response options were exported as 1, 3, 4, 5. Remap them onto the
    four contiguous codes 1, 2, 3, 4. Only act when the observed value set is
    exactly the broken one, so a corrected file is left alone.
    """
    observed = set(pd.unique(df[GAD7_ITEMS].values.ravel()))
    observed = {int(v) for v in observed if pd.notna(v)}
    log(f"GAD-7 observed item values: {sorted(observed)}")
    if observed == GAD7_FIXED_VALUES:
        log("  GAD-7 already corrected; nothing to do")
        return df
    if not observed.issubset(GAD7_BROKEN_VALUES):
        raise SystemExit(f"Unexpected GAD-7 value set {sorted(observed)}; stopping.")
    old_mean = df["gad7"].mean()
    for it in GAD7_ITEMS:
        df[it] = df[it].map(GAD7_GAP_MAP).astype(df[it].dtype)
    df["gad7"] = df[GAD7_ITEMS].sum(axis=1).astype(int)
    log(f"  remapped 1,3,4,5 -> 1,2,3,4; GAD-7 total mean {old_mean:.2f} -> "
        f"{df['gad7'].mean():.2f}, range {df['gad7'].min()}-{df['gad7'].max()}")
    return df


def prorate(df, items, total, label, flags, log):
    """Prorate totals with isolated missing items. Return the all-missing mask."""
    n_items = len(items)
    n_missing = df[items].isna().sum(axis=1)
    partial = (n_missing > 0) & (n_missing < n_items)
    whole = n_missing == n_items
    for idx in df.index[partial]:
        answered = df.loc[idx, items].dropna()
        new = int(round(answered.mean() * n_items))
        log(f"  {label}: row {idx} omitted {int(n_missing[idx])} item(s); "
            f"total {df.loc[idx, total]} -> {new} (prorated)")
        df.loc[idx, total] = new
        flags.setdefault(idx, []).append(
            f"{total}_prorated_from_{n_items - int(n_missing[idx])}_of_{n_items}_items")
    return whole


def fix_missing(df, log):
    """Prorate partial scales, impute the one all-missing scale, and flag both."""
    log("Missing item handling:")
    flags = {}
    whole_sss = prorate(df, SSS8_ITEMS, "sss8", "SSS-8", flags, log)
    whole_pcs = prorate(df, PCS_ITEMS, "pcs", "PCS", flags, log)
    if whole_pcs.any():
        raise SystemExit("A participant omitted the whole PCS; not expected, stopping.")
    if whole_sss.any():
        others = df.loc[~whole_sss, "sss8"]
        fill = int(round(others.mean()))
        for idx in df.index[whole_sss]:
            log(f"  SSS-8: row {idx} omitted all 8 items; total "
                f"{df.loc[idx, 'sss8']} -> {fill} (sample mean of the other {len(others)})")
            df.loc[idx, "sss8"] = fill
            flags.setdefault(idx, []).append("sss8_mean_imputed_all_8_items_missing")
    if (df["sss8"] < len(SSS8_ITEMS)).any():
        raise SystemExit("SSS-8 totals below scale minimum remain; stopping.")

    # Rebuild the flag column from scratch every run, so it stays deterministic.
    df[FLAG_COL] = ""
    for idx, notes in flags.items():
        df.loc[idx, FLAG_COL] = ";".join(notes)
    log(f"  flagged {len(flags)} row(s) in the '{FLAG_COL}' column")
    return df


def drop_duplicate_completion(df, log):
    """Remove second completions by a participant who took the survey twice.

    Two routes, both selecting the same row. Where real Prolific IDs are
    present, any ID occurring more than once keeps only its earliest start
    time. Where the identifiers are anonymised, rows matching a recorded
    (startdate, enddate, duration) signature are removed. Doing nothing when
    neither applies keeps the function idempotent.
    """
    log("Duplicate completions:")
    drop = set()

    ids = None
    if PROLIFIC_ID_COL in df.columns:
        ids = df[PROLIFIC_ID_COL].astype(str).str.strip().str.lower()
        if ids.isin(ANONYMISED).all():
            ids = None
    if ids is not None:
        for pid, rows in df.groupby(ids).groups.items():
            rows = list(rows)
            if len(rows) < 2:
                continue
            keep = df.loc[rows, "startdate"].astype(float).idxmin()
            for r in rows:
                if r != keep:
                    drop.add(r)
                    log(f"  Prolific ID {pid[:8]}...: row {r} is a repeat completion "
                        f"(start {df.loc[r, 'startdate']}, {df.loc[r, 'duration (in seconds)']} s); "
                        f"keeping row {keep} (start {df.loc[keep, 'startdate']})")
    elif {"startdate", "enddate", "duration (in seconds)"} <= set(df.columns):
        for sd, ed, du in DUPLICATE_SECOND_COMPLETIONS:
            hit = df.index[
                np.isclose(df["startdate"].astype(float), sd)
                & np.isclose(df["enddate"].astype(float), ed)
                & (df["duration (in seconds)"].astype(float) == du)
            ]
            for r in hit:
                drop.add(r)
                log(f"  row {r} matches a recorded repeat-completion signature "
                    f"(start {sd}, {du} s)")

    if not drop:
        log("  none found; nothing to do")
        return df
    df = df.drop(index=sorted(drop)).reset_index(drop=True)
    log(f"  removed {len(drop)} row(s); N is now {len(df)}")
    return df


def correct(df, log=print):
    """Apply every correction. Pure: takes a frame, returns a new frame.

    The duplicate completion goes first, so that the row positions reported by
    the missing-data step refer to the final file.
    """
    df = df.copy()
    df = drop_duplicate_completion(df, log)
    df = fix_stai(df, log)
    df = fix_gad7(df, log)
    df = fix_missing(df, log)
    return df


def self_check():
    """Build a synthetic file with all three defects and prove each is fixed."""
    rng = np.random.default_rng(0)
    n = 60
    d = {}
    # STAI: build a coherent trait score, then UN-reverse items 3 and 14 so the
    # defect is present exactly as it was in the deposited file.
    trait = rng.integers(1, 5, n)
    for i in range(1, 21):
        noise = rng.integers(-1, 2, n)
        d[f"stai_{i}"] = np.clip(trait + noise, 1, 4)
    for i in (3, 14):
        d[f"stai_{i}"] = 5 - d[f"stai_{i}"]
    # GAD-7 with the recode gap: only the values 1, 3, 4, 5 ever occur.
    for i in range(1, 8):
        d[f"gad_7_{i}"] = rng.choice([1, 3, 4, 5], n)
    # SSS-8 and PCS, complete for now.
    for i in range(1, 9):
        d[f"sss_8_{i}"] = rng.integers(1, 6, n)
    for i in range(1, 14):
        d[f"pcs_{i}"] = rng.integers(1, 6, n)
    df = pd.DataFrame(d).astype(float)
    # A participant who completed the survey twice: row 5 repeats row 4.
    df["startdate"] = 45000.0 + np.arange(n) / 1000.0
    df["enddate"] = df["startdate"] + 0.01
    df["duration (in seconds)"] = 900
    df[PROLIFIC_ID_COL] = [f"pid{i:03d}" for i in range(n)]
    df.loc[5, PROLIFIC_ID_COL] = df.loc[4, PROLIFIC_ID_COL]
    df["stai"] = df[STAI_ITEMS].sum(axis=1)
    df["gad7"] = df[GAD7_ITEMS].sum(axis=1)
    # Reproduce the na.rm = TRUE defect: missing items counted as zero.
    df.loc[0, SSS8_ITEMS] = np.nan            # whole scale missing
    df.loc[1, "sss_8_3"] = np.nan             # one item missing
    df.loc[1, "pcs_7"] = np.nan               # one item missing
    df["sss8"] = df[SSS8_ITEMS].sum(axis=1, min_count=0)
    df["pcs"] = df[PCS_ITEMS].sum(axis=1, min_count=0)

    quiet = lambda *a, **k: None
    out = correct(df, log=quiet)
    # `out` has the repeat completion removed, so compare against the same
    # frame with that row dropped and the index reset.
    ref = df.drop(index=5).reset_index(drop=True)

    # 1. STAI: items 3 and 14 flipped back, every item positively keyed.
    assert (out["stai_3"] == 5 - ref["stai_3"]).all(), "stai_3 was not reversed"
    assert (out["stai_14"] == 5 - ref["stai_14"]).all(), "stai_14 was not reversed"
    for it in ("stai_3", "stai_14"):
        assert item_total_r(out, it, STAI_ITEMS, "stai") > 0, f"{it} still negatively keyed"
    assert (out["stai"] == out[STAI_ITEMS].sum(axis=1)).all(), "STAI total is not the item sum"
    # A correctly keyed item must be left alone.
    assert (out["stai_2"] == ref["stai_2"]).all(), "stai_2 must not change"

    # 2. GAD-7: gap closed, total recomputed, order preserved.
    vals = set(pd.unique(out[GAD7_ITEMS].values.ravel()))
    assert vals == {1.0, 2.0, 3.0, 4.0}, f"GAD-7 values not contiguous: {sorted(vals)}"
    assert (out["gad7"] == out[GAD7_ITEMS].sum(axis=1)).all(), "GAD-7 total is not the item sum"
    assert out["gad7"].corr(ref["gad7"]) > 0.97, "GAD-7 remap changed the rank order"
    assert out["gad7"].mean() < ref["gad7"].mean(), "GAD-7 mean should fall after the remap"

    # 3. Missing data: prorated, imputed, and both flagged.
    assert out.loc[0, "sss8"] >= len(SSS8_ITEMS), "all-missing SSS-8 left below the scale floor"
    assert "mean_imputed" in out.loc[0, FLAG_COL], "imputed row not flagged"
    exp = int(round(ref.loc[1, SSS8_ITEMS].dropna().mean() * 8))
    assert out.loc[1, "sss8"] == exp, f"SSS-8 proration wrong: {out.loc[1, 'sss8']} != {exp}"
    assert "prorated" in out.loc[1, FLAG_COL], "prorated row not flagged"
    assert (out.loc[2:, FLAG_COL] == "").all(), "complete rows must not be flagged"
    assert (out["sss8"] >= 8).all(), "an SSS-8 total is below the scale floor"

    # 4. Idempotent: a second pass must be a no-op.
    again = correct(out, log=quiet)
    assert again.equals(out), "correct() is not idempotent"

    # 5. Duplicate completion: the later row goes, the earlier one stays.
    assert len(out) == n - 1, f"duplicate row not removed: {len(out)} rows"
    kept = out[PROLIFIC_ID_COL].tolist()
    assert kept.count(ref.loc[4, PROLIFIC_ID_COL]) == 1, "duplicate Prolific ID survived"
    assert out[PROLIFIC_ID_COL].is_unique, "Prolific IDs are not unique after the fix"
    assert (out["startdate"] == ref["startdate"].values).all(), "wrong row was dropped"
    assert out.index.equals(pd.RangeIndex(len(out))), "index was not reset"

    print("SELF-CHECK PASSED (duplicate removal, STAI reversal, GAD-7 recode, "
          "proration, flags, idempotency)")


def main(path):
    before = pd.read_csv(path)
    # An empty flag round-trips through CSV as NaN. Normalise it before the
    # no-op comparison, or every run would claim to have rewritten the file.
    if FLAG_COL in before.columns:
        before[FLAG_COL] = before[FLAG_COL].fillna("").astype(str)
    after = correct(before)

    if after.equals(before):
        print("No changes written.")
    else:
        after.to_csv(path, index=False)
        print(f"Wrote corrected file: {path}")

    print(f"After: STAI range {after['stai'].min()}-{after['stai'].max()}, "
          f"GAD-7 range {after['gad7'].min()}-{after['gad7'].max()} (mean "
          f"{after['gad7'].mean():.2f}), "
          f"SSS-8 range {after['sss8'].min()}-{after['sss8'].max()}, "
          f"PCS range {after['pcs'].min()}-{after['pcs'].max()}, N = {len(after)}")


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--self-check":
        self_check()
    else:
        main(args[0] if args else DEFAULT)
