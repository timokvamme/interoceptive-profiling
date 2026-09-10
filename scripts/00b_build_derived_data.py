"""
================================================================================
00b - BUILD DERIVED DATA FILES
================================================================================

Two files under supplementary_data/ are READ by the pipeline but were never
WRITTEN by it. They were hand-made snapshots, so every time the master data
file changed they went stale and the scripts that read them silently used old
numbers:

  supplementary_data/data_with_ias_iats_ratio.csv
      Read by 18_somatic_factor_structure.R. It is the master file plus six
      derived columns. Before this script existed it still carried the broken
      GAD-7 coding (values 1, 3, 4, 5), so the somatic factor analysis and the
      mental-health composite in script 18 ran on uncorrected data.

  supplementary_data/interoception_item_data.csv
      Read by 02_correlation_ias_iats.R. It holds the 42 IAS and IATS item
      wordings plus each item's correlation with the mental-health composite.
      The wordings are fixed metadata, but the correlations depend on the
      composite, which changed when STAI and GAD-7 were corrected. Script 02
      plots the stored correlation, so a published figure was drawn from stale
      values.

This script rebuilds both from the corrected master, so the pipeline is
self-contained. Item wordings and the two network-derived columns
(edges_to_mh, top_10_predictor) are metadata that no correction can change,
so they are carried over from the existing file rather than invented.

Usage:
  python 00b_build_derived_data.py              # rebuild both derived files
  python 00b_build_derived_data.py --self-check # run the built-in test, write nothing
================================================================================
"""

import os
import sys

import numpy as np
import pandas as pd

MASTER = "dfc_interoception_profiling.csv"

# Resolve the repository root: this script sits beside the data file in the
# working repo, and in scripts/ in the public release. Run from either.
_here = os.path.dirname(os.path.abspath(__file__))
for _c in (_here, os.path.dirname(_here)):
    if os.path.exists(os.path.join(_c, MASTER)):
        os.chdir(_c)
        break
SUPP = "supplementary_data"
RATIO_FILE = f"{SUPP}/data_with_ias_iats_ratio.csv"
ITEM_FILE = f"{SUPP}/interoception_item_data.csv"

IAS_ITEMS = [f"ias_{i}" for i in range(1, 22)]
IATS_ITEMS = [f"iats_{i}" for i in range(1, 22)]

# The mental-health composite, defined exactly as 02_correlation_ias_iats.R
# defines it: the row mean of the six z-scored outcome totals.
OUTCOME_VARS = ["tas", "phq9", "gad7", "stai", "sss8", "pcs"]

# Carried straight over from the previous file. These depend only on variables
# no correction touches (the VVIQ items), so recomputing them is unnecessary
# and would mean guessing at a derivation the repository does not record.
CARRIED_RATIO_COLS = ["vviq_k", "vviq_q"]
CARRIED_ITEM_COLS = ["item", "item_label", "scale", "questionnaire_text",
                     "edges_to_mh", "top_10_predictor"]


def mh_composite(df):
    """Row mean of the z-scored outcome totals. Mirrors rowMeans(scale(...))."""
    z = (df[OUTCOME_VARS] - df[OUTCOME_VARS].mean()) / df[OUTCOME_VARS].std(ddof=1)
    return z.mean(axis=1)


ROW_KEY = ["startdate", "enddate", "duration (in seconds)"]


def align_old_to_master(master, old):
    """Subset the old derived file to the participants still in the master.

    A participant removed from the master (the repeat completion, see
    00_score_corrections.py) leaves the old derived file one row longer. Match
    on the response timestamps, which are stable across every copy of the data,
    and return the old rows in the master's order.
    """
    if len(old) == len(master):
        return old
    if not set(ROW_KEY) <= set(master.columns) | set(old.columns):
        raise SystemExit("Cannot align: the timestamp key columns are missing")
    key = lambda d: list(zip(*[d[c].astype(float).round(6) for c in ROW_KEY]))
    pos = {k: i for i, k in enumerate(key(old))}
    idx = [pos.get(k) for k in key(master)]
    if any(i is None for i in idx):
        raise SystemExit("Cannot align: a master row has no match in the old file")
    print(f"Old derived file has {len(old)} rows for {len(master)} participants; "
          f"aligned on {', '.join(ROW_KEY)} and dropped "
          f"{len(old) - len(master)} row(s) no longer in the master.")
    return old.iloc[idx].reset_index(drop=True)


def check_alignment(master, old):
    """Prove the carried rows belong to the same participants, in the same order.

    The item responses are the fingerprint: if all 64 IAS, IATS and VVIQ item
    columns match row for row, the two frames are the same people in the same
    order. Age and gender are checked as an independent second key.
    """
    if len(master) != len(old):
        raise SystemExit(f"Row count differs: master {len(master)}, old {len(old)}")
    vviq = [c for c in master.columns if c.startswith("vviq_q")]
    keys = IAS_ITEMS + IATS_ITEMS + vviq + ["age", "gender"]
    bad = [c for c in keys if c in old.columns and not old[c].equals(master[c])]
    if bad:
        raise SystemExit(f"Row alignment failed; these columns differ: {bad[:5]}")
    return len(keys)


def build_ratio(master, old):
    """Master file plus the six derived columns, rebuilt from corrected data."""
    out = master.copy()
    out["cohort"] = np.where(out["sass"].notna(), "control", "aphan")
    out["vviq_k"] = old["vviq_k"].values
    out["vviq_q"] = old["vviq_q"].values
    out["ias_total"] = out["ias"] / len(IAS_ITEMS)
    out["iats_total"] = out["iats"] / len(IATS_ITEMS)
    out["ias_iats_ratio"] = out["ias_total"] / out["iats_total"]
    return out


def build_items(master, old):
    """Item metadata kept, correlations recomputed against the fixed composite."""
    mh = mh_composite(master)
    out = old[CARRIED_ITEM_COLS].copy()
    out["correlation_with_mh_composite"] = [master[i].corr(mh) for i in out["item"]]
    out["abs_correlation_with_mh"] = out["correlation_with_mh_composite"].abs()
    # Preserve the original column order of the file the pipeline expects.
    return out[["item", "item_label", "scale", "questionnaire_text",
                "correlation_with_mh_composite", "abs_correlation_with_mh",
                "edges_to_mh", "top_10_predictor"]]


def self_check():
    """Prove the composite and the derived columns behave as the scripts expect."""
    rng = np.random.default_rng(1)
    n = 40
    d = {v: rng.integers(10, 40, n).astype(float) for v in OUTCOME_VARS}
    for i in IAS_ITEMS + IATS_ITEMS:
        d[i] = rng.integers(1, 6, n).astype(float)
    df = pd.DataFrame(d)
    df["ias"] = df[IAS_ITEMS].sum(axis=1)
    df["iats"] = df[IATS_ITEMS].sum(axis=1)
    df["sass"] = [1.0] * 25 + [np.nan] * 15
    df["age"] = rng.integers(18, 70, n)
    df["gender"] = rng.integers(1, 3, n)
    df["vviq_q2_1_1"] = rng.integers(1, 6, n)

    mh = mh_composite(df)
    assert abs(mh.mean()) < 1e-9, "composite should be centred on zero"
    # A composite built from z-scores must be invariant to a constant shift in
    # any component. This is exactly why the 1-indexed coding is harmless here.
    shifted = df.copy()
    shifted["phq9"] = shifted["phq9"] + 9
    assert np.allclose(mh, mh_composite(shifted)), "composite is not shift-invariant"

    old = pd.DataFrame({"vviq_k": rng.integers(1, 6, n), "vviq_q": rng.integers(1, 6, n),
                        "age": df["age"], "gender": df["gender"]})
    for c in IAS_ITEMS + IATS_ITEMS + ["vviq_q2_1_1"]:
        old[c] = df[c]
    n_keys = check_alignment(df, old)
    assert n_keys == 45, f"expected 45 alignment keys, got {n_keys}"

    r = build_ratio(df, old)
    assert np.allclose(r["ias_total"], df["ias"] / 21), "ias_total is not the per-item mean"
    assert np.allclose(r["ias_iats_ratio"], r["ias_total"] / r["iats_total"]), "ratio wrong"
    assert (r.loc[df["sass"].notna(), "cohort"] == "control").all(), "cohort mislabelled"
    assert (r.loc[df["sass"].isna(), "cohort"] == "aphan").all(), "cohort mislabelled"

    # Misaligned input must be refused, not silently carried over.
    broken = old.copy()
    broken.loc[0, "ias_1"] = broken.loc[0, "ias_1"] + 1
    try:
        check_alignment(df, broken)
    except SystemExit:
        pass
    else:
        raise AssertionError("check_alignment did not catch a misaligned frame")

    print("SELF-CHECK PASSED (composite, shift-invariance, derived columns, alignment guard)")


def main():
    master = pd.read_csv(MASTER)

    old_ratio = align_old_to_master(master, pd.read_csv(RATIO_FILE))
    n_keys = check_alignment(master, old_ratio)
    print(f"Row alignment verified against {n_keys} key columns.")

    ratio = build_ratio(master, old_ratio)
    ratio.to_csv(RATIO_FILE, index=False)
    gad = sorted(pd.unique(ratio[[f"gad_7_{i}" for i in range(1, 8)]].values.ravel()))
    print(f"Wrote {RATIO_FILE}: {ratio.shape[0]} rows x {ratio.shape[1]} cols; "
          f"GAD-7 item values now {[int(v) for v in gad]}")

    old_items = pd.read_csv(ITEM_FILE)
    items = build_items(master, old_items)
    items.to_csv(ITEM_FILE, index=False)
    moved = (items["correlation_with_mh_composite"].values
             - old_items["correlation_with_mh_composite"].values)
    print(f"Wrote {ITEM_FILE}: {len(items)} items; correlations with the "
          f"mental-health composite recomputed, largest change "
          f"{np.abs(moved).max():.4f}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--self-check":
        self_check()
    else:
        main()
