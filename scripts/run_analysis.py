"""
================================================================================
RUN ANALYSIS - Main Pipeline Script
================================================================================

Executes all analysis scripts in order. Scripts are numbered to indicate
the logical flow of the analysis pipeline.

Pipeline:
  01. Normality checks (R) - Data assumption validation
  02. IAS-IATS correlation (R) - Basic bivariate relationship (stats only)
  03. Correlation heatmaps (Python) - Broader correlations with MH outcomes
  04. Cluster analysis (R) - K-selection, k-means, robustness
  05. LPA analysis (R) - Latent profile analysis
  06. Simpson's paradox (Python) - Bayes factor analysis
  07. At-risk subcluster (Python) - Subgroup analysis
  08. Cluster visualizations (Python) - 3D density, combined figures
  09. Moderation + demographics (R) - IAS x IATS interaction with controls
  10. Mediation analysis (R) - Dual-pathway mediation (TAS + Somatic)

Usage:
  python run_analysis.py           # Run all scripts
  python run_analysis.py 03 05     # Run only scripts 03 and 05
  python run_analysis.py --list    # List all scripts
"""

import subprocess
import sys
import os
from pathlib import Path

# Define analysis pipeline
PIPELINE = [
    ("01", "01_normality_checks.R", "R", "Data assumption checks"),
    ("02", "02_correlation_ias_iats.R", "R", "IAS-IATS correlation stats"),
    ("03", "03_correlation_heatmaps.py", "Python", "Correlation heatmaps with MH"),
    ("04", "04_cluster_analysis.R", "R", "K-selection + clustering + robustness"),
    ("05", "05_lpa_analysis.R", "R", "Latent profile analysis"),
    ("06", "06_simpsons_paradox_bayes.py", "Python", "Simpson's paradox Bayes factor"),
    ("07", "07_at_risk_subcluster.py", "Python", "At-risk subgroup analysis"),
    ("08", "08_cluster_visualizations.py", "Python", "3D density + combined figures"),
    ("09", "09_moderation_demographics.R", "R", "Moderation + demographic controls"),
    ("10", "10_mediation_analysis.R", "R", "Dual-pathway mediation"),
]

def run_r_script(script_path):
    """Run an R script using Rscript."""
    cmd = ["Rscript", "--vanilla", script_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result

def run_python_script(script_path):
    """Run a Python script."""
    python_exe = r"C:/code/languages/anaconda3/envs/mipype/python.exe"
    cmd = [python_exe, script_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result

def list_scripts():
    """Print list of all scripts in pipeline."""
    print("\n" + "=" * 70)
    print("ANALYSIS PIPELINE")
    print("=" * 70 + "\n")
    print(f"{'#':<4} {'Script':<35} {'Type':<8} {'Description'}")
    print("-" * 70)
    for num, script, lang, desc in PIPELINE:
        print(f"{num:<4} {script:<35} {lang:<8} {desc}")
    print()

def run_script(num, script, lang, desc, script_dir):
    """Run a single script and report result."""
    script_path = script_dir / script

    if not script_path.exists():
        print(f"  WARNING: {script} not found, skipping")
        return False

    print(f"\n[{num}] Running {script}...")
    print(f"     {desc}")

    if lang == "R":
        result = run_r_script(str(script_path))
    else:
        result = run_python_script(str(script_path))

    if result.returncode == 0:
        print(f"     DONE")
        if result.stdout:
            # Print last few lines of output
            lines = result.stdout.strip().split('\n')
            for line in lines[-3:]:
                print(f"     > {line}")
        return True
    else:
        print(f"     ERROR (exit code {result.returncode})")
        if result.stderr:
            print(f"     {result.stderr[:200]}")
        return False

def main():
    script_dir = Path(__file__).parent

    # Parse arguments
    if len(sys.argv) > 1:
        if sys.argv[1] == "--list":
            list_scripts()
            return

        # Run specific scripts
        selected = set(sys.argv[1:])
        scripts_to_run = [(n, s, l, d) for n, s, l, d in PIPELINE if n in selected]
    else:
        scripts_to_run = PIPELINE

    if not scripts_to_run:
        print("No scripts selected. Use --list to see available scripts.")
        return

    print("\n" + "=" * 70)
    print("INTEROCEPTIVE MODERATION ANALYSIS PIPELINE")
    print("=" * 70)
    print(f"\nRunning {len(scripts_to_run)} scripts...")

    success = 0
    failed = 0

    for num, script, lang, desc in scripts_to_run:
        if run_script(num, script, lang, desc, script_dir):
            success += 1
        else:
            failed += 1

    print("\n" + "=" * 70)
    print(f"COMPLETE: {success} succeeded, {failed} failed")
    print("=" * 70 + "\n")

if __name__ == "__main__":
    main()
