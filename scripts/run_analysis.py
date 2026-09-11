"""
================================================================================
RUN ANALYSIS - Main Pipeline Script
================================================================================

Executes all analysis scripts in order. Scripts are numbered to indicate
the logical flow of the analysis pipeline.

Pipeline:
  00b. Build derived data (Python) - rebuilds the two supplementary_data files
       that later scripts read, from the corrected master data file
  01.  Normality checks (R) - data-assumption validation
  02.  IAS-IATS correlation (R) - bivariate relationship (stats only)
  03.  Correlation heatmaps (Python) - correlations with mental-health outcomes
  04.  Cluster analysis (R) - K-selection, k-means, robustness
  05.  LPA analysis (R) - latent profile analysis
  06.  Simpson's paradox (Python) - Bayes factor analysis
  07.  At-risk subcluster (Python) - subgroup analysis
  08.  Cluster visualizations (Python) - 3D density, combined figures
  09.  Item selection (R) - minimal item-set selection
  10.  Network analysis (R) - network community structure
  11.  Combined figure (R) - main figure assembly
  12.  SEM pathway (R) - dual-pathway structural equation model
  12b. Pathway dominance test (R)
  12c. Multi-group omnibus invariance test (R)
  12d. Pathway dominance figure (R)
  12e. Comprehensive moderation figure (R)
  13.  TAS-20 subscale (R) - subscale mediation
  13b. TAS-20 subscale path diagrams (R)
  14.  VVIQ mediation (R) - VVIQ-interoception mediation
  15.  Item-level clustering sensitivity (Python; needs scikit-learn-extra,
       which requires NumPy < 2)
  16.  Alexithymia-somatic pathway (Python)
  17.  Revision sensitivity (R) - age/gender covariates, SSS-8-only mediator
  18.  Somatic composite factor structure (R)
  19.  Exact test statistics for reporting (R)
  s8, s11, s12, s13, s14, s18, s24. Supplementary figure builders, run last
       because they read what the numbered steps write.

  00_score_corrections.py is not part of the pipeline. It was applied once to
  the deposited data file (see its docstring) and is kept for provenance.

Usage:
  python run_analysis.py           # Run all scripts
  python run_analysis.py 03 05     # Run only scripts 03 and 05
  python run_analysis.py s18       # Run one supplementary figure builder
  python run_analysis.py --list    # List all scripts

Interpreters:
  R is called as "Rscript" and Python as the interpreter running this file.
  Neither path is hard-coded. Override either one if it is not what you want:
    set PYTHON_EXE=C:\\path\\to\\python.exe
    set RSCRIPT=C:\\path\\to\\Rscript.exe
  Step 15 needs scikit-learn-extra, which requires NumPy < 2. If that lives in
  a separate environment, PYTHON_EXE is how you point at it.
"""

import subprocess
import sys
import os
from pathlib import Path

# Define analysis pipeline
PIPELINE = [
    ("00b", "00b_build_derived_data.py", "Python", "Rebuild derived data files from the corrected master"),
    ("01", "01_normality_checks.R", "R", "Normality / data-assumption checks"),
    ("02", "02_correlation_ias_iats.R", "R", "IAS-IATS correlation stats"),
    ("03", "03_correlation_heatmaps.py", "Python", "Correlation heatmaps with mental-health outcomes"),
    ("04", "04_cluster_analysis.R", "R", "K-selection, k-means clustering, robustness"),
    ("05", "05_lpa_analysis.R", "R", "Latent profile analysis"),
    ("06", "06_simpsons_paradox_bayes.py", "Python", "Simpson's paradox Bayes factor"),
    ("07", "07_at_risk_subcluster.py", "Python", "At-risk subcluster analysis"),
    ("08", "08_cluster_visualizations.py", "Python", "Cluster visualizations (3D density, combined)"),
    ("09", "09_item_selection_analysis.R", "R", "Minimal item-set selection"),
    ("10", "10_network_analysis.R", "R", "Network community analysis"),
    ("11", "11_figure_combined.R", "R", "Combined main figure assembly"),
    ("19", "19_ert_exact_statistics.R", "R", "Exact test statistics for reporting (writes the values used by the 12 and 12b figures)"),
    ("12", "12_sem_pathway_analysis.R", "R", "Dual-pathway SEM analysis"),
    ("12b", "12b_pathway_dominance_test.R", "R", "Pathway dominance test"),
    ("12c", "12c_multigroup_omnibus_test.R", "R", "Multi-group omnibus invariance test"),
    ("12d", "12d_pathway_dominance_figure.R", "R", "Pathway dominance figure"),
    ("12e", "12e_comprehensive_moderation_figure.R", "R", "Comprehensive moderation figure"),
    ("13", "13_tas_subscale_analysis.R", "R", "TAS-20 subscale mediation"),
    ("13b", "13b_tas_subscale_path_diagram.R", "R", "TAS-20 subscale path diagrams"),
    ("14", "14_vviq_interoception_mediation.R", "R", "VVIQ-interoception mediation"),
    ("15", "15_item_level_clustering_sensitivity.py", "Python", "Item-level clustering sensitivity"),
    ("16", "16_alexithymia_somatic_pathway.py", "Python", "Alexithymia-somatic pathway decomposition"),
    ("17", "17_revision_sensitivity.R", "R", "Revision sensitivity: covariates, SSS-8 only"),
    ("18", "18_somatic_factor_structure.R", "R", "Somatic composite factor structure"),

    # Supplementary figure builders. They run last because they read what the
    # numbered steps write: s8 needs the item ranking from 09, s18 needs the
    # bootstrap table from 19. Keeping them in the pipeline is what stops a
    # deposited figure from surviving a change to the analysis underneath it.
    ("s8", "supplementary_figure_scripts/sem_sensitivity_analysis.R", "R",
     "Supp. Fig. 8: SEM sensitivity, full 42 items vs minimal 21"),
    ("s11", "supplementary_figure_scripts/create_simple_slopes_continuous.R", "R",
     "Supp. Fig. 11: simple slopes for continuous moderation"),
    ("s12", "supplementary_figure_scripts/create_simple_slopes_indirect_continuous.R", "R",
     "Supp. Fig. 12: conditional indirect effects by IATS level"),
    ("s13", "supplementary_figure_scripts/create_3d_combined_wireframe.R", "R",
     "Supp. Fig. 13: three-dimensional interaction surfaces"),
    ("s14", "supplementary_figure_scripts/create_simple_slopes_indirect_profiles.R", "R",
     "Supp. Fig. 14: indirect effects by interoceptive profile"),
    ("s18", "supplementary_figure_scripts/make_subscale_figure.py", "Python",
     "Supp. Fig. 18: TAS-20 subscale decomposition (needs 19)"),
    ("s24", "supplementary_figure_scripts/make_somatic_mediator_figure.py", "Python",
     "Supp. Fig. 24: coefficients by definition of the somatic mediator"),

    # Journal source-data workbooks, one per main figure, built from the
    # supplementary_data files the figure scripts read and write. Needs Rscript
    # for the Figure 3D MDS coordinates.
    ("sd", "make_source_data.py", "Python",
     "Source data workbooks for Figures 1 to 3 (source_data/)"),
]

def run_r_script(script_path):
    """Run an R script using Rscript."""
    rscript = os.environ.get("RSCRIPT", "Rscript")
    cmd = [rscript, "--vanilla", script_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result

def run_python_script(script_path):
    """Run a Python script with the interpreter that is running this file.

    Step 15 needs scikit-learn-extra, which requires NumPy < 2. If that lives in
    a different environment from the one you start the pipeline with, point
    PYTHON_EXE at it. Same idea for RSCRIPT if Rscript is not on PATH.
    """
    python_exe = os.environ.get("PYTHON_EXE", sys.executable)
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

    # Every script resolves its data, output and plot paths relative to the
    # repository root, so run them from there no matter where this was invoked.
    # Two layouts are supported: the working repo, where the scripts sit beside
    # the data file, and the public release, where they sit in scripts/.
    DATA = "dfc_interoception_profiling.csv"
    for candidate in (script_dir, script_dir.parent):
        if (candidate / DATA).exists():
            os.chdir(candidate)
            break
    else:
        sys.exit(f"Cannot find {DATA} in {script_dir} or {script_dir.parent}. "
                 f"Run this from a complete checkout of the repository.")

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
