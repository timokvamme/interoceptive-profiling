# Supplementary figure map

The file names in this folder use a per-analysis prefix (`1_s_*` profiles,
`2_s_*` item selection, `3_s_*` dual-pathway model). Those prefixes record which
analysis produced the figure. They are **not** the numbers used in the paper.

The Supplementary Information numbers the figures S1 to S24 in reading order.
Use this table to go from the paper to the file.

| Supplementary Information | File | Written by |
|---|---|---|
| S1. Cluster validation and visualization | `1_s_2_k_selection_cluster.png` | `08_cluster_visualizations.py` |
| S2. Level-dependent bias in IAS-IATS correlation | `1_s_6_simpsons_paradox_bayes.png` | `06_simpsons_paradox_bayes.py` |
| S3. Latent profile analysis | `1_s_4_lpa_profiles.png` | `05_lpa_analysis.R` |
| S4. Bootstrap stability analysis | `1_s_3_cluster_robustness.png` | `04_cluster_analysis.R` |
| S5. Demographic characteristics by profile | `1_s_5_demographics_by_profile.png` | `04_cluster_analysis.R` |
| S6. Network community structure | `2_s_2_network_communities.png` | `11_figure_combined.R` |
| S7. Cluster preservation with reduced item sets | `2_s_3_cluster_preservation.png` | `09_item_selection_analysis.R` |
| S8. SEM sensitivity analysis | `3_s_1_sem_sensitivity_comparison.png` | `sem_sensitivity_analysis.R` |
| S9. Distribution histograms for interoceptive scales | `1_s_8_normality_histograms.png` | `01_normality_checks.R` |
| S10. Q-Q plots for normality assessment | `1_s_9_normality_qq_plots.png` | `01_normality_checks.R` |
| S11. Simple slopes for continuous moderation | `3_s_3_simple_slopes_continuous.png` | `create_simple_slopes_continuous.R` |
| S12. Conditional indirect effects by IATS level | `3_s_4_simple_slopes_indirect_continuous.png` | `create_simple_slopes_indirect_continuous.R` |
| S13. Three-dimensional interaction surfaces | `3_s_6_3d_interaction_surfaces.png` | `create_3d_combined_wireframe.R` |
| S14. Indirect effects by interoceptive profile | `3_s_5_simple_slopes_indirect_profiles.png` | `create_simple_slopes_indirect_profiles.R` |
| S15. At-risk subcluster validation | `1_s_7_at_risk_cluster_analysis.png` | `07_at_risk_subcluster.py` |
| S16. Item-level clustering sensitivity | `1_s_9_item_level_clustering_sensitivity.png` | `15_item_level_clustering_sensitivity.py` |
| S17. Pathway dominance analysis | `3_s_2_pathway_dominance_statistical.png` | `12b_pathway_dominance_test.R` |
| S18. TAS-20 subscale decomposition of the alexithymia pathway | `3_s_8_tas_subscale_decomposition.png` | `make_subscale_figure.py` |
| S19. Dual-pathway SEM diagrams by TAS-20 subscale | `3_s_9_tas_subscale_path_diagrams.png` | `13b_tas_subscale_path_diagram.R` |
| S20. Comprehensive pathway comparison | `3_s_7_comprehensive_moderated_mediation.png` | `12e_comprehensive_moderation_figure.R` |
| S21. Item-level correlations with mental health | `2_s_1_item_mh_correlations.png` | `02_correlation_ias_iats.R` |
| S22. Profile transition patterns | `2_s_4_transition_heatmap.png` | `09_item_selection_analysis.R` |
| S23. Alexithymia facets and the somatic pathway | `3_s_10_alexithymia_somatic_pathway.png` | `16_alexithymia_somatic_pathway.py` |
| S24. Dual-pathway coefficients by definition of the somatic mediator | `1_s_24_somatic_mediator_definitions.png` | `make_somatic_mediator_figure.py` |

Scripts named above live in `../../scripts/`, except the seven figure builders
in `../../supplementary_figure_scripts/`.

## Files in this folder that the paper does not print

`1_s_1_k_selection_methods.png` is the three-panel k-selection figure produced by
`04_cluster_analysis.R`. Script `08_cluster_visualizations.py` crops three panels
out of it and combines them with the 3D cluster plot to build S1. It is kept here
because the pipeline writes it, and because S1 cannot be rebuilt without it.

## Note on the two files numbered `1_s_9`

`1_s_9_normality_qq_plots.png` and `1_s_9_item_level_clustering_sensitivity.png`
share a prefix number. They are different figures, S10 and S16. The shared number
is an artefact of the file naming, not an error in the paper. Use this table, not
the prefix.

## Figures outside this folder

`../other_plots_vviq_plots/` holds the VVIQ analyses. The paper does not print
them. They are kept because the deposited scripts produce them.
