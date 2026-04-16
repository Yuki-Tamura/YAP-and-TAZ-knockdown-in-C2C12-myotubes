# GitHub Export: siRNA Specificity and RNA-seq Support Scripts

This folder contains copies of the scripts used for the siRNA specificity and related RNA-seq support analyses in this project, renamed with more descriptive file names for easier GitHub publication.

The original scripts remain unchanged in their original locations under `/Users/Yuki/Desktop/siRNA/analysis/`.

## Recommended script order

1. `00_run_full_rnaseq_and_sirna_workflow.R`
   Main driver for the default dataset in the current project directory.

2. `01_run_workflow_on_alternate_count_matrix.R`
   Runs the same workflow on an alternate count matrix and writes results to a user-specified analysis subdirectory.

3. `02_workflow_helper_functions.R`
   Shared helper functions used by the main workflow scripts.

4. `10_build_sirna_seed_and_utr_annotations.py`
   Builds siRNA strand definitions, searches near-perfect transcriptome matches, reconstructs transcript-level 3'UTRs, and annotates seed matches.

5. `11_analyze_single_knockdown_seed_shift_cdfs.R`
   Generates single-knockdown seed-match CDF analyses and related summary tables and figures.

6. `12_analyze_double_knockdown_seed_shift_cdfs.R`
   Generates double-knockdown CDF analyses for YAP/TAZ guide and passenger seed categories, including combined categories such as both guides and both passengers.

7. `13_summarize_cdf_shift_metrics_across_seed_lengths.R`
   Summarizes extended CDF metrics across seed5mer to seed8mer outputs, including delta median, area between CDFs, tail enrichment, and KS D values.

8. `14_compare_observed_seed_effects_against_null_seeds.R`
   Performs null-seed sensitivity analyses by comparing the observed seed effect against matched random control seeds.

9. `15_test_on_target_confounding_of_seed_effects.R`
   Performs sensitivity analyses to test whether apparent seed-associated repression is inflated by likely on-target biology.

10. `16_generate_seed_length_specific_outputs.py`
    Regenerates siRNA specificity outputs for multiple seed-length definitions such as 5mer, 6mer, 7mer, and 8mer.

11. `17_recompute_seed_annotations_for_selected_seed_length.py`
    Recomputes seed annotations for a selected seed-length setting.

## Original-to-export mapping

| Exported file | Original file |
| --- | --- |
| `00_run_full_rnaseq_and_sirna_workflow.R` | `analysis/run_workflow.R` |
| `01_run_workflow_on_alternate_count_matrix.R` | `analysis/run_dataset_workflow.R` |
| `02_workflow_helper_functions.R` | `analysis/R/helpers.R` |
| `10_build_sirna_seed_and_utr_annotations.py` | `analysis/python/sirna_specificity.py` |
| `11_analyze_single_knockdown_seed_shift_cdfs.R` | `analysis/R/sirna_seed_stats.R` |
| `12_analyze_double_knockdown_seed_shift_cdfs.R` | `analysis/R/double_knockdown_seed_cdf.R` |
| `13_summarize_cdf_shift_metrics_across_seed_lengths.R` | `analysis/R/cdf_metric_extensions.R` |
| `14_compare_observed_seed_effects_against_null_seeds.R` | `analysis/R/null_seed_sensitivity.R` |
| `15_test_on_target_confounding_of_seed_effects.R` | `analysis/R/sensitivity_on_target_confounding.R` |
| `16_generate_seed_length_specific_outputs.py` | `analysis/python/generate_seed_length_variants.py` |
| `17_recompute_seed_annotations_for_selected_seed_length.py` | `analysis/python/recompute_seed_annotations.py` |

## Notes

- These are copied export files intended for clearer public sharing.
- Output paths inside the scripts still reflect the current project structure unless you modify them for a standalone repository layout.
- If you want, a next step would be to refactor these into a cleaner GitHub-ready package structure such as `scripts/`, `src/`, and `README` with usage examples.
