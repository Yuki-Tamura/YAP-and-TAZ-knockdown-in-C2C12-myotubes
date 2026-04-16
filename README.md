# YAP-and-TAZ-knockdown-in-C2C12-myotubes

This repository contains analysis scripts associated with the study of YAP and TAZ knockdown in differentiated C2C12 myotubes.

The code is organized into three analysis modules:

- **Calcium kinetics**: analysis of electrically evoked Ca2+ responses in C2C12 myotubes
- **MitobrightRed to Calcein signal ratio**: image-based normalization of mitochondrial signal
- **siRNA specificity**: in silico and RNA-seq-based analyses related to siRNA specificity and seed-associated transcriptome shifts

This repository is intended as a transparent record of the analysis workflows used in the study. It should be regarded as a project code repository rather than as a fully automated end-to-end software package.

---

## Repository structure

```text
.
├── Calcium kinetics/
├── MitobrightRed to Calcein signal ratio/
└── siRNA specificity/
```

---

## Module overview

### 1. Calcium kinetics

This folder contains scripts for Ca2+ imaging analysis in electrically stimulated C2C12 myotubes.

The folder currently includes:

- `README.md`
- `calculate_calcium_kinetics_from_excel.py`

According to the folder-specific documentation, this module is used to calculate replicate-level Ca2+ kinetic parameters from fluorescence time-series data exported to Excel workbooks.

Typical output metrics include:

- mean rise slope
- mean decay slope
- mean time to peak
- mean apparent decay tau
- number of detected events per replicate

Please see the folder-specific `README.md` for expected input layout and usage details.

### 2. MitobrightRed to Calcein signal ratio

This folder contains a simple image-based workflow to estimate normalized mitochondrial signal from TIFF images.

The folder currently includes:

- `README.md`
- `calculate_c003_to_c002_signal_ratios.py`

In this workflow:

- `C002T001` is treated as the Calcein channel
- `C003T001` is treated as the MitoBrightRed channel
- the mean MitoBrightRed intensity is divided by the mean Calcein intensity

The output is intended as a **field-level normalized proxy** for mitochondrial signal, not a segmented per-cell measurement.

### 3. siRNA specificity

This folder contains R and Python scripts for siRNA specificity analysis and RNA-seq support analyses.

The folder currently includes:

- `00_run_full_rnaseq_and_sirna_workflow.R`
- `01_run_workflow_on_alternate_count_matrix.R`
- `02_workflow_helper_functions.R`
- `10_build_sirna_seed_and_utr_annotations.py`
- `11_analyze_single_knockdown_seed_shift_cdfs.R`
- `12_analyze_double_knockdown_seed_shift_cdfs.R`
- `13_summarize_cdf_shift_metrics_across_seed_lengths.R`
- `14_compare_observed_seed_effects_against_null_seeds.R`
- `15_test_on_target_confounding_of_seed_effects.R`
- `16_generate_seed_length_specific_outputs.py`
- `17_recompute_seed_annotations_for_selected_seed_length.py`
- `README.md`

These scripts support analyses such as:

- siRNA seed and 3′UTR annotation
- seed-match-based CDF analyses
- single- and double-knockdown transcriptome shift analysis
- sensitivity analyses using null seeds
- evaluation of possible on-target confounding

Please see the folder-specific `README.md` for the recommended script order and project-specific notes.

---

## Software

This repository contains both **Python** and **R** scripts.

### Python
Different modules use packages such as:

- `numpy`
- `openpyxl`
- `Pillow`
- `imageio`

Some workflows may also use **FFmpeg**, as described in the folder-level documentation.

### R
The siRNA specificity workflow is implemented mainly in **R**, with Python helper scripts used for annotation and preprocessing.

Please check the `README.md` inside each subdirectory for module-specific requirements and usage details.

---

