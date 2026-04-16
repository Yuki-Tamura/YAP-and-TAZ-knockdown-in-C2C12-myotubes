#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys


ROOT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.abspath(".")
ANALYSIS_ROOT = os.path.join(ROOT, "analysis")
PY_DIR = os.path.join(ANALYSIS_ROOT, "python")
R_SCRIPT = os.path.join(ANALYSIS_ROOT, "R", "sirna_seed_stats.R")
RECOMPUTE_SCRIPT = os.path.join(PY_DIR, "recompute_seed_annotations.py")

DATASETS = [
    {
        "analysis_dir": ANALYSIS_ROOT,
        "base_spec_dir": os.path.join(ANALYSIS_ROOT, "04_sirna_specificity"),
    },
    {
        "analysis_dir": os.path.join(ANALYSIS_ROOT, "alt_exclude_TAZ1"),
        "base_spec_dir": os.path.join(ANALYSIS_ROOT, "alt_exclude_TAZ1", "04_sirna_specificity"),
    },
]

STATIC_FILES = [
    "candidate_offtargets_nearperfect.tsv",
    "transcript_utr_reconstruction.tsv",
    "utr_reconstruction_status.txt",
    "required_input_specification.txt",
    "required_inputs_status.txt",
    "sirna_derivation_notes.txt",
]


def run(cmd, env=None):
    subprocess.run(cmd, check=True, env=env)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def prepare_variant(base_spec_dir, variant_dir):
    ensure_dir(variant_dir)
    for name in STATIC_FILES:
        src = os.path.join(base_spec_dir, name)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(variant_dir, name))


def main():
    for dataset in DATASETS:
        base_spec_dir = dataset["base_spec_dir"]
        analysis_dir = dataset["analysis_dir"]
        for seed_len in (5, 6, 7, 8):
            variant_dir = os.path.join(base_spec_dir, f"seed{seed_len}mer")
            prepare_variant(base_spec_dir, variant_dir)
            run(["python", RECOMPUTE_SCRIPT, variant_dir, str(seed_len)])
            run(["Rscript", R_SCRIPT, ROOT, variant_dir, analysis_dir])


if __name__ == "__main__":
    main()
