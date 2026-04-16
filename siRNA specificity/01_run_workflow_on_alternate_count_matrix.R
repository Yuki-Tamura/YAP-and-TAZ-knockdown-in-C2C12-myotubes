#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript analysis/run_dataset_workflow.R <count_matrix_path> <output_analysis_subdir>")
}

project_dir <- normalizePath(".", mustWork = TRUE)
code_dir <- file.path(project_dir, "analysis")
count_file <- normalizePath(args[[1]], mustWork = TRUE)
analysis_dir <- file.path(code_dir, args[[2]])

source(file.path(code_dir, "R", "helpers.R"))

run_workflow(
  project_dir = project_dir,
  analysis_dir = analysis_dir,
  count_file = count_file,
  code_dir = code_dir
)
