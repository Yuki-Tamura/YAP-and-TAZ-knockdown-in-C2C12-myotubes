#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_dir <- normalizePath(if (length(args) >= 1) args[[1]] else ".", mustWork = TRUE)
analysis_dir <- file.path(project_dir, "analysis")

source(file.path(analysis_dir, "R", "helpers.R"))

run_workflow(project_dir = project_dir, analysis_dir = analysis_dir)
