#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- normalizePath(args[[1]], mustWork = TRUE)
spec_root <- normalizePath(args[[2]], mustWork = TRUE)

seed_dirs <- list.dirs(spec_root, recursive = FALSE, full.names = TRUE)
seed_dirs <- seed_dirs[grepl("seed[5-8]mer$", seed_dirs)]

double_res <- read.csv(
  file.path(analysis_dir, "02_deseq2", "siYAPsiTAZ_vs_siScr_all_genes.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

build_category_table <- function(seed_genes) {
  seed_genes %>%
    transmute(
      gene_id = gene_id,
      gene_id_stable = gene_id_stable,
      gene_symbol = ifelse(is.na(gene_name) | gene_name == "", gene_id_stable, gene_name),
      max_utr_length = max_utr_length,
      n_utr_transcripts = n_utr_transcripts,
      YAP_guide = YAP_guide_has_match == 1,
      TAZ_guide = TAZ_guide_has_match == 1,
      YAP_passenger = YAP_passenger_has_match == 1,
      TAZ_passenger = TAZ_passenger_has_match == 1
    ) %>%
    mutate(
      any_guide = YAP_guide | TAZ_guide,
      both_guides = YAP_guide & TAZ_guide,
      any_passenger = YAP_passenger | TAZ_passenger,
      both_passengers = YAP_passenger & TAZ_passenger,
      any_of_four = any_guide | any_passenger
    )
}

compute_cdf_stats <- function(df, match_col, label) {
  matched <- df$log2FoldChange[df[[match_col]]]
  unmatched <- df$log2FoldChange[!df[[match_col]]]
  ks_res <- suppressWarnings(ks.test(matched, unmatched))
  wilcox_res <- suppressWarnings(wilcox.test(matched, unmatched, exact = FALSE))
  tibble(
    category = label,
    match_column = match_col,
    n_matched = sum(df[[match_col]], na.rm = TRUE),
    n_unmatched = sum(!df[[match_col]], na.rm = TRUE),
    median_matched = median(matched, na.rm = TRUE),
    median_unmatched = median(unmatched, na.rm = TRUE),
    delta_median = median(matched, na.rm = TRUE) - median(unmatched, na.rm = TRUE),
    ks_d = unname(ks_res$statistic),
    ks_p = ks_res$p.value,
    wilcox_p = wilcox_res$p.value,
    effect_direction = ifelse(
      median(matched, na.rm = TRUE) < median(unmatched, na.rm = TRUE),
      "Left-shifted in matched genes",
      "No left shift in matched genes"
    )
  )
}

make_cdf_plot <- function(plot_df, title_text) {
  ggplot(plot_df, aes(x = log2FoldChange, color = match_group, linetype = match_group)) +
    stat_ecdf(linewidth = 0.9) +
    scale_color_manual(values = c("Matched" = "#B23A48", "Unmatched" = "grey55")) +
    scale_linetype_manual(values = c("Matched" = "solid", "Unmatched" = "dashed")) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.4) +
    labs(
      title = title_text,
      x = expression(log[2]~fold~change),
      y = "Cumulative fraction of genes",
      color = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
}

category_specs <- tribble(
  ~match_column, ~category_label, ~file_stub,
  "YAP_guide", "YAP guide seed-matched genes", "yap_guide",
  "TAZ_guide", "TAZ guide seed-matched genes", "taz_guide",
  "YAP_passenger", "YAP passenger seed-matched genes", "yap_passenger",
  "TAZ_passenger", "TAZ passenger seed-matched genes", "taz_passenger",
  "any_guide", "Any guide-seed-matched genes", "any_guide",
  "both_guides", "Both guide-seed-matched genes", "both_guides",
  "any_passenger", "Any passenger-seed-matched genes", "any_passenger",
  "both_passengers", "Both passenger-seed-matched genes", "both_passengers",
  "any_of_four", "Any guide or passenger seed-matched genes", "any_of_four"
)

for (seed_dir in seed_dirs) {
  out_dir <- file.path(seed_dir, "double_knockdown_cdf_two_guides_two_passengers")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  seed_genes <- read.delim(
    file.path(seed_dir, "sirna_seed_annotations_genes.tsv"),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  category_df <- build_category_table(seed_genes)
  merged <- double_res %>%
    mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
    inner_join(category_df, by = "gene_id_stable")

  write.table(
    merged,
    file.path(out_dir, "double_knockdown_seed_gene_annotations.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  stats_df <- bind_rows(lapply(seq_len(nrow(category_specs)), function(i) {
    spec <- category_specs[i, ]
    compute_cdf_stats(merged, spec$match_column, spec$category_label)
  })) %>%
    mutate(
      ks_p_bh = p.adjust(ks_p, method = "BH"),
      wilcox_p_bh = p.adjust(wilcox_p, method = "BH")
    )

  write.table(
    stats_df,
    file.path(out_dir, "double_knockdown_seed_cdf_stats.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  plot_source <- bind_rows(lapply(seq_len(nrow(category_specs)), function(i) {
    spec <- category_specs[i, ]
    merged %>%
      transmute(
        category = spec$category_label,
        file_stub = spec$file_stub,
        log2FoldChange = log2FoldChange,
        matched = .data[[spec$match_column]]
      )
  })) %>%
    mutate(
      match_group = ifelse(matched, "Matched", "Unmatched"),
      category = factor(category, levels = category_specs$category_label)
    )

  write.table(
    plot_source,
    file.path(out_dir, "double_knockdown_seed_cdf_source_data.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  panel_plot <- ggplot(plot_source, aes(x = log2FoldChange, color = match_group, linetype = match_group)) +
    stat_ecdf(linewidth = 0.9) +
    facet_wrap(~ category, ncol = 2) +
    scale_color_manual(values = c("Matched" = "#B23A48", "Unmatched" = "grey55")) +
    scale_linetype_manual(values = c("Matched" = "solid", "Unmatched" = "dashed")) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.4) +
    labs(
      title = "Double-knockdown CDFs for guide and passenger seed-match categories",
      subtitle = "siYAP+siTAZ vs siScr",
      x = expression(log[2]~fold~change),
      y = "Cumulative fraction of genes",
      color = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))

  ggsave(file.path(out_dir, "double_knockdown_seed_cdf_panel.png"), panel_plot, width = 11, height = 10, dpi = 300)
  ggsave(file.path(out_dir, "double_knockdown_seed_cdf_panel.pdf"), panel_plot, width = 11, height = 10)

  single_dir <- file.path(out_dir, "single_panels")
  dir.create(single_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(category_specs))) {
    spec <- category_specs[i, ]
    panel_df <- plot_source %>% filter(file_stub == spec$file_stub)
    plot_obj <- make_cdf_plot(panel_df, paste0(spec$category_label, " in double knockdown"))
    ggsave(file.path(single_dir, paste0(spec$file_stub, "_double_knockdown_cdf.png")), plot_obj, width = 6, height = 5, dpi = 300)
    ggsave(file.path(single_dir, paste0(spec$file_stub, "_double_knockdown_cdf.pdf")), plot_obj, width = 6, height = 5)
  }

  memo_lines <- c(
    "Double-knockdown seed CDF memo",
    "",
    "This folder summarizes CDF analyses for siYAP+siTAZ vs siScr using the seed annotations already generated for the current seed-length variant.",
    "Assumption used: because the double knockdown contains two guide candidates and two passenger candidates, both individual-seed categories and combined categories were reported.",
    "",
    "Included categories:",
    "- YAP guide",
    "- TAZ guide",
    "- YAP passenger",
    "- TAZ passenger",
    "- any guide seed",
    "- both guide seeds",
    "- any passenger seed",
    "- both passenger seeds",
    "- any guide or passenger seed",
    "",
    "Interpretation guidance:",
    "- Individual categories show whether one specific seed is associated with the double-knockdown distribution.",
    "- Combined guide categories are most useful for asking whether the double-knockdown shift is concentrated in genes carrying one or both guide-seed annotations.",
    "- These results are supportive only and do not prove that the double-knockdown response is explained by seed-mediated off-target effects."
  )

  writeLines(memo_lines, file.path(out_dir, "double_knockdown_seed_cdf_memo.txt"))
}
