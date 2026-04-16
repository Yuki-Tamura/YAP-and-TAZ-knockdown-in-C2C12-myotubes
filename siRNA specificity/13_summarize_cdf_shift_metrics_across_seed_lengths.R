#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- normalizePath(args[[1]], mustWork = TRUE)
spec_root <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- file.path(spec_root, "cdf_metric_extensions_seed5to8")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed_dirs <- list.dirs(spec_root, recursive = FALSE, full.names = TRUE)
seed_dirs <- seed_dirs[grepl("seed[5-8]mer$", seed_dirs)]

yap_res <- read.csv(file.path(analysis_dir, "02_deseq2", "siYAP_vs_siScr_all_genes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
taz_res <- read.csv(file.path(analysis_dir, "02_deseq2", "siTAZ_vs_siScr_all_genes.csv"), stringsAsFactors = FALSE, check.names = FALSE)

compute_abc <- function(matched, unmatched) {
  x <- sort(unique(c(matched, unmatched)))
  if (length(x) < 2) {
    return(c(signed = NA_real_, absolute = NA_real_))
  }
  Fm <- ecdf(matched)
  Fu <- ecdf(unmatched)
  diff_y <- Fm(x) - Fu(x)
  dx <- diff(x)
  mid_diff <- (diff_y[-length(diff_y)] + diff_y[-1]) / 2
  signed <- sum(mid_diff * dx)
  absolute <- sum(abs(mid_diff) * dx)
  c(signed = signed, absolute = absolute)
}

compute_tail_metrics <- function(matched, unmatched, threshold) {
  m_hit <- sum(matched < threshold, na.rm = TRUE)
  m_non <- sum(matched >= threshold, na.rm = TRUE)
  u_hit <- sum(unmatched < threshold, na.rm = TRUE)
  u_non <- sum(unmatched >= threshold, na.rm = TRUE)
  ft <- fisher.test(matrix(c(m_hit, m_non, u_hit, u_non), nrow = 2), alternative = "greater")
  tibble(
    threshold = threshold,
    matched_tail_fraction = m_hit / length(matched),
    unmatched_tail_fraction = u_hit / length(unmatched),
    tail_fraction_difference = (m_hit / length(matched)) - (u_hit / length(unmatched)),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

compute_metrics <- function(res_df, seed_df, match_col, comparison_label, seed_dir_name) {
  merged <- res_df %>%
    mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
    inner_join(seed_df, by = "gene_id_stable")

  matched <- merged$log2FoldChange[merged[[match_col]] == 1]
  unmatched <- merged$log2FoldChange[merged[[match_col]] == 0]
  ks <- suppressWarnings(ks.test(matched, unmatched))
  abc <- compute_abc(matched, unmatched)

  summary_row <- tibble(
    seed_variant = seed_dir_name,
    comparison = comparison_label,
    match_column = match_col,
    n_matched = length(matched),
    n_unmatched = length(unmatched),
    median_matched = median(matched, na.rm = TRUE),
    median_unmatched = median(unmatched, na.rm = TRUE),
    delta_median = median(matched, na.rm = TRUE) - median(unmatched, na.rm = TRUE),
    abc_signed = unname(abc["signed"]),
    abc_absolute = unname(abc["absolute"]),
    ks_d = unname(ks$statistic),
    ks_p = ks$p.value
  )

  tail_rows <- bind_rows(
    compute_tail_metrics(matched, unmatched, -1.0),
    compute_tail_metrics(matched, unmatched, -0.5),
    compute_tail_metrics(matched, unmatched, -0.3)
  ) %>%
    mutate(
      seed_variant = seed_dir_name,
      comparison = comparison_label,
      match_column = match_col
    ) %>%
    relocate(seed_variant, comparison, match_column)

  list(summary = summary_row, tail = tail_rows)
}

all_summary <- list()
all_tail <- list()

for (seed_dir in seed_dirs) {
  seed_dir_name <- basename(seed_dir)
  seed_df <- read.delim(file.path(seed_dir, "sirna_seed_annotations_genes.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

  runs <- list(
    list(res = yap_res, col = "YAP_guide_has_match", label = "siYAP vs siScr"),
    list(res = yap_res, col = "YAP_passenger_has_match", label = "siYAP vs siScr"),
    list(res = taz_res, col = "TAZ_guide_has_match", label = "siTAZ vs siScr"),
    list(res = taz_res, col = "TAZ_passenger_has_match", label = "siTAZ vs siScr")
  )

  for (run in runs) {
    metrics <- compute_metrics(run$res, seed_df, run$col, run$label, seed_dir_name)
    all_summary[[length(all_summary) + 1]] <- metrics$summary
    all_tail[[length(all_tail) + 1]] <- metrics$tail
  }
}

summary_df <- bind_rows(all_summary) %>%
  mutate(
    seed_length = as.integer(str_match(seed_variant, "seed([0-9]+)mer")[, 2]),
    strand_role = ifelse(grepl("guide", match_column), "guide", "passenger"),
    sirna = ifelse(grepl("^YAP_", match_column), "YAP", "TAZ"),
    left_shift_supported = delta_median < 0
  ) %>%
  relocate(seed_length, .after = seed_variant)

tail_df <- bind_rows(all_tail) %>%
  mutate(
    seed_length = as.integer(str_match(seed_variant, "seed([0-9]+)mer")[, 2]),
    strand_role = ifelse(grepl("guide", match_column), "guide", "passenger"),
    sirna = ifelse(grepl("^YAP_", match_column), "YAP", "TAZ")
  ) %>%
  relocate(seed_length, .after = seed_variant) %>%
  group_by(threshold) %>%
  mutate(p_value_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.table(summary_df, file.path(out_dir, "cdf_extended_metrics_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tail_df, file.path(out_dir, "cdf_tail_enrichment_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

guide_plot_df <- summary_df %>%
  filter(strand_role == "guide") %>%
  select(seed_length, sirna, delta_median, abc_absolute, ks_d) %>%
  pivot_longer(cols = c(delta_median, abc_absolute, ks_d), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("delta_median", "abc_absolute", "ks_d"),
                    labels = c("Delta median", "ABC (absolute area between CDFs)", "D value"))
  )

guide_plot <- ggplot(guide_plot_df, aes(seed_length, value, color = sirna)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = c(5, 6, 7, 8)) +
  labs(
    title = "CDF-derived shift metrics across 5-mer to 8-mer seed definitions",
    subtitle = "Guide-seed comparisons for siYAP and siTAZ",
    x = "Seed length",
    y = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"))

ggsave(file.path(out_dir, "cdf_metric_guide_summary.png"), guide_plot, width = 8, height = 9, dpi = 300)

tail_plot_df <- tail_df %>%
  filter(strand_role == "guide", threshold == -1.0)

tail_plot <- ggplot(tail_plot_df, aes(seed_length, odds_ratio, color = sirna)) +
  geom_hline(yintercept = 1, color = "grey80", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(5, 6, 7, 8)) +
  labs(
    title = "Tail enrichment for strong downregulation",
    subtitle = "Guide-seed matched vs unmatched genes, threshold log2FC < -1",
    x = "Seed length",
    y = "Fisher odds ratio",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "cdf_tail_enrichment_guide_log2fc_lt_minus1.png"), tail_plot, width = 7, height = 4.5, dpi = 300)

memo_lines <- c(
  "CDF metric extensions memo",
  "",
  "This folder provides additional descriptors for seed-matched vs seed-unmatched fold-change distributions across seed5mer to seed8mer outputs.",
  "",
  "Definitions used:",
  "- Delta median: median(log2FC_seed_matched) minus median(log2FC_seed_unmatched). More negative values indicate stronger typical repression in seed-matched genes.",
  "- ABC: numerical area between the two empirical CDFs. Both signed and absolute values are reported; the absolute value summarizes the overall magnitude of distributional separation.",
  "- Tail enrichment: Fisher exact enrichment of strongly downregulated genes among seed-matched genes, evaluated at log2FC thresholds of -1.0, -0.5, and -0.3.",
  "- D value: the Kolmogorov-Smirnov statistic, reported as a descriptive measure of distributional separation.",
  "",
  "Interpretation guidance:",
  "- Delta median is useful for the typical size of repression.",
  "- ABC is useful for the overall shift across the full distribution.",
  "- Tail enrichment is useful for asking whether a small subset of strongly repressed candidate off-targets is present.",
  "- D value is a complementary descriptive statistic and should not be over-interpreted in isolation."
)

writeLines(memo_lines, file.path(out_dir, "cdf_metric_extensions_memo.txt"))
