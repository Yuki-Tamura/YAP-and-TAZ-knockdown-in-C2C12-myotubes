#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- normalizePath(if (length(args) >= 1) args[[1]] else ".", mustWork = TRUE)
base_dir <- file.path(project_dir, "analysis", "04_sirna_specificity")
out_dir <- file.path(base_dir, "null_seed_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv_safe <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
read_tsv_safe <- function(path) read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

strand_df <- read_csv_safe(file.path(base_dir, "sirna_strand_derivation.csv"))
utr_recon <- read_tsv_safe(file.path(base_dir, "transcript_utr_reconstruction.tsv"))
yap_res <- read_csv_safe(file.path(project_dir, "analysis", "02_deseq2", "siYAP_vs_siScr_all_genes.csv"))
taz_res <- read_csv_safe(file.path(project_dir, "analysis", "02_deseq2", "siTAZ_vs_siScr_all_genes.csv"))
yap_go <- read_csv_safe(file.path(project_dir, "analysis", "03_pathways", "siYAP_vs_siScr_go_bp_focused.csv"))

extract_go_genes <- function(df, pattern, top_n = 3) {
  sub <- df %>%
    filter(grepl(pattern, TERM, ignore.case = TRUE)) %>%
    arrange(padj, pvalue)
  if (nrow(sub) == 0) return(character())
  sub <- slice_head(sub, n = min(top_n, nrow(sub)))
  unique(unlist(strsplit(paste(sub$overlap_genes, collapse = ";"), ";", fixed = TRUE)))
}

reverse_complement_rna <- function(seq) {
  chars <- strsplit(seq, "")[[1]]
  comp <- c(A = "U", U = "A", G = "C", C = "G")
  paste(rev(unname(comp[chars])), collapse = "")
}

compute_seed_stats <- function(seed, df, lfc_col) {
  motif <- reverse_complement_rna(seed)
  matched <- grepl(motif, df$utr_seq, fixed = TRUE)
  lfc <- df[[lfc_col]]
  matched_vals <- lfc[matched]
  unmatched_vals <- lfc[!matched]
  if (length(matched_vals) < 10 || length(unmatched_vals) < 10) {
    return(tibble(
      n_seed_matched_genes = length(matched_vals),
      n_seed_unmatched_genes = length(unmatched_vals),
      median_log2fc_seed_matched = median(matched_vals, na.rm = TRUE),
      median_log2fc_seed_unmatched = median(unmatched_vals, na.rm = TRUE),
      median_difference = median(matched_vals, na.rm = TRUE) - median(unmatched_vals, na.rm = TRUE),
      ks_statistic = NA_real_,
      ks_p_value = NA_real_,
      wilcox_statistic = NA_real_,
      wilcox_p_value = NA_real_
    ))
  }
  ks_res <- suppressWarnings(ks.test(matched_vals, unmatched_vals))
  wilcox_res <- suppressWarnings(wilcox.test(matched_vals, unmatched_vals, exact = FALSE))
  tibble(
    n_seed_matched_genes = length(matched_vals),
    n_seed_unmatched_genes = length(unmatched_vals),
    median_log2fc_seed_matched = median(matched_vals, na.rm = TRUE),
    median_log2fc_seed_unmatched = median(unmatched_vals, na.rm = TRUE),
    median_difference = median(matched_vals, na.rm = TRUE) - median(unmatched_vals, na.rm = TRUE),
    ks_statistic = unname(ks_res$statistic),
    ks_p_value = ks_res$p.value,
    wilcox_statistic = unname(wilcox_res$statistic),
    wilcox_p_value = wilcox_res$p.value
  )
}

generate_null_seeds <- function(observed_seed, n = 1000, exclude = character()) {
  set.seed(20260414)
  nucs <- c("A", "U", "G", "C")
  obs_chars <- strsplit(observed_seed, "")[[1]]
  len_target <- nchar(observed_seed)
  gc_target <- sum(obs_chars %in% c("G", "C"))
  obs_comp <- table(factor(obs_chars, levels = nucs))
  candidates <- character()
  while (length(candidates) < n * 10) {
    cand <- paste0(sample(nucs, len_target, replace = TRUE, prob = as.numeric(obs_comp)), collapse = "")
    if (cand %in% c(observed_seed, exclude)) next
    if (sum(strsplit(cand, "")[[1]] %in% c("G", "C")) != gc_target) next
    candidates <- c(candidates, cand)
  }
  candidates <- unique(candidates)
  comp_distance <- vapply(candidates, function(s) {
    cnt <- table(factor(strsplit(s, "")[[1]], levels = nucs))
    sum(abs(cnt - obs_comp))
  }, numeric(1))
  tibble(seed = candidates, composition_distance = comp_distance) %>%
    arrange(composition_distance, seed) %>%
    slice_head(n = min(n, length(candidates))) %>%
    pull(seed)
}

prepare_gene_level_df <- function(yap_res, taz_res, utr_recon) {
  representative_utr <- utr_recon %>%
    filter(utr_reconstructable == 1, !is.na(utr_seq), utr_seq != "") %>%
    mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
    arrange(gene_id_stable, desc(computed_utr_len), transcript_id) %>%
    group_by(gene_id_stable) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      gene_id_stable,
      representative_transcript_id = transcript_id,
      utr_seq = utr_seq,
      utr_source = utr_source,
      representative_utr_length = computed_utr_len
    )

  yap_res %>%
    transmute(gene_id_stable = sub("\\..*$", "", gene_id), gene_id = gene_id, yap_log2FC = log2FoldChange, yap_padj = padj, baseMean = baseMean) %>%
    inner_join(
      taz_res %>% transmute(gene_id_stable = sub("\\..*$", "", gene_id), taz_log2FC = log2FoldChange, taz_padj = padj),
      by = "gene_id_stable"
    ) %>%
    inner_join(representative_utr, by = "gene_id_stable")
}

gene_df <- prepare_gene_level_df(yap_res, taz_res, utr_recon)

yap_seed_6mer <- unique(strand_df$guide_seed_2_7[strand_df$sirna == "YAP"])[1]
taz_seed_6mer <- unique(strand_df$guide_seed_2_7[strand_df$sirna == "TAZ"])[1]
yap_seed_7mer <- substr(unique(strand_df$guide_candidate[strand_df$sirna == "YAP"])[1], 2, 8)
taz_seed_7mer <- substr(unique(strand_df$guide_candidate[strand_df$sirna == "TAZ"])[1], 2, 8)

null_seed_n <- 1000
yap_null_seeds <- generate_null_seeds(yap_seed_6mer, n = null_seed_n, exclude = c(taz_seed_6mer))
taz_null_seeds <- generate_null_seeds(taz_seed_6mer, n = 500, exclude = c(yap_seed_6mer))

run_seed_distribution <- function(observed_seed, null_seeds, df, lfc_col, analysis_name) {
  observed <- compute_seed_stats(observed_seed, df, lfc_col) %>%
    mutate(seed = observed_seed, seed_type = "observed", analysis = analysis_name)
  nulls <- bind_rows(lapply(null_seeds, function(s) {
    compute_seed_stats(s, df, lfc_col) %>%
      mutate(seed = s, seed_type = "null", analysis = analysis_name)
  }))
  bind_rows(observed, nulls) %>%
    dplyr::select(analysis, seed_type, seed, everything())
}

yguide_dist <- run_seed_distribution(yap_seed_6mer, yap_null_seeds, gene_df, "yap_log2FC", "YAP_guide_6mer")
tguide_dist <- run_seed_distribution(taz_seed_6mer, taz_null_seeds, gene_df, "taz_log2FC", "TAZ_guide_6mer")

write.table(yguide_dist, file.path(out_dir, "yguide_null_seed_distribution.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tguide_dist, file.path(out_dir, "tguide_null_seed_distribution.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summarize_observed_vs_null <- function(dist_df, label) {
  obs <- dist_df %>% filter(seed_type == "observed")
  nulls <- dist_df %>% filter(seed_type == "null")
  tibble(
    analysis = label,
    observed_seed = obs$seed,
    n_null_seeds = nrow(nulls),
    observed_n_seed_matched = obs$n_seed_matched_genes,
    observed_median_difference = obs$median_difference,
    observed_ks_statistic = obs$ks_statistic,
    observed_ks_p_value = obs$ks_p_value,
    observed_wilcox_statistic = obs$wilcox_statistic,
    observed_wilcox_p_value = obs$wilcox_p_value,
    empirical_percentile_median_difference = mean(nulls$median_difference <= obs$median_difference, na.rm = TRUE),
    empirical_p_left_tail_median_difference = mean(nulls$median_difference <= obs$median_difference, na.rm = TRUE),
    empirical_percentile_ks_statistic = mean(nulls$ks_statistic <= obs$ks_statistic, na.rm = TRUE),
    empirical_p_right_tail_ks_statistic = mean(nulls$ks_statistic >= obs$ks_statistic, na.rm = TRUE)
  )
}

observed_vs_null <- bind_rows(
  summarize_observed_vs_null(yguide_dist, "YAP_guide_6mer"),
  summarize_observed_vs_null(tguide_dist, "TAZ_guide_6mer")
)
write.table(observed_vs_null, file.path(out_dir, "observed_vs_null_seed_statistics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_null_distribution <- function(dist_df, title_text, out_png, out_pdf = NULL) {
  obs <- dist_df %>% filter(seed_type == "observed")
  nulls <- dist_df %>% filter(seed_type == "null")
  p1 <- ggplot(nulls, aes(x = median_difference)) +
    geom_histogram(bins = 50, fill = "grey75", color = "white") +
    geom_vline(xintercept = obs$median_difference, color = "#B23A48", linewidth = 1) +
    labs(
      title = title_text,
      subtitle = sprintf("Observed median difference = %.3f", obs$median_difference),
      x = "Median log2FC difference (seed-matched minus unmatched)",
      y = "Number of null seeds"
    ) +
    theme_bw(base_size = 11)
  p2 <- ggplot(nulls, aes(x = ks_statistic)) +
    geom_histogram(bins = 50, fill = "grey75", color = "white") +
    geom_vline(xintercept = obs$ks_statistic, color = "#2F5D8A", linewidth = 1) +
    labs(
      subtitle = sprintf("Observed KS statistic = %.3f", obs$ks_statistic),
      x = "KS statistic",
      y = "Number of null seeds"
    ) +
    theme_bw(base_size = 11)
  png(file = out_png, width = 1800, height = 1400, res = 200)
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  print(p1)
  print(p2)
  dev.off()
  if (!is.null(out_pdf)) {
    pdf(out_pdf, width = 8, height = 9)
    print(p1)
    print(p2)
    dev.off()
  }
}

plot_null_distribution(
  yguide_dist,
  "Observed YAP-guide seed effect compared with matched null-seed controls",
  file.path(out_dir, "yguide_null_seed_plot.png"),
  file.path(out_dir, "yguide_null_seed_plot.pdf")
)
plot_null_distribution(
  tguide_dist,
  "Observed TAZ-guide seed effect compared with matched null-seed controls",
  file.path(out_dir, "tguide_null_seed_plot.png"),
  NULL
)

mask_list <- list(
  baseline_none = character(),
  yap1_only = unique(gene_df$gene_id_stable[gene_df$gene_symbol %in% c("Yap1", "YAP1")]),
  top50_yap_down = head(yap_res %>% mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>% filter(!is.na(padj), log2FoldChange < 0) %>% arrange(padj, log2FoldChange) %>% pull(gene_id_stable), 50),
  top100_yap_down = head(yap_res %>% mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>% filter(!is.na(padj), log2FoldChange < 0) %>% arrange(padj, log2FoldChange) %>% pull(gene_id_stable), 100),
  muscle_terms = extract_go_genes(yap_go, "muscle contraction|muscle system process|myofibril|sarcomere", top_n = 4)
)
if (length(mask_list$yap1_only) == 0) {
  mask_list$yap1_only <- unique(gene_df$gene_id_stable[gene_df$gene_symbol %in% c("Yap1", "YAP1")])
}

masked_summary <- bind_rows(lapply(names(mask_list), function(mask_name) {
  sub_df <- gene_df %>% filter(!(gene_id_stable %in% mask_list[[mask_name]]))
  dist <- run_seed_distribution(yap_seed_6mer, yap_null_seeds, sub_df, "yap_log2FC", paste0("YAP_guide_6mer_", mask_name))
  obs <- dist %>% filter(seed_type == "observed")
  nulls <- dist %>% filter(seed_type == "null")
  tibble(
    mask_name = mask_name,
    n_excluded = nrow(gene_df) - nrow(sub_df),
    observed_n_seed_matched = obs$n_seed_matched_genes,
    observed_median_difference = obs$median_difference,
    observed_ks_statistic = obs$ks_statistic,
    observed_wilcox_p_value = obs$wilcox_p_value,
    empirical_percentile_median_difference = mean(nulls$median_difference <= obs$median_difference, na.rm = TRUE),
    empirical_p_left_tail_median_difference = mean(nulls$median_difference <= obs$median_difference, na.rm = TRUE)
  )
}))
write.table(masked_summary, file.path(out_dir, "yguide_null_seed_after_masking_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

masked_plot <- ggplot(masked_summary, aes(x = factor(mask_name, levels = mask_name), y = observed_median_difference)) +
  geom_col(fill = "#B23A48") +
  geom_text(aes(label = sprintf("empirical p=%.3f", empirical_p_left_tail_median_difference)), vjust = -0.5, size = 3) +
  labs(
    title = "Sensitivity of YAP-guide seed-associated repression to exclusion of likely on-target YAP-responsive genes",
    x = NULL,
    y = "Observed median difference after masking"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "yguide_null_seed_masked_plot.png"), masked_plot, width = 9, height = 5, dpi = 300)

methods_lines <- c(
  "Null-seed sensitivity methods notes",
  "",
  sprintf("Observed YAP guide seed used for null comparison: 6-mer positions 2-7 = %s.", yap_seed_6mer),
  sprintf("Observed TAZ guide seed used for comparator null comparison: 6-mer positions 2-7 = %s.", taz_seed_6mer),
  sprintf("For reference, the corresponding 7-mer seeds (positions 2-8) are YAP %s and TAZ %s.", yap_seed_7mer, taz_seed_7mer),
  "The existing CDF framework in this project used 6-mer seed matching, so the null analysis was matched to that same 6-mer definition for like-for-like comparison.",
  sprintf("Null seeds were generated with matched seed length and GC content, sampled with probabilities proportional to the observed seed nucleotide composition, then ranked by composition distance. %d YAP null seeds and %d TAZ null seeds were retained after excluding duplicates, the true observed seed, and the other guide seed.", length(yap_null_seeds), length(taz_null_seeds)),
  "Gene-level analysis used one representative reconstructed 3'UTR per gene: the longest reconstructed 3'UTR available for that gene.",
  "Seed-matched genes were defined as genes whose representative 3'UTR contained at least one reverse-complement seed motif.",
  "Statistics were computed as in the observed analysis: number of matched genes, matched/unmatched median log2FC, median difference, KS statistic and P value, and Wilcoxon statistic and P value.",
  "Masked sensitivity summaries reused likely on-target YAP masks derived from prior project analyses."
)
writeLines(methods_lines, file.path(out_dir, "null_seed_methods_notes.txt"))

y_obs <- observed_vs_null %>% filter(analysis == "YAP_guide_6mer")
t_obs <- observed_vs_null %>% filter(analysis == "TAZ_guide_6mer")
memo_lines <- c(
  "YAP-guide null-seed interpretation memo",
  "",
  sprintf("Observed YAP-guide seed: %s (canonical 6-mer, positions 2-7 of the inferred guide strand).", yap_seed_6mer),
  sprintf("Observed YAP-guide median difference = %.3f, KS statistic = %.3f, Wilcoxon P = %s.", y_obs$observed_median_difference, y_obs$observed_ks_statistic, formatC(y_obs$observed_wilcox_p_value, format = "e", digits = 2)),
  sprintf("Compared with %d matched null seeds, the observed YAP-guide effect lay at the %.1fth percentile for left-shift magnitude and had empirical left-tail P = %.3f.", y_obs$n_null_seeds, 100 * y_obs$empirical_percentile_median_difference, y_obs$empirical_p_left_tail_median_difference),
  sprintf("Observed TAZ-guide median difference = %.3f, with empirical left-tail P = %.3f against %d matched null seeds.", t_obs$observed_median_difference, t_obs$empirical_p_left_tail_median_difference, t_obs$n_null_seeds),
  "",
  sprintf("After masking likely on-target YAP-responsive genes, the observed YAP-guide median difference ranged from %.3f to %.3f across masks.", min(masked_summary$observed_median_difference), max(masked_summary$observed_median_difference)),
  sprintf("The empirical left-tail null P values after masking ranged from %.3f to %.3f.", min(masked_summary$empirical_p_left_tail_median_difference), max(masked_summary$empirical_p_left_tail_median_difference)),
  "",
  "Conservative interpretation:",
  if (y_obs$empirical_p_left_tail_median_difference <= 0.01) {
    "- The observed YAP-guide seed effect appears more extreme than would be expected from matched background 3'UTR sequence composition alone."
  } else {
    "- The observed YAP-guide seed effect falls within the range seen for many matched null seeds, so it may reflect background sequence properties and should not be overinterpreted."
  },
  if (max(masked_summary$empirical_p_left_tail_median_difference) > y_obs$empirical_p_left_tail_median_difference + 0.05) {
    "- Excluding likely on-target YAP-responsive genes weakens the apparent YAP-guide seed effect to some degree, so true YAP biology may inflate part of the original shift."
  } else {
    "- Excluding likely on-target YAP-responsive genes does not substantially weaken the observed-vs-null separation, so the signal does not appear to be mainly driven by those masks."
  },
  "- This analysis supports or weakens the case for meaningful seed-associated repression only relative to matched null seeds; it does not fully separate on-target and off-target biology or exclude all off-target effects."
)
writeLines(memo_lines, file.path(out_dir, "yguide_null_seed_interpretation_memo.txt"))
