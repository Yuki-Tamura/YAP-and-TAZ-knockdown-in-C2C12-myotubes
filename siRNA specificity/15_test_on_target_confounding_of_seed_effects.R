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
out_dir <- file.path(base_dir, "sensitivity_on_target_confounding")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_read_csv <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
safe_read_tsv <- function(path) read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

yap_res <- safe_read_csv(file.path(project_dir, "analysis", "02_deseq2", "siYAP_vs_siScr_all_genes.csv"))
taz_res <- safe_read_csv(file.path(project_dir, "analysis", "02_deseq2", "siTAZ_vs_siScr_all_genes.csv"))
seed_genes <- safe_read_tsv(file.path(base_dir, "sirna_seed_annotations_genes.tsv"))
utr_recon <- safe_read_tsv(file.path(base_dir, "transcript_utr_reconstruction.tsv"))
sirna_seed_tx <- safe_read_tsv(file.path(base_dir, "sirna_seed_annotations.tsv"))
yap_go <- safe_read_csv(file.path(project_dir, "analysis", "03_pathways", "siYAP_vs_siScr_go_bp_focused.csv"))

compute_gc <- function(seq) {
  if (is.na(seq) || seq == "") return(NA_real_)
  chars <- strsplit(seq, "")[[1]]
  mean(chars %in% c("G", "C"))
}

extract_go_genes <- function(df, pattern, top_n = 3) {
  sub <- df %>%
    filter(grepl(pattern, TERM, ignore.case = TRUE)) %>%
    arrange(padj, pvalue)
  if (nrow(sub) == 0) return(character())
  sub <- slice_head(sub, n = min(top_n, nrow(sub)))
  unique(unlist(strsplit(paste(sub$overlap_genes, collapse = ";"), ";", fixed = TRUE)))
}

gene_covariates <- utr_recon %>%
  filter(utr_reconstructable == 1) %>%
  mutate(
    gene_id_stable = sub("\\..*$", "", gene_id),
    utr_gc_content = vapply(utr_seq, compute_gc, numeric(1))
  ) %>%
  group_by(gene_id_stable) %>%
  arrange(desc(computed_utr_len), transcript_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    gene_id_stable,
    representative_transcript_id = transcript_id,
    representative_utr_length = computed_utr_len,
    representative_utr_gc_content = utr_gc_content,
    representative_utr_source = utr_source
  )

seed_tx_longest <- sirna_seed_tx %>%
  mutate(gene_id_stable = sub("\\..*$", "", gene_id_stable)) %>%
  inner_join(gene_covariates %>% dplyr::select(gene_id_stable, representative_transcript_id), by = "gene_id_stable") %>%
  filter(transcript_id == representative_transcript_id) %>%
  transmute(
    gene_id_stable,
    yguide_match_longest = YAP_guide_has_match == 1,
    tguide_match_longest = TAZ_guide_has_match == 1
  )

base_df <- yap_res %>%
  transmute(
    gene_id = gene_id,
    gene_id_stable = sub("\\..*$", "", gene_id),
    baseMean = baseMean,
    yap_log2FC = log2FoldChange,
    yap_padj = padj
  ) %>%
  inner_join(
    taz_res %>%
      transmute(
        gene_id_stable = sub("\\..*$", "", gene_id),
        taz_log2FC = log2FoldChange,
        taz_padj = padj
      ),
    by = "gene_id_stable"
  ) %>%
  left_join(seed_genes, by = "gene_id_stable") %>%
  left_join(gene_covariates, by = "gene_id_stable") %>%
  left_join(seed_tx_longest, by = "gene_id_stable") %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_name) | gene_name == "", gene_id_stable, gene_name),
    yguide_seed_match = ifelse(is.na(YAP_guide_has_match), FALSE, YAP_guide_has_match == 1),
    tguide_seed_match = ifelse(is.na(TAZ_guide_has_match), FALSE, TAZ_guide_has_match == 1),
    ypassenger_seed_match = ifelse(is.na(YAP_passenger_has_match), FALSE, YAP_passenger_has_match == 1),
    tpassenger_seed_match = ifelse(is.na(TAZ_passenger_has_match), FALSE, TAZ_passenger_has_match == 1),
    log_baseMean = log10(baseMean + 1),
    log_utr_length = log10(representative_utr_length + 1)
  )

yap_top_down <- yap_res %>%
  mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
  filter(!is.na(padj), log2FoldChange < 0) %>%
  arrange(padj, log2FoldChange)

taz_top_down <- taz_res %>%
  mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
  filter(!is.na(padj), log2FoldChange < 0) %>%
  arrange(padj, log2FoldChange)

mask_list <- list(
  A_intended_target_only = list(
    description = "Mask A: intended target only (Yap1).",
    genes = unique(base_df$gene_id_stable[base_df$gene_symbol == "Yap1"])
  ),
  B_top50_down = list(
    description = "Mask B: top 50 siYAP downregulated genes ranked by adjusted P value then more negative log2FC.",
    genes = head(yap_top_down$gene_id_stable, 50)
  ),
  C_top100_down = list(
    description = "Mask C: top 100 siYAP downregulated genes ranked by adjusted P value then more negative log2FC.",
    genes = head(yap_top_down$gene_id_stable, 100)
  ),
  D_top200_down = list(
    description = "Mask D: top 200 siYAP downregulated genes ranked by adjusted P value then more negative log2FC.",
    genes = head(yap_top_down$gene_id_stable, 200)
  ),
  E_muscle_contractile_go = list(
    description = "Mask E: genes from the top 3 siYAP-enriched muscle/contractile GO terms matching muscle contraction, striated/skeletal contraction, myofibril, or sarcomere wording.",
    genes = extract_go_genes(yap_go, "muscle contraction|striated muscle contraction|skeletal muscle contraction|myofibril|sarcomere", top_n = 3)
  ),
  F_calcium_muscle_system_go = list(
    description = "Mask F: genes from the top 3 siYAP-enriched calcium-handling or muscle-system GO terms matching calcium or muscle system wording.",
    genes = extract_go_genes(yap_go, "calcium|muscle system process", top_n = 3)
  )
)
mask_list$G_combined_union <- list(
  description = "Mask G: union of mask B, mask E, and mask F.",
  genes = unique(c(mask_list$B_top50_down$genes, mask_list$E_muscle_contractile_go$genes, mask_list$F_calcium_muscle_system_go$genes))
)

mask_metadata <- tibble(
  mask_name = names(mask_list),
  description = vapply(mask_list, `[[`, character(1), "description"),
  n_mask_genes = vapply(mask_list, function(x) length(unique(x$genes)), integer(1))
)

run_masked_cdf <- function(df, seed_col, mask_name, mask_genes, label) {
  sub <- df %>%
    filter(!(gene_id_stable %in% mask_genes)) %>%
    filter(!is.na(yap_log2FC))
  matched <- sub$yap_log2FC[sub[[seed_col]]]
  unmatched <- sub$yap_log2FC[!sub[[seed_col]]]
  ks_res <- suppressWarnings(ks.test(matched, unmatched))
  wilcox_res <- suppressWarnings(wilcox.test(matched, unmatched, exact = FALSE))
  tibble(
    mask_name = mask_name,
    analysis_label = label,
    n_excluded = sum(df$gene_id_stable %in% mask_genes),
    n_seed_matched = length(matched),
    n_seed_unmatched = length(unmatched),
    median_seed_matched = median(matched, na.rm = TRUE),
    median_seed_unmatched = median(unmatched, na.rm = TRUE),
    median_difference = median(matched, na.rm = TRUE) - median(unmatched, na.rm = TRUE),
    ks_statistic = unname(ks_res$statistic),
    ks_p = ks_res$p.value,
    wilcox_p = wilcox_res$p.value
  )
}

yguide_masked <- bind_rows(
  tibble(mask_name = "baseline_none", description = "Baseline: no likely on-target mask excluded.", n_mask_genes = 0),
  mask_metadata
) %>%
  rowwise() %>%
  do({
    mask_name <- .$mask_name
    mask_genes <- if (mask_name == "baseline_none") character() else mask_list[[mask_name]]$genes
    stats <- run_masked_cdf(base_df, "yguide_seed_match", mask_name, mask_genes, "YAP guide seed")
    cbind(., stats[setdiff(names(stats), "mask_name")])
  }) %>%
  ungroup() %>%
  mutate(
    ks_p_bh = p.adjust(ks_p, method = "BH"),
    wilcox_p_bh = p.adjust(wilcox_p, method = "BH")
  )

write.table(yguide_masked, file.path(out_dir, "yguide_masked_cdf_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

yguide_plot_df <- bind_rows(lapply(c("baseline_none", names(mask_list)), function(mask_name) {
  mask_genes <- if (mask_name == "baseline_none") character() else mask_list[[mask_name]]$genes
  label <- if (mask_name == "baseline_none") "Baseline (no exclusion)" else mask_name
  base_df %>%
    filter(!(gene_id_stable %in% mask_genes)) %>%
    transmute(
      mask_name = mask_name,
      panel_label = label,
      log2FC = yap_log2FC,
      match_group = ifelse(yguide_seed_match, "YAP guide seed-matched", "Seed-unmatched")
    )
}))

yguide_plot_labels <- yguide_masked %>%
  mutate(
    panel_label = ifelse(mask_name == "baseline_none", "Baseline (no exclusion)", mask_name),
    label = sprintf("n=%d vs %d\nmedian diff=%.3f\nKS FDR=%s\nWilcoxon FDR=%s",
                    n_seed_matched, n_seed_unmatched, median_difference,
                    formatC(ks_p_bh, format = "e", digits = 2),
                    formatC(wilcox_p_bh, format = "e", digits = 2)),
    x = -4.8,
    y = 0.92
  )

yguide_cdf_plot <- ggplot(yguide_plot_df, aes(x = log2FC, color = match_group, linetype = match_group)) +
  stat_ecdf(linewidth = 0.85) +
  facet_wrap(~ panel_label, ncol = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_text(data = yguide_plot_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3) +
  scale_color_manual(values = c("YAP guide seed-matched" = "#B23A48", "Seed-unmatched" = "grey55")) +
  scale_linetype_manual(values = c("YAP guide seed-matched" = "solid", "Seed-unmatched" = "dashed")) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
  labs(
    title = "Sensitivity analysis of YAP-guide seed effects after excluding likely on-target YAP-responsive genes",
    x = "siYAP vs siScr log2 fold change",
    y = "Cumulative fraction of genes",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "yguide_masked_cdf_plot.png"), yguide_cdf_plot, width = 11, height = 12, dpi = 300)

combined_mask_genes <- mask_list$G_combined_union$genes
match_pool <- base_df %>%
  filter(!(gene_id_stable %in% combined_mask_genes)) %>%
  filter(!is.na(yap_log2FC), !is.na(log_baseMean), !is.na(log_utr_length), !is.na(representative_utr_gc_content))

seed_matched_pool <- match_pool %>% filter(yguide_seed_match)
seed_unmatched_pool <- match_pool %>% filter(!yguide_seed_match)

matched_pairs <- bind_rows(lapply(seq_len(nrow(seed_matched_pool)), function(i) {
  target <- seed_matched_pool[i, ]
  candidates <- seed_unmatched_pool %>%
    mutate(
      distance = (log_baseMean - target$log_baseMean)^2 +
        (log_utr_length - target$log_utr_length)^2 +
        (representative_utr_gc_content - target$representative_utr_gc_content)^2
    ) %>%
    arrange(distance) %>%
    slice_head(n = 1)
  tibble(
    seed_gene_id_stable = target$gene_id_stable,
    seed_gene_symbol = target$gene_symbol,
    seed_gene_log2FC = target$yap_log2FC,
    control_gene_id_stable = candidates$gene_id_stable,
    control_gene_symbol = candidates$gene_symbol,
    control_gene_log2FC = candidates$yap_log2FC,
    matching_distance = candidates$distance
  )
}))

matched_diff <- matched_pairs$seed_gene_log2FC - matched_pairs$control_gene_log2FC
matched_wilcox <- suppressWarnings(wilcox.test(matched_pairs$seed_gene_log2FC, matched_pairs$control_gene_log2FC, paired = TRUE, exact = FALSE))

yguide_matched_summary <- tibble(
  analysis = "YAP_guide_matched_control",
  matching_strategy = "1:1 nearest-neighbor matching with replacement among seed-unmatched genes after excluding combined on-target mask G; distance used log10(baseMean+1), log10(representative 3'UTR length+1), and representative 3'UTR GC content from the longest reconstructed transcript per gene.",
  n_seed_genes_matched = nrow(matched_pairs),
  median_seed_gene_log2FC = median(matched_pairs$seed_gene_log2FC, na.rm = TRUE),
  median_control_gene_log2FC = median(matched_pairs$control_gene_log2FC, na.rm = TRUE),
  median_paired_difference = median(matched_diff, na.rm = TRUE),
  paired_wilcox_p = matched_wilcox$p.value,
  mean_matching_distance = mean(matched_pairs$matching_distance, na.rm = TRUE)
)
write.table(yguide_matched_summary, file.path(out_dir, "yguide_matched_control_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

matched_plot_df <- bind_rows(
  tibble(group = "YAP guide seed-matched genes", log2FC = matched_pairs$seed_gene_log2FC),
  tibble(group = "Matched seed-unmatched controls", log2FC = matched_pairs$control_gene_log2FC)
) %>% mutate(group = factor(group, levels = c("YAP guide seed-matched genes", "Matched seed-unmatched controls")))

yguide_match_plot <- ggplot(matched_plot_df, aes(x = log2FC, color = group, linetype = group)) +
  stat_ecdf(linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  scale_color_manual(values = c("YAP guide seed-matched genes" = "#B23A48", "Matched seed-unmatched controls" = "#4C6A92")) +
  scale_linetype_manual(values = c("YAP guide seed-matched genes" = "solid", "Matched seed-unmatched controls" = "dashed")) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
  labs(
    title = "Matched-control analysis of YAP-guide seed-associated repression",
    subtitle = sprintf("Paired Wilcoxon P = %s | median paired difference = %.3f", formatC(yguide_matched_summary$paired_wilcox_p, format = "e", digits = 2), yguide_matched_summary$median_paired_difference),
    x = "siYAP vs siScr log2 fold change",
    y = "Cumulative fraction of genes",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "yguide_matched_control_plot.png"), yguide_match_plot, width = 8, height = 6, dpi = 300)

for (nm in names(mask_list)) {
  base_df[[paste0("mask_", nm)]] <- base_df$gene_id_stable %in% mask_list[[nm]]$genes
}

regression_df <- base_df %>%
  filter(!is.na(yap_log2FC), !is.na(log_baseMean), !is.na(log_utr_length), !is.na(representative_utr_gc_content))

extract_coef <- function(model, term, model_name, outcome_name) {
  coefs <- coef(summary(model))
  if (!(term %in% rownames(coefs))) {
    return(tibble(model = model_name, outcome = outcome_name, term = term, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_, p_value = NA_real_))
  }
  tibble(
    model = model_name,
    outcome = outcome_name,
    term = term,
    estimate = unname(coefs[term, "Estimate"]),
    std_error = unname(coefs[term, "Std. Error"]),
    statistic = unname(coefs[term, ncol(coefs) - 1]),
    p_value = unname(coefs[term, ncol(coefs)])
  )
}

regression_df <- regression_df %>% mutate(yap_downregulated = !is.na(yap_padj) & yap_padj < 0.05 & yap_log2FC < 0)

lm_base <- lm(yap_log2FC ~ yguide_seed_match + log_baseMean + log_utr_length + representative_utr_gc_content, data = regression_df)
lm_combined <- lm(yap_log2FC ~ yguide_seed_match + log_baseMean + log_utr_length + representative_utr_gc_content + mask_G_combined_union, data = regression_df)
glm_combined <- glm(yap_downregulated ~ yguide_seed_match + log_baseMean + log_utr_length + representative_utr_gc_content + mask_G_combined_union, data = regression_df, family = binomial())

yguide_reg_results <- bind_rows(
  extract_coef(lm_base, "yguide_seed_matchTRUE", "lm_base_covariates", "yap_log2FC"),
  extract_coef(lm_combined, "yguide_seed_matchTRUE", "lm_plus_combined_mask", "yap_log2FC"),
  extract_coef(glm_combined, "yguide_seed_matchTRUE", "glm_plus_combined_mask", "yap_downregulated")
) %>% mutate(fdr_bh = p.adjust(p_value, method = "BH"))

write.table(yguide_reg_results, file.path(out_dir, "yguide_adjusted_regression_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

generate_null_seeds <- function(seed, n = 1000) {
  set.seed(20260414)
  nucs <- c("A", "U", "G", "C")
  gc_target <- sum(strsplit(seed, "")[[1]] %in% c("G", "C"))
  seeds <- character()
  while (length(seeds) < n) {
    cand <- paste0(sample(nucs, nchar(seed), replace = TRUE), collapse = "")
    if (cand != seed && sum(strsplit(cand, "")[[1]] %in% c("G", "C")) == gc_target) {
      seeds <- c(seeds, cand)
    }
  }
  unique(seeds)
}

longest_utr_genes <- base_df %>%
  filter(!is.na(representative_transcript_id), !is.na(representative_utr_length), representative_utr_length > 0) %>%
  dplyr::select(gene_id_stable, gene_symbol, yap_log2FC, taz_log2FC, representative_transcript_id)

utr_lookup <- utr_recon %>%
  dplyr::select(transcript_id, utr_seq) %>%
  distinct()

longest_utr_genes <- longest_utr_genes %>%
  left_join(utr_lookup, by = c("representative_transcript_id" = "transcript_id"))

seed_match_stat <- function(seed, lfc, utr_seqs) {
  motif <- chartr("AUGC", "UACG", seed)
  motif <- paste(rev(strsplit(motif, "")[[1]]), collapse = "")
  matched <- grepl(motif, utr_seqs, fixed = TRUE)
  if (sum(matched, na.rm = TRUE) < 10 || sum(!matched, na.rm = TRUE) < 10) {
    return(c(n_matched = sum(matched, na.rm = TRUE), median_diff = NA_real_, ks = NA_real_, wilcox_p = NA_real_))
  }
  c(
    n_matched = sum(matched, na.rm = TRUE),
    median_diff = median(lfc[matched], na.rm = TRUE) - median(lfc[!matched], na.rm = TRUE),
    ks = suppressWarnings(ks.test(lfc[matched], lfc[!matched])$statistic),
    wilcox_p = suppressWarnings(wilcox.test(lfc[matched], lfc[!matched], exact = FALSE)$p.value)
  )
}

run_null_seed_analysis <- function(seed, lfc, utr_df, label_prefix) {
  null_seeds <- generate_null_seeds(seed, n = 1000)
  null_seeds <- null_seeds[seq_len(min(1000, length(null_seeds)))]
  observed <- seed_match_stat(seed, lfc, utr_df$utr_seq)
  null_rows <- bind_rows(lapply(null_seeds, function(s) {
    vals <- seed_match_stat(s, lfc, utr_df$utr_seq)
    tibble(seed = s, type = "null", n_matched = vals["n_matched"], median_difference = vals["median_diff"], ks_statistic = vals["ks"], wilcox_p = vals["wilcox_p"])
  }))
  obs_row <- tibble(seed = seed, type = "observed", n_matched = observed["n_matched"], median_difference = observed["median_diff"], ks_statistic = observed["ks"], wilcox_p = observed["wilcox_p"])
  bind_rows(obs_row, null_rows) %>% mutate(analysis = label_prefix)
}

yguide_seed <- unique(sirna_seed_tx$YAP_guide_seed)[1]
tguide_seed <- unique(sirna_seed_tx$TAZ_guide_seed)[1]

yguide_null <- run_null_seed_analysis(yguide_seed, longest_utr_genes$yap_log2FC, longest_utr_genes, "YAP_guide")
write.table(yguide_null, file.path(out_dir, "yguide_null_seed_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

yguide_obs_stat <- yguide_null %>% filter(type == "observed") %>% pull(median_difference)
yguide_null_plot <- ggplot(yguide_null %>% filter(type == "null"), aes(x = median_difference)) +
  geom_histogram(bins = 50, fill = "grey75", color = "white") +
  geom_vline(xintercept = yguide_obs_stat, color = "#B23A48", linewidth = 1) +
  labs(
    title = "Observed YAP-guide seed effect compared with matched null-seed controls",
    subtitle = sprintf("Observed median difference = %.3f", yguide_obs_stat),
    x = "Median log2FC difference (seed-matched minus unmatched)",
    y = "Number of null seeds"
  ) +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "yguide_null_seed_plot.png"), yguide_null_plot, width = 8, height = 5.5, dpi = 300)

tguide_masks <- tibble(
  mask_name = c("baseline_none", "top50_taz_down", "top100_taz_down"),
  description = c(
    "Baseline: no mask excluded.",
    "Top 50 siTAZ downregulated genes ranked by adjusted P value then more negative log2FC.",
    "Top 100 siTAZ downregulated genes ranked by adjusted P value then more negative log2FC."
  ),
  genes = I(list(character(), head(taz_top_down$gene_id_stable, 50), head(taz_top_down$gene_id_stable, 100)))
)

tguide_masked <- bind_rows(lapply(seq_len(nrow(tguide_masks)), function(i) {
  row <- tguide_masks[i, ]
  sub <- base_df %>% filter(!(gene_id_stable %in% row$genes[[1]]))
  matched <- sub$taz_log2FC[sub$tguide_seed_match]
  unmatched <- sub$taz_log2FC[!sub$tguide_seed_match]
  ks_res <- suppressWarnings(ks.test(matched, unmatched))
  wilcox_res <- suppressWarnings(wilcox.test(matched, unmatched, exact = FALSE))
  tibble(
    analysis_component = "masked_cdf",
    mask_name = row$mask_name,
    description = row$description,
    n_seed_matched = length(matched),
    n_seed_unmatched = length(unmatched),
    median_difference = median(matched, na.rm = TRUE) - median(unmatched, na.rm = TRUE),
    ks_p = ks_res$p.value,
    wilcox_p = wilcox_res$p.value
  )
}))

t_reg_df <- base_df %>% filter(!is.na(taz_log2FC), !is.na(log_baseMean), !is.na(log_utr_length), !is.na(representative_utr_gc_content)) %>%
  mutate(taz_downregulated = !is.na(taz_padj) & taz_padj < 0.05 & taz_log2FC < 0,
         mask_taz_top100 = gene_id_stable %in% head(taz_top_down$gene_id_stable, 100))
t_lm <- lm(taz_log2FC ~ tguide_seed_match + log_baseMean + log_utr_length + representative_utr_gc_content + mask_taz_top100, data = t_reg_df)
t_glm <- glm(taz_downregulated ~ tguide_seed_match + log_baseMean + log_utr_length + representative_utr_gc_content + mask_taz_top100, data = t_reg_df, family = binomial())

tguide_reg <- bind_rows(
  extract_coef(t_lm, "tguide_seed_matchTRUE", "lm_plus_top100_mask", "taz_log2FC"),
  extract_coef(t_glm, "tguide_seed_matchTRUE", "glm_plus_top100_mask", "taz_downregulated")
) %>% mutate(analysis_component = "adjusted_regression", fdr_bh = p.adjust(p_value, method = "BH"))

tguide_null <- run_null_seed_analysis(tguide_seed, longest_utr_genes$taz_log2FC, longest_utr_genes, "TAZ_guide")
tguide_obs <- tguide_null %>% filter(type == "observed")
tguide_null_summary <- tibble(
  analysis_component = "null_seed",
  mask_name = "null_seed_comparison",
  description = "Observed TAZ-guide statistic compared with GC-matched random 6-mer null seeds using longest reconstructed 3'UTR per gene.",
  n_seed_matched = tguide_obs$n_matched,
  n_seed_unmatched = nrow(longest_utr_genes) - tguide_obs$n_matched,
  median_difference = tguide_obs$median_difference,
  ks_p = tguide_obs$ks_statistic,
  wilcox_p = tguide_obs$wilcox_p
)

tguide_summary <- bind_rows(
  tguide_masked %>% mutate(ks_p_bh = p.adjust(ks_p, method = "BH"), wilcox_p_bh = p.adjust(wilcox_p, method = "BH")),
  tguide_reg %>% transmute(
    analysis_component = analysis_component,
    mask_name = model,
    description = paste(outcome, term),
    n_seed_matched = NA_integer_,
    n_seed_unmatched = NA_integer_,
    median_difference = estimate,
    ks_p = p_value,
    wilcox_p = fdr_bh
  ),
  tguide_null_summary %>% mutate(ks_p_bh = NA_real_, wilcox_p_bh = NA_real_)
)
write.table(tguide_summary, file.path(out_dir, "tguide_sensitivity_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

yguide_baseline <- yguide_masked %>% filter(mask_name == "baseline_none")
yguide_combined <- yguide_masked %>% filter(mask_name == "G_combined_union")
yguide_null_pctl <- mean(yguide_null$median_difference[yguide_null$type == "null"] <= yguide_obs_stat, na.rm = TRUE)
tguide_null_pctl <- mean(tguide_null$median_difference[tguide_null$type == "null"] <= tguide_obs$median_difference, na.rm = TRUE)

memo_lines <- c(
  "Sensitivity analysis memo for on-target confounding of seed effects",
  "",
  "Overview:",
  "- The purpose of this analysis was to test whether the apparent YAP-guide seed-associated repression could be inflated by genuine on-target YAP-responsive biology.",
  "- Multiple exclusion masks, matched-control analyses, adjusted regressions, and null-seed comparisons were used.",
  "",
  "Mask definitions:",
  apply(mask_metadata %>% mutate(line = paste0("- ", mask_name, ": ", description, " (n=", n_mask_genes, ").")) %>% dplyr::select(line), 1, identity),
  "",
  sprintf("Baseline YAP-guide median difference (seed-matched minus unmatched) = %.3f.", yguide_baseline$median_difference),
  sprintf("After combined mask G exclusion, median difference = %.3f, compared with %.3f at baseline.", yguide_combined$median_difference, yguide_baseline$median_difference),
  sprintf("Baseline YAP-guide Wilcoxon FDR = %s; combined mask G Wilcoxon FDR = %s.", formatC(yguide_baseline$wilcox_p_bh, format = "e", digits = 2), formatC(yguide_combined$wilcox_p_bh, format = "e", digits = 2)),
  sprintf("Matched-control analysis after combined mask exclusion gave a median paired difference of %.3f with paired Wilcoxon P = %s.", yguide_matched_summary$median_paired_difference, formatC(yguide_matched_summary$paired_wilcox_p, format = "e", digits = 2)),
  sprintf("Adjusted regression results for the YAP-guide seed term remained %s after covariate and combined-mask adjustment in the linear model (estimate %.3f, FDR %s).",
          ifelse(yguide_reg_results$estimate[yguide_reg_results$model == "lm_plus_combined_mask"] < 0 && yguide_reg_results$fdr_bh[yguide_reg_results$model == "lm_plus_combined_mask"] < 0.05, "negative and statistically detectable", "weak or not clearly detectable"),
          yguide_reg_results$estimate[yguide_reg_results$model == "lm_plus_combined_mask"],
          formatC(yguide_reg_results$fdr_bh[yguide_reg_results$model == "lm_plus_combined_mask"], format = "e", digits = 2)),
  sprintf("The observed YAP-guide null-seed statistic lay at the %.1fth percentile of GC-matched random null seeds (more negative values indicate stronger left shift).", 100 * yguide_null_pctl),
  "",
  sprintf("For TAZ guide, the observed null-seed statistic lay at the %.1fth percentile of matched null seeds.", 100 * tguide_null_pctl),
  sprintf("The TAZ top-100-mask regression seed term was %.3f in the linear model.", tguide_reg$estimate[tguide_reg$model == "lm_plus_top100_mask"]),
  "",
  "Conservative interpretation:",
  if (abs(yguide_combined$median_difference) >= 0.85 * abs(yguide_baseline$median_difference)) {
    "- The YAP-guide seed effect remains very similar after excluding likely on-target biology masks, so it does not appear to be largely explained by these on-target gene sets."
  } else if (abs(yguide_combined$median_difference) >= 0.5 * abs(yguide_baseline$median_difference)) {
    "- The YAP-guide seed effect is partially attenuated after excluding likely on-target biology, suggesting some confounding but not a complete explanation."
  } else {
    "- The YAP-guide seed effect is markedly attenuated after excluding likely on-target biology, suggesting meaningful confounding by true on-target YAP-responsive genes."
  },
  "- The matched-control, regression, and null-seed analyses help reduce confounding but do not fully separate on-target from off-target effects.",
  "- If the YAP-guide effect remains more negative than matched controls and more extreme than most null seeds after adjustment, that supports a more robust seed-associated component.",
  "- The TAZ guide analysis provides a comparator and can be interpreted more lightly because the original TAZ seed effect was weaker."
)

writeLines(memo_lines, file.path(out_dir, "sensitivity_on_target_confounding_memo.txt"))
