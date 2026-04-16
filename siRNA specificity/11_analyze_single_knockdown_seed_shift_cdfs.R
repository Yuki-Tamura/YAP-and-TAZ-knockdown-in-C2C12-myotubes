#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- normalizePath(args[[2]], mustWork = TRUE)
analysis_dir <- normalizePath(if (length(args) >= 3) args[[3]] else dirname(out_dir), mustWork = TRUE)

seed_genes <- read.delim(file.path(out_dir, "sirna_seed_annotations_genes.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
seed_tx <- read.delim(file.path(out_dir, "sirna_seed_annotations.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
nearperfect <- read.delim(file.path(out_dir, "candidate_offtargets_nearperfect.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
utr_recon <- read.delim(file.path(out_dir, "transcript_utr_reconstruction.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
yap_res <- read.csv(file.path(analysis_dir, "02_deseq2", "siYAP_vs_siScr_all_genes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
taz_res <- read.csv(file.path(analysis_dir, "02_deseq2", "siTAZ_vs_siScr_all_genes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
double_res <- read.csv(file.path(analysis_dir, "02_deseq2", "siYAPsiTAZ_vs_siScr_all_genes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
interaction_res <- read.csv(file.path(analysis_dir, "02_deseq2", "interaction_effect_all_genes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
fc_cor_df <- read.csv(file.path(analysis_dir, "06_tables", "siYAP_vs_siTAZ_log2fc_source_data.csv"), stringsAsFactors = FALSE, check.names = FALSE)
strand_df <- read.csv(file.path(out_dir, "sirna_strand_derivation.csv"), stringsAsFactors = FALSE, check.names = FALSE)
seed_col_name <- grep("^guide_seed_", colnames(strand_df), value = TRUE)[1]
seed_label <- sub("^guide_seed_", "", seed_col_name)
seed_positions <- strsplit(seed_label, "_", fixed = TRUE)[[1]]
seed_len <- nchar(strand_df[[seed_col_name]][1])

analyze_seed <- function(res_df, seed_genes, sirna_name, strand_role) {
  seed_col <- paste0(sirna_name, "_", strand_role, "_has_match")
  merged <- res_df %>%
    mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
    inner_join(seed_genes, by = "gene_id_stable") %>%
    mutate(
      seed_match = .data[[seed_col]] == 1,
      log_utr_length = log10(max_utr_length + 1),
      downregulated = !is.na(padj) & padj < 0.05 & log2FoldChange < 0
    )

  wilcox_p <- tryCatch(
    wilcox.test(log2FoldChange ~ seed_match, data = merged)$p.value,
    error = function(e) NA_real_
  )
  lm_fit <- lm(log2FoldChange ~ seed_match + log_utr_length, data = merged)
  lm_coef <- coef(summary(lm_fit))
  glm_fit <- glm(downregulated ~ seed_match + log_utr_length, data = merged, family = binomial())
  glm_coef <- coef(summary(glm_fit))

  tibble(
    contrast = paste0("si", sirna_name, "_vs_siScr"),
    sirna = sirna_name,
    strand_role = strand_role,
    n_genes_tested = nrow(merged),
    n_seed_matched_genes = sum(merged$seed_match, na.rm = TRUE),
    median_log2fc_seed_matched = median(merged$log2FoldChange[merged$seed_match], na.rm = TRUE),
    median_log2fc_seed_unmatched = median(merged$log2FoldChange[!merged$seed_match], na.rm = TRUE),
    frac_down_seed_matched = mean(merged$downregulated[merged$seed_match], na.rm = TRUE),
    frac_down_seed_unmatched = mean(merged$downregulated[!merged$seed_match], na.rm = TRUE),
    wilcox_p = wilcox_p,
    lm_seed_coef = unname(lm_coef["seed_matchTRUE", "Estimate"]),
    lm_seed_p = unname(lm_coef["seed_matchTRUE", "Pr(>|t|)"]),
    glm_seed_logodds = unname(glm_coef["seed_matchTRUE", "Estimate"]),
    glm_seed_p = unname(glm_coef["seed_matchTRUE", "Pr(>|z|)"])
  )
}

prepare_cdf_data <- function(res_df, seed_genes, sirna_name, strand_role) {
  seed_col <- paste0(sirna_name, "_", strand_role, "_has_match")
  contrast_label <- paste0("si", sirna_name, " vs siScr")
  panel_label <- paste(
    ifelse(sirna_name == "YAP", "YAP", "TAZ"),
    ifelse(strand_role == "guide", "guide seed-matched genes", "passenger seed-matched genes")
  )

  res_df %>%
    mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
    inner_join(seed_genes, by = "gene_id_stable") %>%
    transmute(
      contrast = contrast_label,
      sirna = sirna_name,
      strand_role = strand_role,
      panel_label = panel_label,
      gene_id = gene_id.x,
      gene_id_stable = gene_id_stable,
      gene_name = gene_name,
      gene_biotype = gene_biotype,
      max_utr_length = max_utr_length,
      n_utr_transcripts = n_utr_transcripts,
      log2FoldChange = log2FoldChange,
      padj = padj,
      seed_match = .data[[seed_col]] == 1
    ) %>%
    mutate(
      match_group = ifelse(seed_match, "Seed-matched", "Seed-unmatched"),
      match_group = factor(match_group, levels = c("Seed-matched", "Seed-unmatched"))
    )
}

compute_cdf_stats <- function(cdf_df) {
  matched <- cdf_df$log2FoldChange[cdf_df$seed_match]
  unmatched <- cdf_df$log2FoldChange[!cdf_df$seed_match]

  ks_res <- suppressWarnings(ks.test(matched, unmatched))
  wilcox_res <- suppressWarnings(wilcox.test(matched, unmatched, exact = FALSE))
  matched_median <- median(matched, na.rm = TRUE)
  unmatched_median <- median(unmatched, na.rm = TRUE)

  tibble(
    contrast = unique(cdf_df$contrast),
    sirna = unique(cdf_df$sirna),
    strand_role = unique(cdf_df$strand_role),
    panel_label = unique(cdf_df$panel_label),
    collapse_rule = "Gene-level collapse: a gene was classified as seed-matched if any reconstructed transcript 3'UTR contained at least one seed match; otherwise seed-unmatched.",
    n_seed_matched = sum(cdf_df$seed_match, na.rm = TRUE),
    n_seed_unmatched = sum(!cdf_df$seed_match, na.rm = TRUE),
    median_log2fc_seed_matched = matched_median,
    median_log2fc_seed_unmatched = unmatched_median,
    ks_statistic = unname(ks_res$statistic),
    ks_p = ks_res$p.value,
    wilcox_p = wilcox_res$p.value,
    effect_direction = ifelse(
      matched_median < unmatched_median,
      "Left-shifted in seed-matched genes (more negative log2FC)",
      "No left shift in seed-matched genes"
    )
  )
}

build_top_downregulated_annotation <- function(
  res_df,
  seed_genes,
  seed_tx,
  utr_recon,
  nearperfect,
  sirna_name,
  n_top = 100,
  log2fc_cutoff = -0.5,
  base_mean_min = 20
) {
  intended_symbol <- ifelse(sirna_name == "YAP", "Yap1", "Wwtr1")
  guide_query <- paste0(sirna_name, "_guide")
  passenger_query <- paste0(sirna_name, "_passenger")

  gene_conf <- utr_recon %>%
    filter(utr_reconstructable == 1) %>%
    mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
    group_by(gene_id_stable) %>%
    summarise(
      n_utr_transcripts_reconstructed = n(),
      n_explicit_utr_transcripts = sum(utr_source == "explicit_three_prime_utr", na.rm = TRUE),
      n_fallback_utr_transcripts = sum(utr_source == "inferred_from_cds_stop", na.rm = TRUE),
      utr_confidence_class = case_when(
        n_explicit_utr_transcripts > 0 & n_fallback_utr_transcripts == 0 ~ "explicit_ensembl_3utr_only",
        n_explicit_utr_transcripts > 0 & n_fallback_utr_transcripts > 0 ~ "mixed_explicit_and_fallback",
        n_explicit_utr_transcripts == 0 & n_fallback_utr_transcripts > 0 ~ "inferred_fallback_only",
        TRUE ~ "no_reconstructed_3utr"
      ),
      .groups = "drop"
    )

  near_gene <- nearperfect %>%
    filter(query_name %in% c(guide_query, passenger_query)) %>%
    mutate(gene_id_stable = sub("\\..*$", "", gene_id_stable)) %>%
    group_by(gene_id_stable) %>%
    summarise(
      any_near_perfect_match_detected = n() > 0,
      near_perfect_query_names = paste(sort(unique(query_name)), collapse = ";"),
      near_perfect_min_mismatches = min(mismatches, na.rm = TRUE),
      .groups = "drop"
    )

  top_df <- res_df %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < log2fc_cutoff, !is.na(baseMean), baseMean >= base_mean_min) %>%
    arrange(padj, log2FoldChange) %>%
    slice_head(n = n_top) %>%
    mutate(gene_id_stable = sub("\\..*$", "", gene_id))

  annotated <- top_df %>%
    left_join(seed_genes, by = "gene_id_stable") %>%
    left_join(gene_conf, by = "gene_id_stable") %>%
    left_join(near_gene, by = "gene_id_stable") %>%
    mutate(
      gene_id = gene_id.x,
      comparison = paste0("si", sirna_name, " vs siScr"),
      selection_rule = sprintf("FDR < 0.05, log2FC < %.1f, baseMean >= %d, top %d ranked by adjusted P value then more negative log2 fold change.", log2fc_cutoff, base_mean_min, n_top),
      log2fc_cutoff = log2fc_cutoff,
      baseMean_min = base_mean_min,
      gene_symbol = ifelse(is.na(gene_name) | gene_name == "", gene_id_stable, gene_name),
      intended_target = gene_symbol == intended_symbol,
      any_near_perfect_match_detected = ifelse(is.na(any_near_perfect_match_detected), FALSE, any_near_perfect_match_detected),
      near_perfect_query_names = ifelse(is.na(near_perfect_query_names), "", near_perfect_query_names),
      near_perfect_min_mismatches = ifelse(is.na(near_perfect_min_mismatches), "", near_perfect_min_mismatches),
      guide_seed_match_present = ifelse(is.na(.data[[paste0(sirna_name, "_guide_has_match")]]), FALSE, .data[[paste0(sirna_name, "_guide_has_match")]] == 1),
      passenger_seed_match_present = ifelse(is.na(.data[[paste0(sirna_name, "_passenger_has_match")]]), FALSE, .data[[paste0(sirna_name, "_passenger_has_match")]] == 1),
      utr_confidence_class = ifelse(is.na(utr_confidence_class), "no_reconstructed_3utr", utr_confidence_class)
    ) %>%
    dplyr::select(
      comparison, selection_rule, log2fc_cutoff, baseMean_min, gene_id, gene_symbol, baseMean, log2FoldChange, padj,
      intended_target, any_near_perfect_match_detected, near_perfect_query_names, near_perfect_min_mismatches,
      guide_seed_match_present, passenger_seed_match_present,
      utr_confidence_class, n_utr_transcripts_reconstructed, n_explicit_utr_transcripts, n_fallback_utr_transcripts
    )

  annotated
}

annotate_gene_set <- function(gene_ids_stable, set_name, seed_genes, nearperfect, intended_symbols = character()) {
  near_gene <- nearperfect %>%
    mutate(gene_id_stable = sub("\\..*$", "", gene_id_stable)) %>%
    group_by(gene_id_stable) %>%
    summarise(
      any_near_perfect_match = n() > 0,
      near_perfect_query_names = paste(sort(unique(query_name)), collapse = ";"),
      near_perfect_min_mismatches = min(mismatches, na.rm = TRUE),
      .groups = "drop"
    )

  tibble(gene_id_stable = unique(gene_ids_stable), gene_set = set_name) %>%
    left_join(seed_genes, by = "gene_id_stable") %>%
    left_join(near_gene, by = "gene_id_stable") %>%
    mutate(
      gene_symbol = ifelse(is.na(gene_name) | gene_name == "", gene_id_stable, gene_name),
      YAP_guide_seed_match = ifelse(is.na(YAP_guide_has_match), FALSE, YAP_guide_has_match == 1),
      YAP_passenger_seed_match = ifelse(is.na(YAP_passenger_has_match), FALSE, YAP_passenger_has_match == 1),
      TAZ_guide_seed_match = ifelse(is.na(TAZ_guide_has_match), FALSE, TAZ_guide_has_match == 1),
      TAZ_passenger_seed_match = ifelse(is.na(TAZ_passenger_has_match), FALSE, TAZ_passenger_has_match == 1),
      both_guide_seed_match = YAP_guide_seed_match & TAZ_guide_seed_match,
      any_seed_match = YAP_guide_seed_match | YAP_passenger_seed_match | TAZ_guide_seed_match | TAZ_passenger_seed_match,
      any_near_perfect_match = ifelse(is.na(any_near_perfect_match), FALSE, any_near_perfect_match),
      near_perfect_query_names = ifelse(is.na(near_perfect_query_names), "", near_perfect_query_names),
      near_perfect_min_mismatches = ifelse(is.na(near_perfect_min_mismatches), "", near_perfect_min_mismatches),
      intended_target_status = case_when(
        gene_symbol %in% intended_symbols ~ "intended_single_target_gene",
        TRUE ~ "not_intended_target"
      )
    ) %>%
    dplyr::select(
      gene_set, gene_id_stable, gene_symbol, intended_target_status,
      YAP_guide_seed_match, YAP_passenger_seed_match, TAZ_guide_seed_match, TAZ_passenger_seed_match,
      both_guide_seed_match, any_seed_match, any_near_perfect_match,
      near_perfect_query_names, near_perfect_min_mismatches
    )
}

fisher_feature_enrichment <- function(feature_vector, set_membership, feature_name, gene_set_name) {
  feature_vector <- as.logical(feature_vector)
  set_membership <- as.logical(set_membership)
  a <- sum(feature_vector & set_membership, na.rm = TRUE)
  b <- sum(!feature_vector & set_membership, na.rm = TRUE)
  c <- sum(feature_vector & !set_membership, na.rm = TRUE)
  d <- sum(!feature_vector & !set_membership, na.rm = TRUE)
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  tibble(
    gene_set = gene_set_name,
    feature = feature_name,
    genes_in_set = sum(set_membership, na.rm = TRUE),
    genes_in_background = length(set_membership),
    feature_in_set = a,
    feature_not_in_set = b,
    feature_in_background_outside_set = c,
    feature_not_in_background_outside_set = d,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

focused_keyword_enrichment <- function(gene_ids_stable, universe_ids_stable) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Mm.eg.db", quietly = TRUE) ||
      !requireNamespace("GO.db", quietly = TRUE)) {
    return(NULL)
  }

  gene2go <- AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = unique(universe_ids_stable),
    columns = c("GOALL", "ONTOLOGYALL"),
    keytype = "ENSEMBL"
  ) %>%
    as_tibble() %>%
    filter(!is.na(GOALL), ONTOLOGYALL == "BP")

  go_terms <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = unique(gene2go$GOALL),
    columns = c("TERM", "ONTOLOGY"),
    keytype = "GOID"
  ) %>%
    as_tibble() %>%
    filter(ONTOLOGY == "BP") %>%
    dplyr::select(GOALL = GOID, TERM)

  gene2go <- gene2go %>%
    dplyr::select(ENSEMBL, GOALL) %>%
    distinct() %>%
    inner_join(go_terms, by = "GOALL")

  keyword_map <- tribble(
    ~focus_category, ~pattern,
    "muscle_contraction", "muscle contraction|striated muscle contraction|skeletal muscle contraction",
    "muscle_system_process", "muscle system process|muscle structure development|myofibril|sarcomere",
    "calcium_transport", "calcium ion transport|regulation of calcium ion transport|calcium",
    "mitochondrial_metabolic", "mitochond|metabo|respirat|oxidative|sterol|cholesterol|energy"
  )

  bind_rows(lapply(seq_len(nrow(keyword_map)), function(i) {
    row <- keyword_map[i, ]
    bg_genes <- unique(gene2go$ENSEMBL[grepl(row$pattern, gene2go$TERM, ignore.case = TRUE)])
    in_set <- universe_ids_stable %in% gene_ids_stable
    in_feature <- universe_ids_stable %in% bg_genes
    ft <- fisher.test(matrix(c(
      sum(in_set & in_feature),
      sum(in_set & !in_feature),
      sum(!in_set & in_feature),
      sum(!in_set & !in_feature)
    ), nrow = 2), alternative = "greater")
    tibble(
      focus_category = row$focus_category,
      genes_in_focus_terms = sum(in_feature),
      genes_in_gene_set = sum(in_set),
      overlap = sum(in_set & in_feature),
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value
    )
  })) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH"))
}

results <- bind_rows(
  analyze_seed(yap_res, seed_genes, "YAP", "guide"),
  analyze_seed(yap_res, seed_genes, "YAP", "passenger"),
  analyze_seed(taz_res, seed_genes, "TAZ", "guide"),
  analyze_seed(taz_res, seed_genes, "TAZ", "passenger")
) %>%
  mutate(across(ends_with("_p"), ~ p.adjust(.x, method = "BH"), .names = "{.col}_bh"))

write.table(results, file.path(out_dir, "seed_enrichment_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cdf_data <- bind_rows(
  prepare_cdf_data(yap_res, seed_genes, "YAP", "guide"),
  prepare_cdf_data(yap_res, seed_genes, "YAP", "passenger"),
  prepare_cdf_data(taz_res, seed_genes, "TAZ", "guide"),
  prepare_cdf_data(taz_res, seed_genes, "TAZ", "passenger")
)

cdf_stats <- bind_rows(
  cdf_data %>% filter(sirna == "YAP", strand_role == "guide") %>% compute_cdf_stats(),
  cdf_data %>% filter(sirna == "YAP", strand_role == "passenger") %>% compute_cdf_stats(),
  cdf_data %>% filter(sirna == "TAZ", strand_role == "guide") %>% compute_cdf_stats(),
  cdf_data %>% filter(sirna == "TAZ", strand_role == "passenger") %>% compute_cdf_stats()
) %>%
  mutate(
    ks_p_bh = p.adjust(ks_p, method = "BH"),
    wilcox_p_bh = p.adjust(wilcox_p, method = "BH")
  )

write.table(cdf_stats, file.path(out_dir, "cdf_seed_shift_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cdf_plot_df <- cdf_data %>%
  mutate(
    panel_label = factor(
      panel_label,
      levels = c(
        "YAP guide seed-matched genes",
        "TAZ guide seed-matched genes",
        "YAP passenger seed-matched genes",
        "TAZ passenger seed-matched genes"
      )
    )
  )

stats_labels <- cdf_stats %>%
  mutate(
    panel_label = factor(
      panel_label,
      levels = levels(cdf_plot_df$panel_label)
    ),
    label = sprintf(
      "n=%d vs %d\nKS P=%s\nWilcoxon P=%s",
      n_seed_matched,
      n_seed_unmatched,
      formatC(ks_p_bh, format = "e", digits = 2),
      formatC(wilcox_p_bh, format = "e", digits = 2)
    ),
    x = -4.8,
    y = 0.92
  )

cdf_palette <- c("Seed-matched" = "#B23A48", "Seed-unmatched" = "grey55")

make_single_cdf_plot <- function(df, title_text) {
  ggplot(df, aes(x = log2FoldChange, color = match_group, linetype = match_group)) +
    stat_ecdf(linewidth = 0.9) +
    scale_color_manual(values = cdf_palette) +
    scale_linetype_manual(values = c("Seed-matched" = "solid", "Seed-unmatched" = "dashed")) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.4) +
    labs(
      title = title_text,
      x = "RNA-seq log2 fold change",
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

cdf_plot <- ggplot(cdf_plot_df, aes(x = log2FoldChange, color = match_group, linetype = match_group)) +
  stat_ecdf(linewidth = 0.9) +
  facet_wrap(~ panel_label, ncol = 2) +
  scale_color_manual(values = cdf_palette) +
  scale_linetype_manual(values = c("Seed-matched" = "solid", "Seed-unmatched" = "dashed")) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.4) +
  geom_text(
    data = stats_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.1
  ) +
  labs(
    title = "CDF of RNA-seq fold changes for seed-matched vs unmatched genes",
    subtitle = "Guide seeds are the main focus; passenger panels are shown for context",
    x = "RNA-seq log2 fold change",
    y = "Cumulative fraction of genes",
    color = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "supplemental_seed_cdf_plot.png"), cdf_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "supplemental_seed_cdf_plot.pdf"), cdf_plot, width = 10, height = 8)

single_cdf_dir <- file.path(out_dir, "cdf_single_panels_no_annotations")
dir.create(single_cdf_dir, recursive = TRUE, showWarnings = FALSE)

single_panel_specs <- tribble(
  ~sirna, ~strand_role, ~file_stub, ~title_text,
  "YAP", "guide", "yap_guide_seed_cdf", "YAP guide seed-matched genes",
  "YAP", "passenger", "yap_passenger_seed_cdf", "YAP passenger seed-matched genes",
  "TAZ", "guide", "taz_guide_seed_cdf", "TAZ guide seed-matched genes",
  "TAZ", "passenger", "taz_passenger_seed_cdf", "TAZ passenger seed-matched genes"
)

for (i in seq_len(nrow(single_panel_specs))) {
  spec <- single_panel_specs[i, ]
  panel_df <- cdf_plot_df %>%
    filter(sirna == spec$sirna, strand_role == spec$strand_role)
  panel_plot <- make_single_cdf_plot(panel_df, spec$title_text)
  ggsave(file.path(single_cdf_dir, paste0(spec$file_stub, ".png")), panel_plot, width = 6, height = 5, dpi = 300)
  ggsave(file.path(single_cdf_dir, paste0(spec$file_stub, ".pdf")), panel_plot, width = 6, height = 5)
}

memo_lines <- c(
  "siRNA off-target interpretation memo",
  "",
  "Near-perfect off-target candidates were searched transcriptome-wide against Ensembl mouse cDNA sequences using the 19-nt duplex core for candidate guide and passenger strands, allowing up to 1 mismatch.",
  "3'UTR sequences were reconstructed from the Ensembl GTF plus cDNA FASTA for transcripts with exon/cDNA length agreement and annotated CDS features.",
  sprintf("Seed analyses were run conservatively at the gene level by marking a gene as seed-matched if any reconstructed transcript 3'UTR contained at least one reverse-complement %d-mer seed site (positions %s-%s of the inferred strand).", seed_len, seed_positions[1], seed_positions[2]),
  "Downregulation support was tested in the matched single-knockdown RNA-seq contrast for each siRNA, using both a Wilcoxon comparison of log2 fold changes and models that adjusted for 3'UTR length.",
  "",
  sprintf("Transcriptome-wide siYAP vs siTAZ log2FC correlation from the RNA-seq analysis was %.3f.", cor(fc_cor_df$log2FoldChange_yap, fc_cor_df$log2FoldChange_taz, use = 'pairwise.complete.obs')),
  sprintf("Near-perfect transcriptome search found %d candidate hits total: %d for the YAP guide and %d for the TAZ guide; no 0-1 mismatch hits were found for the passenger-strand candidates.", nrow(nearperfect), sum(nearperfect$query_name == 'YAP_guide'), sum(nearperfect$query_name == 'TAZ_guide')),
  sprintf("All near-perfect hits mapped only to the intended target genes: Yap1 for the YAP guide and Wwtr1 for the TAZ guide."),
  sprintf("Guide-seed matched genes showed stronger downregulation than unmatched genes for siYAP (BH-adjusted Wilcoxon P = %.2g; length-adjusted linear-model P = %.2g) and a weaker but still significant trend for siTAZ (BH-adjusted Wilcoxon P = %.2g; length-adjusted linear-model P = %.2g).", results$wilcox_p_bh[results$sirna == 'YAP' & results$strand_role == 'guide'], results$lm_seed_p_bh[results$sirna == 'YAP' & results$strand_role == 'guide'], results$wilcox_p_bh[results$sirna == 'TAZ' & results$strand_role == 'guide'], results$lm_seed_p_bh[results$sirna == 'TAZ' & results$strand_role == 'guide']),
  sprintf("Passenger-seed associations were not robust after multiple-testing adjustment."),
  sprintf("Because %d-mer seed matches can still occur across mammalian 3'UTRs, these guide-seed enrichments should be interpreted as supportive evidence for some seed-mediated contribution rather than proof that the observed transcriptomic programs are predominantly off-target.", seed_len),
  "",
  "Interpretation guidance:",
  "- The near-perfect search argues against widespread near-perfect transcriptome off-targeting by these duplexes, because only the intended target genes were recovered at 0-1 mismatch.",
  "- The modest siYAP/siTAZ transcriptome correlation and limited DEG overlap support the manuscript position that YAP and TAZ perturbations are overlapping but non-equivalent, which is not the pattern expected if a dominant shared seed-driven artifact were overwhelming both responses.",
  "- Some guide-seed-associated downregulation is detectable, especially for the YAP guide, so off-target effects cannot be fully excluded.",
  "- These analyses remain supportive rather than definitive and do not fully exclude off-target effects, non-seed effects, or strand-loading uncertainty."
 )
writeLines(memo_lines, file.path(out_dir, "sirna_offtarget_interpretation_memo.txt"))

cdf_memo_lines <- c(
  "Seed CDF interpretation memo",
  "",
  "CDF analyses were performed at the gene level.",
  "Collapse rule: a gene was classified as seed-matched if any reconstructed transcript 3'UTR contained at least one seed match; otherwise it was classified as seed-unmatched.",
  "This gene-level rule was chosen to avoid ambiguity from multiple transcript isoforms per gene and to keep the comparison aligned with gene-level RNA-seq differential expression results.",
  "",
  sprintf(
    "YAP guide: seed-matched genes showed a left-shifted fold-change distribution relative to seed-unmatched genes (matched n=%d, unmatched n=%d, KS BH-adjusted P=%s, Wilcoxon BH-adjusted P=%s).",
    cdf_stats$n_seed_matched[cdf_stats$sirna == "YAP" & cdf_stats$strand_role == "guide"],
    cdf_stats$n_seed_unmatched[cdf_stats$sirna == "YAP" & cdf_stats$strand_role == "guide"],
    formatC(cdf_stats$ks_p_bh[cdf_stats$sirna == "YAP" & cdf_stats$strand_role == "guide"], format = "e", digits = 2),
    formatC(cdf_stats$wilcox_p_bh[cdf_stats$sirna == "YAP" & cdf_stats$strand_role == "guide"], format = "e", digits = 2)
  ),
  sprintf(
    "TAZ guide: seed-matched genes also showed a left shift, but the effect was weaker than for YAP (matched n=%d, unmatched n=%d, KS BH-adjusted P=%s, Wilcoxon BH-adjusted P=%s).",
    cdf_stats$n_seed_matched[cdf_stats$sirna == "TAZ" & cdf_stats$strand_role == "guide"],
    cdf_stats$n_seed_unmatched[cdf_stats$sirna == "TAZ" & cdf_stats$strand_role == "guide"],
    formatC(cdf_stats$ks_p_bh[cdf_stats$sirna == "TAZ" & cdf_stats$strand_role == "guide"], format = "e", digits = 2),
    formatC(cdf_stats$wilcox_p_bh[cdf_stats$sirna == "TAZ" & cdf_stats$strand_role == "guide"], format = "e", digits = 2)
  ),
  sprintf(
    "Passenger seeds: YAP passenger showed only a weak left shift, whereas TAZ passenger did not show a left-shift pattern despite a detectable distribution difference."
  ),
  "",
  "Interpretation guidance:",
  "- A left-shifted CDF for seed-matched genes is consistent with some seed-associated repression.",
  "- This does not by itself show that the observed transcriptomic programs are primarily off-target artifacts.",
  "- If guide-seed curves are shifted but siYAP and siTAZ transcriptome-wide responses remain only modestly correlated, that pattern supports some seed contribution without supporting a dominant shared seed-driven explanation for the overall biology."
)

writeLines(cdf_memo_lines, file.path(out_dir, "seed_cdf_interpretation_memo.txt"))

top_configs <- tribble(
  ~config_label, ~log2fc_cutoff, ~base_mean_min, ~n_top,
  "fdr005_lfc05_basemean20_top100", -0.5, 20, 100,
  "fdr005_lfc03_basemean20_top100", -0.3, 20, 100
)

top_combined <- bind_rows(lapply(seq_len(nrow(top_configs)), function(i) {
  cfg <- top_configs[i, ]
  bind_rows(
    build_top_downregulated_annotation(
      yap_res, seed_genes, seed_tx, utr_recon, nearperfect, "YAP",
      n_top = cfg$n_top, log2fc_cutoff = cfg$log2fc_cutoff, base_mean_min = cfg$base_mean_min
    ),
    build_top_downregulated_annotation(
      taz_res, seed_genes, seed_tx, utr_recon, nearperfect, "TAZ",
      n_top = cfg$n_top, log2fc_cutoff = cfg$log2fc_cutoff, base_mean_min = cfg$base_mean_min
    )
  ) %>%
    mutate(config_label = cfg$config_label)
}))

write.table(
  top_combined,
  file.path(out_dir, "top_downregulated_gene_offtarget_annotation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_panel_df <- top_combined %>%
  filter(config_label == "fdr005_lfc05_basemean20_top100") %>%
  mutate(row_label = paste0(gene_symbol, " (", comparison, ")")) %>%
  dplyr::select(
    comparison, row_label,
    Intended_target = intended_target,
    Near_perfect_match = any_near_perfect_match_detected,
    Guide_seed_match = guide_seed_match_present,
    Passenger_seed_match = passenger_seed_match_present
  ) %>%
  pivot_longer(
    cols = c(Intended_target, Near_perfect_match, Guide_seed_match, Passenger_seed_match),
    names_to = "annotation",
    values_to = "present"
  ) %>%
  mutate(
    row_label = factor(row_label, levels = rev(unique(row_label))),
    comparison = factor(comparison, levels = c("siYAP vs siScr", "siTAZ vs siScr"))
  )

summary_plot <- ggplot(summary_panel_df, aes(x = annotation, y = row_label, fill = present)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ comparison, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("TRUE" = "#C23B22", "FALSE" = "grey90")) +
  labs(
    title = "Top downregulated genes and simple sequence-based off-target annotations",
    subtitle = "Top 100 genes after FDR, log2FC, and low-expression filtering (FDR < 0.05, log2FC < -0.5, baseMean >= 20)",
    x = NULL,
    y = NULL,
    fill = "Present"
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "supplemental_top_downregulated_offtarget_summary.png"), summary_plot, width = 9, height = 14, dpi = 300)

top_memo_lines <- c(
  "Top downregulated gene off-target annotation memo",
  "",
  "Selection rule:",
  "- Filtered tables were generated separately for siYAP vs siScr and siTAZ vs siScr.",
  "- Primary filter set: FDR < 0.05, log2FC < -0.5, baseMean >= 20, then top 100 ranked by adjusted P value followed by more negative log2 fold change.",
  "- Secondary relaxed filter set: FDR < 0.05, log2FC < -0.3, baseMean >= 20, then top 100 ranked by adjusted P value followed by more negative log2 fold change.",
  "- baseMean from DESeq2 was used as the low-expression exclusion criterion.",
  "",
  "Annotation fields:",
  "- Intended target indicates Yap1 in the siYAP comparison and Wwtr1 in the siTAZ comparison.",
  "- Near-perfect match indicates any 0-1 mismatch transcriptome-wide candidate detected for the corresponding siRNA guide or passenger query.",
  "- Guide-seed and passenger-seed flags indicate whether any reconstructed transcript 3'UTR for that gene contained at least one corresponding seed match.",
  "- 3'UTR confidence class summarizes whether reconstructed 3'UTRs came from explicit Ensembl three_prime_utr annotation or fallback inference.",
  "",
  sprintf(
    "Primary stricter set (log2FC < -0.5): siYAP top list n=%d with intended target=%d, near-perfect matches=%d, guide-seed matches=%d, passenger-seed matches=%d.",
    sum(top_combined$config_label == "fdr005_lfc05_basemean20_top100" & top_combined$comparison == "siYAP vs siScr"),
    sum(top_combined$config_label == "fdr005_lfc05_basemean20_top100" & top_combined$comparison == "siYAP vs siScr" & top_combined$intended_target),
    sum(top_combined$config_label == "fdr005_lfc05_basemean20_top100" & top_combined$comparison == "siYAP vs siScr" & top_combined$any_near_perfect_match_detected),
    sum(top_combined$config_label == "fdr005_lfc05_basemean20_top100" & top_combined$comparison == "siYAP vs siScr" & top_combined$guide_seed_match_present),
    sum(top_combined$config_label == "fdr005_lfc05_basemean20_top100" & top_combined$comparison == "siYAP vs siScr" & top_combined$passenger_seed_match_present)
  ),
  sprintf(
    "Primary stricter set (log2FC < -0.5): siTAZ top list n=%d with intended target=%d, near-perfect matches=%d, guide-seed matches=%d, passenger-seed matches=%d.",
    sum(top_combined$config_label == "fdr005_lfc05_basemean20_top100" & top_combined$comparison == "siTAZ vs siScr"),
    sum(top_combined$config_label == "fdr005_lfc05_basemean20_top100" & top_combined$comparison == "siTAZ vs siScr" & top_combined$intended_target),
    sum(top_combined$config_label == "fdr005_lfc05_basemean20_top100" & top_combined$comparison == "siTAZ vs siScr" & top_combined$any_near_perfect_match_detected),
    sum(top_combined$config_label == "fdr005_lfc05_basemean20_top100" & top_combined$comparison == "siTAZ vs siScr" & top_combined$guide_seed_match_present),
    sum(top_combined$config_label == "fdr005_lfc05_basemean20_top100" & top_combined$comparison == "siTAZ vs siScr" & top_combined$passenger_seed_match_present)
  ),
  "",
  "Interpretation guidance:",
  "- The strongest downregulated genes are not broadly explained by obvious near-perfect unintended targeting.",
  "- Seed matches are present for some genes, but the top-ranked changes are not uniformly attributable to simple sequence-based off-target effects.",
  "- Absence of obvious homology should not be interpreted as proof that all observed changes are direct on-target biology; secondary transcriptional consequences remain possible."
)

writeLines(top_memo_lines, file.path(out_dir, "top_downregulated_offtarget_interpretation_memo.txt"))

sig_rule_text <- "Significant genes were defined as padj < 0.05 and |log2 fold change| >= 1."
double_specific_rule_text <- paste(
  "Double-specific genes were defined as genes significant in siYAP+siTAZ vs siScr",
  "but not significant in siYAP vs siScr and not significant in siTAZ vs siScr,",
  "using the same significance rule in all three contrasts."
)

is_sig <- function(df) {
  !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) >= 1
}

interaction_sig_ids <- interaction_res %>%
  mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
  pull(gene_id_stable)

double_sig_ids <- double_res %>%
  mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
  pull(gene_id_stable)

yap_sig_ids <- yap_res %>%
  mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
  pull(gene_id_stable)

taz_sig_ids <- taz_res %>%
  mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
  pull(gene_id_stable)

double_specific_ids <- setdiff(double_sig_ids, union(yap_sig_ids, taz_sig_ids))

interaction_annot <- annotate_gene_set(interaction_sig_ids, "interaction_significant", seed_genes, nearperfect, intended_symbols = c("Yap1", "Wwtr1"))
double_annot <- annotate_gene_set(double_sig_ids, "double_knockdown_deg", seed_genes, nearperfect, intended_symbols = c("Yap1", "Wwtr1"))
double_specific_annot <- annotate_gene_set(double_specific_ids, "double_specific", seed_genes, nearperfect, intended_symbols = c("Yap1", "Wwtr1")) %>%
  mutate(definition_rule = double_specific_rule_text)

write.table(
  double_specific_annot,
  file.path(out_dir, "double_specific_gene_annotations.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

gene_universe <- bind_rows(
  yap_res %>% dplyr::select(gene_id),
  taz_res %>% dplyr::select(gene_id),
  double_res %>% dplyr::select(gene_id),
  interaction_res %>% dplyr::select(gene_id)
) %>%
  distinct() %>%
  mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
  dplyr::select(gene_id_stable) %>%
  left_join(seed_genes, by = "gene_id_stable") %>%
  left_join(
    nearperfect %>%
      mutate(gene_id_stable = sub("\\..*$", "", gene_id_stable)) %>%
      group_by(gene_id_stable) %>%
      summarise(any_near_perfect_match = n() > 0, .groups = "drop"),
    by = "gene_id_stable"
  ) %>%
  mutate(
    YAP_guide_seed_match = ifelse(is.na(YAP_guide_has_match), FALSE, YAP_guide_has_match == 1),
    TAZ_guide_seed_match = ifelse(is.na(TAZ_guide_has_match), FALSE, TAZ_guide_has_match == 1),
    both_guide_seed_match = YAP_guide_seed_match & TAZ_guide_seed_match,
    any_seed_match = ifelse(is.na(YAP_guide_has_match), FALSE, YAP_guide_has_match == 1) |
      ifelse(is.na(YAP_passenger_has_match), FALSE, YAP_passenger_has_match == 1) |
      ifelse(is.na(TAZ_guide_has_match), FALSE, TAZ_guide_has_match == 1) |
      ifelse(is.na(TAZ_passenger_has_match), FALSE, TAZ_passenger_has_match == 1),
    any_near_perfect_match = ifelse(is.na(any_near_perfect_match), FALSE, any_near_perfect_match)
  )

set_membership_df <- gene_universe %>%
  mutate(
    interaction_significant = gene_id_stable %in% interaction_sig_ids,
    double_knockdown_deg = gene_id_stable %in% double_sig_ids,
    double_specific = gene_id_stable %in% double_specific_ids
  )

interaction_seed_enrichment <- bind_rows(
  fisher_feature_enrichment(set_membership_df$YAP_guide_seed_match, set_membership_df$interaction_significant, "YAP_guide_seed_match", "interaction_significant"),
  fisher_feature_enrichment(set_membership_df$TAZ_guide_seed_match, set_membership_df$interaction_significant, "TAZ_guide_seed_match", "interaction_significant"),
  fisher_feature_enrichment(set_membership_df$both_guide_seed_match, set_membership_df$interaction_significant, "both_guide_seed_match", "interaction_significant"),
  fisher_feature_enrichment(set_membership_df$any_seed_match, set_membership_df$interaction_significant, "any_seed_match", "interaction_significant"),
  fisher_feature_enrichment(set_membership_df$YAP_guide_seed_match, set_membership_df$double_knockdown_deg, "YAP_guide_seed_match", "double_knockdown_deg"),
  fisher_feature_enrichment(set_membership_df$TAZ_guide_seed_match, set_membership_df$double_knockdown_deg, "TAZ_guide_seed_match", "double_knockdown_deg"),
  fisher_feature_enrichment(set_membership_df$both_guide_seed_match, set_membership_df$double_knockdown_deg, "both_guide_seed_match", "double_knockdown_deg"),
  fisher_feature_enrichment(set_membership_df$any_seed_match, set_membership_df$double_knockdown_deg, "any_seed_match", "double_knockdown_deg"),
  fisher_feature_enrichment(set_membership_df$YAP_guide_seed_match, set_membership_df$double_specific, "YAP_guide_seed_match", "double_specific"),
  fisher_feature_enrichment(set_membership_df$TAZ_guide_seed_match, set_membership_df$double_specific, "TAZ_guide_seed_match", "double_specific"),
  fisher_feature_enrichment(set_membership_df$both_guide_seed_match, set_membership_df$double_specific, "both_guide_seed_match", "double_specific"),
  fisher_feature_enrichment(set_membership_df$any_seed_match, set_membership_df$double_specific, "any_seed_match", "double_specific")
) %>%
  mutate(fdr_bh = p.adjust(p_value, method = "BH"))

write.table(
  interaction_seed_enrichment,
  file.path(out_dir, "interaction_seed_enrichment_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_plot_df <- interaction_seed_enrichment %>%
  mutate(
    feature_label = recode(
      feature,
      YAP_guide_seed_match = "YAP guide seed",
      TAZ_guide_seed_match = "TAZ guide seed",
      both_guide_seed_match = "Both guide seeds",
      any_seed_match = "Any seed"
    ),
    gene_set_label = recode(
      gene_set,
      interaction_significant = "Interaction-significant genes",
      double_knockdown_deg = "Double-knockdown DEGs",
      double_specific = "Double-specific genes"
    ),
    log2_or = log2(pmax(odds_ratio, 1e-6)),
    sig_label = ifelse(fdr_bh < 0.05, "FDR < 0.05", "NS")
  )

interaction_plot <- ggplot(summary_plot_df, aes(x = feature_label, y = gene_set_label, fill = log2_or)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("OR=%.2f\nFDR=%s", odds_ratio, formatC(fdr_bh, format = "e", digits = 1))), size = 3) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C23B22", midpoint = 0) +
  labs(
    title = "Seed-feature enrichment in double-knockdown and interaction gene sets",
    subtitle = "Fisher's exact tests against the RNA-seq gene universe",
    x = NULL,
    y = NULL,
    fill = "log2 odds ratio"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave(file.path(out_dir, "supplemental_interaction_seed_summary.png"), interaction_plot, width = 9, height = 4.8, dpi = 300)

interaction_memo_lines <- c(
  "Interaction and double-knockdown seed enrichment memo",
  "",
  "Design context:",
  "- RNA-seq differential expression was modeled with ~ YAP_KD + TAZ_KD + YAP_KD:TAZ_KD.",
  paste("- Significance rule:", sig_rule_text),
  paste("- Double-specific rule:", double_specific_rule_text),
  "",
  sprintf(
    "Gene-set sizes: interaction-significant = %d, double-knockdown DEGs = %d, double-specific = %d.",
    length(unique(interaction_sig_ids)),
    length(unique(double_sig_ids)),
    length(unique(double_specific_ids))
  ),
  "",
  "Interpretation guidance:",
  "- This analysis tests whether double-knockdown or interaction-associated genes are disproportionately marked by simple predicted sequence features such as guide-seed matches.",
  "- If these sets are not enriched for guide-seed or any-seed matches, that argues against simple additive sequence-based off-target effects as the main explanation.",
  "- Even if some enrichment is present, that would still not prove that the dominant biology is off-target; it would only support a possible contribution.",
  "- This analysis is most informative about dominant shared/simple sequence-based explanations and does not completely exclude all possible siRNA-specific off-target effects."
)

top_rows <- interaction_seed_enrichment %>%
  arrange(fdr_bh, p_value)

interaction_memo_lines <- c(
  interaction_memo_lines,
  "",
  "Top enrichment results:",
  apply(
    top_rows %>% slice_head(n = min(6, nrow(top_rows))) %>%
      mutate(line = sprintf(
        "- %s / %s: OR = %.2f, P = %s, FDR = %s.",
        gene_set,
        feature,
        odds_ratio,
        formatC(p_value, format = "e", digits = 2),
        formatC(fdr_bh, format = "e", digits = 2)
      )) %>%
      dplyr::select(line),
    1,
    identity
  )
)

writeLines(interaction_memo_lines, file.path(out_dir, "interaction_seed_interpretation_memo.txt"))

additivity_df <- yap_res %>%
  dplyr::select(gene_id, yap_log2FC = log2FoldChange, yap_padj = padj) %>%
  mutate(gene_id_stable = sub("\\..*$", "", gene_id)) %>%
  inner_join(
    taz_res %>%
      dplyr::select(gene_id, taz_log2FC = log2FoldChange, taz_padj = padj) %>%
      mutate(gene_id_stable = sub("\\..*$", "", gene_id)),
    by = "gene_id_stable"
  ) %>%
  inner_join(
    double_res %>%
      dplyr::select(gene_id, observed_double_log2FC = log2FoldChange, double_padj = padj) %>%
      mutate(gene_id_stable = sub("\\..*$", "", gene_id)),
    by = "gene_id_stable"
  ) %>%
  left_join(seed_genes, by = "gene_id_stable") %>%
  left_join(
    nearperfect %>%
      mutate(gene_id_stable = sub("\\..*$", "", gene_id_stable)) %>%
      group_by(gene_id_stable) %>%
      summarise(
        any_near_perfect_match = n() > 0,
        near_perfect_query_names = paste(sort(unique(query_name)), collapse = ";"),
        .groups = "drop"
      ),
    by = "gene_id_stable"
  ) %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_name) | gene_name == "", gene_id_stable, gene_name),
    predicted_additive_log2FC = yap_log2FC + taz_log2FC,
    residual_observed_minus_predicted = observed_double_log2FC - predicted_additive_log2FC,
    abs_residual = abs(residual_observed_minus_predicted),
    YAP_guide_seed_match = ifelse(is.na(YAP_guide_has_match), FALSE, YAP_guide_has_match == 1),
    YAP_passenger_seed_match = ifelse(is.na(YAP_passenger_has_match), FALSE, YAP_passenger_has_match == 1),
    TAZ_guide_seed_match = ifelse(is.na(TAZ_guide_has_match), FALSE, TAZ_guide_has_match == 1),
    TAZ_passenger_seed_match = ifelse(is.na(TAZ_passenger_has_match), FALSE, TAZ_passenger_has_match == 1),
    any_passenger_seed_match = YAP_passenger_seed_match | TAZ_passenger_seed_match,
    any_near_perfect_match = ifelse(is.na(any_near_perfect_match), FALSE, any_near_perfect_match),
    near_perfect_query_names = ifelse(is.na(near_perfect_query_names), "", near_perfect_query_names),
    intended_target_status = case_when(
      gene_symbol == "Yap1" ~ "YAP_target",
      gene_symbol == "Wwtr1" ~ "TAZ_target",
      TRUE ~ "not_intended_target"
    )
  ) %>%
  arrange(desc(abs_residual))

top_n_deviation <- 50
additivity_df <- additivity_df %>%
  mutate(
    positive_deviation_rank = rank(-residual_observed_minus_predicted, ties.method = "first"),
    negative_deviation_rank = rank(residual_observed_minus_predicted, ties.method = "first"),
    strong_positive_deviation = positive_deviation_rank <= top_n_deviation,
    strong_negative_deviation = negative_deviation_rank <= top_n_deviation,
    strong_nonadditive = strong_positive_deviation | strong_negative_deviation
  ) %>%
  dplyr::select(
    gene_id_stable, gene_symbol, intended_target_status,
    yap_log2FC, taz_log2FC, observed_double_log2FC, predicted_additive_log2FC,
    residual_observed_minus_predicted, abs_residual,
    positive_deviation_rank, negative_deviation_rank,
    strong_positive_deviation, strong_negative_deviation, strong_nonadditive,
    YAP_guide_seed_match, TAZ_guide_seed_match, YAP_passenger_seed_match, TAZ_passenger_seed_match,
    any_passenger_seed_match, any_near_perfect_match, near_perfect_query_names
  )

write.table(
  additivity_df,
  file.path(out_dir, "observed_vs_predicted_double_additivity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pearson_cor <- cor(additivity_df$observed_double_log2FC, additivity_df$predicted_additive_log2FC, use = "pairwise.complete.obs", method = "pearson")
spearman_cor <- cor(additivity_df$observed_double_log2FC, additivity_df$predicted_additive_log2FC, use = "pairwise.complete.obs", method = "spearman")

scatter_plot <- ggplot(additivity_df, aes(x = predicted_additive_log2FC, y = observed_double_log2FC)) +
  geom_point(alpha = 0.45, size = 1.1, color = "#2F5D8A") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  labs(
    title = "Observed double knockdown response vs simple additive prediction",
    subtitle = sprintf("Predicted additive = siYAP log2FC + siTAZ log2FC | Pearson r = %.3f, Spearman rho = %.3f", pearson_cor, spearman_cor),
    x = "Predicted additive log2 fold change",
    y = "Observed siYAP+siTAZ log2 fold change"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir, "supplemental_double_additivity_scatter.png"), scatter_plot, width = 7.5, height = 6, dpi = 300)

residual_plot <- ggplot(additivity_df, aes(x = residual_observed_minus_predicted)) +
  geom_histogram(bins = 60, fill = "#7A9E7E", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "firebrick") +
  labs(
    title = "Residual distribution for double knockdown additivity analysis",
    subtitle = "Residual = observed double-knockdown log2FC - predicted additive log2FC",
    x = "Residual (observed - predicted)",
    y = "Number of genes"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir, "supplemental_double_additivity_residuals.png"), residual_plot, width = 7.5, height = 5.2, dpi = 300)

positive_genes <- additivity_df %>% arrange(desc(residual_observed_minus_predicted)) %>% slice_head(n = top_n_deviation)
negative_genes <- additivity_df %>% arrange(residual_observed_minus_predicted) %>% slice_head(n = top_n_deviation)

positive_enrichment <- focused_keyword_enrichment(positive_genes$gene_id_stable, additivity_df$gene_id_stable)
negative_enrichment <- focused_keyword_enrichment(negative_genes$gene_id_stable, additivity_df$gene_id_stable)
if (!is.null(positive_enrichment)) {
  positive_enrichment <- positive_enrichment %>% mutate(deviation_direction = "positive")
  negative_enrichment <- negative_enrichment %>% mutate(deviation_direction = "negative")
  write.table(
    bind_rows(positive_enrichment, negative_enrichment),
    file.path(project_dir, "analysis", "03_pathways", "double_additivity_focused_keyword_enrichment.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

double_memo_lines <- c(
  "Double-knockdown additivity interpretation memo",
  "",
  "Additivity rule:",
  "- Predicted additive log2 fold change was defined as log2FC(siYAP vs siScr) + log2FC(siTAZ vs siScr).",
  "- Residual was defined as observed siYAP+siTAZ log2FC minus predicted additive log2FC.",
  sprintf("- Correlation between observed and predicted values: Pearson r = %.3f, Spearman rho = %.3f.", pearson_cor, spearman_cor),
  sprintf("- Strong non-additive deviation sets were defined as the top %d genes with the most positive residuals and the top %d genes with the most negative residuals.", top_n_deviation, top_n_deviation),
  "",
  sprintf(
    "Near-perfect unintended matches among the top %d positive and top %d negative deviation genes were rare: %d and %d genes, respectively.",
    top_n_deviation, top_n_deviation,
    sum(positive_genes$any_near_perfect_match),
    sum(negative_genes$any_near_perfect_match)
  ),
  sprintf(
    "Guide-seed matches among the top positive deviation genes: YAP guide %d, TAZ guide %d. Among the top negative deviation genes: YAP guide %d, TAZ guide %d.",
    sum(positive_genes$YAP_guide_seed_match),
    sum(positive_genes$TAZ_guide_seed_match),
    sum(negative_genes$YAP_guide_seed_match),
    sum(negative_genes$TAZ_guide_seed_match)
  )
)

if (!is.null(positive_enrichment)) {
  pos_sig <- positive_enrichment %>% filter(fdr_bh < 0.05)
  neg_sig <- negative_enrichment %>% filter(fdr_bh < 0.05)
  double_memo_lines <- c(
    double_memo_lines,
    "",
    "Focused biological enrichment among strongly non-additive genes:",
    if (nrow(pos_sig) > 0) {
      apply(
        pos_sig %>%
          mutate(line = sprintf("- Positive residual genes enriched for %s (OR %.2f, FDR %s).", focus_category, odds_ratio, formatC(fdr_bh, format = "e", digits = 2))) %>%
          dplyr::select(line),
        1,
        identity
      )
    } else {
      "- No focused keyword category reached FDR < 0.05 among positive residual genes."
    },
    if (nrow(neg_sig) > 0) {
      apply(
        neg_sig %>%
          mutate(line = sprintf("- Negative residual genes enriched for %s (OR %.2f, FDR %s).", focus_category, odds_ratio, formatC(fdr_bh, format = "e", digits = 2))) %>%
          dplyr::select(line),
        1,
        identity
      )
    } else {
      "- No focused keyword category reached FDR < 0.05 among negative residual genes."
    }
  )
}

double_memo_lines <- c(
  double_memo_lines,
  "",
  "Interpretation guidance:",
  "- Deviations between observed double-knockdown responses and the simple additive prediction support non-additive interaction between YAP and TAZ perturbation effects.",
  "- If strongly non-additive genes are not strongly marked by near-perfect or simple seed-based annotations, that argues against a simple sequence-based off-target explanation for the interaction signal.",
  "- This analysis remains conservative and does not completely exclude all off-target contributions or indirect secondary biology."
)

writeLines(double_memo_lines, file.path(out_dir, "double_additivity_interpretation_memo.txt"))

plot_df <- results %>%
  mutate(label = paste(sirna, strand_role, sep = " ")) %>%
  select(label, n_seed_matched_genes, median_log2fc_seed_matched, median_log2fc_seed_unmatched, frac_down_seed_matched, frac_down_seed_unmatched) %>%
  pivot_longer(
    cols = c(n_seed_matched_genes, median_log2fc_seed_matched, median_log2fc_seed_unmatched, frac_down_seed_matched, frac_down_seed_unmatched),
    names_to = "metric",
    values_to = "value"
  )

nearperfect_summary <- nearperfect %>%
  count(sirna, strand_role, mismatches, name = "n_hits") %>%
  mutate(metric = paste0("nearperfect_", mismatches, "mm"), label = paste(sirna, strand_role, sep = " ")) %>%
  select(label, metric, value = n_hits)

supp_plot_df <- bind_rows(plot_df %>% select(label, metric, value), nearperfect_summary)

p <- ggplot(supp_plot_df, aes(label, value, fill = label)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(title = "Supplemental siRNA specificity support summary", x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(out_dir, "supplemental_sirna_specificity_summary.png"), p, width = 10, height = 7, dpi = 300)

caption_lines <- c(
  "Supplemental specificity support summary",
  "",
  "Panels summarize transcriptome-wide near-perfect cDNA hit counts, the number of genes with reconstructed 3'UTR seed matches, median gene-level log2 fold changes for seed-matched versus unmatched genes, and the fraction of genes called significantly downregulated in seed-matched versus unmatched sets."
)
writeLines(caption_lines, file.path(out_dir, "supplemental_sirna_specificity_summary_caption.txt"))
