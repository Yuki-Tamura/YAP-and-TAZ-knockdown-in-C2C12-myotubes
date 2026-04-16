suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

safe_require <- function(pkg) {
  suppressWarnings(requireNamespace(pkg, quietly = TRUE))
}

write_note <- function(path, lines) {
  ensure_dir(dirname(path))
  writeLines(lines, path)
}

find_count_matrix_file <- function(project_dir) {
  candidates <- c(
    list.files(project_dir, pattern = "\\.xlsx$", full.names = TRUE),
    list.files(project_dir, pattern = "\\.csv$", full.names = TRUE)
  )
  if (length(candidates) == 0) {
    stop("No .xlsx or .csv count matrix file found in project directory.")
  }
  candidates[[1]]
}

load_counts_from_excel <- function(path) {
  dat <- readxl::read_excel(path, col_names = TRUE)
  colnames(dat)[1] <- "gene_id"
  counts <- as.data.frame(dat)
  counts$gene_id <- as.character(counts$gene_id)
  sample_cols <- setdiff(colnames(counts), "gene_id")
  for (nm in sample_cols) {
    counts[[nm]] <- as.integer(round(as.numeric(counts[[nm]])))
  }
  counts
}

load_counts_from_csv <- function(path) {
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(dat)[1] <- "gene_id"
  counts <- as.data.frame(dat)
  counts$gene_id <- as.character(counts$gene_id)
  sample_cols <- setdiff(colnames(counts), "gene_id")
  for (nm in sample_cols) {
    counts[[nm]] <- as.integer(round(as.numeric(counts[[nm]])))
  }
  counts
}

load_counts_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    return(load_counts_from_excel(path))
  }
  if (ext == "csv") {
    return(load_counts_from_csv(path))
  }
  stop(sprintf("Unsupported count matrix format: %s", path))
}

parse_metadata <- function(sample_names) {
  condition <- sub("_(\\d+)$", "", sample_names)
  replicate <- as.integer(sub("^.*_(\\d+)$", "\\1", sample_names))
  tibble(
    sample = sample_names,
    condition = condition,
    replicate = replicate,
    YAP_KD = ifelse(grepl("YAP", condition, fixed = TRUE), "yes", "no"),
    TAZ_KD = ifelse(grepl("TAZ", condition, fixed = TRUE), "yes", "no")
  ) %>%
    mutate(
      condition = factor(condition, levels = c("siScr", "siYAP", "siTAZ", "siYAP+siTAZ")),
      YAP_KD = factor(YAP_KD, levels = c("no", "yes")),
      TAZ_KD = factor(TAZ_KD, levels = c("no", "yes"))
    )
}

qc_library_sizes <- function(count_matrix, metadata, out_dir) {
  lib_sizes <- tibble(
    sample = colnames(count_matrix),
    library_size = colSums(count_matrix),
    detected_genes = colSums(count_matrix > 0)
  ) %>%
    left_join(metadata, by = "sample")

  write.csv(lib_sizes, file.path(out_dir, "library_size_summary.csv"), row.names = FALSE)

  p <- ggplot(lib_sizes, aes(x = sample, y = library_size, fill = condition)) +
    geom_col(width = 0.7) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Library size", title = "RNA-seq library sizes") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "library_sizes_barplot.png"), p, width = 8, height = 4, dpi = 300)

  p2 <- ggplot(lib_sizes, aes(x = sample, y = detected_genes, fill = condition)) +
    geom_col(width = 0.7) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = "Detected genes (>0 counts)", title = "Detected genes per sample") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "detected_genes_barplot.png"), p2, width = 8, height = 4, dpi = 300)

  lib_sizes
}

vst_matrix_from_counts <- function(count_matrix) {
  log2(count_matrix + 1)
}

qc_distance_heatmap <- function(vst_mat, metadata, out_dir) {
  d <- dist(t(vst_mat))
  dist_mat <- as.matrix(d)
  write.csv(dist_mat, file.path(out_dir, "sample_distance_matrix.csv"))

  if (safe_require("pheatmap")) {
    annotation_df <- metadata %>% dplyr::select(sample, condition, YAP_KD, TAZ_KD) %>% as.data.frame()
    rownames(annotation_df) <- annotation_df$sample
    annotation_df$sample <- NULL
    pheatmap::pheatmap(
      dist_mat,
      annotation_col = annotation_df,
      annotation_row = annotation_df,
      filename = file.path(out_dir, "sample_distance_heatmap.png"),
      width = 8,
      height = 7
    )
  } else {
    dist_long <- as.data.frame(as.table(dist_mat))
    p <- ggplot(dist_long, aes(Var1, Var2, fill = Freq)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "firebrick") +
      labs(x = NULL, y = NULL, fill = "Distance", title = "Sample-to-sample distance") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out_dir, "sample_distance_heatmap.png"), p, width = 8, height = 7, dpi = 300)
  }

  dist_mat
}

qc_pca <- function(vst_mat, metadata, out_dir) {
  pca <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
  var_explained <- summary(pca)$importance[2, 1:4]
  pca_df <- as_tibble(pca$x[, 1:4], rownames = "sample") %>%
    left_join(metadata, by = "sample")

  write.csv(pca_df, file.path(out_dir, "pca_coordinates.csv"), row.names = FALSE)
  write.csv(
    data.frame(component = names(var_explained), proportion_variance = as.numeric(var_explained)),
    file.path(out_dir, "pca_variance_explained.csv"),
    row.names = FALSE
  )

  p <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.7, size = 3) +
    labs(
      title = "PCA of log2(count + 1)",
      x = sprintf("PC1 (%.1f%%)", 100 * var_explained[[1]]),
      y = sprintf("PC2 (%.1f%%)", 100 * var_explained[[2]])
    ) +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "pca_plot.png"), p, width = 7, height = 5, dpi = 300)

  pca_df
}

run_deseq2 <- function(count_df, metadata, out_dir, notes_dir) {
  if (!safe_require("DESeq2")) {
    write_note(
      file.path(notes_dir, "software_limitations.txt"),
      c(
        "DESeq2 could not be loaded locally, so differential expression analysis could not be executed.",
        "Install Bioconductor packages including DESeq2 before rerunning analysis/run_workflow.R."
      )
    )
    return(NULL)
  }

  tryCatch({
    counts <- as.matrix(count_df[, metadata$sample, drop = FALSE])
    rownames(counts) <- count_df$gene_id

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = as.data.frame(metadata),
      design = ~ YAP_KD + TAZ_KD + YAP_KD:TAZ_KD
    )
    keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
    dds <- dds[keep, ]
    dds <- DESeq2::DESeq(dds)

    norm_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
    norm_counts <- tibble::rownames_to_column(norm_counts, "gene_id")
    write.csv(norm_counts, file.path(out_dir, "normalized_counts.csv"), row.names = FALSE)

    vst <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
    vst_mat <- SummarizedExperiment::assay(vst)
    write.csv(
      tibble::rownames_to_column(as.data.frame(vst_mat), "gene_id"),
      file.path(out_dir, "vst_matrix.csv"),
      row.names = FALSE
    )

    results_list <- list(
      siYAP_vs_siScr = DESeq2::results(dds, contrast = c("YAP_KD", "yes", "no")),
      siTAZ_vs_siScr = DESeq2::results(dds, contrast = c("TAZ_KD", "yes", "no")),
      siYAPsiTAZ_vs_siScr = DESeq2::results(dds, contrast = list(c("YAP_KD_yes_vs_no", "TAZ_KD_yes_vs_no", "YAP_KDyes.TAZ_KDyes"))),
      interaction_effect = DESeq2::results(dds, name = "YAP_KDyes.TAZ_KDyes")
    )

    exported <- list()
    for (nm in names(results_list)) {
      res_df <- as.data.frame(results_list[[nm]]) %>%
        tibble::rownames_to_column("gene_id") %>%
        arrange(padj, pvalue)
      res_df$contrast <- nm
      exported[[nm]] <- res_df
      write.csv(res_df, file.path(out_dir, paste0(nm, "_all_genes.csv")), row.names = FALSE)

      sig_df <- res_df %>%
        filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
      write.csv(sig_df, file.path(out_dir, paste0(nm, "_significant_degs.csv")), row.names = FALSE)
    }

    saveRDS(dds, file.path(out_dir, "dds.rds"))
    saveRDS(vst, file.path(out_dir, "vst.rds"))
    write_note(file.path(notes_dir, "software_limitations.txt"), "DESeq2 loaded and differential expression analysis completed successfully in this run.")

    list(dds = dds, vst = vst, results = exported)
  }, error = function(e) {
    write_note(
      file.path(notes_dir, "software_limitations.txt"),
      c(
        "DESeq2 analysis failed after package load.",
        paste("Error:", conditionMessage(e))
      )
    )
    NULL
  })
}

make_volcano <- function(res_df, title, out_path) {
  plot_df <- res_df %>%
    mutate(
      neglog10padj = -log10(ifelse(is.na(padj) | padj <= 0, 1, padj)),
      significant = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1
    )

  p <- ggplot(plot_df, aes(log2FoldChange, neglog10padj, color = significant)) +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(
      title = title,
      x = "log2 fold change",
      y = "-log10 adjusted P value",
      color = "DEG"
    ) +
    theme_bw(base_size = 11)
  ggsave(out_path, p, width = 6, height = 5, dpi = 300)
}

make_top_deg_tables <- function(results_list, out_dir) {
  top_tables <- lapply(names(results_list), function(nm) {
    results_list[[nm]] %>%
      filter(!is.na(padj)) %>%
      arrange(padj, desc(abs(log2FoldChange))) %>%
      slice_head(n = 25) %>%
      mutate(contrast = nm)
  })
  top_df <- bind_rows(top_tables)
  write.csv(top_df, file.path(out_dir, "top_25_degs_per_contrast.csv"), row.names = FALSE)
  top_df
}

deg_overlap_and_similarity <- function(results_list, out_fig_dir, out_table_dir) {
  yap <- results_list$siYAP_vs_siScr
  taz <- results_list$siTAZ_vs_siScr
  combo <- results_list$siYAPsiTAZ_vs_siScr

  deg_sets <- list(
    siYAP = yap %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>% pull(gene_id),
    siTAZ = taz %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>% pull(gene_id),
    siYAPsiTAZ = combo %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>% pull(gene_id)
  )

  overlap_summary <- tibble(
    comparison = c("siYAP_only", "siTAZ_only", "shared_siYAP_siTAZ", "siYAPsiTAZ"),
    n_genes = c(
      length(setdiff(deg_sets$siYAP, deg_sets$siTAZ)),
      length(setdiff(deg_sets$siTAZ, deg_sets$siYAP)),
      length(intersect(deg_sets$siYAP, deg_sets$siTAZ)),
      length(deg_sets$siYAPsiTAZ)
    )
  )
  write.csv(overlap_summary, file.path(out_table_dir, "deg_overlap_summary.csv"), row.names = FALSE)

  corr_df <- yap %>%
    dplyr::select(gene_id, log2FoldChange_yap = log2FoldChange, padj_yap = padj) %>%
    inner_join(taz %>% dplyr::select(gene_id, log2FoldChange_taz = log2FoldChange, padj_taz = padj), by = "gene_id")

  cor_est <- suppressWarnings(cor(corr_df$log2FoldChange_yap, corr_df$log2FoldChange_taz, use = "pairwise.complete.obs"))
  write.csv(corr_df, file.path(out_table_dir, "siYAP_vs_siTAZ_log2fc_source_data.csv"), row.names = FALSE)

  p <- ggplot(corr_df, aes(log2FoldChange_yap, log2FoldChange_taz)) +
    geom_point(alpha = 0.45, size = 1.1, color = "steelblue4") +
    geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
    labs(
      title = sprintf("siYAP vs siTAZ log2FC correlation (r = %.3f)", cor_est),
      x = "siYAP vs siScr log2FC",
      y = "siTAZ vs siScr log2FC"
    ) +
    theme_bw(base_size = 11)
  ggsave(file.path(out_fig_dir, "siYAP_vs_siTAZ_log2fc_correlation.png"), p, width = 6, height = 5, dpi = 300)

  overlap_plot_df <- tibble(
    category = factor(overlap_summary$comparison, levels = overlap_summary$comparison),
    n_genes = overlap_summary$n_genes
  )
  p2 <- ggplot(overlap_plot_df, aes(category, n_genes, fill = category)) +
    geom_col(show.legend = FALSE) +
    labs(title = "DEG overlap summary", x = NULL, y = "Number of genes") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(out_fig_dir, "deg_overlap_summary.png"), p2, width = 6, height = 4, dpi = 300)

  list(correlation = cor_est, overlap_summary = overlap_summary)
}

run_pathway_analysis <- function(results_list, out_dir, notes_dir) {
  required <- c("org.Mm.eg.db", "GO.db", "AnnotationDbi")
  missing <- required[!vapply(required, safe_require, logical(1))]
  if (length(missing) > 0) {
    write_note(
      file.path(notes_dir, "pathway_analysis_status.txt"),
      c(
        sprintf("Pathway analysis not run because these packages are unavailable locally: %s", paste(missing, collapse = ", ")),
        "Scripts and output locations are prepared; rerun after installing the missing packages."
      )
    )
    return(NULL)
  }

  tryCatch({
    suppressPackageStartupMessages(library(AnnotationDbi))
    universe <- unique(unlist(lapply(results_list, function(df) df$gene_id)))
    gene2go <- AnnotationDbi::select(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = universe,
      columns = c("GOALL", "ONTOLOGYALL", "SYMBOL"),
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
      dplyr::select(ENSEMBL, SYMBOL, GOALL) %>%
      distinct() %>%
      inner_join(go_terms, by = "GOALL")

    pathways_written <- list()
    focus_pattern <- "(muscle|contract|calcium|mitochond|metabo|redox|oxid|energy|respirat|sarcom|myofibr|ATP)"

    for (nm in names(results_list)) {
      res_df <- results_list[[nm]] %>%
        filter(!is.na(padj))
      sig_genes <- unique(res_df %>%
                            filter(padj < 0.05, abs(log2FoldChange) >= 1) %>%
                            pull(gene_id))
      bg_genes <- unique(res_df$gene_id)

      term_stats <- gene2go %>%
        filter(ENSEMBL %in% bg_genes) %>%
        distinct(ENSEMBL, GOALL, TERM) %>%
        group_by(GOALL, TERM) %>%
        summarise(
          term_size = n_distinct(ENSEMBL),
          overlap = n_distinct(ENSEMBL[ENSEMBL %in% sig_genes]),
          overlap_genes = paste(sort(unique(ENSEMBL[ENSEMBL %in% sig_genes])), collapse = ";"),
          .groups = "drop"
        ) %>%
        filter(term_size >= 10, term_size <= 500, overlap >= 2) %>%
        rowwise() %>%
        mutate(
          pvalue = fisher.test(matrix(
            c(
              overlap,
              length(sig_genes) - overlap,
              term_size - overlap,
              length(bg_genes) - term_size - (length(sig_genes) - overlap)
            ),
            nrow = 2
          ), alternative = "greater")$p.value
        ) %>%
        ungroup() %>%
        mutate(
          padj = p.adjust(pvalue, method = "BH"),
          contrast = nm
        ) %>%
        arrange(padj, pvalue, desc(overlap))

      write.csv(term_stats, file.path(out_dir, paste0(nm, "_go_bp.csv")), row.names = FALSE)

      focused <- term_stats %>%
        filter(grepl(focus_pattern, TERM, ignore.case = TRUE))
      write.csv(focused, file.path(out_dir, paste0(nm, "_go_bp_focused.csv")), row.names = FALSE)

      pathways_written[[nm]] <- term_stats
    }

    write_note(
      file.path(notes_dir, "pathway_analysis_status.txt"),
      "GO biological process enrichment completed using a direct Fisher-test over-representation analysis from org.Mm.eg.db and GO.db."
    )
    pathways_written
  }, error = function(e) {
    write_note(
      file.path(notes_dir, "pathway_analysis_status.txt"),
      c(
        "Pathway analysis failed.",
        paste("Error:", conditionMessage(e))
      )
    )
    NULL
  })
}

write_supplemental_captions <- function(out_dir) {
  write_note(
    file.path(out_dir, "supplemental_figure_captions.txt"),
    c(
      "Supplemental Figure S1. RNA-seq quality control in differentiated C2C12 myotubes after siScr, siYAP, siTAZ, or combined siYAP+siTAZ treatment. Panels show library sizes, detected genes per sample, principal component analysis of log2(count+1) values, and the sample-to-sample distance heatmap.",
      "Supplemental Figure S2. Supportive siRNA specificity analysis. Panels summarize inferred guide and passenger strands, availability of local annotation resources for sequence-based off-target assessment, transcriptome-level similarity between siYAP and siTAZ perturbations, and DEG overlap summaries."
    )
  )
}

write_methods_memo <- function(project_dir, analysis_dir, count_file, metadata, de_ran, pathway_ran, sirna_ran) {
  lines <- c(
    "Methods memo",
    "",
    sprintf("Project directory: %s", project_dir),
    sprintf("Input count matrix: %s", basename(count_file)),
    "",
    "Sample metadata derivation:",
    "- Sample names were parsed from the Excel header row.",
    "- Condition labels were inferred from prefixes: siScr, siYAP, siTAZ, and siYAP+siTAZ.",
    "- Binary factors YAP_KD and TAZ_KD were derived from the presence of 'YAP' and 'TAZ' in each condition label.",
    "",
    "QC:",
    "- Library sizes and detected-gene counts were summarized directly from the count matrix.",
    "- PCA and sample-to-sample distances were computed on log2(count + 1) values, or on DESeq2 VST values when DESeq2 completed.",
    "",
    "Differential expression:",
    "- The intended design was ~ YAP_KD + TAZ_KD + YAP_KD:TAZ_KD.",
    "- Requested simple contrasts were siYAP vs siScr, siTAZ vs siScr, and siYAP+siTAZ vs siScr.",
    sprintf("- Differential expression completed: %s", ifelse(de_ran, "yes", "no")),
    "",
    "Pathway analysis:",
    "- GO biological process enrichment was prepared for significant DEGs from each contrast.",
    sprintf("- Pathway analysis completed: %s", ifelse(pathway_ran, "yes", "no")),
    "",
    "siRNA specificity support:",
    "- Sense strands provided by the user were used to derive candidate passenger and reverse-complement guide strands with 3' overhang handling documented in analysis outputs.",
    "- Local transcriptome and 3'UTR resources were checked before attempting sequence-based off-target searches.",
    sprintf("- siRNA specificity support completed: %s", ifelse(sirna_ran, "yes", "partially; placeholder scripts and requirements were generated")),
    "",
    "Metadata table columns:",
    paste(colnames(metadata), collapse = ", ")
  )
  write_note(file.path(analysis_dir, "07_methods_memo", "methods_memo.txt"), lines)
}

write_summary <- function(analysis_dir, de_outputs, similarity_stats) {
  summary_path <- file.path(analysis_dir, "07_methods_memo", "key_findings_summary.txt")
  lines <- c("Key findings summary", "")

  if (is.null(de_outputs)) {
    lines <- c(lines, "DESeq2 was not available locally, so only workflow setup, QC summaries, and siRNA-support scaffolding were completed in this run.")
  } else {
    get_nsig <- function(name) {
      sum(!is.na(de_outputs$results[[name]]$padj) &
            de_outputs$results[[name]]$padj < 0.05 &
            abs(de_outputs$results[[name]]$log2FoldChange) >= 1)
    }
    lines <- c(
      lines,
      sprintf("Significant DEGs (padj < 0.05 and |log2FC| >= 1): siYAP %d, siTAZ %d, siYAP+siTAZ %d, interaction %d.",
              get_nsig("siYAP_vs_siScr"),
              get_nsig("siTAZ_vs_siScr"),
              get_nsig("siYAPsiTAZ_vs_siScr"),
              get_nsig("interaction_effect"))
    )
    if (!is.null(similarity_stats)) {
      lines <- c(
        lines,
        sprintf("Transcriptome-wide siYAP vs siTAZ log2FC correlation: %.3f.", similarity_stats$correlation),
        sprintf("Shared siYAP/siTAZ DEG count: %d.", similarity_stats$overlap_summary$n_genes[similarity_stats$overlap_summary$comparison == "shared_siYAP_siTAZ"]),
        "Interpretation should remain conservative: similarity or overlap alone does not prove on-target equivalence or exclude off-target effects."
      )
    }
  }

  lines <- c(
    lines,
    "",
    "Off-target interpretation framing:",
    "- This workflow provides supportive evidence against dominant shared sequence-based off-target suppression only if siYAP and siTAZ transcriptomic responses are not excessively similar.",
    "- Local annotation resources needed for direct transcriptome/3'UTR off-target searches were absent in this repository at run time."
  )

  write_note(summary_path, lines)
}

run_sirna_specificity <- function(analysis_dir, results_list = NULL, code_dir = analysis_dir, project_dir = dirname(code_dir)) {
  spec_dir <- ensure_dir(file.path(analysis_dir, "04_sirna_specificity"))
  notes_dir <- ensure_dir(file.path(analysis_dir, "00_notes"))
  figure_dir <- ensure_dir(file.path(analysis_dir, "05_figures"))
  table_dir <- ensure_dir(file.path(analysis_dir, "06_tables"))

  script <- file.path(code_dir, "python", "sirna_specificity.py")
  cmd <- sprintf("python %s %s %s", shQuote(script), shQuote(spec_dir), shQuote(project_dir))
  status <- system(cmd)
  if (status != 0) {
    write_note(
      file.path(notes_dir, "sirna_specificity_status.txt"),
      c("siRNA specificity helper script failed to run.", sprintf("Command attempted: %s", cmd))
    )
  }

  similarity_stats <- NULL
  if (!is.null(results_list)) {
    similarity_stats <- deg_overlap_and_similarity(results_list, figure_dir, table_dir)
  }

  if (file.exists(file.path(spec_dir, "sirna_specificity_summary.csv"))) {
    spec_df <- read.csv(file.path(spec_dir, "sirna_specificity_summary.csv"), stringsAsFactors = FALSE)
    p <- ggplot(spec_df, aes(sirna, n_distinct_seed_6mer_matches, fill = strand_role)) +
      geom_col(position = "dodge") +
      labs(
        title = "Seed-match support summary",
        x = NULL,
        y = "Distinct seed-match count (if 3'UTR FASTA available)"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(figure_dir, "sirna_specificity_support_summary.png"), p, width = 7, height = 4, dpi = 300)
  }

  similarity_stats
}

run_workflow <- function(project_dir, analysis_dir, count_file = NULL, code_dir = analysis_dir) {
  dirs <- c(
    "00_notes", "01_qc", "02_deseq2", "03_pathways",
    "04_sirna_specificity", "05_figures", "06_tables", "07_methods_memo"
  )
  invisible(lapply(file.path(analysis_dir, dirs), ensure_dir))

  if (is.null(count_file)) {
    count_file <- find_count_matrix_file(project_dir)
  }
  count_df <- load_counts_table(count_file)
  sample_names <- setdiff(colnames(count_df), "gene_id")
  metadata <- parse_metadata(sample_names)

  write.csv(metadata, file.path(analysis_dir, "06_tables", "sample_metadata.csv"), row.names = FALSE)
  write.csv(count_df, file.path(analysis_dir, "06_tables", "raw_counts_imported.csv"), row.names = FALSE)

  count_matrix <- as.matrix(count_df[, sample_names, drop = FALSE])
  rownames(count_matrix) <- count_df$gene_id

  qc_library_sizes(count_matrix, metadata, file.path(analysis_dir, "01_qc"))
  log_mat <- vst_matrix_from_counts(count_matrix)
  qc_distance_heatmap(log_mat, metadata, file.path(analysis_dir, "01_qc"))
  qc_pca(log_mat, metadata, file.path(analysis_dir, "01_qc"))

  de_outputs <- run_deseq2(
    count_df = count_df,
    metadata = metadata,
    out_dir = file.path(analysis_dir, "02_deseq2"),
    notes_dir = file.path(analysis_dir, "00_notes")
  )

  if (!is.null(de_outputs)) {
    qc_distance_heatmap(SummarizedExperiment::assay(de_outputs$vst), metadata, file.path(analysis_dir, "01_qc"))
    qc_pca(SummarizedExperiment::assay(de_outputs$vst), metadata, file.path(analysis_dir, "01_qc"))

    lapply(
      c("siYAP_vs_siScr", "siTAZ_vs_siScr", "siYAPsiTAZ_vs_siScr"),
      function(nm) {
        make_volcano(
          de_outputs$results[[nm]],
          title = nm,
          out_path = file.path(analysis_dir, "05_figures", paste0(nm, "_volcano.png"))
        )
      }
    )

    make_top_deg_tables(de_outputs$results, file.path(analysis_dir, "06_tables"))
    write.csv(
      de_outputs$results$interaction_effect,
      file.path(analysis_dir, "06_tables", "interaction_effect_genes.csv"),
      row.names = FALSE
    )
  }

  pathway_outputs <- NULL
  if (!is.null(de_outputs)) {
    pathway_outputs <- run_pathway_analysis(
      results_list = de_outputs$results,
      out_dir = file.path(analysis_dir, "03_pathways"),
      notes_dir = file.path(analysis_dir, "00_notes")
    )
  }

  similarity_stats <- run_sirna_specificity(
    analysis_dir = analysis_dir,
    results_list = if (is.null(de_outputs)) NULL else de_outputs$results,
    code_dir = code_dir,
    project_dir = project_dir
  )

  write_supplemental_captions(file.path(analysis_dir, "07_methods_memo"))
  write_methods_memo(
    project_dir = project_dir,
    analysis_dir = analysis_dir,
    count_file = count_file,
    metadata = metadata,
    de_ran = !is.null(de_outputs),
    pathway_ran = !is.null(pathway_outputs),
    sirna_ran = TRUE
  )
  write_summary(analysis_dir, de_outputs, similarity_stats)
}
