#!/usr/bin/env Rscript

# ============================================================
# Professional DESeq2 RNA-seq Analysis
# Arabidopsis thaliana - Heat Stress
#
# Design:
#   Normal:     3 biological replicates
#   HeatStress: 3 biological replicates
#
# Contrast:
#   HeatStress vs Normal
# ============================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

project_dir <- path.expand("~/Arabidopsis_HeatStress_RNAseq")

count_file <- file.path(
    project_dir,
    "results/counts/gene_counts.txt"
)

output_dir <- file.path(
    project_dir,
    "results/deseq2"
)

plot_dir <- file.path(
    output_dir,
    "plots"
)

summary_dir <- file.path(
    output_dir,
    "summary"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Sample metadata
# ------------------------------------------------------------

sample_names <- c(
    "Normal_1",
    "Normal_2",
    "Normal_3",
    "HeatStress_1",
    "HeatStress_2",
    "HeatStress_3"
)

condition <- factor(
    c(
        "Normal",
        "Normal",
        "Normal",
        "HeatStress",
        "HeatStress",
        "HeatStress"
    ),
    levels = c("Normal", "HeatStress")
)

coldata <- data.frame(
    condition = condition,
    row.names = sample_names
)

write.csv(
    coldata,
    file.path(output_dir, "sample_metadata.csv"),
    quote = FALSE
)

# ------------------------------------------------------------
# 3. Check input
# ------------------------------------------------------------

if (!file.exists(count_file)) {
    stop("ERROR: Count file does not exist:\n", count_file)
}

cat("Reading count matrix...\n")

# ------------------------------------------------------------
# 4. Read featureCounts output
# ------------------------------------------------------------

counts <- read.delim(
    count_file,
    comment.char = "#",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (ncol(counts) < 7) {
    stop("ERROR: featureCounts file does not contain expected columns.")
}

# featureCounts structure:
# Geneid Chr Start End Strand Length Sample1 Sample2 ...

count_matrix <- counts[, 7:ncol(counts), drop = FALSE]

# ------------------------------------------------------------
# 5. Force correct sample names
# ------------------------------------------------------------

if (ncol(count_matrix) != length(sample_names)) {
    stop(
        "ERROR: Expected 6 samples but found ",
        ncol(count_matrix)
    )
}

colnames(count_matrix) <- sample_names

rownames(count_matrix) <- counts$Geneid

count_matrix <- as.matrix(count_matrix)

storage.mode(count_matrix) <- "integer"

# ------------------------------------------------------------
# 6. Check for duplicated gene IDs
# ------------------------------------------------------------

if (anyDuplicated(rownames(count_matrix)) > 0) {

    cat("WARNING: duplicated Gene IDs detected.\n")
    cat("Keeping the first occurrence of each Gene ID.\n")

    count_matrix <- count_matrix[
        !duplicated(rownames(count_matrix)),
        ,
        drop = FALSE
    ]
}

# ------------------------------------------------------------
# 7. Initial statistics
# ------------------------------------------------------------

genes_initial <- nrow(count_matrix)

cat("Genes in original count matrix:",
    genes_initial, "\n")

# ------------------------------------------------------------
# 8. Remove genes with zero counts
# ------------------------------------------------------------

count_matrix <- count_matrix[
    rowSums(count_matrix) > 0,
    ,
    drop = FALSE
]

genes_after_zero_filter <- nrow(count_matrix)

cat("Genes after removing all-zero genes:",
    genes_after_zero_filter, "\n")

# ------------------------------------------------------------
# 9. Verify sample names
# ------------------------------------------------------------

if (!identical(
    colnames(count_matrix),
    rownames(coldata)
)) {
    stop(
        "ERROR: Sample names/order in count matrix and metadata do not match."
    )
}

# ------------------------------------------------------------
# 10. Create DESeq2 object
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ condition
)

# ------------------------------------------------------------
# 11. Low-count filtering
#
# Keep genes having at least 10 reads across all samples.
# ------------------------------------------------------------

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep, ]

genes_after_lowcount_filter <- nrow(dds)

cat("Genes after low-count filtering:",
    genes_after_lowcount_filter, "\n")

# ------------------------------------------------------------
# 12. Run DESeq2
# ------------------------------------------------------------

cat("\nRunning DESeq2...\n")

dds <- DESeq(dds)

cat("DESeq2 completed.\n")

# ------------------------------------------------------------
# 13. Normalized counts
# ------------------------------------------------------------

normalized_counts <- counts(
    dds,
    normalized = TRUE
)

normalized_counts_df <- as.data.frame(
    normalized_counts
)

normalized_counts_df$GeneID <- rownames(
    normalized_counts_df
)

normalized_counts_df <- normalized_counts_df[
    ,
    c("GeneID", sample_names)
]

write.csv(
    normalized_counts_df,
    file.path(
        output_dir,
        "normalized_counts.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 14. PCA analysis
# ------------------------------------------------------------

cat("Generating PCA plot...\n")

vsd <- vst(
    dds,
    blind = FALSE
)

pca_data <- plotPCA(
    vsd,
    intgroup = "condition",
    returnData = TRUE
)

percent_var <- round(
    100 * attr(
        pca_data,
        "percentVar"
    )
)

pca_plot <- ggplot(
    pca_data,
    aes(
        x = PC1,
        y = PC2,
        label = name,
        shape = condition
    )
) +
    geom_point(size = 4) +
    geom_text(
        vjust = -1,
        size = 3.5
    ) +
    xlab(
        paste0(
            "PC1: ",
            percent_var[1],
            "% variance"
        )
    ) +
    ylab(
        paste0(
            "PC2: ",
            percent_var[2],
            "% variance"
        )
    ) +
    ggtitle(
        "PCA - Arabidopsis Heat Stress RNA-seq"
    ) +
    theme_bw()

ggsave(
    file.path(
        plot_dir,
        "PCA_plot.pdf"
    ),
    pca_plot,
    width = 8,
    height = 6
)

# ------------------------------------------------------------
# 15. Sample correlation heatmap
# ------------------------------------------------------------

cat("Generating sample correlation heatmap...\n")

vsd_matrix <- assay(vsd)

cor_matrix <- cor(
    vsd_matrix,
    method = "pearson"
)

pdf(
    file.path(
        plot_dir,
        "sample_correlation_heatmap.pdf"
    ),
    width = 8,
    height = 7
)

pheatmap(
    cor_matrix,
    main = "Sample-to-Sample Correlation",
    display_numbers = TRUE,
    number_format = "%.2f"
)

dev.off()

# ------------------------------------------------------------
# 16. Dispersion plot
# ------------------------------------------------------------

cat("Generating dispersion plot...\n")

pdf(
    file.path(
        plot_dir,
        "dispersion_plot.pdf"
    ),
    width = 8,
    height = 6
)

plotDispEsts(dds)

dev.off()

# ------------------------------------------------------------
# 17. Differential expression
#
# HeatStress relative to Normal
# ------------------------------------------------------------

cat("Calculating HeatStress vs Normal...\n")

res <- results(
    dds,
    contrast = c(
        "condition",
        "HeatStress",
        "Normal"
    ),
    alpha = 0.05
)

# ------------------------------------------------------------
# 18. Log2 fold-change shrinkage
# ------------------------------------------------------------

cat("Applying log2 fold-change shrinkage...\n")

if (requireNamespace(
    "apeglm",
    quietly = TRUE
)) {

    res_shrunk <- lfcShrink(
        dds,
        coef = "condition_HeatStress_vs_Normal",
        type = "apeglm"
    )

} else {

    cat(
        "apeglm not installed; using normal shrinkage.\n"
    )

    res_shrunk <- lfcShrink(
        dds,
        coef = "condition_HeatStress_vs_Normal",
        type = "normal"
    )
}

# ------------------------------------------------------------
# 19. Prepare complete results
# ------------------------------------------------------------

res_df <- as.data.frame(res)

res_df$GeneID <- rownames(res_df)

res_df <- res_df[
    ,
    c(
        "GeneID",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj"
    )
]

res_df <- res_df[
    order(res_df$padj, na.last = TRUE),
    ,
    drop = FALSE
]

write.csv(
    res_df,
    file.path(
        output_dir,
        "DESeq2_HeatStress_vs_Normal_all.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 20. Shrunk results
# ------------------------------------------------------------

shrunk_df <- as.data.frame(
    res_shrunk
)

shrunk_df$GeneID <- rownames(
    shrunk_df
)

shrunk_df <- shrunk_df[
    ,
    c(
        "GeneID",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "pvalue",
        "padj"
    )
]

shrunk_df <- shrunk_df[
    order(
        shrunk_df$padj,
        na.last = TRUE
    ),
    ,
    drop = FALSE
]

write.csv(
    shrunk_df,
    file.path(
        output_dir,
        "DESeq2_HeatStress_vs_Normal_LFC_shrunk.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 21. Significant DEGs
#
# padj < 0.05
# |log2FC| >= 1
# ------------------------------------------------------------

sig <- shrunk_df[
    !is.na(shrunk_df$padj) &
    shrunk_df$padj < 0.05 &
    abs(shrunk_df$log2FoldChange) >= 1,
    ,
    drop = FALSE
]

write.csv(
    sig,
    file.path(
        output_dir,
        "DESeq2_significant_DEGs.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 22. Upregulated genes
# ------------------------------------------------------------

up <- sig[
    sig$log2FoldChange >= 1,
    ,
    drop = FALSE
]

write.csv(
    up,
    file.path(
        output_dir,
        "DESeq2_upregulated_DEGs.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 23. Downregulated genes
# ------------------------------------------------------------

down <- sig[
    sig$log2FoldChange <= -1,
    ,
    drop = FALSE
]

write.csv(
    down,
    file.path(
        output_dir,
        "DESeq2_downregulated_DEGs.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 24. MA plot
# ------------------------------------------------------------

cat("Generating MA plot...\n")

pdf(
    file.path(
        plot_dir,
        "MA_plot.pdf"
    ),
    width = 8,
    height = 6
)

plotMA(
    res,
    alpha = 0.05,
    ylim = c(-6, 6),
    main = "MA Plot - HeatStress vs Normal"
)

dev.off()

# ------------------------------------------------------------
# 25. Volcano plot
# ------------------------------------------------------------

cat("Generating volcano plot...\n")

volcano_df <- shrunk_df

volcano_df$negLog10Padj <- -log10(
    volcano_df$padj
)

volcano_df$Status <- "Not significant"

volcano_df$Status[
    !is.na(volcano_df$padj) &
    volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange >= 1
] <- "Upregulated"

volcano_df$Status[
    !is.na(volcano_df$padj) &
    volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange <= -1
] <- "Downregulated"

volcano_df <- volcano_df[
    is.finite(volcano_df$negLog10Padj),
    ,
    drop = FALSE
]

volcano_plot <- ggplot(
    volcano_df,
    aes(
        x = log2FoldChange,
        y = negLog10Padj,
        shape = Status
    )
) +
    geom_point(
        alpha = 0.6,
        size = 1.8
    ) +
    geom_vline(
        xintercept = c(-1, 1),
        linetype = "dashed"
    ) +
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed"
    ) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-value") +
    ggtitle(
        "Volcano Plot - HeatStress vs Normal"
    ) +
    theme_bw()

ggsave(
    file.path(
        plot_dir,
        "volcano_plot.pdf"
    ),
    volcano_plot,
    width = 8,
    height = 6
)

# ------------------------------------------------------------
# 26. Save DESeq2 object
# ------------------------------------------------------------

saveRDS(
    dds,
    file.path(
        output_dir,
        "DESeq2_dataset.rds"
    )
)

# ------------------------------------------------------------
# 27. Analysis summary
# ------------------------------------------------------------

summary_file <- file.path(
    summary_dir,
    "DESeq2_summary.txt"
)

sink(summary_file)

cat("============================================\n")
cat("DESeq2 ANALYSIS SUMMARY\n")
cat("============================================\n\n")

cat("Organism: Arabidopsis thaliana\n")
cat("Experiment: Heat Stress\n\n")

cat("Samples:\n")
cat("Normal: 3 biological replicates\n")
cat("HeatStress: 3 biological replicates\n\n")

cat("Reference condition: Normal\n")
cat("Contrast: HeatStress vs Normal\n\n")

cat("Genes in original count matrix:",
    genes_initial, "\n")

cat("Genes after removing zero-count genes:",
    genes_after_zero_filter, "\n")

cat("Genes after low-count filtering:",
    genes_after_lowcount_filter, "\n\n")

cat("Significant DEG criteria:\n")
cat("Adjusted p-value < 0.05\n")
cat("|log2FC| >= 1\n\n")

cat("Total significant DEGs:",
    nrow(sig), "\n")

cat("Upregulated genes:",
    nrow(up), "\n")

cat("Downregulated genes:",
    nrow(down), "\n\n")

cat("DESeq2 version:",
    as.character(
        packageVersion("DESeq2")
    ),
    "\n")

cat("R version:",
    R.version.string,
    "\n")

cat("\n============================================\n")

sink()

# ------------------------------------------------------------
# 28. Finished
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("DESeq2 ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat("Genes analyzed:",
    genes_after_lowcount_filter,
    "\n")

cat("Significant DEGs:",
    nrow(sig),
    "\n")

cat("Upregulated:",
    nrow(up),
    "\n")

cat("Downregulated:",
    nrow(down),
    "\n\n")

cat("Results:",
    output_dir,
    "\n")

cat("Plots:",
    plot_dir,
    "\n")

cat("Summary:",
    summary_file,
    "\n")

cat("============================================\n")
