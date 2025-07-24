#!/usr/bin/env Rscript

# DESeq2 analysis comparing Salmon and featureCounts quantification
# Usage: Called by Snakemake rule

library(DESeq2)
library(tximport)
library(readr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(htmlwidgets)
library(plotly)
library(DT)
library(rmarkdown)

# Get input files from snakemake object
salmon_files <- snakemake@input[["salmon_files"]]
featurecounts_file <- snakemake@input[["featurecounts"]]
design_file <- snakemake@input[["design"]]
output_file <- snakemake@output[[1]]

cat("Starting DESeq2 analysis...\n")
cat("Salmon files:", length(salmon_files), "\n")
cat("FeatureCounts file:", featurecounts_file, "\n")
cat("Design file:", design_file, "\n")

# Read experimental design
design <- read.csv(design_file, stringsAsFactors = FALSE)
cat("Design dimensions:", dim(design), "\n")
cat("Design columns:", colnames(design), "\n")

# Extract sample names from salmon file paths
sample_names <- basename(dirname(salmon_files))
names(salmon_files) <- sample_names

cat("Sample names from salmon files:\n")
print(sample_names)

# Match design to samples
design_matched <- design[match(sample_names, design$sample), ]
rownames(design_matched) <- sample_names

cat("Matched design:\n")
print(design_matched)

# Check if we have a condition column
if(!"condition" %in% colnames(design_matched)) {
    # Try to infer condition from sample names
    if(any(grepl("GFP", sample_names))) {
        design_matched$condition <- ifelse(grepl("GFP", sample_names), "GFP", "treatment")
    } else {
        # Create a simple condition based on sample groups
        design_matched$condition <- rep(c("condition_A", "condition_B"), length.out = nrow(design_matched))
    }
    cat("Created condition column:\n")
    print(table(design_matched$condition))
}

# Create output directory
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Create R Markdown report
rmd_content <- '
---
title: "RNA-seq Differential Expression Analysis"
output: 
  html_document:
    theme: flatly
    toc: true
    toc_float: true
    code_folding: hide
date: "`r Sys.Date()`"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.width = 10, fig.height = 6)
```

# Analysis Overview

This report compares RNA-seq quantification results from Salmon and featureCounts, followed by differential expression analysis using DESeq2.

## Sample Information

```{r sample_info}
DT::datatable(design_matched, options = list(pageLength = 20))
```

# Salmon Analysis

## Import Salmon Data

```{r salmon_import}
# For transcript-level analysis, we need a transcript-to-gene mapping
# Since we don'\''t have it, we'\''ll treat transcripts as genes for now
tx2gene <- data.frame(
  transcript_id = character(0),
  gene_id = character(0)
)

# Import salmon data at transcript level
txi <- tximport(salmon_files, type = "salmon", txOut = TRUE)

cat("Salmon import summary:\n")
cat("Samples:", ncol(txi$counts), "\n")
cat("Transcripts:", nrow(txi$counts), "\n")
cat("Total counts range:", range(colSums(txi$counts)), "\n")
```

## Salmon Quality Control

```{r salmon_qc}
# Sample correlation heatmap
cor_matrix <- cor(txi$counts + 1, method = "spearman")
pheatmap(cor_matrix, 
         annotation_row = design_matched["condition"],
         annotation_col = design_matched["condition"],
         main = "Sample Correlation (Salmon)")

# Library size distribution
lib_sizes <- data.frame(
  sample = colnames(txi$counts),
  library_size = colSums(txi$counts),
  condition = design_matched$condition
)

ggplot(lib_sizes, aes(x = sample, y = library_size, fill = condition)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Library Sizes (Salmon)", y = "Total Counts") +
  scale_y_continuous(labels = scales::comma)
```

## Salmon DESeq2 Analysis

```{r salmon_deseq2}
# Create DESeq2 dataset
dds_salmon <- DESeqDataSetFromTximport(txi, design_matched, design = ~ condition)

# Filter low count transcripts
keep <- rowSums(counts(dds_salmon)) >= 10
dds_salmon <- dds_salmon[keep,]

cat("Retained transcripts after filtering:", nrow(dds_salmon), "\n")

# Run DESeq2
dds_salmon <- DESeq(dds_salmon)

# Get results
res_salmon <- results(dds_salmon)
res_salmon_df <- as.data.frame(res_salmon)
res_salmon_df$transcript <- rownames(res_salmon_df)

cat("Salmon DESeq2 results:\n")
summary(res_salmon)
```

# FeatureCounts Analysis

## Import FeatureCounts Data

```{r featurecounts_import}
# Read featureCounts output
fc_data <- read.table(featurecounts_file, header = TRUE, sep = "\t", skip = 1)

# Extract count matrix
count_matrix <- as.matrix(fc_data[, 7:ncol(fc_data)])
rownames(count_matrix) <- fc_data$Geneid

# Clean column names to match sample names
colnames(count_matrix) <- gsub(".*\\/([^/]+)_Aligned\\.sortedByCoord\\.out\\.bam", "\\\\1", colnames(count_matrix))

cat("FeatureCounts import summary:\n")
cat("Samples:", ncol(count_matrix), "\n")
cat("Genes:", nrow(count_matrix), "\n")
cat("Total counts range:", range(colSums(count_matrix)), "\n")
```

## FeatureCounts Quality Control

```{r fc_qc}
# Sample correlation heatmap
fc_cor_matrix <- cor(count_matrix + 1, method = "spearman")
pheatmap(fc_cor_matrix, 
         annotation_row = design_matched["condition"],
         annotation_col = design_matched["condition"],
         main = "Sample Correlation (FeatureCounts)")

# Library size comparison
fc_lib_sizes <- data.frame(
  sample = colnames(count_matrix),
  library_size = colSums(count_matrix),
  condition = design_matched$condition
)

ggplot(fc_lib_sizes, aes(x = sample, y = library_size, fill = condition)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Library Sizes (FeatureCounts)", y = "Total Counts") +
  scale_y_continuous(labels = scales::comma)
```

## FeatureCounts DESeq2 Analysis

```{r fc_deseq2}
# Create DESeq2 dataset
dds_fc <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = design_matched,
                                design = ~ condition)

# Filter low count genes
keep_fc <- rowSums(counts(dds_fc)) >= 10
dds_fc <- dds_fc[keep_fc,]

cat("Retained genes after filtering:", nrow(dds_fc), "\n")

# Run DESeq2
dds_fc <- DESeq(dds_fc)

# Get results
res_fc <- results(dds_fc)
res_fc_df <- as.data.frame(res_fc)
res_fc_df$gene <- rownames(res_fc_df)

cat("FeatureCounts DESeq2 results:\n")
summary(res_fc)
```

# Comparison of Methods

## Volcano Plots

```{r volcano_plots}
# Salmon volcano plot
salmon_volcano <- res_salmon_df %>%
  filter(!is.na(padj)) %>%
  mutate(significant = padj < 0.05 & abs(log2FoldChange) > 1)

p1 <- ggplot(salmon_volcano, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Salmon - Volcano Plot") +
  scale_color_manual(values = c("grey", "red"))

# FeatureCounts volcano plot
fc_volcano <- res_fc_df %>%
  filter(!is.na(padj)) %>%
  mutate(significant = padj < 0.05 & abs(log2FoldChange) > 1)

p2 <- ggplot(fc_volcano, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "FeatureCounts - Volcano Plot") +
  scale_color_manual(values = c("grey", "red"))

print(p1)
print(p2)
```

## Summary Statistics

```{r summary_stats}
# Count significant results
salmon_sig <- sum(res_salmon_df$padj < 0.05 & abs(res_salmon_df$log2FoldChange) > 1, na.rm = TRUE)
fc_sig <- sum(res_fc_df$padj < 0.05 & abs(res_fc_df$log2FoldChange) > 1, na.rm = TRUE)

comparison_stats <- data.frame(
  Method = c("Salmon", "FeatureCounts"),
  Total_Features = c(nrow(res_salmon_df), nrow(res_fc_df)),
  Significant_DE = c(salmon_sig, fc_sig),
  Percent_DE = c(round(salmon_sig/nrow(res_salmon_df)*100, 2), 
                round(fc_sig/nrow(res_fc_df)*100, 2))
)

DT::datatable(comparison_stats, options = list(pageLength = 10))
```

# Top Results Tables

## Top Salmon Results

```{r top_salmon}
top_salmon <- res_salmon_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50)

DT::datatable(top_salmon, options = list(pageLength = 20, scrollX = TRUE)) %>%
  formatRound(columns = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), digits = 4)
```

## Top FeatureCounts Results

```{r top_fc}
top_fc <- res_fc_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50)

DT::datatable(top_fc, options = list(pageLength = 20, scrollX = TRUE)) %>%
  formatRound(columns = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), digits = 4)
```

# Session Information

```{r session_info}
sessionInfo()
```
'

# Write the R Markdown file
rmd_file <- sub("\\.html$", ".Rmd", output_file)
writeLines(rmd_content, rmd_file)

# Render the report
rmarkdown::render(rmd_file, output_file = output_file)

cat("Analysis complete! Report saved to:", output_file, "\n")