#!/usr/bin/env Rscript

# PRO-seq Relevant RNA-seq Analysis
# Comparing the actual plasmids used in PRO-seq experiments

# Load libraries
library(tximport)
library(DESeq2)
library(ggplot2)

# Create output directory
dir.create("results/proseq_relevant", recursive = TRUE, showWarnings = FALSE)

# Load data (assuming it's already loaded, if not reload)
if (!exists("txi")) {
  cat("Loading data...\n")
  design <- read.csv("experimental_design.csv")
  salmon_files <- list.files("results/salmon", pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
  names(salmon_files) <- basename(dirname(salmon_files))
  txi <- tximport(salmon_files, type = "salmon", txOut = TRUE)
  
  salmon_samples <- colnames(txi$counts)
  design_matched <- design[match(salmon_samples, design$sample), ]
  rownames(design_matched) <- salmon_samples
}

cat("=== PRO-SEQ RELEVANT RNA-SEQ ANALYSIS ===\n\n")

# Analysis 1: p65 (656) vs 2xTIMs (917) ----
cat("=== ANALYSIS 1: p65 vs 2xTIMs (the main comparison!) ===\n")

design_p65_vs_tims <- design_matched
design_p65_vs_tims$plasmid <- NA
design_p65_vs_tims$plasmid[design_p65_vs_tims$condition == "condition_656"] <- "p65"
design_p65_vs_tims$plasmid[design_p65_vs_tims$condition == "condition_917"] <- "tims"

# Keep only p65 and 2xTIMs samples
plasmid_samples <- !is.na(design_p65_vs_tims$plasmid)
design_plasmids <- design_p65_vs_tims[plasmid_samples,]
txi_plasmids <- list(
  counts = txi$counts[, plasmid_samples],
  abundance = txi$abundance[, plasmid_samples],
  length = txi$length[, plasmid_samples],
  countsFromAbundance = txi$countsFromAbundance
)

cat("Samples per plasmid:\n")
print(table(design_plasmids$plasmid))

# Run DESeq2
dds_plasmids <- DESeqDataSetFromTximport(txi_plasmids, design_plasmids, design = ~ plasmid)
keep_plasmids <- rowSums(counts(dds_plasmids)) >= 10
dds_plasmids <- dds_plasmids[keep_plasmids,]
dds_plasmids <- DESeq(dds_plasmids)

res_p65_vs_tims <- results(dds_plasmids, contrast = c("plasmid", "tims", "p65"))
cat("Results summary (2xTIMs vs p65):\n")
summary(res_p65_vs_tims)

# Save results
write.csv(as.data.frame(res_p65_vs_tims), "results/proseq_relevant/p65_vs_tims_results.csv")

# Create volcano plot
volcano_data_plasmids <- data.frame(
  log2FoldChange = res_p65_vs_tims$log2FoldChange,
  neglog10pval = -log10(res_p65_vs_tims$padj),
  significant = res_p65_vs_tims$padj < 0.05 & abs(res_p65_vs_tims$log2FoldChange) > 1
)
volcano_data_plasmids <- volcano_data_plasmids[complete.cases(volcano_data_plasmids),]

p_volcano_plasmids <- ggplot(volcano_data_plasmids, aes(x = log2FoldChange, y = neglog10pval, color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot: 2xTIMs vs p65",
       x = "Log2 Fold Change (2xTIMs vs p65)", 
       y = "-Log10(Adjusted P-value)") +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlim(-5, 5) +
  ylim(0, 20)

ggsave("results/proseq_relevant/volcano_p65_vs_tims.png", p_volcano_plasmids, width = 10, height = 6)

# Analysis 2: p65 vs GFP control ----
cat("\n=== ANALYSIS 2: p65 vs GFP control ===\n")

design_p65_vs_gfp <- design_matched
design_p65_vs_gfp$comparison <- NA
design_p65_vs_gfp$comparison[design_p65_vs_gfp$condition == "condition_656"] <- "p65"
design_p65_vs_gfp$comparison[design_p65_vs_gfp$condition == "GFP_control"] <- "GFP"

# Keep only p65 and GFP samples
p65_gfp_samples <- !is.na(design_p65_vs_gfp$comparison)
design_p65_gfp <- design_p65_vs_gfp[p65_gfp_samples,]
txi_p65_gfp <- list(
  counts = txi$counts[, p65_gfp_samples],
  abundance = txi$abundance[, p65_gfp_samples],
  length = txi$length[, p65_gfp_samples],
  countsFromAbundance = txi$countsFromAbundance
)

cat("Samples per group:\n")
print(table(design_p65_gfp$comparison))

# Run DESeq2
dds_p65_gfp <- DESeqDataSetFromTximport(txi_p65_gfp, design_p65_gfp, design = ~ comparison)
keep_p65_gfp <- rowSums(counts(dds_p65_gfp)) >= 10
dds_p65_gfp <- dds_p65_gfp[keep_p65_gfp,]
dds_p65_gfp <- DESeq(dds_p65_gfp)

res_p65_vs_gfp <- results(dds_p65_gfp, contrast = c("comparison", "p65", "GFP"))
cat("Results summary (p65 vs GFP):\n")
summary(res_p65_vs_gfp)

# Save results
write.csv(as.data.frame(res_p65_vs_gfp), "results/proseq_relevant/p65_vs_gfp_results.csv")

# Analysis 3: 2xTIMs vs GFP control ----
cat("\n=== ANALYSIS 3: 2xTIMs vs GFP control ===\n")

design_tims_vs_gfp <- design_matched
design_tims_vs_gfp$comparison <- NA
design_tims_vs_gfp$comparison[design_tims_vs_gfp$condition == "condition_917"] <- "tims"
design_tims_vs_gfp$comparison[design_tims_vs_gfp$condition == "GFP_control"] <- "GFP"

# Keep only 2xTIMs and GFP samples
tims_gfp_samples <- !is.na(design_tims_vs_gfp$comparison)
design_tims_gfp <- design_tims_vs_gfp[tims_gfp_samples,]
txi_tims_gfp <- list(
  counts = txi$counts[, tims_gfp_samples],
  abundance = txi$abundance[, tims_gfp_samples],
  length = txi$length[, tims_gfp_samples],
  countsFromAbundance = txi$countsFromAbundance
)

cat("Samples per group:\n")
print(table(design_tims_gfp$comparison))

# Run DESeq2
dds_tims_gfp <- DESeqDataSetFromTximport(txi_tims_gfp, design_tims_gfp, design = ~ comparison)
keep_tims_gfp <- rowSums(counts(dds_tims_gfp)) >= 10
dds_tims_gfp <- dds_tims_gfp[keep_tims_gfp,]
dds_tims_gfp <- DESeq(dds_tims_gfp)

res_tims_vs_gfp <- results(dds_tims_gfp, contrast = c("comparison", "tims", "GFP"))
cat("Results summary (2xTIMs vs GFP):\n")
summary(res_tims_vs_gfp)

# Save results
write.csv(as.data.frame(res_tims_vs_gfp), "results/proseq_relevant/tims_vs_gfp_results.csv")

# Create volcano plots for GFP comparisons
volcano_data_p65_gfp <- data.frame(
  log2FoldChange = res_p65_vs_gfp$log2FoldChange,
  neglog10pval = -log10(res_p65_vs_gfp$padj),
  significant = res_p65_vs_gfp$padj < 0.05 & abs(res_p65_vs_gfp$log2FoldChange) > 1
)
volcano_data_p65_gfp <- volcano_data_p65_gfp[complete.cases(volcano_data_p65_gfp),]

p_volcano_p65_gfp <- ggplot(volcano_data_p65_gfp, aes(x = log2FoldChange, y = neglog10pval, color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot: p65 vs GFP Control",
       x = "Log2 Fold Change (p65 vs GFP)", 
       y = "-Log10(Adjusted P-value)") +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlim(-5, 5) +
  ylim(0, 20)

ggsave("results/proseq_relevant/volcano_p65_vs_gfp.png", p_volcano_p65_gfp, width = 10, height = 6)

volcano_data_tims_gfp <- data.frame(
  log2FoldChange = res_tims_vs_gfp$log2FoldChange,
  neglog10pval = -log10(res_tims_vs_gfp$padj),
  significant = res_tims_vs_gfp$padj < 0.05 & abs(res_tims_vs_gfp$log2FoldChange) > 1
)
volcano_data_tims_gfp <- volcano_data_tims_gfp[complete.cases(volcano_data_tims_gfp),]

p_volcano_tims_gfp <- ggplot(volcano_data_tims_gfp, aes(x = log2FoldChange, y = neglog10pval, color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot: 2xTIMs vs GFP Control",
       x = "Log2 Fold Change (2xTIMs vs GFP)", 
       y = "-Log10(Adjusted P-value)") +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlim(-5, 5) +
  ylim(0, 20)

ggsave("results/proseq_relevant/volcano_tims_vs_gfp.png", p_volcano_tims_gfp, width = 10, height = 6)

# Summary report ----
cat("\n=== PRO-SEQ RELEVANT ANALYSIS SUMMARY ===\n")

sig_p65_vs_tims <- sum(res_p65_vs_tims$padj < 0.05, na.rm = TRUE)
sig_p65_vs_gfp <- sum(res_p65_vs_gfp$padj < 0.05, na.rm = TRUE)
sig_tims_vs_gfp <- sum(res_tims_vs_gfp$padj < 0.05, na.rm = TRUE)

cat("Significant transcripts (p < 0.05):\n")
cat("p65 vs 2xTIMs:", sig_p65_vs_tims, "\n")
cat("p65 vs GFP:", sig_p65_vs_gfp, "\n")
cat("2xTIMs vs GFP:", sig_tims_vs_gfp, "\n")

cat("\nFiles created in results/proseq_relevant/:\n")
cat("- p65_vs_tims_results.csv\n")
cat("- p65_vs_gfp_results.csv\n")
cat("- tims_vs_gfp_results.csv\n")
cat("- volcano_p65_vs_tims.png\n")
cat("- volcano_p65_vs_gfp.png\n")
cat("- volcano_tims_vs_gfp.png\n")

cat("\nThese comparisons are directly relevant to your PRO-seq experiments!\n")
cat("The p65 vs 2xTIMs comparison shows the direct transcriptional differences\n")
cat("between your two main experimental constructs.\n")