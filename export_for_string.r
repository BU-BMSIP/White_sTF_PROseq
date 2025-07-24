#!/usr/bin/env Rscript

# Export differentially expressed genes for STRING analysis
# Converts Ensembl transcript IDs to gene symbols and creates STRING-ready files

# Load required libraries
if (!require("biomaRt", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("biomaRt")
}

library(biomaRt)

# Create output directory
dir.create("results/string_analysis", recursive = TRUE, showWarnings = FALSE)

# Function to convert transcripts to genes and export for STRING
export_for_string <- function(results, comparison_name, padj_cutoff = 0.05, log2fc_cutoff = 1) {
  cat("Processing", comparison_name, "...\n")
  
  # Filter significant results
  sig_results <- results[!is.na(results$padj) & 
                           results$padj < padj_cutoff & 
                           abs(results$log2FoldChange) > log2fc_cutoff, ]
  
  if (nrow(sig_results) == 0) {
    cat("No significant transcripts found for", comparison_name, "\n")
    return(NULL)
  }
  
  cat("Found", nrow(sig_results), "significant transcripts\n")
  
  # Get transcript IDs (remove version numbers for biomaRt)
  transcript_ids <- gsub("\\.\\d+$", "", rownames(sig_results))
  
  # Connect to Ensembl
  cat("Connecting to Ensembl database...\n")
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Convert transcript IDs to gene symbols
  cat("Converting transcript IDs to gene symbols...\n")
  gene_info <- getBM(attributes = c("ensembl_transcript_id", "external_gene_name", "description", "gene_biotype"),
                     filters = "ensembl_transcript_id",
                     values = transcript_ids,
                     mart = mart)
  
  # Merge with results
  sig_results$transcript_id_clean <- gsub("\\.\\d+$", "", rownames(sig_results))
  merged_results <- merge(sig_results, gene_info, 
                          by.x = "transcript_id_clean", 
                          by.y = "ensembl_transcript_id", 
                          all.x = TRUE)
  
  # Create STRING-ready formats
  
  # 1. Simple gene list (just gene symbols, one per line)
  gene_symbols <- unique(merged_results$external_gene_name[merged_results$external_gene_name != ""])
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  
  # Save simple gene list
  simple_filename <- paste0("results/string_analysis/", gsub(" ", "_", comparison_name), "_gene_list.txt")
  writeLines(gene_symbols, simple_filename)
  cat("Saved simple gene list:", simple_filename, "(", length(gene_symbols), "genes )\n")
  
  # 2. Detailed CSV with fold changes and descriptions
  detailed_results <- merged_results[!is.na(merged_results$external_gene_name) & 
                                       merged_results$external_gene_name != "", ]
  detailed_results <- detailed_results[order(detailed_results$padj), ]
  
  # Select and rename columns for clarity
  export_detailed <- data.frame(
    Gene_Symbol = detailed_results$external_gene_name,
    Transcript_ID = rownames(detailed_results),
    Log2FoldChange = round(detailed_results$log2FoldChange, 3),
    AdjustedPvalue = detailed_results$padj,
    BaseMean = round(detailed_results$baseMean, 1),
    Direction = ifelse(detailed_results$log2FoldChange > 0, "UP", "DOWN"),
    Gene_Description = detailed_results$description,
    Gene_Biotype = detailed_results$gene_biotype
  )
  
  detailed_filename <- paste0("results/string_analysis/", gsub(" ", "_", comparison_name), "_detailed.csv")
  write.csv(export_detailed, detailed_filename, row.names = FALSE)
  cat("Saved detailed results:", detailed_filename, "(", nrow(export_detailed), "genes )\n")
  
  # 3. Separate files for up and down regulated genes
  up_genes <- gene_symbols[merged_results$external_gene_name %in% gene_symbols & 
                             merged_results$log2FoldChange > 0]
  down_genes <- gene_symbols[merged_results$external_gene_name %in% gene_symbols & 
                               merged_results$log2FoldChange < 0]
  
  if (length(up_genes) > 0) {
    up_filename <- paste0("results/string_analysis/", gsub(" ", "_", comparison_name), "_upregulated.txt")
    writeLines(unique(up_genes), up_filename)
    cat("Saved upregulated genes:", up_filename, "(", length(unique(up_genes)), "genes )\n")
  }
  
  if (length(down_genes) > 0) {
    down_filename <- paste0("results/string_analysis/", gsub(" ", "_", comparison_name), "_downregulated.txt")
    writeLines(unique(down_genes), down_filename)
    cat("Saved downregulated genes:", down_filename, "(", length(unique(down_genes)), "genes )\n")
  }
  
  # Summary
  cat("Summary for", comparison_name, ":\n")
  cat("- Total unique genes:", length(gene_symbols), "\n")
  cat("- Upregulated:", length(unique(up_genes)), "\n")
  cat("- Downregulated:", length(unique(down_genes)), "\n\n")
  
  return(list(
    genes = gene_symbols,
    detailed = export_detailed,
    up = unique(up_genes),
    down = unique(down_genes)
  ))
}

# Load results if they exist
cat("=== EXPORTING DIFFERENTIAL EXPRESSION RESULTS FOR STRING ===\n\n")

if (exists("res_p65")) {
  p65_export <- export_for_string(res_p65, "P65 Mutated vs Wildtype")
} else if (file.exists("results/r_analysis/p65_mutated_vs_wildtype_results.csv")) {
  cat("Loading P65 results from file...\n")
  res_p65 <- read.csv("results/r_analysis/p65_mutated_vs_wildtype_results.csv", row.names = 1)
  p65_export <- export_for_string(res_p65, "P65 Mutated vs Wildtype")
}

if (exists("res_tims")) {
  tims_export <- export_for_string(res_tims, "2xTIMs Mutated vs Wildtype")
} else if (file.exists("results/r_analysis/tims_mutated_vs_wildtype_results.csv")) {
  cat("Loading 2xTIMs results from file...\n")
  res_tims <- read.csv("results/r_analysis/tims_mutated_vs_wildtype_results.csv", row.names = 1)
  tims_export <- export_for_string(res_tims, "2xTIMs Mutated vs Wildtype")
}

if (exists("res_all")) {
  all_export <- export_for_string(res_all, "All Treatments vs Control")
} else if (file.exists("results/r_analysis/all_vs_control_results.csv")) {
  cat("Loading All vs Control results from file...\n")
  res_all <- read.csv("results/r_analysis/all_vs_control_results.csv", row.names = 1)
  all_export <- export_for_string(res_all, "All Treatments vs Control")
}

cat("=== EXPORT COMPLETE ===\n")
cat("Files saved in results/string_analysis/\n")
cat("\nFor STRING analysis:\n")
cat("1. Go to https://string-db.org\n")
cat("2. Click 'Search' -> 'Multiple proteins'\n")
cat("3. Upload one of the *_gene_list.txt files\n")
cat("4. Set organism to 'Homo sapiens'\n")
cat("5. Analyze networks and pathways!\n")

cat("\nFiles created:\n")
list.files("results/string_analysis/", full.names = FALSE)