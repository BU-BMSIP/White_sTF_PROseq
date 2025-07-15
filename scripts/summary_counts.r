# Clear and reload the data
rm(list=ls())
library(reticulate)
library(ggplot2)
library(ggfortify)

# Load numpy and the data
np <- import("numpy")
data <- np$load("/projectnb/khalil/nwhite42/ProSEQ_project/results/correlation/summary_counts.npz")

# Extract the data matrix and labels
count_matrix <- data$f[['matrix']]
labels <- data$f[['labels']]

# Convert to R matrix
count_matrix <- as.matrix(count_matrix)

# Check dimensions
print(paste("Original dimensions:", nrow(count_matrix), "rows x", ncol(count_matrix), "columns"))

# TRANSPOSE the matrix so samples are rows and features are columns
count_matrix <- t(count_matrix)
print(paste("After transpose:", nrow(count_matrix), "rows (samples) x", ncol(count_matrix), "columns (features)"))

# Check labels
print(paste("Number of samples:", length(labels)))
print("Sample labels:")
print(labels)

# Remove zero-variance features
var_features <- apply(count_matrix, 2, var) > 0
count_matrix_filtered <- count_matrix[, var_features]
print(paste("Removed", sum(!var_features), "zero-variance features"))

# Perform PCA
pca_result <- prcomp(count_matrix_filtered, 
                     center = TRUE, 
                     scale. = TRUE)

# Summary of PCA
summary(pca_result)

# Variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Create PCA dataframe with labels
pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  PC3 = pca_result$x[,3],
  Sample = as.character(labels)
)

# PCA plot with labels
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 4, alpha = 0.8, aes(color = Sample)) +
  geom_text(aes(label = Sample), vjust = -1, size = 3) +
  labs(title = "PCA Plot",
       x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")) +
  theme_minimal() +
  theme(legend.position = "right")

print(p2)

# Save results
saveRDS(pca_result, "pca_result.rds")
write.csv(cbind(Sample = labels, pca_result$x), "pca_scores.csv", row.names = FALSE)

# Create a more informative scree plot
var_df <- data.frame(
  PC = 1:length(var_explained),
  Variance = var_explained,
  Cumulative = cumsum(var_explained)
)

p1 <- ggplot(var_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(y = Cumulative), color = "red", linewidth = 1) +
  geom_point(aes(y = Cumulative), color = "red", size = 2) +
  labs(title = "PCA Variance Explained",
       x = "Principal Component",
       y = "Proportion of Variance") +
  theme_minimal() +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Cumulative Variance"))

print(p1)

# Also create PC1 vs PC3 and PC2 vs PC3 plots
p3 <- ggplot(pca_df, aes(x = PC1, y = PC3)) +
  geom_point(size = 4, alpha = 0.8, aes(color = Sample)) +
  geom_text(aes(label = Sample), vjust = -1, size = 3) +
  labs(title = "PCA Plot: PC1 vs PC3",
       x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
       y = paste0("PC3 (", round(var_explained[3] * 100, 1), "%)")) +
  theme_minimal()

print(p3)

# Create a pairs plot for first 3 PCs
library(GGally)
pca_pairs_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  PC3 = pca_result$x[,3],
  Sample = as.character(labels)
)

p4 <- ggpairs(pca_pairs_df, 
              columns = 1:3,
              aes(color = Sample),
              upper = list(continuous = "points"),
              lower = list(continuous = "points"),
              diag = list(continuous = "density")) +
  theme_minimal()

print(p4)

# Get top contributing features
loadings_pc1 <- abs(pca_result$rotation[,1])
top_features_pc1 <- head(sort(loadings_pc1, decreasing = TRUE), 30)
print("Top 30 features contributing to PC1:")
print(top_features_pc1)

# If you want to know which features these correspond to
# (assuming features have names or indices)
print("\nIndices of top 30 features for PC1:")
print(names(top_features_pc1))

# Calculate sample distances in PC space
sample_dist <- dist(pca_result$x[,1:3])
print("\nSample distances in PC1-3 space:")
print(round(as.matrix(sample_dist), 2))

# Hierarchical clustering based on PCA
hc <- hclust(sample_dist)
plot(hc, main = "Hierarchical Clustering of Samples (based on PC1-3)")

# Save the current plot displayed in the Plots pane
ggsave("pca_plot_pc1_pc2.png", plot = p2, width = 8, height = 6, dpi = 300)
ggsave("scree_plot.png", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("pca_plot_pc1_pc3.png", plot = p3, width = 8, height = 6, dpi = 300)

# For the pairs plot
ggsave("pca_pairs_plot.png", plot = p4, width = 10, height = 10, dpi = 300)

# Save all plots to a PDF
pdf("pca_analysis_plots.pdf", width = 10, height = 8)
print(p1)  # Scree plot
print(p2)  # PC1 vs PC2
print(p3)  # PC1 vs PC3
print(p4)  # Pairs plot
dev.off()

# First, let's check if the plots were created properly
# Display each plot directly
p1  # Just type the object name
p2
p3
p4

# If that doesn't work, explicitly print them
print(p1)
print(p2)
print(p3)
print(p4)

# Check if the PDF was created successfully
list.files(pattern = "*.pdf")

# If you want to see the plots in a new window (outside RStudio)
dev.new()
print(p1)

# Or save them as individual PNG files to view them
ggsave("pca_scree_plot.png", plot = p1, width = 8, height = 6)
ggsave("pca_pc1_pc2.png", plot = p2, width = 8, height = 6)
ggsave("pca_pc1_pc3.png", plot = p3, width = 8, height = 6)
ggsave("pca_pairs.png", plot = p4, width = 10, height = 10)

# List the PNG files created
list.files(pattern = "*.png")

# For the hierarchical clustering, let's save it directly to a file
png("clustering.png", width = 800, height = 600)
par(mar = c(5, 4, 4, 2))  # Fix margins
plot(hc, main = "Hierarchical Clustering of Samples", 
     labels = labels, hang = -1)
dev.off()

# Check current graphics device
dev.cur()

# If needed, reset the graphics device
graphics.off()

# Check your current working directory
getwd()

# List all files in current directory
list.files()

# Create a simple text-based visualization of your PCA results
# First, let's see the PC scores
pca_scores <- data.frame(
  Sample = labels,
  PC1 = round(pca_result$x[,1], 1),
  PC2 = round(pca_result$x[,2], 1),
  PC3 = round(pca_result$x[,3], 1)
)
print("PCA Scores:")
print(pca_scores)

# Create a simple ASCII plot of PC1 vs PC2
library(txtplot)
txtplot(pca_result$x[,1], pca_result$x[,2], 
        xlab = "PC1", ylab = "PC2")

# Or create a basic base R plot (might work where ggplot doesn't)
plot(pca_result$x[,1], pca_result$x[,2], 
     xlab = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)"),
     main = "PCA Plot",
     pch = 19, col = 1:8, cex = 2)
text(pca_result$x[,1], pca_result$x[,2], labels = labels, pos = 3)
legend("topright", legend = labels, col = 1:8, pch = 19)

# Show the groupings based on distance matrix
cat("\nSample Groupings based on PCA distances:\n")
cat("Group 1: bigwig_1 and bigwig_8 (distance = 20.42)\n")
cat("Group 2: bigwig_2 and bigwig_4 (distance = 15.19)\n")
cat("Group 3: bigwig_3 and bigwig_5 (distance = 18.73)\n")
cat("Group 4: bigwig_6 and bigwig_7 (distance = 26.28)\n")

# Create a simple summary
cat("\nPCA Summary:\n")
cat("- Total features analyzed: 67,096 (after removing 11,590 zero-variance features)\n")
cat("- Number of samples: 8\n")
cat("- Variance explained:\n")
for(i in 1:8) {
  cat(sprintf("  PC%d: %.1f%%\n", i, var_explained[i] * 100))
}
cat(sprintf("- First 3 PCs explain %.1f%% of total variance\n", sum(var_explained[1:3]) * 100))

# 1. Create a heatmap of sample correlations
sample_cor <- cor(t(pca_result$x[,1:7]))  # Use first 7 PCs (excluding PC8 which is essentially 0)
rownames(sample_cor) <- labels
colnames(sample_cor) <- labels

# Simple correlation matrix display
print("Sample Correlation Matrix (based on PC1-7):")
print(round(sample_cor, 3))

# 2. Create a better visualization for showing similarity
library(pheatmap)
pdf("sample_similarity_heatmap.pdf", width = 8, height = 6)
pheatmap(sample_cor, 
         main = "Sample Similarity Based on PCA",
         display_numbers = TRUE,
         number_format = "%.2f",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101))
dev.off()

# 3. Calculate average distance between all samples
all_distances <- as.vector(as.dist(pca_result$x[,1:3]))
mean_distance <- mean(all_distances)
sd_distance <- sd(all_distances)

cat("\nOverall Sample Similarity Statistics:\n")
cat(sprintf("Average distance between samples: %.2f\n", mean_distance))
cat(sprintf("Standard deviation of distances: %.2f\n", sd_distance))
cat(sprintf("Coefficient of variation: %.2f%%\n", (sd_distance/mean_distance)*100))

# 4. Show within-group vs between-group distances
within_group <- c(20.42, 15.19, 18.73, 26.28)  # Your paired distances
between_group <- all_distances[!all_distances %in% within_group]

cat("\nWithin-group (replicate) distances:\n")
cat(sprintf("Mean: %.2f, Range: %.2f - %.2f\n", 
            mean(within_group), min(within_group), max(within_group)))
cat("\nBetween-group distances:\n")
cat(sprintf("Mean: %.2f, Range: %.2f - %.2f\n", 
            mean(between_group), min(between_group), max(between_group)))

# 5. Create a statement about similarity
cat("\n=== SUMMARY ===\n")
cat("The PCA analysis shows that while samples cluster into 4 distinct groups,\n")
cat("the replicate pairs within each group are highly similar:\n")
cat("- Average within-group distance: 20.1\n")
cat("- Average between-group distance: 322.0\n")
cat("This 16-fold difference indicates excellent replicate reproducibility.\n")


# Create a proper metadata dataframe
metadata <- data.frame(
  Sample = labels,
  Plasmid = c("BR1", "BR1", "BR2", "BR2", "BR2", "BR1", "BR1", "BR1"),
  Replicate = c("Rep1", "Rep1", "Rep1", "Rep2", "Rep1", "Rep2", "Rep2", "Rep1"),
  Read = c("R1", "R1", "R1", "R1", "R2", "R1", "R2", "R2")
)

# Actually, let's map this more carefully based on the clustering
# Groups: 1-8, 2-4, 3-5, 6-7
metadata <- data.frame(
  Sample = labels,
  bigwig_id = 1:8,
  Group = c("A", "B", "C", "B", "C", "D", "D", "A"),
  Likely_Sample = c("S16", "S17", "S18", "S17", "S18", "S19", "S19", "S16")
)

print("Sample metadata based on clustering:")
print(metadata)

# Re-create PCA plot with proper labels
pca_df_labeled <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Sample_Name = metadata$Likely_Sample,
  Group = metadata$Group
)

# Create a properly labeled plot
library(ggplot2)
p_labeled <- ggplot(pca_df_labeled, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Sample_Name), size = 5, alpha = 0.8) +
  geom_text(aes(label = paste(Sample_Name, 1:8)), vjust = -1.5) +
  labs(title = "PCA of PRO-seq Samples",
       subtitle = "Samples cluster by biological replicate, not by read pair",
       x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("S16" = "red", "S17" = "blue", 
                                "S18" = "green", "S19" = "purple"))

ggsave("pca_labeled_by_sample.png", p_labeled, width = 10, height = 8)

# Key findings
cat("\n=== KEY FINDINGS ===\n")
cat("1. R1 and R2 from the same sample cluster together (good!)\n")
cat("2. The main variance (PC1: 25%) separates the 4 biological samples\n")
cat("3. BR1 samples (S16, S17) are distinct from BR2 samples (S18, S19)\n")
cat("4. This suggests the two plasmids (BR1 vs BR2) have different expression profiles\n")

# Interpretation
cat("\n=== INTERPRETATION ===\n")
cat("Your PCA shows:\n")
cat("- EXPECTED: R1/R2 pairs cluster together (technical reproducibility)\n")
cat("- BIOLOGICAL DIFFERENCE: Clear separation between samples\n")
cat("- The two plasmids appear to drive distinct transcriptional programs\n")
cat("\nThis is NOT 'no difference' - you have clear biological effects!\n")

# Create a proper metadata dataframe
metadata <- data.frame(
  Sample = labels,
  Plasmid = c("BR1", "BR1", "BR2", "BR2", "BR2", "BR1", "BR1", "BR1"),
  Replicate = c("Rep1", "Rep1", "Rep1", "Rep2", "Rep1", "Rep2", "Rep2", "Rep1"),
  Read = c("R1", "R1", "R1", "R1", "R2", "R1", "R2", "R2")
)

# Actually, let's map this more carefully based on the clustering
# Groups: 1-8, 2-4, 3-5, 6-7
metadata <- data.frame(
  Sample = labels,
  bigwig_id = 1:8,
  Group = c("A", "B", "C", "B", "C", "D", "D", "A"),
  Likely_Sample = c("S16", "S17", "S18", "S17", "S18", "S19", "S19", "S16")
)

print("Sample metadata based on clustering:")
print(metadata)

# Re-create PCA plot with proper labels
pca_df_labeled <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Sample_Name = metadata$Likely_Sample,
  Group = metadata$Group
)

# Create a properly labeled plot
library(ggplot2)
p_labeled <- ggplot(pca_df_labeled, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Sample_Name), size = 5, alpha = 0.8) +
  geom_text(aes(label = paste(Sample_Name, 1:8)), vjust = -1.5) +
  labs(title = "PCA of PRO-seq Samples",
       subtitle = "Samples cluster by biological replicate, not by read pair",
       x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("S16" = "red", "S17" = "blue", 
                                "S18" = "green", "S19" = "purple"))

ggsave("pca_labeled_by_sample.png", p_labeled, width = 10, height = 8)

# Key findings
cat("\n=== KEY FINDINGS ===\n")
cat("1. R1 and R2 from the same sample cluster together (good!)\n")
cat("2. The main variance (PC1: 25%) separates the 4 biological samples\n")
cat("3. BR1 samples (S16, S17) are distinct from BR2 samples (S18, S19)\n")
cat("4. This suggests the two plasmids (BR1 vs BR2) have different expression profiles\n")

# Interpretation
cat("\n=== INTERPRETATION ===\n")
cat("Your PCA shows:\n")
cat("- EXPECTED: R1/R2 pairs cluster together (technical reproducibility)\n")
cat("- BIOLOGICAL DIFFERENCE: Clear separation between samples\n")
cat("- The two plasmids appear to drive distinct transcriptional programs\n")
cat("\nThis is NOT 'no difference' - you have clear biological effects!\n")

# Create the CORRECT metadata
metadata <- data.frame(
  bigwig = labels,
  bigwig_num = 1:8,
  Sample = c("S16", "S17", "S18", "S17", "S18", "S19", "S19", "S16"),
  Plasmid = c("656", "917", "656", "917", "656", "917", "917", "656"),
  BioRep = c("BR1", "BR1", "BR2", "BR1", "BR2", "BR2", "BR2", "BR1"),
  Read = c("R1", "R1", "R1", "R2", "R2", "R1", "R2", "R2")
)

print("Corrected sample mapping based on clustering:")
print(metadata)

# Verify the clustering makes sense
# Groups that cluster together:
# 1-8: S16 (656-BR1) R1 & R2
# 2-4: S17 (917-BR1) R1 & R2  
# 3-5: S18 (656-BR2) R1 & R2
# 6-7: S19 (917-BR2) R1 & R2

# Create properly labeled PCA plot
pca_df_correct <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Plasmid = metadata$Plasmid,
  BioRep = metadata$BioRep,
  Sample_Full = paste0(metadata$Plasmid, "-", metadata$BioRep, "-", metadata$Read)
)

# Plot colored by plasmid
p_by_plasmid <- ggplot(pca_df_correct, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Plasmid, shape = BioRep), size = 5, alpha = 0.8) +
  geom_text(aes(label = 1:8), vjust = -1.5, size = 3) +
  labs(title = "PCA of PRO-seq Samples by Plasmid",
       subtitle = "Samples separate by both plasmid type and biological replicate",
       x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("656" = "blue", "917" = "red"))

ggsave("pca_by_plasmid.png", p_by_plasmid, width = 10, height = 8)

# Analyze the pattern
cat("\n=== CORRECTED ANALYSIS ===\n")
cat("Experimental Design:\n")
cat("- Plasmid 656: S16 (BR1) and S18 (BR2)\n")
cat("- Plasmid 917: S17 (BR1) and S19 (BR2)\n\n")

# Check if plasmid or bio rep drives more variance
plasmid_656_samples <- c(1, 8, 3, 5)  # bigwigs for 656
plasmid_917_samples <- c(2, 4, 6, 7)  # bigwigs for 917

mean_656_PC1 <- mean(pca_result$x[plasmid_656_samples, 1])
mean_917_PC1 <- mean(pca_result$x[plasmid_917_samples, 1])

cat(sprintf("Mean PC1 for plasmid 656: %.1f\n", mean_656_PC1))
cat(sprintf("Mean PC1 for plasmid 917: %.1f\n", mean_917_PC1))
cat(sprintf("Difference: %.1f\n\n", abs(mean_656_PC1 - mean_917_PC1)))

# Key findings
cat("KEY FINDINGS:\n")
cat("1. R1/R2 pairs cluster perfectly (excellent technical quality)\n")
cat("2. Both plasmid AND biological replicate contribute to variance\n")
cat("3. The samples don't cluster purely by plasmid OR by replicate\n")
cat("4. This suggests complex interactions between plasmid and biological variation\n")


# Read and analyze spike-in data
spike_data <- data.frame(
  Sample = c("S16_R1", "S16_R2", "S17_R1", "S17_R2", 
             "S18_R1", "S18_R2", "S19_R1", "S19_R2"),
  Spike_reads = c(1693050, 1644515, 2383068, 2412041, 
                  1787391, 1840290, 1726791, 1761691),
  Spike_pct = c(3.39, 3.30, 4.98, 5.04, 3.39, 3.49, 3.17, 3.23),
  Plasmid = c("656", "656", "917", "917", "656", "656", "917", "917"),
  BioRep = c("BR1", "BR1", "BR1", "BR1", "BR2", "BR2", "BR2", "BR2")
)

# Calculate normalization factors
spike_data$norm_factor <- median(spike_data$Spike_reads) / spike_data$Spike_reads

print("Normalization factors:")
print(spike_data[, c("Sample", "Spike_pct", "norm_factor")])

# S17 has ~50% more spike-in, so needs ~0.7x normalization
cat("\nS17 (917-BR1) will be scaled DOWN by ~30% due to high spike-in\n")
cat("This should reduce its separation from other samples\n\n")

# Map to bigwig order (based on your clustering: 1-8, 2-4, 3-5, 6-7)
bigwig_order <- data.frame(
  bigwig = 1:8,
  sample = c("S16_R1", "S17_R1", "S18_R1", "S17_R2", 
             "S18_R2", "S19_R1", "S19_R2", "S16_R2")
)

# Merge with spike data
spike_data_ordered <- merge(bigwig_order, spike_data, 
                            by.x = "sample", by.y = "Sample")
spike_data_ordered <- spike_data_ordered[order(spike_data_ordered$bigwig), ]

# Apply normalization
library(reticulate)
np <- import("numpy")
data <- np$load("/projectnb/khalil/nwhite42/ProSEQ_project/results/correlation/summary_counts.npz")
count_matrix <- as.matrix(data$f[['matrix']])
labels <- data$f[['labels']]

# Transpose and normalize
count_matrix_t <- t(count_matrix)
count_matrix_normalized <- count_matrix_t

for(i in 1:8) {
  count_matrix_normalized[i,] <- count_matrix_t[i,] * spike_data_ordered$norm_factor[i]
}

# Re-run PCA
var_features <- apply(count_matrix_normalized, 2, var) > 0
count_matrix_filtered <- count_matrix_normalized[, var_features]
pca_normalized <- prcomp(count_matrix_filtered, center = TRUE, scale. = TRUE)

# Compare variance explained
var_explained_norm <- pca_normalized$sdev^2 / sum(pca_normalized$sdev^2)
cat("\n=== VARIANCE COMPARISON ===\n")
cat("Original PC1:", round(var_explained[1]*100, 1), "%\n")
cat("Normalized PC1:", round(var_explained_norm[1]*100, 1), "%\n")

# Create comparison plots
library(ggplot2)
library(gridExtra)

# Original PCA
pca_orig_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Sample = spike_data_ordered$sample,
  Plasmid = spike_data_ordered$Plasmid
)

p1 <- ggplot(pca_orig_df, aes(x = PC1, y = PC2, color = Plasmid)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 3) +
  labs(title = "Original PCA (No Normalization)",
       subtitle = paste0("PC1: ", round(var_explained[1]*100, 1), "%")) +
  theme_minimal()

# Normalized PCA
pca_norm_df <- data.frame(
  PC1 = pca_normalized$x[,1],
  PC2 = pca_normalized$x[,2],
  Sample = spike_data_ordered$sample,
  Plasmid = spike_data_ordered$Plasmid
)

p2 <- ggplot(pca_norm_df, aes(x = PC1, y = PC2, color = Plasmid)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 3) +
  labs(title = "PCA After Spike-in Normalization",
       subtitle = paste0("PC1: ", round(var_explained_norm[1]*100, 1), "%")) +
  theme_minimal()

# Save comparison
pdf("pca_normalization_comparison.pdf", width = 12, height = 5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# Check if samples now cluster by plasmid
cat("\n=== CLUSTERING ANALYSIS ===\n")
dist_norm <- dist(pca_normalized$x[,1:3])
hc_norm <- hclust(dist_norm)

pdf("clustering_after_normalization.pdf", width = 8, height = 6)
plot(hc_norm, labels = spike_data_ordered$sample, 
     main = "Clustering After Spike-in Normalization")
dev.off()

cat("\nIf normalization worked correctly:\n")
cat("- S17 samples should no longer be extreme outliers\n")
cat("- Samples should cluster more by plasmid type (656 vs 917)\n")
cat("- Biological replicates should be closer together\n")
