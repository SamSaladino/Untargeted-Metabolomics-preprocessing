# 04_visualization.R
# Exploratory data analysis and visualization
# This script generates PCA plots and other exploratory visualizations

# =============================================================================
# SETUP
# =============================================================================

# Load required libraries
library(notame)
library(tidyverse)
library(Biobase)
library(ggplot2)
library(pcaMethods)
library(RColorBrewer)
library(pheatmap)

# Source configuration
source("config.R")

message("=== Step 4: Visualization & Exploratory Analysis ===\n")

# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading normalized MetaboSet object...")

# Load the normalized data from step 3
load(file.path(OBJECTS_PATH, "metaboset_normalized.RData"))

message("✓ Loaded MetaboSet with ", nrow(metaboset_normalized), " features and ", 
        ncol(metaboset_normalized), " samples\n")

# Create visualization directory
vis_plots_dir <- file.path(PLOTS_PATH, "visualization")
dir.create(vis_plots_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PRINCIPAL COMPONENT ANALYSIS (PCA)
# =============================================================================

message("Performing PCA...")

# Prepare data for PCA (samples as rows, features as columns)
pca_data <- t(exprs(metaboset_normalized))

# Remove features with zero variance
feature_var <- apply(pca_data, 2, var, na.rm = TRUE)
pca_data <- pca_data[, feature_var > 0 & !is.na(feature_var)]

message("  Using ", ncol(pca_data), " features for PCA")

# Perform PCA
pca_result <- pcaMethods::pca(pca_data, 
                               nPcs = min(N_PCA_COMPONENTS, nrow(pca_data) - 1),
                               method = "svd",
                               center = TRUE,
                               scale = "uv")

# Extract scores and loadings
pca_scores <- as.data.frame(scores(pca_result))
pca_loadings <- as.data.frame(loadings(pca_result))

# Get variance explained
variance_explained <- (pca_result@R2 * 100)

message("✓ PCA complete")
message("  Variance explained by PC1: ", round(variance_explained[1], 2), "%")
message("  Variance explained by PC2: ", round(variance_explained[2], 2), "%")

# Add sample metadata to PCA scores
pca_scores <- cbind(pca_scores, pData(metaboset_normalized))

# =============================================================================
# PCA PLOTS
# =============================================================================

message("\nGenerating PCA plots...")

# Get color palette
n_groups <- length(unique(pca_scores[[GROUP_COLUMN]]))
colors <- brewer.pal(min(9, max(3, n_groups)), COLOR_PALETTE)

# 1. PC1 vs PC2 colored by group
p1 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = .data[[GROUP_COLUMN]])) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "PCA: PC1 vs PC2",
    x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
    color = GROUP_COLUMN
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  scale_color_manual(values = colors)

ggsave(file.path(vis_plots_dir, "pca_pc1_pc2.pdf"), 
       p1, width = PLOT_WIDTH, height = PLOT_HEIGHT)
ggsave(file.path(vis_plots_dir, "pca_pc1_pc2.png"), 
       p1, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI)

# 2. PC2 vs PC3 (if available)
if (ncol(pca_scores) >= 5) {  # PC1, PC2, PC3 + metadata columns
  p2 <- ggplot(pca_scores, aes(x = PC2, y = PC3, color = .data[[GROUP_COLUMN]])) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = "PCA: PC2 vs PC3",
      x = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
      y = paste0("PC3 (", round(variance_explained[3], 2), "%)"),
      color = GROUP_COLUMN
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    ) +
    scale_color_manual(values = colors)
  
  ggsave(file.path(vis_plots_dir, "pca_pc2_pc3.pdf"), 
         p2, width = PLOT_WIDTH, height = PLOT_HEIGHT)
}

# 3. Scree plot (variance explained)
scree_data <- data.frame(
  PC = paste0("PC", 1:length(variance_explained)),
  Variance = variance_explained
)
scree_data$PC <- factor(scree_data$PC, levels = scree_data$PC)

p3 <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(group = 1), color = "red", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "PCA Scree Plot",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(vis_plots_dir, "pca_scree_plot.pdf"), 
       p3, width = PLOT_WIDTH, height = PLOT_HEIGHT)

# 4. PCA biplot (if batch information is available)
if (BATCH_COLUMN %in% colnames(pca_scores)) {
  p4 <- ggplot(pca_scores, aes(x = PC1, y = PC2, 
                                color = .data[[GROUP_COLUMN]], 
                                shape = .data[[BATCH_COLUMN]])) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = "PCA: Colored by Group, Shaped by Batch",
      x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
      y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
      color = GROUP_COLUMN,
      shape = BATCH_COLUMN
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    ) +
    scale_color_manual(values = colors)
  
  ggsave(file.path(vis_plots_dir, "pca_batch_effect.pdf"), 
         p4, width = PLOT_WIDTH + 2, height = PLOT_HEIGHT)
}

message("✓ PCA plots saved")

# =============================================================================
# HEATMAP
# =============================================================================

message("\nGenerating heatmap...")

# Select top variable features for heatmap (e.g., top 50)
feature_vars <- apply(exprs(metaboset_normalized), 1, var, na.rm = TRUE)
top_features <- names(sort(feature_vars, decreasing = TRUE)[1:min(50, nrow(metaboset_normalized))])

# Prepare heatmap data
heatmap_data <- exprs(metaboset_normalized)[top_features, ]

# Scale features (z-score)
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# Prepare annotation
annotation_col <- data.frame(
  Group = pData(metaboset_normalized)[[GROUP_COLUMN]],
  row.names = colnames(heatmap_data)
)

# Create heatmap
pdf(file.path(vis_plots_dir, "heatmap_top_features.pdf"), 
    width = PLOT_WIDTH, height = PLOT_HEIGHT)

pheatmap(
  heatmap_data_scaled,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap of Top 50 Variable Features",
  fontsize = 8,
  fontsize_col = 6
)

dev.off()

message("✓ Heatmap saved")

# =============================================================================
# SAMPLE CORRELATION MATRIX
# =============================================================================

message("\nGenerating sample correlation matrix...")

# Calculate correlation between samples
sample_cor <- cor(exprs(metaboset_normalized), use = "pairwise.complete.obs")

# Plot correlation heatmap
pdf(file.path(vis_plots_dir, "sample_correlation.pdf"), 
    width = PLOT_WIDTH, height = PLOT_HEIGHT)

pheatmap(
  sample_cor,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Sample-to-Sample Correlation",
  fontsize = 6
)

dev.off()

message("✓ Correlation matrix saved")

# =============================================================================
# COEFFICIENT OF VARIATION (CV) PLOT
# =============================================================================

message("\nGenerating CV plot...")

# Calculate CV for each feature
feature_means <- rowMeans(exprs(metaboset_normalized), na.rm = TRUE)
feature_sds <- apply(exprs(metaboset_normalized), 1, sd, na.rm = TRUE)
feature_cvs <- (feature_sds / feature_means) * 100

# Create CV distribution plot
cv_data <- data.frame(CV = feature_cvs)

p5 <- ggplot(cv_data, aes(x = CV)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  labs(
    title = "Coefficient of Variation Distribution",
    x = "CV (%)",
    y = "Number of Features"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(vis_plots_dir, "cv_distribution.pdf"), 
       p5, width = PLOT_WIDTH, height = PLOT_HEIGHT)

message("✓ CV plot saved")

# =============================================================================
# SAVE RESULTS
# =============================================================================

message("\nSaving PCA results...")

# Save PCA scores
write_csv(pca_scores, file.path(TABLES_PATH, "pca_scores.csv"))

# Save PCA loadings
loadings_df <- tibble::rownames_to_column(pca_loadings, "Feature_ID")
write_csv(loadings_df, file.path(TABLES_PATH, "pca_loadings.csv"))

# Save variance explained
variance_df <- data.frame(
  PC = paste0("PC", 1:length(variance_explained)),
  Variance_Explained = variance_explained,
  Cumulative_Variance = cumsum(variance_explained)
)
write_csv(variance_df, file.path(TABLES_PATH, "pca_variance.csv"))

message("✓ PCA results saved")

# Save visualization summary
summary_file <- file.path(TABLES_PATH, "04_visualization_summary.txt")
sink(summary_file)
cat("=== Visualization Summary ===\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("PCA Analysis:\n")
cat("  Number of components:", length(variance_explained), "\n")
cat("  Variance explained by PC1:", round(variance_explained[1], 2), "%\n")
cat("  Variance explained by PC2:", round(variance_explained[2], 2), "%\n")
cat("  Cumulative variance (PC1+PC2):", round(sum(variance_explained[1:2]), 2), "%\n\n")
cat("Feature statistics:\n")
cat("  Mean CV:", round(mean(feature_cvs, na.rm = TRUE), 2), "%\n")
cat("  Median CV:", round(median(feature_cvs, na.rm = TRUE), 2), "%\n\n")
cat("Plots generated:\n")
cat("  - PCA score plots\n")
cat("  - Scree plot\n")
cat("  - Heatmap of top variable features\n")
cat("  - Sample correlation matrix\n")
cat("  - CV distribution\n")
sink()

message("✓ Summary saved to: ", summary_file)

message("\n=== Visualization Complete ===\n")
message("All plots saved to: ", vis_plots_dir)
message("Next step: Run 05_statistical_analysis.R")
