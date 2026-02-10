# 03_normalization.R
# Normalization and transformation of metabolomics data
# This script applies normalization methods to correct for technical variation

# =============================================================================
# SETUP
# =============================================================================

# Load required libraries
library(notame)
library(tidyverse)
library(Biobase)

# Source configuration
source("config.R")

message("=== Step 3: Normalization & Transformation ===\n")

# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading filtered MetaboSet object...")

# Load the filtered data from step 2
load(file.path(OBJECTS_PATH, "metaboset_filtered.RData"))

message("✓ Loaded MetaboSet with ", nrow(metaboset_filtered), " features and ", 
        ncol(metaboset_filtered), " samples\n")

# Keep a copy of pre-normalized data for comparison
metaboset_prenorm <- metaboset_filtered

# =============================================================================
# MISSING VALUE IMPUTATION
# =============================================================================

message("Handling missing values...")

# Count missing values
n_missing <- sum(is.na(exprs(metaboset_filtered)))
pct_missing <- n_missing / (nrow(metaboset_filtered) * ncol(metaboset_filtered)) * 100

message("  Missing values: ", n_missing, " (", round(pct_missing, 2), "%)")

if (pct_missing > 0) {
  if (IMPUTATION_METHOD != "none") {
    message("  Performing ", IMPUTATION_METHOD, " imputation...")
    
    # Impute missing values
    # Note: notame provides impute_rf() for random forest imputation
    # Here we use a simple approach - replace NA with half minimum
    # Adjust based on your specific needs
    
    if (IMPUTATION_METHOD == "knn") {
      # Use k-nearest neighbors imputation
      if (requireNamespace("impute", quietly = TRUE)) {
        expr_matrix <- exprs(metaboset_filtered)
        expr_imputed <- impute::impute.knn(expr_matrix, k = K_NEIGHBORS)$data
        exprs(metaboset_filtered) <- expr_imputed
        message("✓ KNN imputation complete")
      } else {
        warning("Package 'impute' not available. Skipping imputation.")
      }
    } else if (IMPUTATION_METHOD == "half_min") {
      # Replace NA with half of minimum detected value per feature
      expr_matrix <- exprs(metaboset_filtered)
      for (i in 1:nrow(expr_matrix)) {
        min_val <- min(expr_matrix[i, ], na.rm = TRUE)
        expr_matrix[i, is.na(expr_matrix[i, ])] <- min_val / 2
      }
      exprs(metaboset_filtered) <- expr_matrix
      message("✓ Half-minimum imputation complete")
    }
    
    # Verify no missing values remain
    n_missing_after <- sum(is.na(exprs(metaboset_filtered)))
    if (n_missing_after > 0) {
      warning("  Still ", n_missing_after, " missing values after imputation")
    }
  } else {
    message("  Skipping imputation (method = 'none')")
  }
}

# =============================================================================
# NORMALIZATION
# =============================================================================

message("\nApplying normalization...")
message("  Method: ", NORMALIZATION_METHOD)

if (NORMALIZATION_METHOD != "none") {
  
  if (NORMALIZATION_METHOD == "pqn") {
    # Probabilistic Quotient Normalization
    # This is recommended for metabolomics data
    
    # Calculate reference spectrum (median across samples)
    expr_matrix <- exprs(metaboset_filtered)
    reference <- apply(expr_matrix, 1, median, na.rm = TRUE)
    
    # Calculate quotients
    quotients <- expr_matrix / reference
    
    # Calculate normalization factors (median of quotients per sample)
    norm_factors <- apply(quotients, 2, median, na.rm = TRUE)
    
    # Normalize
    expr_normalized <- sweep(expr_matrix, 2, norm_factors, "/")
    exprs(metaboset_filtered) <- expr_normalized
    
    message("✓ PQN normalization complete")
    
  } else if (NORMALIZATION_METHOD == "median") {
    # Median normalization
    expr_matrix <- exprs(metaboset_filtered)
    sample_medians <- apply(expr_matrix, 2, median, na.rm = TRUE)
    global_median <- median(sample_medians, na.rm = TRUE)
    norm_factors <- sample_medians / global_median
    expr_normalized <- sweep(expr_matrix, 2, norm_factors, "/")
    exprs(metaboset_filtered) <- expr_normalized
    
    message("✓ Median normalization complete")
    
  } else if (NORMALIZATION_METHOD == "sum") {
    # Total sum normalization
    expr_matrix <- exprs(metaboset_filtered)
    sample_sums <- colSums(expr_matrix, na.rm = TRUE)
    global_mean <- mean(sample_sums, na.rm = TRUE)
    norm_factors <- sample_sums / global_mean
    expr_normalized <- sweep(expr_matrix, 2, norm_factors, "/")
    exprs(metaboset_filtered) <- expr_normalized
    
    message("✓ Sum normalization complete")
    
  } else if (NORMALIZATION_METHOD == "quantile") {
    # Quantile normalization
    if (requireNamespace("preprocessCore", quietly = TRUE)) {
      expr_matrix <- exprs(metaboset_filtered)
      expr_normalized <- preprocessCore::normalize.quantiles(expr_matrix)
      dimnames(expr_normalized) <- dimnames(expr_matrix)
      exprs(metaboset_filtered) <- expr_normalized
      message("✓ Quantile normalization complete")
    } else {
      warning("Package 'preprocessCore' not available. Skipping normalization.")
    }
    
  } else {
    warning("Unknown normalization method: ", NORMALIZATION_METHOD)
  }
  
} else {
  message("  Skipping normalization (method = 'none')")
}

# =============================================================================
# TRANSFORMATION
# =============================================================================

message("\nApplying transformation...")
message("  Method: ", TRANSFORMATION)

if (TRANSFORMATION != "none") {
  
  if (TRANSFORMATION == "log") {
    # Log transformation
    # Add small constant to avoid log(0)
    expr_matrix <- exprs(metaboset_filtered)
    min_nonzero <- min(expr_matrix[expr_matrix > 0], na.rm = TRUE)
    const <- min_nonzero / 10
    
    expr_transformed <- log2(expr_matrix + const)
    exprs(metaboset_filtered) <- expr_transformed
    
    message("✓ Log2 transformation complete")
    
  } else if (TRANSFORMATION == "sqrt") {
    # Square root transformation
    expr_matrix <- exprs(metaboset_filtered)
    expr_transformed <- sqrt(expr_matrix)
    exprs(metaboset_filtered) <- expr_transformed
    
    message("✓ Square root transformation complete")
    
  } else {
    warning("Unknown transformation method: ", TRANSFORMATION)
  }
  
} else {
  message("  Skipping transformation (method = 'none')")
}

# =============================================================================
# VISUALIZATION
# =============================================================================

message("\nGenerating normalization QC plots...")

# Create plots directory
norm_plots_dir <- file.path(PLOTS_PATH, "normalization")
dir.create(norm_plots_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Boxplot comparison: before vs after normalization
pdf(file.path(norm_plots_dir, "boxplot_comparison.pdf"), 
    width = PLOT_WIDTH * 1.5, height = PLOT_HEIGHT)
par(mfrow = c(1, 2))

# Before normalization
boxplot(log2(exprs(metaboset_prenorm) + 1), 
        main = "Before Normalization",
        xlab = "Samples", 
        ylab = "log2(Intensity)",
        las = 2, 
        cex.axis = 0.5,
        outline = FALSE)

# After normalization
boxplot(exprs(metaboset_filtered), 
        main = "After Normalization",
        xlab = "Samples", 
        ylab = if (TRANSFORMATION == "log") "log2(Intensity)" else "Intensity",
        las = 2, 
        cex.axis = 0.5,
        outline = FALSE)

dev.off()

# 2. Density plots
pdf(file.path(norm_plots_dir, "density_comparison.pdf"), 
    width = PLOT_WIDTH * 1.5, height = PLOT_HEIGHT)
par(mfrow = c(1, 2))

# Before normalization
plot(density(log2(exprs(metaboset_prenorm)[, 1] + 1), na.rm = TRUE),
     main = "Before Normalization",
     xlab = "log2(Intensity)",
     ylab = "Density",
     col = 1, lwd = 1)
for (i in 2:min(20, ncol(metaboset_prenorm))) {
  lines(density(log2(exprs(metaboset_prenorm)[, i] + 1), na.rm = TRUE), 
        col = i, lwd = 1)
}

# After normalization
plot(density(exprs(metaboset_filtered)[, 1], na.rm = TRUE),
     main = "After Normalization",
     xlab = if (TRANSFORMATION == "log") "log2(Intensity)" else "Intensity",
     ylab = "Density",
     col = 1, lwd = 1)
for (i in 2:min(20, ncol(metaboset_filtered))) {
  lines(density(exprs(metaboset_filtered)[, i], na.rm = TRUE), 
        col = i, lwd = 1)
}

dev.off()

# 3. Total intensity per sample
pdf(file.path(norm_plots_dir, "total_intensity.pdf"), 
    width = PLOT_WIDTH * 1.5, height = PLOT_HEIGHT)
par(mfrow = c(1, 2))

# Before
barplot(colSums(exprs(metaboset_prenorm), na.rm = TRUE),
        main = "Before Normalization",
        xlab = "Samples",
        ylab = "Total Intensity",
        las = 2,
        cex.names = 0.5)

# After
barplot(colSums(exprs(metaboset_filtered), na.rm = TRUE),
        main = "After Normalization",
        xlab = "Samples",
        ylab = if (TRANSFORMATION == "log") "Total log2(Intensity)" else "Total Intensity",
        las = 2,
        cex.names = 0.5)

dev.off()

message("✓ Normalization plots saved to: ", norm_plots_dir)

# =============================================================================
# SAVE NORMALIZED DATA
# =============================================================================

message("\nSaving normalized MetaboSet...")

# Rename for clarity
metaboset_normalized <- metaboset_filtered

# Save normalized object
save_path <- file.path(OBJECTS_PATH, "metaboset_normalized.RData")
save(metaboset_normalized, file = save_path)

message("✓ Saved to: ", save_path)

# Export normalized data as CSV
csv_path <- file.path(TABLES_PATH, "normalized_data.csv")
expr_df <- as.data.frame(exprs(metaboset_normalized))
expr_df <- tibble::rownames_to_column(expr_df, "Feature_ID")
write_csv(expr_df, csv_path)

message("✓ Exported to CSV: ", csv_path)

# Save normalization summary
summary_file <- file.path(TABLES_PATH, "03_normalization_summary.txt")
sink(summary_file)
cat("=== Normalization Summary ===\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Normalization method:", NORMALIZATION_METHOD, "\n")
cat("Transformation method:", TRANSFORMATION, "\n")
cat("Imputation method:", IMPUTATION_METHOD, "\n\n")
cat("Data dimensions:\n")
cat("  Features:", nrow(metaboset_normalized), "\n")
cat("  Samples:", ncol(metaboset_normalized), "\n\n")
cat("Missing values:\n")
cat("  Before imputation:", n_missing, "(", round(pct_missing, 2), "%)\n")
cat("  After imputation:", sum(is.na(exprs(metaboset_normalized))), "\n")
sink()

message("✓ Summary saved to: ", summary_file)

message("\n=== Normalization Complete ===\n")
message("Next step: Run 04_visualization.R")
