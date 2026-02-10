# 02_qc_filtering.R
# Quality control and filtering of metabolomics data
# This script filters low-quality features and samples

# =============================================================================
# SETUP
# =============================================================================

# Load required libraries
library(notame)
library(tidyverse)
library(Biobase)

# Source configuration
source("config.R")

message("=== Step 2: Quality Control & Filtering ===\n")

# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading MetaboSet object...")

# Load the raw data from step 1
load(file.path(OBJECTS_PATH, "metaboset_raw.RData"))

message("✓ Loaded MetaboSet with ", nrow(metaboset), " features and ", 
        ncol(metaboset), " samples\n")

# Store original dimensions for reporting
n_features_original <- nrow(metaboset)
n_samples_original <- ncol(metaboset)

# =============================================================================
# IDENTIFY QC SAMPLES
# =============================================================================

message("Identifying QC samples...")

# Check if QC column exists
if (QC_COLUMN %in% colnames(pData(metaboset))) {
  qc_samples <- pData(metaboset)[[QC_COLUMN]] == QC_LABEL
  n_qc <- sum(qc_samples, na.rm = TRUE)
  message("✓ Found ", n_qc, " QC samples")
} else {
  warning("QC column '", QC_COLUMN, "' not found in metadata. Skipping QC-based filtering.")
  qc_samples <- rep(FALSE, ncol(metaboset))
  n_qc <- 0
}

# =============================================================================
# FEATURE FILTERING
# =============================================================================

message("\nFiltering features...")

# Calculate detection rate (fraction of non-missing values)
detection_matrix <- !is.na(exprs(metaboset))

# Detection rate in QC samples
if (n_qc > 0) {
  detection_rate_qc <- rowMeans(detection_matrix[, qc_samples, drop = FALSE])
  qc_filter <- detection_rate_qc >= MIN_DETECTION_RATE
  message("  - Features with ≥", MIN_DETECTION_RATE * 100, 
          "% detection in QC: ", sum(qc_filter))
} else {
  qc_filter <- rep(TRUE, nrow(metaboset))
}

# Detection rate in all samples
detection_rate_all <- rowMeans(detection_matrix)
fraction_filter <- detection_rate_all >= MIN_FRACTION_SAMPLES
message("  - Features with ≥", MIN_FRACTION_SAMPLES * 100, 
        "% detection in samples: ", sum(fraction_filter))

# RSD filtering in QC samples
if (n_qc > 0) {
  # Calculate RSD (Relative Standard Deviation) in QC samples
  qc_data <- exprs(metaboset)[, qc_samples, drop = FALSE]
  qc_mean <- rowMeans(qc_data, na.rm = TRUE)
  qc_sd <- apply(qc_data, 1, sd, na.rm = TRUE)
  qc_rsd <- (qc_sd / qc_mean) * 100
  
  rsd_filter <- qc_rsd <= RSD_THRESHOLD | is.na(qc_rsd)
  message("  - Features with RSD ≤", RSD_THRESHOLD, "% in QC: ", sum(rsd_filter))
} else {
  rsd_filter <- rep(TRUE, nrow(metaboset))
}

# Combine filters
feature_filter <- qc_filter & fraction_filter & rsd_filter
n_features_removed <- sum(!feature_filter)

message("\n  Removing ", n_features_removed, " features (", 
        round(n_features_removed / n_features_original * 100, 1), "%)")

# Apply feature filter
metaboset_filtered <- metaboset[feature_filter, ]

message("✓ Retained ", nrow(metaboset_filtered), " features")

# =============================================================================
# SAMPLE FILTERING
# =============================================================================

message("\nFiltering samples...")

# Calculate feature detection rate per sample
sample_detection <- colMeans(!is.na(exprs(metaboset_filtered)))
sample_filter <- sample_detection >= MIN_FEATURES_PER_SAMPLE

n_samples_removed <- sum(!sample_filter)
message("  Removing ", n_samples_removed, " samples with <", 
        MIN_FEATURES_PER_SAMPLE * 100, "% feature detection")

# Identify which samples are being removed
if (n_samples_removed > 0) {
  removed_samples <- colnames(metaboset_filtered)[!sample_filter]
  message("  Removed samples: ", paste(removed_samples, collapse = ", "))
}

# Apply sample filter
metaboset_filtered <- metaboset_filtered[, sample_filter]

message("✓ Retained ", ncol(metaboset_filtered), " samples")

# =============================================================================
# REMOVE BLANK SAMPLES
# =============================================================================

message("\nHandling blank samples...")

# Remove blank samples if they exist
if (GROUP_COLUMN %in% colnames(pData(metaboset_filtered))) {
  blank_samples <- pData(metaboset_filtered)[[GROUP_COLUMN]] == BLANK_LABEL
  n_blanks <- sum(blank_samples, na.rm = TRUE)
  
  if (n_blanks > 0) {
    message("  Removing ", n_blanks, " blank samples")
    metaboset_filtered <- metaboset_filtered[, !blank_samples]
  } else {
    message("  No blank samples found")
  }
}

# =============================================================================
# DRIFT CORRECTION (OPTIONAL)
# =============================================================================

# If you have injection order information, you can perform drift correction
# This is particularly useful for large batches
# Example:
# if ("Injection_Order" %in% colnames(pData(metaboset_filtered))) {
#   message("\nPerforming drift correction...")
#   metaboset_filtered <- correct_drift(metaboset_filtered)
# }

# =============================================================================
# VISUALIZATION
# =============================================================================

message("\nGenerating QC plots...")

# Create plots directory if needed
qc_plots_dir <- file.path(PLOTS_PATH, "qc")
dir.create(qc_plots_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Feature detection distribution
pdf(file.path(qc_plots_dir, "feature_detection_distribution.pdf"), 
    width = PLOT_WIDTH, height = PLOT_HEIGHT)
hist(detection_rate_all, 
     breaks = 50,
     main = "Feature Detection Rate Distribution",
     xlab = "Detection Rate (fraction of samples)",
     ylab = "Number of Features",
     col = "steelblue")
abline(v = MIN_FRACTION_SAMPLES, col = "red", lwd = 2, lty = 2)
legend("topright", 
       legend = paste("Threshold:", MIN_FRACTION_SAMPLES),
       col = "red", lty = 2, lwd = 2)
dev.off()

# 2. Sample detection distribution
pdf(file.path(qc_plots_dir, "sample_detection_distribution.pdf"), 
    width = PLOT_WIDTH, height = PLOT_HEIGHT)
hist(sample_detection, 
     breaks = 50,
     main = "Sample Detection Rate Distribution",
     xlab = "Detection Rate (fraction of features)",
     ylab = "Number of Samples",
     col = "coral")
abline(v = MIN_FEATURES_PER_SAMPLE, col = "red", lwd = 2, lty = 2)
legend("topleft", 
       legend = paste("Threshold:", MIN_FEATURES_PER_SAMPLE),
       col = "red", lty = 2, lwd = 2)
dev.off()

# 3. RSD distribution in QC samples (if applicable)
if (n_qc > 0) {
  pdf(file.path(qc_plots_dir, "qc_rsd_distribution.pdf"), 
      width = PLOT_WIDTH, height = PLOT_HEIGHT)
  hist(qc_rsd[!is.na(qc_rsd)], 
       breaks = 50,
       main = "RSD Distribution in QC Samples",
       xlab = "RSD (%)",
       ylab = "Number of Features",
       col = "lightgreen")
  abline(v = RSD_THRESHOLD, col = "red", lwd = 2, lty = 2)
  legend("topright", 
         legend = paste("Threshold:", RSD_THRESHOLD, "%"),
         col = "red", lty = 2, lwd = 2)
  dev.off()
}

message("✓ QC plots saved to: ", qc_plots_dir)

# =============================================================================
# SAVE FILTERED DATA
# =============================================================================

message("\nSaving filtered MetaboSet...")

# Save filtered object
save_path <- file.path(OBJECTS_PATH, "metaboset_filtered.RData")
save(metaboset_filtered, file = save_path)

message("✓ Saved to: ", save_path)

# Save filtering summary
summary_file <- file.path(TABLES_PATH, "02_filtering_summary.txt")
sink(summary_file)
cat("=== QC Filtering Summary ===\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Original data:\n")
cat("  Features:", n_features_original, "\n")
cat("  Samples:", n_samples_original, "\n\n")
cat("Filtering parameters:\n")
cat("  Min detection rate (QC):", MIN_DETECTION_RATE, "\n")
cat("  Min detection rate (all):", MIN_FRACTION_SAMPLES, "\n")
cat("  Max RSD (QC):", RSD_THRESHOLD, "%\n")
cat("  Min features per sample:", MIN_FEATURES_PER_SAMPLE, "\n\n")
cat("Filtered data:\n")
cat("  Features:", nrow(metaboset_filtered), "\n")
cat("  Samples:", ncol(metaboset_filtered), "\n\n")
cat("Removed:\n")
cat("  Features:", n_features_removed, "(", 
    round(n_features_removed / n_features_original * 100, 1), "%)\n")
cat("  Samples:", n_samples_removed, "(", 
    round(n_samples_removed / n_samples_original * 100, 1), "%)\n")
sink()

message("✓ Summary saved to: ", summary_file)

message("\n=== QC Filtering Complete ===\n")
message("Retained ", nrow(metaboset_filtered), " features and ", 
        ncol(metaboset_filtered), " samples")
message("Next step: Run 03_normalization.R")
