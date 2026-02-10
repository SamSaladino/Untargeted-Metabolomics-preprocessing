# 05_statistical_analysis.R
# Statistical analysis of metabolomics data
# This script performs differential analysis between groups

# =============================================================================
# SETUP
# =============================================================================

# Load required libraries
library(notame)
library(tidyverse)
library(Biobase)
library(ggplot2)

# Source configuration
source("config.R")

message("=== Step 5: Statistical Analysis ===\n")

# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading normalized MetaboSet object...")

# Load the normalized data from step 3
load(file.path(OBJECTS_PATH, "metaboset_normalized.RData"))

message("✓ Loaded MetaboSet with ", nrow(metaboset_normalized), " features and ", 
        ncol(metaboset_normalized), " samples\n")

# Create statistical analysis directory
stats_dir <- file.path(PLOTS_PATH, "statistics")
dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PREPARE DATA FOR ANALYSIS
# =============================================================================

message("Preparing data for statistical analysis...")

# Extract expression data and metadata
expr_data <- exprs(metaboset_normalized)
sample_data <- pData(metaboset_normalized)

# Get unique groups
groups <- unique(sample_data[[GROUP_COLUMN]])
message("  Groups found: ", paste(groups, collapse = ", "))

# Check if we have exactly 2 groups for pairwise comparison
if (length(groups) == 2) {
  message("  Performing pairwise comparison between ", groups[1], " and ", groups[2])
  comparison_mode <- "pairwise"
} else if (length(groups) > 2) {
  message("  Multiple groups detected. Will perform ANOVA.")
  comparison_mode <- "anova"
} else {
  stop("Need at least 2 groups for statistical analysis")
}

# =============================================================================
# PAIRWISE STATISTICAL TESTS (FOR 2 GROUPS)
# =============================================================================

if (comparison_mode == "pairwise") {
  
  message("\nPerforming pairwise statistical tests...")
  
  # Get samples for each group
  group1_samples <- sample_data[[SAMPLE_ID_COLUMN]][sample_data[[GROUP_COLUMN]] == groups[1]]
  group2_samples <- sample_data[[SAMPLE_ID_COLUMN]][sample_data[[GROUP_COLUMN]] == groups[2]]
  
  # Initialize results data frame
  results <- data.frame(
    Feature_ID = rownames(expr_data),
    stringsAsFactors = FALSE
  )
  
  # Calculate fold change and statistics for each feature
  message("  Calculating statistics for ", nrow(expr_data), " features...")
  
  for (i in 1:nrow(expr_data)) {
    
    # Get values for each group
    values_group1 <- expr_data[i, group1_samples]
    values_group2 <- expr_data[i, group2_samples]
    
    # Calculate means
    mean_group1 <- mean(values_group1, na.rm = TRUE)
    mean_group2 <- mean(values_group2, na.rm = TRUE)
    
    # Calculate fold change (group2 / group1)
    fold_change <- mean_group2 / mean_group1
    log2_fc <- log2(fold_change)
    
    # T-test
    if (sum(!is.na(values_group1)) >= 2 && sum(!is.na(values_group2)) >= 2) {
      t_test <- t.test(values_group2, values_group1)
      p_value <- t_test$p.value
    } else {
      p_value <- NA
    }
    
    # Store results
    results$Mean_Group1[i] <- mean_group1
    results$Mean_Group2[i] <- mean_group2
    results$Fold_Change[i] <- fold_change
    results$Log2_FC[i] <- log2_fc
    results$P_Value[i] <- p_value
  }
  
  # Adjust p-values for multiple testing
  results$P_Adjusted <- p.adjust(results$P_Value, method = ADJUST_METHOD)
  
  # Add significance column
  results$Significant <- results$P_Adjusted < ALPHA & abs(results$Log2_FC) > log2(FC_THRESHOLD)
  results$Significant[is.na(results$Significant)] <- FALSE
  
  # Add regulation direction
  results$Regulation <- ifelse(results$Significant,
                                ifelse(results$Log2_FC > 0, "Up", "Down"),
                                "Not Significant")
  
  message("✓ Statistical tests complete")
  message("  Significant features (p < ", ALPHA, ", FC > ", FC_THRESHOLD, "): ", 
          sum(results$Significant, na.rm = TRUE))
  message("    Up-regulated: ", sum(results$Regulation == "Up", na.rm = TRUE))
  message("    Down-regulated: ", sum(results$Regulation == "Down", na.rm = TRUE))
  
  # =============================================================================
  # VOLCANO PLOT
  # =============================================================================
  
  message("\nGenerating volcano plot...")
  
  # Prepare data for volcano plot
  volcano_data <- results %>%
    mutate(
      neg_log10_p = -log10(P_Adjusted),
      Color = case_when(
        Regulation == "Up" ~ "Up-regulated",
        Regulation == "Down" ~ "Down-regulated",
        TRUE ~ "Not Significant"
      )
    )
  
  # Create volcano plot
  p_volcano <- ggplot(volcano_data, aes(x = Log2_FC, y = neg_log10_p, color = Color)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-log2(FC_THRESHOLD), log2(FC_THRESHOLD)), 
               linetype = "dashed", color = "gray40") +
    scale_color_manual(
      values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "gray70"),
      name = "Regulation"
    ) +
    labs(
      title = paste0("Volcano Plot: ", groups[2], " vs ", groups[1]),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(file.path(stats_dir, "volcano_plot.pdf"), 
         p_volcano, width = PLOT_WIDTH, height = PLOT_HEIGHT)
  ggsave(file.path(stats_dir, "volcano_plot.png"), 
         p_volcano, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI)
  
  message("✓ Volcano plot saved")
  
  # =============================================================================
  # MA PLOT
  # =============================================================================
  
  message("Generating MA plot...")
  
  # Prepare data for MA plot
  ma_data <- results %>%
    mutate(
      A = (Mean_Group1 + Mean_Group2) / 2,  # Average expression
      M = Log2_FC,  # Log2 fold change
      Color = Color
    )
  
  # Create MA plot
  p_ma <- ggplot(ma_data, aes(x = A, y = M, color = Color)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_hline(yintercept = c(-log2(FC_THRESHOLD), log2(FC_THRESHOLD)), 
               linetype = "dashed", color = "gray40") +
    scale_color_manual(
      values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "gray70"),
      name = "Regulation"
    ) +
    labs(
      title = paste0("MA Plot: ", groups[2], " vs ", groups[1]),
      x = "Average Expression",
      y = "Log2 Fold Change"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(file.path(stats_dir, "ma_plot.pdf"), 
         p_ma, width = PLOT_WIDTH, height = PLOT_HEIGHT)
  
  message("✓ MA plot saved")
  
  # =============================================================================
  # SAVE RESULTS
  # =============================================================================
  
  message("\nSaving statistical results...")
  
  # Save full results
  results_sorted <- results %>%
    arrange(P_Adjusted, desc(abs(Log2_FC)))
  
  write_csv(results_sorted, file.path(TABLES_PATH, "statistical_results.csv"))
  
  # Save significant features only
  significant_features <- results_sorted %>%
    filter(Significant)
  
  write_csv(significant_features, file.path(TABLES_PATH, "significant_features.csv"))
  
  message("✓ Results saved")
  message("  All features: ", nrow(results_sorted))
  message("  Significant features: ", nrow(significant_features))
  
}

# =============================================================================
# ANOVA (FOR >2 GROUPS)
# =============================================================================

if (comparison_mode == "anova") {
  
  message("\nPerforming ANOVA...")
  
  # Initialize results
  results <- data.frame(
    Feature_ID = rownames(expr_data),
    stringsAsFactors = FALSE
  )
  
  # Perform ANOVA for each feature
  message("  Testing ", nrow(expr_data), " features...")
  
  for (i in 1:nrow(expr_data)) {
    
    # Prepare data for ANOVA
    feature_values <- expr_data[i, ]
    anova_data <- data.frame(
      Value = feature_values,
      Group = sample_data[[GROUP_COLUMN]]
    )
    
    # Perform ANOVA
    if (sum(!is.na(feature_values)) >= length(groups)) {
      anova_result <- aov(Value ~ Group, data = anova_data)
      anova_summary <- summary(anova_result)
      p_value <- anova_summary[[1]]$`Pr(>F)`[1]
    } else {
      p_value <- NA
    }
    
    # Calculate mean per group
    group_means <- tapply(feature_values, sample_data[[GROUP_COLUMN]], mean, na.rm = TRUE)
    
    # Store results
    results$P_Value[i] <- p_value
    for (grp in names(group_means)) {
      results[[paste0("Mean_", grp)]][i] <- group_means[grp]
    }
  }
  
  # Adjust p-values
  results$P_Adjusted <- p.adjust(results$P_Value, method = ADJUST_METHOD)
  results$Significant <- results$P_Adjusted < ALPHA
  results$Significant[is.na(results$Significant)] <- FALSE
  
  message("✓ ANOVA complete")
  message("  Significant features (p < ", ALPHA, "): ", sum(results$Significant, na.rm = TRUE))
  
  # Save results
  results_sorted <- results %>%
    arrange(P_Adjusted)
  
  write_csv(results_sorted, file.path(TABLES_PATH, "anova_results.csv"))
  
  significant_features <- results_sorted %>%
    filter(Significant)
  
  write_csv(significant_features, file.path(TABLES_PATH, "significant_features.csv"))
  
  message("✓ Results saved")
}

# =============================================================================
# SUMMARY REPORT
# =============================================================================

message("\nGenerating summary report...")

summary_file <- file.path(TABLES_PATH, "05_statistical_summary.txt")
sink(summary_file)
cat("=== Statistical Analysis Summary ===\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Analysis type:", comparison_mode, "\n")
cat("Groups:", paste(groups, collapse = ", "), "\n")
cat("Significance threshold (alpha):", ALPHA, "\n")
cat("P-value adjustment method:", ADJUST_METHOD, "\n")

if (comparison_mode == "pairwise") {
  cat("Fold change threshold:", FC_THRESHOLD, "\n\n")
  cat("Results:\n")
  cat("  Total features tested:", nrow(results), "\n")
  cat("  Significant features:", sum(results$Significant, na.rm = TRUE), "\n")
  cat("  Up-regulated:", sum(results$Regulation == "Up", na.rm = TRUE), "\n")
  cat("  Down-regulated:", sum(results$Regulation == "Down", na.rm = TRUE), "\n")
} else {
  cat("\nResults:\n")
  cat("  Total features tested:", nrow(results), "\n")
  cat("  Significant features:", sum(results$Significant, na.rm = TRUE), "\n")
}

sink()

message("✓ Summary saved to: ", summary_file)

message("\n=== Statistical Analysis Complete ===\n")
message("All results saved to: ", TABLES_PATH)
message("\nWorkflow complete! Review results in the '", RESULTS_PATH, "' directory.")
