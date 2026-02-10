# 01_data_import.R
# Import and prepare LC-MS metabolomics data for notame workflow
# This script loads raw data and creates a MetaboSet object

# =============================================================================
# SETUP
# =============================================================================

# Load required libraries
library(notame)
library(tidyverse)
library(Biobase)

# Source configuration
source("config.R")

message("=== Step 1: Data Import ===\n")

# =============================================================================
# LOAD SAMPLE METADATA
# =============================================================================

message("Loading sample metadata...")

# Read metadata file
if (!file.exists(METADATA_PATH)) {
  stop("Metadata file not found: ", METADATA_PATH, 
       "\nPlease create a metadata file with columns: ", 
       paste(c(SAMPLE_ID_COLUMN, GROUP_COLUMN, BATCH_COLUMN, QC_COLUMN), collapse = ", "))
}

metadata <- read_csv(METADATA_PATH, show_col_types = FALSE)

# Verify required columns exist
required_cols <- c(SAMPLE_ID_COLUMN, GROUP_COLUMN)
missing_cols <- required_cols[!required_cols %in% colnames(metadata)]
if (length(missing_cols) > 0) {
  stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

message("✓ Loaded metadata for ", nrow(metadata), " samples")
message("  - Groups: ", paste(unique(metadata[[GROUP_COLUMN]]), collapse = ", "))

# =============================================================================
# LOAD FEATURE DATA
# =============================================================================

message("\nLoading feature abundance data...")

# Option 1: Load from pre-processed feature table (CSV format)
if (file.exists(FEATURE_TABLE)) {
  message("Loading feature table from: ", FEATURE_TABLE)
  
  # Read feature table
  # Expected format: rows = features, columns = samples
  # First column should be feature ID, subsequent columns are sample abundances
  feature_data <- read_csv(FEATURE_TABLE, show_col_types = FALSE)
  
  # Extract feature IDs
  feature_ids <- feature_data[[1]]
  
  # Extract abundance matrix
  abundance_matrix <- as.matrix(feature_data[, -1])
  rownames(abundance_matrix) <- feature_ids
  
  # Create feature metadata (mz, RT, etc. if available)
  # Adjust based on your feature table structure
  feature_meta <- data.frame(
    Feature_ID = feature_ids,
    row.names = feature_ids
  )
  
  message("✓ Loaded ", nrow(abundance_matrix), " features across ", 
          ncol(abundance_matrix), " samples")
  
} else {
  # Option 2: Load from mzML/mzXML files using notame's built-in functions
  # This requires peak-picked data from tools like xcms
  
  message("Feature table not found. Looking for mzML/mzXML files in: ", DATA_PATH)
  
  # List available files
  mzml_files <- list.files(DATA_PATH, pattern = "\\.(mzML|mzXML)$", 
                           full.names = TRUE, ignore.case = TRUE)
  
  if (length(mzml_files) > 0) {
    message("Found ", length(mzml_files), " LC-MS files")
    
    # Use notame's import functions
    # Note: This assumes data has been pre-processed with peak picking
    # Adjust based on your specific data format
    
    stop("Direct import from mzML/mzXML not yet implemented in this script.\n",
         "Please either:\n",
         "1. Provide a pre-processed feature table as '", FEATURE_TABLE, "', or\n",
         "2. Implement mzML import using notame::read_from_excel() or similar")
    
  } else {
    stop("No feature data found. Please provide either:\n",
         "  - A feature table at: ", FEATURE_TABLE, "\n",
         "  - LC-MS files (mzML/mzXML) in: ", DATA_PATH)
  }
}

# =============================================================================
# CREATE METABOSET OBJECT
# =============================================================================

message("\nCreating MetaboSet object...")

# Ensure sample names match between metadata and feature data
sample_names_meta <- metadata[[SAMPLE_ID_COLUMN]]
sample_names_data <- colnames(abundance_matrix)

# Check for mismatches
if (!all(sample_names_data %in% sample_names_meta)) {
  warning("Some samples in feature data not found in metadata")
}

# Align metadata with feature data
metadata_aligned <- metadata %>%
  filter(.data[[SAMPLE_ID_COLUMN]] %in% sample_names_data) %>%
  arrange(match(.data[[SAMPLE_ID_COLUMN]], sample_names_data))

# Convert metadata to data frame with proper row names
pdata <- as.data.frame(metadata_aligned)
rownames(pdata) <- pdata[[SAMPLE_ID_COLUMN]]

# Create feature metadata
fdata <- as.data.frame(feature_meta)

# Create ExpressionSet
eset <- ExpressionSet(
  assayData = abundance_matrix[, pdata[[SAMPLE_ID_COLUMN]]],
  phenoData = AnnotatedDataFrame(pdata),
  featureData = AnnotatedDataFrame(fdata)
)

# Create MetaboSet from ExpressionSet
metaboset <- construct_metabosets(
  exprs = Biobase::exprs(eset),
  pheno_data = Biobase::pData(eset),
  feature_data = Biobase::fData(eset),
  group_col = GROUP_COLUMN
)

message("✓ Created MetaboSet object")
message("  - ", nrow(metaboset), " features")
message("  - ", ncol(metaboset), " samples")

# =============================================================================
# INITIAL QUALITY CHECKS
# =============================================================================

message("\nInitial data summary:")

# Feature statistics
n_features <- nrow(metaboset)
n_samples <- ncol(metaboset)

message("  Features: ", n_features)
message("  Samples: ", n_samples)

# Missing value statistics
missing_pct <- sum(is.na(exprs(metaboset))) / (n_features * n_samples) * 100
message("  Missing values: ", round(missing_pct, 2), "%")

# Sample group distribution
group_counts <- table(pData(metaboset)[[GROUP_COLUMN]])
message("\n  Sample distribution by group:")
for (grp in names(group_counts)) {
  message("    ", grp, ": ", group_counts[grp])
}

# =============================================================================
# SAVE METABOSET OBJECT
# =============================================================================

message("\nSaving MetaboSet object...")

# Save as RData file
save_path <- file.path(OBJECTS_PATH, "metaboset_raw.RData")
save(metaboset, file = save_path)

message("✓ Saved to: ", save_path)

# Also save a summary report
summary_file <- file.path(TABLES_PATH, "01_import_summary.txt")
sink(summary_file)
cat("=== Data Import Summary ===\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Input files:\n")
cat("  Metadata:", METADATA_PATH, "\n")
cat("  Features:", FEATURE_TABLE, "\n\n")
cat("Data dimensions:\n")
cat("  Features:", n_features, "\n")
cat("  Samples:", n_samples, "\n")
cat("  Missing values:", round(missing_pct, 2), "%\n\n")
cat("Sample groups:\n")
print(group_counts)
cat("\n")
sink()

message("✓ Summary saved to: ", summary_file)

message("\n=== Data Import Complete ===\n")
message("Next step: Run 02_qc_filtering.R")
