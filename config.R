# config.R
# Configuration file for notame metabolomics workflow
# Modify these parameters according to your experiment

# =============================================================================
# FILE PATHS
# =============================================================================

# Input data paths
DATA_PATH <- "data/raw"                    # Path to raw LC-MS data
METADATA_PATH <- "data/metadata/sample_metadata.csv"  # Sample metadata file
FEATURE_TABLE <- "data/raw/feature_table.csv"         # Feature abundance table (if pre-processed)

# Output paths
RESULTS_PATH <- "results"                  # Main results directory
PLOTS_PATH <- file.path(RESULTS_PATH, "plots")
TABLES_PATH <- file.path(RESULTS_PATH, "tables")
OBJECTS_PATH <- file.path(RESULTS_PATH, "objects")

# Create output directories if they don't exist
dir.create(RESULTS_PATH, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_PATH, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLES_PATH, showWarnings = FALSE, recursive = TRUE)
dir.create(OBJECTS_PATH, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# EXPERIMENT PARAMETERS
# =============================================================================

# Sample information
SAMPLE_ID_COLUMN <- "Sample_ID"           # Column name for sample identifiers
GROUP_COLUMN <- "Group"                    # Column name for experimental groups
BATCH_COLUMN <- "Batch"                    # Column name for batch information (optional)
QC_COLUMN <- "QC"                          # Column indicating QC samples

# QC sample labels
QC_LABEL <- "QC"                           # Label for quality control samples
BLANK_LABEL <- "Blank"                     # Label for blank samples

# =============================================================================
# QUALITY CONTROL PARAMETERS
# =============================================================================

# Feature filtering
MIN_DETECTION_RATE <- 0.8                  # Minimum detection rate in QC samples (0-1)
MIN_FRACTION_SAMPLES <- 0.5                # Minimum fraction of samples with detection (0-1)
RSD_THRESHOLD <- 30                        # Maximum RSD% in QC samples

# Sample filtering
MIN_FEATURES_PER_SAMPLE <- 0.5             # Minimum fraction of features detected per sample

# =============================================================================
# NORMALIZATION PARAMETERS
# =============================================================================

# Normalization method
# Options: "pqn", "median", "sum", "quantile", "vsn", "none"
NORMALIZATION_METHOD <- "pqn"              # Probabilistic Quotient Normalization

# Transformation
# Options: "log", "sqrt", "none"
TRANSFORMATION <- "log"                    # Log transformation

# =============================================================================
# STATISTICAL ANALYSIS PARAMETERS
# =============================================================================

# Statistical tests
ALPHA <- 0.05                              # Significance threshold
ADJUST_METHOD <- "fdr"                     # P-value adjustment method (fdr, bonferroni, holm, etc.)
FC_THRESHOLD <- 1.5                        # Fold change threshold for volcano plots

# PCA parameters
N_PCA_COMPONENTS <- 5                      # Number of principal components to calculate

# =============================================================================
# VISUALIZATION PARAMETERS
# =============================================================================

# Plot dimensions
PLOT_WIDTH <- 10                           # inches
PLOT_HEIGHT <- 8                           # inches
DPI <- 300                                 # Resolution for saved plots

# Color schemes
COLOR_PALETTE <- "Set1"                    # RColorBrewer palette name

# Heatmap parameters
N_TOP_FEATURES_HEATMAP <- 50               # Number of top variable features for heatmap

# =============================================================================
# ADVANCED PARAMETERS
# =============================================================================

# Missing value imputation
IMPUTATION_METHOD <- "knn"                 # k-nearest neighbors imputation
K_NEIGHBORS <- 5                           # Number of neighbors for KNN

# Parallel processing
N_CORES <- 1                               # Number of CPU cores to use (1 = no parallelization)

# Random seed for reproducibility
SEED <- 42
set.seed(SEED)

# =============================================================================
# PRINT CONFIGURATION
# =============================================================================

message("=== Workflow Configuration Loaded ===")
message("Data path: ", DATA_PATH)
message("Results path: ", RESULTS_PATH)
message("Normalization: ", NORMALIZATION_METHOD)
message("Transformation: ", TRANSFORMATION)
message("Alpha: ", ALPHA)
message("=====================================")
