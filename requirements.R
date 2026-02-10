# requirements.R
# Install required packages for notame metabolomics workflow
# Run this script once to set up your environment

# Function to install packages if not already installed
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing", pkg, "..."))
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
  } else {
    message(paste(pkg, "already installed"))
  }
}

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Core packages from CRAN
cran_packages <- c(
  "tidyverse",    # Data manipulation and visualization
  "readxl",       # Read Excel files
  "ggplot2",      # Plotting
  "cowplot",      # Plot arrangement
  "RColorBrewer", # Color palettes
  "pheatmap",     # Heatmaps
  "car",          # Statistical tests
  "Hmisc"         # Data analysis tools
)

# Core packages from Bioconductor
bioc_packages <- c(
  "notame",       # Main workflow package
  "Biobase",      # Base structures for bioinformatics
  "pcaMethods",   # PCA analysis
  "impute"        # Missing value imputation
)

# Install CRAN packages
message("\n=== Installing CRAN packages ===")
for (pkg in cran_packages) {
  install_if_missing(pkg, bioc = FALSE)
}

# Install Bioconductor packages
message("\n=== Installing Bioconductor packages ===")
for (pkg in bioc_packages) {
  install_if_missing(pkg, bioc = TRUE)
}

# Verify installation
message("\n=== Verifying installation ===")
all_packages <- c(cran_packages, bioc_packages)
missing_packages <- all_packages[!sapply(all_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) == 0) {
  message("\n✓ All packages installed successfully!")
} else {
  warning("\n✗ The following packages failed to install:")
  print(missing_packages)
}

# Print session info
message("\n=== Session Info ===")
sessionInfo()
