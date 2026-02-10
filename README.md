# Untargeted Metabolomics Preprocessing Workflow

A comprehensive R-based workflow for processing untargeted metabolomics data using the [notame](https://bioconductor.org/packages/notame/) package from Bioconductor.

## Overview

This repository provides a complete, step-by-step workflow for preprocessing LC-MS metabolomics data, including:
- Data import and quality control
- Feature and sample filtering
- Normalization and transformation
- Exploratory data analysis (PCA, heatmaps)
- Statistical analysis and differential metabolite detection

## Features

✅ **Modular workflow** - Five sequential R scripts covering the complete preprocessing pipeline  
✅ **Flexible configuration** - Centralized config file for easy parameter adjustment  
✅ **Quality control** - Automated QC filtering based on detection rates and RSD  
✅ **Multiple normalization methods** - PQN, median, sum, quantile normalization  
✅ **Comprehensive visualization** - PCA plots, heatmaps, volcano plots, and more  
✅ **Statistical analysis** - T-tests, ANOVA, multiple testing correction  
✅ **Well-documented** - Clear code comments and detailed output summaries  

## Installation

### Prerequisites

- R (version ≥ 4.0.0)
- RStudio (recommended)

### Install Required Packages

Run the installation script to automatically install all required packages:

```r
source("requirements.R")
```

This will install:
- **notame** - Main metabolomics preprocessing package
- **tidyverse** - Data manipulation and visualization
- **Biobase** - Bioconductor core functionality
- **pcaMethods** - PCA analysis
- And several other dependencies

## Quick Start

### 1. Prepare Your Data

Place your data files in the appropriate directories:

```
data/
├── raw/
│   └── feature_table.csv          # Your LC-MS feature abundance table
└── metadata/
    └── sample_metadata.csv        # Sample information (see template)
```

**Feature table format:**
- Rows = features (metabolites)
- Columns = samples
- First column = feature IDs
- Remaining columns = sample abundances

**Metadata format:**
- Required columns: `Sample_ID`, `Group`
- Optional: `Batch`, `QC` (for quality control samples)
- See `data/metadata/sample_metadata.csv` for an example

### 2. Configure Parameters

Edit `config.R` to set your analysis parameters:

```r
# Key parameters to adjust:
GROUP_COLUMN <- "Group"              # Column name for experimental groups
MIN_DETECTION_RATE <- 0.8            # Minimum detection rate in QC samples
NORMALIZATION_METHOD <- "pqn"        # Normalization method
ALPHA <- 0.05                        # Significance threshold
```

### 3. Run the Workflow

Execute the scripts in order:

```r
# Step 1: Import data
source("scripts/01_data_import.R")

# Step 2: Quality control and filtering
source("scripts/02_qc_filtering.R")

# Step 3: Normalization
source("scripts/03_normalization.R")

# Step 4: Visualization and exploratory analysis
source("scripts/04_visualization.R")

# Step 5: Statistical analysis
source("scripts/05_statistical_analysis.R")
```

### 4. Review Results

All results are saved in the `results/` directory:

```
results/
├── plots/              # Visualizations (PDF and PNG)
│   ├── qc/            # Quality control plots
│   ├── normalization/ # Before/after normalization
│   ├── visualization/ # PCA, heatmaps
│   └── statistics/    # Volcano plots, MA plots
├── tables/            # CSV files with analysis results
│   ├── statistical_results.csv
│   ├── significant_features.csv
│   └── pca_scores.csv
└── objects/           # Saved R objects (.RData)
    ├── metaboset_raw.RData
    ├── metaboset_filtered.RData
    └── metaboset_normalized.RData
```

## Workflow Steps

### Script 1: Data Import (`01_data_import.R`)

- Loads LC-MS feature data and sample metadata
- Creates a MetaboSet object
- Performs initial data quality checks
- **Output:** `metaboset_raw.RData`

### Script 2: Quality Control (`02_qc_filtering.R`)

- Filters low-quality features based on:
  - Detection rate in QC samples
  - Detection rate across all samples  
  - RSD (Relative Standard Deviation) in QC samples
- Removes poor-quality samples
- Generates QC visualizations
- **Output:** `metaboset_filtered.RData`, QC plots

### Script 3: Normalization (`03_normalization.R`)

- Imputes missing values (KNN, half-minimum)
- Applies normalization (PQN, median, sum, quantile)
- Performs data transformation (log, sqrt)
- Generates before/after comparison plots
- **Output:** `metaboset_normalized.RData`, normalized data table

### Script 4: Visualization (`04_visualization.R`)

- Principal Component Analysis (PCA)
- Scree plots and score plots
- Heatmaps of top variable features
- Sample correlation matrices
- Coefficient of variation (CV) analysis
- **Output:** PCA plots, heatmaps, correlation matrices

### Script 5: Statistical Analysis (`05_statistical_analysis.R`)

- Differential analysis (t-test or ANOVA)
- Multiple testing correction (FDR, Bonferroni)
- Fold change calculation
- Volcano plots and MA plots
- **Output:** Statistical results tables, significance plots

## Configuration Options

Key parameters in `config.R`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NORMALIZATION_METHOD` | `"pqn"` | Normalization method (pqn, median, sum, quantile) |
| `TRANSFORMATION` | `"log"` | Data transformation (log, sqrt, none) |
| `MIN_DETECTION_RATE` | `0.8` | Minimum detection rate in QC samples (0-1) |
| `RSD_THRESHOLD` | `30` | Maximum RSD% in QC samples |
| `ALPHA` | `0.05` | Statistical significance threshold |
| `FC_THRESHOLD` | `1.5` | Fold change threshold for volcano plots |
| `ADJUST_METHOD` | `"fdr"` | P-value adjustment method |

## Input Data Requirements

### Feature Table

Your feature table should be a CSV file with:
- **First column:** Feature IDs (unique identifiers)
- **Remaining columns:** Sample abundances (one column per sample)
- Feature IDs can be in formats like: `M123T456`, `Compound_001`, etc.
- Missing values should be represented as `NA` or empty cells

Example:
```csv
Feature_ID,Sample_001,Sample_002,Sample_003
M123.45T67.8,15000,18000,12000
M234.56T78.9,25000,NA,23000
```

### Sample Metadata

Required columns:
- `Sample_ID`: Must match column names in feature table
- `Group`: Experimental groups for comparison (e.g., "Control", "Treatment")

Optional columns:
- `Batch`: Batch information for visualization
- `QC`: Quality control sample indicator (e.g., "QC", "Sample", "Blank")

## Normalization Methods

The workflow supports multiple normalization methods:

- **PQN (Probabilistic Quotient Normalization)** ⭐ Recommended
  - Robust to dilution effects
  - Good for biological samples
  
- **Median Normalization**
  - Simple and fast
  - Assumes similar overall abundance
  
- **Sum Normalization**
  - Normalizes to total intensity
  - Sensitive to high-abundance features
  
- **Quantile Normalization**
  - Forces same distribution across samples
  - May over-normalize biological variation

## Troubleshooting

### "Metadata file not found"
- Ensure `data/metadata/sample_metadata.csv` exists
- Check that the path in `config.R` is correct

### "Feature table not found"
- Place your feature table at `data/raw/feature_table.csv`
- Or update `FEATURE_TABLE` path in `config.R`

### "Sample names don't match"
- Verify that sample IDs in metadata match column names in feature table
- Check for typos or extra spaces

### Missing packages
- Re-run `source("requirements.R")`
- Manually install missing packages: `BiocManager::install("package_name")`

## Citation

If you use this workflow, please cite the notame package:

```
Klåvus A, Kokla M, Noerman S, Koistinen VM, Tuomainen M, Zarei I, Meuronen T, 
Häkkinen MR, Rummukainen S, Babu AF, Gomez-Gallego C, Schwab U, Auriola S, 
Hanhineva K, Koistinen V (2020). "notame: Workflow for Non-Targeted LC–MS 
Metabolic Profiling." Metabolites, 10(4), 135. doi: 10.3390/metabo10040135
```

## Contributing

Contributions are welcome! Please feel free to:
- Report issues
- Suggest improvements
- Submit pull requests

## License

This workflow is distributed under the terms included in the LICENSE file.

## Resources

- [notame package documentation](https://bioconductor.org/packages/notame/)
- [Bioconductor workflows](https://bioconductor.org/packages/release/BiocViews.html#___Workflow)
- [Metabolomics data analysis guide](https://www.metaboanalyst.ca/)

## Contact

For questions or issues, please open an issue on GitHub.
