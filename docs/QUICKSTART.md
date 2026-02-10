# Quick Start Guide

Get started with the metabolomics preprocessing workflow in 5 minutes.

## 1. Install Dependencies

Open R or RStudio and run:

```r
source("requirements.R")
```

This will install all necessary packages from CRAN and Bioconductor.

## 2. Prepare Your Data

Place your files in the correct locations:

```
data/
├── raw/
│   └── feature_table.csv          # Your LC-MS data
└── metadata/
    └── sample_metadata.csv        # Your sample info
```

**Feature table format:**
- First column: Feature IDs
- Other columns: Sample abundances
- Samples as columns, features as rows

**Metadata format:**
- Required columns: `Sample_ID`, `Group`
- Sample_ID must match column names in feature table

See `data/README.md` for detailed format requirements.

## 3. Configure Settings (Optional)

Edit `config.R` to adjust parameters:

```r
# Most important settings:
GROUP_COLUMN <- "Group"              # Your group column name
NORMALIZATION_METHOD <- "pqn"        # Normalization method
ALPHA <- 0.05                        # Significance threshold
```

For first-time users: **keep default settings**.

## 4. Run the Workflow

Execute scripts in order:

```r
# Step 1: Import data
source("scripts/01_data_import.R")

# Step 2: Quality control
source("scripts/02_qc_filtering.R")

# Step 3: Normalization
source("scripts/03_normalization.R")

# Step 4: Visualization
source("scripts/04_visualization.R")

# Step 5: Statistical analysis
source("scripts/05_statistical_analysis.R")
```

**Or run all at once:**

```r
# Run complete workflow
source("scripts/01_data_import.R")
source("scripts/02_qc_filtering.R")
source("scripts/03_normalization.R")
source("scripts/04_visualization.R")
source("scripts/05_statistical_analysis.R")
```

## 5. Review Results

Check the `results/` directory:

**Key files to examine:**
1. `plots/visualization/pca_pc1_pc2.pdf` - Do samples cluster by group?
2. `plots/statistics/volcano_plot.pdf` - Significant features visualization
3. `tables/significant_features.csv` - List of differential metabolites
4. `tables/05_statistical_summary.txt` - Summary statistics

## Expected Runtime

For typical dataset (5,000 features, 50 samples):
- Step 1: ~10 seconds
- Step 2: ~30 seconds
- Step 3: ~20 seconds
- Step 4: ~1 minute
- Step 5: ~30 seconds

**Total: ~3 minutes**

Larger datasets may take longer.

## Troubleshooting

### Installation fails
- Update R to latest version (>= 4.0)
- Update Bioconductor: `BiocManager::install(version = "3.18")`

### "File not found" error
- Check file paths in `config.R`
- Verify files exist in `data/raw/` and `data/metadata/`

### "Sample names don't match"
- Ensure Sample_ID in metadata exactly matches feature table columns
- Check for spaces, special characters, or typos

### No significant features
- Normal if biological effect is weak
- Try lowering FC_THRESHOLD in `config.R`
- Increase sample size if possible

## Next Steps

After running the workflow:

1. **Examine PCA plots** - Understand your data structure
2. **Review significant features** - Identify candidates for validation
3. **Adjust parameters** - Re-run if needed with different settings
4. **Annotate features** - Link feature IDs to metabolite names
5. **Pathway analysis** - Use MetaboAnalyst or similar tools

## Learn More

- **Full documentation:** See `README.md`
- **Detailed workflow guide:** See `docs/WORKFLOW.md`
- **Data format help:** See `data/README.md`
- **notame package:** https://bioconductor.org/packages/notame/

## Example Workflow Session

```r
# Complete example session
# 1. Set working directory to project root
setwd("/path/to/Untargeted-Metabolomics-preprocessing")

# 2. Install packages (first time only)
source("requirements.R")

# 3. Run workflow
source("scripts/01_data_import.R")
source("scripts/02_qc_filtering.R")
source("scripts/03_normalization.R")
source("scripts/04_visualization.R")
source("scripts/05_statistical_analysis.R")

# 4. Check results
list.files("results/tables")        # List result files
list.files("results/plots")         # List plot files

# 5. Load specific results
results <- read.csv("results/tables/significant_features.csv")
head(results)

# 6. View PCA scores
pca <- read.csv("results/tables/pca_scores.csv")
head(pca)
```

---

**Ready to start?** Make sure you have your data files ready and run `source("requirements.R")`!
