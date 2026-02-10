# Detailed Workflow Guide

This document provides detailed information about each step of the metabolomics preprocessing workflow.

## Table of Contents

1. [Workflow Overview](#workflow-overview)
2. [Detailed Step Descriptions](#detailed-step-descriptions)
3. [Understanding the Output](#understanding-the-output)
4. [Best Practices](#best-practices)
5. [Advanced Usage](#advanced-usage)

## Workflow Overview

The metabolomics preprocessing workflow consists of five sequential steps:

```
Raw Data → Import → QC Filtering → Normalization → Visualization → Statistics → Results
```

Each step builds upon the previous one, transforming raw LC-MS data into meaningful biological insights.

## Detailed Step Descriptions

### Step 1: Data Import

**Purpose:** Load raw LC-MS data and create a structured MetaboSet object.

**What it does:**
1. Reads feature abundance table (CSV format)
2. Reads sample metadata (CSV format)
3. Validates data consistency
4. Creates a MetaboSet object (notame's core data structure)
5. Generates initial data summary

**Key checks:**
- Sample IDs match between data and metadata
- Required metadata columns are present
- No duplicate feature IDs
- Basic statistics (missing values, dimensions)

**Output files:**
- `metaboset_raw.RData` - Raw MetaboSet object
- `01_import_summary.txt` - Import statistics

**Common issues:**
- **Sample mismatch:** Ensure sample IDs are identical in both files
- **Missing columns:** Check that metadata has required columns
- **File encoding:** Use UTF-8 encoding for CSV files

---

### Step 2: Quality Control and Filtering

**Purpose:** Remove low-quality features and samples to improve data reliability.

**Feature filtering criteria:**

1. **Detection rate in QC samples** (`MIN_DETECTION_RATE`)
   - Removes features not consistently detected in QC samples
   - Default: 80% detection required
   - Rationale: Features with low QC detection are unreliable

2. **Detection rate across all samples** (`MIN_FRACTION_SAMPLES`)
   - Removes features detected in too few samples
   - Default: 50% detection required
   - Rationale: Rare features may be noise or contaminants

3. **Relative Standard Deviation (RSD) in QC** (`RSD_THRESHOLD`)
   - Removes features with high technical variation
   - Default: RSD ≤ 30%
   - Rationale: High RSD indicates poor measurement reproducibility

**Sample filtering:**

- Removes samples with low feature detection rates
- Default: Samples must detect ≥50% of features
- Identifies outliers or failed injections

**Visualization:**
- Feature detection distribution histogram
- Sample detection distribution histogram
- RSD distribution in QC samples

**Output files:**
- `metaboset_filtered.RData` - Filtered MetaboSet object
- `02_filtering_summary.txt` - Filtering statistics
- QC plots in `results/plots/qc/`

**Tuning parameters:**
- **Strict filtering:** Increase thresholds (e.g., 90% detection, RSD ≤ 20%)
- **Lenient filtering:** Decrease thresholds (e.g., 60% detection, RSD ≤ 40%)
- Balance between data quality and feature retention

---

### Step 3: Normalization and Transformation

**Purpose:** Correct for technical variation and prepare data for statistical analysis.

**Missing value imputation:**

Options:
- `knn` - K-nearest neighbors (default)
  - Uses similar samples to estimate missing values
  - Preserves data structure
  
- `half_min` - Half minimum
  - Replaces NA with half the minimum detected value
  - Conservative approach, assumes below detection limit
  
- `none` - No imputation
  - Keeps missing values (some analyses may fail)

**Normalization methods:**

1. **PQN - Probabilistic Quotient Normalization** ⭐ Recommended
   - Robust to dilution effects
   - Calculates quotients relative to reference (median spectrum)
   - Normalizes using median quotient per sample
   - Best for biological samples with variable dilution

2. **Median Normalization**
   - Scales samples to have equal median intensity
   - Simple and fast
   - Assumes similar overall metabolite abundance

3. **Sum Normalization**
   - Scales to equal total intensity
   - Sensitive to high-abundance features
   - Use with caution if few features dominate

4. **Quantile Normalization**
   - Forces identical intensity distributions
   - May over-normalize biological variation
   - Better for technical replicates than biological samples

**Transformation:**

- **Log transformation** (recommended)
  - Reduces heteroscedasticity
  - Makes data more normally distributed
  - Appropriate for most metabolomics data
  
- **Square root transformation**
  - Less aggressive than log
  - Alternative for count-like data
  
- **No transformation**
  - Keep original scale
  - May violate statistical test assumptions

**Visualization:**
- Boxplots before/after normalization
- Density plots showing distribution changes
- Total intensity per sample

**Output files:**
- `metaboset_normalized.RData` - Normalized MetaboSet
- `normalized_data.csv` - Normalized feature matrix
- `03_normalization_summary.txt` - Normalization details
- Normalization plots in `results/plots/normalization/`

**Choosing normalization:**
- Start with PQN + log transformation (default)
- Check normalization plots - distributions should align
- If batch effects persist, consider batch correction methods

---

### Step 4: Visualization and Exploratory Analysis

**Purpose:** Explore data structure, identify patterns, and assess quality.

**Principal Component Analysis (PCA):**

- Reduces dimensionality while preserving variance
- Identifies major sources of variation
- Detects outliers and batch effects
- Shows sample clustering by experimental groups

**Key PCA outputs:**

1. **Score plots** (PC1 vs PC2, PC2 vs PC3)
   - Each point = one sample
   - Colored by experimental group
   - Close points = similar metabolic profiles
   - Expected: Samples cluster by biological group

2. **Scree plot**
   - Shows variance explained by each PC
   - First 2-3 PCs typically capture most variation
   - Expected: PC1 > PC2 > PC3...

3. **Loadings**
   - Shows feature contributions to PCs
   - Identifies features driving separation
   - Saved in `pca_loadings.csv`

**Interpretation:**

Good quality data:
- Samples cluster by experimental group
- QC samples cluster tightly together
- First 2 PCs explain >30% variance
- No strong batch effects

Red flags:
- Samples don't cluster by group (weak biological effect)
- QC samples spread out (technical variation)
- Batch effects dominate PC1 (need batch correction)
- Single outlier samples (potential failed injections)

**Other visualizations:**

1. **Heatmap of top variable features**
   - Shows hierarchical clustering
   - Features (rows) × Samples (columns)
   - Color intensity = scaled abundance
   - Identifies feature patterns

2. **Sample correlation matrix**
   - Pearson correlation between all sample pairs
   - Expected: High correlation within groups
   - Identifies outliers (low correlation to others)

3. **Coefficient of Variation (CV) distribution**
   - Shows feature variability
   - Lower CV = more reproducible
   - Typical range: 10-50% for metabolomics

**Output files:**
- PCA plots (PDF and PNG)
- Heatmaps and correlation matrices
- `pca_scores.csv`, `pca_loadings.csv`, `pca_variance.csv`
- `04_visualization_summary.txt`

---

### Step 5: Statistical Analysis

**Purpose:** Identify significantly different metabolites between groups.

**Analysis types:**

1. **Pairwise comparison (2 groups)**
   - T-test for each feature
   - Calculates fold change (FC)
   - Generates volcano and MA plots
   
2. **Multiple groups (>2 groups)**
   - One-way ANOVA
   - Identifies features varying across groups
   - Post-hoc tests for pairwise comparisons

**Statistical workflow:**

1. For each feature:
   - Calculate group means
   - Perform statistical test (t-test or ANOVA)
   - Calculate fold change (2 groups only)

2. Multiple testing correction:
   - Adjusts p-values for many comparisons
   - Methods: FDR (Benjamini-Hochberg), Bonferroni, Holm
   - Default: FDR (less conservative, higher power)

3. Apply significance criteria:
   - Adjusted p-value < α (default: 0.05)
   - Absolute fold change > threshold (default: 1.5)
   - Both criteria must be met

**Volcano Plot:**

- X-axis: Log2 fold change (effect size)
- Y-axis: -Log10 adjusted p-value (significance)
- Red points: Up-regulated (higher in group 2)
- Blue points: Down-regulated (lower in group 2)
- Gray points: Not significant

**MA Plot:**

- X-axis: Average expression (A)
- Y-axis: Log2 fold change (M)
- Similar information to volcano plot
- Better shows relationship between expression level and fold change

**Output files:**
- `statistical_results.csv` - All features with statistics
- `significant_features.csv` - Only significant features
- Volcano and MA plots
- `05_statistical_summary.txt`

**Interpreting results:**

Columns in results table:
- `Feature_ID` - Feature identifier
- `Mean_Group1`, `Mean_Group2` - Average abundance per group
- `Fold_Change` - Ratio of group means (Group2/Group1)
- `Log2_FC` - Log2-transformed fold change
- `P_Value` - Raw p-value from statistical test
- `P_Adjusted` - Multiple testing corrected p-value
- `Significant` - Meets both p-value and FC criteria
- `Regulation` - Up/Down/Not Significant

---

## Understanding the Output

### Directory Structure

```
results/
├── plots/                      # All visualizations
│   ├── qc/                    # Quality control plots
│   │   ├── feature_detection_distribution.pdf
│   │   ├── sample_detection_distribution.pdf
│   │   └── qc_rsd_distribution.pdf
│   ├── normalization/          # Normalization assessment
│   │   ├── boxplot_comparison.pdf
│   │   ├── density_comparison.pdf
│   │   └── total_intensity.pdf
│   ├── visualization/          # Exploratory analysis
│   │   ├── pca_pc1_pc2.pdf
│   │   ├── pca_scree_plot.pdf
│   │   ├── heatmap_top_features.pdf
│   │   └── sample_correlation.pdf
│   └── statistics/             # Statistical results
│       ├── volcano_plot.pdf
│       └── ma_plot.pdf
├── tables/                     # CSV result files
│   ├── 01_import_summary.txt
│   ├── 02_filtering_summary.txt
│   ├── 03_normalization_summary.txt
│   ├── 04_visualization_summary.txt
│   ├── 05_statistical_summary.txt
│   ├── normalized_data.csv
│   ├── pca_scores.csv
│   ├── pca_loadings.csv
│   ├── statistical_results.csv
│   └── significant_features.csv
└── objects/                    # Saved R objects
    ├── metaboset_raw.RData
    ├── metaboset_filtered.RData
    └── metaboset_normalized.RData
```

### Key Result Files

1. **normalized_data.csv**
   - Final processed feature intensities
   - Ready for downstream analysis
   - Can be imported into other tools (MetaboAnalyst, etc.)

2. **statistical_results.csv**
   - Complete statistical analysis results
   - All features, sorted by significance
   - Use for comprehensive review

3. **significant_features.csv**
   - Only features meeting significance criteria
   - Primary candidates for biological interpretation
   - Use for pathway analysis, annotation

4. **pca_scores.csv**
   - Sample coordinates in PCA space
   - Can be used for custom plotting
   - Includes metadata for coloring

---

## Best Practices

### Data Preparation

1. **Quality control samples**
   - Include QC samples (pooled samples) throughout acquisition
   - Minimum: 1 QC per 10-15 samples
   - Use for RSD filtering and normalization assessment

2. **Blank samples**
   - Include blank samples (extraction without biological material)
   - Identify background contamination
   - Remove or subtract from biological signals

3. **Randomization**
   - Randomize sample injection order
   - Reduces systematic drift effects
   - Record injection order in metadata

### Parameter Selection

1. **Start with defaults**
   - Default parameters work well for most datasets
   - Adjust based on your data characteristics

2. **QC filtering**
   - Stricter for high-quality instruments (RSD ≤ 20%)
   - More lenient for older instruments (RSD ≤ 30-40%)
   - Consider your sample size (strict filtering reduces features)

3. **Statistical thresholds**
   - α = 0.05 is standard
   - For discovery studies: FDR correction (less conservative)
   - For targeted validation: Bonferroni correction (more conservative)

### Workflow Execution

1. **Run sequentially**
   - Scripts must be run in order (01 → 05)
   - Each depends on previous outputs

2. **Check intermediate results**
   - Review summary files after each step
   - Examine QC plots before proceeding
   - Verify expected number of features/samples retained

3. **Save intermediate objects**
   - MetaboSet objects are saved at each step
   - Can restart from any point without re-running earlier steps

### Data Interpretation

1. **PCA before statistics**
   - Always examine PCA first
   - If groups don't separate, statistical tests may find little
   - Check for outliers, batch effects, QC clustering

2. **Volcano plot assessment**
   - Expect roughly symmetric distribution (equal up/down)
   - Skewed distribution may indicate normalization issues
   - Check for outlier features with extreme fold changes

3. **Significance criteria**
   - Use both p-value AND fold change
   - P-value alone: statistically significant but biologically irrelevant
   - Fold change alone: large effect but unreliable

---

## Advanced Usage

### Custom Normalization

You can implement custom normalization methods in `03_normalization.R`:

```r
# Example: Internal standard normalization
is_feature <- "M123.45T67.8"  # Your internal standard
is_values <- expr_matrix[is_feature, ]
norm_factors <- is_values / median(is_values)
expr_normalized <- sweep(expr_matrix, 2, norm_factors, "/")
```

### Batch Correction

For datasets with batch effects, add drift correction after filtering:

```r
# In 02_qc_filtering.R, before saving
if ("Injection_Order" %in% colnames(pData(metaboset_filtered))) {
  metaboset_filtered <- correct_drift(metaboset_filtered)
}
```

### Post-hoc Tests (ANOVA)

For >2 groups, add pairwise comparisons after ANOVA:

```r
# In 05_statistical_analysis.R
for (i in 1:nrow(expr_data)) {
  if (results$Significant[i]) {
    posthoc <- pairwise.t.test(
      expr_data[i, ], 
      sample_data[[GROUP_COLUMN]], 
      p.adjust.method = ADJUST_METHOD
    )
    # Save posthoc results
  }
}
```

### Annotation Integration

Link features to metabolite identities:

```r
# After statistical analysis
# Merge with annotation database
annotation <- read_csv("annotation_database.csv")
results_annotated <- left_join(
  results, 
  annotation, 
  by = "Feature_ID"
)
```

### Pathway Analysis

Export significant features for pathway analysis:

```r
# Extract significant feature list
significant_ids <- significant_features$Feature_ID

# Export for MetaboAnalyst
write.table(
  significant_ids, 
  "metaboanalyst_input.txt", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)
```

---

## Troubleshooting

### Issue: Too many features removed

**Cause:** Filtering thresholds too strict

**Solution:**
- Lower `MIN_DETECTION_RATE` (e.g., 0.7 instead of 0.8)
- Lower `MIN_FRACTION_SAMPLES` (e.g., 0.4 instead of 0.5)
- Increase `RSD_THRESHOLD` (e.g., 40 instead of 30)

### Issue: PCA shows no group separation

**Cause:** Weak biological effect or technical issues

**Solution:**
1. Check if technical variation dominates (batch effects)
2. Verify sample labeling is correct
3. Consider stricter QC filtering
4. May indicate true biological similarity (not always a problem)

### Issue: Volcano plot heavily skewed

**Cause:** Normalization problem or outlier samples

**Solution:**
1. Review normalization plots
2. Try different normalization method
3. Check for outlier samples in PCA
4. Verify correct group assignments

### Issue: No significant features

**Cause:** 
- Insufficient sample size
- Weak biological effect
- Too stringent thresholds

**Solution:**
- Check if p-values are generally low (power issue) vs high (no effect)
- Increase sample size if possible
- Lower FC threshold (e.g., 1.2 instead of 1.5)
- Use FDR instead of Bonferroni

---

## Additional Resources

- [notame vignettes](https://bioconductor.org/packages/release/bioc/vignettes/notame/inst/doc/notame.html)
- [Metabolomics Society guidelines](http://metabolomicssociety.org/)
- [MetaboAnalyst tutorial](https://www.metaboanalyst.ca/docs/Tutorial.xhtml)
- [PCA interpretation guide](https://towardsdatascience.com/pca-clearly-explained-how-when-why-to-use-it-and-feature-importance-a-guide-in-python-7c274582c37e)

---

**Last updated:** 2026-02-10
