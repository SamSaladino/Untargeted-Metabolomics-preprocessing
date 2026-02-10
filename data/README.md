# Data Directory

This directory contains input data for the metabolomics preprocessing workflow.

## Directory Structure

```
data/
├── raw/                    # Raw LC-MS data files
│   └── feature_table.csv  # Feature abundance matrix (required)
└── metadata/              # Sample metadata
    └── sample_metadata.csv  # Sample information (required)
```

## Required Files

### 1. Feature Table (`raw/feature_table.csv`)

This file contains the abundance matrix from your LC-MS analysis.

**Format:**
- **Rows:** Features (metabolites/compounds)
- **Columns:** Samples
- **First column:** Feature identifiers
- **Remaining columns:** Sample intensities

**Example:**

```csv
Feature_ID,Sample_001,Sample_002,Sample_003,Sample_004,QC_001
M100.0523T45.2,15234,18456,12789,16234,15123
M150.0841T67.8,45678,52341,48901,49876,48234
M200.1234T89.3,98234,87654,92345,95678,91234
```

**Important:**
- Column names (sample IDs) must match exactly with `Sample_ID` in metadata
- Missing values should be represented as `NA` or left empty
- Feature IDs should be unique
- Common feature ID formats:
  - `M{mz}T{rt}` (e.g., M123.4567T89.12)
  - `Compound_{number}` (e.g., Compound_001)
  - HMDB or other database IDs

**How to generate:**
- Export from your LC-MS processing software (e.g., XCMS, MZmine, Progenesis)
- Ensure proper peak picking and alignment have been performed
- Include only peak areas/heights, not metadata columns

---

### 2. Sample Metadata (`metadata/sample_metadata.csv`)

This file contains information about each sample.

**Required columns:**
- `Sample_ID` - Unique identifier for each sample (must match feature table columns)
- `Group` - Experimental group/condition (e.g., "Control", "Treatment")

**Optional columns:**
- `Batch` - Batch or run number (for batch effect assessment)
- `QC` - Sample type indicator (e.g., "Sample", "QC", "Blank")
- `Injection_Order` - Order of sample injection (for drift correction)
- Any other metadata you want to track (age, gender, time point, etc.)

**Example:**

```csv
Sample_ID,Group,Batch,QC,Injection_Order,Time_Point
Sample_001,Control,1,Sample,1,0h
Sample_002,Control,1,Sample,3,0h
Sample_003,Control,1,Sample,5,0h
Sample_004,Treatment,1,Sample,7,24h
Sample_005,Treatment,1,Sample,9,24h
Sample_006,Treatment,1,Sample,11,24h
QC_001,QC,1,QC,6,NA
QC_002,QC,1,QC,12,NA
Blank_001,Blank,1,Blank,2,NA
```

**Important:**
- Sample_ID values must be unique
- Sample_ID must match exactly with column names in feature table
- Group should have at least 2 different values (for comparisons)
- QC samples are highly recommended (for quality assessment)
- Include blank samples if available (for background subtraction)

---

## Data Preparation Checklist

Before running the workflow, ensure:

- [ ] Feature table has samples as columns (not rows)
- [ ] First column of feature table contains unique feature IDs
- [ ] Sample IDs match exactly between feature table and metadata
- [ ] No special characters in sample IDs (use letters, numbers, underscore)
- [ ] Metadata has required columns: `Sample_ID` and `Group`
- [ ] Group column contains descriptive names (not just 1/2)
- [ ] QC samples are labeled consistently (recommended)
- [ ] Files are saved as CSV with UTF-8 encoding
- [ ] No empty rows or columns at the end of files

---

## Tips for Data Quality

### Quality Control Samples

**What are QC samples?**
- Pooled samples created by mixing aliquots from all study samples
- Injected repeatedly throughout the analytical batch
- Used to assess technical reproducibility

**Recommendations:**
- Include 1 QC sample per 10-15 study samples
- Inject QC samples regularly throughout the run
- Inject 3-5 QC samples at the beginning (for column conditioning)

### Blank Samples

**What are blank samples?**
- Samples processed identically to study samples but without biological material
- Used to identify background contamination

**Recommendations:**
- Include at least 2-3 blank samples per batch
- Process blanks at the beginning, middle, and end of batch
- Use same solvents and extraction procedure as study samples

### Randomization

**Why randomize?**
- Prevents systematic bias
- Distributes technical drift across all groups
- Improves statistical validity

**Recommendations:**
- Randomize injection order within batches
- If studying time points, randomize within each time point
- Record injection order in metadata for drift correction

---

## Common Data Formats

### From XCMS (R)

```r
# Export feature table from XCMS
write.csv(featureValues(xset), "feature_table.csv")

# Export metadata
write.csv(phenoData(xset), "sample_metadata.csv")
```

### From MZmine

1. Export → CSV file
2. Choose "Export for nontargeted" or similar
3. Transpose if needed (samples should be columns)

### From Progenesis QI

1. Export → Compound measurements
2. Choose "All compounds" and "Normalized abundance"
3. Format: CSV with samples as columns

### From MetaboAnalyst

If you processed data in MetaboAnalyst:
1. Download the normalized data table
2. Format should already be compatible (features × samples)
3. Create metadata file separately

---

## Example Datasets

A template metadata file is provided at `data/metadata/sample_metadata.csv`.

For testing the workflow with real data:
1. Replace the template with your actual metadata
2. Add your feature table to `data/raw/feature_table.csv`
3. Update `config.R` if your column names differ

---

## File Size Considerations

**Large datasets:**
- For >10,000 features or >500 samples, processing may be slow
- Consider increasing memory allocation in R: `options(java.parameters = "-Xmx8g")`
- Use high-performance computing if available

**Storage:**
- Raw data files are not tracked in git (see `.gitignore`)
- Store original data separately
- Backup processed results (`results/` directory)

---

## Troubleshooting

### "Sample names don't match"

Check for:
- Extra spaces in sample IDs
- Inconsistent capitalization (Sample_001 vs sample_001)
- Special characters (-, ., / may cause issues)
- Excel auto-formatting (dates, scientific notation)

**Solution:** Use plain text editor to verify exact matches

### "Cannot find feature_table.csv"

- Ensure file is in `data/raw/` directory
- Check file extension (.csv not .CSV)
- Verify path in `config.R` matches actual location

### "Missing required columns"

- Metadata must have `Sample_ID` and `Group` columns
- Column names are case-sensitive
- Check for typos in column headers

---

## Getting Help

If you encounter issues with data formatting:

1. Check that files match the examples above
2. Verify with a small test dataset first
3. Review error messages for specific issues
4. Consult the main README.md for general troubleshooting

---

**Note:** The workflow expects pre-processed data (peak picking and alignment completed). It does not perform raw spectrum processing from mzML/mzXML files.
