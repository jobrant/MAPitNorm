# MAPitNorm
Normalization of accessibility methylation data from genome-wide, single molecule, MAPit sequence data

# MAPitNorm: Acessibility Methylation Data Normalization Package

## Overview
MAPitNorm is a specialized R package designed for normalizing methylation data from genome-wide MAPit experiments, such as MAPit-CpGiant or MAPit-RRBS. It implements a sophisticated two-step normalization process that addresses both technical variation between replicates and systematic differences between experimental groups while preserving biological signals.

## Core Functions

### Data Loading and Preparation
- `load_data()`: Efficiently loads and processes methylation call files; optionally with multiple cores
  - Handles both GCH and HCG methylation contexts
  - Performs initial data validation and organization
  - Returns a structured list of samples grouped by type

### Normalization Process

#### 1. Coverage Normalization (`normalize_coverage`)
- Normalizes coverage depth across samples
- Features:
  - Within-group normalization for technical replicates
  - Optional between-group normalization
  - Scaling factor calculation based on total coverage
  - Preserves relative coverage patterns

#### 2. Methylation Rate Normalization (`normalize_methylation_rates`)
- Implements advanced quantile-based normalization with structure preservation
- Key features:
  - Dynamic quantile binning based on dataset size
  - Structure preservation using alpha parameter (default = 0.3)
  - Two-step process:
    - Within-group normalization (technical replicates)
    - Optional between-group normalization (experimental conditions)
  - Preserves biological differences while reducing technical variation

### Visualization and Diagnostics

#### Within-Group Report (`plot_within_group_report`)
- Focused on technical replicate normalization
- Includes:
  - QQ plots comparing replicates
  - Density distributions
  - Multiple visualizations per experimental group

#### Between-Group Report (`plot_between_group_report`)
- Evaluates group-level normalization effects
- Features:
  - Group mean comparisons
  - Distribution plots
  - Box plots showing group-level changes

### Key Features
1. **Efficient Processing**
   - Optimized for large datasets
   - Uses data.table for fast operations
   - Memory-efficient visualization using hexbin plots

2. **Flexible Normalization**
   - Adjustable parameters for different experimental needs
   - Structure preservation to maintain biological signals
   - Independent control of within and between-group normalization

3. **Comprehensive Visualization**
   - Dense data visualization using hexbin plots
   - Separate reports for within and between-group effects
   - Customizable plotting parameters

4. **Quality Control**
   - Diagnostic information during normalization
   - Visual assessment of normalization effects
   - Multiple metrics for evaluating normalization success

## Technical Details

### Normalization Parameters
- `sites_per_quantile`: Controls quantile bin size (default = 1000)
- `within_alpha`: Structure preservation for technical replicates (default = 0.3)
- `between_alpha`: Structure preservation for group differences (default = 0.5)
- Dynamic quantile calculation based on dataset size

### Visualization Options
- Customizable color schemes (viridis color palettes)
- Adjustable bin sizes for density plots
- Optional data sampling for quick visualization
- PDF export capability for reports

## Future Developments
- Additional normalization methods
- Enhanced visualization options
- Performance optimizations for very large datasets
- Extended diagnostic capabilities

## Usage
The package is designed for methylation data analysis workflows, particularly suited for:
- MAPit experimental data
- Multi-group experimental designs
- Datasets with technical replicates
- Large-scale methylation studies

This implementation provides a robust framework for normalizing methylation data while maintaining biological significance and offering comprehensive diagnostic tools for quality assessment.
