# PGS-Compare

PGS-Compare is a Python package for analyzing and comparing Polygenic Scores (PGS) across ancestry groups. It uses the PGS Catalog and 1000 Genomes data to help researchers evaluate the stability of PGS scores across different ancestry groups.

## Features

- Download necessary data from the 1000 Genomes project and reference panels
- Fetch and calculate PGS scores for specific traits using the PGS Catalog
- Analyze PGS score distributions across different ancestry groups
- Compare consistency of PGS scores from different studies
- Visualize results with publication-ready plots
- Assess stability of PGS predictions across and within ancestry groups

## Prerequisites

The package relies on the following external tools:

1. [PLINK 2](https://www.cog-genomics.org/plink/2.0/) - For genetic data processing
2. [Nextflow](https://www.nextflow.io/) - For running the PGS Catalog Calculator

Make sure these tools are installed and available in your PATH before using PGS-Compare.

## Installation

Install the package from PyPI:

```bash
pip install pgs-compare
```

Or install directly from GitHub:

```bash
pip install git+https://github.com/yourusername/pgs-compare.git
```

## Getting Started

### Basic Usage

```python
from pgs_compare import PGSCompare

# Initialize with automatic dependency checking and data downloading
pgs = PGSCompare()

# Run the full pipeline for a specific trait
# Example: Parkinson's disease (MONDO_0005180)
results = pgs.run_pipeline("MONDO_0005180")

# The results include:
# - Calculation results (PGS scores)
# - Analysis results (summary statistics, correlations, etc.)
# - Visualization results (paths to plots)
```

### Command-line Interface

PGS-Compare also provides a command-line interface:

```bash
# Run calculations for Parkinson's disease
pgs-compare calculate MONDO_0005180

# Run calculations with custom options
pgs-compare calculate MONDO_0005180 --exclude-child-pgs --max-variants 1000000 --run-ancestry

# Analyze the results
pgs-compare analyze --trait-id MONDO_0005180

# Analyze with a specific scores file
pgs-compare analyze --trait-id MONDO_0005180 --scores-file path/to/scores.txt.gz

# Generate visualizations
pgs-compare visualize --trait-id MONDO_0005180 --show-error-bars

# Or run the full pipeline
pgs-compare pipeline MONDO_0005180 --show-error-bars

# Run the pipeline with custom options
pgs-compare pipeline MONDO_0005180 --exclude-child-pgs --max-variants 1000000 --run-ancestry --skip-visualize
```

## API Reference

### PGSCompare Class

The main class for interacting with the package.

```python
from pgs_compare import PGSCompare

pgs = PGSCompare(data_dir=None, download_data=True)
```

Parameters:

- `data_dir` (str, optional): Directory to store data. Default is "data" in the current directory.
- `download_data` (bool): Whether to download missing data during initialization. Defaults to True.
  If set to False, will still check for dependencies but won't download missing data.

Methods:

#### calculate

```python
pgs.calculate(trait_id, include_child_pgs=True, max_variants=None,
              run_ancestry=False, reference_panel=None)
```

Run PGS calculations for a specific trait.

Parameters:

- `trait_id` (str): Trait ID (e.g., "MONDO_0005180" for Parkinson's disease)
- `include_child_pgs` (bool): Whether to include child-associated PGS IDs
- `max_variants` (int, optional): Maximum number of variants to include in PGS
- `run_ancestry` (bool): Whether to run ancestry analysis
- `reference_panel` (str, optional): Path to reference panel for ancestry analysis.

Returns:

- dict: Information about the calculation including success status and output path

#### analyze

```python
pgs.analyze(trait_id=None, scores_file=None)
```

Analyze PGS scores across ancestry groups.

Parameters:

- `trait_id` (str, optional): Trait ID. Used for organizing output if provided.
- `scores_file` (str, optional): Path to the scores file (aggregated_scores.txt.gz).
  If None, will look in the standard location based on trait_id.

Returns:

- dict: Analysis results including summary statistics, correlations, and variance metrics

#### visualize

```python
pgs.visualize(trait_id=None, analysis_results=None, show_error_bars=False)
```

Visualize PGS analysis results.

Parameters:

- `trait_id` (str, optional): Trait ID. Used for organizing output if provided.
- `analysis_results` (dict, optional): Analysis results from analyze().
  If None, will try to load from the standard location based on trait_id.
- `show_error_bars` (bool): Whether to display error bars on plots. Default is False.

Returns:

- dict: Dictionary with paths to the generated plots

#### run_pipeline

```python
pgs.run_pipeline(trait_id, include_child_pgs=True, max_variants=None,
                run_ancestry=False, visualize=True, show_error_bars=False)
```

Run the full pipeline (calculate, analyze, visualize) for a specific trait.

Parameters:

- `trait_id` (str): Trait ID (e.g., "MONDO_0005180" for Parkinson's disease)
- `include_child_pgs` (bool): Whether to include child-associated PGS IDs
- `max_variants` (int, optional): Maximum number of variants to include in PGS
- `run_ancestry` (bool): Whether to run ancestry analysis
- `visualize` (bool): Whether to generate visualization plots
- `show_error_bars` (bool): Whether to display error bars on plots. Default is False.

Returns:

- dict: Pipeline results containing calculation, analysis, and visualization results

## Finding Trait IDs

You can find trait IDs by searching the [PGS Catalog](https://www.pgscatalog.org/). Some common traits:

- Parkinson's disease: `MONDO_0005180`
- Coronary artery disease: `EFO_0001645`
- Body height: `OBA_VT0001253`
- Breast cancer: `MONDO_0007254`
- Alzheimer disease: `MONDO_0004975`

## Understanding Results

The analysis results include:

1. **Summary Statistics**: Basic statistics of PGS scores by ancestry group and PGS study
2. **Correlations**: Correlation matrices showing how different PGS studies relate to each other
3. **Individual Variance**: Measurement of how consistently different PGS studies rank individuals within each ancestry group
4. **PGS Variance**: Measurement of how each PGS deviates from the consensus prediction across individuals

Visualizations include:

1. **Distribution plots** by ancestry group for each PGS
2. **Standardized score distributions** (z-scores)
3. **Correlation heatmaps** showing relationships between different PGS models
4. **Individual Variance plots** showing the stability of PGS predictions across ancestry groups
5. **PGS Variance plots** showing how each PGS deviates from the consensus prediction
6. **Average PGS Variance plots** showing the average variance by ancestry group

### Individual Variance Metric

The individual variance metric quantifies the stability of PGS predictions across different studies:

- For each individual, we calculate the variance of their z-scores across all PGS studies
- These individual variances are then averaged within each ancestry group
- Lower variance indicates more stable predictions (i.e., different PGS models consistently rank individuals similarly)
- Higher variance suggests less consistency across different PGS models

This metric is particularly useful for comparing prediction stability between European and non-European ancestry groups, as PGS studies typically show higher variance in non-European populations due to training bias.

### PGS Variance Metric

The PGS variance metric quantifies how each PGS deviates from the consensus prediction:

- For each individual, we calculate the average z-score across all PGS studies, which serves as the "true" z-score
- For each PGS, we calculate the variance of its predictions from these "true" z-scores
- PGS with higher variance are those that tend to deviate more from the consensus prediction
- This metric helps identify which PGS are outliers compared to the overall consensus

By comparing PGS variance across ancestry groups, researchers can identify which models are more consistent or divergent for different populations, potentially highlighting PGS that might be more biased for certain ancestry groups.

## Citing PGS-Compare

If you use PGS-Compare in your research, please cite:

```
PGS-Compare: A tool for analyzing the stability of polygenic scores across ancestry groups
```

And also cite the underlying tools:

- The PGS Catalog: https://doi.org/10.1038/s41588-021-00783-5
- The 1000 Genomes Project: https://doi.org/10.1038/nature15393
- PLINK 2: https://www.cog-genomics.org/plink/2.0/

## License

This project is licensed under the MIT License - see the LICENSE file for details.
