# PGS-Compare

PGS-Compare is a Python package for analyzing and comparing Polygenic Scores (PGS) across ancestry groups. It uses the PGS Catalog and 1000 Genomes data to help researchers evaluate the stability of PGS scores across different ancestry groups.

## Features

- Download necessary data from the 1000 Genomes project and reference panels
- Fetch and calculate PGS scores for specific traits using the PGS Catalog
- Analyze PGS score distributions across different ancestry groups
- Compare consistency of PGS scores from different studies
- Visualize results with publication-ready plots

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

# Analyze the results
pgs-compare analyze --trait-id MONDO_0005180

# Generate visualizations
pgs-compare visualize --trait-id MONDO_0005180

# Or run the full pipeline
pgs-compare pipeline MONDO_0005180
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
- dict: Information about the calculation

#### analyze

```python
pgs.analyze(trait_id=None, scores_file=None)
```

Analyze PGS scores across ancestry groups.

Parameters:
- `trait_id` (str, optional): Trait ID. Used for organizing output if provided.
- `scores_file` (str, optional): Path to the scores file (aggregated_scores.txt.gz).

Returns:
- dict: Analysis results

#### visualize

```python
pgs.visualize(trait_id=None, analysis_results=None)
```

Visualize PGS analysis results.

Parameters:
- `trait_id` (str, optional): Trait ID. Used for organizing output if provided.
- `analysis_results` (dict, optional): Analysis results from analyze().

Returns:
- dict: Dictionary with paths to the generated plots

#### run_pipeline

```python
pgs.run_pipeline(trait_id, include_child_pgs=True, max_variants=None,
                run_ancestry=False, visualize=True)
```

Run the full pipeline (calculate, analyze, visualize) for a specific trait.

Parameters:
- `trait_id` (str): Trait ID (e.g., "MONDO_0005180" for Parkinson's disease)
- `include_child_pgs` (bool): Whether to include child-associated PGS IDs
- `max_variants` (int, optional): Maximum number of variants to include in PGS
- `run_ancestry` (bool): Whether to run ancestry analysis
- `visualize` (bool): Whether to generate visualization plots

Returns:
- dict: Pipeline results

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
3. **Deviations**: Analysis of how each PGS study deviates from the average across ancestry groups

Visualizations include:

1. Distribution plots by ancestry group for each PGS
2. Standardized score distributions (z-scores)
3. Correlation heatmaps
4. Deviation plots showing the stability of PGS scores across ancestry groups

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