# TensorQTL Mapping Script

This script is designed to map cis-QTLs using TensorQTL. It primarily follows the example provided in the [TensorQTL Examples](https://github.com/broadinstitute/tensorqtl/blob/master/example/tensorqtl_examples.ipynb).

## Features

- Load phenotypes and covariates
- Read genotypes using PLINK
- Perform cis-QTL mapping with nominal p-value calculation
- Map cis-QTLs with empirical p-values for genome-wide FDR
- Format output and extract significant QTL pairs

## Prerequisites

Before running this script, ensure you have the following installed:

- Python 3
- Pandas
- TensorQTL
- Torch

## Usage

This script accepts command-line arguments for processing. The main parameters include:

- `-output`: Directory path to save output files (required)
- `-prefix`: Prefix for output file names (default: "cis_qtls")
- `-pheno`: Phenotype input file (required)
- `-covariates`: Covariates input file (required)
- `-plink`: Prefix for PLINK .bim/.fam/.bed genotype files (required)
- `-window`: Window size for cis-QTL mapping (default: 1e6)
- `-run_example`: Flag to run the example data (default: False)

Run the script from the command line as follows:

```bash
python 01_run_tensorqtl.py -output [output directory] -pheno [phenotype file] -covariates [covariates file] -plink [PLINK prefix]
```

Replace `[output directory]`, `[phenotype file]`, `[covariates file]`, and `[PLINK prefix]` with your specific file paths.

## Example

To run the script with example parameters:

```bash
python 01_run_tensorqtl.py -output /path/to/output -pheno /path/to/pheno.bed -covariates /path/to/covariates.txt -plink /path/to/plink_prefix
```

Run “Python3 01_qtl_mapping.py -h” for information on the required input.

See example script `01_qtl_mapping/01_map_qtls.sh` for an example on how to run this script.

## Output

The script outputs several files in the specified output directory, including:

- parquet files containing QTL mapping results for each chromosome
- `cis_qtl_signif_pairs.csv`
    - FDR corrected significant QTLs
- `cis_qtls.time_to_run.txt`
    - runtime of script
- `cis_qtl.cis_qtl.txt.gz`
    - summary statistics output from TensorQTL
- `cis_qtl_summary_stats.csv`
    - summary statistics output from TensorQTL saved as CSV



See examples in `example_output/tensorqtl`

