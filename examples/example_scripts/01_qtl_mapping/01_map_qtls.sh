#!/usr/bin/env bash

# map qtls using qtl pipeline script

DATADIR="../../example_data"

PHENO_FILE="${DATADIR}/tensorqtl/GEUVADIS.445_samples.expression.bed.gz"
COVS_FILE="${DATADIR}/example_data/tensorqtl/GEUVADIS.445_samples.covariates.txt"
PLINK_PATH="${DATADIR}/tensorqtl/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18"

OUTPUT_DIR="../../example_output/tensorqtl/"
mkdir -p $OUTPUT_DIR

# make sure conda environment is activated
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate tensorqtl

# select the cis-window
CIS_WINDOW=1000000

python3 ../../01_qtl_mapping/01_run_tensorqtl.py \
    -output $OUTPUT_DIR \
    -pheno $PHENO_FILE \
    -covariates $COVS_FILE \
    -plink $PLINK_PATH \
    -window $CIS_WINDOW \
    -run_example True
