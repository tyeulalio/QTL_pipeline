#!/usr/bin/env bash

# original data was downloaded from the TensorQTL Github's example data
# this script converts Plink2 format to Plink1 format to work with the current scripts in this pipeline
plink2 --pfile ../example_data/tensorqtl/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18 \
    --make-bed \
    --output-chr chrM \
    --out ../example_data/tensorqtl/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18
