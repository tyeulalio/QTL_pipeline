#!/usr/bin/env bash

# create the genotype counts matrix for DAP-G using plink2
# make sure that you are counting the correct ref/alt for this!

plink2 \
    --pfile "../example_data/tensorqtl/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18" \
    --out "../example_data/tensorqtl/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18" \
    --export A-transpose
