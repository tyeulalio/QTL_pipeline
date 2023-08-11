#!/usr/bin/env bash

# generate genotype VCFs using plink

plink2 --bfile "../example_data/tensorqtl/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18" \
    --recode vcf \
    --out "../example_data/tensorqtl/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18"
