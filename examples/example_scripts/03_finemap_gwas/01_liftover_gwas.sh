#!/usr/bin/env bash

# liftover the VCFs from hg37 to hg38
# using GATK liftoverVCF

# make sure to load the GATK module
module load gatk/4.0.10.0

# GWAS data downloaded from the wightman et al 2021 Alzheimer's GWAS
DATAFILE="../../example_data/gwas/wightman_etal_2021_results.vcf"

# output the lifted file
OUTPUT_FILE="../../example_data/gwas/wightman_etal_2021_results_hg38.vcf"

# file to store rejected entries
REJECT_FILE="../../example_output/03_gwas/01_liftover/rejected_records.vcf"

# CHAIN file can be downloaded from UCSC
# REFERENCE_SEQUENCE is a genome reference file for the appropriate build

gatk LiftoverVcf \
    --INPUT="${DATAFILE}" \
    --OUTPUT="${OUTPUT_FILE}" \
    --VERBOSITY=DEBUG \
    --CHAIN="/home/eulalio/shared/liftOver/chains/hg19ToHg38.over.chain.gz" \
    --REFERENCE_SEQUENCE="/reference/RefGenomes/GATK_Resource_Bundle/hg38/hg38.fa" \
    --REJECT="${REJECT_FILE}" \
    --MAX_RECORDS_IN_RAM 10000 \
    --RECOVER_SWAPPED_REF_ALT
