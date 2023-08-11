#!/usr/bin/env bash

# run fastenloc

NUM_GWAS_VARS=12688339

OUTDIR="../../example_output/04_colocalization/01_run_fastenloc"
mkdir -p $OUTDIR

# -total_variants = number of GWAS variants tested

# paths to fine-mapped qtls and gwas
QTL_FILE="../../example_output/02_finemapping/05_fastqtl_qtl_annotations/fastenloc.qtl.annotation.vcf.gz"
GWAS_FILE="../../example_output/03_gwas/06_formatted_finemapped_gwas/formatted_fastqtl_gwas_input_pips.vcf.gz"

../../../03_colocalization/01_run_fastenloc.sh \
    $OUTDIR \
    $QTL_FILE \
    $GWAS_FILE \
    4 \
    $NUM_GWAS_VARS
