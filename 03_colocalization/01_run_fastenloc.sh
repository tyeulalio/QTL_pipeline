#!/usr/bin/env bash

# run fastenloc

OUTDIR=$1
QTL_FILE=$2
GWAS_FILE=$3
THREADS=$4
NUM_GWAS_VARS=$5


# create the output directory if it doesn't exist
mkdir -p $OUTDIR

# run fastenloc 
fastenloc \
    -total_variants $NUM_GWAS_VARS \
    -thread $THREADS \
    -prefix "${OUTDIR}/fastqtl_" \
    -eqtl $QTL_FILE \
    -gwas $GWAS_FILE
