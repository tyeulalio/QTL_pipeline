#!/usr/bin/env bash

# run torus to estimate fine-mapping priors


# input file for torus
INPUT_FILE=$1

# output directory
# if directory doesn't exist, create it
OUTDIR=$2
mkdir -p $OUTDIR

time torus -d $INPUT_FILE \
    --fastqtl \
    -dump_prior "${OUTDIR}/torus" \
    > >(tee -a "${OUTDIR}/stdout.log") \
    2> >(tee -a "${OUTDIR}/stderr.log" >&2)
