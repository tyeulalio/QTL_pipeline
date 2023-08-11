#!/usr/bin/env bash


#OUTDIR="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/02_concatenated_weights/${PROPORTION_TYPE}"
OUTDIR="../../example_output/05_ptwas/02_concatenated_weights/"
mkdir -p $OUTDIR

# directory containing the ptwas weights
#WEIGHTSDIR="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/01_ptwas_weights/${PROPORTION_TYPE}/${RUNNAME}/ptwas_weights"
WEIGHTSDIR="../../example_output/05_ptwas/01_ptwas_weights/ptwas_weights"

## --- run script here
time ../../../04_ptwas/02_concat_weights.sh \
    $WEIGHTSDIR \
    $OUTDIR


