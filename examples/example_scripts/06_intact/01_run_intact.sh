#!/usr/bin/env bash

# load conda environment
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate r


SAVEDIR="../../example_output/06_intact/"
PTWAS_FILE="../../example_output/05_ptwas/04_ptwas_scan/all_chroms_ptwas_scan.stratified_out.txt"
FASTENLOC_FILE="../../example_output/04_colocalization/01_run_fastenloc/fastqtl_.enloc.gene.out"

Rscript ../../../05_intact/01_run_intact.R \
    $SAVEDIR \
    $PTWAS_FILE \
    $FASTENLOC_FILE
