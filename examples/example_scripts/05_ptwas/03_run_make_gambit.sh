#!/usr/bin/env bash

# load conda environment
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate r

# output directory
OUTDIR="../../example_output/05_ptwas/03_gambit_db"
mkdir -p $OUTDIR

# concatenated ptwas weights file from previous step
PTWAS_WEIGHTS_FILE="../../example_output/05_ptwas/02_concatenated_weights/all_gene.ptwas_weights.txt.gz"

# output file
OUTPUT_FILE="${OUTDIR}/all_gene.ptwas_weights.gambit.vcf"
rm -f $OUTPUT_FILE
touch $OUTPUT_FILE

# get command line input
# run the r script
Rscript ../../../04_ptwas/make_GAMBIT.DB.R \
    -d ${PTWAS_WEIGHTS_FILE} \
    -o ${OUTPUT_FILE}


# format the output file
mamba activate colocalization

../../../04_ptwas/03_run_make_gambit.sh \
    ${PTWAS_WEIGHTS_FILE} \
    $OUTDIR \
    $OUTPUT_FILE
