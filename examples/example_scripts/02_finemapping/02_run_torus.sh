#!/usr/bin/env bash

# run torus to get SNP PIPs


# create the main output directory
OUTPUT_DIR="../../example_output/finemapping/02_run_torus"
mkdir -p $OUTPUT_DIR

# output from previous step 01_create_torus_input.py
INPUT_FILE="../../example_output/finemapping/01_create_torus_input/GEUVADIS_fastqtl_single_snp_output.tsv.gz"


# call the run_torus script from the qtl pipeline
../../../02_finemapping/finemap_qtls/02_run_torus.sh \
    $INPUT_FILE \
    $OUTPUT_DIR
