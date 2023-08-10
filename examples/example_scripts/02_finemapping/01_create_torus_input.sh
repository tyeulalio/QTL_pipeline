#!/usr/bin/env bash

# create torus input from TensorQTL output

# output directory
OUTPUT_DIR="../../example_output/finemapping/01_create_torus_input"
mkdir -p $OUTPUT_DIR

# prefix for the output file
PREFIX="GEUVADIS"

# input qtl files
QTL_PREFIX="${CELL_TYPE}_${REGION_TYPE}_${SUMMARY_TYPE}"
set -f #disable star expansion
QTL_FILES="../../example_output/tensorqtl/cis_qtls.cis_qtl_pairs.chr*.parquet"

python3 ../../../02_finemapping/finemap_qtls/01_create_torus_input.py \
    -o $OUTPUT_DIR \
    -p $PREFIX \
    -qtl $QTL_FILES  \
    -run_example True
