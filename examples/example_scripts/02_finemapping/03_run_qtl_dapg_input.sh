#!/usr/bin/env bash


OUTPUT_DIR="../../example_output/finemapping/03_qtl_dapg_input"
mkdir -p $OUTPUT_DIR


# TORUS PRIORS DIRECTORY
PRIORS_DIR="../../example_output/finemapping/02_run_torus/torus"

# phenotype directory
#PHENO_DIR="/home/eulalio/deconvolution/new_rosmap/output/05_qtl_analysis/11_formatted_qtl_input/${PROPORTION_TYPE}"
PHENO_FILE="../../example_data/tensorqtl/GEUVADIS.445_samples.expression.bed.gz"

# genotype counts file, this file can be created using plink
GENO_FILE="../../example_data/tensorqtl/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18.traw"

#GENOTYPE_DIR="${HOME}/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/sbam_genotype_files"
#GENOTYPE_DIR="../../example_data/tensorqtl/"


python3 ../../../02_finemapping/finemap_qtls/03_create_qtl_dapg_input_parallel.py \
    -output_dir $OUTPUT_DIR \
    -pheno_file $PHENO_FILE \
    -geno_file $GENO_FILE
    #$SLURM_ARRAY_TASK_ID \
    #$OUTPUT_DIR \
    #$NUM_JOBS \
    #$PRIORS_DIR \
    #$PHENO_DIR \
    #$GENOTYPE_DIR
