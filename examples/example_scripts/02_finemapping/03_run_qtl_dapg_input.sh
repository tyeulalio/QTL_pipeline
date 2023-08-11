#!/usr/bin/env bash


OUTPUT_DIR="../../example_output/02_finemapping/03_qtl_dapg_input"
mkdir -p $OUTPUT_DIR


# TORUS PRIORS DIRECTORY
# priors computed in the previous script 02
PRIORS_DIR="../../example_output/02_finemapping/02_run_torus/torus"

# controlled covariates file
# these are the covariates used in QTL mapping
#COVARIATES_FILE="/home/eulalio/deconvolution/new_rosmap/output/05_qtl_analysis/11_formatted_qtl_input/${PROPORTION_TYPE}"
COVARIATES_FILE="../../example_data/tensorqtl/GEUVADIS.445_samples.covariates.txt"

# molecular phenotype file used in QTL mapping
PHENO_FILE="../../example_data/tensorqtl/GEUVADIS.445_samples.expression.bed.gz"

# genotype counts file, this file can be created using plink
GENO_FILE="../../example_data/tensorqtl/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18.traw"

#GENOTYPE_DIR="${HOME}/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/sbam_genotype_files"
#GENOTYPE_DIR="../../example_data/tensorqtl/"


python3 ../../../02_finemapping/finemap_qtls/03_create_qtl_dapg_input_parallel.py \
    -output_dir $OUTPUT_DIR \
    -pheno_file $PHENO_FILE \
    -geno_file $GENO_FILE \
    -priors_dir $PRIORS_DIR \
    -covariates_file $COVARIATES_FILE
    #$SLURM_ARRAY_TASK_ID \
    #$OUTPUT_DIR \
    #$NUM_JOBS \
    #$PRIORS_DIR \
    #$PHENO_DIR \
    #$GENOTYPE_DIR
