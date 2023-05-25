#!/bin/bash

# run the deterministic approximation of posteriors (DAP)
# to get pobabilistic GWAS for fastENLOC
# following GTEx instructions here: https://github.com/xqwen/dap/tree/master/gtex_v8_analysis

# -d_z z-score file
# -d_ld ld file
# --all forces output of all SNPs, not just noteworthy ones
# -t number of parallel threads
# -ld_control lowest LD threshold (R^2) to admit a SNP into a signal cluster

# command line arguments
# REGIONID = region id
# ZSCORE_FILE = string, file containing zscores for SNPs tested with this gene (window defined in QTL analysis, eg 1Mb from TSS)
# LD_FILE = LD matrix for SNPs tested within the window of this gene
# GWAS = gwas file
REGIONID=$1
ZSCORE_FILE=$2
LD_FILE=$3
OUTDIR=$4
THREADS=$5

# create the output directory
mkdir -p $OUTDIR

ld_cutoff=0.25

# name of the output file
OUTFILE="${OUTDIR}${REGIONID}.dapg" 

echo $ZSCORE_FILE
echo $LD_FILE

# run dap-g
dap-g -d_z ${ZSCORE_FILE} \
    -d_ld ${LD_FILE} \
    -ld_control ${ld_cutoff} \
    --all \
    -t $THREADS \
    -n ${REGIONID} \
    -d_n 762917 \
    > $OUTFILE 

    #-d_n 762917 \
