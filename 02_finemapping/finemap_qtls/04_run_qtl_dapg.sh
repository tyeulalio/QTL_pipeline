#!/usr/bin/env bash

# run the deterministic approximation of posteriors (DAP)
# to get pobabilistic QTL annotations for fastENLOC
# following GTEx instructions here: https://github.com/xqwen/dap/tree/master/gtex_v8_analysis

# -d sbam file
# -p prior file from torus
# --all forces output of all SNPs, not just noteworthy ones
# -t number of parallel threads
# -ld_control lowest LD threshold (R^2) to admit a SNP into a signal cluster


SBAM_FILE=$1
PRIOR_FILE=$2
LD_CONTROL=$3
GENE=$4
OUTFILE=$5
THREADS=$6

# run dap-g
time dap-g -d $SBAM_FILE \
    -p $PRIOR_FILE \
    -ld_control $LD_CONTROL \
    --all \
    -t $THREADS \
    -n ${GENE} \
    > $OUTFILE
