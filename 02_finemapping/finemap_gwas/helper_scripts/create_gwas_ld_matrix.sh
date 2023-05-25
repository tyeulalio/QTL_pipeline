#!/usr/bin/env bash

# using plink2 on scg
#module load plink/1.90b6.25.lua
#module load plink2

# load the anaconda environment
. ${HOME}/micromamba/etc/profile.d/conda.sh
. ${HOME}/micromamba/etc/profile.d/mamba.sh
mamba activate base
mamba activate colocalization


# this script runs plink to generate LD matrices
# call this script for each gene and the SNPs within the QTL window

# get command line input
# REGIONID = region id
# SNPFILE = text file with SNPs that fall within the tested QTL window for this gene (eg 1Mb of TSS)
# OUTDIR = string, output directory to save files
# BFILE = plink formatted genotype data
# THREADS = number of threads to use
REGIONID=$1
SNPFILE=$2
OUTDIR=$3
BFILE=$4
THREADS=$5

# create a directory to save the log files
LOGDIR="${OUTDIR}/logs"
mkdir -p ${LOGDIR}

echo "region $REGIONID"
echo "snpfile $SNPFILE"

# make LD matrix
plink \
    --bfile $BFILE \
    --r square\
    --extract "${SNPFILE}" \
    --threads $THREADS \
    --write-snplist \
    --out "${OUTDIR}/rosmap_ld_${REGIONID}" \
    > >(tee -a "${LOGDIR}/${REGIONID}_stdout.log") \
    2> >(tee -a "${LOGDIR}/${REGIONID}_stderr.log" >&2)
