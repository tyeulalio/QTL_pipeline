#!/usr/bin/env bash

# use fastqtl to create annoations from the fine-mapped QTL data
# following the tutorial on the github https://github.com/xqwen/fastenloc/tree/master/tutorial/
# formats the QTL DAP-G output for FastQTL
# the summarize_dap2enlco.pl script come from the FastEnloc github
# this script requires a SNP annotation file in VCF form - created by script 05/02 using plink

# arguments
# dir = directory with DAP fine-mapped qtl results
# vcf = file to annotate all SNP's positions
# outdir = output directory

DAP_DIR=$1
VCF=$2
OUTDIR=$3

time ${HOME}/programs/cloned_repo/fastenloc/src/summarize_dap2enloc.pl \
    -dir $DAP_DIR \
    -vcf $VCF \
    | gzip - > "${OUTDIR}/fastenloc.qtl.annotation.vcf.gz"
