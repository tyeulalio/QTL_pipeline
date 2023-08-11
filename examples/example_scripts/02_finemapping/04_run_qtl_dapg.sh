#!/usr/bin/env bash

# run the deterministic approximation of posteriors (DAP)
# to get pobabilistic QTL annotations for fastENLOC
# following GTEx instructions here: https://github.com/xqwen/dap/tree/master/gtex_v8_analysis

# -d sbam file
# -p prior file from torus
# --all forces output of all SNPs, not just noteworthy ones
# -t number of parallel threads
# -ld_control lowest LD threshold (R^2) to admit a SNP into a signal cluster

INPUT_DIR="../../example_output/02_finemapping/03_qtl_dapg_input"

# get file with genes to process
GENESFILE="../../example_output/02_finemapping/03_qtl_dapg_input/genes_list.txt"
NUM_GENES=$(wc -l < $GENESFILE)

echo "Processing ${NUM_GENES} genes"

# -- input and output directory
OUTDIR="../../example_output/02_finemapping/04_dapg"
mkdir -p $OUTDIR

for ((i=1; i<=$NUM_GENES; i++))
do
    # get the gene from the file
    GENE=$(sed "${i}q;d" $GENESFILE)

    echo "$GENE $i"

    # format gene properly
    FORMATTED_GENE=$(echo "$GENE" | tr \\. - )

    OUTFILE="${OUTDIR}/${FORMATTED_GENE}.dapg" 

    # skip ifoutput file exists
    #if [ -f "$OUTFILE" ]; then
        #echo "output file exists"
        #exit
    #fi

    SBAM_FILE="${INPUT_DIR}/sbams/${GENE}_fastqtl_singl_snp_output.sbam" 
    PRIORS_FILE="${INPUT_DIR}/priors/${GENE}_priors.txt" 


    # run dap-g
    dap-g -d $SBAM_FILE \
        -p $PRIORS_FILE \
        -ld_control 0.25 \
        --all \
        -t 4 \
        -n ${GENE} \
        > $OUTFILE

done
