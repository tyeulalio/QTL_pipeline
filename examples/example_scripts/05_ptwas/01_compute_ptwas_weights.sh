#!/usr/bin/env bash

# compute weights for PTWAS


# get list of genes 
DAPDIR="../../example_output/02_finemapping/04_dapg"
GENES_LIST=$(ls $DAPDIR | sed -r 's/.dapg//g')
FILE_COUNT=$(echo $GENES_LIST | wc -w)
echo "Total files: $FILE_COUNT"

# SBAMS used for DAPG finemapping of QTLs
SBAMSDIR="../../example_output/02_finemapping/03_qtl_dapg_input/sbams"

# output directory
OUTDIR="../../example_output/05_ptwas/01_ptwas_weights"
mkdir -p $OUTDIR

# create array of files
readarray -t GENES_ARRAY <<<"$GENES_LIST"

time for i in $(seq 1 $FILE_COUNT); do
    echo "starting $i"
    GENE_NUM=$i

    # select the gene that we're working with
    GENE=${GENES_ARRAY[$GENE_NUM]}

    # run ptwas builder for each gene
    # file containing sbam info
    SBAMGENE=$(echo $GENE | sed 's/-/./1')
    SBAMFILE="${SBAMSDIR}/${SBAMGENE}_fastqtl_singl_snp_output.sbam"
    echo "SBAM FILE: ${SBAMFILE}"

    # file containing dap-g output
    DAPGENE=$(echo $GENE | tr \\. -)
    DAPFILE="${DAPDIR}/${DAPGENE}.dapg"
    echo "DAP FILE: ${DAPFILE}"


    echo "processing $GENE $GENE_NUM out of $FILE_COUNT"
    ((i+=1))

    # format output file
    # put genes in same format
    GENE=$(echo $GENE | tr - .)

    ../../../04_ptwas/01_compute_ptwas_weights.sh \
        $DAPFILE \
        $SBAMFILE \
        $GENE \
        $OUTDIR

done
