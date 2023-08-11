#!/usr/bin/env bash

PTWAS_WEIGHTS_FILE=$1
OUTDIR=$2
OUTPUT_FILE=$3

# format the output file by chromosome
# to save time in scan step
for CHROM in {1..22};
do
    echo "Processing chromosome ${CHROM}"

    # get the header from gabmit file
    FORMATTED_OUTFILE="${OUTDIR}/chrom${CHROM}_all_gene.ptwas_weights.gambit.vcf"
    touch $FORMATTED_OUTFILE

    echo "Writing output to ${FORMATTED_OUTFILE}"

    echo "printing header"
    cat $OUTPUT_FILE | head | grep "^#" > $FORMATTED_OUTFILE

    # remove chr from chromosome and from ID
    echo "printing remaining formatted lines"
    cat $OUTPUT_FILE | grep -v "^#" | grep "^chr${CHROM}\s" | sed -e "s/chr//g" >> $FORMATTED_OUTFILE

    # bgzip the output file
    echo "zipping output file"
    bgzip -f $FORMATTED_OUTFILE

    # tabix the output
    echo "running tabix"
    tabix -f -p 'vcf' "${FORMATTED_OUTFILE}.gz"
done
