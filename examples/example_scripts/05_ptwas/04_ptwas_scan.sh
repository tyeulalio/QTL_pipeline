#!/usr/bin/env bash

# output directory
OUTDIR="../../example_output/05_ptwas/04_ptwas_scan"
mkdir -p $OUTDIR

# GWAS FILE - should be formatted appropriately
GWAS_FILE="../../example_data/gwas/wightman_gwas_file.vcf.gz"


for CHROM in {1..22}
do
    PTWAS_WEIGHT_FILE="../../example_output/05_ptwas/03_gambit_db/chrom${CHROM}_all_gene.ptwas_weights.gambit.vcf.gz"
    echo $PTWAS_WEIGHT_FILE

    NUM_LINES=$(zcat ${PTWAS_WEIGHT_FILE} | wc -l)
    echo "num lines ${NUM_LINES}"

    if [ $NUM_LINES -lt 3 ]; then
        echo "Less than 3 lines in input file. Skipping chromosome ${CHROM}"
        continue
    fi

    # ld panel from the GAMBIT github, updated for hg38
    LD_PANEL_FILES="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/06_swapped_ld_panels/$GWAS_TYPE/chr*.vcf.gz"
    set -f

    OUTPUT_PREFIX="${RUN_OUTDIR}/chr${CHROM}_ptwas_scan"

    ## --- run script here
    ../../../04_ptwas/04_ptwas_scan.sh \
        $GWAS_FILE \
        $PTWAS_WEIGHT_FILE \
        $LD_PANEL_FILES \
        $OUTPUT_PREFIX
done
