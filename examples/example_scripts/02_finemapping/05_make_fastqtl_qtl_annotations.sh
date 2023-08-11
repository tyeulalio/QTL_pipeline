
# use fastqtl to create annoations from the fine-mapped QTL data
# following the tutorial on the github https://github.com/xqwen/fastenloc/tree/master/tutorial/
# formats the QTL DAP-G output for FastQTL
# the summarize_dap2enlco.pl script come from the FastEnloc github
# this script requires a SNP annotation file in VCF form - created by script 05/02 using plink

# arguments
# dir = directory with DAP fine-mapped qtl results
# vcf = file to annotate all SNP's positions


# create output directory
OUTDIR="../../example_output/02_finemapping/05_fastqtl_qtl_annotations"
mkdir -p ${OUTDIR}

# input DAP-G directory
DAPDIR="../../example_output/02_finemapping/04_dapg"

# need a genotype VCF (can use PLINK to generate this)
# see example in generate_data directory
GENO_VCF="../../example_data/tensorqtl/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.chr18.vcf.gz"

time /home/eulalio/programs/cloned_repo/fastenloc/src/summarize_dap2enloc.pl \
    -dir $DAPDIR \
    -vcf $GENO_VCF \
    | gzip - > "${OUTDIR}/fastenloc.qtl.annotation.vcf.gz"
