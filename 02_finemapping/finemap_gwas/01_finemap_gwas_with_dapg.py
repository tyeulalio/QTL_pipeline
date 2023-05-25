# format the GWAS LD matrices for DAPG input
# create an LD matrix for each gene
# calls dap-g for each gene also

import os
import pandas as pd
import gzip
import sys
import numpy as np
import subprocess
import argparse
from datetime import datetime


def get_args():
    # get command line arguments
    parser = argparse.ArgumentParser()
    # required parameters
    parser.add_argument('-m', '--map', type=str, required=True, help="Map of SNPs to regions (e.g. LD blocks)")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output directory")
    parser.add_argument('-z', '--zscore', type=str, required=True, help="Z-score file")
    parser.add_argument('-p', '--plink', type=str, required=True, help="PLINK formatted genotype data")
    # optional parameters
    parser.add_argument('-t', '--tmpdir', type=str, required=False, help="Temporary directory to store intermediate temp files")
    parser.add_argument('-chk', '--check', type=bool, required=False, help="Check setting for debugging")
    parser.add_argument('-r', '--region', type=int, required=False, help="Region number to fine-map")
    parser.add_argument('-th', '--threads', type=int, required=False, help="Number of threads to use for parallel steps")
    parser.add_argument('-s', '--savetmp', type=bool, required=False, help="Save temporary files", default=False)

    args = parser.parse_args()

    return(args)




def get_snp_map(snp_map_file):
    # get the overlap of ld blocks with the snps on this chromosome

    print(snp_map_file)
    ld_snps = pd.read_csv(snp_map_file, sep='\t')
    # print(ld_snps.head())

    # select the columns we need
    qtl_res = ld_snps.loc[:, ['phenotype_id', 'variant_id']]
    print(qtl_res.head())

    return(qtl_res)



def create_region_ld_matrix(region_id, ordered_snps, tmpdir, plink_file, threads):
    # create LD matrix for each gene separately

    # skip if we already have the ld matrix saved
    ldfile = "{}/{}_ld_matrix.txt".format(tmpdir, region_id) 
    # if (os.path.isfile(ldfile)):
        # return()


    # sort the snp list by position
    ordered_snps['pos'] = ordered_snps['variant_id'].str.split('_', expand=True)[1]
    ordered_snps = ordered_snps.sort_values('pos')
    # print(ordered_snps.head(10))

    # create a snp file for this gene
    snpfile = "{}/{}_snp_list.txt".format(tmpdir, region_id) 
    snps_overlap = ordered_snps.variant_id.values
    pd.DataFrame(snps_overlap).to_csv(snpfile, header=False, index=False)

    # correlation file for this gene
    corr_file = "{}/rosmap_ld_{}.ld".format(tmpdir, region_id)

    # create the LD matrix if it doesn't exit
    # if not os.path.isfile(corr_file):
    # get the correlations using plink script
    ld_helper_script = "/home/eulalio/qtl_pipeline/02_finemapping/finemap_gwas/helper_scripts/create_gwas_ld_matrix.sh"
    cmd = [ld_helper_script, region_id, snpfile, tmpdir, plink_file, str(threads)]
    p = subprocess.Popen(cmd)
    p.wait()

    return()

    
def create_zscore_file(region_id, snp_map, tmpdir, zscore_file):
    # create the zscore for this region

    # get the snps for this region
    region_snps = snp_map.loc[snp_map.phenotype_id == region_id,].variant_id.drop_duplicates()
    print("region snps length", len(region_snps))


    # load the zscores
    zscores = pd.read_csv(zscore_file, sep=r'\s+', names=['snp_id', 'zscore'])

    # create a list
    region_snps_list = list(region_snps)

    # overlap the snps in gwas and qtls
    snps_overlap = list(set(region_snps_list).intersection(set(zscores.snp_id.tolist())))
    # put snps overlap in a dataframe already
    ordered_snps = pd.DataFrame(snps_overlap, columns=['variant_id'])
    print(ordered_snps.head())

    print('snp overlap length', len(snps_overlap))

    # get the gwas zscores for these SNPs
    zscores = zscores.set_index('snp_id', drop=False)

    region_zscores = zscores.loc[ordered_snps.variant_id.values,]
    print(region_zscores.head())

    zscores_savefile = "{}/{}_zscores.txt".format(tmpdir, region_id)
    print("Saving tmp zscores to", zscores_savefile)
    region_zscores.to_csv(zscores_savefile, sep=' ', header=False, index=False)

    return(ordered_snps)
    

def remove_tmp_files(region_id, tmpdir):
    cmd = "rm -f {}/*{}*".format(tmpdir, region_id)
    subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()


def run_dapg(region_id, tmpdir, output_dir, threads):
    # run dapg for this gene
    zscore_file = "{}/{}_ordered_zscores.txt".format(tmpdir, region_id)
    ld_file = "{}/rosmap_ld_reformatted_{}.ld".format(tmpdir, region_id) 

    dapg_script = "/home/eulalio/qtl_pipeline/02_finemapping/finemap_gwas/helper_scripts/run_gwas_dapg.sh"
    dapg_outdir = "{}/gwas_dapg/".format(output_dir)
    cmd = [dapg_script, region_id, zscore_file, ld_file, dapg_outdir, str(threads)]
    subprocess.Popen(cmd).communicate()


def check_output(final_dir, genes):
    missing_genes = []

    # check if all of the genes were processed
    for i in range(len(genes)):
        gene_id = genes[i]
        final_file = "{}/{}.dapg".format(final_dir, gene_id)

        if not os.path.isfile(final_file):
            missing_genes.append(gene_id)

    print("Found {} missing genes out of {} total genes".format(len(missing_genes), len(genes)))

    return(missing_genes)

def order_zscores(region_id, tmpdir):
    # reorder the zscores file to match the ld matrix

    # read in ld matrix snplist
    snpfile = "{}/rosmap_ld_{}.snplist".format(tmpdir, region_id)
    snplist = pd.read_csv(snpfile, names=['snp'])
    # print(snplist.head())
    
    # read in zscores
    zfile = "{}/{}_zscores.txt".format(tmpdir, region_id)
    zscores = pd.read_csv(zfile, names=['snp', 'zscore'], delimiter=' ')
    # print(zscores.head())

    # select zscore in order of snplist
    ordered_zscores = zscores.set_index('snp')
    ordered_zscores = ordered_zscores.reindex(index = snplist['snp'])
    ordered_zscores = ordered_zscores.reset_index()
    # print(ordered_zscores.head())
    # print(snplist.shape)
    # print(ordered_zscores.shape)

    # write to output file
    zscores_savefile = "{}/{}_ordered_zscores.txt".format(tmpdir, region_id)
    ordered_zscores.to_csv(zscores_savefile, sep=' ', header=False, index=False)
    

def format_ld_matrix(region_id, tmpdir):
    # format the ld matrix
    # change tabs to spaces

    ldfile = "{}/rosmap_ld_{}.ld".format(tmpdir, region_id)
    outfile = "{}/rosmap_ld_reformatted_{}.ld".format(tmpdir, region_id)
    out = open(outfile, 'w+')
    
    i = 0
    with open(ldfile, "r") as f:
        for line in f:
            line = line.split()

            outline = " ".join(line)
            out.write("{}\n".format(outline))

    out.close()

def process_region(region_num, region_id, regions, snp_map, tmpdir, start, zscore_file, plink_file, threads, savetmp, output_dir):
    print("Starting on region {} {} out of {} {}".format(region_num, region_id, len(regions), datetime.now()))

    # create gene z score file
    print("Creating Zscore file {}".format(datetime.now()))
    ordered_snps = create_zscore_file(region_id, snp_map, tmpdir, zscore_file)

    print("Creating LD matrix {}".format(datetime.now()))
    create_region_ld_matrix(region_id, ordered_snps, tmpdir, plink_file, threads)

    # rewrite the zscores to match ld matrix cpg order
    print("Ordering Z-scores")
    order_zscores(region_id, tmpdir)

    # change format of ld matrix
    print("Formatting LD matrix")
    format_ld_matrix(region_id, tmpdir)

    # run dap-g
    print("Run DAP-G {}".format(datetime.now()))
    run_dapg(region_id, tmpdir, output_dir, threads)

    # remove all temp files
    # don't need to do this when submitting sbatch script
    if (not savetmp):
        print("Removing temp files {}".format(datetime.now()))
        remove_tmp_files(region_id, tmpdir)

    print("Finished gene {} out of {} {}".format(region_num, len(regions), datetime.now()))

    print("Finished {}".format(datetime.now()))
    endtime = datetime.now()

    duration = endtime - start
    print("Duration {}".format(duration))



def main():
    # get the command line arguments
    args = get_args()

    # set all of the arguments to variables
    snp_map_file = args.map
    check = args.check
    region_num = args.region
    output_dir = args.output
    zscore_file = args.zscore
    tmpdir = args.tmpdir
    plink_file = args.plink
    threads = args.threads
    savetmp = args.savetmp

    if tmpdir == None:
        tmpdir = output_dir


    start = datetime.now()

    # get genes mapped to SNPs tested in QTL analysis
    print("Loading SNP map {}".format(datetime.now()))
    snp_map = get_snp_map(snp_map_file)

    # create LD matrix for each gene separately
    regions = snp_map.phenotype_id.unique()
    print("Regions to fine-map:", len(regions))

    # if no gene id exists, do not continue with LD and DAPG
    if region_num is None:
        print("Include a region number between 1 and {} to fine-map".format(len(regions)))
        return()

    print("Processing region number", region_num)

    # check that the number is within bounds
    region_num <- region_num - 1  # adjust for zero-index
    if (region_num < 0) | (region_num > len(regions)):
        print("Region number is out of bounds")
        print("Include a region number between 1 and {} to fine-map".format(len(regions)))
        return()

    # assign the gene id 
    region_id = regions[region_num]
    print(region_id)



    # set up output file
    final_file = "{}/{}.dapg".format(output_dir, region_id)

    # process this region
    process_region(region_num, region_id, regions, snp_map, tmpdir, start, zscore_file, plink_file, threads, savetmp, output_dir)

if __name__ == '__main__':
    main()
