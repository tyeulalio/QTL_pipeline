# create the input file for torus using the output from fastqtl
import os
import pandas as pd
import gzip
import sys
from datetime import datetime
import argparse


def get_args():
    # get command line arguments
    parser = argparse.ArgumentParser()
    # required parameters
    parser.add_argument('-qtl', '--qtl_file', type=str, required=True, help="QTL parquet files")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output directory")

    # optional parameters
    parser.add_argument('-p', '--prefix', type=str, required=False, default="torus", help="Prefix for output file name")
    # parser.add_argument('-chroms', '--chromosomes', type=list, required=False, default=range(1,23), help="List of chromosomes to process")
    parser.add_argument('-run_example', type=bool, required=False, help="Boolean, indicating whether to run the example data")

    args = parser.parse_args()

    return(args)


def process_chromosome(chrom, savefile, qtl_file):
    # use results of all snp-phenotype pairs
    datafile = qtl_file.replace("*", str(chrom))

    print("Reading file", datafile)

    if not os.path.isfile(datafile):
        print("not a file")

    # read in the data file
    qtls = pd.read_parquet(datafile)

    # subset to the colums that we want to keep
    qtls = qtls.loc[:,['phenotype_id', 'variant_id', 'tss_distance', 'pval_nominal', 'slope', 'slope_se']]

    return(qtls)


    


def main():
    start = datetime.now()

    # get command line arguments
    args = get_args()

    savedir = args.output
    if not os.path.exists(savedir):
        print("creating output directory", savedir)
        os.makedirs(savedir)

    # create the savefile for output
    savefile = "{}/{}_fastqtl_single_snp_output.tsv".format(savedir, args.prefix)
    print("output file:", savefile)

    # clear the file after testing
    with open(savefile, "w+"):
        next

    # process each chromosome
    chromosomes = range(1,23)
    if args.run_example:
        print("Running example")
        chromosomes = [18]

    full_qtls = pd.DataFrame()
    for chrom in chromosomes:
        chrom_qtls = process_chromosome(chrom, savefile, args.qtl_file)

        # concatenate all of the chromosome files together
        if (full_qtls.empty):
            full_qtls = chrom_qtls
        else:
            full_qtls = pd.concat([full_qtls, chrom_qtls])
        print(chrom, chrom_qtls.shape, full_qtls.shape)

        print(full_qtls.shape)

    # save data to output file
    print("Writing output to file")
    full_qtls.to_csv(savefile, mode='a', index=False, header=False, sep=' ')

    # gzip the file
    print("Zipping output file")
    zipfile = "{}.gz".format(savefile)
    with open(savefile, 'rb') as f_in, gzip.open(zipfile, 'wb') as f_out:
        f_out.writelines(f_in)


    print("--> Finished mapping", datetime.now())
    end = datetime.now()
    duration = end - start
    print(duration)


if __name__ == '__main__':
    main()
