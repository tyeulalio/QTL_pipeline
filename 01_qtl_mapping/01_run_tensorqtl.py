#!/usr/bin/env python3

## This script maps cis-QTLs using TensorQTL
# mainly follows the example found here: https://github.com/broadinstitute/tensorqtl/blob/master/example/tensorqtl_examples.ipynb

import pandas as pd
import torch
from datetime import datetime
import subprocess
import os
import pickle
import sys
import argparse
import tensorqtl

from tensorqtl import genotypeio, cis, trans, susie, post
from post import *


def get_args():
    # get command line arguments
    print("Loading command line arguments")

    # create argument parser
    parser = argparse.ArgumentParser(prog="01_run_tensorqtl",
                                     description="Map QTLs using TensorQTL")
    # output directory
    parser.add_argument('-output',
                        type=str, help="Directory path to save output files",
                        required=True)
    # output filename prefix
    parser.add_argument('-prefix',
                        type=str, help="Prefix for output file names",
                        required=False, default = "cis_qtls")
    # input files
    parser.add_argument('-pheno',
                        type=str, help="Phenotype input file",
                        required=True)
    parser.add_argument('-covariates',
                        type=str, help="Covariates input file",
                        required=True)
    parser.add_argument('-plink',
                        type=str, help="Prefix for PLINK .bim/.fam/.bed genotype files",
                        required=True)
    # cis window size
    parser.add_argument('-window',
                        type=int, help="Window size for cis-QTL mapping",
                        required=False,
                        default=1e6
                        )

    # get arguments
    args = parser.parse_args()

    # use a dictionary
    args_dict = {"output_dir": args.output,
                 "pheno_file": args.pheno,
                 "covs_file" : args.covariates,
                 "plink_prefix_path": args.plink,
                 "window" : args.window,
                 "output_prefix" : "{}/{}".format(args.output, args.prefix)
                 }

    print("Saving output files to", args_dict['output_dir'])

    # print(args_dict)
    return(args_dict)


def load_input(phenotype_bed_file, covariates_file, plink_prefix_path):
    # --- load phenotypes and covariates
    print("--> loading phenotypes", datetime.now())
    print(phenotype_bed_file)
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)

    print("phenotype_df:", phenotype_df.shape)
    print(phenotype_df.head())

    print("phenotype_pos_df:", phenotype_pos_df.shape)
    print(phenotype_pos_df.head())


    print("--> loading covariates", datetime.now())
    print(covariates_file)
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)

    print("covariates_df:", covariates_df.shape)
    print(covariates_df.head())


    # --- PLINK reader for genotypes
    print("--> reading genotypes", datetime.now())
    print(plink_prefix_path)
    pr = genotypeio.PlinkReader(plink_prefix_path)
    print("--> loading genotypes", datetime.now())
    genotype_df = pr.load_genotypes()
    print("--> loading variants", datetime.now())
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]


    data_dict = {'phenotype_df':phenotype_df,
            'phenotype_pos_df':phenotype_pos_df,
            'covariates_df':covariates_df,
            'genotype_df':genotype_df,
            'variant_df':variant_df
            }

    return(data_dict)



def cis_mapping(data_dict, cis_window, output_prefix, output_dir):
    # ---  cis-QTL: nominal p-value for all variant-phenotype pairs
    print("--> mapping cis-QTLs", datetime.now())

    genotype_df = data_dict['genotype_df']
    variant_df = data_dict['variant_df']
    phenotype_df = data_dict['phenotype_df']
    phenotype_pos_df = data_dict['phenotype_pos_df']
    covariates_df = data_dict['covariates_df']

    # map all cis-associations (results for each chrom written to file)

    print("cis window: {}".format(cis_window))

    # get summary statistics for all variant-phenotype pairs
    # this outputs to the parquet files for each chromosome
    # maps nominal associations for all variant-phenotype pairs
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, output_prefix, covariates_df, window=cis_window)


    # --- cis-QTL: empirical p-values for phenotypes
    # permutes data to get empirical p-values,
    # enabling calculation of genome-wide FDR
    # output contains top variant for each outcome phenotype
    print("--> mapping cis-QTLs empirical p-values", datetime.now())
    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, seed=123456, window=cis_window)

    # get the nominal associations for all variant-phenotype pairs
    tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85, fdr=0.05)
    out_file = os.path.join("{}.cis_qtl.txt.gz".format(output_prefix))
    cis_df.to_csv(out_file, sep='\t', float_format='%.6g')


    # save the summary stats for the permutation runs
    cis_df.to_csv("{}/cis_qtl_summary_stats.csv".format(output_dir))



def format_output(output_prefix, output_dir):
    # read the cis qtls file
    datafile = "{}.cis_qtl.txt.gz".format(output_prefix)
    print("Reading file", datafile)

    # read in the data file
    # qtls = pd.read_parquet(datafile)
    qtls = pd.read_csv(datafile, sep='\t', index_col=0)
    print(qtls.head())
    # get significant qtl pairs
    # nominal_files = "{}/all_chroms_hg38.cis_qtl_pairs".format(output_dir)
    # nominal_files = "{}/cis_qtls.cis_qtl_pairs".format(output_dir)
    # print(nominal_files)
    signif_df = get_significant_pairs(qtls, output_prefix)
    print(signif_df.head())

    signif_df.to_csv("{}/cis_qtl.signif_pairs.csv".format(output_dir))




def susie_finemapping(data_dict):
    # this function is not fully implemented
    # use DAPG for fine-mapping instead 
    print("--> mapping cis-QTLs with Susie", datetime.now())

    genotype_df = data_dict['genotype_df']
    variant_df = data_dict['variant_df']
    phenotype_df = data_dict['phenotype_df']
    phenotype_pos_df = data_dict['phenotype_pos_df']
    covariates_df = data_dict['covariates_df']

    # susie needs the output from running the cis mapping with q-values
    signif_df = pd.read_csv("{}/cis_qtls.csv".format(output_dir))
    print(signif_df.head())

    ix = phenotype_df.index[phenotype_df.index.isin(signif_df['phenotype_id'].unique())]
    summary_dict = susie.map(genotype_df, variant_df,
            phenotype_df.loc[ix], phenotype_pos_df.loc[ix],
            covariates_df, max_iter=500) # the summary_only parameter doesn't exist??

    # convert output dictionary to data fame
    summary_df = susie.get_summary(summary_dict)

    # save as binary file
    summary_df.to_parquet(os.path.join(output_dir, 'susie_summary.parquet'))

    # save as csv - don't do this if the dataset if really huge
    summary_df.to_csv("{}/susie_summary.csv".format(output_dir), index=False)



def map_trans():
    # --- trans-QTL mapping
    # skipping trans QTL mapping for now
    print("--> mapping trans-QTLs empirical p-values", datetime.now())
    # run mapping
    # to limit output size, only associations with p-value <= 1e-5 are returned
    trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df, batch_size=10000,
                               return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)

    # remove cis-associations
    trans_df = trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df, window=5000000)

    trans_df.to_csv("{}/trans_qtls.csv".format(output_dir))


def main():
    start = datetime.now()

    # get command line arguments
    args_dict = get_args()

    # load the input data
    data_dict = load_input(args_dict['pheno_file'], args_dict['covs_file'], args_dict['plink_prefix_path'])

    # map cis qtls
    cis_mapping(data_dict, args_dict['window'], args_dict['output_prefix'], args_dict['output_dir'])

    # format output
    format_output(args_dict['output_prefix'], args_dict['output_dir'])

    # fine-mapping with susie
    # susie_finemapping(data_dict)

    print("--> Finished mapping", datetime.now())
    end = datetime.now()
    duration = end - start
    print(duration)

    # print the time to a file
    timefile = "{}.time_to_run.txt".format(args_dict['output_prefix'])
    with open(timefile, 'w+') as f:
        f.write("{}\n".format(duration))
    

if __name__ == '__main__':
    main()
