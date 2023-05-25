# create input files for ENLOC

import os
import gzip
import pandas as pd
from datetime import datetime
import sys
import argparse
import subprocess

limit = -1



def get_args():
    # get command line arguments
    parser = argparse.ArgumentParser()
    # required parameters
    parser.add_argument('-o', '--output', type=str, required=True, help="Output directory")
    parser.add_argument('-t', '--torus_dir', type=str, required=True, help="Directory containing Torus priors for all genes")
    parser.add_argument('-g', '--geno_file', type=str, required=True, help="Genotype reference counts file")
    parser.add_argument('-p', '--pheno_file', type=str, required=True, help="Phenotype file")
    parser.add_argument('-c', '--covs_file', type=str, required=True, help="File with formatted covariant data")

    # optional parameters
    parser.add_argument('-i', '--index', type=int, required=False, default=-1, help="Index of gene to process")

    args = parser.parse_args()

    return(args)

def process_phenotype(savefile, gene, phenofile):
    # print("Processing phenotype data", datetime.now())
    # read in the phenotype data

    # grab lines that we care about
    cmd = "zcat {} | grep -v '##' | grep '#|{}'".format(phenofile, gene)
    process = subprocess.Popen(cmd.split(),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    lines = stdout.decode('utf-8').split('\n')
    header = lines[0].split('\t')
    line = lines[1].split('\t')

    # get the fields for output
    field1 = "pheno"
    pheno_id = line[3]
    group_id = "group"

    # create output
    outlist = [field1, pheno_id, group_id] + line[4:]
    outline = " ".join(outlist)

    # output to savefile in appropriate format
    out = open(savefile, "a+")
    out.write("{}\n".format(outline)) 
    out.close()

    # use this to match order of samples for geno and controlled pheno
    ordered_samples = header[4:]

    return(ordered_samples)


def get_genes(savedir, phenofile):
    # get a list of all the genes
    print("Getting gene data", datetime.now())

    genesdf = pd.read_csv(phenofile, sep='\t')

    # get the genes
    geneslist = genesdf.gene_id.values

    # output genes list
    genesfile = "{}/genes_list.txt".format(savedir)
    genesonly = genesdf.gene_id
    genesonly.to_csv(genesfile, index=False, header=False)

    return(geneslist)
    

def process_genotype(sbamfile, priorsfile, ordered_samples, gene, genotypes_dict, missing_priors_file, torus_dir):
    # print("Processing genotype data", datetime.now())

    # get list of snps associated with this gene
    snpfile = "{}/{}.prior".format(torus_dir, gene)

    if not os.path.exists(snpfile):
        print("no priors file")
        with open(missing_priors_file, 'a+') as f:
            f.write("{}\n".format(gene))
        return()


    snpdf = pd.read_csv(snpfile, delim_whitespace=True, names=['pheno', 'pip'])
    snplist = snpdf.pheno.values


    # get the total number of spns
    numsnps = len(snplist)

    # print("num snps", numsnps)
    # print(snplist[1:100])

    out = open(sbamfile, "a+")

    header = []
    samples_idx = []

    notfound=0
    snpcount = 0
    for snp in snplist:
        # get the key for genotypes dict
        snpcount += 1

        if not snp in genotypes_dict:
            notfound+=1
            continue

        outline = genotypes_dict[snp]
        out.write("{}\n".format(outline)) 


    out.close()

    # print("not found", notfound)

    # write a new snp pip file
    # originally used to remove SNPs but don't need to do that
    # need to rewrite file for some reason or else error in dap-g occurs
    # snpdf.index = snpdf.pheno.values
    # keepdf = snpdf.loc[keepsnps,:]
    snpdf.to_csv(priorsfile, sep=' ', index=False, header=False)


def process_variables(savefile, ordered_samples, covs_file):
    # print("Processing covariates", datetime.now())
    # read in the variable data
    # output to savefile in appropriate format

    out = open(savefile, "a+")

    # need to transpose data so read in with pandas
    covs = pd.read_csv(covs_file, sep='\t')

    # samples alraedy match ordering --
    
    # need to save each column to output file
    cols = covs.columns
    for i in range(1, covs.shape[1]):
        # grab the column
        col = covs.iloc[:, i].values
        colname = cols[i]

        # create the output
        outlist = ['controlled', colname, 'group'] + list(col)

        outline = " ".join([str(x) for x in outlist])

        out.write("{}\n".format(outline))


    out.close()

def load_genotypes(genofile):
    # load the genotypes into a dictionary 
    # for quicker look-ups

    # read in the genotype data
    # output to savefile in appropriate format
    header = []
    samples_idx = []

    genotypes_dict = {}

    round = 0
    with open(genofile, 'r') as f:
        for line in f:
            # split on tab
            line = line.strip().split('\t')

            # store the header line
            if round == 0:
                header = line

                # remove 0_ from the sample names
                header = [x.replace('0_', '') for x in header]

                # get the index of the ordered samples in header
                # samples_idx = [header.index(x) for x in ordered_samples]
                round += 1
                continue



            # get the fields for output
            field1 = "geno"
            pheno_id = line[1]
            group_id = "group"

            # switch counts from ref to alt by subtracting from 2
            genos = line[6:]
            genos = [str(2-int(x)) if x != 'NA' else x for x in genos]


            # create output
            outlist = [field1, pheno_id, group_id] + genos
            outline = " ".join(outlist)

            genotypes_dict[pheno_id] = outline

            if round % 1000000 == 0:
                print("Loaded {} genotypes {}".format(round, datetime.now()))

            # stop for testing
            # if round == 10:
                # print("USING TESTING MODE FOR ONLY 10 GENES")
                # break
            round += 1

    return(genotypes_dict)


def main():
    start = datetime.now()
    # get command line args
    args = get_args()

    # create the output directories
    # savedir = "/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/03_dapg_input"
    savedir = args.output
    if not os.path.exists(savedir):
        print("creating output directory", savedir)
        os.makedirs(savedir)

    priorsdir = "{}/priors/".format(savedir)
    if not os.path.exists(priorsdir):
        print("creating output directory", priorsdir)
        os.makedirs(priorsdir)

    sbamdir = "{}/sbams/".format(savedir)
    if not os.path.exists(sbamdir):
        print("creating output directory", sbamdir)
        os.makedirs(sbamdir)


    # create this file to keep track of gene with missing priors
    missing_priors_file = "{}/missing_priors_file.txt".format(savedir)
    # with open(missing_priors_file, 'w+') as f:
        # None

    # need to output each gene separately - get list of genes
    geneslist = get_genes(savedir, args.pheno_file)
    print("Num genes:", len(geneslist))

    if args.index == -1:
        print("Enter a gene number")

    # select the gene to process
    gene = geneslist[args.index]
    print("Processing gene", args.index, gene)


    print("Loading genotype {}".format(datetime.now()))
    genotypes_dict = load_genotypes(args.geno_file)

    # write an output file for each gene

    # create the savefile for output
    sbamfile = "{}/{}_fastqtl_singl_snp_output.sbam".format(sbamdir, gene)
    # use this to not overwrite any saved files
    # if os.path.isfile(sbamfile):
        # return

    # clear the file after testing
    with open(sbamfile, "w+"):
        None

    # create priors output file
    priorsfile = "{}/{}_priors.txt".format(priorsdir, gene)
    with open(sbamfile, "w+"):
        None


    # process pheno
    # print("Processing pheno")
    ordered_samples = process_phenotype(sbamfile, gene, args.pheno_file)

    # process geno
    # print("Processing geno")
    process_genotype(sbamfile, priorsfile, ordered_samples, gene, genotypes_dict, missing_priors_file, args.torus_dir)

    # process controlled variables
    # print("Processing controlled variables")
    process_variables(sbamfile, ordered_samples, args.covs_file)

    print("--> Finished mapping", datetime.now())
    end = datetime.now()
    duration = end - start
    print(duration)


if __name__ == '__main__':
    main()
