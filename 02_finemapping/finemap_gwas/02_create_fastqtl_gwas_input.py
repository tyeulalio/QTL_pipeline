# create the gwas input for fastqtl

import os
from datetime import datetime
import pandas as pd
import gzip
import argparse


def get_args():
    # get command line arguments
    parser = argparse.ArgumentParser()
    # required parameters
    parser.add_argument('-o', '--output', type=str, required=True, help="Output directory")
    parser.add_argument('-d', '--dap_dir', type=str, required=True, help="Directory containing fine-mapped DAP-G output")

    args = parser.parse_args()

    return(args)



def get_genes(dap_dir):
    # read in the filenames from the dap-g output
    print("Reading fine-mapped files from {}".format(dap_dir))
    gene_files = os.listdir(dap_dir)
    gene_files = [f for f in gene_files if os.path.isfile("{}/{}".format(dap_dir, f))]
    ld_files = [f for f in gene_files if f.startswith('ld')]

    genes = [x.replace('.dapg', '') for x in ld_files]

    return(genes)

def write_gene_output(gene, out, dap_dir):
    # read in the gwas pips data
    # output in proper format to file
    dap_file = "{}/{}.dapg".format(dap_dir, gene)

    gene = gene.replace(".", "-")

    i = 0
    with open(dap_file, 'r') as f:
        for line in f:
            line = line.strip()

            # keep only the SNP data
            if not '((' in line:
                continue

            # split to the line to a list
            # grab the fields we need
            line = line.split()
            # print(line)

            _,snp,pip,_,signal_id,_,_=line
            # print(snp, pip, signal_id)

            # add gene name before SNP
            gene_snp = "{}_{}".format(gene, snp)

            # add gene name before signal id
            if signal_id == "-1":
                continue

            gene_signal_id = "g{}:{}".format(gene, signal_id)

            # create the output line
            outline = "\t".join([snp, gene_signal_id, pip])

            # write to output file
            out.write("{}\n".format(outline))

            i+=1


def main():

    # get command line arguments
    args = get_args()

    # create output directory if it doesn't exit
    print("Save output to", args.output)
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # get list of genes
    print("Loading GWAS genes {}".format(datetime.now()))
    genes = get_genes(args.dap_dir)

    # output file
    savefile = "{}/fastqtl_gwas_input_pips.vcf".format(args.output)
    out = open(savefile, 'w+')

    print("Processing {} genes {}".format(len(genes), datetime.now()))

    # write each gene's output to file
    j = 0
    for gene in genes:
        write_gene_output(gene, out, args.dap_dir)

        if j % 1000 == 0:
            print("Processed {} genes {}".format(j, datetime.now()))
        j += 1

    out.close()

    
    # gzip the file
    zipfile = "{}.gz".format(savefile)
    with open(savefile, 'rb') as f_in, gzip.open(zipfile, 'wb') as f_out:
        f_out.writelines(f_in)


if __name__ == '__main__':
    main()
