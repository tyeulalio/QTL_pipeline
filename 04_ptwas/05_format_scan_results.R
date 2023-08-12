library(tidyverse)

# check the ptwas scan results
# this can be done by concatenating the separate chromosome files
# but you need to make sure the variants are sorted by position

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

ptwas_scan_dir <- args[1]
ptwas_weights_dir <- args[2]

print(paste("Reading PTWAS scan files from", ptwas_scan_dir))

testing=FALSE
if (testing){
    # directory containing ptwas scan results
    ptwas_scan_dir <- "../examples/example_output/05_ptwas/04_ptwas_scan"
    ptwas_weights_dir <- "../examples/example_output/05_ptwas/01_ptwas_weights/ptwas_weights/"
    chrom <- 18
}


process_chrom <- function(chrom){
    # read in ptwas results for each chromosomes

    # summary file
    datafile <- paste0(ptwas_scan_dir, "/chr", chrom, "_ptwas_scan.summary_out.txt")
    datafile
    if (!file.exists(datafile)){
        print(paste("No output for chromosome", chrom))
        return(NA)
    }

    datafile
    summary_res <- read_tsv(datafile, col_names=c('chrom', 'pos', 'gene', 'n_snps', 'n_classes', 'top_class', 'top_subclass',
                                                  'min_unadj_pval', 'naive_pval', 'pval', 'info'), skip=1,
                            show_col_types =FALSE
    )
    head(summary_res)

    # stratified file
    datafile <- paste0(ptwas_scan_dir, "/chr", chrom, "_ptwas_scan.stratified_out.txt")
    strat_res <- read_tsv(datafile, col_names=c('chrom', 'pos', 'gene', 'class', 'subclass', 'n_snps', 'stat', 'pval', 'info'),
                          skip=1,
                          show_col_types =FALSE
    )
    head(strat_res)

    strat_res %>%
        arrange(pval)

    head(summary_res)
    head(strat_res)

    dim(summary_res)
    dim(strat_res)

    return(list(strat_res=strat_res,
                summary_res=summary_res))
}

chrom_res <- lapply(1:22, process_chrom)
head(chrom_res)

# handle missing chromosomes
if (all(is.na(chrom_res))){
    print(paste("No output for all chromosomes"))
    quit()
}
full_res <- chrom_res[!is.na(chrom_res)]

# combine the data together
combined_strat <- do.call(rbind, lapply(full_res, function(x) x$strat_res))
head(combined_strat)

combined_summary <- do.call(rbind, lapply(full_res, function(x) x$summary_res))
head(combined_summary)


nrow(combined_summary)
table(combined_summary$pval < 0.05)

# turn this on to check missing genes
check_missing = FALSE
if (check_missing){
    # load genes from previous step
    genes <- list.files(ptwas_weights_dir) %>%
        str_remove('_ptwas_weights.txt')
    head(genes)
    
    # check which genes are missing
    head(combined_summary)
    missing_genes <- setdiff(genes, unique(combined_summary$gene))
    length(missing_genes)
    head(missing_genes)
}

# save to output file
(savefile <- paste0(ptwas_scan_dir, "/all_chroms_ptwas_scan.summary_out.txt"))
write_tsv(combined_summary, savefile)

(savefile <- paste0(ptwas_scan_dir, "/all_chroms_ptwas_scan.stratified_out.txt"))
write_tsv(combined_strat, savefile)

print(paste("Num genes:", length(unique(combined_summary$gene))))
