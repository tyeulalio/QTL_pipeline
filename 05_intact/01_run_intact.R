library(INTACT)
#library(biomaRt)
library(tidyverse)

# run intact on the ptwas and colocalization output


run_intact <- function(ptwas, fastenloc, savedir){
    # run intact
    head(ptwas) 
    head(fastenloc)

    summary(ptwas$stat)

    # join the columns that we want
    sub_ptwas <- ptwas %>%
        select(gene, zscore=stat)
    head(sub_ptwas)

    # match the gene names 
    sub_fastenloc <- fastenloc %>%
        select(gene=Gene, GLCP) %>%
        mutate(gene = str_replace_all(gene, "-", "\\."))
    head(sub_fastenloc)

    ptwas_fastenloc <- sub_ptwas %>%
        inner_join(sub_fastenloc)
    head(ptwas_fastenloc)
    dim(ptwas_fastenloc)

    res <- intact(GLCP_vec=ptwas_fastenloc$GLCP,
                  z_vec=ptwas_fastenloc$zscore,
                  prior_fun=linear
    )
    head(res)
    summary(res)

    # combine results
    intact_res <- cbind(ptwas_fastenloc, intact_pip=res) %>%
        arrange(-intact_pip)
    head(intact_res)

    (savefile <- paste0(savedir, "all_intact_results.tsv"))
    print(paste("Saving output to", savefile))
    saveRDS(intact_res, savefile)
    1
}

main <- function(){
    # get command line arguments
    args <- commandArgs(trailingOnly = TRUE)

    # arguments
    savedir <- args[1]
    # ptwas stratefied file
    ptwas_file <- args[2]
    # fastenloc output
    fastenloc_file <- args[3]

    for_testing <- FALSE
    if (for_testing){
        savedir <- "../examples/example_output/06_intact/"

        ptwas_file <- paste0("../examples/example_output/05_ptwas/04_ptwas_scan/all_chroms_ptwas_scan.stratified_out.txt")
        fastenloc_file <- paste0("../examples/example_output/04_colocalization/01_run_fastenloc/fastqtl_.enloc.gene.out")
    }

    # create the output directory if it doesn't exist
    dir.create(savedir, showWarnings = FALSE)

    # load ptwas data
    ptwas <- read_table(ptwas_file)
    head(ptwas)

    # load fastenloc GLCP data
    fastenloc <- read_table(fastenloc_file)
    head(fastenloc)

    res <- run_intact(ptwas, fastenloc, savedir)
    print("Finished INTACT")
    1
}

main()
