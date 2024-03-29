#!/home/eulalio/micromamba/bin/r

# this script was taken from the PTWAS github

require(data.table)

###################################################################
# OPTIONS
###################################################################

# output file name
OUT_FILE <- "PTWAS_weights.vcf"

# input directory containing $tissue.twb.gz files for each tissue
INPUT_DIR <- "UNDEF"
INPUT_FILE <- "UNDEF"

# no. sig figs for weights
n_sig_digits <- 4

# compress the combined weight file?
BGZIP <- FALSE

# create tabix index for combined weight file?
INDEX <- FALSE

OUT_FILE <- "PTWAS.db"

args <- commandArgs(trailingOnly = TRUE)
for  (i in 1:length(args)){
    
    if(args[i] == "-d"){
        i = i+1
        INPUT_FILE = args[i]
    }
    
    if(args[i] == "-D"){
        i = i+1
        INPUT_DIR = args[i]
    }
    
    if(args[i] == "-o"){
        i= i+1
        OUT_FILE = args[i]
    }

}

OUT_FILE_US = paste0(OUT_FILE,".us")

get_tissue<-function(in_file){
    return(fread( cmd = paste0("zcat ", in_file, " | head -1 | awk '{print $1 }' "), header=FALSE)$V1)
}


if(INPUT_FILE == "UNDEF" & INPUT_DIR == "UNDEF"){
    cat("Error: no input file or directory ...\n")
    q()
}

if(INPUT_FILE != "UNDEF" & INPUT_DIR == "UNDEF"){
    db_files = c(INPUT_FILE)
}

if(INPUT_FILE == "UNDEF" & INPUT_DIR != "UNDEF"){
    db_files <- fread(cmd = paste0('ls ',INPUT_DIR,'/*.gz'), header = FALSE)$V1
}




###################################################################
# READ WEIGHT FILES 
###################################################################


# get unique tissues
cat('Found ',length(db_files),' tissue-specific weight files ... \n')
tissues <- sapply(db_files, get_tissue);


names(db_files) <- tissues 

cat("Processing input weight files ... \n")
d <- rbindlist(lapply(tissues, function(tissue){
        d <- fread( cmd = paste0("zcat ", db_files[tissue], " | sed -e 's/ /\\t/g' | cut -f2- | sed -e 's/ \\+/_/g' | tr '_' '\\t' | cut -f1-5,7-8 | awk '{if( $6 * $6 > 0 ) print $0 }' "), header = FALSE )
	setnames(d, 
		c('gene', 'chr', 'pos', 'ref', 'alt', 'beta', 'z')
	)
	d$tissue <- tissue
	setcolorder(d, c('chr', 'pos', 'ref', 'alt', 'gene', 'tissue', 'beta', 'z'))
	
	cat('Processed ',tissue, ' ... \r')
	
	d
}))

# convert tissue IDs to integer keys (optional)
d$tissue_f <- as.integer(
	factor(d$tissue, levels = tissues, ordered = TRUE)
)-1

# round to desired significant figures
d[,bb:=signif(beta,n_sig_digits),]

# collapse weights across tissues for each eGene:eVariant
cat('Aggregating tissues for each eGene:eVariant ... \n')
dd <- d[,list(
	'beta'=paste0(gene[1],'=',paste(bb,tissue_f,sep='@',collapse=';'))
),by=list(chr,pos,ref,alt,gene)]

# collapse weights across eGenes for each eVariant
cat('Aggregating eGenes for each variant ... \n')
dd <- dd[,list(
	'beta'=paste(beta, collapse = '|'), 
	'id'=paste(chr[1],pos[1],ref[1],alt[1],sep=':')
),by=list(chr,pos,ref,alt)][order(chr,pos,ref,alt)]

# set order of columns
setcolorder(dd, c('chr','pos','id','ref','alt','beta'))

# header linking tissue names to integer keys
header_idx <- paste0('## SUBCLASS_IDS=', paste(sort(unique(d$tissue_f)), tissues, sep = ':', collapse = ','))
header_vcf <- c("#chr\tpos\tid\tref\talt\tbeta")
# write header to output file
writeLines(text = header_idx, OUT_FILE_US)
writeLines(text = c(header_idx, header_vcf), OUT_FILE)
# add hash before column IDs
setnames(dd, 'chr', '#chr')

# write output to file
cat('Writing output ... \n')

fwrite(dd, OUT_FILE_US, append = TRUE, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

system(paste("sort -V -k1,1 -k2,2n ", OUT_FILE_US, " >> ",OUT_FILE))
system(paste("rm ", OUT_FILE_US))

# bgzip output file
if( BGZIP ){
	cat('Compressing output file ... \n')
	system(paste("bgzip", OUT_FILE))
}

# create tabix index
if( INDEX ){
	cat('Making tabix index ... \n')
	system(paste0("tabix -p 'vcf' ", OUT_FILE, ".gz"))
}

