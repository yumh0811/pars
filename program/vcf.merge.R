#!/usr/bin/env Rscript

library(plyr)
library(getopt)
library(ape)

spec = matrix(
    c(
        "help",
        "h",
        0,
        "logical",
        "brief help message",

        "name",
        "n",
        1,
        "character",
        "input name",
        
        "area",
        "a",
        1,
        "character",
        "input area eg. cds/utr/syn/nsy",
        
        "outfile",
        "o",
        1,
        "character",
        "output filename"
    ),
    byrow = TRUE,
    ncol = 5
)
opt = getopt(spec)

name <- opt$name
path <- paste0("~/data/mrna-structure/vcf/1011Matrix.gvcf/", name, collapse = NULL)
setwd(path)
file_pars <- paste0("data_SNPs_PARS_",opt$area,".pars.tsv", collapse = NULL)
snp_pars <- read.csv(file_pars, head= TRUE, sep = '\t')
snp_ext <- read.csv("1011Matrix.ext.tsv", head= TRUE, sep = '\t')

snp_merge <- merge(snp_pars, snp_ext, by = "name", all = F)
write.table (snp_merge, file = paste0(opt$area,"_snp.merge.tsv", collapse = NULL), row.names = F,sep = "\t")
