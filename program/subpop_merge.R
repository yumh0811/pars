#!/usr/bin/env Rscript

library(plyr)
library(getopt)
library(ape)

library(plyr)
library(getopt)
library(ape)
library(reshape2)
library(ggplot2)
library(scales)
library(reshape)
library(pander)
library(gridExtra)
library(dplyr)
library(proto)
library(gsubfn)
library(RSQLite)
library(sqldf)

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
        
        "infile",
        "i",
        1,
        "character",
        "input filename",
       
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
path <- paste0("~/data/mrna-structure/result/", name, "/subpop", collapse = NULL)
setwd(path)

file_syn_snp <- paste0("~/data/mrna-structure/result/", name, "/data_SNPs_PARS_syn.update_codon.csv", collapse = NULL)
data_syn_snp <- read.csv2(file_syn_snp,header=T,sep=",")
data_subpop_info <- read.csv2("subpop.csv",header=T,sep=",")
C <- merge(data_syn_snp, data_subpop_info, by="name")
write.csv(C, paste0("~/data/mrna-structure/result/", name, "/subpop/SNPs_syn_mt.csv", collapse = NULL), row.names = F)
