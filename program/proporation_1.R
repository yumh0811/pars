#!/usr/bin/env Rscript

setwd("~/data/mrna-structure/gene_filiter/")
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

        "infile",
        "i",
        1,
        "character",
        "input filename",
        
        "range",
        "r",
        1,
        "character",
        "protein coding range chr",
        
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

data_gene_range = read.csv(opt$infile,head=TRUE)
data_cds_range_chr = read.table(opt$range,head=FALSE,sep = ",")
colnames(data_cds_range_chr) = c("gene","chr","start","end")
data_gene_range = merge(data_cds_range_chr,data_gene_range,by="gene")
data_out <- arrange(data_gene_range,data_gene_range$chr,data_gene_range$start)
data_out_1 = subset(data_out, data_out$proporation==1)
write.csv(data_out_1,file=opt$outfile,row.names = FALSE)
