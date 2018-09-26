#!/usr/bin/env Rscript

library(plyr)
library(getopt)
library(ape)

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
area <- opt$area
path <- paste0("~/data/mrna-structure/result/", name, collapse = NULL)
setwd(path)
file_pars <- paste0("data_SNPs_PARS_",area,".csv", collapse = NULL)
snp_pars <- read.csv(file_pars, head= TRUE, sep = ',')

vcf.list <- paste0("~/data/mrna-structure/vcf/1011Matrix.gvcf/", name, opt$outfile, "/", name, opt$outfile, ".", area, "_snp.merge.pro.txt", collapse = NULL)
snp_vcf <- read.csv(vcf.list, head= TRUE, sep = ',')

data_update <- merge(snp_pars, snp_vcf, by="name")
data_update_non <- sqldf('SELECT * FROM [snp_pars] EXCEPT SELECT * FROM [data_update]')

dd_SNPs <- data.frame(name = "update", SNPs = c(nrow(data_update)))
dd_SNPs <- rbind(dd_SNPs, data.frame(name = "update_non", SNPs = nrow(data_update_non)))

write.csv(data_update,file=paste0(path, "/data_SNPs_PARS_", area, ".update", opt$outfile, ".csv", collapse = NULL),row.names = FALSE)
write.csv(data_update_non,file=paste0(path, "/data_SNPs_PARS_", area, ".update_non", opt$outfile, ".csv", collapse = NULL),row.names = FALSE)
write.csv(dd_SNPs,file=paste0(path, "/stat_", area, ".update", opt$outfile, ".csv", collapse = NULL),row.names = FALSE)

