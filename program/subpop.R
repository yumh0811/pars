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

data_snp_list <- read.csv("genelist.csv",header = T, sep = ",")

cds_snp <- read.csv("/Users/yumh/data/mrna-structure/vcf/1011Matrix.gvcf/Scer_n128_Spar.wild/Scer_n128_Spar.wild.cds_snp.merge.pro.tsv",header = T, sep = "\t")
merge_snp <- merge(cds_snp,data_snp_list,by="gene")
merge_snp_ATCG <- subset(merge_snp, structure=="stem")
merge_snp_ATCG <- subset(merge_snp_ATCG, mutant_to_pars == "A->G"|mutant_to_pars == "T->G"|mutant_to_pars == "A->C"|mutant_to_pars == "T->C")

write.csv(merge_snp_ATCG,file = "merge_snp.csv",row.names = F)

snp <- read.csv("total_snp.csv",header=T,sep=",")
mvar <- read.csv("/Users/yumh/data/mrna-structure/xlsx/Scer_n128_Spar.mvar.gene_list.csv",header=T,sep=",")
a <- merge(mvar,snp,by="snp_id")

b <- subset(a,select = c("name","snp_freq","snp_outgroup_base","snp_all_bases","snp_mutant_to"))

c <- merge(merge_snp_ATCG,b,by="name")

write.csv(c,file = "filiter_snp.csv",row.names = F)
