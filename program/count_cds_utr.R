#!/usr/bin/env Rscript

library(getopt)
library(ape)
library(ggplot2)
library(scales)
library(reshape)
library(pander)
library(gridExtra)
library(plyr)
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

#name <- "Scer_n128_Spar"
name <- opt$name
path <- paste0("~/data/mrna-structure/result/", name, collapse = NULL)
setwd(path)

#snp
file_cds_per_gene_ATGC <- paste0(path,"/data_SNPs_PARS_cds.update_per_gene_ATGC.csv",collapse = NULL)
data_cds_per_gene_ATGC <- read.csv(file_cds_per_gene_ATGC, header = TRUE, sep = ",")
file_utr_per_gene_ATGC <- paste0(path,"/data_SNPs_PARS_utr.update_per_gene_ATGC.csv",collapse = NULL)
data_utr_per_gene_ATGC <- read.csv(file_utr_per_gene_ATGC, header = TRUE, sep = ",")
colnames(data_cds_per_gene_ATGC) <- c("gene","cds_stem_AT_GC","cds_stem_GC_AT","cds_stem_total","cds_loop_AT_GC","cds_loop_GC_AT","cds_loop_total")
colnames(data_utr_per_gene_ATGC) <- c("gene","utr_stem_AT_GC","utr_stem_GC_AT","utr_stem_total","utr_loop_AT_GC","utr_loop_GC_AT","utr_loop_total")

#length
file_cds_per_gene_length <- paste0(path,"/stem_loop_cds_length.csv",collapse = NULL)
data_cds_per_gene_length <- read.csv(file_cds_per_gene_length, header = TRUE, sep = ",")
file_utr_per_gene_length <- paste0(path,"/stem_loop_utr_length.csv",collapse = NULL)
data_utr_per_gene_length <- read.csv(file_utr_per_gene_length, header = TRUE, sep = ",")
colnames(data_cds_per_gene_length) <- c("gene","chr","cds_stem_length","cds_loop_length")
colnames(data_utr_per_gene_length) <- c("gene","chr","utr_stem_length","utr_loop_length")

data_cds_snp_length <- merge(data_cds_per_gene_length,data_cds_per_gene_ATGC,by="gene")
data_utr_snp_length <- merge(data_utr_per_gene_length,data_utr_per_gene_ATGC,by="gene")
data_utr_snp_length<- subset(data_utr_snp_length,select = c("gene","utr_stem_length","utr_loop_length","utr_stem_AT_GC","utr_stem_GC_AT","utr_stem_total","utr_loop_AT_GC","utr_loop_GC_AT","utr_loop_total"))
data_snp_length <- merge(data_cds_snp_length,data_utr_snp_length,by="gene",all=T)

file_info <- paste0(path,"/",name,".gene_variation.fold_class.csv",collapse = NULL)
data_info <- read.csv(file_info , header = TRUE, sep = ",")

data_info <- merge(data_info ,data_snp_length,by="gene")
write.csv(data_info,file ="data_info_cds_utr.update.csv",row.names = FALSE)
data_info<-na.omit(data_info)
write.csv(data_info,file ="data_info_cds_utr.update_overlap.csv",row.names = FALSE)


