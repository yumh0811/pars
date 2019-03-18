#!/usr/bin/env Rscript

library(plyr)
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
# name <- "Scer_n7_Spar"
name <- opt$name
path <- paste0("~/data/mrna-structure/result/", name, collapse = NULL)
setwd(path)

# snp所有信息
data_SNPs_total_info_vep_non_overlapped <- read.csv(paste0(path, "/", name, "_SNPs_total_info_vep_non_overlapped.tsv",collapse = NULL),  header = T, sep = "\t")

## 统计有PARS数据的转录本中的SNPs
file_SNPs_fold_info <- paste0("~/data/mrna-structure/process/", name, ".gene_variation.var_pars.tsv", collapse = NULL)
data_SNPs_fold_info <- read.csv(file_SNPs_fold_info,  header = TRUE, sep = "\t")
data_SNPs_fold_info <- subset(data_SNPs_fold_info, select=-c(gene))
colnames(data_SNPs_fold_info) [1] <- "location"

## 得到每一个有PARS信息的转录本的茎环长度、GC含量等信息
file_fold_class <- paste0("~/data/mrna-structure/result/",name,"/",name, ".gene_variation.fold_class.csv", collapse = NULL)
data_fold_class <- read.csv(file_fold_class,  header = TRUE, sep = ",")

data_SNPs_total <- merge(data_SNPs_total_info_vep_non_overlapped, data_fold_class, by="gene") # check in process/fail_pos.txt "genes don't match length"
data_SNPs_total <- merge(data_SNPs_fold_info, data_SNPs_total,  by="location") # check in process/fail_pos.txt "SNPs in overlapped Gene"

rm(data_SNPs_total_info_vep_non_overlapped,  data_SNPs_fold_info, data_fold_class)

#consequence类型
consequence <- unique(data_SNPs_total["consequence"])
rownames(consequence) <- NULL

# 求出mRNA中SNPs数量和mRNA数量
dd_SNPs <- data.frame(name="mRNA", SNPs=nrow(data_SNPs_total))
dd_gene <- data.frame(name="mRNA", gene=nrow(unique(data_SNPs_total["gene"],fromLast=TRUE)))

# 取出mRNA中CDS，求出SNPs数量和mRNA数量
data_SNPs_cds <- sqldf('SELECT * FROM [data_SNPs_total] where CDS_position != "-"')
dd_SNPs <- rbind(dd_SNPs, data.frame(name="CDS", SNPs=nrow(data_SNPs_cds)))
dd_gene <- rbind(dd_gene, data.frame(name="CDS", gene=nrow(unique(data_SNPs_cds["gene"],fromLast=TRUE))))

# 取出mRNA中UTR，求出SNPs数量和mRNA数量
data_SNPs_utr <- sqldf('SELECT * FROM [data_SNPs_total] where CDS_position = "-"')
dd_SNPs <- rbind(dd_SNPs, data.frame(name="UTR", SNPs=nrow(data_SNPs_utr)))
dd_gene <- rbind(dd_gene, data.frame(name="UTR", gene=nrow(unique(data_SNPs_utr["gene"],fromLast=TRUE))))

# 取出syn，求出SNPs数量和mRNA数量
data_SNPs_syn <- sqldf('SELECT * FROM [data_SNPs_total] where Consequence == "stop_retained_variant" OR Consequence == "synonymous_variant"')
dd_SNPs <- rbind(dd_SNPs, data.frame(name="SYN", SNPs=nrow(data_SNPs_syn)))
dd_gene <- rbind(dd_gene, data.frame(name="SYN", gene=nrow(unique(data_SNPs_syn["gene"],fromLast=TRUE))))

# 取出nsy，求出SNPs数量和mRNA数量
data_SNPs_nsy <- sqldf('SELECT * FROM [data_SNPs_total] where Consequence == "missense_variant" OR Consequence == "start_lost" OR Consequence == "stop_gained" OR Consequence == "stop_lost"')
dd_SNPs <- rbind(dd_SNPs, data.frame(name="NSY", SNPs=nrow(data_SNPs_nsy)))
dd_gene <- rbind(dd_gene, data.frame(name="NSY", gene=nrow(unique(data_SNPs_nsy["gene"],fromLast=TRUE))))

write.csv(data_SNPs_total, file=paste0(path, "/data_SNPs_PARS_mRNA.csv",collapse = NULL),row.names = FALSE)
write.csv(data_SNPs_cds, file=paste0(path, "/data_SNPs_PARS_cds.csv",collapse = NULL),row.names = FALSE)
write.csv(data_SNPs_utr, file=paste0(path, "/data_SNPs_PARS_utr.csv",collapse = NULL),row.names = FALSE)
write.csv(data_SNPs_syn, file=paste0(path, "/data_SNPs_PARS_syn.csv",collapse = NULL),row.names = FALSE)
write.csv(data_SNPs_nsy, file=paste0(path, "/data_SNPs_PARS_nsy.csv",collapse = NULL),row.names = FALSE)
write.csv(dd_SNPs, file=paste0(path, "/dd_SNPs.csv",collapse = NULL),row.names = FALSE)
write.csv(dd_gene, file=paste0(path, "/dd_gene.csv",collapse = NULL),row.names = FALSE)
