name <- "Scer_n157_nonMosaic_Spar" 
path <- paste0("~/data/mrna-structure/result/", name, collapse = NULL)
setwd(path)

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
library(sm)
library(vioplot)

#snp
file_cds_per_gene_ATGC <- paste0(path,"/data_SNPs_PARS_cds_per_gene_ATGC.csv",collapse = NULL)
data_cds_per_gene_ATGC <- read.csv(file_cds_per_gene_ATGC, header = TRUE, sep = ",")
colnames(data_cds_per_gene_ATGC) <- c("gene","cds_stem_AT_GC","cds_stem_GC_AT","cds_stem_total","cds_loop_AT_GC","cds_loop_GC_AT","cds_loop_total")

#length
file_cds_per_gene_length <- "~/data/mrna-structure/phylogeny/stem_loop_cds_length.csv"
data_cds_per_gene_length <- read.csv(file_cds_per_gene_length, header = TRUE, sep = ",")
colnames(data_cds_per_gene_length) <- c("gene","chr","cds_stem_length","cds_loop_length")

data_cds_snp_length <- merge(data_cds_per_gene_length,data_cds_per_gene_ATGC,by="gene")

file_info <- paste0(path,"/",name,".gene_variation.fold_class_cds.csv",collapse = NULL)
data_info <- read.csv(file_info , header = TRUE, sep = ",")

data_info <- merge(data_info ,data_cds_snp_length,by="gene")
write.csv(data_info,file ="data_info_cds.csv",row.names = FALSE)


#snp
file_syn_per_gene_ATGC <- paste0(path,"/data_SNPs_PARS_syn_per_gene_ATGC.csv",collapse = NULL)
data_syn_per_gene_ATGC <- read.csv(file_syn_per_gene_ATGC, header = TRUE, sep = ",")
colnames(data_syn_per_gene_ATGC) <- c("gene","syn_stem_AT_GC","syn_stem_GC_AT","syn_stem_total","syn_loop_AT_GC","syn_loop_GC_AT","syn_loop_total")

#length
file_cds_per_gene_length <- "~/data/mrna-structure/phylogeny/stem_loop_cds_length.csv"
data_cds_per_gene_length <- read.csv(file_cds_per_gene_length, header = TRUE, sep = ",")
colnames(data_cds_per_gene_length) <- c("gene","chr","syn_stem_length","syn_loop_length")

data_syn_snp_length <- merge(data_cds_per_gene_length,data_syn_per_gene_ATGC,by="gene")

file_info <- paste0(path,"/",name,".gene_variation.fold_class_cds.csv",collapse = NULL)
data_info <- read.csv(file_info , header = TRUE, sep = ",")

data_info <- merge(data_info ,data_syn_snp_length,by="gene")
write.csv(data_info,file ="data_info_syn.csv",row.names = FALSE)

#snp
file_nsy_per_gene_ATGC <- paste0(path,"/data_SNPs_PARS_nsy_per_gene_ATGC.csv",collapse = NULL)
data_nsy_per_gene_ATGC <- read.csv(file_nsy_per_gene_ATGC, header = TRUE, sep = ",")
colnames(data_nsy_per_gene_ATGC) <- c("gene","nsy_stem_AT_GC","nsy_stem_GC_AT","nsy_stem_total","nsy_loop_AT_GC","nsy_loop_GC_AT","nsy_loop_total")

#length
file_cds_per_gene_length <- "~/data/mrna-structure/phylogeny/stem_loop_cds_length.csv"
data_cds_per_gene_length <- read.csv(file_cds_per_gene_length, header = TRUE, sep = ",")
colnames(data_cds_per_gene_length) <- c("gene","chr","nsy_stem_length","nsy_loop_length")

data_nsy_snp_length <- merge(data_cds_per_gene_length,data_nsy_per_gene_ATGC,by="gene")

file_info <- paste0(path,"/",name,".gene_variation.fold_class_cds.csv",collapse = NULL)
data_info <- read.csv(file_info , header = TRUE, sep = ",")

data_info <- merge(data_info ,data_nsy_snp_length,by="gene")
write.csv(data_info,file ="data_info_nsy.csv",row.names = FALSE)