name <- "Scer_n128_Seub" 
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

# 统计所有snp的数量
file_SNPs_all <- paste0("~/data/mrna-structure/xlsx/", name, ".mvar.gene_list.csv", collapse = NULL)
data_SNPs_all <- read.csv(file_SNPs_all,  header = TRUE, sep = ",")
dd_SNPs <- data.frame(name = c("all"),SNPs=c(nrow(data_SNPs_all)))
rownames(data_SNPs_all) <- NULL # suppress rownames

# 统计intergenic snp的数量
file_SNPs_intergenic <- paste0("~/data/mrna-structure/process/", name, ".snp.intergenic.pos.txt", collapse = NULL)
data_intergenic <- read.csv(file_SNPs_intergenic,  header = FALSE, sep = "\t")
names(data_intergenic) =c("name")
data_intergenic <- merge(data_SNPs_all, data_intergenic, by="name")
dd_SNPs <- rbind(dd_SNPs, data.frame(name="intergenic", SNPs=nrow(data_intergenic) ))
rownames(data_intergenic) <- NULL # suppress rownames

# 统计有PARS数据的转录本中的SNPs
file_SNPs_PARS_transcripts <- paste0("~/data/mrna-structure/process/", name, ".gene_variation.var_pars.tsv", collapse = NULL)
data_SNPs_PARS_transcripts <- read.csv(file_SNPs_PARS_transcripts,  header = TRUE, sep = "\t")

# 合并data_SNPs_all和data_SNPs_PARS_transcripts
data_SNPs_PARS_transcripts <- merge(data_SNPs_PARS_transcripts, data_SNPs_all, by="name")
dd_SNPs <- rbind(dd_SNPs, data.frame(name="PARS_transcripts", SNPs=nrow(data_SNPs_PARS_transcripts) ))

# 得到每一个有PARS信息的转录本的茎环长度、GC含量等信息
file_fold_class <- paste0("~/data/mrna-structure/result/",name,"/",name, ".gene_variation.fold_class.csv", collapse = NULL)
data_fold_class <- read.csv(file_fold_class,  header = TRUE, sep = ",")

# 得到protein coding gene的list
file_protein_coding_list <- "~/data/mrna-structure/phylogeny/protein_coding_list.csv"
data_protein_coding_list <- read.csv(file_protein_coding_list,  header = FALSE, sep = ",")
colnames(data_protein_coding_list) <- c("gene")

# 取出cds alignment proportation = 1 的基因
file_proporation_1_gene <- paste0("~/data/mrna-structure/phylogeny/",name,"_distance_processed_pro_1.csv", collapse = NULL)
data_proporation_1_gene <- read.csv(file_proporation_1_gene,  header = TRUE, sep = ",")
data_proporation_1_gene <- subset(data_proporation_1_gene,select = gene)
data_protein_coding_list <- merge(data_protein_coding_list,data_proporation_1_gene,by="gene")

# 将有PARS信息的转录本分为mRNA和非mRNA
data_fold_class_mRNA <- merge(data_fold_class, data_protein_coding_list , by="gene")
data_fold_class_non_mRNA <- sqldf('SELECT * FROM [data_fold_class] EXCEPT SELECT * FROM [data_fold_class_mRNA]') # 备用

# 得到有PARS信息的mRNA和非mRNA的SNPs
data_SNPs_PARS_mRNA <- merge(data_SNPs_PARS_transcripts , data_fold_class_mRNA , by="gene")
data_SNPs_PARS_non_mRNA <- merge(data_SNPs_PARS_transcripts , data_fold_class_non_mRNA , by="gene") # 备用

# 统计有PARS信息的mRNA中SNPs的数量
dd_SNPs <- rbind(dd_SNPs, data.frame(name="PARS_mRNA", SNPs=nrow(data_SNPs_PARS_mRNA) ))
# 统计有PARS信息的mRNA的数量
data_gene_process <- data_SNPs_PARS_mRNA["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
colnames(data_gene_process) <- c("PARS_mRNA_gene")
write.csv(data_gene_process,file=paste0(path,"/PARS_mRNA_gene.csv",collapse = NULL),row.names = FALSE)
dd_gene <- data.frame(name = c("PARS_mRNA"),gene=c(nrow(data_gene_process)))

# 去除有PARS信息的mRNA的SNPs中complex，求出SNPs数量和mRNA数量
data_SNPs_PARS_mRNA <- subset(data_SNPs_PARS_mRNA, data_SNPs_PARS_mRNA$mutant_to != "Complex")
write.csv(data_SNPs_PARS_mRNA,file=paste0(path,"/data_SNPs_PARS_mRNA.csv",collapse = NULL),row.names = FALSE)
dd_SNPs <- rbind(dd_SNPs, data.frame(name="PARS_mRNA_non_complex", SNPs=nrow(data_SNPs_PARS_mRNA) ))
data_gene_process <- data_SNPs_PARS_mRNA["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
colnames(data_gene_process) <- c("PARS_mRNA_non_complex_gene")
write.csv(data_gene_process,file=paste0(path,"/PARS_mRNA_non_complex_gene.csv",collapse = NULL),row.names = FALSE)
dd_gene <- rbind(dd_gene, data.frame(name="PARS_mRNA_non_complex", gene=nrow(data_gene_process) ))

# 取出mRNA中CDS，求出SNPs数量和mRNA数量
file_SNPs_cds <- paste0("~/data/mrna-structure/process/", name, ".snp.cds.pos.txt", collapse = NULL)
data_cds <- read.csv(file_SNPs_cds,  header = FALSE, sep = "\t")
names(data_cds) =c("name")
data_cds <- merge(data_SNPs_PARS_mRNA, data_cds, by="name")
dd_SNPs <- rbind(dd_SNPs, data.frame(name="PARS_cds", SNPs=nrow(data_cds) ))
rownames(data_cds) <- NULL # suppress rownames
data_gene_process <- data_cds["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
colnames(data_gene_process) <- c("PARS_cds_gene")
write.csv(data_gene_process,file=paste0(path,"/PARS_cds_gene.csv",collapse = NULL),row.names = FALSE)
dd_gene <- rbind(dd_gene, data.frame(name="PARS_cds", gene=nrow(data_gene_process) ))
write.csv(data_cds,file=paste0(path,"/data_SNPs_PARS_cds.csv",collapse = NULL),row.names = FALSE)

# 取出mRNA中UTR，求出SNPs数量和mRNA数量
file_SNPs_utr <- paste0("~/data/mrna-structure/process/", name, ".snp.utr.pos.txt", collapse = NULL)
data_utr <- read.csv(file_SNPs_utr,  header = FALSE, sep = "\t")
names(data_utr) =c("name")
data_utr <- merge(data_SNPs_PARS_mRNA, data_utr, by="name")
dd_SNPs <- rbind(dd_SNPs, data.frame(name="PARS_utr", SNPs=nrow(data_utr) ))
rownames(data_utr) <- NULL # suppress rownames
data_gene_process <- data_utr["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
colnames(data_gene_process) <- c("PARS_utr_gene")
write.csv(data_gene_process,file=paste0(path,"/PARS_utr_gene.csv",collapse = NULL),row.names = FALSE)
dd_gene <- rbind(dd_gene, data.frame(name="PARS_utr", gene=nrow(data_gene_process) ))
write.csv(data_utr,file=paste0(path,"/data_SNPs_PARS_utr.csv",collapse = NULL),row.names = FALSE)

# 取出syn，求出SNPs数量和mRNA数量
data_SNPs_PARS_syn <- sqldf('SELECT * FROM [data_SNPs_PARS_mRNA] where syn > 0 AND nsy == 0' )
dd_SNPs <- rbind(dd_SNPs, data.frame(name="PARS_syn", SNPs=nrow(data_SNPs_PARS_syn) ))
data_gene_process <- data_SNPs_PARS_syn["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
colnames(data_gene_process) <- c("PARS_syn_gene")
write.csv(data_gene_process,file=paste0(path,"/PARS_syn_gene.csv",collapse = NULL),row.names = FALSE)
dd_gene <- rbind(dd_gene, data.frame(name="PARS_syn", gene=nrow(data_gene_process) ))
write.csv(data_SNPs_PARS_syn,file=paste0(path,"/data_SNPs_PARS_syn.csv",collapse = NULL),row.names = FALSE)

# 取出nsy，求出SNPs数量和mRNA数量
data_SNPs_PARS_nsy <- sqldf('SELECT * FROM [data_SNPs_PARS_mRNA] where nsy > 0 AND syn == 0' )
dd_SNPs <- rbind(dd_SNPs, data.frame(name="PARS_nsy", SNPs=nrow(data_SNPs_PARS_nsy) ))
data_gene_process <- data_SNPs_PARS_nsy["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
colnames(data_gene_process) <- c("PARS_nsy_gene")
write.csv(data_gene_process,file=paste0(path,"/PARS_nsy_gene.csv",collapse = NULL),row.names = FALSE)
dd_gene <- rbind(dd_gene, data.frame(name="PARS_nsy", gene=nrow(data_gene_process) ))
write.csv(data_SNPs_PARS_nsy,file=paste0(path,"/data_SNPs_PARS_nsy.csv",collapse = NULL),row.names = FALSE)

write.csv(dd_SNPs,file=paste0(path,"/dd_SNPs.csv",collapse = NULL),row.names = FALSE)
write.csv(dd_gene,file=paste0(path,"/dd_gene.csv",collapse = NULL),row.names = FALSE)