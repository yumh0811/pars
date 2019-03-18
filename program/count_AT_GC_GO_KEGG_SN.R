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
#name <- "Scer_n128_Spar.update"
name <- opt$name
path <- paste0("~/data/mrna-structure/result/", name, collapse = NULL)
setwd(path)

#输入csv
file_SNPs_PARS_mRNA <- paste0(path,'/data_SNPs_PARS_mRNA.tsv',collapse = NULL)
data_SNPs_PARS_mRNA <- read.csv(file_SNPs_PARS_mRNA,header = TRUE,sep = "\t")
dd_SNP <- data.frame(name = "mRNA",SNP = c(nrow(data_SNPs_PARS_mRNA)))
data_gene_process <- data_SNPs_PARS_mRNA["gene"]
dd_gene<- data.frame(name = "mRNA", gene = c(nrow(data_gene_process)))

# 取出syn，求出SNPs数量和mRNA数量
data_SNPs_PARS_syn <- sqldf('SELECT * FROM [data_SNPs_PARS_mRNA] where Consequence == "stop_retained_variant" OR Consequence == "synonymous_variant"')
dd_SNP <- rbind(dd_SNP, data.frame(name="SYN", SNP=nrow(data_SNPs_PARS_syn)))
dd_gene <- rbind(dd_gene, data.frame(name="SYN", gene=nrow(unique(data_SNPs_PARS_syn["gene"],fromLast=TRUE))))

# 取出nsy，求出SNPs数量和mRNA数量
data_SNPs_PARS_nsy <- sqldf('SELECT * FROM [data_SNPs_PARS_mRNA] where Consequence == "missense_variant" OR Consequence == "start_lost" OR Consequence == "stop_gained" OR Consequence == "stop_lost"')
dd_SNP <- rbind(dd_SNP, data.frame(name="NSY", SNP=nrow(data_SNPs_PARS_nsy)))
dd_gene <- rbind(dd_gene, data.frame(name="NSY", gene=nrow(unique(data_SNPs_PARS_nsy["gene"],fromLast=TRUE))))

#go_kegg
file_gene_go_kegg <- "Scer_n128_Spar_go_kegg.csv" 
data_gene_go_kegg <- read.csv(file_gene_go_kegg,header = FALSE ,sep = ",")
colnames(data_gene_go_kegg) <- c(1:41)
for (go_kegg in 1:41){
  gene <- data.frame(data_gene_go_kegg[3:nrow(data_gene_go_kegg),go_kegg],row.names = NULL)
  colnames(gene) <- "gene"
  snp <- assign(paste0("data_SNPs_PARS_syn_go_kegg_",go_kegg),merge(data_SNPs_PARS_syn,gene,by="gene"))
  dd_SNP <- rbind(dd_SNP,data.frame(name=paste0("go_kegg_",go_kegg,"_",data_gene_go_kegg[2,go_kegg]),SNP = c(nrow(snp))))
  data_gene_process <- snp["gene"]
  data_gene_process <- unique(data_gene_process,fromLast=TRUE)
  dd_gene <- rbind(dd_gene, data.frame(name = data_gene_go_kegg[2,go_kegg], gene = c(nrow(data_gene_process))))
  write.csv(data_gene_process , file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_gene_',go_kegg,"_",data_gene_go_kegg[2,go_kegg],".csv"), row.names = FALSE)
  write.csv(snp , file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_SNP_',go_kegg,"_",data_gene_go_kegg[2,go_kegg],".csv"), row.names = FALSE)
  
    if(max(snp$freq)>=10){
      for(i in 1:10){
        # 统计每个freq的总SNPs和总gene的情况
        n <- assign(paste0('snp_',i,collapse = NULL),subset(snp,freq <= max(freq)*(i/10) & freq > max(freq)*((i-1)/10)))
        if(i==1){
          dd_SNPs_freq <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(n)))
        }else{
          dd_SNPs_freq <- rbind(dd_SNPs_freq, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(n)))
        }
        data_gene_process <- n["gene"]
        data_gene_process <- unique(data_gene_process,fromLast=TRUE)
        if(i==1){
          dd_gene_freq <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), gene = c(nrow(data_gene_process)))
        }else{
          dd_gene_freq <- rbind(dd_gene_freq, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), gene = nrow(data_gene_process) ))
        }  
        # 统计每个freq的stem-loop的SNPs和gene的情况 
        'stem'
        
        data_stem <- subset(n,(structure == "stem"))
        if(i==1){
          dd_SNPs_freq_stem <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_stem)))
        }else{
          dd_SNPs_freq_stem <- rbind(dd_SNPs_freq_stem, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_stem)))
        }
        
        data_stem_AT_GC <- sqldf('SELECT * FROM [data_stem] where mutant_to == "A->G" OR mutant_to == "A->C" OR mutant_to == "T->G" OR mutant_to == "T->C"' )
        if(i==1){
          dd_SNPs_freq_stem_AT_GC <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_stem_AT_GC)))
        }else{
          dd_SNPs_freq_stem_AT_GC <- rbind(dd_SNPs_freq_stem_AT_GC, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_stem_AT_GC)))
        }
        
        data_stem_GC_AT <- sqldf('SELECT * FROM [data_stem] where mutant_to == "G->A" OR mutant_to == "C->A" OR mutant_to == "G->T" OR mutant_to == "C->T"')  
        if(i==1){
          dd_SNPs_freq_stem_GC_AT <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_stem_GC_AT)))
        }else{
          dd_SNPs_freq_stem_GC_AT <- rbind(dd_SNPs_freq_stem_GC_AT, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_stem_GC_AT)))
        }
        
        'loop'
        
        data_loop <- subset(n,(structure == "loop"))
        if(i==1){
          dd_SNPs_freq_loop <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_loop)))
        }else{
          dd_SNPs_freq_loop <- rbind(dd_SNPs_freq_loop, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_loop)))
        }
        
        data_loop_AT_GC <- sqldf('SELECT * FROM [data_loop] where mutant_to == "A->G" OR mutant_to == "A->C" OR mutant_to == "T->G" OR mutant_to == "T->C"' )
        if(i==1){
          dd_SNPs_freq_loop_AT_GC <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_loop_AT_GC)))
        }else{
          dd_SNPs_freq_loop_AT_GC <- rbind(dd_SNPs_freq_loop_AT_GC, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_loop_AT_GC)))
        }
        
        data_loop_GC_AT <- sqldf('SELECT * FROM [data_loop] where mutant_to == "G->A" OR mutant_to == "C->A" OR mutant_to == "G->T" OR mutant_to == "C->T"')  
        if(i==1){
          dd_SNPs_freq_loop_GC_AT <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_loop_GC_AT)))
        }else{
          dd_SNPs_freq_loop_GC_AT <- rbind(dd_SNPs_freq_loop_GC_AT, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_loop_GC_AT)))
        }
        
        assign(paste0('data_stat_',i),data.frame(structure=c("stem","loop"),AT_GC=c(nrow(data_stem_AT_GC),nrow(data_loop_AT_GC)),GC_AT=c(nrow(data_stem_GC_AT),nrow(data_loop_GC_AT))))
        
      }
      # 合并多个数据框
      data <- lapply(paste0('data_stat_',1:10), function(data_stat_) eval(as.name(data_stat_)))
      data_stat <- do.call("rbind", data)
      
      write.csv(data_stat, file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_stat_',go_kegg,"_",data_gene_go_kegg[2,go_kegg],'_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq, file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_',go_kegg,'_stat_SNPs_freq_10.csv'), row.names = FALSE)
      write.csv(dd_gene_freq, file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_',go_kegg,'_stat_gene_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_stem, file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_',go_kegg,'_stat_stem_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_loop, file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_',go_kegg,'_stat_loop_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_stem_AT_GC, file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_',go_kegg,'_stat_stem_AT_GC_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_loop_AT_GC, file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_',go_kegg,'_stat_loop_AT_GC_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_stem_GC_AT, file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_',go_kegg,'_stat_stem_GC_AT_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_loop_GC_AT, file=paste0(path,'/freq_10/go_kegg/syn/go_kegg_',go_kegg,'_stat_loop_GC_AT_freq_10.csv'), row.names = FALSE)
    }
    
}
  
write.csv(dd_gene, file=paste0(path,'/freq_10/go_kegg/syn/dd_gene.csv'), row.names = FALSE)
write.csv(dd_SNP, file=paste0(path,'/freq_10/go_kegg/syn/dd_SNP.csv'), row.names = FALSE)

rm (dd_gene,dd_SNP)

#输入csv
file_SNPs_PARS_mRNA <- paste0(path,'/data_SNPs_PARS_mRNA.tsv',collapse = NULL)
data_SNPs_PARS_mRNA <- read.csv(file_SNPs_PARS_mRNA,header = TRUE,sep = "\t")
dd_SNP <- data.frame(name = "mRNA",SNP = c(nrow(data_SNPs_PARS_mRNA)))
data_gene_process <- data_SNPs_PARS_mRNA["gene"]
dd_gene<- data.frame(name = "mRNA", gene = c(nrow(data_gene_process)))

# 取出syn，求出SNPs数量和mRNA数量
data_SNPs_PARS_syn <- sqldf('SELECT * FROM [data_SNPs_PARS_mRNA] where Consequence == "stop_retained_variant" OR Consequence == "synonymous_variant"')
dd_SNP <- rbind(dd_SNP, data.frame(name="SYN", SNP=nrow(data_SNPs_PARS_syn)))
dd_gene <- rbind(dd_gene, data.frame(name="SYN", gene=nrow(unique(data_SNPs_PARS_syn["gene"],fromLast=TRUE))))

# 取出nsy，求出SNPs数量和mRNA数量
data_SNPs_PARS_nsy <- sqldf('SELECT * FROM [data_SNPs_PARS_mRNA] where Consequence == "missense_variant" OR Consequence == "start_lost" OR Consequence == "stop_gained" OR Consequence == "stop_lost"')
dd_SNP <- rbind(dd_SNP, data.frame(name="NSY", SNP=nrow(data_SNPs_PARS_nsy)))
dd_gene <- rbind(dd_gene, data.frame(name="NSY", gene=nrow(unique(data_SNPs_PARS_nsy["gene"],fromLast=TRUE))))

#go_kegg
file_gene_go_kegg <- "Scer_n128_Spar_go_kegg.csv" 
data_gene_go_kegg <- read.csv(file_gene_go_kegg,header = FALSE ,sep = ",")
colnames(data_gene_go_kegg) <- c(1:41)
for (go_kegg in 1:41){
  gene <- data.frame(data_gene_go_kegg[3:nrow(data_gene_go_kegg),go_kegg],row.names = NULL)
  colnames(gene) <- "gene"
  snp <- assign(paste0("data_SNPs_PARS_nsy_go_kegg_",go_kegg),merge(data_SNPs_PARS_nsy,gene,by="gene"))
  dd_SNP <- rbind(dd_SNP,data.frame(name=paste0("go_kegg_",go_kegg,"_",data_gene_go_kegg[2,go_kegg]),SNP = c(nrow(snp))))
  data_gene_process <- snp["gene"]
  data_gene_process <- unique(data_gene_process,fromLast=TRUE)
  dd_gene <- rbind(dd_gene, data.frame(name = data_gene_go_kegg[2,go_kegg], gene = c(nrow(data_gene_process))))
  write.csv(data_gene_process , file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_gene_',go_kegg,"_",data_gene_go_kegg[2,go_kegg],".csv"), row.names = FALSE)
  write.csv(snp , file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_SNP_',go_kegg,"_",data_gene_go_kegg[2,go_kegg],".csv"), row.names = FALSE)
  
    if(max(snp$freq)>=10){
      for(i in 1:10){
        # 统计每个freq的总SNPs和总gene的情况
        n <- assign(paste0('snp_',i,collapse = NULL),subset(snp,freq <= max(freq)*(i/10) & freq > max(freq)*((i-1)/10)))
        if(i==1){
          dd_SNPs_freq <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(n)))
        }else{
          dd_SNPs_freq <- rbind(dd_SNPs_freq, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(n)))
        }
        data_gene_process <- n["gene"]
        data_gene_process <- unique(data_gene_process,fromLast=TRUE)
        if(i==1){
          dd_gene_freq <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), gene = c(nrow(data_gene_process)))
        }else{
          dd_gene_freq <- rbind(dd_gene_freq, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), gene = nrow(data_gene_process) ))
        }  
        # 统计每个freq的stem-loop的SNPs和gene的情况 
        'stem'
        
        data_stem <- subset(n,(structure == "stem"))
        if(i==1){
          dd_SNPs_freq_stem <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_stem)))
        }else{
          dd_SNPs_freq_stem <- rbind(dd_SNPs_freq_stem, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_stem)))
        }
        
        data_stem_AT_GC <- sqldf('SELECT * FROM [data_stem] where mutant_to == "A->G" OR mutant_to == "A->C" OR mutant_to == "T->G" OR mutant_to == "T->C"' )
        if(i==1){
          dd_SNPs_freq_stem_AT_GC <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_stem_AT_GC)))
        }else{
          dd_SNPs_freq_stem_AT_GC <- rbind(dd_SNPs_freq_stem_AT_GC, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_stem_AT_GC)))
        }
        
        data_stem_GC_AT <- sqldf('SELECT * FROM [data_stem] where mutant_to == "G->A" OR mutant_to == "C->A" OR mutant_to == "G->T" OR mutant_to == "C->T"')  
        if(i==1){
          dd_SNPs_freq_stem_GC_AT <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_stem_GC_AT)))
        }else{
          dd_SNPs_freq_stem_GC_AT <- rbind(dd_SNPs_freq_stem_GC_AT, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_stem_GC_AT)))
        }
        
        'loop'
        
        data_loop <- subset(n,(structure == "loop"))
        if(i==1){
          dd_SNPs_freq_loop <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_loop)))
        }else{
          dd_SNPs_freq_loop <- rbind(dd_SNPs_freq_loop, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_loop)))
        }
        
        data_loop_AT_GC <- sqldf('SELECT * FROM [data_loop] where mutant_to == "A->G" OR mutant_to == "A->C" OR mutant_to == "T->G" OR mutant_to == "T->C"' )
        if(i==1){
          dd_SNPs_freq_loop_AT_GC <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_loop_AT_GC)))
        }else{
          dd_SNPs_freq_loop_AT_GC <- rbind(dd_SNPs_freq_loop_AT_GC, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_loop_AT_GC)))
        }
        
        data_loop_GC_AT <- sqldf('SELECT * FROM [data_loop] where mutant_to == "G->A" OR mutant_to == "C->A" OR mutant_to == "G->T" OR mutant_to == "C->T"')  
        if(i==1){
          dd_SNPs_freq_loop_GC_AT <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_loop_GC_AT)))
        }else{
          dd_SNPs_freq_loop_GC_AT <- rbind(dd_SNPs_freq_loop_GC_AT, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_loop_GC_AT)))
        }
        
        assign(paste0('data_stat_',i),data.frame(structure=c("stem","loop"),AT_GC=c(nrow(data_stem_AT_GC),nrow(data_loop_AT_GC)),GC_AT=c(nrow(data_stem_GC_AT),nrow(data_loop_GC_AT))))
        
      }
      # 合并多个数据框
      data <- lapply(paste0('data_stat_',1:10), function(data_stat_) eval(as.name(data_stat_)))
      data_stat <- do.call("rbind", data)
      
      write.csv(data_stat, file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_stat_',go_kegg,"_",data_gene_go_kegg[2,go_kegg],'_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq, file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_',go_kegg,'_stat_SNPs_freq_10.csv'), row.names = FALSE)
      write.csv(dd_gene_freq, file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_',go_kegg,'_stat_gene_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_stem, file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_',go_kegg,'_stat_stem_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_loop, file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_',go_kegg,'_stat_loop_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_stem_AT_GC, file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_',go_kegg,'_stat_stem_AT_GC_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_loop_AT_GC, file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_',go_kegg,'_stat_loop_AT_GC_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_stem_GC_AT, file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_',go_kegg,'_stat_stem_GC_AT_freq_10.csv'), row.names = FALSE)
      write.csv(dd_SNPs_freq_loop_GC_AT, file=paste0(path,'/freq_10/go_kegg/nsy/go_kegg_',go_kegg,'_stat_loop_GC_AT_freq_10.csv'), row.names = FALSE)
    }
    
}
  
write.csv(dd_gene, file=paste0(path,'/freq_10/go_kegg/nsy/dd_gene.csv'), row.names = FALSE)
write.csv(dd_SNP, file=paste0(path,'/freq_10/go_kegg/nsy/dd_SNP.csv'), row.names = FALSE)
