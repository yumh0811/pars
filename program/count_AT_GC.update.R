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
path <- paste0("~/data/mrna-structure/result/", name,".update", collapse = NULL)
setwd(path)

#输入csv
file_SNPs_PARS_mRNA <- paste0("~/data/mrna-structure/result/",name,"/data_SNPs_PARS_mRNA.csv",collapse = NULL)
data_SNPs_PARS_mRNA <- read.csv(file_SNPs_PARS_mRNA,header = TRUE,sep = ",")
file_SNPs_update <- paste0("~/data/mrna-structure/vcf/1011Matrix.gvcf/",name,".mRNA.snp.update.txt",collapse = NULL)
data_SNPs_update <- read.csv(file_SNPs_update,header = TRUE,sep = ",")

data_SNPs_PARS_mRNA <- merge(data_SNPs_PARS_mRNA,data_SNPs_update,by="Location")
data_SNPs_PARS_mRNA_gene <- unique(data_SNPs_PARS_mRNA["Gene"])
write.csv(data_SNPs_PARS_mRNA, file="data_SNPs_PARS_mRNA.update.csv", row.names = FALSE)
write.csv(data_SNPs_PARS_mRNA_gene, file="mRNA.gene.list.update.csv", row.names = FALSE)

group = c("mRNA")
for (g in group){
  t <- get(paste0('data_SNPs_PARS_',g,collapse = NULL))
  for(i in 1:max(t$Freq)){
    n <- assign(paste0('data_SNPs_PARS_',g,'_',i,collapse = NULL),subset(t, Freq == i))
    if(i==1){
      dd_SNPs_Freq <- data.frame(name = i, SNPs = c(nrow(n)))
    }else{
      dd_SNPs_Freq <- rbind(dd_SNPs_Freq, data.frame(name = i, SNPs = nrow(n)))
    }
    data_gene_process <- n["Gene"]
    data_gene_process <- unique(data_gene_process,fromLast=TRUE)
    if(i==1){
      dd_gene_Freq <- data.frame(name = i, gene=c(nrow(data_gene_process)))
    }else{
      dd_gene_Freq <- rbind(dd_gene_Freq, data.frame(name = i, gene = nrow(data_gene_process) ))
    }
    # 统计每个Freq的stem-loop的SNPs和gene的情况 
    'stem'
    
    data_stem <- subset(n,(structure == "stem"))
    if(i==1){
      dd_SNPs_Freq_stem <- data.frame(name = i, SNPs = c(nrow(data_stem)))
    }else{
      dd_SNPs_Freq_stem <- rbind(dd_SNPs_Freq_stem , data.frame(name = i, SNPs = nrow(data_stem)))
    }
    
    data_stem_AT_GC <- sqldf('SELECT * FROM [data_stem] where mutant_to == "A->G" OR mutant_to == "A->C" OR mutant_to == "T->G" OR mutant_to == "T->C"' )
    if(i==1){
      dd_SNPs_Freq_stem_AT_GC <- data.frame(name = i, SNPs = c(nrow(data_stem_AT_GC)))
    }else{
      dd_SNPs_Freq_stem_AT_GC <- rbind(dd_SNPs_Freq_stem_AT_GC , data.frame(name = i, SNPs = nrow(data_stem_AT_GC)))
    }
    
    data_stem_GC_AT <- sqldf('SELECT * FROM [data_stem] where mutant_to == "G->A" OR mutant_to == "C->A" OR mutant_to == "G->T" OR mutant_to == "C->T"') 
    if(i==1){
      dd_SNPs_Freq_stem_GC_AT <- data.frame(name = i, SNPs = c(nrow(data_stem_GC_AT)))
    }else{
      dd_SNPs_Freq_stem_GC_AT <- rbind(dd_SNPs_Freq_stem_GC_AT , data.frame(name = i, SNPs = nrow(data_stem_GC_AT)))
    }
    
    'loop'
    
    data_loop <- subset(n,(structure == "loop"))
    if(i==1){
      dd_SNPs_Freq_loop <- data.frame(name = i, SNPs = c(nrow(data_loop)))
    }else{
      dd_SNPs_Freq_loop <- rbind(dd_SNPs_Freq_loop , data.frame(name = i, SNPs = nrow(data_loop)))
    }
    
    data_loop_AT_GC <- sqldf('SELECT * FROM [data_loop] where mutant_to == "A->G" OR mutant_to == "A->C" OR mutant_to == "T->G" OR mutant_to == "T->C"' ) 
    if(i==1){
      dd_SNPs_Freq_loop_AT_GC <- data.frame(name = i, SNPs = c(nrow(data_loop_AT_GC)))
    }else{
      dd_SNPs_Freq_loop_AT_GC <- rbind(dd_SNPs_Freq_loop_AT_GC , data.frame(name = i, SNPs = nrow(data_loop_AT_GC)))
    }
    
    data_loop_GC_AT <- sqldf('SELECT * FROM [data_loop] where mutant_to == "G->A" OR mutant_to == "C->A" OR mutant_to == "G->T" OR mutant_to == "C->T"')    
    if(i==1){
      dd_SNPs_Freq_loop_GC_AT <- data.frame(name = i, SNPs = c(nrow(data_loop_GC_AT)))
    }else{
      dd_SNPs_Freq_loop_GC_AT <- rbind(dd_SNPs_Freq_loop_GC_AT , data.frame(name = i, SNPs = nrow(data_loop_GC_AT)))
    }
    
    assign(paste0('data_stat_',i),data.frame(structure=c("stem","loop"),AT_GC=c(nrow(data_stem_AT_GC),nrow(data_loop_AT_GC)),GC_AT=c(nrow(data_stem_GC_AT),nrow(data_loop_GC_AT))))
  }
  # 合并多个数据框
  data <- lapply(paste0('data_stat_',1:max(t$Freq)), function(data_stat_) eval(as.name(data_stat_)))
  data_stat <- do.call("rbind", data)

  write.csv(data_stat, file=paste0(path,'/freq_each/PARS_',g,'_stat.csv'), row.names = FALSE)
  
  write.csv(dd_SNPs_Freq, file=paste0(path,'/freq_each/PARS_',g,'_stat_SNPs.csv'), row.names = FALSE)
  write.csv(dd_gene_Freq, file=paste0(path,'/freq_each/PARS_',g,'_stat_gene.csv'), row.names = FALSE)
  write.csv(dd_SNPs_Freq_stem, file=paste0(path,'/freq_each/PARS_',g,'_stat_stem.csv'), row.names = FALSE)
  write.csv(dd_SNPs_Freq_loop, file=paste0(path,'/freq_each/PARS_',g,'_stat_loop.csv'), row.names = FALSE)
  write.csv(dd_SNPs_Freq_stem_AT_GC, file=paste0(path,'/freq_each/PARS_',g,'_stat_stem_AT_GC.csv'), row.names = FALSE)
  write.csv(dd_SNPs_Freq_loop_AT_GC, file=paste0(path,'/freq_each/PARS_',g,'_stat_loop_AT_GC.csv'), row.names = FALSE)
  write.csv(dd_SNPs_Freq_stem_GC_AT, file=paste0(path,'/freq_each/PARS_',g,'_stat_stem_GC_AT.csv'), row.names = FALSE)
  write.csv(dd_SNPs_Freq_loop_GC_AT, file=paste0(path,'/freq_each/PARS_',g,'_stat_loop_GC_AT.csv'), row.names = FALSE)
}



# 把Freq分组（先判断Freq是否>=10)
group = c("mRNA")
for (g in group){
  t <- get(paste0('data_SNPs_PARS_',g,collapse = NULL))
  if(max(t$Freq)>=10){
    for(i in 1:10){
      # 统计每个Freq的总SNPs和总gene的情况
      n <- assign(paste0('data_SNPs_PARS_',g,'_',i,collapse = NULL),subset(t,Freq <= max(Freq)*(i/10) & Freq > max(Freq)*((i-1)/10)))
      if(i==1){
        dd_SNPs_Freq <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(n)))
      }else{
        dd_SNPs_Freq <- rbind(dd_SNPs_Freq, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(n)))
      }
      data_gene_process <- n["Gene"]
      data_gene_process <- unique(data_gene_process,fromLast=TRUE)
      if(i==1){
        dd_gene_Freq <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), gene = c(nrow(data_gene_process)))
      }else{
        dd_gene_Freq <- rbind(dd_gene_Freq, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), gene = nrow(data_gene_process) ))
      }  
      # 统计每个Freq的stem-loop的SNPs和gene的情况 
      'stem'
      
      data_stem <- subset(n,(structure == "stem"))
      if(i==1){
        dd_SNPs_Freq_stem <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_stem)))
      }else{
        dd_SNPs_Freq_stem <- rbind(dd_SNPs_Freq_stem, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_stem)))
      }
      
      data_stem_AT_GC <- sqldf('SELECT * FROM [data_stem] where mutant_to == "A->G" OR mutant_to == "A->C" OR mutant_to == "T->G" OR mutant_to == "T->C"' )
      if(i==1){
        dd_SNPs_Freq_stem_AT_GC <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_stem_AT_GC)))
      }else{
        dd_SNPs_Freq_stem_AT_GC <- rbind(dd_SNPs_Freq_stem_AT_GC, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_stem_AT_GC)))
      }
      
      data_stem_GC_AT <- sqldf('SELECT * FROM [data_stem] where mutant_to == "G->A" OR mutant_to == "C->A" OR mutant_to == "G->T" OR mutant_to == "C->T"')  
      if(i==1){
        dd_SNPs_Freq_stem_GC_AT <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_stem_GC_AT)))
      }else{
        dd_SNPs_Freq_stem_GC_AT <- rbind(dd_SNPs_Freq_stem_GC_AT, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_stem_GC_AT)))
      }
      
      'loop'
      
      data_loop <- subset(n,(structure == "loop"))
      if(i==1){
        dd_SNPs_Freq_loop <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_loop)))
      }else{
        dd_SNPs_Freq_loop <- rbind(dd_SNPs_Freq_loop, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_loop)))
      }
      
      data_loop_AT_GC <- sqldf('SELECT * FROM [data_loop] where mutant_to == "A->G" OR mutant_to == "A->C" OR mutant_to == "T->G" OR mutant_to == "T->C"' )
      if(i==1){
        dd_SNPs_Freq_loop_AT_GC <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_loop_AT_GC)))
      }else{
        dd_SNPs_Freq_loop_AT_GC <- rbind(dd_SNPs_Freq_loop_AT_GC, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_loop_AT_GC)))
      }
      
      data_loop_GC_AT <- sqldf('SELECT * FROM [data_loop] where mutant_to == "G->A" OR mutant_to == "C->A" OR mutant_to == "G->T" OR mutant_to == "C->T"')  
      if(i==1){
        dd_SNPs_Freq_loop_GC_AT <- data.frame(name = paste0("0-",i,"0%",collapse = NULL), SNPs = c(nrow(data_loop_GC_AT)))
      }else{
        dd_SNPs_Freq_loop_GC_AT <- rbind(dd_SNPs_Freq_loop_GC_AT, data.frame(name = paste0(i-1,"0","-",i,"0%",collapse = NULL), SNPs = nrow(data_loop_GC_AT)))
      }
      
      assign(paste0('data_stat_',i),data.frame(structure=c("stem","loop"),AT_GC=c(nrow(data_stem_AT_GC),nrow(data_loop_AT_GC)),GC_AT=c(nrow(data_stem_GC_AT),nrow(data_loop_GC_AT))))
      
    }
    # 合并多个数据框
    data <- lapply(paste0('data_stat_',1:10), function(data_stat_) eval(as.name(data_stat_)))
    data_stat <- do.call("rbind", data)
    
    write.csv(data_stat, file=paste0(path,'/freq_10/PARS_',g,'_stat_freq_10.csv'), row.names = FALSE)
    
    write.csv(dd_SNPs_Freq, file=paste0(path,'/freq_10/PARS_',g,'_stat_SNPs_freq_10.csv'), row.names = FALSE)
    write.csv(dd_gene_Freq, file=paste0(path,'/freq_10/PARS_',g,'_stat_gene_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_Freq_stem, file=paste0(path,'/freq_10/PARS_',g,'_stat_stem_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_Freq_loop, file=paste0(path,'/freq_10/PARS_',g,'_stat_loop_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_Freq_stem_AT_GC, file=paste0(path,'/freq_10/PARS_',g,'_stat_stem_AT_GC_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_Freq_loop_AT_GC, file=paste0(path,'/freq_10/PARS_',g,'_stat_loop_AT_GC_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_Freq_stem_GC_AT, file=paste0(path,'/freq_10/PARS_',g,'_stat_stem_GC_AT_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_Freq_loop_GC_AT, file=paste0(path,'/freq_10/PARS_',g,'_stat_loop_GC_AT_freq_10.csv'), row.names = FALSE)
  }

}
