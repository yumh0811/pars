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

#输入csv
file_SNPs_PARS_cds <- paste0(path,'/data_SNPs_PARS_cds.update_codon.csv',collapse = NULL)
data_SNPs_PARS_cds <- read.csv(file_SNPs_PARS_cds,header = TRUE,sep = ",")

#codon
data_SNPs_PARS_4D <- sqldf('SELECT * FROM [data_SNPs_PARS_cds] where snp_codons ==  
                           "CCT->CCC" OR snp_codons == "CCT->CCA" OR snp_codons == "CCT->CCG" OR snp_codons == "CCC->CCT" OR snp_codons == "CCC->CCA" OR snp_codons == "CCC->CCG" OR snp_codons == "CCA->CCC" OR snp_codons == "CCA->CCT" OR snp_codons == "CCA->CCG" OR snp_codons == "CCG->CCC" OR snp_codons == "CCG->CCT" OR snp_codons == "CCG->CCA" OR snp_codons == 
                           "GTT->GTC" OR snp_codons == "GTT->GTA" OR snp_codons == "GTT->GTG" OR snp_codons == "GTC->GTT" OR snp_codons == "GTC->GTA" OR snp_codons == "GTC->GTG" OR snp_codons == "GTA->GTC" OR snp_codons == "GTA->GTT" OR snp_codons == "GTA->GTG" OR snp_codons == "GTG->GTC" OR snp_codons == "GTG->GTT" OR snp_codons == "GTG->GTA" OR snp_codons == 
                           "ACT->ACC" OR snp_codons == "ACT->ACA" OR snp_codons == "ACT->ACG" OR snp_codons == "ACC->ACT" OR snp_codons == "ACC->ACA" OR snp_codons == "ACC->ACG" OR snp_codons == "ACA->ACC" OR snp_codons == "ACA->ACT" OR snp_codons == "ACA->ACG" OR snp_codons == "ACG->ACC" OR snp_codons == "ACG->ACT" OR snp_codons == "ACG->ACA" OR snp_codons == 
                           "GCT->GCC" OR snp_codons == "GCT->GCA" OR snp_codons == "GCT->GCG" OR snp_codons == "GCC->GCT" OR snp_codons == "GCC->GCA" OR snp_codons == "GCC->GCG" OR snp_codons == "GCA->GCC" OR snp_codons == "GCA->GCT" OR snp_codons == "GCA->GCG" OR snp_codons == "GCG->GCC" OR snp_codons == "GCG->GCT" OR snp_codons == "GCG->GCA" OR snp_codons == 
                           "GGT->GGC" OR snp_codons == "GGT->GGA" OR snp_codons == "GGT->GGG" OR snp_codons == "GGC->GGT" OR snp_codons == "GGC->GGA" OR snp_codons == "GGC->GGG" OR snp_codons == "GGA->GGC" OR snp_codons == "GGA->GGT" OR snp_codons == "GGA->GGG" OR snp_codons == "GGG->GGC" OR snp_codons == "GGG->GGT" OR snp_codons == "GGG->GGA"
                           ')
write.csv(data_SNPs_PARS_4D,file=paste0(path,"/data_SNPs_PARS_4D.csv",collapse = NULL),row.names = FALSE)

#tRNA
data_SNPs_PARS_tRNA <- sqldf('SELECT * FROM [data_SNPs_PARS_cds] where 
                             snp_codons == "CTA->CTG" OR snp_codons == "CTG->CTA" OR 
                             snp_codons == "GTT->GTC" OR snp_codons == "GTT->GTA" OR snp_codons == "GTC->GTT" OR snp_codons == "GTC->GTA" OR snp_codons == "GTA->GTT" OR snp_codons == "GTA->GTC" OR
                             snp_codons == "GGT->GGC" OR snp_codons == "GGC->GGT" OR
                             snp_codons == "GCT->GCC" OR snp_codons == "GCT->GCA" OR snp_codons == "GCC->GCT" OR snp_codons == "GCC->GCA" OR snp_codons == "GCA->GCT" OR snp_codons == "GCA->GCC" OR
                             snp_codons == "AGA->AGG" OR snp_codons == "AGG->AGA" OR
                             snp_codons == "CGT->CGC" OR snp_codons == "CGT->CGA" OR snp_codons == "CGC->CGT" OR snp_codons == "CGC->CGA" OR snp_codons == "CGA->CGT" OR snp_codons == "CGA->CGC" OR
                             snp_codons == "AAA->AAG" OR snp_codons == "AAG->AAA" OR
                             snp_codons == "GAA->GAG" OR snp_codons == "GAG->GAA" OR
                             snp_codons == "GAT->GAC" OR snp_codons == "GAC->GAT" OR
                             snp_codons == "ACT->ACC" OR snp_codons == "ACT->ACA" OR snp_codons == "ACC->ACT" OR snp_codons == "ACC->ACA" OR snp_codons == "ACA->ACT" OR snp_codons == "ACA->ACC" OR
                             snp_codons == "TAT->TAC" OR snp_codons == "TAC->TAT" OR
                             snp_codons == "TCT->TCC" OR snp_codons == "TCT->TCA" OR snp_codons == "TCC->TCT" OR snp_codons == "TCC->TCA" OR snp_codons == "TCA->TCT" OR snp_codons == "TCA->TCC" OR
                             snp_codons == "TCA->TCG" OR snp_codons == "TCG->TCA" OR
                             snp_codons == "CAT->CAC" OR snp_codons == "CAC->CAT" OR
                             snp_codons == "TTT->TTC" OR snp_codons == "TTC->TTT" OR
                             snp_codons == "TGT->TGC" OR snp_codons == "TGC->TGT" ')
write.csv(data_SNPs_PARS_tRNA,file=paste0(path,"/data_SNPs_PARS_tRNA.csv",collapse = NULL),row.names = FALSE)

group = c("4D","tRNA")
for (g in group){
  t <- get(paste0('data_SNPs_PARS_',g,collapse = NULL))
  for(i in 1:max(t$freq)){
    n <- assign(paste0('data_SNPs_PARS_',g,'_',i,collapse = NULL),subset(t, freq == i))
    if(i==1){
      dd_SNPs_freq <- data.frame(name = i, SNPs = c(nrow(n)))
    }else{
      dd_SNPs_freq <- rbind(dd_SNPs_freq, data.frame(name = i, SNPs = nrow(n)))
    }
    data_gene_process <- n["gene"]
    data_gene_process <- unique(data_gene_process,fromLast=TRUE)
    if(i==1){
      dd_gene_freq <- data.frame(name = i, gene=c(nrow(data_gene_process)))
    }else{
      dd_gene_freq <- rbind(dd_gene_freq, data.frame(name = i, gene = nrow(data_gene_process) ))
    }
    # 统计每个freq的stem-loop的SNPs和gene的情况 
    'stem'
    
    data_stem <- subset(n,(structure == "stem"))
    if(i==1){
      dd_SNPs_freq_stem <- data.frame(name = i, SNPs = c(nrow(data_stem)))
    }else{
      dd_SNPs_freq_stem <- rbind(dd_SNPs_freq_stem , data.frame(name = i, SNPs = nrow(data_stem)))
    }
    
    data_stem_AT_GC <- sqldf('SELECT * FROM [data_stem] where mutant_to == "A->G" OR mutant_to == "A->C" OR mutant_to == "T->G" OR mutant_to == "T->C"' )
    if(i==1){
      dd_SNPs_freq_stem_AT_GC <- data.frame(name = i, SNPs = c(nrow(data_stem_AT_GC)))
    }else{
      dd_SNPs_freq_stem_AT_GC <- rbind(dd_SNPs_freq_stem_AT_GC , data.frame(name = i, SNPs = nrow(data_stem_AT_GC)))
    }
    
    data_stem_GC_AT <- sqldf('SELECT * FROM [data_stem] where mutant_to == "G->A" OR mutant_to == "C->A" OR mutant_to == "G->T" OR mutant_to == "C->T"') 
    if(i==1){
      dd_SNPs_freq_stem_GC_AT <- data.frame(name = i, SNPs = c(nrow(data_stem_GC_AT)))
    }else{
      dd_SNPs_freq_stem_GC_AT <- rbind(dd_SNPs_freq_stem_GC_AT , data.frame(name = i, SNPs = nrow(data_stem_GC_AT)))
    }
    
    'loop'
    
    data_loop <- subset(n,(structure == "loop"))
    if(i==1){
      dd_SNPs_freq_loop <- data.frame(name = i, SNPs = c(nrow(data_loop)))
    }else{
      dd_SNPs_freq_loop <- rbind(dd_SNPs_freq_loop , data.frame(name = i, SNPs = nrow(data_loop)))
    }
    
    data_loop_AT_GC <- sqldf('SELECT * FROM [data_loop] where mutant_to == "A->G" OR mutant_to == "A->C" OR mutant_to == "T->G" OR mutant_to == "T->C"' ) 
    if(i==1){
      dd_SNPs_freq_loop_AT_GC <- data.frame(name = i, SNPs = c(nrow(data_loop_AT_GC)))
    }else{
      dd_SNPs_freq_loop_AT_GC <- rbind(dd_SNPs_freq_loop_AT_GC , data.frame(name = i, SNPs = nrow(data_loop_AT_GC)))
    }
    
    data_loop_GC_AT <- sqldf('SELECT * FROM [data_loop] where mutant_to == "G->A" OR mutant_to == "C->A" OR mutant_to == "G->T" OR mutant_to == "C->T"')    
    if(i==1){
      dd_SNPs_freq_loop_GC_AT <- data.frame(name = i, SNPs = c(nrow(data_loop_GC_AT)))
    }else{
      dd_SNPs_freq_loop_GC_AT <- rbind(dd_SNPs_freq_loop_GC_AT , data.frame(name = i, SNPs = nrow(data_loop_GC_AT)))
    }
    
    assign(paste0('data_stat_',i),data.frame(structure=c("stem","loop"),AT_GC=c(nrow(data_stem_AT_GC),nrow(data_loop_AT_GC)),GC_AT=c(nrow(data_stem_GC_AT),nrow(data_loop_GC_AT))))
  }
  # 合并多个数据框
  data <- lapply(paste0('data_stat_',1:max(t$freq)), function(data_stat_) eval(as.name(data_stat_)))
  data_stat <- do.call("rbind", data)
  
  write.csv(data_stat, file=paste0(path,'/freq_each/PARS_',g,'_stat.csv'), row.names = FALSE)
  
  write.csv(dd_SNPs_freq, file=paste0(path,'/freq_each/PARS_',g,'_stat_SNPs.csv'), row.names = FALSE)
  write.csv(dd_gene_freq, file=paste0(path,'/freq_each/PARS_',g,'_stat_gene.csv'), row.names = FALSE)
  write.csv(dd_SNPs_freq_stem, file=paste0(path,'/freq_each/PARS_',g,'_stat_stem.csv'), row.names = FALSE)
  write.csv(dd_SNPs_freq_loop, file=paste0(path,'/freq_each/PARS_',g,'_stat_loop.csv'), row.names = FALSE)
  write.csv(dd_SNPs_freq_stem_AT_GC, file=paste0(path,'/freq_each/PARS_',g,'_stat_stem_AT_GC.csv'), row.names = FALSE)
  write.csv(dd_SNPs_freq_loop_AT_GC, file=paste0(path,'/freq_each/PARS_',g,'_stat_loop_AT_GC.csv'), row.names = FALSE)
  write.csv(dd_SNPs_freq_stem_GC_AT, file=paste0(path,'/freq_each/PARS_',g,'_stat_stem_GC_AT.csv'), row.names = FALSE)
  write.csv(dd_SNPs_freq_loop_GC_AT, file=paste0(path,'/freq_each/PARS_',g,'_stat_loop_GC_AT.csv'), row.names = FALSE)
}



# 把freq分组（先判断freq是否>=10)
group = c("4D","tRNA")
for (g in group){
  t <- get(paste0('data_SNPs_PARS_',g,collapse = NULL))
  if(max(t$freq)>=10){
    for(i in 1:10){
      # 统计每个freq的总SNPs和总gene的情况
      n <- assign(paste0('data_SNPs_PARS_',g,'_',i,collapse = NULL),subset(t,freq <= max(freq)*(i/10) & freq > max(freq)*((i-1)/10)))
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
    
    write.csv(data_stat, file=paste0(path,'/freq_10/PARS_',g,'_stat_freq_10.csv'), row.names = FALSE)
    
    write.csv(dd_SNPs_freq, file=paste0(path,'/freq_10/PARS_',g,'_stat_SNPs_freq_10.csv'), row.names = FALSE)
    write.csv(dd_gene_freq, file=paste0(path,'/freq_10/PARS_',g,'_stat_gene_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_stem, file=paste0(path,'/freq_10/PARS_',g,'_stat_stem_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_loop, file=paste0(path,'/freq_10/PARS_',g,'_stat_loop_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_stem_AT_GC, file=paste0(path,'/freq_10/PARS_',g,'_stat_stem_AT_GC_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_loop_AT_GC, file=paste0(path,'/freq_10/PARS_',g,'_stat_loop_AT_GC_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_stem_GC_AT, file=paste0(path,'/freq_10/PARS_',g,'_stat_stem_GC_AT_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_loop_GC_AT, file=paste0(path,'/freq_10/PARS_',g,'_stat_loop_GC_AT_freq_10.csv'), row.names = FALSE)
  }
  
}
