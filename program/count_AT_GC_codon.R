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

# name <- "Scer_n128_Spar"
name <- opt$name
area <- opt$area
path <- paste0("~/data/mrna-structure/result/", name, collapse = NULL)
setwd(path)

#输入csv
file_SNPs_PARS_cds <- paste0(path,'/data_SNPs_PARS_syn_codon.csv',collapse = NULL)
data_SNPs_PARS_cds <- read.csv(file_SNPs_PARS_cds,header = TRUE,sep = ",")

#codon
data_SNPs_PARS_4D <- sqldf('SELECT * FROM [data_SNPs_PARS_cds] where Codons ==  
                           "CCT->CCC" OR Codons == "CCT->CCA" OR Codons == "CCT->CCG" OR Codons == "CCC->CCT" OR Codons == "CCC->CCA" OR Codons == "CCC->CCG" OR Codons == "CCA->CCC" OR Codons == "CCA->CCT" OR Codons == "CCA->CCG" OR Codons == "CCG->CCC" OR Codons == "CCG->CCT" OR Codons == "CCG->CCA" OR Codons == 
                           "GTT->GTC" OR Codons == "GTT->GTA" OR Codons == "GTT->GTG" OR Codons == "GTC->GTT" OR Codons == "GTC->GTA" OR Codons == "GTC->GTG" OR Codons == "GTA->GTC" OR Codons == "GTA->GTT" OR Codons == "GTA->GTG" OR Codons == "GTG->GTC" OR Codons == "GTG->GTT" OR Codons == "GTG->GTA" OR Codons == 
                           "ACT->ACC" OR Codons == "ACT->ACA" OR Codons == "ACT->ACG" OR Codons == "ACC->ACT" OR Codons == "ACC->ACA" OR Codons == "ACC->ACG" OR Codons == "ACA->ACC" OR Codons == "ACA->ACT" OR Codons == "ACA->ACG" OR Codons == "ACG->ACC" OR Codons == "ACG->ACT" OR Codons == "ACG->ACA" OR Codons == 
                           "GCT->GCC" OR Codons == "GCT->GCA" OR Codons == "GCT->GCG" OR Codons == "GCC->GCT" OR Codons == "GCC->GCA" OR Codons == "GCC->GCG" OR Codons == "GCA->GCC" OR Codons == "GCA->GCT" OR Codons == "GCA->GCG" OR Codons == "GCG->GCC" OR Codons == "GCG->GCT" OR Codons == "GCG->GCA" OR Codons == 
                           "GGT->GGC" OR Codons == "GGT->GGA" OR Codons == "GGT->GGG" OR Codons == "GGC->GGT" OR Codons == "GGC->GGA" OR Codons == "GGC->GGG" OR Codons == "GGA->GGC" OR Codons == "GGA->GGT" OR Codons == "GGA->GGG" OR Codons == "GGG->GGC" OR Codons == "GGG->GGT" OR Codons == "GGG->GGA"
                           ')
write.csv(data_SNPs_PARS_4D,file=paste0(path,"/data_SNPs_PARS_4D.csv",collapse = NULL),row.names = FALSE)

#tRNA
data_SNPs_PARS_tRNA <- sqldf('SELECT * FROM [data_SNPs_PARS_cds] where 
                             Codons == "CTA->CTG" OR Codons == "CTG->CTA" OR 
                             Codons == "GTT->GTC" OR Codons == "GTT->GTA" OR Codons == "GTC->GTT" OR Codons == "GTC->GTA" OR Codons == "GTA->GTT" OR Codons == "GTA->GTC" OR
                             Codons == "GGT->GGC" OR Codons == "GGC->GGT" OR
                             Codons == "GCT->GCC" OR Codons == "GCT->GCA" OR Codons == "GCC->GCT" OR Codons == "GCC->GCA" OR Codons == "GCA->GCT" OR Codons == "GCA->GCC" OR
                             Codons == "AGA->AGG" OR Codons == "AGG->AGA" OR
                             Codons == "CGT->CGC" OR Codons == "CGT->CGA" OR Codons == "CGC->CGT" OR Codons == "CGC->CGA" OR Codons == "CGA->CGT" OR Codons == "CGA->CGC" OR
                             Codons == "AAA->AAG" OR Codons == "AAG->AAA" OR
                             Codons == "GAA->GAG" OR Codons == "GAG->GAA" OR
                             Codons == "GAT->GAC" OR Codons == "GAC->GAT" OR
                             Codons == "ACT->ACC" OR Codons == "ACT->ACA" OR Codons == "ACC->ACT" OR Codons == "ACC->ACA" OR Codons == "ACA->ACT" OR Codons == "ACA->ACC" OR
                             Codons == "TAT->TAC" OR Codons == "TAC->TAT" OR
                             Codons == "TCT->TCC" OR Codons == "TCT->TCA" OR Codons == "TCC->TCT" OR Codons == "TCC->TCA" OR Codons == "TCA->TCT" OR Codons == "TCA->TCC" OR
                             Codons == "TCA->TCG" OR Codons == "TCG->TCA" OR
                             Codons == "CAT->CAC" OR Codons == "CAC->CAT" OR
                             Codons == "TTT->TTC" OR Codons == "TTC->TTT" OR
                             Codons == "TGT->TGC" OR Codons == "TGC->TGT" ')
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
