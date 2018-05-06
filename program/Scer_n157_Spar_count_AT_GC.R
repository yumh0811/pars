name <- "Scer_n157_Spar" 
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

#输入csv
file_SNPs_PARS_cds <- paste0(path,'/data_SNPs_PARS_cds.csv',collapse = NULL)
data_SNPs_PARS_cds <- read.csv(file_SNPs_PARS_cds,header = TRUE,sep = ",")
file_SNPs_PARS_utr <- paste0(path,'/data_SNPs_PARS_utr.csv',collapse = NULL)
data_SNPs_PARS_utr <- read.csv(file_SNPs_PARS_utr,header = TRUE,sep = ",")
file_SNPs_PARS_syn <- paste0(path,'/data_SNPs_PARS_syn.csv',collapse = NULL)
data_SNPs_PARS_syn <- read.csv(file_SNPs_PARS_syn,header = TRUE,sep = ",")
file_SNPs_PARS_nsy <- paste0(path,'/data_SNPs_PARS_nsy.csv',collapse = NULL)
data_SNPs_PARS_nsy <- read.csv(file_SNPs_PARS_nsy,header = TRUE,sep = ",")

group = c("cds","utr","syn","nsy")
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
group = c("cds","utr","syn","nsy")
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
