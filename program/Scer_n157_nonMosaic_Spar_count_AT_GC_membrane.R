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

#输入csv
file_SNPs_PARS_cds <- paste0(path,'/data_SNPs_PARS_cds.csv',collapse = NULL)
data_SNPs_PARS_cds <- read.csv(file_SNPs_PARS_cds,header = TRUE,sep = ",")
data_gene_process <- data_SNPs_PARS_cds["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
dd_gene<- data.frame(name = "cds", gene = c(nrow(data_gene_process)))

file_gene_nucle <- "~/data/mrna-structure/phylogeny/membrane/Scer_n157_nonMosaic_gene_nucleoplasm.list"
data_gene_nucle <- read.csv(file_gene_nucle,header = TRUE,sep = ",")
data_SNPs_PARS_cds_nucle <- merge(data_SNPs_PARS_cds,data_gene_nucle,by="gene")
data_gene_process <- data_SNPs_PARS_cds_nucle["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
dd_gene <- rbind(dd_gene, data.frame(name = "cds_nucle", gene = c(nrow(data_gene_process))))
write.csv(data_gene_process , file=paste0(path,'/freq_10/membrane/PARS_cds_gene_nucle.csv'), row.names = FALSE)
write.csv(data_SNPs_PARS_cds_nucle , file=paste0(path,'/freq_10/membrane/data_SNPs_PARS_cds_nucle.csv'), row.names = FALSE)

file_gene_golgi <- "~/data/mrna-structure/phylogeny/membrane/Scer_n157_nonMosaic_gene_Golgi_membrane.list"
data_gene_golgi <- read.csv(file_gene_golgi,header = TRUE,sep = ",")
data_SNPs_PARS_cds_golgi <- merge(data_SNPs_PARS_cds,data_gene_golgi,by="gene")
data_gene_process <- data_SNPs_PARS_cds_golgi["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
dd_gene <- rbind(dd_gene, data.frame(name = "cds_golgi", gene = c(nrow(data_gene_process))))
write.csv(data_gene_process , file=paste0(path,'/freq_10/membrane/PARS_cds_gene_golgi.csv'), row.names = FALSE)
write.csv(data_SNPs_PARS_cds_golgi , file=paste0(path,'/freq_10/membrane/data_SNPs_PARS_cds_golgi.csv'), row.names = FALSE)

file_gene_mitout <- "~/data/mrna-structure/phylogeny/membrane/Scer_n157_nonMosaic_gene_mitochondrial_outer_membrane.list"
data_gene_mitout <- read.csv(file_gene_mitout,header = TRUE,sep = ",")
data_SNPs_PARS_cds_mitout <- merge(data_SNPs_PARS_cds,data_gene_mitout,by="gene")
data_gene_process <- data_SNPs_PARS_cds_mitout["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
dd_gene <- rbind(dd_gene, data.frame(name = "cds_mitout", gene = c(nrow(data_gene_process))))
write.csv(data_gene_process , file=paste0(path,'/freq_10/membrane/PARS_cds_gene_mitout.csv'), row.names = FALSE)
write.csv(data_SNPs_PARS_cds_mitout , file=paste0(path,'/freq_10/membrane/data_SNPs_PARS_cds_mitout.csv'), row.names = FALSE)

file_gene_mitin <- "~/data/mrna-structure/phylogeny/membrane/Scer_n157_nonMosaic_gene_mitochondrial_inner_membrane.list"
data_gene_mitin <- read.csv(file_gene_mitin,header = TRUE,sep = ",")
data_SNPs_PARS_cds_mitin <- merge(data_SNPs_PARS_cds,data_gene_mitin,by="gene")
data_gene_process <- data_SNPs_PARS_cds_mitin["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
dd_gene <- rbind(dd_gene, data.frame(name = "cds_mitin", gene = c(nrow(data_gene_process))))
write.csv(data_gene_process , file=paste0(path,'/freq_10/membrane/PARS_cds_gene_mitin.csv'), row.names = FALSE)
write.csv(data_SNPs_PARS_cds_mitin , file=paste0(path,'/freq_10/membrane/data_SNPs_PARS_cds_mitin.csv'), row.names = FALSE)

file_gene_ftv <- "~/data/mrna-structure/phylogeny/membrane/Scer_n157_nonMosaic_gene_fungal_type_vacuole_membrane.list"
data_gene_ftv <- read.csv(file_gene_ftv,header = TRUE,sep = ",")
data_SNPs_PARS_cds_ftv <- merge(data_SNPs_PARS_cds,data_gene_ftv,by="gene")
data_gene_process <- data_SNPs_PARS_cds_ftv["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
dd_gene <- rbind(dd_gene, data.frame(name = "cds_ftv", gene = c(nrow(data_gene_process))))
write.csv(data_gene_process , file=paste0(path,'/freq_10/membrane/PARS_cds_gene_ftv.csv'), row.names = FALSE)
write.csv(data_SNPs_PARS_cds_ftv , file=paste0(path,'/freq_10/membrane/data_SNPs_PARS_cds_ftv.csv'), row.names = FALSE)

file_gene_ER <- "~/data/mrna-structure/phylogeny/membrane/Scer_n157_nonMosaic_gene_endoplasmic_reticulum_membrane.list"
data_gene_ER <- read.csv(file_gene_ER,header = TRUE,sep = ",")
data_SNPs_PARS_cds_ER <- merge(data_SNPs_PARS_cds,data_gene_ER,by="gene")
data_gene_process <- data_SNPs_PARS_cds_ER["gene"]
data_gene_process <- unique(data_gene_process,fromLast=TRUE)
dd_gene <- rbind(dd_gene, data.frame(name = "cds_ER", gene = c(nrow(data_gene_process))))
write.csv(data_gene_process , file=paste0(path,'/freq_10/membrane/PARS_cds_gene_ER.csv'), row.names = FALSE)
write.csv(data_SNPs_PARS_cds_ER , file=paste0(path,'/freq_10/membrane/data_SNPs_PARS_cds_ER.csv'), row.names = FALSE)

write.csv(dd_gene , file=paste0(path,'/freq_10/membrane/PARS_cds_gene_membrane.csv'), row.names = FALSE)


# 把freq分组（先判断freq是否>=10)
group = c("cds_nucle","cds_golgi","cds_mitout","cds_mitin","cds_ftv","cds_ER")
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
    
    write.csv(data_stat, file=paste0(path,'/freq_10/membrane/PARS_',g,'_stat_freq_10.csv'), row.names = FALSE)
    
    write.csv(dd_SNPs_freq, file=paste0(path,'/freq_10/membrane/PARS_',g,'_stat_SNPs_freq_10.csv'), row.names = FALSE)
    write.csv(dd_gene_freq, file=paste0(path,'/freq_10/membrane/PARS_',g,'_stat_gene_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_stem, file=paste0(path,'/freq_10/membrane/PARS_',g,'_stat_stem_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_loop, file=paste0(path,'/freq_10/membrane/PARS_',g,'_stat_loop_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_stem_AT_GC, file=paste0(path,'/freq_10/membrane/PARS_',g,'_stat_stem_AT_GC_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_loop_AT_GC, file=paste0(path,'/freq_10/membrane/PARS_',g,'_stat_loop_AT_GC_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_stem_GC_AT, file=paste0(path,'/freq_10/membrane/PARS_',g,'_stat_stem_GC_AT_freq_10.csv'), row.names = FALSE)
    write.csv(dd_SNPs_freq_loop_GC_AT, file=paste0(path,'/freq_10/membrane/PARS_',g,'_stat_loop_GC_AT_freq_10.csv'), row.names = FALSE)
  }

}
