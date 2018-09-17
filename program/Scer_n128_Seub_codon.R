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

file_Scer_n128_snp <- paste0(path,"/Scer_n128_Seub.snp.csv")
data_Scer_n128_snp <- read.csv(file_Scer_n128_snp, header = TRUE, sep = ",")

data_Scer_n128_snp <- sqldf('SELECT * FROM [data_Scer_n128_snp] where snp_codon_pos == "2"' )
data_Scer_n128_snp <- sqldf('SELECT * FROM [data_Scer_n128_snp] where snp_syn >0 AND snp_nsy == 0' )
data_Scer_n128_snp <- sqldf('SELECT * FROM [data_Scer_n128_snp] where snp_mutant_to != "Complex"' )

data_Scer_n128_snp <- subset(data_Scer_n128_snp,select=c("snp_id","snp_pos","snp_target_base","snp_query_base","snp_outgroup_base","snp_mutant_to","snp_freq","snp_codon_pos","snp_codons","snp_syn","snp_nsy","snp_stop"))
write.csv(data_Scer_n128_snp,file = paste0(path,"/Scer_n128_Seub.snp.codon.csv"),row.names = F)
