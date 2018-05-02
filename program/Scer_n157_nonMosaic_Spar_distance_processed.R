setwd("~/data/mrna-structure/phylogeny/")

data_gene_range = read.csv("~/data/mrna-structure/phylogeny/Scer_n157_nonMosaic_Spar_gene_range.csv",head=TRUE)
data_mean_distance = read.csv("~/data/mrna-structure/phylogeny/Scer_n157_nonMosaic_Spar_mean_distance.csv",head=TRUE)

#data_gene_range = subset(data_gene_range, data_gene_range$proporation==1)
data_mean_distance = merge(data_gene_range,data_mean_distance,by="gene")

data_cds_range_chr = read.table("~/data/mrna-structure/phylogeny/protein_coding_list_range_chr.csv",head=FALSE,sep = ",")
colnames(data_cds_range_chr) = c("gene","chr","start","end")

data_out = merge(data_cds_range_chr,data_mean_distance,by="gene")

write.csv(data_out,file="~/data/mrna-structure/phylogeny/Scer_n157_nonMosaic_Spar_distance_processed.csv",row.names = FALSE)