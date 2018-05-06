setwd("~/data/mrna-structure/phylogeny/")
library(plyr)

data_gene_range = read.csv("~/data/mrna-structure/phylogeny/Scer_n157_nonMosaic_Spar_gene_range.csv",head=TRUE)
data_cds_range_chr = read.table("~/data/mrna-structure/phylogeny/protein_coding_list_range_chr.csv",head=FALSE,sep = ",")
colnames(data_cds_range_chr) = c("gene","chr","start","end")
data_gene_range = merge(data_cds_range_chr,data_gene_range,by="gene")
data_out <- arrange(data_gene_range,data_gene_range$chr,data_gene_range$start)
data_out_1 = subset(data_out, data_out$proporation==1)
write.csv(data_out_1,file="~/data/mrna-structure/phylogeny/Scer_n157_nonMosaic_Spar_distance_processed_pro_1.csv",row.names = FALSE)

data_mean_distance = read.csv("~/data/mrna-structure/phylogeny/Scer_n157_nonMosaic_Spar_mean_distance.csv",head=TRUE)
data_out_1 = merge(data_out_1,data_mean_distance,by="gene")
data_out_1 = subset(data_out_1, select = c(gene,chr,start,end,alignment_cds_length,sgd_cds_length,intersection_length,proporation,asian,wine,beer1,beer2,mixed))

data_gene_nodom <- data.frame()
data_gene_limdom <- data.frame()
data_gene_strdom <- data.frame()

for (gene in data_out_1$gene){
  data_gene = data_out_1[which(data_out_1$gene == gene),]
  
  # compare asian,wine,beer1,beer2
  min = min(c(data_gene$asian,data_gene$wine,data_gene$beer1,data_gene$beer2))
  
  # "min" represents the shortest distance, namely closer evolution. 
  if (data_gene$asian == min){
    data_gene_nodom <- rbind(data_gene_nodom, data.frame(gene = data_gene$gene))
  }else if (data_gene$wine == min){
    data_gene_limdom <- rbind(data_gene_limdom, data.frame(gene = data_gene$gene))
  }else if (data_gene$beer1 == min | data_gene$beer2 == min){
    data_gene_strdom <- rbind(data_gene_strdom, data.frame(gene = data_gene$gene))
  }
}

write.table(data_gene_nodom,file="~/data/mrna-structure/phylogeny/Scer_n157_nonMosaic_gene_nodom.list",row.names = FALSE)
write.table(data_gene_limdom,file="~/data/mrna-structure/phylogeny/Scer_n157_nonMosaic_gene_limdom.list",row.names = FALSE)
write.table(data_gene_strdom,file="~/data/mrna-structure/phylogeny/Scer_n157_nonMosaic_gene_strdom.list",row.names = FALSE)
