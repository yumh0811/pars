setwd("/Volumes/Backup/yumh/data/vcf/1011Matrix.gvcf/")

chr1_pars <- read.csv("chr1.pars.tsv", head= TRUE, sep = '\t')
chr1_ext <- read.csv("chr1.subset.ext.tsv", head= TRUE, sep = '\t')

chr1_merge <- merge(chr1_pars, chr1_ext, by = "name", all = F)
write.table (chr1_merge, file = "chr1.merge.subset.tsv", row.names = F,sep = "\t")
