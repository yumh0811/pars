export NAME=Scer_n157_nonMosaic_Spar
cd ~/data/mrna-structure/result/$NAME 
perl ~/Scripts/pars/program/count_per_gene_ACGT_percent.pl --file data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds_per_gene_ATGC.csv
perl ~/Scripts/pars/program/count_per_gene_ACGT_percent.pl --file data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr_per_gene_ATGC.csv
unset NAME
 