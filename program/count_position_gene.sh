export NAME=Scer_n157_nonMosaic_Spar
cd ~/data/mrna-structure/result/$NAME 
#perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --origin data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds_pos.csv
perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --origin data_SNPs_PARS_mRNA.csv --output data_SNPs_PARS_mRNA_pos.csv
unset NAME
