export NAME=Scer_n157_nonMosaic_Spar
cd ~/data/mrna-structure/result/$NAME
perl ~/Scripts/pars/program/count_structure_length_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --name ~/data/mrna-structure/result/Scer_n157_nonMosaic_Spar/PARS_cds_gene.csv --structure stem --output stem_length_cds.csv
perl ~/Scripts/pars/program/count_structure_length_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --name ~/data/mrna-structure/result/Scer_n157_nonMosaic_Spar/PARS_cds_gene.csv --structure loop --output loop_length_cds.csv
unset NAME
