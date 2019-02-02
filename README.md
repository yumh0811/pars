# Processing Yeast PARS Data

[TOC level=1-3]: # " "
- [Processing Yeast PARS Data](#processing-yeast-pars-data)
- [Install needed softwares](#install-needed-softwares)
- [Download reference data](#download-reference-data)
    - [Download PARS10 full site.](#download-pars10-full-site)
    - [SGD](#sgd)
- [Download genomes of strains](#download-genomes-of-strains)
    - [Download S288c (soft-masked) from Ensembl](#download-s288c-soft-masked-from-ensembl)
    - [Download strains from NCBI assembly](#download-strains-from-ncbi-assembly)
    - [Download strains from NCBI WGS](#download-strains-from-ncbi-wgs)
    - [Download strains from 1002genomes project](#download-strains-from-1002genomes-project)
- [Prepare sequences (RepeatMasker)](#prepare-sequences-repeatmasker)
- [Align](#align)
    - [Sanger](#sanger)
    - [PacBio](#pacbio)
    - [Illumina](#illumina)
- [Blast](#blast)
- [Gene_filiter](#gene_filiter)
    - [create protein coding gene list](#create-protein-coding-gene-list)
    - [cut mRNA alignment](#cut-mrna-alignment)
        - [create mRNA_yml](#create-mrna_yml)
        - [cut alignment by mRNA_yml](#cut-alignment-by-mrna_yml)
        - [count mRNA_alignment proporation in sgd](#count-mrna_alignment-proporation-in-sgd)
- [Extract SNP-list](#extract-snp-list)
- [Features](#features)
- [Real Processing n7](#real-processing-n7)
- [Real Processing n7p](#real-processing-n7p)
- [Real Processing n128](#real-processing-n128)
- [SNP](#snp)
    - [count per gene GC content](#count-per-gene-gc-content)
    - [count SNPs and gene](#count-snps-and-gene)
    - [vcf](#vcf)
    - [update](#update)
    - [count A/T <-> G/C](#count-at---gc)
    - [count stem length selection](#count-stem-length-selection)
    - [count_codon_gene](#count_codon_gene)
    - [count per gene cds_utr](#count-per-gene-cds_utr)
    - [count GO KEGG](#count-go-kegg)
    - [stat subpopulation SNPs frequency](#stat-subpopulation-snps-frequency)


# Install needed softwares

```bash
brew tap wang-q/tap
brew install multiz faops

```

# Download reference data

## Download PARS10 full site.

```bash
mkdir -p ~/data/mrna-structure/PARS10
cd ~/data/mrna-structure/PARS10

perl ~/Scripts/download/list.pl -u http://genie.weizmann.ac.il/pubs/PARS10/
perl ~/Scripts/download/download.pl -i pubs_PARS10.yml

find . -name "*.gz" |
    parallel -j 1 'gzip -dcf {} > {.}'

```

## SGD

```bash
mkdir -p ~/data/mrna-structure/sgd
cd ~/data/mrna-structure/sgd

aria2c -c http://downloads.yeastgenome.org/sequence/S288C_reference/intergenic/NotFeature.fasta.gz
aria2c -c http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz
aria2c -c http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_genomic_all.fasta.gz
aria2c -c http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff

find . -name "*.gz" |
    parallel -j 1 'gzip -dcf {} > {.}'

```

# Download genomes of strains

## Download S288c (soft-masked) from Ensembl

```bash
mkdir -p ~/data/mrna-structure/ensembl/
cd ~/data/mrna-structure/ensembl/

aria2c -x 6 -s 3 -c ftp://ftp.ensembl.org/pub/release-94/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz
aria2c -x 6 -s 3 -c ftp://ftp.ensembl.org/pub/release-94/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.94.gff3.gz

find . -name "*.gz" | xargs gzip -t

```

## Download strains from NCBI assembly

```bash
cd ~/data/mrna-structure/

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/pars/scer.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/scer.assembly.rsync.sh

bash ASSEMBLY/scer.assembly.collect.sh

```

## Download strains from NCBI WGS

```bash
cd ~/data/mrna-structure/

perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/pars/scer.wgs.tsv \
    --fix \
    -o WGS

bash WGS/scer.wgs.rsync.sh

```

## Download strains from 1002genomes project

```bash
mkdir -p ~/data/mrna-structure/download/
cd ~/data/mrna-structure/download/

wget -c http://1002genomes.u-strasbg.fr/files/1011Assemblies.tar.gz

#tar -zxvf 1011Assemblies.tar.gz

```

# Prepare sequences (RepeatMasker)

```bash
cd ~/data/mrna-structure/

# reference
egaz prepseq \
    ensembl/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz \
    --repeatmasker "--species Fungi --gff --parallel 12" \
    --min 1000 --gi -v \
    -o GENOMES/S288c

gzip -dcf ensembl/Saccharomyces_cerevisiae.R64-1-1.94.gff3.gz > GENOMES/S288c/chr.gff

# prep assembly
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    --min 1000 --about 1_000_000 \
    -v --repeatmasker "--species Fungi --parallel 12"

bash GENOMES/0_prep.sh

# prep wgs
egaz template \
    WGS \
    --prep -o GENOMES \
    --min 1000 --about 1_000_000 \
    -v --repeatmasker "--species Fungi --parallel 12"

bash GENOMES/0_prep.sh

```

# Align

## Sanger

```bash
mkdir -p ~/data/mrna-structure/alignment
cd ~/data/mrna-structure/alignment

ln -s ~/data/mrna-structure/GENOMES .

egaz template \
    GENOMES/S288c GENOMES/EC1118 GENOMES/Kyokai_no_7 GENOMES/RM11_1a \
    GENOMES/Sigma1278b GENOMES/T7 GENOMES/YJM789 GENOMES/Spar \
    --multi -o n7 \
    --multiname Scer_n7_Spar --outgroup Spar \
    --vcf --aligndb \
    --order -v --parallel 12

bash n7/1_pair.sh
bash n7/3_multi.sh
bash n7/6_chr_length.sh
bash n7/7_multi_aligndb.sh

# clean
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw"   | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

## PacBio

```bash
cd ~/data/mrna-structure/alignment

egaz template \
    GENOMES/S288c GENOMES/DBVPG6044 GENOMES/UWOPS03_461_4 GENOMES/Y12 \
    GENOMES/SK1 GENOMES/YPS128 GENOMES/DBVPG6765 GENOMES/Spar \
    --multi -o n7p \
    --multiname Scer_n7p_Spar --outgroup Spar \
    --vcf --aligndb \
    --order -v --parallel 12

bash n7p/1_pair.sh
bash n7p/3_multi.sh
bash n7p/6_chr_length.sh
bash n7p/7_multi_aligndb.sh

```

## Illumina

```bash
cd ~/data/mrna-structure/alignment

egaz template \
    GENOMES/S288c \
    $(
        cat ~/Scripts/pars/group_phylo.tsv |
            grep -v "^#" |
            cut -f 2 |
            tr "," "\n" |
            sed 's/^/GENOMES\//'
    ) \
    GENOMES/Spar GENOMES/Seub \
    --multi -o n128 \
    -v --parallel 16

bash n128/1_pair.sh

egaz template \
    GENOMES/S288c \
    $(
        cat ~/Scripts/pars/group_phylo.tsv |
            grep -v "^#" |
            cut -f 2 |
            tr "," "\n" |
            sed 's/^/GENOMES\//'
    ) \
    GENOMES/Spar \
    --multi -o n128 \
    --multiname Scer_n128_Spar --outgroup Spar \
    --vcf --aligndb \
    --order -v --parallel 16

bash n128/3_multi.sh
bash n128/6_chr_length.sh
bash n128/7_multi_aligndb.sh

egaz template \
    GENOMES/S288c \
    $(
        cat ~/Scripts/pars/group_phylo.tsv |
            grep -v "^#" |
            cut -f 2 |
            tr "," "\n" |
            sed 's/^/GENOMES\//'
    ) \
    GENOMES/Seub \
    --multi -o n128 \
    --multiname Scer_n128_Seub --outgroup Seub \
    --vcf --aligndb \
    --order -v --parallel 16

bash n128/3_multi.sh
bash n128/6_chr_length.sh
bash n128/7_multi_aligndb.sh

```

# Blast

Prepare a combined fasta file of yeast genome and blast genes against the genome.

```bash
mkdir -p ~/data/mrna-structure/blast
cd ~/data/mrna-structure/blast

cat ~/data/alignment/egaz/S288c/{I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI,Mito}.fa \
    > S288c.fa

perl -nl -i -e '/^>/ or $_ = uc $_; print'  S288c.fa
faops size S288c.fa > S288c.sizes

# formatdb
~/share/blast/bin/formatdb -p F -o T -i S288c.fa

# blast every transcripts against genome
~/share/blast/bin/blastall -p blastn -F "m D" -m 0 -b 10 -v 10 -e 1e-3 -a 4 \
    -i ../PARS10/pubs/PARS10/data/sce_genes.fasta -d S288C.fa -o sce_genes.blast
    
# parse blastn output
perl ~/Scripts/pars/blastn_transcript.pl -f sce_genes.blast -m 0

```

# Gene_filiter

## create protein coding gene list

```bash
mkdir -p ~/data/mrna-structure/gene_filiter
cd ~/data/mrna-structure/gene_filiter

# sgd/saccharomyces_cerevisiae.gff → protein coding gene list
perl ~/Scripts/pars/program/protein_coding_list.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list.csv
perl ~/Scripts/pars/program/protein_coding_list_range.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list_range.csv
perl ~/Scripts/pars/program/protein_coding_list_range_chr.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list_range_chr.csv
perl ~/Scripts/pars/program/protein_coding_list_range_chr_strand.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list_range_chr_strand.csv

# sgd non-overlap
perl ~/Scripts/pars/program/protein_coding_overlap.pl --file ~/data/mrna-structure/gene_filiter/protein_coding_list_range_chr_strand.csv --output ~/data/mrna-structure/gene_filiter/protein_coding_overlap.csv
perl ~/Scripts/pars/program/protein_coding_overlap_yml.pl --file ~/data/mrna-structure/gene_filiter/protein_coding_overlap.csv --output ~/data/mrna-structure/gene_filiter/protein_coding_overlap.yml
cat ~/data/mrna-structure/gene_filiter/protein_coding_overlap.csv | perl -nl -a -F"," -e 'print qq{$F[0]};' | sort | uniq | sed -e "/target_gene/d"  > protein_coding_overlap_unique.csv
cat protein_coding_list.csv protein_coding_overlap_unique.csv | sort | uniq -u > protein_coding_non_overlap_unique.csv

# PARS genes
perl ~/Scripts/pars/program/PARS_genes_list_range.pl --file ~/data/mrna-structure/blast/sce_genes.blast.tsv --output ~/data/mrna-structure/gene_filiter/PARS_genes_list_range.csv
perl ~/Scripts/pars/program/PARS_genes_list_range_chr.pl --file ~/data/mrna-structure/blast/sce_genes.blast.tsv --output ~/data/mrna-structure/gene_filiter/PARS_genes_list_range_chr.csv
perl ~/Scripts/pars/program/PARS_genes_list_range_chr_strand.pl --file ~/data/mrna-structure/blast/sce_genes.blast.tsv --output ~/data/mrna-structure/gene_filiter/PARS_genes_list_range_chr_strand.csv
cat ~/data/mrna-structure/gene_filiter/PARS_genes_list_range_chr_strand.csv | perl -nl -a -F"," -e 'print qq{$F[0]};'  > PARS_genes_list.csv

#merge non-overlap
cat protein_coding_non_overlap_unique.csv PARS_genes_list.csv | sort | uniq -d > PARS_genes_non_overlap_unique_protein_coding.csv

```

## cut mRNA alignment

### create mRNA_yml

```bash
cd ~/data/mrna-structure/gene_filiter

#cut mRNA in PARS_Blast
mkdir -p ~/data/mrna-structure/gene_filiter/gene_mRNA_yml
perl ~/Scripts/pars/program/cut_mRNA_yml.pl --file protein_coding_list_range_chr.csv --output gene_mRNA_yml
```

### cut alignment by mRNA_yml

```bash

export NAME=Scer_n7_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/multi8/${NAME}_refined ~/data/mrna-structure/gene_filiter/${NAME}_refined
cd ~/data/mrna-structure/gene_filiter/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/gene_filiter/${NAME}_gene_alignment_mRNA
cd ~/data/mrna-structure/gene_filiter/${NAME}_gene_alignment_mRNA
cat ../PARS_genes_non_overlap_unique_protein_coding.csv |
   parallel --line-buffer -j 16 '
       fasops slice ../${NAME}_refined/species.fas ../gene_mRNA_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n7p_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/multi8p/${NAME}_refined ~/data/mrna-structure/gene_filiter/${NAME}_refined
cd ~/data/mrna-structure/gene_filiter/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/gene_filiter/${NAME}_gene_alignment_mRNA
cd ~/data/mrna-structure/gene_filiter/${NAME}_gene_alignment_mRNA
cat ../PARS_genes_non_overlap_unique_protein_coding.csv |
   parallel --line-buffer -j 16 '
       fasops slice ../${NAME}_refined/species.fas ../gene_mRNA_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n128_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/multi128_Spar/${NAME}_refined ~/data/mrna-structure/gene_filiter/${NAME}_refined
cd ~/data/mrna-structure/gene_filiter/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/gene_filiter/${NAME}_gene_alignment_mRNA
cd ~/data/mrna-structure/gene_filiter/${NAME}_gene_alignment_mRNA
cat ../PARS_genes_non_overlap_unique_protein_coding.csv |
   parallel --line-buffer -j 16 '
       fasops slice ../${NAME}_refined/species.fas ../gene_mRNA_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n128_Seub
cp -rf ~/data/mrna-structure/alignment/scer_wgs/multi128_Seub/${NAME}_refined ~/data/mrna-structure/gene_filiter/${NAME}_refined
cd ~/data/mrna-structure/gene_filiter/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/gene_filiter/${NAME}_gene_alignment_mRNA
cd ~/data/mrna-structure/gene_filiter/${NAME}_gene_alignment_mRNA
cat ../PARS_genes_non_overlap_unique_protein_coding.csv |
   parallel --line-buffer -j 16 '
       fasops slice ../${NAME}_refined/species.fas ../gene_mRNA_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME
```

### count mRNA_alignment proporation in sgd

```bash

export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/gene_filiter
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_mRNA --output ${NAME}_protein_coding_range.csv
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/gene_filiter
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_mRNA --output ${NAME}_protein_coding_range.csv
unset NAME

export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/gene_filiter
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_mRNA --output ${NAME}_protein_coding_range.csv
unset NAME

export NAME=Scer_n128_Seub
cd ~/data/mrna-structure/gene_filiter
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_mRNA --output ${NAME}_protein_coding_range.csv
unset NAME
```

```bash
#生成alignment_proporation_1.list

export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/gene_filiter
Rscript ~/Scripts/pars/program/proporation_1.R -i ${NAME}_protein_coding_range.csv -r protein_coding_list_range_chr.csv -o ${NAME}_non-overlap_pro_1.csv
cat ~/data/mrna-structure/gene_filiter/${NAME}_non-overlap_pro_1.csv | perl -nl -a -F"," -e 'print qq{$F[0]};' | sed "s/\"//g" | sed -e "/gene/d" > ${NAME}_final_genes.csv
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/gene_filiter
Rscript ~/Scripts/pars/program/proporation_1.R -i ${NAME}_protein_coding_range.csv -r protein_coding_list_range_chr.csv -o ${NAME}_non-overlap_pro_1.csv
cat ~/data/mrna-structure/gene_filiter/${NAME}_non-overlap_pro_1.csv | perl -nl -a -F"," -e 'print qq{$F[0]};' | sed "s/\"//g" | sed -e "/gene/d" > ${NAME}_final_genes.csv
unset NAME

export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/gene_filiter
Rscript ~/Scripts/pars/program/proporation_1.R -i ${NAME}_protein_coding_range.csv -r protein_coding_list_range_chr.csv -o ${NAME}_non-overlap_pro_1.csv
cat ~/data/mrna-structure/gene_filiter/${NAME}_non-overlap_pro_1.csv | perl -nl -a -F"," -e 'print qq{$F[0]};' | sed "s/\"//g" | sed -e "/gene/d" > ${NAME}_final_genes.csv
unset NAME

export NAME=Scer_n128_Seub
cd ~/data/mrna-structure/gene_filiter
Rscript ~/Scripts/pars/program/proporation_1.R -i ${NAME}_protein_coding_range.csv -r protein_coding_list_range_chr.csv -o ${NAME}_non-overlap_pro_1.csv
cat ~/data/mrna-structure/gene_filiter/${NAME}_non-overlap_pro_1.csv | perl -nl -a -F"," -e 'print qq{$F[0]};' | sed "s/\"//g" | sed -e "/gene/d" > ${NAME}_final_genes.csv
unset NAME
```

# Extract SNP-list

```bash

mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

export NAME=Scer_n7_Spar
mkdir -p ${NAME}_snp
cat ../gene_filiter/${NAME}_final_genes.csv |
   parallel --line-buffer -j 16 '
       fasops vars --nosingle --outgroup --nocomplex ../gene_filiter/${NAME}_gene_alignment_mRNA/{}.fas.fas -o ${NAME}_snp/{}.SNPs.tsv
   '
cat ${NAME}_snp/*.SNPs.tsv | perl -nl -a -F"\t" -e 'print qq{$F[1]\t$F[3]\t$F[3]\t$F[5]/$F[6]};' > ${NAME}.total.SNPs.tsv
unset NAME

export NAME=Scer_n7p_Spar
mkdir -p ${NAME}_snp
cat ../gene_filiter/${NAME}_final_genes.csv |
   parallel --line-buffer -j 16 '
       fasops vars --nosingle --outgroup --nocomplex ../gene_filiter/${NAME}_gene_alignment_mRNA/{}.fas.fas -o ${NAME}_snp/{}.SNPs.tsv
   '
cat ${NAME}_snp/*.SNPs.tsv | perl -nl -a -F"\t" -e 'print qq{$F[1]\t$F[3]\t$F[3]\t$F[5]/$F[6]};' > ${NAME}.total.SNPs.tsv
unset NAME

export NAME=Scer_n128_Spar
mkdir -p ${NAME}_snp
cat ../gene_filiter/${NAME}_final_genes.csv |
   parallel --line-buffer -j 16 '
       fasops vars --nosingle --outgroup --nocomplex ../gene_filiter/${NAME}_gene_alignment_mRNA/{}.fas.fas -o ${NAME}_snp/{}.SNPs.tsv
   '
cat ${NAME}_snp/*.SNPs.tsv | perl -nl -a -F"\t" -e 'print qq{$F[1]\t$F[3]\t$F[3]\t$F[5]/$F[6]};' > ${NAME}.total.SNPs.tsv
unset NAME

export NAME=Scer_n128_Seub
mkdir -p ${NAME}_snp
cat ../gene_filiter/${NAME}_final_genes.csv |
   parallel --line-buffer -j 16 '
       fasops vars --nosingle --outgroup --nocomplex ../gene_filiter/${NAME}_gene_alignment_mRNA/{}.fas.fas -o ${NAME}_snp/{}.SNPs.tsv
   '
cat ${NAME}_snp/*.SNPs.tsv | perl -nl -a -F"\t" -e 'print qq{$F[1]\t$F[3]\t$F[3]\t$F[5]/$F[6]};' > ${NAME}.total.SNPs.tsv
unset NAME


```

# Features

```bash
mkdir -p ~/data/mrna-structure/process
cd ~/data/mrna-structure/process

#----------------------------------------------------------#
# gene
#----------------------------------------------------------#
# parse blastn output
perl ~/Scripts/pars/blastn_transcript.pl -f ../blast/sce_genes.blast -m 0

# produce transcript set
# YLR167W	568	chrXII	498888	499455	+
cat sce_genes.blast.tsv \
    | perl -nla -e 'print qq{$F[2]:$F[3]-$F[4]}' \
    | sort \
    > sce_genes.pos.txt
jrunlist cover sce_genes.pos.txt -o sce_genes.yml

#----------------------------------------------------------#
# intergenic
#----------------------------------------------------------#
cat ../sgd/NotFeature.fasta \
    | perl -n -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+(\d+)\-(\d+)/ or next;
        $1 eq "Mito" and next;
        print qq{$1:$2-$3\n};
    ' \
    > sce_intergenic.pos.txt
jrunlist cover sce_intergenic.pos.txt -o sce_intergenic.yml

#----------------------------------------------------------#
# intron
#----------------------------------------------------------#
cat ../sgd/orf_coding_all.fasta \
    | perl -n -MAlignDB::IntSpan -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+([\d,-]+)/ or next;
        $1 eq "Mito" and next;

        my $chr = $1;
        my $range = $2;
        my @ranges = sort { $a <=> $b } grep {/^\d+$/} split /,|\-/, $range;
        my $intspan = AlignDB::IntSpan->new()->add_range(@ranges);
        my $hole = $intspan->holes;

        printf qq{%s:%s\n}, $chr, $hole->as_string if $hole->is_not_empty;
    ' \
    > sce_intron.pos.txt
jrunlist cover sce_intron.pos.txt -o sce_intron.yml

#----------------------------------------------------------#
# utr
#----------------------------------------------------------#
# produce orf_genomic set
cat ../sgd/orf_genomic_all.fasta \
    | perl -n -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+(\d+)\-(\d+)/ or next;
        $1 eq "Mito" and next;

        if ($2 == $3) {
            print qq{$1:$2\n};
        }
        elsif ($2 < $3) {
            print qq{$1:$2-$3\n};
        }
        else {
            print qq{$1:$3-$2\n};
        }
    ' \
    > sce_orf_genomic.pos.txt
jrunlist cover sce_orf_genomic.pos.txt -o sce_orf_genomic.yml

jrunlist compare --op diff sce_genes.yml sce_orf_genomic.yml -o sce_utr.yml
runlist convert sce_utr.yml -o sce_utr.pos.txt

jrunlist compare --op diff sce_genes.yml sce_intron.yml -o sce_mRNA.yml
runlist convert sce_mRNA.yml -o sce_mRNA.pos.txt

jrunlist compare --op diff sce_mRNA.yml sce_utr.yml -o sce_cds.yml
runlist convert sce_cds.yml -o sce_cds.pos.txt

# Stats
printf "| %s | %s | %s | %s |\n" \
    "Name" "chrLength" "size" "coverage" \
    > coverage.stat.md
printf "|:--|--:|--:|--:|\n" >> coverage.stat.md

for f in genes intergenic intron orf_genomic utr mRNA cds; do
    printf "| %s | %s | %s | %s |\n" \
        ${f} \
        $(
            jrunlist stat ../blast/S288c.sizes sce_${f}.yml --all -o stdout \
            | grep -v coverage \
            | sed "s/,/ /g"
        )
done >> coverage.stat.md

cat coverage.stat.md
```

| Name        | chrLength |    size | coverage |
|:------------|----------:|--------:|---------:|
| genes       |  12071326 | 4235405 |   0.3509 |
| intergenic  |  12071326 | 2864170 |   0.2373 |
| intron      |  12071326 |   65144 |   0.0054 |
| orf_genomic |  12071326 | 8895737 |   0.7369 |
| utr         |  12071326 |  516569 |   0.0428 |
| mRNA        |  12071326 | 4233361 |   0.3507 |
| cds         |  12071326 | 3716792 |   0.3079 |

# Real Processing n7

```bash
export NAME=Scer_n7_Spar

cd ~/data/mrna-structure/process

# SNPs within transcripts
runlist position --op superset \
    sce_genes.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.gene.pos.txt

# read gene and snp info file
# produce ${NAME}.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl \
    --pars ../PARS10/pubs/PARS10/data \
    --gene sce_genes.blast.tsv \
    --pos  ${NAME}.snp.gene.pos.txt \
    > fail_pos.txt

# review fail_pos.txt to find SNPs located in overlapped genes

# process ${NAME}.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml

# SNPs within orf_genomic regions
runlist position --op superset \
    sce_orf_genomic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within cds
runlist position --op superset \
    sce_cds.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.cds.pos.txt
unset NAME
```

# Real Processing n7p

```bash
export NAME=Scer_n7p_Spar

cd ~/data/mrna-structure/process

# SNPs within transcripts
runlist position --op superset \
    sce_genes.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.gene.pos.txt

# read gene and snp info file
# produce ${NAME}.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl \
    --pars ../PARS10/pubs/PARS10/data \
    --gene sce_genes.blast.tsv \
    --pos  ${NAME}.snp.gene.pos.txt \
    > fail_pos.txt

# review fail_pos.txt to find SNPs located in overlapped genes

# process ${NAME}.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml

# SNPs within orf_genomic regions
runlist position --op superset \
    sce_orf_genomic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within cds
runlist position --op superset \
    sce_cds.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.cds.pos.txt
unset NAME
```

# Real Processing n128

```bash
#Spar
export NAME=Scer_n128_Spar

cd ~/data/mrna-structure/process

# SNPs within transcripts
runlist position --op superset \
    sce_genes.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.gene.pos.txt

# read gene and snp info file
# produce ${NAME}.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl \
    --pars ../PARS10/pubs/PARS10/data \
    --gene sce_genes.blast.tsv \
    --pos  ${NAME}.snp.gene.pos.txt \
    > fail_pos.txt

# review fail_pos.txt to find SNPs located in overlapped genes

# process ${NAME}.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml

# SNPs within orf_genomic regions
runlist position --op superset \
    sce_orf_genomic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within cds
runlist position --op superset \
    sce_cds.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.cds.pos.txt
unset NAME

#Seub
export NAME=Scer_n128_Seub

cd ~/data/mrna-structure/process

# SNPs within transcripts
runlist position --op superset \
    sce_genes.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.gene.pos.txt

# read gene and snp info file
# produce ${NAME}.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl \
    --pars ../PARS10/pubs/PARS10/data \
    --gene sce_genes.blast.tsv \
    --pos  ${NAME}.snp.gene.pos.txt \
    > fail_pos.txt

# review fail_pos.txt to find SNPs located in overlapped genes

# process ${NAME}.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml

# SNPs within orf_genomic regions
runlist position --op superset \
    sce_orf_genomic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within cds
runlist position --op superset \
    sce_cds.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.cds.pos.txt
unset NAME
```

# SNP

## count per gene GC content

```bash

export NAME=Scer_n7_Spar
mkdir -p ~/data/mrna-structure/result/${NAME}
cd ~/data/mrna-structure/result/${NAME}
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --varfold ~/data/mrna-structure/process/${NAME}.gene_variation.fold_class.tsv --output ${NAME}.gene_variation.fold_class.csv
unset NAME

export NAME=Scer_n7p_Spar
mkdir -p ~/data/mrna-structure/result/${NAME}
cd ~/data/mrna-structure/result/${NAME}
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --varfold ~/data/mrna-structure/process/${NAME}.gene_variation.fold_class.tsv --output ${NAME}.gene_variation.fold_class.csv
unset NAME

export NAME=Scer_n128_Spar
mkdir -p ~/data/mrna-structure/result/${NAME}
cd ~/data/mrna-structure/result/${NAME}
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --varfold ~/data/mrna-structure/process/${NAME}.gene_variation.fold_class.tsv --output ${NAME}.gene_variation.fold_class.csv
unset NAME

export NAME=Scer_n128_Seub
mkdir -p ~/data/mrna-structure/result/${NAME}
cd ~/data/mrna-structure/result/${NAME}
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --varfold ~/data/mrna-structure/process/${NAME}.gene_variation.fold_class.tsv --output ${NAME}.gene_variation.fold_class.csv
unset NAME
```

## count SNPs and gene

```bash

Rscript -e 'install.packages("getopt", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("ape", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("ggplot2", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("scales", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("reshape", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("pander", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("gridExtra", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("plyr", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("dplyr", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("proto", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("gsubfn", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("RSQLite", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("sqldf", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("sm", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'

export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/result/${NAME}
Rscript ~/Scripts/pars/program/stat_SNPs.R -n ${NAME}
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/result/${NAME}
Rscript ~/Scripts/pars/program/stat_SNPs.R -n ${NAME}
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME}
Rscript ~/Scripts/pars/program/stat_SNPs.R -n ${NAME}
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

export NAME=Scer_n128_Seub
cd ~/data/mrna-structure/result/${NAME}
Rscript ~/Scripts/pars/program/stat_SNPs.R -n ${NAME}
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME
```

## vcf

```bash
mkdir -p ~/data/mrna-structure/vcf
cd ~/data/mrna-structure/vcf
wget -c http://1002genomes.u-strasbg.fr/files/1011Matrix.gvcf.gz
gzip -d 1011Matrix.gvcf.gz

#rsync -av ymh@wq.nju.edu.cn:~/data/vcf/1011Matrix.gvcf /Volumes/Backup/yumh/data/vcf/1011Matrix.gvcf/
#ln -s /Volumes/Backup/yumh/data/vcf/1011Matrix.gvcf.gz .
#mkdir -p ~/data/mrna-structure/vcf/1011Matrix.gvcf
#cd 1011Matrix.gvcf/
#ln -s /Volumes/Backup/yumh/data/vcf/1011Matrix.gvcf/1011Matrix.gvcf .

# 1011
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf
perl ~/Scripts/pars/program/vcf.cut.pl --file 1011Matrix.gvcf --output 1011Matrix.tsv

export NAME=Scer_n7_Spar
mkdir -p ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}

perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds.pars.tsv
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr.pars.tsv
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_syn.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_syn.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file data_SNPs_PARS_syn.csv.bak --output data_SNPs_PARS_syn.pars.tsv
rm -rf data_SNPs_PARS_syn.csv.bak
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_nsy.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_nsy.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file data_SNPs_PARS_nsy.csv.bak --output data_SNPs_PARS_nsy.pars.tsv
rm -rf data_SNPs_PARS_nsy.csv.bak

perl ~/Scripts/pars/program/vcf.extract.pl --file ../1011Matrix.tsv --output 1011Matrix.ext.tsv

Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a cds
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a utr
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a syn
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a nsy

perl ~/Scripts/pars/program/vcf.merge.pro.pl --file cds_snp.merge.tsv --output ${NAME}.cds_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file utr_snp.merge.tsv --output ${NAME}.utr_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file syn_snp.merge.tsv --output ${NAME}.syn_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file nsy_snp.merge.tsv --output ${NAME}.nsy_snp.merge.pro.tsv

cat ${NAME}.cds_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.cds_snp.merge.pro.txt
cat ${NAME}.utr_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.utr_snp.merge.pro.txt
cat ${NAME}.syn_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.syn_snp.merge.pro.txt
cat ${NAME}.nsy_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.nsy_snp.merge.pro.txt
unset NAME

export NAME=Scer_n7p_Spar
mkdir -p ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}

perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds.pars.tsv
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr.pars.tsv
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_syn.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_syn.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file data_SNPs_PARS_syn.csv.bak --output data_SNPs_PARS_syn.pars.tsv
rm -rf data_SNPs_PARS_syn.csv.bak
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_nsy.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_nsy.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file data_SNPs_PARS_nsy.csv.bak --output data_SNPs_PARS_nsy.pars.tsv
rm -rf data_SNPs_PARS_nsy.csv.bak

perl ~/Scripts/pars/program/vcf.extract.pl --file ../1011Matrix.tsv --output 1011Matrix.ext.tsv

Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a cds
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a utr
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a syn
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a nsy

perl ~/Scripts/pars/program/vcf.merge.pro.pl --file cds_snp.merge.tsv --output ${NAME}.cds_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file utr_snp.merge.tsv --output ${NAME}.utr_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file syn_snp.merge.tsv --output ${NAME}.syn_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file nsy_snp.merge.tsv --output ${NAME}.nsy_snp.merge.pro.tsv

cat ${NAME}.cds_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.cds_snp.merge.pro.txt
cat ${NAME}.utr_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.utr_snp.merge.pro.txt
cat ${NAME}.syn_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.syn_snp.merge.pro.txt
cat ${NAME}.nsy_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.nsy_snp.merge.pro.txt
unset NAME

export NAME=Scer_n128_Spar
mkdir -p ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}

perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds.pars.tsv
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr.pars.tsv
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_syn.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_syn.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file data_SNPs_PARS_syn.csv.bak --output data_SNPs_PARS_syn.pars.tsv
rm -rf data_SNPs_PARS_syn.csv.bak
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_nsy.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_nsy.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file data_SNPs_PARS_nsy.csv.bak --output data_SNPs_PARS_nsy.pars.tsv
rm -rf data_SNPs_PARS_nsy.csv.bak

perl ~/Scripts/pars/program/vcf.extract.pl --file ../1011Matrix.tsv --output 1011Matrix.ext.tsv

Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a cds
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a utr
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a syn
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a nsy

perl ~/Scripts/pars/program/vcf.merge.pro.pl --file cds_snp.merge.tsv --output ${NAME}.cds_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file utr_snp.merge.tsv --output ${NAME}.utr_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file syn_snp.merge.tsv --output ${NAME}.syn_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file nsy_snp.merge.tsv --output ${NAME}.nsy_snp.merge.pro.tsv

cat ${NAME}.cds_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.cds_snp.merge.pro.txt
cat ${NAME}.utr_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.utr_snp.merge.pro.txt
cat ${NAME}.syn_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.syn_snp.merge.pro.txt
cat ${NAME}.nsy_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.nsy_snp.merge.pro.txt
unset NAME

export NAME=Scer_n128_Seub
mkdir -p ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}

perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds.pars.tsv
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr.pars.tsv
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_syn.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_syn.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file data_SNPs_PARS_syn.csv.bak --output data_SNPs_PARS_syn.pars.tsv
rm -rf data_SNPs_PARS_syn.csv.bak
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_nsy.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_nsy.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file data_SNPs_PARS_nsy.csv.bak --output data_SNPs_PARS_nsy.pars.tsv
rm -rf data_SNPs_PARS_nsy.csv.bak

perl ~/Scripts/pars/program/vcf.extract.pl --file ../1011Matrix.tsv --output 1011Matrix.ext.tsv

Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a cds
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a utr
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a syn
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME} -a nsy

perl ~/Scripts/pars/program/vcf.merge.pro.pl --file cds_snp.merge.tsv --output ${NAME}.cds_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file utr_snp.merge.tsv --output ${NAME}.utr_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file syn_snp.merge.tsv --output ${NAME}.syn_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file nsy_snp.merge.tsv --output ${NAME}.nsy_snp.merge.pro.tsv

cat ${NAME}.cds_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.cds_snp.merge.pro.txt
cat ${NAME}.utr_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.utr_snp.merge.pro.txt
cat ${NAME}.syn_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.syn_snp.merge.pro.txt
cat ${NAME}.nsy_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.nsy_snp.merge.pro.txt
unset NAME

# wild.strains in 1011
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf

perl -i -pe 's/chromosome4\t193242.*\n//g;s/chromosome4\t193246.*\n//g;s/chromosome4\t88:2:49\..*\n//g;s/chromosome4\t88:268.*\n//g;' 1011Matrix.gvcf

bcftools view 1011Matrix.gvcf -s CCL,BBQ,BBS,BFP,BTG,CLC,CLB,CLD,BAM,BAQ,BAG,BAH,BAL,AMH,CEG,CEI,CCQ,CCR,CCS,BAK,BAI,ACQ,CCN,CDL,SACE_YCR,BMA,AKM,BMB,BMC,SACE_MAL,SACE_YCY,BAN,BAP,CMP,CCH,ACC,CCC,CCD,CCE,CCF,CCG,CCI,CMQ,CDF,CDG,CDH,CDI,AVI,ACD,ANF,ANH,ANC,ANE,ANG,AND,ANK,ANI,AKN,SACE_YBS,SACE_YCU | bcftools +fill-tags -o 1011Matrix.wild.gvcf

perl ~/Scripts/pars/program/vcf.cut.pl --file 1011Matrix.wild.gvcf --output 1011Matrix.wild.tsv

export NAME=Scer_n7_Spar
mkdir -p ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}.wild
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}.wild

perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds.pars.tsv
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr.pars.tsv
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_syn.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_syn.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file data_SNPs_PARS_syn.csv.bak --output data_SNPs_PARS_syn.pars.tsv
rm -rf data_SNPs_PARS_syn.csv.bak
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_nsy.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_nsy.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file data_SNPs_PARS_nsy.csv.bak --output data_SNPs_PARS_nsy.pars.tsv
rm -rf data_SNPs_PARS_nsy.csv.bak

perl ~/Scripts/pars/program/vcf.extract.pl --file ../1011Matrix.wild.tsv --output 1011Matrix.ext.tsv

Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a cds
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a utr
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a syn
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a nsy

perl ~/Scripts/pars/program/vcf.merge.pro.pl --file cds_snp.merge.tsv --output ${NAME}.wild.cds_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file utr_snp.merge.tsv --output ${NAME}.wild.utr_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file syn_snp.merge.tsv --output ${NAME}.wild.syn_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file nsy_snp.merge.tsv --output ${NAME}.wild.nsy_snp.merge.pro.tsv

cat ${NAME}.wild.cds_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.cds_snp.merge.pro.txt
cat ${NAME}.wild.utr_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.utr_snp.merge.pro.txt
cat ${NAME}.wild.syn_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.syn_snp.merge.pro.txt
cat ${NAME}.wild.nsy_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.nsy_snp.merge.pro.txt
unset NAME

export NAME=Scer_n7p_Spar
mkdir -p ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}.wild
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}.wild

perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds.pars.tsv
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr.pars.tsv
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_syn.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_syn.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file data_SNPs_PARS_syn.csv.bak --output data_SNPs_PARS_syn.pars.tsv
rm -rf data_SNPs_PARS_syn.csv.bak
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_nsy.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_nsy.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_7.pl --file data_SNPs_PARS_nsy.csv.bak --output data_SNPs_PARS_nsy.pars.tsv
rm -rf data_SNPs_PARS_nsy.csv.bak

perl ~/Scripts/pars/program/vcf.extract.pl --file ../1011Matrix.wild.tsv --output 1011Matrix.ext.tsv

Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a cds
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a utr
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a syn
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a nsy

perl ~/Scripts/pars/program/vcf.merge.pro.pl --file cds_snp.merge.tsv --output ${NAME}.wild.cds_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file utr_snp.merge.tsv --output ${NAME}.wild.utr_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file syn_snp.merge.tsv --output ${NAME}.wild.syn_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file nsy_snp.merge.tsv --output ${NAME}.wild.nsy_snp.merge.pro.tsv

cat ${NAME}.wild.cds_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.cds_snp.merge.pro.txt
cat ${NAME}.wild.utr_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.utr_snp.merge.pro.txt
cat ${NAME}.wild.syn_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.syn_snp.merge.pro.txt
cat ${NAME}.wild.nsy_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.nsy_snp.merge.pro.txt
unset NAME

export NAME=Scer_n128_Spar
mkdir -p ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}.wild
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}.wild

perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds.pars.tsv
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr.pars.tsv
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_syn.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_syn.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file data_SNPs_PARS_syn.csv.bak --output data_SNPs_PARS_syn.pars.tsv
rm -rf data_SNPs_PARS_syn.csv.bak
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_nsy.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_nsy.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file data_SNPs_PARS_nsy.csv.bak --output data_SNPs_PARS_nsy.pars.tsv
rm -rf data_SNPs_PARS_nsy.csv.bak

perl ~/Scripts/pars/program/vcf.extract.pl --file ../1011Matrix.wild.tsv --output 1011Matrix.ext.tsv

Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a cds
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a utr
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a syn
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a nsy

perl ~/Scripts/pars/program/vcf.merge.pro.pl --file cds_snp.merge.tsv --output ${NAME}.wild.cds_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file utr_snp.merge.tsv --output ${NAME}.wild.utr_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file syn_snp.merge.tsv --output ${NAME}.wild.syn_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file nsy_snp.merge.tsv --output ${NAME}.wild.nsy_snp.merge.pro.tsv

cat ${NAME}.wild.cds_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.cds_snp.merge.pro.txt
cat ${NAME}.wild.utr_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.utr_snp.merge.pro.txt
cat ${NAME}.wild.syn_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.syn_snp.merge.pro.txt
cat ${NAME}.wild.nsy_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.nsy_snp.merge.pro.txt
unset NAME

export NAME=Scer_n128_Seub
mkdir -p ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}.wild
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}.wild

perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds.pars.tsv
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr.pars.tsv
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_syn.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_syn.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file data_SNPs_PARS_syn.csv.bak --output data_SNPs_PARS_syn.pars.tsv
rm -rf data_SNPs_PARS_syn.csv.bak
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_nsy.csv | perl -p -e 's/^(.+?),(.+?),/$2,$1,/;' > data_SNPs_PARS_nsy.csv.bak
perl ~/Scripts/pars/program/vcf.merge.pre_128.pl --file data_SNPs_PARS_nsy.csv.bak --output data_SNPs_PARS_nsy.pars.tsv
rm -rf data_SNPs_PARS_nsy.csv.bak

perl ~/Scripts/pars/program/vcf.extract.pl --file ../1011Matrix.wild.tsv --output 1011Matrix.ext.tsv

Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a cds
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a utr
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a syn
Rscript ~/Scripts/pars/program/vcf.merge.R -n ${NAME}.wild -a nsy

perl ~/Scripts/pars/program/vcf.merge.pro.pl --file cds_snp.merge.tsv --output ${NAME}.wild.cds_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file utr_snp.merge.tsv --output ${NAME}.wild.utr_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file syn_snp.merge.tsv --output ${NAME}.wild.syn_snp.merge.pro.tsv
perl ~/Scripts/pars/program/vcf.merge.pro.pl --file nsy_snp.merge.tsv --output ${NAME}.wild.nsy_snp.merge.pro.tsv

cat ${NAME}.wild.cds_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.cds_snp.merge.pro.txt
cat ${NAME}.wild.utr_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.utr_snp.merge.pro.txt
cat ${NAME}.wild.syn_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.syn_snp.merge.pro.txt
cat ${NAME}.wild.nsy_snp.merge.pro.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.wild.nsy_snp.merge.pro.txt
unset NAME

```

## update

```bash
export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/result/${NAME}
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a cds
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a utr
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a syn
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a nsy
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a cds -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a utr -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a syn -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a nsy -o .wild
cat data_SNPs_PARS_cds.update.csv data_SNPs_PARS_utr.update.csv | sort | uniq | perl -e 'print reverse <>' > data_SNPs_PARS_mRNA.update.csv
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/result/${NAME}
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a cds
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a utr
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a syn
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a nsy
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a cds -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a utr -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a syn -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a nsy -o .wild
cat data_SNPs_PARS_cds.update.csv data_SNPs_PARS_utr.update.csv | sort | uniq | perl -e 'print reverse <>' > data_SNPs_PARS_mRNA.update.csv
unset NAME

export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME}
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a cds
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a utr
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a syn
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a nsy
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a cds -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a utr -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a syn -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a nsy -o .wild
cat data_SNPs_PARS_cds.update.csv data_SNPs_PARS_utr.update.csv | sort | uniq | perl -e 'print reverse <>' > data_SNPs_PARS_mRNA.update.csv
unset NAME

export NAME=Scer_n128_Seub
cd ~/data/mrna-structure/result/${NAME}
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a cds
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a utr
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a syn
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a nsy
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a cds -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a utr -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a syn -o .wild
Rscript ~/Scripts/pars/program/update_SNPs.R -n ${NAME} -a nsy -o .wild
cat data_SNPs_PARS_cds.update.csv data_SNPs_PARS_utr.update.csv | sort | uniq | perl -e 'print reverse <>' > data_SNPs_PARS_mRNA.update.csv
unset NAME
```

## count A/T <-> G/C

```bash

export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/result/${NAME}
mkdir -p ~/data/mrna-structure/result/${NAME}/freq_each
Rscript ~/Scripts/pars/program/count_AT_GC.R -n ${NAME} 
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_cds_stat.csv --output freq_each/PARS_cds_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_utr_stat.csv --output freq_each/PARS_utr_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_syn_stat.csv --output freq_each/PARS_syn_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_nsy_stat.csv --output freq_each/PARS_nsy_stat_chi_square.csv
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/result/${NAME}
mkdir -p ~/data/mrna-structure/result/${NAME}/freq_each
Rscript ~/Scripts/pars/program/count_AT_GC.R -n ${NAME} 
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_cds_stat.csv --output freq_each/PARS_cds_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_utr_stat.csv --output freq_each/PARS_utr_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_syn_stat.csv --output freq_each/PARS_syn_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_nsy_stat.csv --output freq_each/PARS_nsy_stat_chi_square.csv
unset NAME

export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME}
mkdir -p ~/data/mrna-structure/result/${NAME}/freq_each
mkdir -p ~/data/mrna-structure/result/${NAME}/freq_10
Rscript ~/Scripts/pars/program/count_AT_GC.R -n ${NAME} 
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_cds_stat.csv --output freq_each/PARS_cds_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_utr_stat.csv --output freq_each/PARS_utr_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_syn_stat.csv --output freq_each/PARS_syn_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_nsy_stat.csv --output freq_each/PARS_nsy_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_mRNA_stat.csv --output freq_each/PARS_mRNA_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_cds_stat_freq_10.csv --output freq_10/PARS_cds_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_utr_stat_freq_10.csv --output freq_10/PARS_utr_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_syn_stat_freq_10.csv --output freq_10/PARS_syn_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_nsy_stat_freq_10.csv --output freq_10/PARS_nsy_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_mRNA_stat_freq_10.csv --output freq_10/PARS_mRNA_stat_freq_10_chi_square.csv
unset NAME

export NAME=Scer_n128_Seub
cd ~/data/mrna-structure/result/${NAME}
mkdir -p ~/data/mrna-structure/result/${NAME}/freq_each
mkdir -p ~/data/mrna-structure/result/${NAME}/freq_10
Rscript ~/Scripts/pars/program/count_AT_GC.R -n ${NAME} 
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_cds_stat.csv --output freq_each/PARS_cds_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_utr_stat.csv --output freq_each/PARS_utr_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_syn_stat.csv --output freq_each/PARS_syn_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_nsy_stat.csv --output freq_each/PARS_nsy_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_cds_stat_freq_10.csv --output freq_10/PARS_cds_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_utr_stat_freq_10.csv --output freq_10/PARS_utr_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_syn_stat_freq_10.csv --output freq_10/PARS_syn_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_nsy_stat_freq_10.csv --output freq_10/PARS_nsy_stat_freq_10_chi_square.csv
unset NAME

```

## count stem length selection

```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME} 
mkdir -p freq_10/stem_length

perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --origin data_SNPs_PARS_cds.update.csv --output data_SNPs_PARS_cds.update_pos.csv
perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --origin data_SNPs_PARS_utr.update.csv --output data_SNPs_PARS_utr.update_pos.csv

cat data_SNPs_PARS_cds.update_pos.csv data_SNPs_PARS_utr.update_pos.csv | sort | uniq | perl -e 'print reverse <>' > data_SNPs_PARS_mRNA.update_pos.csv 

Rscript ~/Scripts/pars/program/count_AT_GC_gene_trait.R -n ${NAME}

cat data_SNPs_PARS_mRNA.update_pos.csv | perl -nl -a -F"," -e 'print qq{$F[1]};' | sort | uniq | perl -e 'print reverse <>' > mRNA.gene.list.csv

perl ~/Scripts/pars/program/count_structure_length_gene.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --name ~/data/mrna-structure/result/${NAME}/mRNA.gene.list.csv --structure stem --output stem_length_mRNA.csv
perl ~/Scripts/pars/program/count_structure_length_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --name ~/data/mrna-structure/result/${NAME}/mRNA.gene.list.csv --structure loop --output loop_length_mRNA.csv

#perl ~/Scripts/pars/program/count_structure_length_gene.update.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.update.yml --name ~/data/mrna-structure/result/${NAME}/mRNA.gene.list.csv --structure stem --output stem_length_cds.update.csv
#perl ~/Scripts/pars/program/count_structure_length_gene.update.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.update.yml --name ~/data/mrna-structure/result/${NAME}/mRNA.gene.list.csv --structure loop --output loop_length_cds.update.csv

unset NAME

```

## count_codon_gene

```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME} 
perl ~/Scripts/pars/program/count_codon_gene.pl --origin data_SNPs_PARS_cds.update.csv --output data_SNPs_PARS_cds.update_codon.csv
perl ~/Scripts/pars/program/count_codon_gene.pl --origin data_SNPs_PARS_syn.update.csv --output data_SNPs_PARS_syn.update_codon.csv
Rscript ~/Scripts/pars/program/count_AT_GC_codon.R -n ${NAME}
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_tRNA_stat.csv --output freq_each/PARS_tRNA_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_4D_stat.csv --output freq_each/PARS_4D_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_tRNA_stat_freq_10.csv --output freq_10/PARS_tRNA_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_4D_stat_freq_10.csv --output freq_10/PARS_4D_stat_freq_10_chi_square.csv
unset NAME

```

## count per gene cds_utr

```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME} 
perl ~/Scripts/pars/program/count_cut_range.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --cut ~/data/mrna-structure/process/sce_cds.yml --output stem_loop_cds_length.csv 
perl ~/Scripts/pars/program/count_cut_range.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --cut ~/data/mrna-structure/process/sce_utr.yml --output stem_loop_utr_length.csv
perl ~/Scripts/pars/program/count_per_gene_ACGT_percent.pl --file data_SNPs_PARS_cds.update.csv --output data_SNPs_PARS_cds.update_per_gene_ATGC.csv
perl ~/Scripts/pars/program/count_per_gene_ACGT_percent.pl --file data_SNPs_PARS_utr.update.csv --output data_SNPs_PARS_utr.update_per_gene_ATGC.csv
Rscript ~/Scripts/pars/program/count_cds_utr.R -n ${NAME}
unset NAME
```

## count GO KEGG

```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME}

#筛选snp（cds stem中A/T->G/C）
cat ~/data/mrna-structure/vcf/1011Matrix.gvcf/${NAME}/${NAME}.cds_snp.merge.pro.tsv \
    | perl -nla -F"\t" -e '
        if ( ( $F[2] eq "stem" ) && ( ($F[3] eq "A->G") || ($F[3] eq "A->C") || ($F[3] eq "T->C") || ($F[3] eq "T->G") ) ){
            print qq{$F[1]};
        }
    ' \
    | sort | uniq > ${NAME}.cds.filtrate.txt

#将filtrate后的gene list输入 https://david.ncifcrf.gov/ 中，得到GO，KEGG信息
mkdir -p freq_10/GO
mkdir -p freq_10/KEGG
Rscript ~/Scripts/pars/program/count_AT_GC_GO.R -n ${NAME}
Rscript ~/Scripts/pars/program/count_AT_GC_KEGG.R -n ${NAME}

#得到Scer_n128_Spar_go_kegg.csv

mkdir -p freq_10/go_kegg
mkdir -p freq_10/go_kegg/syn
mkdir -p freq_10/go_kegg/nsy
Rscript ~/Scripts/pars/program/count_AT_GC_go_kegg.R -n ${NAME}

unset NAME

```

## stat subpopulation SNPs frequency

```bash
export NAME=Scer_n128_Spar
mkdir -p ~/data/mrna-structure/result/${NAME}/subpop
cd ~/data/mrna-structure/result/${NAME}/subpop

#generate strainlist, order same as egaz template
echo -e "S288c\nbeer001\nbeer003\nbeer004\nbeer005\nbeer006\nbeer007\nbeer008\nbeer009\nbeer010\nbeer011\nbeer012\nbeer013\nbeer014\nbeer015\nbeer016\nbeer020\nbeer021\nbeer022\nbeer023\nbeer024\nbeer025\nbeer026\nbeer027\nbeer028\nbeer029\nbeer030\nbeer031\nbeer032\nbeer033\nbeer034\nbeer036\nbeer037\nbeer038\nbeer040\nbeer041\nbeer043\nbeer044\nbeer045\nbeer046\nbeer047\nbeer048\nbeer049\nbeer050\nbeer051\nbeer052\nbeer053\nbeer054\nbeer055\nbeer056\nbeer059\nbeer061\nbeer062\nbeer063\nbeer064\nbeer065\nbeer066\nbeer067\nbeer068\nbeer069\nbeer070\nbeer071\nbeer073\nbeer075\nbeer076\nbeer077\nbeer078\nbeer079\nbeer080\nbeer081\nbeer082\nbeer083\nbeer084\nbeer085\nbeer086\nbeer087\nbeer088\nbeer089\nbeer090\nbeer091\nbeer092\nbeer094\nbeer095\nbeer096\nbeer097\nbeer098\nbeer099\nbeer100\nbeer101\nbeer102\nbioethanol001\nbioethanol003\nbioethanol004\nbread001\nbread002\nbread003\nbread004\nsake001\nsake002\nsake003\nsake004\nsake005\nsake006\nsake007\nspirits001\nspirits002\nspirits003\nspirits004\nspirits005\nspirits011\nwine001\nwine003\nwine004\nwine005\nwine006\nwine007\nwine009\nwine010\nwine011\nwine012\nwine013\nwine014\nwine015\nwine017\nwine018\nwild005\nwild006\nwild007" > strainlist.csv
 

#download total SNPs from MySQL → total_snp.csv
cat total_snp.csv | wc -l #check number

#get genelist by filitering strong selection from GO/KEGG annotation and deleting repeating item

echo -e "gene\nYDL174C\nYKR066C\nYNR001C\nYDR487C\nYOR065W\nYCL057W\nYIR037W\nYMR203W\nYPL132W\nYLR304C\nYJL054W\nYNL055C\nYJR104C\nYKR071C\nYKL150W\nYBR056W\nYDR226W\nYDR375C\nYKL053C-A\nYKL087C\nYGL187C\nYLR259C\nYFR033C\nYJR121W\nYMR145C\nYAL039C\nYKL067W\nYDL120W\nYDR353W\nYLL009C\nYDR511W\nYPR140W\nYJL143W\nYHR116W\nYNR022C\nYDR296W\nYHR147C\nYBR268W\nYDR116C\nYKR006C\nYLR439W\nYPL183W-A\nYMR286W\nYDL202W\nYKL138C\nYJL063C\nYMR024W\nYNL005C\nYBR122C\nYJL096W\nYOR150W\nYDR237W\nYBL038W\nYLR312W-A\nYKL167C\nYKR085C\nYLR008C\nYKL016C\nYDR204W\nYMR241W\nYKL148C\nYGR096W\nYJL166W\nYOR065W\nYOR297C\nYMR089C\nYGR257C\nYPL134C\nYCL044C\nYPR024W\nYJR045C\nYKR087C\nYGL187C\nYBR039W\nYEL052W\nYGR222W\nYER017C\nYKR065C\nYLR259C\nYFR033C\nYMR035W\nYBR185C\nYMR256C\nYOR176W\nYGR033C\nYPR140W\nYML110C\nYKL141W\nYLR203C\nYHR024C\nYBR085W\nYDL174C\nYMR301C\nYER078C\nYML120C\nYPL270W\nYLR253W\nYLL041C\nYLR164W\nYNL003C\nYER141W\nYPR191W\nYGR235C\nYPL132W\nYJL054W\nYBL030C\nYGR062C\nYOL008W\nYPL063W\nYDR298C\nYLR188W\nYDR236C\nYDR375C\nYKR052C\nYKL087C\nYEL024W\nYGR101W\nYHR037W\nYDR377W\nYPL078C\nYPL271W\nYJR121W\nYOR232W\nYOR356W\nYBR291C\nYNL100W\nYAL039C\nYBR003W\nYDL120W\nYDL004W\nYPL189C-A\nYOR125C\nYCL057C-A\nYKL120W\nYIL134W\nYIL022W\nYOR222W\nYJL143W\nYOL027C\nYGR082W\nYNL026W\nYLR099W-A\nYLR090W\nYML086C\nYNL055C\nYKL150W\nYNL131W\nYMR110C\nYNL121C\nYER019W\nYPR140W\nYNL070W\nYHR117W\nYLR008C\nYGR082W\nYER017C\nYKR065C\nYOR232W\nYMR203W\nYLR090W\nYPL063W\nYPR024W\nYNL131W\nYJR045C\nYNL121C\nYGR033C\nYIL022W\nYNL070W\nYJL143W\nYHR024C\nYLR008C\nYDL174C\nYNR001C\nYOR065W\nYMR203W\nYOR297C\nYLR304C\nYJL054W\nYNL055C\nYPL063W\nYDR375C\nYKL053C-A\nYGL187C\nYNL121C\nYHR117W\nYNL070W\nYGR082W\nYKR065C\nYNL026W\nYLR259C\nYOR027W\nYJR121W\nYOR232W\nYHR083W\nYGR028W\nYNL064C\nYDL120W\nYNL131W\nYLL009C\nYPR140W\nYGR033C\nYIL022W\nYJL143W\nYMR301C\nYIL003W\nYPL270W\nYMR312W\nYCR011C\nYHR169W\nYGL008C\nYMR089C\nYDL007W\nYPR173C\nYLR397C\nYOR259C\nYDR091C\nYLR188W\nYDL166C\nYNL329C\nYHR187W\nYJR045C\nYGR262C\nYDR375C\nYBR039W\nYGR210C\nYBL022C\nYDL126C\nYFL028C\nYER017C\nYER036C\nYDR377W\nYLR259C\nYDR061W\nYPL271W\nYJR121W\nYJR072C\nYNL290W\nYCL047C\nYOR291W\nYEL031W\nYOR117W\nYGR028W\nYDL100C\nYLR249W\nYDL004W\nYPL226W\nYOR153W\nYKL073W\nYGL048C\nYBR080C\nYDR011W\nYOR278W\nYAL039C\nYDR047W\nYDL120W\nYDR232W\nYKL087C\nYER141W\nYGL040C\nYDR044W\nYGL245W\nYPL172C\nYOR176W\nYDL174C\nYDL078C\nYPL061W\nYDR272W\nYML004C\nYBL015W\nYBR218C\nYOR374W\nYOR347C\nYPL262W\nYKL085W\nYPL028W\nYER178W\nYMR110C\nYLR153C\nYNL071W\nYBR221C\nYNR016C" | sort | uniq | perl -e 'print reverse <>' > genelist.csv

#filiter SNPs
RScript ~/Scripts/pars/program/subpop.R -n ${NAME} -i genelist.csv -o filiter_snp.csv

#delete "double quotation marks" and "blank" in filiter_snp.csv

#calculate subpopulation SNPs proporation
perl ~/Scripts/pars/program/subpop.pl filiter_snp.csv strainlist.csv > subpop.csv

RScript ~/Scripts/pars/program/subpop_merge.R -n ${NAME}

unset NAME

```

