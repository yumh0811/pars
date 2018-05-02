# Processing Yeast PARS Data

[TOC level=1-3]: # " "
- [Processing Yeast PARS Data](#processing-yeast-pars-data)
- [Download reference data](#download-reference-data)
    - [Download PARS10 full site.](#download-pars10-full-site)
    - [Download S288c annotation data from ensembl by rsync](#download-s288c-annotation-data-from-ensembl-by-rsync)
    - [SGD](#sgd)
- [Download strains and outgroups](#download-strains-and-outgroups)
    - [Sanger (NCBI WGS)](#sanger-ncbi-wgs)
    - [Illumina (NCBI ASSEMBLY)](#illumina-ncbi-assembly)
- [Plans of alignments](#plans-of-alignments)
- [AlignDB](#aligndb)
    - [Build alignDB for multiple genomes](#build-aligndb-for-multiple-genomes-n7)
    - [Extract gene-list and snp-codon-list](#extract-gene-list-and-snp-codon-list-n7)
    - [SNPs and indels](#snps-and-indels-n7)
- [Blast](#blast)
- [Features](#features)
- [Real Processing](#real-processing)
- [Download other reference data](#download-other-reference-data)
    - [mRNA levels](#mrna-levels)
    - [ess, rich/minimal and chem](#ess-richminimal-and-chem)
    - [Recombination rates](#recombination-rates)
    - [Protein-protein interactions](#protein-protein-interactions)


# Download reference data

## Download PARS10 full site.

```bash
mkdir -p ~/data/mrna-structure/PARS10
cd ~/data/mrna-structure/PARS10

perl ~/Scripts/download/list.pl -u http://genie.weizmann.ac.il/pubs/PARS10/
perl ~/Scripts/download/download.pl -i pubs_PARS10.yml

find . -name "*.gz" | xargs gzip -d
```

## Download S288c annotation data from ensembl by rsync

http://www.ensembl.org/info/data/ftp/rsync.html?redirect=no

S288c assembly version is not changed since 2011, R64-1-1
(GCA_000146045.2).

```bash
mkdir -p ~/data/mrna-structure/ensembl82/mysql
cd ~/data/mrna-structure/ensembl82/mysql
rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-82/mysql/saccharomyces_cerevisiae_core_82_4 .

mkdir -p ~/data/mrna-structure/ensembl82/fasta
cd ~/data/mrna-structure/ensembl82/fasta
rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-82/fasta/saccharomyces_cerevisiae .

perl ~/Scripts/withncbi/ensembl/build_ensembl.pl -e ~/data/mrna-structure/ensembl82/mysql/saccharomyces_cerevisiae_core_82_4 --checksum
perl ~/Scripts/withncbi/ensembl/build_ensembl.pl -e ~/data/mrna-structure/ensembl82/mysql/saccharomyces_cerevisiae_core_82_4 --initdb --db saccharomyces_cerevisiae_core_29_82_4
```

## SGD

```bash
mkdir -p ~/data/mrna-structure/sgd
cd ~/data/mrna-structure/sgd

aria2c -c http://downloads.yeastgenome.org/sequence/S288C_reference/intergenic/NotFeature.fasta.gz
aria2c -c http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz
aria2c -c http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_genomic_all.fasta.gz
aria2c -c http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff

find . -name "*.gz" \
    | parallel -j 1 "
        echo {};
        gzip -d -c {} > {.};
    "
```


# Download strains and outgroups

## Sanger (NCBI WGS)

```bash
mkdir -p ~/data/mrna-structure/GENOMES
cd ~/data/mrna-structure/GENOMES

perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/pars/scer_wgs.tsv \
    --fix -a \
    -o WGS

aria2c -UWget -x 6 -s 3 -c -i WGS/scer_wgs.url.txt

find WGS -name "*.gz" | xargs gzip -t
```

## Illumina (NCBI ASSEMBLY)

```bash
mkdir -p ~/data/mrna-structure/GENOMES/ASSEMBLIES
cd ~/data/mrna-structure/GENOMES/ASSEMBLIES

# Download, rename files and change fasta headers
perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
    -p -f ~/Scripts/pars/scer_assembly.csv

```

# Plans of alignments

```bash
# create downloaded genome list
cat ~/Scripts/pars/scer_assembly.csv \
    | grep -v "^#" \
    | cut -d',' -f1,3 \
    | uniq \
    | perl -nl -a -F"," -e 'printf qq{    --download "name=%s;taxon=%s" \\\n}, $F[0], $F[1];'

mkdir -p ~/data/mrna-structure/alignment/scer_wgs
cd ~/data/mrna-structure/alignment/scer_wgs

perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
    -i ~/data/mrna-structure/GENOMES/WGS/scer_wgs.data.yml \
    -o scer_wgs.plan.yml \
    -d ~/data/mrna-structure/GENOMES/WGS \
    -m prefix \
    -r '*.fsa_nt.gz' \
    --opt group_name=scer_wgs \
    --opt base_dir='~/data/mrna-structure/alignment' \
    --opt data_dir="~/data/mrna-structure/alignment/scer_wgs" \
    --opt rm_species=Fungi \
    --dd ~/data/mrna-structure/GENOMES/ASSEMBLIES \
    --download "name=S288c;taxon=559292" \
    --download "name=DBVPG6044;taxon=100000000001" \
    --download "name=UWOPS03_461_4;taxon=100000000002" \
    --download "name=Y12;taxon=100000000003" \
    --download "name=SK1;taxon=100000000004" \
    --download "name=YPS128;taxon=100000000005" \
    --download "name=DBVPG6765;taxon=100000000006" \
    --download "name=EC1118;taxon=643680" \
    --plan 'name=Scer_n7_Spar;t=S288c;qs=EC1118,Kyokai_no_7,RM11_1a,Sigma1278b,T7,YJM789,Spar;o=Spar' \
    --plan 'name=Scer_n7p_Spar;t=S288c;qs=DBVPG6044,UWOPS03_461_4,Y12,SK1,YPS128,DBVPG6765,Spar;o=Spar' \
    --plan 'name=Scer_n157_Spar;t=S288c;qs=beer001,beer002,beer003,beer004,beer005,beer006,beer007,beer008,beer009,beer010,beer011,beer012,beer013,beer014,beer015,beer016,beer017,beer018,beer019,beer020,beer021,beer022,beer023,beer024,beer025,beer026,beer027,beer028,beer029,beer030,beer031,beer032,beer033,beer034,beer035,beer036,beer037,beer038,beer039,beer040,beer041,beer042,beer043,beer044,beer045,beer046,beer047,beer048,beer049,beer050,beer051,beer052,beer053,beer054,beer055,beer056,beer057,beer058,beer059,beer060,beer061,beer062,beer063,beer064,beer065,beer066,beer067,beer068,beer069,beer070,beer071,beer072,beer073,beer074,beer075,beer076,beer077,beer078,beer079,beer080,beer081,beer082,beer083,beer084,beer085,beer086,beer087,beer088,beer089,beer090,beer091,beer092,beer093,beer094,beer095,beer096,beer097,beer098,beer099,beer100,beer101,beer102,bioethanol001,bioethanol002,bioethanol003,bioethanol004,bioethanol005,bread001,bread002,bread003,bread004,laboratory001,laboratory002,sake001,sake002,sake003,sake004,sake005,sake006,sake007,spirits001,spirits002,spirits003,spirits004,spirits005,spirits006,spirits007,spirits008,spirits009,spirits010,spirits011,wine001,wine002,wine003,wine004,wine005,wine006,wine007,wine008,wine009,wine010,wine011,wine012,wine013,wine014,wine015,wine016,wine017,wine018,wine019,wild001,wild002,wild003,wild004,wild005,wild006,wild007,Spar;o=Spar' \
    --plan 'name=Scer_n157_nonMosaic_Spar;t=S288c;qs=beer001,beer002,beer003,beer004,beer005,beer006,beer007,beer008,beer009,beer010,beer011,beer012,beer013,beer014,beer015,beer016,beer020,beer021,beer022,beer023,beer024,beer025,beer026,beer027,beer028,beer029,beer030,beer031,beer032,beer033,beer034,beer036,beer037,beer038,beer039,beer040,beer041,beer043,beer044,beer045,beer046,beer047,beer048,beer049,beer050,beer051,beer052,beer053,beer054,beer055,beer056,beer059,beer061,beer062,beer063,beer064,beer065,beer066,beer067,beer068,beer069,beer070,beer071,beer073,beer075,beer076,beer077,beer078,beer079,beer080,beer081,beer082,beer083,beer084,beer085,beer086,beer087,beer088,beer089,beer090,beer091,beer092,beer094,beer095,beer096,beer097,beer098,beer099,beer100,beer101,beer102,bioethanol001,bioethanol003,bioethanol004,bread001,bread002,bread003,bread004,sake001,sake002,sake003,sake004,sake005,sake006,sake007,spirits001,spirits002,spirits003,spirits004,spirits005,spirits011,wine001,wine003,wine004,wine005,wine006,wine007,wine009,wine010,wine011,wine012,wine013,wine014,wine015,wine017,wine018,wild004,wild005,wild006,wild007,Spar;o=Spar' \
    -y

# pop_prep.pl
perl ~/Scripts/withncbi/pop/pop_prep.pl -p 16 -i scer_wgs.plan.yml

bash 01_file.sh
bash 02_rm.sh
bash 03_strain_info.sh

# plan_ALL.sh
bash plan_ALL.sh

bash 1_real_chr.sh
bash 3_pair_cmd.sh
bash 4_rawphylo.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

# other plans
bash plan_Scer_n7_Spar.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

# other plans
bash plan_Scer_n7p_Spar.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

# other plans
bash plan_Scer_n157_Spar.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

# other plans
bash plan_Scer_n157_nonMosaic_Spar.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh
```

# AlignDB

## Build alignDB for multiple genomes n7

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n7_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n7_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n7_Spar -r 1-60

```

## Build alignDB for multiple genomes n7p

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n7p_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n7p_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n7p_Spar -r 1-60

```

## Build alignDB for multiple genomes n157

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n157_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_Spar -r 1-60

```

## Build alignDB for multiple genomes n157_nonMosaic

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_nonMosaic_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n157_nonMosaic_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_nonMosaic_Spar -r 1-60

```

## Extract gene-list and snp-codon-list n7

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n7_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n7_Spar.mvar.gene_list.csv
```

## Extract gene-list and snp-codon-list n7p

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n7p_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n7p_Spar.mvar.gene_list.csv
```

## Extract gene-list and snp-codon-list n157

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_Spar.mvar.gene_list.csv
```

## Extract gene-list and snp-codon-list n157_nonMosaic

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_nonMosaic_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_nonMosaic_Spar.mvar.gene_list.csv
```

## SNPs and indels n7

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n7_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n7_Spar.indel.pos.txt

```

## SNPs and indels n7p

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n7p_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n7p_Spar.indel.pos.txt

```

## SNPs and indels n157

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_Spar.indel.pos.txt

```

## SNPs and indels n157_nonMosaic

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_nonMosaic_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_nonMosaic_Spar.indel.pos.txt

```

# Blast

Prepare a combined fasta file of yeast genome and blast genes against
the genome.

```bash
mkdir -p ~/data/mrna-structure/blast
cd ~/data/mrna-structure/blast

cat ~/data/mrna-structure/alignment/scer_wgs/Genomes/S288c/{I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI,Mito}.fa \
    > S288c.fa

perl -nl -i -e '/^>/ or $_ = uc $_; print'  S288c.fa
faops size S288c.fa > S288c.sizes

# formatdb
~/share/blast/bin/formatdb -p F -o T -i S288c.fa

# blast every transcripts against genome
~/share/blast/bin/blastall -p blastn -F "m D" -m 0 -b 10 -v 10 -e 1e-3 -a 4 \
    -i ../PARS10/pubs/PARS10/data/sce_genes.fasta -d S288C.fa -o sce_genes.blast
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
| Name | chrLength | size | coverage |
|:--|--:|--:|--:|
| genes | 12071326 | 4235405 | 0.3509 |
| intergenic | 12071326 | 2864170 | 0.2373 |
| intron | 12071326 | 65144 | 0.0054 |
| orf_genomic | 12071326 | 8895737 | 0.7369 |
| utr | 12071326 | 516569 | 0.0428 |
| mRNA | 12071326 | 4233361 | 0.3507 |
| cds | 12071326 | 3716792 | 0.3079 |


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

# SNPs within utr regions
runlist position --op superset \
    sce_utr.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within orf_genomic
runlist position --op superset \
    sce_orf_genomic.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

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

# SNPs within utr regions
runlist position --op superset \
    sce_utr.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within orf_genomic
runlist position --op superset \
    sce_orf_genomic.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

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

# Real Processing n157

```bash
export NAME=Scer_n157_Spar

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

# SNPs within utr regions
runlist position --op superset \
    sce_utr.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within orf_genomic
runlist position --op superset \
    sce_orf_genomic.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

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

# Real Processing n157_nonMosaic

```bash
export NAME=Scer_n157_nonMosaic_Spar

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

# SNPs within utr regions
runlist position --op superset \
    sce_utr.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within orf_genomic
runlist position --op superset \
    sce_orf_genomic.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

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

# Download other reference data

## mRNA levels

* Data source: Quantification of the yeast transcriptome by
  single-molecule sequencing. Lipson, D. et al. Nature Biotechnology 27,
  652-658 (2009) doi:10.1038/nbt.1551

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://www.nature.com/nbt/journal/v27/n7/extref/nbt.1551-S2.xls

perl ~/Scripts/fig_table/xlsx2csv.pl -f nbt.1551-S2.xls --sheet 'counts' \
    | perl -nla -F/,/ -e '
    if ( /^#/ ) {
        print qq{#ORF\tGene\tAvg};
    }
    elsif ( /^\d/ ) {
        print qq{$F[1]\t$F[2]\t$F[9]};
    }
    ' \
    > mrna_levels.tsv
```

## ess, rich/minimal and chem

* ess: Giaever, G., et al. Functional Profiling of theSaccharomyces
  cerevisiae Genome. Nature 418, 387-391. (2002)
* rich/minimal: Mechanisms of Haploinsufficiency Revealed by Genome-Wide
  Profiling in Yeast Deutschbauer, AM. et al. GENETICS April 1, 2005
  vol. 169 no. 4 1915-1925; 10.1534/genetics.104.036871
* chem: The Chemical Genomic Portrait of Yeast: Uncovering a Phenotype
  for All Genes. Hillenmeyer, M.E. et al. Science 18 Apr 2008: Vol. 320,
  Issue 5874, pp. 362-365 DOI: 10.1126/science.1150021

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt

cat Essential_ORFs.txt \
    | perl -nl -e '
        next unless /^\d/;
        my $orf = ( split /\s+/ )[1];
        print qq{$orf};
    ' \
    > ess_orf.tsv

wget -N http://www-sequence.stanford.edu/group/research/HIP_HOP/supplements/01yfh/files/OrfGeneData.txt

cat OrfGeneData.txt \
    | perl -nla -F"\t" -e '
        printf q{#} if /^orf/;
        print qq{$F[0]\t$F[1]\t$F[5]\t$F[13]\t$F[17]};
    ' \
    > rich_orf.tsv

# http://chemogenomics.stanford.edu/supplements/global/

wget -N http://chemogenomics.stanford.edu/supplements/global/download/data/hom.z_tdist_pval_nm.counts.smallmol.cutoff.01.xls

perl ~/Scripts/fig_table/xlsx2csv.pl -f hom.z_tdist_pval_nm.counts.smallmol.cutoff.01.xls --sheet 'hom.z_tdist_pval_nm.smallmol.co' \
    | perl -nla -F"," -MText::CSV_XS -e '
    BEGIN { 
        our $csv = Text::CSV_XS->new(); 
        print qq{#ORF\tGene\tCount};
    }
    
    if ( /^Y/ ) {
        if ($csv->parse($_)) {
            my @fields = $csv->fields();
            print qq{$fields[0]\t$fields[1]\t$fields[6]};
        }
    }
    ' \
    > chem_orf.tsv

```

## Recombination rates

* Data source: Global mapping of meiotic recombination hotspots and
  coldspots in the yeast Saccharomyces cerevisiae. vol. 97 no. 21 PNAS
  Jennifer L. Gerton, 11383–11390

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://derisilab.ucsf.edu/data/hotspots/forWebORFs.txt

cat forWebORFs.txt \
    | perl -nla -F"\t" -MStatistics::Lite -e '
        next unless /^Y/;    # ORF stable id start with a "Y"
        next if @F < 2;
        my $rec_rate  = Statistics::Lite::median(grep {defined} @F[1 .. 7]);
        print qq{$F[0]\t$rec_rate};
    ' \
    > rec_rate.tsv
```

## Protein-protein interactions

* Data source:

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://drygin.ccbr.utoronto.ca/%7Ecostanzo2009/sgadata_costanzo2009_stringentCutoff_101120.txt.gz

gzip -d -c sgadata_costanzo2009_stringentCutoff_101120.txt.gz \
    | perl -nla -F"\t" -e '
        BEGIN { our %interact_of; }
        
        next unless /^Y/;    # ORF stable id start with a "Y"
        next if @F != 7;
        $interact_of{$F[0]}++;
        
        END {
            for my $key ( sort keys %interact_of ) {
                print qq{$key\t$interact_of{$key}};
            }
        }
    ' \
    > interact_count.tsv
```

# Phylogeny

## sgd/saccharomyces_cerevisiae.gff → protein coding gene list

```bash
mkdir -p ~/data/mrna-structure/phylogeny
cd ~/data/mrna-structure/phylogeny

perl ~/Scripts/pars/program/protein_coding_list.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list.csv

perl ~/Scripts/pars/program/protein_coding_list_range.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list_range.csv

perl ~/Scripts/pars/program/protein_coding_list_range_chr.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list_range_chr.csv

```

## cut cds alignment

### create cds_yml

```bash
cd ~/data/mrna-structure/phylogeny
mkdir -p ~/data/mrna-structure/phylogeny/gene_cds_yml

perl ~/Scripts/pars/program/cut_cds_yml.pl --file protein_coding_list_range_chr.csv --output gene_cds_yml

```

### cut cds_alignment by cds_yml (n157_nonMosaic)

```bash
export NAME=Scer_n157_nonMosaic_Spar

cp -rf ~/data/mrna-structure/alignment/scer_wgs/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas

mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds

cat ../protein_coding_list.csv |
   parallel --line-buffer -j 8 '
   	   fasops slice ../${NAME}_refined/species.fas ../gene_cds_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

```

### count cds_alignment proporation in sgd (n157_nonMosaic)

```bash
export NAME=Scer_n157_nonMosaic_Spar

cd ~/data/mrna-structure/phylogeny

perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_cds --output ${NAME}_gene_range.csv

unset NAME
```

## create gene_phylogeny (n157_nonMosaic)

```bash
export NAME=Scer_n157_nonMosaic_Spar

mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_newick
cd ~/data/mrna-structure/phylogeny/${NAME}_newick

cat ../protein_coding_list.csv |
  parallel --line-buffer -j 1 '
      if [ -e {}.nwk ]; then
          echo >&2 '    {}.nwk already presents  '
          exit;
      fi
      egaz raxml ../${NAME}_gene_alignment_cds/{}.fas.fas --seed 999 --tmp . --parallel 8 --outgroup Spar -o {}.nwk
  '
unset NAME
```

## count distance

```bash
export NAME=Scer_n157_nonMosaic_Spar

cd ~/data/mrna-structure/phylogeny
echo 'S288c,beer001,beer002,beer003,beer004,beer005,beer006,beer007,beer008,beer009,beer010,beer011,beer012,beer013,beer014,beer015,beer016,beer020,beer021,beer022,beer023,beer024,beer025,beer026,beer027,beer028,beer029,beer030,beer031,beer032,beer033,beer034,beer036,beer037,beer038,beer039,beer040,beer041,beer043,beer044,beer045,beer046,beer047,beer048,beer049,beer050,beer051,beer052,beer053,beer054,beer055,beer056,beer059,beer061,beer062,beer063,beer064,beer065,beer066,beer067,beer068,beer069,beer070,beer071,beer073,beer075,beer076,beer077,beer078,beer079,beer080,beer081,beer082,beer083,beer084,beer085,beer086,beer087,beer088,beer089,beer090,beer091,beer092,beer094,beer095,beer096,beer097,beer098,beer099,beer100,beer101,beer102,bioethanol001,bioethanol003,bioethanol004,bread001,bread002,bread003,bread004,sake001,sake002,sake003,sake004,sake005,sake006,sake007,spirits001,spirits002,spirits003,spirits004,spirits005,spirits011,wine001,wine003,wine004,wine005,wine006,wine007,wine009,wine010,wine011,wine012,wine013,wine014,wine015,wine017,wine018,wild004,wild005,wild006,wild007,Spar' \
| tr "," "\n" \
> ${NAME}_strain_name.list

mkdir ~/data/mrna-structure/phylogeny/${NAME}_distance
cat protein_coding_list.csv |
   parallel --line-buffer -j 8 '
    if [ -e "${NAME}_newick/{}.nwk" ]; then
       perl ~/Scripts/pars/program/count_distance.pl --file ~/data/mrna-structure/phylogeny/${NAME}_newick/{}.nwk --list ~/data/mrna-structure/phylogeny/${NAME}_strain_name.list --output ${NAME}_distance/{}.csv
    fi
    '

unset NAME
```
