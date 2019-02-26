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
- [Gene filter](#gene-filter)
    - [create protein coding gene list](#create-protein-coding-gene-list)
    - [Intact mRNAs](#intact-mrnas)
    - [Cut mRNA alignments and extract SNP list](#cut-mrna-alignments-and-extract-snp-list)
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

# clean
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw"   | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

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

# clean
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw"   | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

# Blast

Prepare a combined fasta file of yeast genome and blast genes against the genome.

```bash
mkdir -p ~/data/mrna-structure/blast
cd ~/data/mrna-structure/blast

cat ~/data/mrna-structure/GENOMES/S288c/{I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI,Mito}.fa \
    > S288c.fa

perl -nl -i -e '/^>/ or $_ = uc $_; print'  S288c.fa
faops size S288c.fa > S288c.sizes

# formatdb
makeblastdb -dbtype nucl -in S288c.fa -parse_seqids 

# blast every transcripts against genome
blastn -task blastn -evalue 1e-3 -num_threads 4 -num_descriptions 10 -num_alignments 10 -outfmt 0 \
    -dust yes -soft_masking true \
    -db S288c.fa -query ../PARS10/pubs/PARS10/data/sce_genes.fasta -out sce_genes.blast

# parse blastn output
perl ~/Scripts/pars/blastn_transcript.pl -f sce_genes.blast -m 0

```

# Gene filter

## create protein coding gene list

```bash
mkdir -p ~/data/mrna-structure/gene-filter
cd ~/data/mrna-structure/gene-filter

# sgd/saccharomyces_cerevisiae.gff → protein coding gene list
cat ../sgd/saccharomyces_cerevisiae.gff |
    perl -nla -e '
        next if /^#/;
        next unless $F[2] eq q{mRNA};
        my $annotation = $F[8];
        $annotation =~ /ID=(.*)_mRNA;Name=/;
        my $ID = $1;
        my $chr = $F[0];
        $chr =~ s/^chr//i;
        next if $chr eq q{mt}; # Skip genes on mitochondria
        print join qq{,}, $ID, qq{$chr($F[6]):$F[3]-$F[4]};
    ' \
    > protein_coding_list.csv

mkdir -p mRNAs
cat protein_coding_list.csv |
    parallel --colsep ',' --no-run-if-empty --linebuffer -k -j 12 '
        echo {1}
        echo {2} | jrunlist cover stdin -o mRNAs/{1}.yml
    ' 
jrunlist merge mRNAs/*.yml -o mRNAs.merge.yml
rm -fr mRNAs

# overlapped regions
cut -d, -f 2 protein_coding_list.csv |
    jrunlist cover -c 2 stdin -o overlapped.yml

jrunlist statop \
    ../blast/S288c.sizes \
    mRNAs.merge.yml overlapped.yml \
    --op intersect --all -o stdout |
    grep -v "^key" |
    perl -nla -F, -e '
        $F[4] == 0 and print $F[0];
    ' \
    > non-overlapped.lst

# PARS genes
cat non-overlapped.lst |
    grep -Fx -f <(cut -f 1 ../blast/sce_genes.blast.tsv) \
    > PARS-non-overlapped.lst

cat ../blast/sce_genes.blast.tsv |
    perl -nla -e '
        next if /^#/;
        my $ID = $F[0];
        my $chr = $F[2];
        next if $chr eq q{mt}; # Skip genes on mitochondria
        print join qq{,}, $ID, qq{$chr($F[5]):$F[3]-$F[4]};
    ' \
    > PARS_gene_list.csv

mkdir -p PARS
cat PARS_gene_list.csv |
    parallel --colsep ',' --no-run-if-empty --linebuffer -k -j 12 '
        echo {1}
        echo {2} | jrunlist cover stdin -o PARS/{1}.yml
    ' 
jrunlist merge PARS/*.yml -o PARS.merge.yml
rm -fr PARS

jrunlist some mRNAs.merge.yml PARS-non-overlapped.lst -o mRNAs.non-overlapped.yml
#jrunlist split mRNAs.non-overlapped.yml -o mRNAs

jrunlist some PARS.merge.yml PARS-non-overlapped.lst -o PARS.non-overlapped.yml
jrunlist split PARS.non-overlapped.yml -o PARS

```

## Intact mRNAs

```bash
cd ~/data/mrna-structure/gene-filter

gzip -dcf ../alignment/n7/Scer_n7_Spar_refined/*.gz |
    pigz > Scer_n7_Spar.fas.gz
gzip -dcf ../alignment/n7p/Scer_n7p_Spar_refined/*.gz |
    pigz > Scer_n7p_Spar.fas.gz

gzip -dcf ../alignment/n128/Scer_n128_Spar_refined/*.gz |
    pigz > Scer_n128_Spar.fas.gz
gzip -dcf ../alignment/n128/Scer_n128_Seub_refined/*.gz |
    pigz > Scer_n128_Seub.fas.gz

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    echo "==> ${NAME}"
    fasops covers -n S288c ${NAME}.fas.gz -o ${NAME}.yml
done

# intact mRNAs
for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    echo "==> ${NAME}"
    jrunlist statop \
        ../blast/S288c.sizes \
        mRNAs.non-overlapped.yml ${NAME}.yml \
        --op intersect --all -o stdout |
        grep -v "^key" |
        perl -nla -F, -e '
            $F[2] == $F[4] and print $F[0];
        ' \
        > ${NAME}.intact.lst
done

wc -l *.lst ../blast/*.tsv* |
    grep -v "total$" |
    datamash reverse -W |
    (echo -e "File\tCount" && cat) |
    mlr --itsv --omd cat

```

| File                              | Count |
|:----------------------------------|------:|
| non-overlapped.lst                |  5344 |
| PARS-non-overlapped.lst           |  2494 |
| Scer_n128_Seub.intact.lst         |  1500 |
| Scer_n128_Spar.intact.lst         |  1997 |
| Scer_n7p_Spar.intact.lst          |  2270 |
| Scer_n7_Spar.intact.lst           |  2234 |
| ../blast/sce_genes.blast.tsv      |  2980 |
| ../blast/sce_genes.blast.tsv.skip |   216 |

## Cut mRNA alignments and extract SNP list

```bash
cd ~/data/mrna-structure/gene-filter

# PARS slices
for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    echo "==> ${NAME}"
    mkdir -p PARS_${NAME}
    
    cat ${NAME}.intact.lst |
        parallel --no-run-if-empty --linebuffer -k -j 12 "
           fasops slice ${NAME}.fas.gz PARS/{}.yml -n S288c -o PARS_${NAME}/{}.fas
        "
done

# SNP list
for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    echo "==> ${NAME}"
    mkdir -p SNP_${NAME}
    
    cat ${NAME}.intact.lst |
        parallel --no-run-if-empty --linebuffer -k -j 12 "
           fasops vars --outgroup --nocomplex PARS_${NAME}/{}.fas -o SNP_${NAME}/{}.SNPs.tsv
        "

    cat ${NAME}.intact.lst | while read i
    do
        file=SNP_${NAME}/${i}.SNPs.tsv
        export prefix=$(echo ${file} | xargs basename | perl -p -e 's/\.SNPs\.tsv//')
        cat ${file} | perl -nl -e 'print "$_\t$ENV{prefix}"' > SNP_${NAME}/${i}.tsv
        unset prefix
        unset file
    done
    rm -fr SNP_${NAME}/*.SNPs.tsv
    
    cat SNP_${NAME}/*.tsv |
        perl -nla -F"\t" -e 'print qq{$F[4]\t$F[5]\t$F[6]\t$F[8]\t$F[9]\t$F[7]\t$F[13]};' > ${NAME}.total.SNPs.info.tsv #loccation,REF,ALT,mutant_to,freq,occured,gene

done

wc -l *.total.SNPs.info.tsv |
    grep -v "total$" |
    datamash reverse -W |
    (echo -e "File\tCount" && cat) |
    mlr --itsv --omd cat

```

| File | Count |
| --- | --- |
| Scer_n128_Seub.total.SNPs.info.tsv | 30834 |
| Scer_n128_Spar.total.SNPs.info.tsv | 50268 |
| Scer_n7p_Spar.total.SNPs.info.tsv | 38481 |
| Scer_n7_Spar.total.SNPs.info.tsv | 30038 |

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
cat 1011Matrix.tsv | 
    perl -nla -F"\t" -e '
        next if /^#/;
        my $loca = $F[0];
        $loca =~ /^chromosome(\d+)/;
        $chr = &trans($1);
        my $R = length $F[3];
        my $A = length $F[4];
        my @info = split /;/, $F[7];
        my @AF = split /=/, $info[1];
        my $Freq_vcf = $AF[1];
        my @AC = split /=/, $info[0];
        my $REF_vcf = $AC[1];
        my @AN = split /=/, $info[2];
        my $ALT_vcf = $AN[1]-$AC[1];
        if ($R == 1 && $A == 1){
            print qq{$chr:$F[1]\t$F[3]\t$F[4]\t$Freq_vcf\t$REF_vcf\t$ALT_vcf};
        } 
        BEGIN{
            print qq{Location\tREF\tALT\tFreq_vcf\tREF_vcf\tALT_vcf};
            sub trans {
		            my %roman = (16=>"XVI",15=>"XV",14=>"XIV",13=>"XIII",12=>"XII",11=>"XI",10=>"X",9=>"IX",8=>"VIII",7=>"VII",6=>"VI",5=>"V",4=>"IV",3=>"III",2=>"II",1=>"I");
		            $roman{"$_[0]"};
            }
        }
    ' \
> 1011Matrix.ext.tsv
cat 1011Matrix.ext.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > 1011Matrix.ext.txt

# wili in 1011
cd ~/data/mrna-structure/vcf/1011Matrix.gvcf

perl -i -pe 's/chromosome4\t193242.*\n//g;s/chromosome4\t193246.*\n//g;s/chromosome4\t88:2:49\..*\n//g;s/chromosome4\t88:268.*\n//g;' 1011Matrix.gvcf
bcftools view 1011Matrix.gvcf -s CCL,BBQ,BBS,BFP,BTG,CLC,CLB,CLD,BAM,BAQ,BAG,BAH,BAL,AMH,CEG,CEI,CCQ,CCR,CCS,BAK,BAI,ACQ,CCN,CDL,SACE_YCR,BMA,AKM,BMB,BMC,SACE_MAL,SACE_YCY,BAN,BAP,CMP,CCH,ACC,CCC,CCD,CCE,CCF,CCG,CCI,CMQ,CDF,CDG,CDH,CDI,AVI,ACD,ANF,ANH,ANC,ANE,ANG,AND,ANK,ANI,AKN,SACE_YBS,SACE_YCU | bcftools +fill-tags -o 1011Matrix.wild.gvcf
perl ~/Scripts/pars/program/vcf.cut.pl --file 1011Matrix.wild.gvcf --output 1011Matrix.wild.tsv

cat 1011Matrix.wild.tsv | 
    perl -nla -F"\t" -e '
        next if /^#/;
        my $loca = $F[0];
        $loca =~ /^chromosome(\d+)/;
        $chr = &trans($1);
        my $R = length $F[3];
        my $A = length $F[4];
        my @info = split /;/, $F[7];
        my @AF = split /=/, $info[1];
        my $Freq_vcf = $AF[1];
        my @AC = split /=/, $info[0];
        my $REF_vcf = $AC[1];
        my @AN = split /=/, $info[2];
        my $ALT_vcf = $AN[1]-$AC[1];
        if ($R == 1 && $A == 1){
            print qq{$chr:$F[1]\t$F[3]\t$F[4]\t$Freq_vcf\t$REF_vcf\t$ALT_vcf};
        } 
        BEGIN{
            print qq{Location\tREF\tALT\tFreq_vcf\tREF_vcf\tALT_vcf};
            sub trans {
		            my %roman = (16=>"XVI",15=>"XV",14=>"XIV",13=>"XIII",12=>"XII",11=>"XI",10=>"X",9=>"IX",8=>"VIII",7=>"VII",6=>"VI",5=>"V",4=>"IV",3=>"III",2=>"II",1=>"I");
		            $roman{"$_[0]"};
            }
        }
    ' \
> 1011Matrix.ext.wild.tsv
cat 1011Matrix.ext.wild.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > 1011Matrix.ext.wild.txt

```
| File                          | Count |
|:------------------------------|------:|
| 1011Matrix.ext.tsv | 1544488 |
| 1011Matrix.ext.wild.tsv | 1544485 |

# VEP
```bash
mkdir -p ~/data/mrna-structure/vep
cd ~/data/mrna-structure/vep

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    tsv-join --z \
        ~/data/mrna-structure/vcf/1011Matrix.gvcf/1011Matrix.ext.txt \
        -f ../gene-filter/${NAME}.total.SNPs.info.tsv \
        --key-fields 1 \
        --append-fields 2,3,4,5,6,7 \
    > ${NAME}.total.SNPs.info.update.tsv
    cat ${NAME}.total.SNPs.info.update.tsv | datamash check
    cat ${NAME}.total.SNPs.info.update.tsv | 
        perl -nla -F"\t" -e '
            my $loca = $F[0];
            $loca =~ /^(.*):(.*)/;
            my $Chr = $1;
            my $position = $2;
            print qq{$Chr\t$position\t$position\t$F[1]\t$F[2]}; 
        ' \
    > ${NAME}.total.SNPs.update.tsv
done

wc -l *.total.SNPs.update.tsv |
    grep -v "total$" |
    datamash reverse -W |
    (echo -e "File\tCount" && cat) |
    mlr --itsv --omd cat

```
| File | Count |
| --- | --- |
| Scer_n128_Seub.total.SNPs.update.tsv | 27323 |
| Scer_n128_Spar.total.SNPs.update.tsv | 44251 |
| Scer_n7_Spar.total.SNPs.update.tsv | 26814 |
| Scer_n7p_Spar.total.SNPs.update.tsv | 35078 |

upload ${NAME}.total.SNPs.update.tsv in https://asia.ensembl.org/Tools/VEP

Species: Saccharomyces cerevisiae(Saccharomyces cerevisiae)
Additional configurations: 
    Additional_annotations: 
        Upstream/Downstream distance (bp): 1
Download VEP format profiles to `vep`, and rename ${NAME}.total.SNPs.vep.txt


```bash
cd ~/data/mrna-structure/vep

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    cat ${NAME}.total.SNPs.vep.txt | 
        perl -nla -F"\t" -e '
            next if /^#/;
            my $loca = $F[1];
            $loca =~ /^(.*)-[0-9]+/;
            my $ID = $1;
            print qq{$ID\t$F[2]\t$F[3]\t$F[6]\t$F[8]\t$F[10]\t$F[11]\t$F[12]}; #location,allele,gene,consequence,CDS_position,amino_acids,codons,existing_variation
        ' \
    > ${NAME}.total.SNPs.vep.tsv
done

wc -l *.total.SNPs.vep.tsv |
    grep -v "total$" |
    datamash reverse -W |
    (echo -e "File\tCount" && cat) |
    mlr --itsv --omd cat

```
| File | Count |
| --- | --- |
| Scer_n128_Seub.total.SNPs.vep.tsv | 29125 |
| Scer_n128_Spar.total.SNPs.vep.tsv | 47578 |
| Scer_n7_Spar.total.SNPs.vep.tsv | 29109 |
| Scer_n7p_Spar.total.SNPs.vep.tsv | 37964 |

# Process PARS data

```bash
mkdir -p ~/data/mrna-structure/process
cd ~/data/mrna-structure/process

perl ~/Scripts/pars/blastn_transcript.pl -f ../blast/sce_genes.blast -m 0

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    echo "==> ${NAME}"
    
    cat ../gene-filter/${NAME}.total.SNPs.info.csv | 
        perl -nla -F"," -e 'print  qq{$F[0]};' | sort -u  > ${NAME}.snp.gene.pos.txt

    perl ~/Scripts/pars/read_fold.pl \
        --pars ../PARS10/pubs/PARS10/data \
        --gene sce_genes.blast.tsv \
        --pos  ${NAME}.snp.gene.pos.txt \
    > ${NAME}_fail_pos.txt # review fail_pos.txt to find SNPs located in overlapped genes
    
    perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml
done

cd ~/data/mrna-structure/process

#----------------------------------------------------------#
# gene
#----------------------------------------------------------#

# produce transcript set
# YLR167W	568	chrXII	498888	499455	+
cat sce_genes.blast.tsv \
    | perl -nla -e 'print qq{$F[2]:$F[3]-$F[4]}' \
    | sort \
    > sce_genes.pos.txt
jrunlist cover sce_genes.pos.txt -o sce_genes.yml

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

#----------------------------------------------------------#
# utr
#----------------------------------------------------------#
jrunlist compare --op diff sce_genes.yml sce_orf_genomic.yml -o sce_utr.yml
runlist convert sce_utr.yml -o sce_utr.pos.txt

#----------------------------------------------------------#
# mRNA
#----------------------------------------------------------#
jrunlist compare --op diff sce_genes.yml sce_intron.yml -o sce_mRNA.yml
runlist convert sce_mRNA.yml -o sce_mRNA.pos.txt

#----------------------------------------------------------#
# cds
#----------------------------------------------------------#
jrunlist compare --op diff sce_mRNA.yml sce_utr.yml -o sce_cds.yml
runlist convert sce_cds.yml -o sce_cds.pos.txt

# Stats
printf "| %s | %s | %s | %s |\n" \
    "Name" "chrLength" "size" "coverage" \
    > coverage.stat.md
printf "|:--|--:|--:|--:|\n" >> coverage.stat.md

for f in genes intron orf_genomic utr mRNA cds; do
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
| genes | 12071326 | 4236728 | 0.3510 |
| intron | 12071326 | 65144 | 0.0054 |
| orf_genomic | 12071326 | 8895737 | 0.7369 |
| utr | 12071326 | 516701 | 0.0428 |
| mRNA | 12071326 | 4234684 | 0.3508 |
| cds | 12071326 | 3717983 | 0.3080 |


# SNP

## count per gene GC content

```bash
for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    mkdir -p ~/data/mrna-structure/result/${NAME}
    cd ~/data/mrna-structure/result/${NAME}
    perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --varfold ~/data/mrna-structure/process/${NAME}.gene_variation.fold_class.tsv --output ${NAME}.gene_variation.fold_class.csv
done
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

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    cd ~/data/mrna-structure/result/${NAME}
    tsv-join --z \
        ~/data/mrna-structure/vep/${NAME}.total.SNPs.info.update.tsv \
        -f ~/data/mrna-structure/vep/${NAME}.total.SNPs.vep.tsv \
        --key-fields 1 \
        --append-fields 2-8 \
    >${NAME}.tsv
    cat ${NAME}.tsv | 
        perl -nla -F"\t" -e '
            if ($F[8] eq "-" || $F[6] eq $F[8]){
                splice @F, 8, 1 ;
                my $F = join ("\t",@F);
                print qq{$F};
            }
            BEGIN{
                print qq{location\tREF\tALT\tmutant_to\tfreq\toccured\tgene\tallele\tconsequence\tCDS_position\tamino_acids\tcodons\texisting_variation};
            }
        ' \
    > ${NAME}_SNPs_total_info_vep_non_overlapped.tsv
    rm ${NAME}.tsv
    Rscript ~/Scripts/pars/program/stat_SNPs.R -n ${NAME}  
done
unset NAME

```

## count A/T <-> G/C

```bash

for NAME in Scer_n7_Spar Scer_n7p_Spar; do
    cd ~/data/mrna-structure/result/${NAME}
    mkdir -p ~/data/mrna-structure/result/${NAME}/freq_each
    
    Rscript ~/Scripts/pars/program/count_AT_GC.R -n ${NAME} 
    for AREA  in mRNA cds utr syn nsy; do
        perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_${AREA}_stat.csv --output freq_each/PARS_${AREA}_stat_chi_square.csv
    done
    unset AREA
done
unset NAME

for NAME in Scer_n128_Spar Scer_n128_Seub; do
    cd ~/data/mrna-structure/result/${NAME}
    mkdir -p ~/data/mrna-structure/result/${NAME}/freq_each
    mkdir -p ~/data/mrna-structure/result/${NAME}/freq_10
    
    Rscript ~/Scripts/pars/program/count_AT_GC.R -n ${NAME} 
    for AREA  in mRNA cds utr syn nsy; do
        perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_${AREA}_stat.csv --output freq_each/PARS_${AREA}_stat_chi_square.csv
        perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_${AREA}_stat_freq_10.csv --output freq_10/PARS_${AREA}_stat_freq_10_chi_square.csv
    done
    unset AREA
done
unset NAME

```

## count stem length selection

```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME} 
mkdir -p freq_10/stem_length

perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --origin data_SNPs_PARS_mRNA.csv --output data_SNPs_PARS_mRNA_pos.csv

Rscript ~/Scripts/pars/program/count_AT_GC_gene_trait.R -n ${NAME}

cat data_SNPs_PARS_mRNA.csv | perl -nl -a -F"," -e 'print qq{$F[8]};' | sort -u | perl -nl -a -F"," -e 'next if /"gene"/; print qq{$F[0]}; BEGIN{print qq{gene};}' > mRNA.gene.list.csv

for STRUCTURE in stem loop; do
    perl ~/Scripts/pars/program/count_structure_length_gene.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --name ~/data/mrna-structure/result/${NAME}/mRNA.gene.list.csv --structure ${STRUCTURE} --output ${STRUCTURE}_length_mRNA.csv
done
unset STRUCTURE
unset NAME

```

## count_codon_gene

```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME} 
perl ~/Scripts/pars/program/count_codon_gene.pl --origin data_SNPs_PARS_syn.csv --output data_SNPs_PARS_syn_codon.csv
Rscript ~/Scripts/pars/program/count_AT_GC_codon.R -n ${NAME}
for AREA  in tRNA 4D; do
        perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_${AREA}_stat.csv --output freq_each/PARS_${AREA}_stat_chi_square.csv
        perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_${AREA}_stat_freq_10.csv --output freq_10/PARS_${AREA}_stat_freq_10_chi_square.csv
done
unset AREA
unset NAME

```

## count per gene cds_utr

```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME} 
for AREA in cds utr; do
    perl ~/Scripts/pars/program/count_cut_range.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --cut ~/data/mrna-structure/process/sce_${AREA}.yml --output stem_loop_${AREA}_length.csv
    perl ~/Scripts/pars/program/count_per_gene_ACGT_percent.pl --file data_SNPs_PARS_${AREA}.csv --output data_SNPs_PARS_${AREA}_per_gene_ATGC.csv
done

Rscript ~/Scripts/pars/program/count_cds_utr.R -n ${NAME}
unset NAME
```


# GO/KEGG
```bash
cd ~/data/mrna-structure/vcf

export NAME=Scer_n128_Spar
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_mRNA.csv | sed 's/\"//g' |
    perl -nla -F"," -e '
        next if /^location/;
        my $Freq = $F[12]/128;
        my %count;
        my @occured = split //, $F[13];
        my @uniq = grep { ++$count{$_} < 2; } @occured;
        my $REF_pars = $count{ $occured[0] };
        my $REF = $occured[0];
        my $ALT_pars = 128 - $count{ $occured[0] };
        my $ALT;
        if ( $uniq[0] eq $REF ){
            $ALT = $uniq[1];
        }else{
            $ALT = $uniq[0];
        }
        print qq{$F[0]\t$F[8]\t$F[6]\t$F[11]\t$REF\t$ALT\t$Freq\t$REF_pars\t$ALT_pars}; 
        BEGIN{
            print qq{location\tgene\tstructure\tmutant_to\tREF\tALT\tfreq_pars\tREF_pars\tALT_pars};
        }
    ' \
> ${NAME}_data_SNPs_PARS_mRNA.pars.tsv

tsv-join --z \
    ${NAME}_data_SNPs_PARS_mRNA.pars.tsv \
    -f 1011Matrix.gvcf/1011Matrix.ext.wild.tsv  \
    --key-fields 1 \
    --append-fields 2-6\
> ${NAME}_data_SNPs_PARS_mRNA.merge.wild.tsv

cat ${NAME}_data_SNPs_PARS_mRNA.merge.wild.tsv | sed 's/\"//g' |
    perl -nla -F"\t" -e '
        next if /^location/;
        if ($F[4]eq$F[9] && $F[5]eq$F[10]){
            my $minus = $F[6] - $F[11];
            my $obs = [ [ $F[7], $F[8] ], [ $F[12], $F[13] ] ];
            my $chi = new Statistics::ChisqIndep;
            $chi->load_data($obs);
            #$chi->print_summary();
            $Chi = ${$chi}{'chisq_statistic'};
            $P = ${$chi}{'p_value'};
            my $chi = 
            print qq{$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[11]\t$minus\t$Chi\t$P}; 
        }
        BEGIN{
            print qq{location\tgene\tstructure\tmutant_to\tREF\tALT\tfreq_pars\tfreq_vcf\tfreq_minus\tchi\tp};
            use Statistics::ChisqIndep;
        }
    ' \
> ${NAME}_data_SNPs_PARS_mRNA.merge.wild.Chi.tsv

## extract SNP list, 1011=1011_wild
cat ${NAME}_data_SNPs_PARS_mRNA.merge.wild.Chi.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]};' > ${NAME}.mRNA.wild.snp.update.txt

## extract freq_minus>0,p<0.05 in 1011.wild
cat ${NAME}_data_SNPs_PARS_mRNA.merge.wild.Chi.tsv | 
    perl -nla -F"\t" -e '
        next if /^location/;
        if ($F[8]>0 && $F[10]<0.05){
            print qq{$F[0]};
        }
        BEGIN{
            print qq{location};
        }
    ' \
> ${NAME}.mRNA.wild.snp.update.filter.txt

unset NAME

```

## update wild

```bash
export NAME=Scer_n128_Spar
mkdir -p ~/data/mrna-structure/result/${NAME}.update
mkdir -p ~/data/mrna-structure/result/${NAME}.update/freq_each
mkdir -p ~/data/mrna-structure/result/${NAME}.update/freq_10

cd ~/data/mrna-structure/result/${NAME}.update

tsv-join --z \
    ~/data/mrna-structure/vcf/${NAME}.mRNA.wild.snp.update.filter.txt \
    -f <(cat ../${NAME}/data_SNPs_PARS_mRNA.csv | sed 's/\"//g' | csv2tsv) \
    --key-fields 1 \
    --append-fields 2-63 \
> data_SNPs_PARS_mRNA.tsv

Rscript ~/Scripts/pars/program/count_AT_GC.update.R -n ${NAME}

perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_mRNA_stat.csv --output freq_each/PARS_mRNA_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_mRNA_stat_freq_10.csv --output freq_10/PARS_mRNA_stat_freq_10_chi_square.csv

```

upload Scer_n128_Spar.update/mRNA.gene.list.update.csv in https://david.ncifcrf.gov/ , get GO/KEGG information

```bash
mkdir -p freq_10/GO
mkdir -p freq_10/KEGG
Rscript ~/Scripts/pars/program/count_AT_GC_GO_KEGG.R -n ${NAME}.update

#obtain Scer_n128_Spar_go_kegg.csv

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

#get genelist by filtering strong selection from GO/KEGG annotation and deleting repeating item

echo -e "gene\nYDL174C\nYKR066C\nYNR001C\nYDR487C\nYOR065W\nYCL057W\nYIR037W\nYMR203W\nYPL132W\nYLR304C\nYJL054W\nYNL055C\nYJR104C\nYKR071C\nYKL150W\nYBR056W\nYDR226W\nYDR375C\nYKL053C-A\nYKL087C\nYGL187C\nYLR259C\nYFR033C\nYJR121W\nYMR145C\nYAL039C\nYKL067W\nYDL120W\nYDR353W\nYLL009C\nYDR511W\nYPR140W\nYJL143W\nYHR116W\nYNR022C\nYDR296W\nYHR147C\nYBR268W\nYDR116C\nYKR006C\nYLR439W\nYPL183W-A\nYMR286W\nYDL202W\nYKL138C\nYJL063C\nYMR024W\nYNL005C\nYBR122C\nYJL096W\nYOR150W\nYDR237W\nYBL038W\nYLR312W-A\nYKL167C\nYKR085C\nYLR008C\nYKL016C\nYDR204W\nYMR241W\nYKL148C\nYGR096W\nYJL166W\nYOR065W\nYOR297C\nYMR089C\nYGR257C\nYPL134C\nYCL044C\nYPR024W\nYJR045C\nYKR087C\nYGL187C\nYBR039W\nYEL052W\nYGR222W\nYER017C\nYKR065C\nYLR259C\nYFR033C\nYMR035W\nYBR185C\nYMR256C\nYOR176W\nYGR033C\nYPR140W\nYML110C\nYKL141W\nYLR203C\nYHR024C\nYBR085W\nYDL174C\nYMR301C\nYER078C\nYML120C\nYPL270W\nYLR253W\nYLL041C\nYLR164W\nYNL003C\nYER141W\nYPR191W\nYGR235C\nYPL132W\nYJL054W\nYBL030C\nYGR062C\nYOL008W\nYPL063W\nYDR298C\nYLR188W\nYDR236C\nYDR375C\nYKR052C\nYKL087C\nYEL024W\nYGR101W\nYHR037W\nYDR377W\nYPL078C\nYPL271W\nYJR121W\nYOR232W\nYOR356W\nYBR291C\nYNL100W\nYAL039C\nYBR003W\nYDL120W\nYDL004W\nYPL189C-A\nYOR125C\nYCL057C-A\nYKL120W\nYIL134W\nYIL022W\nYOR222W\nYJL143W\nYOL027C\nYGR082W\nYNL026W\nYLR099W-A\nYLR090W\nYML086C\nYNL055C\nYKL150W\nYNL131W\nYMR110C\nYNL121C\nYER019W\nYPR140W\nYNL070W\nYHR117W\nYLR008C\nYGR082W\nYER017C\nYKR065C\nYOR232W\nYMR203W\nYLR090W\nYPL063W\nYPR024W\nYNL131W\nYJR045C\nYNL121C\nYGR033C\nYIL022W\nYNL070W\nYJL143W\nYHR024C\nYLR008C\nYDL174C\nYNR001C\nYOR065W\nYMR203W\nYOR297C\nYLR304C\nYJL054W\nYNL055C\nYPL063W\nYDR375C\nYKL053C-A\nYGL187C\nYNL121C\nYHR117W\nYNL070W\nYGR082W\nYKR065C\nYNL026W\nYLR259C\nYOR027W\nYJR121W\nYOR232W\nYHR083W\nYGR028W\nYNL064C\nYDL120W\nYNL131W\nYLL009C\nYPR140W\nYGR033C\nYIL022W\nYJL143W\nYMR301C\nYIL003W\nYPL270W\nYMR312W\nYCR011C\nYHR169W\nYGL008C\nYMR089C\nYDL007W\nYPR173C\nYLR397C\nYOR259C\nYDR091C\nYLR188W\nYDL166C\nYNL329C\nYHR187W\nYJR045C\nYGR262C\nYDR375C\nYBR039W\nYGR210C\nYBL022C\nYDL126C\nYFL028C\nYER017C\nYER036C\nYDR377W\nYLR259C\nYDR061W\nYPL271W\nYJR121W\nYJR072C\nYNL290W\nYCL047C\nYOR291W\nYEL031W\nYOR117W\nYGR028W\nYDL100C\nYLR249W\nYDL004W\nYPL226W\nYOR153W\nYKL073W\nYGL048C\nYBR080C\nYDR011W\nYOR278W\nYAL039C\nYDR047W\nYDL120W\nYDR232W\nYKL087C\nYER141W\nYGL040C\nYDR044W\nYGL245W\nYPL172C\nYOR176W\nYDL174C\nYDL078C\nYPL061W\nYDR272W\nYML004C\nYBL015W\nYBR218C\nYOR374W\nYOR347C\nYPL262W\nYKL085W\nYPL028W\nYER178W\nYMR110C\nYLR153C\nYNL071W\nYBR221C\nYNR016C" | sort | uniq | perl -e 'print reverse <>' > genelist.csv

#filter SNPs
RScript ~/Scripts/pars/program/subpop.R -n ${NAME} -i genelist.csv -o filter_snp.csv

#delete "double quotation marks" and "blank" in filter_snp.csv

#calculate subpopulation SNPs proporation
perl ~/Scripts/pars/program/subpop.pl filter_snp.csv strainlist.csv > subpop.csv

RScript ~/Scripts/pars/program/subpop_merge.R -n ${NAME}

unset NAME

```

## stat subpopulation SNPs frequency

```bash
export NAME=Scer_n128_Spar
mkdir -p ~/data/mrna-structure/result/${NAME}.update/subpop
cd ~/data/mrna-structure/result/${NAME}.update/subpop

#generate strainlist, order same as egaz template
cat ~/Scripts/pars/group_phylo.tsv | grep -v "^#" | cut -f 2 | tr "," "\n" | perl -nl -a -F"\t" -e 'print qq{$F[0]};BEGIN{print qq{S288c};}' > strainlist.tsv

tsv-join --z \
    <(cat ~/data/mrna-structure/vcf/${NAME}_data_SNPs_PARS_mRNA.merge.wild.Chi.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\t$F[11]};') \
    -f ../data_SNPs_PARS_mRNA.tsv \
    --key-fields 1 \
    --append-fields 2-63 \
> data_SNPs_PARS_mRNA_all.tsv

perl ~/Scripts/pars/program/subpop.pl data_SNPs_PARS_mRNA_all.tsv strainlist.tsv > subpop.csv
rm data_SNPs_PARS_mRNA_all.tsv

tsv-join --z \
    <(cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_syn.csv | sed 's/\"//g' \
        | perl -nl -a -F"," -e '
              next if /^location/;
              if (($F[6]eq"stem")&&($F[11]eq"A->G"||$F[11]eq"A->C"||$F[11]eq"T->G"||$F[11]eq"T->C")) {
                  my $F = join ("\t",@F);
                  print qq{$F};
              }
              BEGIN{
                  print qq{location\tfold_length\tgene_base\tgene_pos\tpars\tpair_base\tstructure\tstrand\tgene\tREF\tALT\tmutant_to\tfreq\toccured\tallele\tconsequence\tCDS_position\tamino_acids\tcodons\texisting_variation\tlength\tmF\tfold_dot_length\tfold_dot_vars\tfold_left_length\tfold_left_vars\tfold_right_length\tfold_right_vars\tstem_A_num\tstem_A_per\tstem_C_num\tstem_C_per\tstem_G_num\tstem_G_per\tstem_U_num\tstem_U_per\tloop_A_num\tloop_A_per\tloop_C_num\tloop_C_per\tloop_G_num\tloop_G_per\tloop_U_num\tloop_U_per\tA_num\tA_per\tC_num\tC_per\tG_num\tG_per\tU_num\tU_per\tstem_AU_num\tstem_CG_num\tloop_AU_num\tloop_CG_num\tAU_num\tCG_num\tstem_CG_content\tloop_CG_content\tCG_content\tX2\tP};
              }
          ') \
    -f <(cat subpop.csv | csv2tsv) \
    --key-fields 1 \
    --append-fields 3-11 \
> subpop.syn.tsv












