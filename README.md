# Processing Yeast PARS Data

[TOC level=1-3]: # " "
- [Processing Yeast PARS Data](#processing-yeast-pars-data)
- [Install needed softwares](#install-needed-softwares)
    - [Homebrew packages](#homebrew-packages)
    - [Perl modules](#perl-modules)
    - [R packages](#r-packages)
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
- [VCF of 1011 project](#vcf-of-1011-project)
- [VEP](#vep)
- [Process PARS data](#process-pars-data)
- [SNP](#snp)
    - [count per gene GC content](#count-per-gene-gc-content)
    - [count SNPs and gene](#count-snps-and-gene)
    - [count A/T <-> G/C](#count-at---gc)
    - [count stem length selection](#count-stem-length-selection)
    - [count_codon_gene](#count_codon_gene)
    - [count per gene cds_utr](#count-per-gene-cds_utr)
- [GO/KEGG](#gokegg)
    - [update wild](#update-wild)
- [subpop](#subpop)
    - [stat subpopulation SNPs frequency](#stat-subpopulation-snps-frequency)


# Install needed softwares

## Homebrew packages

```bash
brew install bcftools blast pigz samtools   

brew tap brewsci/bio
brew install mafft raxml

brew tap wang-q/tap
brew install faops lastz multiz sparsemem

curl -fsSL https://raw.githubusercontent.com/wang-q/App-Egaz/master/share/check_dep.sh | bash

```

## Perl modules

```bash
cpanm App::Fasops App::Rangeops App::Egaz

```

## R packages

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

pigz -dcf ../alignment/n7/Scer_n7_Spar_refined/*.gz |
    pigz > Scer_n7_Spar.fas.gz
pigz -dcf ../alignment/n7p/Scer_n7p_Spar_refined/*.gz |
    pigz > Scer_n7p_Spar.fas.gz

pigz -dcf ../alignment/n128/Scer_n128_Spar_refined/*.gz |
    pigz > Scer_n128_Spar.fas.gz
pigz -dcf ../alignment/n128/Scer_n128_Seub_refined/*.gz |
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

| File                               | Count |
|:-----------------------------------|:------|
| Scer_n128_Seub.total.SNPs.info.tsv | 30834 |
| Scer_n128_Spar.total.SNPs.info.tsv | 50268 |
| Scer_n7p_Spar.total.SNPs.info.tsv  | 38481 |
| Scer_n7_Spar.total.SNPs.info.tsv   | 30038 |

# VCF of 1011 project

```bash
mkdir -p ~/data/mrna-structure/vcf
cd ~/data/mrna-structure/vcf
aria2c -c http://1002genomes.u-strasbg.fr/files/1011Matrix.gvcf.gz
pigz -dcf 1011Matrix.gvcf.gz > 1011Matrix.gvcf

# 1011
cat 1011Matrix.gvcf |
    perl -nla -F"\t" -e '
        /^\#\#/ and next;
        splice @F, 8;
        print join qq{\t}, @F;
    ' \
    > 1011Matrix.tsv

cat 1011Matrix.tsv |
    perl -nla -F"\t" -e '
        BEGIN {
            print qq{Location\tREF\tALT\tFreq_vcf\tREF_vcf\tALT_vcf};
            our %roman = (
                16 => "XVI",
                15 => "XV",
                14 => "XIV",
                13 => "XIII",
                12 => "XII",
                11 => "XI",
                10 => "X",
                9  => "IX",
                8  => "VIII",
                7  => "VII",
                6  => "VI",
                5  => "V",
                4  => "IV",
                3  => "III",
                2  => "II",
                1  => "I"
            );
        }
        next if /^#/;
        my $loca = $F[0];
        $loca =~ /^chromosome(\d+)/;
        $chr = $roman{$1};
        my $R        = length $F[3];
        my $A        = length $F[4];
        my @info     = split /;/, $F[7];
        my @AF       = split /=/, $info[1];
        my $Freq_vcf = $AF[1];
        my @AC       = split /=/, $info[0];
        my @AN       = split /=/, $info[2];
        my $ALT_vcf  = $AC[1];
        my $REF_vcf  = $AN[1] - $AC[1];

        if ( $R == 1 && $A == 1 ) {
            print qq{$chr:$F[1]\t$F[3]\t$F[4]\t$Freq_vcf\t$REF_vcf\t$ALT_vcf};
        }
    ' \
    > 1011Matrix.ext.tsv
cut -f 1 1011Matrix.ext.tsv > 1011Matrix.ext.txt

# wild strains in 1011
#perl -pi -e '
#    s/chromosome4\t193242.*\n//g;
#    s/chromosome4\t193246.*\n//g;
#    s/chromosome4\t88:2:49\..*\n//g;
#    s/chromosome4\t88:268.*\n//g;
#    ' 1011Matrix.gvcf
bcftools view 1011Matrix.gvcf -s CCL,BBQ,BBS,BFP,BTG,CLC,CLB,CLD,BAM,BAQ,BAG,BAH,BAL,AMH,CEG,CEI,CCQ,CCR,CCS,BAK,BAI,ACQ,CCN,CDL,SACE_YCR,BMA,AKM,BMB,BMC,SACE_MAL,SACE_YCY,BAN,BAP,CMP,CCH,ACC,CCC,CCD,CCE,CCF,CCG,CCI,CMQ,CDF,CDG,CDH,CDI,AVI,ACD,ANF,ANH,ANC,ANE,ANG,AND,ANK,ANI,AKN,SACE_YBS,SACE_YCU |
    bcftools +fill-tags -o 1011Matrix.wild.gvcf
cat 1011Matrix.wild.gvcf |
    perl -nla -F"\t" -e '
        /^\#\#/ and next;
        splice @F, 8;
        print join qq{\t}, @F;
    ' \
    > 1011Matrix.wild.tsv

cat 1011Matrix.wild.tsv | 
    perl -nla -F"\t" -e '
        BEGIN {
            print qq{Location\tREF\tALT\tFreq_vcf\tREF_vcf\tALT_vcf};
            our %roman = (
                16 => "XVI",
                15 => "XV",
                14 => "XIV",
                13 => "XIII",
                12 => "XII",
                11 => "XI",
                10 => "X",
                9  => "IX",
                8  => "VIII",
                7  => "VII",
                6  => "VI",
                5  => "V",
                4  => "IV",
                3  => "III",
                2  => "II",
                1  => "I"
            );
        }
        next if /^#/;
        my $loca = $F[0];
        $loca =~ /^chromosome(\d+)/;
        $chr = $roman{$1};
        my $R        = length $F[3];
        my $A        = length $F[4];
        my @info     = split /;/, $F[7];
        my @AF       = split /=/, $info[1];
        my $Freq_vcf = $AF[1];
        my @AC       = split /=/, $info[0];
        my @AN       = split /=/, $info[2];
        my $ALT_vcf  = $AC[1];
        my $REF_vcf  = $AN[1] - $AC[1];

        if ( $R == 1 && $A == 1 ) {
            print qq{$chr:$F[1]\t$F[3]\t$F[4]\t$Freq_vcf\t$REF_vcf\t$ALT_vcf};
        }
    ' \
    > 1011Matrix.ext.wild.tsv
cut -f 1 1011Matrix.ext.wild.tsv > 1011Matrix.ext.wild.txt

rm 1011Matrix.gvcf 1011Matrix.wild.gvcf

```

| File                    |   Count |
|:------------------------|--------:|
| 1011Matrix.ext.tsv      | 1544488 |
| 1011Matrix.ext.wild.tsv | 1544485 |

# VEP

```bash
mkdir -p ~/data/mrna-structure/vep
cd ~/data/mrna-structure/vep

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    tsv-join --z \
        ~/data/mrna-structure/vcf/1011Matrix.ext.txt \
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

| File                                 | Count |
|:-------------------------------------|:------|
| Scer_n128_Seub.total.SNPs.update.tsv | 27323 |
| Scer_n128_Spar.total.SNPs.update.tsv | 44251 |
| Scer_n7_Spar.total.SNPs.update.tsv   | 26814 |
| Scer_n7p_Spar.total.SNPs.update.tsv  | 35078 |

Upload ${NAME}.total.SNPs.update.tsv to https://asia.ensembl.org/Tools/VEP

* Species: Saccharomyces cerevisiae (Saccharomyces cerevisiae)
* Additional_annotations:
    * Upstream/Downstream distance (bp): 1

* Download VEP format profiles to `vep/`, and rename it to `${NAME}.total.SNPs.vep.txt`


```bash
cd ~/data/mrna-structure/vep

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    cat ${NAME}.total.SNPs.vep.txt | 
        perl -nla -F"\t" -e '
            next if /^#/;
            my $loca = $F[1];
            $loca =~ /^(.*)-[0-9]+/;
            my $ID = $1;
            #location,allele,gene,consequence,CDS_position,amino_acids,codons,existing_variation
            print qq{$ID\t$F[2]\t$F[3]\t$F[6]\t$F[8]\t$F[10]\t$F[11]\t$F[12]};
        ' \
    > ${NAME}.total.SNPs.vep.tsv
done

wc -l *.total.SNPs.vep.tsv |
    grep -v "total$" |
    datamash reverse -W |
    (echo -e "File\tCount" && cat) |
    mlr --itsv --omd cat

```

| File                              | Count |
|:----------------------------------|:------|
| Scer_n128_Seub.total.SNPs.vep.tsv | 27324 |
| Scer_n128_Spar.total.SNPs.vep.tsv | 44252 |
| Scer_n7_Spar.total.SNPs.vep.tsv   | 26818 |
| Scer_n7p_Spar.total.SNPs.vep.tsv  | 35082 |

# Process PARS data

```bash
mkdir -p ~/data/mrna-structure/process
cd ~/data/mrna-structure/process

perl ~/Scripts/pars/blastn_transcript.pl -f ../blast/sce_genes.blast -m 0

for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    echo "==> ${NAME}"
    
    cat ../vep/${NAME}.total.SNPs.info.update.tsv | 
        perl -nla -F"\t" -e 'print  qq{$F[0]};' |
        sort -u \
        > ${NAME}.snp.gene.pos.txt

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
# YLR167W 568 chrXII 498888 499455 +
cat sce_genes.blast.tsv |
    perl -nla -e 'print qq{$F[2]:$F[3]-$F[4]}' |
    sort \
    > sce_genes.pos.txt
jrunlist cover sce_genes.pos.txt -o sce_genes.yml

#----------------------------------------------------------#
# intron
#----------------------------------------------------------#
cat ../sgd/orf_coding_all.fasta |
    perl -n -MAlignDB::IntSpan -e '
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
cat ../sgd/orf_genomic_all.fasta |
    perl -n -e '
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
            jrunlist stat ../blast/S288c.sizes sce_${f}.yml --all -o stdout |
                grep -v coverage |
                sed "s/,/ /g"
        )
done >> coverage.stat.md

cat coverage.stat.md

```

| Name        | chrLength |    size | coverage |
|:------------|----------:|--------:|---------:|
| genes       |  12071326 | 4236728 |   0.3510 |
| intron      |  12071326 |   65144 |   0.0054 |
| orf_genomic |  12071326 | 8895737 |   0.7369 |
| utr         |  12071326 |  516701 |   0.0428 |
| mRNA        |  12071326 | 4234684 |   0.3508 |
| cds         |  12071326 | 3717983 |   0.3080 |

# SNP

## count per gene GC content

```bash
for NAME in Scer_n7_Spar Scer_n7p_Spar Scer_n128_Spar Scer_n128_Seub; do
    mkdir -p ~/data/mrna-structure/result/${NAME}
    cd ~/data/mrna-structure/result/${NAME}
    perl ~/Scripts/pars/program/count_ACGT_percent.pl \
        --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml \
        --varfold ~/data/mrna-structure/process/${NAME}.gene_variation.fold_class.tsv \
        --output ${NAME}.gene_variation.fold_class.csv
done

```

## count SNPs and gene

```bash
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
done

for NAME in Scer_n128_Spar Scer_n128_Seub; do
    cd ~/data/mrna-structure/result/${NAME}
    mkdir -p ~/data/mrna-structure/result/${NAME}/freq_each
    mkdir -p ~/data/mrna-structure/result/${NAME}/freq_10
    
    Rscript ~/Scripts/pars/program/count_AT_GC.R -n ${NAME} 
    for AREA  in mRNA cds utr syn nsy; do
        perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_${AREA}_stat.csv --output freq_each/PARS_${AREA}_stat_chi_square.csv
        perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_${AREA}_stat_freq_10.csv --output freq_10/PARS_${AREA}_stat_freq_10_chi_square.csv
    done
done

```

## count stem length selection

```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/${NAME} 
mkdir -p freq_10/stem_length

perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml --origin data_SNPs_PARS_mRNA.csv --output data_SNPs_PARS_mRNA_pos.csv

Rscript ~/Scripts/pars/program/count_AT_GC_gene_trait.R -n ${NAME}

for CATEGORY in stem_AT_GC stem_GC_AT loop_AT_GC loop_GC_AT; do
    cat freq_10/stem_length/PARS_mRNA_1_stat_${CATEGORY}_freq_10.csv | csv2tsv > list.tmp
    for ((i=2; i<=15; i++)) 
    do     
        tsv-join \
            list.tmp \
            -f <(cat freq_10/stem_length/PARS_mRNA_${i}_stat_${CATEGORY}_freq_10.csv | csv2tsv) \
            --key-fields 1 \
            --append-fields 2 \
        > list.tmp.bak
        cat list.tmp.bak > list.tmp
    done
    cat list.tmp > stem_length_PARS_mRNA_stat_${CATEGORY}_freq_10.tsv
    rm list.tmp
    rm list.tmp.bak
done

cat data_SNPs_PARS_mRNA.csv |
    perl -nl -a -F"," -e 'print qq{$F[8]};' |
    sort -u |
    perl -nl -a -F"," -e 'next if /"gene"/; print qq{$F[0]}; BEGIN{print qq{gene};}' \
    > mRNA.gene.list.csv

for STRUCTURE in stem loop; do
    perl ~/Scripts/pars/program/count_structure_length_gene.pl \
        --file ~/data/mrna-structure/process/${NAME}.gene_variation.process.yml \
        --name ~/data/mrna-structure/result/${NAME}/mRNA.gene.list.csv \
        --structure ${STRUCTURE} \
        --output ${STRUCTURE}_length_mRNA.csv
done

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
cat ~/data/mrna-structure/result/${NAME}/data_SNPs_PARS_mRNA.csv |
    sed 's/\"//g' |
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
    --append-fields 2-6 |
    perl -nla -F"\t" -e '
        my $mutant_to = $F[3];
        $mutant_to =~ /^(.*)->/;
        my $REF = $1;
        if ($REF ne $F[9]){
            $F[11] = 1 - $F[11];
        }
        my $F = join("\t",@F);
        print qq{$F};
    ' \
    > ${NAME}_data_SNPs_PARS_mRNA.merge.wild.tsv

cat ${NAME}_data_SNPs_PARS_mRNA.merge.wild.tsv |
    sed 's/\"//g' |
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
cat ${NAME}_data_SNPs_PARS_mRNA.merge.wild.Chi.tsv |
    perl -nl -a -F"\t" -e 'print qq{$F[0]};' \
    > ${NAME}.mRNA.wild.snp.update.txt

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

unset NAME

```

upload Scer_n128_Spar.update/mRNA.gene.list.update.csv in https://david.ncifcrf.gov/ , get GO/KEGG
information

```bash
export NAME=Scer_n128_Spar

cd ~/data/mrna-structure/result/${NAME}.update

mkdir -p freq_10/GO
mkdir -p freq_10/KEGG
Rscript ~/Scripts/pars/program/count_AT_GC_GO_KEGG.R -n ${NAME}.update

# process
##BP
for CATEGORY in stem_AT_GC; do
    cat freq_10/GO/BP_1_stat_${CATEGORY}_freq_10.csv | csv2tsv > list.tmp
    for ((i=2; i<=123; i++)) 
    do     
        tsv-join \
            list.tmp \
            -f <(cat freq_10/GO/BP_${i}_stat_${CATEGORY}_freq_10.csv | csv2tsv) \
            --key-fields 1 \
            --append-fields 2 \
        > list.tmp.bak
        cat list.tmp.bak > list.tmp
    done
    cat list.tmp > BP_stat_${CATEGORY}_freq_10.tsv
    rm list.tmp
    rm list.tmp.bak
done

##CC
for CATEGORY in stem_AT_GC; do
    cat freq_10/GO/CC_1_stat_${CATEGORY}_freq_10.csv | csv2tsv > list.tmp
    for ((i=2; i<=72; i++)) 
    do     
        tsv-join \
            list.tmp \
            -f <(cat freq_10/GO/CC_${i}_stat_${CATEGORY}_freq_10.csv | csv2tsv) \
            --key-fields 1 \
            --append-fields 2 \
        > list.tmp.bak
        cat list.tmp.bak > list.tmp
    done
    cat list.tmp > CC_stat_${CATEGORY}_freq_10.tsv
    rm list.tmp
    rm list.tmp.bak
done

##MF
for CATEGORY in stem_AT_GC; do
    cat freq_10/GO/MF_1_stat_${CATEGORY}_freq_10.csv | csv2tsv > list.tmp
    for ((i=2; i<=52; i++)) 
    do     
        tsv-join \
            list.tmp \
            -f <(cat freq_10/GO/MF_${i}_stat_${CATEGORY}_freq_10.csv | csv2tsv) \
            --key-fields 1 \
            --append-fields 2 \
        > list.tmp.bak
        cat list.tmp.bak > list.tmp
    done
    cat list.tmp > MF_stat_${CATEGORY}_freq_10.tsv
    rm list.tmp
    rm list.tmp.bak
done

##KEGG
for CATEGORY in stem_AT_GC; do
    cat freq_10/KEGG/KEGG_1_stat_${CATEGORY}_freq_10.csv | csv2tsv > list.tmp
    for ((i=2; i<=30; i++)) 
    do     
        tsv-join \
            list.tmp \
            -f <(cat freq_10/KEGG/KEGG_${i}_stat_${CATEGORY}_freq_10.csv | csv2tsv) \
            --key-fields 1 \
            --append-fields 2 \
        > list.tmp.bak
        cat list.tmp.bak > list.tmp
    done
    cat list.tmp > KEGG_stat_${CATEGORY}_freq_10.tsv
    rm list.tmp
    rm list.tmp.bak
done

#evaluate γ (Matlab), obtain Scer_n128_Spar_go_kegg.csv

mkdir -p freq_10/go_kegg
mkdir -p freq_10/go_kegg/syn
mkdir -p freq_10/go_kegg/nsy
Rscript ~/Scripts/pars/program/count_AT_GC_GO_KEGG_SN.R -n ${NAME}.update

# process
##go_kegg
for AREA in syn nsy; do
    for CATEGORY in stem_AT_GC; do
        cat freq_10/go_kegg/${AREA}/go_kegg_1_stat_${CATEGORY}_freq_10.csv | csv2tsv > list.tmp
        for ((i=2; i<=41; i++)) 
        do     
            tsv-join \
            list.tmp \
                -f <(cat freq_10/go_kegg/${AREA}/go_kegg_${i}_stat_${CATEGORY}_freq_10.csv | csv2tsv) \
                --key-fields 1 \
                --append-fields 2 \
            > list.tmp.bak
            cat list.tmp.bak > list.tmp
        done
        cat list.tmp > go_kegg_stat_${CATEGORY}_freq_10_${AREA}.tsv
        rm list.tmp
        rm list.tmp.bak
    done
done
unset NAME

```

# subpop

## stat subpopulation SNPs frequency

```bash
export NAME=Scer_n128_Spar
mkdir -p ~/data/mrna-structure/result/${NAME}.update/subpop
cd ~/data/mrna-structure/result/${NAME}.update/subpop

#get genelist by filtering strong selection from GO/KEGG annotation (mt) and deleting repeating item
echo -e "YLR305C\nYHR037W\nYHR063C\nYML110C\nYDR232W\nYER183C\nYMR205C\nYER017C\nYDL126C\nYBR111C\nYBR085W\nYBL030C\nYDR477W\nYMR301C\nYOR153W\nYLR188W\nYPL270W\nYDR011W\nYBL022C\nYDR194C\nYDL181W\nYGR008C\nYDL203C\nYEL052W\nYHL021C\nYHR198C\nYIL087C\nYMR157C\nYOR215C\nYPR004C\nYHR199C\nYLR370C\nYJL178C\nYLR356W\nYNL315C\nYER155C\nYNL305C\nYPR113W\nYER026C\nYJR060W\nYPL215W\nYIL043C\nYLR390W-A\nYGR207C\nYKL137W\nYPL189C-A\nYPL132W\nYML129C\nYER141W\nYHR116W\nYLR380W\nYBR291C\nYBR037C\nYMR264W\nYML086C\nYDL174C\nYDL178W\nYGL120C\nYOR151C\nYLL006W\nYLR390W\nYOR246C\nYAR002C-A\nYGL002W\nYNL125C\nYPL078C\nYDR298C\nYJR121W\nYKL016C\nYDL004W\nYPL271W\nYDR377W\nYPR020W\nYBR039W\nYOL077W-A\nYKR052C\nYGL091C\nYLL027W\nYIL065C\nYLR077W\nYFL046W\nYNL168C\nYER004W\nYOR271C\nYGL225W\nYGL020C\nYBR004C\nYOR355W\nYLR091W\nYDL198C\nYGL008C\nYKL084W\nYOR020C\nYJR045C\nYPR067W\nYER048W-A\nYJL094C\nYCL064C\nYCL005W\nYIL070C\nYBR185C\nYOR298C-A\nYLR253W\nYKL053C-A\nYLR069C\nYOR232W\nYJR144W\nYCL044C\nYDR296W\nYKL195W\nYCL057C-A\nYGR235C\nYNL100W\nYMR002W\nYLR190W\nYML128C\nYOR354C\nYGR028W\nYLR203C\nYGR257C\nYCL033C\nYDR493W\nYMR145C\nYML120C\nYNL200C\nYHR179W\nYPR048W\nYER125W\nYHR086W\nYGR147C\nYPR149W\nYHR042W\nYPL226W\nYKL120W\nYBR129C\nYKR065C\nYLR008C\nYGR178C\nYGR193C\nYGR222W\nYNL003C\nYBR106W\nYGL023C\nYMR129W\nYFL005W\nYOR089C\nYFL038C\nYER031C\nYML001W\nYEL037C\nYKL205W\nYLR084C\nYNR018W\nYLR059C\nYPR165W\nYDR233C\nYNL063W\nYHR083W\nYNL026W\nYKR042W\nYDR511W\nYBR269C\nYGL228W\nYIL016W\nYEL059C-A\nYML052W\nYPL105C\nYML072C\nYBR091C\nYOR297C\nYGR033C\nYJL054W\nYKL056C\nYGR082W\nYNL131W\nYMR203W\nYPR133W-A\nYOR045W\nYNL070W\nYCR083W\nYLR193C\nYOL129W\nYLR090W\nYMR241W\nYMR152W\nYBR054W\nYHR017W\nYDR326C\nYOL109W\nYNL310C\nYMR243C\nYOR221C\nYMR108W\nYNR016C\nYBL015W\nYOL140W\nYLR304C\nYJL200C\nYKL192C\nYKL094W\nYIL124W\nYJL005W\nYDR226W\nYLR089C\nYMR083W\nYGL256W\nYOR374W\nYPL061W\nYIL125W\nYDR148C\nYJR062C\nYDR046C\nYHR189W\nYER078C\nYEL063C\nYKL106W\nYLR027C\nYDR375C\nYOL066C\nYDL168W\nYJL130C\nYHR190W\nYIR038C\nYJR073C\nYJR057W\nYHR208W\nYLR259C\nYNR001C\nYCR005C\nYLL009C\nYGR155W\nYJR048W\nYGL187C\nYMR256C\nYLR395C\nYDL067C\nYGL191W\nYLR038C\nYNL052W\nYKL087C\nYKL150W\nYKR066C\nYFL018C\nYNL071W\nYJR016C\nYPR183W\nYKR001C\nYLL001W\nYKR071C\nYOL021C\nYOR176W\nYDL120W\nYIL134W\nYCR004C\nYKL060C\nYPL262W\nYEL047C\nYBR229C\nYBR196C\nYMR062C\nYDL215C\nYDL171C\nYGL245W\nYHL032C\nYOL059W\nYAL044C\nYMR189W\nYBR263W\nYER020W\nYMR110C\nYGL253W\nYPR033C\nYPL148C\nYAL039C\nYDR234W\nYGL221C\nYDL157C\nYDR061W\nYGR266W\nYJL127C-B\nYJL133C-A\nYJR085C\nYKR018C\nYKR070W\nYMR252C\nYNL122C\nYNL320W\nYNR040W\nYPR024W\nYJL104W\nYMR267W\nYER019W\nYEL031W\nYDL066W\nYBL098W\nYJL060W\nYPL004C\nYGR086C\nYOR317W\nYNL073W\nYPR140W\nYMR089C\nYOR206W\nYKL085W\nYGR062C\nYER154W\nYOR316C\nYFR044C\nYKR087C\nYCL057W\nYLR389C\nYPL134C\nYOR222W\nYDR347W\nYPL118W\nYPL013C\nYDR337W\nYBR251W\nYMR158W\nYBR146W\nYNL137C\nYER050C\nYNR037C\nYIL093C\nYFR049W\nYDL045W-A\nYKL003C\nYNL306W\nYCR046C\nYKL167C\nYDR116C\nYNR022C\nYDL202W\nYKR006C\nYLR312W-A\nYHR147C\nYNL005C\nYKR085C\nYOR150W\nYDR462W\nYMR024W\nYKL138C\nYMR286W\nYBR122C\nYBR268W\nYML009C\nYLR439W\nYBL038W\nYJL096W\nYDR237W\nYJL063C\nYHR162W\nYHR024C\nYBL016W\nYKL113C\nYAL029C\nYKL067W\nYBR104W\nYKL196C\nYDL040C\nYML078W\nYNL135C\nYIR037W\nYPL206C\nYKL212W\nYBR200W\nYGR254W\nYNL055C\nYGR231C\nYGL011C\nYOL038W\nYFR004W\nYLR164W\nYNL121C\nYHR117W\nYIL022W\nYPL063W\nYHR005C-A\nYJL143W\nYNR017W\nYJR135W-A\nYPL154C\nYPL172C\nYNL292W\nYDR384C\nYGR012W\nYDR196C\nYOR356W\nYLR351C\nYER057C\nYOR196C\nYNR036C\nYNL081C\nYPL183W-A\nYOR125C\nYIL114C\nYKR049C\nYER178W\nYBR221C\nYOR347C\nYGL080W\nYOR004W\nYGR101W\nYDR236C\nYOL027C\nYCR008W\nYHR135C\nYJL109C\nYML008C\nYKL141W\nYKL148C\nYLL041C\nYDR178W\nYGR244C\nYJR010W\nYJR104C\nYHR008C\nYDR120C\nYER168C\nYLR105C\nYPL083C\nYHR003C\nYKR079C\nYKL182W\nYGR096W\nYLR043C\nYBL064C\nYDR353W\nYHR106W\nYOR286W\nYER086W\nYDL015C\nYBR003W\nYOR187W\nYLR291C\nYGR162W\nYAL035W\nYDL217C\nYDR074W\nYOR086C\nYNL256W\nYPL231W\nYBR084W\nYDR050C\nYDR268W\nYDL230W\nYOR065W\nYEL024W\nYPR191W\nYFR033C\nYDR529C\nYJL166W\nYLR290C\nYDR204W\nYOL008W\nYML021C\nYGR094W\nYGR285C" | sort -u | perl -nl -a -F"\t" -e 'print qq{$F[0]};BEGIN{print qq{gene};}'> genelist.csv

#generate strainlist, order same as egaz template
cat ~/Scripts/pars/group_phylo.tsv | grep -v "^#" | cut -f 2 | tr "," "\n" | perl -nl -a -F"\t" -e 'print qq{$F[0]};BEGIN{print qq{S288c};}' > strainlist.tsv

tsv-join --z \
    <(cat ~/data/mrna-structure/vcf/${NAME}_data_SNPs_PARS_mRNA.merge.wild.Chi.tsv | perl -nl -a -F"\t" -e 'print qq{$F[0]\t$F[6]\t$F[7]\t$F[8]\t$F[10]};') \
    -f ../data_SNPs_PARS_mRNA.tsv \
    --key-fields 1 \
    --append-fields 2-63 \
> data_SNPs_PARS_mRNA_all.tsv

rm data_SNPs_PARS_mRNA_filiter.tsv
cat genelist.csv | while read i
    do
        export GENE=${i}
        cat data_SNPs_PARS_mRNA_all.tsv | 
            perl -nl -a -F"\t" -e '
                if ($F[12] eq $ENV{GENE}){
                    my $F = join("\t",@F);
                    print qq{$F};
                }
            ' >> data_SNPs_PARS_mRNA_filiter.tsv
    done

perl ~/Scripts/pars/program/subpop.pl data_SNPs_PARS_mRNA_filiter.tsv strainlist.tsv > subpop.csv
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

```


