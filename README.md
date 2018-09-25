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
- [Phylogeny](#phylogeny)
    - [create protein coding gene list](#create-protein-coding-gene-list)
    - [cut cds alignment](#cut-cds-alignment)
        - [create cds_yml](#create-cds_yml)
        - [cut cds_alignment by cds_yml](#cut-cds_alignment-by-cds_yml)
        - [count cds_alignment proporation in sgd](#count-cds_alignment-proporation-in-sgd)
    - [create gene_phylogeny (n157_nonMosaic)](#create-gene_phylogeny-n157_nonmosaic)
    - [count distance (n157_nonMosaic)](#count-distance-n157_nonmosaic)
- [SNP](#snp)

# Download reference data

## Download PARS10 full site.

```bash
mkdir -p ~/data/mrna-structure/PARS10
cd ~/data/mrna-structure/PARS10

perl ~/Scripts/download/list.pl -u http://genie.weizmann.ac.il/pubs/PARS10/
perl ~/Scripts/download/download.pl -i pubs_PARS10.yml

find . -name "*.gz" | xargs gzip -d
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

# Download strains genomes

## Download S288c (soft-masked) from Ensembl

```bash
mkdir -p ~/data/alignment/egaz/download/S288c
cd ~/data/alignment/egaz/download/S288c

aria2c -x 6 -s 3 -c ftp://ftp.ensembl.org/pub/release-82/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz
aria2c -x 6 -s 3 -c ftp://ftp.ensembl.org/pub/release-82/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.82.gff3.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz S288c.fa
```

## Download strains from NCBI assembly

**DBVPG6044**
```bash
mkdir ~/data/alignment/egaz/download/DBVPG6044
cd ~/data/alignment/egaz/download/DBVPG6044
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/025/GCA_002079025.1_ASM207902v1/GCA_002079025.1_ASM207902v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079025.1_ASM207902v1_genomic.fna.gz DBVPG6044.fa
```

**Y12**
```bash
mkdir ~/data/alignment/egaz/download/Y12
cd ~/data/alignment/egaz/download/Y12
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/058/645/GCA_002058645.1_ASM205864v1/GCA_002058645.1_ASM205864v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002058645.1_ASM205864v1_genomic.fna.gz Y12.fa
```

**SK1**
```bash
mkdir ~/data/alignment/egaz/download/SK1
cd ~/data/alignment/egaz/download/SK1
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/057/885/GCA_002057885.1_ASM205788v1/GCA_002057885.1_ASM205788v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002057885.1_ASM205788v1_genomic.fna.gz SK1.fa
```

**UWOPS03-461.4**
```bash
mkdir ~/data/alignment/egaz/download/UWOPS03_461_4
cd ~/data/alignment/egaz/download/UWOPS03_461_4
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/058/095/GCA_002058095.1_ASM205809v1/GCA_002058095.1_ASM205809v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002058095.1_ASM205809v1_genomic.fna.gz UWOPS03_461_4.fa
```

**YPS128**
```bash
mkdir ~/data/alignment/egaz/download/YPS128
cd ~/data/alignment/egaz/download/YPS128
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/057/995/GCA_002057995.1_ASM205799v1/GCA_002057995.1_ASM205799v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002057995.1_ASM205799v1_genomic.fna.gz YPS128.fa
```

**DBVPG6765**
```bash
mkdir ~/data/alignment/egaz/download/DBVPG6765
cd ~/data/alignment/egaz/download/DBVPG6765
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/057/805/GCA_002057805.1_ASM205780v1/GCA_002057805.1_ASM205780v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002057805.1_ASM205780v1_genomic.fna.gz DBVPG6765.fa
```

**CBS432**
```bash
mkdir ~/data/alignment/egaz/download/CBS432
cd ~/data/alignment/egaz/download/CBS432
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/055/GCA_002079055.1_ASM207905v1/GCA_002079055.1_ASM207905v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079055.1_ASM207905v1_genomic.fna.gz CBS432.fa
```

**N44**
```bash
mkdir ~/data/alignment/egaz/download/N44
cd ~/data/alignment/egaz/download/N44
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/085/GCA_002079085.1_ASM207908v1/GCA_002079085.1_ASM207908v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079085.1_ASM207908v1_genomic.fna.gz N44.fa
```

**YPS138**
```bash
mkdir ~/data/alignment/egaz/download/YPS138
cd ~/data/alignment/egaz/download/YPS138
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/115/GCA_002079115.1_ASM207911v1/GCA_002079115.1_ASM207911v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079115.1_ASM207911v1_genomic.fna.gz YPS138.fa
```
**UFRJ50816**
```bash
mkdir ~/data/alignment/egaz/download/UFRJ50816
cd ~/data/alignment/egaz/download/UFRJ50816
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/145/GCA_002079145.1_ASM207914v1/GCA_002079145.1_ASM207914v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079145.1_ASM207914v1_genomic.fna.gz UFRJ50816.fa
```

**UWOPS91-917.1**
```bash
mkdir ~/data/alignment/egaz/download/UWOPS91_917_1
cd ~/data/alignment/egaz/download/UWOPS91_917_1
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/175/GCA_002079175.1_ASM207917v1/GCA_002079175.1_ASM207917v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079175.1_ASM207917v1_genomic.fna.gz UWOPS91_917_1.fa
```

**EC1118**
```bash
mkdir ~/data/alignment/egaz/download/EC1118
cd ~/data/alignment/egaz/download/EC1118
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/218/975/GCA_000218975.1_ASM21897v1/GCA_000218975.1_ASM21897v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_000218975.1_ASM21897v1_genomic.fna.gz EC1118.fa
```

**Seub**
```bash
mkdir ~/data/alignment/egaz/download/Seub
cd ~/data/alignment/egaz/download/Seub
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/298/625/GCA_001298625.1_SEUB3.0/GCA_001298625.1_SEUB3.0_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_001298625.1_SEUB3.0_genomic.fna.gz Seub.fa
```

## Download strains from NCBI WGS

**scer**
```bash
cd ~/data/alignment/egaz/download
perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/pars/scer_wgs.tsv \
    --fix
bash scer_wgs.rsync.sh
```

**spar**
```bash
cd ~/data/alignment/egaz/download
perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/pars/spar_wgs.tsv \
    --fix
bash spar_wgs.rsync.sh
```

## Download strains from 1002genomes project

```bash
cd ~/data/alignment/egaz/download
wget -c http://1002genomes.u-strasbg.fr/files/1011Assemblies.tar.gz
tar -zxvf 1011Assemblies.tar.gz
```

# RepeatMasker

```bash
cd ~/data/alignment/egaz

egaz prepseq download/S288c/S288c.fa -o S288c -v
gzip -d -c download/S288c/Saccharomyces_cerevisiae.R64-1-1.82.gff3.gz > S288c/chr.gff
egaz masked S288c/*.fa -o S288c/repeat.yml

egaz prepseq \
    download/DBVPG6044/DBVPG6044.fa -o DBVPG6044 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/Y12/Y12.fa -o Y12 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/SK1/SK1.fa -o SK1 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/UWOPS03_461_4/UWOPS03_461_4.fa -o UWOPS03_461_4 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/YPS128/YPS128.fa -o YPS128 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/DBVPG6765/DBVPG6765.fa -o DBVPG6765 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/CBS432/CBS432.fa -o CBS432 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/N44/N44.fa -o N44 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/YPS138/YPS138.fa -o YPS138 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/UFRJ50816/UFRJ50816.fa -o UFRJ50816 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/UWOPS91_917_1/UWOPS91_917_1.fa -o UWOPS91_917_1 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/EC1118/EC1118.fa -o EC1118 \
    --repeatmasker '--species Fungi --parallel 8' -v

mkdir -p ~/data/alignment/egaz/Seub
cd ~/data/alignment/egaz/Seub
RepeatMasker --species Fungi --parallel 16 -xsmall ../download/Seub/Seub.fa
egaz prepseq \
    ../download/Seub/Seub.fa.masked -v

cd ~/data/alignment/egaz
cat download/scer_wgs.csv \
    | grep -v "^prefix" \
    | cut -d',' -f1,3 \
    | uniq \
    | perl -nl -a -F"," -e 'printf qq{egaz prepseq \\\n   download/%s/%s.*.fsa_nt.gz -o %s \\\n   --about 2000000 --repeatmasker " --species Fungi --parallel 8" --min 1000 --gi -v \n}, $F[1], $F[0], $F[1];' > rm_scer_wgs.sh   
bash rm_scer_wgs.sh

cat download/spar_wgs.csv \
    | grep -v "^prefix" \
    | cut -d',' -f1,3 \
    | uniq \
    | perl -nl -a -F"," -e 'printf qq{egaz prepseq \\\n   download/%s/%s.*.fsa_nt.gz -o %s \\\n   --about 2000000 --repeatmasker " --species Fungi --parallel 8" --min 1000 --gi -v \n}, $F[1], $F[0], $F[1];' > rm_spar_wgs.sh
bash rm_spar_wgs.sh

# 1011
cd ~/data/alignment/egaz/download
for file in $(ls ~/data/alignment/egaz/download/GENOMES_ASSEMBLED/*.re.fa);
do

filename=$(basename $file)
dir=$(echo $filename | perl -p -e 's/^([A-Za-z]+).+/$1/;')

if [ -d ../$dir ];
then echo -n;
else
cd ~/data/alignment/egaz/download
cat $file | perl -nl -e '

if (m/^>([A-Za-z]+)/){

my $dir = $1;
mkdir ("../$dir") unless (-d "../$dir");
last;

}
'
cp -rf $file ../$dir
sed -i".bak" "s/-/_/g" ../$dir/*.re.fa
sed -i".bak" "s/\./_/g" ../$dir/*.re.fa

faops filter -a 1000 ../$dir/*.re.fa ../$dir/$dir.fasta
rm -rf ../$dir/*.re.fa
RepeatMasker --species Fungi --parallel 16 -xsmall ../$dir/$dir.fasta

cd ~/data/alignment/egaz/$dir
egaz prepseq \
    ../$dir/$dir.fasta.masked -v
fi
done
```

# Align

## Sanger
```bash
mkdir -p ~/data/mrna-structure/alignment/scer_wgs
cd ~/data/mrna-structure/alignment/scer_wgs
ln -s ~/data/alignment/egaz .

egaz template \
    egaz/S288c egaz/EC1118 egaz/Kyokai_no_7 egaz/RM11_1a egaz/Sigma1278b egaz/T7 egaz/YJM789 egaz/Spar \
    --multi -o multi8/ \
    --rawphylo --order --parallel 16 -v
bash multi8/1_pair.sh
bash multi8/2_rawphylo.sh
bash multi8/3_multi.sh

egaz template \
    egaz/S288c egaz/EC1118 egaz/Kyokai_no_7 egaz/RM11_1a egaz/Sigma1278b egaz/T7 egaz/YJM789 egaz/Spar \
    --multi -o multi8/ \
    --multiname Scer_n7_Spar --order --tree multi8/Results/multi8.nwk --outgroup Spar \
    --vcf --aligndb \
    --parallel 16 -v

bash multi8/3_multi.sh
bash multi8/6_chr_length.sh
bash multi8/7_multi_aligndb.sh

mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n7_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/multi8/Scer_n7_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/multi8/Results/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/multi8/Results/chr_length.csv \
    -lt 1000 --parallel 16 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n7_Spar -r 1-60
```

## PacBio
```bash
mkdir -p ~/data/mrna-structure/alignment/scer_wgs
cd ~/data/mrna-structure/alignment/scer_wgs
ln -s ~/data/alignment/egaz .

egaz template \
    egaz/S288c egaz/DBVPG6044 egaz/UWOPS03_461_4 egaz/Y12 egaz/SK1 egaz/YPS128 egaz/DBVPG6765 egaz/Spar \
    --multi -o multi8p/ \
    --rawphylo --order --parallel 16 -v
bash multi8p/1_pair.sh
bash multi8p/2_rawphylo.sh
bash multi8p/3_multi.sh

egaz template \
    egaz/S288c egaz/DBVPG6044 egaz/UWOPS03_461_4 egaz/Y12 egaz/SK1 egaz/YPS128 egaz/DBVPG6765 egaz/Spar \
    --multi -o multi8p/ \
    --multiname Scer_n7p_Spar --order --tree multi8p/Results/multi8p.nwk --outgroup Spar \
    --vcf --aligndb \
    --parallel 16 -v

bash multi8p/3_multi.sh
bash multi8p/6_chr_length.sh
bash multi8p/7_multi_aligndb.sh

mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n7p_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/multi8p/Scer_n7p_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/multi8p/Results/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/multi8p/Results/chr_length.csv \
    -lt 1000 --parallel 16 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n7p_Spar -r 1-60
```

## Illumina 
```bash
# n128_Spar

mkdir -p ~/data/mrna-structure/alignment/scer_wgs
cd ~/data/mrna-structure/alignment/scer_wgs
ln -s ~/data/alignment/egaz .

egaz template \
    egaz/S288c egaz/beer001 egaz/beer003 egaz/beer004 egaz/beer005 egaz/beer006 egaz/beer007 egaz/beer008 egaz/beer009 egaz/beer010 egaz/beer011 egaz/beer012 egaz/beer013 egaz/beer014 egaz/beer015 egaz/beer016 egaz/beer020 egaz/beer021 egaz/beer022 egaz/beer023 egaz/beer024 egaz/beer025 egaz/beer026 egaz/beer027 egaz/beer028 egaz/beer029 egaz/beer030 egaz/beer031 egaz/beer032 egaz/beer033 egaz/beer034 egaz/beer036 egaz/beer037 egaz/beer038 egaz/beer040 egaz/beer041 egaz/beer043 egaz/beer044 egaz/beer045 egaz/beer046 egaz/beer047 egaz/beer048 egaz/beer049 egaz/beer050 egaz/beer051 egaz/beer052 egaz/beer053 egaz/beer054 egaz/beer055 egaz/beer056 egaz/beer059 egaz/beer061 egaz/beer062 egaz/beer063 egaz/beer064 egaz/beer065 egaz/beer066 egaz/beer067 egaz/beer068 egaz/beer069 egaz/beer070 egaz/beer071 egaz/beer073 egaz/beer075 egaz/beer076 egaz/beer077 egaz/beer078 egaz/beer079 egaz/beer080 egaz/beer081 egaz/beer082 egaz/beer083 egaz/beer084 egaz/beer085 egaz/beer086 egaz/beer087 egaz/beer088 egaz/beer089 egaz/beer090 egaz/beer091 egaz/beer092 egaz/beer094 egaz/beer095 egaz/beer096 egaz/beer097 egaz/beer098 egaz/beer099 egaz/beer100 egaz/beer101 egaz/beer102 egaz/bioethanol001 egaz/bioethanol003 egaz/bioethanol004 egaz/bread001 egaz/bread002 egaz/bread003 egaz/bread004 egaz/sake001 egaz/sake002 egaz/sake003 egaz/sake004 egaz/sake005 egaz/sake006 egaz/sake007 egaz/spirits001 egaz/spirits002 egaz/spirits003 egaz/spirits004 egaz/spirits005 egaz/spirits011 egaz/wine001 egaz/wine003 egaz/wine004 egaz/wine005 egaz/wine006 egaz/wine007 egaz/wine009 egaz/wine010 egaz/wine011 egaz/wine012 egaz/wine013 egaz/wine014 egaz/wine015 egaz/wine017 egaz/wine018 egaz/wild005 egaz/wild006 egaz/wild007 egaz/Spar \
    --multi -o multi128_Spar/ \
    --rawphylo --order --parallel 16 -v
bash multi128_Spar/1_pair.sh
bash multi128_Spar/2_rawphylo.sh
bash multi128_Spar/3_multi.sh

egaz template \
    egaz/S288c egaz/beer001 egaz/beer003 egaz/beer004 egaz/beer005 egaz/beer006 egaz/beer007 egaz/beer008 egaz/beer009 egaz/beer010 egaz/beer011 egaz/beer012 egaz/beer013 egaz/beer014 egaz/beer015 egaz/beer016 egaz/beer020 egaz/beer021 egaz/beer022 egaz/beer023 egaz/beer024 egaz/beer025 egaz/beer026 egaz/beer027 egaz/beer028 egaz/beer029 egaz/beer030 egaz/beer031 egaz/beer032 egaz/beer033 egaz/beer034 egaz/beer036 egaz/beer037 egaz/beer038 egaz/beer040 egaz/beer041 egaz/beer043 egaz/beer044 egaz/beer045 egaz/beer046 egaz/beer047 egaz/beer048 egaz/beer049 egaz/beer050 egaz/beer051 egaz/beer052 egaz/beer053 egaz/beer054 egaz/beer055 egaz/beer056 egaz/beer059 egaz/beer061 egaz/beer062 egaz/beer063 egaz/beer064 egaz/beer065 egaz/beer066 egaz/beer067 egaz/beer068 egaz/beer069 egaz/beer070 egaz/beer071 egaz/beer073 egaz/beer075 egaz/beer076 egaz/beer077 egaz/beer078 egaz/beer079 egaz/beer080 egaz/beer081 egaz/beer082 egaz/beer083 egaz/beer084 egaz/beer085 egaz/beer086 egaz/beer087 egaz/beer088 egaz/beer089 egaz/beer090 egaz/beer091 egaz/beer092 egaz/beer094 egaz/beer095 egaz/beer096 egaz/beer097 egaz/beer098 egaz/beer099 egaz/beer100 egaz/beer101 egaz/beer102 egaz/bioethanol001 egaz/bioethanol003 egaz/bioethanol004 egaz/bread001 egaz/bread002 egaz/bread003 egaz/bread004 egaz/sake001 egaz/sake002 egaz/sake003 egaz/sake004 egaz/sake005 egaz/sake006 egaz/sake007 egaz/spirits001 egaz/spirits002 egaz/spirits003 egaz/spirits004 egaz/spirits005 egaz/spirits011 egaz/wine001 egaz/wine003 egaz/wine004 egaz/wine005 egaz/wine006 egaz/wine007 egaz/wine009 egaz/wine010 egaz/wine011 egaz/wine012 egaz/wine013 egaz/wine014 egaz/wine015 egaz/wine017 egaz/wine018 egaz/wild005 egaz/wild006 egaz/wild007 egaz/Spar \
    --multi -o multi128_Spar/ \
    --multiname Scer_n128_Spar --order --tree multi128_Spar/Results/multi128_Spar.nwk --outgroup Spar \
    --vcf --aligndb \
    --parallel 16 -v

bash multi128_Spar/3_multi.sh
bash multi128_Spar/6_chr_length.sh
bash multi128_Spar/7_multi_aligndb.sh

mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n128_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/multi128_Spar/Scer_n128_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/multi128_Spar/Results/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/multi128_Spar/Results/chr_length.csv \
    -lt 1000 --parallel 16 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n128_Spar -r 1-60

# n128_Seub

mkdir -p ~/data/mrna-structure/alignment/scer_wgs
cd ~/data/mrna-structure/alignment/scer_wgs
ln -s ~/data/alignment/egaz .

egaz template \
    egaz/S288c egaz/beer001 egaz/beer003 egaz/beer004 egaz/beer005 egaz/beer006 egaz/beer007 egaz/beer008 egaz/beer009 egaz/beer010 egaz/beer011 egaz/beer012 egaz/beer013 egaz/beer014 egaz/beer015 egaz/beer016 egaz/beer020 egaz/beer021 egaz/beer022 egaz/beer023 egaz/beer024 egaz/beer025 egaz/beer026 egaz/beer027 egaz/beer028 egaz/beer029 egaz/beer030 egaz/beer031 egaz/beer032 egaz/beer033 egaz/beer034 egaz/beer036 egaz/beer037 egaz/beer038 egaz/beer040 egaz/beer041 egaz/beer043 egaz/beer044 egaz/beer045 egaz/beer046 egaz/beer047 egaz/beer048 egaz/beer049 egaz/beer050 egaz/beer051 egaz/beer052 egaz/beer053 egaz/beer054 egaz/beer055 egaz/beer056 egaz/beer059 egaz/beer061 egaz/beer062 egaz/beer063 egaz/beer064 egaz/beer065 egaz/beer066 egaz/beer067 egaz/beer068 egaz/beer069 egaz/beer070 egaz/beer071 egaz/beer073 egaz/beer075 egaz/beer076 egaz/beer077 egaz/beer078 egaz/beer079 egaz/beer080 egaz/beer081 egaz/beer082 egaz/beer083 egaz/beer084 egaz/beer085 egaz/beer086 egaz/beer087 egaz/beer088 egaz/beer089 egaz/beer090 egaz/beer091 egaz/beer092 egaz/beer094 egaz/beer095 egaz/beer096 egaz/beer097 egaz/beer098 egaz/beer099 egaz/beer100 egaz/beer101 egaz/beer102 egaz/bioethanol001 egaz/bioethanol003 egaz/bioethanol004 egaz/bread001 egaz/bread002 egaz/bread003 egaz/bread004 egaz/sake001 egaz/sake002 egaz/sake003 egaz/sake004 egaz/sake005 egaz/sake006 egaz/sake007 egaz/spirits001 egaz/spirits002 egaz/spirits003 egaz/spirits004 egaz/spirits005 egaz/spirits011 egaz/wine001 egaz/wine003 egaz/wine004 egaz/wine005 egaz/wine006 egaz/wine007 egaz/wine009 egaz/wine010 egaz/wine011 egaz/wine012 egaz/wine013 egaz/wine014 egaz/wine015 egaz/wine017 egaz/wine018 egaz/wild005 egaz/wild006 egaz/wild007 egaz/Seub \
    --multi -o multi128_Seub/ \
    --rawphylo --order --parallel 16 -v
bash multi128_Seub/1_pair.sh
bash multi128_Seub/2_rawphylo.sh
bash multi128_Seub/3_multi.sh

egaz template \
    egaz/S288c egaz/beer001 egaz/beer003 egaz/beer004 egaz/beer005 egaz/beer006 egaz/beer007 egaz/beer008 egaz/beer009 egaz/beer010 egaz/beer011 egaz/beer012 egaz/beer013 egaz/beer014 egaz/beer015 egaz/beer016 egaz/beer020 egaz/beer021 egaz/beer022 egaz/beer023 egaz/beer024 egaz/beer025 egaz/beer026 egaz/beer027 egaz/beer028 egaz/beer029 egaz/beer030 egaz/beer031 egaz/beer032 egaz/beer033 egaz/beer034 egaz/beer036 egaz/beer037 egaz/beer038 egaz/beer040 egaz/beer041 egaz/beer043 egaz/beer044 egaz/beer045 egaz/beer046 egaz/beer047 egaz/beer048 egaz/beer049 egaz/beer050 egaz/beer051 egaz/beer052 egaz/beer053 egaz/beer054 egaz/beer055 egaz/beer056 egaz/beer059 egaz/beer061 egaz/beer062 egaz/beer063 egaz/beer064 egaz/beer065 egaz/beer066 egaz/beer067 egaz/beer068 egaz/beer069 egaz/beer070 egaz/beer071 egaz/beer073 egaz/beer075 egaz/beer076 egaz/beer077 egaz/beer078 egaz/beer079 egaz/beer080 egaz/beer081 egaz/beer082 egaz/beer083 egaz/beer084 egaz/beer085 egaz/beer086 egaz/beer087 egaz/beer088 egaz/beer089 egaz/beer090 egaz/beer091 egaz/beer092 egaz/beer094 egaz/beer095 egaz/beer096 egaz/beer097 egaz/beer098 egaz/beer099 egaz/beer100 egaz/beer101 egaz/beer102 egaz/bioethanol001 egaz/bioethanol003 egaz/bioethanol004 egaz/bread001 egaz/bread002 egaz/bread003 egaz/bread004 egaz/sake001 egaz/sake002 egaz/sake003 egaz/sake004 egaz/sake005 egaz/sake006 egaz/sake007 egaz/spirits001 egaz/spirits002 egaz/spirits003 egaz/spirits004 egaz/spirits005 egaz/spirits011 egaz/wine001 egaz/wine003 egaz/wine004 egaz/wine005 egaz/wine006 egaz/wine007 egaz/wine009 egaz/wine010 egaz/wine011 egaz/wine012 egaz/wine013 egaz/wine014 egaz/wine015 egaz/wine017 egaz/wine018 egaz/wild005 egaz/wild006 egaz/wild007 egaz/Seub \
    --multi -o multi128_Seub/ \
    --multiname Scer_n128_Seub --order --tree multi128_Seub/Results/multi128_Seub.nwk --outgroup Seub \
    --vcf --aligndb \
    --parallel 16 -v

bash multi128_Seub/3_multi.sh
bash multi128_Seub/6_chr_length.sh
bash multi128_Seub/7_multi_aligndb.sh

mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n128_Seub \
    -da ~/data/mrna-structure/alignment/scer_wgs/multi128_Seub/Scer_n128_Seub_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/multi128_Seub/Results/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/multi128_Seub/Results/chr_length.csv \
    -lt 1000 --parallel 16 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n128_Seub -r 1-60
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

## Extract gene-list and snp-codon-list n128

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n128_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n128_Spar.mvar.gene_list.csv


cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Seub.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n128_Seub.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Seub.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n128_Seub.mvar.gene_list.csv
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

## SNPs and indels n128

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n128_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n128_Spar.indel.pos.txt


Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Seub.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n128_Seub.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Seub.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n128_Seub.indel.pos.txt
```

# Blast

Prepare a combined fasta file of yeast genome and blast genes against
the genome.

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

# Phylogeny

## create protein coding gene list

```bash
mkdir -p ~/data/mrna-structure/phylogeny
cd ~/data/mrna-structure/phylogeny

# sgd/saccharomyces_cerevisiae.gff → protein coding gene list

perl ~/Scripts/pars/program/protein_coding_list.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list.csv

perl ~/Scripts/pars/program/protein_coding_list_range.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list_range.csv

perl ~/Scripts/pars/program/protein_coding_list_range_chr.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list_range_chr.csv
```

## cut mRNA alignment

### create mRNA_yml

```bash
cd ~/data/mrna-structure/phylogeny
mkdir -p ~/data/mrna-structure/phylogeny/gene_mRNA_yml

#cut mRNA in Scer.gff
perl ~/Scripts/pars/program/cut_mRNA_yml.pl --file protein_coding_list_range_chr.csv --output gene_mRNA_yml
```

### cut alignment by mRNA_yml

```bash

export NAME=Scer_n7_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/multi8/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_mRNA
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_mRNA
cat ../protein_coding_list.csv |
   parallel --line-buffer -j 8 '
       fasops slice ../${NAME}_refined/species.fas ../gene_mRNA_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n7p_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/multi8p/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_mRNA
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_mRNA
cat ../protein_coding_list.csv |
   parallel --line-buffer -j 8 '
       fasops slice ../${NAME}_refined/species.fas ../gene_mRNA_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n128_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/multi128_Spar/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_mRNA
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_mRNA
cat ../protein_coding_list.csv |
   parallel --line-buffer -j 8 '
   	   fasops slice ../${NAME}_refined/species.fas ../gene_mRNA_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n128_Seub
cp -rf ~/data/mrna-structure/alignment/scer_wgs/multi128_Seub/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_mRNA
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_mRNA
cat ../protein_coding_list.csv |
   parallel --line-buffer -j 8 '
   	   fasops slice ../${NAME}_refined/species.fas ../gene_mRNA_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME
```

### count mRNA_alignment proporation in sgd

```bash

export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/phylogeny
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_mRNA --output ${NAME}_gene_range.csv
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/phylogeny
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_mRNA --output ${NAME}_gene_range.csv
unset NAME
 
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/phylogeny
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_mRNA --output ${NAME}_gene_range.csv
unset NAME

export NAME=Scer_n128_Seub
cd ~/data/mrna-structure/phylogeny
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_mRNA --output ${NAME}_gene_range.csv
unset NAME
```

```bash
#生成alignment_proporation_1.list

export NAME=Scer_n7_Spar
Rscript ~/Scripts/pars/program/${NAME}_distance_processed.R
unset NAME

export NAME=Scer_n7p_Spar
Rscript ~/Scripts/pars/program/${NAME}_distance_processed.R
unset NAME

export NAME=Scer_n128_Spar
Rscript ~/Scripts/pars/program/${NAME}_distance_processed.R
unset NAME

export NAME=Scer_n128_Seub
Rscript ~/Scripts/pars/program/${NAME}_distance_processed.R
unset NAME
```

# SNP

## count per gene GC content

```bash

export NAME=Scer_n7_Spar
mkdir -p ~/data/mrna-structure/result/$NAME
cd ~/data/mrna-structure/result/$NAME
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --varfold ~/data/mrna-structure/process/$NAME.gene_variation.fold_class.tsv --output $NAME.gene_variation.fold_class.csv
unset NAME

export NAME=Scer_n7p_Spar
mkdir -p ~/data/mrna-structure/result/$NAME
cd ~/data/mrna-structure/result/$NAME
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --varfold ~/data/mrna-structure/process/$NAME.gene_variation.fold_class.tsv --output $NAME.gene_variation.fold_class.csv
unset NAME

export NAME=Scer_n128_Spar
mkdir -p ~/data/mrna-structure/result/$NAME
cd ~/data/mrna-structure/result/$NAME
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --varfold ~/data/mrna-structure/process/$NAME.gene_variation.fold_class.tsv --output $NAME.gene_variation.fold_class.csv
unset NAME

export NAME=Scer_n128_Seub
mkdir -p ~/data/mrna-structure/result/$NAME
cd ~/data/mrna-structure/result/$NAME
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --varfold ~/data/mrna-structure/process/$NAME.gene_variation.fold_class.tsv --output $NAME.gene_variation.fold_class.csv
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

export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/result/$NAME
Rscript ~/Scripts/pars/program/${NAME}_stat_SNPs.R
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/result/$NAME
Rscript ~/Scripts/pars/program/${NAME}_stat_SNPs.R
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/$NAME
Rscript ~/Scripts/pars/program/${NAME}_stat_SNPs.R
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

export NAME=Scer_n128_Seub
cd ~/data/mrna-structure/result/$NAME
Rscript ~/Scripts/pars/program/${NAME}_stat_SNPs.R
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME
```

## count A/T <->G/C

```bash

export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/result/$NAME
mkdir -p ~/data/mrna-structure/result/$NAME/freq_each
Rscript ~/Scripts/pars/program/${NAME}_count_AT_GC.R
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_cds_stat.csv --output freq_each/PARS_cds_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_utr_stat.csv --output freq_each/PARS_utr_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_syn_stat.csv --output freq_each/PARS_syn_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_nsy_stat.csv --output freq_each/PARS_nsy_stat_chi_square.csv
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/result/$NAME
mkdir -p ~/data/mrna-structure/result/$NAME/freq_each
Rscript ~/Scripts/pars/program/${NAME}_count_AT_GC.R
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_cds_stat.csv --output freq_each/PARS_cds_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_utr_stat.csv --output freq_each/PARS_utr_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_syn_stat.csv --output freq_each/PARS_syn_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_nsy_stat.csv --output freq_each/PARS_nsy_stat_chi_square.csv
unset NAME

export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/$NAME
mkdir -p ~/data/mrna-structure/result/$NAME/freq_each
mkdir -p ~/data/mrna-structure/result/$NAME/freq_10
Rscript ~/Scripts/pars/program/${NAME}_count_AT_GC.R
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_cds_stat_freq_10.csv --output freq_10/PARS_cds_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_utr_stat_freq_10.csv --output freq_10/PARS_utr_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_syn_stat_freq_10.csv --output freq_10/PARS_syn_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_nsy_stat_freq_10.csv --output freq_10/PARS_nsy_stat_freq_10_chi_square.csv
unset NAME

export NAME=Scer_n128_Seub
cd ~/data/mrna-structure/result/$NAME
mkdir -p ~/data/mrna-structure/result/$NAME/freq_each
mkdir -p ~/data/mrna-structure/result/$NAME/freq_10
Rscript ~/Scripts/pars/program/${NAME}_count_AT_GC.R
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_cds_stat_freq_10.csv --output freq_10/PARS_cds_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_utr_stat_freq_10.csv --output freq_10/PARS_utr_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_syn_stat_freq_10.csv --output freq_10/PARS_syn_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_nsy_stat_freq_10.csv --output freq_10/PARS_nsy_stat_freq_10_chi_square.csv
unset NAME

```

## count stem length selection cds
```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/$NAME 
mkdir -p freq_10/stem_length
perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --origin data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds_pos.csv
Rscript ~/Scripts/pars/program/$NAME_count_AT_GC_gene_trait.R
unset NAME

```

## count stem length selection mRNA
```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/$NAME 
mkdir -p freq_10/stem_length
perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --origin data_SNPs_PARS_mRNA.csv --output data_SNPs_PARS_mRNA_pos.csv
Rscript ~/Scripts/pars/program/$NAME_count_AT_GC_gene_trait.R
unset NAME

```

```bash
export NAME=Scer_n128_Seub
cd ~/data/mrna-structure/result/$NAME 
mkdir -p freq_10/stem_length
perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --origin data_SNPs_PARS_mRNA.csv --output data_SNPs_PARS_mRNA_pos.csv
Rscript ~/Scripts/pars/program/$NAME_count_AT_GC_gene_trait.R
unset NAME

```

## count GO KEGG
```bash
export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/$NAME 
mkdir -p freq_10/GO
mkdir -p freq_10/KEGG
Rscript ~/Scripts/pars/program/$NAME_count_AT_GC_GO.R
Rscript ~/Scripts/pars/program/$NAME_count_AT_GC_KEGG.R
unset NAME

```

```bash
export NAME=Scer_n128_Seub
cd ~/data/mrna-structure/result/$NAME 
mkdir -p freq_10/GO
mkdir -p freq_10/KEGG
Rscript ~/Scripts/pars/program/Scer_n128_Seub_count_AT_GC_GO.R 
Rscript ~/Scripts/pars/program/Scer_n128_Seub_count_AT_GC_KEGG.R
unset NAME

```

## count per gene cds_utr
```bash
export NAME=Scer_n128_Spar
perl ~/Scripts/pars/program/count_cut_range.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --cut ~/data/mrna-structure/process/sce_cds.yml --output stem_loop_cds_length.csv 
perl ~/Scripts/pars/program/count_cut_range.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --cut ~/data/mrna-structure/process/sce_utr.yml --output stem_loop_utr_length.csv
perl ~/Scripts/pars/program/count_per_gene_ACGT_percent.pl --file data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds_per_gene_ATGC.csv
perl ~/Scripts/pars/program/count_per_gene_ACGT_percent.pl --file data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr_per_gene_ATGC.csv
Rscript ~/Scripts/pars/program/$NAME_cds_utr.R
unset NAME
```

# vcf
```bash
mkdir -p ~/data/mrna-structure/vcf
cd ~/data/mrna-structure/vcf
wget -c http://1002genomes.u-strasbg.fr/files/1011Matrix.gvcf.gz
gzip -d 1011Matrix.gvcf.gz

#ln -s /Volumes/Backup/yumh/data/vcf/1011Matrix.gvcf.gz .
#mkdir -p ~/data/mrna-structure/vcf/1011Matrix.gvcf
#cd 1011Matrix.gvcf/
#ln -s /Volumes/Backup/yumh/data/vcf/1011Matrix.gvcf/1011Matrix.gvcf .

# 1011
cd ~/data/vcf/1011Matrix.gvcf
perl ~/Scripts/pars/program/vcf.cut.pl --file 1011Matrix.gvcf --output 1011Matrix.tsv
perl ~/Scripts/pars/program/vcf.merge_pre.pl --file ~/data/mrna-structure/result/Scer_n128_Spar/data_SNPs_PARS_cds.csv --output chr1.pars.tsv
perl extract.pl --file chr1.tsv --output chr1.ext.tsv
Rscript merge.R
perl merge_pro.pl --file chr1.merge.tsv --output chr1.merge.pro.tsv


# wild.strains in 1011
cd ~/yumh/data/vcf/1011Matrix.gvcf
bcftools view chr1.gvcf -s CCL,BBQ,BBS,BFP,BTG,CLC,CLB,CLD,BAM,BAQ,BAG,BAH,BAL,AMH,CEG,CEI,CCQ,CCR,CCS,BAK,BAI,ACQ,CCN,CDL,SACE_YCR,BMA,AKM,BMB,BMC,SACE_MAL,SACE_YCY,BAN,BAP,CMP,CCH,ACC,CCC,CCD,CCE,CCF,CCG,CCI,CMQ,CDF,CDG,CDH,CDI,AVI,ACD,ANF,ANH,ANC,ANE,ANG,AND,ANK,ANI,AKN,SACE_YBS,SACE_YCU | bcftools +fill-tags -o chr1.subset.vcf
perl cut.pl --file chr1.subset.vcf --output chr1.subset.tsv
perl merge_pre.pl --file ~/data/mrna-structure/result/Scer_n128_Spar/data_SNPs_PARS_cds.csv --output chr1.pars.tsv
perl extract.pl --file chr1.subset.tsv --output chr1.subset.ext.tsv
Rscript merge.subset.R
perl merge_pro.pl --file chr1.merge.subset.tsv --output chr1.merge.subset.pro.tsv
```