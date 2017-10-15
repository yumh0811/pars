# Processing Yeast PARS Data

[TOC level=1-3]: # " "
- [Processing Yeast PARS Data](#processing-yeast-pars-data)
- [Download reference data](#download-reference-data)
    - [Download PARS10 full site.](#download-pars10-full-site)
    - [Download S288c annotation data from ensembl by rsync](#download-s288c-annotation-data-from-ensembl-by-rsync)
    - [SGD](#sgd)
    - [mRNA levels](#mrna-levels)
    - [ess, rich/minimal and chem](#ess-richminimal-and-chem)
    - [Recombination rates](#recombination-rates)
    - [Protein-protein interactions](#protein-protein-interactions)
- [Other strains and outgroups](#other-strains-and-outgroups)
    - [Sanger (NCBI WGS)](#sanger-ncbi-wgs)
    - [Illumina (NCBI ASSEMBLY)](#illumina-ncbi-assembly)
- [Plans of alignments](#plans-of-alignments)
- [Build alignDB for multiple genomes](#build-aligndb-for-multiple-genomes)
    - [Extract `gene_list` and `snp_codon_list`](#extract-gene-list-and-snp-codon-list)
    - [SNPs and indels](#snps-and-indels)
- [Blast](#blast)
- [Features](#features)
- [Real Processing](#real-processing)
- [Stats](#stats)
- [Pack all things up](#pack-all-things-up)


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

# Other strains and outgroups

* `withncbi/db/`: taxonomy database
* `withncbi/pop/`: scer_wgs alignments

## Sanger (NCBI WGS)

| Strain                       | Taxonomy ID | Sequencing Technology | Total length (bp) |
|:-----------------------------|:------------|:----------------------|:------------------|
| *S. cerevisiae* S288c        | 559292      | Sanger                | 12,157,105        |
| *S. cerevisiae* EC1118       | 643680      | 6x Sanger; 17.6x 454  | 11,659,512        |
| *S. cerevisiae* Kyokai no. 7 | 721032      | 9.1x Sanger           | 12,370,866        |
| *S. cerevisiae* RM11-1a      | 285006      | 10x Sanger            | 11,675,031        |
| *S. cerevisiae* Sigma1278b   | 658763      | 45x Sanger/Illumina   | 11,906,055        |
| *S. cerevisiae* T7           | 929585      | 25.4x Sanger/454      | 11,758,843        |
| *S. cerevisiae* YJM789       | 307796      | 10x Sanger            | 11,990,995        |
| *S. paradoxus* NRRL Y-17217  | 226125      | 7.7x Sanger           | 11,872,617        |

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
    -p -f ~/Scripts/pars/scer_100.seq.csv

```

## Illumina (NCBI WGS_157)

```bash
mkdir -p ~/data/mrna-structure/GENOMES
cd ~/data/mrna-structure/GENOMES

perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/pars/scer_wgs_157.tsv \
    --fix -a \
    -o WGS

aria2c -UWget -x 6 -s 3 -c -i WGS/scer_wgs_157.url.txt

find WGS -name "*.gz" | xargs gzip -t
```


# Plans of alignments

## Illumina (NCBI ASSEMBLY)

```bash
# create downloaded genome list
cat ~/Scripts/pars/scer_100.seq.csv \
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
    --download "name=EC1118;taxon=643680" \
    --download "name=YJM993;taxon=1294331" \
    --download "name=YJM195;taxon=1294305" \
    --download "name=YJM270;taxon=1294308" \
    --download "name=YJM470;taxon=1294313" \
    --download "name=YJM683;taxon=1294320" \
    --download "name=YJM689;taxon=1294321" \
    --download "name=YJM693;taxon=1294322" \
    --download "name=YJM1248;taxon=1294340" \
    --download "name=YJM1273;taxon=1294343" \
    --download "name=YJM1385;taxon=1294357" \
    --download "name=YJM1388;taxon=1294360" \
    --download "name=YJM1389;taxon=1294361" \
    --download "name=YJM1399;taxon=1294362" \
    --download "name=YJM1402;taxon=1294365" \
    --download "name=YJM1418;taxon=1294368" \
    --download "name=YJM1439;taxon=1294372" \
    --download "name=YJM1443;taxon=1294373" \
    --download "name=YJM1444;taxon=1294374" \
    --download "name=YJM1447;taxon=1294375" \
    --download "name=YJM1460;taxon=1294377" \
    --download "name=YJM1549;taxon=1294384" \
    --download "name=YJM1573;taxon=1294385" \
    --download "name=YJM1592;taxon=1294387" \
    --download "name=YJM244;taxon=1294306" \
    --download "name=YJM1083;taxon=1292971" \
    --download "name=YJM1129;taxon=1293430" \
    --download "name=YJM193;taxon=1294304" \
    --download "name=YJM248;taxon=1294307" \
    --download "name=YJM320;taxon=947042" \
    --download "name=YJM326;taxon=468558" \
    --download "name=YJM428;taxon=947044" \
    --download "name=YJM451;taxon=502869" \
    --download "name=YJM453;taxon=1294311" \
    --download "name=YJM456;taxon=1294312" \
    --download "name=YJM541;taxon=1294314" \
    --download "name=YJM555;taxon=1294316" \
    --download "name=YJM627;taxon=1294317" \
    --download "name=YJM681;taxon=1294318" \
    --download "name=YJM969;taxon=1294323" \
    --download "name=YJM972;taxon=1294324" \
    --download "name=YJM978;taxon=1294326" \
    --download "name=YJM981;taxon=1294327" \
    --download "name=YJM984;taxon=1294328" \
    --download "name=YJM996;taxon=1294332" \
    --download "name=YJM1190;taxon=1294334" \
    --download "name=YJM1199;taxon=1294335" \
    --download "name=YJM1202;taxon=1294336" \
    --download "name=YJM1250;taxon=1294341" \
    --download "name=YJM1311;taxon=1294346" \
    --download "name=YJM1332;taxon=1294348" \
    --download "name=YJM1338;taxon=1294350" \
    --download "name=YJM1341;taxon=1294351" \
    --download "name=YJM1355;taxon=1294353" \
    --download "name=YJM1356;taxon=1294354" \
    --download "name=YJM1383;taxon=1294356" \
    --download "name=YJM1386;taxon=1294358" \
    --download "name=YJM1400;taxon=1294363" \
    --download "name=YJM1401;taxon=1294364" \
    --download "name=YJM1415;taxon=1294366" \
    --download "name=YJM1433;taxon=1294370" \
    --download "name=YJM1450;taxon=1294376" \
    --download "name=YJM1463;taxon=1294378" \
    --download "name=YJM1478;taxon=1294380" \
    --download "name=YJM1479;taxon=1294381" \
    --download "name=YJM1526;taxon=1294382" \
    --download "name=YJM1527;taxon=1294383" \
    --download "name=YJM1574;taxon=1294386" \
    --download "name=YJM1615;taxon=1294388" \
    --download "name=YJM1304;taxon=1294344" \
    --download "name=YJM1434;taxon=1294371" \
    --download "name=YJM1078;taxon=1296266" \
    --download "name=YJM450;taxon=1294310" \
    --download "name=YJM990;taxon=1294330" \
    --download "name=YJM1242;taxon=1294338" \
    --download "name=YJM1244;taxon=1294339" \
    --download "name=YJM1307;taxon=1294345" \
    --download "name=YJM1336;taxon=1294349" \
    --download "name=YJM1381;taxon=1294355" \
    --download "name=YJM1387;taxon=1294359" \
    --download "name=YJM1419;taxon=1294369" \
    --download "name=YJM1477;taxon=1294379" \
    --download "name=YJM189;taxon=1294303" \
    --download "name=YJM271;taxon=1294309" \
    --download "name=YJM554;taxon=1294315" \
    --download "name=YJM682;taxon=1294319" \
    --download "name=YJM975;taxon=1294325" \
    --download "name=YJM987;taxon=1294329" \
    --download "name=YJM1133;taxon=1294333" \
    --download "name=YJM1208;taxon=1294337" \
    --download "name=YJM1252;taxon=1294342" \
    --download "name=YJM1326;taxon=1294347" \
    --download "name=YJM1342;taxon=1294352" \
    --download "name=YJM1417;taxon=1294367" \
    --plan 'name=Scer_n7_pop;t=S288c;qs=EC1118,Kyokai_no_7,RM11_1a,Sigma1278b,T7,YJM789' \
    --plan 'name=Scer_n7_Spar;t=S288c;qs=EC1118,Kyokai_no_7,RM11_1a,Sigma1278b,T7,YJM789,Spar;o=Spar' \
    --plan 'name=Scer_n94_pop;t=S288c;qs=YJM993,YJM1078,YJM195,YJM270,YJM470,YJM683,YJM689,YJM693,YJM1248,YJM1252,YJM1273,YJM1342,YJM1385,YJM1387,YJM1388,YJM1389,YJM1399,YJM1402,YJM1418,YJM1439,YJM1443,YJM1444,YJM1447,YJM1460,YJM1549,YJM1573,YJM1592,YJM244,YJM1083,YJM1129,YJM189,YJM193,YJM248,YJM271,YJM320,YJM326,YJM428,YJM450,YJM451,YJM453,YJM456,YJM541,YJM554,YJM555,YJM627,YJM681,YJM682,YJM969,YJM972,YJM975,YJM978,YJM981,YJM984,YJM987,YJM990,YJM996,YJM1133,YJM1190,YJM1199,YJM1202,YJM1208,YJM1242,YJM1244,YJM1250,YJM1307,YJM1311,YJM1326,YJM1332,YJM1336,YJM1338,YJM1341,YJM1355,YJM1356,YJM1381,YJM1383,YJM1386,YJM1400,YJM1401,YJM1415,YJM1417,YJM1419,YJM1433,YJM1450,YJM1463,YJM1477,YJM1478,YJM1479,YJM1526,YJM1527,YJM1574,YJM1615' \
    --plan 'name=Scer_n94_Spar;t=S288c;qs=YJM993,YJM1078,YJM195,YJM270,YJM470,YJM683,YJM689,YJM693,YJM1248,YJM1252,YJM1273,YJM1342,YJM1385,YJM1387,YJM1388,YJM1389,YJM1399,YJM1402,YJM1418,YJM1439,YJM1443,YJM1444,YJM1447,YJM1460,YJM1549,YJM1573,YJM1592,YJM244,YJM1083,YJM1129,YJM189,YJM193,YJM248,YJM271,YJM320,YJM326,YJM428,YJM450,YJM451,YJM453,YJM456,YJM541,YJM554,YJM555,YJM627,YJM681,YJM682,YJM969,YJM972,YJM975,YJM978,YJM981,YJM984,YJM987,YJM990,YJM996,YJM1133,YJM1190,YJM1199,YJM1202,YJM1208,YJM1242,YJM1244,YJM1250,YJM1307,YJM1311,YJM1326,YJM1332,YJM1336,YJM1338,YJM1341,YJM1355,YJM1356,YJM1381,YJM1383,YJM1386,YJM1400,YJM1401,YJM1415,YJM1417,YJM1419,YJM1433,YJM1450,YJM1463,YJM1477,YJM1478,YJM1479,YJM1526,YJM1527,YJM1574,YJM1615,Spar;o=Spar' \
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
bash plan_Scer_n7_pop.sh
bash 5_multi_cmd.sh

# other plans
bash plan_Scer_n7_Spar.sh
bash 5_multi_cmd.sh

# other plans
bash plan_Scer_n94_pop.sh
bash 5_multi_cmd.sh

# other plans
bash plan_Scer_n94_Spar.sh
bash 5_multi_cmd.sh

```

## Illumina (NCBI WGS_157)

```bash

mkdir -p ~/data/mrna-structure/alignment/scer_wgs
cd ~/data/mrna-structure/alignment/scer_wgs

perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
    -i ~/data/mrna-structure/GENOMES/WGS/scer_wgs_157.data.yml \
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
    --plan 'name=Scer_n157_pop;t=S288c;qs=wild007,wild006,wild005,wild004,wild003,wild002,wild001,wine019,wine018,wine017,wine016,wine015,wine014,wine013,wine012,wine011,wine010,wine009,wine008,wine007,wine006,wine005,wine004,wine003,wine002,wine001,spirits011,spirits010,spirits009,spirits008,spirits007,spirits006,spirits005,spirits004,spirits003,spirits002,spirits001,sake007,sake006,sake005,sake004,sake003,sake002,sake001,laboratory002,laboratory001,bread004,bread003,bread002,bread001,bioethanol005,bioethanol004,bioethanol003,bioethanol002,bioethanol001,beer102,beer101,beer100,beer099,beer098,beer097,beer096,beer095,beer094,beer093,beer092,beer091,beer090,beer089,beer088,beer087,beer086,beer085,beer084,beer083,beer082,beer081,beer080,beer079,beer078,beer077,beer076,beer075,beer074,beer073,beer072,beer071,beer070,beer069,beer068,beer067,beer066,beer065,beer064,beer063,beer062,beer061,beer060,beer059,beer058,beer057,beer056,beer055,beer054,beer053,beer052,beer051,beer050,beer049,beer048,beer047,beer046,beer045,beer044,beer043,beer042,beer041,beer040,beer039,beer038,beer037,beer036,beer035,beer034,beer033,beer032,beer031,beer030,beer029,beer028,beer027,beer026,beer025,beer024,beer023,beer022,beer021,beer020,beer019,beer018,beer017,beer016,beer015,beer014,beer013,beer012,beer011,beer010,beer009,beer008,beer007,beer006,beer005,beer004,beer003,beer002,beer001' \
    --plan 'name=Scer_n157_Spar;t=S288c;qs=wild007,wild006,wild005,wild004,wild003,wild002,wild001,wine019,wine018,wine017,wine016,wine015,wine014,wine013,wine012,wine011,wine010,wine009,wine008,wine007,wine006,wine005,wine004,wine003,wine002,wine001,spirits011,spirits010,spirits009,spirits008,spirits007,spirits006,spirits005,spirits004,spirits003,spirits002,spirits001,sake007,sake006,sake005,sake004,sake003,sake002,sake001,laboratory002,laboratory001,bread004,bread003,bread002,bread001,bioethanol005,bioethanol004,bioethanol003,bioethanol002,bioethanol001,beer102,beer101,beer100,beer099,beer098,beer097,beer096,beer095,beer094,beer093,beer092,beer091,beer090,beer089,beer088,beer087,beer086,beer085,beer084,beer083,beer082,beer081,beer080,beer079,beer078,beer077,beer076,beer075,beer074,beer073,beer072,beer071,beer070,beer069,beer068,beer067,beer066,beer065,beer064,beer063,beer062,beer061,beer060,beer059,beer058,beer057,beer056,beer055,beer054,beer053,beer052,beer051,beer050,beer049,beer048,beer047,beer046,beer045,beer044,beer043,beer042,beer041,beer040,beer039,beer038,beer037,beer036,beer035,beer034,beer033,beer032,beer031,beer030,beer029,beer028,beer027,beer026,beer025,beer024,beer023,beer022,beer021,beer020,beer019,beer018,beer017,beer016,beer015,beer014,beer013,beer012,beer011,beer010,beer009,beer008,beer007,beer006,beer005,beer004,beer003,beer002,beer001,Spar;o=Spar' \
    --plan 'name=Scer_n157_Beer1_Spar;t=S288c;qs=beer072,beer093,beer074,beer017,beer046,beer008,beer041,beer073,beer012,beer087,beer015,beer043,beer081,beer016,spirits005,beer022,beer026,beer076,beer077,beer075,beer031,beer078,beer060,beer035,beer051,beer065,beer071,beer067,beer007,beer102,beer098,beer082,beer044,beer068,beer058,beer042,beer052,beer097,beer095,beer096,beer019,beer057,beer066,beer094,beer045,beer050,beer090,beer099,beer001,beer069,beer055,beer056,beer100,beer064,beer053,beer070,beer047,wine012,beer009,beer089,beer079,beer010,beer036,beer048,beer037,beer101,beer049,beer054,Spar;o=Spar' \
    --plan 'name=Scer_n157_Beer2_Spar;t=S288c;qs=beer002,wine019,beer003,beer063,beer086,beer092,beer021,beer080,beer059,beer032,beer062,beer091,beer013,beer084,beer083,beer034,beer004,beer027,beer085,beer040,beer011,beer039,Spar;o=Spar' \
    --plan 'name=Scer_n157_Mixed_Spar;t=S288c;qs=spirits009,beer005,wine008,beer029,beer028,beer006,spirits001,beer023,bread004,beer038,spirits003,beer061,bread002,bread001,beer025,bread003,wild006,wild007,wild005,Spar;o=Spar' \
    --plan 'name=Scer_n157_Wine_Spar;t=S288c;qs=beer014,beer030,wine014,wine007,wine018,wine010,wine003,wine017,wine011,wine015,sake002,wine006,beer020,spirits011,wine004,wine005,beer024,beer088,wine009,wine001,beer033,wine013,spirits002,Spar;o=Spar' \
    --plan 'name=Scer_n157_Wild_Spar;t=S288c;qs=wine002,wild004,bioethanol004,bioethanol003,bioethanol001,sake006,sake001,sake004,sake007,sake005,sake003,wild002,wild003,wild001,Spar;o=Spar' \
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
bash plan_Scer_n157_pop.sh
bash 5_multi_cmd.sh

# other plans
bash plan_Scer_n157_Spar.sh
bash 5_multi_cmd.sh

# other plans
bash plan_Scer_n157_Beer1_Spar.sh
bash 5_multi_cmd.sh
bash plan_Scer_n157_Beer2_Spar.sh
bash 5_multi_cmd.sh
bash plan_Scer_n157_Mixed_Spar.sh
bash 5_multi_cmd.sh
bash plan_Scer_n157_Wine_Spar.sh
bash 5_multi_cmd.sh
bash plan_Scer_n157_Wild_Spar.sh
bash 5_multi_cmd.sh

```

# Build alignDB for multiple genomes n7

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n7_pop \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n7_pop_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n7_pop -r 1-60

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

# Build alignDB for multiple genomes n94

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n94_pop \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n94_pop_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n94_pop -r 1-60

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n94_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n94_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n94_Spar -r 1-60

```

# Build alignDB for multiple genomes n7p

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n7p_pop \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n7p_pop_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n7p_pop -r 1-60

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

# Build alignDB for multiple genomes n157

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_pop \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n157_pop_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_pop -r 1-60

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

# Build alignDB for multiple genomes n157_Beer1

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_Beer1_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n157_Beer1_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_Beer1_Spar -r 1-60

```

# Build alignDB for multiple genomes n157_Beer2

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_Beer2_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n157_Beer2_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_Beer2_Spar -r 1-60

```

# Build alignDB for multiple genomes n157_Mixed

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_Mixed_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n157_Mixed_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_Mixed_Spar -r 1-60

```

# Build alignDB for multiple genomes n157_Wine

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_Wine_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n157_Wine_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_Wine_Spar -r 1-60

```

# Build alignDB for multiple genomes n157_Wild

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_Wild_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n157_Wild_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_Wild_Spar -r 1-60

```

##  Extract `gene_list` and `snp_codon_list` n7_Spar

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n7_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n7_Spar.mvar.gene_list.csv
```


##  Extract `gene_list` and `snp_codon_list` n7_pop

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_pop.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n7_pop.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_pop.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n7_pop.mvar.gene_list.csv
```


##  Extract `gene_list` and `snp_codon_list` n94_Spar

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n94_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n94_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n94_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n94_Spar.mvar.gene_list.csv
```

##  Extract `gene_list` and `snp_codon_list` n94_pop

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n94_pop.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n94_pop.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n94_pop.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n94_pop.mvar.gene_list.csv
```

##  Extract `gene_list` and `snp_codon_list` n7p_Spar

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n7p_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n7p_Spar.mvar.gene_list.csv
```


##  Extract `gene_list` and `snp_codon_list` n7p_pop

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_pop.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n7p_pop.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_pop.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n7p_pop.mvar.gene_list.csv
```

##  Extract `gene_list` and `snp_codon_list` n157_Spar

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_Spar.mvar.gene_list.csv
```


##  Extract `gene_list` and `snp_codon_list` n157_pop

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_pop.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_pop.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_pop.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_pop.mvar.gene_list.csv
```

##  Extract `gene_list` and `snp_codon_list` n157_Beer1_Spar

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Beer1_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_Beer1_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Beer1_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_Beer1_Spar.mvar.gene_list.csv
```

##  Extract `gene_list` and `snp_codon_list` n157_Beer2_Spar

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Beer2_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_Beer2_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Beer2_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_Beer2_Spar.mvar.gene_list.csv
```

##  Extract `gene_list` and `snp_codon_list` n157_Mixed_Spar

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Mixed_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_Mixed_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Mixed_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_Mixed_Spar.mvar.gene_list.csv
```

##  Extract `gene_list` and `snp_codon_list` n157_Wine_Spar

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Wine_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_Wine_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Wine_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_Wine_Spar.mvar.gene_list.csv
```

##  Extract `gene_list` and `snp_codon_list` n157_Wild_Spar

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Wild_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_Wild_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Wild_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_Wild_Spar.mvar.gene_list.csv
```

## SNPs and indels n7_Spar

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


## SNPs and indels n7_pop

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_pop.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n7_pop.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_pop.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n7_pop.indel.pos.txt

```


## SNPs and indels n94_Spar

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n94_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n94_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n94_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n94_Spar.indel.pos.txt

```


## SNPs and indels n94_pop

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n94_pop.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n94_pop.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n94_pop.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n94_pop.indel.pos.txt

```

## SNPs and indels n7p_Spar

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


## SNPs and indels n7p_pop

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_pop.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n7p_pop.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_pop.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n7p_pop.indel.pos.txt

```

## SNPs and indels n157_Spar

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


## SNPs and indels n157_pop

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_pop.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_pop.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_pop.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_pop.indel.pos.txt

```

## SNPs and indels n157_Beer1_Spar

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Beer1_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_Beer1_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Beer1_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_Beer1_Spar.indel.pos.txt

```

## SNPs and indels n157_Beer2_Spar

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Beer2_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_Beer2_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Beer2_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_Beer2_Spar.indel.pos.txt

```

## SNPs and indels n157_Mixed_Spar

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Mixed_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_Mixed_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Mixed_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_Mixed_Spar.indel.pos.txt

```

## SNPs and indels n157_Wine_Spar

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Wine_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_Wine_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Wine_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_Wine_Spar.indel.pos.txt

```

## SNPs and indels n157_Wild_Spar

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Wild_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_Wild_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Wild_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_Wild_Spar.indel.pos.txt

```

List all valid genes.

```mysql
SELECT *
FROM gene g, window w
where g.window_id  = w.window_id
and gene_is_full = 1
and gene_biotype = 'protein_coding'
and w.window_differences > 2
and g.gene_description not like "Putative%"
and g.gene_description not like "Dubious%"
and g.gene_description not like "Identified%"
order by w.window_length
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
# utr (5' and 3')
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

# utr 5' and 3'
jrunlist compare --op diff sce_genes.yml sce_orf_genomic.yml -o sce_utr.yml
runlist convert sce_utr.yml -o sce_utr.pos.txt

# 5' of orf_genomic
cat ../sgd/orf_genomic_all.fasta \
    | perl -n -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+(\d+)\-(\d+)/ or next;
        $1 eq "Mito" and next;

        if ($2 < $3) {
            print qq{$1:$2\n};
        }
        else {
            print qq{$1:$3\n};
        }
    ' \
    > sce_5prime.pos.txt
jrunlist cover sce_5prime.pos.txt -o sce_5prime.yml
jrunlist span --op pad -n 5 sce_5prime.yml -o sce_5prime_pad5.yml

# split 5' and 3' utr
runlist position --op overlap \
    sce_5prime_pad5.yml sce_utr.pos.txt \
    -o sce_utr5.pos.txt
jrunlist cover sce_utr5.pos.txt -o sce_utr5.yml

runlist position --op non-overlap \
    sce_5prime_pad5.yml sce_utr.pos.txt \
    -o sce_utr3.pos.txt
jrunlist cover sce_utr3.pos.txt -o sce_utr3.yml

# Stats
printf "| %s | %s | %s | %s |\n" \
    "Name" "chrLength" "size" "coverage" \
    > coverage.stat.md
printf "|:--|--:|--:|--:|\n" >> coverage.stat.md

for f in genes intergenic intron orf_genomic utr utr5 utr3; do
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
| utr5        |  12071326 |  275817 |   0.0228 |
| utr3        |  12071326 |  240752 |   0.0199 |

# Real Processing n7_Spar

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

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n7_pop

```bash
export NAME=Scer_n7_pop

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

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n94_Spar

```bash
export NAME=Scer_n94_Spar

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

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

# SNPs within orf_genomic
runlist position --op superset \
    sce_orf_genomic.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

unset NAME
```

# Real Processing n94_pop

```bash
export NAME=Scer_n94_pop

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

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n7p_Spar

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

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n7p_pop

```bash
export NAME=Scer_n7p_pop

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

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n157_Spar

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

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n157_pop

```bash
export NAME=Scer_n157_pop

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

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n157_Beer1_Spar

```bash
export NAME=Scer_n157_Beer1_Spar

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

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n157_Beer2_Spar

```bash
export NAME=Scer_n157_Beer2_Spar

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

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n157_Mixed_Spar

```bash
export NAME=Scer_n157_Mixed_Spar

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

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n157_Wine_Spar

```bash
export NAME=Scer_n157_Wine_Spar

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

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Real Processing n157_Wild_Spar

```bash
export NAME=Scer_n157_Wild_Spar

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

# SNPs within 5' and 3' utr
runlist position --op superset \
    sce_utr5.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr5.pos.txt

runlist position --op superset \
    sce_utr3.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr3.pos.txt

unset NAME
```

# Stats

Switch to RStudio, let R do its jobs.

```bash
open -a RStudio ~/data/mrna-structure

```

# Pack all things up

```bash
cd  ~/data
tar -cf - mrna-structure/ | xz -9 -c - > mrna-structure.tar.xz

```

