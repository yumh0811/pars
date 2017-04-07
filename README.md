# Processing Yeast PARS Data

[TOC]: # " "

- [Download reference data](#download-reference-data)
- [Other strains and outgroups](#other-strains-and-outgroups)
    - [Sanger (WGS)](#sanger-wgs)
- [Build alignDB for multiple genomes](#build-aligndb-for-multiple-genomes)
- [Build self alignDB for gene information](#build-self-aligndb-for-gene-information)
- [Blast](#blast)
- [SNPs and indels](#snps-and-indels)
- [Real Processing](#real-processing)
- [Pack all things up](#pack-all-things-up)
- [Stats](#stats)


## Download reference data

Download PARS10 full site.

```bash
mkdir -p ~/data/mrna-structure/PARS10
cd ~/data/mrna-structure/PARS10

perl ~/Scripts/download/list.pl -u http://genie.weizmann.ac.il/pubs/PARS10/
perl ~/Scripts/download/download.pl -i pubs_PARS10.yml

find . -name "*.gz" | xargs gzip -d
```

Download S288c annotation data from ensembl by
[rsync](http://www.ensembl.org/info/data/ftp/rsync.html?redirect=no).

S288c assembly version is not changed from 2011, R64-1-1 (GCA_000146045.2).

```bash
mkdir -p ~/data/mrna-structure/ensembl82/mysql
cd ~/data/mrna-structure/ensembl82/mysql
rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-82/mysql/saccharomyces_cerevisiae_core_82_4 .

mkdir -p ~/data/mrna-structure/ensembl82/fasta
cd ~/data/mrna-structure/ensembl82/fasta
rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-82/fasta/saccharomyces_cerevisiae .

perl ~/Scripts/withncbi/ensembl/build_ensembl.pl -e ~/data/mrna-structure/ensembl82/mysql/saccharomyces_cerevisiae_core_82_4 --checksum
perl ~/Scripts/withncbi/ensembl/build_ensembl.pl -e ~/data/mrna-structure/ensembl82/mysql/saccharomyces_cerevisiae_core_82_4 --initdb --db yeast_82
```

SGD.

```bash
mkdir -p ~/data/mrna-structure/sgd
cd ~/data/mrna-structure/sgd
aria2c -c http://downloads.yeastgenome.org/sequence/S288C_reference/intergenic/NotFeature.fasta.gz
aria2c -c http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
aria2c -c http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_genomic_all.fasta.gz

find . -name "*.gz" | xargs gzip -d
```

## Other strains and outgroups

* `withncbi/db/`: taxonomy database
* `withncbi/pop/`: scer_wgs alignments

### Sanger (NCBI WGS)

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

### Illumina (NCBI ASSEMBLY)

```bash
mkdir -p ~/data/mrna-structure/GENOMES/ASSEMBLIES
cd ~/data/mrna-structure/GENOMES/ASSEMBLIES

# Download, rename files and change fasta headers
perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
    -p -f ~/Scripts/pars/scer_100.seq.csv

```

## Plans of alignments

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

## Build alignDB for multiple genomes

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Scer_n8_Spar \
    -da ~/data/alignment/Fungi/scer_wgs/Scer_n8_Spar_refined \
    --ensembl yeast_82 \
    --block \
    --outgroup \
    --chr ~/data/alignment/Fungi/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run all
```

Extract `gene_list` from `Scer_n8_Spar.mvar.xlsx`.

Extract `snp_codon_list` from `Scer_n8_Spar.mvar.xlsx`.

```bat
cd /d d:\data\mrna-structure\xlsx\

perl d:\Scripts\fig_table\collect_excel.pl -f Scer_n8_Spar.mvar.xlsx -s gene_list -n gene_list -o Scer_n8_Spar.mvar.gene_list.xlsx
perl d:\Scripts\fig_table\xlsx2xls.pl -d Scer_n8_Spar.mvar.gene_list.xlsx --csv
rm Scer_n8_Spar.mvar.gene_list.xlsx

perl d:\Scripts\fig_table\collect_excel.pl -f Scer_n8_Spar.mvar.xlsx -s snp_codon_list -n snp_codon_list -o Scer_n8_Spar.mvar.snp_codon_list.xlsx
perl d:\Scripts\fig_table\xlsx2xls.pl -d Scer_n8_Spar.mvar.snp_codon_list.xlsx --csv
rm Scer_n8_Spar.mvar.snp_codon_list.xlsx
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

## Build self alignDB for gene information

```bash
perl ~/Scripts/alignDB/init/init_alignDB.pl \
    -d S288Cvsself_gene \
    -taxon ~/data/alignment/Fungi/scer_wgs/taxon.csv \
    -chr ~/data/alignment/Fungi/scer_wgs/chr_length.csv

perl ~/Scripts/alignDB/init/gen_alignDB_genome.pl \
    -d S288Cvsself_gene \
    -t "559292,S288c" \
    -dir ~/data/alignment/Fungi/scer_wgs/Genomes/S288c \
    --length 1_000_000 \
    --parallel 8

perl ~/Scripts/alignDB/gene/insert_gene.pl -d S288Cvsself_gene -e yeast_82 --batch 1 --parallel 8

perl ~/Scripts/alignDB/init/update_sw_cv.pl -d S288Cvsself_gene --batch 1 --parallel 8
perl ~/Scripts/alignDB/init/update_feature.pl -d S288Cvsself_gene -e yeast_82 --batch 1 --parallel 8

perl ~/Scripts/alignDB/gene/update_gene_yeast_ess.pl -d S288Cvsself_gene
perl ~/Scripts/alignDB/gene/update_gene_yeast_quan.pl -d S288Cvsself_gene

cat <<EOF > query.tmp
SELECT
    g.gene_stable_id gene,
    g.gene_feature10 quan,
    g.gene_feature4 ess,
    g.gene_feature7 interact,
    g.gene_feature6 rec,
    AVG(sw.codingsw_cv) avg_cv,
    AVG(sw.codingsw_intra_cv) avg_intra_cv
FROM
    gene g,
    exon e,
    codingsw sw,
    window w
WHERE
    1 = 1
        AND g.gene_id = e.gene_id
        AND e.exon_id = sw.exon_id
        AND sw.window_id = w.window_id
        AND sw.codingsw_distance < 0
        AND g.gene_feature10 IS NOT NULL
GROUP BY g.gene_id
EOF

perl ~/Scripts/alignDB/util/query_sql.pl -d S288Cvsself_gene -f query.tmp -o ~/data/mrna-structure/xlsx/S288Cvsself_gene.csv
rm query.tmp

```

## Blast

Prepare a combined fasta file of yeast genome.

```bash
cd ~/data/alignment/Fungi/scer_wgs/Genomes/S288c

cat \
I.fa     \
II.fa    \
III.fa   \
IV.fa    \
V.fa     \
VI.fa    \
VII.fa   \
VIII.fa  \
IX.fa    \
X.fa     \
XI.fa    \
XII.fa   \
XIII.fa  \
XIV.fa   \
XV.fa    \
XVI.fa   \
> S288c.fa
```

Move `S288c.fa` to the blast folder and blast genes against the genome.

```bash
mkdir -p ~/data/mrna-structure/blast
cd ~/data/mrna-structure/blast

mv ~/data/alignment/Fungi/scer_wgs/Genomes/S288c/S288c.fa .

perl -nl -i -e '/^>/ or $_ = uc $_; print'  S288c.fa

# formatdb
~/share/blast/bin/formatdb -p F -o T -i S288c.fa

# blast every transcripts against genome
~/share/blast/bin/blastall -p blastn -F "m D" -m 0 -b 10 -v 10 -e 1e-3 -a 4 -i ~/data/mrna-structure/PARS10/pubs/PARS10/data/sce_genes.fasta -d S288C.fa -o sce_genes.blast
```

## SNPs and indels

Select columns `chr_name	snp_pos snp_pos` and manually create snp bed file
`~/data/mrna-structure/process/Scer_n8_Spar.snp.bed` from
`~/data/mrna-structure/xlsx/Scer_n8_Spar.mvar.xlsx`

Select columns `chr_name	indel_start indel_end` and manually create indel bed file
`~/data/mrna-structure/process/Scer_n8_Spar.indel.bed`.

## Real Processing

```bash
mkdir -p ~/data/mrna-structure/process
export NAME=Scer_n8_Spar

#----------------------------------------------------------#
# gene
#----------------------------------------------------------#
cd ~/data/mrna-structure/process

# parse blastn output
perl ~/Scripts/pars/blastn_transcript.pl -f ~/data/mrna-structure/blast/sce_genes.blast -m 0

# produce transcript set
# YLR167W	568	chrXII	498888	499455	+
perl -an -e 'print qq{$F[2]\t$F[3]\t$F[4]\n}' \
    ~/data/mrna-structure/process/sce_genes.blast.tsv \
    | sort \
    > ~/data/mrna-structure/process/sce_genes.bed

# snps within transcripts
perl ~/Scripts/alignDB/ofg/bed_op.pl --op bed_intersect --file ~/data/mrna-structure/process/$NAME.snp.bed --name ~/data/mrna-structure/process/sce_genes.bed
mv $NAME.snp.bed.bed_intersect $NAME.snp.gene.bed

# read gene and snp info file
# produce $NAME.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl --pars ~/data/mrna-structure/PARS10/pubs/PARS10/data --gene ~/data/mrna-structure/process/sce_genes.blast.tsv --pos $NAME.snp.gene.bed  > fail_pos.txt

# check fail_pos.txt to find snps located in overlap genes
# produce fail_pos.xlsx

# process $NAME.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.yml

#----------------------------------------------------------#
# intergenic
#----------------------------------------------------------#
cd ~/data/mrna-structure/process

# produce intergenic set
perl -n -e '/>/ or next; /Chr\s+(\w+)\s+from\s+(\d+)\-(\d+)/ and print qq{$1\t$2\t$3\n}' ~/data/mrna-structure/sgd/NotFeature.fasta > ~/data/mrna-structure/process/intergenic.bed

# snps within intergenic
perl ~/Scripts/alignDB/ofg/bed_op.pl --op bed_intersect --file ~/data/mrna-structure/process/$NAME.snp.bed --name ~/data/mrna-structure/process/intergenic.bed
mv $NAME.snp.bed.bed_intersect $NAME.snp.intergenic.bed

# convert to unique snp name
perl -an -e '
    BEGIN{print qq{name\n}};
    print qq{$F[0]:$F[1]\n}
    ' \
    ~/data/mrna-structure/process/$NAME.snp.intergenic.bed \
    > ~/data/mrna-structure/process/$NAME.intergenic.snp.tsv

#----------------------------------------------------------#
# utr (5' and 3')
#----------------------------------------------------------#
# FIXME: seperate 5' and 3' utrs
cd ~/data/mrna-structure/process

# produce orf_genomic set
perl -n -e '/>/ or next; /Chr\s+(\w+)\s+from\s+(\d+)\-(\d+)/ and print qq{$1\t$2\t$3\n}' ~/data/mrna-structure/sgd/orf_genomic_all.fasta > ~/data/mrna-structure/process/orf_genomic.bed

# remove in orf_genomic set from snp.gene, so the left is snps.utr
perl ~/Scripts/alignDB/ofg/bed_op.pl --op bed_diff --file ~/data/mrna-structure/process/$NAME.snp.gene.bed --name  ~/data/mrna-structure/process/orf_genomic.bed
mv $NAME.snp.gene.bed.bed_diff $NAME.snp.utr.bed

# convert to unique snp name
perl -an -e '
    BEGIN{print qq{name\n}};
    print qq{$F[0]:$F[1]\n};
    ' \
    ~/data/mrna-structure/process/$NAME.snp.utr.bed \
    > ~/data/mrna-structure/process/$NAME.utr.snp.tsv

unset NAME
```

## Pack all things up

```bash
cd  ~/data
tar -cf - mrna-structure/ | xz -9 -c - > mrna-structure.tar.xz

```

## Stats

Switch to RStudio, let R do its jobs.

```bash
open -a RStudio ~/data/mrna-structure

```

