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

### Sanger (WGS)

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

```bash
mkdir -p ~/data/mrna-structure/GENOMES/ASSEMBLIES
cd ~/data/mrna-structure/GENOMES/ASSEMBLIES

# Download S288c, EC1118 separately
perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
    -f ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_assembly_report.txt \
    --nuclear -name S288c \
    > S288c.seq.csv

perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
    -f ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/218/975/GCA_000218975.1_ASM21897v1/GCA_000218975.1_ASM21897v1_assembly_report.txt \
    --nuclear --genbank --scaffold -name EC1118 \
    > EC1118.seq.csv

echo "#strain_name,accession,strain_taxon_id,seq_name" > scer_wgs.seq.csv
cat S288c.seq.csv EC1118.seq.csv \
    | perl -nl -e '/^#/ and next; /^\s*$/ and next; print;' \
    >> scer_wgs.seq.csv

# Download, rename files and change fasta headers
perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
    -p -f scer_wgs.seq.csv

```

```bash
# create downloaded genome list
cat ~/data/mrna-structure/GENOMES/ASSEMBLIES/scer_wgs.seq.csv \
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
    --plan 'name=Scer_n7_pop;t=S288c;qs=EC1118,Kyokai_no_7,RM11_1a,Sigma1278b,T7,YJM789' \
    --plan 'name=Scer_n7_Spar;t=S288c;qs=EC1118,Kyokai_no_7,RM11_1a,Sigma1278b,T7,YJM789,Spar;o=Spar' \
    -y

# pop_prep.pl
perl ~/Scripts/withncbi/pop/pop_prep.pl -p 8 -i scer_wgs.plan.yml

sh 01_file.sh
sh 02_rm.sh
sh 03_strain_info.sh

# plan_ALL.sh
sh plan_ALL.sh

sh 1_real_chr.sh
sh 3_pair_cmd.sh
sh 4_rawphylo.sh
sh 5_multi_cmd.sh

# other plans
sh plan_Scer_n7_pop.sh
sh 5_multi_cmd.sh
sh 7_multi_db_only.sh

# other plans
sh plan_Scer_n7_Spar.sh
sh 5_multi_cmd.sh
sh 7_multi_db_only.sh

find . -name "*.gz" | xargs gunzip
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

```sql
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

Select columns `chr_name	snp_pos  snp_pos` and manually create snp bed file `~/data/mrna-structure/process/Scer_n8_Spar.snp.bed` from `~/data/mrna-structure/xlsx/Scer_n8_Spar.mvar.xlsx`

Select columns `chr_name	indel_start    indel_end` and manually create indel bed file `~/data/mrna-structure/process/Scer_n8_Spar.indel.bed`.

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
