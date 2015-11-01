# Processing Yeast PARS Data

## Preparing steps

* `withncbi/db`: taxonomy database
* `withncbi/pop/`: scer_wgs alignments

## Download data

Download PARS10 full site.

```bash
mkdir -p ~/data/mrna-structure/PARS10
cd ~/data/mrna-structure/PARS10

perl ~/Scripts/download/list.pl -u http://genie.weizmann.ac.il/pubs/PARS10/
perl ~/Scripts/download/download.pl -i pubs_PARS10.yml

find . -name "*.gz" | xargs gunzip
```

Download S288c annotation data from ensembl by [rsync](http://www.ensembl.org/info/data/ftp/rsync.html?redirect=no).

S288c assembly version is not changed from 2011, R64-1-1 (GCA_000146045.2).

```bash
mkdir -p ~/data/ensembl82/mysql
cd ~/data/ensembl82/mysql
rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-82/mysql/saccharomyces_cerevisiae_core_82_4 .

rsync -avP rsync://ftp.ensembl.org/ensembl/pub/grch37/release-82/mysql/homo_sapiens_core_82_37 .
rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-82/mysql/mus_musculus_core_82_38 .

mkdir -p ~/data/ensembl82/fasta
cd ~/data/ensembl82/fasta
rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-82/fasta/saccharomyces_cerevisiae .

perl ~/Scripts/alignDB/util/build_ensembl.pl -e ~/data/ensembl82/mysql/saccharomyces_cerevisiae_core_82_4 --checksum
perl ~/Scripts/alignDB/util/build_ensembl.pl -e ~/data/ensembl82/mysql/saccharomyces_cerevisiae_core_82_4 --initdb --db yeast_82
```

SGD.

```bash
mkdir -p ~/data/mrna-structure/sgd
cd ~/data/mrna-structure/sgd
wget -N http://downloads.yeastgenome.org/sequence/S288C_reference/intergenic/NotFeature.fasta.gz
wget -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
wget -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_genomic_all.fasta.gz

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
    --id ~/data/alignment/Fungi/scer_wgs/id2name.csv \
    --outgroup \
    -taxon ~/data/alignment/Fungi/scer_wgs/taxon.csv \
    -chr ~/data/alignment/Fungi/scer_wgs/chr_length.csv \
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
