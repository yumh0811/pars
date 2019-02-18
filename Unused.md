
## Download S288c annotation data from ensembl by rsync

http://www.ensembl.org/info/data/ftp/rsync.html?redirect=no

S288c assembly version is not changed since 2011, R64-1-1 (GCA_000146045.2).

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


## mRNA levels

* Data source: Quantification of the yeast transcriptome by single-molecule sequencing. Lipson, D.
  et al. Nature Biotechnology 27, 652-658 (2009) doi:10.1038/nbt.1551

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://www.nature.com/nbt/journal/v27/n7/extref/nbt.1551-S2.xls

perl ~/Scripts/fig_table/xlsx2csv.pl -f nbt.1551-S2.xls --sheet 'counts' |
    perl -nla -F/,/ -e '
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

* ess: Giaever, G., et al. Functional Profiling of theSaccharomyces cerevisiae Genome. Nature 418,
  387-391. (2002)
* rich/minimal: Mechanisms of Haploinsufficiency Revealed by Genome-Wide Profiling in Yeast
  Deutschbauer, AM. et al. GENETICS April 1, 2005 vol. 169 no. 4 1915-1925;
  10.1534/genetics.104.036871
* chem: The Chemical Genomic Portrait of Yeast: Uncovering a Phenotype for All Genes. Hillenmeyer,
  M.E. et al. Science 18 Apr 2008: Vol. 320, Issue 5874, pp. 362-365 DOI: 10.1126/science.1150021

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt

cat Essential_ORFs.txt |
    perl -nl -e '
        next unless /^\d/;
        my $orf = ( split /\s+/ )[1];
        print qq{$orf};
    ' \
    > ess_orf.tsv

wget -N http://www-sequence.stanford.edu/group/research/HIP_HOP/supplements/01yfh/files/OrfGeneData.txt

cat OrfGeneData.txt |
    perl -nla -F"\t" -e '
        printf q{#} if /^orf/;
        print qq{$F[0]\t$F[1]\t$F[5]\t$F[13]\t$F[17]};
    ' \
    > rich_orf.tsv

# http://chemogenomics.stanford.edu/supplements/global/

wget -N http://chemogenomics.stanford.edu/supplements/global/download/data/hom.z_tdist_pval_nm.counts.smallmol.cutoff.01.xls

perl ~/Scripts/fig_table/xlsx2csv.pl \
    -f hom.z_tdist_pval_nm.counts.smallmol.cutoff.01.xls \
    --sheet 'hom.z_tdist_pval_nm.smallmol.co' |
    perl -nla -F"," -MText::CSV_XS -e '
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

* Data source: Global mapping of meiotic recombination hotspots and coldspots in the yeast
  Saccharomyces cerevisiae. vol. 97 no. 21 PNAS Jennifer L. Gerton, 11383–11390

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://derisilab.ucsf.edu/data/hotspots/forWebORFs.txt

cat forWebORFs.txt |
    perl -nla -F"\t" -MStatistics::Lite -e '
        next unless /^Y/;    # ORF stable id start with a "Y"
        next if @F < 2;
        my $rec_rate  = Statistics::Lite::median(grep {defined} @F[1 .. 7]);
        print qq{$F[0]\t$rec_rate};
    ' \
    > rec_rate.tsv
```


# HKA test

```bash
cd ~/Scripts/pars/
cat HKA_prepare_chisq.txt |
    grep -v '^#' |
    parallel --line-buffer -j 16 '
        echo >&2 {}
        HKA_NAME=$(echo {} | cut -f 1)
        HKA_SIZE=$(echo {} | cut -f 2)
        HKA_NUMBERS=$(echo {} | cut -f3-10)
        # echo $HKA_NUMBERS

        python hka_test.py --name ${HKA_NAME} --size ${HKA_SIZE} $HKA_NUMBERS
    ' \
    > HKA_chisq.txt

```

## Protein-protein interactions

* Data source:

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://drygin.ccbr.utoronto.ca/%7Ecostanzo2009/sgadata_costanzo2009_stringentCutoff_101120.txt.gz

gzip -d -c sgadata_costanzo2009_stringentCutoff_101120.txt.gz |
    perl -nla -F"\t" -e '
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
