#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use File::Find::Rule;
use AlignDB::IntSpan;
use YAML::Syck qw();

use Getopt::Long qw();
use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/count_gene_range.pl --file ~/data/mrna-structure/phylogeny/protein_coding_list_range.csv --dir ~/data/mrna-structure/phylogeny/Scer_n128_Spar_gene_alignment_cds --output ~/data/mrna-structure/phylogeny/gene_range.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'dir|d=s'   => \my $dir,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

my @files = File::Find::Rule->file()->name('*.fas.fas')
  ->in($dir);
  
my $csv_fh;
open $csv_fh, '<', $file;
open OUT, ">",$output or die $!;
open OUT1, '>>', $output;
print OUT1 "gene" . ","
  . "alignment_length" . ","
  . "reference_length" . ","
  . "intersection_length" . ","
  . "proporation" . "\n";

my %hash1;
my %hash2;
foreach my $fas (@files) {
    my @data = split "/", $fas;
    my @data1 = split ".fas.fas", $data[-1], 2;

    #print $data1[0]."\n";

    open my $file_fh, $fas or die $!;

    my $range;
    my $range2 = AlignDB::IntSpan->new();
    while (<$file_fh>) {
        $_ =~ /^>S288c.*:(.*)\n/;
        $range = $1;
        $range2->add($range);

        #print $range."\n";
    }

    $hash1{"$data1[0]"} = $range2;

    #print $range2."\n";sleep 1;

}

while (<$csv_fh>) {
    chomp;

    #print $_, "\n";
    my @content = split /,/, $_;
    my $range2 = AlignDB::IntSpan->new;
    $range2->add_pair( $content[1], $content[2] );

    #print $range2."\n";
    $hash2{"$content[0]"} = $range2;
}

#print YAML::Syck::Dump(\%hash1),"\n";
#print $hash1{"YAL007C"},"\n";
#print keys %hash1,"\n";

foreach $_ ( sort keys %hash2 ) {
    my $gene = $_;
    print $gene,"\n";
    my @test = keys %hash1;
    if ( grep {$_ eq $gene} @test ) {
        my $value1 = $hash1{"$gene"};
        my $value2 = $hash2{"$gene"};
        my $intersection =
          AlignDB::IntSpan::intersect( $hash2{"$gene"}, $hash1{"$gene"} );
        #print $intersection, "\n";
        my $range_cds_length    = $value1->cardinality();
        my $range2_cds_length    = $value2->cardinality();
        my $intersection_length = $intersection->cardinality();
        my $proporation;

        if ( $range2_cds_length != 0 ) {
            $proporation = $intersection_length / $range2_cds_length;
        }
        #print $proporation,"\n";
        open OUT1,'>>',$output;
        print OUT1 "$gene" . ","
          . "$range_cds_length" . ","
          . "$range2_cds_length" . ","
          . "$intersection_length" . ","
          . "$proporation" . "\n";
    }
}

close OUT;
close OUT1;

__END__
