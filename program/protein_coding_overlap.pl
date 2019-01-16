#!/usr/bin/perl
use strict;
#use warnings;
use autodie;

use Getopt::Long qw();
use FindBin;
use YAML::Syck qw();

use Path::Tiny;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use List::Util qw/min/;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/protein_coding_overlap.pl --file ~/data/mrna-structure/phylogeny/protein_coding_list_range_chr_strand.csv --output ~/data/mrna-structure/phylogeny/protein_coding_overlap.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Process protein_coding info...");

#----------------------------------------------------------#
# start process
#----------------------------------------------------------#

open OUT, ">",$output;
print OUT "target_gene" . ",". "target_chr\(strand\)" . ",". "target_range" . "," . "target_length" . "," .  "query_gene" . "," ."query_chr\(strand\)" . "," . "query_range" . "," . "query_length" . "," ."overlap" . "," . "overlap_length" . "\n";

my $csv_fh;
open $csv_fh,'<',$file;

my %gene;

while (<$csv_fh>) {
    chomp;
    my @gene_info = split /,/, $_;
    my $gene_range = AlignDB::IntSpan->new;
    $gene_range->add("$gene_info[3]-$gene_info[4]");
    my @info = ($gene_info[1],$gene_info[2],$gene_range);
    $gene{"$gene_info[0]"} = \@info;
    #print YAML::Syck::Dump(\%gene),"\n";sleep 1;
}

foreach my $target_gene ( sort keys %gene) {
    print $target_gene,"\n";
    foreach my $query_gene ( sort keys %gene) {
    	  print "pairwise $target_gene - $query_gene"."\n";
    	  if ($target_gene ne $query_gene){
    	  	  
    	  	  
            if ($gene{"$target_gene"}->[0] eq $gene{"$query_gene"}->[0]){
            	  
            	  my $target_length = $gene{"$target_gene"}->[2]->cardinality();
            	  my $query_length = $gene{"$query_gene"}->[2]->cardinality();
            	  
            		my $intersection = AlignDB::IntSpan::intersect( $gene{"$target_gene"}->[2], $gene{"$query_gene"}->[2] );	
            	  my $intersection_length = $intersection->cardinality();
            	  
            	  #print YAML::Syck::Dump($target_length),"\n";sleep 1;
            	  #print YAML::Syck::Dump($query_length),"\n";sleep 1;
            	  #print YAML::Syck::Dump($intersection),"\n";sleep 1;
            	  #print YAML::Syck::Dump($intersection_length),"\n";sleep 1;
            	  
            	  if ( $intersection_length != 0 ){
            	      
            	      my @inf;
            	      $inf[0] = $target_gene;
            	      my $chr1 = $gene{"$target_gene"}->[0];
            	      my $strand1 = $gene{"$target_gene"}->[1];
            	      $inf[1] = "$chr1".'('."$strand1".')';
            	      $inf[2] = $gene{"$target_gene"}->[2];
            	      $inf[3] = $target_length;
            	      $inf[4] = $query_gene;
            	      my $chr2 = $gene{"$query_gene"}->[0];
            	      my $strand2 = $gene{"$query_gene"}->[1];
            	      $inf[5] = "$chr2".'('."$strand2".')';
            	      $inf[6] = $gene{"$query_gene"}->[2];
            	      $inf[7] = $query_length;
            	      $inf[8] = $intersection;
            	      $inf[9] = $intersection_length;
            	      my $inf = join ",", @inf;
            	      
            	      open OUT, ">>",$output;
            	      print OUT $inf , "\n";
            	      
                }
            }
    	  }
    }
}

close OUT;

__END__
