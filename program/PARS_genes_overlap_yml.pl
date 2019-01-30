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

    perl ~/Scripts/pars/program/PARS_genes_overlap_yml.pl --file ~/data/mrna-structure/phylogeny/PARS_genes_overlap.csv --output ~/data/mrna-structure/phylogeny/overlap.yml

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
$stopwatch->start_message("Process PARS info...");

#----------------------------------------------------------#
# start process
#----------------------------------------------------------#

my $csv_fh;
open $csv_fh,'<',$file;

my %gene_info;
my %gene_info_new;

while (<$csv_fh>) {
    chomp;
    my @gene_info = split /,/, $_ ;
    
    if($gene_info[0] ne "target_gene"){
    	  my $chr = $gene_info[1];
    	  $chr =~ /([a-zA-Z]+)\(/;
    	  $chr = $1;

        push @{$gene_info{"$chr"}} , $gene_info[8];
        
        #print YAML::Syck::Dump(\%gene_info),"\n";sleep 1;   
    }
}
foreach my $chr ( sort keys %gene_info ) {
    my @join = @{$gene_info{"$chr"}};
		my $join = join (",",@join);
		my $gene_range = AlignDB::IntSpan->new("$join");
		$gene_info_new{"$chr"} = $gene_range -> as_string;
}

#print YAML::Syck::Dump(\%gene_info),"\n";sleep 5;

open OUT, ">",$output;
print OUT YAML::Syck::Dump(\%gene_info_new),"\n";

close OUT;

__END__
