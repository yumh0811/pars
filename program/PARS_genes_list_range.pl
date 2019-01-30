#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use FindBin;
use YAML::Syck qw();

use Path::Tiny;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use App::Rangeops;

#use List::Util;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/PARS_genes_list_range_chr.pl --file ~/data/mrna-structure/process/sce_genes.blast.tsv --output ~/data/mrna-structure/phylogeny/PARS_genes_list_range_chr.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);



my $tsv_fh;
open $tsv_fh, '<', $file
  or die "Could not open sce_genes.blast.tsv:$! ";

open OUT, '>', $output
  || die "$!";

while (<$tsv_fh>) {
    chomp;
    my @content = split /	/, $_;
    my @content_new = ();
    $content_new[0] = $content[0];
    $content_new[1] = $content[3];
    $content_new[2] = $content[4];
    my $content_new = join ",", @content_new;
    open OUT, '>>', $output;
    print OUT $content_new, "\n";
}

close $tsv_fh;
close OUT;

__END__
