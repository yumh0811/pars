#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use FindBin;
use YAML::Syck qw();

use Path::Tiny;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/cut_cds_yml.pl --file ~/data/mrna-structure/phylogeny/protein_coding_list_range_chr.csv --output gene_cds_all_yml

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

open my $in_fh,"<",$file;

while (<$in_fh>) {
     chomp;
     my @data = split ",", $_;
     #print $_,"\n";
     #print $data[0],"\n";
     open OUT,"> $output/$data[0].yml"|| die "$!";
     print OUT "--- ","\n",$data[1],": ",$data[2],"-",$data[3];
        }
close OUT;