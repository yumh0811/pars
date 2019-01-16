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

    perl ~/Scripts/pars/program/protein_coding_list_range_chr.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output ~/data/mrna-structure/phylogeny/protein_coding_list_range_chr.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);



my $tsv_fh;
open $tsv_fh, '<', $file
  or die "Could not open saccharomyces_cerevisiae.gff:$! ";

open OUT, '>', $output
  || die "$!";

while (<$tsv_fh>) {
    chomp;
    my @content = split /	/, $_;
    my @content_new = ();
    if ( $content[0] ne "#" ) {
        if (m/^chr/) {
            if ( $content[2] eq 'mRNA' ) {
                my $process = $content[8];
                $process =~ /ID=(.*)_mRNA;Name=/;
                $content_new[0] = $1;
                my $chr = $content[0];
                $chr =~ /^chr(.*)/;
                $content_new[1] = $1;
                $content_new[2] = $content[6];
                $content_new[3] = $content[3];
                $content_new[4] = $content[4];
                my $content_new = join ",", @content_new;
                open OUT, '>>', $output;
                print OUT $content_new, "\n";
            }
        }
    }
}

close $tsv_fh;
close OUT;

__END__
