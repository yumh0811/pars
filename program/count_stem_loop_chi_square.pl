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

use Statistics::ChisqIndep;
use POSIX;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file ~/data/mrna-structure/result/$NAMEfreq_10/PARS_nsy_stat_freq_10.csv --output ~/data/mrna-structure/result/$NAMEfreq_10/PARS_nsy_stat_freq_10_chi_square.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

  
my $csv_fh;
open $csv_fh, '<',$file;;

open OUT,'>',$output;

my ($stem_AT_GC,$stem_GC_AT,$loop_AT_GC,$loop_GC_AT);

while (<$csv_fh>) {
    chomp;
    my @content = split /,/, $_;
    if ( $content[0] ne '"structure"' ) {
        if ( $content[0] eq '"stem"' ) {
            $stem_AT_GC = $content[1];
            $stem_GC_AT = $content[2];
            $content[3] = $content[1] / ( $content[1] + $content[2] );
            my $content = join ",", @content;
            open OUT,'>>',$output;
            print OUT $content, "\n";
            next;
        }
        else {
            $loop_AT_GC = $content[1];
            $loop_GC_AT = $content[2];
            $content[3] = $content[1] / ( $content[1] + $content[2] );
            my $obs =
              [ [ $stem_AT_GC, $stem_GC_AT ], [ $loop_AT_GC, $loop_GC_AT ] ];
            my $chi = new Statistics::ChisqIndep;
            $chi->load_data($obs);
            $chi->print_summary();
            my $chisq = ${$chi}{'chisq_statistic'};
            my $P     = ${$chi}{'p_value'};
            $content[4] = $chisq;
            $content[5] = $P;
            my $content = join ",", @content;
            open OUT,'>>',$output;
            print OUT $content, "\n";
            next;

        }
    }
    else {
        $content[3] = "AT_GC_ratio";
        $content[4] = "X2";
        $content[5] = "p-value";
        my $content = join ",", @content;
        open OUT,'>>',$output;
        print OUT $content, "\n";

    }
}

close $csv_fh;
close OUT;

__END__
