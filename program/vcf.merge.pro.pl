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

use POSIX;
use Data::Dumper "Dumper";

use Statistics::ChisqIndep;
use POSIX;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS
    
    cd /Volumes/Backup/yumh/data/vcf/1011Matrix.gvcf
    perl merge_pro.pl --file chr1.merge.tsv --output chr1.merge.pro.tsv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

my $tsv_fh;
open $tsv_fh, '<', $file;

open OUT, '>', $output;

while (<$tsv_fh>) {
    chomp;
    s/\"//g;
    if (m/^name/) {
        my @content = split /\t/, $_;
        splice @content, 10, 1;
        $content[9]  = "mutant_to_vcf";
        $content[11] = "REF_vcf";
        $content[13] = "ALT_vcf";
        push @content, "freq_vcf";
        push @content, "minus";
        push @content, "chi";
        push @content, "p-valve";
        my $content = join "\t", @content;
        open OUT, '>>', $output;
        print OUT $content, "\n";
    }
    else {
        my @content = split /\t/, $_;
        if ( $content[8] == 0 ) {
            $content[9] = "$content[9]->$content[10]";
        }
        else {
            $content[9] = "$content[10]->$content[9]";
        }
        splice @content, 10, 1;
        if ( $content[8] == 0 ) {
            $content[14] = $content[12] if ( $content[3] eq $content[9] );
            $content[11] = $content[13] - $content[11] if ( $content[14] );
        	  $content[13] = $content[13] - $content[11] if ( $content[14] );
        }
        else {
            $content[14] = 1 - $content[12] if ( $content[3] eq $content[9] );
            $content[13] = $content[13] - $content[11] if ( $content[14] );
        }
        
        $content[15] = $content[5] - $content[14]  if ( $content[14] );
        if ( $content[14] ) {
            my $obs =
              [ [ $content[6], $content[7] ], [ $content[11], $content[13] ] ];
            my $chi = new Statistics::ChisqIndep;
            $chi->load_data($obs);
            $chi->print_summary();
            $content[16] = ${$chi}{'chisq_statistic'};
            $content[17] = ${$chi}{'p_value'};
        }
        my $content = join "\t", @content;
        open OUT, '>>', $output;
        print OUT $content, "\n";
    }
}

close $tsv_fh;
close OUT;

__END__
