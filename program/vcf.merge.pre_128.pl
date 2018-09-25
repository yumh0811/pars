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

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS
    
    cd /Volumes/Backup/yumh/data/vcf/1011Matrix.gvcf/
    perl vcf.merge.pre.pl --file /Users/yumh/data/mrna-structure/result/Scer_n128_Spar/data_SNPs_PARS_cds.csv --output chr1.pars.tsv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

my $csv_fh;
open $csv_fh, '<', $file;

open OUT, '>', $output;

while (<$csv_fh>) {
    chomp;
    s/\"//g;
    if (m/^name/) {
        my @content_new;
        $content_new[0] = "name";
        $content_new[1] = "gene";
        $content_new[2] = "structure";
        $content_new[3] = "mutant_to_pars";
        $content_new[4] = "freq";
        $content_new[5] = "freq_pars";
        $content_new[6] = "REF_pars";
        $content_new[7] = "ALT_pars";
        $content_new[8] = "target";
        my $content_new = join "\t", @content_new;
        open OUT,'>>',$output;;
        print OUT $content_new, "\n";
    }else{
        		my @content = split /,/, $_;
        		my @content_new;
        		$content_new[0] = $content[0];
        		$content_new[1] = $content[1];
        		$content_new[2] = $content[7];
        		#my @mutant = split "->", $content[10];
        		#if ($content[13] == 0){
        		#		$content_new[3] = "$mutant[0]<->$mutant[1]";
        		#}else{
        		#		$content_new[3] = "$mutant[1]<->$mutant[0]";
        		#}
        		$content_new[3] = $content[10];
        		$content_new[4] = $content[11];
        		$content_new[5] = $content[11]/128;
        		if ($content[13] == 0){
        				$content_new[6] = $content[11];
        		}else{
        				$content_new[6] = 128-$content[11];
        		}
        		$content_new[7] = 128-$content_new[6];
        		$content_new[8] = $content[13];
        		my $content_new = join "\t", @content_new;
        		my $content = join "\t", @content;
        		open OUT,'>>',$output;;
        		print OUT $content_new, "\n";
    }
}



close $csv_fh;
close OUT;

__END__



