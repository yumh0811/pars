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
    perl ~/Scripts/pars/program/vcf.cut.pl --file 1011Matrix.gvcf --output 1011Matrix.tsv

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
    if (m/^\#\#/) {
        next;
    }else{
        my @content = split /\t/, $_;
        splice @content, 8;
        my $content = join "\t", @content;
        open OUT,'>>',$output;;
        print OUT $content, "\n";
    }
}

close $tsv_fh;
close OUT;

__END__