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

    perl ~/Scripts/pars/program/count_codon_gene.pl --origin data_SNPs_PARS_cds.update.csv --output data_SNPs_PARS_cds.update_codon.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'origin|or=s'   => \my $origin,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Process folding info...");

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

my $csv_fh;
open $csv_fh,'<',$origin;

open OUT,'>',$output;

while (<$csv_fh>) {
    chomp;
    s/"//g;
    my @snp = split /,/, $_;

    if ( $snp[1] eq "gene" ) {
    
        my $snp = join ",", @snp;

        open OUT,'>>',$output;
        print OUT $snp, "\n";

    }
    else {
        if ( $snp[8] eq '-' ) {             
						$snp[10] =~ tr/ATGC/TACG/;
        }
				
				my @codon_occured = split "\\|", $snp[17] if defined( $snp[17] );
				my %count;
        my @unicodon = grep { ++$count{$_} < 2; } @codon_occured;
        my $codon_num = scalar( keys %count );
				my $codon_to;

            if ( $codon_num != 2 ) {
                $codon_to = 'o1';
            }
            else {
                if ( $snp[10] eq 'Complex' ) {
                    $codon_to = 'o2';
                }
                else {
                    #print YAML::Syck::Dump( \@unicodon );
                    $codon_to = "$unicodon[0]\|$unicodon[1]"
                }
            }

        splice @snp, 17, 1, $codon_to if defined( $snp[17] );
        my $snp = join ",", @snp;
        
        
        open OUT,'>>',$output;
        print OUT $snp, "\n";


    }
}

close $csv_fh;
close OUT;

__END__
