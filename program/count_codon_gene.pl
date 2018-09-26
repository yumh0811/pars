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
				
				my @codon_occured = split "\\|", $snp[17] if defined( $snp[17] );
				my %count;
        my @unicodon = grep { ++$count{$_} < 2; } @codon_occured;
        my $codon_num = scalar( keys %count );
				my $codon_to;

            if ( $codon_num != 2 ) {
                $codon_to = join( '|', keys %count );
            }
            else {
                if ( $snp[10] eq 'Complex' ) {
                    $codon_to = join( '|', keys %count );
                }
                else {
                    #print YAML::Syck::Dump( \@unicodon );
                    if ( $codon_num != 0 ) {
                        if ( $count{ $unicodon[1] } ne $count{ $unicodon[0] } )
                        {
                            if ( $count{ $unicodon[1] } == $snp[11] ) {
                                $codon_to = "$unicodon[0]->$unicodon[1]";
                            }
                            else {
                                $codon_to = "$unicodon[1]->$unicodon[0]";
                            }
                        }else{
                    	  		if ( $snp[8] eq '+'){
                    	  				$codon_to = "$unicodon[0]->$unicodon[1]" if substr($unicodon[0], $snp[14], 1) eq substr($snp[10], 0 ,1);
                    	  				$codon_to = "$unicodon[1]->$unicodon[0]" if substr($unicodon[1], $snp[14], 1) eq substr($snp[10], 0 ,1);
                    	  		}else{
                    	  				my $mutant_to = $snp[10];              
                    	  		    $mutant_to =~ tr/ATGC/TACG/;
                    	  				$codon_to = "$unicodon[0]->$unicodon[1]" if substr($unicodon[0], $snp[14], 1) eq substr($mutant_to, 0 ,1);
                    	  				$codon_to = "$unicodon[1]->$unicodon[0]" if substr($unicodon[1], $snp[14], 1) eq substr($mutant_to, 0 ,1);
                    	  		}
                    	  }
                    }
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
