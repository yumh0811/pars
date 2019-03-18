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
$stopwatch->start_message("Process codon info...");

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

    if ( $snp[8] eq "gene" ) {
    
        splice @snp, 19, 0, "Codon_pos";
        my $snp = join ",", @snp;

        open OUT,'>>',$output;
        print OUT $snp, "\n";

    }
    else {
				
				my @codon = split "/", $snp[18] if defined( $snp[18] );
				my @mutant = split "->", $snp[11] if defined( $snp[11] );
				my $codon_to;
				my $codon_pos;
            
        if ( $snp[7] eq '+' ){
            if ($codon[0] =~ m/$mutant[0]/){
                $codon_to = "$codon[0]->$codon[1]";
                $codon_pos=index($codon[0],$mutant[0]);
            }else{
                $codon_to = "$codon[1]->$codon[0]";
                $codon_pos=index($codon[1],$mutant[0]);
            }
        }else{
        	  $mutant[0] =~ tr/ATGC/TACG/;
        	  $mutant[1] =~ tr/ATGC/TACG/;
            if ($codon[0] =~ m/$mutant[0]/){
                $codon_to = "$codon[0]->$codon[1]";
                $codon_pos=index($codon[0],$mutant[0]);
            }else{
                $codon_to = "$codon[1]->$codon[0]";
                $codon_pos=index($codon[1],$mutant[0]);
            }        	                 	
        }
        $codon_to = uc($codon_to);
        
        splice @snp, 18, 1, $codon_to if defined( $snp[18] );
        splice @snp, 19, 0, $codon_pos;
        my $snp = join ",", @snp;
        
        
        open OUT,'>>',$output;
        print OUT $snp, "\n";


    }
}

close $csv_fh;
close OUT;

__END__
