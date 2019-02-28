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
use Data::Dumper "Dumper";

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/count_per_gene_ACGT_percent.pl --file data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds_per_gene_ATGC.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

my $csv_fh;
open $csv_fh, '<', $file;

open OUT, '>', $output;

my %mutation_AT_GC =
  ( 'A->C' => 0, 'A->G' => 0, 'T->C' => 0, 'T->G' => 0 );    #A/T->G/C
my %mutation_GC_AT =
  ( 'G->A' => 0, 'G->T' => 0, 'C->A' => 0, 'C->T' => 0 );    #G/C->A/T
my %mutation_other =
  ( 'G->C' => 0, 'C->G' => 0, 'T->A' => 0, 'A->T' => 0 );    #G<->C,A<->T
my %gene_name;
my %count_SNPs = (
    'stem_AT_GC' => 0,
    'stem_GC_AT' => 0,
    'loop_AT_GC' => 0,
    'loop_GC_AT' => 0
);

while (<$csv_fh>) {
    chomp;
    s/\"//g;
    my @content = split /,/, $_;
    my @content_new;
    if ( $content[8] eq "gene" ) {
        $content_new[0] = "gene";
        $content_new[1] = "stem_AT_GC";
        $content_new[2] = "stem_GC_AT";
        $content_new[3] = "stem_total";
        $content_new[4] = "loop_AT_GC";
        $content_new[5] = "loop_GC_AT";
        $content_new[6] = "loop_total";
        #$content_new[7] = "stem_AT_GC_ratio";
        #$content_new[8] = "loop_AT_GC_ratio";
        #$content_new[9] = "X2";
        #$content_new[10] = "p_value";
        my $content_new = join ",", @content_new;
        open OUT,'>>',$output;
        print OUT $content_new, "\n";

    }
    else {

        if ( exists $gene_name{ $content[8] } ) {

            if ( $content[6] eq "stem" ) {
                if ( exists $mutation_AT_GC{ $content[11] } ) {
                    $count_SNPs{'stem_AT_GC'} += 1;
                }
                elsif ( exists $mutation_GC_AT{ $content[11] } ) {
                    $count_SNPs{'stem_GC_AT'} += 1;
                }elsif ( exists $mutation_other{ $content[11] } ) {
                	  $count_SNPs{'stem_other'} += 1;
                }
            }
            else {
                if ( exists $mutation_AT_GC{ $content[11] } ) {
                    $count_SNPs{'loop_AT_GC'} += 1;
                }
                elsif ( exists $mutation_GC_AT{ $content[11] } ) {
                    $count_SNPs{'loop_GC_AT'} += 1;
                }elsif ( exists $mutation_other{ $content[11] } ) {
                	  $count_SNPs{'loop_other'} += 1;
                }
            }

            #my %hash = %count_SNPs;
            $gene_name{ $content[8] } = {%count_SNPs};

            #print $content[1],"\n";
            #print Dumper(\%gene_name);sleep 2;
        }
        else {

            %count_SNPs = (
                'stem_AT_GC' => 0,
                'stem_GC_AT' => 0,
                'loop_AT_GC' => 0,
                'loop_GC_AT' => 0,
                'stem_other' => 0,
                'loop_other' => 0
            );

            if ( $content[6] eq "stem" ) {
                if ( exists $mutation_AT_GC{ $content[11] } ) {
                    $count_SNPs{'stem_AT_GC'} += 1;
                }
                elsif ( exists $mutation_GC_AT{ $content[11] } ) {
                    $count_SNPs{'stem_GC_AT'} += 1;
                }elsif ( exists $mutation_other{ $content[11] } ) {
                	  $count_SNPs{'stem_other'} += 1;
                }
            }
            elsif ( $content[6] eq "loop" ) {
                if ( exists $mutation_AT_GC{ $content[11] } ) {
                    $count_SNPs{'loop_AT_GC'} += 1;
                }
                elsif ( exists $mutation_GC_AT{ $content[11] } ) {
                    $count_SNPs{'loop_GC_AT'} += 1;
                }elsif ( exists $mutation_other{ $content[11] } ) {
                	  $count_SNPs{'loop_other'} += 1;
                }
            }

            #my %hash = %count_SNPs;
            $gene_name{ $content[8] } = {%count_SNPs};

            #print $content[1],"\n";
            #print Dumper(\%gene_name);sleep 2;

        }
    }

}

#print "=============";
#print  scalar(keys %gene_name);

foreach my $keys ( sort keys %gene_name ) {
    my @content_new;
    $content_new[0] = $keys;

    my $stem_AT_GC = $gene_name{$keys}{'stem_AT_GC'};
    my $stem_GC_AT = $gene_name{$keys}{'stem_GC_AT'};
    my $loop_AT_GC = $gene_name{$keys}{'loop_AT_GC'};
    my $loop_GC_AT = $gene_name{$keys}{'loop_GC_AT'};
    my $stem_other = $gene_name{$keys}{'stem_other'};
    my $loop_other = $gene_name{$keys}{'loop_other'};
    my $stem_total = $stem_AT_GC + $stem_GC_AT + $stem_other;
    my $loop_total = $loop_AT_GC + $loop_GC_AT + $loop_other;
    my $stem_AT_GC_ratio;
    my $loop_AT_GC_ratio;

#    if ( $stem_AT_GC + $stem_GC_AT + $stem_other!= 0 ) {
#        $stem_AT_GC_ratio = $stem_AT_GC / ( $stem_AT_GC + $stem_GC_AT + $stem_other );
#    }
#    if ( $loop_AT_GC + $loop_GC_AT + $loop_other!= 0 ) {
#        $loop_AT_GC_ratio = $loop_AT_GC / ( $loop_AT_GC + $loop_GC_AT + $loop_other);
#    }
    
    if ( $stem_AT_GC + $loop_AT_GC != 0 ) {
        $stem_AT_GC_ratio = $stem_AT_GC / ( $stem_AT_GC + $loop_AT_GC );
    }
    if ( $stem_AT_GC + $loop_AT_GC != 0 ) {
        $loop_AT_GC_ratio = $loop_AT_GC / ( $stem_AT_GC + $loop_AT_GC );
    }

    my $obs = [ [ $stem_AT_GC, $stem_GC_AT ], [ $loop_AT_GC, $loop_GC_AT ] ];
    my $chi = new Statistics::ChisqIndep;
    $chi->load_data($obs);
    $chi->print_summary();
    my $chisq = ${$chi}{'chisq_statistic'};
    my $P     = ${$chi}{'p_value'};

    $content_new[1] = $stem_AT_GC;
    $content_new[2] = $stem_GC_AT;
    $content_new[3] = $stem_total;
    $content_new[4] = $loop_AT_GC;
    $content_new[5] = $loop_GC_AT;
    $content_new[6] = $loop_total;
    #$content_new[7] = $stem_AT_GC_ratio;
    #$content_new[8] = $loop_AT_GC_ratio;
    #$content_new[9] = $chisq;
    #$content_new[10] = $P;

    my $content_new = join ",", @content_new;
    open OUT,'>>',$output;;
    print OUT $content_new, "\n";
}

close $csv_fh;
close OUT;

__END__
