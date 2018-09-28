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

use List::Util qw/min max/;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/count_structure_length_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --name ~/data/mrna-structure/result/Scer_n157_nonMosaic_Spar/PARS_cds_gene.csv --structure stem --output stem_length_cds.csv
    perl ~/Scripts/pars/program/count_structure_length_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --name ~/data/mrna-structure/result/Scer_n157_nonMosaic_Spar/PARS_cds_gene.csv --structure loop --output loop_length_cds.csv

=cut

Getopt::Long::GetOptions(
    'help|?'        => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'      => \my $file,
    'name|n=s'      => \my $name,
    'structure|n=s' => \my $structure,
    'output|o=s'    => \my $output,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Process folding info...");

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
print "Load $file\n";
my $gene_info_of = YAML::Syck::LoadFile($file);

my $csv_fh;
open $csv_fh, '<', $name;

open OUT, '>', $output;

while (<$csv_fh>) {
    chomp;
    s/"//g;
    my @snp = split /,/, $_;
    
    if ( $snp[0] eq "gene" ) {
        my @snp_new;
        $snp_new[0] = "gene";
        $snp_new[1] = "max";
        for ( my $i = 1 ; $i <= 200 ; $i = $i + 1 ) {
            $snp_new[$i+1] = $i,;
        }
        my $snp = join ",", @snp_new;
        open OUT, '>>', $output;
        print OUT $snp, "\n";

    }
    else {
        if ( $structure eq "stem" ) {
            my @snp_new;
            $snp_new[0] = $snp[0];
             
            if ( grep ( $snp[0], keys %{$gene_info_of} ) )
            {    # $snp[0] represents snp in which gene

                my $info = $gene_info_of->{ $snp[0] };

                my $stem = AlignDB::IntSpan->new;

                $stem->add( $info->{fold_left} );
                $stem->add( $info->{fold_right} );

                my @stem_length   = ();
                my @stem_position = $stem->sets();
                foreach my $ranges1 (@stem_position) {
                    my $max1    = $ranges1->max();
                    my $min1    = $ranges1->min();
                    my $length1 = $max1 - $min1 + 1;
                    push( @stem_length, $length1 );
                }
                my %stem_length;
                foreach (@stem_length) {
                    $stem_length{$_}++;
                }

                foreach ( sort { $a <=> $b } keys %stem_length ) {
                    $snp_new[$_+1] = $stem_length{$_};
                }
                my $max = max(keys %stem_length);
                $snp_new[1] = $max;
                #print $max."\n";
                for ( my $i = 1 ; $i <= 200 ; $i = $i + 1 ){
                    $snp_new[$i+1] = 0 unless (defined($snp_new[$i+1]));                    
                }
            }
           
            my $snp_processed = join ",", @snp_new;

            open OUT, '>>', $output;
            print OUT $snp_processed, "\n";

        }
        else {
            my @snp_new;
            $snp_new[0] = $snp[0];
            
            if ( grep ( $snp[0], keys %{$gene_info_of} ) )
            {    # $snp[0] represents snp in which gene

                my $info = $gene_info_of->{ $snp[0] };

                my $loop = AlignDB::IntSpan->new( $info->{fold_dot} );

                my @loop_length   = ();
                my @loop_position = $loop->sets();
                foreach my $ranges2 (@loop_position) {
                    my $max2    = $ranges2->max();
                    my $min2    = $ranges2->min();
                    my $length2 = $max2 - $min2 + 1;
                    push( @loop_length, $length2 );
                }
                my %loop_length;
                foreach (@loop_length) {
                    $loop_length{$_}++;
                }
                foreach ( sort { $a <=> $b } keys %loop_length ) {
                    $snp_new[$_+1] = $loop_length{$_};
                }
                my $max = max(keys %loop_length);
                $snp_new[1] = $max;
                #print $max."\n";
                for ( my $i = 1 ; $i <= 200 ; $i = $i + 1 ){
                    $snp_new[$i+1] = 0 unless (defined($snp_new[$i+1]));                    
                }
            }           

            my $snp_processed = join ",", @snp_new;

            open OUT, '>>', $output;
            print OUT $snp_processed, "\n";

        }
    }
}

close $csv_fh;
close OUT;

__END__
