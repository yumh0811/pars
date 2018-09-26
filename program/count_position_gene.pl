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

use List::Util qw/min/;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --origin data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds_pos.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
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
print "Load $file\n";
my $gene_info_of = YAML::Syck::LoadFile($file);

my $csv_fh;
open $csv_fh,'<',$origin;

open OUT,'>',$output;

while (<$csv_fh>) {
    chomp;
    s/"//g;
    my @snp = split /,/, $_;

    if ( $snp[1] eq "gene" ) {
        $snp[61] = "snp_pos";
        $snp[62] = "island_length";
        my $snp = join ",", @snp;

        open OUT,'>>',$output;
        print OUT $snp, "\n";

    }
    else {

        if ( $snp[7] eq "stem" ) {

            if ( grep ( $snp[1], keys %{$gene_info_of} ) )
            {    # $snp[1] represents snp in which gene

                my $info = $gene_info_of->{ $snp[1] };

                my $stem = AlignDB::IntSpan->new;

                $stem->add( $info->{fold_left} );
                $stem->add( $info->{fold_right} );

                #print YAML::Syck::Dump($stem);
                #my $loop = AlignDB::IntSpan -> new($info->{fold_dot});

                my @stem_junction = $stem->ranges();
                my @positon       = ();

                foreach my $junction (@stem_junction) {
                    push( @positon, abs( $junction - $snp[4] ) )
                      ;    # $snp[4] represents snp absolute position
                    
                    
                }

                my $min = List::Util::min @positon;

                push( @snp, $min );
                
                my $island = $stem ->find_islands( $snp[4] );
                my $length = $island -> cardinality();
                push( @snp, $length );
                
            }

            my $snp_processed = join ",", @snp;

            open OUT,'>>',$output;
            print OUT $snp_processed, "\n";
        }
        else {
            if ( grep ( $snp[1], keys %{$gene_info_of} ) )
            {    # $snp[1] represents snp in which gene

                my $info = $gene_info_of->{ $snp[1] };

                my $stem = AlignDB::IntSpan->new;

                $stem->add( $info->{fold_left} );
                $stem->add( $info->{fold_right} );


                my $loop = AlignDB::IntSpan -> new($info->{fold_dot});

                my @stem_junction = $stem->ranges();
                my @positon       = ();

                foreach my $junction (@stem_junction) {
                    push( @positon, abs( $junction - $snp[4] ) )
                      ;    # $snp[4] represents snp absolute position
                }

                my $min = List::Util::min @positon;

                push( @snp, -$min );

                my $island = $loop -> find_islands( $snp[4] );
                my $length = $island -> cardinality();
                push( @snp, $length );

            }

            my $snp_processed = join ",", @snp;

            open OUT,'>>',$output;
            print OUT $snp_processed, "\n";
        }

    }
}

close $csv_fh;
close OUT;

__END__
