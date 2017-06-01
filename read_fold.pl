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
use App::RL::Common;
use App::Fasops::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/read_fold.pl \
        --pars ~/data/mrna-structure/PARS10/pubs/PARS10/data \
        --gene ~/data/mrna-structure/process/sce_genes.blast.tsv \
        --pos $NAME.snp.gene.pos.txt  \
        > fail_pos.txt

=cut

Getopt::Long::GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'pars=s' => \my $dir_pars,
    'gene=s' => \( my $file_gene = "sce_genes.blast.tsv" ),
    'pos=s'  => \( my $file_pos  = ".snp.gene.pos.txt" ),
    'pos2=s'     => \my $file_pos2,    # for indel
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Process folding data...");

if ( !$output ) {
    $output = path($file_pos)->basename;
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.gene_variation";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $gene_info_of = {};

# read gene chr pos
# YKL152C	872	chrXI	163525	164396	-
{
    print "Read $file_gene\n";
    open my $fh_gene, '<', $file_gene;
    while ( my $line = <$fh_gene> ) {
        chomp $line;
        my ( $gene, undef, $chr, $start, $end, $strand ) = split /\t/, $line;
        $gene_info_of->{$gene} = {
            chr    => $chr,
            start  => $start,
            end    => $end,
            strand => $strand,
            length => $end - $start + 1,
            vars   => [],                  # prepare for variation section
            vars2  => [],                  # prepare for variation section
        };
    }
    close $fh_gene;
    print "\n";
}

# read gene fold
{
    my $file_fold = path( $dir_pars, "sce_genes_folded.tab" )->absolute->stringify;
    print "Read $file_fold\n";

    my ( @non_match, @non_exist );

    open my $fh_fold, '<', $file_fold;
    while ( my $line = <$fh_fold> ) {
        chomp $line;
        my ( $gene, $seq, $fold ) = split /\t/, $line;
        if ( exists $gene_info_of->{$gene} ) {
            if ( length $seq != $gene_info_of->{$gene}{length} ) {
                delete $gene_info_of->{$gene};
                push @non_match, $gene;
            }
            else {
                $gene_info_of->{$gene}{seq}  = $seq;
                $gene_info_of->{$gene}{fold} = $fold;
            }
        }
        else {
            push @non_exist, $gene;
        }
    }
    close $fh_fold;

    if (@non_match) {
        print scalar @non_match, " genes don't match length\n";
        print join( " ", @non_match, "\n" );
    }

    if (@non_exist) {
        print scalar @non_exist, " genes don't exist in $file_gene\n";
        print join( " ", @non_exist, "\n" );
    }
    print "\n";
}

# read gene score
{
    my $file_score = path( $dir_pars, "sce_Score.tab" )->absolute->stringify;
    print "Read $file_score\n";

    my (@non_exist);

    open my $fh, '<', $file_score;
    while ( my $line = <$fh> ) {
        chomp $line;
        my ( $gene, undef, $score ) = split /\t/, $line;

        if ( exists $gene_info_of->{$gene} ) {
            my @scores = split /\;/, $score;
            $gene_info_of->{$gene}{mF_score}    = App::Fasops::Common::mean(@scores);
            $gene_info_of->{$gene}{pars_scores} = [@scores];
        }
        else {
            push @non_exist, $gene;
        }
    }
    close $fh;

    if (@non_exist) {
        print scalar @non_exist, " genes don't exist in $file_gene\n";
        print join( " ", @non_exist, "\n" );
    }
    print "\n";
}

# read variation pos, handle snps
# chrI	35070	35070
{
    print "Read $file_pos\n";

    open my $fh_pos, '<', $file_pos;
    while ( my $line = <$fh_pos> ) {
        next if substr( $line, 0, 1 ) eq "#";
        chomp $line;

        my $info = App::RL::Common::decode_header($line);
        next unless App::RL::Common::info_is_valid($info);

        my @genes = grep {
                    $gene_info_of->{$_}{start} <= $info->{start}
                and $info->{end} <= $gene_info_of->{$_}{end}
            }
            grep { $gene_info_of->{$_}{chr} eq $info->{chr} } keys %{$gene_info_of};

        my $count = scalar @genes;
        if ( $count != 1 ) {
            print join "\t", $count, $line, @genes, "\n";
        }
        else {
            my $gene = $genes[0];
            my $name = $gene_info_of->{$gene}{chr} . ":" . $info->{start};
            push @{ $gene_info_of->{$gene}{vars} }, { chr_pos => $info->{start}, name => $name, };
        }
    }
    close $fh_pos;
    print "\n";
}

# read variation pos
if ($file_pos2) {
    print "Read $file_pos2\n";

    open my $fh_pos, '<', $file_pos2;
    while ( my $line = <$fh_pos> ) {
        next if substr( $line, 0, 1 ) eq "#";
        chomp $line;

        my $info = App::RL::Common::decode_header($line);
        next unless App::RL::Common::info_is_valid($info);

        my @genes = grep {
                    $gene_info_of->{$_}{start} <= $info->{start}
                and $info->{end} <= $gene_info_of->{$_}{end}
            }
            grep { $gene_info_of->{$_}{chr} eq $info->{chr} } keys %{$gene_info_of};

        my $count = scalar @genes;
        if ( $count != 1 ) {
            print join "\t", $count, $line, @genes, "\n";
        }
        else {
            my $gene = $genes[0];
            my $name = $gene_info_of->{$gene}{chr} . ":" . $info->{start};
            push @{ $gene_info_of->{$gene}{vars2} },
                { chr_start => $info->{start}, chr_end => $info->{end}, name => $name, };
        }
    }
    close $fh_pos;
    print "\n";
}

YAML::Syck::DumpFile( "$output.yml", $gene_info_of );

$stopwatch->end_message;

exit;

__END__
