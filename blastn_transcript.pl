#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use Bio::SearchIO;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

blastn_transcript.pl - get transcripts' genome location

=head1 SYNOPSIS

    perl blastn_transcript.pl [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     blast result file
        --view          -m  INT     blast output format, same as `blastall -m`
                                    0 => Pairwise
                                    7 => BLAST XML
                                    9 => Hit Table
        --identityi     -i  INT     default is [99]
        --coverage      -c  INT     default is [1]

=cut

GetOptions(
    'help|?'   => sub { HelpMessage(0) },
    'file|f=s' => \my $file,
    'view|m=s'     => \( my $alignment_view = 0 ),
    'identity|i=i' => \( my $identity       = 99 ),
    'coverage|c=i' => \( my $coverage       = 1 ),
) or HelpMessage(1);

$|++;
my $output;

my $view_name = {
    0 => "blast",         # Pairwise
    7 => "blastxml",      # BLAST XML
    9 => "blasttable",    # Hit Table
};

my $result_format = $view_name->{$alignment_view};

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Find transcripts...");

if ( !$output ) {
    $output = path($file)->basename;
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.blast.tsv";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
open my $blast_fh, '<', $file;
open my $out_fh,   ">", $output;
open my $skip_fh,  ">", "$output.skip";

my $searchio = Bio::SearchIO->new(
    -format => $result_format,
    -fh     => $blast_fh,
);

QUERY: while ( my $result = $searchio->next_result ) {
    my $query_name   = $result->query_name;
    my $query_length = $result->query_length;
    print "name $query_name\tlength $query_length\n";
    while ( my $hit = $result->next_hit ) {
        my $hit_name = $hit->name;

        my $hit_length = $hit->length;
        my $query_set  = AlignDB::IntSpan->new;
        while ( my $hsp = $hit->next_hsp ) {

            # process the Bio::Search::HSP::HSPI object
            my $hsp_length = $hsp->length( ['query'] );

            # use "+" for default strand
            # -1 = Minus strand, +1 = Plus strand
            my ( $query_strand, $hit_strand ) = $hsp->strand("list");
            my $hsp_strand = "+";
            if ( $query_strand + $hit_strand == 0 and $query_strand != 0 ) {
                $hsp_strand = "-";
            }

            my $align_obj   = $hsp->get_aln;                             # a Bio::SimpleAlign object
            my ($query_obj) = $align_obj->each_seq_with_id($query_name);
            my ($hit_obj)   = $align_obj->each_seq_with_id($hit_name);

            my $q_start = $query_obj->start;
            my $q_end   = $query_obj->end;
            if ( $q_start > $q_end ) {
                ( $q_start, $q_end ) = ( $q_end, $q_start );
            }
            $query_set->add_range( $q_start, $q_end );

            my $h_start = $hit_obj->start;
            my $h_end   = $hit_obj->end;
            if ( $h_start > $h_end ) {
                ( $h_start, $h_end ) = ( $h_end, $h_start );
            }

            my $query_coverage = $query_set->size / $query_length;
            my $hsp_identity   = $hsp->percent_identity;

            if ( $query_coverage < $coverage or $hsp_identity < $identity ) {
                print {$skip_fh} join "\t", $query_name, $query_length,
                    $hit_name,
                    $query_coverage, $hsp_identity;
                print {$skip_fh} "\n";
            }
            else {
                print {$out_fh} join "\t", $query_name, $query_length,
                    $hit_name,
                    $h_start, $h_end, $hsp_strand;
                print {$out_fh} "\n";
            }

            next QUERY;    # only get the first hsp of the first hit
        }
    }
}

close $blast_fh;
close $out_fh;
close $skip_fh;

$stopwatch->end_message;

exit;

__END__
