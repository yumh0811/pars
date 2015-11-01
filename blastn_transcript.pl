#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Basename;
use Bio::SearchIO;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;

#use lib "$FindBin::Bin/../lib";

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
);

my $file;
my $alignment_view = 0;    # blastall -m

my $identity = 99;
my $coverage = 1;

my $output;

my $man  = 0;
my $help = 0;

$|++;

GetOptions(
    'help|?'       => \$help,
    'man|m'        => \$man,
    'f|file=s'     => \$file,
    'm|view=s'     => \$alignment_view,
    'i|identity=i' => \$identity,
    'c|coverage=i' => \$coverage,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

my $view_name = {
    0 => "blast",         # Pairwise
    7 => "blastxml",      # BLAST XML
    9 => "blasttable",    # Hit Table
};

my $result_format = $view_name->{$alignment_view};

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Find paralog...");

if ( !$output ) {
    $output = basename($file);
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

            my $align_obj = $hsp->get_aln;    # a Bio::SimpleAlign object
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

=head1 NAME

    update_align_paralog.pl - Add additional paralog info to alignDB

=head1 SYNOPSIS

    update_align_paralog.pl [options]
      Options:
        --help               brief help message
        --man                full documentation
        --server             MySQL server IP/Domain name
        --db                 database name
        --username           username
        --password           password
        --datalib|da         blast database
        --megablast|mega     use megablast or not
        --view|v             blast output format

    update_align_paralog.pl -d=Nipvs9311 -da=nip_chro --mega=1 -v=9

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut
