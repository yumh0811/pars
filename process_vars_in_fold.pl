#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/process_vars_in_fold.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.yml

=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'file|f=s'   => \my $file,
    'output|o=s' => \my $output,
) or HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Process folding info...");

if ( !$output ) {
    $output = path($file)->basename;
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.gene_variation";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
print "Load $file\n";
my $gene_info_of = LoadFile($file);

my $fold_class = {
    fold_dot   => '.',
    fold_left  => '(',
    fold_right => ')',
};

print "Process...\n";
for my $gene ( sort keys %{$gene_info_of} ) {
    my $info = $gene_info_of->{$gene};    # shortcut

    # parse fold string
    for my $class ( keys %{$fold_class} ) {
        my $sites = char_sites( $info->{fold}, $fold_class->{$class} );
        my $set = AlignDB::IntSpan->new;
        for ( @{$sites} ) {
            $set = $set->add( $_->{start} .. $_->{end} );
        }
        $info->{$class}               = $set;
        $info->{ $class . "_length" } = $set->size;
        $info->{ $class . "_vars" }   = 0;
    }

    # position of the paired base
    # for each '(', find forwrod balanced ')'
    # and record for each side of pair
    my $pair_pos_of = {};
    for my $pos ( 1 .. $info->{length} ) {
        my $class = substr $info->{fold}, $pos - 1, 1;

        if ( $class eq '.' ) {
            $pair_pos_of->{$pos} = undef;
        }
        elsif ( $class eq '(' ) {
            my $pair_pos = find_pair( $info->{fold}, $pos );
            if ($pair_pos) {
                $pair_pos_of->{$pos}      = $pair_pos;
                $pair_pos_of->{$pair_pos} = $pos;
            }
        }
    }
    $info->{pair_pos} = $pair_pos_of;

    my $chr_set
        = AlignDB::IntSpan->new( $info->{start} . "-" . $info->{end} );    # gene set in chr, no gap
    $info->{chr_set} = $chr_set;

    # vars count
    my $count = scalar @{ $info->{vars} };
    $info->{var_count} = $count;
    next if $count == 0;

    for my $i ( 0 .. $count - 1 ) {
        my $var = $info->{vars}[$i];                                       # shortcut, too

        # transform vars coordinate to local
        my $gene_pos = $chr_set->index( $var->{chr_pos} );
        if ( $info->{strand} eq '-' ) {
            $gene_pos = $info->{length} - $gene_pos + 1;
        }
        $var->{gene_pos} = $gene_pos;

        # var fold class
        for my $class ( keys %{$fold_class} ) {
            if ( $info->{$class}->contains($gene_pos) ) {
                $var->{fold_class} = $class;
                my $island = $info->{$class}->find_islands($gene_pos);
                $var->{fold_length} = $island->size;

                # count vars in each class
                $info->{ $class . "_vars" }++;
                last;
            }
        }

        # var base in gene
        $var->{gene_base} = substr $info->{seq}, $gene_pos - 1, 1;

        # var PARS score
        $var->{pars} = $info->{pars_scores}[ $gene_pos - 1 ];

        my $paired_pos = $info->{pair_pos}{$gene_pos};
        $var->{pair_base} = undef;    # for '.', remain undef
        if ($paired_pos) {
            $var->{pair_base} = substr $info->{seq}, $paired_pos - 1, 1;
        }
    }
}

#----------------------------#
# outputs
#----------------------------#
# gene
{
    open my $out_fh, '>', "$output.fold_class.tsv";
    my @heads
        = map { ( $_ . "_length", $_ . "_vars" ) } sort keys %{$fold_class};
    print {$out_fh} join( "\t", "gene", "length", "mF", @heads ), "\n";
    for my $gene ( sort keys %{$gene_info_of} ) {
        my $info = $gene_info_of->{$gene};    # shortcut

        # transform AlignDB::IntSpan object to runlist
        $info->{chr_set} = $info->{chr_set}->runlist;
        for my $class ( keys %{$fold_class} ) {
            $info->{$class} = $info->{$class}->runlist;
        }

        print {$out_fh} join "\t", $gene, $info->{length}, $info->{mF_score};
        print {$out_fh} "\t" . $info->{$_} for @heads;
        print {$out_fh} "\n";

        # this array is too large
        delete $info->{pars_scores};
    }
    close $out_fh;
}

# vars
{
    open my $out_fh, '>', "$output.var_pars.tsv";
    my @heads = qw{ name fold_length gene_base gene_pos pars pair_base};
    print {$out_fh} join( "\t", @heads, "gene", "structure", "strand" ), "\n";
    for my $gene ( sort keys %{$gene_info_of} ) {
        my $info = $gene_info_of->{$gene};    # shortcut

        for my $var ( @{ $info->{vars} } ) {
            for (@heads) {
                if ( defined $var->{$_} ) {
                    print {$out_fh} $var->{$_} . "\t";
                }
                else {
                    #print " " x 4, $var->{name}, " $_ not defined\n";
                    print {$out_fh} "\t";
                }
            }
            print {$out_fh} $gene . "\t";
            print {$out_fh} $var->{fold_class} eq "fold_dot" ? "loop" : "stem", "\t";
            print {$out_fh} $info->{strand};
            print {$out_fh} "\n";
        }
    }
    close $out_fh;
}

DumpFile( "$output.process.yml", $gene_info_of );

$stopwatch->end_message;

exit;

sub char_sites {
    my $string = shift;
    my $char   = shift;

    my $legnth = length $string;
    my @sites;

    my $offset = 0;
    my $start  = 0;
    my $end    = 0;
    for my $pos ( 1 .. $legnth ) {
        my $pos_char = substr( $string, $pos - 1, 1 );
        if ( $pos_char eq $char ) {
            if ( $offset == 0 ) {
                $start = $pos;
            }
            $offset++;
        }
        else {
            if ( $offset != 0 ) {
                $end = $pos - 1;
                my $length = $end - $start + 1;

                push @sites,
                    {
                    length => $length,
                    start  => $start,
                    end    => $end,
                    };
            }
            $offset = 0;
        }
    }
    if ( $offset != 0 ) {
        $end = $legnth;
        my $length = $end - $start + 1;

        push @sites,
            {
            length => $length,
            start  => $start,
            end    => $end,
            };
    }
    @sites = sort { $a->{start} <=> $b->{start} } @sites;

    return \@sites;
}

sub find_pair {
    my $string = shift;
    my $pos    = shift;

    my $depth = 1;    # starting depth
    my $pair_pos;

    while ( $pos++ ) {
        my $char = substr $string, $pos - 1, 1;
        if ( $char eq '(' ) {
            $depth++;
            next;
        }
        elsif ( $char eq ')' ) {
            $depth--;
            if ( $depth == 0 ) {
                $pair_pos = $pos;
                last;
            }
        }
    }

    return $pair_pos;
}

__END__
