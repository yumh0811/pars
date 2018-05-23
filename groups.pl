#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive qw();
use FindBin;
use YAML::Syck qw();
use Path::Tiny qw();

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
(
    #@type Getopt::Long::Descriptive::Opts
    my $opt,

    #@type Getopt::Long::Descriptive::Usage
    my $usage,
    )
    = Getopt::Long::Descriptive::describe_options(
    'perl %c [options] <snp.csv>',
    [ 'help|h', 'print usage message and exit' ],
    [],
    [ 'strain|s=s', 'list of strains', { default => "$FindBin::RealBin/strain_non_mosaic.csv" }, ],
    [ 'group|g=s',  'list of groups',  { default => "$FindBin::RealBin/group_phylo.tsv" }, ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

#----------------------------------------------------------#
# Run
#----------------------------------------------------------#
my @strains = Path::Tiny::path( $opt->{strain} )->lines( { chomp => 1, } );

my @groups;
my $strain_of = {};
for my $line ( Path::Tiny::path( $opt->{group} )->lines( { chomp => 1, } ) ) {
    my @fields = split /\t/, $line;
    next if $fields[0] =~ /^#/;
    next if scalar @fields < 2;

    my ( $group, $list ) = @fields;
    push @groups, $group;
    $strain_of->{$group} = [ split /,/, $list ];
}

#warn YAML::Syck::Dump $strain_of;

print join( ",", ( "name", @groups ) ), "\n";

for my $line ( Path::Tiny::path( $ARGV[0] )->lines( { chomp => 1, } ) ) {
    my @fields = split /,/, $line;
    next if $fields[0] eq "name";

    my %count = map { $_ => 0 } @strains;

    my @bases = split //, $fields[3];
    if ( @bases != @strains ) {
        warn "bases and strains don't match each others\n";
        warn YAML::Syck::Dump {
            bases   => scalar @bases,
            strains => scalar @strains,
        };

        next;
    }

    for my $i ( 0 .. $#strains ) {
        if ( $bases[$i] ne $fields[2] ) {
            $count{ $strains[$i] }++;
        }
    }

    my %total = map { $_ => 0 } @groups;
    my %mean  = map { $_ => 0 } @groups;

    for my $group (@groups) {
        for my $strain ( @{ $strain_of->{$group} } ) {
            $total{$group} += $count{$strain};
        }
        $mean{$group} = $total{$group} / scalar( @{ $strain_of->{$group} } );
    }

    print join( ",", ( $fields[0], @mean{@groups} ) ), "\n";
}
