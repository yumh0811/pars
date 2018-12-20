#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use YAML::Syck qw();

use Path::Tiny qw();

use AlignDB::IntSpan;

use List::Util qw/min/;


my $file_snp    = shift || "filiter_snp.csv";
my $file_strain = shift || "strainlist.csv";

my @asian = (
    "bioethanol001", "bioethanol003", "bioethanol004", "sake001",
    "sake003",       "sake004",       "sake005",       "sake006",
    "sake007",       
);
my @wine = (
    "beer014", "beer020",    "beer024",    "beer030",    "beer033", "beer088",
    "sake002", "spirits002", "spirits004", "spirits011", "wine001", "wine003",
    "wine004", "wine005",    "wine006",    "wine007",    "wine009", "wine010",
    "wine011", "wine013",    "wine014",    "wine015",    "wine017", "wine018"
);
my @beer1 = (
    "beer001", "beer007", "beer008",    "beer009", "beer010", "beer012",
    "beer015", "beer016", "beer022",    "beer026", "beer031", "beer036",
    "beer037", "beer041", "beer043",    "beer044", "beer045", "beer046",
    "beer047", "beer048", "beer049",    "beer050", "beer051", "beer052",
    "beer053", "beer054", "beer055",    "beer056", "beer064", "beer065",
    "beer066", "beer067", "beer068",    "beer069", "beer070", "beer071",
    "beer073", "beer075", "beer076",    "beer077", "beer078", "beer079",
    "beer081", "beer082", "beer087",    "beer089", "beer090", "beer094",
    "beer095", "beer096", "beer097",    "beer098", "beer099", "beer100",
    "beer101", "beer102", "spirits005", "wine012"
);
my @beer2 = (
    "beer003", "beer004", "beer011", "beer013", "beer021", "beer027", "beer032",
    "beer034", "beer040", "beer059", "beer062", "beer063", "beer080", "beer083",
    "beer084", "beer085", "beer086", "beer091", "beer092"
);
my @mixed = (
    "beer005",    "beer006",    "beer023",  "beer025",  "beer028",  "beer029",
    "beer038",    "beer061",    "bread001", "bread002", "bread003", "bread004",
    "spirits001", "spirits003", "wild005",  "wild006",  "wild007"
);

my @strains = Path::Tiny::path($file_strain)->lines( { chomp => 1, } );

#warn YAML::Syck::Dump \@strains;

print join( ",", ( "name", "gene", "asian", "wine", "beer1", "beer2", "mixed", "vcf", "minus", "p-value" ) ), "\n";

for my $line ( Path::Tiny::path($file_snp)->lines( { chomp => 1, } ) ) {
    my @fields = split /,/, $line;
    next if ( $fields[0] eq "name" );

    my %count = map { $_ => 0 } @strains;

    my @bases = split //, $fields[20];
    if ( @bases != @strains ) {
        warn "bases and strains don't match each others\n";
        warn YAML::Syck::Dump {
            bases   => scalar @bases,
            strains => scalar @strains,
        };

        next;
    }

    for my $i ( 0 .. $#strains ) {
        if ( $bases[$i] ne $fields[19] ) {
            $count{ $strains[$i] }++;
        }
    }

    my $asian_total = 0;
    foreach (@asian) {
        $asian_total += $count{$_};
    }
    my $asian_mean = $asian_total / scalar(@asian);

    my $wine_total = 0;
    foreach (@wine) {
        $wine_total += $count{$_};
    }
    my $wine_mean = $wine_total / scalar(@wine);

    my $beer1_total = 0;
    foreach (@beer1) {
        $beer1_total += $count{$_};
    }
    my $beer1_mean = $beer1_total / scalar(@beer1);

    my $beer2_total = 0;
    foreach (@beer2) {
        $beer2_total += $count{$_};
    }
    my $beer2_mean = $beer2_total / scalar(@beer2);

    my $mixed_total = 0;
    foreach (@mixed) {
        $mixed_total += $count{$_};
    }
    my $mixed_mean = $mixed_total / scalar(@mixed);

    print
        join( ",", ( $fields[0], $fields[1], $asian_mean, $wine_mean, $beer1_mean, $beer2_mean, $mixed_mean ,$fields[14], $fields[15], $fields[17])),
        "\n";
}