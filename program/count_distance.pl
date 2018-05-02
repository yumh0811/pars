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

use Bio::Phylo::IO qw(parse);
#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl count_distance.pl --file $NAME.nwk  --list $NAME_strain_name.list

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'list|l=s'   => \my $list,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);


my $strain_name_fh;
open $strain_name_fh, "<",$list;

if ( !$output ) {
    $output = path($file)->basename;
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.csv";
}

  
my $tree = parse( -format => 'newick', -file => "$file" )->first;

# Bio::Phylo::Forest::Node
my $target = $tree->get_by_name("S288c");
my @content;
while (<$strain_name_fh>) {
    chomp;
    if($_ ne "S288c"){
        my $query = $tree->get_by_name("$_");
        my $dis = $query->calc_patristic_distance($target);
        $content[0] = $_; 
        $content[1] = $dis;

        open OUT,
           ">>", $output;
        printf OUT "$content[0]".','."%.9f"."\n",$content[1];
    }

}

close $strain_name_fh;
#close OUT1;
close OUT;

__END__
