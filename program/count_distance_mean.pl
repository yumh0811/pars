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

    perl ~/Scripts/pars/program/count_distance_mean.pl --file ${NAME}_distance/{}.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
) or Getopt::Long::HelpMessage(1);
my $gene_name = path($file)->basename('.csv');

my @asian = ("bioethanol001","bioethanol003","bioethanol004","sake001","sake003","sake004","sake005","sake006","sake007","wild004");
my @wine = ("beer014","beer020","beer024","beer030","beer033","beer088","sake002","spirits002","spirits004","spirits011","wine001","wine003","wine004","wine005","wine006","wine007","wine009","wine010","wine011","wine013","wine014","wine015","wine017","wine018");
my @beer1 = ("beer001","beer007","beer008","beer009","beer010","beer012","beer015","beer016","beer022","beer026","beer031","beer036","beer037","beer041","beer043","beer044","beer045","beer046","beer047","beer048","beer049","beer050","beer051","beer052","beer053","beer054","beer055","beer056","beer064","beer065","beer066","beer067","beer068","beer069","beer070","beer071","beer073","beer075","beer076","beer077","beer078","beer079","beer081","beer082","beer087","beer089","beer090","beer094","beer095","beer096","beer097","beer098","beer099","beer100","beer101","beer102","spirits005","wine012");
my @beer2 = ("beer002","beer003","beer004","beer011","beer013","beer021","beer027","beer032","beer034","beer039","beer040","beer059","beer062","beer063","beer080","beer083","beer084","beer085","beer086","beer091","beer092");
my @mixed = ("beer005","beer006","beer023","beer025","beer028","beer029","beer038","beer061","bread001","bread002","bread003","bread004","spirits001","spirits003","wild005","wild006","wild007");
my @beer1_beer2 = ("beer001","beer007","beer008","beer009","beer010","beer012","beer015","beer016","beer022","beer026","beer031","beer036","beer037","beer041","beer043","beer044","beer045","beer046","beer047","beer048","beer049","beer050","beer051","beer052","beer053","beer054","beer055","beer056","beer064","beer065","beer066","beer067","beer068","beer069","beer070","beer071","beer073","beer075","beer076","beer077","beer078","beer079","beer081","beer082","beer087","beer089","beer090","beer094","beer095","beer096","beer097","beer098","beer099","beer100","beer101","beer102","spirits005","wine012","beer002","beer003","beer004","beer011","beer013","beer021","beer027","beer032","beer034","beer039","beer040","beer059","beer062","beer063","beer080","beer083","beer084","beer085","beer086","beer091","beer092");

my @beer1_beer2_mixed = ("beer001","beer007","beer008","beer009","beer010","beer012","beer015","beer016","beer022","beer026","beer031","beer036","beer037","beer041","beer043","beer044","beer045","beer046","beer047","beer048","beer049","beer050","beer051","beer052","beer053","beer054","beer055","beer056","beer064","beer065","beer066","beer067","beer068","beer069","beer070","beer071","beer073","beer075","beer076","beer077","beer078","beer079","beer081","beer082","beer087","beer089","beer090","beer094","beer095","beer096","beer097","beer098","beer099","beer100","beer101","beer102","spirits005","wine012","beer002","beer003","beer004","beer011","beer013","beer021","beer027","beer032","beer034","beer039","beer040","beer059","beer062","beer063","beer080","beer083","beer084","beer085","beer086","beer091","beer092","beer005","beer006","beer023","beer025","beer028","beer029","beer038","beer061","bread001","bread002","bread003","bread004","spirits001","spirits003","wild005","wild006","wild007");

my @wine_beer1_beer2_mixed = ("beer014","beer020","beer024","beer030","beer033","beer088","sake002","spirits002","spirits004","spirits011","wine001","wine003","wine004","wine005","wine006","wine007","wine009","wine010","wine011","wine013","wine014","wine015","wine017","wine018","beer001","beer007","beer008","beer009","beer010","beer012","beer015","beer016","beer022","beer026","beer031","beer036","beer037","beer041","beer043","beer044","beer045","beer046","beer047","beer048","beer049","beer050","beer051","beer052","beer053","beer054","beer055","beer056","beer064","beer065","beer066","beer067","beer068","beer069","beer070","beer071","beer073","beer075","beer076","beer077","beer078","beer079","beer081","beer082","beer087","beer089","beer090","beer094","beer095","beer096","beer097","beer098","beer099","beer100","beer101","beer102","spirits005","wine012","beer002","beer003","beer004","beer011","beer013","beer021","beer027","beer032","beer034","beer039","beer040","beer059","beer062","beer063","beer080","beer083","beer084","beer085","beer086","beer091","beer092","beer005","beer006","beer023","beer025","beer028","beer029","beer038","beer061","bread001","bread002","bread003","bread004","spirits001","spirits003","wild005","wild006","wild007");

open my $in_fh , "<" ,$file;

my %dis;
while (<$in_fh>) {
    chomp;
    my @content = split /,/, $_;
    $dis{$content[0]} = $content[1];
}

#print YAML::Syck::Dump(\%dis);

my $asian_dis;
foreach (@asian){
    $asian_dis += $dis{$_};
}
my $asian_dis_mean = $asian_dis/scalar(@asian);

my $wine_dis;
foreach (@wine){
    $wine_dis += $dis{$_}; 
}
my $wine_dis_mean = $wine_dis/scalar(@wine);

my $beer1_dis;
foreach (@beer1){
    $beer1_dis += $dis{$_}; 
}
my $beer1_dis_mean = $beer1_dis/scalar(@beer1);

my $beer2_dis;
foreach (@beer2){
    $beer2_dis += $dis{$_}; 
}
my $beer2_dis_mean = $beer2_dis/scalar(@beer2);

my $mixed_dis;
foreach (@mixed){
    $mixed_dis += $dis{$_}; 
}
my $mixed_dis_mean = $mixed_dis/scalar(@mixed);

my $beer1_beer2_dis;
foreach (@beer1_beer2){
    $beer1_beer2_dis += $dis{$_}; 
}
my $beer1_beer2_dis_mean = $beer1_beer2_dis/scalar(@beer1_beer2);

my $beer1_beer2_mixed_dis;
foreach (@beer1_beer2_mixed){
    $beer1_beer2_mixed_dis += $dis{$_}; 
}
my $beer1_beer2_mixed_dis_mean = $beer1_beer2_mixed_dis/scalar(@beer1_beer2_mixed);

my $wine_beer1_beer2_mixed_dis;
foreach (@wine_beer1_beer2_mixed){
    $wine_beer1_beer2_mixed_dis += $dis{$_}; 
}
my $wine_beer1_beer2_mixed_dis_mean = $wine_beer1_beer2_mixed_dis/scalar(@wine_beer1_beer2_mixed);

#print "$asian_dis","\n";
#print "$wine_dis","\n";
#print "$beer1_dis","\n";
#print "$beer2_dis","\n";
#print "$mixed_dis","\n";


print "$gene_name".','."$asian_dis_mean".','."$wine_dis_mean".','."$beer1_dis_mean".','."$beer2_dis_mean".','."$mixed_dis_mean".','."$beer1_beer2_dis_mean".','."$beer1_beer2_mixed_dis_mean".','."$wine_beer1_beer2_mixed_dis_mean","\n";


__END__
