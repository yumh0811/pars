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

    perl ~/Scripts/pars/program/PARS_score.process.pl --file ~/data/mrna-structure/PARS10/pubs/PARS10/data/sce_Score.tab --output ~/data/mrna-structure/process/${NAME}.gene_variation.process.update.yml

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

my $csv_fh;
open $csv_fh,'<',$file;

open OUT,'>',$output;

my %info;
while (<$csv_fh>) {
		chomp;
    my @info = split /\t/, $_;
    my @pars = split /;/, $info[2];
    
    my %gene_length;
  	$gene_length{"length"} = $info[1];
  	#print YAML::Syck::Dump(\%gene_length),"\n";sleep 2;
  	my @stem = ();
  	my @loop = ();
  	for( my $i=0; $i<=$#pars; $i++ ){
  			if ($pars[$i]>=0){
						push @stem, $i+1;
				}else{
						push @loop, $i+1;
				}
	  }
	  my $stem = join (",",@stem);
		my $loop = join (",",@loop);
	  my $stem_range = AlignDB::IntSpan->new("$stem")->as_string;
	  my $loop_range = AlignDB::IntSpan->new("$loop")->as_string;
	  my %structure = ( "stem" => $stem_range , "loop" => $loop_range );
	  
	  my %in = (%gene_length,%structure);
    $info{"$info[0]"}= \%in;

    #print YAML::Syck::Dump(\%info),"\n";sleep 2;
}
open OUT, ">",$output;
print OUT YAML::Syck::Dump(\%info),"\n";

close OUT;

__END__
