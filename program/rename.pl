#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use FindBin;
use YAML::Syck qw();

use Path::Tiny;
use Bio::SearchIO;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

$^I =".bak";
my $n = 0;
while (<>){
	if(m/^>/){
	   $n = $n + 1;
		 s/^>.*/>seq$n/;
	   print;
	}else{
	   print;
	}	
}
