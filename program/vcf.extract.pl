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
use App::Rangeops;

use POSIX;
use Data::Dumper "Dumper";

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS
    
    cd /Volumes/Backup/yumh/data/vcf/1011Matrix.gvcf/
    perl extract.pl --file chr1.tsv --output chr1.ext.tsv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

my $tsv_fh;
open $tsv_fh, '<', $file;

open OUT, '>', $output;

while (<$tsv_fh>) {
    chomp;
    if (m/^\#/) {
        my @content;
        $content[0] = "name";
        $content[1] = "REF";
        $content[2] = "ALT";
        $content[3] = "QUAL";
        $content[4] = "AC";
        $content[5] = "AF";
        $content[6] = "AN";
        my $content = join "\t", @content;
        open OUT,'>>',$output;;
        print OUT $content, "\n";
    }else{
        my @content = split /\t/, $_;
        my @content_new;
        my $chr;
        if ($content[0] =~ /^chromosome(\d+)/){
        		$chr = &trans($1);
        }
        
        $content_new[0] = "$chr:$content[1]";
        $content_new[1] = $content[3];
        $content_new[2] = $content[4];
        $content_new[3] = $content[5];
        my @info = split /;/, $content[7];
        my @AC = split /=/, $info[0];
        $content_new[4] = $AC[1];
        my @AF = split /=/, $info[1];
        $content_new[5] = $AF[1];
        my @AN = split /=/, $info[2];
        $content_new[6] = $AN[1];
        my $content = join "\t", @content_new;
        open OUT,'>>',$output;;
        print OUT $content, "\n";
    }
}

sub trans {
		my %roman = (16=>"XVI",15=>"XV",14=>"XIV",13=>"XIII",12=>"XII",11=>"XI",10=>"X",9=>"IX",8=>"VIII",7=>"VII",6=>"VI",5=>"V",4=>"IV",3=>"III",2=>"II",1=>"I");
		$roman{"$_[0]"};
}

close $tsv_fh;
close OUT;

__END__



