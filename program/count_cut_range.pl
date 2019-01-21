#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use File::Find::Rule;
use AlignDB::IntSpan;
use YAML::Syck qw();

use Getopt::Long qw();
use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 SYNOPSIS

    perl ~/Scripts/pars/program/count_cut_range.pl --file ~/data/mrna-structure/process/Scer_n128_Spar.gene_variation.process.yml --cut ~/data/mrna-structure/process/sce_cds.yml --output stem_loop_cds_length.csv 
    perl ~/Scripts/pars/program/count_cut_range.pl --file ~/data/mrna-structure/process/Scer_n128_Spar.gene_variation.process.yml --cut ~/data/mrna-structure/process/sce_utr.yml --output stem_loop_utr_length.csv

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'cut|c=s'   => \my $cut,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

  
open OUT, ">",$output;
print OUT "gene" . ",". "chr" . ",". "stem_length" .",". "loop_length"."\n";

print "Load $file\n";
my $gene_info_of = YAML::Syck::LoadFile($file);

foreach my $gene ( sort keys %{$gene_info_of} ){
	print $gene."\n";
	my $info = $gene_info_of->{ $gene };
  my $stem = AlignDB::IntSpan->new;
  $stem->add( $info->{fold_left} );
  $stem->add( $info->{fold_right} );
  #print $stem."\n"."\n";sleep 1;
  my $loop = AlignDB::IntSpan -> new($info->{fold_dot});
  #print $loop."\n";sleep 1;
  
  my $chr = $info->{chr};
  my $chr_set = AlignDB::IntSpan->new($info->{chr_set});
  my @chr_set = $chr_set->ranges();
  #print $chr_set[0]."\n";

  #transfer relative position to absolute position
  my $stem2 = $stem->map_set(sub{$_+$chr_set[0]-1});
  #print $stem2."\n"."\n";sleep 1;
  my $loop2 = $loop->map_set(sub{$_+$chr_set[0]-1});
  #print $loop2."\n"."\n";sleep 1;
  
  print "Load $cut\n";
  my $gene_info2_of = YAML::Syck::LoadFile($cut);
  
  my $info2 = $gene_info2_of->{ $chr };

  
  my $cut_range = AlignDB::IntSpan->new;
  $cut_range->add( $info2 );
  #print $cut_range."\n";sleep 1;
  
  my $stem_intersection = AlignDB::IntSpan::intersect( $stem2, $cut_range );
  print $stem_intersection."\n";
  my $stem_intersection_length = $stem_intersection->cardinality();           
  my $loop_intersection = AlignDB::IntSpan::intersect( $loop2, $cut_range );
  print $loop_intersection."\n";
  my $loop_intersection_length = $loop_intersection->cardinality();             
                
  open OUT,">>",$output;
  print OUT "$gene" . ",". "$chr" . ",". "$stem_intersection_length" . ",". "$loop_intersection_length" . "\n";     
}	
	


close OUT;
__END__
