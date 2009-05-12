package GetDiscriminantScorePlot;
require Exporter;

#####################################################################################################
#        Copyright (C) 2009, Jennifer Siepen, University of Manchester                              #
#                                                                                                   #
#    This program is free software: you can redistribute it and/or modify                           #
#    it under the terms of the GNU General Public License as published by                           #
#    the Free Software Foundation, either version 3 of the License, or                              #
#    (at your option) any later version.                                                            #
#                                                                                                   #
#    This program is distributed in the hope that it will be useful,                                #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                                 #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                  #
#    GNU General Public License for more details.                                                   #
#                                                                                                   #
#    You should have received a copy of the GNU General Public License                              #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.                          #
#                                                                                                   #
#                                                                                                   #
#####################################################################################################


our @ISA = qw(Exporter);
our $VERSION = 1.0;;
our @EXPORT = qw(GetScoreDistribution);

use strict;
use URI::Escape;

use lib qw(~/LIVERPOOL/bin/perl_modules/);
use strict;
use MascotParser;
use CreatePNG;


sub GetScoreDistribution
{
my $imagefile = shift;
my $scoretype = shift;
my $setype = shift;
my @results = @_;

my @graph;
my %distribution;
my $min=1;
my $max = 0;
my @scoresorter;

my $engine;
 if($setype eq "M")
 {
 $engine = "Mascot";
 }
 elsif($setype eq "O")
 {
 $engine = "Omssa";
 }
 elsif($setype eq "T")
 {
 $engine = "XTandem";
 }


 for(my $r=1 ; $r<scalar(@results) ; $r++)
 {
  #sometimes mascot doesn't cover all sequence

  if($results[$r][1]{'sequence'} && $results[$r][1]{'sequence'} ne "NULL")
  {
   #only rank 1
   for(my $rank=1 ; $rank<2 ; $rank++)
   {

   #get the identity score
   my $score = $results[$r][$rank]{'ionscore'};
   my $disc;
    #discriminant score required?
    if($scoretype eq "D")
    {
    my $identity = GetIdentity($results[$r][$rank]{'qmatch'});
    $disc = $score - $identity;
    }
    else
    {
    $disc = $score;
    }

   $disc = int($disc);

    if($disc>$max)
    {
    $max = $disc;
    }

   my $rev = 0;
    if($results[$r][$rank]{'protein'} =~ m/REV\_/)
    {
    $rev = 1;
    }   

   

    if($distribution{$disc}[$rev] && $distribution{$disc}[$rev]<1000)
    {
    $distribution{$disc}[$rev]++;
    }
    elsif(!$distribution{$disc}[$rev])
    {
    $distribution{$disc}[$rev] = 1;
    } 

   }
  }
 }

my @plotdata;
my $count = 0;
 foreach my $score (keys %distribution)
 {
push(@scoresorter,$score);
  $plotdata[0][$count] = $score;
   if($distribution{$score}[0])
   {
   $plotdata[1][$count] = $distribution{$score}[0];
   }
   if($distribution{$score}[1])
   {
   $plotdata[2][$count] = $distribution{$score}[1];

   }
  #print "\n"; 
  $count++;
 }


#my @sorted = sort { $a <=> $b } @{$plotdata[0]};
my @sorted = sort {$a <=> $b} (@scoresorter);

my @sorted_plot_data;
my $count = 0;
 for(my $s=0 ; $s<scalar(@sorted) ; $s++)
 {
  if($sorted[$s-1] && $sorted[$s-1] != $sorted[$s])
  {
 $sorted_plot_data[0][$count] = $sorted[$s];
 $sorted_plot_data[1][$count] = $distribution{$sorted[$s]}[0];
 $sorted_plot_data[2][$count] = $distribution{$sorted[$s]}[1];;
 $count++;
 }
}
 $scoretype .= " score";
 my $title = $engine . " Score Distribution";
 my $image = GetPNGBars($title,$scoretype,'number spectra',$min,$max,@sorted_plot_data);


 my $mainimage = new GD::Image((450), (300));
 my $white = $mainimage->colorAllocate(255,255,255);

 $mainimage->copy($image,0,0,0,0,450,300);
 open (CHART, ">$imagefile") or print "unable to open the $imagefile file ";
 print CHART $mainimage->png;
 close CHART;

return 1;

} 



return 1;


 
sub GetDbPath
{
my $file = shift;
my $db;
open(FILE,"<$file") or die "unable to open the mascot results, $file\n";
 while(my $line = <FILE>)
 {
  if($line =~ m /^fastafile\=/)
  {
  my @split = split/\=/,$line;
  $db = $split[1];
  close FILE;
  }  
 }
close FILE;

$db =~ s/\/usr\/local\/mascot\//\/fs\/msct\//;


return $db;

}

sub GetIdentity
{
my $qmatch = shift;
my $id = 0;

 if($qmatch<1)
 {
 $id = -10.0 * log(1.0/10.0)/log(10);
 }
 else
 {
 $id = -10 * log(1.0/(1.0*$qmatch))/log(10);
 }

return $id;

}
