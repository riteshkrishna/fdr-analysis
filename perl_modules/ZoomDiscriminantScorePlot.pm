package ZoomDiscriminantScorePlot;
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
our @EXPORT = qw(ZoomScoreDistribution);

use strict;
use URI::Escape;

use lib qw(~/LIVERPOOL/bin/perl_modules/);
use strict;
use MascotParser;
use CreatePNG;


sub ZoomScoreDistribution
{
my $imagefile = shift;
my $scoretype = shift;
my $setype = shift;
my @results = @_;

my @graph;
my %tmpdistribution;
my %distribution;
my $min=1;
my $max = 0;

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

  if($results[$r][1]{'sequence'})
  {
   #only rank 1
   for(my $rank=1 ; $rank<2 ; $rank++)
   {

   #get the discriminant score
   my $identity = GetIdentity($results[$r][$rank]{'qmatch'});
   #get the identity score
   my $score = $results[$r][$rank]{'ionscore'};
   my $disc;

    if($scoretype eq "D")
    {
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

   

    if($distribution{$disc}[$rev] && $distribution{$disc}[$rev]<100)
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

#put the distribution into bins
foreach my $score (keys %distribution)
{
my $newscore = $score + (-$score % 5);
 if($tmpdistribution{$newscore}[0])
 {
 $tmpdistribution{$newscore}[0] += $distribution{$score}[0];
 }
 else
 {
 $tmpdistribution{$newscore}[0] = $distribution{$score}[0];
 }
 
 if($tmpdistribution{$newscore}[1])
 {
 $tmpdistribution{$newscore}[1] += $distribution{$score}[1];
 }
 else
 {
 $tmpdistribution{$newscore}[1] = $distribution{$score}[1];
 }
}

my %plotdata;
my $count = 0;
my @scorearray;

 foreach my $score (keys %tmpdistribution)
 {
  #print "$score,";
  $plotdata{$score}[0] = $score;
push(@scorearray,$score);
   if($tmpdistribution{$score}[0] && $tmpdistribution{$score}[1])
   {
   #this will be the estimated correct
   $plotdata{$score}[1] = ($tmpdistribution{$score}[0]-$tmpdistribution{$score}[1]);
   }
   elsif($tmpdistribution{$score}[0])
   {
   $plotdata{$score}[1] = $tmpdistribution{$score}[0];
   }
   if($tmpdistribution{$score}[1])
   {
   #and the estimated incorrect
   $plotdata{$score}[2] = $tmpdistribution{$score}[1];
   }
  $count++;
 }


#my @sorted = sort { $a <=> $b } @{$plotdata[0]};
my @sorted = sort {$a <=> $b} (@scorearray);
my @sorted_plot_data;
my $count = 0;

 for(my $s=0 ; $s<scalar(@sorted) ; $s++)
 {
 $sorted_plot_data[0][$count] = $sorted[$s];
 $sorted_plot_data[1][$count] = $plotdata{$sorted[$s]}[1];
 $sorted_plot_data[2][$count] = $plotdata{$sorted[$s]}[2];
 $count++;
 }

 $scoretype .= " score";
 my $title = $engine . " Estimation of Correct/Incorrect";
 my $image = GetPNGLines($title,$scoretype,'number spectra',$min,$max,'EC',@sorted_plot_data);

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
