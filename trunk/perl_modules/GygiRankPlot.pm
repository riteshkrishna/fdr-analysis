package GygiRankPlot;
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
our @EXPORT = qw(GetGygiRankPlot);

use strict;
use URI::Escape;

use lib qw(~/LIVERPOOL/bin/perl_modules/);
use strict;
use CreatePNG;


sub GetGygiRankPlot
{
my $se = shift;
my $imagefile = shift;
my $cutoff = shift;
my @results = @_;

my $engine;
 if($se eq "M")
 {
 $engine = "Mascot";
 }
 elsif($se eq "O")
 {
 $engine = "Omssa";
 }
 elsif($se eq "T")
 {
 $engine = "XTandem";
 } 

my @graph;

my $max_rank = 0;
#only include reverse if rank 1 is > $cutoff expect


 for(my $r=1 ; $r<scalar(@results) ; $r++)
 {
  #sometimes mascot doesn't cover all sequence

  if($results[$r][1]{'sequence'} && $results[$r][1]{'protein'} && $results[$r][1]{'expect'}<$cutoff)
  {
   for(my $rank=1 ; $rank<(scalar(@{$results[$r]})) ; $rank++)
   {

    if($rank>$max_rank)
    {
    $max_rank = $rank;
    }

    if($results[$r][$rank]{'protein'} =~ m/^REV\_/ || $results[$r][$rank]{'protein'} =~ m/\_r$/)
    {  

     if($graph[1][$rank])
     {
     $graph[1][$rank]++;
     }
     else
     {
     $graph[1][$rank] = 1;
     }
    }
    elsif($results[$r][$rank]{'protein'})
    {
     if($graph[0][$rank])
     {
     $graph[0][$rank]++;
     }
     else
     {
     $graph[0][$rank] = 1;
     }
    }
   }
  }
 }

#for each rank calculate the total forward and reverse
my @plotdata;
my $min = 0;
my $max = 10;

$max_rank++;

 for(my $rank=1 ; $rank<$max_rank ; $rank++)
 {
 my $for;
 my $rev;
  if(!$graph[0][$rank] && !$graph[1][$rank])
  {
  $for = 0;
  $rev = 0;
  }
  elsif(!$graph[0][$rank] && $graph[1][$rank])
  {
  $for = 0;
  $rev = 100; 
  }
  elsif($graph[0][$rank] && !$graph[1][$rank])
  {
  $for = 100;
  $rev = 0;
  }
  else
  {
  my $total = $graph[0][$rank] + $graph[1][$rank];
  $for = $graph[0][$rank]/$total*100;
  $rev = $graph[1][$rank]/$total*100;
  }


 $plotdata[0][$rank-1] = $rank;
 $plotdata[1][$rank-1] = $for;
 $plotdata[2][$rank-1] = $rev;
 }

#and now the png
my $image;

#my $title = $engine . " Gygi Rank Plot";
my $title = "$engine";

$image = GetPNGBars($title,'rank','percent',$min,$max,@plotdata);

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


