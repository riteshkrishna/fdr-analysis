package DeltaMassPlot;
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
our @EXPORT = qw(GetDeltaMassPlot);

use strict;
use URI::Escape;

use lib qw(~/LIVERPOOL/bin/perl_modules/);
use strict;
use MascotParser;
use CreatePNG;

sub GetDeltaMassPlot
{
my $imagefile = shift;
my $scoretype = shift;
my $setype = shift;
my @results = @_;

my @graph;

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


my %delta;
my %delta1;
my $max_tmp=0;
 for(my $r=1 ; $r<scalar(@results) ; $r++)
 {
  #sometimes mascot doesn't cover all sequence

  if($results[$r][1]{'sequence'} && $results[$r][1]{'sequence'} ne "NULL")
  {
   #only rank 1
   for(my $rank=1 ; $rank<2 ; $rank++)
   {

   #get the discriminant score
   my $identity = GetIdentity($results[$r][$rank]{'qmatch'});
   #get the identity score
   my $score = $results[$r][$rank]{'ionscore'};
   my $disc;

 
    #discriminant score required?
    if($scoretype eq "D")
    {
    $disc = $score - $identity;
    }
    else
    {
    $disc = $score;
    }

   if($setype ne "O")
   {
   $disc = int($disc);
   }
   else
   {
   $disc = sprintf("%.5f",$disc);
   }



    #make the delat to 2 dp
    $results[$r][$rank]{'delta'} = sprintf("%.2f",$results[$r][$rank]{'delta'});


    #for the reverses
    if($results[$r][$rank]{'protein'} =~ m/^REV\_/ || $results[$r][$rank]{'protein'} =~ m/\_r$/)
    {  
    #reverse hits
    $delta{$results[$r][$rank]{'delta'}} = $disc;
    }
    elsif($results[$r][$rank]{'protein'})
    {
    #forward hits
    $delta1{$results[$r][$rank]{'delta'}} = $disc;     
    }
   }
  }
 }

 #forwards

 foreach my $d (keys %delta1)
 {
 #print "$d,$delta1{$d}, \n";
 }
 #reverses
 foreach my $d (keys %delta)
 {
 #print "$d, ,$delta{$d}\n";
 }



my $image;
my %plotdata;
my $count = 0;
my $min = 1;
my $max = 0;
my @all;

#get the min and max score
 foreach my $del (keys %delta1)
 {
  if($del>$max)
  {
  $max = $del;
  } 
  if($del<$min)
  {
  $min = $del;
  }
 }
 
 foreach my $del (keys %delta)
 {
  if($del>$max)
  {
  $max = $del;
  }
  if($del<$min)
  {
  $min = $del;
  }
 }

#push the deltas into mass bins of 0.1 and average the score
 for(my $s=$min ; $s<$max ; $s=$s+0.025)
 {

 my $total_score = 0;
 my $count = 0; 
  foreach my $del (keys %delta1) 
  {
   if($del>=$s && $del<$s+0.025)
   {
   $total_score += $delta1{$del};
   $count++;
   }
  }

 if($total_score ==0 || $count ==0)
 {
 $total_score = 0;
 }
 else
 {
 $total_score = $total_score/$count;
 }

 $plotdata{$s}[1] = $total_score;
  $total_score = 0;
  $count = 0;

  foreach my $del (keys %delta)
  {
   if($del>=$s && $del<$s+0.1)
   {
   $total_score = $delta{$del};
   $count++;
   }
  }

 if($total_score ==0 || $count ==0)
 {
 $total_score = 0;
 }
 else
 {
 $total_score = $total_score/$count;
 }
  push(@all,$s);
  $plotdata{$s}[2] = $total_score;
  } 

# #for the forwards
# foreach my $del (keys %delta1)
# {
#  if($del>$max)
#  {
#  $max = $del;
#  }
#push(@all,$del);
#
# $plotdata{$del}[1] = $delta1{$del};
# $count++;
# }
#
# #and the reverse
# foreach my $del (keys %delta)
# {
#push(@all,$del);
#  if($del>$max)
#  {
#  $max = $del;
#  }

# $plotdata{$del}[2] = $delta{$del};
# $count++;
# }


#i now want to sort the plot data
#my @sorted = sort { $a <=> $b } @{$plotdata[0]};
my @sorted = sort {$a <=> $b} (@all);
my @sorted_plot_data;
my $count = 0;

 for(my $s=0 ; $s<scalar(@sorted) ; $s++)
 {
 $sorted_plot_data[0][$count] = $sorted[$s];
 $sorted_plot_data[1][$count] = $plotdata{$sorted[$s]}[1];
 $sorted_plot_data[2][$count] = $plotdata{$sorted[$s]}[2];
 $count++;
 }

 $scoretype = "score";
# my $title = $engine . " Delta Mass Plot";
 my $title = "$engine";
 $image = GetPNGPoints($title,'Delta Mass',$scoretype,$min,$max,@sorted_plot_data);
 
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
