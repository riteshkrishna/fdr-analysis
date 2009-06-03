package GetNterminaldistribution;
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
our @EXPORT = qw(GetNtermPlot);

use strict;
use URI::Escape;

use lib qw(~/LIVERPOOL/bin/perl_modules/);
use strict;
use MascotParser;
use CreatePNG;

sub GetNtermPlot
{
my $se = shift;
my $imagefile = shift;
my $REVTAG = shift;
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

my %nterm;


 for(my $r=1 ; $r<scalar(@results) ; $r++)
 {
  #sometimes mascot doesn't cover all sequence

  if($results[$r][1]{'sequence'} && $results[$r][1]{'sequence'} ne "NULL")
  {
   #only rank 1
   for(my $rank=1 ; $rank<2 ; $rank++)
   {

   my $rev = 0;
    if($results[$r][$rank]{'protein'} =~ m/$REVTAG/)
    {
    $rev = 1;
    }   

    my $start = $results[$r][$rank]{'start'};
    if($results[$r][$rank]{'expect'}<0.05)
    {
     if($nterm{$start}[$rev])
     {
     $nterm{$start}[$rev]++;
     }
     else
     {
     $nterm{$start}[$rev] = 1;
     } 
    }
   }
  }
 }

my @plotdata;
my $count = 0;
my $min = 1;
my $max = 1;

my @total;
$total[0] = 0;
$total[1] = 0;

 foreach my $position (keys %nterm)
 {
#  if($nterm{$position}[0])
#  {
  $total[0] += $nterm{$position}[0];
#  }
#  if($nterm{$position}[1])
#  {
  $total[1] += $nterm{$position}[1];
#  }
 }
 
 foreach my $position (keys %nterm)
 {
  if($position>$max)
  {
  $max = $position;
  }
  $plotdata[0][$position] = $position;
#   if($nterm{$position}[0])
#   {
   $plotdata[1][$position] = ($nterm{$position}[0]/$total[0])*100;
#   }
   if($nterm{$position}[1] && $total[1])
   {
   $plotdata[2][$position] = ($nterm{$position}[1]/$total[1])*100; 
   }
  $count++;
 }

my @sorted = sort { $a <=> $b } @{$plotdata[0]};

my @sorted_plot_data;
my $count = 0;
 for(my $s=0 ; $s<scalar(@sorted) ; $s++)
 {
  if($sorted[$s])
  {
  $sorted_plot_data[0][$count] = $sorted[$s];
  $sorted_plot_data[1][$count] = $nterm{$sorted[$s]}[0];
  $sorted_plot_data[2][$count] = $nterm{$sorted[$s]}[1];
#$sorted_plot_data[1][$count] = $plotdata[1]{$sorted[$s]};
#$sorted_plot_data[2][$count] = $plotdata[2]{$sorted[$s]};
  $count++;
  }
 }
 my $title = $engine . " Nterminal Distribution";
print "in Nter score plot about to call bars\n";
my $image = GetPNGBars($title,'sequence position','number spectra',$min,$max,@sorted_plot_data);

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
