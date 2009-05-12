package GetFDRValues;
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
our @EXPORT = qw(GetFDRValues);

use URI::Escape;

use lib qw(~/LIVERPOOL/bin/perl_modules/);
use MascotParser;

sub GetFDRValues
{
my $max_fdr = shift;
my $fdrtype = shift;
my $setype = shift;
my $nterminal = shift;
my $rev_tag = shift;
my @results = @_;

my @graph;
my %fdr_gygi;
my %fdr_jones;
my $min = 1;
my $max = 0;

my %score_count;
my %gygi_obs;
my %jones_obs;
my %nter_count;
my @scorerecord;
my $identity_count = 0;
my $identity;
my %unique_peptides;

my $expect_count = 0;
my %expect_peptides;

my $minimum_score = 10000;


 #for all the peptides
 for(my $s=0 ; $s<scalar(@results) ; $s++)
 {
 my $r = 1;
  if($results[$s][$r]{'sequence'} && $results[$s][$r]{'sequence'} ne "NULL")
  {
  my $disc = $results[$s][$r]{'expect'};
   #get the sum of for and negative hits above this score

   if(!$fdr_gygi{$disc} && $disc<$minimum_score)
#if(!$fdr_gygi{$disc})
   {

   push(@scorerecord,$disc);
   my @sum = GetSums($disc,'D',$rev_tag,@results);
   #forward is sum[0] and reverse sum[1]
   #Gygi
   $fdr = 0; 
   $fdr = ($sum[1])+($sum[0]);
   $fdr = (2*$sum[1])/$fdr;
   $fdr_gygi{$disc} = $fdr;

   #Jones
   $fdr = 0;
   $fdr = $sum[0];
   $fdr = $sum[1]/$fdr;
   $fdr_jones{$disc} = $fdr;


   #which FDR is lower?
   if($fdr>$tmp_gygi_fdr)
   {
    if($fdr>=0.6 && $disc<$minimum_score)
    {
    $minimum_score = $disc;
    }
   }
   else
   {
    if($tmp_gygi_fdr>=0.6 && $disc<$minimum_score)
    {
    $minimum_score = $disc;
    }
   }
  } 


   if($score_count{$disc} && !$gygi_obs{$results[$s][$r]{'sequence'}})
   {
   $score_count{$disc}++;
   $gygi_obs{$results[$s][$r]{'sequence'}} = $disc;
    if($results[$s][$r]{'start'}<3 && $nter_count{$disc})
    {
    $nter_count{$disc}++;
    }
    elsif($results[$s][$r]{'start'}<3)
    {
    $nter_count{$disc} = 1;
    }
   }
   #otherwise we have seen the sequence before but with sa worse score
   elsif($score_count{$disc} && $disc<$gygi_obs{$results[$s][$r]{'sequence'}})
   {

   $score_count{$gygi_obs{$results[$s][$r]{'sequence'}}}--;
   $score_count{$disc}++;
    if($results[$s][$r]{'start'}<3 && $nter_count{$disc})
    {

    $nter_count{$gygi_obs{$results[$s][$r]{'sequence'}}}--;
    $nter_count{$disc}++;
    }
    elsif($results[$s][$r]{'start'}<3)
    {

    $nter_count{$gygi_obs{$results[$s][$r]{'sequence'}}}--;
    $nter_count{$disc} = 1;
    }
    $gygi_obs{$results[$s][$r]{'sequence'}} = $disc;
   }
   elsif(!$gygi_obs{$results[$s][$r]{'sequence'}})
   {
   $score_count{$disc} = 1;
   $gygi_obs{$results[$s][$r]{'sequence'}} = $disc;
    if($results[$s][$r]{'start'}<3 && $nter_count{$disc})
    {
    $nter_count{$disc}++;
    }
    elsif($results[$s][$r]{'start'}<3)
    {
    $nter_count{$disc} = 1;
    }
   }
   elsif($disc<$gygi_obs{$results[$s][$r]{'sequence'}})
   {
   $score_count{$gygi_obs{$results[$s][$r]{'sequence'}}}--;
   $score_count{$disc} = 1;
    if($results[$s][$r]{'start'}<3 && $nter_count{$disc})
    {
    $nter_count{$gygi_obs{$results[$s][$r]{'sequence'}}}--;
    $nter_count{$disc}++;
    }
    elsif($results[$s][$r]{'start'}<3)
    {
    $nter_count{$gygi_obs{$results[$s][$r]{'sequence'}}}--;
    $nter_count{$disc} = 1;
    }
    $gygi_obs{$results[$s][$r]{'sequence'}} = $disc;
    }
  }
 }


my %counter;
my %ntercounter;
my $fdr;
my $max_05 = 0;

 $counter{$max_fdr} = 0;
 $ntercounter{$max_fdr} = 0;
  foreach my $score (keys %fdr_gygi)
  {
   #if the fdr at this score is better than threshold
   if($fdr_gygi{$score} <= $max_fdr)
   { 
   $counter{$max_fdr} += $score_count{$score};
   $ntercounter{$max_fdr} += $nter_count{$score}; 

    if($fdr_gygi{$score}<=$max_fdr && $score>$max_05)
    {
    $max_05 = $score;
    }
   }
  }
 

#get the peptides
my %pep;
 for(my $s=0 ; $s<scalar(@results) ; $s++)
 {
 my $r = 1;
  if($results[$s][$r]{'sequence'} && $results[$s][$r]{'sequence'} ne "NULL" && $results[$s][$r]{'expect'}<=$max_05)
  {
  $pep{$max_fdr} .= "$results[$s][$r]{'protein'}\#\#\#$results[$s][$r]{'sequence'}\#\#\#$results[$s][$r]{'start'}\#\#\#$results[$s][$r]{'ionscore'}\#\#\#$results[$s][$r]{'expect'}\*\*\*";
  }
 } 

 my %res;
 #for all the fdr scores
 foreach my $fdr (keys %counter)
 {
 $res{$fdr}[0] = $counter{$fdr};
 $res{$fdr}[1] = $ntercounter{$fdr};
 $res{$fdr}[4] = $pep{$fdr};
 }


#repeat for Jones
%counter = ();
%ntercounter = ();
$fdr = "";
my $max_05 = 0;
 $counter{$max_fdr} = 0;
 $ntercounter{$max_fdr} = 0;
  foreach my $score (keys %fdr_jones)
  {
   #if the fdr at this score is better than threshold
   if($fdr_jones{$score} <= $max_fdr)
   {
   $counter{$max_fdr} += $score_count{$score};
   $ntercounter{$max_fdr} += $nter_count{$score};
   }
   if($fdr_jones{$score}<=$max_fdr && $score>$max_05)
   {
   $max_05 = $score;
   }
  }


#my %pep;
 for(my $s=0 ; $s<scalar(@results) ; $s++)
 {
 my $r = 1;
  if($results[$s][$r]{'sequence'} && $results[$s][$r]{'sequence'} ne "NULL" && $results[$s][$r]{'expect'}<=$max_05)
  {
  $pep{$max_fdr} .= "$results[$s][$r]{'protein'}\#\#\#$results[$s][$r]{'sequence'}\#\#\#$results[$s][$r]{'start'}\#\#\#$results[$s][$r]{'ionscore'}\#\#\#$results[$s][$r]{'expect'}\*\*\*";
  }
 }


 foreach my $fdr (keys %counter)
 {
 $res{$fdr}[2] = $counter{$fdr};
 $res{$fdr}[3] = $ntercounter{$fdr};
 $res{$fdr}[5] = $pep{$fdr};
 }
 


return %res;
} 


sub InitColors {
    my($im) = $_[0];
    # Allocate colors
    $white = $im->colorAllocate(255,255,255);
    $black = $im->colorAllocate(0,0,0);
    $red = $im->colorAllocate(255,0,0);
    $blue = $im->colorAllocate(0,0,255);
    $green = $im->colorAllocate(0, 255, 0);

    $brown = $im->colorAllocate(255, 0x99, 0);
    $violet = $im->colorAllocate(255, 0, 255);
    $yellow = $im->colorAllocate(255, 255, 0);
}

sub CountSpectra
{
my $threshold = shift;
my %counts = %{$_[0]};
my $result = 0;

 #for all the results
 foreach my $score (keys %counts)
 {
  if($score >= $threshold)
  {
  $result += $counts{$score};
  }
 }

return $result;
}


return 1;


sub GetSums
{
my $threshold = shift;
my $scoretype = shift;
my $reverse_tag = shift;
my @results = @_;
my @sum;

#forward
$sum[0] = 0;
#reverse
$sum[1] = 0;


 for(my $s=0 ; $s<scalar(@results) ; $s++)
 {
 my $r=1;
  if($results[$s][$r]{'sequence'} && $results[$s][$r]{'sequence'} ne "NULL")
  {
   my $disc = $results[$s][$r]{'expect'};
   if($results[$s][$r]{'protein'} !~ m/^$reverse_tag/)
   {

    if($disc <= $threshold)
    {
    $sum[0]++;
    }
   }
   else
   {
    if($disc <= $threshold)
    { 
    $sum[1]++;
    }
   }
  }
 }


return @sum;
}



 
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
