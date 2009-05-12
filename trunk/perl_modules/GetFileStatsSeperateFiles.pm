package GetFileStatsSeperateFiles;
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
our @EXPORT = qw(GetSeperateStats);

use URI::Escape;

use lib qw(~/LIVERPOOL/bin/perl_modules/);
use MascotParser;
use CreatePNG;
use GD;
use GD::Image;

sub GetSeperateStats
{
my $scoretype = shift;
my $fdrtype = shift;
my $setype = shift;
my @results = @_;


my @graph;
my %fdr_all;
my $min = 1;
my $max = 0;

my %score_count;
my %nter_count;
my @scorerecord;
my $identity_count = 0;
my $identity;
my %unique_peptides;

my @expect_count;
$expect_count[0] = 0;
$expect_count[1] = 0;
my %expect_peptides;
my %expect_proteins;
#calculate a global identity score
my $gidentity;
my $gidentity_count = 0;
my %g_unique_peptides;
my $count = 0;
 for(my $r=1 ; $r<scalar(@results) ; $r++)
 {
  for(my $rank=1 ; $rank<2 ; $rank++)
  {
   if($results[$r][$rank]{'qmatch'} && $results[$r][$rank]{'qmatch'}>0)
   {
   $gidentity += GetIdentity($results[$r][$rank]{'qmatch'});
   $count++;
   }
  }
 }
  if($count>0)
  {
  $gidentity = $gidentity/$count;
  }
  else
  {
  $gidentity = 0;
    }
 for(my $r=1 ; $r<scalar(@results) ; $r++)
 {
  #sometimes mascot doesn't cover all sequence

  if($results[$r][1]{'sequence'} && $results[$r][1]{'sequence'} ne "NULL")
  {
   #only rank 1
   for(my $rank=1 ; $rank<2 ; $rank++)
   {

   $identity = 0;
   #get the discriminant score

   $identity = GetIdentity($results[$r][$rank]{'qmatch'});
   #get the identity score
   my $score = $results[$r][$rank]{'ionscore'};
   my $disc;


    #individual identity scores
    if($identity>0 && $results[$r][$rank]{'ionscore'} >= $identity && !$unique_peptides{$results[$r][$rank]{'sequence'}} && $results[$r][$rank]{'protein'} !~ m/^REV\_/)
    {
    $identity_count++;
    $unique_peptides{$results[$r][$rank]{'sequence'}} = 1;
    }
    elsif($identity == 0 && $results[$r][$rank]{'expect'}<=0.05 && !$unique_peptides{$results[$r][$rank]{'sequence'}} && $results[$r][$rank]{'protein'} !~ m/^REV\_/)
    {
    $identity_count++;
    $unique_peptides{$results[$r][$rank]{'sequence'}} = 1;
    }
    #average identity
    if($results[$r][$rank]{'ionscore'} >= $gidentity && !$g_unique_peptides{$results[$r][$rank]{'sequence'}} && $results[$r][$rank]{'protein'} !~ m/^REV\_/)
    {
    $gidentity_count++;
    $g_unique_peptides{$results[$r][$rank]{'sequence'}} = 1;
    }
    #expect matches
    if($results[$r][$rank]{'expect'} <= 0.05 && !$expect_peptides{$results[$r][$rank]{'sequence'}} && $results[$r][$rank]{'protein'} !~ m/^REV\_/)
    {
     if($results[$r][$rank]{'start'}<3)
     {
     $expect_count[0]++;
     }
     else
     {
     $expect_count[1]++;
     }
    $expect_peptides{$results[$r][$rank]{'sequence'}} = 1;
    $expect_proteins{$results[$r][$rank]{'protein'}} = 1;
    }


    #discriminant score required?
    if($setype eq "M")
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

    if($disc>$max)
    {
    $max = $disc;
    }

   #get the sum of for and negative hits above this score
    if(!$fdr_all{$disc})
    {
    push(@scorerecord,$disc);
    my @sum = GetSums($disc,$scoretype,$setype,@results);

    #forward is sum[0] and reverse sum[1]
     if($fdrtype == '0')
     {
     #forward is sum[0] and reverse sum[1]
     $fdr = ($sum[1])+($sum[0]);
     $fdr = ($sum[1])/$fdr;
     $fdr_all{$disc} = $fdr;
     }
     else
     {
     $fdr = $sum[0];
      if($fdr>0)
      {
      $fdr = $sum[1]/$fdr;
      }
     $fdr_all{$disc} = $fdr;
     }


    $fdr_all{$disc} = $fdr;
    }
    if($results[$r][$rank]{'protein'} !~ m/^REV\_/ && $score_count{$disc}[0] && $score_count{$disc}[1] !~ m/$results[$r][$rank]{'sequence'}/)
    {
    $score_count{$disc}[0]++;
    $score_count{$disc}[1] .= $results[$r][$rank]{'sequence'};
    }
    elsif($results[$r][$rank]{'protein'} !~ m/^REV\_/ && $score_count{$disc}[1] !~ m/$results[$r][$rank]{'sequence'}/)
    {
    $score_count{$disc}[0] = 1;
    $score_count{$disc}[1] .= $results[$r][$rank]{'sequence'};
    }
    if($results[$r][$rank]{'start'}<3 && $nter_count{$disc}[0] && $results[$r][$rank]{'protein'} !~ m/^REV\_/ && $nter_count{$disc}[1] !~ m/$results[$r][$rank]{'sequence'}/)
    {
    $nter_count{$disc}[0]++;
    $nter_count{$disc}[1] .= $results[$r][$rank]{'sequence'};
    }
    elsif($results[$r][$rank]{'start'}<3 && $results[$r][$rank]{'protein'} !~ m/^REV\_/ && !$nter_count{$disc}[1] !~ m/$results[$r][$rank]{'sequence'}/)
    {
    $nter_count{$disc}[0] = 1;
    $nter_count{$disc}[1] .= $results[$r][$rank]{'sequence'};
    }

   }
  }
 }



my %counter;
my %ntercounter;
my $fdr;
 #for a set of different FDR levels
 for($fdr=0.01 ; $fdr<0.54 ; $fdr = $fdr+0.01)
 {
 my $min_score=1000;

 $counter{$fdr} = 0;
 $ntercounter{$fdr} = 0;
  foreach my $score (keys %fdr_all)
  {
   #if the fdr at this score is better than threshold
   if($fdr_all{$score} <= $fdr)
   {
    if($score<$min_score)
    {
    $min_score = $score;
    }
   $counter{$fdr} += $score_count{$score}[0];
   $ntercounter{$fdr} += $nter_count{$score}[0]; 
   }
  }
#print "FDR $fdr count $counter{$fdr} nter $ntercounter{$fdr}\n";
 #print "FDR is $fdr and score $min_score\n";
$counter{$fdr} = $counter{$fdr} - $ntercounter{$fdr};
 }

print "FDR 0.01 $counter{'0.01'} nter $ntercounter{'0.01'}\n";
print "FDR 0.05 $counter{'0.05'} nter $ntercounter{'0.05'}\n";
print "FDR 0.15 $counter{'0.15'} nter $ntercounter{'0.15'}\n";
print "FDR 0.2  $counter{'0.2'}  nter $ntercounter{'0.2'}\n"; 

$im = new GD::Image(450, 300);

# Allocate some colors
&InitColors($im);

# Make the background transparent and interlaced
$im->transparent($white);
$im->interlaced('true');

# Draw text in small font
$im->string(gdLargeFont, 2, 0 + 40, "Results Statistics", $blue);
$im->string(gdSmallFont, 2, 20 + 40, "Individual identity threshold $identity ", $black);
$im->string(gdSmallFont, 2, 35 + 40, "unique peptide matches is $identity_count ", $black);

$im->string(gdSmallFont, 2, 55 + 40, "Average identity threshold $gidentity ", $blue);
$im->string(gdSmallFont, 2, 70 + 40, "unique peptide matches is $gidentity_count ", $blue);

$im->string(gdSmallFont, 2, 90 + 40, "$expect_count[0] unique N-ter and $expect_count[1] unique NonNter expect<0.05", $black);

 #for all the fdrs
 my $pos = 110;
 foreach my $fdr (sort keys %counter)
 {
 $im->string(gdSmallFont, 2, $pos + 40, "There are $counter{$fdr} id's at FDR $fdr ($ntercounter{$fdr} nterminal)", $blue);
 $pos += 15;
 }

return $im;
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
my $setype = shift;
my @results = @_;
my @sum;

#forward
$sum[0] = 0;
#reverse
$sum[1] = 0;

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
    if($setype eq "M")
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

    if($results[$r][$rank]{'protein'} =~ m/REV\_/ && $disc >= $threshold)
    {
    $sum[1]++;
    }    
    elsif($disc >= $threshold)
    {
    $sum[0]++;
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
