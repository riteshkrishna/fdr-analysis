package ZoomDiscriminantScorePlot;
require Exporter;


#####################################################################################################
#        Copyright (C) 2009, Jennifer Siepen and David Wedge, University of Manchester              #
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

#use lib "/var/www/localhost/cgi-bin/FDRAnalysis/perl_modules/";
use strict;
use CreatePNG;

sub ZoomScoreDistribution
{
	my $imagefile = shift;
	my $scoretype = shift;
	my $setype = shift;
	my $revtag = shift;
	my $decoy_size = shift;
	my @results = @_;

	my @graph;
	my %tmpdistribution;
	my %distribution;
	my $min=1000;
	my $max = 0;

	my $binWidth = 5;

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
	elsif($setype eq "C")
	{
		$engine = "Consensus";
		$binWidth = 0.1;
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

				if($setype eq "O")
				{
					#DCW - 211209
					$disc = $results[$r][$rank]{'expect'};
					#$disc = -10*log($disc)/log(10);#convert e-value to score
					$disc = -log($disc);#convert e-value to score
					$disc = sprintf("%.1f",$disc);
				}
				elsif($setype eq "C")
				{
					$disc = $results[$r][$rank]{'FDR_score'};
				}
				else
				{
					$disc = int($disc);
				}

				my $rev = 0;
				if($results[$r][$rank]{'protein'} =~ m/$revtag/)#DCW - variable revtag
				{
					$rev = 1;
				} 

				if($results[$r][$rank]{'protein'} =~ m/\S/)#protein must contain non-whitespace
				{
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
	}

	#put the distribution into bins
	my @totals;
	foreach my $score (keys %distribution)
	{
		my $newscore;
		if($binWidth>=1)#condition added 030310
		{
			if($score<0)
			{
				$newscore =  -int(-$score/$binWidth+0.5) * $binWidth;
			}
			else
			{
				$newscore =  int($score/$binWidth+0.5) * $binWidth;
			}
		}
		else
		{
			$newscore = sprintf("%.1f",$score);
		}

		if($newscore>$max)
		{
			$max = $newscore;
		}
		if($newscore<$min)
		{
			$min = $newscore;
		}

		if($tmpdistribution{$newscore}[0])
		{
			$tmpdistribution{$newscore}[0] += $distribution{$score}[0];
			$totals[0]+= $distribution{$score}[0];
		}
 		else
		{
			$tmpdistribution{$newscore}[0] = $distribution{$score}[0];
		}
 
		if($tmpdistribution{$newscore}[1])
		{
			$tmpdistribution{$newscore}[1] += $distribution{$score}[1]/$decoy_size; #scale according to decoy size
			$totals[1]+= $distribution{$score}[1];
		}
		else
		{
			$tmpdistribution{$newscore}[1] = $distribution{$score}[1]/$decoy_size; #scale according to decoy size
		}
	}

	$min = $min-$binWidth;
	$max = $max+$binWidth;

	if($setype eq "C")
	{
		$min = 0.0;
		$max = 1.0;
	}

	my %plotdata;
	my @scorearray;
	#add points with zero value
	if($binWidth>=1)#condition added 020310
	{
		for(my $i = $min; $i<=$max; $i += $binWidth)
		{
			$plotdata{$i}[0] = $i;
			$plotdata{$i}[1] = 0;
			$plotdata{$i}[2] = 0;
			push(@scorearray,$i);
		}
	}
	else
	{
		for(my $i=int($min*10);$i<=int($max*10);$i++)
		{
			$plotdata{$i/10.0}[0] = $i/10.0;
			$plotdata{$i/10.0}[1] = 0;
			$plotdata{$i/10.0}[2] = 0;
			push(@scorearray,$i/10.0);
		}
	}

	foreach my $score (keys %tmpdistribution)
	{
		$plotdata{$score}[0] = $score;
	   	if($tmpdistribution{$score}[0] && $tmpdistribution{$score}[1])
   		{
 			#this will be the estimated correct
			my $estimatedCorrect = $tmpdistribution{$score}[0]-$tmpdistribution{$score}[1];
			#print("score=".$score."est correct=".$estimatedCorrect."\n");
			#DCW - don't use values less than zero
			if($estimatedCorrect > 0)
			{
				$plotdata{$score}[1] = $estimatedCorrect;
			}
			else
			{
				$plotdata{$score}[1] = 0;
			}
		}
		elsif($tmpdistribution{$score}[0])
		{
			#print("score=".$score."forward dist=".$tmpdistribution{$score}[0]."\n");
			$plotdata{$score}[1] = $tmpdistribution{$score}[0];
		}
		if($tmpdistribution{$score}[1])
		{
			#print("score=".$score."decoy dist=".$tmpdistribution{$score}[0]."\n");
			$plotdata{$score}[2] = $tmpdistribution{$score}[1];
		}
	}

	#my @sorted = sort { $a <=> $b } @{$plotdata[0]};
	my @sorted = sort {$a <=> $b} (@scorearray);
	my @sorted_plot_data;

	my $count = 0;
	#my $count = 1;#DCW - add extra point at the start with zero value

	for(my $s=0 ; $s<scalar(@sorted) ; $s++)
	{
		$sorted_plot_data[0][$count] = $sorted[$s];
		$sorted_plot_data[1][$count] = $plotdata{$sorted[$s]}[1];
		$sorted_plot_data[2][$count] = $plotdata{$sorted[$s]}[2];
		$count++;
	}

	if($setype eq "O")
	{
		#$scoretype = "-10log(expect)";
		$scoretype = "-log(expect)";
	}
	elsif($setype eq "C")
	{
		$scoretype = "FDR";
	}
	elsif($setype eq "T")
	{
		$scoretype = "hyperscore";
	}
	else
	{
		#$scoretype .= " score";
		$scoretype = "score";
	}

	my $title = $engine . " Estimation of Correct/Incorrect";
	my $image = GetPNGLines($title,$scoretype,'number spectra',@sorted-1,$min,$max,'EC',@sorted_plot_data);
	#my $image = GetPNGLines($title,$scoretype,'number spectra',0,$max,'EC',@sorted_plot_data);

	my $mainimage = new GD::Image((450), (300));
	my $white = $mainimage->colorAllocate(255,255,255);

	$mainimage->copy($image,0,0,0,0,450,300);
	open (CHART, ">$imagefile") or print "unable to open the $imagefile file ";
	print CHART $mainimage->png;
	close CHART;

	return 1;
} 

return 1;

#DCW 190110- not sure what this is for!
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
