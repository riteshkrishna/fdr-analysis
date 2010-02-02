package GetDiscriminantScorePlot;
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
our @EXPORT = qw(GetScoreDistribution);

use strict;
use URI::Escape;

#use lib "/var/www/localhost/cgi-bin/FDRAnalysis/perl_modules/";
use strict;
use CreatePNG;


sub GetScoreDistribution
{
	my $imagefile = shift;
	my $scoretype = shift;
	my $setype = shift;
	my $REVTAG = shift;
	my @results = @_;

	my @graph;
	my %distribution;
	#my $min=1;
	my $min=1000;
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

				if($setype ne "O")
				{
					#$disc = int($disc);
				}
				else
				{
					#DCW - 211209
					$disc = $results[$r][$rank]{'expect'};
					$disc = -log($disc);

				}

				#DCW 190110 - convert to integer value, allowing for negative numbers
				if($disc<0)
				{
					$disc=-int(-$disc-0.000001)-1;
				}
				else
				{
					$disc=int($disc);
				}

				if($disc>$max)
				{
					$max = $disc;
				}
				#DCW
				if($disc<$min)
				{
					$min = $disc;
				}

				my $rev = 0;

				if($results[$r][$rank]{'protein'} =~ m/$REVTAG/)
				{
					$rev = 1;
				}   

				if($results[$r][$rank]{'protein'} =~ m/\S/)#protein must contain non-whitespace
				{
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
	}

	#print("discriminant score plot: max=".$max.", min=".$min."\n");

	#DCW 190110 - floor/ceil no nearest $labelInterval
	my $labelInterval = 5;
	if($min<0)
	{
		$min = -(int((-$min-0.000001)/$labelInterval)+1)*$labelInterval;
	}
	else
	{
		$min = (int(($min-0.000001)/$labelInterval)+1)*$labelInterval;
	}
	if($max<0)
	{
		$max = -(int((-$max-0.000001)/$labelInterval)+1)*$labelInterval;
	}
	else
	{
		$max = (int(($max-0.000001)/$labelInterval)+1)*$labelInterval;
	}

	my @plotdata;
	#DCW 190110
	for(my $i=$min;$i<=$max;$i++)
	{
		if($i % $labelInterval == 0)
		{
		  	$plotdata[0][$i-$min] = $i;
		}
		else
		{
			$plotdata[0][$i-$min] = '';
		}
	}

	my @totals;
	foreach my $score (keys %distribution)
	{
		if($distribution{$score}[0])
		{
			$totals[0] += $distribution{$score}[0];
		}
		if($distribution{$score}[1])
		{
			$totals[1] += $distribution{$score}[1];
		}
	}

	#DCW 190110
	foreach my $score (keys %distribution)
	{
		if($distribution{$score}[0])
			{
			#$plotdata[1][$score-$min] = $distribution{$score}[0];
			$plotdata[1][$score-$min] = $distribution{$score}[0]/$totals[0]*100;#convert to %
		}
		if($distribution{$score}[1])
		{
			#$plotdata[2][$score-$min] = $distribution{$score}[1];
			$plotdata[2][$score-$min] = $distribution{$score}[1]/$totals[1]*100;#convert to %
		}
	}

	if($setype ne "O")
	{
		#$scoretype .= " score";
		$scoretype = "score";
	}
	else
	{
		$scoretype = "-log(expect)";
	}
	my $title = $engine . " Score Distribution";

	#my $image = GetPNGBars($title,$scoretype,'number spectra',$min,$max,@sorted_plot_data);
	#my $image = GetPNGBars($title,$scoretype,'number spectra',0,$max,@sorted_plot_data);
	#my $image = GetPNGBars($title,$scoretype,'% of spectra',$min,$max,@plotdata);#DCW - data are already sorted
	my $image = GetPNGBars($title,$scoretype,'% of spectra',0,$max,@plotdata);#DCW - min value messes up graph

	my $mainimage = new GD::Image((450), (300));
	my $white = $mainimage->colorAllocate(255,255,255);

	$mainimage->copy($image,0,0,0,0,450,300);
	open (CHART, ">$imagefile") or print "unable to open the $imagefile file ";
	print CHART $mainimage->png;
	close CHART;
	return 1;
} 

return 1;

#DCW -not sure what this is for
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
