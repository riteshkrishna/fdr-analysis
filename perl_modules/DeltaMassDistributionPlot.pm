package DeltaMassDistributionPlot;
require Exporter;

#####################################################################################################
#        Copyright (C) 2010, David Wedge, University of Manchester                                  #
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
our @EXPORT = qw(GetDeltaMassDistributionPlot);

use strict;
use URI::Escape;

use strict;
use CreatePNG;

sub GetDeltaMassDistributionPlot
{
	my $imagefile = shift;
	my $scoretype = shift;
	my $setype = shift;
	my $REVTAG = shift;
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

	#print("DeltaMass results size=".scalar(@results)."\n");

	my %delta;
	my %delta1;
	my $max_tmp=0;
	for(my $r=1 ; $r<scalar(@results) ; $r++)
	{
		#print("delta mass sequence: ".$results[$r][1]{'sequence'}."\n");

		#sometimes mascot doesn't cover all sequence
		if($results[$r][1]{'sequence'} && $results[$r][1]{'sequence'} ne "NULL")
		{
			#only rank 1
			for(my $rank=1 ; $rank<2 ; $rank++)
			{

				#make the delta to 2 dp
				$results[$r][$rank]{'delta'} = sprintf("%.2f",$results[$r][$rank]{'delta'});

				#for the reverses
				if($results[$r][$rank]{'protein'} =~ m/$REVTAG/)#DCW - revtag doesn't have to be at the start
				{
					#reverse hits
					if($delta{$results[$r][$rank]{'delta'}})
					{
						$delta{$results[$r][$rank]{'delta'}} += 1;
					}
					else
					{
						$delta{$results[$r][$rank]{'delta'}} = 1;
					}
				}
				#elsif($results[$r][$rank]{'protein'})
				elsif($results[$r][$rank]{'protein'} =~ m/\S/)#protein must contain non-whitespace
				{
					#forward hits
					if($delta1{$results[$r][$rank]{'delta'}})
					{
						$delta1{$results[$r][$rank]{'delta'}} += 1;
					}
					else
					{
						$delta1{$results[$r][$rank]{'delta'}} = 1;
					}  
				}
			}
		}
	}



	my $image;
	my %plotdata;
	my $count = 0;
	#my $min = 1;
	#my $max = 0;
	#DCW
	my $min = 1000;
	my $max = -1000;
	my %all;

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


	#for the forwards
	foreach my $del (keys %delta1)
	{
		$all{$del} = 1;
		$count++;
		foreach my $del2 (keys %delta1)
		{
			if(abs($del-$del2)<=$tolerance)
			{
				if($plotdata{$del}[1])
				{
					$plotdata{$del}[1] ++;
				}
				else
				{
					$plotdata{$del}[1] = 1;
				}
				
			}
		}
	}

	#and the reverse
	foreach my $del (keys %delta)
	{
		$all{$del} = 1;
		$count++;
		foreach my $del2 (keys %delta)
		{
			if(abs($del-$del2)<=$tolerance)
			{
				if($plotdata{$del}[2])
				{
					$plotdata{$del}[2] ++;
				}
				else
				{
					$plotdata{$del}[2] = 1;
				}
				
			}
		}
	}
####################################################

	#i now want to sort the plot data
	#my @sorted = sort { $a <=> $b } @{$plotdata[0]};
	my @sorted = sort {$a <=> $b} (keys %all);
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

	#DCW add a small amount to the width, to make sure that the points are displayed
	#$max += 0.05;
	#$min -= 0.05;
	$max += 0.02;
	$min -= 0.02;
	#floor/ceil to nearest 0.1
	my $tmpmax = int($max*10)+1;
	my $tmpmin = int($min*10)-1;
	my $noLabels = $tmpmax-$tmpmin;
	if($noLabels>50)
	{
		#floor/ceil to nearest 1
		$tmpmax = int($max)+1;
		$tmpmin = int($min)-1;
		$noLabels = $tmpmax-$tmpmin;
		$max = $tmpmax;
		$min = $tmpmin;
	}
	elsif($noLabels>10)
	{
		#floor/ceil to nearest 0.5
		$tmpmax = int($max*2)+1;
		$tmpmin = int($min*2)-1;
		$noLabels = $tmpmax-$tmpmin;
		$max = $tmpmax/2.0;
		$min = $tmpmin/2.0;
	}
	else
	{
		$max = $tmpmax/10.0;
		$min = $tmpmin/10.0;
	}

	#$image = GetPNGPoints($title,'Delta Mass',$scoretype,$min,$max,@sorted_plot_data);
	$image = GetPNGLines($title,'Delta Mass Distribution',$scoretype,$noLabels,$min,$max,@sorted_plot_data); #extra input for setaxis (noLabels)
 
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
