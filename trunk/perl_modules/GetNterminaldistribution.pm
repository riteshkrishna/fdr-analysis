package GetNterminaldistribution;
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
our $VERSION = 1.0;
our @EXPORT = qw(GetNtermPlot);

use strict;
use URI::Escape;

#use lib "/var/www/localhost/cgi-bin/FDRAnalysis/perl_modules/";
use strict;
use CreatePNG;

sub GetNtermPlot
{
	my $se = shift;
	my $imagefile = shift;
	my $REVTAG = shift;
	my $cutoff = shift;#DCW
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

	my $min = 0;#DCW
	#my $min = 1;
	my $max = 1;

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

				print("Nterm expect value=".$results[$r][$rank]{'expect'}.", REV=".$rev."\n");

				#if($results[$r][$rank]{'expect'}<0.05)
				if($results[$r][$rank]{'expect'}<$cutoff)#DCW - allow user to choose cutoff
				{
					#DCW
					if($start>$max)
					{
						$max = $start;
					}

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

	my @total;
	$total[0] = 0;
	$total[1] = 0;

	foreach my $position (keys %nterm)
	{
		$total[0] += $nterm{$position}[0];
		$total[1] += $nterm{$position}[1];
	}

	print("Nterm max=".$max."\n");

	#DCW - split into compartments, plus N-terminal
	my $noSplits = 8;
	my %groupedData;
	#set all bar heights to zero, otherwise they will not be displayed
	#treat start of 1 or 2 specially
	$groupedData{"1-2"}[0] = 0;
	$groupedData{"1-2"}[1] = 0;
	my $tmpMax = $max-2;
	if($tmpMax<$noSplits+2)
	{
		$tmpMax=$noSplits+2;
	}
	for(my $i=0;$i<$noSplits;$i++)
	{
		#my $minPos = int($i*$max/$noSplits)+1;
		#my $maxPos = int(($i+1)*$max/$noSplits);
		my $minPos = int($i*$tmpMax/$noSplits)+3;
		my $maxPos = int(($i+1)*$tmpMax/$noSplits)+2;
		my $groupedPos = $minPos . "-" . $maxPos;
		$groupedData{$groupedPos}[0] = 0;
		$groupedData{$groupedPos}[1] = 0;
	}

	foreach my $position (keys %nterm)
	{
		#my $groupedPos = int(($position-1)*10/$max)*$max/10+$max/20;
		my $minPos;
		my $maxPos;
		if($position<=2)
		{
			$minPos = 1;
			$maxPos = 2;
		}
		else
		{
			#$minPos = int(int(($position-1)*$noSplits/$max)*$max/10)+1;
			#$maxPos = int(int((($position-1)*10/$max)+1)*$max/10);
			$minPos = int(int(($position-3)*$noSplits/$tmpMax)*$tmpMax/$noSplits)+3;
			$maxPos = int(int((($position-3)*$noSplits/$tmpMax)+1)*$tmpMax/$noSplits)+2;
		}
		my $groupedPos = $minPos . "-" . $maxPos;

		print("groupedPos=".$groupedPos.", nTerm0=".$nterm{$position}[0].", nTerm1=".$nterm{$position}[1]."\n");
		if($groupedData{$groupedPos}[0])
		{
			$groupedData{$groupedPos}[0] += $nterm{$position}[0];
		}
		else
		{
			$groupedData{$groupedPos}[0] = $nterm{$position}[0];
		}
		if($groupedData{$groupedPos}[1])
		{
			$groupedData{$groupedPos}[1] += $nterm{$position}[1];
		}
		else
		{
			$groupedData{$groupedPos}[1] = $nterm{$position}[1];
		}
	}
 
 #foreach my $position (keys %nterm)
 #{
  #if($position>$max)
  #{
  #$max = $position;
  #}
  #$plotdata[0][$position] = $position;
  # $plotdata[1][$position] = ($nterm{$position}[0]/$total[0])*100;
  # if($nterm{$position}[1] && $total[1])
  # {
  # $plotdata[2][$position] = ($nterm{$position}[1]/$total[1])*100; 
  # }
  #$count++;
 #}

	#DCW
	foreach my $position (keys %groupedData)
	{
		$plotdata[0][$position] = $position;
		$plotdata[1][$position] = ($groupedData{$position}[0]/$total[0])*100;
		if($groupedData{$position}[1] && $total[1])
		{
			$plotdata[2][$position] = ($groupedData{$position}[1]/$total[1])*100;
		}
		$count++;
	}

	my @sorted = sort { $a <=> $b } @{$plotdata[0]};

	my @sorted_plot_data;
	my $count = 0;
# for(my $s=0 ; $s<scalar(@sorted) ; $s++)
# {
#  if($sorted[$s])
#  {
#  $sorted_plot_data[0][$count] = $sorted[$s];
#  $sorted_plot_data[1][$count] = $nterm{$sorted[$s]}[0];
#  $sorted_plot_data[2][$count] = $nterm{$sorted[$s]}[1];
#  #$sorted_plot_data[1][$count] = $plotdata[1]{$sorted[$s]};
#  #$sorted_plot_data[2][$count] = $plotdata[2]{$sorted[$s]};
#  $count++;
#  }
# }

	#DCW
	for(my $s=0 ; $s<scalar(@sorted) ; $s++)
	{
		if($sorted[$s])
		{
			print("PLOT DATA:".$sorted[$s].",".$groupedData{$sorted[$s]}[0].",".$groupedData{$sorted[$s]}[1].",".$plotdata[1][$sorted[$s]].",".$plotdata[2][$sorted[$s]]."\n");

			$sorted_plot_data[0][$count] = $sorted[$s];
			#$sorted_plot_data[1][$count] = $groupedData{$sorted[$s]}[0];
			#$sorted_plot_data[2][$count] = $groupedData{$sorted[$s]}[1];
			$sorted_plot_data[1][$count] = $plotdata[1][$sorted[$s]];
			$sorted_plot_data[2][$count] = $plotdata[2][$sorted[$s]];

			$count++;
		}
	}

	#my $title = $engine . " Nterminal Distribution";
	my $title = $engine . " Peptide Start Position Distribution";
	print "in Nter score plot about to call bars\n";
	#my $image = GetPNGBars($title,'sequence position','number spectra',10,$min,$max,@sorted_plot_data);
	#my $image = GetPNGBars($title,'sequence position','% of spectra',10,$min,$max,@sorted_plot_data);
	my $image = GetPNGBars($title,'sequence position','% of spectra',$min,$max,@sorted_plot_data);#no longer need to set axis

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
