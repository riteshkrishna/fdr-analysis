#!/usr/bin/perl


####################################################################################################
#	 Copyright (C) 2009, Andrew Jones, University of Liverpool									   #
#    If you use this software, please cite the following paper:									   #
#    Jones, A. R., Siepen, J. A., Hubbard, S. J., Paton, N. W., Improving sensitivity in proteome  #
#    studies by analysis of false discovery rates for multiple search engines. 					   #
#    PROTEOMICS 2009, 9, 1220-1229.		   														   #	
#   																							   #
#    																							   #
#    This program is free software: you can redistribute it and/or modify						   #
#    it under the terms of the GNU General Public License as published by						   #
#    the Free Software Foundation, either version 3 of the License, or							   #
#    (at your option) any later version.														   #
#						 																		   #
#    This program is distributed in the hope that it will be useful,							   #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of								   #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the								   #	
#    GNU General Public License for more details.												   #
#																								   #
#    You should have received a copy of the GNU General Public License							   #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.						   #
#   																							   #
#   																							   #
####################################################################################################



use strict;

## This code is for calculating the point on a curve when the points are closest to linear and reporting the intercept

sub calculateGradients{

	my $tmp0 = shift;
	my $tmp1 = shift;
	my @x_values = @{$tmp0};
	my @y_values = @{$tmp1};
	
	my %slopes;
	my %mean_vals;
	#my %yvalues;
	my %std_y;

	for(my $i =0; $i<@x_values; $i++){
	
		#print $x_values[$i] . "\t" . $y_values[$i] . "\n";
	
		 
		my $tmp_xval = @x_values[$i];
		if($i > 2 && $i < @x_values-2 && $tmp_xval != @x_values[$i-1] && $tmp_xval !=  @x_values[$i-2]&& $tmp_xval != @x_values[$i+1] && $tmp_xval !=  @x_values[$i+2]){
			my %tmp = {};
			$tmp{@x_values[$i-2]} = @y_values[$i-2];
			$tmp{@x_values[$i-1]} = @y_values[$i-1];
			$tmp{@x_values[$i]} = @y_values[$i];
			$tmp{@x_values[$i+1]} = @y_values[$i+1];
			$tmp{@x_values[$i+2]} = @y_values[$i+2];
			
			
			my ($slope, $intercept, $rsq, $mean_y,$stddev_y) = linearRegression({%tmp});
			
			my $currentXval = @x_values[$i];
			#print "Current $currentXval $slope $intercept\n";
			
			
			#$evalSlopes{$currEval} = $slope;
			#$evalIntercepts{$currEval} = $intercept;
			
			$slopes{$currentXval} = $slope;
			#$yvalues{$currentXval} = $intercept;
			$mean_vals{$currentXval} = $mean_y;
			$std_y{$currentXval} = $stddev_y;
			
		}
		
	}

	#Now calculate linear point
	#print "Slopes\n";
	
	my $prevSlope = 0;
	
	my $foundXVal;
	my $foundSlope = 10000;
	
	my $mean_y=0;
	
	for my $xval ( keys %slopes ) {
	
		my $currSlope = $slopes{$xval};	
		$mean_y = $mean_vals{$xval};
		#print $xval . "\t" . $currSlope . "\t" . $intercept . "\t".  $prevSlope. "\n";
		#print "x:$xval\tslope:$currSlope\tmeanx:$mean_y\n";
	
		$prevSlope = $currSlope;
		
		if(abs($currSlope) < abs($foundSlope) ){
			$foundSlope = $currSlope;
			$foundXVal = $xval;
			#print "Found: intercept: " . $intercept . " for xval: $xval \n\n";
			# !!!  Note: not exactly the intercept that I want... For this point on the graph, want to know the y-value !!
			
		}		
	}
	
	#print "intercept for slope closest to zero: $foundIntercept";
	$mean_y = $mean_vals{$foundXVal};
	my $std_y = $std_y{$foundXVal};
	#print "Xval at closest to linear: $foundXVal , slope: $foundSlope , mean:$mean_y stdev: $std_y \n";
		
	return $mean_y,$foundSlope,$std_y;
}




#return slope and intercept of best fit line for values of x and y
#return slope and intercept of best fit line for values of x and y
sub linearRegression{

	my $temp = shift;
	
	#my @xvalues = @{$tmpX};
	#my @yvalues = @{$tmpY};
	
	my %values = %{$temp};
	
	my $newslope = 0;
	my $newintercept = 0;
	
	my $sumProductXY = 0;
	my $sumX = 0;
	my $sumY = 0;
	my $sumSquareX = 0;
	my $sumSquareY = 0;
	
		
	#if(@xvalues != @yvalues){
	#	print "Error, supplied X and Y value arrays are not the same length!\n";
	#	exit;
	#}
	
	#my $n = keys %values;
	my $n = 0;
	for my $x ( keys %values ) {
		
		my $y = $values{$x};
		
		#print "x: $x y: $y \n"; 
		
		if($y ne "INF" && $y){
			$sumX = $sumX + $x;
			$sumY = $sumY + $y;
			$sumSquareX = $sumSquareX + ($x * $x);
			$sumSquareY = $sumSquareY + ($y * $y);
			$sumProductXY = $sumProductXY + ($x * $y);
			$n++;
		}
    }	
		
	#print "SumX $sumX\n";
	#print "SumY $sumY\n";
	#print "Sum squareX: $sumSquareX\n";
	#print "Sum productXY: $sumProductXY\n";
	#print "N: $n \n";
			
	my $mean_ProductXY = 	$sumProductXY / $n;
	my $mean_x = $sumX / $n;
	my $mean_y = $sumY / $n;
	my $mean_X_Cross_mean_Y = $mean_x * $mean_y;
	
	my $stddev_x = sqrt( ($sumSquareX/$n) - (($sumX/$n)*($sumX/$n))  );
	my $stddev_y = sqrt( ($sumSquareY/$n) - (($sumY/$n)*($sumY/$n))  );
		
	my $rsquared;
	if($stddev_x * $stddev_y != 0){
		$rsquared	= ($mean_ProductXY - $mean_X_Cross_mean_Y) / ($stddev_x * $stddev_y );
	}
	else{
		$rsquared	= 0;
	}
		
		
	$newslope = (($n * $sumProductXY) - ($sumX * $sumY)) / (($n * $sumSquareX) - ($sumX * $sumX));
	$newintercept = ($sumY - ($newslope * $sumX)) / $n;
	
	#print "M: $mean_y\t";
	#print "S: $newslope\t";
	#print "I: $newintercept\t";
	#print "R: $rsquared\t";
	#print "STD: $stddev_y\n";
	
	#So var(x) is by definition <x^2>-<x>^2. In other words, first square 
	#the x's and average them. Then average the x's and square the result.  
	#Take the square root of the difference to get the standard deviation. 
	#I get 1.118 for this.
	
	 #<xy>-<x><y>
   #-------------- = r = correlation coefficient
    #std(x)*std(y
		
	
	return ($newslope, $newintercept, $rsquared, $mean_y, $stddev_y);
}

1;