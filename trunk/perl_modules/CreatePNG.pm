package CreatePNG;
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
our @EXPORT = qw(GetPNGBars GetPNGLines GetPNGPoints GetFinalPNG GetFinalPNGNter GetPNGVenn);

use strict;
use URI::Escape;
use GD::Graph::lines;
use GD::Graph::bars;
use GD::Graph::points;
use GD;

use lib qw(~/LIVERPOOL/bin/perl_modules/);
use strict;

my $width="450";
my $height="300";
my $fontpath = "/usr/share/fonts/truetype/comic.ttf";
my $x_max;


sub GetPNGVenn
{
my $imagefile = shift;
my @res = @_;

#labels are
#res[0][1] label 1
#res[0][2] label 2
#res[0][3] label 3

#Note the numbers to go in the venn diagram are as follows
#res[1][0] = number category 1
#res[2][0] = number category 2
#res[3][0] = number category 3
#res[0][0] = number category 1,2 & 3
#res[1][2] = number category 1 & 2
#res[1][3] = number category 1 & 3
#res[2][3] = number category 2 & 3


my $im = new GD::Image(100,100);

# allocate some colors
my $white = $im->colorAllocate(255,255,255);
my $black = $im->colorAllocate(0,0,0);

# make the background transparent and interlaced
$im->transparent($white);
$im->interlaced('true');

my $mainimage = new GD::Image((450), (300));
my $white = $mainimage->colorAllocate(255,255,255);

# Draw the first circle
$im->arc(50,50,100,100,0,360,$black);
$im->fill(50,50,$white);
$mainimage->copy($im,50,50,0,0,100,100);
$mainimage->string(gdSmallFont,10,25,$res[0][1],$black);

#draw the second circle
 if($res[0][2])
 {
 $im->arc(50,50,100,100,0,360,$black);
 $mainimage->copy($im,100,50,0,0,100,100);
 $mainimage->string(gdSmallFont,200,25,$res[0][2],$black);
 }

##draw the third
 if($res[0][3])
 {
 $im->arc(50,50,100,100,0,360,$black);
 $mainimage->copy($im,75,100,0,0,100,100);
 $mainimage->string(gdSmallFont,110,225,$res[0][3],$black);
 }

#and now for the values
 if($res[1][0])
 {
 $mainimage->string(gdSmallFont,75,75,$res[1][0],$black);
 }
 if($res[1][2])
 {
 $mainimage->string(gdSmallFont,115,75,$res[1][2],$black);
 }
 if($res[2][0])
 {
 $mainimage->string(gdSmallFont,165,75,$res[2][0],$black);
 }
 if($res[1][3])
 {
 $mainimage->string(gdSmallFont,90,125,$res[1][3],$black);
 }
 if($res[0][0]) 
 {
 $mainimage->string(gdSmallFont,110,110,$res[0][0],$black);
 }
 if($res[2][3])
 {
 $mainimage->string(gdSmallFont,150,125,$res[2][3],$black);
 }
 if($res[3][0])
 {
 $mainimage->string(gdSmallFont,125,160,$res[3][0],$black);
 }

 #leftovers
 if($res[3][3])
 {
 $mainimage->string(gdSmallFont,150,160,$res[3][3],$black);
 }

open (CHART, ">$imagefile") or return 0;
print CHART $mainimage->png;
close CHART;

return 1;
}

sub GetPNGBars
{

my $chartname = shift;
my $xlabel = shift;
my $ylabel = shift;
my $setaxis = shift;
$x_max = shift;
$x_max = $x_max + ($x_max/100*10);
my @data = @_;



my $myimage;
my $mygraph;
 #plot the values
 $mygraph = GD::Graph::bars->new($width, $height);


  $mygraph->set(
               x_label     => $xlabel,
               y_label     => $ylabel,
               title       => $chartname,
               x_max_value=>$x_max,
               bar_spacing=>1, 
               ) or warn $mygraph->error;


$mygraph = SetConstants($mygraph);

if($setaxis)
{
$mygraph->set(x_tick_number=>'10');
} 




   $mygraph->set_legend('Forward', 'Reverse');
   $myimage = $mygraph->plot(\@data) or die $mygraph->error;

return $myimage;
}


sub GetPNGLines
{
my $chartname = shift;
my $xlabel = shift;
my $ylabel = shift;
my $setaxis = shift;
$x_max = shift;
$x_max = $x_max + ($x_max/100*10);
my $type = shift;
my @data = @_;


my $myimage;
my $mygraph;

$mygraph = GD::Graph::lines->new($width, $height);


  $mygraph->set(
               x_label     => $xlabel,
               y_label     => $ylabel,
               x_max_value => $x_max,
               title       => $chartname,
               line_types  => [1, 2],
               line_width  => 1,
               ) or warn $mygraph->error;
               
$mygraph = SetConstants($mygraph);

if($setaxis)
{
$mygraph->set(x_tick_number=>'10');
}

 if($type eq "EC")
 {
 $mygraph->set_legend('estimated correct', 'estimated incorrect');
 $mygraph->set(y_label=> $ylabel);
 }
 elsif($type eq "FDR")
 {
 $mygraph->set(two_axes=>1);
 $mygraph->set(y1_label=> $ylabel);
 $mygraph->set(y2_label=> "spectra count");
 $mygraph->set_legend('FDR', 'forward spectra count');
 }


   $myimage = $mygraph->plot(\@data) or die $mygraph->error;

return $myimage;


}


sub GetPNGPoints
{
my $chartname = shift;
my $xlabel = shift;
my $ylabel = shift;
my $setaxis = shift;
$x_max = shift;
$x_max = $x_max + ($x_max/100*10);
my @data = @_;


my $myimage;
my $mygraph;

$mygraph = GD::Graph::points->new($width, $height);
  $mygraph->set(
               x_label     => $xlabel,
               y_label     => $ylabel,
               title       => $chartname,
               markers  => [3, 4],
               marker_size  => 2,
               x_max_value => $x_max,
               ) or warn $mygraph->error;

if($setaxis)
{
$mygraph->set(x_tick_number=>'10');
}


$mygraph = SetConstants($mygraph);
$mygraph->set_legend('forward', 'reverse');
$myimage = $mygraph->plot(\@data) or die $mygraph->error;

return $myimage;


}



return 1;

sub SetConstants
{
my $mygraph = shift;

#font
$mygraph->set_title_font($fontpath,12);
$mygraph->set_y_axis_font($fontpath,10);
$mygraph->set_y_label_font($fontpath,12);
$mygraph->set_x_axis_font($fontpath,10);
$mygraph->set_x_label_font($fontpath,12);
$mygraph->set_legend_font($fontpath,12);

#legend 
$mygraph->set(legend_placement=>'CB');

#x axis
$mygraph->set(x_label_position=>0.5);
$mygraph->set(zero_axis_only=>1);

#y axis
$mygraph->set( 'x_number_format' => \&x_format );


#dataset colours
$mygraph->set(dclrs=>['#336633', '#888888', 'cyan','green']);


return $mygraph;
}

sub x_format
{
my $value = shift;
my $ret;
if($x_max>5)
{ 
$ret = int($value);
}
else
{ 
$ret = sprintf("%.1f", $value);
}
return $ret;
}


sub GetFinalPNG
{

my $chartname = shift;
my @myimage = @_;

my $space = 10;

my $mainimage = new GD::Image(((2*$width)+(4*$space)), ((4*$height)+(5*$space)));
my $white = $mainimage->colorAllocate(255,255,255);


#A1
$mainimage->copy($myimage[0],$space,$space,0,0,$width,$height);
#B1
$mainimage->copy($myimage[1],($width+(2*$space)),$space,0,0,$width,$height);
#A2
$mainimage->copy($myimage[2],$space,($height+(2*$space)),0,0,$width,$height);
#B2
$mainimage->copy($myimage[3],($width+(2*$space)),($height+(2*$space)),0,0,$width,$height);
#A3
$mainimage->copy($myimage[4],$space,((2*$height)+(3*$space)),0,0,$width,$height);
#B3
$mainimage->copy($myimage[5],($width+(2*$space)),((2*$height)+(3*$space)),0,0,$width,$height);
#A4
#$mainimage->copy($myimage[6],$space,((3*$height)+(4*$space)),0,0,$width,$height);

open (CHART, ">$chartname") or print "unable to open the $chartname file ";
print CHART $mainimage->png;
close CHART;

}


sub GetFinalPNGNter
{

my $chartname = shift;
my @myimage = @_;

my $space = 10;

my $mainimage = new GD::Image(((2*$width)+(4*$space)), ((4*$height)+(5*$space)));
my $white = $mainimage->colorAllocate(255,255,255);


#A1
$mainimage->copy($myimage[0],$space,$space,0,0,$width,$height);
#B1
$mainimage->copy($myimage[1],($width+(2*$space)),$space,0,0,$width,$height);
#A2
$mainimage->copy($myimage[2],$space,($height+(2*$space)),0,0,$width,$height);
#B2
$mainimage->copy($myimage[3],($width+(2*$space)),($height+(2*$space)),0,0,$width,$height);
#A3
$mainimage->copy($myimage[4],$space,((2*$height)+(3*$space)),0,0,$width,$height);
#B3
$mainimage->copy($myimage[5],($width+(2*$space)),((2*$height)+(3*$space)),0,0,$width,$height);
#A4
$mainimage->copy($myimage[6],$space,((3*$height)+(4*$space)),0,0,$width,$height);

open (CHART, ">$chartname") or print "unable to open the $chartname file ";
print CHART $mainimage->png;
close CHART;

}
return 1;




