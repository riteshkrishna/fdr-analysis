#!/usr/bin/perl -w

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

#####################################################################################################
#        WARNING: This script should be used with care for security reasons                         #
#####################################################################################################

use strict;

my $web_cgipath = "http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/";#the web location of the perl scripts

my $send_to = shift;
my $sid = shift;
my $files = shift;
my @file_array = split(",",$files);
#print("file list length=".@file_array."\n");

my $sendmail = "/usr/sbin/sendmail -t";
my $from = "From: noreply\@manchester.ac.uk\n";
my $subject  = "Subject: FDR analysis report\n";
my $content  = "Dear FDR Analysis web user,\nThe results of your FDR analysis search are available at\n".$web_cgipath."FDR_analysis_static.pl?".$sid."\nThe results will available for a limited period (up to a week).\n\nBest wishes,\nThe FDR Analysis Team";
#my $content  = "Dear FDR Analysis web user,\nThe results of your FDR analysis search will be available at\n".$web_cgipath."FDR_analysis_static.pl?sid=".$sid.".\nIt will take a few minutes for the analysis to complete and, once complete, the results will available for a limited period (up to a week).\n\nBest wishes,\nThe FDR Analysis Team";

foreach my $file(@file_array)
{
	if(! -e $file || ! -s $file)
	{
		$content="Dear FDR Analysis web user,\nThere was a problem running your FDR analysis. Please try again.\n\nBest wishes,\nThe FDR Analysis Team";
		last;
	}
}

$send_to = "To: " . $send_to."\n";

open(SENDMAIL, "|$sendmail") or die "Cannot open $sendmail: $!";
print SENDMAIL $from;
print SENDMAIL $send_to;
print SENDMAIL $subject;
print SENDMAIL "Content-type: text/plain\n\n";
print SENDMAIL $content;
close(SENDMAIL);