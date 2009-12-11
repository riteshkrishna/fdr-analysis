#!/usr/bin/perl -w

#####################################################################################################
#        Copyright (C) 2009, Jennifer Siepen & David Wedge, University of Manchester                              #
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


use strict;
use DBI;
use CGI;
use CGI::Session;
use FileHandle;
use Benchmark;

my $QSUB = "/opt/sge/bin/lx26-amd64/qsub";
$ENV{'SGE_ROOT'} = "/opt/sge";
$ENV{'SGE_CELL'} = "ispider";


my $cgi = new CGI;
my $sid;
my $debugger = 0;


if($cgi->cookie("CGISESSID"))
{
$sid = $cgi->cookie("CGISESSID");
$debugger = 2;
}
else
{
$sid = undef;
$debugger = 1;
}

my $session = new CGI::Session(undef, $sid, {Directory=>'/var/www/tmp/'});
my $cookie = $cgi->cookie(CGISESSID => $session->id );
print $cgi->header(-cookie => $cookie );

if ($cgi->param("submit_clicked"))
{
	$session->clear();
	$session->flush();
}

$sid = $session->id();

main_content($cgi->param("refresh"));


footer();

sub main_content()
{
	my $refresh = shift;
	my $results_ready = 1;
	#print("refresh:".$refresh);
	my $file1;
	my $file2;
	my $file3;
	my $combine_param;
	if($refresh eq "true")
	{
		$file1=$cgi->param("input1");
		$file2=$cgi->param("input2");
		$file3=$cgi->param("input3");
		$combine_param = $cgi->param("combine");
		$session->param("paramsset",0);
		$results_ready=0;
	}

	my $refresh_browser;
	my $true = 0;
	my @FDR_image_file;

	if($file1 ||($refresh ne "true" && $session->param("mascot") == 1))
	{
		#print("checking file1");

		$FDR_image_file[0] = "/var/www/localhost/htdocs/FDRAnalysis/tmp/mascot_fdranalysis_" . $session->id . ".png";
		if($refresh eq "true")
		{
			#$FDR_image_file[0] = "/var/www/localhost/htdocs/FDRAnalysis/tmp/mascot_fdranalysis_" . $session->id . ".png";
			my $tmp = $FDR_image_file[0];
			$tmp =~ s/\.png/GygiRank\.png/;
 			if(-e $tmp)
			{
				#REMOVE IMAGE FILE
				my $cmd = "rm $tmp";
				system($cmd);
			}
		}
		else
		{
			my $tmp = $FDR_image_file[0];
			$tmp =~ s/\.png/GygiRank\.png/;
 			if(!-e $tmp)
 			{
 				$results_ready = 0;
				$refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl\">";
			}
		}
		$session->param('mascot', 1);
		$session->flush();
	}
	else
	{
		$session->param('mascot', 0);
		$session->flush();
	}
	
	if($file2 ||($refresh ne "true" && $session->param("omssa") == 1))
	{
		$FDR_image_file[1] = "/var/www/localhost/htdocs/FDRAnalysis/tmp/omssa_fdranalysis_" . $session->id . ".png";
		if($refresh eq "true")
		{
			#$FDR_image_file[1] = "/var/www/localhost/htdocs/FDRAnalysis/tmp/omssa_fdranalysis_" . $session->id . ".png";
			my $tmp = $FDR_image_file[1];
			$tmp =~ s/\.png/GygiRank\.png/;
 			if(-e $tmp)
			{
				#REMOVE IMAGE FILE
				my $cmd = "rm $tmp";
				system($cmd);
			}
		}
		else
		{
			my $tmp = $FDR_image_file[1];
			$tmp =~ s/\.png/GygiRank\.png/;
 			if(!-e $tmp)
 			{
 				$results_ready = 0;
				$refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl\">";
			}
		}
		$session->param('omssa', 1);
		$session->flush();
	}
	else
	{
		$session->param('omssa', 0);
		$session->flush();
	}
	
	if($file3 || ($refresh ne "true" && $session->param("tandem") == 1))
	{
		$FDR_image_file[2] = "/var/www/localhost/htdocs/FDRAnalysis/tmp/tandem_fdranalysis_" . $session->id . ".png";
		if($refresh eq "true")
		{
			#$FDR_image_file[2] = "/var/www/localhost/htdocs/FDRAnalysis/tmp/tandem_fdranalysis_" . $session->id . ".png";
			my $tmp = $FDR_image_file[2];
			$tmp =~ s/\.png/GygiRank\.png/;
 			if(-e $tmp)
			{
				#REMOVE IMAGE FILE
				my $cmd = "rm $tmp";
				system($cmd);
			}
		}
		else
		{
			my $tmp = $FDR_image_file[2];
			$tmp =~ s/\.png/GygiRank\.png/;
 			if(!-e $tmp)
 			{
 				$results_ready = 0;
				$refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl\">";
			}
		}
		$session->param('tandem', 1);
		$session->flush();
	}
	else
	{
		$session->param('tandem', 0);
		$session->flush();
	}

	#DCW - tmp2 taken out of 'if'
	my $tmp2 = "/var/www/tmp/" . $session->id . "_summary.txt";
	if($refresh eq "true")
	{
		my $cmd = "rm $tmp2";
		system($cmd);
	}
   	elsif(!-e $tmp2)
   	{
		$results_ready = 0;
   		$refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl\">";
  	}

	#if($combine_param || $session->param("combine") == 1)
	if($combine_param)
	{
		my $tmp = "/var/www/tmp/" . $session->id . "combined_peptides.out";
		if($refresh eq "true")
		{
			my $cmd = "rm $tmp";
			system($cmd);
		}
   		elsif(!-e $tmp || -z $tmp)
   		{
			$results_ready = 0;
   			$refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl\">";
  		}
  		$session->param('combine', 1);
  		$session->flush();
  	}
  	else
  	{
  		$session->param('combine', 0);
  		$session->flush();
	}

	$refresh="false";
	
	if(!$session->param("paramsset"))
	{
		#print "setting params\n";

		#moved to here 091209
		if(!$session->param("tandem") && !$session->param("omssa") && !$session->param("mascot"))
  		{
			ErrorMsg("You need to run a search before you can view the results.\n");
			exit(1);
		}

		$session->param("paramsset",1);
		$results_ready = 0;
		$refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl\">";

		my $params = "/var/www/tmp/fdranalysis_param_" . $session->id . ".txt";
		my $cmd = "rm $params";
		system($cmd);

		##Mascot
		if($file1)
		{
			#deal with the rev_file
			#$cgi->param('rev_file', DealWithFileUpload($cgi->param("input1"),"mascot_rev"));
			my $tmpFile = DealWithFileUpload($cgi->param("input1"),"mascot_rev");#091209
			#write the Mascot parameters to the param file
			open(PARAM,">$params") or die "unable to open the file $params to write to\n";
			#print PARAM "rev_file\t" . $cgi->param("rev_file") . "\n";
			print PARAM "rev_file\t" . $tmpFile . "\n";#091209
			print PARAM "mascot_search_type\tconcatenated forward\n";
			print PARAM "image\t" . $FDR_image_file[0] . "\n";
			close PARAM;
		}

		##Omssa
		if($file2)
		{
			#deal with the rev_file
			#$cgi->param('omssa_rev_file', DealWithFileUpload($cgi->param("input2"),"omssa_rev"));
			my $tmpFile = DealWithFileUpload($cgi->param("input2"),"omssa_rev");#091209
			#write the Omssa parameters to the param file
			open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
			#print PARAM "omssa_rev_file\t" . $cgi->param("omssa_rev_file") . "\n";
			print PARAM "omssa_rev_file\t" . $tmpFile . "\n";#091209
			print PARAM "omssa_search_type\tconcatenated forward\n";
			print PARAM "omssa_image\t" . $FDR_image_file[1] . "\n";
			close PARAM;
		}

		##Tandem
		if($file3)
		{
			#deal with the rev_file
			my $tmpFile = DealWithFileUpload($cgi->param("input3"),"tandem_rev");#091209
			#write the Tandem parameters to the param file
			open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
			#print PARAM "tandem_rev_file\t" . $cgi->param("tandem_rev_file") . "\n";
			print PARAM "tandem_rev_file\t" . $tmpFile . "\n";#091209
			print PARAM "tandem_search_type\tconcatenated forward\n";
			print PARAM "tandem_image\t" . $FDR_image_file[2] . "\n";
			close PARAM;
		}		
		
		my $tmpParam;#091209

		$session->param("nterminal", $cgi->param("nter"));
		open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
		#print PARAM "analysis_type\t" . $cgi->param("nter") . "\n";
		print PARAM "analysis_type\t" . $session->param("nterminal") . "\n";#091209
		close PARAM;

		#$tmpParam=$cgi->param("combine");
		if($combine_param)
		{
			open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
			#print PARAM "oversize\t" . $cgi->param("oversize") . "\n";
			#print PARAM "combine\t" . $cgi->param("combine") . "\n";
			print PARAM "combine\t" . $combine_param . "\n";
			close PARAM;
		}
		open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
		print PARAM "oversize\t" . $cgi->param("oversize") . "\n";
		close PARAM;

		$tmpParam = $cgi->param("max_expect");
	   	$session->param('max_expect', $tmpParam);
		$tmpParam = $cgi->param("fdr_value");
	   	$session->param('fdr_value', $tmpParam);
		$tmpParam = $cgi->param("rev_tag");
		$session->param('rev_tag', $cgi->param("rev_tag"));
	   	open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
	   	print PARAM "max_expect\t" . $session->param("max_expect") . "\n";
	   	print PARAM "fdr_value\t" . $session->param("fdr_value") . "\n";
	   	print PARAM "rev_tag\t" . $session->param("rev_tag") . "\n";
	   	close PARAM;

		if(!$session->param("max_expect"))
		{
			ErrorMsg("You must set a maximum expectation value to use\n");
			exit(1);
		}
		if(!$session->param("fdr_value"))
		{
			ErrorMsg("You must set a maximum FDR value to use\n");
			exit(1);
		}
		if($session->param("fdr_value")>0.5)
		{
			ErrorMsg("The maximum FDR value is 0.5.\n");
			exit(1);
		}
		
		$session->flush();
		RunAnalysis($params);		
	}
	

	start_html($refresh_browser);

	if($results_ready == 0)
	{
		print  $cgi->div({id=>"main_result"},
			$cgi->div({id=>"content"},
                     $cgi->br(),
			$cgi->p("The analysis is running, please be patient this shouldn't take too long."),
			$cgi->p("This page will automatically refresh, if it doesn't please click ",$cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl"},"here")),
 			$cgi->br(),
			$cgi->img({src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/Waiting.gif"}),
			$cgi->br(),
             		)
       	);
 	} 
	else
	{
		#print $cgi->h2("initial image filename= $FDR_image_file[0]");
		#change all image files to a url
		for(my $i=0 ; $i<scalar(@FDR_image_file) ; $i++)
		{ 
			$FDR_image_file[$i] =~ s/\/var\/www\/localhost\/htdocs/http:\/\/ispider.smith.man.ac.uk/;
		}
		my $final_image_line;
		my $blurb;
		my $title;
		if(!$cgi->param("result_view") || $cgi->param("result_view") eq "summary")
		{
			#only one summary file
			my $summaryfile = "/var/www/tmp/" . $session->id . "_summary.txt";
			open(FILE,"<$summaryfile") or print "can't open the summary $summaryfile\n";
			while(my $line = <FILE>)
			{
				$blurb .= $line;
				$blurb .= "<BR><BR>";
			}
			close FILE;

			my $newpeplist = "/var/www/localhost/htdocs/FDRAnalysis/tmp/" . $session->id . "_peptidelist.txt";
			my $cmd = "cp /var/www/tmp/" . $session->id . "_peptidelist.txt " .  $newpeplist;
			system($cmd);
			$newpeplist =~ s/\/var\/www\/localhost\/htdocs\//http:\/\/ispider\.smith\.man\.ac\.uk\//;
			$blurb .= "<BR><font size=-2><I><B>E</B> FDR method (Elias & Gygi (2007) Nat. Methods 4: 207-214) <BR> <B>K</B> FDR method (Kall et al. (2008) J. Proteome Res. 7:29-34).</I></FONT>"; 
			$blurb .= "<BR><font size=-2><I><B>*</B> Combined method as described in Jones et al. (2009) Proteomics 9: 1220-9.</I></FONT>"; 
			$blurb .=  "<BR><BR><P>To download a list of identified peptides/proteins please click <a href=\"$newpeplist\" target=\"top\">here</a><BR>";

			if($session->param("combine") == 1)
			{
				my $combinedlist = "/var/www\/localhost\/htdocs\/FDRAnalysis\/tmp\//" . $session->id . "combined_peptides.out";
				$combinedlist =~ s/\/var\/www\/localhost\/htdocs\//http:\/\/ispider\.smith\.man\.ac\.uk\//;
				$blurb .=  "<BR><P>To download a list of identified peptides from the combined analysis please click <a href=\"$combinedlist\" target=\"top\">here</a><BR>";
			}

			$title = "Results Summary";
		}
		elsif($cgi->param("result_view") eq "rankplot")
		{
			#for all results sets
			for(my $i=0 ; $i<scalar(@FDR_image_file) ; $i++)
			{
				$FDR_image_file[$i] =~ s/\.png/GygiRank\.png/;
				$final_image_line .= "<IMG SRC=\"$FDR_image_file[$i]\" title=\"Rank plot image\">";
			}
			$blurb = "<BR><P>This plot is the percentage of forward and reverse hits at different identified ranks where the rank 1 peptide has an expect value of " . $session->param("max_expect") . " or less.  Rank 1 should have very few reverse hits whereas ranks 2-> should have approximately 50% forward and 50% reverse (as these represent false positives)(based on that of Elias & Gygi (2007) Nat. Methods 4: 207-214).";
			$title = "Rank Plot";
		}
		elsif($cgi->param("result_view") eq "deltamass")
		{
			for(my $i=0 ; $i<scalar(@FDR_image_file) ; $i++) 
			{
				$FDR_image_file[$i] =~ s/\.png/DeltaMass\.png/;
				$final_image_line .=  "<IMG SRC=\"$FDR_image_file[$i]\" title=\"Rank plot image\">";
			}
			$blurb = "<BR><P>This image represents a plot of the mass differences of all identifications vs. score between experimental and theoretical masses.  A normal distribution around the 0 axis is expected";
			$title = "Delta Mass";
		}
		elsif($cgi->param("result_view") eq "nterdist")
		{
			for(my $i=0 ; $i<scalar(@FDR_image_file) ; $i++)
                     {
				$FDR_image_file[$i] =~ s/\.png/NterDist\.png/;
				$final_image_line .=  "<IMG SRC=\"$FDR_image_file[$i]\" title=\"Start position  plot image\">"; 
			}
			$blurb = "<BR><P>This plot is part of the N-terminal proteome analysis pipeline.  This analysis specifically targets the N-terminal peptides and as a result the forward hits should be located predominantly at positions 1 and 2 and the reverse scattered randomly amongst all other positions. (Note max spectra shown on y-axis is 100 for a clearer view).";
			$title = "N-terminal distribution";
		}
		elsif($cgi->param("result_view") eq "scoredist")
		{
			for(my $i=0 ; $i<scalar(@FDR_image_file) ; $i++)
                     {
				$FDR_image_file[$i] =~ s/\.png/ScoreDist\.png/;
				$final_image_line .=  "<IMG SRC=\"$FDR_image_file[$i]\" title=\"Score distribution plot image\">";
			}
			$blurb = "<BR><P>This plot is the distribution of scores for both the forward and reverse hits.  Two normal distributions should be evident, a distribution for reverse scores centered around a poor score and a seperate distribution for the forward hits, centered around a better score.";
			$title = "Score distribution";
		}
		elsif($cgi->param("result_view") eq "zoomscore")
		{
			for(my $i=0 ; $i<scalar(@FDR_image_file) ; $i++)
			{
				$FDR_image_file[$i] =~ s/\.png/ZoomScore\.png/;
				$final_image_line .=  "<IMG SRC=\"$FDR_image_file[$i]\" title=\"Zoomed Score distribution plot image\">";
			}
			$blurb = "<BR><P>This plot shows the estimation of the distribution of correct and incorrect spectra at different scores.  The two distributions should be seperate, with the incorrect distribution centered with the poorer scores and the correct at the more confident scores.";
			$title = "Correct / Incorrect";
		} 
		elsif($cgi->param("result_view") eq "combined")
		{
			my $imagefile = "http:\/\/ispider\.smith\.man\.ac\.uk\/FDRAnalysis\/tmp\/" . $session->id . "_VennDiagram.png";

			if($session->param("combine") == 1)
			{
				my $combinedimage = "http:\/\/ispider\.smith\.man\.ac\.uk\/FDRAnalysis\/tmp\/" . $session->id . "_CombinedVennDiagram.png";
				$final_image_line .=  "<IMG SRC=\"$imagefile\" title=\"Venn diagram of Individual Identification results\"><P>This plot represents the overlap of peptide identifications after combining the results using the FDRScore* method at FDR** " . $session->param("fdr_value") . " or better.<BR><P><FONT size=-1><I>*Jones et al. (2009) Proteomics 9: 1220-9.</I><P><I>**Kall et al. (2008) J. Proteome Res. 7:29-34</I><BR><IMG SRC=\"$combinedimage\" title=\"Venn diagram of Peptides identified from the combined analysis\"><BR>";
			}
			else
			{
				$final_image_line .=  "<IMG SRC=\"$imagefile\" title=\"Venn diagram of Individual Identification results\">";
			}
			$blurb = "<BR><P>This plot represents the overlap of peptide identification from the different search engines at FDR* " . $session->param("fdr_value") . " or better.";
			$blurb .= "<BR><font size=-2><I><b>*</b> Elias and Gygi (2007) Nat. Methods 4:207-214</i></FONT>";
			$title = "Overlap of Peptides from different search engines";
		}

		#only if the user requested Nter do we display the nter option
		my $nterminal = " ";
		if($session->param("nterminal") == 1)
		{
			$nterminal = qq{<li><a href="http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl?result_view=nterdist">Nterminal Distribution</a></li>};
		}

		#  #delta
		#   for(my $i=0 ; $i<scalar(@FDR_image_file) ; $i++)
		#   {
		#   $FDR_image_file[$i] =~ s/\.png/DeltaMass\.png/;
		#   $final_image_line .=  "<IMG SRC=\"$FDR_image_file[$i]\" title=\"Delta mass plot image\">";
		#   }

		print  $cgi->div({id=>"main_result"},
				$cgi->div({id=>"content"},
				$cgi->br(),
				$cgi->h2($title),
 				$cgi->br(),
				$cgi->table(
					$cgi->Tr(
						$cgi->td({valign=>"top"},
						$cgi->div({id=>"nav"},
                                   		$cgi->table(
								$cgi->th($cgi->p("Select display format")),
                                          		$cgi->Tr( 
									$cgi->td(
										$cgi->ul(
                                							$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl?result_view=summary"},"Summary")),
											$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl?result_view=combined"},"Peptide Overlap")),
											$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl?result_view=rankplot"},"Rank plot")),
                               							$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl?result_view=deltamass"},"Delta Mass")),
                                							$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl?result_view=scoredist"},"Score Distribution")),
                                							$nterminal,
											$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDR_analysis_static.pl?result_view=zoomscore"},"Correct/Incorrect"))
										)
									)
								)
							)
						)
					),
                            	$cgi->td(
						$cgi->p($blurb)
						)
					),
					$cgi->Tr(
						$cgi->td($cgi->br()),
						$cgi->td($final_image_line)
                			)
				)
			)
		);
	}
	#backButton();
}

#DCW 101209 - not currently used
sub backButton()
{
	print  $cgi->div({id=>"main_result"},
		$cgi->div({id=>"content"},
			$cgi->br(),
			$cgi->start_form(-name=>'fdr_search',-action=>"http://www.ispider.manchester.ac.uk/FDRAnalysis/FDR_analysis_search.html", -method=>"post"),
			$cgi->submit(-value=>"New Analysis", id=>'textarea_border'),
			$cgi->input({-type=>"hidden", -name=>"submit_clicked"}),
			$cgi->input({-type=>"hidden", -name=>"result_view"}),
			$cgi->end_form(),
			$cgi->br()
		)
        );
}

sub DealWithFileUpload
{
	my $file = shift;
	my $search = shift;

	#print("dealing with file $file\n");
	
 	if(length($file)<1)
 	{
 		ErrorMsg("Sorry, could not understand the filename.\n");
 		exit(1);
 	}
 	elsif(length($file)>150)
 	{
 		ErrorMsg("Sorry the filename is too large\n");
 		exit(1);
 	} 
 	elsif($file =~ m/\// || $file =~ m/^\./)
 	{
 		ErrorMsg("Sorry the filename is not valid\n");
 		exit(1);
 	}

	#is the file too big?
	if($ENV{'CONTENT_LENGTH'} > 1024*1024*200)
	{
		ErrorMsg("Sorry, the file ($file) you have tried to upload is too big.\n");
		exit(1);
	}

	my $tmp = "/var/www/tmp/" . $search . "_" . $session->id . ".txt";
	my $cmd = "rm $tmp";
	system($cmd);

	my $exists = -e $file;
	my $size = -s $file;

	open(LOCAL,">$tmp") or die "unable to upload file $file";
 	while(<$file>)
 	{
 		print LOCAL $_;
	}
	close LOCAL;

	return $tmp;
}

sub footer()
{
	print $cgi->div({id=>"footer"},
	$cgi->div({id=>"path"},
	$cgi->p("Copyright &copy; 2008",$cgi->a({href=>"mailto:jennifer.siepen\@manchester.ac.uk"},"Jennifer Siepen"))));
}


sub navigation_bar()
{
	print $cgi->div({id=>"navigation"},
		$cgi->table({width=>"100%"},
 			$cgi->td($cgi->a({href=>"http://www.ispider.manchester.ac.uk/FDRAnalysis/FDR_analysis_home.html"},$cgi->img({width=>"60",height=>"50",alt=>"link to home",title=>"HOME",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/home.jpg"}))),
              	$cgi->td($cgi->a({href=>"http://www.ispider.manchester.ac.uk/FDRAnalysis/FDR_analysis_search.html"},$cgi->img({title=>"SEARCH", width=>"60",height=>"50",alt=>"link to search",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/search.jpg"}))),
			$cgi->td($cgi->a({href=>"FDR_analysis_static.pl"},$cgi->img({title=>"ANALYSE",width=>"60",height=>"50",alt=>"link to FDR analysis",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/analysis.jpg"}))),
			$cgi->td($cgi->a({href=>"http://www.ispider.manchester.ac.uk/FDRAnalysis/FDR_analysis_help.html"},$cgi->img({title=>"HELP", width=>"60",height=>"50",alt=>"link to help",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg"}))),
			$cgi->td($cgi->a({href=>"mailto:david.wedge\@manchester.ac.uk"},$cgi->img({title=>"CONTACT",width=>"60",height=>"50",alt=>"link to contact",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/contact.jpg"})))
		)
 	);
}

sub ErrorMsg()
{
	my $line = shift;

	start_html();
	print  $cgi->div({id=>"main_result"},
                $cgi->div({id=>"content"},
                	$cgi->br(),
			$cgi->p("Error:", $line),
                )
        );
}


sub start_html()
{
	my $refresh = shift;
	header($refresh);
	navigation_bar();
}

sub header()
{
	my $refresh = shift;
	if($refresh)
	{
		print $cgi->start_html(
			-head=>$cgi->meta({-http_equiv=>'refresh',-content=>"3;FDR_analysis_static.pl"}),
			-title=>'FDR Analysis',
			-style=>{'src'=>'/FDRAnalysis/fdranalysis.css','media'=>'screen, tv, projection','title'=>'Default'},
		);
	}
	else
	{
		print $cgi->start_html(
			-title=>'FDR Analysis',
			-style=>{'src'=>'/FDRAnalysis/fdranalysis.css','media'=>'screen, tv, projection','title'=>'Default'},
		);
	}
}


sub RunAnalysis
{
	my $param_file = shift;

	my $cmd = "/var/www/localhost/cgi-bin/FDRAnalysis/TestAllDecoySearchesInOne.pl ";

	my $rev;
	my $for;
	my $omssa_for;
	my $omssa_rev;
	my $tandem_for;
	my $tandem_rev;
	my $oversize = 1;

	#get the parameters
	open(FILE,"<$param_file") or die "RunAnalysis unable to open the params file $param_file\n";
	while(my $line = <FILE>)
	{
		my @split = split/\t/,$line;
		if($split[0] eq "analysis_type")
		{
			if($split[1] =~ m/on/)
			{
				$cmd .= "-n 1 ";
			}
		}
		if($split[0] eq "fdr_type")
		{
 			if($split[1] =~ m/jones/)
 			{
				$cmd .= "-q 1 ";
			}
		}
  		#deal with the mascot files
		elsif($split[0] eq "mascot_search_type")
		{
			if($split[1] =~ m/Mascot\s+decoy/)
			{
				$cmd .= "-d $rev ";
			}   
			elsif($split[1] =~ m/concatenated\s+forward/)
			{
				$cmd .= "-c $rev ";
			}
			elsif($split[1] =~ m/separate\s+forward/)
			{
				$cmd .= "-d $rev -f $for ";
			}
		}
		#deal with omssa
		elsif($split[0] eq "omssa_search_type")
		{
 			if($split[1] =~ m/concatenated\s+forward/)
			{
				$cmd .= " -r $omssa_rev ";
			}
			elsif($split[1] =~ m/separate\s+forward/)
			{
				$cmd .= " -p $omssa_for -i $omssa_rev ";
			}
		}
		#deal with tandem
		elsif($split[0] eq "tandem_search_type")
		{
			if($split[1] =~ m/concatenated\s+forward/)
			{
				$cmd .= " -x $tandem_rev ";
			}
			elsif($split[1] =~ m/separate\s+forward/)
			{
				$cmd .= " -u $tandem_for -v $tandem_rev ";
			}
		}
		elsif($split[0] =~ m/^for\_file/)
		{
			$for = $split[1];
			$for =~ s/\n//g;
		}
		elsif($split[0] =~ m/omssa\_for\_file/)
		{
			$omssa_for = $split[1];
			$omssa_for =~ s/\n//g;
		}
		elsif($split[0] =~ m/tandem\_for\_file/)
		{
			$tandem_for = $split[1];
			$tandem_for =~ s/\n//g;
		}
		elsif($split[0] =~ m/^rev\_file/) 
		{
			$rev = $split[1];
			$rev =~ s/\n//g;
		}
		elsif($split[0] =~ m/omssa\_rev\_file/)
		{
			$omssa_rev = $split[1];
			$omssa_rev =~ s/\n//g;
		}
		elsif($split[0] =~ m/tandem\_rev\_file/)
		{
			$tandem_rev = $split[1];
			$tandem_rev =~ s/\n//g;
		}
		elsif($split[0] =~ m/^image/)
		{
			my $image = $split[1];
			$image =~ s/\n//g;
			$cmd .= " -o $image ";
		}
		elsif($split[0] =~ m/omssa\_image/)
		{
			my $image = $split[1];
			$image =~ s/\n//g;
			$cmd .= " -j $image ";
		}
		elsif($split[0] =~ m/tandem\_image/)
		{
			my $image = $split[1];
			$image =~ s/\n//g;
			$cmd .= " -e $image ";
		} 
		elsif($split[0] =~ m/oversize/)
		{
			$oversize = $split[1];
		}
		elsif($split[0] =~ m/combine/) 
		{
			if($split[1] =~ /on/)
			{
				$oversize =~ s/\n//g;
				$cmd .= " -z $oversize";
				$cmd .= " -D $oversize"; #DCW
			}
		}
		elsif($split[0] =~ m/max\_expect/)
		{
			my $value = $split[1];
			$value =~ s/\n//g;
			$cmd .= " -k $value ";
		}
		elsif($split[0] =~ m/fdr\_value/)
		{
			my $value = $split[1];
			$value =~ s/\n//g;
			$cmd .= " -w $value ";
		}
		elsif($split[0] =~ m/rev\_tag/)
		{
			my $value = $split[1];
			$value =~ s/\n//g;
			if($value)
			{
				$cmd .= " -t $value ";
			}
		}
	}
	close FILE;

	$cmd .= "  -s /var/www/tmp/" . $session->id . " ";

	$cmd .= " -I"; #DCW - produce images
	#DCW - produces combined search info
	#if($cgi->param("combine"))
	#if($combine_param)
	if($session->param("combine") == 1)
	{
		$cmd .= " -z 1";
	}

	#$cmd .= " &"; #this MUST be commented out if using QSUB!!!
	#system($cmd);
	#print "cmd is $cmd\n";
	my $shell = "/tmp/RunFDRAnalysis_" . $session->id . ".sh";

	open(SHELL,">$shell") or die "unable to create a shell file, $shell\n";
	print SHELL "#!/bin/sh \n";
	print SHELL $cmd;
	print SHELL "\n\n";

	my $sys_call = $QSUB." -q www -o /var/www/tmp -e /var/www/tmp " . $shell . " &";
	#print "syscall is $sys_call\n";
	system($sys_call);

	return 1;
}



