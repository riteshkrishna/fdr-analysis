#!/usr/bin/perl

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

#the display
#header();
#navigation_bar();
main_content($cgi->param("content"));
footer();



sub main_content()
{
my $page = shift;


 #HOME
 if(!$page ||$page eq "home")
 {
 start_html();
 print	$cgi->div({id=>"main_result"},
		$cgi->div({id=>"content"},
 				$cgi->h1("FDRAnalysis"),
             			$cgi->p("These pages should allow you to an FDR analysis on your data.  It can accept Mascot, X!Tandem and OMSSA results.  It also provides the option to run an X!Tandem and/or OMSSA search"),
            			$cgi->p("Please use the menu buttons above to navigate through the analysis")
 	      	  		)
			);
 
 }

 #upload
 elsif($page eq "upload")
 {
 start_html();

 #determine the file upload required - it depends upon the file(s) chosen
 my $file_upload_text;
  if($cgi->param("mascot_search_type") =~ m/Mascot\s+decoy/)
  { 
  $file_upload_text = qq{$cgi->p("please enter the Mascot results file (from the Mascot-decoy search)"),$cgi->input({type=>"file", name=>"mas_result"}),};

  }
  elsif($cgi->param("mascot_search_type") =~ m/concatanated\s+forward\/reverse/)
  {
   $file_upload_text = "$cgi->p(\"please enter the Mascot results file (from the concatanated forward/reverse Mascot search)\"),$cgi->input({type=>\"file\", name=>\"fr_result\"})";
  }
  elsif($cgi->param("mascot_search_type") =~ m/seperate\s+forward\/reverse/)
  {
  $file_upload_text = qq{<p>please enter the Mascot results file (from the 'forward' search)</p><input type="file" name="for_result">};
  $file_upload_text .= qq{<p>please enter the Mascot results file (from the 'reverse' search</p><input type="file" name="rev_result">};
  }

 print  $cgi->div({id=>"main_result"},
                $cgi->div({id=>"content"},
                                $cgi->br(),
				$cgi->h3("Have you set the ",$cgi->a({href=>"FDRanalysis.pl?content=search"},"parameters?")),
  				$cgi->br(),
				$cgi->p("Please Upload your file(s)"),
 				$file_upload_text,	
          			$cgi->br(),
				$cgi->submit(-value=>"Perform Analysis"),		
		)
	);
 } 


 #search
 elsif($page eq "search")
 {
 start_html();

 $session->clear();
 $session->flush();


 my $cmd = "rm /var/www/localhost/htdocs/FDRAnalysis/tmp/*fdranalysis_" . $cgi->cookie("CGISESSID") . "*.png";
 system($cmd); 

 print  $cgi->div({id=>"main_result"},
                $cgi->div({id=>"content"},
			$cgi->start_form(-name=>'fdr_search',-action=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse", -method=>"post"),
			$cgi->br(),
 			$cgi->h2("Parameters for the FDR Analysis"),
			$cgi->br(),
 			$cgi->div({id=>"innerboxtable"},
  			$cgi->div({id=>"innerboxheader"},
                        $cgi->p("Analysis Details")),
                        $cgi->div({id=>"innerbox",align=>"left"},
 			$cgi->br(),
                        $cgi->input({-type=>'checkbox',-name=>'nter',id=>'textarea_border'},'   Include N terminal analysis'),
			$cgi->br(),
                        $cgi->input({-type=>'text',-name=>'rev_tag',value=>"REV_",id=>'textarea_border',size=>"5"},'   Tag used in Decoy search'),
			$cgi->br(),
			$cgi->input({-type=>'checkbox',-name=>'combine',id=>'textarea_border'},'   Combine results from different search engines'),
                        $cgi->br(),
			$cgi->input({-type=>'text',-name=>'max_expect',value=>"0.05",id=>'textarea_border',size=>"5"},'   Maximum expectation value to use for the Rank plot'),
                        $cgi->br(),
                        $cgi->input({-type=>'text',-name=>'fdr_value',value=>"0.05",id=>'textarea_border',size=>"5"},'   FDR rate'),
                        $cgi->input({-type=>'hidden',-name=>'analysis_type'}),
 			$cgi->br(),
                        $cgi->br(),	
			),
			$cgi->br(),	
			$cgi->br(),
			$cgi->div({id=>"innerboxheader"},
 			$cgi->p("Mascot File(s) Upload")),
                        $cgi->div({id=>"innerbox",align=>"left"},
			$cgi->br(),
			$cgi->popup_menu({type=>"popup_menu",id=>'textarea_border', name=>"mascot_search_type",width=>"50",values=>[" ","Mascot decoy","concatanated forward/reverse","seperate forward/reverse"]}),
 			"   Select search type",
			$cgi->br(),
			$cgi->input({type=>"file", id=>'textarea_border', name=>"rev_file",size=>"50",value=>"/fs/msct/data/20080704/F291465405.dat"}," Decoy/Concatanated/Reverse"),
                        $cgi->br(),
			$cgi->input({type=>"file", id=>'textarea_border', name=>"for_file",size=>"50"}," forward"),
			$cgi->br(),
                        $cgi->br(),
			),
                        $cgi->div({id=>"innerboxheader"},
                        $cgi->p("Omssa File(s) Upload")),
                        $cgi->div({id=>"innerbox",align=>"left"},
                        $cgi->br(),
  			$cgi->popup_menu({type=>"popup_menu",id=>'textarea_border', name=>"omssa_search_type",width=>"50",values=>[" ","concatanated forward/reverse","seperate forward/reverse"]}),
			"   Select search type",
                        $cgi->br(),
                        $cgi->input({type=>"file", id=>'textarea_border', name=>"omssa_rev_file",size=>"50",value=>""}," Concatanated/Reverse"),
                        $cgi->br(),
                        $cgi->input({type=>"file", id=>'textarea_border', name=>"omssa_for_file",size=>"50"}," forward"),
                        $cgi->br(),
			$cgi->br(),
                        ),
                        $cgi->div({id=>"innerboxheader"},
                        $cgi->p("X!Tandem File(s) Upload")),
                        $cgi->div({id=>"innerbox",align=>"left"},
                        $cgi->br(),
                        $cgi->popup_menu({type=>"popup_menu",id=>'textarea_border', name=>"tandem_search_type",width=>"50",values=>[" ","concatanated forward/reverse","seperate forward/reverse"]}),
 			"   Select search type",
                        $cgi->br(),
                        $cgi->input({type=>"file", id=>'textarea_border', name=>"tandem_rev_file",size=>"50",value=>""}," Concatanated/Reverse"),
                        $cgi->br(),
                        $cgi->input({type=>"file", id=>'textarea_border', name=>"tandem_for_file",size=>"50"}," forward"),
                        $cgi->br(),
			$cgi->br(),
                        $cgi->input({-type=>"hidden",-name=>"content",value=>"analyse"}),
                        ),
                        $cgi->div({id=>"innerboxheader"},
                        $cgi->p("")),
                        $cgi->br(),
 			$cgi->submit(-value=>"Perform Analysis"),
			$cgi->br(),
			$cgi->br(),
			)
                )
        );


 }

 #analyse
 elsif($page eq "analyse")
 {

 my $results_ready = 1;
 my $refresh_browser;
 my $true = 0;
 my @FDR_image_file;
  # mascot?
  if($cgi->param("mascot_search_type") =~ m/[Mcs]/ || $session->param("mascot") == 1)
  {
  $FDR_image_file[0] = "/var/www/localhost/htdocs/FDRAnalysis/tmp/mascot_fdranalysis_" . $session->id . ".png";
  my $tmp = $FDR_image_file[0];
  $tmp =~ s/\.png/GygiRank\.png/;
   if(!-e $tmp)
   {
   $results_ready = 0;
   $refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse\">";
   }
   $session->param('mascot', 1);
   $session->flush();
  }
  if($cgi->param("omssa_search_type")  =~ m/^[cs]/ || $session->param("omssa") == 1)
  {
  $FDR_image_file[1] = "/var/www/localhost/htdocs/FDRAnalysis/tmp/omssa_fdranalysis_" . $session->id . ".png";
  my $tmp = $FDR_image_file[1];
  $tmp =~ s/\.png/GygiRank\.png/;
   if(!-e $tmp)
   {
   $results_ready = 0;
   $refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse\">";
   }
  $session->param('omssa', 1);
  $session->flush();
  }
  else
  {
  $session->param('omssa', 0);
  $session->flush();
  }
  if($cgi->param("combine") || $session->param("combine") == 1)
  {
my $tmp = "/var/www/tmp/" . $session->id . "combined_peptides.out";
my $tmp2 = "/var/www/tmp/" . $session->id . "_summary.txt";
   if(!-e $tmp || -z $tmp || !-e $tmp2)
   {
$results_ready = 0;
   $refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse\">";

  }

  $session->param('combine', 1);
  $session->flush();
  }
  else
  {
  $session->param('combine', 0);
  $session->flush();
  }
  if($cgi->param("tandem_search_type") =~ m/^[cs]/ || $session->param("tandem") == 1)
  {
  $FDR_image_file[2] = "/var/www/localhost/htdocs/FDRAnalysis/tmp/tandem_fdranalysis_" . $session->id . ".png"; 
  my $tmp = $FDR_image_file[2];
   $tmp =~ s/\.png/GygiRank\.png/;
   if(!-e $tmp)
   {
   $results_ready = 0;
   $refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse\">";
   }
  $session->param("tandem",1);
  $session->flush();
  }
  else
  {
  $session->param("tandem",0);
  $session->flush();
  }

  if(!$session->param("tandem") && !$session->param("omssa") && !$session->param("mascot"))
  {
  ErrorMsg("You need to run a search before you can view the results.\n");
  exit(1);
  }


  if(!$session->param("paramsset"))
  {
$session->param("paramsset",1);
$results_ready = 0;
  $refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse\">";

  my $params = "/var/www/tmp/fdranalysis_param_" . $session->id . ".txt";
  my $cmd = "rm $params";
  system($cmd);
  ##Mascot
   if($cgi->param("mascot_search_type") =~ m/[Mcs]/)
   {
   #was a file uploaded? 
    if(!$cgi->param("rev_file"))
    {
    ErrorMsg("You have not uploaded a Mascot file to analyse\n");
    exit(1);
    }
    elsif($cgi->param('mascot_search_type') =~ m/seperate\s+forward\/reverse/ && !$cgi->param("for_file"))
    {
    ErrorMsg("You requested the analysis of two seperate files yet you have only uploaded one file.\n");
    exit(1);
    }
   #deal with the rev_file
   $cgi->param('rev_file', DealWithFileUpload($cgi->param("rev_file"),"mascot_rev"));
    if($cgi->param("for_file"))
    {
    $cgi->param('for_file', DealWithFileUpload($cgi->param("for_file"),"mascot_for"));
    }
   #write the Mascot parameters to the param file
   open(PARAM,">$params") or die "unable to open the file $params to write to\n";
   print PARAM "rev_file\t" . $cgi->param("rev_file") . "\n";
   print PARAM "for_file\t" . $cgi->param("for_file") . "\n";
   print PARAM "mascot_search_type\t" . $cgi->param("mascot_search_type") . "\n";
   print PARAM "image\t" . $FDR_image_file[0] . "\n";
   close PARAM;

   }
   elsif($cgi->param("rev_file") && $cgi->param("mascot_search_type") !~ m/[Mcs]/)
   {
   ErrorMsg("You have not chosen the type of Mascot search\n");
   exit(1);
   }

   ##Omssa
   if($cgi->param("omssa_search_type")  =~ m/^[cs]/)
   {
   #was a file uploaded?
    if(!$cgi->param("omssa_rev_file"))
    {
    ErrorMsg("You have not uploaded an OMSSA file to analyse \n");
    exit(1);
    }
    elsif($cgi->param('omssa_search_type') =~ m/seperate\s+forward\/reverse/ && !$cgi->param("for_file"))
    {
    ErrorMsg("You requested the analysis of two seperate files yet you have only uploaded one omssa file.\n");
    exit(1);
    }

   $cgi->param('omssa_rev_file', DealWithFileUpload($cgi->param("omssa_rev_file"),"omssa_rev"));
    if($cgi->param("omssa_for_file"))
    { 
    $cgi->param('omssa_for_file', DealWithFileUpload($cgi->param("omssa_for_file"),"omssa_for"));
    }
   #write the Omssa parameters to the param file
   open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
   print PARAM "omssa_rev_file\t" . $cgi->param("omssa_rev_file") . "\n";
   print PARAM "omssa_for_file\t" . $cgi->param("omssa_for_file") . "\n";
   print PARAM "omssa_search_type\t" . $cgi->param("omssa_search_type") . "\n";
   print PARAM "omssa_image\t" . $FDR_image_file[1] . "\n";
   close PARAM;
   }

   if($cgi->param("nter"))
   {
   $session->param('nterminal', 1);
   open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
   print PARAM "analysis_type\t" . $cgi->param("nter") . "\n";
   close PARAM;
   }
   else
   {
   $session->param('nterminal', 0);
   }
   if($cgi->param("combine"))
   {
   open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
   print PARAM "combine\t" . $cgi->param("combine") . "\n";
   close PARAM;
   }

   #add the max_expect
   open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
   print PARAM "max_expect\t" . $cgi->param("max_expect") . "\n";
   print PARAM "fdr_value\t" . $cgi->param("fdr_value") . "\n";
   close PARAM;
   $session->param('max_expect', $cgi->param("max_expect"));
   $session->param('fdr_value', $cgi->param("fdr_value"));
   ##X!Tandem
   if($cgi->param("tandem_search_type") =~ m/^[cs]/)
   {
   #was a file uploaded?
    if(!$cgi->param("tandem_rev_file"))
    {
    ErrorMsg("You have not uploaded an X!TANDEM file to analyse\n");
    exit(1);
    }
    elsif($cgi->param('tandem_search_type') =~ m/seperate\s+forward\/reverse/ && !$cgi->param("for_file"))
    {
    ErrorMsg("You requested the analysis of two seperate files yet you have only uploaded one tandem file.\n");
    exit(1);
    }
   $cgi->param('tandem_rev_file', DealWithFileUpload($cgi->param("tandem_rev_file"),"tandem_rev"));
    if($cgi->param("tandem_for_file"))
    {
    $cgi->param('tandem_for_file', DealWithFileUpload($cgi->param("tandem_for_file"),"tandem_for"));
    }
   #write the X!Tandem parameters to the param file
   open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
   print PARAM "tandem_rev_file\t" . $cgi->param("tandem_rev_file") . "\n";
   print PARAM "tandem_for_file\t" . $cgi->param("tandem_for_file") . "\n";
   print PARAM "tandem_search_type\t" . $cgi->param("tandem_search_type") . "\n";
   print PARAM "tandem_image\t" . $FDR_image_file[2] . "\n";
   close PARAM;
   } 


  if(!$cgi->param("max_expect"))
  {
  ErrorMsg("You must set a maximum expectation value to use\n");
  exit(1);
  }
  if(!$cgi->param("fdr_value"))
  {
  ErrorMsg("You must set a maximum FDR value to use\n");
  exit(1);
  }
  if($cgi->param("fdr_value")>0.5)
  {
  ErrorMsg("The maximum FDR value is 0.5.\n");
  exit(1);
  }


  $session->flush();
  #Run the analysis
  RunAnalysis($params);
  }#true!=0

start_html($refresh_browser);
 if($results_ready == 0)
 {
 print  $cgi->div({id=>"main_result"},
               $cgi->div({id=>"content"},
                               $cgi->br(),
			$cgi->p("The analysis is running, please be patient this shouldn't take too long."),
			$cgi->p("This page will automatically refresh, if it doesn't please click ",$cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse"},"here")),
 			$cgi->br(),
			$cgi->img({src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/Waiting.gif"}),
			$cgi->br(),
               )
       );

 } 
 else
 {
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
   #actually a blacnk image as the results are now in the blurb
#   $final_image_line = "<IMG SRC=\"http://ispider.smith.man.ac.uk/FDRAnalysis/summary.png\" title=\"blank image\" width=\"5%\">\n";
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
    for(my $i=0 ; $i<scalar(@FDR_image_file) ; $i++)                                                                                           {
    $FDR_image_file[$i] =~ s/\.png/NterDist\.png/;
    $final_image_line .=  "<IMG SRC=\"$FDR_image_file[$i]\" title=\"Start position  plot image\">"; 
    }
   $blurb = "<BR><P>This plot is part of the N-terminal proteome analysis pipeline.  This analysis specifically targets the N-terminal peptides and as a result the forward hits should be located predominantly at positions 1 and 2 and the reverse scattered randomly amongst all other positions. (Note max spectra shown on y-axis is 100 for a clearer view).";
   $title = "N-terminal distribution";
   }
   elsif($cgi->param("result_view") eq "scoredist")
   {
    for(my $i=0 ; $i<scalar(@FDR_image_file) ; $i++)                                                                                           {
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
   $final_image_line .=  "<IMG SRC=\"$imagefile\" title=\"Venn diagram of combined results\">";
   $blurb = "<BR><P>This plot represents the overlap of peptide identification from the different search engines at FDR* " . $session->param("fdr_value") . " or better.";
   $blurb .= "<BR><font size=-2><I><b>*</b> Elias and Gygi (2007) Nat. Methods 4:207-214</i></FONT>";
   $title = "Overlap of Peptides from different search engines";

   }

  #only if the user requested Nter do we display the nter option
  my $nterminal = " ";
  if($session->param("nterminal") == 1)
  {
   $nterminal = qq{<li><a href="http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse;result_view=nterdist">Nterminal Distribution</a></li>};
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
                                		 $cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse;result_view=summary"},"Summary")),
                                			$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse;result_view=rankplot"},"Rank plot")),
                               				 $cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse;result_view=deltamass"},"Delta Mass")),
                                			$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse;result_view=scoredist"},"Score Distribution")),
                                			$nterminal,
							$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse;result_view=zoomscore"},"Correct/Incorrect")),
							$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis.pl?content=analyse;result_view=combined"},"Peptide Overlap"))
							)

							))))),
                                $cgi->td(
				$cgi->p($blurb))),
				$cgi->Tr(
				$cgi->td($cgi->br()),
				$cgi->td($final_image_line)	
	
                ))));
        


  }

 }

 #HELP
 elsif($page eq "help")
 {

start_html();



print  $cgi->div({id=>"main_result"},
                $cgi->div({id=>"content"},
                                $cgi->br(),
				$cgi->div({id=>"innerboxheader"},
				$cgi->p("Help Page - Contents")),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->a({href=>"#intro"},"Introduction"),$cgi->br(),
				$cgi->a({href=>"#analysis"},"Search Form"),$cgi->br(),$cgi->br(),
				$cgi->a({href=>"#analysis"},"&nbsp;&nbsp;&nbsp;Analysis"),$cgi->br(),
				$cgi->a({href=>"#mascot_files"},"&nbsp;&nbsp;&nbsp;Mascot Files"),$cgi->br(),
				$cgi->a({href=>"#omssa_files"},"&nbsp;&nbsp;&nbsp;OMSSA/X!Tandem Files"),$cgi->br(),$cgi->br(),
				$cgi->a({href=>"#results_summary"},"Results"),$cgi->br(),
				$cgi->a({href=>"#results_rank"},"&nbsp;&nbsp;&nbsp;Results - Rank Plot"),$cgi->br(),
				$cgi->a({href=>"#results_delta"},"&nbsp;&nbsp;&nbsp;Results - Delta Mass Plot"),$cgi->br(),
				$cgi->a({href=>"#results_score"},"&nbsp;&nbsp;&nbsp;Results - Score Distribution"),$cgi->br(),
				$cgi->a({href=>"#results_correct"},"&nbsp;&nbsp;&nbsp;Results - Correct/Incorrect"),$cgi->br(),
				$cgi->a({href=>"#results_nter"},"&nbsp;&nbsp;&nbsp;Results - Nterminal Distribution"),$cgi->br(),
				),
				$cgi->br(),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"intro"},$cgi->p("FDR Analysis Help Page"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p("This software enables the user to analyse the results of a Mascot, Omssa and /or X!Tandem database search, particuarly an in depth look at the results from a decoy search")),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"analysis"},$cgi->p("Analysis Details"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p($cgi->b("N terminal analysis"),"if checked will include an N-terminal analysis.  This included looking at the fraction of identifications that are the N-terminal most peptide of the protein"),
				$cgi->p($cgi->b("Decoy tag"),"We are dealing with decoy database searching.  Except for a Mascot decoy search the tag used to mark the protein accession as a decoy protein is required.  For example if the target database has accession >P12345 and the equivalent decoy protein is >REVERSE_P12345 then the tag is 'REVERSE_'"),
				$cgi->p($cgi->b("Combine Results"),"This employs the algorithm described in Jones et al. Proteomics (2009) 9: 1220-9 to combine the results of two or more searches from the different search engines."),
  				),	
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"mascot_files"},$cgi->p("Mascot File(s) Upload"))),
				$cgi->div({id=>"innerbox",align=>"left"},
                                $cgi->p($cgi->b("Search type")," There are 3 types of decoy search a concatanated forward/decoy whereby the target and decoy databases were joined and one search was performed.  a Mascot decoy search whereby the 'decoy' checkbox was selected in the Mascot search and finally two searches, one on the target database and a second on the decoy database."),
				$cgi->p("The first file upload box is the main upload box and is for any file except the search run soley on the target database, which the second box will contain if applicable."),
				),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"omssa_files"},$cgi->p("OMSSA/X!Tandem File(s) Upload"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p($cgi->b("Search type")," There are 2 types of decoy search for both search engines, a concatanated forward/decoy whereby the target and decoy databases were joined and one search was performed or where two searches were performed one on the target database and a second on the decoy database."), 
				$cgi->p("The first file upload box is the main upload box and is for any file except the search run soley on the target database, which the second box will contain if applicable."),
				),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_summary"},$cgi->p("Results Pages - summary"))),
			        $cgi->div({id=>"innerbox",align=>"left"},
			        $cgi->p("A simple table displays the peptide frequencies identified from the search at FDR 0.05.  There are two different FDR calculation methods .... see ....."),
				),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_rank"},$cgi->p("Results Pages - Rank Plot"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p("Rank plots are shown for the different search engines used in the query form.  Note that X!Tandem provides only the top ranked peptide identifications and OMSSA doesn't always report all ranks.  The description of the plot can be seen to the left of the plot."),
				),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_delta"},$cgi->p("Results Pages - Delta Mass Plot"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p(" Delta Mass plots are shown for the different search engines used in the query form.  Note that X!Tandem provides only the top ranked peptide identifications and OMSSA doesn't always report all ranks.  The description of the plot can be seen to the left of the plot."),
                                ),
                                $cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_score"},$cgi->p("Results Pages - Score Distribution"))),
				$cgi->div({id=>"innerbox",align=>"left"},
			        $cgi->p("Score distribution plots are shown for the different search engines used in the query form. "),
	                        ),	
                                $cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_correct"},$cgi->p("Results Pages - Correct/Incorrect"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p("Correct/Incorrect plots are shown for the different search engines used in the query form. They should give an indication of the score distributions of the forward and reverse hits."),
                                ),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_nter"},$cgi->p("Results Pages - Nterminal Distribution"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p("Nterminal distributions are shown for the different search engines used in the query form only if the Nter analysis opetion was checked. They show the start position of the identified peptides from the protein sequences."),
				),
                )
        );
 }

}


sub DealWithFileUpload
{
my $file = shift;
my $search = shift;

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
open(LOCAL,">$tmp") or die "$!";
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
print $cgi->div({id=>'header'},$cgi->a({href=>"FDRanalysis.pl?content=home"}), 

  	$cgi->ul({id=>"navigation"},
 			$cgi->li($cgi->a({href=>"FDRanalysis.pl?content=home"},$cgi->img({width=>"60",height=>"50",alt=>"link to home",title=>"HOME",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/home.jpg"}))),
                        $cgi->li($cgi->a({href=>"FDRanalysis.pl?content=search"},$cgi->img({title=>"SEARCH", width=>"60",height=>"50",alt=>"link to search",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/search.jpg"}))),
			$cgi->li($cgi->a({href=>"FDRanalysis.pl?content=analyse"},$cgi->img({title=>"ANALYSE",width=>"60",height=>"50",alt=>"link to FDR analysis",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/analysis.jpg"}))),
			$cgi->li($cgi->a({href=>"FDRanalysis.pl?content=help"},$cgi->img({title=>"HELP", width=>"60",height=>"50",alt=>"link to help",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg"}))),
			$cgi->li($cgi->a({href=>"mailto:jennifer.siepen\@manchester.ac.uk"},$cgi->img({title=>"CONTACT",width=>"60",height=>"50",alt=>"link to contact",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/contact.jpg"})))

 			));
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
  -head=>$cgi->meta({-http_equiv=>'refresh',-content=>"3;FDRanalysis.pl?content=analyse"}),
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

sub news()
{

#data source name
my $username = "ispider_web";
my $password = "1sp1der";
my $db = "ispiderRegistry";
my $dsn = "DBI:mysql:" . $db . ":node12";
#connect to the database
my $dbh = DBI->connect($dsn, $username,$password, {RaiseError =>1});

my @date;
my @text;
my $count = 0;

my $sth = $dbh->prepare(qq{select date,text from news order by id desc limit 2});
$sth->execute();
 while(my @array = $sth->fetchrow_array ())
 {
 $date[$count] = $array[0];
 $text[$count] = $array[1];
 $count++;
 }

$sth->finish;

print $cgi->div({class=>"main_right_box"},$cgi->div({id=>"nav"},$cgi->h3("News"),$cgi->p("$date[0]",$cgi->br(),"$text[0]"),$cgi->p("$date[1]",$cgi->br(),"$text[1]"),$cgi->p($cgi->a({href=>"ispider.pl?content=news"},"More news..."))));

}

sub services()
{
 print $cgi->div({class=>"main_right_box"},
  $cgi->div({id=>"nav"},
  $cgi->h3("Services Available"),
  ));

}

sub RunAnalysis
{
my $param_file = shift;

my $cmd = "/var/www/localhost/cgi-bin/FDRAnalysis/TestMascotAllDecoySearchesInOne.pl ";

my $rev;
my $for;
my $omssa_for;
my $omssa_rev;
my $tandem_for;
my $tandem_rev;

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
   elsif($split[1] =~ m/concatanated\s+forward/)
   {
   $cmd .= "-c $rev ";
   }
   elsif($split[1] =~ m/seperate\s+forward/)
   {
   $cmd .= "-d $rev -f $for ";
   }
  }
  #deal with omssa
  elsif($split[0] eq "omssa_search_type")
  {
   if($split[1] =~ m/concatanated\s+forward/)
   {
   $cmd .= " -r $omssa_rev ";
   }
   elsif($split[1] =~ m/seperate\s+forward/)
   {
   $cmd .= " -p $omssa_for -i $omssa_rev ";
   }
  }
  #deal with tandem
  elsif($split[0] eq "tandem_search_type")
  {
   if($split[1] =~ m/concatanated\s+forward/)
   {
   $cmd .= " -x $tandem_rev ";
   }
   elsif($split[1] =~ m/seperate\s+forward/)
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
  elsif($split[0] =~ m/combine/) 
  {
   if($split[1] =~ /on/)
   {
   $cmd .= " -z "; 
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
 }
close FILE;

$cmd .= "  -s /var/www/tmp/" . $session->id . " ";


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



