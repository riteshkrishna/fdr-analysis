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
             			$cgi->p("This application enables an analysis of false discovery rates (FDR) of Mass spectrometry identification results.  Currently it can accept Mascot, X!Tandem and OMSSA results and can analyse all the searches separately but can also combine the search from different search engines, enhancing confidence in the results."),
            			$cgi->p("Please use the menu buttons above to run the analysis")
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
  elsif($cgi->param("mascot_search_type") =~ m/concatenated\s+forward\/reverse/)
  {
   $file_upload_text = "$cgi->p(\"please enter the Mascot results file (from the concatenated forward/reverse Mascot search)\"),$cgi->input({type=>\"file\", name=>\"fr_result\"})";
  }
  elsif($cgi->param("mascot_search_type") =~ m/separate\s+forward\/reverse/)
  {
  $file_upload_text = qq{<p>please enter the Mascot results file (from the 'forward' search)</p><input type="file" name="for_result">};
  $file_upload_text .= qq{<p>please enter the Mascot results file (from the 'reverse' search</p><input type="file" name="rev_result">};
  }

 print  $cgi->div({id=>"main_result"},
                $cgi->div({id=>"content"},
                                $cgi->br(),
				$cgi->h3("Have you set the ",$cgi->a({href=>"FDRanalysis_NEW.pl?content=search"},"parameters?")),
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
			#$cgi->start_form(-name=>'fdr_search', -method=>"post"),
			$cgi->start_form(-name=>'fdr_search',-action=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse", -method=>"post"),
			$cgi->br(),
 			$cgi->h2("Parameters for the FDR Analysis"),
			$cgi->br(),
 			$cgi->div({id=>"innerboxtable"},
  			$cgi->div({id=>"innerboxheader"},
                        $cgi->p("Analysis Details")),
                        $cgi->div({id=>"innerbox",align=>"left"},
 			$cgi->br(),
                        $cgi->input({-type=>'checkbox',-name=>'nter',id=>'textarea_border'},'   Include N terminal analysis'),
			$cgi->img({height => "20", alt=>"N-terminal help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=N_terminal_analysis','helpwindow','width=450,height=250');"}),
			$cgi->br(),
                        $cgi->input({-type=>'text',-name=>'rev_tag',value=>"REV_",id=>'textarea_border',size=>"5"},'   Tag used in Decoy search'),
			$cgi->img({height => "20", alt=>"Decoy tag help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=decoy_tag','helpwindow','width=450,height=250');"}),
			$cgi->br(),
			$cgi->input({-type=>'text',-name=>'max_expect',value=>"0.05",id=>'textarea_border',size=>"5"},'   Maximum expectation value to use for the Rank plot'),
			$cgi->img({height => "20", alt=>"Expectation value help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=expectation_value','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
                        $cgi->input({-type=>'text',-name=>'fdr_value',value=>"0.05",id=>'textarea_border',size=>"5"},'   FDR rate'),
                        $cgi->input({-type=>'hidden',-name=>'analysis_type'}),
			$cgi->img({height => "20", alt=>"FDR help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=FDR','helpwindow','width=450,height=250');"}),
 			$cgi->br(),
			$cgi->br(),
                        $cgi->input({-type=>'checkbox',-name=>'combine',id=>'textarea_border'},'   Combine results from different search engines'),
			$cgi->img({height => "20", alt=>"Combine results help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=combine_results','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
                        $cgi->input({-type=>'text',-name=>'oversize',id=>'textarea_border',value=>'1',size=>'2'},'  If using an oversized decoy enter size increase here (specific to combined analysis)'),
                        $cgi->br(),
			$cgi->br(),
                        $cgi->input({-type=>'text',-name=>'missed_cleavages',value=>"1",id=>'textarea_border',size=>"2"},'   Number of missed cleavages'),
			$cgi->img({height => "20", alt=>"Missed cleavages help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=missed_cleavages','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
                        $cgi->br(),
			"Parent mass type:  ",$cgi->radio_group({type=>"radio_group",name=>"parent_mass_type", values=>["monoisotopic","average"], default=>"monoisotopic"}),
			$cgi->img({height => "20", alt=>"Parent mass type help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=parent_mass_type','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
			"Parent tolerance ± ",$cgi->input({-type=>'text',-name=>'parent_tolerance',value=>"1.0",id=>'textarea_border',size=>"2",align=>"right"}),
			$cgi->popup_menu({type=>"popup_menu",id=>'textarea_border', name=>"parent_units",width=>"50",values=>["Da","ppm"]}),
                        $cgi->br(),
			"Fragment mass type:  ",$cgi->radio_group({type=>"radio_group",name=>"fragment_mass_type", values=>["monoisotopic","average"], default=>"monoisotopic"}),
			$cgi->img({height => "20", alt=>"Fragment mass type help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=fragment_mass_type','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
			"Fragment tolerance ± ",$cgi->input({-type=>'text',-name=>'fragment_tolerance',value=>"0.5",id=>'textarea_border',size=>"2",align=>"right"}),
			$cgi->popup_menu({type=>"popup_menu",id=>'textarea_border', name=>"fragment_units",width=>"50",values=>["Da","ppm"]}),
                        $cgi->br(),
                        $cgi->br(),
			$cgi->popup_menu({type=>"popup_menu",id=>'textarea_border', name=>"enzyme",width=>"50",values=>["Trypsin","Arg-C","Asp-N","Asp-N_ambic","Chymotrypsin","CNBr","CNBr+Trypsin","Formic_acid","Lys-C","Lys-C/P","PepsinA","Tryp-CNBr","TrypChymo","Trypsin/P","V8-DE","V8-E","semiTrypsin","LysC+AspN","None","No_cleavage"]}),
 			"   Select enzyme",
			$cgi->img({height => "20", alt=>"Select enzyme help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=enzyme','helpwindow','width=450,height=250');"}),
			$cgi->br(),
                     $cgi->br(),
			$cgi->popup_menu({id=>'textarea_border', name=>"fixed_mods",width=>"50",selectmode =>"multiple",multiple=>"5",values=>["Acetyl(K)","Acetyl(N-term)","Acetyl(ProteinN-term)","Amidated(C-term)","Amidated(ProteinC-term)","Ammonia-loss(N-termC)","Biotin(K)","Biotin(N-term)","Carbamidomethyl(C)","Carbamyl(K)","Carbamyl(N-term)","Carboxymethyl(C)","Cation:Na(C-term)","Cation:Na(DE)","Deamidated(NQ)","Dehydrated(N-termC)","Dehydro(C)","Dioxidation(M)","Ethanolyl(C)","ExacTagAmine(K)","ExacTagThiol(C)","Formyl(N-term)","Formyl(ProteinN-term)","Gln-&gt;pyro-Glu(N-termQ)","Glu-&gt;pyro-Glu(N-termE)","Guanidinyl(K)","ICAT-C(C)","ICAT-C:13C(9)(C)","ICPL(K)","ICPL(ProteinN-term)","ICPL:13C(6)(K)","ICPL:13C(6)(ProteinN-term)","ICPL:2H(4)(K)","ICPL:2H(4)(ProteinN-term)","iTRAQ4plex(K)","iTRAQ4plex(N-term)","iTRAQ4plex(Y)","iTRAQ8plex(K)","iTRAQ8plex(N-term)","iTRAQ8plex(Y)","Label:18O(1)(C-term)","Label:18O(2)(C-term)","Met-&gt;Hse(C-termM)","Met-&gt;Hsl(C-termM)","Methyl(C-term)","Methyl(DE)","Methylthio(C)","NIPCAM(C)","Oxidation(HW)","Oxidation(M)","Phospho(ST)","Phospho(Y)","Propionamide(C)","Pyridylethyl(C)","Pyro-carbamidomethyl(N-termC)","Sulfo(S)","Sulfo(T)","Sulfo(Y)","TMT(K)","TMT(N-term)","TMT2plex(K)","TMT2plex(N-term)","TMT6plex(K)","TMT6plex(N-term)"]}),
 			"   Select fixed mods",
			$cgi->img({height => "20", alt=>"Select fixed mods help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=fixed_mods','helpwindow','width=450,height=250');"}),
			$cgi->br(),
			$cgi->br(),
			$cgi->popup_menu({id=>'textarea_border', name=>"variable_mods",width=>"50",selectmode =>"multiple",multiple=>"5",values=>["Acetyl(K)","Acetyl(N-term)","Acetyl(ProteinN-term)","Amidated(C-term)","Amidated(ProteinC-term)","Ammonia-loss(N-termC)","Biotin(K)","Biotin(N-term)","Carbamidomethyl(C)","Carbamyl(K)","Carbamyl(N-term)","Carboxymethyl(C)","Cation:Na(C-term)","Cation:Na(DE)","Deamidated(NQ)","Dehydrated(N-termC)","Dehydro(C)","Dioxidation(M)","Ethanolyl(C)","ExacTagAmine(K)","ExacTagThiol(C)","Formyl(N-term)","Formyl(ProteinN-term)","Gln-&gt;pyro-Glu(N-termQ)","Glu-&gt;pyro-Glu(N-termE)","Guanidinyl(K)","ICAT-C(C)","ICAT-C:13C(9)(C)","ICPL(K)","ICPL(ProteinN-term)","ICPL:13C(6)(K)","ICPL:13C(6)(ProteinN-term)","ICPL:2H(4)(K)","ICPL:2H(4)(ProteinN-term)","iTRAQ4plex(K)","iTRAQ4plex(N-term)","iTRAQ4plex(Y)","iTRAQ8plex(K)","iTRAQ8plex(N-term)","iTRAQ8plex(Y)","Label:18O(1)(C-term)","Label:18O(2)(C-term)","Met-&gt;Hse(C-termM)","Met-&gt;Hsl(C-termM)","Methyl(C-term)","Methyl(DE)","Methylthio(C)","NIPCAM(C)","Oxidation(HW)","Oxidation(M)","Phospho(ST)","Phospho(Y)","Propionamide(C)","Pyridylethyl(C)","Pyro-carbamidomethyl(N-termC)","Sulfo(S)","Sulfo(T)","Sulfo(Y)","TMT(K)","TMT(N-term)","TMT2plex(K)","TMT2plex(N-term)","TMT6plex(K)","TMT6plex(N-term)"]}), 			"   Select variable mods",
			$cgi->img({height => "20", alt=>"Select variable mods help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=variable_mods','helpwindow','width=450,height=250');"}),
			$cgi->br(),
			$cgi->br()
			),

			$cgi->br(),	
			$cgi->br(),
			$cgi->div({id=>"innerboxheader"},
 			$cgi->p("Mascot File(s) Upload")),
                        $cgi->div({id=>"innerbox",align=>"left"},
			$cgi->br(),
			$cgi->popup_menu({type=>"popup_menu",id=>'textarea_border', name=>"mascot_search_type",width=>"50",values=>[" ","Mascot decoy","concatenated forward/reverse","separate forward/reverse"]}),
 			"   Select search type",
			$cgi->img({height => "20", alt=>"Mascot search help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=mascot_search','helpwindow','width=450,height=250');"}),
			$cgi->br(),
			$cgi->input({type=>"file", id=>'textarea_border', name=>"rev_file",size=>"50"}," Decoy/Concatenated/Reverse"),
			$cgi->img({height => "20", alt=>"Mascot decoy files help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=mascot_decoy_files','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
			$cgi->input({type=>"file", id=>'textarea_border', name=>"for_file",size=>"50"}," Forward"),
			$cgi->img({height => "20", alt=>"Mascot forward files help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=mascot_forward_files','helpwindow','width=450,height=250');"}),
			$cgi->br(),
                        $cgi->br(),
			),
                        $cgi->div({id=>"innerboxheader"},
                        $cgi->p("Omssa File(s) Upload")),
                        $cgi->div({id=>"innerbox",align=>"left"},
                        $cgi->br(),
 			$cgi->input({type=>"text",id=>'textarea_border', name=>"omssa_software_version",size=>"5",value=>"2.14"}),
			"   OMSSA software version",
			$cgi->img({height => "20", alt=>"OMSSA version help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=OMSSA_version','helpwindow','width=450,height=250');"}),
                       $cgi->br(),
  			$cgi->popup_menu({type=>"popup_menu",id=>'textarea_border', name=>"omssa_search_type",width=>"50",values=>[" ","concatenated forward/reverse","separate forward/reverse"]}),
			"   Select search type",
			$cgi->img({height => "20", alt=>"OMSSA search help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=OMSSA_search','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
                        $cgi->input({type=>"file", id=>'textarea_border', name=>"omssa_rev_file",size=>"50",value=>""}," Concatenated/Reverse"),
			$cgi->img({height => "20", alt=>"OMSSA reverse files help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=OMSSA_decoy_files','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
                        $cgi->input({type=>"file", id=>'textarea_border', name=>"omssa_for_file",size=>"50"}," Forward"),
			$cgi->img({height => "20", alt=>"OMSSA forward files help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=OMSSA_forward_files','helpwindow','width=450,height=250');"}),
                        $cgi->br(),

			$cgi->br(),
                        ),
                        $cgi->div({id=>"innerboxheader"},
                        $cgi->p("X!Tandem File(s) Upload")),
                        $cgi->div({id=>"innerbox",align=>"left"},
                        $cgi->br(),
			$cgi->input({type=>"text",id=>'textarea_border', name=>"omssa_software_version",size=>"15",value=>"2008.12.1.1"}),
			"   X!Tandem software version",
			$cgi->img({height => "20", alt=>"X!Tandem version help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=Tandem_version','helpwindow','width=450,height=250');"}),
                       $cgi->br(),
                        $cgi->popup_menu({type=>"popup_menu",id=>'textarea_border', name=>"tandem_search_type",width=>"50",values=>[" ","concatenated forward/reverse","separate forward/reverse"]}),
 			"   Select search type",
			$cgi->img({height => "20", alt=>"X!Tandem search help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=Tandem_search','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
                        $cgi->input({type=>"file", id=>'textarea_border', name=>"tandem_rev_file",size=>"50",value=>""}," Concatenated/Reverse"),
			$cgi->img({height => "20", alt=>"X!Tandem reverse files help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=Tandem_decoy_files','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
                        $cgi->input({type=>"file", id=>'textarea_border', name=>"tandem_for_file",size=>"50"}," Forward"),
			$cgi->img({height => "20", alt=>"X!Tandem forward files help", src => "http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg", onclick =>"window.open('http://ispider.smith.man.ac.uk/cgi-bin/FDRAnalysis/help_page.pl?content=Tandem_forward_files','helpwindow','width=450,height=250');"}),
                        $cgi->br(),
			$cgi->br(),
                        $cgi->input({-type=>"hidden",-name=>"content",value=>"analyse"}),
                        ),
                        $cgi->div({id=>"innerboxheader"},
                        $cgi->p("")),
                        $cgi->br(),
 			$cgi->submit(-value=>"Perform Analysis", id=>'textarea_border'),
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
   $refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse\">";
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
   $refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse\">";
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
   $refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse\">";

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
   $refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse\">";
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
  $refresh_browser = "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"2\;URL=http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse\">";

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
    elsif($cgi->param('mascot_search_type') =~ m/separate\s+forward\/reverse/ && !$cgi->param("for_file"))
    {
    ErrorMsg("You requested the analysis of two separate files yet you have only uploaded one file.\n");
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
    elsif($cgi->param('omssa_search_type') =~ m/separate\s+forward\/reverse/ && !$cgi->param("for_file"))
    {
    ErrorMsg("You requested the analysis of two separate files yet you have only uploaded one omssa file.\n");
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
   print PARAM "oversize\t" . $cgi->param("oversize") . "\n";
   print PARAM "combine\t" . $cgi->param("combine") . "\n";
   close PARAM;
   }

   #add the max_expect
   open(PARAM,">>$params") or die "unable to open the file $params to write to\n";
   print PARAM "max_expect\t" . $cgi->param("max_expect") . "\n";
   print PARAM "fdr_value\t" . $cgi->param("fdr_value") . "\n";
   print PARAM "rev_tag\t" . $cgi->param("rev_tag") . "\n";
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
    elsif($cgi->param('tandem_search_type') =~ m/separate\s+forward\/reverse/ && !$cgi->param("for_file"))
    {
    ErrorMsg("You requested the analysis of two separate files yet you have only uploaded one tandem file.\n");
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
  if(!$cgi->param("rev_tag") && ($cgi->param('tandem_search_type') =~ m/concat/ || $cgi->param('mascot_search_type') =~ m/concat/ || $cgi->param('omssa_search_type') =~ m/concat/))
  {
  ErrorMsg("You must define the reverse tag used in a concatenated search otherwise I cannot determine which are decoy proteins!\n");
  exit(1);
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
			$cgi->p("This page will automatically refresh, if it doesn't please click ",$cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse"},"here")),
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
   #actually a blank image as the results are now in the blurb
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
   $blurb = "<BR><P>This plot is the distribution of scores for both the forward and reverse hits.  Two normal distributions should be evident, a distribution for reverse scores centered around a poor score and a separate distribution for the forward hits, centered around a better score.";
   $title = "Score distribution";
   }
   elsif($cgi->param("result_view") eq "zoomscore")
   {
    for(my $i=0 ; $i<scalar(@FDR_image_file) ; $i++)
    {
    $FDR_image_file[$i] =~ s/\.png/ZoomScore\.png/;
    $final_image_line .=  "<IMG SRC=\"$FDR_image_file[$i]\" title=\"Zoomed Score distribution plot image\">";
    }
   $blurb = "<BR><P>This plot shows the estimation of the distribution of correct and incorrect spectra at different scores.  The two distributions should be separate, with the incorrect distribution centered with the poorer scores and the correct at the more confident scores.";
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
   $nterminal = qq{<li><a href="http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse;result_view=nterdist">Nterminal Distribution</a></li>};
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
                                		 $cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse;result_view=summary"},"Summary")),
							$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse;result_view=combined"},"Peptide Overlap")),
$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse;result_view=rankplot"},"Rank plot")),
                               				 $cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse;result_view=deltamass"},"Delta Mass")),
                                			$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse;result_view=scoredist"},"Score Distribution")),
                                			$nterminal,
							$cgi->li($cgi->a({href=>"http://www.ispider.manchester.ac.uk/cgi-bin/FDRAnalysis/FDRanalysis_NEW.pl?content=analyse;result_view=zoomscore"},"Correct/Incorrect"))
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
				$cgi->p("This software enables the user to analyse the results of a Mascot, Omssa and /or X!Tandem database search.  It is specifically targetted towards searches including a decoy database, either searched separately or concatenated to the target database.")),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"analysis"},$cgi->p("Analysis Details"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p($cgi->b("N terminal analysis"),"if checked will include an N-terminal analysis, which will determine the fraction of identifications that are the N-terminal most peptide of the protein (based on the start position reported by each of the different search engines)"),
				$cgi->p($cgi->b("Decoy tag")," This is the 'tag' used to mark which proteins are from the decoy search.  This is required when a concatenated search was performed.  It can be included if a tag was used in the decoy database of a search done separately against target and decoy databases.  It is not required for a Mascot decoy search.  For example if the target database has accession >P12345 and the equivalent decoy protein is >REVERSE_P12345 then the tag is 'REVERSE_'.  <B>Please note that the tag must be at the start of the protein accession and have no spaces between it and the protein</B>."),
				$cgi->p($cgi->b("Combine Results"),"This employs the algorithm described in Jones et al. Proteomics (2009) 9: 1220-9 to combine the results of two or more searches from the different search engines.  This search can take a few extra minutes (depending on file sizes) and so please be patient."),
				$cgi->p($cgi->b("Maximum Expectation value"),"Related to one of the resulting plots, useful for investigating the results."),
 				$cgi->p($cgi->b("FDR rate"),"The false discovery rate (FDR) you would like the results for."),
  				),	
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"mascot_files"},$cgi->p("Mascot File(s) Upload"))),
				$cgi->div({id=>"innerbox",align=>"left"},
                            $cgi->p($cgi->b("Search type")," There are 3 types of decoy search a concatenated target/decoy whereby the target and decoy databases were joined and one search was performed.  A Mascot decoy search whereby the 'decoy' checkbox was selected in the Mascot search and finally two searches, one on the target database and a second on the decoy database. Mascot '.dat' files are expected here."),
				$cgi->p("The first file upload box is the main upload box and is for any file except the search run solely on the target database, which the second box will contain if applicable."),
				),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"omssa_files"},$cgi->p("OMSSA/X!Tandem File(s) Upload"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p($cgi->b("Search type")," There are 2 types of decoy search for both search engines, a concatenated target/decoy whereby the target and decoy databases were joined and one search was performed or where two searches were performed one on the target database and a second on the decoy database. OMSSA '.csv' result files and X!Tandem '.xml' result files are expected here."), 
				$cgi->p("The first file upload box is the main upload box and is for any file except the search run solely on the target database, which the second box will contain if applicable."),
				),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_summary"},$cgi->p("Results Pages - summary"))),
			        $cgi->div({id=>"innerbox",align=>"left"},
			        $cgi->p("A simple table displays the peptide frequencies identified from the search at the chosen FDR.  There are two different FDR calculation methods the first based on a method described by Elias & Gygi (2007) Nat. Methods 4: 207-214 and the second a method described by Kall et al. (2008) J. Proteome Res. 7:29-34.  If a combined search was selected these results will be displayed in the table on this page also.  This page also enables the user to download a text file containing a tab separated file of all the peptides and associated proteins identified by the different search engines at the chosen FDR. "),
				),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_rank"},$cgi->p("Results Pages - Rank Plot"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p(" If target and decoy peptides are equally likey to be selected as matches, lower ranked peptides (which are assumed to be incorrect) should be evenly distributed between the decoy and database portions.  The top ranking identifications on the otherhand are assumed to be correct and so should have a large bias towards target peptides.  The e-value cut-off for the top ranking peptides can be changed on the query for, the default is 0.05.  Rank plots are shown for the different search engines used in the query form.  Note that X!Tandem provides only the top ranked peptide identifications and OMSSA doesn't always report all ranks.  The description of the plot can be seen to the left of the plot."),
				),
				$cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_delta"},$cgi->p("Results Pages - Delta Mass Plot"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p("Delta Mass plots are shown for the different search engines used in the query form.  These plots represent the difference in mass between the calculated and experimental masses at different scores.  They are based on Elias & Gygi (2007) Nat. Methods 4: 207-214.  Note that X!Tandem provides only the top ranked peptide identifications and OMSSA doesn't always report all ranks.  The description of the plot can be seen to the left of the plot."),
                                ),
                                $cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_score"},$cgi->p("Results Pages - Score Distribution"))),
				$cgi->div({id=>"innerbox",align=>"left"},
			        $cgi->p("Score distribution plots are shown for the different search engines used in the query form and can be used to look at the distribution of decoy and target identifications at different scores.  Peptides with lower scores are less likely to be correct and so these plots should show an approximate equal districution of target/decoy peptides at low scores, but a clear difference in the two at higher scores. "),
	                        ),	
                                $cgi->div({id=>"innerboxheader"},
				$cgi->a({name=>"results_correct"},$cgi->p("Results Pages - Correct/Incorrect"))),
				$cgi->div({id=>"innerbox",align=>"left"},
				$cgi->p("Correct/Incorrect plots are shown for the different search engines used in the query form. They are based on those described in Elias & Gygi (2007) Nat. Methods 4: 207-214.  For separate searches estimated correct identifications were calculated by subtracting  the number of decoy peptides from target peptides at a given score threshold, and incorrect identifications were estimated as the minumum number of target or decoy peptides returned at a given score.  For concatenated search correct identifications were estimated by subtracting twice the number of decoy peptides at a given score threshold and incorrect estimated by doubling the number of decoy peptides at a given threshold."),
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
#print $cgi->div({id=>'header'},#$cgi->a({href=>"FDRanalysis_NEW.pl?content=home"}), 
#print 	$cgi->ul({id=>"navigation"},
#			$cgi->li($cgi->a({href=>"FDRanalysis_NEW.pl?content=home"},$cgi->img({width=>"60",height=>"50",alt=>"link to home",title=>"HOME",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/home.jpg"}))),
 #                    $cgi->li($cgi->a({href=>"FDRanalysis_NEW.pl?content=search"},$cgi->img({title=>"SEARCH", width=>"60",height=>"50",alt=>"link to search",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/search.jpg"}))),
#			$cgi->li($cgi->a({href=>"FDRanalysis_NEW.pl?content=analyse"},$cgi->img({title=>"ANALYSE",width=>"60",height=>"50",alt=>"link to FDR analysis",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/analysis.jpg"}))),
#			$cgi->li($cgi->a({href=>"FDRanalysis_NEW.pl?content=help"},$cgi->img({title=>"HELP", width=>"60",height=>"50",alt=>"link to help",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg"}))),
#			$cgi->li($cgi->a({href=>"mailto:david.wedge\@manchester.ac.uk"},$cgi->img({title=>"CONTACT",width=>"60",height=>"50",alt=>"link to contact",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/contact.jpg"})))
#		));

 print $cgi->div({id=>"navigation"},	$cgi->table({width=>"100%"},
		#$cgi=>tr(
 			$cgi->td($cgi->a({href=>"FDRanalysis_NEW.pl?content=home"},$cgi->img({width=>"60",height=>"50",alt=>"link to home",title=>"HOME",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/home.jpg"}))),
                    $cgi->td($cgi->a({href=>"FDRanalysis_NEW.pl?content=search"},$cgi->img({title=>"SEARCH", width=>"60",height=>"50",alt=>"link to search",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/search.jpg"}))),
			$cgi->td($cgi->a({href=>"FDRanalysis_NEW.pl?content=analyse"},$cgi->img({title=>"ANALYSE",width=>"60",height=>"50",alt=>"link to FDR analysis",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/analysis.jpg"}))),
			$cgi->td($cgi->a({href=>"FDRanalysis_NEW.pl?content=help"},$cgi->img({title=>"HELP", width=>"60",height=>"50",alt=>"link to help",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg"}))),
			$cgi->td($cgi->a({href=>"mailto:david.wedge\@manchester.ac.uk"},$cgi->img({title=>"CONTACT",width=>"60",height=>"50",alt=>"link to contact",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/contact.jpg"})))
		)
 			);

#print $cgi->div({id=>'navigation'},
#$cgi->table({align=>"center"},
#	$cgi->tr(
#			$cgi->a({href=>"FDRanalysis_NEW.pl?content=home"},$cgi->img({width=>"60",height=>"50",alt=>"link to home",title=>"HOME",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/home.jpg"})),
#	              $cgi->a({href=>"FDRanalysis_NEW.pl?content=search"},$cgi->img({title=>"SEARCH", width=>"60",height=>"50",alt=>"link to search",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/search.jpg"})),
#			$cgi->a({href=>"FDRanalysis_NEW.pl?content=analyse"},$cgi->img({title=>"ANALYSE",width=>"60",height=>"50",alt=>"link to FDR analysis",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/analysis.jpg"})),
#			$cgi->a({href=>"FDRanalysis_NEW.pl?content=help"},$cgi->img({title=>"HELP", width=>"60",height=>"50",alt=>"link to help",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/help.jpg"})),
#			$cgi->a({href=>"mailto:david.wedge\@manchester.ac.uk"},$cgi->img({title=>"CONTACT",width=>"60",height=>"50",alt=>"link to contact",src=>"http://ispider.smith.man.ac.uk/FDRAnalysis/contact.jpg"}))
#		)));
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
  -head=>$cgi->meta({-http_equiv=>'refresh',-content=>"3;FDRanalysis_NEW.pl?content=analyse"}),
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


#$cmd .= " &"; #this MUST be commented out if using QSUB!!!
#system($cmd);
#print "cmd is $cmd\n";
my $shell = "/tmp/RunFDRAnalysis_" . $session->id . ".sh";

open(SHELL,">$shell") or die "unable to create a shell file, $shell\n";
print SHELL "#!/bin/sh \n";
print SHELL $cmd;
print SHELL "\n\n";

my $sys_call = $QSUB." -q www -o /var/www/tmp -e /var/www/tmp " . $shell . " &";
print "syscall is $sys_call\n";
system($sys_call);


return 1;
}



