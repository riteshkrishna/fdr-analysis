#!/usr/bin/perl

#####################################################################################################
#        Copyright (C) 2009, David Wedge, University of Manchester                                  #
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

use CGI qw/:standard/;

my $cgi = new CGI;
my $subject = ($cgi->param("content"));

if($subject eq "N_terminal_analysis")
{
	print header,
	start_html({-title=>'Help: N-terminal analysis','-onload' => 'focus();'}),
	h2('Help: N-terminal analysis'),
	"If the N-terminal analysis box is checked, the analysis will determine the fraction of identifications that are the N-terminal most peptide of the protein (based on the start position reported by each of the different search engines).";
}
elsif($subject eq "decoy_tag")
{
	print header,
	start_html('Help: decoy tags'),
	h2('Help: decoy tags'),
	"This is the 'tag' used to mark which proteins are from the decoy search.  This is required when a concatenated search was performed.  It can be included if a tag was used in the decoy database of a search done separately against target and decoy databases.  It is not required for a Mascot decoy search.  For example if the target database has accession >P12345 and the equivalent decoy protein is >REVERSE_P12345 then the tag is 'REVERSE_'.  <B>Please note that the tag must be at the start of the protein accession and have no spaces between it and the protein</B>.";
}
elsif($subject eq "expectation_value")
{
	print header,
	start_html('Help: expectation value'),
	h2('Help: expectation value'),
	"The maximum expectation value is the cut-off for one of the resulting plots. This plot is useful for investigating the results.";
}
elsif($subject eq "FDR")
{
	print header,
	start_html('Help: FDR'),
	h2('Help: false discovery rate'),
	"The false discovery rate (FDR) you would like the results for.";
}
elsif($subject eq "combine_results")
{
	print header,
	start_html('Help: combining results'),
	h2('Help: combining results'),
	"If 'combine results' is true, the algorithm described in <i>Jones et al. Proteomics (2009) 9: 1220-9</i> is employed to combine the results of two or more searches from the different search engines.  This search can take a few extra minutes (depending on file sizes), so please be patient.";
}
elsif($subject eq "mascot_search")
{
	print header,
	start_html('Help: Mascot search types'),
	h2('Help: Mascot search types'),
	"There are 3 types of decoy search: <ul><li>A concatenated target/decoy whereby the target and decoy databases were joined and one search was performed.<li>A Mascot decoy search whereby the 'decoy' checkbox was selected in the Mascot search.<li>Two separate searches, one on the target database and a second on the decoy database.";
}
elsif($subject eq "mascot_decoy_files")
{
	print header,
	start_html('Help: Mascot decoy and concatenated files'),
	h2('Help: Mascot decoy and concatenated files'),
	"A Mascot '.dat' file is expected here. This upload box should be used to upload the results of searching aganist a decoy, concatenated or reverse database. In the last case, the results of a forward search should be uploaded in the box below.";
}
elsif($subject eq "mascot_forward_files")
{
	print header,
	start_html('Help: Mascot forward files'),
	h2('Help: Mascot forward files'),
	"A Mascot '.dat' file is expected here. This upload box should only be used to upload the results of a search run solely on the target database, without decoys. The results of a decoy (or reverse) search should be uploaded in the box above.";
}
elsif($subject eq "OMSSA_search")
{
	print header,
	start_html('Help: OMSSA search types'),
	h2('Help: OMSSA search types'),
	"There are 2 types of decoy search for OMSSA:<ul><li>A single search was performed on a concatenated target/decoy database containing both targets and decoys.<li>Two separate searches were performed, one on the target database and one on the decoy database.</ul>";
}
elsif($subject eq "OMSSA_version")
{
	print header,
	start_html('Help: OMSSA version'),
	h2('Help: OMSSA version'),
	"Enter the OMSSA version number used to perform the searches. Run 'omssacl -version' to check your software version.";
}
elsif($subject eq "OMSSA_decoy_files")
{
	print header,
	start_html('Help: OMSSA reverse and concatenated files'),
	h2('Help: OMSSA reverse and concatenated files'),
	"An OMSSA '.csv' results file is expected here. This upload box should be used to upload the results of searching aganist a concatenated or reverse database. In the second case, the results of a forward search should be uploaded in the box below.";
}
elsif($subject eq "OMSSA_forward_files")
{
	print header,
	start_html('Help: OMSSA forward files'),
	h2('Help: OMSSA forward files'),
	"An OMSSA '.csv' results file is expected here. This upload box should only be used to upload the results of a search run solely on the target database, without decoys. The results of a decoy (or reverse) search should be uploaded in the box above.";
}
elsif($subject eq "Tandem_search")
{
	print header,
	start_html('Help: X!Tandem search types'),
	h2('Help: X!Tandem search types'),
	"There are 2 types of decoy search for X!Tandem:<ul><li>A single search was performed on a concatenated target/decoy database containing both targets and decoys.<li>Two separate searches were performed, one on the target database and one on the decoy database.</ul>";


}
elsif($subject eq "Tandem_version")
{
	print header,
	start_html('Help: X!Tandem version'),
	h2('Help: X!Tandem version'),
	"Enter the X!Tandem version number used to perform the searches. Run 'tandem' without any command line parameters to check your software version.";
}
elsif($subject eq "Tandem_decoy_files")
{
	print header,
	start_html('Help: X!Tandem reverse and concatenated files'),
	h2('Help: X!Tandem reverse and concatenated files'),
	"An X!Tandem '.csv' results file is expected here. This upload box should be used to upload the results of searching aganist a concatenated or reverse database. In the second case, the results of a forward search should be uploaded in the box below.";
}
elsif($subject eq "Tandem_forward_files")
{
	print header,
	start_html('Help: X!Tandem forward files'),
	h2('Help: X!Tandem forward files'),
	"An X!Tandem '.xml' results file is expected here. This upload box should only be used to upload the results of a search run solely on the target database, without decoys. The results of a decoy (or reverse) search should be uploaded in the box above.";
}
elsif($subject eq "parent_mass_type")
{
	print header,
	start_html('Help: mass types'),
	h2('Help: mass types'),
	"There are 2 possible atomic mass types used during search. The monoisotopic mass is the mass of the most commonly occurring isotope, the average mass is averaged across all isotopes (weighted by distribution).";
}
elsif($subject eq "fragment_mass_type")
{
	print header,
	start_html('Help: mass types'),
	h2('Help: mass types'),
	"There are 2 possible atomic mass types used during search. The monoisotopic mass is the mass of the most commonly occurring isotope, the average mass is averaged across all isotopes (weighted by distribution).";
}
elsif($subject eq "select_enzyme")
{
	print header,
	start_html('Help: select enzyme'),
	h2('Help: select enzyme'),
	"Enter the enzyme used to digest proteins prior to mass spectrometry.";
}
elsif($subject eq "fixed_mods")
{
	print header,
	start_html('Help: fixed mods'),
	h2('Help: fixed modifications (mods)'),
	"Select the chemical modifications that were applied to every peptide during the search.";
}
elsif($subject eq "variable_mods")
{
	print header,
	start_html('Help: variable mods'),
	h2('Help: variable modifications (mods)'),
	"Select the variable chemical modifications that were applied during the search. Variable mods may not apply to every peptide found.";
}
else
{
	print header,
	start_html('Error'),
	h2('Error: unknown error topic'),
}
