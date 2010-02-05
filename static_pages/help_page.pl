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
	"If the N-terminal analysis box is checked, the analysis will include an N-terminal analysis, which will determine the number of identifications that are the N-terminal most peptide of the protein (based on the start position reported by each of the different search engines).";
}
elsif($subject eq "decoy_tag")
{
	print header,
	start_html('Help: decoy tags'),
	h2('Help: decoy tags'),
	"This is the 'tag' used to mark the proteins which are decoys. For example, if the target database has a protein with accession number >P12345 and the equivalent decoy protein is >REVERSE_P12345 then the tag is 'REVERSE_'. It is very important for the FDR calculations that this tag is entered correctly.";
}
elsif($subject eq "decoy_ratio")
{
	print header,
	start_html('Help: decoy size'),
	h2('Help: decoy size'),
	"This is the size of the decoy database, relative to the target database. If the 2 databases are of equal size this value should be 1. It is important for the FDR calculations that this value is entered is correctly.";
}

elsif($subject eq "expectation_value")
{
	print header,
	start_html('Help: expectation value'),
	h2('Help: expectation value'),
	"The maximum expectation value determines which identifications will be included in the rank plot and start position plot.";
}
elsif($subject eq "FDR")
{
	print header,
	start_html('Help: FDR'),
	h2('Help: false discovery rate'),
	"Results will be reported up to the false discovery rate (FDR) entered here.";
}
elsif($subject eq "combine_results")
{
	print header,
	start_html('Help: combining results'),
	h2('Help: combining results'),
	"If 'combine results' is true, the algorithm described in <i>Jones et al. Proteomics (2009) 9: 1220-9</i> is employed to combine the results of two or more searches from the different search engines.  This search can take a few extra minutes (depending on file sizes), so please be patient.";
}
elsif($subject eq "input_files")
{
	print header,
	start_html('Help: Input files'),
	h2('Help: Input files'),
	"Uploaded files should be the results of a concatenated target/decoy search, in which the target and decoy databases were joined and one search was performed.";
}
elsif($subject eq "file_types")
{
	print header,
	start_html('Help: Input file types'),
	h2('Help: Input file types'),
	"The accepted files types are:  Mascot .dat files; OMSSA .csv files; <span style='white-space:nowrap'>X!Tandem</span> .xml files; .mzident files from either Mascot, OMSSA or <span style='white-space:nowrap'>X!Tandem</span>.";
}
elsif($subject eq "Tandem_version")
{
	print header,
	start_html('Help: X!Tandem version'),
	h2('Help: X!Tandem version'),
	"Enter the X!Tandem version number used to perform the searches. Run 'tandem' without any command line parameters to check your software version.";
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
elsif($subject eq "tolerance")
{
	print header,
	start_html('Help: peptide mass tolerance'),
	h2('Help: peptide mass tolerance'),
	"A peptide mass identified by mass spectometry may not exactly match the theoretical mass. The mass tolerance sets the level of 'error' allowed between the theoretical and observed masses when making a positive identification. Mass tolerances may be set in Daltons or in parts per million, and separate tolerances are set for parent and fragment ions.";
}

elsif($subject eq "select_enzyme")
{
	print header,
	start_html('Help: select enzyme'),
	h2('Help: select enzyme'),
	"Enter the enzyme used to digest proteins prior to mass spectrometry.";
}
elsif($subject eq "missed_cleavages")
{
	print header,
	start_html('Help: missed cleavages'),
	h2('Help: missed cleavages'),
	"Enzymes cleave proteins at specific amino acid residues. However, the cleavage of proteins is often incomplete. Enter a maximum number of missed cleavages to be allowed in the resultant peptides.";
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
