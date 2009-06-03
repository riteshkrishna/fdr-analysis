package DecoyMascotParser;
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


####################################################################
#Conatins all the pasers for the different identification software!#
####################################################################


our @ISA = qw(Exporter);
our $VERSION = 1.0;;
our @EXPORT = qw(ParseDecoyMascot);

use strict;
use URI::Escape;

use lib qw(~/LIVERPOOL/bin/perl_modules/);




sub ParseDecoyMascot
{

#Written by Jennifer Siepen June 2006, adapted February 2008
#parse a dat file


my $mascot_file = shift;
my $mgf_file = shift;
my $db_path = shift;
my $output_type = shift;
my $cutoff = shift;
my $REVTAG = shift;
my @spectrum_number;


 if(!$mascot_file && !$mgf_file && !$db_path)
 {
 print "ERROR: No dat file and/or mgf file passed to the Mascot Parser\n";
 exit(1)
 }

 #if there is no output type use machine as default
 if(!$output_type)
 {
 $output_type = "M";
 }


#get the unimod modifications!
#JS These are now actually in the Mascot dat file - need to parse out!
my %unimod;



my $params = 0;
my $summary = 0;
my $peptides = 0;
my $proteins = 0;
my $query = 0;
my %search;
my %ions;
my @results;
my $enzyme = 0;
my $unimodxml = 0;
my $unimodtxt;

#Firstly get the hash of mgf titles - to correspond to the spectrum number
my %mgf_titles;

 if($mgf_file)
 {
 %mgf_titles  = GetMGFDetails($mgf_file);
 }
 else
 {
 #print "no MGF file supplied!\n";
 }

#open the mascot file
open (MASCOTFILE,"<$mascot_file") or die "unable to open the mascot file, $mascot_file\n";
 while(my $line = <MASCOTFILE>)
 {
  #starts at 'name="parameters"'
  if($line =~ m/name=\"parameters\"/)
  {
  $params = 1;
  }
  elsif($line =~ m/name=\"unimod\"/)
  {
  $params = 0;
  $unimodxml = 1;
  }
  elsif($line =~ m/name=\"enzyme\"/)
  {
  $unimodxml = 0;
   if($unimodtxt)
   {
   #%unimod = GetUniModInfo($unimodtxt);
   }
  $enzyme = 1;
  } 
  elsif($unimodxml == 1)
  {
  $unimodtxt .= $line;
  }
  #end the enzyme bit and start the summary
  elsif($line =~ m/name=\"decoy\_summary\"/)
  {
  $enzyme = 0;
  $summary = 1;
  }
  elsif($line =~ m/name=\"peptides\"/)
  {
  $summary = 0;
  }
  #end the summary bit and start the peptides
  elsif($line =~ m/name=\"decoy\_peptides\"/)
  {
  $summary = 0;
  $peptides = 1;
  }
  #end the peptides bit and start the proteins
  elsif($line =~ m/name=\"proteins\"/)
  {
  $peptides = 0;
  $proteins = 1;
  }
  #end the proteins and get the query stuff
  elsif($line =~ m/name=\"query\d+\"/)
  {
  $proteins = 0;
  #get the query number
  my @tmp = split/name\=/, $line;
  $query = $tmp[1];
  $query =~ s/\D+//g;
  }
  #end the query stuff
  elsif($line =~ m/name=\"index"/)
  {
  $query = 0;
  }
  #use the query mass & charge to map id's back to MGF file
  elsif($summary>0 && $line =~ m/qexp\d+\=/)
  {
  my @tmp = split/\=/,$line;
  #the query number is tmp[0]
  $tmp[0] =~ s/qexp//;

  #the mass & intensity are tmp1
  my @get_values = split/\,/,$tmp[1];
  #round the mass to 4dp to match the mgf file
  $get_values[0] = sprintf "%.4f", $get_values[0];
  #and the charge (remove the sign)
  $get_values[1] =~ s/\+//;
  $get_values[1] =~ s/\r\n//;
  $get_values[1] =~ s/\n//;
  #now save the full name
  my $match = $get_values[0] . "-" . $get_values[1];
   if($mgf_file)
   { 
   $spectrum_number[$tmp[0]] = $mgf_titles{$match};
   }
   else
   {
   $spectrum_number[$tmp[0]] = $tmp[0];
   }
  }	  
  #if it is the results summary - get the precursor information
  elsif($summary>0 && $line =~ m/qmatch\d+\=/)
  {
   if($line =~ m/^qmatch/)
   {
   my @tmp = split/\=/, $line;
   #what is the query number?
   $tmp[0] =~ s/\D+//;
   #what is the intensity
   $tmp[1] =~ s/\r\n//;
   $ions{'match'}[$tmp[0]] = $tmp[1];
   #and store the qmatch for this query to calculate the mascot Expect later
   $tmp[1] =~ s/\n//;
   $results[$spectrum_number[$tmp[0]]][1]{'qmatch'} = $tmp[1];
   }
  }
  #if it is the peptide results
  elsif($peptides == 1 && $line =~ m /^q/)
  {
   #I am only interested in all ranking peptides
   if($line =~ m/q\d+\_p\d+\=/)
   {
   my @tmp = split/\=/, $line;
    #as long as it wasn't a '-1'
    if($tmp[1] !~ m/^\-1/)
    {
    #get the q and p numbers
    my @num = split/\_/,$tmp[0];
    $num[0] =~ s/q//;
    $num[1] =~ s/p//;
    $tmp[1] =~ s/\n//;
      
    #split the string
    my @tmp2 = split/\,/,$tmp[1];

    #and now get the start residue
    my @split1 = split/\;/,$tmp2[10];
     #really I should be doing this for all the proteins - use the first unless another is REV
     my $protein;
     my $start = 0;
     for(my $p=1 ; $p<scalar(@split1) ; $p++)
     {
     my @split2 = split/\:/,$split1[$p];
     #and the protein name (is it a reverse?)
     $split2[0] =~ s/\"//g;
      #if there is no protein already
      if(!$protein)
      {
      $protein = $split2[0];
      $start = $split2[2];
      }
      if($protein && $split2[0]=~ m/^$REVTAG/)
      {
      $protein = $split2[0];
      $start = $split2[2];
      }   
     }
    $results[$spectrum_number[$num[0]]][$num[1]]{'sequence'} = $tmp2[4];
    $results[$spectrum_number[$num[0]]][$num[1]]{'ionscore'} = $tmp2[7];
    $results[$spectrum_number[$num[0]]][$num[1]]{'start'} = $start;
    $results[$spectrum_number[$num[0]]][$num[1]]{'protein'} = $protein;
    $results[$spectrum_number[$num[0]]][$num[1]]{'expect'} = GetExpect($tmp2[7],$results[$spectrum_number[$num[0]]][1]{'qmatch'});
    $results[$spectrum_number[$num[0]]][$num[1]]{'delta'} = $tmp2[2];
    }
}
    #need to put blanks if nothing there
    #else
    #{
    #my @tmp = split/\=/, $line;
    ##get the q and p numbers
    #my @num = split/\_/,$tmp[0];
   # #$num[0] =~ s/q//;
   # $num[1] =~ s/p//;
   #
   # $results[$spectrum_number[$num[0]]][$num[1]]{'sequence'} = "NULL";

    #$results[$spectrum_number[$num[0]]][$num[1]]{'ionscore'} = "-1";
    #$results[$spectrum_number[$num[0]]][$num[1]]{'expect'} = "100";
    #}
   
  }
  elsif($params == 1 && $line =~ m/\=/)
  {
  my @tmp = split/\=/, $line;
  $tmp[1] =~ s/\n//;
  $search{$tmp[0]} = $tmp[1];
  }

 }
close MASCOTFILE; 


#get the nter peps
#print "getting the N-terminal peptides\n";
#my %nter = GetNterPeps($db_path);
my %nter;


#what do we have in results?
 for(my $r=1 ; $r<scalar(@results) ; $r++)
 {

  if($results[$r][1]{'sequence'})
  {
   #for all the ranks
   for(my $rank=1 ; $rank<(scalar(@{$results[$r]})) ; $rank++)
   {
   my $pos = 0;
 
 
    if($results[$r][$rank]{'start'}<3)
    {
    $pos = 1;
    }
     $results[$r][$rank]{'nterminal'} = $pos; 
   }
  }
 }

 #if they want machine output
 if($output_type eq "M" || $output_type eq "m")
 {
 return(\@results);
 }

 #otherwise human readable
 else
 { 
 my %unique_peptides; 
 my $nter = 0;
 my $non_nter = 0; 
 my %unique_rev;
 my $rev = 0;
  #if there was not cutoff score set, set it to 0
  if(!$cutoff)
  {
  $cutoff = 0;
  }


  for(my $r=1 ; $r<scalar(@results) ; $r++)
  {
   if($results[$r][1]{'sequence'})
   {
    for(my $rank=1 ; $rank<2 ; $rank++)
    {
     #does it meet the search criteria?
     if($results[$r][$rank]{'expect'} < $cutoff)
     {
      if($results[$r][$rank]{'protein'} =~ m/^$REVTAG/ && !$unique_rev{$results[$r][$rank]{'sequence'}})
      {
      $rev++;
      $unique_rev{$results[$r][$rank]{'sequence'}} = 1;
      }
      elsif(!$unique_peptides{$results[$r][$rank]{'sequence'}} && $results[$r][$rank]{'nterminal'} == 1 && $results[$r][$rank]{'protein'} !~ m/^$REVTAG/)
      {
      $nter++;
      $unique_peptides{$results[$r][$rank]{'sequence'}} = 1;
      }
      elsif(!$unique_peptides{$results[$r][$rank]{'sequence'}} && $results[$r][$rank]{'nterminal'} == 0 && $results[$r][$rank]{'protein'} !~ m/^REV\_/)
      {
      $non_nter++;
      $unique_peptides{$results[$r][$rank]{'sequence'}} = 1; 
      }
     }
    }
   }
  }
 return(\@results);
 }

return 1;

}

return 1;


sub GetNterPeps
{

my $db_path = shift;
my %nterpeps;
my $seq = 0;

#open the sequence database
open(DB,"<$db_path") or die "unable to open the database file $db_path\n";
 while(my $line = <DB>)
 {
  if($line =~ m/^\>/)
  {
   $seq =~ s/\s+//g;
   $seq =~ s/\n//g;
   $seq =~ s/RP/\*/g;
   $seq =~ s/KP/\#/g;
   
   #split at all K
   my @ksplit = split/K/,$seq; 

   my $count = 0;

   #for the first two peps split at R
   for(my $i=0 ; $i<2 ; $i++)
   {
   my @rsplit = split/R/,$ksplit[$i];
    
    for(my $j=0 ; $j<scalar(@rsplit) ; $j++)
    {
     if($count<2)
     {
     $rsplit[$j] =~ s/\*/RP/g;
     $rsplit[$j] =~ s/\#/KP/g;
     $nterpeps{$rsplit[$j]} = 1;
     $count++; 

      if($j == 0)
      {
      $rsplit[$j] =~ s/^M//;
      $nterpeps{$rsplit[$j]} = 1;
      } 
     }
    }
   }
   if($count<2)
   {
$ksplit[1] =~ s/\*/RP/g;
$ksplit[1] =~ s/\#/KP/g;
   $nterpeps{$ksplit[1]} = 1;
   } 
  $seq = "";
  }
  else
  {
  $seq .= $line;
  }
 }   
return %nterpeps;
}


sub GetMGFDetails ()
{
my $mgf_file = shift;
my %titles;
my $count=1;

#get all the title
my $name="";
my $mass;
my $charge;
open(MGF,"<$mgf_file") or die "unable to open the MGF file in Mascot parser, $mgf_file\n";
 while(my $line = <MGF>)
 {
  if($line =~ m/^PEPMASS/)
  {
  my @tmp = split/\=/,$line;
  $tmp[1] =~ s/\n//;
  #round to 4dp to ease matching with mascot!
  my @tmp2 = split/\s+/,$tmp[1];
  $tmp[1] = $tmp2[0];
  $tmp[1] = sprintf "%.4f", $tmp[1];
  $mass = $tmp[1];
  }
  elsif($line =~ m/^CHARGE/)
  {
  my @tmp = split/\=/,$line;




  $tmp[1] =~ s/\r\n//;
  #remove the sign
  $tmp[1] =~ s/\+//;
  $charge = $tmp[1];
  }
  elsif($line =~ m/^END\s+IONS/)
  { 
  my $match = $mass . "-" . $charge;
  $titles{$match} = $count;

  $mass = "";
  $charge = "";
  $count++;
  }
 }
close MGF;
return %titles;
}

#sub GetUniModInfo
#{
#my $xml = shift;
#my %unimod;
#
##split the xml at ever new line
#my @lines = split/\n/,$xml;
#
#
# #for all the lines
# for(my $l=0 ; $l<scalar(@lines) ; $l++)
# { 
#  #only if it is a modification entry
#  if($lines[$l] =~ m/^\s+\<modifications_row/)
#  {
#  $lines[$l] =~ s/^\s+//;
#  my $record = 0;
#  my $avge = 0;
#  my $mono = 0;
# 
#  #get the info
#  my @tmp = split/[\=\"]/,$lines[$l];
 #  #for each of these get the name values
 #  for(my $i=1 ; $i<scalar(@tmp) ; $i++)
 #  {
 #  $tmp[$i] =~ s/^\s+//;
 #  #get the record id
 #   if($tmp[$i] eq "record_id")
 #   {
 #   $record = $tmp[$i+2];
 #   }
 #   elsif($tmp[$i] eq "full_name")
 #   {
 #   #lowercase
 #   $tmp[$i+2] =~ tr/[A-Z]/[a-z]/;
 #   $unimod{$tmp[$i+2]}[0] = $record;
 #   $unimod{$tmp[$i+2]}[1] = $avge;
 #   $unimod{$tmp[$i+2]}[2] = $mono;
 #   }
 ##   elsif($tmp[$i] eq "code_name")
 #   {
 #   #lowercase
 #   $tmp[$i+2] =~ tr/[A-Z]/[a-z]/;
 #   $unimod{$tmp[$i+2]}[0] = $record;
  #  $unimod{$tmp[$i+2]}[1] = $avge;
  #  $unimod{$tmp[$i+2]}[2] = $mono;
  #  }
  #  elsif($tmp[$i] eq "avge_mass")
  #  {
  #  $avge = $tmp[$i+2];
  #  }
  #  elsif($tmp[$i] eq "mono_mass")
  #  {
  #  $mono = $tmp[$i+2];
  #  }
  # }
 # }
# }
#
#return %unimod;


#}


sub GetExpect()
{

my $score = shift;
my $qmatch = shift;


my $entries = 0;
my $evalue= 10;

 if($qmatch && $qmatch > 0)
 {
 $evalue = 0.05 * $qmatch * (10**(-$score/10));
 }
return $evalue;
}



