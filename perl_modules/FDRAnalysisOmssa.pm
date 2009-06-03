package FDRAnalysisOmssa;
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
our @EXPORT = qw(RunOmssa ParseOmssa);


use strict;
use IPC::Open2;
use strict;
use DBI;
use LWP;
use HTML::LinkExtor;
use HTML::Form;
use URI;
use Getopt::Long;
use Data::Dumper;
use MIME::Base64;
use Math::Complex;
use POSIX qw(log10);
use Parser;



my $cgi;

#a globl string to hold the mod info
my %mod_params;
my %params;
my %search;
my @mod_order;


sub RunOmssa
{

my $commandline = shift;
my @arguments = split/\s+/,$commandline;


my $result_file = $arguments[0];

 for(my $i=1 ; $i<scalar(@arguments) ; $i++)
 {
 my @tmp = split/\:\:\:/,$arguments[$i];
 #Here get rid of wether it a guess or not
 $tmp[1] =~ s/\:\d+$//;
 $tmp[1] =~ s/\:\d+\,$//;
 $params{$tmp[0]} = $tmp[1];
 }

#create the input file to  Omssa 
 
#run tandem
RunOmssaMSMS($result_file);

return $result_file;
}

return 1;

sub ParseOmssa
{

my $OMSSA = shift;
my $database_searched = shift;
my $file;
my $html_stuff;
my $hits = 0;

my $found = 0;
my %protein_details;
my @results;
my %peptide;
my $pep_count=1;
my $db_source;
my $mod_name = $OMSSA . ".params";

#the results for each pep id
my $protein_hit; 
my %pep_result;
my %mgf_titles;
my $title = 0;
my @spectrum_number;
my $query = 0;
my %protein_group;

my $pep_count = 1;
my $query = 0;

my $NIST = 0;

open(RES,"<$OMSSA") or return "ERROR: There is a problem opening the OMSSA results, $OMSSA\n";
 while(my $line = <RES>)
 {
  #ignore the first line
  if($line =~m /^Spectrum\s+number/)
  {
   if($line =~ m/NIST\s+score/)
   {
   $NIST = 1;
   }
  next;
  }
  else
  {

  my @tmp = split/\,/,$line;
  #is it the same query?
   if($tmp[0] == $query)
   {
   $pep_count++;
   } 
   else
   {
   $query = $tmp[0];
   $pep_count = 1;
   }
  my $seq = $tmp[2];
  #omssa makes the residues with mods lower case - change them to upper
  $seq = uc($seq);
  my $mh = $tmp[4];
  my $size  = scalar(@tmp);
  #BUG FIX for large omssa description lines
  $tmp[12] = $tmp[$size-3];

  my $delta = $tmp[12] - $tmp[4];
  my $eval = $tmp[3];
  #if there is more than one mod in the variable sequence omssa puts commas therefore get the p val as the last element of the array
  my $last = scalar(@tmp) - 1;
   #BUG FIX for omssa results that have a NIST result
   $last = $last - $NIST;
  my $pval = $tmp[$last];
  $pval =~ s/\n//;

#   my @split_prot = split/\|/,$tmp[9];
#  my $prot = $split_prot[0];
  my $prot = $tmp[6];


#BUG FIX - Needs looking at
   #if the protein accession in 6 is all numbers use the start of the description line??
   if($prot !~ m/[ABCDEFGHIJKLMNOPQRSTUVQXYZ]/)
   {
   my @split_prot = split/\s+/,$tmp[9];
   $prot = $split_prot[0];
   }
  my $start = $tmp[7];
  my $end = $tmp[8];
  my $title = $tmp[1];
  my $mods = $tmp[10];

  #beware of commas in multiple mods!
  my $t = scalar(@tmp);
  my $theo_mass = $tmp[$t-2];
  my $charge = $tmp[$t-3];
  my $mod;
   #and the mods
   for(my $i=10 ; $i<($t-3) ; $i++)
   {
    if($tmp[$i])
    {
    $mod .= "," . $tmp[$i];
    }
   }
  $mod =~ s/^\,//;
  #if it is zero or
  $protein_details{$prot}[0] = "nan";

  #deal with the protein_details 
  #must account for protein accessions that contains a comma!
  my $quote_counter = 0;
  my $moveon = 0;
  my $seen = 0;
   for(my $element = 9 ; $element<scalar(@tmp) ; $element++)
   {
   $quote_counter = 0;
    while($tmp[$element] =~ m/\"/g)
    {
    $quote_counter++;
    }
    if($quote_counter == 0 && $seen == 0)
    { 
    $moveon = $element+1;
    last;
    }
    elsif($quote_counter == 0 && $seen == 1)
    {
    $protein_details{$prot}[1] . $tmp[$element];
    }
    elsif($quote_counter == 1 && $seen == 0)
    {
    $protein_details{$prot}[1] . $tmp[$element];
    $seen = 1;
    }
    elsif($quote_counter == 1 && $seen == 1)
    {
    $protein_details{$prot}[1] . $tmp[$element];
    $moveon = $element+1;
    last;
    }
    elsif($quote_counter == 2)
    {
    $protein_details{$prot}[1] . $tmp[$element];
    $moveon = $element+1;
    last;
    }
   }
  #and now group the results in proteins
  push(@{$protein_group{$prot}},$query . "," . $pep_count);

  $results[$query][$pep_count]{'sequence'} = $seq;
  $results[$query][$pep_count]{'expect'} = $eval;
  $results[$query][$pep_count]{'pvalue'} = $pval;
  my $ionscore = 0;
#   if($pval>0)
#   { 
#   $ionscore = log($pval);
#   }
#  $ionscore = $ionscore/log(10);
#  $ionscore = $ionscore*-10;
$ionscore = $eval;

$results[$query][$pep_count]{'ionscore'} = $ionscore;
  $results[$query][$pep_count]{'protein'} = $prot;
  $results[$query][$pep_count]{'start'} = $start;
  $results[$query][$pep_count]{'end'} = $end;
  $results[$query][$pep_count]{'delta'} = $delta;
  $results[$query][$pep_count]{'title'} = $title;
  $results[$query][$pep_count]{'charge'} = $charge;
  $results[$query][$pep_count]{'prec_masstocharge'} = $mh;
  $results[$query][$pep_count]{'theo_mass'} = $theo_mass;
  $results[$query][$pep_count]{'mods'} = $mod;
  #and the spectrum number
  $spectrum_number[$query] = $mgf_titles{$tmp[1]};
  }
 }
 for(my $m=0 ; $m<scalar(@mod_order) ; $m++)
 {
 my $t = $m+1;
 my $d = "delta" . $t;
 $search{$d} = "?,$mod_order[$m]";

 } 

close FILE;
my %nter;

if($database_searched)
{
%nter = GetNterPeps($database_searched);
}
 #what do we have in results?
 for(my $r=0 ; $r<scalar(@results) ; $r++)
 {
 my $pos = 0;

  if($nter{$results[$r][1]{'sequence'}})
  {
  $pos = 1;
  }
  #for all the keys
  foreach my $attr (keys %{$results[$r][1]})
  {

  $results[$r][1]{'nterminal'} = $pos;
  }
 }
return(\@results);

}


sub GetModString()
{

my $sequence = shift;
my $mod = shift;
my $length = length($sequence);

 my $string; 

 for(my $a=0 ; $a<($length+2) ; $a++)
 {
 $string .= "0";
 }

 if($mod)
 {
 #split the mod to get the position
 my @tmp = split/\:/,$mod;

 #get rid of the residue at the end
 $tmp[0] =~ s/\s+\w$//;

 my $position = $tmp[1];
 my $number;
 my $found = 0;
  #do we have this mod in the order thing?
  for(my $m=0 ; $m<scalar(@mod_order) ; $m++)
  {
   if($mod_order[$m] eq $tmp[0])
   {
   $number = $m+1;
   $found = 1;
   }
  } 
  if($found == 0)
  {
  push(@mod_order,$tmp[0]);
  $number = scalar(@mod_order);
  }
 
  #if it is the N or C term
  if($position =~ m/N/i)
  {
  $position = 0; 
  }
  elsif($position =~ m/C/i)
  {
  $position = $length+1;
  }
  #now redo the string
  $string = "";
  for(my $a=0 ; $a<($length+2) ; $a++)
  {
   if($a == $position)
   {
   $string .= $number;
   }
   else
   {
   $string .= "0";
   }
  }
 }
return $string;
}


sub GetMGFTitles ()
{
my $mgf_file = shift;
my %titles;
my $count=1;

#get all the title
my $name="";
my $mass;
my $charge;
open(MGF,"<$mgf_file") or print "**WARNING: unable to open the MGF file in Omssa, $mgf_file\n";
 while(my $line = <MGF>)
 {
  if($line =~ m/^TITLE/)
  {
  my @tmp = split/\=/,$line;
  $tmp[1] =~ s/\n//;
  $titles{$tmp[1]} = $count;
  }
  elsif($line =~ m/^END\s+IONS/)
  {
  $count++;
  }
 }
close MGF;
return %titles;
}



sub FillSearchHash()
{
my $command = shift;
my %search;

my @arguments = split/\s+/,$command;

my $result_file = $arguments[0];

 for(my $i=1 ; $i<scalar(@arguments) ; $i++)
 {
 my @tmp = split/\:\:\:/,$arguments[$i];
 #Here get rid of wether it a guess or not
 $tmp[1] =~ s/\:\d+$//;
 $tmp[1] =~ s/\:\d+\,$//;
 $search{$tmp[0]} = $tmp[1];
 }



return %search;
}


sub RunOmssaMSMS()
{

my $result_file = shift;

my $sys_call;

my $in_file = $params{'FILE'};

my $database = GetDbPath($params{'DB'});

#now I must copy the db to a local dir and format db it!
my $call = "/fs/ls1/home/mjfssjs2/LIVERPOOL/bin/perl_scripts/GetDbInfo.pl " . $params{'DB'} . " omssa";
system($call);

my @tmp = split/\//,$database;
my $size = scalar(@tmp);
$database = "/fs/ls1/home/mjfssjs2/LIVERPOOL/databases/" . $tmp[$size-1];


#deal with enzyme, if typsin then leave as default
my $enzyme = GetEnzymeRule($params{'CLE'});

my $tol;
my $itol;
 #and the error tolerance - needs to be in Daltons
 if($params{'ITOLU'} eq "Da")
 {
 $itol = $params{'ITOL'};
 }
 else
 {
 print "ERROR: OMSSA Does not support anything but Daltons! setting to default of 0.8 Da\n";
 $itol = "0.8";
 $params{'ITOL'} = "0.8";
 $params{'ITOLU'} = "Da";
 }
 if($params{'TOLU'} eq "Da")
 {
 $tol = $params{'TOL'};
 }
 else
 {
 print "ERROR: OMSSA Does not support anything but Daltons, setting to default!\n";
 $tol = "2.0";
 $params{'TOL'} = "2.0";
 $params{'TOLU'} = "Da"; 
 }

 #and the modifications
 my @varmods = split/\,/,$params{"IT_MODS"};
 my @fixmods = split/\,/,$params{"MODS"};
 my $vmods = GetMods(@varmods);
 my $fmods = GetMods(@fixmods);


  if($vmods)
  {
  $vmods = " -mv " . $vmods;
  }
  if($fmods)
  {
  $fmods = " -mf " . $fmods;
  }

#run the omssa program
my $system = "/usr/local/software/omssa/omssacl -fm " . $in_file . " -d " . $database . " -to " . $itol . " -te " . $tol . $fmods . " " . $vmods . " -oc " . $result_file;
print "system is $system\n";

system($system);

}


sub GetDBVersion()
{
my $db_name = shift;
my $version = "unknown" ;
my $found = 0;

 #open the tax file
 open(TAX,"</tmp/taxonomy.xml") or print "****WARNING: Unable to open the taxonomy file fo th db version - /tmp/taxonomy.xml\n";
  while(my $line = <TAX>)
  {
   if($line =~ m/label=\"$db_name\"/)
   {
   $found = 1;
   }
   elsif($line =~ m/\<taxon\s+label/)
   {
   $found = 0;
   }
   elsif($found == 1 && $line =~ /\<file\s+format/)
   {
   my @tmp = split/\s+/,$line;

    for(my $t=0 ; $t<scalar(@tmp) ; $t++)
    {
     if($tmp[$t] =~ m/^URL\=/)
     {
     $tmp[$t] =~ s/\"//g;
     my @tmp2 = split/\=/,$tmp[$t];
     $version = $tmp2[1];
     }
    }
   }
  }
 close TAX;
return $version;

}


sub GetMods()
{
my @mods = @_;
my $string="";
#get the mods from the mod_file
my %mascot_mods = GetModsFromModFile();

#the first mod is the mass type - i.e. monoisotopic or average
 for(my $m=0 ; $m<scalar(@mods) ; $m++)
 {

$mods[$m] =~ s/\*/\ /g;
$mods[$m] =~ s/\n//;

 #store the name of the mod
 #make sure there is no residue info at the end
  #ignore terminus for now
  if($mods[$m] !~ m/term/)
  {
   if($mascot_mods{$mods[$m]})
   {
   $string .=  $mascot_mods{$mods[$m]} . ",";
   }
  }
 } 

 if($string)
 {
 $string =~ s/\,$//;
 }

$string =~ s/\n//;
return $string;

}

sub GetEnzymeRule()
{

my $enzyme = shift;

 if($enzyme =~/^Trypsin$/i)
 {
 return;
 }
 elsif($enzyme eq "Arg-C")
 {
 return "-v 1";
 }
 elsif($enzyme eq "Asp-N")
 {
 return "-v 12";
 }
 elsif($enzyme eq "Chymotrypsin")
 {
 return "-v 3";
 }
 elsif($enzyme eq "CNBr")
 {
 return "-v 2";
 }
 elsif($enzyme eq "CNBr+Trypsin")
 {
 return "-v 8";
 }
 elsif($enzyme eq "Formic_acid")
 {
 return "-v 4";
 }
 elsif($enzyme eq "Lys-C")
 {
 return "-v 5";
 }
 elsif($enzyme eq "Lys-Ci/P")
 {
 return "-v 6";
 }
 elsif($enzyme eq "PepsinA")
 {
 return "-v 7";
 }
 elsif($enzyme eq "TrypChymo")
 {
 return "-v 9";
 }
 elsif($enzyme eq "TrypChymo")
 {
 return "-v 9";
 }
 elsif($enzyme eq "Trypsin/P")
 {
 return "-v 10";
 }
 elsif($enzyme eq "semiTrypsin")
 {
 return "-v 16";
 }
 elsif($enzyme eq "none")
 {
 return "-v 17";
 }
 else
 {
 return "";
 }


}

sub GetModsFromModFile()
{

my %mods;

my $title;

open(MOD,"</var/www/localhost/cgi-bin/pipeline/REQD_FILES/mod_file_mascot_vs_omssa");

 while(my $line = <MOD>)
 {
 
  if($line =~ m/^Title/)
  {
  my @tmp = split/\:/,$line;
  $title = $tmp[1];
  $title =~ s/\s+$//;
  }
  elsif($line =~ m/^OM/)
  {
  my @tmp = split/\:/,$line;
   if($title)
   {
   $mods{$title} = $tmp[1];
   }
  }
  elsif($line =~ m/^\*/)
  {
  $title = "";
  }
 }

close MOD;


return %mods;
}

sub UpdateTaxonomy()
{

my $file = shift;

my $db_counter = 0;
my $start = 0;
my @split_line;
open(TAX,">$file") or print "unable to open the taxonomy file to write to ($file)\n";;
print TAX "<?xml version=\"1.0\"?>\n<bioml label=\"x! taxon-to-file matching list\">\n";

open(DBFILE,"</tmp/mascot.dat") or die "unable to open the mascot.dat file\n";
 while(my $line = <DBFILE>)
 {
  #is it the start of a new entry?
  if($line =~ m/^Databases/)
  {
  $start = 1;
  }
  elsif($line =~ m/^end/)
  {
  $start = 0;
  }
  elsif($start == 1)
  {
  $line =~ s/^\#//;
  $line =~ s/^\s+//;
  my @split_line = split/\s+/,$line;
   #does this location already exist?
   if($split_line[2] =~ m/AA/)
   {
   $split_line[2] = "aminoacid";
   }
   else
   {
   $split_line[2] = "nucleic";
   }
  #get rid of trailing white space
  $split_line[1] =~ s/\s+//;
   
  #need to get full location, number of sequences and the system date
  my $filename;

  #first full filename
  #does file name contain a *?
   if($split_line[1] =~ m/\*/)
   {
   #if we are on ispider change the db path
   $split_line[1] =~ s/\/usr\/share\//\/db\//;
   $filename = GetFileName($split_line[1]);
   }
   else
   {
   $filename = $split_line[1];
   $filename =~ s/\/usr\/share\//\/db\//;
   }
   #check it isn't the line with my descrition in it!
   if($filename !~ m/by/)
   {

    if($split_line[2] eq "aminoacid")
    {
    print TAX "<taxon label=\"" . $split_line[0] . "\">\n";
    print TAX "\t<file format=\"peptide\" URL=\"" . $filename . "\" />\n";
    print TAX "</taxon>\n";
    }
   }
  }
 } 
close DBFILE;


print TAX "</bioml>\n";

close TAX;

}


sub GetFileName()
{

my $oldfile = shift;
my $newfile = $oldfile;

#get the directory name
my @tmp = split/\//, $oldfile;
my $dir = "";

 for(my $x=0 ; $x<scalar(@tmp)-1 ; $x++)
 {
 $dir = $dir . $tmp[$x] . "/";
 }
my $length = scalar(@tmp);
#open the directory and get the filename

opendir(DIR, "$dir") or return $newfile;
my @files = readdir DIR;
closedir DIR;

#find the filename I am interested in
my @name = split/\*/,$tmp[$length-1];
 #for all the files find 'the one'
 for(my $x=0 ; $x<scalar(@files) ; $x++)
 {
  if( $files[$x] =~ m/^$name[0]/ && $files[$x] =~ m/$name[1]/)
  {
  $newfile = $dir . $files[$x];
  }
 }

return $newfile;



}

sub GetDbPath
{

my $dbname = shift;
my $path;
#open the mascot.dat file
open(FILE,"</fs/ls1/home/mjfssjs2/LIVERPOOL/tmp_files/mascot.dat") or die "unable to open the file, /fs/ls1/home/mjfssjs2/LIVERPOOL/tmp_files/mascot.dat\n";

 while(my $line = <FILE>)
 {
  if($line =~ m/^$dbname/)
  {
  my @split = split/\s+/,$line;
  $path = $split[1];
  $path =~ s/\/usr\/local\/mascot\//\/fs\/msct\//;

  #now get the real name of the file
  my $call = "ls $path |";
  open(PIPE,$call);
   while(my $l = <PIPE>)
   {
   $l =~ s/\n//;
    if($l)
    {
    $path = $l;
    }
   }
  close PIPE;
  close FILE;
  }

 }
close FILE;

return $path;

}


sub GetNterPeps
{

my $db_path = shift;
my %nterpeps;
my $seq = 0;

#open the sequence database
open(DB,"<$db_path") or print "unable to open the database file $db_path\n";
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

