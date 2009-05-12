package FDRAnalysisXTandem;
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
our @EXPORT = qw(RunTandem ParseTandem);


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
use lib qw(~/LIVERPOOL/bin/perl_modules/);
use Parser;


my %aa;

$aa{'A'} = "71.03711";
$aa{'C'} = "103.00919";
$aa{'D'} = "115.02694";
$aa{'E'} = "129.04259";
$aa{'F'} = "147.06841";
$aa{'G'} = "57.02146";
$aa{'H'} = "137.05891";
$aa{'I'} = "113.08406";
$aa{'K'} = "128.09496";
$aa{'J'} = "113.08406";
$aa{'L'} = "113.08406";
$aa{'M'} = "131.04049";
$aa{'N'} = "114.04293";
$aa{'P'} = "97.05276";
$aa{'Q'} = "128.05858";
$aa{'R'} = "156.10111";
$aa{'S'} = "87.03203";
$aa{'T'} = "131.04768";
$aa{'V'} = "99.06841";
$aa{'W'} = "186.07931";
$aa{'Y'} = "163.06333";
$aa{'nterm'} = "0.00";
$aa{'cterm'} = "0.00";




#a globl string to hold the mod info
my %mod_params;
my %params;
my %search;

sub RunTandem
{

my $commandline = shift;

my @arguments = split/\s+/,$commandline;
my $taxonomy_file;

my $result_file = $arguments[0];
#default score is the X!tandem score
my $scoretype = "tandem";

 for(my $i=1 ; $i<scalar(@arguments) ; $i++)
 {
  
 my @tmp = split/\:\:\:/,$arguments[$i];
 #Here get rid of wether it a guess or not
 $tmp[1] =~ s/\:\d+$//;
 $tmp[1] =~ s/\:\d+\,$//;

  #if there is a taxonomy file
  if($tmp[0] eq "TAXONOMYFILE")
  {
  $taxonomy_file = $tmp[1];
  }
  elsif($tmp[0] eq "SCORE")
  {
  $scoretype = $tmp[1];
  }
  else
  { 
  $params{$tmp[0]} = $tmp[1];
  }

 }

 


#create the input file to x!tandem 
print "in the XTandem function and the result file is $result_file\n"; 
#run tandem
RunTandemMSMS($result_file,$taxonomy_file,$scoretype);

return $result_file;
}

return 1;

sub ParseTandem
{
my $TANDEM = shift;
my $input = shift;
my $taxonomy = shift;
my $pride_command = shift;

my $file;
my $html_stuff;
my $hits = 0;

my $found = 0;
my %protein_details;
my @results;
my %peptide;
my $pep_count=1;
my $db_source;
my $mod_name =$TANDEM . ".params";

#fill a hash with the search params!
my %search = FillSearchHash($pride_command);

#the results for each pep id
my $protein_hit; 
my %pep_result;
my %mgf_titles;
my $title = 0;
my @spectrum_number;
my $query = 0;
my %protein_group;
my @mod_pattern;

my $spectrum;

print "about to open the tandem result $TANDEM\n";
open(RES,"<$TANDEM") or print "WARNING: There is a problem opening the X!tandem results, $TANDEM\n";
 while(my $line = <RES>)
 {
  #get the mgf file name 
  if($line =~ m/\<bioml\s+xmlns\:GAML\=/)
  {
  my @tmp = split/label/,$line;
  #the file name should be after
  #the filename is in ''
  my @tmp2 = split/\'/,$tmp[1];
  #get the spectra numbers!
  #use the module from Parser.pm
  %mgf_titles = GetMGFTitles($tmp2[1]);
  }
  #is it a new query/group
  if($line =~ m/\<group\s+id/)
  {
  $query++;
 
  #find the spectrum number
  my @split1 = split/\s+/,$line;
  #only interested in the second element of the split1 array
  my @split2 = split/\=/,$split1[1];  
  $split2[1] =~ s/\"//g;
  $spectrum = $split2[1];

  $pep_count = 1;
  }
  #is it a protein?
  elsif($line =~ m/^\<protein/)
  {
  $found = 1;
  #get the protein information
  #split at whitespace
  my @prot = split/\s+/,$line;
   #for all the different fields, save protein info in a hash
   for(my $p=0 ; $p<scalar(@prot) ; $p++)
   {
   my @tmp = split/\=/,$prot[$p];
   #remove the "
   $tmp[1] =~ s/\"//g;
    if($tmp[0] eq "label")
    {
    #protein description
    $protein_details{$tmp[1]}[0] = "nan";
    $protein_details{$tmp[1]}[1] = $tmp[1];
    $protein_hit = $tmp[1];
    #there is no mass info
    #and now store the protein_group info
    push(@{$protein_group{$tmp[1]}},$query . "," . $pep_count);
    }
   }
  }
  elsif($line =~ m/^\<\/protein/)
  {
  $protein_hit = "";
  $found = 0;
  }
  elsif($line =~ m/^\<domain/)
  {
  #get the peptide info
  #split at the whitespace
  my @pep = split/\s+/,$line;
   #for all the different fields, save the peptide info to the hash
   for(my $p=0 ; $p<scalar(@pep) ; $p++)
   {
   my @tmp = split/\=/,$pep[$p];
   #remove the "
   $tmp[1] =~ s/\"//g;
   $tmp[1] =~ s/\>//g;
   $pep_result{$tmp[0]} = $tmp[1];
   }
  }
  elsif($line =~ m/^\<\/domain/)
  {
  my $modification_string;
   for(my $m=0 ; $m<scalar(@mod_pattern) ; $m++)
   {
   $modification_string .= $mod_pattern[$m];
   }
  my $ions_count = $pep_result{'y_ions'} + $pep_result{'b_ions'};

 $results[$spectrum][1]{'sequence'} = $pep_result{'seq'};
 $results[$spectrum][1]{'ionscore'} = $pep_result{'hyperscore'};
 $results[$spectrum][1]{'delta'} = $pep_result{'delta'};
 $results[$spectrum][1]{'start'} = $pep_result{'start'};
 $results[$spectrum][1]{'expect'} = $pep_result{'expect'};
 $results[$spectrum][1]{'protein'} = $protein_hit;
  %pep_result=();

  $spectrum="";
  $pep_count++;

  }

  elsif($line=~ m/\<group\s+type\=\"support\"\s+label\=\"fragment\s+ion\s+mass\s+spectrum\"/)
  {
  $title = 1;
  }
  elsif($line =~ m/\<note\s+label\=\"Description\"/ && $title == 1)
  {
  $title = 0;
  my @tmp = split/[<>]/,$line;
  $spectrum_number[$query] = $mgf_titles{$tmp[1]};
  }
 }
close FILE;

#and determine if they are n ter
my $db_path = GetDbPath($search{'DB'},$taxonomy);
my %nter = GetNterPeps($db_path);


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

sub GetMGFTitles ()
{
my $mgf_file = shift;
my %titles;
my $count=1;

#get all the title
my $name="";
my $mass;
my $charge;
open(MGF,"<$mgf_file") or print "unable to open the MGF file in Xtandem, $mgf_file\n";
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


sub GetDelta()
{
my $mass = shift;
my $number = 0;
my $found = 0;
my $delta = 0;


 foreach my $key (keys %search)
 {
  if($key =~ m/^delta\d/)
  {
   if($search{$key} =~m/\?\,$mass/)
   {
   $number = $key;
   $number =~ s/delta//;
   $found = 1;
    if($number>$delta)
    {
    $delta = $number;
    }
   }
  my $tmp = $key;
  $tmp =~ s/delta//;
   if($tmp>$delta)
   {
   $delta = $tmp;
   }  
  }
 }

 if($found == 0)
 {
 $delta++;
 my $tmp = "delta" . $delta;
 $number = $delta;
 }

return $number;

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


#and get the database version
my $db_version = GetDBVersion($search{"DB"});
$search{'release'} = $db_version;

return %search;
}


sub RunTandemMSMS()
{

my %data;
my $form;
my $result_file = shift;
my $tax_file = shift;
my $scoretype = shift;

 my $call = "/home/jensiepen/LIVERPOOL/bin/perl_scripts/GetDbInfo.pl " . $params{"DB"} . " tandem"; 
 system($call);
# $tax_file = "/fs/ls1/home/mjfssjs2/LIVERPOOL/tmp_files/taxonomy.xml";
my $tax_file = "/home/jensiepen/LIVERPOOL/taxonomy.xml";

#my $input_file = "/fs/ls1/home/mjfssjs2/LIVERPOOL/tmp_files/input.xml";
my $input_file = "/home/jensiepen/LIVERPOOL/input.xml"; 

#$params{'FILE'} = "/fs/ls1/home/mjfssjs2/Tandem_K/tandem-linux-08-02-01-1/bin/test_spectra.mgf";

#open the defaults  file (note will need to make this a unique file!)
open(NEW_INPUT,">$input_file") or print "WARNING:Unable to open the file $input_file  to write to\n";

#create the input file
print NEW_INPUT "<?xml version=\"1.0\"?>\n<bioml>\n";
print NEW_INPUT "<note type=\"input\" label=\"list path, default parameters\">/var/www/localhost/cgi-bin/pipeline/default_input.xml</note>\n";
print NEW_INPUT "<note type=\"input\" label=\"list path, taxonomy information\">$tax_file</note>\n";
#also add the file name and the output filename
print NEW_INPUT "<note type=\"input\" label=\"spectrum, path\">".$params{'FILE'} ."</note>\n";
print NEW_INPUT "<note type=\"input\" label=\"output, maximum valid expectation value\">0.01</note>\n";
#get the var mods - its an array
my @varmods = split/\,/,$params{"IT_MODS"};
my $var_mod_line = GetMods('0',$params{'MASS'},@varmods);
print NEW_INPUT "<note type=\"input\" label=\"refine, potential modification mass\">" . $var_mod_line . "</note>\n";
#nter mods
#get the fix mods - its an array
my @fixmods = split/\,/,$params{"MODS"};
my $nter_mod = GetMods('1',$params{'MASS'},@fixmods);
 if($nter_mod)
 {
 print NEW_INPUT "<note type=\"input\" label=\"protein, N-terminal residue modification mass\">" . $nter_mod . "</note>\n";
 }

my $fix_mod_line = GetMods('0',$params{'MASS'},@fixmods);
print NEW_INPUT "<note type=\"input\" label=\"residue, modification mass\">" . $fix_mod_line . "</note>\n";

my $iu;
my $pu;
 if($params{'ITOLU'} eq "Da")
 {
 $iu = "Daltons";
 }
 elsif($params{'ITOLU'} eq "ppm")
 {
 $iu = "ppm";
 }
 else
 {
 open(FILE,">$result_file");
 print FILE "Sorry you cannot use ". $params{'ITOLU'} . "(ion mass tolerance) with x!tandem, Daltons or ppm only<BR>";
 print "Sorry you cannot use ". $params{'ITOLU'} . "(ion mass tolerance) with x!tandem, Daltons or ppm only<BR>";
 close FILE;
 exit(1);
 }
 if($params{'TOLU'} eq "Da")
 {
 $pu = "Daltons";
 }
 elsif($params{'TOLU'} eq "ppm")
 {
 $pu = "ppm";
 }
 else
 {
 open(FILE,">$result_file");
 print FILE "Sorry you cannot use ". $params{'TOLU'} . " (peptide tolerance) with x!tandem, Daltons or ppm only<BR>";
 close FILE;
 exit(1);
 }

#make the mass type all lower case
 $params{'MASS'} = lc $params{'MASS'};

 print NEW_INPUT "<note type=\"input\" label=\"spectrum, fragment monoisotopic mass error\">" . $params{'ITOL'} . "</note>\n";
 print NEW_INPUT "<note type=\"input\" label=\"spectrum, parent " . $params{'MASS'} . " mass error plus\">" . $params{'TOL'} . "</note>\n";
 print NEW_INPUT "<note type=\"input\" label=\"spectrum, parent " . $params{'MASS'} . " mass error minus\">" . $params{'TOL'} . "</note>\n";
 print NEW_INPUT "<note type=\"input\" label=\"refine, maximum valid expectation value\">10</note>\n";
 print NEW_INPUT "<note type=\"input\" label=\"spectrum, parent " . $params{'MASS'} . " mass isotope error\">yes</note>\n";
 print NEW_INPUT "<note type=\"input\" label=\"spectrum, fragment " . $params{'MASS'} . " mass error units\">" . $params{'ITOLU'} . "</note>\n";
 print NEW_INPUT "<note type=\"input\" label=\"spectrum, parent " . $params{'MASS'} . " mass error units\">" . $params{'TOLU'} . "</note>\n";
 print NEW_INPUT "<note type=\"input\" label=\"spectrum, fragment mass type\">" . $params{'MASS'} . "</note>\n";
 #database stuff
 print NEW_INPUT "<note type=\"input\" label=\"protein, taxon\">" . $params{"DB"} . "<\/note>\n";
 #Enzyme stuff
 my $cle_string = GetEnzymeRule($params{'CLE'});
 print NEW_INPUT "<note type=\"input\" label=\"protein, cleavage site\">" . $cle_string . "<\/note>\n"; 
 print NEW_INPUT "<note type=\"input\" label=\"output, path\">" . $result_file . "</note>\n";
 print NEW_INPUT "<note type=\"input\" label=\"output, path hashing\">no</note>\n";
 print NEW_INPUT "<note type=\"input\" label=\"output, sort results by\">spectrum</note>\n";

 #and now if it is the kscore scoring add the following lines
 if($scoretype eq "kscore")
 {
 print NEW_INPUT "<note label=\"scoring, algorithm\" type=\"input\">k-score</note>\n";
 print NEW_INPUT "<note label=\"spectrum, use conditioning\" type=\"input\">no</note>\n";
 print NEW_INPUT "<note label=\"scoring, minimum ion count\" type=\"input\">1</note>\n";
 }

print NEW_INPUT "</bioml>\n\n";
close NEW_INPUT;

#run the tandem program
my $sys_call = "/fs/ls1/home/mjfssjs2/Tandem_K/tandem-linux-08-02-01-1/bin/tandem.exe " . $input_file ;


print "system call is $sys_call\n";
print "result_file is $result_file\n";

$result_file = $result_file . ".params";
system($sys_call);



print "resuklt_file $result_file created\n";

}


sub GetDBVersion()
{
my $db_name = shift;
my $taxonomy = shift;
my $version = "unknown" ;
my $found = 0;

 #open the tax file
 open(TAX,"<$taxonomy") or print "WARNING: Unable to open the taxonomy file for the db version - $taxonomy\n";
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
my $nter = shift;
my @mods = @_;
my $string="";

#get the mods from the mod_file
my %mascot_mods = GetModsFromModFile();

#the first mod is the mass type - i.e. monoisotopic or average
 for(my $m=1 ; $m<scalar(@mods) ; $m++)
 {

 #store the name of the mod
 $mod_params{'name'}[$m-1] = $mods[$m];

 #make sure there is no residue info at the end
  #ignore terminus for now
  #if not the nter bit
  if($nter == 0 && $mods[$m] !~ m/term/)
  {
   if($mascot_mods{$mods[$m]})
   {
   my @tmp = split/\s+/,$mascot_mods{$mods[$m]};
    if($mods[0] =~ m/monoisotopic/i)
    {
    #put a plus if it a positive number
    #minus the value of the amino acid
    $tmp[1] = $tmp[1] - $aa{$tmp[0]};

    $string = $string . ",+" . $tmp[1] . "@" . $tmp[0];
    #store the mass and residue
    $mod_params{'mass'}[$m-1] = $tmp[1];
    $mod_params{'residue'}[$m-1] =  $tmp[0];
    
    #and store this information in the search hash
    my $delta = "delta" . $m;
    $search{$delta} = $tmp[1] . "," . $mods[$m];
    }
    else
    {
    $tmp[2] = $tmp[2] - $aa{$tmp[0]};
    $string = $string . ",+" . $tmp[2] . "@" . $tmp[0];
    $mod_params{'mass'}[$m-1] = $tmp[2];
    $mod_params{'residue'}[$m-1] =  $tmp[0];

    #and store this information in the search hash
    my $delta = "delta" . $m;
    $search{$delta} = $tmp[2] . "," . $mods[$m];
    } 
   }
  } 
  #if it is for ntermods
  elsif($nter == 1 && $mods[$m] =~ m/N\-term/)
  {
  if($mascot_mods{$mods[$m]})
   {
   my @tmp = split/\s+/,$mascot_mods{$mods[$m]};
    if($mods[0] =~ m/monoisotopic/i)
    {
    $tmp[1] = $tmp[1]-$aa{$tmp[0]};

    $string = $string . ",+" . $tmp[1] . "@[";
    $mod_params{'mass'}[$m-1] = $tmp[1];
    $mod_params{'residue'}[$m-1] =  $tmp[0];
    my $delta = "delta" . $m;
    $search{$delta} = $tmp[1] . "," . $mods[$m];
    }
    else
    {
    $tmp[2] = $tmp[2] - $aa{$tmp[0]};
    $string = $string . ",+" . $tmp[2] . "@[";
    $mod_params{'mass'}[$m-1] = $tmp[2];
    $mod_params{'residue'}[$m-1] =  $tmp[0];
    my $delta = "delta" . $m;
    $search{$delta} = $tmp[2] . "," . $mods[$m];
    }
   }
  }
 }

#remove leading ,
$string =~ s/^\,//;

return $string;

}




sub GetEnzymeRule()
{

my $enzyme = shift;
my $string = "";

my $cle;
my $restrict;
my $term;
my $found = 0;

 if(!$enzyme)
 {
 $enzyme = "Trypsin";
 }
$enzyme =~ s/\+/\\\+/;
$enzyme =~ s/\-/\\\-/;
$enzyme =~ s/\_/\\\_/;
$enzyme =~ s/\//\\\//;

open(CLE,"</var/www/localhost/cgi-bin/pipeline/enzymes") or print "WARNING: unable to open the enzymes file\n";;
 while(my $line = <CLE>)
 {
 
  if($enzyme && $line =~ m/^Title:$enzyme$/)
  {
  $found = 1;
  }
  elsif($line =~ m/^Title/)
  {
  $found = 0;
  }
  elsif($found==1 && $line =~ m/^Cleavage/)
  {
  my @tmp = split/\:/,$line;
  $cle = $tmp[1];
  }
  elsif($found==1 && $line =~ m/^Restrict/)
  {
  my @tmp = split/\:/,$line;
  $restrict = $tmp[1];
  }
  elsif($found==1 && $line =~ m/^[NC]term/)
  {
  $term = $line;
  }
  elsif($found==1 && $line =~ m/^\*/)
  {
  $found = 0;
  }
 } 


#get rid of line break
$cle =~ s/\n$//;
$restrict =~ s/\n$//;
 #if there is a term
 #did we get where it cleaves?
 if($cle)
 {
  if($term)
  { 
   if($term =~ /Cterm/)
   {
   $string = $string . "[" . $cle . "]";
   }
  }
  else
  {
  $string = "[" . $cle . "]";
  }
 }
$string = $string . "|";

 #is there a restrict?
 if($restrict)
 {
 $string = $string . "{" . $restrict . "}";
 }



return $string;
}

sub GetModsFromModFile()
{

my %mods;

my $title;

open(MOD,"</var/www/localhost/cgi-bin/pipeline/mod_file") or print "WARNING: cannot open the mod file!!!\n";;

 while(my $line = <MOD>)
 {
 
  if($line =~ m/^Title/)
  {
  my @tmp = split/\:/,$line;
  $title = $tmp[1];
  $title =~ s/\s+$//;
  $title =~ s/\s/\*/;
  }
  elsif($line =~ m/^Residues/)
  {
  my @tmp = split/\:/,$line;
   if($title)
   {
   $mods{$title} = $tmp[1];
   }
 
  }
  elsif($line =~ m/term/)
  {
   $line =~ s/[\:\s]/\ /;
   if($title)
   {
   $mods{$title} = $line;
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

sub CheckDbForProtein()
{
my $db = shift;
my $type = 0;
open(FILE,"</home/jensiepen/LIVERPOOL/tandem_files/mascot.dat") or print "unable to open the file /home/jensiepen/LIVERPOOL/tandem_files/mascot.dat\n";

 while(my $line = <FILE>)
 {
  if($line =~ m/^$db\s+/)
  { 
  my @split = split/\s+/,$line; 
   if($split[2] eq "AA")
   {
   $type = 1;
   }
   else
   {
   $type = 0;
   }
  }
 }
close FILE;

return $type;

}

sub UpdateTaxonomy()
{
my $file = shift;

# $command = "cp /var/www/localhost/cgi-bin/pipeline/ISPIDER_DEFAULT_TAXONOMY_FILE.xml " . $file;

#system($command); 

}

#I have set a crojob to get the current file names - so this is old
sub OLDUpdateTaxonomy()
{

my $file = shift;

my $db_counter = 0;
my $start = 0;
my @split_line;
open(TAX,">$file") or print "WARNING: unable to open the taxonomy file to write to ($file)\n";;
print TAX "<?xml version=\"1.0\"?>\n<bioml label=\"x! taxon-to-file matching list\">\n";

#my $sys = "scp msct:/usr/local/share/mascot/config/mascot.dat /tmp/mascot.dat";
#system($sys);

if(! -e "/home/jensiepen/LIVERPOOL/tandem_files/mascot.dat")
{
my $sys = "cp /var/www/localhost/cgi-bin/pipeline/mascot.dat /home/jensiepen/LIVERPOOL/tandem_files/mascot.dat";
system($sys);
}
open(DBFILE,"</home/jensiepen/LIVERPOOL/tandem_files/mascot.dat") or print "unable to open the mascot.dat file\n";
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
   #$split_line[1] =~ s/\/usr\/share\//\/db\//;
   $split_line[1] =~ s/\/usr\/local\/share\/mascot\/sequence/\/fs\/msct\/sequence/;
   $filename = GetFileName($split_line[1]);
   }
   else
   {
   $filename = $split_line[1];
   #$filename =~ s/\/usr\/share\//\/db\//;
   $split_line[1] =~ s/\/usr\/local\/share\/mascot\/sequence/\/fs\/msct\/sequence/;
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

#add the reverse IPI database 
print TAX "<taxon label=\"IPI_arath\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_arath/current/ipi.ARATH.v3_23.fasta\" />\n";
print TAX "</taxon>\n";
print TAX "<taxon label=\"IPI_arathr\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_arath/current/reverse/ipi.ARATH.v3_23.reverse\" />\n";
print TAX "</taxon>\n";

print TAX "<taxon label=\"IPI_bovin\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_bovin/current/ipi.BOVIN.v3_11.fasta\" />\n";
print TAX "</taxon>\n";
print TAX "<taxon label=\"IPI_bovinr\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_bovin/current/reverse/ipi.BOVIN.v3_11.reverse\" />\n";
print TAX "</taxon>\n";

print TAX "<taxon label=\"IPI_brare\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_brare/current/ipi.BRARE.v3_24.fasta\" />\n";
print TAX "</taxon>\n";
print TAX "<taxon label=\"IPI_brarer\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_brare/current/reverse/ipi.BRARE.v3_24.reverse\" />\n";
print TAX "</taxon>\n";

print TAX "<taxon label=\"IPI_chick\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_chick/current/ipi.CHICK.v3_19.fasta\" />\n";
print TAX "</taxon>\n";
print TAX "<taxon label=\"IPI_chickr\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_chick/current/reverse/ipi.CHICK.v3_19.reverse\" />\n";
print TAX "</taxon>\n";

print TAX "<taxon label=\"IPI_human\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_human/current/ipi.HUMAN.fasta\" />\n";
print TAX "</taxon>\n";
print TAX "<taxon label=\"IPI_humanr\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_human/current/reverse/ipi.HUMAN.reverse\" />\n";
print TAX "</taxon>\n";

print TAX "<taxon label=\"IPI_mouse\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_mouse/current/ipi.MOUSE.v3_25.fasta\" />\n";
print TAX "</taxon>\n";
print TAX "<taxon label=\"IPI_mouser\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_mouse/current/reverse/ipi.MOUSE.v3_25.reverse\" />\n";
print TAX "</taxon>\n";

print TAX "<taxon label=\"IPI_rat\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_rat/current/ipi.RAT.v3_25.fasta\" />\n";
print TAX "</taxon>\n";
print TAX "<taxon label=\"IPI_ratr\">\n";
print TAX "\t<file format=\"peptide\" URL=\"/fs/msct/sequence/IPI_rat/current/reverse/ipi.RAT.v3_25.reverse\" />\n";
print TAX "</taxon>\n";


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


sub GetDbPath
{

my $dbname = shift;
my $taxonomy = shift;
my $path;
#open the taxonomy xml file
open(FILE,"<$taxonomy") or print "unable to open the file, $taxonomy\n";

 while(my $line = <FILE>)
 {
  if($line =~ m/\<file\s+format/)
  {
  my @split = split/\=/,$line;
  $path = $split[2];
  $path =~ s/\"//g;
  $path =~ s/\>//g;
  $path =~ s/\s+$//g;
  $path =~ s/\/$//;
  }
 }
close FILE;

return $path;

}



