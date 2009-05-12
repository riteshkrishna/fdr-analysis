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
use lib qw(/var/www/localhost/cgi-bin/FDRAnalysis/perl_modules/);
use RunMascot;
use FDRAnalysisMascotParser;
use FDRAnalysisXTandem;
use FDRAnalysisOmssa;
use DecoyMascotParser;
use DeltaMassPlot;
use GygiRankPlot;
use GetFileStats;
use GetFDRPlot;
use GetFDRValues;
use GetDiscriminantScorePlot;
use ZoomDiscriminantScorePlot;
use GetNterminaldistribution;
use CreatePNG;
use Getopt::Std;
use GetFileStatsSeperateFiles;
use List::Util qw[min max];


########################################################################
####  An uber Mascot Parse to parse all types of decoy mascot searces ##
########################################################################

our($opt_w,$opt_k,$opt_c, $opt_m, $opt_f, $opt_d, $opt_r, $opt_q, $opt_s, $opt_o, $opt_n, $opt_h, $opt_l, $opt_t, $opt_p, $opt_b, $opt_x, $opt_e, $opt_j, $opt_y, $opt_a, $opt_z, $opt_g, $opt_u, $opt_v, $opt_i);
getopts('c:m:f:d:r:q:s:o:n:ht:l:p:b:x:e:j:y:a:zg:u:v:i:k:w:');

   if($opt_h)
   {
   ErrorMsg('2');
   exit(1);
   }
   elsif(!$opt_c && !$opt_f && !$opt_d && !$opt_x && !$opt_r && !$opt_l)
   {
   ErrorMsg('1');
   exit(1);
   }
   

 #if a decoy method what is reverse tag?
 if(!$opt_t)
 {
 $opt_t = "REV_";
 }

#global structures
my @images;
my %gygi_results;
my %gygi_peptides;
my %jones_results;
my %jones_peptides;
my $jones_combined_file = $opt_s . "tmp_combined.txt";
my $cmd = "rm $jones_combined_file";
my %investigate_proteins;
my %investigate_proteins_jones;
my %summary_table;
my %peptidelist;
system($cmd);


 # Mascot Concatanated?              
 if($opt_c)
 {
 MascotConcat($opt_c);
 }

 # Mascot seperate forward/reverse
 if($opt_f && $opt_d)
 {
 AnalyseSeperate('M');
 }
 elsif($opt_d)
 {
 MascotDecoy();
 }

 #list of concats
 # best way to handle this?
 # run each file individually and generate a file for each and at the same time get a unique peptide list? or just look at overlap?
 if($opt_l)
 {
 ConcatList('M',$opt_l);
 }


 #X!Tandem
 if($opt_x)
 {
 TandemConcat($opt_x);
 }
 if($opt_u && $opt_v)
 {
#
 AnalyseSeperate('T');
 } 

 #Omssa
 if($opt_r)
 {
 OmssaConcat($opt_r)
 }
 if($opt_p && $opt_i)
 {
 AnalyseSeperate('O');
 }

if($opt_g)
{
AnalyseProteins('K');
AnalyseProteins('E');
}

RunAndysStuff();





sub AnalyseProteins
{
my $type = shift;

my %proteins;
 if($type eq "K")
 {
 %proteins = %investigate_proteins_jones;
 }
 else
 { 
 %proteins = %investigate_proteins;
 }


my $count = 0;
my $max = 0;;
my $max2 = 0;


my $fdr = 0.01;

 #for all the proteins
 foreach my $prot (keys %proteins)
 {

 my %peps;
  #for all potential matches 
  if($proteins{$prot}{$fdr}[0])
  {
  for(my $i=0 ; $i<scalar(@{$proteins{$prot}{$fdr}}) ; $i++)
   { 
   $peps{$proteins{$prot}{$fdr}[$i]{'sequence'}}[0] = $proteins{$prot}{$fdr}[$i]{'start'};
   $peps{$proteins{$prot}{$fdr}[$i]{'sequence'}}[1] = $proteins{$prot}{$fdr}[$i]{'ionscore'};
   $peps{$proteins{$prot}{$fdr}[$i]{'sequence'}}[2] = $proteins{$prot}{$fdr}[$i]{'original'};
   }
  }

 # now for all the peptides if there is more than 1
 my $count = 0;
 my $diff = 0;
 my %old;
  foreach my $p (keys %peps)
  {  
  $count++;
   if(%old)
   {    
    if($old{$p} && $old{$p} != $peps{$p}[0])
    {
    $diff++;
    }  
   $old{$p} = $peps{$p}[0];
   }
  }
  if($diff>0)
  {
  #print "$type method at $fdr fdr $prot has both Nter and internal peps\n";
   foreach my $p (keys %peps)
   {
  # print "$p\t$peps{$p}[0]\t$peps{$p}[1]\t$peps{$p}[2]\n";
   }
  }
  elsif($count>1)
  {
  #print "$type at $fdr fdr $prot has more than 1 peptide identified\n\n";
   foreach my $p (keys %peps)
   {
   #print "$p\t$peps{$p}[0]\t$peps{$p}[1]\t$peps{$p}[2]\n";
   } 
  }
 }

}



sub AddToPeptideList
{
my $se = shift;
my $file = shift;
my $ptr = shift;
my %res = %{$ptr};


$peptidelist{'file'}{'file'}{'file'} = $opt_s . "_peptidelist.txt";
 foreach my $fdr (keys %res)
 {
 $peptidelist{$se}{'E'}{$fdr} = $res{$fdr}[4];
 $peptidelist{$se}{'K'}{$fdr} = $res{$fdr}[5];
 }
}

sub AddToCreateSummaryTable
{

my $se = shift;
my $file = shift;
my $ptr = shift;
my %res = %{$ptr};
my $ptr2 = shift;
my %jones = %{$ptr2};

$summary_table{'file'}{'file'}{'0'}[0] = $opt_s . "_summary.txt";


 #for all the gygi results 
 foreach my $fdr (sort keys %res)
 {
  if($fdr<=0.5)
  {
  $summary_table{$se}{'E'}{$fdr}[0] = $res{$fdr}[0];

   if($opt_n)
   {
   $summary_table{$se}{'E'}{$fdr}[1] = $res{$fdr}[1];
   }
  }
 }
 foreach my $fdr (sort keys %jones)
 {
  if($fdr<=0.5)
  {
  $summary_table{$se}{'K'}{$fdr}[0] = $jones{$fdr}[0];
   if($opt_n)
   {
   $summary_table{$se}{'K'}{$fdr}[1] = $jones{$fdr}[1];
   }
  }
 }


}


sub CreatePeptideList
{

#open the peptide file
open(NEW,">$peptidelist{'file'}{'file'}{'file'}") or die "unable to open the file to write to $peptidelist{'file'}{'file'} to write to\n";
print NEW "Search Engine\tProtein\tPeptide\tStart position\tScore\tExpect\n";
 foreach my $se (keys %peptidelist)
 {
 my $engine;
  if($se eq "M")
  {
  $engine = "Mascot";
  }
  elsif($se eq "T")
  {
  $engine = "XTandem";
  }
  elsif($se eq "O")
  {
  $engine = "Omssa";
  }
  elsif($se eq "C")
  {
  $engine = "Combined Search engines";
  }

  foreach my $fdr (keys %{$peptidelist{$se}{'E'}})
  {
   if($fdr<=$opt_w)
   {
   #get the peptides
   my @split = split/\*\*\*/,$peptidelist{$se}{'E'}{$fdr};
   #for all the peptides
   for(my $i=0 ; $i<scalar(@split) ; $i++)
   {
   my @list = split/\#\#\#/,$split[$i];
   print NEW "$engine\t$list[0]\t$list[1]\t$list[2]\t$list[3]\t$list[4]\n";
   }
  } 
 }
}

close NEW;

}

sub CreateSummary
{


open(NEW,">$summary_table{'file'}{'file'}{'0'}[0]") or die "**unable to open the file $summary_table{'file'}{'file'}{'0'}[0] to write to\n";
#results for each file

print NEW "<TABLE BORDER='1'>";
print NEW "<TR><TH>Search Engine</TH><TH>FDR method</TH><TH>FDR value</TH><TH>Total unique peptide count</TH>";
 if($opt_n)
 {
 print NEW "<TH>unique Nter count</TH>";
 }
print NEW "</TR>";

 #for all search engines
 foreach my $se ( keys %summary_table)
 {

 my $engine;
  if($se eq "M")
  {
  $engine = "Mascot";
  }
  elsif($se eq "T")
  {
  $engine = "XTandem";
  }
  elsif($se eq "O")
  {
  $engine = "Omssa";
  }
  elsif($se eq "C")
  {
  $engine = "**Combined Search engines";
  }

  if($se ne "file")
  {
   #for the different FDR types
   foreach my $fdrtype (keys %{$summary_table{$se}})
   {
   my $bgcolor = " ";
    if($fdrtype eq "E")
    {
    $bgcolor = qq{ bgcolor="#cccccc" };
    }
    elsif($se eq "C")
    {
    $bgcolor = qq{ bgcolor="#669966" };
    print NEW "<TR><TD><BR></TD><TD><BR></TD><TD><BR></TD><TD><BR></TD>";
      if($opt_n)
      {
      print NEW "<TD><BR></TD>";
      }

    }

    #for the different FDR values
    foreach my $fdr (sort keys %{$summary_table{$se}{$fdrtype}})
    {
     if($fdr<=$opt_w)
     {
      if($summary_table{$se}{$fdrtype}{$fdr}[0])
      {
      print NEW "<TR><TD $bgcolor>$engine</TD><TD $bgcolor>$fdrtype</TD><TD $bgcolor>$fdr</TD><TD $bgcolor>$summary_table{$se}{$fdrtype}{$fdr}[0]</TD>";
      }
      else
      {
      print NEW "<TR><TD $bgcolor>$engine</TD><TD $bgcolor>$fdrtype</TD><TD $bgcolor>$fdr</TD><TD $bgcolor>0</TD>";
      }
      if($opt_n)
      {
       if($summary_table{$se}{$fdrtype}{$fdr}[1])
       {
       print NEW "<TD $bgcolor>$summary_table{$se}{$fdrtype}{$fdr}[1]</TD>";
       }
       else
       {
       print NEW "<TD $bgcolor>0</TD>";
       }
      }
     print NEW "</TR>";
     }
    }

   }
  } 
 }

print NEW "</TABLE>";

print "file written ... \n";

}



sub RunAndysStuff
{
 #remove the old files
 my $cmd = "rm /var/www/tmp/" . $opt_s . "combined*";
 system($cmd);


 #run andys!
 if($opt_z && $jones_combined_file)
 {
 #remove the files
 my $cmd = "rm " . $opt_s . "*combined*out";
system($cmd);

print "about to run the command\n";


#my $cmd = "/fs/san/home/mjfssjs2/LIVERPOOL/bin/perl_scripts/andys_MS-Weighted.pl ";
my $cmd = "perl /fs/san/home/mjfssjs2/FDRScore_March09/MultipleSearch/web_FDRScore.pl ";
$cmd .= "-F $jones_combined_file -R 1 -I $opt_t -O ". $opt_s . "combined_results.out -P " . $opt_s . "combined_proteins.out -Q " . $opt_s . "combined_peptides.out -T " . $opt_w; 

system($cmd);

print "andys command is $cmd\n";

 ParseCombined();
 CreateSummary();
 CreatePeptideList();

 #so the combined peptides file is downloadable move it to the other tmp dir
 my $cmd = "cp " . $opt_s . "combined_peptides.out /var/www/localhost/htdocs/FDRAnalysis/tmp/";
 system($cmd);


 my $imagefile = $opt_s . "_VennDiagram.png";
 $imagefile =~ s/\/var\/www\/tmp\//\/var\/www\/localhost\/htdocs\/FDRAnalysis\/tmp\//;

 } 
 else
 {
print "About to create the peptidelist and summary\n";
 CreateSummary();
 CreatePeptideList();
 } 

#also want to create the venn diagram for the peptide
my @values = GetVennNumbers();
my $imagefile = $opt_s . "_VennDiagram.png";
$imagefile =~ s/\/var\/www\/tmp\//\/var\/www\/localhost\/htdocs\/FDRAnalysis\/tmp\//;

 if(!GetPNGVenn($imagefile,@values))
 {
 print "Error trying to run the GetPNGVenn\n";
 exit(1);
 }

print "Venn diagram should have been created\n";


}

sub GetVennNumbers
{


#peptidelist
my $newpeplist = $opt_s . "_peptidelist.txt";

my %peps;
my @res;
my %se;

 open(FILE,"<$newpeplist") or die "unable to open the file to read $newpeplist\n";
  while(my $line = <FILE>)
  {
   if($line !~ m/^Search/)
   {
   my @split = split/\t/,$line;
   $peps{$split[2]}{$split[0]} = 1;
   $se{$split[0]} = 1;
   }
 }
 close FILE;

my %tmplabels;

 my $count = 1;
 my @labels;

 foreach my $engine (keys %se)
 {
 $labels[$count] = $engine;
 $res[0][$count] = $engine;
 $count++;
 }


 #and now get the results
 foreach my $seq (keys %peps)
 {
  #are all 3 present?
  if($peps{$seq}{$labels[1]} && $peps{$seq}{$labels[2]} && $peps{$seq}{$labels[3]})
  {
   if($res[0][0])
   {
   $res[0][0]++;
   }
   else
   {
   $res[0][0] = 1;
   }
  }
  #otherwise if 1 and 2 present
  elsif($peps{$seq}{$labels[1]} && $peps{$seq}{$labels[2]})
  {
   if($res[1][2])
   {
   $res[1][2]++;
   }
   else
   {
   $res[1][2] = 1;
   }
  }
  #otherwise 2 and 3 presant
  elsif($peps{$seq}{$labels[2]} && $peps{$seq}{$labels[3]})
  {
   if($res[2][3])
   {
   $res[2][3]++;
   }
   else
   {
   $res[2][3] = 1;
   }
  }
  #otherwise 1 and 3
  elsif($peps{$seq}{$labels[1]} && $peps{$seq}{$labels[3]})
  {
   if($res[1][3])
   {
   $res[1][3]++;
   }
   else
   {
   $res[1][3] =1;
   }
  }
  #otherwise just 1
  elsif($peps{$seq}{$labels[1]})
  {
   if($res[1][0])
   {
   $res[1][0]++;
   }
   else
   {
   $res[1][0] = 1;
   }
  }
  #otherwise just 2
  elsif($peps{$seq}{$labels[2]})
  {
   if($res[2][0])
   {
   $res[2][0]++;
   }
   else
   {
   $res[2][0] = 1;
   }
  }
  #otherwise just 3
  elsif($peps{$seq}{$labels[3]})
  {
   if($res[3][0])
   {
   $res[3][0]++;
   }
   else
   {
   $res[3][0] = 1;
   }
  }
  else
  {
print "missed sequence is $seq\n";
   if($res[3][3])
   {
   $res[3][3]++;
   }
   else
   {
   $res[3][3] = 1;
   }
  }
 }

#remove the file

return @res;

}


 sub ParseCombined
 {

 my %res;
 my $file = $opt_s . "combined_peptides.out";
print "looking at file $file\n";

#open the peptide results
 open(FILE,"<$file") or print "problem with Combined results ($file)\n";
  while(my $line = <FILE>)
  {
  my @split = split/\t/,$line;
  $res{$split[0]}{$split[1]} = $split[2];
  print "adding res{$split[0]}{$split[1]} = $split[2]\n";
  }
 close FILE;

my %observed;
my $count = 0;
my $nter = 0;
#for the different FDRs
 foreach my $fdr (keys %res)
 {
 print "fdr is $fdr and about to see if $fdr >= $opt_w\n";
  if($fdr>=$opt_w)
  {

  print "The FDR $fdr matched\n";
   foreach my $seq (keys %{$res{$fdr}})
   {
    if(!$observed{$seq})
    {
print "LOOK $seq\n";
$observed{$seq} = 1;
    $count++;
     if($res{$fdr}{$seq}<3)
     {
     $nter++;
     }
    } 
   }
  } 
 } 
 print "count is $count\n";
 $summary_table{'C'}{'K'}{$opt_w}[0] = $count;
 $summary_table{'C'}{'K'}{$opt_w}[1] = $nter;

 }

 sub MascotDecoy
 {

 my $mascotfile = $opt_d;
 my $db = GetDbPath($mascotfile);
 my $mgf = "";
 my $ptr = ParseDecoyMascot($mascotfile,$mgf,$db,'H','0.05');
 my @resultsTwo = @{$ptr};

 my $ptr = ParseMascot($mascotfile,$mgf,$db,'H','0.05');
 my @results = @{$ptr};


 my %peptide_hits;

 my $fdrtype = 0;
 my $setype = "M";
 my $max  = max(scalar(@results),scalar(@resultsTwo));

 my @new_results = ProduceUber($max,\@results,\@resultsTwo);


my $n = 0;
my $in = 0;
my $rn = 0;
my $rin = 0;

for(my $s = 0 ; $s<scalar(@new_results) ; $s++)
{
 for(my $r=2 ; $r<11 ; $r++) 
 {
  if($new_results[$s][$r]{'sequence'})
  {
   if($new_results[$s][$r]{'protein'} =~ m/^REV\_/)
   {
    if($new_results[$s][$r]{'start'}<3)
    {
    $rn++;
    }
    else
    {
    $rin++;
    }
   }
   elsif($new_results[$s][$r]{'start'}<3)
   {
   $n++;
   }
   else
   {
   $in++;
   }
  }
 }
}
print "for nter is $n for internal is $in revese internal is $rin reverse nter is $rn\n";

 my @image;
my $fdrtype = 0;
  if($opt_q)
  {
  $fdrtype = $opt_q;
  }

 my $setype = "M";

 my $scoretype = "S";

 my $imagefile = $opt_f;
 $imagefile =~ s/\.dat//;
 $imagefile = $imagefile . "_" . $scoretype . ".png";
  if($opt_o)
  {
  $imagefile = $opt_o;
  }

 #Get the actual FDR values for both Gyig and Jones method
 my %all_results = GetFDRValues($opt_w,'0','M',$opt_n,$opt_t,@new_results);
 
 my %gygi_results;
 my %jones_results;
 
 #put the different results into the respective hashes
  foreach my $fdr (keys %all_results)
  {
  $gygi_results{$fdr}[0] = $all_results{$fdr}[0];
  $gygi_results{$fdr}[1] = $all_results{$fdr}[1];
  $jones_results{$fdr}[0] = $all_results{$fdr}[2];
  $jones_results{$fdr}[1] = $all_results{$fdr}[3];
  }



 #summary file
 my $plot = $imagefile;
 $plot =~ s/\.png/Summary\.txt/;
 AddToCreateSummaryTable('M',$plot,\%gygi_results,\%jones_results);
 AddToPeptideList('M',$plot,\%all_results);

 #deltamass plot
 my $plot = $imagefile;
 $plot =~ s/\.png/DeltaMass\.png/;
 GetDeltaMassPlot($plot,$scoretype,$setype,@new_results);

 #gygi rank plot
 my $plot = $imagefile;
 $plot =~ s/\.png/GygiRank\.png/;
 GetGygiRankPlot($setype,$plot,$opt_k,@new_results);

 #Score distribution
 my $plot = $imagefile;
 $plot =~ s/\.png/ScoreDist\.png/;
 GetScoreDistribution($plot,$scoretype,$setype,@new_results);

 #Zoom score plot
 my $plot = $imagefile;
 $plot =~ s/\.png/ZoomScore\.png/;
 ZoomScoreDistribution($plot,$scoretype,$setype,@new_results);


  if($opt_n)
  {
  #nter plot
  my $plot = $imagefile;
  $plot =~ s/\.png/NterDist\.png/;
  GetNtermPlot($setype,$plot,@new_results);
  }

 }


 sub AnalyseSeperate
 {
 my $setype = shift;
 my @results;
 my @resultsTwo;

 my $imagefile;

  if($setype eq "M")
  {
  my $mascotfile = $opt_f;
  my $filetwo = $opt_d;
  my $mgf = $opt_m;
  my $db = GetDbPath($mascotfile);

  my $ptr = ParseMascot($mascotfile,$mgf,$db,'H','0.05');
  @results = @{$ptr};

  $db = GetDbPath($filetwo);
  #get the second result file
  my $ptr2 = ParseMascot($filetwo,$mgf,$db,'H','0.05');
  @resultsTwo = @{$ptr2};
 
  $imagefile = $opt_o;
  }

  if($setype eq "T")
  {
  my $forward = $opt_u;
  my $reverse = $opt_v;
  my $ptr = ParseTandem($forward,'/fs/san/home/mjfssjs2/input.xml','/fs/san/home/mjfssjs2/taxonomy.xml');
  @results = @{$ptr};
 
  my $ptr2 = ParseTandem($reverse,'/fs/san/home/mjfssjs2/input.xml','/fs/san/home/mjfssjs2/taxonomy.xml');
  @resultsTwo = @{$ptr2};

  $imagefile = $opt_e;
  }
 
  if($setype eq "O")
  {
  my $forward = $opt_p;
  my $reverse = $opt_i;

  my $ptr = ParseOmssa($forward,"");
  @results = @{$ptr};

  my $ptr2 = ParseOmssa($reverse,"");
  @resultsTwo = @{$ptr2};
 
  $imagefile = $opt_j;
  }

 my %peptide_hits;

 my $fdrtype = 0;
 my $max  = max(scalar(@results),scalar(@resultsTwo));
 my @new_results = ProduceUber($max,\@results,\@resultsTwo);


my $n = 0;
my $in = 0;
my $rn = 0;
my $rin = 0;


  #if combining the result
  if($opt_z && $setype eq "O")
  {
  SimpleCount('OMSSA',@new_results);
  }
  elsif($opt_z && $setype eq "M")
  {
  SimpleCount('Mascot',@new_results);
  } 
  elsif($opt_z && $setype eq "T")
  {
  SimpleCount('Tandem',@new_results);
  }


 my @image;
my $fdrtype = 0;
  if($opt_q)
  {
  $fdrtype = $opt_q;
  }


 my $scoretype = "S";

 #Get the actual FDR values for both Gyig and Jones method
 my %all_results = GetFDRValues($opt_w,'0','M',$opt_n,$opt_t,@new_results);
 
 my %gygi_results;
 my %jones_results;

 #put the different results into the respective hashes
  foreach my $fdr (keys %all_results)
  {
  $gygi_results{$fdr}[0] = $all_results{$fdr}[0];
  $gygi_results{$fdr}[1] = $all_results{$fdr}[1];
  $jones_results{$fdr}[0] = $all_results{$fdr}[2];
  $jones_results{$fdr}[1] = $all_results{$fdr}[3];
  }



#summary file
my $plot = $imagefile;
$plot =~ s/\.png/Summary\.txt/;
AddToCreateSummaryTable($setype,$plot,\%gygi_results,\%jones_results);
AddToPeptideList($setype,$plot,\%all_results);

#deltamass plot
my $plot = $imagefile;
$plot =~ s/\.png/DeltaMass\.png/;
GetDeltaMassPlot($plot,$scoretype,$setype,@new_results);

#gygi rank plot
my $plot = $imagefile;
$plot =~ s/\.png/GygiRank\.png/;
GetGygiRankPlot($setype,$plot,$opt_k,$opt_k,@new_results);

#Score distribution
my $plot = $imagefile;
$plot =~ s/\.png/ScoreDist\.png/;
GetScoreDistribution($plot,$scoretype,$setype,@new_results);

#Zoom score plot
my $plot = $imagefile;
$plot =~ s/\.png/ZoomScore\.png/;
ZoomScoreDistribution($plot,$scoretype,$setype,@new_results);

  if($opt_n)
  {
#nter plot
my $plot = $imagefile;
$plot =~ s/\.png/NterDist\.png/;
GetNtermPlot($setype,$plot,@new_results);
  }

 }


 sub ProduceUber
 {
 my $max = shift;
 my $ptr = shift;
 my $ptr2 = shift;
 my @results = @{$ptr};
 my @resultsTwo = @{$ptr2};

 my @new_results;

my $count = 0;
 for(my $r=0 ; $r<$max ; $r++)
 {
 #find the maximum hit
 my $limit = 10000;
 my $rank = 0;
 #rev get best hit (do rev first)
  for(my $k=1 ; $k<11 ; $k++)
  {
   #is there a reverse hit?
   if($resultsTwo[$r][$k]{'sequence'})
   {
    if(!$new_results[$r][$rank+1]{'ionscore'})
    {
    $rank++;
     if($rank<11)
     {
      foreach my $attr (keys %{$resultsTwo[$r][$k]})
      {
      $new_results[$r][$rank]{$attr} = $resultsTwo[$r][$k]{$attr};
      }
     }

    }
    elsif($resultsTwo[$r][$k]{'ionscore'}>$new_results[$r][$rank+1]{'ionscore'})
    {
    $rank++;
    #push reverse down the list
     if($rank<10)
     {
      foreach my $attr (keys %{$new_results[$r][$rank]})
      {
      $new_results[$r][$rank+1]{$attr} = $new_results[$r][$rank]{$attr};
      }
     }
     if($rank<11)
     {
      foreach my $attr (keys %{$resultsTwo[$r][$k]})
      {
      $new_results[$r][$rank]{$attr} = $resultsTwo[$r][$k]{$attr};
      }
     }
    }

    elsif($resultsTwo[$r][$k]{'ionscore'}<$new_results[$r][$rank+1]{'ionscore'})
    {
    $rank++;
    my $limit = GetMax(\@new_results,$r);
     for(my $position = $rank+1 ; $position<=$limit ; $position++)
     {
      if($resultsTwo[$r][$k]{'ionscore'}>$new_results[$r][$rank+1]{'ionscore'})
      {
       for(my $push= 9 ; $push>=$position ; $push--)
       {
        foreach my $attr (keys %{$new_results[$r][$push]})
        {
        $new_results[$r][$push+1]{$attr} = $new_results[$r][$push]{$attr};
        }
       }
       foreach my $attr (keys %{$results[$r][$k]})
       {
       $new_results[$r][$position]{$attr} = $resultsTwo[$r][$k]{$attr};
       }
      $position = 11;
      }
      else
      {
       foreach my $attr (keys %{$results[$r][$k]})
       {
       $new_results[$r][$rank+2]{$attr} = $resultsTwo[$r][$k]{$attr};
       }
      }
     } 
    
    }
   }#end of is there a reverse hit

   if($results[$r][$k]{'sequence'} && $results[$r][$k]{'sequence'} ne "NULL")
   {
   #if the forward hit does better than reverse for this particular rank
    if($results[$r][$k]{'ionscore'}>$new_results[$r][$rank]{'ionscore'} && $new_results[$r][$rank]{'ionscore'})
    {
   #push reverse down the list
     for(my $push= 9 ; $push>($rank-1) ; $push--)
     {
      if($new_results[$r][$push]{'sequence'})
      {
       foreach my $attr (keys %{$new_results[$r][$push]})
       {
       $new_results[$r][$push+1]{$attr} = $new_results[$r][$push]{$attr};
       }
      }
     }
    my @split = split/\|/,$results[$r][$k]{'protein'};
    $results[$r][$k]{'protein'} = $split[0];
     if($rank<11)
     {
      foreach my $attr (keys %{$results[$r][$k]})
      {
      $new_results[$r][$rank]{$attr} = $results[$r][$k]{$attr};
      }
     }
    }   
    elsif($results[$r][$k]{'ionscore'}>$new_results[$r][$rank]{'ionscore'})
    {
    $rank++;
    my @split = split/\|/,$results[$r][$k]{'protein'};
    $results[$r][$k]{'protein'} = $split[0];
     if($rank<11)
     {   
      foreach my $attr (keys %{$results[$r][$k]})
      {
      $new_results[$r][$rank]{$attr} = $results[$r][$k]{$attr};
      }
     }

    }
    elsif($results[$r][$k]{'ionscore'} == $new_results[$r][$rank]{'ionscore'})
    {
    $rank++;
     if($rank<11)
     {
      foreach my $attr (keys %{$results[$r][$k]})
      {
      $new_results[$r][$rank]{$attr} = $results[$r][$k]{$attr};

      }
     }
    }
    #otherwise push this into the next position
    else
    {
    #where does this score fit?
     #for all the positions
     for(my $position = $rank+1 ; $position<11 ; $position++)
     {
      if($results[$r][$k]{'ionscore'}>$new_results[$r][$position]{'ionscore'} || !$new_results[$r][$position]{'ionscore'} )
      {
       for(my $push= 9 ; $push>=$position ; $push--)
       {
        foreach my $attr (keys %{$results[$r][$k]})
        {
        $new_results[$r][$push+1]{$attr} = $new_results[$r][$push]{$attr};
        }
       }
       foreach my $attr (keys %{$results[$r][$k]})
       {
       $new_results[$r][$position]{$attr} = $results[$r][$k]{$attr};
       }
       $position = 11;
       }
      }

     }
    }
   } #end all rev ranks
  } #end all spectra



 return @new_results;
 }


 sub OmssaConcat
 {
 my $omssafile = shift;

 my($ptr1) = ParseOmssa($omssafile,$opt_a);
 my @results = @{$ptr1};
  #if combining the result
  if($opt_z)
  { 
  SimpleCount('OMSSA',@results);
  }

 my $scoretype = 'S';
 my $fdrtype = 5;
 my $setype = "O";

 my $imagefile = $opt_j;

 #Get the actual FDR values for both Gyig and Jones method
 my %all_results = GetFDRValues($opt_w,'0','O',$opt_n,$opt_t,@results);
 
 my %gygi_results;
 my %jones_results;

 #put the different results into the respective hashes
  foreach my $fdr (keys %all_results)
  {
  $gygi_results{$fdr}[0] = $all_results{$fdr}[0];
  $gygi_results{$fdr}[1] = $all_results{$fdr}[1];
  $jones_results{$fdr}[0] = $all_results{$fdr}[2];
  $jones_results{$fdr}[1] = $all_results{$fdr}[3];
  }

 

 #summary file
 my $plot = $imagefile;
 $plot =~ s/\.png/Summary\.txt/;
 AddToCreateSummaryTable('O',$plot,\%gygi_results,\%jones_results);
 AddToPeptideList('O',$plot,\%all_results);

 #deltamass plot
 my $plot = $imagefile;
 $plot =~ s/\.png/DeltaMass\.png/;
 GetDeltaMassPlot($plot,$scoretype,$setype,@results);

 #gygi rank plot
 my $plot = $imagefile;
 $plot =~ s/\.png/GygiRank\.png/;
 GetGygiRankPlot($setype,$plot,$opt_k,@results);

 #Score distribution
 my $plot = $imagefile;
 $plot =~ s/\.png/ScoreDist\.png/;
 GetScoreDistribution($plot,$scoretype,$setype,@results);

 #Zoom score plot
 my $plot = $imagefile;
 $plot =~ s/\.png/ZoomScore\.png/;
 ZoomScoreDistribution($plot,$scoretype,$setype,@results);


  if($opt_n)
  {
  #nter plot
  my $plot = $imagefile;
  $plot =~ s/\.png/NterDist\.png/;
  GetNtermPlot($setype,$plot,@results);
  }


 }

 sub TandemConcat
 {
 my $tandemfile = shift;

 my $taxonomy= $opt_y;
 my $input;
 my($ptr1) = ParseTandem($tandemfile,$input,$taxonomy);
 my @results = @{$ptr1};

  if($opt_z)
  {
  SimpleCount('Tandem',@results);
  }
 my $scoretype = 'S';
 my $fdrtype = 5;
 my $setype = "T";
 
 my $imagefile = $opt_e;

 #Get the actual FDR values for both Gyig and Jones method
 my %all_results = GetFDRValues($opt_w,'0','T',$opt_n,$opt_t,@results);
 
 my %gygi_results;
 my %jones_results;

 #put the different results into the respective hashes
  foreach my $fdr (keys %all_results)
  {
  $gygi_results{$fdr}[0] = $all_results{$fdr}[0];
  $gygi_results{$fdr}[1] = $all_results{$fdr}[1];
  $jones_results{$fdr}[0] = $all_results{$fdr}[2];
  $jones_results{$fdr}[1] = $all_results{$fdr}[3];
  }



 #summary file
 my $plot = $imagefile;
 $plot =~ s/\.png/Summary\.txt/;
 AddToCreateSummaryTable('T',$plot,\%gygi_results,\%jones_results);
 AddToPeptideList('T',$plot,\%all_results);

 #deltamass plot
 my $plot = $imagefile;
 $plot =~ s/\.png/DeltaMass\.png/;
 GetDeltaMassPlot($plot,$scoretype,$setype,@results);

 #gygi rank plot
 my $plot = $imagefile;
 $plot =~ s/\.png/GygiRank\.png/;
 GetGygiRankPlot($setype,$plot,$opt_k,@results);

 #Score distribution
 my $plot = $imagefile;
 $plot =~ s/\.png/ScoreDist\.png/;
 GetScoreDistribution($plot,$scoretype,$setype,@results);

 #Zoom score plot
 my $plot = $imagefile;
 $plot =~ s/\.png/ZoomScore\.png/;
 ZoomScoreDistribution($plot,$scoretype,$setype,@results);


  if($opt_n)
  {
  #nter plot
  my $plot = $imagefile;
  $plot =~ s/\.png/NterDist\.png/;
  GetNtermPlot($setype,$plot,@results);
  }

 }



 sub MascotConcat
 {
 my $mascotfile = shift;
 my $mgf = $opt_m;
 my $fdrtype = 0;
  if($opt_q)
  {
  $fdrtype = $opt_q;
  }
  
 my $setype = "M";

 my $scoretype = "S";

 my $imagefile = $opt_o;
 #get the database
 my $db = GetDbPath($mascotfile);

 my $ptr = ParseMascot($mascotfile,$mgf,$db,'H','0.05');
 my @results = @{$ptr};

 if($opt_z)
 {
 SimpleCount('Mascot',@results);
 }

 #Get the actual FDR values for both Gyig and Jones method
 my %all_results = GetFDRValues($opt_w,'0','M',$opt_n,$opt_t,@results);

my %gygi_results;
my %jones_results;

 #put the different results into the respective hashes
 foreach my $fdr (keys %all_results)
 {
 $gygi_results{$fdr}[0] = $all_results{$fdr}[0];
 $gygi_results{$fdr}[1] = $all_results{$fdr}[1];
 $jones_results{$fdr}[0] = $all_results{$fdr}[2];
 $jones_results{$fdr}[1] = $all_results{$fdr}[3];
 }
    

 #summary file
 my $plot = $imagefile;
 $plot =~ s/\.png/Summary\.txt/;
 AddToCreateSummaryTable('M',$plot,\%gygi_results,\%jones_results);
 AddToPeptideList('M',$plot,\%all_results);

 #deltamass plot
 my $plot = $imagefile;
 $plot =~ s/\.png/DeltaMass\.png/;
 GetDeltaMassPlot($plot,$scoretype,$setype,@results);

 #gygi rank plot
 my $plot = $imagefile;
 $plot =~ s/\.png/GygiRank\.png/;
 GetGygiRankPlot($setype,$plot,$opt_k,@results);

 #Score distribution
 my $plot = $imagefile;
 $plot =~ s/\.png/ScoreDist\.png/;
 GetScoreDistribution($plot,$scoretype,$setype,@results);

 #Zoom score plot
 my $plot = $imagefile;
 $plot =~ s/\.png/ZoomScore\.png/;
 ZoomScoreDistribution($plot,$scoretype,$setype,@results);

  if($opt_n)
  {
  #nter plot
  my $plot = $imagefile;
  $plot =~ s/\.png/NterDist\.png/;
  GetNtermPlot($setype,$plot,@results);
  }

 }


sub CreateConcatChart
{

my $original_file = shift;
my $imagefile = shift;
my $scoretype = shift;
my $setype = shift;
my $fdrtype = shift;
my @results = @_;


my @image;


 #Get the actual FDR values for both Gyig and Jones method
 my %all_results = GetFDRValues($opt_w,'0','M',$opt_n,$opt_t,@results);
 
 my %gygi_results;
 my %jones_results;

 #put the different results into the respective hashes
  foreach my $fdr (keys %all_results)
  {
  $gygi_results{$fdr}[0] = $all_results{$fdr}[0];
  $gygi_results{$fdr}[1] = $all_results{$fdr}[1];
  $jones_results{$fdr}[0] = $all_results{$fdr}[2];
  $jones_results{$fdr}[1] = $all_results{$fdr}[3];
  }


#summary file
my $plot = $imagefile;
$plot =~ s/\.png/Summary\.txt/;
AddToCreateSummaryTable($setype,$plot,\%gygi_results,\%jones_results);
AddToPeptideList($setype,$plot,\%all_results);

#deltamass plot
my $plot = $imagefile;
$plot =~ s/\.png/DeltaMass\.png/;
#print "whats this $plot?\n";
GetDeltaMassPlot($plot,$scoretype,$setype,@results);

#gygi rank plot
my $plot = $imagefile;
$plot =~ s/\.png/GygiRank\.png/;
GetGygiRankPlot($setype,$plot,$opt_k,@results);

#Score distribution
my $plot = $imagefile;
$plot =~ s/\.png/ScoreDist\.png/;
GetScoreDistribution($plot,$scoretype,$setype,@results);

#Zoom score plot
my $plot = $imagefile;
$plot =~ s/\.png/ZoomScore\.png/;
ZoomScoreDistribution($plot,$scoretype,$setype,@results);



# push(@image,GetStats($scoretype,$fdrtype,$setype,$opt_n,@results));

# push(@image,GetDeltaMassPlot($scoretype,$setype,@results));
# push(@image,GetGygiRankPlot(@results));
# push(@image,GetFDR($scoretype,$fdrtype,$setype,@results));
# push(@image,GetScoreDistribution($scoretype,$setype,@results));
# push(@image,ZoomScoreDistribution($scoretype,$setype,@results));

  if($opt_n)
  {
#  push(@image,GetNtermPlot(@results));
my $plot = $imagefile;
$plot =~ s/\.png/NterDist\.png/;
GetNtermPlot($setype,$plot,@results);
#  GetFinalPNGNter($imagefile,@image);
  }
#  else
#  {
#  GetFinalPNG($imagefile,@image);
#  }

print "imagefile is $imagefile\n";

 return $imagefile;
}


sub StoreGygiResults
{
my $file = shift;
my $ptr = shift;
my %fdr = %{$ptr};
my $ptr2 = shift;
my @results = @{$ptr2};


 #for all the different fdrs
 foreach my $f (keys %fdr)
 {
  #for all the results
  for(my $s=0 ; $s<scalar(@results) ; $s++)
  {
  my $r = 1;
   if($results[$s][$r]{'sequence'} && $results[$s][$r]{'sequence'} ne "NULL")
   {
    #is the score bigger than the threshold for this fdr and forward?
    if($results[$s][$r]{'ionscore'} >= $fdr{$f} && $results[$s][$r]{'proteinlist'} !~ m/^REV/ && $results[$s][$r]{'proteinlist'} !~ m/\#\*\#REV/)
    {
    #is it nter?
     if($results[$s][$r]{'start'}<3 && !$gygi_results{$file}{$f}[0])
     { 
     $gygi_results{$file}{$f}[0] = 1;
     }
     elsif($results[$s][$r]{'start'}<3)
     {
     $gygi_results{$file}{$f}[0]++;
     }

     #store all results 
     if(!$gygi_results{$file}{$f}[1])
     {
     $gygi_results{$file}{$f}[1] = 1;          
     }
     else
     {
     $gygi_results{$file}{$f}[1]++;
     }

     #and store all the peptides
     if($results[$s][$r]{'start'}<3)
     {
     $gygi_peptides{$f}[0]{$results[$s][$r]{'sequence'}} = $results[$s][$r]{'ionscore'};
     }
    $gygi_peptides{$f}[1]{$results[$s][$r]{'sequence'}} = $results[$s][$r]{'ionscore'};

    #store the values for investigate peptides - make sure this protein/peptide not already stored for this file
    DealWithProteins($results[$s][$r]{'sequence'},$results[$s][$r]{'protein'},$file,$results[$s][$r]{'start'},$results[$s][$r]{'ionscore'},$f);

    }
   }
  }
 }

}

sub DealWithJonesProteins
{

my $sequence = shift;
my $protein = shift;
my $file = shift;
my $start = shift;
my $score = shift;
my $f = shift;

#store the values for investigate peptides - make sure this protein/peptide not already stored for this file
my $seen = 0;
my $original = $protein;
$protein =~ s/\_SIG$//;
$protein =~ s/\_\_/\_/g;
$protein =~ s/SR\_//g;
$protein =~ s/\_SS\d+//;
$protein =~ s/\_SR\_//;
if($protein =~ m/^REV\_/)
{
my @split = split/\_/,$protein;
$protein = $split[0] . "_" . $split[1];
}
else
{
my @split = split/\_/,$protein;
$protein = $split[0];
}

 if($investigate_proteins_jones{$protein}{$f}[0])
 {
  for(my $i=0 ; $i<scalar(@{$investigate_proteins_jones{$protein}{$f}}) ; $i++)
  {
   if($investigate_proteins_jones{$protein}{$f}[$i]{'sequence'} eq $sequence && $investigate_proteins_jones{$protein}{$f}[$i]{'file'} eq $file)
   {
   $seen = 1;
   }
  }
 }
 #if we haven't already stored this protein/peptide for this file
 if($seen == 0)
 {
 my $i = 0;
  if(%investigate_proteins_jones)
  {
  $i = scalar(@{$investigate_proteins_jones{$protein}{$f}});
  }


 $investigate_proteins_jones{$protein}{$f}[$i]{'sequence'} = $sequence;
 $investigate_proteins_jones{$protein}{$f}[$i]{'start'} = $start;
 $investigate_proteins_jones{$protein}{$f}[$i]{'ionscore'} = $score;
 $investigate_proteins_jones{$protein}{$f}[$i]{'file'} = $file;
 $investigate_proteins_jones{$protein}{$f}[$i]{'original'} = $original;
 }

}



sub DealWithProteins
{

my $sequence = shift;
my $protein = shift;
my $file = shift;
my $start = shift;
my $score = shift;
my $f = shift;

my $original = $protein;
$protein =~ s/\_SIG$//;
$protein =~ s/\_\_/\_/g;
$protein =~ s/SR\_//g;
$protein =~ s/\_SS\d+//;
$protein =~ s/\_SR\_//;
if($protein =~ m/^REV\_/)
{
my @split = split/\_/,$protein;
$protein = $split[0] . "_" . $split[1];
}
else
{
my @split = split/\_/,$protein;
$protein = $split[0];
}


my $seen = 0; 
 if($investigate_proteins{$protein}{$f}[0])
 {
 #store the values for investigate peptides - make sure this protein/peptide not already stored for this file
  for(my $i=0 ; $i<scalar(@{$investigate_proteins{$protein}{$f}}) ; $i++)
  {
   if($investigate_proteins{$protein}{$f}[$i]{'sequence'} eq $sequence && $investigate_proteins{$protein}{$f}[$i]{'file'} eq $file)
   {
   $seen = 1;
   }
  }
 }
 #if we haven't already stored this protein/peptide for this file
 if($seen == 0)
 {
 my $i =0 ;


  if(%investigate_proteins)
  {
  $i = scalar(@{$investigate_proteins{$protein}{$f}});
  }
 $investigate_proteins{$protein}{$f}[$i]{'sequence'} = $sequence;
 $investigate_proteins{$protein}{$f}[$i]{'start'} = $start;
 $investigate_proteins{$protein}{$f}[$i]{'ionscore'} = $score;
 $investigate_proteins{$protein}{$f}[$i]{'file'} = $file;
 $investigate_proteins{$protein}{$f}[$i]{'original'} = $original;
 }

}



sub StoreJonesResults
{
my $file = shift;
my $ptr = shift;
my %fdr = %{$ptr};
my $ptr2 = shift;
my @results = @{$ptr2};

 #for all the different fdrs
 foreach my $f (keys %fdr)
 {

  #for all the results
  for(my $s=0 ; $s<scalar(@results) ; $s++)
  {
  my $r = 1;
   if($results[$s][$r]{'sequence'} && $results[$s][$r]{'sequence'} ne "NULL")
   {
    #is the score bigger than the threshold for this fdr and forward?
    if($results[$s][$r]{'ionscore'} >= $fdr{$f} && $results[$s][$r]{'proteinlist'} !~ m/^REV/ && $results[$s][$r]{'proteinlist'} !~ m/\#\*\
#REV/)
    {

    #is it nter?
     if($results[$s][$r]{'start'}<3 && !$jones_results{$file}{$f}[0])
     {
     $jones_results{$file}{$f}[0] = 1;
     }
     elsif($results[$s][$r]{'start'}<3)
     {
     $jones_results{$file}{$f}[0]++;
     }

     #store all results
     if(!$jones_results{$file}{$f}[1])
     {
     $jones_results{$file}{$f}[1] = 1;
     }
     else
     {
     $jones_results{$file}{$f}[1]++;
     }

     #and store all the peptides
     if($results[$s][$r]{'start'}<3)
     {
     $jones_peptides{$f}[0]{$results[$s][$r]{'sequence'}} = $results[$s][$r]{'ionscore'};
     }
    $jones_peptides{$f}[1]{$results[$s][$r]{'sequence'}} = $results[$s][$r]{'ionscore'};

    #store the values for investigate peptides - make sure this protein/peptide not already stored for this file
    DealWithJonesProteins($results[$s][$r]{'sequence'},$results[$s][$r]{'protein'},$file,$results[$s][$r]{'start'},$results[$s][$r]{'ionscore'},$f);


    }
   }
  }
 }

}


sub SimpleCount
{
my $se = shift;
my @results = @_;

 if(!$opt_z)
 {
 return;
 }

if($se eq "Mascot")
{
$se = "mascot";
}
elsif($se eq "OMSSA")
{
$se = "omssa";
}
else
{
$se = "X!Tandem";
}
#print "se is now $se\n"; 


open(RES,">>$jones_combined_file") or die "unable to open the file $jones_combined_file to append to\n";

 for(my $s=0 ; $s<scalar(@results) ; $s++)
 {
  for(my $rank=1 ; $rank<11 ; $rank++)
  {
   if($results[$s][$rank]{'sequence'})
   {
    #if engine is omssa, score is \n
    my $score;
    if($se eq "omssa")
    {
    $score = '\\N';
    }
    else
    {
    $score = $results[$s][$rank]{'ionscore'};
    }
   #try and remove any fancy protein naming
   #signal peps
   $results[$s][$rank]{'protein'} =~ s/\_SIG$//;
   $results[$s][$rank]{'protein'} =~ s/\_\_/\_/g;
   $results[$s][$rank]{'protein'} =~ s/SR\_//g;
   $results[$s][$rank]{'protein'} =~ s/\_SS\d+//;
   $results[$s][$rank]{'protein'} =~ s/\_SR\_//;
    #what about the frayed proteins?
    if($results[$s][$rank]{'protein'} =~ m/^REV\_/)
    {
    my @split = split/\_/,$results[$s][$rank]{'protein'};
    $results[$s][$rank]{'protein'} = $split[0] . "_" . $split[1];
    } 
    else
    { 
    my @split = split/\_/,$results[$s][$rank]{'protein'};
    $results[$s][$rank]{'protein'} = $split[0];
    }


   my @split = split/\|/,$results[$s][$rank]{'protein'};
   $results[$s][$rank]{'protein'} = $split[0];
   print RES "$results[$s][$rank]{'protein'}\t$results[$s][$rank]{'sequence'}\t$score\t$results[$s][$rank]{'expect'}\t$rank\t$s\t$se\tNterminalaccession1\tNterminal study\t$results[$s][$rank]{'start'}\n";


   }
  }
 }
close RES;





}

 
sub GetDbPath
{
my $file = shift;
my $db;
open(FILE,"<$file") or die "unable to open the mascot results, $file\n";
 while(my $line = <FILE>)
 {
  if($line =~ m /^fastafile\=/)
  {
  my @split = split/\=/,$line;
  $db = $split[1];
  close FILE;
  }  
 }
close FILE;

$db =~ s/\/usr\/local\/mascot\//\/fs\/msct\//;


return $db;

}

sub GetMax
{
my $ptr = shift;
my $s = shift;
my @res = @{$ptr};

my $count = 1;
 for(my $r=1 ; $r<11 ; $r++)
 {
  if($res[$s][$r]{'sequence'})
  {
  $count++;
  }
 }

return $count;
}



sub ErrorMsg
{

my $response = shift;

my $msg;
 if($response == 1)
 {
 $msg = "You need to enter a result file";
 }  
 if($response == 3)
 {
 $msg = "You need to include the OMSSA database file\n";
 }
 if($response == 4)
 {
 $msg = "The list of omssa files should contain the omssa file,db like the following\n";
 $msg .= "my_omssa_file.csv,my_database.fasta\n";
 }


 if($msg)
 {
 print "\n\tIncorrect Usage: \n";
 print "\n\t$msg\n\n";
 }
 else
 {
 print "\n\tUsage:\n";
 }

print "\tInput files\n";
print "\t\t-m\tMGF file\n";
print "\t\t-c\tConcatanated Mascot file\n";
print "\t\t-b\tFrayed database search\n";
print "\t\t-f\tMascot .dat file (forward only)\n";
print "\t\t-d\tMascot .dat file (Decoy only)\n";
print "\t\t-l\tname of file containing a list of mascot files - analyse over all\n";
print "\n\t\t-x\tConcatanated Xtandem file\n";
print "\t\t-u\tXTandem forward only file\n";
print "\t\t-v\tXTandem reverse only file\n";
print "\t\t-y\tTaxonomy file associated with the X!Tandem search\n";
print "\t\t-e\tTandem output image name\n";
print "\n\t\t-r\tConcatanated Omssa result file\n";
print "\t\t-j\tOmssa image name\n";
print "\t\t-p\tOmssa forward only file\n";
print "\t\t-i\tOmssa reverse only file\n";
print "\t\t-a\tname of db single omssa searched\n";
print "\n\tAnalysis Settings\n";
print "\t\t-z\tCombine SearchEngine Results using Kall method\n";
print "\t\t-t\tDecoy tag (default REV_)\n";
print "\t\t-q\tFDR type (0 (default) =  Gygi, 1 = Kall)\n";
print "\t\t-s\tSession id\n";
print "\t\t-n\tN-terminal analysis (yes = 1, no = 0 (default)\n";
print "\t\t-g\tProteinAnalysis\n";
print "\n\tOutputs\n";
print "\t\t-o\tMascot Output Image file name (default is png of dat filename)\n";
print "\t\t-h\tThis help page\n";


print "\n";

}







