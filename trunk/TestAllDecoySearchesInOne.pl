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
#use lib qw(/nfs/san_home/mbdssdw6/LIVERPOOL/bin/perl_modules/);
use lib qw(/var/www/localhost/cgi-bin/FDRAnalysis/perl_modules/);
#use lib '../bin/perl_modules';
use RunMascot;
use FDRAnalysisMascotParser;
use FDRAnalysisXTandem;
use FDRAnalysisOmssa;
use DecoyMascotParser;
use DeltaMassPlot;
use GygiRankPlot;
##use GetFileStats;
##use GetFDRPlot;
use GetFDRValues;
use GetDiscriminantScorePlot;
use ZoomDiscriminantScorePlot;
use GetNterminaldistribution;
use CreatePNG;
use Getopt::Std;
use GetFileStatsSeperateFiles;
use List::Util qw[min max];


########################################################################
####  An uber Mascot Parse to parse all types of decoy mascot searches ##
########################################################################

#our($opt_I,$opt_w,$opt_k,$opt_c, $opt_m, $opt_f, $opt_d, $opt_r, $opt_q, $opt_s, $opt_o, $opt_n, $opt_h, $opt_l, $opt_t, $opt_p, $opt_b, $opt_x, $opt_e, $opt_j, $opt_y, $opt_a, $opt_z, $opt_g, $opt_u, $opt_v, $opt_i);
#getopts('c:m:f:d:r:q:s:o:n:ht:l:p:b:x:e:j:y:a:zg:u:v:i:kw:I');
#DCW - opt_D added for decoy ratio
our($opt_D,$opt_I,$opt_w,$opt_k,$opt_c, $opt_m, $opt_f, $opt_d, $opt_r, $opt_q, $opt_s, $opt_o, $opt_n, $opt_h, $opt_l, $opt_t, $opt_p, $opt_b, $opt_x, $opt_e, $opt_j, $opt_y, $opt_a, $opt_z, $opt_g, $opt_u, $opt_v, $opt_i);
getopts('c:m:f:d:r:q:s:o:n:ht:l:p:b:x:e:j:y:a:z:g:u:v:i:k:w:ID:');

	print "TADSIO opt_h is ".$opt_h."\n";
	print "TADSIO opt_s is ".$opt_s."\n";
	print "TADSIO opt_c is ".$opt_c."\n";

   if($opt_h)
   {
   ErrorMsg('2');
   exit(1);
   }
   elsif(!$opt_c && !$opt_f && !$opt_d && !$opt_x && !$opt_r && !$opt_l && !$opt_q && !$opt_g)
   {
   ErrorMsg('1');
   exit(1);
   }

 if(!$opt_s)
 {
 print "You have not entered a file name for the output files\n";
 exit(1);
 }


 #if a decoy method what is reverse tag?
 if(!$opt_t)
 {
 $opt_t = "REV_";
 }

print "rev_tag=".$opt_t."\n";

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

print "OK1\n";

 # Mascot Concatanated?              
 if($opt_c)
 {
 MascotConcat($opt_c,'1');
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
 my @files = split/\,/,$opt_l;
 my $spectra_count = 1;
  for(my $f=0 ; $f<scalar(@files) ; $f++)
  {
  $files[$f] =~ s/\s+//g;
  my $db = GetDbPath($files[$f]);
  $spectra_count = MascotConcat($files[$f],$spectra_count);
  }				        
 }

 if($opt_q)
 {
 my @files = split/\,/,$opt_q;
 my $spectra_count = 1;
  for(my $f=0 ; $f<scalar(@files) ; $f++)
  {
  $files[$f] =~ s/\s+//g;
  $spectra_count = TandemConcat($files[$f],$spectra_count);
  }

 }

print "OK2";

 #X!Tandem
 if($opt_x)
 {
 TandemConcat($opt_x,'1');
 }
 if($opt_u && $opt_v)
 {
 AnalyseSeperate('T');
 } 

 if($opt_g)
 {
 my @files = split/\,/,$opt_g;
 my $spectra_count = 1;
  for(my $f=0 ; $f<scalar(@files) ; $f++)
  {
  $files[$f] =~ s/\s+//g;
  $spectra_count = OmssaConcat($files[$f],$spectra_count);
  }

 }

print "OK3";

 #Omssa
 if($opt_r)
 {
OmssaConcat($opt_r,'1')
 }
 if($opt_p && $opt_i)
 {
 AnalyseSeperate('O');
 }

print "OK4";

RunAndysStuff();

print "OK5";


sub AddToPeptideList
{
my $se = shift;
my $file = shift;
my $ptr = shift;
my %res = %{$ptr};

print "peptidelist going into " . $opt_s . "_peptidelist.txt\n";
$peptidelist{'file'}{'file'}{'file'}[0] = $opt_s . "_peptidelist.txt";
 foreach my $fdr (keys %res)
 {
 push(@{$peptidelist{$se}{'E'}{$fdr}},$res{$fdr}[4]);
 push(@{$peptidelist{$se}{'K'}{$fdr}},$res{$fdr}[5]);
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

$summary_table{'file'}{'file'}{'0'}[0][0] = $opt_s . "_summary.txt";


 #for all the gygi results 
 foreach my $fdr (sort keys %res)
 {
  #print("res fdr: $fdr\n");
  if($fdr<=0.5)
  {
  push(@{$summary_table{$se}{'E'}{$fdr}[0]},$res{$fdr}[0]);

   if($opt_n)
   {
   push(@{$summary_table{$se}{'E'}{$fdr}[1]},$res{$fdr}[1]);
   }
  }
 }
 foreach my $fdr (sort keys %jones)
 {
  #print("jones fdr: $fdr\n");
  if($fdr<=0.5)
  {
  push(@{$summary_table{$se}{'K'}{$fdr}[0]},$jones{$fdr}[0]);
   if($opt_n)
   {
   push(@{$summary_table{$se}{'K'}{$fdr}[1]},$jones{$fdr}[1]);
   }
  }
 }


}


sub CreatePeptideList
{

my %seq;
my %prot;
my %total_seq;
my %total_prot;
my %nter;
my %total_nter;


#open(TMP,">THISISTMP2.txt") or die "unable to open THISISTMP2.txt\n";

#open the peptide file
open(NEW,">$peptidelist{'file'}{'file'}{'file'}[0]") or die "unable to open the file $peptidelist{'file'}{'file'}{'file'}[0] to write to\n";
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
   #for all the files
    for(my $file=0 ; $file<scalar(@{$peptidelist{$se}{'K'}{$fdr}}) ; $file++)
    {
    my @split = split/\*\*\*/,$peptidelist{$se}{'K'}{$fdr}[$file];
     #for all the peptides
     for(my $i=0 ; $i<scalar(@split) ; $i++)
     {
     my @list = split/\#\#\#/,$split[$i];
     $seq{$engine}{$list[1]} = 1;
     $prot{$engine}{$list[0]} = 1;
     $total_seq{$list[1]} = 1;
     $total_prot{$list[0]} =1;
#print TMP "$list[1]\t$list[0]\n";
     print NEW "$engine\t$list[0]\t$list[1]\t$list[2]\t$list[3]\t$list[4]\n";
      if($list[2]<3)
      {
      $total_nter{$list[1]} = 1;
      $nter{$engine}{$list[1]} = 1;
      }
     }
    } 
   }  
  }
 }

close NEW;

#close TMP;
my $count = 0;

 foreach my $s (keys %seq)
 {
  foreach my $a (keys %{$seq{$s}})
  {
  $count++;
  }
 $count = 0;
 }
 
 if($opt_n)
 {
 my $count = 0;

  foreach my $s (keys %nter)
  {
   foreach my $a (keys %{$nter{$s}})
   {
   $count++;
   }
  $count = 0;
  }
 }

$count = 0; 
 foreach my $s (keys %prot)
 {
  foreach my $a (keys %{$prot{$s}})
  {
  $count++;
  }
 $count = 0;
 }

$count = 0;
 foreach my $s (keys %total_seq)
 {
 $count++;
 }

 if($opt_n)
 {
 $count = 0;
  foreach my $s (keys %total_nter)
  {
  $count++;
  }
 }




$count = 0;
 foreach my $s (keys %total_prot)
 {
 $count++;
 }




}

sub CreateSummary
{


#results for each file

print "Search Engine\tFDR method\tFDR value\tTotal unique peptide count";
 if($opt_n)
 {
 print "\tunique Nter count";
 }
print "\n\n";

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
    print " \t\t \t\t \t\t ";
      if($opt_n)
      {
      print "\t\t \t\t";
      }
    print "\n";
    }

    #print("writing summary\n");
    #for the different FDR values
    foreach my $fdr (sort keys %{$summary_table{$se}{$fdrtype}})
    {
     if($fdr<=$opt_w)
     {
      if($summary_table{$se}{$fdrtype}{$fdr}[0])
      {
       #print("***$fdr!!!\n");

       for(my $file=0 ; $file<scalar(@{$summary_table{$se}{$fdrtype}{$fdr}[0]}) ; $file++)
       {
        print "(file $file)$engine\t\t$fdrtype\t\t$fdr\t\t$summary_table{$se}{$fdrtype}{$fdr}[0][$file]";
        if($opt_n && $summary_table{$se}{$fdrtype}{$fdr}[1][$file])
	 {
	  print "\t\t$summary_table{$se}{$fdrtype}{$fdr}[1][$file]";
	 }
	 elsif($opt_n)
	 {
	  print "\t\t0";
	 }
        print "\n";	
       }
      }
      else
      {
       #print("else branch ***$fdr!!!\n");
       if($opt_n && $summary_table{$se}{$fdrtype}{$fdr}[1])
       {
        for(my $file=0 ; $file<scalar(@{$summary_table{$se}{$fdrtype}{$fdr}[0]}) ; $file++)
        {
	  print "(file $file) $engine\t\t$fdrtype\t\t$fdr\t\t0\t\t$summary_table{$se}{$fdrtype}{$fdr}[1][$file]\n";
	 }
       }
       else
       {
        print "$engine\t\t$fdrtype\t\t$fdr\t\t0";
        if($opt_n)
	 {
	  print "\t\t0";
	 }
        print "\n";
       }
      }
     }
    } 
   }
  } 
 }

#DCW - WRITE TO HTML (from TestMascotAllDecoySearchesInOne)
open(NEW,">$summary_table{'file'}{'file'}{'0'}[0][0]") or die "**unable to open the file $summary_table{'file'}{'file'}{'0'}[0] to write to\n";
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
      print NEW "<TR><TD $bgcolor>$engine</TD><TD $bgcolor>$fdrtype</TD><TD $bgcolor>$fdr</TD><TD $bgcolor>$summary_table{$se}{$fdrtype}{$fdr}[0][0]</TD>";
      }
      else
      {
      print NEW "<TR><TD $bgcolor>$engine</TD><TD $bgcolor>$fdrtype</TD><TD $bgcolor>$fdr</TD><TD $bgcolor>0</TD>";
      }
      if($opt_n)
      {
       if($summary_table{$se}{$fdrtype}{$fdr}[1])
       {
       print NEW "<TD $bgcolor>$summary_table{$se}{$fdrtype}{$fdr}[1][0]</TD>";
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
 #my $cmd = "rm " . $opt_s . "combined*";
 #system($cmd);


 #run andys!
 #
 print "opt_z is $opt_z and jones_combined file is *  $jones_combined_file *\n";
 if($opt_z && $jones_combined_file)
 {
 #remove the files
 #my $cmd = "rm " . $opt_s . "*combined*out";
 #my $cmd = "rm " . $opt_s . "*combined*"; #DCW
 #system($cmd);

 #my $cmd = "perl LIVERPOOL/bin/perl_scripts/Andys/andys_scripts/FDRScore_March09/MultipleSearch/web_FDRScore.pl ";
my $cmd = "perl /nfs/san_home/mbdssdw6/LIVERPOOL/bin/perl_scripts/Andys/andys_scripts/FDRScore_March09/MultipleSearch/web_FDRScore.pl ";
 $cmd .= "-D $opt_D -F $jones_combined_file -R 1 -I $opt_t -O ". $opt_s . "combined_results.out -P " . $opt_s . "combined_proteins.out -Q " . $opt_s . "combined_peptides.out -T " . $opt_w; 
print "andys command is $cmd\n";
system($cmd);

 ParseCombined();
 CreateSummary();
 CreatePeptideList();

print("copying peptide file");

 #so the combined peptides file is downloadable move it to the other tmp dir
 my $cmd = "cp " . $opt_s . "combined_peptides.out /var/www/localhost/htdocs/FDRAnalysis/tmp/";
 system($cmd);

  #print("before images\n");

  my $peplist = $opt_s . "_peptidelist.txt";
  if($opt_I)
  {


 my $comblist = $opt_s . "combined_peptides.out";
 my $imagefile = $opt_s . "_CombinedVennDiagram.png";
 $imagefile =~ s/\/var\/www\/tmp\//\/var\/www\/localhost\/htdocs\/FDRAnalysis\/tmp\//;
 my @values = GetVennNumbersFromCombine($comblist);

   if(!GetPNGVenn($imagefile,@values))
   {
   print "Error trying to run the GetPNGVenn\n";
   exit(1);
   }
  }
  #print("after images\n");
 #CreatePeptideList();#commented out DCW

 }
 else
 {
 CreateSummary();
 CreatePeptideList();
 } 

print("before images2\n");

 if($opt_I)
 {
 #also want to create the venn diagram for the peptide
 my $peplist = $opt_s . "_peptidelist.txt";
 my @values = GetVennNumbers($peplist);

print("copying peptide list");
#DCW - copy peptidelist
 $cmd = "cp " . $opt_s . "_peptidelist.txt /var/www/localhost/htdocs/FDRAnalysis/tmp/";
 system($cmd);

 my $imagefile = $opt_s . "_VennDiagram.png";
 $imagefile =~ s/\/var\/www\/tmp\//\/var\/www\/localhost\/htdocs\/FDRAnalysis\/tmp\//;

  if(!GetPNGVenn($imagefile,@values))
  {
  print "Error trying to run the GetPNGVenn\n";
  exit(1);
  }
 }

print("after images2\n");

 #if any duplicates are involved
 if($opt_g || $opt_q || $opt_l)
 {
 my $peplist = $opt_s . "_peptidelist.txt";
 my %seqs;
 my %nter;
 open(FILE,"<$peplist") or die "*unable to open the peplist file $peplist*\n";
 while(my $line = <FILE>)
 {
  if($line !~ m/^Search\s+Engine/)
  {
  my @split = split/\t/,$line;
  $seqs{$split[0]}{$split[2]} = 1;
   if($split[3]<3)
   {
   $nter{$split[0]}{$split[2]} = 1;
   }
  }
 }
 close FILE;

print("before overlap\n");

#Look at overlap
my %seq_overall;
my %nter_overall;

 foreach my $engine (keys %seqs)
 {
 my $count = 0;
 my $nter_count = 0;

  #print("OVERLAP optk $opt_k\n");

  #if the full overlap required
  if($opt_k)
  {
  GetFullOverlap($engine,\%{$seqs{$engine}},\%{$nter{$engine}});
  }

  foreach my $sequence (keys %{$seqs{$engine}})
  {
  $seq_overall{$sequence} = 1;
  $count++;
  }

  foreach my $sequence (keys %{$nter{$engine}})
  {
  $nter_overall{$sequence} = 1;
  $nter_count++;
  }

 print "$engine $count unique sequences (of which $nter_count are nter)\n";
 }

my $count = 0 ;
my $nter_count = 0;
 foreach my $sequence (keys %seq_overall)
 {
 $count++;
 }
 foreach my $sequence (keys %nter_overall)
 {
 $nter_count++;
 }
print "Over all search engines $count unique sequences, $nter_count of which are Nterminal\n"; 


 }

}


sub GetFullOverlap
{
print("Running GetFullOverlap\n");
my $engine = shift;
my $ptr = shift;
my $ptr2 = shift;

my %seq = %{$ptr};
my %nter = %{$ptr2};


my @files;
my @sequences;
 #the files are
 if($engine =~ m/^M/)
 {
 @files = split/\,/,$opt_l;
  for(my $f=0 ; $f<scalar(@files) ; $f++)
  {
  my $ptr = ParseMascot($files[$f],'','','H','0.05');
  my @results = @{$ptr};
  my %all_results = GetFDRValues('S',$opt_w,'0','M',$opt_n,$opt_t,@results);
   foreach my $fdr (keys %all_results)
   {
   my @split1 = split/\*\*\*/,$all_results{$fdr}[5];
    #for each protein
    for(my $prot=0 ; $prot<scalar(@split1) ; $prot++)
    {
    my @split2 = split/\#\#\#/,$split1[$prot];
    $sequences[$f]{$split2[1]} = 1;
    }
   }
  } 
 }
 elsif($engine =~ m/^O/)
 {
 @files = split/\,/,$opt_g;

  for(my $f=0 ; $f<scalar(@files) ; $f++)
  {
  my $ptr = ParseOmssa($files[$f],$opt_a);
  my @results = @{$ptr};
  my %all_results = GetFDRValues('S',$opt_w,'0','O',$opt_n,$opt_t,@results);
   foreach my $fdr (keys %all_results)
   {
   my @split1 = split/\*\*\*/,$all_results{$fdr}[5];
    #for each protein
    for(my $prot=0 ; $prot<scalar(@split1) ; $prot++)
    {
    my @split2 = split/\#\#\#/,$split1[$prot];
    $sequences[$f]{$split2[1]} = 1;
    }
   }
  }
 }
 elsif($engine =~ m/^[TX]/)
 {

 @files = split/\,/,$opt_q;
  for(my $f=0 ; $f<scalar(@files) ; $f++)
  {
  my $taxonomy= $opt_y;
  my $input;
print "ParseTandem($files[$f],$input,$taxonomy)\n";
  my($ptr1) = ParseTandem($files[$f],$input,$taxonomy);
  my @results = @{$ptr1};

  my %all_results = GetFDRValues('S',$opt_w,'0','T',$opt_n,$opt_t,@results);#code from this file
  #my %all_results = GetFDRValues('C',$opt_w,'0','T',$opt_n,$opt_t,@results);#code from this file
  #my %all_results = GetFDRValues('S',$opt_w,'0','O',$opt_n,$opt_t,@results);#JENNY'S ORIGINAL CODE
   foreach my $fdr (keys %all_results)
   {
   my @split1 = split/\*\*\*/,$all_results{$fdr}[5];
    #for each protein
    for(my $prot=0 ; $prot<scalar(@split1) ; $prot++)
    {
    my @split2 = split/\#\#\#/,$split1[$prot];
    $sequences[$f]{$split2[1]} = 1;
    }
   }
  } 
 }

#obs and nter_obs look like this: "0.", "0.1.", "0.1.2.", etc.
my %obs;
my %nter_obs;
 #for the sequences I am interested in
 foreach my $pep (keys %seq)
 {

  #for all the resplicates
  for(my $f=0 ; $f<scalar(@sequences) ; $f++)
  {
   if($sequences[$f]{$pep})
   {
   $obs{$pep} .= "$f.";
   }
  }
 }
 foreach my $pep (keys %nter)
 {

  for(my $f=0 ; $f<scalar(@sequences) ; $f++)
  {
   if($sequences[$f]{$pep})
   {
   $nter_obs{$pep} .= "$f.";
   }
  }
 }
my %counter;
my %nter_counter;
 foreach my $pep (keys %obs)
 {
  if($counter{$obs{$pep}})
  {
  $counter{$obs{$pep}}++;
  }
  else
  {
  $counter{$obs{$pep}} = 1;
  }
 }
 foreach my $pep (keys %nter_obs)
 {

  if($nter_counter{$nter_obs{$pep}})
  {
  $nter_counter{$nter_obs{$pep}}++;
  }
  else
  {
  $nter_counter{$nter_obs{$pep}} = 1;
  }
 }

print "For search engine $engine (based on Kall FDR method) the overlap between files are as follows\n";
print "Combination (i.e files seen in separated by '.'\tpeps seen\tof which nterminal\n";

foreach my $combination (keys %counter)
{
print "$combination\t$counter{$combination}\t$nter_counter{$combination}\n";
}
foreach my $combination (keys %nter_counter)
{
 if(!$counter{$combination})
 {
 print "$combination\t0\t$nter_counter{$combination}\n";
 }
}


}
sub GetVennNumbersFromCombine
{
my $combine = shift;

my @res;
my %peps;
my @seen;

open(FILE,"<$combine") or die "unable to open the results file $combine\n";
 while(my $line = <FILE>)
 {
 my @split = split/\t/,$line;
 my $seq = $split[1];
 my $m = $split[4];
 my $o = $split[5];
 my $t = $split[6];

  if($m)
  {
  $seen[1] = 1;
  $peps{$seq}[1] = 1;
  }
  if($o)
  {
  $seen[2] = 1;
  $peps{$seq}[2] = 1;
  }
  if($t)
  {
  $seen[3] = 1;
  $peps{$seq}[3] = 1;
  }
 }
close FILE;

 if($seen[1])
 {
 $res[0][1] = "Mascot";
 }
 if($seen[2])
 {
 $res[0][2] = "Omssa";
 }
 if($seen[3])
 {
 $res[0][3] = "X!Tandem";
 }


  #for all the peptides in the file
  foreach my $seq (keys %peps)
  {
  my $m = $peps{$seq}[1];
  my $o = $peps{$seq}[2];
  my $t = $peps{$seq}[3];

   if($m && $o && $t)
   {
    if($res[0][0])
    {
    $res[0][0]++;
    }
    else
    {
    $res[0][0]=1;
    }
   }
   elsif($m && $o)
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
   elsif($m && $t)
   {
    if($res[1][3])
    {
    $res[1][3]++;
    }
    else
    {
    $res[1][3]=1;
    }
   }
   elsif($o && $t)
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
   elsif($m)
   {
    if($res[1][0])
    {
    $res[1][0]++;
    }
    else
    {
    $res[1][0]=1;
    }
   }
   elsif($o)
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
   elsif($t)
   {
    if($res[3][0])
    {
    $res[3][0]++;
    }
    else
    {
    $res[3][0]=1;
    }
   } 
  }

return @res;

}


sub GetVennNumbers
{


#peptidelist
my $newpeplist = shift;

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
 print("ParseCombined running\n");
 my %res;
 my $file = $opt_s . "combined_peptides.out";

#open the peptide results
 open(FILE2,"<$file") or print "problem with Combined results ($file)\n";
  while(my $line = <FILE2>)
  {
  my @split = split/\t/,$line;
  $res{$split[0]}{$split[1]} = $split[3];
  }
 close FILE2;

my %observed;
my $count = 0;
my $nter = 0;
#for the different FDRs
 foreach my $fdr (keys %res)
 {
  print("COMBINED $fdr,$opt_w\n");

  if($fdr>=$opt_w)
  {

   foreach my $seq (keys %{$res{$fdr}})
   {
    if(!$observed{$seq})
    {
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
 $summary_table{'C'}{'K'}{$opt_w}[0][0] = $count;
 $summary_table{'C'}{'K'}{$opt_w}[1][0] = $nter;

 }

 sub MascotDecoy
 {

 my $mascotfile = $opt_d;
 my $db = GetDbPath($mascotfile);
 my $mgf = "";
 my $ptr = ParseDecoyMascot($mascotfile,$mgf,$db,'H','0.05',$opt_t);
 my @resultsTwo = @{$ptr};
 my $ptr = ParseMascot($mascotfile,$mgf,$db,'H','0.05');
 my @results = @{$ptr};


 my %peptide_hits;

 my $fdrtype = 0;
 my $setype = "M";
 my $max  = max(scalar(@results),scalar(@resultsTwo));

 my @new_results = ProduceUber($max,\@results,\@resultsTwo);
 my @reranked_results = ProduceUberReRanked($max,\@results,\@resultsTwo);


 my @image;
my $fdrtype = 0;

  #DCW - from TestMascotAllDecoySearchesInOne
  if($opt_q)
  {
  $fdrtype = $opt_q;
  }

 my $scoretype = "S";

 #my $imagefile = "MascotDecoy.png";
 my $imagefile = $opt_f; #DCW - from TestMascotAllDecoySearchesInOne

 $imagefile =~ s/\.dat//;
 $imagefile = $imagefile . "_" . $scoretype . ".png";

  #DCW - from TestMascotAllDecoySearchesInOne
  if($opt_o)
  {
  $imagefile = $opt_o;
  }

 #print("MascotDecoy calling with search type S\n");

 #Get the actual FDR values for both Gygi and Jones method
 my %all_results = GetFDRValues('S',$opt_w,'0','M',$opt_n,$opt_t,@new_results);
 
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

  if($opt_I)
  {
  #deltamass plot
  my $plot = $imagefile;
  $plot =~ s/\.png/DeltaMass\.png/;
  GetDeltaMassPlot($plot,$scoretype,$setype,$opt_t,@new_results);
 
  #gygi rank plot
  my $plot = $imagefile;
  $plot =~ s/\.png/GygiRank\.png/;

  GetGygiRankPlot($opt_t,$setype,$plot,$all_results{$opt_w}[7],@reranked_results);#DCW - use FDR threshold, NOT expect value, must add a small amount
  # GetGygiRankPlot($opt_t,$setype,$plot,'0.05',@new_results);
  #GetGygiRankPlot($opt_t,$setype,$plot,'0.05',@reranked_results);
  #Score distribution
  my $plot = $imagefile;
  $plot =~ s/\.png/ScoreDist\.png/;
  GetScoreDistribution($plot,$scoretype,$setype,$opt_t,@new_results);

  #Zoom score plot
  my $plot = $imagefile;
  $plot =~ s/\.png/ZoomScore\.png/;
  ZoomScoreDistribution($plot,$scoretype,$setype,@new_results);


   if($opt_n)
   {
   #nter plot
   my $plot = $imagefile;
   $plot =~ s/\.png/NterDist\.png/;
   GetNtermPlot($setype,$plot,$opt_t,@new_results);
   }
  }

  #DCW 290709
  if($opt_z)
  {
   SimpleCount('Mascot','1',@new_results);
  }
 }


 sub AnalyseSeperate
 {
 my $setype = shift;
 my $inputFile = 'LIVERPPOL/test/input.xml';
 my $taxonomyFile = 'LIVERPPOL/test/taxonomy.xml';
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
 
  #$imagefile = "MascotSeperate.png";
  $imagefile = $opt_o; #DCW -  from TestMascotAllDecoySearchesInOne
  }

  if($setype eq "T")
  {
  my $forward = $opt_u;
  my $reverse = $opt_v;
  my $ptr = ParseTandem($forward,$inputFile,$taxonomyFile);
  @results = @{$ptr};
 
  my $ptr2 = ParseTandem($reverse,$inputFile,$taxonomyFile);
  @resultsTwo = @{$ptr2};

  #$imagefile = "TandemSeperate.png";
  $imagefile = $opt_e; #DCW -  from TestMascotAllDecoySearchesInOne

  }
 
  if($setype eq "O")
  {
  my $forward = $opt_p;
  my $reverse = $opt_i;

  my $ptr = ParseOmssa($forward,"");
  @results = @{$ptr};

  my $ptr2 = ParseOmssa($reverse,"");
  @resultsTwo = @{$ptr2};
 
  #$imagefile = "OmssaSeperate.png";
  $imagefile = $opt_j; #DCW - from TestMascotAllDecoySearchesInOne
  }

 my %peptide_hits;

 my $fdrtype = 0;
 my $max  = max(scalar(@results),scalar(@resultsTwo));
 my @new_results = ProduceUber($max,\@results,\@resultsTwo);


  #if combining the result
  if($opt_z && $setype eq "O")
  {
  SimpleCount('OMSSA','1',@new_results);
  }
  elsif($opt_z && $setype eq "M")
  {
  SimpleCount('Mascot','1',@new_results);
  } 
  elsif($opt_z && $setype eq "T")
  {
  SimpleCount('Tandem','1',@new_results);
  }


 my @image;
my $fdrtype = 0;


 my $scoretype = "S";

 #Get the actual FDR values for both Gyig and Jones method
 my %all_results = GetFDRValues('S',$opt_w,'0','M',$opt_n,$opt_t,@new_results);
 
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

 if($opt_I)
 {
 #deltamass plot
 my $plot = $imagefile;
 $plot =~ s/\.png/DeltaMass\.png/;
 GetDeltaMassPlot($plot,$scoretype,$setype,$opt_t,@new_results);

 # gygi rank plot
 my $plot = $imagefile;
 $plot =~ s/\.png/GygiRank\.png/;
 #GetGygiRankPlot($opt_t,$setype,$plot,'0.05',@new_results);
 GetGygiRankPlot($opt_t,$setype,$plot,$all_results{$opt_w}[7],@new_results);#DCW - use FDR threshold, NOT expect value, must add a small amount

 #Score distribution
 my $plot = $imagefile;
 $plot =~ s/\.png/ScoreDist\.png/;
 GetScoreDistribution($plot,$scoretype,$setype,$opt_t,@new_results);

 #Zoom score plot
 my $plot = $imagefile;
 $plot =~ s/\.png/ZoomScore\.png/;
 ZoomScoreDistribution($plot,$scoretype,$setype,@new_results);

  if($opt_n)
  {
  #nter plot
  my $plot = $imagefile;
  $plot =~ s/\.png/NterDist\.png/;
  GetNtermPlot($setype,$plot,$opt_t,@new_results);
  }
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
my $size = scalar(@results);
 for(my $r=0 ; $r<$max ; $r++)
 {
 #find the maximum hit
  for(my $k=1 ; $k<11 ; $k++)
  {
   if($resultsTwo[$r][$k]{'sequence'})
   {
   $size++;
   foreach my $attr (keys %{$resultsTwo[$r][$k]})
    {
     if($attr eq "protein")
     {
      if($resultsTwo[$r][$k]{$attr} !~ m/^$opt_t/)
      {
      $resultsTwo[$r][$k]{$attr} = $opt_t . $resultsTwo[$r][$k]{$attr};
      }
     }
    $new_results[$size][$k]{$attr} = $resultsTwo[$r][$k]{$attr};
    }
   }
  
   if($results[$r][$k]{'sequence'})
   {
    foreach my $attr (keys %{$results[$r][$k]})
    {
    $new_results[$r][$k]{$attr} = $results[$r][$k]{$attr};
    }
   }
  } 
 }
return @new_results;
}

 sub ProduceUberReRanked
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
    #no target hit
    if(!$new_results[$r][$rank+1]{'ionscore'})
    {
    $rank++;
     if($rank<11)
     {
      foreach my $attr (keys %{$resultsTwo[$r][$k]})
      {
       #make sure there is a tag on the prot - if not add one! 
       if($attr eq "protein")
       {
        if($resultsTwo[$r][$k]{$attr} !~ m/^$opt_t/)
        {
        $resultsTwo[$r][$k]{$attr} = $opt_t . $resultsTwo[$r][$k]{$attr};
        }
       }	
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

       #make sure there is a tag on the prot - if not add one!
       if($attr eq "protein")
       {
        if($resultsTwo[$r][$k]{$attr} !~ m/^$opt_t/)
        {
        $resultsTwo[$r][$k]{$attr} = $opt_t . $resultsTwo[$r][$k]{$attr};
        }
       }

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
        #make sure there is a tag on the prot - if not add one!
        if($attr eq "protein")
        {
         if($resultsTwo[$r][$k]{$attr} !~ m/^$opt_t/)
         {
         $resultsTwo[$r][$k]{$attr} = $opt_t . $resultsTwo[$r][$k]{$attr};
         }
        }
       $new_results[$r][$position]{$attr} = $resultsTwo[$r][$k]{$attr};
       }
      $position = 11;
      }
      else
      {
       foreach my $attr (keys %{$results[$r][$k]})
       {

       #make sure there is a tag on the prot - if not add one!
       if($attr eq "protein")
       { 
        if($resultsTwo[$r][$k]{$attr} !~ m/^$opt_t/)
        {
        $resultsTwo[$r][$k]{$attr} = $opt_t . $resultsTwo[$r][$k]{$attr};
        }
       }
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
 my $spectra_count = shift;

 my($ptr1) = ParseOmssa($omssafile,$opt_a);
 my @results = @{$ptr1};
  #if combining the result
  if($opt_z)
  { 
  SimpleCount('OMSSA',$spectra_count,@results);
  }

 my $scoretype = 'S';
 my $fdrtype = 5;
 my $setype = "O";

 #my $imagefile = "OmssaConcat" . $spectra_count . ".png";
 my $imagefile = $opt_j; #DCW - from TestMascotAllDecoySearchesInOne

 #Get the actual FDR values for both Gyig and Jones method
 my %all_results = GetFDRValues('C',$opt_w,'0','O',$opt_n,$opt_t,@results);
 
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


  if($opt_I)
  {
  #deltamass plot
  my $plot = $imagefile;
  $plot =~ s/\.png/DeltaMass\.png/;

print "In omssa GetDeltaMassPlot($plot,$scoretype,$setype,$opt_t,\n";
  GetDeltaMassPlot($plot,$scoretype,$setype,$opt_t,@results);

  #gygi rank plot
  my $plot = $imagefile;
  $plot =~ s/\.png/GygiRank\.png/;
  #GetGygiRankPlot($opt_t,$setype,$plot,'0.05',@results);
  GetGygiRankPlot($opt_t,$setype,$plot,$all_results{$opt_w}[7],@results);#DCW - use FDR threshold, NOT expect value, must add a small amount

  #Score distribution
  my $plot = $imagefile;
  $plot =~ s/\.png/ScoreDist\.png/;
  GetScoreDistribution($plot,$scoretype,$setype,$opt_t,@results);

  #Zoom score plot
  my $plot = $imagefile;
  $plot =~ s/\.png/ZoomScore\.png/;
  ZoomScoreDistribution($plot,$scoretype,$setype,@results);


   if($opt_n)
   {
   #nter plot
   my $plot = $imagefile;
   $plot =~ s/\.png/NterDist\.png/;
   GetNtermPlot($setype,$plot,$opt_t,@results);
   }  
  }

my $size = scalar(@results);
$spectra_count = $size + $spectra_count;
return $spectra_count;

 }

 sub TandemConcat
 {
 my $tandemfile = shift;
 my $spectra_count = shift;

 my $taxonomy= $opt_y;
 my $input;
print("TandemConcat ParseTandem($tandemfile,$input,$taxonomy)\n");
 my($ptr1) = ParseTandem($tandemfile,$input,$taxonomy);
 my @results = @{$ptr1};


  if($opt_z)
  {
  SimpleCount('Tandem',$spectra_count,@results);
  }
 my $scoretype = 'S';
 my $fdrtype = 5;
 my $setype = "T";
 
 #my $imagefile = "TandemConcat" . $spectra_count . ".png";
 my $imagefile = $opt_e; #DCW - from TestMascotAllDecoySearchesInOne

 #Get the actual FDR values for both Gyig and Jones method
 my %all_results = GetFDRValues('C',$opt_w,'0','T',$opt_n,$opt_t,@results);
 
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


  if($opt_I)
  {
  #deltamass plot
  my $plot = $imagefile;
  $plot =~ s/\.png/DeltaMass\.png/;
  GetDeltaMassPlot($plot,$scoretype,$setype,$opt_t,@results);

  #gygi rank plot
  my $plot = $imagefile;
  $plot =~ s/\.png/GygiRank\.png/;
  #GetGygiRankPlot($opt_t,$setype,$plot,'0.05',@results);
  GetGygiRankPlot($opt_t,$setype,$plot,$all_results{$opt_w}[7],@results);#DCW - use FDR threshold, NOT expect value, must add a small amount

  #Score distribution
  my $plot = $imagefile;
  $plot =~ s/\.png/ScoreDist\.png/;
  GetScoreDistribution($plot,$scoretype,$setype,$opt_t,@results);

  #Zoom score plot
  my $plot = $imagefile;
  $plot =~ s/\.png/ZoomScore\.png/;
  ZoomScoreDistribution($plot,$scoretype,$setype,@results);

   if($opt_n)
   {
   #nter plot
   my $plot = $imagefile;
   $plot =~ s/\.png/NterDist\.png/;
   GetNtermPlot($setype,$plot,$opt_t,@results);
   }
  }  

my $size = scalar(@results);
my $spectra_count = $size + $spectra_count;
return $spectra_count;

 }



 sub MascotConcat
 {
print("Starting MascotConcat\n");

 my $mascotfile = shift;
 my $spectra_count = shift;
 my $mgf = $opt_m;


 my $setype = "M";

 my $scoretype = "S";

 #my $imagefile = "MascotConcat" . $spectra_count . ".png";
 my $imagefile = $opt_o; #DCW - from TestMascotAllDecoySearchesInOne

 #get the database
 my $db = GetDbPath($mascotfile);

print("parsing mascot\n");

 my $ptr = ParseMascot($mascotfile,$mgf,$db,'H','0.05');
 my @results = @{$ptr};

print("results size: " . scalar(@results) . "\n");

print("finished parsing mascot\n");

 if($opt_z)
 {
 SimpleCount('Mascot',$spectra_count,@results);
 }

print("Getting FDR values\n");

 #Get the actual FDR values for both Gygi and Jones method
 my %all_results = GetFDRValues('C',$opt_w,'0','M',$opt_n,$opt_t,@results);

my %gygi_results;
my %jones_results;

print("creating hashes\n");

print("all_results size: " . keys(%all_results) . "\n");

 #put the different results into the respective hashes
 foreach my $fdr (keys %all_results)
 {
 $gygi_results{$fdr}[0] = $all_results{$fdr}[0];
 $gygi_results{$fdr}[1] = $all_results{$fdr}[1];
 $jones_results{$fdr}[0] = $all_results{$fdr}[2];
 $jones_results{$fdr}[1] = $all_results{$fdr}[3];
 }
    
print("creating summary file\n");

 #summary file
 my $plot = $imagefile;
 $plot =~ s/\.png/Summary\.txt/;

 AddToCreateSummaryTable('M',$plot,\%gygi_results,\%jones_results);
 AddToPeptideList('M',$plot,\%all_results);

print("creating images\n");

  if($opt_I)
  {
	print("***IMAGE FILE = $imagefile");
  #deltamass plot
  my $plot = $imagefile;
  $plot =~ s/\.png/DeltaMass\.png/;
  GetDeltaMassPlot($plot,$scoretype,$setype,$opt_t,@results);

#print("max_05: ".$all_results{$opt_w}[6].",".$all_results{$opt_w}[7]."\n");
  #gygi rank plot
  my $plot = $imagefile;
  $plot =~ s/\.png/GygiRank\.png/;
 #print "GetGygiRankPlot($opt_t,$setype,$plot,".$all_results{"max_05"}[1].",\@results)\n";
  #GetGygiRankPlot($opt_t,$setype,$plot,'0.05',@results);
  GetGygiRankPlot($opt_t,$setype,$plot,$all_results{$opt_w}[7],@results);#DCW - use FDR threshold, NOT expect value, must add a small amount
  #GetGygiRankPlot($opt_t,$setype,$plot,0.169,@results);#DCW - use FDR threshold, NOT expect value

  #Score distribution
  my $plot = $imagefile;
  $plot =~ s/\.png/ScoreDist\.png/;
  GetScoreDistribution($plot,$scoretype,$setype,$opt_t,@results);

  #Zoom score plot
  my $plot = $imagefile;
  $plot =~ s/\.png/ZoomScore\.png/;
  ZoomScoreDistribution($plot,$scoretype,$setype,@results);

   if($opt_n)
   {
   #nter plot
   my $plot = $imagefile;
   $plot =~ s/\.png/NterDist\.png/;
   GetNtermPlot($setype,$plot,$opt_t,@results);
   }
  }

my $size = scalar(@results);
$spectra_count = $spectra_count + $size;

print("Finishing MascotConcat\n");

return $spectra_count;

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
 my %all_results = GetFDRValues('C',$opt_w,'0','M',$opt_n,$opt_t,@results);
 
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

  if($opt_I)
  {
  #deltamass plot
  my $plot = $imagefile;
  $plot =~ s/\.png/DeltaMass\.png/;
  GetDeltaMassPlot($plot,$scoretype,$setype,$opt_t,@results);

  #gygi rank plot
  my $plot = $imagefile;
  $plot =~ s/\.png/GygiRank\.png/;
  #GetGygiRankPlot($opt_t,$setype,$plot,'0.05',@results);
  GetGygiRankPlot($opt_t,$setype,$plot,$all_results{$opt_w}[7],@results);#DCW - use FDR threshold, NOT expect value, must add a small amount

  #Score distribution
  my $plot = $imagefile;
  $plot =~ s/\.png/ScoreDist\.png/;
  GetScoreDistribution($plot,$scoretype,$setype,$opt_t,@results);

  #Zoom score plot
  my $plot = $imagefile;
  $plot =~ s/\.png/ZoomScore\.png/;
  ZoomScoreDistribution($plot,$scoretype,$setype,@results);

   if($opt_n)
   {
   my $plot = $imagefile;
   $plot =~ s/\.png/NterDist\.png/;
   GetNtermPlot($setype,$plot,$opt_t,@results);
   }
  }

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
    if($results[$s][$r]{'ionscore'} >= $fdr{$f} && $results[$s][$r]{'proteinlist'} !~ m/^$opt_t/ && $results[$s][$r]{'proteinlist'} !~ m/\#\*\#REV/)
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
if($protein =~ m/^$opt_t/)
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
if($protein =~ m/^$opt_t/)
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
    if($results[$s][$r]{'ionscore'} >= $fdr{$f} && $results[$s][$r]{'proteinlist'} !~ m/^$opt_t/ && $results[$s][$r]{'proteinlist'} !~ m/\#\*\
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
my $spectra_count = shift;
$spectra_count--;
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
    #DCW get correct 'start' value for frayed proteins
    my $frayedStart=0;
    if($results[$s][$rank]{'protein'} =~ m/^$opt_t/)
    {
    my @split = split/\_/,$results[$s][$rank]{'protein'};
    $results[$s][$rank]{'protein'} = $split[0] . "_" . $split[1];
    } 
    else
    { 
    my @split = split/\_/,$results[$s][$rank]{'protein'};
    $results[$s][$rank]{'protein'} = $split[0];
	$frayedStart = $split[2]; #DCW
    }

   #DCW
   my $startPos = $results[$s][$rank]{'start'}+$frayedStart;
   $results[$s][$rank]{'start'}+=$frayedStart;
   #print "startPos: " . $startPos . ", fray: ".$frayedStart."\n";

   my @split = split/\|/,$results[$s][$rank]{'protein'};
   $results[$s][$rank]{'protein'} = $split[0];
   my $spec_number = $s + $spectra_count;

   #DCW Extra tab inserted, because mods are expected before 'start'
   #print RES "$results[$s][$rank]{'protein'}\t$results[$s][$rank]{'sequence'}\t$score\t$results[$s][$rank]{'expect'}\t$rank\t$spec_number\t$se\tNterminalaccession1\tNterminal study\t$results[$s][$rank]{'start'}\n";
   #print RES "$results[$s][$rank]{'protein'}\t$results[$s][$rank]{'sequence'}\t$score\t$results[$s][$rank]{'expect'}\t$rank\t$spec_number\t$se\tNterminalaccession1\tNterminal study\t\t$results[$s][$rank]{'start'}\n";
   #mods added
   print RES "$results[$s][$rank]{'protein'}\t$results[$s][$rank]{'sequence'}\t$score\t$results[$s][$rank]{'expect'}\t$rank\t$spec_number\t$se\tNterminalaccession1\tNterminal study\t$results[$s][$rank]{'mods'}\t$startPos\n";




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
print "\t\t-l\tname of comma seperated mascot files - analyse over all\n";
print "\n\t\t-x\tConcatanated Xtandem file\n";
print "\t\t-u\tXTandem forward only file\n";
print "\t\t-v\tXTandem reverse only file\n";
print "\t\t-y\tTaxonomy file associated with the X!Tandem search\n";
print "\t\t-q\tcomma seperated name of tandem files - analyse over all\n";
print "\n\t\t-r\tConcatanated Omssa result file\n";
print "\t\t-p\tOmssa forward only file\n";
print "\t\t-i\tOmssa reverse only file\n";
print "\t\t-a\tname of db single omssa searched\n";
print "\t\t-g\tcomma speparated names of omssa files - analyse over all\n";
print "\n\tAnalysis Settings\n";
print "\t\t-w\tFDR threshold to use\n";
print "\t\t-z\tCombine SearchEngine Results using Kall method\n";
print "\t\t-t\tDecoy tag (default REV_)\n";
print "\t\t-s\toutputfile name\n";
print "\t\t-n\tN-terminal analysis (yes = 1, no = 0 (default)\n";
print "\n\tOutputs\n";
print "\t\t-k\tDo full overlap analysis of replicates and searchengines\n";
print "\t\t-h\tThis help page\n";
print "\t\t-I\tProduceImages\n";

print "\n";

}







