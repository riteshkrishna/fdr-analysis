#!/usr/bin/perl


####################################################################################################
#	 Copyright (C) 2009, Andrew Jones, University of Liverpool									   #
#    If you use this software, please cite the following paper:									   #
#    Jones, A. R., Siepen, J. A., Hubbard, S. J., Paton, N. W., Improving sensitivity in proteome  #
#    studies by analysis of false discovery rates for multiple search engines. 					   #
#    PROTEOMICS 2009, 9, 1220-1229.		   														   #	
#   																							   #
#    																							   #
#    This program is free software: you can redistribute it and/or modify						   #
#    it under the terms of the GNU General Public License as published by						   #
#    the Free Software Foundation, either version 3 of the License, or							   #
#    (at your option) any later version.														   #
#						 																		   #
#    This program is distributed in the hope that it will be useful,							   #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of								   #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the								   #	
#    GNU General Public License for more details.												   #
#																								   #
#    You should have received a copy of the GNU General Public License							   #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.						   #
#   																							   #
#   																							   #
####################################################################################################


use strict;

use FindBin qw($Bin);
use lib "$Bin";	#add current directory to library, allows script to be run from anywhere

use Getopt::Std;

use ExperimentTitle;
use PeptideIdent;
use ConsensusPeptide;
use ConsensusProtein;
use SpectrumIdent;
use ProteinAmbiguityGroup;

		
our $exp_acc;
our($opt_F,$opt_R,$opt_I,$opt_O, $opt_P, $opt_S, $opt_D, $opt_A, $opt_V, $opt_Q, $opt_T);

getopts('F:R:I:O:P:S:D:A:V:Q:T:');


if(!$opt_F || !$opt_O  || !$opt_I){
	print "\n\nCorrect usage:\n";
	print "FDRScore.pl -F [input_file] -O [output_file] -I [string_to_identify_decoy_accessions_eg_Rev]\n";
	print "\nNon mandatory options:\n";
	print "-A [filepath] Local fasta database to extract protein description lines\n";
	print "-D [integer] Ratio of decoys to targets  (default = 1 i.e. [1]decoy:1target)\n";
	print "-P [filepath] output proteins to file \n";
	print "-R [integer] Ranked peptides to accept (default = 1 i.e. only top ranked peptide-spectrum match) \n";
	print "-S [mascot,omssa,X!Tandem] include search engines (default = mascot,omssa,X!Tandem)\n";
	print "-V [filepath] verbose output to file\n";
	print "-Q [filepath] output peptides to a file\n";
	print "-T [float] FDR threshold (default 0.05)\n";

	exit;
}



our $pep_cutoff = 0.5;	# when window FDR for consensus peptides reaches > 0.5 (i.e. at least half are false positives), do not contribute to protein scoring

my $input_file = $opt_F;

our $test_rank = 1;
if($opt_R){
	$test_rank = $opt_R;
}
our $test_se = $opt_S;
our $decoy_ratio = 1;
our $desc_db = $opt_A;		#This is a database to retrieve accession lines from
our %db_seq;



if(!$desc_db){
	print "No database set, no description lines will be reported\n";
}
else{
	
	open(desc_in, "$desc_db") or die "cannot open $desc_db\n";
	while (my $ln = <desc_in>) {
		if($ln =~ s/\>//){
			$ln =~ s/\n//;
			my @tmp = split(/ /, $ln);
			
			my $desc = substr($ln,index($ln," ") + 1);			

			if(index($desc,"product=") != -1){
				
				my $start_pos = index($desc,"product=")+8;
				my $str_len = index($desc,"|",$start_pos) - $start_pos;
				$desc = substr($desc,$start_pos, $str_len);
			}
			
			$db_seq{$tmp[0]} = $desc; 				
		}
	}
}


if($opt_D){
	$decoy_ratio = $opt_D;
}
print "Ratio of targets to decoys is 1:$decoy_ratio\n";
print "Only accepting identifications with rank <= $test_rank\n";

if($test_se){
	print "Only accepting identifications from: $test_se\n";
}

our $rev_string = $opt_I;
my $out_file = $opt_O;

our $prot_out = $opt_P;

our $pep_out = $opt_Q; #JAS added

#JAS
if(!$opt_T) {
$opt_T = 0.05;
}

if(!$test_rank){
	print "Error, no test rank entered!";
	exit;
}
if(!$rev_string){
	print "Error, no reverse string entered!";
	exit;
}
else{
	print "Decoy string: $rev_string\n";
}

if(!$out_file){
	print "Error, no output file entered\n";
	exit;
}



open(input, "<$input_file") or die "cannot open $input_file\n";

our %exp_title_prot_cons;

my $current_etitle = "XXXX";
my $current_se = "XXXX";

our $pepCounter = 0;

my (%peptides_se,%peptides_cons,%peaklists_se,%prots_se,%prots_cons);
our %conflict_prots;

our %pag_mappings;	#mappings from each protein accession to the protein ambiguity group it belongs to
our @prot_ambiguity_groups;	#These two values get assigned to each exp_title_obj



my @temp_results;
my $temp_counter = 0;

open(outfile, ">$out_file");

our $verbose_out = $opt_V;

if($verbose_out){
	open(verbose, ">$verbose_out");
}



if($prot_out){
	open(prot_out, ">$prot_out");
}

#JAS and the peptide file
if($pep_out){
	open(pep_out, ">$pep_out");
}

my $exp_title_obj =  eval { new ExperimentTitle(); }  or die ($@);

while (my $ln = <input>) {
	
	$ln =~ s/\n//;	#remove end of line char
	my @tmp = split(/\t/, $ln);
	my $exp_title = $tmp[8];
	
	my $rank = $tmp[4];
	my $engine = $tmp[6];
	
	my $se_passes = 1;	
	
	if($test_se){
		if(index($test_se,$engine) == -1){
			$se_passes = 0;
		}		
	}
	
	if(@temp_results && $engine ne $current_se){
		#print "New search engine at line $pepCounter \n";
		
		processSearchEngine(\@temp_results,$current_se);		
		
		@temp_results = ();		
		$temp_counter = 0;		
	}

	if($rank <= $test_rank && $se_passes){	#e.g. can set to 1 to only allow top-ranking peptides
		$temp_results[$temp_counter] = $ln;
		$temp_counter++;
	}
	
							
	if($pepCounter != 0 && $exp_title ne $current_etitle){	
		
		print "Processed experiment: $current_etitle\n";
		
		$exp_title_obj->expTitle($current_etitle);
		$exp_title_prot_cons{$current_etitle} = $exp_title_obj;	

		$exp_title_obj->pagMappings(\%pag_mappings);
		$exp_title_obj->pags(\@prot_ambiguity_groups);
		
		$exp_title_obj =  eval { new ExperimentTitle(); }  or die ($@);
		
		%pag_mappings = ();	
		@prot_ambiguity_groups = ();

		print "Processing experiment: $exp_title\n";			
	
	}		
	
	$current_etitle = $exp_title;	
	$pepCounter++;	
	$current_se =$engine;
	
}

#Send the last set of results
if(@temp_results){
	processSearchEngine(\@temp_results,$current_se);
}
print "Processed last experiment: $current_etitle\n";

$exp_title_obj->pagMappings(\%pag_mappings);
$exp_title_obj->pags(\@prot_ambiguity_groups);		
$exp_title_obj->expTitle($current_etitle);
$exp_title_prot_cons{$current_etitle} = $exp_title_obj;

makeConsensus();

close pep_out;


sub makeConsensus(){	
	
	combineAcrossSE();
	
	#Only converts peptides to proteins if a protein file is specified
	if($prot_out){
		resolveConsensusConflicts();
		assignProteinQValues();
	}
	
}

sub processSearchEngine{
	
	my $temp = shift;
	my $search_engine = shift;
	my @results = @{$temp};
	
	print "\tLoading $search_engine results \n";
	
	my (%all_idents, $tmp_count, $exp_title);
	
	foreach my $result (@results){
		
		my @tmp = split(/\t/, $result);
		
		my $peaklist = trim($tmp[5]);
		my $prot_acc = trim($tmp[0]);
		
		if(!$prot_acc){
			print "Error! No protein accession, for result: $result\n";
			exit;
		}
			
		$exp_acc = $tmp[7];	#is the same for all entries	
		my $rank = $tmp[4];
		my $eval = $tmp[3];	
		$exp_title = $tmp[8];		
		my $mod = $tmp[9];		
		my $pep_seq = uc($tmp[1]);	#Omssa puts mods in lower case, convert all peptides to lower case
                my $nter_start = $tmp[10]; #JAS start residue
                $nter_start =~ s/\n//g;  

		
		my $conflict = $peaklist . "_". $rank;					
		my $pep_hit;
		
		#If a PSM has been seen before with the same spectrum number and rank i.e. it is a different peptide2protein mapping
		if($all_idents{$conflict}){
			
			my $prev_hit = $all_idents{$conflict};
			my $prev_seq = $prev_hit->pepSeq();
			my @prev_prot_accs = $prev_hit->getProtAccs();
			my $prev_acc = $prev_prot_accs[0];
			
			#It is the same peptide, mapped to a different protein
			if($prev_hit->evalue() == $eval && $prev_seq eq $pep_seq && $prev_hit->searchEngine() eq $search_engine){				
				
				if($prev_acc eq $prot_acc){
					#print "Identical entries $prev_acc $prot_acc ($pep_seq (e: vs $eval".$prev_hit->pepSeq() . $prev_hit->evalue() .") with same conflict number for search engine: $search_engine\n";
				}
				else{
					my $pag_id = $pag_mappings{$prot_acc};					
									
					if(!defined($pag_id)){	#create new PAG (Protein Ambiguity Group)								
					
						if(defined($pag_mappings{$prev_acc})){
								#Only add new prot acc to PAG mapping if it is not already present
							if(index($prot_ambiguity_groups[$pag_mappings{$prev_acc}],$prot_acc) == -1){ 
								$prot_ambiguity_groups[$pag_mappings{$prev_acc}] .= ",". $prot_acc;	#Add accession to end of PAG string		
								#print "Added to PAG: . " . $prot_ambiguity_groups[$pag_mappings{$prev_acc}] . "\n";
							}
						}
						else{						
							
							my $new_pag = $prev_acc . ",". $prot_acc; 
							$pag_mappings{$prot_acc} = @prot_ambiguity_groups;
							$pag_mappings{$prev_acc} = @prot_ambiguity_groups;
							push(@prot_ambiguity_groups, $new_pag);
							#print "Added new PAG: $new_pag\n";
						}
					}				
					else{
						my $pag_string = $prot_ambiguity_groups[$pag_id];
						if(substr($pag_string, $prot_acc) == -1 || substr($pag_string, $prev_acc) == -1 ){
							die "Fatal error in code, adding proteins to ambiguity group $prot_acc prev: $prev_acc, PAG: $pag_string \n";
						}
					}				
					$prev_hit->addProtAcc($prot_acc);	#add new protein accession
				}
			}
			else{
				die "Fatal error, different peptide ($pep_seq (e: vs $eval".$prev_hit->pepSeq() . $prev_hit->evalue() .") with same PSM number for search engine: $search_engine\n";
			}		
		}
		else{
		
			$pep_hit =  eval { new PeptideIdent(); }  or die ($@);
			$pep_hit->expTitle($exp_title);
			$pep_hit->score($tmp[2]);
			$pep_hit->evalue($eval);
			$pep_hit->pepSeq($pep_seq);
			$pep_hit->rank($rank);		
			$pep_hit->searchEngine($search_engine);
			$pep_hit->peaklist($peaklist);			
			$pep_hit->mods($mod);
			$pep_hit->addProtAcc($prot_acc);
 			$pep_hit->addStart($nter_start);		
	
			#This adds a single peptide to protein mapping
			$all_idents{$conflict} = $pep_hit;
		}		
	}
	
	if($search_engine eq "mascot"){
	
		#Order all idents
		print verbose "Raw peptide results from Mascot from $exp_title\n";
		print verbose "Is Decoy\tSequence\tE-value\tLocal FDR\tQ-value\tCombined FDR Score\tWindow FDR\n";
		my @ordered_idents = @{orderIdents(\%all_idents,1)};	
		
		$exp_title_obj->mascotOrderedIdents(\@ordered_idents);
		$exp_title_obj->mascotIdents(\%all_idents);
	}
	elsif($search_engine eq "omssa"){
	
		print verbose "Raw peptide results from OMSSA from $exp_title\n";
		print verbose "Is Decoy\tSequence\tE-value\tLocal FDR\tQ-value\tCombined FDR Score\tWindow FDR\n";
		my @ordered_idents = @{orderIdents(\%all_idents,1)};
		$exp_title_obj->omssaOrderedIdents(\@ordered_idents);
		$exp_title_obj->omssaIdents(\%all_idents);
	}
	elsif($search_engine  eq "X!Tandem"){
		print verbose "Raw peptide results from X!Tandem from $exp_title\n";
		print verbose "Is Decoy\tSequence\tE-value\tLocal FDR\tQ-value\tCombined FDR Score\tWindow FDR\n";
		my @ordered_idents = @{orderIdents(\%all_idents,1)};
		$exp_title_obj->tandemOrderedIdents(\@ordered_idents);
		$exp_title_obj->tandemIdents(\%all_idents);
	}
	else{
		print "Search engine not recognized: $search_engine, error\n";
		exit;
	}	
	
}


#1. Orders identifications by either e-value (from individual search engines) or by FDR Score (consensus peptides)
#2. Calculates local FDR, q-value and FDR Score
#3. Also calculates the FDR in moving windows - only consensus peptides within a window less than a specified value contribute to proteins
sub orderIdents{

	my $tmp = shift;	
	my $option = shift;	
	my @ordered;
		
	if($option == 1){
		my %idents = %{$tmp};		#put hash into array
		@ordered = ();		
		for my $ident_key ( keys %idents ) {				
			my $pep_hit = $idents{$ident_key};
			push(@ordered,$pep_hit);	
		}
	}
	else{	# Array sent as param
		@ordered = @{$tmp};		
	}	
	
	my $sorted = 0;	
	while(!$sorted){
		
		$sorted = 1;	
		for (my $i=0; $i<(@ordered-1); $i++){			
			my ($first_metric,$second_metric);
			
			my $first_ident =  $ordered[$i];
			my $second_ident = $ordered[$i+1];
			
			if($option == 1){
				$first_metric = $first_ident->evalue();
				$second_metric = $second_ident->evalue();
			}
			else{
				$first_metric = $first_ident->geoMeanFDR();
				$second_metric = $second_ident->geoMeanFDR();
			}					
			
			if($second_metric < $first_metric && $second_ident){
				$ordered[$i] = $second_ident;
				$ordered[$i+1] = $first_ident;
				$sorted = 0;
			}		
		}
	}
	
	my $decoy_count = 0;
	my $all_idents = 0;
	my $prev_fdr = 0;	
	my $last_metric;
	my $window_decoys = 0;
	
	my $window_size;
	if(@ordered / $decoy_ratio < 30){		# This is a very small window so FDR values would have large errors
		$window_size = ((5 * $decoy_ratio) + 5)/2;
	}
	elsif(@ordered / $decoy_ratio < 100){
		$window_size = ((10 * $decoy_ratio) + 10)/2;
	}
	elsif(@ordered / $decoy_ratio < 300){
		$window_size = ((20 * $decoy_ratio) + 20)/2;
	}
	else{
		$window_size = ((30 * $decoy_ratio) + 30)/2;
	}
	
	#print "Using window size $window_size\n";	

	my $window_fdr = 0;
	
	#First run through the loop assigns local FDR values for each ident	
	for(my $i=0; $i<@ordered; $i++){	

		$all_idents++;
		my $ident =  $ordered[$i];
				
		$ident->orderRank($i);	#set the ordered rank of the identification to ensure it is maintained in later re-ordering stage
		
		my @prot_accs = $ident->getProtAccs();		
		
		my $is_decoy = 0;
		my $is_target = 0;		
		
		foreach my $acc (@prot_accs){		
			if(index($acc,$rev_string) == -1){
				$is_target = 1;
			}
			else{
				$is_decoy = 1;
			}
		}		
		
		#This code just calculates a window of FDR, not used in FDR Score algorithm		
		if(@ordered >  $window_size * 2){		
			
			if($i < $window_size-1){#count backwards for half-window size
				if($is_decoy){	$window_decoys++; 	}
			}
			elsif($i == $window_size-1){	#special case, count forward for window size / 2
				for(my $j=$i; $j<=($window_size+$i); $j++){
					my $tmp_ident =  $ordered[$j];
					my @tmp_prot_accs = $tmp_ident->getProtAccs();
					my $acc = $tmp_prot_accs[0];				
					if(index($acc,$rev_string)!=-1){					
						$window_decoys++;
					}				
				}
			}
			elsif($i >= (@ordered - $window_size)){
				#do nothing at the moment, last values 
			}
			else{			#check current-window_size and current+window_size			
				
				#First remove decoys from prev start position, then add from positions + 1
				my $startWin_ident =  $ordered[$i-$window_size];
				my @tmp_prot_accs = $startWin_ident->getProtAccs();
				my $acc = $tmp_prot_accs[0];				
				if(index($acc,$rev_string)!=-1){					
					$window_decoys--;
				}
			
				my $endWin_ident = $ordered[$i+$window_size];
				@tmp_prot_accs = $endWin_ident->getProtAccs();
				$acc = $tmp_prot_accs[0];				
				if(index($acc,$rev_string)!=-1){					
					$window_decoys++;
				}
				
				my $win_targets = ($window_size*2) - $window_decoys;
				my $win_fp  = $window_decoys / $decoy_ratio;
				my $win_tp = $win_targets - $win_fp;
				if($win_tp < 0){$win_tp = 0};	# if TP estimate less than zero, set to zero
				$window_fdr = $win_fp / ($win_tp + $win_fp);		

				#print "$i\t$window_decoys\t$win_targets\t$win_fp\t$win_tp\t$window_fdr\n";				
			}
		}
		
		#If a peptide is matched in both the target and decoy database, it has been set as a decoy, can print out here if needed		
		if($is_decoy && $is_target){
			#print "$pepSeq, pep matched to decoy and target, $pepSeq\n";	
		}		
	
		if($is_decoy){
			$ident->isDecoy($is_decoy);
			$decoy_count++;
		}
	
		my $all_targets = $all_idents - $decoy_count;

		my $fp = $decoy_count / $decoy_ratio;
		my $tp = $all_targets - $fp;		

		if($tp < 0){
			$tp = 0;
		}
		my $fdr = $fp / ( $tp + $fp);
		
		$ident->localFDR($fdr);
		$ident->qValue($fdr);
		$ident->winFDR($window_fdr);		
		
		if($i==(@ordered-1) && $option == 2){				
			$last_metric = $ident->geoMeanFDR();
		}
	}		

			
	my @step_points;	
	my $step_count = 0;
	
	#An artificial decoy hit is added at the end of the series in the combined FDR Score calculation to ensure that no ident as combined FDR Score =0
	if($option == 2){	
		$decoy_count++;
		$all_idents++;
		my $all_targets = $all_idents - $decoy_count;

		my $fp = $decoy_count / $decoy_ratio;
		my $tp = $all_targets - $fp;	
		
		my $fdr = 1;
		
		if(($tp + $fp) != 0){
			$fdr = $fp / ( $tp + $fp);
		}
		
		if($fdr > 1){
			$fdr = 1;
		}

		my $new_ident = eval { new ConsensusPeptide(); }  or die ($@);	
		$new_ident->isDecoy(1);							#Force last ident to be a decoy
		$new_ident->orderRank(@ordered);
		
		$new_ident->geoMeanFDR($last_metric);	#Same geoMeanFDR as worst scoring ident
		$new_ident->qValue($fdr);
		push(@ordered,$new_ident);
		#$ordered[@ordered] = $new_ident;
		my $array_size = @ordered-1;
		$step_points[0] = $array_size;
		$step_count++;
		
	}	
	
	my $current_lowest = 1;
	
	my $prev_q = 100;
	my $prev_metric = 10000;
	
	#Second run through converts FDR to q-values by storing current lowest at all stages	
	for(my $i=(@ordered-2); $i>=0; $i--){

		$all_idents++;
		my $ident =  $ordered[$i];
		my $identFDR = $ident->qValue();		
		my $metric;
		
		if($option ==1){
			$metric = $ident->evalue();
		}
		else{
			$metric = $ident->geoMeanFDR();
		}
				
		if($identFDR < $current_lowest){
			$current_lowest = $identFDR;
		}
		else{
			$ident->qValue($current_lowest);
			$identFDR = $current_lowest;	
		}

		my ($slope, $intercept);
		my $step = 0;
		
		#Definition of a step point is that 
		# ordered - 1 has lower FDR, and ordered + 1 has higher metric
		
		my $step_ident = $ordered[$i+1];
		my $higher_ident =$ordered[$i+2];
		my $is_step = 0;
		
		if($higher_ident){
			if($option == 1){
				
				my $tmp1 = $higher_ident->evalue();
				my $tmp2 = $step_ident->evalue();
				my $tmp3 = 1;
				
				if($tmp2 != 0){
					$tmp3 = $tmp1 / $tmp2;
				}
				
				my $tmp4 = abs($tmp3 - 1);
				
				if($identFDR != $step_ident->qValue() && $tmp4 > 0.000000001){
					$step_points[$step_count] = $i+1;
					$is_step = 1;
				}
			}
			else{
				
				my $tmp1 = $higher_ident->geoMeanFDR();
				my $tmp2 = $step_ident->geoMeanFDR();
				
				my $tmp3 = 1;
				
				if($tmp2 != 0 && $tmp2 ne "0" && $tmp2 && $tmp2 > 0.0000000001){
					$tmp3	= $tmp1 / $tmp2;
				}
				
				my $tmp4 = abs($tmp3 - 1);
				
				#only do regression if there is a more than miniscule difference between the metrics, sometimes rounding errors cause inequality comparison to fail
				if($identFDR != $step_ident->qValue() && $tmp4 > 0.00000001){					
					$step_points[$step_count] = $i+1;
					$is_step = 1;
				}			
			}
		}
		
		if($step_count>0 && $is_step){
			my %vals_to_regress;
			
			my ($curr_step_metric,$prev_step_metric);
			my $curr_step_ident = $ordered[$i+1];
			my $prev_ident = $ordered[$step_points[$step_count-1]];
			
			if($option ==1){
				$curr_step_metric = $curr_step_ident->evalue();
				$prev_step_metric = $prev_ident->evalue();
			}
			else{
				$curr_step_metric = $curr_step_ident->geoMeanFDR();
				$prev_step_metric = $prev_ident->geoMeanFDR();
			}
			
			my $curr_step_q = $curr_step_ident->qValue();
			my $prev_step_q = $prev_ident->qValue();					
								
			$vals_to_regress{$curr_step_metric} = $curr_step_q;
			$vals_to_regress{$prev_step_metric} = $prev_step_q;
		
			($slope, $intercept) = linearRegression(\%vals_to_regress);
			
			my $tmp = $i + 1;
			for(my $j = $i+1; $j < $step_points[$step_count-1]; $j++){
				my $tmp_ident = $ordered[$j];
				my $seq;						
				my $qval = $tmp_ident->qValue();
				my $fdrScore;
				
				if($option ==1){
					$metric =  $tmp_ident->evalue();
					$fdrScore = ($metric *  $slope ) + $intercept;
					$tmp_ident->fdrScore($fdrScore);
				}
				else{
					$metric = $tmp_ident->geoMeanFDR();
					$fdrScore = ($metric *  $slope ) + $intercept;
					$tmp_ident->weightedFDR($fdrScore);
				}						
			}					
		}
		
		if($is_step){
			$step_count++;
		}
		
		if($i==0){	#special case, put 0,0 in the array
			my %vals_to_regress;
			my $curr_step_metric = 0;
			my $curr_step_q = 0;
			
			my $prev_ident = $ordered[$step_points[$step_count-1]];
			
			my $prev_step_metric;
			my $prev_step_q = $prev_ident->qValue();
			
			if($option ==1){
				$prev_step_metric = $prev_ident->evalue();
			}
			else{
				$prev_step_metric = $prev_ident->geoMeanFDR();
			}			
				
			$vals_to_regress{$curr_step_metric} = $curr_step_q;
			$vals_to_regress{$prev_step_metric} = $prev_step_q;			
			
			($slope, $intercept) = linearRegression(\%vals_to_regress);		

			my $end = $step_points[@step_points -1];
		
			for(my $j = 0; $j < $end; $j++){
				my $tmp_ident = $ordered[$j];
				my $qval = $tmp_ident->qValue();
				my ($seq,$fdrScore); 
				
				my $metric;
				if($option ==1){
					$metric = $tmp_ident->evalue();
					$seq = $tmp_ident->pepSeq();

					$fdrScore = ($metric *  $slope );	#no intercept needed this goes thru origin 
					$tmp_ident->fdrScore($fdrScore);
				}
				else{
					$metric = $tmp_ident->geoMeanFDR();
					$seq = $tmp_ident->getPepSeq();

					$fdrScore = ($metric *  $slope );	#no intercept needed this goes thru origin 
					$tmp_ident->weightedFDR($fdrScore);
				}
				
				$ordered[$j] = 	$tmp_ident;					
				#print "$seq\t$metric\t$qval\t$fdrScore\n";
			}
			$step_count++;
		}		
		
		$prev_q = $identFDR;
		$prev_metric = $metric;
	}	
	
	#Special case to assign q-values equal to FDR scores for last values on first run through, on second run through these are be set okay by regression method
	if($option ==1){
		for(my $j = $step_points[0]; $j < @ordered; $j++){
			my $tmp_ident = $ordered[$j];
			my ($seq,$metric,$fdrScore); 
			my $qval = $tmp_ident->qValue();
			
			$metric = $tmp_ident->evalue();
			$seq = $tmp_ident->pepSeq();
			$fdrScore = $qval;
			$tmp_ident->fdrScore($fdrScore);
		}	
	}
	
	for(my $i=0; $i<@ordered; $i++){
		my $tmp_ident = $ordered[$i];
		my ($seq,$localFDR,$metric,$qval,$fdrScore);
		
		my $is_decoy = $tmp_ident->isDecoy();
		if(!$is_decoy){
			$is_decoy = 0;
		}
		
		my $win_fdr = 0;
		
		if($option == 1){
			$seq = $tmp_ident->pepSeq();
			$localFDR = $tmp_ident->localFDR();
			$metric = $tmp_ident->evalue();
			$qval = $tmp_ident->qValue();
			$fdrScore =$tmp_ident->fdrScore();
			$win_fdr = $tmp_ident->winFDR();
		}
		else{
			$seq = $tmp_ident->getPepSeq();
			$localFDR = $tmp_ident->localFDR();
			$metric = $tmp_ident->geoMeanFDR();
			$qval = $tmp_ident->qValue();
			$fdrScore =$tmp_ident->weightedFDR();
			$win_fdr = $tmp_ident->winFDR();
		}	
	
		print verbose "$is_decoy\t$seq\t$metric\t$localFDR\t$qval\t$fdrScore\t$win_fdr\n";
	}
		
	print verbose "\n\n";	
	
	pop(@ordered);			#remove last item that is an artificial decoy added to normalise sets with no decoys in...	
	
	return \@ordered;
	#Third run through for printing and performing regression
	#print "\n\n\nOrdered peptides\n";
}

sub orderFinalPeptides{

	my $tmp = shift;
	my @ordered = @{$tmp};		
	
	my $sorted = 0;
	
	while(!$sorted){
		
		$sorted = 1;
	
		for (my $i=0; $i<(@ordered-1); $i++){
			
			my ($first_metric,$second_metric);
			
			my $first_ident =  $ordered[$i];
			my $second_ident = $ordered[$i+1];
			
			$first_metric = $first_ident->weightedFDR();
			$second_metric = $second_ident->weightedFDR();				
			
			if($second_metric < $first_metric && $second_ident){

				$ordered[$i] = $second_ident;
				$ordered[$i+1] = $first_ident;
				$sorted = 0;
			}		
		}
	}

	return \@ordered;
}



#return slope and intercept of best fit line for values of x and y
sub linearRegression{

	my $tmp = shift;	

	
	my %values = %{$tmp};
	
	my $slope = 0;
	my $intercept = 0;
	
	my $sumProductXY = 0;
	my $sumX = 0;
	my $sumY = 0;
	my $sumSquareX = 0;
	
	#for my $key  (keys %values ){
	#	my $e = $values{$key};				
	#	print verbose "$e\t$key\n";				
	#}	

	if(keys %values > 1){
	
		my $n = 0;
		for my $x ( keys %values ) {
			
			my $y = $values{$x};
			
			if($y ne "INF"){
				$sumX = $sumX + $x;
				$sumY = $sumY + $y;
				$sumSquareX = $sumSquareX + ($x * $x);
				$sumProductXY = $sumProductXY + ($x * $y);
				$n++;
			}
	    }	
			
		#print verbose "SumX $sumX\n";
		#print verbose "SumY $sumY\n";
		#print verbose "Sum squareX: $sumSquareX\n";
		#print verbose "Sum productXY: $sumProductXY\n";
		
		my $numerator = ($n * $sumProductXY) - ($sumX * $sumY);
		my $denominator = ($n * $sumSquareX) - ($sumX * $sumX);

		if($denominator == 0){
			if($numerator ==0){			
				$slope = 0;
			}
			elsif($numerator ==1){
				$slope = 1;
			}
			else{
				print "Error, numerator: $numerator\n";
				
				for my $x ( keys %values ) {
					my $y = $values{$x};
					print "$x --> $y\n";
				}
				
				$slope = 1;
			}
		}
		else{
			$slope = $numerator / $denominator;
		}

		$intercept = ($sumY - ($slope * $sumX)) / $n;
	}
	else{
		print verbose "Error: insufficient values for linear regression\n";
		exit;
	}
	
	#print "Slope: $slope\n";
	#print "Intercept: $intercept\n";	
	
	return ($slope, $intercept);
}




sub combineAcrossSE{

	my ($total_tp_001,$total_tp_005,$total_tp_010,$total_tp_025);
	my $total_fp_001 = 0;
	my $total_fp_005 = 0;
	my $total_fp_010 = 0;
	my $total_fp_025 = 0;
	
	my ($mascot_total_tp05,$mascot_total_fp05 ,$omssa_total_tp05,$omssa_total_fp05,$tandem_total_tp05,$tandem_total_fp05);
	my ($mascot_total_tp01,$mascot_total_fp01 ,$omssa_total_tp01,$omssa_total_fp01,$tandem_total_tp01,$tandem_total_fp01);

	for my $exp_title ( keys %exp_title_prot_cons ) {	
	
		print "Combining results from experiment: $exp_title \n";
		print verbose "\n\n\n********************************************\nExperiment: $exp_title\n";
	
		my $expt_res = $exp_title_prot_cons{$exp_title};	
	
		my (%mascot_idents, %omssa_idents, %tandem_idents);
	
		if($expt_res->mascotIdents()){
			%mascot_idents = %{$expt_res->mascotIdents()};
		}
		else{
			print "No Mascot results for $exp_title\n";
		}
		
		if($expt_res->omssaIdents()){
			%omssa_idents = %{$expt_res->omssaIdents()};
		}
		else{
			print "No Omssa results for $exp_title\n";
		}
		
		if($expt_res->tandemIdents()){
			%tandem_idents = %{$expt_res->tandemIdents()};
		}
		else{
			print "No Tandem results for $exp_title\n";
		}		
	
		my %spectrum_idents;
				
		for my $m_ident_key ( keys %mascot_idents ) {
			my $s_ident;						
						
			if($spectrum_idents{$m_ident_key}){
				$s_ident = $spectrum_idents{$m_ident_key};
			}
			else{
				$s_ident = eval { new SpectrumIdent(); }  or die ($@);
			}
			
			$s_ident->mascotPeptide($mascot_idents{$m_ident_key});
			$spectrum_idents{$m_ident_key} = $s_ident;
		}
				
		for my $o_ident_key ( keys %omssa_idents ) {
			my $s_ident;
						
			if($spectrum_idents{$o_ident_key}){
				$s_ident = $spectrum_idents{$o_ident_key};
			}
			else{
				$s_ident = eval { new SpectrumIdent(); }  or die ($@);
			}
			
			$s_ident->omssaPeptide($omssa_idents{$o_ident_key});
			$spectrum_idents{$o_ident_key} = $s_ident;
		}
		
		for my $t_ident_key ( keys %tandem_idents ) {
			my $s_ident;
						
			if($spectrum_idents{$t_ident_key}){
				$s_ident = $spectrum_idents{$t_ident_key};
			}
			else{
				$s_ident = eval { new SpectrumIdent(); }  or die ($@);
			}
			
			$s_ident->tandemPeptide($tandem_idents{$t_ident_key});		
			$spectrum_idents{$t_ident_key} = $s_ident;
		}		
		
		my (@m_idents,@o_idents,@t_idents,@mo_idents,@mt_idents,@ot_idents,@mot_idents);
		
		for my $s_ident_key ( keys %spectrum_idents ) {
		
			my $s_ident = $spectrum_idents{$s_ident_key};
			my $m_ident = $s_ident->mascotPeptide();
			my $o_ident = $s_ident->omssaPeptide();
			my $t_ident = $s_ident->tandemPeptide();			
			
			my($m_pep,$o_pep,$t_pep,$geomean_fdrScore,$mean_rank);
			
			my $aggFDRScore = 1;
			my $se_count = 0;
			
			if($m_ident){
				$m_pep = $m_ident->pepSeq();
			}
			if($o_ident){
				$o_pep = $o_ident->pepSeq();
			}
			if($t_ident){
				$t_pep = $t_ident->pepSeq();
			}			
			
			if($m_pep && $m_pep eq $o_pep && $o_pep eq $t_pep){
				
				$aggFDRScore = ($m_ident->fdrScore()) * ($o_ident->fdrScore()) * ($t_ident->fdrScore());
				$geomean_fdrScore = $aggFDRScore ** (1/3);				
				$mean_rank = ($m_ident->orderRank() + $o_ident->orderRank() +  $t_ident->orderRank())/3;			
				
				my $cons_ident = eval { new ConsensusPeptide(); }  or die ($@);			
				$cons_ident->geoMeanFDR($geomean_fdrScore); 
				$cons_ident->meanRank($mean_rank);				
				$cons_ident->assignPeptide($m_ident);
				$cons_ident->assignPeptide($o_ident);
				$cons_ident->assignPeptide($t_ident);
				$cons_ident->spectrumID($s_ident_key);
				$cons_ident->isDecoy($m_ident->isDecoy());

                                my $start_res = $m_ident->addStart();
 				$cons_ident->addStart($start_res);

				push(@mot_idents, $cons_ident);

					
			}
			
			if($m_pep && $m_pep eq $o_pep && $o_pep ne $t_pep){
			
				$aggFDRScore = ($m_ident->fdrScore()) * ($o_ident->fdrScore());
				$geomean_fdrScore = $aggFDRScore ** (1/2);				
				$mean_rank = ($m_ident->orderRank() + $o_ident->orderRank())/2;
				
				
				my $cons_ident = eval { new ConsensusPeptide(); }  or die ($@);			
				$cons_ident->geoMeanFDR($geomean_fdrScore);
				$cons_ident->assignPeptide($m_ident);
				$cons_ident->assignPeptide($o_ident);
				$cons_ident->isDecoy($m_ident->isDecoy());
				$cons_ident->spectrumID($s_ident_key);				
				$cons_ident->meanRank($mean_rank);
                                
      				my $start_res = $m_ident->addStart();
				$cons_ident->addStart($start_res);
			
				push(@mo_idents, $cons_ident);				
			}
			
			if($m_pep && $m_pep eq $t_pep && $t_pep ne $o_pep){
				$aggFDRScore = ($m_ident->fdrScore()) * ($t_ident->fdrScore());			
				$geomean_fdrScore = $aggFDRScore ** (1/2);				
				$mean_rank = ($m_ident->orderRank() + $t_ident->orderRank())/2;				
				
				my $cons_ident = eval { new ConsensusPeptide(); }  or die ($@);			
				$cons_ident->geoMeanFDR($geomean_fdrScore);
				$cons_ident->meanRank($mean_rank);
								
				$cons_ident->assignPeptide($m_ident);
				$cons_ident->assignPeptide($t_ident);
				$cons_ident->spectrumID($s_ident_key);
				$cons_ident->isDecoy($m_ident->isDecoy());	

				my $start_res = $m_ident->addStart();
                                $cons_ident->addStart($start_res);

				push(@mt_idents, $cons_ident);
							
			}
			
			if($o_pep && $o_pep eq $t_pep && $o_pep ne $m_pep){
				$aggFDRScore = ($o_ident->fdrScore()) * ($t_ident->fdrScore());
				$geomean_fdrScore = $aggFDRScore ** (1/2);
				$mean_rank = ($o_ident->orderRank() + $t_ident->orderRank())/2;			

				
				my $cons_ident = eval { new ConsensusPeptide(); }  or die ($@);			
				$cons_ident->geoMeanFDR($geomean_fdrScore);
				$cons_ident->meanRank($mean_rank);
				$cons_ident->assignPeptide($o_ident);
				$cons_ident->assignPeptide($t_ident);
				$cons_ident->spectrumID($s_ident_key);
				$cons_ident->isDecoy($o_ident->isDecoy());	
				
				my $start_res = $o_ident->addStart();
                                $cons_ident->addStart($start_res);

				push(@ot_idents, $cons_ident);
			
			}
			
			if($m_pep && $m_pep ne $o_pep &&  $m_pep ne $t_pep){
				$geomean_fdrScore = $m_ident->fdrScore();
				$mean_rank = $m_ident->orderRank();

				if(!$geomean_fdrScore){
					my $pep = $m_ident->pepSeq();
					$geomean_fdrScore = 1;			

				}		
				
				my $cons_ident = eval { new ConsensusPeptide(); }  or die ($@);			
				$cons_ident->geoMeanFDR($geomean_fdrScore);
				$cons_ident->meanRank($mean_rank);
				
				$cons_ident->assignPeptide($m_ident);
				$cons_ident->isDecoy($m_ident->isDecoy());	
				$cons_ident->spectrumID($s_ident_key);
                                my $start_res = $m_ident->addStart();
                                $cons_ident->addStart($start_res);

				push(@m_idents, $cons_ident);


								
			}
			
			if($o_pep && $o_pep ne $m_pep &&  $o_pep ne $t_pep){
				$geomean_fdrScore = $o_ident->fdrScore();
				$mean_rank = $o_ident->orderRank();						

				if(!$geomean_fdrScore){
					my $pep = $o_ident->pepSeq();
					$geomean_fdrScore = 1;					
				}
				
				my $cons_ident = eval { new ConsensusPeptide(); }  or die ($@);			
				$cons_ident->geoMeanFDR($geomean_fdrScore);
				$cons_ident->meanRank($mean_rank);				
				$cons_ident->assignPeptide($o_ident);
				$cons_ident->isDecoy($o_ident->isDecoy());				
				$cons_ident->spectrumID($s_ident_key);
				
                                my $start_res = $o_ident->addStart();
                                $cons_ident->addStart($start_res);

				push(@o_idents, $cons_ident);				
			}
			
			if($t_pep && $t_pep ne $m_pep &&  $t_pep ne $o_pep){

				$geomean_fdrScore = $t_ident->fdrScore();
		
				if(!$geomean_fdrScore){
					my $pep = $t_ident->pepSeq();
					$geomean_fdrScore = 1;				
				}
			
				$mean_rank = $t_ident->orderRank();		
				
				my $cons_ident = eval { new ConsensusPeptide(); }  or die ($@);			
				$cons_ident->geoMeanFDR($geomean_fdrScore);
				$cons_ident->meanRank($mean_rank);				
				$cons_ident->assignPeptide($t_ident);
				$cons_ident->spectrumID($s_ident_key);
				$cons_ident->isDecoy($t_ident->isDecoy());
				
				my $start_res = $t_ident->addStart();
                                $cons_ident->addStart($start_res);

				push(@t_idents, $cons_ident);
		

			}
						
		}
		
		print verbose "\n\n\nTandem only\n";
		print verbose "Is Decoy\tSequence\tMean FDR Score\tLocal FDR\tQ-value\tCombined FDR Score\tWindow FDR\n";
		@t_idents = @{orderIdents(\@t_idents,2)};
		print verbose "\n\n\nMascot only\n";
		print verbose "Is Decoy\tSequence\tMean FDR Score\tLocal FDR\tQ-value\tCombined FDR Score\tWindow FDR\n";		
		@m_idents = @{orderIdents(\@m_idents,2)};
		print verbose "\n\n\nOmssa only\n";
		print verbose "Is Decoy\tSequence\tMean FDR Score\tLocal FDR\tQ-value\tCombined FDR Score\tWindow FDR\n";		
		@o_idents = @{orderIdents(\@o_idents,2)};
		print verbose "\n\n\nOmssa Tandem \n";
		print verbose "Is Decoy\tSequence\tMean FDR Score\tLocal FDR\tQ-value\tCombined FDR Score\tWindow FDR\n";		
		@ot_idents = @{orderIdents(\@ot_idents,2)};
		print verbose "\n\n\nMascot Tandem \n";
		print verbose "Is Decoy\tSequence\tMean FDR Score\tLocal FDR\tQ-value\tCombined FDR Score\tWindow FDR\n";		
		@mt_idents = @{orderIdents(\@mt_idents,2)};
		print verbose "\n\n\nMascot Omssa \n";	
		print verbose "Is Decoy\tSequence\tMean FDR Score\tLocal FDR\tQ-value\tCombined FDR Score\tWindow FDR\n";	
		@mo_idents = @{orderIdents(\@mo_idents,2)};
		print verbose "\n\n\nMascot Omssa Tandem \n";
		print verbose "Is Decoy\tSequence\tMean FDR Score\tLocal FDR\tQ-value\tCombined FDR Score\tWindow FDR\n";		
		@mot_idents = @{orderIdents(\@mot_idents,2)};
		
		print outfile "**************************************\n";
		print outfile "**************************************\n";
		print outfile "Peptide results for assay $exp_title \n";
		
		my (@mascot_idents,@omssa_idents,@tandem_idents);
		
		if($expt_res->mascotOrderedIdents()){
			@mascot_idents = @{$expt_res->mascotOrderedIdents()};
		}
		if($expt_res->omssaOrderedIdents()){
			@omssa_idents = @{$expt_res->omssaOrderedIdents()};
		}
		if($expt_res->tandemOrderedIdents()){
			@tandem_idents = @{$expt_res->tandemOrderedIdents()};
		}
		
		print outfile "\n\n********************\nPeptide Summary Stats (see verbose output for sequences)\n\n";
		print outfile "FDR Threshold\tTrue positives\tFalse Positives\tCalculate FDR\tSearch engine\n";
		my ($mascot_tp01,$mascot_fp01) = printPeptideResults(\@mascot_idents,$opt_T,"Mascot");
		#my ($mascot_tp05,$mascot_fp05) = printPeptideResults(\@mascot_idents,0.05,"Mascot");
		my ($omssa_tp01,$omssa_fp01) = printPeptideResults(\@omssa_idents,$opt_T,"Omssa");
		#my ($omssa_tp05,$omssa_fp05) = printPeptideResults(\@omssa_idents,0.05,"Omssa");
		my ($tandem_tp01,$tandem_fp01) = printPeptideResults(\@tandem_idents,$opt_T,"X!Tandem");
		#my ($tandem_tp05,$tandem_fp05) = printPeptideResults(\@tandem_idents,0.05,"X!Tandem");
		
		$mascot_total_tp01 +=$mascot_tp01;
		$mascot_total_fp01 +=$mascot_fp01;
		$omssa_total_tp01 +=$omssa_tp01;
		$omssa_total_fp01 +=$omssa_fp01;
		$tandem_total_tp01 +=$tandem_tp01;
		$tandem_total_fp01 +=$tandem_fp01;
		
		#$mascot_total_tp05 +=$mascot_tp05;
		#$mascot_total_fp05 +=$mascot_fp05;
		#$omssa_total_tp05 +=$omssa_tp05;
		#$omssa_total_fp05 +=$omssa_fp05;
		#$tandem_total_tp05 +=$tandem_tp05;
		#$tandem_total_fp05 +=$tandem_fp05;		
		
		my @final_peptides;
		push(@final_peptides,@mot_idents);
		push(@final_peptides,@mo_idents);
		push(@final_peptides,@ot_idents);
		push(@final_peptides,@mt_idents);
		push(@final_peptides,@m_idents);
		push(@final_peptides,@o_idents);
		push(@final_peptides,@t_idents);

		@final_peptides = @{orderFinalPeptides(\@final_peptides)};	
		my ($tmp1,$tmp2,$tmp3) = printFinalPeptides(\@final_peptides);

		my @tp_results = @{$tmp1};
		my @fp_results = @{$tmp2};
		my %all_proteins = %{$tmp3};	
		
		$total_tp_001 += $tp_results[0];
		$total_fp_001 += $fp_results[0];
		$total_tp_005 += $tp_results[1];
		$total_fp_005 += $fp_results[1];
		$total_tp_010 += $tp_results[2];
		$total_fp_010 += $fp_results[2];
		$total_tp_025 += $tp_results[3];
		$total_fp_025 += $fp_results[3];		
	
		my $num = keys %all_proteins;
		$expt_res->consensusProteins(\%all_proteins);		
	}
	
	
	my $fdr_001 = $total_fp_001 / ($total_fp_001 + $total_tp_001);
	my $fdr_005 = $total_fp_005 / ($total_fp_005 + $total_tp_005);
	
	my $fdr_010 = $total_fp_010 / ($total_fp_010 + $total_tp_010);
	my $fdr_025 = $total_fp_025 / ($total_fp_005 + $total_tp_025);
	
	my ($mascot_fdr01,$omssa_fdr01,$tandem_fdr01,$mascot_fdr05,$omssa_fdr05,$tandem_fdr05);
	
	if($mascot_total_fp01 + $mascot_total_tp01 != 0){
		$mascot_fdr01 = $mascot_total_fp01 / ($mascot_total_fp01 + $mascot_total_tp01);
	}
	if($omssa_total_fp01 + $omssa_total_tp01 != 0){
		$omssa_fdr01 = $omssa_total_fp01 / ($omssa_total_fp01 + $omssa_total_tp01);
	}
	if($tandem_total_fp01 + $tandem_total_tp01 != 0){
		$tandem_fdr01 = $tandem_total_fp01 / ($tandem_total_fp01 + $tandem_total_tp01);	
	}

	if($mascot_total_fp05 + $mascot_total_tp05 != 0){
		$mascot_fdr05 = $mascot_total_fp05 / ($mascot_total_fp05 + $mascot_total_tp05);
	}
	if($omssa_total_fp05 + $omssa_total_tp05 != 0){
		$omssa_fdr05 = $omssa_total_fp05 / ($omssa_total_fp05 + $omssa_total_tp05);
	}
	if($tandem_total_fp05 + $tandem_total_tp05 != 0){
		$tandem_fdr05 = $tandem_total_fp05 / ($tandem_total_fp05 + $tandem_total_tp05);
	}
	
	
	
	print outfile "\n\n*****************************\n*****************************\n*****************************\nSummary Results across all files\n";
	print outfile "Threshold\tTP\tFP\tFDR\tSearch engine\n";
	
	print outfile "0.01\t$mascot_total_tp01\t$mascot_total_fp01\t$mascot_fdr01\tMascot\n";
	print outfile "0.01\t$omssa_total_tp01\t$omssa_total_fp01\t$omssa_fdr01\tOmssa\n";
	print outfile "0.01\t$tandem_total_tp01\t$tandem_total_fp01\t$tandem_fdr01\tTandem\n";
	
	print outfile "0.05\t$mascot_total_tp05\t$mascot_total_fp05\t$mascot_fdr05\tMascot\n";
	print outfile "0.05\t$omssa_total_tp05\t$omssa_total_fp05\t$omssa_fdr05\tOmssa\n";
	print outfile "0.05\t$tandem_total_tp05\t$tandem_total_fp05\t$tandem_fdr05\tTandem\n";
	
	print outfile "\n";
	print outfile "Consensus peptide summary:\n";
	print outfile "Threshold\tTP\tFP\tFDR\n";
	print outfile "0.01\t$total_tp_001\t$total_fp_001\t$fdr_001\n";
	print outfile "0.05\t$total_tp_005\t$total_fp_005\t$fdr_005\n";
	print outfile "0.1\t$total_tp_010\t$total_fp_010\t$fdr_010\n";
	print outfile "0.25\t$total_tp_025\t$total_fp_025\t$fdr_025\n";
	
	
	print outfile "\nPeptideSummary\nExp acc\tThreshold\tMascot TP\tOmssa TP\tTandem TP\tCombined TP\n";
	print outfile "$exp_acc\t0.01\t$mascot_total_tp01\t$omssa_total_tp01\t$tandem_total_tp01\t$total_tp_001\n";
	print outfile "$exp_acc\t0.05\t$mascot_total_tp05\t$omssa_total_tp05\t$tandem_total_tp05\t$total_tp_005\n";	

	#print outfile "\n\nProteinResults\n";
	
}


sub addProtein{

	my ($pep_hit,$tmp_all_proteins,$conflict) = @_;
	
	my %all_proteins = %{$tmp_all_proteins};
	
	
	my @prot_accs = $pep_hit->getProtAccs();		
	my $pep_seq = $pep_hit->getPepSeq();	
	
	foreach my $prot_acc (@prot_accs){
	
	#$all_idents{$conflict} = $pep_hit;
	
		#print "addprotein: $prot_acc for: $conflict \n";
	
		my $cons_peptide;
		my %prot_pep_idents;				
		
		if($all_proteins{$prot_acc}){	
			%prot_pep_idents = %{$all_proteins{$prot_acc}};			
		}
		else{
			%prot_pep_idents = ();
		}


		#if($prot_pep_idents{$conflict}){			
		#	$cons_peptide = $prot_pep_idents{$conflict};
		#}
		#else{
		#	$cons_peptide = eval { new ConsensusPeptide(); }  or die ($@);	
		#}		

		$prot_pep_idents{$conflict}=$pep_hit;				
		$all_proteins{$prot_acc} = \%prot_pep_idents;
		
		#print "$pep_seq\t$conflict\t$prot_acc\n";
	}

	return \%all_proteins;	
	
}


sub printPeptideResults{

	my ($tmp1,$threshold,$se) = @_;
	my @results = @{$tmp1};
	
	my ($all_idents,$decoy_count,$tp,$fp,$fdr,$prev_tp,$prev_fp,$prev_fdr);		

	for (my $i = 0; $i<@results; $i++){	

		my $ident = $results[$i];			
		my $is_decoy = $ident->isDecoy();			
		$all_idents++;
		
		if($is_decoy){
			$decoy_count++;
		}			
		my $all_targets = $all_idents - $decoy_count;
		$fp = $decoy_count / $decoy_ratio;
		$tp = $all_targets - $fp;	
		
		if($tp < 0){
			$tp = 0;
		}
		
		if($tp == 0 && $fp ==0){
			$fdr = 0;
		}
		else{	
			$fdr = $fp / ( $tp + $fp);
		}			
		
		if($fdr >= $threshold){
			$i= @results;					
		}
		else{
			$prev_tp = $tp;
			$prev_fp = $fp;
			$prev_fdr = $fdr;
		}		
	}

	print outfile "$threshold\t$prev_tp\t$prev_fp\t$prev_fdr\t$se\n";	
	return ($prev_tp,$prev_fp);	
}


#All peptides have now been ordered according to FDRScore
sub printFinalPeptides{

	my ($tmp1) = @_;	
	my @result_set = @{$tmp1};
	my %all_proteins = ();
	
	print outfile "\nConsensus peptides (FDR < 0.25)\n********************\n\n";
	print verbose "\nConsensus peptides\n********************\n\n";
	print outfile "Is decoy\tSequence\tSearch engines\tProtein accessions\tFDR Score\tTrue positives\tFalse positives\tCumulative FDR\tWindow FDR\n";

	
	my @tp_results;
	my @fp_results;	
	
	my $total_tp = 0;
	my $total_fp = 0;
	my $total_fdr = 0;
	
	#my @test_thresholds = (0.01,0.05,0.1,0.25,1);
	my @test_thresholds = (0.01,0.05,0.1,0.25);		
		
	my $all_idents = 0;
	my $decoy_count = 0;
	
	my $count = 0;
	my $fdr = 0;
	my $prev_fdr = 0;
	my ($prev_tp,$prev_fp);
	my $t = 0;
	my $current_threshold = $test_thresholds[$t];
	
	my $window_decoys = 0;
	
	my $window_size;
	
	if(@result_set / $decoy_ratio < 30){		# This is a very small window so FDR values would have large errors
		$window_size = ((5 * $decoy_ratio) + 5)/2;
	}
	elsif(@result_set / $decoy_ratio < 100){
		$window_size = ((10 * $decoy_ratio) + 10)/2;
	}
	elsif(@result_set / $decoy_ratio < 300){
		$window_size = ((20 * $decoy_ratio) + 20)/2;
	}
	else{
		$window_size = ((30 * $decoy_ratio) + 30)/2;
	}
	
	#print "Using window size $window_size\n\n";
	
	#my $window_size = ((10 * $decoy_ratio) + 10)/2;		#real window is actually double this value, this is used for calculating backwards and forward positions
		
	my $window_fdr = 0;
	
	for (my $i = 0; $i<@result_set; $i++){
		
		my $ident = $result_set[$i];
		#my $geo_score = $ident->geoMeanFDR();
		#my $qval = $ident->qValue();	#These values are irrelevant once all 7 sets have been combined
		my $weightFDR = $ident->weightedFDR();
		my $is_decoy = $ident->isDecoy();		
		my $specID = $ident->spectrumID();		
		my @prot_accs = $ident->getProtAccs();
		
		$all_idents++;
		if($is_decoy){
			$decoy_count++;
		}
   		  if($weightFDR<=$opt_T && $pep_out)
                  {
 		  my $pep_seq = $ident->getPepSeq();
                  print pep_out "$opt_T\t$pep_seq\t$weightFDR\t" . $ident->addStart() . "\t\n"; 
                  }
		
		my $all_targets = $all_idents - $decoy_count;
		my $fp = $decoy_count / $decoy_ratio;
		my $tp = $all_targets - $fp;	
		
		if($tp < 0){
			$tp = 0;
		}
		if($tp + $fp != 0){
			$fdr = $fp / ( $tp + $fp);
		}
		
		if($i < $window_size-1){#count backwards for half-window size
			if($is_decoy){		
				$window_decoys++;
			}
		}
		elsif($i == $window_size-1){	#special case, count forward for window size / 2
			for(my $j=$i; $j<=($window_size+$i); $j++){
				my $tmp_ident =  $result_set[$j];
				my @tmp_prot_accs = $tmp_ident->getProtAccs();
				my $acc = $tmp_prot_accs[0];				
				if(index($acc,$rev_string)!=-1){					
					$window_decoys++;
				}				
			}
		}
		elsif($i >= (@result_set - $window_size)){
			#do nothing at the moment, last values 
		}
		else{			#check current-window_size and current+window_size			
			
			#First remove decoys from prev start position, then add from positions + 1
			my $startWin_ident =  $result_set[$i-$window_size];
			my @tmp_prot_accs = $startWin_ident->getProtAccs();
			my $acc = $tmp_prot_accs[0];				
			if(index($acc,$rev_string)!=-1){					
				$window_decoys--;
			}
		
			my $endWin_ident = $result_set[$i+$window_size];
			@tmp_prot_accs = $endWin_ident->getProtAccs();
			$acc = $tmp_prot_accs[0];				
			if(index($acc,$rev_string)!=-1){					
				$window_decoys++;
			}
			
			my $win_targets = ($window_size*2) - $window_decoys;
			my $win_fp  = $window_decoys / $decoy_ratio;
			my $win_tp = $win_targets - $win_fp;
			if($win_tp < 0){$win_tp = 0};	# if TP estimate less than zero, set to zero
			$window_fdr = $win_fp / ($win_tp + $win_fp);
			
			#print "$i\t$window_decoys\t$win_targets\t$win_fp\t$win_tp\t$window_fdr\n";				
		}
		
		
		if($window_fdr < $pep_cutoff){	# only add to protein array if below cutoff
			%all_proteins = %{addProtein($ident,\%all_proteins,$specID)};	
		}
			
		my $pep_seq = $ident->getPepSeq();
		my $win_fdr1 = $ident->winFDR();	#This value is fairly useless due to the small size of some results sets
		my $win_fdr2 = $window_fdr;
		
		if(!$is_decoy){$is_decoy = 0;}	#for printing, is not decoy, set to 0

		my $se_idents = "";
		if($ident->mascotPeptide()){
			$se_idents .= "m";
		}
		if($ident->omssaPeptide()){
			$se_idents .= "o";
		}
		if($ident->tandemPeptide()){
			$se_idents .= "t";
		}
		
		if($weightFDR < 0.25){
			print outfile "$is_decoy\t$pep_seq\t$se_idents\t@prot_accs\t$weightFDR\t$tp\t$fp\t$fdr\t$win_fdr2\n";
		}
		
		print verbose "$is_decoy\t$pep_seq\t$se_idents\t@prot_accs\t$weightFDR\t$tp\t$fp\t$fdr\t$win_fdr2\n";			
		
		if($weightFDR >= $current_threshold || $i==@result_set-1){ 	
			push(@tp_results,$prev_tp);
			push(@fp_results,$prev_fp);
			$t++;
			$current_threshold = $test_thresholds[$t];
		}
		$prev_tp = $tp;
		$prev_fp =$fp;

	}
	
	print outfile "\n\n********************\nConsensus Peptide Summary Stats\n\n";
	print outfile "FDR Threshold\tTrue positives\tFalse Positives\tCalculate FDR\tSearch engine\n";
	
	
	for (my $k = 0; $k < @test_thresholds; $k++){
		my $threshold = $test_thresholds[$k];
		my $tp = $tp_results[$k];
		my $fp = $fp_results[$k];
		my $fdr = 0;
		if($tp + $fp != 0){
			$fdr = $fp / ($tp + $fp);
		}
			
		print outfile "$threshold\t$tp\t$fp\t$fdr\tCombined\n";
	}

	print outfile "\n\n";
	
	return (\@tp_results,\@fp_results,\%all_proteins);
}



sub resolveConsensusConflicts{

			
	print "Resolving conflicting peptide to protein matches for consensus proteins...\n";		
		
	for my $exp_title ( keys %exp_title_prot_cons ) {				
		my $expt_res = $exp_title_prot_cons{$exp_title};
		
		my %consensus_proteins = %{$expt_res->consensusProteins()};		
		my %resolved_proteins = ();	#put passing proteins in here	
		my %non_unique_proteins = ();	#put discarded proteins in here		
		
		my %pag_mappings = %{$exp_title_obj->pagMappings()};
		my @pags = @{$exp_title_obj->pags()};
		
		my %final_pags = ();			

		for my $pag_string (@pags){						
	
			my @accs = split(/,/,$pag_string);
			
			my @main_peps;
			my $main_acc= "";
			my %main_pep_idents;
			
			my $new_pag = eval { new ProteinAmbiguityGroup(); }  or die ($@);	
					
			
			foreach my $acc (@accs){
				
				if($consensus_proteins{$acc}){	#there will be no consensus protein entry for this accession if the peptides fall below a set cutoff
					my %prot_pep_idents = %{$consensus_proteins{$acc}};				
				
					if(!@main_peps){	#first ident
						$main_acc = $acc;
						for my $pl_rank (keys %prot_pep_idents){		
					
							my $ident = $prot_pep_idents{$pl_rank};						
							push(@main_peps, $ident->getPepSeq());
						}
						@main_peps = sort(@main_peps); 
						%main_pep_idents = %prot_pep_idents;
					}
					else{
						my @new_peps;
						for my $pl_rank (keys %prot_pep_idents){		
					
							my $ident = $prot_pep_idents{$pl_rank};						
							push(@new_peps, $ident->getPepSeq());
						}
						@new_peps = sort(@new_peps);
											
						my $pep_seq = join(',', @main_peps);
						my $new_pep_seq = join(',', @new_peps);
						
						if($pep_seq eq $new_pep_seq){						
							#print "Equal idents, main: $main_acc , same prot: $acc \n";
							
							$new_pag->addSameProteins($acc);
							$non_unique_proteins{$acc} = \%prot_pep_idents;
							
						}
						else{	# need to test if one contains the other
												
							if(@new_peps < @main_peps){
								my $is_subpep = 1;
								foreach my $npep (@new_peps){
									if(index($pep_seq,$npep) == -1){
										$is_subpep = 0;										
									}
								}
								
								if($is_subpep){
									#print "$main_acc resolves pep: $acc \n";
									$new_pag->addSubProteins($acc);
									$non_unique_proteins{$acc} = \%prot_pep_idents;									
								}
								else{
									#print "Conflict between $main_acc and $acc \n";
									$new_pag->addConflicts($acc);
									
									#Also create a conflict entry for other one
									my $tmp_pag = $final_pags{$acc};
									if(!$tmp_pag){
										$tmp_pag = eval { new ProteinAmbiguityGroup(); }  or die ($@);
									}
									else{
										$tmp_pag = $final_pags{$acc};
									}
									$tmp_pag->addConflicts($main_acc);
									$final_pags{$acc} = $tmp_pag;
									}
							}
							elsif(@new_peps == @main_peps){
								#print "both must have unique peptides, same length but not equal: $acc ($new_pep_seq) and $main_acc ($pep_seq)\n";
								#print "Conflict between $main_acc and $acc \n";
								$new_pag->addConflicts($acc);
																	#Also create a conflict entry for other one
								#Also create a conflict entry for other one
								my $tmp_pag = $final_pags{$acc};
								if(!$tmp_pag){
									$tmp_pag = eval { new ProteinAmbiguityGroup(); }  or die ($@);
								}
								else{
									$tmp_pag = $final_pags{$acc};
								}
								$tmp_pag->addConflicts($main_acc);
								$final_pags{$acc} = $tmp_pag;
							}
							else{
								my $is_subpep = 1;
								foreach my $opep (@main_peps){
									if(index($new_pep_seq,$opep) == -1){
										$is_subpep = 0;
										#print "$opep ($main_acc, $pep_seq) not contained within $new_pep_seq ($acc) \n";
										#print "Conflict between  $acc and $main_acc\n";
										$new_pag->addConflicts($main_acc);

										#Also create a conflict entry for other one
										my $tmp_pag = $final_pags{$main_acc};
										if(!$tmp_pag){
											$tmp_pag = eval { new ProteinAmbiguityGroup(); }  or die ($@);
										}
										else{
											$tmp_pag = $final_pags{$main_acc};
										}
										$tmp_pag->addConflicts($acc);
										$final_pags{$main_acc} = $tmp_pag;
																				
										$main_acc = $acc;
										@main_peps = @new_peps;
										%main_pep_idents = %prot_pep_idents;
									}
								}
								
								if($is_subpep){
									#print "$acc resolves pep: $main_acc \n";

									$non_unique_proteins{$main_acc} = \%main_pep_idents;
									$new_pag->addSubProteins($main_acc);
									$main_acc = $acc;
									@main_peps = @new_peps;
									%main_pep_idents = %prot_pep_idents;
								}
							}					
						}						
					}
				}
			}
						
			$final_pags{$main_acc} = $new_pag;
			#$resolved_proteins{$main_acc} = \%main_pep_idents;	#I don't think this is needed
		}		

		#$expt_res->resolvedProteins(\%resolved_proteins);
		$expt_res->nonUniqueProteins(\%non_unique_proteins);
		$expt_res->finalPags(\%final_pags);
	}	
	
}


#For consensus proteins 
sub assignProteinQValues{


	my ($o1_tp_at_01,$o1_tp_at_05,$o1_tp_at_10,$o1_tp_at_20);
	my ($o2_tp_at_01,$o2_tp_at_05,$o2_tp_at_10,$o2_tp_at_20);
	my ($o3_tp_at_01,$o3_tp_at_05,$o3_tp_at_10,$o3_tp_at_20);
	
	for my $exp_title ( keys %exp_title_prot_cons ) {

		my $expt_res = $exp_title_prot_cons{$exp_title};	
		my %consensus_proteins = %{$expt_res->consensusProteins()};	
		
		my %non_unique_proteins = %{$expt_res->nonUniqueProteins()};
		my %final_pags = %{$expt_res->finalPags()};

		my $num = keys %consensus_proteins;	

		my ($tmp2,$tmp_subProts) = assignProteins(\%consensus_proteins,\%non_unique_proteins);
	
		#my @ordered1 = @{$tmp1};	#test metric 1 first
		my @ordered2 = @{$tmp2};
		#my @ordered3 = @{$tmp3};
		
		my @sub_proteins = @{$tmp_subProts};
		
		#May want to create running totals of TP at q <= 0.01, 0.05, 0.1		
		
		#my ($temp1,$temp2,$temp3,$temp4) =  assignProtQ(\@ordered1,1,$exp_title);
		
		#$o1_tp_at_01 = $o1_tp_at_01 + $temp1;
		#$o1_tp_at_05 = $o1_tp_at_05 + $temp2;
		#$o1_tp_at_10 = $o1_tp_at_10 + $temp3;
		#$o1_tp_at_20 = $o1_tp_at_20 + $temp4;
		
		my ($temp1,$temp2,$temp3,$temp4) = assignProtQ(\@ordered2,2,$exp_title,\%final_pags);
		
		$o2_tp_at_01 = $o2_tp_at_01 + $temp1;
		$o2_tp_at_05 = $o2_tp_at_05 + $temp2;
		$o2_tp_at_10 = $o2_tp_at_10 + $temp3;
		$o2_tp_at_20 = $o2_tp_at_20 + $temp4;
		
		#($temp1,$temp2,$temp3,$temp4) = assignProtQ(\@ordered3,3,$exp_title);
		
		#$o3_tp_at_01 = $o3_tp_at_01 + $temp1;
		$o3_tp_at_05 = $o3_tp_at_05 + $temp2;
		$o3_tp_at_10 = $o3_tp_at_10 + $temp3;
		$o3_tp_at_20 = $o3_tp_at_20 + $temp4;	

		print prot_out "\n*****************\n\n";		
		
		print prot_out "\nProteins with no unique peptide assignments from $exp_title(referenced above):\n\n";
		print prot_out "Acc\tProteinScore\tPeptides\n";
		
		foreach my $prot (@sub_proteins){

			my $prot_acc = $prot->protAcc();		
			my %pep_idents = %{$prot->idents()};
			my $peps;
			
			for my $pl_rank (keys %pep_idents){
				my $ident = $pep_idents{$pl_rank};
				my $pep = $ident->getPepSeq();
				my $q = $ident->weightedFDR();
				$peps .= "$pep ($q), ";				
			}
			
			my $metric = $prot->aggNRQ();
			print prot_out "$prot_acc\t$metric\t$peps\n";
		}				
		
		print prot_out "\n\n";
	}

	print prot_out "**************************\n**************************\n**************************\n";
	print prot_out "\n\nSummary of protein true positives\n";
	print prot_out "Order scheme\tTP at 1%\tTP at 5%\tTP at 10%\tTP at 20%\n";
	
	print prot_out "Aggregate NR FDR Score\t$o2_tp_at_01\t$o2_tp_at_05\t$o2_tp_at_10\t$o2_tp_at_20\n";
}



sub assignProteins{

	my ($tmp1,$tmp2) = @_;
	my %proteins = %{$tmp1};
	my %sub_proteins = %{$tmp2};	#these are proteins that do not have unique peptide assignments and are excluded from the calculation
	
	my $num = keys %proteins;	
	#print "Assigning proteins for " . $num ." values ...\n";
	
	my (@cons_proteins, @sub_proteins);
	
	#First create consensus protein entries 
	for my $prot_acc (keys %proteins){	

		my $cons_prot = eval { new ConsensusProtein(); }  or die ($@);			
		my %prot_pep_idents = %{$proteins{$prot_acc}};				
		
		my $agg_redundat_qval =1;	
		my %nr_pep_seqs = ();
		my $is_decoy = 0;
		
		my $tmp = keys %prot_pep_idents;
	
		for my $pl_rank (keys %prot_pep_idents){
			
			my $ident = $prot_pep_idents{$pl_rank};	
			
			if(index($prot_acc,$rev_string) != -1){
				$is_decoy = 1;
			}
			
			my $pepSeq = $ident->getPepSeq();
			my $q = $ident->weightedFDR();
			
			$agg_redundat_qval = $agg_redundat_qval * $q;
			
			if($nr_pep_seqs{$pepSeq}){
				my $test_ident = $nr_pep_seqs{$pepSeq};
				my $test_q = $test_ident->weightedFDR();
				
				if($q < $test_q){	#lower q-value for this ident
					$nr_pep_seqs{$pepSeq} = $ident;
				}				
			}
			else{
				$nr_pep_seqs{$pepSeq} = $ident;
			}			
		}
		
		$cons_prot->aggRedundantQ($agg_redundat_qval);
		$cons_prot->nonRedPeps(\%nr_pep_seqs);
		$cons_prot->isDecoy($is_decoy);
		$cons_prot->protAcc($prot_acc);
		$cons_prot->idents(\%prot_pep_idents);
		
		my $count = 0;
		my $agg_nr_q = 1;
		my $geo_mean_q = 1;
		

		#calculate aggregate NR q and geomean q
		for my $pep_seq (keys %nr_pep_seqs){
			my $ident = $nr_pep_seqs{$pep_seq};	
			my $q = $ident->weightedFDR();			
			
			#$nr_peps .= "$pep_seq ($q),";
			
			$agg_nr_q = $agg_nr_q * $q;
			$count++;
		}		
		
		$geo_mean_q = $agg_nr_q ** (1/$count);		
		
		$cons_prot->aggNRQ($agg_nr_q);
		$cons_prot->geoMeanQ($geo_mean_q);	
		
		#print "$prot_acc\t$agg_redundat_qval\t$agg_nr_q\t$geo_mean_q\tnr:$nr_peps\n";
		if(!$sub_proteins{$prot_acc}){
			push(@cons_proteins,$cons_prot);
		}
		else{
			push(@sub_proteins,$cons_prot);
		}
		
	}

	#print "Total consensus proteins: " . @cons_proteins . "\n";
	
	#my @ordered1 = @{orderProteins(\@cons_proteins,1)};
	my @ordered2 = @{orderProteins(\@cons_proteins,2)};
	#my @ordered3 = @{orderProteins(\@cons_proteins,3)};
	@sub_proteins = @{orderProteins(\@sub_proteins,2)};
	
	return (\@ordered2,\@sub_proteins);
	#return (\@ordered1,\@ordered2,\@ordered3);

}

#orders proteins according to one of three criteria
	#Various options for how to order proteins
	
	#Option 1 - order by aggregate redundant peptide q-value
	#Option 2 - order by aggregate non-redundant peptide q-value
	#Option 3 - geometric mean of NR peptide q-value
	
sub orderProteins{

	my ($tmp1,$option) = @_;
	#my %proteins = %{$tmp1};
	my @ordered = @{$tmp1}; 

	#print "Order proteins for " . @ordered . " values\n";
		
	my $sorted = 0;
	
	while(!$sorted){		
		$sorted = 1;
		
		for (my $i=0; $i<(@ordered-1); $i++){
					
			my $first_ident =  $ordered[$i];
			my $second_ident = $ordered[$i+1];
			
			my ($first_metric,$second_metric);
			
			if($option==1){
				$first_metric = $first_ident->aggRedundantQ();
				$second_metric = $second_ident->aggRedundantQ();
			}
			elsif($option==2){
				$first_metric = $first_ident->aggNRQ();
				$second_metric = $second_ident->aggNRQ();
			}
			elsif($option==3){
				$first_metric = $first_ident->geoMeanQ();
				$second_metric = $second_ident->geoMeanQ();
			}			
		
			if($second_metric < $first_metric){
				$ordered[$i] = $second_ident;
				$ordered[$i+1] = $first_ident;
				$sorted = 0;
			}		
		}
	}	
	
	return \@ordered;
}


sub assignProtQ{
	
	my ($tmp,$option,$exp_title,$tmp1) = @_;
	my @ordered = @{$tmp};
	my %final_pags = %{$tmp1};
	
	my $all_idents = 0;
	my $decoy_count = 0;
	
	#print "Assign protein Q-values for " . @ordered ." values according to aggregate FDR Score for NR peptides \n";
	
	print prot_out "********************************\n";
	print prot_out "********************************\n";		
	print prot_out "Proteins for assay $exp_title\n\n";	

	
	for(my $i=0; $i<@ordered; $i++){
	
		$all_idents++;
		my $ident =  $ordered[$i];
		my $is_decoy = $ident->isDecoy();
		
		if($is_decoy){
			$decoy_count++;
		}	
	
		my $all_targets = $all_idents - $decoy_count;

		my $fp = $decoy_count / $decoy_ratio;
		my $tp = $all_targets - $fp;		

		if($tp < 0){
			$tp = 0;
		}
		my $fdr = $fp / ( $tp + $fp);
		
		if($option == 1){
			$ident->protQ1($fdr);
		}
		elsif($option == 2){
			$ident->protQ2($fdr);
		}
		elsif($option == 3){
			$ident->protQ3($fdr);
		}
		else{
			print "Error, no option given\n";
			exit;
		}		
	}	
	

	my $current_lowest = 1;
	
	#Second run through converts FDR to q-values by storing current lowest at all stages	
	for(my $i=@ordered-1; $i>=0; $i--){	

		$all_idents++;
		my $ident =  $ordered[$i];
		my $identFDR; 
		
		if($option == 1){
			$identFDR = $ident->protQ1();
		}
		elsif($option == 2){
			$identFDR = $ident->protQ2();
		}
		elsif($option == 3){
			$identFDR = $ident->protQ3();
		}		
		
		if($identFDR < $current_lowest){
			$current_lowest = $identFDR;
		}
		else{
			if($option == 1){
				$ident->protQ1($current_lowest);
			}
			elsif($option == 2){
				$ident->protQ2($current_lowest);
			}
			elsif($option == 3){
				$ident->protQ3($current_lowest);
			}
		}		
	}	
	
	my $prev_q = 0;
	my $tp_count = 0;
	
	my ($tp_at_01FDR,$tp_at_05FDR,$tp_at_10FDR,$tp_at_20FDR);

	$decoy_count = 0;
	my $tp = 0;
	my $target_count = 0;
	
	
	print prot_out "Acc\tDescription\tTP\tT/D\tProteinScore\tFDR\tSame peps\tSubsequence\tConflict\tPeptides\n";
	
	for(my $i=0; $i<@ordered; $i++){
		my $ident =  $ordered[$i];
		my $is_decoy = $ident->isDecoy();
		
		my $identFDR;
		my $metric;
		if($option == 1){
			$identFDR = $ident->protQ1();
			$metric = $ident->aggRedundantQ();
		}
		elsif($option == 2){
			$identFDR = $ident->protQ2();
			$metric = $ident->aggNRQ();
		}
		elsif($option == 3){
			$identFDR = $ident->protQ3();
			$metric = $ident->geoMeanQ();			
		}		
		
		my $prot_acc = $ident->protAcc();		
		my %idents = %{$ident->idents()};		
		my $peps = "";
		
		for my $pl_rank (keys %idents){
			my $ident = $idents{$pl_rank};
			my $pep = $ident->getPepSeq();
			my $q = $ident->weightedFDR();
			$peps .= "$pep ($q) [";
			
			my $se_ident;
			my $mods = "";
			if($ident->mascotPeptide()){
				$se_ident = $ident->mascotPeptide();
				$mods .= ",m";
				my $tmp_mods = $se_ident->mods();
				if($tmp_mods){
					$mods .= ":$tmp_mods";
				}
			}
			if($ident->omssaPeptide()){
				$se_ident = $ident->omssaPeptide();
				$mods .= ",o";
				my $tmp_mods = $se_ident->mods();
				if($tmp_mods){
					$mods .= ":$tmp_mods";
				}
			}
			if($ident->tandemPeptide()){
				$se_ident = $ident->tandemPeptide();
				$mods .= ",t";
				my $tmp_mods = $se_ident->mods();
				if($tmp_mods){				
					$mods .= ":$tmp_mods";
				}
			}
			
			$mods =~ s/\,//;			
			$peps .= $mods."],";
		}	

		my ($same_prots,$sub_prots,$conflicts);
		if($final_pags{$prot_acc}){
			my $pag = $final_pags{$prot_acc};
			my @sameProts =  $pag->getSameProteins();
			my @subProts =$pag->getSubProteins();
			my @conflicts =$pag->getConflicts();
			
			foreach my $prot (@sameProts){
				$same_prots .= $prot . ",";
			}
			foreach my $prot (@subProts){
				$sub_prots .= $prot . ",";
			}
			foreach my $prot (@conflicts){
				$conflicts .= $prot . ",";
			}
			chop($same_prots);
			chop($sub_prots);
			chop($conflicts);
		}
		
		
		my $desc_line = "";	
		
		if(%db_seq){			
			$desc_line = $db_seq{$prot_acc};
		}
		
		if($is_decoy){
			$decoy_count++;
			my $fp = $decoy_count / $decoy_ratio;			
			$tp = $target_count - $fp;
			
			#print "$prot_acc\t$tp\tDecoy\t$metric\t$identFDR\t$same_prots\t$sub_prots\t$conflicts\t$peps\n";
			if($prot_out && $option == 2){
				print prot_out "$prot_acc\t$desc_line\t$tp\tDecoy\t$metric\t$identFDR\t$same_prots\t$sub_prots\t$conflicts\t$peps\n";			
			}
		}
		else{
			$target_count++;
			my $fp = $decoy_count / $decoy_ratio;
			$tp = $target_count - $fp;
			
			#print "$prot_acc\t$tp\tTarget\t$metric\t$identFDR\t$same_prots\t$sub_prots\t$conflicts\t$peps\n";
			if($prot_out && $option == 2){
				print prot_out "$prot_acc\t$desc_line\t$tp\tTarget\t$metric\t$identFDR\t$same_prots\t$sub_prots\t$conflicts\t$peps\n";		
			}
		}		
		
		if($tp < 0){
			$tp = 0;
		}	
		
		if($identFDR > 0.01 && !$tp_at_01FDR){
			$tp_at_01FDR = $tp;
		}
		if($identFDR > 0.05 && !$tp_at_05FDR){
			$tp_at_05FDR= $tp;
		}
		if($identFDR > 0.1 && !$tp_at_10FDR){
			$tp_at_10FDR= $tp;
		}
		if($identFDR > 0.2 && !$tp_at_20FDR){
			$tp_at_20FDR= $tp;
		}
		
		#Uncomment this part to stop printing output
		#if($identFDR > 0.25){
		#	$i = @ordered;
		#}		
	}
	
	
	print prot_out "\n\n\n*************************************\n";
	print prot_out "Protein stats for Experiment: $exp_title\n";
	print prot_out "TP at 1%FDR\tTP at 5%FDR\tTP at 10%FDR\n";
	print prot_out "$tp_at_01FDR\t$tp_at_05FDR\t$tp_at_10FDR\n";
	
	print prot_out "\n\n";	
	
	return ($tp_at_01FDR,$tp_at_05FDR,$tp_at_10FDR,$tp_at_20FDR);
}


sub trim{
	my $string = shift;	
	
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	
	return $string;
}


