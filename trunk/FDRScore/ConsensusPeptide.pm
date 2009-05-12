
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




package ConsensusPeptide;
use strict;

#constructor
sub new {
    my ($class) = @_;
	
	#This is a hash of key, value pairs, it may not be possible to put an array in a hash!
	
    my $self = {
	    
		#what kind of data structures are needed for a protein consensus...
			# also multi-dimensional...
		_omssa_peptide => undef,		
		_mascot_peptide => undef,
		_tandem_peptide => undef,
		_spectrum_ID => undef,
		_pep_seq => undef,
		_prot_accs=>[],
		_geomean_fdr => undef,
		_qvalue => undef,
		_win_fdr => undef,		#FDR calculated for a window, added Jan 2009
		_local_fdr => undef,
		_weighted_fdr => undef,
		_mean_rank => undef,
		_order_rank => undef,
		_is_decoy => undef
		
		#_peaklist => undef	#the peaklist is the key used to find this...
		
    };
	bless $self, $class;
    return $self;	
}

sub meanRank{
    my ( $self, $rank ) = @_;
    $self->{_mean_rank} = $rank if defined($rank);
    return $self->{_mean_rank};
}

sub orderRank{
    my ( $self, $rank ) = @_;
    $self->{_order_rank} = $rank if defined($rank);
    return $self->{_order_rank};
}


sub isDecoy{
    my ( $self, $decoy ) = @_;
    $self->{_is_decoy} = $decoy if defined($decoy);
    return $self->{_is_decoy};
}

#added Apr 2009
sub addStart{
my ( $self, $start ) = @_;
$self->{_start} = $start if defined($start);
return $self->{_start};
}


sub localFDR {
    my ( $self, $q ) = @_;
    $self->{_local_fdr} = $q if defined($q);
    return $self->{_local_fdr};
}


sub assignPeptide{
    my ( $self, $pep ) = @_;
	
	my $se = $pep->searchEngine();
	my $pepSeq = $pep->pepSeq();
	
	if($se eq "omssa"){
		$self->{_omssa_peptide} = $pep if defined($pep);
	}
	elsif($se eq "mascot"){
		$self->{_mascot_peptide} = $pep if defined($pep);
	}
	elsif($se eq "X!Tandem"){
		$self->{_tandem_peptide} = $pep if defined($pep);
	}
	else{
		print "No search engine, exit";
		exit;
	}
	$self->{_pep_seq} = $pepSeq if defined($pepSeq);
}

sub getPepSeq{
	my ( $self) = @_;
	return $self->{_pep_seq};
}

sub spectrumID{
	   my ( $self, $id ) = @_;
    $self->{_spectrum_ID} = $id if defined($id);
    return $self->{_spectrum_ID};
}


sub getProtAccs{
	my ($self) = @_;
	
	
	my $pep;
	
	if($self->{_mascot_peptide}){
		$pep = $self->{_mascot_peptide};
	}
	elsif($self->{_tandem_peptide}){
		$pep = $self->{_tandem_peptide};
	}
	elsif($self->{_omssa_peptide}){
		$pep = $self->{_omssa_peptide};
	}
	else{
		print "Error, no search peptide stored \n";
		exit;
	}
	
	my @prot_accs = $pep->getProtAccs();
	
	return @prot_accs;
}



sub omssaPeptide{
    my ( $self, $pep ) = @_;
    $self->{_omssa_peptide} = $pep if defined($pep);
    return $self->{_omssa_peptide};
}

sub tandemPeptide{
    my ( $self, $pep ) = @_;
    $self->{_tandem_peptide} = $pep if defined($pep);
    return $self->{_tandem_peptide};
}

sub mascotPeptide{
    my ( $self, $pep ) = @_;
    $self->{_mascot_peptide} = $pep if defined($pep);
    return $self->{_mascot_peptide};
}

sub peaklist {
    my ( $self, $pl ) = @_;
    $self->{_peaklist} = $pl if defined($pl);
    return $self->{_peaklist};
}

sub geoMeanFDR{
    my ( $self, $score ) = @_;
    $self->{_geomean_fdr} = $score if defined($score);
    return $self->{_geomean_fdr};
}

sub qValue {
    my ( $self, $q ) = @_;
    $self->{_qvalue} = $q if defined($q);
    return $self->{_qvalue};
}

sub winFDR {
    my ( $self, $wf ) = @_;
    $self->{_win_fdr} = $wf if defined($wf);
    return $self->{_win_fdr};
}


sub weightedFDR{
    my ( $self, $score ) = @_;
    $self->{_weighted_fdr} = $score if defined($score);
    return $self->{_weighted_fdr};
}


1;
