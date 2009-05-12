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



package PeptideIdent;
use strict;

#constructor
sub new {
    my ($class) = @_;
	
	#This is a hash of key, value pairs, it may not be possible to put an array in a hash!
	
    my $self = {
	    _prot_accs => [],
        _pep_seq => undef,
		_search_engine => undef,
		_peaklist => undef,
		_rank => undef,
		_order_rank => undef,
		_qvalue => undef,
		_fdr_score => undef,
		_win_fdr => undef,		#FDR calculated for a window, added Jan 2009
		_evalue => undef,
		_score => undef,
		_local_fdr => undef,
		#_exp_acc => undef,			#can't see a reason to use this,...
		_exp_title => undef,
		_is_redundant => undef,	#set if there is a stronger peptide hit within the same protein
		_mods => undef,		# Stores the mod string for each peptide
		_is_decoy => undef
    };
	bless $self, $class;
    return $self;	
}


sub isDecoy{
    my ( $self, $decoy ) = @_;
    $self->{_is_decoy} = $decoy if defined($decoy);
    return $self->{_is_decoy};
}


sub addProtAcc{
    my ( $self, $protAcc ) = @_;
	push ( @{ $self->{_prot_accs} } ,$protAcc);		
}

sub getProtAccs{
	my ($self) = @_;
	return @{$self->{_prot_accs}};
}


sub expAcc{
    my ( $self, $exp_acc ) = @_;
    $self->{_exp_acc} = $exp_acc if defined($exp_acc);
    return $self->{_exp_acc};
}

#added Apr 2009 JAS
sub addStart{
my ( $self, $start ) = @_;
$self->{_start} = $start if defined($start);
return $self->{_start};
}



sub expTitle{
    my ( $self, $exp_title ) = @_;
    $self->{_exp_title} = $exp_title if defined($exp_title);
    return $self->{_exp_title};
}


sub pepSeq {
    my ( $self, $pep_seq ) = @_;
    $self->{_pep_seq} = $pep_seq if defined($pep_seq);
    return $self->{_pep_seq};
}

sub score {
    my ( $self, $score ) = @_;
    $self->{_score} = $score if defined($score);
    return $self->{_score};
}

sub localFDR {
    my ( $self, $q ) = @_;
    $self->{_local_fdr} = $q if defined($q);
    return $self->{_local_fdr};
}

sub qValue {
    my ( $self, $q ) = @_;
    $self->{_qvalue} = $q if defined($q);
    return $self->{_qvalue};
}

sub fdrScore {
    my ( $self, $score ) = @_;
    $self->{_fdr_score} = $score if defined($score);
    return $self->{_fdr_score};
}


sub searchEngine {
    my ( $self, $se ) = @_;
    $self->{_search_engine} = $se if defined($se);
    return $self->{_search_engine};
}

sub peaklist {
    my ( $self, $pl ) = @_;
    $self->{_peaklist} = $pl if defined($pl);
    return $self->{_peaklist};
}
sub rank {
    my ( $self, $rank ) = @_;
    $self->{_rank} = $rank if defined($rank);
    return $self->{_rank};
}
sub orderRank {
    my ( $self, $rank ) = @_;
    $self->{_order_rank} = $rank if defined($rank);
    return $self->{_order_rank};
}
sub evalue {
    my ( $self, $eval ) = @_;
    $self->{_evalue} = $eval if defined($eval);
    return $self->{_evalue};
}
sub winFDR {
    my ( $self, $wf ) = @_;
    $self->{_win_fdr} = $wf if defined($wf);
    return $self->{_win_fdr};
}

sub isRedundant {
    my ( $self, $eval ) = @_;
    $self->{_is_redundant} = $eval if defined($eval);
    return $self->{_is_redundant};
}

sub mods {
    my ( $self, $mods ) = @_;
    $self->{_mods} = $mods if defined($mods);
    return $self->{_mods};
}



1;
