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



package ExperimentTitle;
use strict;

#constructor
sub new {
    my ($class) = @_;
	
	#This is a hash of key, value pairs, it may not be possible to put an array in a hash!
	
    my $self = {
        _exp_title => undef,
		_mascot_proteins => undef,
        _mascot_peaklists => undef,
		_mascot_conflicts => undef,
		_mascot_resolved => undef,
		_tandem_resolved => undef,
		_omssa_resolved => undef,
		_tandem_conflicts => undef,
		_omsss_conflicts => undef,
		_omssa_proteins => undef,
		_omssa_peaklists => undef,		
		_tandem_proteins => undef,
		_tandem_peaklists => undef,		
		_consensus_proteins => undef,
		_consensus_proteins_nr => undef,
		_resolved_proteins_nr => undef,
		_non_unique_proteins => undef, #Added by ARJ 23 Jan to handle the proteins discarded as not having a unique peptide
		
		_pags => undef,				#Added by ARJ 26 Jan to handle protein ambiguity groups
		_pag_Mappings => undef,				#Added by ARJ 26 Jan to handle protein ambiguity groups
		_final_pags => undef,				#Added by ARJ 26 Jan to handle protein ambiguity groups
		
		_resolving_proteins_nr => undef,
		_all_unique_hits => undef, 
		_consensus_peaklists => [],

		_mascot_idents => undef,
		_mascot_ordered_idents => undef,
		_omssa_idents => undef,
		_omssa_ordered_idents => undef,
		_tandem_idents => undef,
		_tandem_ordered_idents => undef

    };
	bless $self, $class;
    return $self;	
}

sub expTitle {
    my ( $self, $exp_title ) = @_;
    $self->{_exp_title} = $exp_title if defined($exp_title);
    return $self->{_exp_title};
}

sub pags {
    my ( $self, $pags ) = @_;
    $self->{_pags} = $pags if defined($pags);
    return $self->{_pags};
}

sub finalPags {
    my ( $self, $pags ) = @_;
    $self->{_final_pags} = $pags if defined($pags);
    return $self->{_final_pags};
}


sub pagMappings {
    my ( $self, $pags ) = @_;
    $self->{_pag_Mappings} = $pags if defined($pags);
    return $self->{_pag_Mappings};
}


sub omssaProteins {
	my ( $self, $prots ) = @_;
    $self->{_omssa_proteins} = $prots if defined($prots);
    return $self->{_omssa_proteins};	
}
sub mascotProteins {
	my ( $self, $prots ) = @_;
    $self->{_mascot_proteins} = $prots if defined($prots);
    return $self->{_mascot_proteins};	
}

sub tandemProteins {
	my ( $self, $prots ) = @_;
    $self->{_tandem_proteins} = $prots if defined($prots);
    return $self->{_tandem_proteins};	
}


sub omssaIdents {
	my ( $self, $idents ) = @_;
    $self->{_omssa_idents} = $idents if defined($idents);
    return $self->{_omssa_idents};	
}
sub mascotIdents {
	my ( $self, $idents ) = @_;
    $self->{_mascot_idents} = $idents if defined($idents);
    return $self->{_mascot_idents};	
}
sub tandemIdents {
	my ( $self, $idents ) = @_;
    $self->{_tandem_idents} = $idents if defined($idents);
    return $self->{_tandem_idents};	
}

sub omssaOrderedIdents {
	my ( $self, $idents ) = @_;
    $self->{_omssa_ordered_idents} = $idents if defined($idents);
    return $self->{_omssa_ordered_idents};	
}
sub mascotOrderedIdents {
	my ( $self, $idents ) = @_;
    $self->{_mascot_ordered_idents} = $idents if defined($idents);
    return $self->{_mascot_ordered_idents};	
}
sub tandemOrderedIdents {
	my ( $self, $idents ) = @_;
    $self->{_tandem_ordered_idents} = $idents if defined($idents);
    return $self->{_tandem_ordered_idents};	
}


sub consensusProteins {
	my ( $self, $prots ) = @_;
    $self->{_consensus_proteins} = $prots if defined($prots);
    return $self->{_consensus_proteins};	
}

sub consensusProteins_nr {
	my ( $self, $prots ) = @_;
    $self->{_consensus_proteins_nr} = $prots if defined($prots);
    return $self->{_consensus_proteins_nr};	
}

sub resolvedProteins {
	my ( $self, $prots ) = @_;
    $self->{_resolved_proteins_nr} = $prots if defined($prots);
    return $self->{_resolved_proteins_nr};	
}

sub nonUniqueProteins {
	my ( $self, $prots ) = @_;
    $self->{_non_unique_proteins} = $prots if defined($prots);
    return $self->{_non_unique_proteins};	
}


sub resolvingProteins{
	my ( $self, $prots ) = @_;
    $self->{_resolving_proteins_nr} = $prots if defined($prots);
    return $self->{_resolving_proteins_nr};	
}


sub omssaPeaklists {
	my ( $self, $prots ) = @_;
    $self->{_omssa_peaklists} = $prots if defined($prots);
    return $self->{_omssa_peaklists};	
}
sub mascotPeaklists{
	my ( $self, $prots ) = @_;
    $self->{_mascot_peaklists} = $prots if defined($prots);
    return $self->{_mascot_peaklists};	
}
sub tandemPeaklists {
	my ( $self, $prots ) = @_;
    $self->{_tandem_peaklists} = $prots if defined($prots);
    return $self->{_tandem_peaklists};	
}

sub mascotConflicts {
	my ( $self, $prots ) = @_;
    $self->{_mascot_conflicts} = $prots if defined($prots);
    return $self->{_mascot_conflicts};	
}

sub tandemConflicts {
	my ( $self, $prots ) = @_;
    $self->{_tandem_conflicts} = $prots if defined($prots);
    return $self->{_tandem_conflicts};	
}

sub omssaConflicts {
	my ( $self, $prots ) = @_;
    $self->{_omssa_conflicts} = $prots if defined($prots);
    return $self->{_omssa_conflicts};	
}

sub mascotResolved {
	my ( $self, $prots ) = @_;
    $self->{_mascot_resolved} = $prots if defined($prots);
    return $self->{_mascot_resolved};	
}

sub tandemResolved {
	my ( $self, $prots ) = @_;
    $self->{_tandem_resolved} = $prots if defined($prots);
    return $self->{_tandem_resolved};	
}

sub omssaResolved {
	my ( $self, $prots ) = @_;
    $self->{_omssa_resolved} = $prots if defined($prots);
    return $self->{_omssa_resolved};	
}

sub allUniqueHits{
	my ( $self, $prots ) = @_;
    $self->{_all_unique_hits} = $prots if defined($prots);
    return $self->{_all_unique_hits};
}

1;

