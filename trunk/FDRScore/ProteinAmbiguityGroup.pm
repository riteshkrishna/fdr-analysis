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



package ProteinAmbiguityGroup;
use strict;

#constructor
sub new {
    my ($class) = @_;
	
	#This is a hash of key, value pairs, it may not be possible to put an array in a hash!
	
    my $self = {
	    _same_prots => [],
		_conflict => [],
        _sub_prots => []

    };
	bless $self, $class;
    return $self;	
}


sub addSameProteins {
    my ( $self, $sameHits ) = @_;
	push ( @{ $self->{_same_prots} } ,$sameHits);	
}

sub addSubProteins {
    my ( $self, $subHits ) = @_;
	push ( @{ $self->{_sub_prots} } ,$subHits);
}

sub getSameProteins{
	my ($self) = @_;
	return @{$self->{_same_prots}};
}
sub getSubProteins{
	my ($self) = @_;
	return @{$self->{_sub_prots}};
}

sub addConflicts {
    my ( $self, $subHits ) = @_;
	push ( @{ $self->{_conflict} } ,$subHits);
}

sub getConflicts{
	my ($self) = @_;
	return @{$self->{_conflict}};
}

1;