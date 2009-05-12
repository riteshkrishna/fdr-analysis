
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



package ConsensusProtein;
use strict;

#constructor
sub new {
    my ($class) = @_;
	
	#This is a hash of key, value pairs, it may not be possible to put an array in a hash!
	
    my $self = {
	    
		#what kind of data structures are needed for a protein consensus...
			# also multi-dimensional...
		_idents => undef,
		_prot_acc => undef,			
		_aggregate_red_q => undef,
		_aggregate_nr_q => undef,
		_nr_peps => undef,
		_prot_q1 => undef,
		_prot_q2 => undef,
		_prot_q3 => undef,		
		_geomean_nr_q => undef,
		_is_decoy => undef

		
		#_peaklist => undef	#the peaklist is the key used to find this...
		
    };
	bless $self, $class;
    return $self;	
}

sub isDecoy{
    my ( $self, $decoy ) = @_;
    $self->{_is_decoy} = $decoy if defined($decoy);
    return $self->{_is_decoy};
}

sub protAcc{
    my ( $self, $acc ) = @_;
    $self->{_prot_acc} = $acc if defined($acc);
    return $self->{_prot_acc};
}



sub idents{
    my ( $self, $idents ) = @_;
    $self->{_idents} = $idents if defined($idents);
    return $self->{_idents};
}

sub aggRedundantQ{
    my ( $self, $agg_q ) = @_;
    $self->{_aggregate_red_q} = $agg_q if defined($agg_q);
    return $self->{_aggregate_red_q};
}

sub aggNRQ{
    my ( $self, $agg_q ) = @_;
    $self->{_aggregate_nr_q} = $agg_q if defined($agg_q);
    return $self->{_aggregate_nr_q};
}

sub geoMeanQ{
    my ( $self, $agg_q ) = @_;
    $self->{_geomean_nr_q} = $agg_q if defined($agg_q);
    return $self->{_geomean_nr_q};
}


sub nonRedPeps{
    my ( $self, $peps ) = @_;
    $self->{_nr_peps} = $peps if defined($peps);
    return $self->{_nr_peps};
}

sub protQ1{
   my ( $self, $q ) = @_;
    $self->{_prot_q1} = $q if defined($q);
    return $self->{_prot_q1};
}

sub protQ2{
   my ( $self, $q ) = @_;
    $self->{_prot_q2} = $q if defined($q);
    return $self->{_prot_q2};
}

sub protQ3{
   my ( $self, $q ) = @_;
    $self->{_prot_q3} = $q if defined($q);
    return $self->{_prot_q3};
}




1;