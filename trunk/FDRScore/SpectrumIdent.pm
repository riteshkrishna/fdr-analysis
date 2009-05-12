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



package SpectrumIdent;
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
		_tandem_peptide => undef

		#_peaklist => undef	#the peaklist is the key used to find this...
		
    };
	bless $self, $class;
    return $self;	
}




sub assignPeptide{
    my ( $self, $pep ) = @_;
	
	my $se = $pep->searchEngine();
	
	if($se eq "omssa"){
		$self->{_omssa_peptide} = $pep if defined($pep);
	}
	elsif($se eq "mascot"){
		$self->{_mascot_peptide} = $pep if defined($pep);
	}
	elsif($se eq "X!Tandem"){
		$self->{_tandem_peptide} = $pep if defined($pep);
	} 

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



1;