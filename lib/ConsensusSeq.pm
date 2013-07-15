package ConsensusSeq;

use warnings;
use strict;
use Carp;

#################### main pod documentation begin ###################

=head1 NAME

BoundaryFinder - Exon/Intron Boundary Finder

=head1 SYNOPSIS

  use BoundaryFinder;

=head1 DESCRIPTION

The ConsensusSeq object that gives methods for accessing and scoring using a consensus sequence.

=head1 USAGE

#################### main pod documentation end ###################

=head1 METHODS

=head2 new

 Usage     : $consensus = ConsensusSeq->new

 Purpose   : Gets a new, empty consensus sequence

 Returns   : ConsensusSeq

 Args      : 

 Throws    : 

 Comment   : 

See Also   : ConsensusSeqIO

=cut

#
# new
#
sub new {
    my $class = shift;    

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    # Defaults
    my $self = {
	-title => "",
	-missing_index_score => 0,
	-boundary => 0
    };

    bless ($self, $class);
    
    # Check Params
    # No Params to check

    # Required default values
    $self->{-length} = 0;
    $self->{-consensus_seq} = [];
    $self->{-seq_keys} = [];

    return $self;
}


=head2 add_element

 Usage     : $consensus->add_element('a' => 0.25, 'c' => 0.25, 't' => 0.25, 'g' => 0.25);

 Purpose   : Adds an element to the end of the consensus sequence

 Returns   : 

 Argument  : Turns named args into a hash in which determines score for the element.

 Throws    : 

 Comment   : Use single characte keys only. Multi characters will not be evauluated.

 See Also   : 

=cut

sub add_element {
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;

    my %params = @_;

    my %vals = map { lc $_ => $params{$_} } keys %params;
    push @{$self->{-consensus_seq}}, \%vals;

    # Add keys to seq_keys
    foreach my $k (keys %vals) {
	if (not grep /^$k$/, @{$self->{-seq_keys}}) {
	    push @{ $self->{-seq_keys} }, $k;
	}
    }

    $self->{-length}++;
}

=head2 boundary

 Usage     : my $boundary = $consensus->boundary;
           : $consensus->missing_index_score(-boundary => 7);

 Purpose   : gets/sets the boundary index

 Returns   : The current boundary index

 Argument  : [-boundary] The a new boundary index

 Throws    : 

 Comment   : The boundary index is the index from 0 to self->length-1 that
           : indicates the boundary of the exon and intron in a sequence.
           : This is important because the boundaries is usually in the middle
           : of the ConsensusSeq.

 See Also   : 

=cut

sub boundary {
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    $self->{-boundary} = $params{-boundary} if defined ($params{-boundary});
    return $self->{-boundary};
}


=head2 length

 Usage     : $consensus->length

 Purpose   : Gets the lenth of the consensus sequence

 Returns   : Length of the consensus sequence

 Argument  : 

 Throws    : 

 Comment   : 

See Also   : 

=cut

sub length {
    my ($self) = shift;
    return $self->{-length};
}

=head2 missing_index_score

 Usage     : my $mis = $consensus->missing_index_score;
           : $consensus->missing_index_score(-score => 0.01);
 Purpose   : gets/sets the lenth of the consensus sequence
 Returns   : The current missing index score
 Argument  : The new missing index score
 Throws    : 
 Comment   : The missing index score is the score when an indext is not found. Defualt: 0
           : E.g. If you did not set a score for an "n" in a sequence element, the score for the "n" would be the missing index score. 

 See Also   : 

=cut

sub missing_index_score {
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    $self->{-missing_index_score} = $params{-score} if defined ($params{-score});
    return $self->{-missing_index_score};
}

=head2 score_element

 Usage     : my $score = $consensus->score_element(nucl => "a", index => $index); 

 Purpose   : Gets the score for a NT at a given index

 Returns   : The score of an NT with index

 Argument  : nucl :  A nucleotide
           : index : The index of the nt
 Throws    : 

 Comment   : Returns the missing index score if not found
           : Case insensitive

See Also   : 

=cut

sub score_element {
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    for my $required ( qw{ nucl index } ) {
	croak "Required parameter '$required' not passed to score_element"
            unless exists $params{$required}; 
    }

    # Check if the value exists before accessing it
    if ( exists ($self->{-consensus_seq}[$params{index}]->{$params{nucl}}) ){
	return $self->{-consensus_seq}[$params{index}]->{$params{nucl}};
    } 
    # Else Missing idx score
    return $self->missing_index_score;
}

=head2 score_seq

 Usage     : my $score = $consensus->score_seq(seq => $seqence);

 Purpose   : Gets the score for a NT sequence

 Returns   : The score of an NT with index

 Argument  : seq:  A nucleotide

 Comment   : Returns the score of given sequence. If the sequence are not the same length, it uses the missing_index_score for a sequence that is shorter that the consensus, and the first _length for sequences that are longer.

=cut

sub score_seq {
    my $self = shift;
    my $score;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    for my $required (qw{ seq }) {
	croak "Required parameter '$required' not passed to score_seq"
            unless exists $params{$required}; 
    }

    my $sequence = lc($params{seq});

    $score = $self->score_element(nucl => substr($sequence, 0, 1), index => 0);
    for (my $i = 1; $i < $self->length; $i++) {

	#
	# Verify that $i is not out of the $sequence bounds
	# 
	if ($i < CORE::length ($sequence)) {
	    $score *= $self->score_element(nucl => substr($sequence, $i, 1), index => $i);
	}
	else {
	    $score = 0;
	    last;
	}
    }

    return $score;
}

=head2 title

 Usage     : my $title = $consensus->title;
           : $consensus->title(title => "Some title");
 Purpose   : Simple title for sequence to help differeniate between sequences eg, 5' vs 3'

 Returns   : The title of the consensus sequence

 Argument  : title : The new title

 Throws    : 
 
 Comment   : 

 See Also   : 

=cut

sub title {
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;
    
    $self->{-title} = $params{title} if defined ($params{title});
    return $self->{-title};
}

=head2 to_string

 Usage     : print $conSeq->toString

 Purpose   : returns the consensus sequence as a readable string.

 Returns   : string

 Argument  : 

 Throws    : 

 Comment   : 

See Also   : 

=cut

sub to_string {
    my $self = shift;
    my $string;

    $string = '>' . $self->title . ',' . $self->missing_index_score . "\n";
    $string .= "keys: " . join (",", @{ $self->{-seq_keys} }) . "\n";

    for (my $i = 0; $i < $self->length; $i++) {
	my @tmp_vals;
	foreach (@{ $self->{-seq_keys} }) {
	    push @tmp_vals, $self->{-consensus_seq}[$i]->{$_};
	}
	$string .= join(",", @tmp_vals) . "\n";
    }

    return $string;

}

=head1 BUGS


=head1 SUPPORT


=head1 AUTHOR

    Nathan Elmore
    nate@elmoren.com
    http://www.elmoren.com

=head1 COPYRIGHT

=head1 SEE ALSO

=cut

1;
