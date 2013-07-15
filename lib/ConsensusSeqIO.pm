package ConsensusSeqIO;

use warnings;
use strict;
use Carp;

#################### main pod documentation begin ###################

=head1 NAME

BoundaryFinder - ConsensusSeqIO

=head1 SYNOPSIS

  use ConsensusSeq;
  use ConsensusSeqIO;

  my $consSeqIO = ConsensusSeqIO->new ( path => "path/to/consensus" );
  my $consensus = ConsensusSeqIO->next_seq

=head1 DESCRIPTION

This is an object that gives methods for reading and/or writing using a ConsensusSeq.

Since there is not a defined format to store a consensus sequence, we use a CSV format. Each
row contains values and represents the a position in the sequence, each column represents 
a NT. 

NOTE: Empty values in row default to the ConsensusSeq->missing_index_score. This is important
because an empty value may not be zero, and you may want 0 values.

Like the fasta format, each sequence starts with a description line starting with a '>'
character. It contains the title, and a missing index score. Seperated by commas.

Example consensus sequence:

> 5' Drosophila Mel., 0.0
a,t,c,g,n
0.1,0.4,0,0.5,0.1
0.1,0.5,0.1,0.5,
0.2,0.3,0,0.5,0.1

This is a consensus sequence with length 3. It has 5 values that it will score, any
other values in the search sequence will be the missing index score.

Note that in the second row, the sequence omits the last value. This will default to the
missing index score when there is an 'n' in the second position of the sequence.

#################### main pod documentation end ###################

=head1 METHODS

=head2 new
    Usage       : $consensusIO = ConsensusSeqIO->new(path => "path/to/file");
    Purpose     : Gets a ConsensusSeqIO object for reading/writed stored consensus sequences
    Returns     : ConsensusSeqIO
    Args        : path => "path/to/file"
                : clobber => 0/1 # clobber determines whether to open for append or write
=cut

sub new {
    my $class = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;
    
    my $self = {
	-path => "",
	-clobber => 0
    };
    
    bless ($self, $class);
    
    # Validate params by checking if they exist is self
    for my $attr (keys %params) {	
        croak "Invalid parameter '$attr' passed to '$class' constructor"
            unless defined ($self->{$attr});

	$self->{$attr} = $params{$attr};
    }
    
    if ($self->{-path} eq "") {
	carp "Warning: ConsensusSeqIO: Path not defined\n";
    }
    
    return $self;
}

=head2 next_seq
    Usage       : $consensusSeq = ConsensusSeqIO->next_seq
    Purpose     : Gets a ConsensusSeqIO object for reading/writed stored consensus sequences
    Returns     : ConsensusSeq object
    Args        : 
=cut

sub next_seq {
    my $self = shift;

    if (!$self->{_fh}) {
	open my $tmpFH, "<", $self->{-path} or 
	    croak "ConsensusSeqIO: Error, cannot open file: $!";
	$self->{_fh} = $tmpFH;
    }

# Read the next sequence from file
    my $entry = $self->_read_seq;
    return unless $entry;

# Try and create a ConsensusSeq
    my $consensusSeq = ConsensusSeq->new();

    #
    # The first line is the header.
    #
    # carp "ConsensusSeqIO:: Parsling line: '@{$entry}[0]'";
    my @header = split /,/, shift @{$entry};

    $header[0] =~ s/^>\s*//g;

    $consensusSeq->title(title => $header[0]) if $header[0];
    $consensusSeq->missing_index_score(-score => $header[1]) if $header[1];
    $consensusSeq->boundary(-boundary => $header[2]) if $header[2];

    #
    # The second line is the list of possible nucleotides
    #
    my @nucleotides = split /,/, shift @{$entry};
    my %seqHash = map {$_  => 0} @nucleotides;

    foreach (@{$entry}) {
	
        # carp "ConsensusSeqIO:: Parsling line: '$_'";
	
	my @vals = split /,/, $_;

	if (@vals != @nucleotides) {
	    croak "Invalid sequence format, skipping rest of sequence. There must be equal numbers of scores as there are nucleotides.\n";
	    last;
	}

	for (my $i = 0; $i < @nucleotides; $i++) {
	    $seqHash{$nucleotides[$i]} = $vals[$i];
	}

	$consensusSeq->add_element(%seqHash);
    }

    return $consensusSeq;
}

=head2 path

 Usage     : my $title = $SeqIO->path;
           : $SeqIO->path(path => "Some Path")
 Purpose   : Simple title for sequence to help differeniate between sequences eg, 5' vs 3'
 Returns   : The title of the consensus sequence
 Argument  : A new string title
 Throws    : 
 Comment   : 

See Also   : 

=cut

sub path {
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;
    
    $self->{-path} = $params{path} if defined ($params{path});
    return $self->{-path};
}

=head2 _read_seq
    Purpose     : Attempts to read the next sequence from fh
    Returns     : array reference that contains each line in a sequence.
    Args        : 
=cut

sub _read_seq {
    
    my $self = shift;
    
    my @seq;  # The Sequence
    my $pos;  # Keep track of position in file so we can seek back to it.
    
    return if !defined($self->{_fh});

    my $fh = $self->{_fh};
#
# Read until a ">" is read, a ">" indicates a new ConSequenceIO
#
#    carp "ConsSeqIO:: Reading until a > is found...";
    while (<$fh>) {	
#	carp "ConsSeqIO:: Next Line: '$_'";
	chomp;
	next if length($_) == 0;

	# Start with a header line?
	if ($_ =~ m/^>/) {
#	    carp "ConsSeqIO:: Found the next sequence!";
	    push @seq, $_;
	    last;
	}
    }
    
    return unless @seq;

    #
    # Read until the next ">" is found, while storing all consecutive lines.
    #
    # NOTE: cant use while(<fh>) because you need to use tell() before the read
    #
    while ( not eof($fh) ) {
	$pos = tell();
	$_ = <$fh>;
	chomp;

	# Skip blank lines
	next if length($_) == 0;

	# If we found another ConSequence, back up and return the @seq
	if ($_ =~ m/^>/) {
#	    carp "ConsSeqIO:: Found a secon sequence. Using seek...";
	    seek $self->{_fh}, $pos, 0 or die "Couldn't seek to $pos: $!\n";
	    last;
	}
	else {
	    push @seq, $_;
	}
	
    }
    
    return \@seq;    
}

=head2 write_seq
    Usage       : $consensusSeq->write_seq($seq)
    Purpose     : Writes a sequence 
    Returns     : 1 for success, 0 for error
    Args        : ConsensusSeq object
=cut

sub write_seq {
    my ($self, $seq) = @_;

    if (!$self->_fh) {
	if ($self->{clobber}) {
	    open $self->{_fh} , ">", $self->{path} or die "Error opening fh: $!\n";
	} 
	else {
	    open $self->{_fh} , ">>", $self->{path} or die "Error opening fh: $!\n";
	}
    }
    
    print {$self->{_fh}} $seq->to_string;
}

1;

