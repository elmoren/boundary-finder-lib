package BoundaryFinder;

use warnings;
use strict;
use Carp;

BEGIN {
     require Exporter;
     use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
     our $VERSION     = '0.01';
     our @ISA         = qw(Exporter);

     our %EXPORT_TAGS = (
	 'all' => [ qw (
                         score_seq
                         score_hit
                  ) ]
     );
     our @EXPORT_OK   = qw( @{ $EXPORT_TAGS{'all'} } );
     our @EXPORT      = qw( score_seq );
}

#################### main pod documentation begin ###################

=head1 NAME

BoundaryFinder - Exon/Intron Boundary Finder

=head1 SYNOPSIS

  use BoundaryFinder;

  Coming soon...

=head1 DESCRIPTION

The boundary finder uses consensus sequences to find and score boundaries in a 
sequence. This was developed for finding exon boundaries, but it can be used 
with any consensus sequence.


=cut

#################### main pod documentation end ###################

=head1 METHODS

=cut

#################### subroutine header begin #######################

=head2 score_seq

 Usage     : my $scoresRef = BoundaryFinder::score_seq( 
                                               -seq => $seq, 
                                               -cons_seq => $cons_seq
                                              );

 Purpose   : Attempts to score a sequence hit
             BoundaryFinder::ConsensusSeq object.

 Argument  : -seq : A sequence to score
             -cons_seq : initialized BoundaryFinder::ConsensusSeq object
             [-min_score] : The minimum score to save (Defualt: 0)

 Returns   : A hash reference in which the key is the index of the start of
             the match, and the value is the score.

             $scores = {
               17 => 1.23e-08,
               50 => 7.23e-04,
               ...
             }

 Comment   : 

=cut

sub score_seq
{

    croak "Illegal parameter list has odd number of values" if @_ % 2;

    my %params = (
	-min_score => 0,
	@_
    );

    my %results;

    # Check required parameters
    for my $required (qw{ -seq -cons_seq }) {
	croak "Required parameter '$required' not passed to score_seq"
            unless exists $params{$required};
    }
    
    # Use the cons seq to score small chunks and keep track of the hits above 
    # the min_score
    my $cons = $params{-cons_seq};
    for (my $i = 0; $i < length($params{-seq}); $i++) {
	my $sub_seq = substr($params{-seq}, $i, $cons->length);
	my $score = $cons->score_seq(seq => $sub_seq);

	# Add if it meets the minimum score
	if ($score > $params{-min_score}) {
	    $results{$i} = $score;
	}
    }

    return \%results;
}

=head2 score_hit

 Usage     : my $ref_to_hit = $blastResultsObj->get(-idx => 0);
             my $scoresRef = BoundaryFinder::score_hit ( 
                                               -hit => $ref_to_hit, 
                                               -db => $db_name,
                                               -cons_seq_5 => $cons_seq_5_prime,
                                               -cons_seq_3 => $cons_seq_3_prime,
                                               -padding => $padding_len,
                                               -min_sore => $min_score
                                              );

 Purpose   : Searches for Boundaries given a BlastResults hit. It attempts to 
             pull the sequence from the database to score it. 

             It uses two consensus sequences, one for the 5' and 3' end. The 
             padding value is optionally set to adjust the 'n' NTs upstream and 
             downstream of 'sstart' and 'send'. Otherwise it searches the whole 
             blast hit, with a default 50 NT upstream and downstream.

             [sstart - 50] to [send + 50]

             The boundary finder will ignore scores less than the optional 
             '-min_score' argument.

 Argument  : -hit
             -db : the database name
              -cons_seq_5 : The 5' consensus sequence (ConsensusSeq object)
             -cons_seq_3 : The 3' consensus sequence (ConsensusSeq object)
             [-padding] : The number of nucleotides to pad the hit on both the
                          5' and 3' ends of the blast hit. Otherwise it defaults
                          to 50 nucleotides.
             [-min_score] : Minimum score to keep in the results

 Returns   : A hash reference in which the key is the index of the start of
             the match, and the value is the score.

             $scores = {
               5_prime_title => { 17 => 1.23e-08,
                                  50 => 7.23e-04,
                                   ... },
               3_prime_title => { 300 => 2.93e-10,
                                  327 => 6.03e-05,
                                   ... }
             }

 Comment   : This class is a wrapper class for the blast command
             line interface.

=cut

sub score_hit
{

    croak "Illegal parameter list has odd number of values" if @_ % 2;

    my %params = (
	-padding => 50,
	-min_score => 0,
	@_
    );
    
    #
    # Check the required parameters
    # 
    for my $required ( qw { -db -hit -cons_seq_5 -cons_seq_3 } ) {
	croak "Required parameter '$required' not passed to score_seq"
            unless exists $params{$required};
    }

    # Variables required for the Boundary Finder
    my ($seq_5_prime, $seq_3_prime, $start, $stop, $cons_seq);
    my $BM  = BlastManager->new(-db => $params{-db});
    my $hit = $params{-hit};
    my %results;
    
    #
    # Check the required values in the blast hit
    #
    for my $required ( qw { sstart send sseqid } ) {
	croak "Required parameter '$required' not passed to score_seq"
            unless exists $hit->{$required};
    }

    #
    # Aquire the 5' and 3 prime sequences
    #

    # 5' Sequence
    $start       = $hit->{sstart} - $params{-padding};
    $start       = 0 if $start < 0;
    $stop        = $hit->{send}; # + $params{-padding};
    $seq_5_prime = $BM->get_sequence(-db     => $params{-db},
				     -seq_id => $hit->{sseqid},
				     -start  => $start,
				     -end    => $stop
    );
    
    # 3' Sequence
    $start       = $hit->{sstart}; # - $params{-padding};
    $start       = 0 if $start < 0;
    $stop        = $hit->{send} + $params{-padding};
    $seq_3_prime = $BM->get_sequence(-db     => $params{-db},
				     -seq_id => $hit->{sseqid},
				     -start  => $start,
				     -end    => $stop
    );
    
    
    $cons_seq = $params{-cons_seq_5};
    $results{$cons_seq->title} = score_seq(-seq       => $seq_5_prime, 
					   -cons_seq  => $cons_seq, 
					   -min_score => $params{-min_score}
    );

    $cons_seq = $params{-cons_seq_3};
    $results{$cons_seq->title} = score_seq(-seq       => $seq_3_prime, 
					   -cons_seq  => $cons_seq, 
					   -min_score => $params{-min_score}
    );

    return \%results;
}

#################### footer pod doc ###############################

=head1 BUGS


=head1 AUTHOR

    Nathan Elmore
    nate@elmoren.com
    http://www.elmoren.com

=head1 COPYRIGHT

=head1 SEE ALSO

=cut

1;
