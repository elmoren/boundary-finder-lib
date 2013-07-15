# -*- perl -*-

# Test ConsensusSeq

use Test::More tests => 7;



BEGIN { 
    use_ok('ConsensusSeq');
    use_ok('ConsensusSeqIO');
}

my $seqIO = ConsensusSeqIO->new(-path => "t/demo.cons");

isa_ok ($seqIO, 'ConsensusSeqIO');

my $seq = $seqIO->next_seq;

isa_ok ($seq, 'ConsensusSeq');

ok ($seq->length == 3, "Check Length");
ok ($seq->score_seq(seq => "caG") == 0.25, "Score Sequence");
ok ($seq->score_seq(seq => "aag") == 0.0, "Score Sequence");
