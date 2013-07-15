# -*- perl -*-

# Test ConsensusSeq

use Test::More tests => 7;

BEGIN { 
    use_ok( 'ConsensusSeq' ); 
}

my $consensus = ConsensusSeq->new();

isa_ok ($consensus, 'ConsensusSeq');

# add a few elements
$consensus->add_element('a' => 0.5, 't' => 0.5, 'c' => 0.0, 'g' => 0.0);
$consensus->add_element('a' => 0.0, 't' => 0.0, 'c' => 0.5, 'g' => 0.5);

# check the scores
ok ($consensus->score_element(nucl => 'a', index => 0) == 0.5, "Score Element");
ok ($consensus->score_element(nucl => 'c', index =>0) == 0.0, "Score Element");
ok ($consensus->score_element(nucl => 'c', index => 1) == 0.5, "Score Element");

# Check a sequence
ok ($consensus->score_seq(seq => "ac") == 0.25, "Score Sequence");

# verify length is correct
ok ($consensus->length == 2, "Check Length");

print "Testing to_string:\n";
print $consensus->to_string, "\n ...done \n";

