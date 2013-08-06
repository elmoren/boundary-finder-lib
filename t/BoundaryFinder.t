use strict;
use warnings;
use Test::More tests => 24;
use Carp;
use Data::Dumper;

BEGIN {
      use_ok('BlastManager');
      use_ok('BlastResults');
      use_ok('BoundaryFinder');
      use_ok('ConsensusSeq');
      use_ok('ConsensusSeqIO');
}

my $seq;
my $cons_seq;

#
# Consesnus Sequence
#
my $seqIO = ConsensusSeqIO->new(-path => "t/demo.cons");

isa_ok($seqIO, 'ConsensusSeqIO');
$cons_seq = $seqIO->next_seq;

ok ($cons_seq->length == 3, "Check Length");
ok ($cons_seq->score_seq(seq => "caG") == 0.25, "Score Sequence");
ok ($cons_seq->score_seq(seq => "aag") == 0.0, "Score Sequence");
ok ($cons_seq->missing_index_score == 0.1, "Missing Index Score");

# =============================================================================
# Test BoundaryFinder::score_seq with quick and simple sequence
# =============================================================================

$seq = "aaactgcaaggt";

my $hits = BoundaryFinder::score_seq(-seq => $seq, -cons_seq => $cons_seq);

ok (keys( %{ $hits } ) == 2, "Num Keys");
ok ($hits->{3} == 0.25, "Hit 1");
ok ($hits->{6} == 0.125, "Hit 2");

# =============================================================================
# Test score_hit with the full process from blast to BoundaryFinder::score_hit
# =============================================================================

#
# First run blasts
#

my $query = "t/blastTest/DMel_prot_a.faa";
my $db = "t/blastTest/dmel";
my $bm = BlastManager->new(-blast_db => $db);
$bm->evalue(-e => 1e-1);
my $br = $bm->tblastn(-query => $query);

#
# Load the 5' and 3' consensus sequences for Drosophila
#

$seqIO = new ConsensusSeqIO(-path => "t/drosophila.cons");
my $cons_seq_5;
my $cons_seq_3;

$cons_seq_5 = $seqIO->next_seq;
$cons_seq_3 = $seqIO->next_seq;

#
# Verify loading two sequences from the same file work
#
ok ($cons_seq_5->missing_index_score == 0, "5' Missing Index Score");
ok ($cons_seq_3->missing_index_score == 0, "3' Missing Index Score");
ok ($cons_seq_5->boundary ==  5, "5' Boundary Index");
ok ($cons_seq_3->boundary == 14, "3' Boundary Index");


#
# Not the best but...
# Testing by asking if its marginally better than the top score
# Because perl goes to more decimal places than my personal calculator
#
ok($cons_seq_5->score_seq(seq => "aaaaggtaagtat") > 0.00045436327, "Test max 5' score");
ok($cons_seq_3->score_seq(seq => "ttttttttttacaggtc") >  4.966e-6, "Test max 3' score");

#
# Go through each blast hit
#
#

# Verify size
ok($br->size == 8, "Assert BlastResults Size");

#
# Run the Boundary Finder score_hit function
#
#

my $hit = $br->get(-idx => 3);
ok($hit->{sseqid} eq "2R:4840277-4848161", "Testing Sub Seq ID After Load");

#carp Dumper($hit);

#
# The 5' and 3' are switched because that is the way the consensus sequences
# are from the literature
#
my $bf_results = BoundaryFinder::score_hit(-hit => $hit, 
				    -db => $db,
				    -cons_seq_5 => $cons_seq_3,
				    -cons_seq_3 => $cons_seq_5
    );

#
# Now you can process the results to determine which sequences you want. It
# is up to the user to determine how they want to examine these results.
#

# carp Dumper($bf_results);

#
# Get the upstream hit
#

# We know the hit is boundary results 40 and 314
# Use a hit + the consensus sequence boundary to get the actual exon

my $up_index = 318 + $cons_seq_3->boundary;
my $dn_index = 642 + $cons_seq_5->boundary - 1;

my $exon = $bm->get_sequence(-db => $db, 
			     -seq_id => $hit->{sseqid},
			     -start => $up_index,
			     -end => $dn_index
    );

#carp $exon;

my $expected = "cctcccccattgagttcgtaatggacaccagcctaaatggctctcgatctgatcctgcgactgcaacccaccctggcaaatggccaccgacgaccaaagcgccggcgctcagggctccagctggaacagctggacatgcttatcagtcgccatcgtcgtcgctggccgctgataacagaagccacgacaacaacaacgctagcgccgtgtccatgctgctgccgcaagacggcgacgcatctggagctgtggcccctgctgttactccccagctgcccatctacattgctcagcctagcgcgaaaaagccagaaa";
ok (lc($exon) eq $expected, "Check exon against original actuall exon");

#
# Quick test of pulling a boundary area
#

#
# Upstream, bad float comparison
#
$exon = $bm->get_sequence(-db => $db, 
			  -seq_id => $hit->{sseqid},
			  -start => 318,
			  -end => 318 + $cons_seq_3->length -1
    );

ok ($cons_seq_3->score_seq(seq => $exon) > 2.1336e-10, "Upstream pull and score");

#
# Downstream, yet another bad float comparison
#
$exon = $bm->get_sequence(-db => $db, 
			  -seq_id => $hit->{sseqid},
			  -start => 642,
			  -end => 642 + $cons_seq_5->length -1
    );

ok ($cons_seq_5->score_seq(seq => $exon) > 1.9012e-05, "Downstream pull and score");
