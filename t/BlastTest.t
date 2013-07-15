use Test::More tests => 9;
use Carp;

BEGIN {
    use_ok('BlastManager');
    use_ok('BlastResults');
}     

my $query = "t/blastTest/DMel_prot_a.faa";
my $db = "t/blastTest/dmel";

my $bm = BlastManager->new(-blast_db => $db);
my $br = $bm->tblastn(-query => $query);

isa_ok($br, 'BlastResults');
ok($br->size == 9, "Check Size");

$bm->evalue(-e => 1e-1);
$br = $bm->tblastn(-query => $query);
ok($br->size == 8, "Check Size");

$br->save(-file => "t/blastTest/dmel.csv");

my $seq = $bm->get_sequence(-seq_id => "2R:4840277-4848161", -start => 1, -end => 10);

$seq = lc ($seq);
ok($seq eq "gccgctcctc", "Test Get Sequence");

$seq = $bm->get_sequence(-seq_id => "2R:4840277-4848161", -start => 10, -end => 1);
$seq = lc ($seq);
ok($seq eq "gaggagcggc", "Test Reverse Complement Sequence");

$br->load_file(-file => "t/blastTest/dmel.csv");
ok($br->size == 8, "Check Size After Load");
ok($br->get(-idx => 0)->{sseqid} eq "2R:4840277-4848161", "Testing Sub Seq ID After Load");
