# -*- perl -*-

# t/001_load.t - check module loading

use Test::More tests => 5;

BEGIN { 

    use_ok( 'BlastManager' ); 
    use_ok( 'BlastResults' ); 
    use_ok( 'BoundaryFinder' ); 
    use_ok( 'ConsensusSeq' ); 
    use_ok( 'ConsensusSeqIO' ); 

}
