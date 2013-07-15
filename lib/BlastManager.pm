package BlastManager;

use warnings;
use strict;

# For determining blast+ tools availability
use Carp;
use File::Which;

#################### main pod documentation begin ###################

=head1 NAME

BlastManager - Blast Tools Wrapper

=head1 SYNOPSIS

  use BlastManager;

  my $BM = BlastManager->new(-blast-db => "path/to/blast/db");

  # FASTA formatted query file
  my $query_path = "path/to/query.fna";

  my $blast_results = $BM->blastn(-query => $query_path);
  ... or
  my $blast_results = $BM->tblastn(-query => $query_path);
  ... or
  my $blast_results = $BM->blastx(-query => $query_path);

=head1 DESCRIPTION

The BlastManager is a wrapper class for the blast command line tools. You specify the source
database, the blast type (blastn, tblastn, etc) and the query. It then gives you a BlastResults
object from the result. 

The BlastResults object is a wrapper to the outputs of a blast query, providing methods to get
sequences of the blast, the score, and the human readable output. See BlastResults.

The BlastManager was made specifically for the BoundaryFinder. It can be used for other purposes,
but it is not as flexible as teh command line tools.

=cut

#################### main pod documentation end ###################

=head1 METHODS

=cut

#################### subroutine header begin #######################

=head2 new

 Usage     : my $BM = BlastManager->new();

 Purpose   : Creates new BlastManager

 Returns   : A BlastManager

 Argument  : 

 Comment   : This class is a wrapper class for the blast command
             line interface.

=cut

sub new
{
    my $class = shift;
    
    my $self = {
	-blast_db => "",
	-evalue => 1,
	@_
    };

    bless ($self, ref ($class) || $class);

    # Test if Blast Exists on the system.
    if ( which('blastn') ) {
	print "Warning: Blast+ tools not found in PATH\n";
    }

    return $self;
}

=head2 blast_db

 Usage     : my $blast_db = $BM->blast_db(-db => "some/path/blastdb_name");
           : my $blast_db = $BM->blast_db();

 Purpose   : Getter/Setter for the current blastdb parms

 Returns   : The blastdb path, 

 Argument  : -db: database path and name as string

 Comment   : 

=cut

sub blast_db
{
    my $self = shift;
    
    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;
    
    $self->{-blast_db} = $params{-db} if defined ($params{-db});
    return $self->{-blast_db}; 
}

=head2 blastn

 Usage     : my $blast_results = $BM->blastn({-query => "queryPath"});

 Purpose   : Runs a blastn query with the given query path or sequence

 Returns   : A BlastResults object

 Argument  : -query: Path to a blast query Database
             -query-seq: a seqence of NT

 Comment   : If both arguments are provided, it uses the query path method and gives a warning

=cut

sub blastn
{
    
}

=head2 blastx

 Usage     : my $blast_results = $BM->blastx(-query => "queryPath");

 Purpose   : Runs a blastx query with the given query path or sequence

 Returns   : A BlastResults object

 Argument  : -query: Path to a blast query Database
             -query-seq: a seqence of NT

 Comment   : If both arguments are provided, it uses the query path method and gives a warning

=cut

sub blastx
{
    croak "BlastManager Error: Not supported yet";
}

=head2 evalue

 Usage     : my $evalue = $BM->evalue(-e => 10e-1);
             my $evalue = $BM->evalue;

 Purpose   : Gets and sets the evalue for the queries

 Returns   : A BlastResults object

 Argument  : -query: Path to a blast query Database
             -query-seq: a seqence of NT

 Comment   : The expected value for saving hits (evalue) is useful when want to remove
             hits that do not have enough matches. Or, it can be used to try and acquire
             more hits that have less matches. 

             See http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#expect

=cut

sub evalue
{
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    $self->{-evalue} = $params{-e} if defined ($params{-e});
    return $self->{-evalue};
}


=head2 get_sequence

 Usage     : my $sequence = $BM->get_sequence(-db => "dmel", -seq_id => "seq_id", -start => 1, -end => 10);

 Purpose   : gets a sequence from the database

 Returns   : A fasta sequence from the blast database

 Argument  : [-db] : The database to get it from (not required if -blast_db set)
             -seq_id : The sequence id of the database
             -start : The start index of the sequence
             -end : then end index

 Comment   :


=cut

sub get_sequence
{
    my $self = shift;

    croak "BlastManager->get_sequence(): Illegal parameter list has odd number of values" if @_ % 2;
    
    my %args = @_;

    for my $required (qw{ -seq_id -start -end }) {
	croak "Required parameter '$required' not passed to score_seq"
            unless exists $args{$required}; 
    }
    
    my ($seq, $db, $range, $seq_id);
    
    if (exists $args{-db}) {
	$db = $args{-db};
    } 
    elsif (exists $self->{-blast_db} && $self->{-blast_db} ne "")
    {
	$db = $self->{-blast_db};
    }
    else
    {
	croak "BlastManager->get_sequence(): No blast db to query";
    }

    # Regular order
    $seq_id = $args{-seq_id};
    if ($args{-start} <= $args{-end}) {
	
	$range = "$args{-start}-$args{-end}";
	$seq = `blastdbcmd -db $db -range $range -entry $seq_id`;

	my @lines = split /\n/, $seq;
        shift @lines;
	$seq = join "", @lines;
    }
    # Reversed Complement Order
    else {

	$range = "$args{-end}-$args{-start}";
	$seq = `blastdbcmd -db $db -range $range -entry $args{-seq_id}`;
	
	my @lines = split /\n/, $seq;
        shift @lines;
	$seq = join "", @lines;

	$seq = reverse($seq);    
	# complement the reversed DNA sequence
	$seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    }

    return $seq;
}

=head2 tblastn

 Usage     : my $blast_results = $BM->tblastn(-query => "queryPath");

 Purpose   : Runs a tblastn query with the given path to query

 Returns   : A BlastResults object, or croaks if error in query

 Argument  : -query: Path to a blast query file
             [-db] : Path to subject database. (Required if -blast_db not set)

 Comment   : Other args to come, maybe. Custom formats not supported, yet.

=cut

sub tblastn
{
    my $self = shift;

    croak "BlastManager->tblastn(): Illegal parameter list has odd number of values" if @_ % 2;
    
    my $query = "";
    my $db = "";
    my $parms .= "-evalue " . $self->evalue;
    my $fmt = ""; # Not currently supported
    my %args =	@_;
 
    for my $required (qw{ -query }) {
	croak "Required parameter '$required' not passed to score_seq"
            unless exists $args{$required}; 
    }

    $query = $args{-query};

    if (exists $args{-db}) {
	$db .= $args{-db};
    }
    elsif ($self->{-blast_db} ne "") {
	$db .= $self->{-blast_db};
    }
    else {
	croak "BlastManager->tblastn(): Error: No subject database";
    }

    my $csv   = `tblastn -query $query -db $db $parms -outfmt "10 $fmt"`;
    my $blast = `tblastn -query $query -db $db $parms`;

    # Check results, and if okay, Create BlastResults mand return it
    
    if ($csv =~ m/error/) {
	croak "BlastManager->tblastn(): Query error: \n $csv";
    }
    
    # Get the DB Source...?
    # my $source;

    # @TODO: Figure out how tot get the db source.
    # Currently assumes/requires the source to be $db.fas

    return BlastResults->new( -raw_csv => $csv, 
			      -raw_blast => $blast, 
			      -blast_type => 'tblastn',
			      -db_source => $db . '.fas',
			      -db_name => $db
	);

}

#################### subroutine header end ########################

#################### footer pod doc ###############################

=head1 BUGS

1. Database source file naming

Assumes the database fasta source file is the same as the database name with
a .fas extension. Not true because the user can name the database whatever
they want with makeblastdb. Looking for solution. It may have to force naming
convention.

=head1 AUTHOR

    Nathan Elmore
    nate@elmoren.com
    http://www.elmoren.com

=head1 COPYRIGHT

=cut

1;
