package BlastResults;

use warnings;
use strict;

use Carp;
use File::Basename;
use Data::Dumper;
#################### main pod documentation begin ###################

=head1 NAME

BlastResults

=head1 SYNOPSIS

my $results = BlastResults->new();
$results->load_file(-file => "results.csv");

for (my $i = 0; $i < $results->size; $i++) {

    my $hit = get(-idx => $i);
    print $hit->{bitscore};

}

=head1 DESCRIPTION

The BlastResults object is a wrapper to the outputs of a blast query, providing 
methods to get data from the blast, the score, and the human readable output.
It also has methods to save and load previous blast results.

The BlastResults object contains a list of blast hits returned from the blast 
query. Iterate through the the results with and index and size.

The user may remove hits and save the results if they wish, but need to be careful if
they edit the data in the hits.

=head2 SAVE/LOAD RESULTS

The BlastResults object provides methods for saving and loading the results with
save_results and load_results. The user can remove unwanted hits before saving 
with the remove method. When the user calls "save". The save method expects a 
name and a directory.

 The BlastResults writes:

"<directory>/<name>.csv." - CSV file containing the output and other meta data

=cut

#################### main pod documentation end ###################

=head1 METHODS

=cut

#################### subroutine header begin #######################

#
# Defualt blast format used for parsing csv files when formats not specified.
#
my $BLAST_FMT =  [ 'qseqid', 'sseqid', 'pident', 'length', 'mismatch',
          	   'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 
                   'bitscore' ];



=head2 db_source

 Usage     : my $src = $BR->db_source();
             my $src = $BM->db_source(-file => "path/to/fasta");

 Purpose   : Gets and sets the blast source file. This is the path the FASTA 
             file that is used to create the subject database. 
             
             Setting this is not recommended unless you are creating a 
             BlastResults object without using a BlastManager instances.
             
 Returns   : Path to a Blast database source.

 Argument  : [-file] : Path to set

 Comment   : This is stored so that the Boundary Finder (or user) can extract a
             sequence from the database. It is a path to a FASTA file

=cut

sub db_source
{
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    $self->{-db_source} = $params{-file} if defined ( $params{-file} );
    return $self->{-db_source};
}

=head2 blast_type

 Usage     : my $source = $BR->blast_type();
             my $evalue = $BR->blast_type(-type => "blastn");
 Purpose   : Gets and sets the blast type
             
 Returns   : Path to a Blast database source.

 Argument  : [-type] : Path to set

 Comment   : Settign this field not recommended if you use BlastManager to get
             the BlastResults object.

=cut

sub blast_type
{
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    $self->{-blast_type} = $params{-type} if defined ( $params{-type} );
    return $self->{-blast_type};
}

=head2 csv_file

 Usage     : my $file = $BR->csv_file();
             my $file = $BR->csv_file(-file => "some/path/csv");

 Purpose   : Gets and sets the blast type             

 Returns   : Path to the a BlastResults CSV

 Argument  : [-file] : Path to set

 Comment   : Settign this field not recommended if you use BlastManager to get
             the BlastResults object.

=cut

sub csv_file
{
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    $self->{-csv_file} = $params{-file} if defined ($params{-file});
    return $self->{-csv_file};
}

=head2 db_name

 Usage     : my $source = $BR->db_name();
             my $evalue = $BR->db_name(-name => "dmel");

 Purpose   : Gets and sets the blast database name

 Returns   : Path to a Blast database name

 Argument  : [-type] : Type to set

 Comment   : This is the name of the blast database. Which is not necessarily 
             the same as the database source name. You can have a fasta file,
             dmel.fas and name the database 'drosophila'. You can tell the
             name of the blase database by the name of the .nsq, .nhr, and .nin

=cut

sub db_name
{
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    $self->{-db_name} = $params{-name} if defined ($params{-name});
    return $self->{-db_name};
}

=head2 get

 Usage     : my $hit = get(-idx => $i);

 Purpose   : Gets a hit at -idx

 Returns   : Hash containing the blast hit data

 Argument  : -idx : index of hit

 Comment   : 

=cut

sub get
{
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    croak "Required parameter '-idx' not passed to score_seq"
	if not exists $params{-idx};

    if ($params{-idx} < 0 || $params{-idx} >= $self->size) {
	croak "BlastResults->get(): Parameter [-idx] is out of bounds";
    }
    
    return $self->{_hits}->[$params{-idx}];
}


=head2 new

 Usage     : my $br = BlastManager::BlastResults->new(-csv_file => "path/to/results/csv");

 Purpose   : Creates a new blast results wrapper to access the data.

 Returns   : A BlastManager::BlastResults Object

 Argument  : [-raw_csv] : The output of a blast with '-outfmt 10'
           : [-format] : Reference to and array of outfmt csv format. Or it uses 
                         default
           : [-raw_blast] : The raw blast output (Default -outfmt). Available 
                            when using BlastManager
           : [-blast_type] : The blast type type used (blastx, tblastn, etc)
           : [-db_source] : Path to fasta file. Cannot extract sequences if this
                            is not set.
           : [-db_name] : Name of the blast database

 Comment   : There are multiple ways to instantiate this object:
           1. Use the blast manager to run a blast, and it will return a BlastResults 
              object
           2. Create a new BlastResults object and pass it raw csv data by adding 
              -outfmt 10 to your query. Use -raw_csv parameter or loat_text()
           3. Create a new BlastResults object and calling the load_file function.

=cut

sub new 
{
    my $class = shift;

    croak "BlastResults->new(): Illegal parameter list has odd number of values" if @_ % 2;

    my %params = @_;
    
    # default params
    my $self = {
	-raw_csv => "",
	-csv_file => "",
	-format => undef,
	-raw_blast => "",
	-blast_type => "",
	-db_source => "",
	-db_name  => ""
    };
    
    # Validate params by checking if they exist is self
    for my $attr (keys %params) {	
        croak "Invalid parameter '$attr' passed to '$class' constructor"
            unless defined ($self->{$attr});

	$self->{$attr} = $params{$attr};
    }


    # Init attributes user is not allowed to specify
    $self->{_size} = 0;
    $self->{_hits} = ();
    
    #
    # Init object if raw csv, or csv_file provided. CANNOT DO BOTH
    #    
    if ( $self->{-raw_csv} ne "" && $self->{-csv_file} ne "" ) {
	croak "BlastResults: does not support initialization of -raw_csv and -csv_file simultaneously" 
    }

    bless ($self, ref ($class) || $class);
    #
    # Init the blast data for the object
    #
    if ( $self->{-raw_csv} ne "" ) {
	$self->load_text( -text => $self->{-raw_csv} );
    }
    elsif ( $self->{-csv_file} ne "" ) {
	$self->load_file( -file => $self->{-csv_file} );
    }
    else {
	
    }


    return $self;
}


=head2 load_text

 Usage     : $blast_results->load_text(-text => $str);

 Purpose   : Loads resutls from raw csv text

 Returns   : Nothing

 Argument  : -text : CSV String
             [-format] : format (if not defualt)

 Comment   : 

=cut

sub load_text
{
    my $self = shift;
    
    croak "BlastResults->load_text(): Illegal parameter list has odd number of values" if @_ % 2;
    
    my $args = {
	-format => undef,
	@_
    };

    for my $required (qw{ -text }) {
	croak "Required parameter '$required' not passed to score_seq"
            unless exists $args->{$required}; 
    }
    
    #
    # First clear out old data.
    #
    $self->{_size} = 0;
    $self->{_hits} = ();
    $self->{_file} = 0; # Not set by default
    
    if ( $args->{-format}) {
	$self->{-format} = $args->{-format};
    } 
    else {
	$self->{-format} = $BLAST_FMT;
    }

    #
    # Parse The CSV Data
    #

    # carp Dumper($self->{-format});

    foreach my $line (split (/\n/, $args->{-text})) {
	chomp $line;
	my $i = 0;
	my %hit = map { $self->{-format}->[$i++] => $_ } split (",", $line);
	
	push @{ $self->{_hits} }, \%hit;
	$self->{_size}++;
    }
    
    # carp Dumper($self->{_hits});

}

=head2 load_file

 Usage     : $blast_results->load_file(-file => $file);

 Purpose   : Loads blast results from a file

 Returns   : Nothing

 Argument  : [-file] : The path to csv file (If its not specified, will try
                       $self->{-csv_file} before giving an angry message

 Comment   : Checks the first line of the file for the format, if not present,
             uses the defualt of blasts 'outfmt 10' option

=cut

sub load_file
{
    my $self = shift;
    
    croak "BlastResults->load_file(): Illegal parameter list has odd number of values" if @_ % 2;
    
    my $args = {
	@_
    };
    
    # Check Optional File argument
    if (exists $args->{-file}) {
	$self->csv_file( -file => $args->{-file} );
    } 

    #
    # Look for .csv (required)
    #

    open CSV, "<", $self->csv_file or 
	croak "BlastResults::load_results Error: File not found: $!";

    #
    # Check for headers
    #
    
    my $line;
    my $csv_data = "";
    my $format;

    foreach my $line (<CSV>) {
	chomp $line;
	if ($line =~ m/DB_SOURCE=/) {
	    $self->{-db_source} = (split "=", $line)[1];
	}
	elsif ($line =~ m/DB_NAME=/) {
	    $self->{-db_name} = (split "=", $line)[1];
	}
	elsif ($line =~ m/FORMAT=/) {
	    $format = [ split /,/, $line ];
	}
	elsif ($line =~ m/BLAST_TYPE=/) {
	    $self->{-blast_type} = (split "=", $line)[1];
	}
	else {
	    $csv_data .= $line . "\n";
	}
    }


    close CSV;

    # Use defualt format if one was not specified in the file.
    $format = $BLAST_FMT if not $format;

    $self->load_text(-text => $csv_data, -format => $format);
}

=head2 raw_blast

 Usage     : my $file = $BR->raw_blast;
             my $file = $BR->raw_blast(-data => $data);

 Purpose   : Gets and sets the raw blast output

 Returns   : Raw Blast as string

 Argument  : [-data] : Raw blast string

 Comment   : 

=cut

sub raw_blast
{
    my $self = shift;

    croak "Illegal parameter list has odd number of values" if @_ % 2;
    my %params = @_;

    $self->{-raw_blast} = $params{-data} if defined ($params{-data});
    return $self->{-raw_blast};
}

=head2 save

 Usage     : $blast_results->save(-file => $file);
 Purpose   : Saves the results for future reference.
 Returns   : Nothing
 Argument  : -file : Desired file to write to. If you want to rewrite the file,
                     use $self->csv_file. Clobbers the old file, if it exists.
 Comment   : 

=cut

sub save
{

    my $self = shift;
    
    croak "BlastResults->save(): Illegal parameter list has odd number of values" if @_ % 2;
    
    my $args = {
	@_
    };

    # See if there is a file to save to
    $self->{_file} = $args->{-file} if exists $args->{-file};
    
    # open the file or croak
    if (exists $self->{_file}) {
	open OUT, ">$self->{_file}" or croak "BlastResults->save_results(): $!";
    }
    else {
	croak "BlastResults->save_results(): no file to write to.\n";
    }

    #
    # Write the headers
    #
 
    print OUT "DB_SOURCE=" . $self->db_source . "\n";
    print OUT "DB_NAME=" . $self->db_name . "\n";
    print OUT "BLAST_TYPE=" . $self->blast_type . "\n";
    
    my $fmt_string = join(",", @{ $self->{-format} }) if exists $self->{-format};
    print OUT "FORMAT=" . $fmt_string . "\n";
    
    #
    # Write the data
    #
    for (my $i = 0; $i < $self->size; $i++) {
	my @hit;
	foreach my $key ( @{ $self->{-format} } ) {
	    push @hit, $self->{_hits}->[$i]->{$key};
	}
	print OUT join(",", @hit) . "\n";
    }
    close OUT;
}

=head2 size

 Usage     : $num_result = $blast_results->size
 Purpose   : Gets number of blast hits
 Returns   : The total number of blast hits contained in the BlastResults
 Argument  : None

 Comment   : 

=cut

sub size
{
    my $self = shift;
    return $self->{_size};
}

#################### subroutine header end ########################

#################### footer pod doc ###############################

=head1 BUGS


=head1 AUTHOR

    Nathan Elmore
    nate@elmoren.com
    http://www.elmoren.com

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

Blast Command Line Tools, BlastResults

=cut

1;

# The preceding line will help the module return a true value
