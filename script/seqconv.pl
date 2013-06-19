#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Bio::Phylo::Util::Logger;

# process command line arguments
my ( $deserializer, $serializer, $verbosity );
my $infile  = '-';
my $outfile = '-';
GetOptions(
	'infile=s'       => \$infile,
	'outfile=s'      => \$outfile,
	'serializer=s'   => \$serializer,
	'deserializer=s' => \$deserializer,
	'verbose'        => \$verbosity,
);

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new(
	'-class' => 'main',
	'-level' => $verbosity,
);

# create in handle
my $reader;
if ( $infile eq '-' ) {
	$reader = Bio::SeqIO->new(
		'-fh'     => \*STDIN,
		'-format' => $deserializer,
	);
}
else {
	$reader = Bio::SeqIO->new(
		'-file'   => $infile,
		'-format' => $deserializer,
	);
}

# create out handle
my $writer;
if ( $outfile eq '-' ) {
	$reader = Bio::SeqIO->new(
		'-fh'     => \*STDOUT,
		'-format' => $serializer,
	);
}
else {
	$reader = Bio::SeqIO->new(
		'-file'   => $outfile,
		'-format' => $serializer,
	);
}

# now iterate over the records
while( my $seq = $reader->next_seq ) {
	$writer->write_seq($seq);
}

