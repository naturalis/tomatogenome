#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $in  = \*STDIN;
my $out = \*STDOUT;

my $reader = Bio::SeqIO->new(
	'-format' => 'genbank',
	'-fh'     => $in,	
);

my $writer = Bio::SeqIO->new(
	'-format' => 'fasta',
	'-fh'     => $out,
);

while ( my $seq = $reader->next_seq ) {
	$writer->write_seq( $seq );
}
