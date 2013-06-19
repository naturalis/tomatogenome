#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_matrix';

# process command line arguments
my ( $infile, $outfile, $skip );
GetOptions(
	'infile=s'  => \$infile,
	'outfile=s' => \$outfile,
	'skip=s'    => \$skip,
);

# read input matrix
my $matrix = parse_matrix(
	'-format' => 'fasta',
	'-type'   => 'dna',
	'-file'   => $infile,
	'-as_project' => 1,
);

# delete outgroup row
my $row = $matrix->get_by_name($skip);
$matrix->delete($row) if $row;

# convert n/N to ?
$matrix->visit(sub{
	my $row = shift;
	my $char = $row->get_char;
	$char =~ s/[nN]/?/g;
	$row->set_char($char);
});

# compute proportion of invariant sites
my $pinvar = $matrix->calc_prop_invar;

# write output
open my $fh, '>>', $outfile or die $!;
print $fh $infile, "\t", $pinvar, "\n";