#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::IO qw'parse_tree parse_matrix';

# process command line arguments
my ( $newickfile, $phylipfile );
GetOptions(
	'tree=s'      => \$newickfile,
	'alignment=s' => \$phylipfile,
);

# read tree file
my $tree = parse_tree(
	'-format'     => 'newick',
	'-file'       => $newickfile,
	'-as_project' => 1,
);

# read matrix file
my $matrix = parse_matrix(
	'-format'     => 'phylip',
	'-type'       => 'dna',
	'-file'       => $phylipfile,
	'-as_project' => 1,
);

# pad matrix, if need be
my $nchar = $matrix->get_nchar;
if ( my $padding = $nchar % 3 ) {
	$matrix->visit(sub{
		my $row = shift;
		my $char = $row->get_char;
		$char .= '?' x ( $padding + 1 );
		$row->set_char($char);
	
	});
}

# do the object juggling to make a nexus file with referential integrity
my $fac = Bio::Phylo::Factory->new;
my $proj   = $fac->create_project;
my $forest = $fac->create_forest;
$forest->insert($tree);

my $ftaxa = $forest->make_taxa;
my $mtaxa = $matrix->make_taxa;
my $merged = $ftaxa->merge_by_name($mtaxa);

$forest->set_taxa($merged);
$matrix->set_taxa($merged);

$proj->insert($merged);
$proj->insert($matrix);
$proj->insert($forest);

print $proj->to_nexus(
	'-make_translate' => 0,
);