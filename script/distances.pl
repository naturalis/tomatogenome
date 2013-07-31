#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util 'looks_like_number';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':simple';

sub check_args {

	# process command line arguments
	my $verbosity = WARN;
	my ( $infile, $from, $to );
	GetOptions(
		'verbose+' => \$verbosity,
		'infile=s' => \$infile,
		'from=s'   => \$from,
		'to=s'     => \$to,
	);
	VERBOSE( '-level' => $verbosity, '-class' => 'main' );
	INFO "infile is $infile";
	INFO "from is $from";
	INFO "to is $to";
	my @from = grep { /\S/ } split /,/, $from;
	my @to   = grep { /\S/ } split /,/, $to;
	return 'infile' => $infile, 'from' => \@from, 'to' => \@to;
}

sub main {
	my %args = check_args();
	
	# parse the input tree
	my $tree = parse_tree(
		'-file'       => $args{'infile'},
		'-format'     => 'newick',
		'-as_project' => 1,
	);
	
	# print the table header
	print "GENE\tFROM\t", join( "\t", map { "TO:$_" } @{ $args{'to'} } ), "\n";
	my $gene = $args{'infile'};
	$gene =~ s/.+\/(.+)\..+$/$1/;
	
	# iterate over input names
	for my $from_name ( @{ $args{'from'} } ) {
		my @result = ( $gene, $from_name );
		my $from = $tree->get_by_name( $from_name );
		
		# the "to" taxon might not exist in all trees!
		for my $to_name ( @{ $args{'to'} } ) {
			my $to = $tree->get_by_name( $to_name );
			if ( $to ) {
				push @result, sprintf('%.9f', $from->calc_patristic_distance($to) );
			}
			else {
				push @result, 'n/a';
			}
		}
		my ($nearest) = sort { $a <=> $b } grep { looks_like_number $_ } @result;
		print join( "\t", map { $_ eq $nearest ? "*$_" : $_ } @result ), "\n";
	}
}
main();