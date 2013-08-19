#!/usr/bin/perl
use strict;
use warnings;
use Bio::Phylo::IO 'parse';

my $project = parse(
	'-format' => 'nexus',
	'-file'   => shift(@ARGV),
	'-as_project' => 1,
);

my ($taxa)   = @{ $project->get_taxa };
my ($matrix) = @{ $project->get_matrices };
my ($forest) = @{ $project->get_forests };

my %keep = (
	'046'        => q['S. pimpinellifolium'],
	'051'        => q['S. chiemliewskii'],
	'054'        => q['S. cheesmaniae'],
	'073'        => q['S. pennellii'],
	'ITAG2'      => q['S. lycopersicum "Heinz"'],
	'U0015716'   => q['S. lycopersicum U0015716, 1891'],
	'U0015717'   => q['S. lycopersicum U0015717, 1830'],
	'WAG0463703' => q['S. peruvianum WAG0463703, 1958'],
);

for my $block ( $taxa, $matrix ) {
	my @things = grep { $keep{$_->get_name} } sort { $a->get_name cmp $b->get_name } @{ $block->get_entities };
	$block->clear;
	$block->insert($_->set_name($keep{$_->get_name})) for @things;
}

$project->delete($forest);
print $project->to_nexus;

